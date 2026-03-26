import express from 'express';
import cors from 'cors';
import path from 'path';
import fs from 'fs';
import { fileURLToPath } from 'url';
import { Client } from '@notionhq/client';
import Anthropic from '@anthropic-ai/sdk';
import dotenv from 'dotenv';

dotenv.config();

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

const app = express();
app.use(cors());
app.use(express.json());

// Serve built frontend files
const distPath = path.join(__dirname, '..', 'dist');
app.use(express.static(distPath));

const notion = new Client({ auth: process.env.NOTION_API_KEY });
const anthropic = new Anthropic({ apiKey: process.env.ANTHROPIC_API_KEY });

const DATABASE_ID = process.env.NOTION_DATABASE_ID || '7b5d9b8c069c4f3bbe40b25f55afdf96';

// Local JSON cache file
const CACHE_FILE = path.join(__dirname, 'cache.json');

// In-memory cache
let cachedPages = [];
let lastFetchTime = 0;
let cachedDataSourceId = null;

// Auto-sync: poll for new/updated pages every 2 minutes
const SYNC_INTERVAL = 2 * 60 * 1000;
let lastSyncEditTime = null;

// --- Local file cache ---
function loadCacheFromDisk() {
  try {
    if (fs.existsSync(CACHE_FILE)) {
      const data = JSON.parse(fs.readFileSync(CACHE_FILE, 'utf-8'));
      cachedPages = data.pages || [];
      lastFetchTime = data.lastFetchTime || 0;
      lastSyncEditTime = data.lastSyncEditTime || null;
      console.log(`[Cache] Loaded ${cachedPages.length} pages from disk cache.`);
      return true;
    }
  } catch (e) {
    console.log('[Cache] Failed to load disk cache:', e.message);
  }
  return false;
}

function saveCacheToDisk() {
  try {
    fs.writeFileSync(CACHE_FILE, JSON.stringify({
      pages: cachedPages,
      lastFetchTime,
      lastSyncEditTime,
      savedAt: new Date().toISOString(),
    }));
    console.log(`[Cache] Saved ${cachedPages.length} pages to disk.`);
  } catch (e) {
    console.log('[Cache] Failed to save disk cache:', e.message);
  }
}

// --- Notion API ---
async function getPageMarkdown(pageId) {
  try {
    const response = await notion.pages.retrieveMarkdown({ page_id: pageId });
    return response.markdown || '';
  } catch (e) {
    return '';
  }
}

function extractPageProperties(page) {
  const props = page.properties;
  const title = props['Name']?.title?.map(t => t.plain_text).join('') || '';
  const description = props['설명']?.rich_text?.map(t => t.plain_text).join('') || '';
  const category = props['논문/미팅 분류']?.select?.name || '';
  const importance = props['중요도']?.select?.name || '';
  const link = props['링크 페이지']?.rich_text?.map(t => t.plain_text).join('') || '';
  const lastEdited = props['생성 일자']?.last_edited_time || '';
  const createdTime = props['첫 생성']?.created_time || '';

  return { title, description, category, importance, link, lastEdited, createdTime };
}

async function getDataSourceId() {
  if (cachedDataSourceId) return cachedDataSourceId;
  const db = await notion.databases.retrieve({ database_id: DATABASE_ID });
  if (db.data_sources && db.data_sources.length > 0) {
    const ds = db.data_sources[0];
    cachedDataSourceId = ds.data_source_id || ds.id;
  } else {
    cachedDataSourceId = db.data_source_id || db.id || DATABASE_ID;
  }
  return cachedDataSourceId;
}

async function fetchPageList() {
  const dataSourceId = await getDataSourceId();
  const pages = [];
  let cursor;
  do {
    const response = await notion.dataSources.query({
      data_source_id: dataSourceId,
      start_cursor: cursor,
      page_size: 100,
    });
    pages.push(...response.results);
    cursor = response.has_more ? response.next_cursor : undefined;
  } while (cursor);
  return pages;
}

async function fetchAllPages(force = false) {
  if (!force && cachedPages.length > 0) {
    return cachedPages;
  }

  console.log('[Fetch] Full fetch from Notion...');
  const pages = await fetchPageList();
  console.log(`[Fetch] Found ${pages.length} pages. Fetching content...`);

  const enrichedPages = [];
  let count = 0;
  for (const page of pages) {
    const props = extractPageProperties(page);
    let content = '';
    try {
      content = await getPageMarkdown(page.id);
    } catch (e) {
      // skip
    }

    enrichedPages.push({
      id: page.id,
      url: page.url,
      ...props,
      content: content.slice(0, 3000),
    });

    count++;
    if (count % 50 === 0) console.log(`[Fetch] Progress: ${count}/${pages.length}`);

    // Rate limit
    await new Promise(r => setTimeout(r, 50));
  }

  cachedPages = enrichedPages;
  lastFetchTime = Date.now();
  console.log(`[Fetch] Done! ${enrichedPages.length} pages loaded.`);

  saveCacheToDisk();
  return enrichedPages;
}

// Incremental sync: only fetch new/updated pages
async function incrementalSync() {
  try {
    const dataSourceId = await getDataSourceId();
    const response = await notion.dataSources.query({
      data_source_id: dataSourceId,
      page_size: 100,
      sorts: [{ property: '생성 일자', direction: 'descending' }],
    });

    if (response.results.length === 0) return;

    const latestEditTime = response.results[0].last_edited_time;

    if (lastSyncEditTime && latestEditTime === lastSyncEditTime) {
      return; // No changes
    }

    console.log('[Sync] Changes detected! Updating...');

    // Find pages edited after last sync
    let updatedCount = 0;
    for (const page of response.results) {
      const editTime = page.last_edited_time;
      if (lastSyncEditTime && editTime <= lastSyncEditTime) break;

      const props = extractPageProperties(page);
      let content = '';
      try {
        content = await getPageMarkdown(page.id);
      } catch (e) {
        // skip
      }

      const enriched = {
        id: page.id,
        url: page.url,
        ...props,
        content: content.slice(0, 3000),
      };

      // Update or add
      const existingIdx = cachedPages.findIndex(p => p.id === page.id);
      if (existingIdx >= 0) {
        cachedPages[existingIdx] = enriched;
      } else {
        cachedPages.push(enriched);
      }
      updatedCount++;
    }

    lastSyncEditTime = latestEditTime;
    lastFetchTime = Date.now();
    saveCacheToDisk();
    console.log(`[Sync] Updated ${updatedCount} pages. Total: ${cachedPages.length}`);
  } catch (e) {
    console.error('[Sync] Error:', e.message);
  }
}

// --- Search ---
function searchPages(pages, query) {
  const queryLower = query.toLowerCase();
  const queryTerms = queryLower.split(/\s+/).filter(t => t.length > 1);

  const scored = pages.map(page => {
    const searchText = `${page.title} ${page.description} ${page.category} ${page.content}`.toLowerCase();
    let score = 0;

    if (searchText.includes(queryLower)) score += 10;

    for (const term of queryTerms) {
      const regex = new RegExp(term, 'gi');
      const matches = searchText.match(regex);
      if (matches) score += matches.length;
      if (page.title.toLowerCase().includes(term)) score += 5;
      if (page.description.toLowerCase().includes(term)) score += 3;
    }

    return { ...page, score };
  });

  return scored
    .filter(p => p.score > 0)
    .sort((a, b) => b.score - a.score)
    .slice(0, 8);
}

// --- API Routes ---
app.post('/api/chat', async (req, res) => {
  try {
    const { message, history = [] } = req.body;
    if (!message) return res.status(400).json({ error: 'Message is required' });

    if (cachedPages.length === 0) {
      return res.status(503).json({ error: 'Server is still loading data. Please wait a moment.' });
    }

    const relevantPages = searchPages(cachedPages, message);

    const context = relevantPages.map((page, i) => {
      return `[Page ${i + 1}] Title: ${page.title}
Category: ${page.category}
Importance: ${page.importance}
Description: ${page.description}
Content: ${page.content.slice(0, 1500)}
URL: ${page.url}`;
    }).join('\n\n---\n\n');

    const messages = [];
    for (const h of history.slice(-6)) {
      messages.push({ role: h.role, content: h.content });
    }
    messages.push({ role: 'user', content: message });

    const systemPrompt = `You are a knowledgeable research assistant that answers questions based on a Notion database containing notes about:
- All-solid-state batteries (ASSB)
- NCM cathode materials
- Solid electrolytes (LPSCl, sulfide-based)
- DEM (Discrete Element Method) simulation
- Electrode fabrication and characterization
- Meeting notes, seminar content, and paper summaries

IMPORTANT RULES:
1. Answer in the same language as the user's question (Korean or English).
2. Base your answers ONLY on the provided database content below.
3. If the database doesn't contain relevant information, say so honestly.
4. Always cite which pages you referenced using their titles.
5. Be concise but thorough.

--- DATABASE CONTENT ---
${context || 'No relevant pages found for this query.'}
--- END DATABASE CONTENT ---`;

    const response = await anthropic.messages.create({
      model: 'claude-sonnet-4-20250514',
      max_tokens: 2000,
      system: systemPrompt,
      messages,
    });

    const answer = response.content[0].text;

    const references = relevantPages.map(p => ({
      id: p.id,
      title: p.title,
      url: p.url,
      category: p.category,
      importance: p.importance,
      description: p.description,
    }));

    res.json({ answer, references });
  } catch (error) {
    console.error('Chat error:', error);
    res.status(500).json({ error: error.message || 'Internal server error' });
  }
});

app.get('/api/pages', (req, res) => {
  const summaries = cachedPages.map(p => ({
    id: p.id,
    title: p.title,
    url: p.url,
    category: p.category,
    importance: p.importance,
    description: p.description,
    lastEdited: p.lastEdited,
  }));
  res.json({ pages: summaries, total: summaries.length });
});

app.get('/api/search', (req, res) => {
  const { q } = req.query;
  if (!q) return res.json({ results: [] });
  const results = searchPages(cachedPages, q);
  res.json({
    results: results.map(p => ({
      id: p.id,
      title: p.title,
      url: p.url,
      category: p.category,
      importance: p.importance,
      description: p.description,
      score: p.score,
    })),
  });
});

app.post('/api/refresh', async (req, res) => {
  try {
    const pages = await fetchAllPages(true);
    res.json({ message: `Refreshed ${pages.length} pages` });
  } catch (error) {
    res.status(500).json({ error: error.message });
  }
});

app.get('/api/status', (req, res) => {
  res.json({
    totalPages: cachedPages.length,
    lastFetchTime: lastFetchTime ? new Date(lastFetchTime).toISOString() : null,
    cacheFileExists: fs.existsSync(CACHE_FILE),
    autoSyncEnabled: true,
    syncIntervalMs: SYNC_INTERVAL,
  });
});

// SPA fallback
app.use((req, res, next) => {
  if (req.path.startsWith('/api')) return next();
  res.sendFile(path.join(distPath, 'index.html'));
});

const PORT = process.env.PORT || 3001;
app.listen(PORT, '0.0.0.0', () => {
  console.log(`Server running on http://localhost:${PORT}`);

  // 1. Try loading from disk cache (instant)
  const hasCache = loadCacheFromDisk();

  if (hasCache) {
    // Start incremental sync immediately
    console.log('[Startup] Using disk cache. Starting incremental sync...');
    incrementalSync().catch(console.error);
  } else {
    // First time: full fetch
    console.log('[Startup] No cache found. Doing full fetch...');
    fetchAllPages().catch(console.error);
  }

  // Start periodic sync
  setInterval(incrementalSync, SYNC_INTERVAL);
  console.log(`[Sync] Auto-sync every ${SYNC_INTERVAL / 1000}s`);
});
