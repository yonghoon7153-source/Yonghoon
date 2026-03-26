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

const distPath = path.join(__dirname, '..', 'dist');
app.use(express.static(distPath));

const notion = new Client({ auth: process.env.NOTION_API_KEY });
const anthropic = new Anthropic({ apiKey: process.env.ANTHROPIC_API_KEY });

const DATABASE_ID = process.env.NOTION_DATABASE_ID || '7b5d9b8c069c4f3bbe40b25f55afdf96';
const CACHE_FILE = path.join(__dirname, 'cache.json');

let cachedPages = [];
let lastFetchTime = 0;
let cachedDataSourceId = null;
const SYNC_INTERVAL = 2 * 60 * 1000;
let lastSyncEditTime = null;

function loadCacheFromDisk() {
  try {
    if (fs.existsSync(CACHE_FILE)) {
      const data = JSON.parse(fs.readFileSync(CACHE_FILE, 'utf-8'));
      cachedPages = data.pages || [];
      lastFetchTime = data.lastFetchTime || 0;
      lastSyncEditTime = data.lastSyncEditTime || null;
      console.log('[Cache] Loaded ' + cachedPages.length + ' pages from disk cache.');
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
      pages: cachedPages, lastFetchTime, lastSyncEditTime,
      savedAt: new Date().toISOString(),
    }));
    console.log('[Cache] Saved ' + cachedPages.length + ' pages to disk.');
  } catch (e) {
    console.log('[Cache] Failed to save disk cache:', e.message);
  }
}

async function getPageContent(pageId) {
  try {
    const blocks = [];
    let cursor;
    do {
      const response = await notion.blocks.children.list({ block_id: pageId, start_cursor: cursor, page_size: 100 });
      blocks.push(...response.results);
      cursor = response.has_more ? response.next_cursor : undefined;
    } while (cursor);
    const textParts = [];
    let firstImage = '';
    for (const block of blocks) {
      const type = block.type;
      const content = block[type];
      if (!content) continue;
      if (content.rich_text) { const text = content.rich_text.map(t => t.plain_text).join(''); if (text) textParts.push(text); }
      if (content.caption) { const caption = content.caption.map(t => t.plain_text).join(''); if (caption) textParts.push(caption); }
      if (!firstImage && type === 'image') { firstImage = content.file?.url || content.external?.url || ''; }
    }
    return { text: textParts.join('\\n'), firstImage };
  } catch (e) { return { text: '', firstImage: '' }; }
}

function extractPageProperties(page) {
  const props = page.properties;
  const title = props['Name']?.title?.map(t => t.plain_text).join('') || '';
  const description = props['\uc124\uba85']?.rich_text?.map(t => t.plain_text).join('') || '';
  const category = props['\ub17c\ubb38/\ubbf8\ud305 \ubd84\ub958']?.select?.name || '';
  const importance = props['\uc911\uc694\ub3c4']?.select?.name || '';
  const link = props['\ub9c1\ud06c \ud398\uc774\uc9c0']?.rich_text?.map(t => t.plain_text).join('') || '';
  const lastEdited = props['\uc0dd\uc131 \uc77c\uc790']?.last_edited_time || '';
  const createdTime = props['\uccab \uc0dd\uc131']?.created_time || '';
  let imageUrl = '';
  const photoFiles = props['\uc0ac\uc9c4']?.files;
  if (photoFiles && photoFiles.length > 0) { const file = photoFiles[0]; imageUrl = file.file?.url || file.external?.url || ''; }
  return { title, description, category, importance, link, lastEdited, createdTime, imageUrl };
}

async function getDataSourceId() {
  if (cachedDataSourceId) return cachedDataSourceId;
  const db = await notion.databases.retrieve({ database_id: DATABASE_ID });
  if (db.data_sources && db.data_sources.length > 0) { const ds = db.data_sources[0]; cachedDataSourceId = ds.data_source_id || ds.id; }
  else { cachedDataSourceId = db.data_source_id || db.id || DATABASE_ID; }
  return cachedDataSourceId;
}

async function fetchPageList() {
  const dataSourceId = await getDataSourceId();
  const pages = []; let cursor;
  do {
    const response = await notion.dataSources.query({ data_source_id: dataSourceId, start_cursor: cursor, page_size: 100 });
    pages.push(...response.results);
    cursor = response.has_more ? response.next_cursor : undefined;
  } while (cursor);
  return pages;
}

async function fetchAllPages(force = false) {
  if (!force && cachedPages.length > 0) return cachedPages;
  console.log('[Fetch] Full fetch from Notion...');
  const pages = await fetchPageList();
  console.log('[Fetch] Found ' + pages.length + ' pages. Fetching content...');
  const enrichedPages = []; let count = 0;
  for (const page of pages) {
    const props = extractPageProperties(page);
    let pageData;
    try { pageData = await getPageContent(page.id); } catch (e) { pageData = { text: '', firstImage: '' }; }
    const finalImageUrl = props.imageUrl || pageData.firstImage;
    enrichedPages.push({ id: page.id, url: page.url, ...props, imageUrl: finalImageUrl, content: pageData.text.slice(0, 3000) });
    count++;
    if (count % 50 === 0) console.log('[Fetch] Progress: ' + count + '/' + pages.length);
    await new Promise(r => setTimeout(r, 50));
  }
  cachedPages = enrichedPages; lastFetchTime = Date.now();
  console.log('[Fetch] Done! ' + enrichedPages.length + ' pages loaded.');
  saveCacheToDisk(); return enrichedPages;
}

async function incrementalSync() {
  try {
    const dataSourceId = await getDataSourceId();
    const response = await notion.dataSources.query({ data_source_id: dataSourceId, page_size: 100, sorts: [{ property: '\uc0dd\uc131 \uc77c\uc790', direction: 'descending' }] });
    if (response.results.length === 0) return;
    const latestEditTime = response.results[0].last_edited_time;
    if (lastSyncEditTime && latestEditTime === lastSyncEditTime) return;
    console.log('[Sync] Changes detected! Updating...');
    let updatedCount = 0;
    for (const page of response.results) {
      const editTime = page.last_edited_time;
      if (lastSyncEditTime && editTime <= lastSyncEditTime) break;
      const props = extractPageProperties(page);
      let pageData;
      try { pageData = await getPageContent(page.id); } catch (e) { pageData = { text: '', firstImage: '' }; }
      const finalImageUrl = props.imageUrl || pageData.firstImage;
      const enriched = { id: page.id, url: page.url, ...props, imageUrl: finalImageUrl, content: pageData.text.slice(0, 3000) };
      const existingIdx = cachedPages.findIndex(p => p.id === page.id);
      if (existingIdx >= 0) { cachedPages[existingIdx] = enriched; } else { cachedPages.push(enriched); }
      updatedCount++;
    }
    lastSyncEditTime = latestEditTime; lastFetchTime = Date.now(); saveCacheToDisk();
    console.log('[Sync] Updated ' + updatedCount + ' pages. Total: ' + cachedPages.length);
  } catch (e) { console.error('[Sync] Error:', e.message); }
}

function searchPages(pages, query) {
  const queryLower = query.toLowerCase();
  const queryTerms = queryLower.split(/\\s+/).filter(t => t.length > 1);
  const scored = pages.map(page => {
    const searchText = (page.title + ' ' + page.description + ' ' + page.category + ' ' + page.content).toLowerCase();
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
  return scored.filter(p => p.score > 0).sort((a, b) => b.score - a.score).slice(0, 20);
}

app.post('/api/chat', async (req, res) => {
  try {
    const { message, history = [] } = req.body;
    if (!message) return res.status(400).json({ error: 'Message is required' });
    if (cachedPages.length === 0) return res.status(503).json({ error: 'Server is still loading data. Please wait a moment.' });
    const relevantPages = searchPages(cachedPages, message);
    const context = relevantPages.map((page, i) => {
      return '[Page ' + (i+1) + '] Title: ' + page.title + '\\nCategory: ' + page.category + '\\nImportance: ' + page.importance + '\\nDescription: ' + page.description + '\\nContent: ' + page.content.slice(0, 1500) + '\\nURL: ' + page.url;
    }).join('\\n\\n---\\n\\n');
    const messages = [];
    for (const h of history.slice(-6)) { messages.push({ role: h.role, content: h.content }); }
    messages.push({ role: 'user', content: message });
    const systemPrompt = 'You are a knowledgeable research assistant that answers questions based on a Notion database containing notes about:\\n- All-solid-state batteries (ASSB)\\n- NCM cathode materials\\n- Solid electrolytes (LPSCl, sulfide-based)\\n- DEM (Discrete Element Method) simulation\\n- Electrode fabrication and characterization\\n- Meeting notes, seminar content, and paper summaries\\n\\nIMPORTANT RULES:\\n1. Answer in the same language as the user\\'s question (Korean or English).\\n2. Base your answers ONLY on the provided database content below.\\n3. If the database doesn\\'t contain relevant information, say so honestly.\\n4. Always cite which pages you referenced using their titles.\\n5. Be concise but thorough.\\n\\n--- DATABASE CONTENT ---\\n' + (context || 'No relevant pages found for this query.') + '\\n--- END DATABASE CONTENT ---';
    const response = await anthropic.messages.create({ model: 'claude-sonnet-4-20250514', max_tokens: 2000, system: systemPrompt, messages });
    const answer = response.content[0].text;
    const references = relevantPages.map(p => ({ id: p.id, title: p.title, url: p.url, category: p.category, importance: p.importance, description: p.description, imageUrl: p.imageUrl }));
    res.json({ answer, references });
  } catch (error) { console.error('Chat error:', error); res.status(500).json({ error: error.message || 'Internal server error' }); }
});

app.get('/api/pages', (req, res) => {
  const summaries = cachedPages.map(p => ({ id: p.id, title: p.title, url: p.url, category: p.category, importance: p.importance, description: p.description, lastEdited: p.lastEdited }));
  res.json({ pages: summaries, total: summaries.length });
});

app.get('/api/search', (req, res) => {
  const { q } = req.query;
  if (!q) return res.json({ results: [] });
  const results = searchPages(cachedPages, q);
  res.json({ results: results.map(p => ({ id: p.id, title: p.title, url: p.url, category: p.category, importance: p.importance, description: p.description, score: p.score })) });
});

app.post('/api/refresh', async (req, res) => {
  try { const pages = await fetchAllPages(true); res.json({ message: 'Refreshed ' + pages.length + ' pages' }); }
  catch (error) { res.status(500).json({ error: error.message }); }
});

const SAVE_DATABASE_ID = '32f53577a0a58048a64fe440a87bd17d';

app.post('/api/save-to-notion', async (req, res) => {
  try {
    const { title, content } = req.body;
    if (!title || !content) return res.status(400).json({ error: 'Title and content are required' });
    const page = await notion.pages.create({
      parent: { database_id: SAVE_DATABASE_ID },
      properties: {
        '\uc774\ub984': { title: [{ text: { content: title } }] },
        '\ud14d\uc2a4\ud2b8': { rich_text: [{ text: { content: content.slice(0, 2000) } }] },
      },
    });
    res.json({ success: true, url: page.url, id: page.id });
  } catch (error) { console.error('Save to Notion error:', error); res.status(500).json({ error: error.message }); }
});

app.get('/api/status', (req, res) => {
  res.json({ totalPages: cachedPages.length, lastFetchTime: lastFetchTime ? new Date(lastFetchTime).toISOString() : null, cacheFileExists: fs.existsSync(CACHE_FILE), autoSyncEnabled: true, syncIntervalMs: SYNC_INTERVAL });
});

app.use((req, res, next) => {
  if (req.path.startsWith('/api')) return next();
  res.sendFile(path.join(distPath, 'index.html'));
});

const PORT = process.env.PORT || 3001;
app.listen(PORT, '0.0.0.0', () => {
  console.log('Server running on http://localhost:' + PORT);
  const hasCache = loadCacheFromDisk();
  if (hasCache) { console.log('[Startup] Using disk cache. Starting incremental sync...'); incrementalSync().catch(console.error); }
  else { console.log('[Startup] No cache found. Doing full fetch...'); fetchAllPages().catch(console.error); }
  setInterval(incrementalSync, SYNC_INTERVAL);
  console.log('[Sync] Auto-sync every ' + (SYNC_INTERVAL / 1000) + 's');
});
