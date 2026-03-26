import express from 'express';
import cors from 'cors';
import { Client } from '@notionhq/client';
import Anthropic from '@anthropic-ai/sdk';
import dotenv from 'dotenv';

dotenv.config();

const app = express();
app.use(cors());
app.use(express.json());

const notion = new Client({ auth: process.env.NOTION_API_KEY });
const anthropic = new Anthropic({ apiKey: process.env.ANTHROPIC_API_KEY });

const DATABASE_ID = process.env.NOTION_DATABASE_ID || '7b5d9b8c069c4f3bbe40b25f55afdf96';

// Cache for Notion pages
let cachedPages = [];
let lastFetchTime = 0;
const CACHE_TTL = 5 * 60 * 1000; // 5 minutes

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
  const db = await notion.databases.retrieve({ database_id: DATABASE_ID });
  console.log('Database keys:', Object.keys(db));
  // Try multiple possible locations for data source ID
  if (db.data_sources && db.data_sources.length > 0) {
    const ds = db.data_sources[0];
    return ds.data_source_id || ds.id;
  }
  if (db.data_source_id) return db.data_source_id;
  if (db.id) return db.id;
  return DATABASE_ID;
}

async function fetchAllPages() {
  const now = Date.now();
  if (cachedPages.length > 0 && now - lastFetchTime < CACHE_TTL) {
    return cachedPages;
  }

  console.log('Fetching all pages from Notion database...');

  // First, get the data source ID from the database
  const dataSourceId = await getDataSourceId();
  console.log(`Using data source ID: ${dataSourceId}`);

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

  console.log(`Found ${pages.length} pages. Fetching content...`);

  // Fetch content for each page (with rate limiting)
  const enrichedPages = [];
  for (const page of pages) {
    const props = extractPageProperties(page);
    let content = '';
    try {
      content = await getPageMarkdown(page.id);
    } catch (e) {
      console.log(`Failed to fetch content for: ${props.title}`);
    }

    enrichedPages.push({
      id: page.id,
      url: page.url,
      ...props,
      content: content.slice(0, 3000),
    });

    // Rate limit: small delay between requests
    await new Promise(r => setTimeout(r, 100));
  }

  cachedPages = enrichedPages;
  lastFetchTime = now;
  console.log(`Fetched ${enrichedPages.length} pages from Notion.`);
  return enrichedPages;
}

function searchPages(pages, query) {
  const queryLower = query.toLowerCase();
  const queryTerms = queryLower.split(/\s+/).filter(t => t.length > 1);

  const scored = pages.map(page => {
    const searchText = `${page.title} ${page.description} ${page.category} ${page.content}`.toLowerCase();
    let score = 0;

    // Full query match
    if (searchText.includes(queryLower)) score += 10;

    // Individual term matches
    for (const term of queryTerms) {
      const regex = new RegExp(term, 'gi');
      const matches = searchText.match(regex);
      if (matches) score += matches.length;

      // Title match gets bonus
      if (page.title.toLowerCase().includes(term)) score += 5;
      // Description match gets bonus
      if (page.description.toLowerCase().includes(term)) score += 3;
    }

    return { ...page, score };
  });

  return scored
    .filter(p => p.score > 0)
    .sort((a, b) => b.score - a.score)
    .slice(0, 8);
}

// API: Chat endpoint
app.post('/api/chat', async (req, res) => {
  try {
    const { message, history = [] } = req.body;

    if (!message) {
      return res.status(400).json({ error: 'Message is required' });
    }

    // Fetch and search pages
    const allPages = await fetchAllPages();
    const relevantPages = searchPages(allPages, message);

    // Build context from relevant pages
    const context = relevantPages.map((page, i) => {
      return `[Page ${i + 1}] Title: ${page.title}
Category: ${page.category}
Importance: ${page.importance}
Description: ${page.description}
Content: ${page.content.slice(0, 1500)}
URL: ${page.url}`;
    }).join('\n\n---\n\n');

    // Build conversation history
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

    // Return answer with referenced pages
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

// API: Get all page summaries (for browsing)
app.get('/api/pages', async (req, res) => {
  try {
    const pages = await fetchAllPages();
    const summaries = pages.map(p => ({
      id: p.id,
      title: p.title,
      url: p.url,
      category: p.category,
      importance: p.importance,
      description: p.description,
      lastEdited: p.lastEdited,
    }));
    res.json({ pages: summaries, total: summaries.length });
  } catch (error) {
    console.error('Pages error:', error);
    res.status(500).json({ error: error.message });
  }
});

// API: Search pages
app.get('/api/search', async (req, res) => {
  try {
    const { q } = req.query;
    if (!q) return res.json({ results: [] });

    const pages = await fetchAllPages();
    const results = searchPages(pages, q);
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
  } catch (error) {
    res.status(500).json({ error: error.message });
  }
});

const PORT = process.env.PORT || 3001;
app.listen(PORT, () => {
  console.log(`Server running on http://localhost:${PORT}`);
  // Pre-fetch pages on startup
  fetchAllPages().catch(console.error);
});
