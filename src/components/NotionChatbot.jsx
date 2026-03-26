import { useState, useRef, useEffect } from 'react';
import './NotionChatbot.css';

const API_BASE = '';

function NotionChatbot() {
  const [messages, setMessages] = useState([]);
  const [input, setInput] = useState('');
  const [isLoading, setIsLoading] = useState(false);
  const [expandedRefs, setExpandedRefs] = useState({});
  const [savingIdx, setSavingIdx] = useState(null);
  const [savedIdx, setSavedIdx] = useState(new Set());
  const messagesEndRef = useRef(null);
  const inputRef = useRef(null);

  useEffect(() => {
    messagesEndRef.current?.scrollIntoView({ behavior: 'smooth' });
  }, [messages]);

  const toggleRef = (msgIndex) => {
    setExpandedRefs(prev => ({ ...prev, [msgIndex]: !prev[msgIndex] }));
  };

  const sendMessage = async () => {
    const trimmed = input.trim();
    if (!trimmed || isLoading) return;

    const userMsg = { role: 'user', content: trimmed };
    setMessages(prev => [...prev, userMsg]);
    setInput('');
    setIsLoading(true);

    try {
      const history = messages.map(m => ({
        role: m.role,
        content: m.content,
      }));

      const res = await fetch(`${API_BASE}/api/chat`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ message: trimmed, history }),
      });

      if (!res.ok) {
        const err = await res.json();
        throw new Error(err.error || 'Server error');
      }

      const data = await res.json();
      const assistantMsg = {
        role: 'assistant',
        content: data.answer,
        references: data.references || [],
      };
      setMessages(prev => [...prev, assistantMsg]);
    } catch (error) {
      setMessages(prev => [
        ...prev,
        {
          role: 'assistant',
          content: `Error: ${error.message}. Please check that the server is running.`,
          references: [],
        },
      ]);
    } finally {
      setIsLoading(false);
      inputRef.current?.focus();
    }
  };

  const handleKeyDown = (e) => {
    if (e.key === 'Enter' && !e.shiftKey) {
      e.preventDefault();
      sendMessage();
    }
  };

  const saveToNotion = async (msgIdx) => {
    const msg = messages[msgIdx];
    if (!msg || msg.role !== 'assistant') return;

    setSavingIdx(msgIdx);
    try {
      const userMsg = messages[msgIdx - 1];
      const question = userMsg ? userMsg.content : '';
      const title = question.length > 50 ? question.slice(0, 50) + '...' : question || 'Notion Chat Answer';
      const content = 'Q: ' + question + '\n\nA: ' + msg.content;

      const res = await fetch(API_BASE + '/api/save-to-notion', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ title, content }),
      });

      if (!res.ok) {
        const err = await res.json();
        throw new Error(err.error || 'Failed to save');
      }

      const data = await res.json();
      setSavedIdx(prev => new Set([...prev, msgIdx]));
      window.open(data.url, '_blank');
    } catch (error) {
      alert('Notion save failed: ' + error.message);
    } finally {
      setSavingIdx(null);
    }
  };

  const importanceStars = (importance) => {
    if (!importance) return null;
    return <span className="ref-importance">{importance}</span>;
  };

  const renderMarkdown = (text) => {
    if (!text) return '';
    let html = text
      .replace(/^### (.+)$/gm, '<h4>$1</h4>')
      .replace(/^## (.+)$/gm, '<h3>$1</h3>')
      .replace(/^# (.+)$/gm, '<h2>$1</h2>')
      .replace(/\*\*(.+?)\*\*/g, '<strong>$1</strong>')
      .replace(/\*(.+?)\*/g, '<em>$1</em>')
      .replace(/`(.+?)`/g, '<code>$1</code>')
      .replace(/^- (.+)$/gm, '<li>$1</li>')
      .replace(/^\d+\. (.+)$/gm, '<li>$1</li>')
      .replace(/\n\n/g, '</p><p>')
      .replace(/\n/g, '<br/>');
    html = html.replace(/((?:<li>.*?<\/li>\s*)+)/g, '<ul>$1</ul>');
    return '<p>' + html + '</p>';
  };

  const suggestions = [
    '\uc804\uace0\uccb4 \ubc30\ud130\ub9ac\uc5d0\uc11c pressure\uac00 \uc911\uc694\ud55c \uc774\uc720\ub294?',
    'Ag-C anode interlayer\uc5d0 \ub300\ud574 \uc124\uba85\ud574\uc918',
    'DEM simulation\uc774\ub780 \ubb34\uc5c7\uc778\uac00?',
    'NCM cathode\uc758 degradation \uba54\ucee4\ub2c8\uc998\uc740?',
  ];

  return (
    <div className="chatbot-container">
      <div className="chatbot-header">
        <div className="chatbot-header-icon">
          <svg width="24" height="24" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
            <path d="M21 15a2 2 0 0 1-2 2H7l-4 4V5a2 2 0 0 1 2-2h14a2 2 0 0 1 2 2z" />
          </svg>
        </div>
        <div>
          <h2 className="chatbot-title">Notion Knowledge Base Chat</h2>
          <p className="chatbot-subtitle">Ask questions about your research database</p>
        </div>
      </div>

      <div className="chatbot-messages">
        {messages.length === 0 && (
          <div className="chatbot-welcome">
            <div className="welcome-icon">
              <svg width="48" height="48" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="1.5">
                <path d="M12 6.253v13m0-13C10.832 5.477 9.246 5 7.5 5S4.168 5.477 3 6.253v13C4.168 18.477 5.754 18 7.5 18s3.332.477 4.5 1.253m0-13C13.168 5.477 14.754 5 16.5 5c1.747 0 3.332.477 4.5 1.253v13C19.832 18.477 18.247 18 16.5 18c-1.746 0-3.332.477-4.5 1.253" />
              </svg>
            </div>
            <h3>Welcome!</h3>
            <p>Ask me anything about your Notion research database.</p>
            <div className="welcome-suggestions">
              {suggestions.map((s, i) => (
                <button key={i} onClick={() => setInput(s)}>{s}</button>
              ))}
            </div>
          </div>
        )}

        {messages.map((msg, idx) => (
          <div key={idx} className={`chat-message ${msg.role}`}>
            <div className="message-avatar">
              {msg.role === 'user' ? (
                <svg width="20" height="20" viewBox="0 0 24 24" fill="currentColor">
                  <path d="M12 12c2.21 0 4-1.79 4-4s-1.79-4-4-4-4 1.79-4 4 1.79 4 4 4zm0 2c-2.67 0-8 1.34-8 4v2h16v-2c0-2.66-5.33-4-8-4z" />
                </svg>
              ) : (
                <svg width="20" height="20" viewBox="0 0 24 24" fill="currentColor">
                  <path d="M12 2C6.48 2 2 6.48 2 12s4.48 10 10 10 10-4.48 10-10S17.52 2 12 2zm-1 17.93c-3.95-.49-7-3.85-7-7.93 0-.62.08-1.21.21-1.79L9 15v1c0 1.1.9 2 2 2v1.93zm6.9-2.54c-.26-.81-1-1.39-1.9-1.39h-1v-3c0-.55-.45-1-1-1H8v-2h2c.55 0 1-.45 1-1V7h2c1.1 0 2-.9 2-2v-.41c2.93 1.19 5 4.06 5 7.41 0 2.08-.8 3.97-2.1 5.39z" />
                </svg>
              )}
            </div>
            <div className="message-content">
              {msg.role === 'assistant' ? (
                <>
                  <div className="message-text" dangerouslySetInnerHTML={{ __html: renderMarkdown(msg.content) }} />
                  <button
                    className={'save-to-notion-btn' + (savedIdx.has(idx) ? ' saved' : '')}
                    onClick={() => saveToNotion(idx)}
                    disabled={savingIdx === idx || savedIdx.has(idx)}
                    title="Save to Notion"
                  >
                    {savingIdx === idx ? (
                      <span className="save-spinner"></span>
                    ) : savedIdx.has(idx) ? (
                      <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2.5"><polyline points="20 6 9 17 4 12" /></svg>
                    ) : (
                      <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2"><path d="M19 21H5a2 2 0 0 1-2-2V5a2 2 0 0 1 2-2h11l5 5v11a2 2 0 0 1-2 2z" /><polyline points="17 21 17 13 7 13 7 21" /><polyline points="7 3 7 8 15 8" /></svg>
                    )}
                    {savedIdx.has(idx) ? 'Saved' : 'Save to Notion'}
                  </button>
                </>
              ) : (
                <div className="message-text">{msg.content}</div>
              )}
              {msg.references && msg.references.length > 0 && (
                <div className="message-references">
                  <button
                    className="ref-toggle"
                    onClick={() => toggleRef(idx)}
                  >
                    <svg
                      width="16"
                      height="16"
                      viewBox="0 0 24 24"
                      fill="none"
                      stroke="currentColor"
                      strokeWidth="2"
                      className={expandedRefs[idx] ? 'ref-arrow expanded' : 'ref-arrow'}
                    >
                      <polyline points="6 9 12 15 18 9" />
                    </svg>
                    Referenced Pages ({msg.references.length})
                  </button>
                  {expandedRefs[idx] && (
                    <div className="ref-list">
                      {msg.references.map((ref, refIdx) => (
                        <a
                          key={refIdx}
                          href={ref.url}
                          target="_blank"
                          rel="noopener noreferrer"
                          className="ref-card"
                        >
                          <div className="ref-card-header">
                            <span className="ref-title">{ref.title || 'Untitled'}</span>
                            {importanceStars(ref.importance)}
                          </div>
                          {ref.category && (
                            <span className="ref-category">{ref.category}</span>
                          )}
                          {ref.description && (
                            <p className="ref-description">{ref.description}</p>
                          )}
                        </a>
                      ))}
                    </div>
                  )}
                </div>
              )}
            </div>
          </div>
        ))}

        {isLoading && (
          <div className="chat-message assistant">
            <div className="message-avatar">
              <svg width="20" height="20" viewBox="0 0 24 24" fill="currentColor">
                <path d="M12 2C6.48 2 2 6.48 2 12s4.48 10 10 10 10-4.48 10-10S17.52 2 12 2zm-1 17.93c-3.95-.49-7-3.85-7-7.93 0-.62.08-1.21.21-1.79L9 15v1c0 1.1.9 2 2 2v1.93zm6.9-2.54c-.26-.81-1-1.39-1.9-1.39h-1v-3c0-.55-.45-1-1-1H8v-2h2c.55 0 1-.45 1-1V7h2c1.1 0 2-.9 2-2v-.41c2.93 1.19 5 4.06 5 7.41 0 2.08-.8 3.97-2.1 5.39z" />
              </svg>
            </div>
            <div className="message-content">
              <div className="typing-indicator">
                <span></span><span></span><span></span>
              </div>
            </div>
          </div>
        )}

        <div ref={messagesEndRef} />
      </div>

      <div className="chatbot-input-area">
        <textarea
          ref={inputRef}
          value={input}
          onChange={(e) => setInput(e.target.value)}
          onKeyDown={handleKeyDown}
          placeholder="Ask a question about your research..."
          rows={1}
          disabled={isLoading}
        />
        <button
          onClick={sendMessage}
          disabled={!input.trim() || isLoading}
          className="send-button"
        >
          <svg width="20" height="20" viewBox="0 0 24 24" fill="currentColor">
            <path d="M2.01 21L23 12 2.01 3 2 10l15 2-15 2z" />
          </svg>
        </button>
      </div>
    </div>
  );
}

export default NotionChatbot;
