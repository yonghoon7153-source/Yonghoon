import { useState } from 'react';
import { POST_PROCESSING } from '../config/calculations';

function Tooltip({ text }) {
  const [pos, setPos] = useState(null);
  const handleEnter = (e) => {
    const rect = e.currentTarget.getBoundingClientRect();
    setPos({ top: rect.top - 8, left: rect.left + rect.width / 2 });
  };
  return (
    <span
      className="tooltip-wrapper"
      onMouseEnter={handleEnter}
      onMouseLeave={() => setPos(null)}
    >
      <span className="tooltip-icon">?</span>
      {pos && (
        <span
          className="tooltip-text"
          style={{ top: pos.top, left: pos.left, transform: 'translate(-50%, -100%)' }}
        >
          {text}
        </span>
      )}
    </span>
  );
}

export default function PostProcessing({ selectedCode, postProcessing, setPostProcessing, analysis }) {
  const [collapsed, setCollapsed] = useState({});

  const toggle = (id) => {
    setPostProcessing((prev) =>
      prev.includes(id) ? prev.filter((p) => p !== id) : [...prev, id]
    );
  };

  const toggleCategory = (catIndex) => {
    setCollapsed((prev) => ({ ...prev, [catIndex]: !prev[catIndex] }));
  };

  const selectAll = (items) => {
    const availableIds = items.filter(it => it.codes.includes(selectedCode)).map(it => it.id);
    const allSelected = availableIds.every(id => postProcessing.includes(id));
    if (allSelected) {
      setPostProcessing(prev => prev.filter(id => !availableIds.includes(id)));
    } else {
      setPostProcessing(prev => [...new Set([...prev, ...availableIds])]);
    }
  };

  return (
    <section className="card">
      <h2>5. Post-Processing</h2>
      <p className="hint">Select calculations to perform. Dependencies are auto-resolved.</p>
      <div className="pp-categories">
        {POST_PROCESSING.map((cat, catIdx) => {
          const availableItems = cat.items.filter((it) => it.codes.includes(selectedCode));
          if (availableItems.length === 0) return null;

          const isCollapsed = collapsed[catIdx];
          const selectedCount = availableItems.filter(it => postProcessing.includes(it.id)).length;

          return (
            <div key={catIdx} className="pp-category">
              <div className="pp-category-header" onClick={() => toggleCategory(catIdx)}>
                <span className={`collapse-arrow ${isCollapsed ? '' : 'open'}`}>&#9654;</span>
                <h3>{cat.category}</h3>
                {selectedCount > 0 && (
                  <span className="selected-badge">{selectedCount}/{availableItems.length}</span>
                )}
                <button
                  className="btn-select-all"
                  onClick={(e) => {
                    e.stopPropagation();
                    selectAll(cat.items);
                  }}
                >
                  {availableItems.every(it => postProcessing.includes(it.id)) ? 'Deselect All' : 'Select All'}
                </button>
              </div>
              {!isCollapsed && (
                <div className="pp-items">
                  {availableItems.map((item) => {
                    const rec = analysis?.postRecommendations?.[item.id];
                    const isRec = rec?.recommended;
                    return (
                      <label
                        key={item.id}
                        className={`pp-item ${postProcessing.includes(item.id) ? 'active' : ''} ${isRec ? 'recommended' : ''}`}
                      >
                        <input
                          type="checkbox"
                          checked={postProcessing.includes(item.id)}
                          onChange={() => toggle(item.id)}
                        />
                        <span className="pp-item-label">
                          {item.label}
                          {isRec && <span className="rec-dot" title={rec.reason}>&#x25CF;</span>}
                        </span>
                        <Tooltip text={isRec ? `${item.description}\n★ ${rec.reason}` : item.description} />
                      </label>
                    );
                  })}
                </div>
              )}
            </div>
          );
        })}
      </div>
    </section>
  );
}
