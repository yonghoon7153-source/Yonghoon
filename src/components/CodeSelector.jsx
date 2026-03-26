import { CODES } from '../config/calculations';

export default function CodeSelector({ selectedCode, setSelectedCode }) {
  return (
    <section className="card">
      <h2>3. Calculation Code</h2>
      <div className="code-grid">
        {CODES.map((code) => (
          <button
            key={code.id}
            className={`code-btn ${selectedCode === code.id ? 'active' : ''}`}
            onClick={() => setSelectedCode(code.id)}
          >
            <span className="code-name">{code.label}</span>
            {code.id === 'vasp' && <span className="code-tag">DFT</span>}
            {code.id === 'qe' && <span className="code-tag">DFT</span>}
            {code.id === 'orca' && <span className="code-tag">DFT/QC</span>}
            {(code.id === 'uma' || code.id === 'mace') && <span className="code-tag">MLIP</span>}
          </button>
        ))}
      </div>
    </section>
  );
}
