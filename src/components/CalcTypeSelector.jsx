import { CALC_TYPES } from '../config/calculations';

export default function CalcTypeSelector({ calcTypes, setCalcTypes }) {
  const toggle = (id) => {
    setCalcTypes((prev) =>
      prev.includes(id) ? prev.filter((c) => c !== id) : [...prev, id]
    );
  };

  return (
    <section className="card">
      <h2>4. Basic Calculation Type</h2>
      <div className="calc-type-grid">
        {CALC_TYPES.map((ct) => (
          <label key={ct.id} className={`calc-type-item ${calcTypes.includes(ct.id) ? 'active' : ''}`}>
            <input
              type="checkbox"
              checked={calcTypes.includes(ct.id)}
              onChange={() => toggle(ct.id)}
            />
            <div>
              <span className="calc-type-label">{ct.label}</span>
              <span className="calc-type-desc">{ct.description}</span>
            </div>
          </label>
        ))}
      </div>
    </section>
  );
}
