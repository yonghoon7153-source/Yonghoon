import { PP_LIBRARIES } from '../config/calculations';

export default function PseudopotentialSelector({ params, setParams, selectedCode }) {
  if (selectedCode !== 'qe') return null;

  return (
    <section className="card">
      <h2>2. Pseudopotential Library</h2>
      <div className="pp-grid">
        {PP_LIBRARIES.map((pp) => (
          <button
            key={pp.id}
            className={`pp-btn ${params.ppLibrary === pp.id ? 'active' : ''}`}
            onClick={() => setParams({ ...params, ppLibrary: pp.id })}
          >
            {pp.label}
          </button>
        ))}
      </div>
      <div className="input-group">
        <label>Pseudopotential directory path:</label>
        <input
          type="text"
          value={params.ppPath}
          onChange={(e) => setParams({ ...params, ppPath: e.target.value })}
          placeholder="./pseudo"
        />
      </div>
    </section>
  );
}
