import { FUNCTIONALS } from '../config/calculations';

export default function ParameterSettings({ params, setParams, selectedCode }) {
  const update = (key, value) => setParams((prev) => ({ ...prev, [key]: value }));

  return (
    <section className="card">
      <h2>6. Calculation Parameters</h2>
      <div className="params-grid">
        {/* Common */}
        <div className="param-group">
          <h4>General</h4>
          <div className="input-group">
            <label>Functional:</label>
            <select value={params.functional} onChange={(e) => update('functional', e.target.value)}>
              {FUNCTIONALS.map((f) => (
                <option key={f.id} value={f.id}>{f.label}</option>
              ))}
            </select>
          </div>
          <div className="input-group">
            <label>Processors (MPI):</label>
            <input
              type="number"
              min={1}
              value={params.nprocs}
              onChange={(e) => update('nprocs', parseInt(e.target.value) || 1)}
            />
          </div>
          <div className="input-group inline">
            <label>
              <input
                type="checkbox"
                checked={params.dftd3}
                onChange={(e) => update('dftd3', e.target.checked)}
              />
              DFT-D3 dispersion correction
            </label>
          </div>
          <div className="input-group inline">
            <label>
              <input
                type="checkbox"
                checked={params.spinPolarized}
                onChange={(e) => update('spinPolarized', e.target.checked)}
              />
              Spin-polarized calculation
            </label>
          </div>
        </div>

        {/* VASP specific */}
        {selectedCode === 'vasp' && (
          <div className="param-group">
            <h4>VASP</h4>
            <div className="input-group">
              <label>ENCUT (eV) — auto from ecutwfc:</label>
              <input
                type="number"
                min={100}
                value={Math.round((params.ecutwfc || 60) * 13.6057)}
                onChange={(e) => update('ecutwfc', Math.round((parseInt(e.target.value) || 500) / 13.6057))}
              />
            </div>
            <div className="input-row">
              <div className="input-group">
                <label>k-points (x):</label>
                <input type="number" min={1} value={params.kx} onChange={(e) => update('kx', parseInt(e.target.value) || 4)} />
              </div>
              <div className="input-group">
                <label>k-points (y):</label>
                <input type="number" min={1} value={params.ky} onChange={(e) => update('ky', parseInt(e.target.value) || 4)} />
              </div>
              <div className="input-group">
                <label>k-points (z):</label>
                <input type="number" min={1} value={params.kz} onChange={(e) => update('kz', parseInt(e.target.value) || 4)} />
              </div>
            </div>
            <div className="input-group">
              <label>POTCAR directory path:</label>
              <input
                type="text"
                value={params.ppPath}
                onChange={(e) => update('ppPath', e.target.value)}
                placeholder="/opt/vasp/potentials/PBE"
              />
            </div>
            <div className="input-group">
              <label>Hubbard U (eV, leave empty if none):</label>
              <input
                type="text"
                value={params.hubbardU}
                onChange={(e) => update('hubbardU', e.target.value)}
                placeholder="e.g., 4.0"
              />
            </div>
          </div>
        )}

        {/* QE specific */}
        {selectedCode === 'qe' && (
          <div className="param-group">
            <h4>Quantum ESPRESSO</h4>
            <div className="input-group">
              <label>ecutwfc (Ry):</label>
              <input
                type="number"
                min={10}
                value={params.ecutwfc}
                onChange={(e) => update('ecutwfc', parseInt(e.target.value) || 60)}
              />
            </div>
            <div className="input-group">
              <label>ecutrho (Ry):</label>
              <input
                type="number"
                min={40}
                value={params.ecutrho}
                onChange={(e) => update('ecutrho', parseInt(e.target.value) || 480)}
              />
            </div>
            <div className="input-row">
              <div className="input-group">
                <label>k-points (x):</label>
                <input
                  type="number"
                  min={1}
                  value={params.kx}
                  onChange={(e) => update('kx', parseInt(e.target.value) || 4)}
                />
              </div>
              <div className="input-group">
                <label>k-points (y):</label>
                <input
                  type="number"
                  min={1}
                  value={params.ky}
                  onChange={(e) => update('ky', parseInt(e.target.value) || 4)}
                />
              </div>
              <div className="input-group">
                <label>k-points (z):</label>
                <input
                  type="number"
                  min={1}
                  value={params.kz}
                  onChange={(e) => update('kz', parseInt(e.target.value) || 4)}
                />
              </div>
            </div>
            <div className="input-group">
              <label>Smearing:</label>
              <select value={params.smearing} onChange={(e) => update('smearing', e.target.value)}>
                <option value="gaussian">Gaussian</option>
                <option value="mv">Marzari-Vanderbilt (cold)</option>
                <option value="mp">Methfessel-Paxton</option>
                <option value="fd">Fermi-Dirac</option>
              </select>
            </div>
            <div className="input-group">
              <label>degauss (Ry):</label>
              <input
                type="number"
                step={0.005}
                min={0.001}
                value={params.degauss}
                onChange={(e) => update('degauss', parseFloat(e.target.value) || 0.02)}
              />
            </div>
            <div className="input-group">
              <label>Hubbard U (eV, leave empty if none):</label>
              <input
                type="text"
                value={params.hubbardU}
                onChange={(e) => update('hubbardU', e.target.value)}
                placeholder="e.g., 4.0"
              />
            </div>
          </div>
        )}

        {/* ORCA specific */}
        {selectedCode === 'orca' && (
          <div className="param-group">
            <h4>ORCA</h4>
            <div className="input-group">
              <label>Basis set:</label>
              <select value={params.basisSet} onChange={(e) => update('basisSet', e.target.value)}>
                <option value="def2-SVP">def2-SVP</option>
                <option value="def2-TZVP">def2-TZVP</option>
                <option value="def2-TZVPP">def2-TZVPP</option>
                <option value="def2-QZVPP">def2-QZVPP</option>
                <option value="cc-pVDZ">cc-pVDZ</option>
                <option value="cc-pVTZ">cc-pVTZ</option>
                <option value="6-31G*">6-31G*</option>
                <option value="6-311G**">6-311G**</option>
              </select>
            </div>
            <div className="input-group">
              <label>Charge:</label>
              <input
                type="number"
                value={params.charge}
                onChange={(e) => update('charge', parseInt(e.target.value) || 0)}
              />
            </div>
            <div className="input-group">
              <label>Multiplicity:</label>
              <input
                type="number"
                min={1}
                value={params.multiplicity}
                onChange={(e) => update('multiplicity', parseInt(e.target.value) || 1)}
              />
            </div>
          </div>
        )}

        {/* MLIP specific */}
        {(selectedCode === 'uma' || selectedCode === 'mace') && (
          <div className="param-group">
            <h4>{selectedCode === 'mace' ? 'MACE' : 'UMA'} Settings</h4>
            <div className="input-group">
              <label>Model path:</label>
              <input
                type="text"
                value={params.mlipModel}
                onChange={(e) => update('mlipModel', e.target.value)}
                placeholder={selectedCode === 'mace' ? 'mace_mp_medium.model' : 'uma_model.pt'}
              />
            </div>
            <div className="input-group">
              <label>MD Temperature (K):</label>
              <input
                type="number"
                min={1}
                value={params.mdTemp}
                onChange={(e) => update('mdTemp', parseInt(e.target.value) || 300)}
              />
            </div>
            <div className="input-group">
              <label>MD Steps:</label>
              <input
                type="number"
                min={100}
                value={params.mdSteps}
                onChange={(e) => update('mdSteps', parseInt(e.target.value) || 5000)}
              />
            </div>
            <div className="input-group">
              <label>MD Timestep (fs):</label>
              <input
                type="number"
                step={0.5}
                min={0.1}
                value={params.mdTimestep}
                onChange={(e) => update('mdTimestep', parseFloat(e.target.value) || 1.0)}
              />
            </div>
          </div>
        )}
      </div>
    </section>
  );
}
