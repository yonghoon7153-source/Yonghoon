import { DEPENDENCIES, CALC_TYPES } from '../config/calculations';

export default function DependencyInfo({ postProcessing, calcTypes }) {
  // Compute auto-included dependencies
  const autoDeps = new Set();
  for (const pp of postProcessing) {
    const deps = DEPENDENCIES[pp] || [];
    for (const dep of deps) {
      if (CALC_TYPES.find(ct => ct.id === dep) && !calcTypes.includes(dep)) {
        autoDeps.add(dep);
      }
    }
  }

  if (autoDeps.size === 0) return null;

  const depLabels = {
    scf: 'SCF',
    nscf: 'NSCF',
    relax: 'Relaxation',
    'vc-relax': 'VC-Relax',
    bands: 'Bands',
    phonon: 'Phonon',
    aimd: 'AIMD',
  };

  return (
    <div className="dependency-info">
      <span className="dep-icon">&#9432;</span>
      <span>
        Auto-included prerequisites:{' '}
        {[...autoDeps].map((d) => (
          <span key={d} className="dep-tag">{depLabels[d] || d}</span>
        ))}
      </span>
    </div>
  );
}
