// Analyze parsed CIF structure and recommend calculations + parameters

const TRANSITION_METALS = ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
  'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
  'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg'];

const MAGNETIC_ELEMENTS = ['Mn', 'Fe', 'Co', 'Ni', 'Cr', 'V', 'Cu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er'];

const HEAVY_ELEMENTS = ['Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi',
  'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu'];

const ALKALI = ['Li', 'Na', 'K', 'Rb', 'Cs'];
const HALIDES = ['F', 'Cl', 'Br', 'I'];
const CHALCOGENIDES = ['O', 'S', 'Se', 'Te'];

// Recommended ecutwfc by element type (Ry)
const ECUTWFC_MAP = {
  'O': 70, 'N': 70, 'F': 70, 'C': 60, 'H': 50,
  'Li': 50, 'Na': 50, 'K': 50, 'S': 50, 'P': 50,
  'Cl': 55, 'Br': 45, 'I': 45, 'Se': 45, 'Te': 45,
  'Si': 45, 'Ge': 45, 'Al': 40, 'Ga': 45,
  'Fe': 70, 'Co': 65, 'Ni': 65, 'Mn': 70, 'Cr': 65,
  'Ti': 55, 'V': 55, 'Cu': 55, 'Zn': 50,
  'default': 60,
};

export function analyzeStructure(parsed) {
  const elements = parsed.elements || [];
  const nat = parsed.nat || 0;
  const a = parsed.a || 0;
  const b = parsed.b || 0;
  const c = parsed.c || 0;

  const result = {
    calcTypeRecommendations: {},  // id -> { recommended, reason }
    postRecommendations: {},      // id -> { recommended, reason }
    paramRecommendations: {},     // key -> value
    systemType: '',
    description: '',
  };

  const hasAlkali = elements.some(e => ALKALI.includes(e));
  const hasHalide = elements.some(e => HALIDES.includes(e));
  const hasChalcogenide = elements.some(e => CHALCOGENIDES.includes(e));
  const hasTM = elements.some(e => TRANSITION_METALS.includes(e));
  const hasMagnetic = elements.some(e => MAGNETIC_ELEMENTS.includes(e));
  const hasHeavy = elements.some(e => HEAVY_ELEMENTS.includes(e));
  const hasP = elements.includes('P');
  const hasS = elements.includes('S');
  const hasO = elements.includes('O');
  const hasLi = elements.includes('Li');
  const hasNa = elements.includes('Na');

  // Detect system type
  const isIonicConductor = (hasLi || hasNa) && (hasHalide || (hasP && hasS) || (hasP && hasO));
  const isBattery = isIonicConductor;
  const isPerovskite = nat >= 5 && (hasO || hasHalide) && elements.length >= 3;
  const isMagnetic = hasMagnetic;
  const isSemiconductor = !hasTM && (elements.includes('Si') || elements.includes('Ge') ||
    (elements.includes('Ga') && (elements.includes('As') || elements.includes('N'))));

  // System description
  if (isIonicConductor) {
    result.systemType = 'ionic_conductor';
    result.description = `Ionic conductor (${hasLi ? 'Li' : 'Na'}-ion) — ${elements.join('-')} system, ${nat} atoms`;
  } else if (isMagnetic) {
    result.systemType = 'magnetic';
    result.description = `Magnetic system — ${elements.join('-')}, ${nat} atoms`;
  } else if (isSemiconductor) {
    result.systemType = 'semiconductor';
    result.description = `Semiconductor — ${elements.join('-')}, ${nat} atoms`;
  } else {
    result.systemType = 'general';
    result.description = `${elements.join('-')} system, ${nat} atoms, ${parsed.ntyp} species`;
  }

  // === Calc type recommendations ===
  result.calcTypeRecommendations['scf'] = { recommended: true, reason: 'Basic energy calculation (always needed)' };
  result.calcTypeRecommendations['relax'] = { recommended: true, reason: 'Optimize atomic positions' };
  result.calcTypeRecommendations['vc-relax'] = { recommended: true, reason: 'Optimize cell + positions for accurate lattice' };
  result.calcTypeRecommendations['nscf'] = { recommended: false, reason: '' };

  // === Post-processing recommendations ===

  // Always recommend
  result.postRecommendations['dos'] = { recommended: true, reason: 'Essential electronic structure analysis' };
  result.postRecommendations['pdos'] = { recommended: true, reason: 'Orbital-resolved electronic structure' };
  result.postRecommendations['bands'] = { recommended: true, reason: 'Band structure visualization' };
  result.postRecommendations['bandgap'] = { recommended: true, reason: 'Determine metallic/insulating character' };
  result.postRecommendations['bond_length'] = { recommended: true, reason: 'Basic structural analysis' };

  // Ionic conductor specific
  if (isIonicConductor) {
    result.postRecommendations['aimd'] = { recommended: true, reason: `${hasLi ? 'Li' : 'Na'}-ion diffusion study` };
    result.postRecommendations['msd'] = { recommended: true, reason: 'Diffusion coefficient & ionic conductivity' };
    result.postRecommendations['rdf'] = { recommended: true, reason: 'Local structure analysis during MD' };
    result.postRecommendations['neb'] = { recommended: true, reason: `${hasLi ? 'Li' : 'Na'} migration barrier` };
    result.postRecommendations['phonon'] = { recommended: true, reason: 'Dynamical stability check' };
    result.postRecommendations['bader'] = { recommended: true, reason: 'Charge transfer analysis' };
    result.postRecommendations['elastic'] = { recommended: true, reason: 'Mechanical stability' };
  }

  // Magnetic systems
  if (isMagnetic) {
    result.postRecommendations['magnetic_moments'] = { recommended: true, reason: 'Magnetic ordering analysis' };
    result.postRecommendations['mae'] = { recommended: true, reason: 'Magnetic anisotropy' };
    result.postRecommendations['cohp'] = { recommended: true, reason: 'Bonding analysis for TM compounds' };
  }

  // Semiconductor
  if (isSemiconductor) {
    result.postRecommendations['effmass'] = { recommended: true, reason: 'Carrier effective mass' };
    result.postRecommendations['dielectric'] = { recommended: true, reason: 'Dielectric response' };
    result.postRecommendations['optical_absorption'] = { recommended: true, reason: 'Optical properties' };
    result.postRecommendations['boltztrap'] = { recommended: true, reason: 'Thermoelectric properties' };
  }

  // Heavy elements → SOC relevant
  if (hasHeavy) {
    result.postRecommendations['mae'] = { recommended: true, reason: 'Strong spin-orbit coupling effects' };
  }

  // Phonon for thermal properties
  if (nat <= 50) {
    result.postRecommendations['phonon'] = { recommended: true, reason: 'Dynamical stability & thermal properties' };
    result.postRecommendations['thermo'] = { recommended: true, reason: 'Thermodynamic functions (Cv, entropy, F)' };
  }

  // === Parameter recommendations ===

  // ecutwfc: take max of all elements
  let maxEcut = 0;
  for (const el of elements) {
    const ecut = ECUTWFC_MAP[el] || ECUTWFC_MAP['default'];
    if (ecut > maxEcut) maxEcut = ecut;
  }
  result.paramRecommendations.ecutwfc = maxEcut;
  result.paramRecommendations.ecutrho = maxEcut * 10;

  // k-points: based on cell size (a, b, c)
  // Rule of thumb: k * a ≈ 30-40 Å for good convergence
  const kTarget = nat > 50 ? 25 : 35;
  result.paramRecommendations.kx = a > 0.1 ? Math.max(2, Math.round(kTarget / a)) : 4;
  result.paramRecommendations.ky = b > 0.1 ? Math.max(2, Math.round(kTarget / b)) : 4;
  result.paramRecommendations.kz = c > 0.1 ? Math.max(1, Math.round(kTarget / c)) : 4;

  // Smearing
  if (isMagnetic || hasTM) {
    result.paramRecommendations.smearing = 'mv';
    result.paramRecommendations.degauss = 0.02;
  } else {
    result.paramRecommendations.smearing = 'gaussian';
    result.paramRecommendations.degauss = 0.02;
  }

  // Spin polarized for magnetic systems
  result.paramRecommendations.spinPolarized = isMagnetic;

  // DFT-D3 for layered/molecular systems or systems with vdW
  const needsVdW = hasS || elements.includes('Se') || elements.includes('Te') ||
    (c / a > 2.5) || (c / b > 2.5);
  result.paramRecommendations.dftd3 = needsVdW;

  // MD temperature for ionic conductors
  if (isIonicConductor) {
    result.paramRecommendations.mdTemp = 600;
    result.paramRecommendations.mdSteps = 10000;
  }

  return result;
}
