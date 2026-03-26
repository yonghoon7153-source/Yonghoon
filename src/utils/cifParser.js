import { getPseudopotentialFile } from '../config/pseudopotentials';

// Simple CIF parser to extract lattice parameters and atomic positions
export function parseCIF(cifText) {
  const lines = cifText.split('\n').map(l => l.trim()).filter(l => l && !l.startsWith('#'));

  const result = {
    a: 0, b: 0, c: 0,
    alpha: 90, beta: 90, gamma: 90,
    spaceGroup: 'P 1',
    atoms: [],
    cellParams: {},
  };

  // Extract cell parameters
  for (const line of lines) {
    if (line.startsWith('_cell_length_a')) result.a = parseFloat(line.split(/\s+/)[1]) || 0;
    if (line.startsWith('_cell_length_b')) result.b = parseFloat(line.split(/\s+/)[1]) || 0;
    if (line.startsWith('_cell_length_c')) result.c = parseFloat(line.split(/\s+/)[1]) || 0;
    if (line.startsWith('_cell_angle_alpha')) result.alpha = parseFloat(line.split(/\s+/)[1]) || 90;
    if (line.startsWith('_cell_angle_beta')) result.beta = parseFloat(line.split(/\s+/)[1]) || 90;
    if (line.startsWith('_cell_angle_gamma')) result.gamma = parseFloat(line.split(/\s+/)[1]) || 90;
    if (line.startsWith('_symmetry_space_group_name_H-M') || line.startsWith('_space_group_name_H-M_alt')) {
      const match = line.match(/'([^']+)'|"([^"]+)"|(\S+)$/);
      if (match) result.spaceGroup = match[1] || match[2] || match[3];
    }
  }

  // Find atom site loop
  let inAtomLoop = false;
  let atomHeaders = [];
  let headerIndex = {};

  for (let i = 0; i < lines.length; i++) {
    const line = lines[i];

    if (line === 'loop_') {
      // Check if next lines are atom_site headers
      const nextLine = lines[i + 1] || '';
      if (nextLine.startsWith('_atom_site')) {
        inAtomLoop = true;
        atomHeaders = [];
        headerIndex = {};
        continue;
      }
    }

    if (inAtomLoop && line.startsWith('_atom_site')) {
      const header = line.trim();
      headerIndex[header] = atomHeaders.length;
      atomHeaders.push(header);
      continue;
    }

    if (inAtomLoop && !line.startsWith('_') && line !== 'loop_') {
      const parts = line.split(/\s+/);
      if (parts.length >= atomHeaders.length) {
        const getVal = (key) => {
          const idx = headerIndex[key];
          return idx !== undefined ? parts[idx] : undefined;
        };

        const label = getVal('_atom_site_label') || getVal('_atom_site_type_symbol') || 'X';
        const symbol = getVal('_atom_site_type_symbol') || label.replace(/[0-9+\-]/g, '');
        const x = parseFloat(getVal('_atom_site_fract_x')) || 0;
        const y = parseFloat(getVal('_atom_site_fract_y')) || 0;
        const z = parseFloat(getVal('_atom_site_fract_z')) || 0;

        result.atoms.push({ label, symbol: symbol.replace(/[0-9+\-]/g, ''), x, y, z });
      } else if (parts.length < 3) {
        inAtomLoop = false;
      }
    }

    if (inAtomLoop && (line === 'loop_' || line.startsWith('_') && !line.startsWith('_atom_site'))) {
      inAtomLoop = false;
    }
  }

  // Compute cell vectors from parameters
  result.cellParams = computeCellVectors(result.a, result.b, result.c, result.alpha, result.beta, result.gamma);

  // Get unique elements
  result.elements = [...new Set(result.atoms.map(a => a.symbol))];
  result.nat = result.atoms.length;
  result.ntyp = result.elements.length;

  return result;
}

function computeCellVectors(a, b, c, alpha, beta, gamma) {
  const toRad = (deg) => deg * Math.PI / 180;
  const cosA = Math.cos(toRad(alpha));
  const cosB = Math.cos(toRad(beta));
  const cosG = Math.cos(toRad(gamma));
  const sinG = Math.sin(toRad(gamma));

  const v1 = [a, 0, 0];
  const v2 = [b * cosG, b * sinG, 0];
  const cx = c * cosB;
  const cy = c * (cosA - cosB * cosG) / sinG;
  const cz = Math.sqrt(c * c - cx * cx - cy * cy);
  const v3 = [cx, cy, cz];

  return { v1, v2, v3 };
}

// Convert cell vectors to Bohr (QE uses Bohr internally, but we use angstrom with ibrav=0)
export function formatCellParameters(parsed) {
  const { v1, v2, v3 } = parsed.cellParams;
  return `CELL_PARAMETERS angstrom
  ${v1.map(v => v.toFixed(10)).join('  ')}
  ${v2.map(v => v.toFixed(10)).join('  ')}
  ${v3.map(v => v.toFixed(10)).join('  ')}`;
}

export function formatAtomicPositions(parsed) {
  let out = 'ATOMIC_POSITIONS crystal\n';
  for (const atom of parsed.atoms) {
    out += `  ${atom.symbol.padEnd(4)} ${atom.x.toFixed(10)}  ${atom.y.toFixed(10)}  ${atom.z.toFixed(10)}\n`;
  }
  return out.trimEnd();
}

export function formatAtomicSpecies(parsed, ppPath, ppLibrary) {
  const massMap = {
    H: 1.008, He: 4.003, Li: 6.941, Be: 9.012, B: 10.81, C: 12.011, N: 14.007, O: 15.999,
    F: 18.998, Ne: 20.180, Na: 22.990, Mg: 24.305, Al: 26.982, Si: 28.086, P: 30.974,
    S: 32.065, Cl: 35.453, Ar: 39.948, K: 39.098, Ca: 40.078, Sc: 44.956, Ti: 47.867,
    V: 50.942, Cr: 51.996, Mn: 54.938, Fe: 55.845, Co: 58.933, Ni: 58.693, Cu: 63.546,
    Zn: 65.38, Ga: 69.723, Ge: 72.64, As: 74.922, Se: 78.96, Br: 79.904, Kr: 83.798,
    Rb: 85.468, Sr: 87.62, Y: 88.906, Zr: 91.224, Nb: 92.906, Mo: 95.96, Tc: 98.0,
    Ru: 101.07, Rh: 102.906, Pd: 106.42, Ag: 107.868, Cd: 112.411, In: 114.818,
    Sn: 118.710, Sb: 121.760, Te: 127.60, I: 126.904, Xe: 131.293, Cs: 132.905,
    Ba: 137.327, La: 138.905, Ce: 140.116, Pr: 140.908, Nd: 144.242, Pm: 145.0,
    Sm: 150.36, Eu: 151.964, Gd: 157.25, Tb: 158.925, Dy: 162.500, Ho: 164.930,
    Er: 167.259, Tm: 168.934, Yb: 173.054, Lu: 174.967, Hf: 178.49, Ta: 180.948,
    W: 183.84, Re: 186.207, Os: 190.23, Ir: 192.217, Pt: 195.084, Au: 196.967,
    Hg: 200.59, Tl: 204.383, Pb: 207.2, Bi: 208.980, Po: 209.0, At: 210.0,
    Rn: 222.0, Fr: 223.0, Ra: 226.0, Ac: 227.0, Th: 232.038, Pa: 231.036,
    U: 238.029, Np: 237.0, Pu: 244.0
  };

  let out = 'ATOMIC_SPECIES\n';
  for (const el of parsed.elements) {
    const mass = massMap[el] || 1.0;
    const ppFile = getPseudopotentialFile(el, ppLibrary);
    out += `  ${el.padEnd(4)} ${mass.toFixed(3)}  ${ppFile}\n`;
  }
  return out.trimEnd();
}
