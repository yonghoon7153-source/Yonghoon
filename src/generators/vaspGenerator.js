// VASP input file generator: INCAR, POSCAR, KPOINTS, POTCAR setup, job script

export function generateVASPScripts(parsed, params, calcTypes, postProcessing) {
  const files = {};

  const kx = params.kx ?? 4;
  const ky = params.ky ?? 4;
  const kz = params.kz ?? 4;
  const encut = Math.round((params.ecutwfc ?? 60) * 13.605693); // Ry → eV
  const functional = params.functional || 'pbe';
  const dftd3 = params.dftd3 || false;
  const spinPolarized = params.spinPolarized || false;
  const hubbardU = params.hubbardU || '';
  const ppPath = params.ppPath || '/opt/vasp/potentials/PBE';
  const nprocs = params.nprocs || 4;

  const needSCF = calcTypes.includes('scf') || postProcessing.length > 0;
  const needRelax = calcTypes.includes('relax') || ['adhesion_energy', 'surface_energy', 'vacancy_energy', 'neb'].some(p => postProcessing.includes(p));
  const needVCRelax = calcTypes.includes('vc-relax');
  const needNSCF = calcTypes.includes('nscf') || ['dos', 'pdos', 'bands', 'bandgap', 'effmass', 'fermi_surface', 'wannier', 'boltztrap', 'optical_absorption', 'dielectric', 'refractive_index'].some(p => postProcessing.includes(p));
  const needBands = postProcessing.includes('bands') || postProcessing.includes('effmass');
  const mixingBeta = parsed.nat > 50 ? 0.3 : parsed.nat > 20 ? 0.4 : 0.7;

  // --- POSCAR ---
  files['POSCAR'] = generatePOSCAR(parsed);

  // --- POTCAR setup script ---
  files['setup_potcar.sh'] = generatePOTCARScript(parsed, ppPath);

  // --- SCF (static) ---
  if (needSCF) {
    files['INCAR_scf'] = generateINCAR('scf', {
      encut, functional, dftd3, spinPolarized, hubbardU, parsed, mixingBeta,
      extra: {},
    });
    files['KPOINTS_scf'] = generateKPOINTS(kx, ky, kz, 'Monkhorst-Pack');
  }

  // --- Relaxation ---
  if (needRelax) {
    files['INCAR_relax'] = generateINCAR('relax', {
      encut, functional, dftd3, spinPolarized, hubbardU, parsed, mixingBeta,
      extra: { ISIF: 2, NSW: 200, IBRION: 2, EDIFFG: -0.01 },
    });
    files['KPOINTS_relax'] = generateKPOINTS(kx, ky, kz, 'Monkhorst-Pack');
  }

  // --- VC-Relax ---
  if (needVCRelax) {
    files['INCAR_vc-relax'] = generateINCAR('vc-relax', {
      encut, functional, dftd3, spinPolarized, hubbardU, parsed, mixingBeta,
      extra: { ISIF: 3, NSW: 200, IBRION: 2, EDIFFG: -0.01 },
    });
    files['KPOINTS_vc-relax'] = generateKPOINTS(kx, ky, kz, 'Monkhorst-Pack');
  }

  // --- DOS ---
  if (postProcessing.includes('dos') || postProcessing.includes('pdos')) {
    files['INCAR_dos'] = generateINCAR('dos', {
      encut, functional, dftd3, spinPolarized, hubbardU, parsed, mixingBeta,
      extra: {
        ICHARG: 11, ISMEAR: -5, SIGMA: 0.05,
        LORBIT: 11, NEDOS: 3001,
        EMIN: -20, EMAX: 20,
      },
    });
    const dkx = Math.min(kx * 2, 12);
    const dky = Math.min(ky * 2, 12);
    const dkz = Math.min(kz * 2, 12);
    files['KPOINTS_dos'] = generateKPOINTS(dkx, dky, dkz, 'Monkhorst-Pack');
  }

  // --- Band structure ---
  if (needBands) {
    files['INCAR_bands'] = generateINCAR('bands', {
      encut, functional, dftd3, spinPolarized, hubbardU, parsed, mixingBeta,
      extra: { ICHARG: 11, ISMEAR: 0, SIGMA: 0.05, LORBIT: 11 },
    });
    files['KPOINTS_bands'] = `k-points along high-symmetry lines
40
Line-mode
Reciprocal
  0.0000  0.0000  0.0000  ! GAMMA
  0.5000  0.0000  0.0000  ! X

  0.5000  0.0000  0.0000  ! X
  0.5000  0.5000  0.0000  ! M

  0.5000  0.5000  0.0000  ! M
  0.0000  0.0000  0.0000  ! GAMMA

  0.0000  0.0000  0.0000  ! GAMMA
  0.0000  0.0000  0.5000  ! Z
`;
  }

  // --- Bader ---
  if (postProcessing.includes('bader')) {
    files['INCAR_bader'] = generateINCAR('bader', {
      encut, functional, dftd3, spinPolarized, hubbardU, parsed, mixingBeta,
      extra: { LAECHG: '.TRUE.', LCHARG: '.TRUE.', NSW: 0 },
    });
    files['KPOINTS_bader'] = generateKPOINTS(kx, ky, kz, 'Monkhorst-Pack');
  }

  // --- CDD ---
  if (postProcessing.includes('cdd')) {
    files['INCAR_cdd'] = generateINCAR('cdd', {
      encut, functional, dftd3, spinPolarized, hubbardU, parsed, mixingBeta,
      extra: { LCHARG: '.TRUE.', LPARD: '.TRUE.', NSW: 0 },
    });
  }

  // --- ELF ---
  if (postProcessing.includes('elf')) {
    files['INCAR_elf'] = generateINCAR('elf', {
      encut, functional, dftd3, spinPolarized, hubbardU, parsed, mixingBeta,
      extra: { LELF: '.TRUE.', NSW: 0 },
    });
  }

  // --- Work function ---
  if (postProcessing.includes('work_function')) {
    files['INCAR_workfunc'] = generateINCAR('work_function', {
      encut, functional, dftd3, spinPolarized, hubbardU, parsed, mixingBeta,
      extra: { LVHAR: '.TRUE.', NSW: 0, IDIPOL: 3, LDIPOL: '.TRUE.' },
    });
  }

  // --- STM ---
  if (postProcessing.includes('stm')) {
    files['INCAR_stm'] = generateINCAR('stm', {
      encut, functional, dftd3, spinPolarized, hubbardU, parsed, mixingBeta,
      extra: { LPARD: '.TRUE.', NBMOD: -3, EINT: '-0.5 0.0', NSW: 0 },
    });
  }

  // --- Optical / dielectric ---
  if (postProcessing.includes('optical_absorption') || postProcessing.includes('dielectric') || postProcessing.includes('refractive_index')) {
    files['INCAR_optics'] = generateINCAR('optics', {
      encut, functional, dftd3, spinPolarized, hubbardU, parsed, mixingBeta,
      extra: {
        LOPTICS: '.TRUE.', CSHIFT: 0.1, NEDOS: 2000,
        NBANDS: Math.max(parsed.nat * 4, 100),
      },
    });
    const okx = Math.min(kx * 2, 12);
    const oky = Math.min(ky * 2, 12);
    const okz = Math.min(kz * 2, 12);
    files['KPOINTS_optics'] = generateKPOINTS(okx, oky, okz, 'Gamma');
  }

  // --- Elastic ---
  if (postProcessing.includes('elastic')) {
    files['INCAR_elastic'] = generateINCAR('elastic', {
      encut: Math.round(encut * 1.3), functional, dftd3, spinPolarized, hubbardU, parsed, mixingBeta,
      extra: { IBRION: 6, ISIF: 3, NSW: 1, NFREE: 4, EDIFF: 1e-8 },
    });
    files['KPOINTS_elastic'] = generateKPOINTS(kx, ky, kz, 'Monkhorst-Pack');
  }

  // --- Piezoelectric ---
  if (postProcessing.includes('piezoelectric')) {
    files['INCAR_piezo'] = generateINCAR('piezo', {
      encut, functional, dftd3, spinPolarized, hubbardU, parsed, mixingBeta,
      extra: { IBRION: 6, ISIF: 3, NSW: 1, NFREE: 4, LEPSILON: '.TRUE.', LCALCEPS: '.TRUE.' },
    });
  }

  // --- Born effective charge ---
  if (postProcessing.includes('born')) {
    files['INCAR_born'] = generateINCAR('born', {
      encut, functional, dftd3, spinPolarized, hubbardU, parsed, mixingBeta,
      extra: { LEPSILON: '.TRUE.', IBRION: 8, NSW: 1 },
    });
  }

  // --- Phonon (VASP+phonopy) ---
  if (postProcessing.includes('phonon') || postProcessing.includes('thermo') || postProcessing.includes('thermal_conductivity')) {
    files['INCAR_phonon'] = generateINCAR('phonon', {
      encut, functional, dftd3, spinPolarized, hubbardU, parsed, mixingBeta,
      extra: { IBRION: -1, NSW: 0, EDIFF: 1e-8, PREC: 'Accurate', LWAVE: '.FALSE.', LCHARG: '.FALSE.' },
    });
    files['KPOINTS_phonon'] = generateKPOINTS(kx, ky, kz, 'Monkhorst-Pack');
    files['phonon_workflow.sh'] = `#!/bin/bash
# Phonon calculation with VASP + phonopy
set -e

echo "=== Step 1: Create phonopy displacement files ==="
phonopy -d --dim="2 2 2" --vasp

echo "=== Step 2: Run VASP for each displacement ==="
for dir in POSCAR-*; do
  num=\${dir#POSCAR-}
  mkdir -p disp-\${num}
  cp \${dir} disp-\${num}/POSCAR
  cp INCAR_phonon disp-\${num}/INCAR
  cp KPOINTS_phonon disp-\${num}/KPOINTS
  cp POTCAR disp-\${num}/POTCAR
  cd disp-\${num}
  mpirun -np ${nprocs} vasp_std
  cd ..
  echo "  Displacement \${num} done"
done

echo "=== Step 3: Collect forces ==="
phonopy -f disp-*/vasprun.xml

echo "=== Step 4: Post-processing ==="
phonopy -p -s band.conf
phonopy -t mesh.conf

echo "Phonon calculation completed!"
`;
    files['band.conf'] = `ATOM_NAME = ${parsed.elements.join(' ')}
DIM = 2 2 2
BAND = 0 0 0  1/2 0 0  1/2 1/2 0  0 0 0  0 0 1/2
BAND_LABELS = \\Gamma X M \\Gamma Z
`;
    files['mesh.conf'] = `ATOM_NAME = ${parsed.elements.join(' ')}
DIM = 2 2 2
MP = 20 20 20
TPROP = .TRUE.
TMIN = 0
TMAX = 1000
TSTEP = 10
`;
  }

  // --- AIMD ---
  if (postProcessing.includes('aimd') || postProcessing.includes('rdf') || postProcessing.includes('msd') || postProcessing.includes('stress_fluct')) {
    const mdTemp = params.mdTemp || 300;
    const mdSteps = params.mdSteps || 5000;
    const mdTimestep = params.mdTimestep || 1.0;
    files['INCAR_aimd'] = generateINCAR('aimd', {
      encut, functional, dftd3, spinPolarized, hubbardU, parsed, mixingBeta,
      extra: {
        IBRION: 0, NSW: mdSteps, POTIM: mdTimestep,
        SMASS: 0, TEBEG: mdTemp, TEEND: mdTemp,
        ISYM: 0, PREC: 'Normal',
        LWAVE: '.FALSE.', LCHARG: '.FALSE.',
        ISIF: 2, EDIFF: 1e-5,
      },
    });
    const mkx = Math.max(1, Math.floor(kx / 2));
    const mky = Math.max(1, Math.floor(ky / 2));
    const mkz = Math.max(1, Math.floor(kz / 2));
    files['KPOINTS_aimd'] = generateKPOINTS(mkx, mky, mkz, 'Gamma');
  }

  // --- NEB ---
  if (postProcessing.includes('neb')) {
    files['INCAR_neb'] = generateINCAR('neb', {
      encut, functional, dftd3, spinPolarized, hubbardU, parsed, mixingBeta,
      extra: {
        IBRION: 3, POTIM: 0.0,
        IOPT: 1, ICHAIN: 0,
        IMAGES: 5, SPRING: -5,
        LCLIMB: '.TRUE.',
        NSW: 200, EDIFFG: -0.03,
        LWAVE: '.FALSE.', LCHARG: '.FALSE.',
      },
    });
    files['KPOINTS_neb'] = generateKPOINTS(kx, ky, kz, 'Monkhorst-Pack');
    files['neb_note.txt'] = `# VASP NEB calculation
# Workflow:
# 1. Optimize initial and final structures
# 2. Create directories: 00/ 01/ 02/ 03/ 04/ 05/ 06/
#    - 00/POSCAR = initial, 06/POSCAR = final
#    - 01-05/ = interpolated images (use nebmake.pl or VTST tools)
# 3. Copy INCAR_neb → INCAR, KPOINTS_neb → KPOINTS, POTCAR to run directory
# 4. Run: mpirun -np ${nprocs} vasp_std
# 5. Use nebbarrier.pl or ASE to extract barrier
`;
  }

  // --- Magnetic moments ---
  if (postProcessing.includes('magnetic_moments')) {
    files['INCAR_mag'] = generateINCAR('magnetic', {
      encut, functional, dftd3, spinPolarized: true, hubbardU, parsed, mixingBeta,
      extra: { LORBIT: 11 },
    });
  }

  // --- MAE ---
  if (postProcessing.includes('mae')) {
    files['INCAR_mae_x'] = generateINCAR('mae_x', {
      encut, functional, dftd3, spinPolarized: true, hubbardU, parsed, mixingBeta,
      extra: {
        LSORBIT: '.TRUE.', LNONCOLLINEAR: '.TRUE.',
        SAXIS: '1 0 0', LORBMOM: '.TRUE.',
        EDIFF: 1e-8, NBANDS: Math.max(parsed.nat * 4, 100),
      },
    });
    files['INCAR_mae_z'] = generateINCAR('mae_z', {
      encut, functional, dftd3, spinPolarized: true, hubbardU, parsed, mixingBeta,
      extra: {
        LSORBIT: '.TRUE.', LNONCOLLINEAR: '.TRUE.',
        SAXIS: '0 0 1', LORBMOM: '.TRUE.',
        EDIFF: 1e-8, NBANDS: Math.max(parsed.nat * 4, 100),
      },
    });
  }

  // --- BoltzTraP ---
  if (postProcessing.includes('boltztrap')) {
    files['INCAR_boltztrap'] = generateINCAR('boltztrap', {
      encut, functional, dftd3, spinPolarized, hubbardU, parsed, mixingBeta,
      extra: { ICHARG: 11, ISMEAR: 0, SIGMA: 0.01, LORBIT: 11 },
    });
    const bkx = Math.min(kx * 3, 16);
    const bky = Math.min(ky * 3, 16);
    const bkz = Math.min(kz * 3, 16);
    files['KPOINTS_boltztrap'] = generateKPOINTS(bkx, bky, bkz, 'Monkhorst-Pack');
  }

  // --- Run script ---
  files['run_all.sh'] = generateVASPRunScript(
    needSCF, needRelax, needVCRelax, needNSCF, needBands, postProcessing, nprocs
  );

  return files;
}

// ---- POSCAR ----
function generatePOSCAR(parsed) {
  const { v1, v2, v3 } = parsed.cellParams;
  const elementCounts = {};
  for (const atom of parsed.atoms) {
    elementCounts[atom.symbol] = (elementCounts[atom.symbol] || 0) + 1;
  }
  const elements = Object.keys(elementCounts);
  const counts = elements.map(e => elementCounts[e]);

  // Sort atoms by element order
  const sortedAtoms = [];
  for (const el of elements) {
    sortedAtoms.push(...parsed.atoms.filter(a => a.symbol === el));
  }

  let out = `${parsed.elements.join('-')} structure\n`;
  out += `1.0\n`;
  out += `  ${v1.map(v => v.toFixed(10)).join('  ')}\n`;
  out += `  ${v2.map(v => v.toFixed(10)).join('  ')}\n`;
  out += `  ${v3.map(v => v.toFixed(10)).join('  ')}\n`;
  out += `  ${elements.join('  ')}\n`;
  out += `  ${counts.join('  ')}\n`;
  out += `Direct\n`;
  for (const atom of sortedAtoms) {
    out += `  ${atom.x.toFixed(10)}  ${atom.y.toFixed(10)}  ${atom.z.toFixed(10)}\n`;
  }
  return out;
}

// ---- KPOINTS ----
function generateKPOINTS(kx, ky, kz, scheme = 'Monkhorst-Pack') {
  return `Automatic mesh
0
${scheme}
  ${kx}  ${ky}  ${kz}
  0  0  0
`;
}

// ---- POTCAR setup ----
function generatePOTCARScript(parsed, ppPath) {
  const elements = parsed.elements;
  // Common recommended POTCAR versions
  const potcarMap = {
    Li: 'Li_sv', Na: 'Na_pv', K: 'K_sv', Rb: 'Rb_sv', Cs: 'Cs_sv',
    Ca: 'Ca_sv', Sr: 'Sr_sv', Ba: 'Ba_sv',
    Sc: 'Sc_sv', Ti: 'Ti_sv', V: 'V_sv', Cr: 'Cr_pv', Mn: 'Mn_pv',
    Fe: 'Fe_pv', Co: 'Co', Ni: 'Ni_pv', Cu: 'Cu_pv', Zn: 'Zn',
    Y: 'Y_sv', Zr: 'Zr_sv', Nb: 'Nb_sv', Mo: 'Mo_sv',
    Ga: 'Ga_d', Ge: 'Ge_d', In: 'In_d', Sn: 'Sn_d',
    La: 'La', Ce: 'Ce', Pr: 'Pr_3', Nd: 'Nd_3',
    Hf: 'Hf_pv', Ta: 'Ta_pv', W: 'W_sv', Re: 'Re_pv',
    Os: 'Os_pv', Ir: 'Ir', Pt: 'Pt', Au: 'Au',
    Pb: 'Pb_d', Bi: 'Bi_d',
  };

  let script = `#!/bin/bash
# Generate POTCAR by concatenating pseudopotentials
# POTCAR path: ${ppPath}
# Elements: ${elements.join(', ')}

POTCAR_DIR="${ppPath}"
cat `;

  const potcars = elements.map(el => {
    const pp = potcarMap[el] || el;
    return `"\${POTCAR_DIR}/${pp}/POTCAR"`;
  });

  script += potcars.join(' \\\n    ');
  script += ` > POTCAR

echo "POTCAR generated for: ${elements.join(' ')}"
echo "Using: ${elements.map(el => potcarMap[el] || el).join(' ')}"

# Verify
grep TITEL POTCAR
`;

  return script;
}

// ---- INCAR ----
function generateINCAR(calcType, opts) {
  const { encut, functional, dftd3, spinPolarized, hubbardU, parsed, mixingBeta, extra } = opts;

  let lines = [];
  lines.push(`# VASP INCAR - ${calcType}`);
  lines.push(`# System: ${parsed.elements.join('-')}, ${parsed.nat} atoms`);
  lines.push('');

  // General
  lines.push('# General');
  lines.push(`SYSTEM = ${parsed.elements.join('-')}_${calcType}`);
  lines.push(`ENCUT = ${encut}`);
  lines.push(`PREC = ${extra.PREC || 'Accurate'}`);
  lines.push(`EDIFF = ${extra.EDIFF || 1e-6}`);
  lines.push(`NELM = 200`);
  lines.push(`AMIX = ${mixingBeta}`);
  lines.push('');

  // Electronic
  lines.push('# Electronic');
  if (extra.ISMEAR !== undefined) {
    lines.push(`ISMEAR = ${extra.ISMEAR}`);
    lines.push(`SIGMA = ${extra.SIGMA || 0.05}`);
  } else {
    lines.push('ISMEAR = 0');
    lines.push('SIGMA = 0.05');
  }
  lines.push('LREAL = Auto');
  lines.push('');

  // Functional
  if (functional === 'hse06') {
    lines.push('# HSE06 hybrid functional');
    lines.push('LHFCALC = .TRUE.');
    lines.push('HFSCREEN = 0.2');
    lines.push('AEXX = 0.25');
    lines.push('ALGO = Damped');
    lines.push('TIME = 0.4');
    lines.push('');
  } else if (functional === 'scan') {
    lines.push('# SCAN meta-GGA');
    lines.push('METAGGA = SCAN');
    lines.push('LASPH = .TRUE.');
    lines.push('');
  } else if (functional === 'pbesol') {
    lines.push('GGA = PS');
    lines.push('');
  }

  // Spin
  if (spinPolarized) {
    lines.push('# Spin polarization');
    lines.push('ISPIN = 2');
    const magmom = parsed.elements.map(el => {
      const count = parsed.atoms.filter(a => a.symbol === el).length;
      const mag = ['Fe', 'Co', 'Ni', 'Mn', 'Cr', 'V'].includes(el) ? 5.0 : 0.6;
      return `${count}*${mag}`;
    }).join(' ');
    lines.push(`MAGMOM = ${magmom}`);
    lines.push('');
  }

  // DFT-D3
  if (dftd3) {
    lines.push('# DFT-D3 dispersion');
    lines.push('IVDW = 12');
    lines.push('');
  }

  // Hubbard U
  if (hubbardU) {
    lines.push('# DFT+U');
    lines.push('LDAU = .TRUE.');
    lines.push('LDAUTYPE = 2');
    const ldul = parsed.elements.map(el =>
      ['Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn'].includes(el) ? 2 : -1
    ).join(' ');
    const lduu = parsed.elements.map(el =>
      ['Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn'].includes(el) ? parseFloat(hubbardU) : 0
    ).join(' ');
    const lduj = parsed.elements.map(() => 0).join(' ');
    lines.push(`LDAUL = ${ldul}`);
    lines.push(`LDAUU = ${lduu}`);
    lines.push(`LDAUJ = ${lduj}`);
    lines.push('');
  }

  // Ionic relaxation
  if (extra.NSW !== undefined && extra.NSW > 0) {
    lines.push('# Ionic relaxation');
    lines.push(`NSW = ${extra.NSW}`);
    lines.push(`IBRION = ${extra.IBRION ?? 2}`);
    if (extra.ISIF !== undefined) lines.push(`ISIF = ${extra.ISIF}`);
    if (extra.EDIFFG) lines.push(`EDIFFG = ${extra.EDIFFG}`);
    if (extra.POTIM !== undefined && extra.POTIM > 0) lines.push(`POTIM = ${extra.POTIM}`);
    lines.push('');
  } else if (extra.NSW === 0 || (!extra.NSW && !extra.IBRION)) {
    lines.push('NSW = 0');
    lines.push('');
  }

  // AIMD specific
  if (extra.SMASS !== undefined) {
    lines.push('# AIMD settings');
    lines.push(`SMASS = ${extra.SMASS}`);
    lines.push(`TEBEG = ${extra.TEBEG}`);
    lines.push(`TEEND = ${extra.TEEND}`);
    if (extra.ISYM !== undefined) lines.push(`ISYM = ${extra.ISYM}`);
    lines.push('');
  }

  // NEB specific
  if (extra.IMAGES !== undefined) {
    lines.push('# NEB settings (VTST)');
    lines.push(`IMAGES = ${extra.IMAGES}`);
    lines.push(`SPRING = ${extra.SPRING}`);
    lines.push(`LCLIMB = ${extra.LCLIMB}`);
    if (extra.IOPT !== undefined) lines.push(`IOPT = ${extra.IOPT}`);
    if (extra.ICHAIN !== undefined) lines.push(`ICHAIN = ${extra.ICHAIN}`);
    lines.push('');
  }

  // Extra flags
  const handled = ['PREC', 'EDIFF', 'ISMEAR', 'SIGMA', 'NSW', 'IBRION', 'ISIF', 'EDIFFG',
    'POTIM', 'SMASS', 'TEBEG', 'TEEND', 'ISYM', 'IMAGES', 'SPRING', 'LCLIMB', 'IOPT', 'ICHAIN'];
  const remaining = Object.entries(extra).filter(([k]) => !handled.includes(k));
  if (remaining.length > 0) {
    lines.push('# Additional settings');
    for (const [key, val] of remaining) {
      lines.push(`${key} = ${val}`);
    }
    lines.push('');
  }

  // Output
  lines.push('# Output');
  if (extra.LWAVE === undefined) lines.push('LWAVE = .FALSE.');
  if (extra.LCHARG === undefined) lines.push('LCHARG = .TRUE.');
  lines.push('');

  return lines.join('\n');
}

// ---- Run script ----
function generateVASPRunScript(needSCF, needRelax, needVCRelax, needNSCF, needBands, postProcessing, nprocs) {
  let s = `#!/bin/bash
# ==============================================
# VASP Workflow Script
# Auto-generated
# ==============================================
set -e
NPROCS=${nprocs}
MPI="mpirun -np \${NPROCS}"
VASP="vasp_std"

# Generate POTCAR first
echo "Setting up POTCAR..."
bash setup_potcar.sh

run_vasp() {
  local step=\$1
  local incar=\$2
  local kpoints=\$3

  echo "[\${step}] Running VASP..."
  cp \${incar} INCAR
  cp \${kpoints} KPOINTS
  \${MPI} \${VASP}
  mkdir -p results_\${step}
  cp OUTCAR OSZICAR CONTCAR vasprun.xml results_\${step}/ 2>/dev/null || true
  [ -f DOSCAR ] && cp DOSCAR results_\${step}/ || true
  [ -f EIGENVAL ] && cp EIGENVAL results_\${step}/ || true
  [ -f CHGCAR ] && cp CHGCAR results_\${step}/ || true
  echo "  \${step} completed."
}
`;

  if (needVCRelax) {
    s += `\nrun_vasp "vc-relax" INCAR_vc-relax KPOINTS_vc-relax\ncp CONTCAR POSCAR\n`;
  }
  if (needRelax) {
    s += `\nrun_vasp "relax" INCAR_relax KPOINTS_relax\ncp CONTCAR POSCAR\n`;
  }
  if (needSCF) {
    s += `\nrun_vasp "scf" INCAR_scf KPOINTS_scf\n`;
  }
  if (postProcessing.includes('dos') || postProcessing.includes('pdos')) {
    s += `\n# Copy CHGCAR from SCF for DOS\ncp results_scf/CHGCAR .\nrun_vasp "dos" INCAR_dos KPOINTS_dos\n`;
  }
  if (postProcessing.includes('bands') || postProcessing.includes('effmass')) {
    s += `\ncp results_scf/CHGCAR .\nrun_vasp "bands" INCAR_bands KPOINTS_bands\n`;
  }
  if (postProcessing.includes('bader')) {
    s += `\nrun_vasp "bader" INCAR_bader KPOINTS_bader\necho "Run: chgsum.pl AECCAR0 AECCAR2 && bader CHGCAR -ref CHGCAR_sum"\n`;
  }
  if (postProcessing.includes('elf')) {
    s += `\nrun_vasp "elf" INCAR_elf KPOINTS_scf\n`;
  }
  if (postProcessing.includes('work_function')) {
    s += `\nrun_vasp "workfunc" INCAR_workfunc KPOINTS_scf\n`;
  }
  if (postProcessing.includes('optical_absorption') || postProcessing.includes('dielectric') || postProcessing.includes('refractive_index')) {
    s += `\nrun_vasp "optics" INCAR_optics KPOINTS_optics\n`;
  }
  if (postProcessing.includes('elastic')) {
    s += `\nrun_vasp "elastic" INCAR_elastic KPOINTS_elastic\necho "Parse OUTCAR for elastic tensor (Cij)"\n`;
  }
  if (postProcessing.includes('born')) {
    s += `\nrun_vasp "born" INCAR_born KPOINTS_scf\n`;
  }
  if (postProcessing.includes('aimd') || postProcessing.includes('rdf') || postProcessing.includes('msd') || postProcessing.includes('stress_fluct')) {
    s += `\nrun_vasp "aimd" INCAR_aimd KPOINTS_aimd\n`;
  }
  if (postProcessing.includes('phonon') || postProcessing.includes('thermo') || postProcessing.includes('thermal_conductivity')) {
    s += `\necho "Running phonon workflow..."\nbash phonon_workflow.sh\n`;
  }
  if (postProcessing.includes('mae')) {
    s += `\nrun_vasp "mae_x" INCAR_mae_x KPOINTS_scf\nrun_vasp "mae_z" INCAR_mae_z KPOINTS_scf\necho "MAE = E(x) - E(z)"\n`;
  }
  if (postProcessing.includes('neb')) {
    s += `\necho "NEB: See neb_note.txt for setup instructions"\n`;
  }

  s += `\necho "=========================================="\necho "All VASP calculations completed!"\necho "=========================================="\n`;
  return s;
}
