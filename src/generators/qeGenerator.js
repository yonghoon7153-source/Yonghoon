import { formatCellParameters, formatAtomicPositions, formatAtomicSpecies } from '../utils/cifParser';

export function generateQEScripts(parsed, params, calcTypes, postProcessing) {
  const files = {};
  const ppPath = params.ppPath || './pseudo';

  // Resolve all needed calculations
  const needSCF = calcTypes.includes('scf') || postProcessing.length > 0;
  const needNSCF = calcTypes.includes('nscf') || ['dos', 'pdos', 'bands', 'bandgap', 'effmass', 'fermi_surface', 'wannier', 'cohp', 'coop', 'cobi', 'stm', 'optical_absorption', 'refractive_index', 'dielectric', 'boltztrap'].some(p => postProcessing.includes(p));
  const needRelax = calcTypes.includes('relax') || ['adhesion_energy', 'surface_energy', 'vacancy_energy', 'neb'].some(p => postProcessing.includes(p));
  const needVCRelax = calcTypes.includes('vc-relax');
  const needBands = postProcessing.includes('bands') || postProcessing.includes('effmass');

  const kx = params.kx ?? 4;
  const ky = params.ky ?? 4;
  const kz = params.kz ?? 4;
  const ecutwfc = params.ecutwfc ?? 60;
  const ecutrho = params.ecutrho ?? 480;
  const smearing = params.smearing || 'gaussian';
  const degauss = params.degauss || 0.02;
  const functional = params.functional || 'pbe';
  const hubbardU = params.hubbardU || '';
  const dftd3 = params.dftd3 || false;
  const spinPolarized = params.spinPolarized || false;

  const inputFunctional = functional === 'pbe' ? 'pbe' : functional === 'pbesol' ? 'pbesol' : functional === 'lda' ? 'pz' : 'pbe';
  const isHybrid = ['hse06', 'b3lyp'].includes(functional);
  const mixingBeta = parsed.nat > 50 ? 0.3 : parsed.nat > 20 ? 0.4 : 0.7;

  const commonSystem = `  ibrav = 0,
  nat = ${parsed.nat},
  ntyp = ${parsed.ntyp},
  ecutwfc = ${ecutwfc},
  ecutrho = ${ecutrho},${isHybrid ? `\n  input_dft = '${functional}',` : ''}${dftd3 ? `\n  vdw_corr = 'dft-d3',` : ''}${spinPolarized ? `\n  nspin = 2,\n  starting_magnetization(1) = 0.5,` : ''}${hubbardU ? `\n  lda_plus_u = .true.,\n  Hubbard_U(1) = ${hubbardU},` : ''}
  occupations = 'smearing',
  smearing = '${smearing}',
  degauss = ${degauss}`;

  const atomicSection = `${formatAtomicSpecies(parsed, ppPath, params.ppLibrary)}

${formatAtomicPositions(parsed)}

${formatCellParameters(parsed)}

K_POINTS automatic
  ${kx} ${ky} ${kz} 0 0 0`;

  // --- SCF ---
  if (needSCF) {
    files['scf.in'] = `&CONTROL
  calculation = 'scf',
  prefix = 'calc',
  outdir = './tmp',
  pseudo_dir = '${ppPath}',
  tprnfor = .true.,
  tstress = .true.,
/

&SYSTEM
${commonSystem}
/

&ELECTRONS
  conv_thr = 1.0d-8,
  mixing_beta = ${mixingBeta},
/

${formatAtomicSpecies(parsed, ppPath, params.ppLibrary)}

${formatAtomicPositions(parsed)}

${formatCellParameters(parsed)}

K_POINTS automatic
  ${kx} ${ky} ${kz} 0 0 0
`;
  }

  // --- NSCF ---
  if (needNSCF) {
    const nscfKx = Math.min(kx * 2, 12);
    const nscfKy = Math.min(ky * 2, 12);
    const nscfKz = Math.min(kz * 2, 12);
    files['nscf.in'] = `&CONTROL
  calculation = 'nscf',
  prefix = 'calc',
  outdir = './tmp',
  pseudo_dir = '${ppPath}',
/

&SYSTEM
${commonSystem}
/

&ELECTRONS
  conv_thr = 1.0d-8,
/

${formatAtomicSpecies(parsed, ppPath, params.ppLibrary)}

${formatAtomicPositions(parsed)}

${formatCellParameters(parsed)}

K_POINTS automatic
  ${nscfKx} ${nscfKy} ${nscfKz} 0 0 0
`;
  }

  // --- Relax ---
  if (needRelax) {
    files['relax.in'] = `&CONTROL
  calculation = 'relax',
  prefix = 'calc',
  outdir = './tmp',
  pseudo_dir = '${ppPath}',
  tprnfor = .true.,
  tstress = .true.,
/

&SYSTEM
${commonSystem}
/

&ELECTRONS
  conv_thr = 1.0d-8,
  mixing_beta = ${mixingBeta},
/

&IONS
  ion_dynamics = 'bfgs',
/

${formatAtomicSpecies(parsed, ppPath, params.ppLibrary)}

${formatAtomicPositions(parsed)}

${formatCellParameters(parsed)}

K_POINTS automatic
  ${kx} ${ky} ${kz} 0 0 0
`;
  }

  // --- VC-Relax ---
  if (needVCRelax) {
    files['vc-relax.in'] = `&CONTROL
  calculation = 'vc-relax',
  prefix = 'calc',
  outdir = './tmp',
  pseudo_dir = '${ppPath}',
  tprnfor = .true.,
  tstress = .true.,
/

&SYSTEM
${commonSystem}
/

&ELECTRONS
  conv_thr = 1.0d-8,
  mixing_beta = ${mixingBeta},
/

&IONS
  ion_dynamics = 'bfgs',
/

&CELL
  cell_dynamics = 'bfgs',
  press_conv_thr = 0.5,
/

${formatAtomicSpecies(parsed, ppPath, params.ppLibrary)}

${formatAtomicPositions(parsed)}

${formatCellParameters(parsed)}

K_POINTS automatic
  ${kx} ${ky} ${kz} 0 0 0
`;
  }

  // --- Bands ---
  if (needBands) {
    files['bands.in'] = `&CONTROL
  calculation = 'bands',
  prefix = 'calc',
  outdir = './tmp',
  pseudo_dir = '${ppPath}',
/

&SYSTEM
${commonSystem}
/

&ELECTRONS
  conv_thr = 1.0d-8,
/

${formatAtomicSpecies(parsed, ppPath, params.ppLibrary)}

${formatAtomicPositions(parsed)}

${formatCellParameters(parsed)}

K_POINTS crystal_b
4
  0.0000 0.0000 0.0000  30  ! Gamma
  0.5000 0.0000 0.0000  30  ! X
  0.5000 0.5000 0.0000  30  ! M
  0.0000 0.0000 0.0000  30  ! Gamma
`;

    files['bands_pp.in'] = `&BANDS
  prefix = 'calc',
  outdir = './tmp',
  filband = 'bands.dat',
/
`;
  }

  // --- DOS ---
  if (postProcessing.includes('dos')) {
    files['dos.in'] = `&DOS
  prefix = 'calc',
  outdir = './tmp',
  fildos = 'dos.dat',
  Emin = -20.0,
  Emax = 20.0,
  DeltaE = 0.01,
/
`;
  }

  // --- PDOS ---
  if (postProcessing.includes('pdos')) {
    files['pdos.in'] = `&PROJWFC
  prefix = 'calc',
  outdir = './tmp',
  filpdos = 'pdos',
  Emin = -20.0,
  Emax = 20.0,
  DeltaE = 0.01,
/
`;
  }

  // --- Bader ---
  if (postProcessing.includes('bader')) {
    files['pp_charge.in'] = `&INPUTPP
  prefix = 'calc',
  outdir = './tmp',
  filplot = 'charge.dat',
  plot_num = 0,
/

&PLOT
  iflag = 3,
  output_format = 6,
  fileout = 'charge.cube',
/
`;
  }

  // --- CDD ---
  if (postProcessing.includes('cdd')) {
    files['pp_cdd.in'] = `&INPUTPP
  prefix = 'calc',
  outdir = './tmp',
  filplot = 'charge_diff.dat',
  plot_num = 9,
/

&PLOT
  iflag = 3,
  output_format = 6,
  fileout = 'charge_diff.cube',
/
`;
  }

  // --- ELF ---
  if (postProcessing.includes('elf')) {
    files['pp_elf.in'] = `&INPUTPP
  prefix = 'calc',
  outdir = './tmp',
  filplot = 'elf.dat',
  plot_num = 8,
/

&PLOT
  iflag = 3,
  output_format = 6,
  fileout = 'elf.cube',
/
`;
  }

  // --- Work function ---
  if (postProcessing.includes('work_function')) {
    files['pp_potential.in'] = `&INPUTPP
  prefix = 'calc',
  outdir = './tmp',
  filplot = 'potential.dat',
  plot_num = 11,
/

&PLOT
  iflag = 3,
  output_format = 6,
  fileout = 'potential.cube',
/
`;
    files['avg_potential.in'] = `1
potential.dat
1.D0
1000
3
3.000000
`;
  }

  // --- STM ---
  if (postProcessing.includes('stm')) {
    files['pp_stm.in'] = `&INPUTPP
  prefix = 'calc',
  outdir = './tmp',
  filplot = 'stm.dat',
  plot_num = 5,
  sample_bias = -0.5,
/

&PLOT
  iflag = 2,
  output_format = 6,
  fileout = 'stm.cube',
/
`;
  }

  // --- Dielectric / optical ---
  if (postProcessing.includes('dielectric') || postProcessing.includes('optical_absorption') || postProcessing.includes('refractive_index')) {
    files['epsilon.in'] = `&INPUTPP
  prefix = 'calc',
  outdir = './tmp',
  calculation = 'eps',
  smeartype = 'gauss',
  intersmear = 0.136,
  nw = 1000,
  wmax = 30.0,
/
`;
  }

  // --- Phonon ---
  if (postProcessing.includes('phonon') || postProcessing.includes('thermo') || postProcessing.includes('thermal_conductivity') || postProcessing.includes('born') || postProcessing.includes('raman') || postProcessing.includes('ir')) {
    files['ph.in'] = `Phonon calculation
&INPUTPH
  prefix = 'calc',
  outdir = './tmp',
  fildyn = 'dyn',
  tr2_ph = 1.0d-14,
  ldisp = .true.,
  nq1 = ${Math.max(Math.floor(kx / 2), 2)},
  nq2 = ${Math.max(Math.floor(ky / 2), 2)},
  nq3 = ${Math.max(Math.floor(kz / 2), 2)},${postProcessing.includes('born') ? '\n  epsil = .true.,' : ''}${postProcessing.includes('raman') ? '\n  lraman = .true.,\n  epsil = .true.,' : ''}
/
`;

    files['q2r.in'] = `&INPUT
  fildyn = 'dyn',
  zasr = 'simple',
  flfrc = 'force_constants.dat',
/
`;

    files['matdyn_dos.in'] = `&INPUT
  asr = 'simple',
  flfrc = 'force_constants.dat',
  flfrq = 'phonon_freq.dat',
  fldos = 'phonon_dos.dat',
  dos = .true.,
  nk1 = ${kx * 2}, nk2 = ${ky * 2}, nk3 = ${kz * 2},
/
`;

    files['matdyn_disp.in'] = `&INPUT
  asr = 'simple',
  flfrc = 'force_constants.dat',
  flfrq = 'phonon_disp.dat',
  q_in_band_form = .true.,
/
4
  0.0000 0.0000 0.0000  30
  0.5000 0.0000 0.0000  30
  0.5000 0.5000 0.0000  30
  0.0000 0.0000 0.0000  1
`;
  }

  // --- Elastic ---
  if (postProcessing.includes('elastic')) {
    files['elastic_note.txt'] = `# Elastic tensor calculation
# Use ElaStic or thermo_pw package for automated elastic constant calculation.
# Workflow:
# 1. Optimize structure (vc-relax)
# 2. Apply strain deformations (6 independent strains for cubic, more for lower symmetry)
# 3. Calculate stress for each deformation
# 4. Fit stress-strain to obtain Cij
# Recommended tool: ElaStic (https://github.com/exciting/ElaStic)
`;

    // Generate a sample strained SCF for reference
    files['elastic_scf_template.in'] = `&CONTROL
  calculation = 'scf',
  prefix = 'elastic',
  outdir = './tmp',
  pseudo_dir = '${ppPath}',
  tprnfor = .true.,
  tstress = .true.,
/

&SYSTEM
${commonSystem}
/

&ELECTRONS
  conv_thr = 1.0d-10,
  mixing_beta = ${mixingBeta},
/

${formatAtomicSpecies(parsed, ppPath, params.ppLibrary)}

${formatAtomicPositions(parsed)}

${formatCellParameters(parsed)}

K_POINTS automatic
  ${kx} ${ky} ${kz} 0 0 0
`;
  }

  // --- EOS ---
  if (postProcessing.includes('eos')) {
    files['eos_note.txt'] = `# Equation of State calculation
# Workflow:
# 1. Run vc-relax to get equilibrium structure
# 2. Scale lattice parameter by factors (0.96, 0.98, 1.00, 1.02, 1.04)
# 3. Run SCF for each volume
# 4. Fit E(V) curve to Birch-Murnaghan EOS
# 5. Extract B₀ (bulk modulus) and V₀ (equilibrium volume)
`;

    const scales = [0.96, 0.97, 0.98, 0.99, 1.00, 1.01, 1.02, 1.03, 1.04];
    for (const scale of scales) {
      const tag = `eos_${scale.toFixed(2)}`;
      files[`${tag}.in`] = `&CONTROL
  calculation = 'scf',
  prefix = '${tag}',
  outdir = './tmp',
  pseudo_dir = '${ppPath}',
  tprnfor = .true.,
  tstress = .true.,
/

&SYSTEM
${commonSystem}
/

&ELECTRONS
  conv_thr = 1.0d-8,
/

${formatAtomicSpecies(parsed, ppPath, params.ppLibrary)}

${formatAtomicPositions(parsed)}

CELL_PARAMETERS angstrom
  ${(parsed.cellParams.v1[0] * scale).toFixed(10)}  ${(parsed.cellParams.v1[1] * scale).toFixed(10)}  ${(parsed.cellParams.v1[2] * scale).toFixed(10)}
  ${(parsed.cellParams.v2[0] * scale).toFixed(10)}  ${(parsed.cellParams.v2[1] * scale).toFixed(10)}  ${(parsed.cellParams.v2[2] * scale).toFixed(10)}
  ${(parsed.cellParams.v3[0] * scale).toFixed(10)}  ${(parsed.cellParams.v3[1] * scale).toFixed(10)}  ${(parsed.cellParams.v3[2] * scale).toFixed(10)}

K_POINTS automatic
  ${kx} ${ky} ${kz} 0 0 0
`;
    }
  }

  // --- AIMD ---
  if (postProcessing.includes('aimd') || postProcessing.includes('rdf') || postProcessing.includes('msd') || postProcessing.includes('stress_fluct')) {
    files['aimd.in'] = `&CONTROL
  calculation = 'md',
  prefix = 'aimd',
  outdir = './tmp',
  pseudo_dir = '${ppPath}',
  dt = ${Math.round((params.mdTimestep || 1.0) / 0.04837)},
  nstep = ${params.mdSteps || 5000},
  tprnfor = .true.,
  tstress = .true.,
/

&SYSTEM
${commonSystem}
  nosym = .true.,
/

&ELECTRONS
  conv_thr = 1.0d-6,
  mixing_beta = ${mixingBeta},
/

&IONS
  ion_dynamics = 'verlet',
  ion_temperature = 'rescaling',
  tempw = ${params.mdTemp || 300},
  tolp = 50.0,
/

${formatAtomicSpecies(parsed, ppPath, params.ppLibrary)}

${formatAtomicPositions(parsed)}

${formatCellParameters(parsed)}

K_POINTS automatic
  ${Math.max(Math.floor(kx / 2), 1)} ${Math.max(Math.floor(ky / 2), 1)} ${Math.max(Math.floor(kz / 2), 1)} 0 0 0
`;
  }

  // --- NEB ---
  if (postProcessing.includes('neb')) {
    files['neb.in'] = `BEGIN
BEGIN_PATH_INPUT
&PATH
  restart_mode = 'from_scratch',
  string_method = 'neb',
  nstep_path = 100,
  num_of_images = 7,
  opt_scheme = 'broyden',
  CI_scheme = 'auto',
/
END_PATH_INPUT

BEGIN_ENGINE_INPUT
&CONTROL
  prefix = 'neb',
  outdir = './tmp',
  pseudo_dir = '${ppPath}',
/

&SYSTEM
${commonSystem}
/

&ELECTRONS
  conv_thr = 1.0d-8,
/

BEGIN_POSITIONS
FIRST_IMAGE
ATOMIC_POSITIONS crystal
! === Insert initial positions here ===
${parsed.atoms.map(a => `  ${a.symbol.padEnd(4)} ${a.x.toFixed(10)}  ${a.y.toFixed(10)}  ${a.z.toFixed(10)}`).join('\n')}
LAST_IMAGE
ATOMIC_POSITIONS crystal
! === Insert final positions here (modify atoms that move) ===
${parsed.atoms.map(a => `  ${a.symbol.padEnd(4)} ${a.x.toFixed(10)}  ${a.y.toFixed(10)}  ${a.z.toFixed(10)}`).join('\n')}
END_POSITIONS

${formatCellParameters(parsed)}

K_POINTS automatic
  ${kx} ${ky} ${kz} 0 0 0

${formatAtomicSpecies(parsed, ppPath, params.ppLibrary)}

END_ENGINE_INPUT
END
`;
  }

  // --- Wannier90 ---
  if (postProcessing.includes('wannier')) {
    files['wannier90.win'] = `num_wann = ${parsed.nat * 4}
num_iter = 200
dis_num_iter = 200

begin unit_cell_cart
Ang
${parsed.cellParams.v1.map(v => v.toFixed(10)).join('  ')}
${parsed.cellParams.v2.map(v => v.toFixed(10)).join('  ')}
${parsed.cellParams.v3.map(v => v.toFixed(10)).join('  ')}
end unit_cell_cart

begin atoms_frac
${parsed.atoms.map(a => `${a.symbol}  ${a.x.toFixed(10)}  ${a.y.toFixed(10)}  ${a.z.toFixed(10)}`).join('\n')}
end atoms_frac

begin kpoint_path
G  0.000  0.000  0.000  X  0.500  0.000  0.000
X  0.500  0.000  0.000  M  0.500  0.500  0.000
M  0.500  0.500  0.000  G  0.000  0.000  0.000
end kpoint_path

mp_grid = ${kx} ${ky} ${kz}

begin kpoints
! k-points will be generated by pw2wannier90
end kpoints

bands_plot = .true.
wannier_plot = .true.
`;

    files['pw2wan.in'] = `&INPUTPP
  prefix = 'calc',
  outdir = './tmp',
  seedname = 'wannier90',
  write_mmn = .true.,
  write_amn = .true.,
  write_unk = .false.,
/
`;
  }

  // --- BoltzTraP ---
  if (postProcessing.includes('boltztrap')) {
    files['boltztrap_note.txt'] = `# BoltzTraP thermoelectric calculation
# Workflow:
# 1. Run SCF → NSCF with dense k-grid
# 2. Use pw2boltztrap or BoltzTraP2 interface
# 3. Outputs: Seebeck coefficient, electrical conductivity, power factor, ZT
# Recommended: BoltzTraP2 (pip install BoltzTraP2)
# Command: btp2 -vv interpolate .
`;
  }

  // --- XPS ---
  if (postProcessing.includes('xps')) {
    files['xps_note.txt'] = `# XPS core-level calculation
# Workflow:
# 1. Run SCF for ground state
# 2. Use core-hole approach: remove one core electron from target atom
# 3. Compare total energies: binding energy = E(core-hole) - E(ground-state)
# 4. For initial state approximation: use Kohn-Sham eigenvalues
`;
  }

  // --- MAE ---
  if (postProcessing.includes('mae')) {
    files['mae_soc_x.in'] = `&CONTROL
  calculation = 'scf',
  prefix = 'mae',
  outdir = './tmp',
  pseudo_dir = '${ppPath}',
/

&SYSTEM
${commonSystem}
  noncolin = .true.,
  lspinorb = .true.,
  angle1(1) = 90.0,
  angle2(1) = 0.0,
/

&ELECTRONS
  conv_thr = 1.0d-10,
/

${formatAtomicSpecies(parsed, ppPath, params.ppLibrary)}

${formatAtomicPositions(parsed)}

${formatCellParameters(parsed)}

K_POINTS automatic
  ${kx * 2} ${ky * 2} ${kz * 2} 0 0 0
`;

    files['mae_soc_z.in'] = files['mae_soc_x.in']
      .replace("prefix = 'mae'", "prefix = 'mae_z'")
      .replace('angle1(1) = 90.0', 'angle1(1) = 0.0');
  }

  // --- Piezoelectric ---
  if (postProcessing.includes('piezoelectric')) {
    files['piezoelectric_note.txt'] = `# Piezoelectric coefficient calculation
# Use DFPT (ph.x with epsil=.true.) or finite field method
# Workflow:
# 1. Optimize structure
# 2. Run ph.x with epsil=.true. to get Born effective charges and dielectric tensor
# 3. Use the Berry phase approach for proper piezoelectric constants
`;
  }

  // --- Generate run script ---
  files['run_all.sh'] = generateRunScript(needSCF, needNSCF, needRelax, needVCRelax, needBands, postProcessing, params);

  return files;
}

function generateRunScript(needSCF, needNSCF, needRelax, needVCRelax, needBands, postProcessing, params) {
  const np = params.nprocs || 4;
  const mpi = `mpirun -np ${np}`;

  let script = `#!/bin/bash
# ==============================================
# DFT Workflow Script - Quantum ESPRESSO
# Auto-generated
# ==============================================

set -e
NPROCS=${np}
MPI="mpirun -np \${NPROCS}"

echo "=========================================="
echo "Starting DFT workflow"
echo "=========================================="
`;

  if (needVCRelax) {
    script += `
echo "[1] Running VC-Relax..."
\${MPI} pw.x < vc-relax.in > vc-relax.out
echo "    VC-Relax completed."
`;
  }

  if (needRelax) {
    script += `
echo "[2] Running Relaxation..."
\${MPI} pw.x < relax.in > relax.out
echo "    Relaxation completed."
`;
  }

  if (needSCF) {
    script += `
echo "[3] Running SCF..."
\${MPI} pw.x < scf.in > scf.out
echo "    SCF completed."
`;
  }

  if (needNSCF) {
    script += `
echo "[4] Running NSCF..."
\${MPI} pw.x < nscf.in > nscf.out
echo "    NSCF completed."
`;
  }

  if (needBands) {
    script += `
echo "[5] Running Bands calculation..."
\${MPI} pw.x < bands.in > bands.out
bands.x < bands_pp.in > bands_pp.out
echo "    Bands completed."
`;
  }

  if (postProcessing.includes('dos')) {
    script += `
echo "[PP] Computing DOS..."
dos.x < dos.in > dos.out
echo "    DOS completed."
`;
  }

  if (postProcessing.includes('pdos')) {
    script += `
echo "[PP] Computing PDOS..."
projwfc.x < pdos.in > pdos.out
echo "    PDOS completed."
`;
  }

  if (postProcessing.includes('bader')) {
    script += `
echo "[PP] Computing charge density for Bader analysis..."
pp.x < pp_charge.in > pp_charge.out
echo "    Run 'bader charge.cube' for Bader analysis."
`;
  }

  if (postProcessing.includes('cdd')) {
    script += `
echo "[PP] Computing charge density difference..."
pp.x < pp_cdd.in > pp_cdd.out
echo "    CDD completed."
`;
  }

  if (postProcessing.includes('elf')) {
    script += `
echo "[PP] Computing ELF..."
pp.x < pp_elf.in > pp_elf.out
echo "    ELF completed."
`;
  }

  if (postProcessing.includes('work_function')) {
    script += `
echo "[PP] Computing work function..."
pp.x < pp_potential.in > pp_potential.out
average.x < avg_potential.in > avg_potential.out
echo "    Work function completed."
`;
  }

  if (postProcessing.includes('stm')) {
    script += `
echo "[PP] Computing STM image..."
pp.x < pp_stm.in > pp_stm.out
echo "    STM completed."
`;
  }

  if (postProcessing.includes('dielectric') || postProcessing.includes('optical_absorption') || postProcessing.includes('refractive_index')) {
    script += `
echo "[PP] Computing dielectric function..."
epsilon.x < epsilon.in > epsilon.out
echo "    Dielectric/optical completed."
`;
  }

  if (postProcessing.includes('phonon') || postProcessing.includes('thermo') || postProcessing.includes('thermal_conductivity') || postProcessing.includes('born') || postProcessing.includes('raman') || postProcessing.includes('ir')) {
    script += `
echo "[PP] Running phonon calculation..."
\${MPI} ph.x < ph.in > ph.out
q2r.x < q2r.in > q2r.out
matdyn.x < matdyn_dos.in > matdyn_dos.out
matdyn.x < matdyn_disp.in > matdyn_disp.out
echo "    Phonon completed."
`;
  }

  if (postProcessing.includes('eos')) {
    script += `
echo "[PP] Running EOS calculations..."
for scale in 0.96 0.97 0.98 0.99 1.00 1.01 1.02 1.03 1.04; do
  echo "  Volume scale: \${scale}"
  \${MPI} pw.x < eos_\${scale}.in > eos_\${scale}.out
done
echo "    EOS completed. Fit E(V) data with ev.x or python script."
`;
  }

  if (postProcessing.includes('aimd') || postProcessing.includes('rdf') || postProcessing.includes('msd') || postProcessing.includes('stress_fluct')) {
    script += `
echo "[PP] Running AIMD..."
\${MPI} pw.x < aimd.in > aimd.out
echo "    AIMD completed."
`;
  }

  if (postProcessing.includes('wannier')) {
    script += `
echo "[PP] Running Wannier90..."
wannier90.x -pp wannier90
\${MPI} pw2wannier90.x < pw2wan.in > pw2wan.out
wannier90.x wannier90
echo "    Wannier90 completed."
`;
  }

  if (postProcessing.includes('neb')) {
    script += `
echo "[PP] Running NEB..."
\${MPI} neb.x -inp neb.in > neb.out
echo "    NEB completed."
`;
  }

  if (postProcessing.includes('elastic')) {
    script += `
echo "[PP] Elastic tensor:"
echo "    Use ElaStic package for automated elastic constant calculation."
echo "    See elastic_note.txt for instructions."
`;
  }

  if (postProcessing.includes('mae')) {
    script += `
echo "[PP] Running MAE (SOC) calculations..."
\${MPI} pw.x < mae_soc_x.in > mae_soc_x.out
\${MPI} pw.x < mae_soc_z.in > mae_soc_z.out
echo "    MAE = E(x) - E(z). Extract total energies from outputs."
`;
  }

  script += `
echo "=========================================="
echo "All calculations completed!"
echo "=========================================="
`;

  return script;
}
