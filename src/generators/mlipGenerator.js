// Generator for UMA and MACE MLIP scripts (Python-based)

export function generateUMAScripts(parsed, params, calcTypes, postProcessing) {
  return generateMLIPScripts(parsed, params, calcTypes, postProcessing, 'uma');
}

export function generateMACEScripts(parsed, params, calcTypes, postProcessing) {
  return generateMLIPScripts(parsed, params, calcTypes, postProcessing, 'mace');
}

function generateMLIPScripts(parsed, params, calcTypes, postProcessing, code) {
  const files = {};
  const modelPath = params.mlipModel || (code === 'mace' ? 'mace_mp_medium.model' : 'uma_model.pt');
  const mdTemp = params.mdTemp || 300;
  const mdSteps = params.mdSteps || 5000;
  const mdTimestep = params.mdTimestep || 1.0;

  // CIF content for ASE to read
  files['structure.cif'] = '# Paste your CIF file content here or use the uploaded CIF\n' +
    `# Structure: ${parsed.elements.join('-')}, ${parsed.nat} atoms\n`;

  const { v1, v2, v3 } = parsed.cellParams;
  const atomLines = parsed.atoms.map(a =>
    `    ['${a.symbol}', [${a.x.toFixed(8)}, ${a.y.toFixed(8)}, ${a.z.toFixed(8)}]]`
  ).join(',\n');

  const setupCode = code === 'mace'
    ? `from mace.calculators import MACECalculator
calc = MACECalculator(model_paths='${modelPath}', device='cuda')`
    : `from uma.calculator import UMACalculator
calc = UMACalculator(model='${modelPath}', device='cuda')`;

  const structureSetup = `import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.build import bulk

# Read structure from CIF
atoms = read('structure.cif')

# Or build manually:
# cell = [
#     [${v1.map(v => v.toFixed(6)).join(', ')}],
#     [${v2.map(v => v.toFixed(6)).join(', ')}],
#     [${v3.map(v => v.toFixed(6)).join(', ')}],
# ]
# atoms = Atoms(
#     symbols=[${parsed.atoms.map(a => `'${a.symbol}'`).join(', ')}],
#     scaled_positions=[
#         ${parsed.atoms.map(a => `[${a.x.toFixed(8)}, ${a.y.toFixed(8)}, ${a.z.toFixed(8)}]`).join(',\n#         ')}
#     ],
#     cell=cell,
#     pbc=True,
# )

${setupCode}
atoms.calc = calc`;

  // --- Single point / SCF ---
  if (calcTypes.includes('scf') || postProcessing.length > 0) {
    files['single_point.py'] = `"""Single-point energy calculation with ${code.toUpperCase()}"""
${structureSetup}

energy = atoms.get_potential_energy()
forces = atoms.get_forces()
stress = atoms.get_stress()

print(f"Total energy: {energy:.6f} eV")
print(f"Max force: {np.max(np.abs(forces)):.6f} eV/Å")
print(f"Stress (Voigt): {stress}")
`;
  }

  // --- Relaxation ---
  if (calcTypes.includes('relax') || calcTypes.includes('vc-relax')) {
    const vcRelax = calcTypes.includes('vc-relax');
    files['optimize.py'] = `"""Geometry optimization with ${code.toUpperCase()}"""
${structureSetup}

from ase.optimize import BFGS${vcRelax ? '\nfrom ase.constraints import ExpCellFilter' : ''}

${vcRelax ? 'ecf = ExpCellFilter(atoms)\nopt = BFGS(ecf, trajectory="opt.traj", logfile="opt.log")' : 'opt = BFGS(atoms, trajectory="opt.traj", logfile="opt.log")'}
opt.run(fmax=0.01)

print(f"Optimized energy: {atoms.get_potential_energy():.6f} eV")
write('optimized.cif', atoms)
print("Optimized structure saved to optimized.cif")
`;
  }

  // --- EOS ---
  if (postProcessing.includes('eos')) {
    files['eos.py'] = `"""Equation of State calculation with ${code.toUpperCase()}"""
${structureSetup}

from ase.eos import EquationOfState

cell0 = atoms.get_cell().copy()
volumes = []
energies = []

for scale in np.linspace(0.94, 1.06, 13):
    a = atoms.copy()
    a.set_cell(cell0 * scale**(1/3), scale_atoms=True)
    a.calc = calc
    e = a.get_potential_energy()
    v = a.get_volume()
    volumes.append(v)
    energies.append(e)
    print(f"  scale={scale:.3f}, V={v:.3f} Å³, E={e:.6f} eV")

eos = EquationOfState(volumes, energies, eos='birchmurnaghan')
v0, e0, B = eos.fit()
print(f"\\nEquilibrium volume: {v0:.3f} Å³")
print(f"Bulk modulus: {B / 1.602176634e-19 * 1e-9:.2f} GPa")
eos.plot('eos_fit.png')
print("Plot saved to eos_fit.png")
`;
  }

  // --- Phonon ---
  if (postProcessing.includes('phonon') || postProcessing.includes('thermo') || postProcessing.includes('thermal_conductivity')) {
    files['phonon.py'] = `"""Phonon calculation with ${code.toUpperCase()} + phonopy"""
${structureSetup}

from phonopy import Phonopy
from phonopy.interface.ase import get_phonopy_structure

ph_atoms = get_phonopy_structure(atoms)
phonon = Phonopy(ph_atoms, supercell_matrix=[[2, 0, 0], [0, 2, 0], [0, 0, 2]])
phonon.generate_displacements(distance=0.01)

supercells = phonon.get_supercells_with_displacements()
forces_sets = []

print(f"Computing forces for {len(supercells)} displacements...")
for i, sc in enumerate(supercells):
    from ase import Atoms as AseAtoms
    a = AseAtoms(
        symbols=sc.symbols,
        positions=sc.positions,
        cell=sc.cell,
        pbc=True,
    )
    a.calc = calc
    f = a.get_forces()
    forces_sets.append(f)
    print(f"  Displacement {i+1}/{len(supercells)} done")

phonon.forces = forces_sets
phonon.produce_force_constants()

# Phonon band structure
path = [[[0, 0, 0], [0.5, 0, 0], [0.5, 0.5, 0], [0, 0, 0]]]
labels = ['Γ', 'X', 'M', 'Γ']
phonon.set_band_structure(path, labels=labels)
phonon.plot_band_structure().savefig('phonon_bands.png')
print("Phonon band structure saved to phonon_bands.png")

# Phonon DOS
phonon.set_mesh([20, 20, 20])
phonon.set_total_DOS()
phonon.plot_total_DOS().savefig('phonon_dos.png')
print("Phonon DOS saved to phonon_dos.png")

# Thermodynamic properties
phonon.set_thermal_properties(t_min=0, t_max=1000, t_step=10)
tp = phonon.get_thermal_properties_dict()
print("\\nThermodynamic properties at 300K:")
idx = list(tp['temperatures']).index(300.0) if 300.0 in tp['temperatures'] else 30
print(f"  Cv = {tp['heat_capacity'][idx]:.3f} J/mol/K")
print(f"  Entropy = {tp['entropy'][idx]:.3f} J/mol/K")
print(f"  Free energy = {tp['free_energy'][idx]:.3f} kJ/mol")
`;
  }

  // --- MD ---
  if (postProcessing.includes('aimd') || postProcessing.includes('rdf') || postProcessing.includes('msd') || postProcessing.includes('stress_fluct')) {
    files['md.py'] = `"""Molecular dynamics with ${code.toUpperCase()}"""
${structureSetup}

from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase import units

# Set up MD
MaxwellBoltzmannDistribution(atoms, temperature_K=${mdTemp})

dyn = Langevin(
    atoms,
    timestep=${mdTimestep} * units.fs,
    temperature_K=${mdTemp},
    friction=0.01,
    trajectory='md.traj',
    logfile='md.log',
)

print(f"Running MD at ${mdTemp} K for ${mdSteps} steps...")
dyn.run(${mdSteps})
print("MD completed. Trajectory saved to md.traj")
`;
  }

  // --- RDF ---
  if (postProcessing.includes('rdf')) {
    files['rdf.py'] = `"""Radial Distribution Function from MD trajectory"""
from ase.io import read
import numpy as np

traj = read('md.traj', index=':')

from ase.geometry.analysis import Analysis
ana = Analysis(traj)

# All-all RDF
rdf = ana.get_rdf(rmax=8.0, nbins=200, return_dists=True)
r, g = rdf[0]

import matplotlib.pyplot as plt
plt.figure(figsize=(8, 5))
plt.plot(r, g)
plt.xlabel('r (Å)')
plt.ylabel('g(r)')
plt.title('Radial Distribution Function')
plt.savefig('rdf.png', dpi=150)
print("RDF saved to rdf.png")
`;
  }

  // --- MSD ---
  if (postProcessing.includes('msd')) {
    files['msd.py'] = `"""Mean Square Displacement & Diffusion Coefficient"""
from ase.io import read
import numpy as np

traj = read('md.traj', index=':')
positions = np.array([a.get_positions() for a in traj])

# Compute MSD
ref = positions[0]
msd = np.mean(np.sum((positions - ref)**2, axis=-1), axis=-1)
time = np.arange(len(msd)) * ${mdTimestep}  # fs

# Fit linear region for diffusion coefficient
from numpy.polynomial import polynomial as P
# Use last 60% of data
start = int(0.4 * len(time))
coeff = P.polyfit(time[start:] * 1e-15, msd[start:] * 1e-20, 1)
D = coeff[1] / 6  # m²/s
print(f"Diffusion coefficient: {D:.4e} m²/s")

# Ionic conductivity (Nernst-Einstein)
from ase import units
T = ${mdTemp}
n_ions = len(traj[0])
V = traj[0].get_volume() * 1e-30  # m³
z = 1  # charge number
e = 1.602e-19
kB = 1.381e-23
sigma = n_ions * z**2 * e**2 * D / (V * kB * T)
print(f"Ionic conductivity: {sigma:.4e} S/m")

import matplotlib.pyplot as plt
plt.figure(figsize=(8, 5))
plt.plot(time, msd)
plt.xlabel('Time (fs)')
plt.ylabel('MSD (Å²)')
plt.title('Mean Square Displacement')
plt.savefig('msd.png', dpi=150)
print("MSD plot saved to msd.png")
`;
  }

  // --- Elastic ---
  if (postProcessing.includes('elastic')) {
    files['elastic.py'] = `"""Elastic constants with ${code.toUpperCase()}"""
${structureSetup}

from ase.optimize import BFGS
from ase.constraints import ExpCellFilter

# First optimize
ecf = ExpCellFilter(atoms)
opt = BFGS(ecf)
opt.run(fmax=0.005)

# Compute elastic constants via stress-strain
from ase.units import GPa

strains = np.linspace(-0.02, 0.02, 5)
cell0 = atoms.get_cell().copy()
C = np.zeros((6, 6))

voigt_map = [(0,0), (1,1), (2,2), (1,2), (0,2), (0,1)]

for i, (a, b) in enumerate(voigt_map):
    stresses = []
    for eps in strains:
        a_copy = atoms.copy()
        strain_matrix = np.eye(3)
        strain_matrix[a, b] += eps / 2
        strain_matrix[b, a] += eps / 2
        new_cell = cell0 @ strain_matrix
        a_copy.set_cell(new_cell, scale_atoms=True)
        a_copy.calc = calc
        s = a_copy.get_stress(voigt=True)
        stresses.append(s)

    stresses = np.array(stresses)
    for j in range(6):
        C[i, j] = -np.polyfit(strains, stresses[:, j], 1)[0] / GPa

print("Elastic tensor Cij (GPa):")
print(np.array2string(C, precision=2, suppress_small=True))

# Voigt averages
B_V = (C[0,0] + C[1,1] + C[2,2] + 2*(C[0,1] + C[0,2] + C[1,2])) / 9
G_V = ((C[0,0] + C[1,1] + C[2,2]) - (C[0,1] + C[0,2] + C[1,2]) + 3*(C[3,3] + C[4,4] + C[5,5])) / 15
E = 9 * B_V * G_V / (3 * B_V + G_V)
nu = (3 * B_V - 2 * G_V) / (2 * (3 * B_V + G_V))

print(f"\\nBulk modulus (Voigt): {B_V:.2f} GPa")
print(f"Shear modulus (Voigt): {G_V:.2f} GPa")
print(f"Young's modulus: {E:.2f} GPa")
print(f"Poisson ratio: {nu:.4f}")
`;
  }

  // --- NEB ---
  if (postProcessing.includes('neb')) {
    files['neb.py'] = `"""NEB calculation with ${code.toUpperCase()}"""
${structureSetup}

from ase.neb import NEB
from ase.optimize import BFGS
from ase.io import read, write

# Define initial and final states
initial = atoms.copy()
final = atoms.copy()
# === Modify final state positions for the reaction path ===
# final.positions[atom_index] = [x, y, z]

n_images = 7
images = [initial.copy() for _ in range(n_images + 2)]
images[-1] = final.copy()

for img in images:
    img.calc = calc

neb = NEB(images, climb=True)
neb.interpolate()

opt = BFGS(neb, trajectory='neb.traj')
opt.run(fmax=0.05)

# Extract barrier
energies = [img.get_potential_energy() for img in images]
barrier = max(energies) - energies[0]
print(f"\\nActivation barrier: {barrier:.4f} eV ({barrier * 96.485:.2f} kJ/mol)")

import matplotlib.pyplot as plt
plt.figure(figsize=(8, 5))
plt.plot(range(len(energies)), [e - energies[0] for e in energies], 'o-')
plt.xlabel('Image')
plt.ylabel('Energy (eV)')
plt.title('NEB Energy Profile')
plt.savefig('neb_profile.png', dpi=150)
print("NEB profile saved to neb_profile.png")
`;
  }

  // --- Formation/surface/vacancy energies ---
  if (postProcessing.includes('formation_energy') || postProcessing.includes('surface_energy') || postProcessing.includes('vacancy_energy') || postProcessing.includes('adhesion_energy') || postProcessing.includes('bonding_energy')) {
    files['energy_analysis.py'] = `"""Energy analysis with ${code.toUpperCase()}"""
${structureSetup}

from ase.io import read
from ase.optimize import BFGS

# Optimize structure
opt = BFGS(atoms, trajectory='opt.traj')
opt.run(fmax=0.01)
E_bulk = atoms.get_potential_energy()
n_atoms = len(atoms)
print(f"Bulk energy: {E_bulk:.6f} eV ({E_bulk/n_atoms:.6f} eV/atom)")

# Formation energy: E_form = E_compound - sum(n_i * E_i)
# You need reference energies for each element
print("\\n--- Formation Energy ---")
print("E_formation = E_total - Σ(n_i × μ_i)")
print("Provide reference chemical potentials (μ_i) for each element.")

# Surface energy: gamma = (E_slab - n * E_bulk/N) / (2 * A)
print("\\n--- Surface Energy ---")
print("γ = (E_slab - n × E_bulk_per_atom) / (2 × A)")
print("Create a slab model and compute E_slab.")

# Vacancy energy: E_vac = E(N-1) + E_atom - E(N)
print("\\n--- Vacancy Formation Energy ---")
print("E_vac = E(N-1 atoms) + μ_removed - E(N atoms)")
print("Remove an atom and recompute energy.")
`;
  }

  // --- Bond length analysis ---
  if (postProcessing.includes('bond_length')) {
    files['bond_analysis.py'] = `"""Bond length analysis"""
from ase.io import read
from ase.geometry.analysis import Analysis

atoms = read('structure.cif')

ana = Analysis(atoms)
print("Bond lengths:")
for i in range(len(atoms)):
    for j in range(i+1, len(atoms)):
        d = atoms.get_distance(i, j, mic=True)
        if d < 3.5:  # cutoff
            print(f"  {atoms[i].symbol}({i}) - {atoms[j].symbol}({j}): {d:.4f} Å")
`;
  }

  // --- Stress fluctuation elastic constants ---
  if (postProcessing.includes('stress_fluct')) {
    files['stress_fluctuation.py'] = `"""Finite-temperature elastic constants via stress fluctuation method"""
from ase.io import read
import numpy as np

traj = read('md.traj', index='500:')  # skip equilibration
V = np.mean([a.get_volume() for a in traj])
T = ${mdTemp}
kB = 8.617e-5  # eV/K

stresses = np.array([a.get_stress(voigt=True) for a in traj])  # 6-component Voigt
# C_ijkl = V/(kB*T) * (<sigma_ij * sigma_kl> - <sigma_ij><sigma_kl>)

C = np.zeros((6, 6))
for i in range(6):
    for j in range(6):
        C[i, j] = V / (kB * T) * (
            np.mean(stresses[:, i] * stresses[:, j]) -
            np.mean(stresses[:, i]) * np.mean(stresses[:, j])
        )

print("Elastic tensor from stress fluctuations (eV/Å³):")
print(np.array2string(C, precision=4))

# Convert to GPa
C_GPa = C * 160.2176634
print("\\nElastic tensor (GPa):")
print(np.array2string(C_GPa, precision=2))
`;
  }

  // --- Run script ---
  files['run_all.sh'] = `#!/bin/bash
# ${code.toUpperCase()} MLIP Workflow Script
set -e

${Object.keys(files).filter(f => f.endsWith('.py')).map(f => {
  const name = f.replace('.py', '');
  return `echo "Running ${name}..."
python ${f} 2>&1 | tee ${name}.log
echo "  ${name} completed."
echo`;
}).join('\n\n')}

echo "All ${code.toUpperCase()} calculations completed!"
`;

  return files;
}
