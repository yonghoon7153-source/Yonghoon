"""
Adhesion v6: 1L NCM + crystalline SE + NCM-frozen SE annealing
Protocol:
  1. xy-shift SE
  2. Stack on NCM (gap 2.5A)
  3. Fix ALL NCM → 500K SE anneal (5ps) → 100K quench (2ps)
  4. Unfreeze NCM → LBFGS relax
  5. Separation + single point → Wad

Quick test: 5 seeds × 5 compositions
Compare with v5 (no anneal) to see if scatter reduces.
"""
import numpy as np
from pymatgen.core import Structure, Lattice
from pymatgen.io.ase import AseAtomsAdaptor
from ase import Atoms, units
from ase.constraints import FixAtoms
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.optimize import LBFGS
from fairchem.core import pretrained_mlip
from fairchem.core.calculate.ase_calculator import FAIRChemCalculator
import json, time

predictor = pretrained_mlip.get_predict_unit("uma-s-1p1", device="cuda")
adaptor = AseAtomsAdaptor()

def new_calc():
    return FAIRChemCalculator(predictor, task_name="omat")

# ═══════════════════════════════════════════════════════
# CONFIG — UPDATE THESE FOR YOUR SYSTEM!
# ═══════════════════════════════════════════════════════
SE_DIR = "/home/ubuntu/runyourai/bml/manuscript_support/structures"

COMPS = {
    'comp3':  {'cif': f'{SE_DIR}/comp3_V0.cif',  'repeat': [2,2,1], 'ncm_nx': 5},
    'comp4':  {'cif': f'{SE_DIR}/comp4_V0.cif',  'repeat': [2,2,1], 'ncm_nx': 5},
    'comp5':  {'cif': f'{SE_DIR}/comp5_V0.cif',  'repeat': [2,2,1], 'ncm_nx': 5},
    'comp1':  {'cif': f'{SE_DIR}/comp1_V0.cif',  'repeat': [2,2,3], 'ncm_nx': 7},
    'comp2B': {'cif': f'{SE_DIR}/comp2B_V0.cif', 'repeat': [2,2,3], 'ncm_nx': 7},
}

seeds = list(range(42, 47))  # 5 seeds for quick test
ANNEAL_T = 500    # K
ANNEAL_PS = 5     # ps
QUENCH_PS = 2     # ps
GAP = 2.5         # Angstrom
VACUUM = 30.0     # Angstrom

# ═══════════════════════════════════════════════════════
# NCM builder: LiNiO2 R-3m hexagonal
# ═══════════════════════════════════════════════════════
def build_ncm_1L(nx):
    a, c = 2.878, 14.19
    lat = Lattice.hexagonal(a, c)
    unit = Structure(lat,
        ["Li", "Ni", "O", "O"],
        [[0,0,0.5], [0,0,0], [0,0,0.2584], [0,0,0.7416]])
    unit.make_supercell([nx, nx, 1])
    atoms = adaptor.get_atoms(unit)
    # Center NCM: shift so z_min = 0
    pos = atoms.get_positions()
    pos[:, 2] -= pos[:, 2].min()
    atoms.set_positions(pos)
    return atoms

# ═══════════════════════════════════════════════════════
# SE loader
# ═══════════════════════════════════════════════════════
def load_se(comp_id):
    cfg = COMPS[comp_id]
    se_struct = Structure.from_file(cfg['cif'])
    se_struct.make_supercell(cfg['repeat'])
    atoms = adaptor.get_atoms(se_struct)
    return atoms

# ═══════════════════════════════════════════════════════
# Interface builder: stack SE on NCM with xy-shift
# ═══════════════════════════════════════════════════════
def build_interface(ncm, se, dx, dy):
    """
    Strain SE to NCM cell (xy), place SE above NCM with gap.
    Returns: (interface_atoms, n_ncm, area_A2)
    """
    ncm_cell = ncm.cell.array.copy()
    se_cart = se.get_positions().copy()

    # Convert SE Cartesian → fractional in NCM cell (applies xy strain)
    ncm_inv = np.linalg.inv(ncm_cell)
    se_frac = se_cart @ ncm_inv

    # xy random shift (PBC-safe)
    se_frac[:, 0] = (se_frac[:, 0] + dx) % 1.0
    se_frac[:, 1] = (se_frac[:, 1] + dy) % 1.0

    # Back to Cartesian in NCM cell
    se_pos = se_frac @ ncm_cell

    # Place SE above NCM
    ncm_pos = ncm.get_positions().copy()
    ncm_zmax = ncm_pos[:, 2].max()
    se_zmin = se_pos[:, 2].min()
    se_pos[:, 2] += (ncm_zmax + GAP - se_zmin)
    se_zmax = se_pos[:, 2].max()

    # Combined cell: NCM xy, z = total + vacuum
    total_z = se_zmax + VACUUM
    combined_cell = ncm_cell.copy()
    combined_cell[2] = [0, 0, total_z]

    # Combine atoms
    symbols = ncm.get_chemical_symbols() + se.get_chemical_symbols()
    positions = np.vstack([ncm_pos, se_pos])
    interface = Atoms(symbols=symbols, positions=positions,
                      cell=combined_cell, pbc=True)

    # Area
    area = np.linalg.norm(np.cross(ncm_cell[0], ncm_cell[1]))

    return interface, len(ncm), area

# ═══════════════════════════════════════════════════════
# v6 protocol: freeze NCM → anneal SE → unfreeze → relax
# ═══════════════════════════════════════════════════════
def anneal_se(interface, n_ncm):
    """Anneal SE surface while NCM is frozen"""
    atoms = interface.copy()

    # Fix ALL NCM atoms
    atoms.set_constraint(FixAtoms(indices=list(range(n_ncm))))
    atoms.calc = new_calc()

    # Set velocities (zero for fixed NCM)
    MaxwellBoltzmannDistribution(atoms, temperature_K=ANNEAL_T)
    vel = atoms.get_velocities()
    vel[:n_ncm] = 0.0
    atoms.set_velocities(vel)

    # 500K anneal — only SE moves
    dyn = Langevin(atoms, 1*units.fs, temperature_K=ANNEAL_T, friction=0.01)
    dyn.run(int(ANNEAL_PS * 1000))

    # Quench to 100K
    dyn2 = Langevin(atoms, 1*units.fs, temperature_K=100, friction=0.05)
    dyn2.run(int(QUENCH_PS * 1000))

    # Remove constraints
    atoms.set_constraint()
    return atoms

# ═══════════════════════════════════════════════════════
# Wad calculation (separation method)
# ═══════════════════════════════════════════════════════
def calc_wad(interface, n_ncm, area):
    """Wad = (E_sep - E_int) / A × 16.0218"""
    # E_int
    atoms_int = interface.copy()
    atoms_int.calc = new_calc()
    E_int = atoms_int.get_potential_energy()

    # E_sep: move SE +30A, expand cell z +30A
    atoms_sep = interface.copy()
    pos = atoms_sep.get_positions()
    pos[n_ncm:, 2] += 30.0
    atoms_sep.set_positions(pos)
    cell = atoms_sep.cell.array.copy()
    cell[2, 2] += 30.0
    atoms_sep.set_cell(cell)
    atoms_sep.calc = new_calc()
    E_sep = atoms_sep.get_potential_energy()

    Wad = (E_sep - E_int) / area * 16.0218
    return Wad, E_int, E_sep

# ═══════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════
print("="*60)
print("Adhesion v6: 1L NCM + NCM-frozen SE annealing")
print(f"Anneal: {ANNEAL_T}K {ANNEAL_PS}ps → quench 100K {QUENCH_PS}ps")
print(f"Seeds: {seeds}")
print("="*60)

# v5 reference (no anneal) for comparison
v5_ref = {
    'comp3': 2.103, 'comp4': 1.970, 'comp5': 1.651,
    'comp1': 1.277, 'comp2B': 1.183
}

all_results = {}
ncm_cache = {}

for comp_id in ['comp3', 'comp4', 'comp5', 'comp1', 'comp2B']:
    cfg = COMPS[comp_id]
    nx = cfg['ncm_nx']

    # Build/cache NCM
    if nx not in ncm_cache:
        ncm_cache[nx] = build_ncm_1L(nx)
        print(f"\nNCM {nx}x{nx}x1: {len(ncm_cache[nx])} atoms")
    ncm = ncm_cache[nx]

    # Load SE
    try:
        se = load_se(comp_id)
    except Exception as e:
        print(f"\n=== {comp_id}: SE load FAILED ({e}) ===")
        print(f"    Check path: {cfg['cif']}")
        continue

    print(f"\n=== {comp_id} (SE {cfg['repeat']}, {len(se)} at) ===")

    wads = []
    rng = np.random.RandomState(42)
    xy_shifts = [(rng.random(), rng.random()) for _ in range(max(seeds) - 41)]

    for s in seeds:
        dx, dy = xy_shifts[s - 42]
        t0 = time.time()

        # Build interface
        interface, n_ncm, area = build_interface(ncm, se, dx, dy)

        # v6: anneal SE
        interface = anneal_se(interface, n_ncm)

        # LBFGS relax (full system, no constraints)
        interface.calc = new_calc()
        try:
            LBFGS(interface, logfile=None).run(fmax=0.01, steps=200)
        except:
            pass

        # Wad
        Wad, E_int, E_sep = calc_wad(interface, n_ncm, area)
        dt = time.time() - t0

        wads.append(Wad)
        print(f"    seed={s}: Wad={Wad:.3f} ({dx:.2f},{dy:.2f}) [{dt:.0f}s]", flush=True)

    mean_w = np.mean(wads)
    std_w = np.std(wads)
    print(f"  v6 mean={mean_w:.3f}+/-{std_w:.3f} (n={len(wads)})")
    print(f"  v5 ref ={v5_ref.get(comp_id, 'N/A')}")

    all_results[comp_id] = {
        'wads': wads, 'mean': mean_w, 'std': std_w,
        'v5_ref': v5_ref.get(comp_id),
        'seeds': seeds
    }

# ═══════════════════════════════════════════════════════
# SUMMARY
# ═══════════════════════════════════════════════════════
print(f"\n{'='*60}")
print("SUMMARY: v6 (anneal) vs v5 (no anneal)")
print(f"{'='*60}")
print(f"{'Comp':<8} {'v6 mean':>10} {'v6 std':>8} {'v5 ref':>8} {'delta':>8}")
for comp_id in ['comp3', 'comp4', 'comp5', 'comp1', 'comp2B']:
    if comp_id in all_results:
        r = all_results[comp_id]
        v5 = r['v5_ref'] or 0
        delta = r['mean'] - v5
        print(f"{comp_id:<8} {r['mean']:>10.3f} {r['std']:>8.3f} {v5:>8.3f} {delta:>+8.3f}")

# Check trends
comps = ['comp3', 'comp4', 'comp5', 'comp1', 'comp2B']
avail = [c for c in comps if c in all_results]
if len(avail) == 5:
    m = {c: all_results[c]['mean'] for c in avail}
    li54_ok = m['comp3'] > m['comp4'] > m['comp5']
    li6_ok = m['comp1'] > m['comp2B']
    cross_ok = m['comp3'] > m['comp1']
    print(f"\nWithin Li5.4: comp3>comp4>comp5 = {'YES' if li54_ok else 'NO'}")
    print(f"Within Li6:   comp1>comp2B     = {'YES' if li6_ok else 'NO'}")
    print(f"Cross-family: comp3>comp1      = {'YES' if cross_ok else 'NO'}")

with open('adhesion_v6_1L_test.json', 'w') as f:
    json.dump(all_results, f, indent=2, default=lambda x: x.tolist() if hasattr(x, 'tolist') else x)

print(f"\nResults saved to adhesion_v6_1L_test.json")
