"""
Model C v2 — Step 3 continuation (ranks 1-4 only)
Step 1-2 done. Rank 0 champion: -258.2051 eV
Top 5 ALL from H0 (S=[1,4]), Li spread = 0.1 meV
"""
import numpy as np
from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor
from ase import units
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.optimize import LBFGS
from fairchem.core import pretrained_mlip
from fairchem.core.calculate.ase_calculator import FAIRChemCalculator
from ase.io import write
import json

predictor = pretrained_mlip.get_predict_unit("uma-s-1p1", device="cuda")
adaptor = AseAtomsAdaptor()

def new_calc():
    return FAIRChemCalculator(predictor, task_name="omat")

# Load reference
ref = Structure.from_file('/scratch/x3430a02/kgy/manuscript_support/comp3_lpsc10b06/configs/comp3_lpsc10b06_config_000.cif')

li_sites = []; p_sites = []; s_framework = []; free_sites = []
for i, site in enumerate(ref):
    sp = str(site.specie)
    fc = site.frac_coords
    if sp == 'Li':
        li_sites.append(fc)
    elif sp == 'P':
        p_sites.append(fc)
    else:
        p_coords = np.array(p_sites) if p_sites else np.zeros((1,3))
        dists = ref.lattice.get_all_distances(fc.reshape(1,-1), p_coords)[0]
        if sp == 'S' and len(p_sites) > 0 and min(dists) < 2.5:
            s_framework.append(fc)
        else:
            free_sites.append(fc)

n_li = len(li_sites)
n_li_select = 27

# Step 2 results (hardcoded from completed run)
rng = np.random.RandomState(42)
Li_configs = [sorted(rng.choice(n_li, n_li_select, replace=False).tolist()) for _ in range(20)]

# Top 5 from Step 2 (all H0, S=[1,4])
top5 = [
    {'h_rank': 0, 'li_trial': 15, 's_indices': [1, 4]},
    {'h_rank': 0, 'li_trial': 0,  's_indices': [1, 4]},
    {'h_rank': 0, 'li_trial': 7,  's_indices': [1, 4]},
    {'h_rank': 0, 'li_trial': 12, 's_indices': [1, 4]},
    {'h_rank': 0, 'li_trial': 16, 's_indices': [1, 4]},
]

# Rank 0 already done: -258.2051
print("Step 3 continuation: ranks 1-4")
print(f"Rank 0 (H0_Li15): already done → -258.2051 eV")
print("="*60)

anneal_results = [{'rank': 0, 'h_rank': 0, 'li_trial': 15, 'e_before': -255.1826, 'e_after': -258.2051}]

for rank in range(1, 5):
    r = top5[rank]
    s_idx = r['s_indices']
    li_select = Li_configs[r['li_trial']]

    species = []; coords = []
    for i in li_select: species.append('Li'); coords.append(li_sites[i])
    for c in p_sites: species.append('P'); coords.append(c)
    for c in s_framework: species.append('S'); coords.append(c)
    for i, c in enumerate(free_sites):
        species.append('S' if i in s_idx else 'Cl'); coords.append(c)

    struct = Structure(ref.lattice, species, coords)
    atoms = adaptor.get_atoms(struct)
    atoms.calc = new_calc()

    e_before = atoms.get_potential_energy()
    MaxwellBoltzmannDistribution(atoms, temperature_K=500)
    Langevin(atoms, 1*units.fs, temperature_K=500, friction=0.01).run(50000)
    Langevin(atoms, 1*units.fs, temperature_K=300, friction=0.05).run(10000)
    try: LBFGS(atoms, logfile=None).run(fmax=0.005, steps=300)
    except: pass
    e_after = atoms.get_potential_energy()

    anneal_results.append({
        'rank': rank, 'h_rank': r['h_rank'], 'li_trial': r['li_trial'],
        'e_before': e_before, 'e_after': e_after
    })
    write(f'modelC_v2_anneal_rank{rank}.xyz', atoms)
    print(f"  Rank {rank} (H0_Li{r['li_trial']}): "
          f"{e_before:.4f} → {e_after:.4f} (dE={(e_after-e_before)*1000:.0f} meV)", flush=True)

anneal_results.sort(key=lambda x: x['e_after'])
champ = anneal_results[0]
print(f"\n★ Champion: rank {champ['rank']} (H0_Li{champ['li_trial']})")
print(f"  E = {champ['e_after']:.4f} eV")

# Save
with open('modelC_v2_results.json', 'w') as f:
    json.dump({
        'step1_best_halogen': [1, 4],
        'step2_top5_li_spread_meV': 0.1,
        'step2_overall_spread_meV': 157.1,
        'anneal': anneal_results
    }, f, indent=2)

print(f"\n{'='*60}")
print(f"★★★ Model C v2 COMPLETE ★★★")
print(f"Champion E = {champ['e_after']:.4f} eV")
print(f"Next: MLIP EOS → DFT EOS")
print(f"{'='*60}")
