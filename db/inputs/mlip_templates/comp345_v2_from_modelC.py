"""
Comp3/4/5 v2 Pipeline — Br doping from Model C champion
Model C champion halogen: S at free_sites [1, 4]
Option A: fix S sites, enumerate Br on remaining 8 Cl sites

comp3: Li5.4PS4.4Cl1.0Br0.6 → 3 Br → C(8,3) = 56
comp4: Li5.4PS4.4Cl0.8Br0.8 → 4 Br → C(8,4) = 70
comp5: Li5.4PS4.4Cl0.6Br1.0 → 5 Br → C(8,5) = 56
"""
import numpy as np
from pymatgen.core import Structure
from pymatgen.io.ase import AseAtomsAdaptor
from itertools import combinations
from ase import units
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution
from ase.optimize import LBFGS
from fairchem.core import pretrained_mlip
from fairchem.core.calculate.ase_calculator import FAIRChemCalculator
from ase.io import write
import json, os

predictor = pretrained_mlip.get_predict_unit("uma-s-1p1", device="cuda")
adaptor = AseAtomsAdaptor()

def new_calc():
    return FAIRChemCalculator(predictor, task_name="omat")

# ═══════════════════════════════════════════
# 1. Load reference structure + classify sites
# ═══════════════════════════════════════════
ref = Structure.from_file('/scratch/x3430a02/kgy/manuscript_support/comp3_lpsc10b06/configs/comp3_lpsc10b06_config_000.cif')
print(f"Reference: {ref.composition}, {len(ref)} atoms")

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
n_li_select = 27  # Li5.4 × 5 f.u.
print(f"Li: {n_li}, P: {len(p_sites)}, S_fw: {len(s_framework)}, Free: {len(free_sites)}")

# Model C champion: S at free_sites [1, 4]
S_INDICES = [1, 4]  # Fixed from Model C v2 champion
print(f"Model C champion S positions: free_sites {S_INDICES}")

# Remaining 8 free sites for Cl/Br assignment
cl_br_site_indices = [i for i in range(len(free_sites)) if i not in S_INDICES]
print(f"Available Cl/Br sites: {len(cl_br_site_indices)} sites at indices {cl_br_site_indices}")

# ═══════════════════════════════════════════
# 2. Composition definitions
# ═══════════════════════════════════════════
compositions = {
    'comp3': {'name': 'Li5.4PS4.4Cl1.0Br0.6', 'n_br': 3, 'n_configs': 56},
    'comp4': {'name': 'Li5.4PS4.4Cl0.8Br0.8', 'n_br': 4, 'n_configs': 70},
    'comp5': {'name': 'Li5.4PS4.4Cl0.6Br1.0', 'n_br': 5, 'n_configs': 56},
}

all_results = {}

for comp_id, comp_info in compositions.items():
    n_br = comp_info['n_br']
    print(f"\n{'='*60}")
    print(f"  {comp_id}: {comp_info['name']}")
    print(f"  Br sites: {n_br}, Cl sites: {8-n_br}")
    print(f"  Br configs: C(8,{n_br}) = {comp_info['n_configs']}")
    print(f"{'='*60}")

    # ═══════════════════════════════════════
    # Step 2a: Br screening (Li fixed from reference)
    # ═══════════════════════════════════════
    print(f"\n--- Step 2a: Br screening ({comp_info['n_configs']} configs, Li fixed) ---")

    br_combos = list(combinations(range(8), n_br))  # Choose n_br from 8 Cl/Br sites
    br_results = []

    for idx, br_local in enumerate(br_combos):
        # br_local: indices within cl_br_site_indices (0-7)
        # Map to actual free_site indices
        br_free_indices = [cl_br_site_indices[j] for j in br_local]

        species = []; coords = []
        # Li (all from reference)
        for c in li_sites[:n_li_select]:
            species.append('Li'); coords.append(c)
        # P
        for c in p_sites:
            species.append('P'); coords.append(c)
        # S framework
        for c in s_framework:
            species.append('S'); coords.append(c)
        # Free sites: S at S_INDICES, Br at br_free_indices, Cl at rest
        for i, c in enumerate(free_sites):
            if i in S_INDICES:
                species.append('S')
            elif i in br_free_indices:
                species.append('Br')
            else:
                species.append('Cl')
            coords.append(c)

        struct = Structure(ref.lattice, species, coords)
        atoms = adaptor.get_atoms(struct)
        atoms.calc = new_calc()
        try: LBFGS(atoms, logfile=None).run(fmax=0.01, steps=200)
        except: pass
        e = atoms.get_potential_energy()
        br_results.append({'idx': idx, 'energy': e, 'br_local': list(br_local),
                           'br_free': br_free_indices})
        print(f"  {idx}/{len(br_combos)}: E={e:.4f}, Br at local {list(br_local)}", flush=True)

    br_results.sort(key=lambda x: x['energy'])
    print(f"\n  Top 5 Br configs:")
    for r in br_results[:5]:
        print(f"    #{r['idx']}: E={r['energy']:.4f}, Br at {r['br_local']}")

    # ═══════════════════════════════════════
    # Step 2b: Top 5 Br × 20 Li = 100 configs
    # ═══════════════════════════════════════
    print(f"\n--- Step 2b: Top 5 Br × 20 Li = 100 configs ---")

    rng = np.random.RandomState(42)
    Li_configs = [sorted(rng.choice(n_li, n_li_select, replace=False).tolist()) for _ in range(20)]

    li_results = []
    for br_rank, br_res in enumerate(br_results[:5]):
        br_free_indices = br_res['br_free']
        print(f"\n  --- Br #{br_res['idx']} (rank {br_rank}) ---", flush=True)

        for li_trial, li_select in enumerate(Li_configs):
            species = []; coords = []
            for i in li_select: species.append('Li'); coords.append(li_sites[i])
            for c in p_sites: species.append('P'); coords.append(c)
            for c in s_framework: species.append('S'); coords.append(c)
            for i, c in enumerate(free_sites):
                if i in S_INDICES:
                    species.append('S')
                elif i in br_free_indices:
                    species.append('Br')
                else:
                    species.append('Cl')
                coords.append(c)

            struct = Structure(ref.lattice, species, coords)
            atoms = adaptor.get_atoms(struct)
            atoms.calc = new_calc()
            try: LBFGS(atoms, logfile=None).run(fmax=0.01, steps=200)
            except: pass
            e = atoms.get_potential_energy()

            li_results.append({
                'br_rank': br_rank, 'br_idx': br_res['idx'],
                'br_local': br_res['br_local'], 'br_free': br_free_indices,
                'li_trial': li_trial, 'li_indices': li_select, 'energy': e
            })
            print(f"    Br{br_rank}_Li{li_trial}: E={e:.4f}", flush=True)

    li_results.sort(key=lambda x: x['energy'])
    E_spread = (li_results[-1]['energy'] - li_results[0]['energy']) * 1000
    print(f"\n  Energy spread: {E_spread:.1f} meV")
    print(f"  Top 5:")
    for r in li_results[:5]:
        print(f"    Br{r['br_rank']}_Li{r['li_trial']}: E={r['energy']:.4f}, Br={r['br_local']}")

    # ═══════════════════════════════════════
    # Step 3: Top 5 → annealing
    # ═══════════════════════════════════════
    print(f"\n--- Step 3: Top 5 annealing (500K, 50ps) ---")

    anneal_results = []
    for rank, r in enumerate(li_results[:5]):
        br_free_indices = r['br_free']
        li_select = r['li_indices']

        species = []; coords = []
        for i in li_select: species.append('Li'); coords.append(li_sites[i])
        for c in p_sites: species.append('P'); coords.append(c)
        for c in s_framework: species.append('S'); coords.append(c)
        for i, c in enumerate(free_sites):
            if i in S_INDICES:
                species.append('S')
            elif i in br_free_indices:
                species.append('Br')
            else:
                species.append('Cl')
            coords.append(c)

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
            'rank': rank, 'br_rank': r['br_rank'], 'li_trial': r['li_trial'],
            'br_local': r['br_local'], 'e_before': e_before, 'e_after': e_after
        })
        write(f'{comp_id}_v2_anneal_rank{rank}.xyz', atoms)
        print(f"  Rank {rank} (Br{r['br_rank']}_Li{r['li_trial']}): "
              f"{e_before:.4f} → {e_after:.4f} (dE={(e_after-e_before)*1000:.0f} meV)", flush=True)

    anneal_results.sort(key=lambda x: x['e_after'])
    champ = anneal_results[0]
    print(f"\n  ★ {comp_id} Champion: rank {champ['rank']} "
          f"(Br{champ['br_rank']}_Li{champ['li_trial']})")
    print(f"    E = {champ['e_after']:.4f} eV")
    print(f"    Br config: {champ['br_local']}")

    all_results[comp_id] = {
        'br_screening_top5': [{'idx': r['idx'], 'energy': r['energy'],
                               'br_local': r['br_local']} for r in br_results[:5]],
        'li_spread_meV': E_spread,
        'li_top5': [{'br_rank': r['br_rank'], 'li_trial': r['li_trial'],
                     'energy': r['energy']} for r in li_results[:5]],
        'anneal': anneal_results,
        'champion_energy': champ['e_after'],
        'champion_br': champ['br_local']
    }

# ═══════════════════════════════════════════
# Save all results
# ═══════════════════════════════════════════
with open('comp345_v2_results.json', 'w') as f:
    json.dump(all_results, f, indent=2, default=str)

print(f"\n{'='*60}")
print(f"★★★ ALL COMPOSITIONS COMPLETE ★★★")
for comp_id in ['comp3', 'comp4', 'comp5']:
    r = all_results[comp_id]
    print(f"  {comp_id}: E={r['champion_energy']:.4f}, Br={r['champion_br']}")
print(f"{'='*60}")
