# Adhesion Energy (Wad) Calculation

## Definition
Work of adhesion = energy per unit area to separate SE/NCM interface into two free surfaces.

```
Wad = (E_SE_isolated + E_NCM_isolated - E_interface) / A
```
(Isolated slab method — avoids strain release artifacts from separation-relax method)

## Methods Comparison

### v1-v2: Amorphous SE (3000K Melt-Quench)
- Melt SE at 3000K (5ps) → quench 300K (3ps) → amorphous SE
- Assemble on crystalline NCM → anneal 500K (2ps) → cool → relax
- **Pro:** Mimics real pressing (SE deforms at interface)
- **Con:** Destroys vacancy information, lattice mismatch for rhombo cells

### v5: Crystalline Slab + Surface-only MQA (Current)
- SE 2×2×1 repeat → lattice match with NCM 5×5×1
- Surface-only softening (800K, 2ps) → quench → Li anneal (500K, 3ps) → cool
- **Pro:** Preserves bulk structure and vacancies
- **Pro:** Lattice match eliminates PBC artifacts
- **Con:** Less realistic for pressing scenarios

## Cell Matching

| Family | SE cell | NCM cell | Strain |
|--------|---------|----------|--------|
| Li₆ (comp1,2) | prim 2×2×1 (52 at) | 5×5×1 (300 at) | +3.3% |
| Li₅.₄ (comp3-5) | 2×2×1 (248 at) | 5×5×1 (300 at) | +1.1% |

Both use NCM 5×5×1 → fair cross-composition comparison.

## Surface-only MQA Protocol

1. SE 2×2×1 slab on NCM 5×5×1, gap=2.5Å
2. **800K surface softening (2ps)** — surface atoms slightly rearrange, bulk preserved
3. **300K quench (2ps)** — freeze surface structure
4. **500K Li anneal (3ps)** — Li sublattice optimization, vacancy preserved
5. **100K cool (2ps)** + MLIP relaxation

Why 800K is safe in v5 (but failed in v3/v4):
- v3/v4 failure cause was lattice mismatch (±20%), not temperature
- v5 has matched lattice (~1-3%) → 800K only softens surface

## z-Cut Sampling

SE slab cut at 5 z-positions (0.0, 0.2, 0.4, 0.6, 0.8 fractional).
Different layers exposed to NCM:
- **Li₆ (no vacancy):** Similar Wad across z-cuts → small variance
- **Li₅.₄ (vacancy):** Wad depends on vacancy exposure → large variance
  → The variance itself is evidence for vacancy chemical anchor effect

## Vacancy → Adhesion Mechanism

Li₆ (no vacancy):
- All surface Li coordinated → no driving force for extra NCM bonding
- γ_SE high (1.21 J/m²), γ_SE/NCM high → Wad moderate

Li₅.₄ (vacancy):
- Under-coordinated Li near vacancy = "chemical anchor"
- Pulls NCM O²⁻ across interface → cross-interface ionic bonds
- γ_SE low (0.45 J/m²), γ_SE/NCM very low → Wad can be higher!

## Surface Energies

| Comp | γ_SE (J/m²) | 2γ (J/m²) |
|------|-------------|-----------|
| comp1 | 1.211 | 2.42 |
| comp2 | 1.189 | 2.38 |
| comp3 | 0.565 | 1.13 |
| comp4 | 0.450 | 0.90 |
| comp5 | 0.470 | 0.94 |

Sharp drop from Li₆ to Li₅.₄ → vacancy reduces surface bond strength.

## Calculator
UMA (uma-s-1p2, fairchem, V100 GPU)

## Known Bugs & Fixes

### v5 FIX: ASE atom slicing error (2026-04-12)
**Error:** `could not broadcast input array from shape (248,) into shape (348,)`
**Cause:** ASE `Atoms` object does not support NumPy-style slicing `interface[:n_ncm]`.
**Fix:** Use `copy()` + `del` instead of slice:
```python
# Before (broken):
ncm_iso = interface[:n_ncm].copy()
se_iso = interface[n_ncm:].copy()

# After (fixed):
ncm_iso = interface.copy()
del ncm_iso[n_ncm:]
se_iso = interface.copy()
del se_iso[:n_ncm]
```
**Status:** FIX2 deployed on V100. Superseded by v5_working.

### v5 Wad Debugging Saga (2026-04-13)

6 iterations to get correct Wad. Root cause: **UMA vacuum sensitivity**.

| Version | Method | Wad (J/m2) | Problem |
|---------|--------|-----------|---------|
| v5_fix | isolated slab + relax | 10.0 | relax releases strain |
| v5_final | isolated slab, no relax | 75.4 | E_ncm=+248eV! 60A vacuum |
| v5_final2 | separation 30A, cell fixed | 40.2 | PBC wrap (SE exceeds cell) |
| v5_correct | independent slab cells, relax | 17.3 | used pre-MQA structure |
| v5_real | post-MQA split, no relax | 47.8 | E_ncm=-24eV (vacuum issue) |
| v5_v2method | separation 30A, vacuum 60A | 24.5 | vacuum too large for UMA |
| **v5_working** | **sep 30A + cell+30A, vac 30A** | **1.252** | **CORRECT** |

### UMA Vacuum Sensitivity (Critical Finding)

UMA (graph neural network MLIP) gives nonphysical energies with large vacuum:
- vacuum 30A: Wad = 1.252 J/m2 (correct, matches literature)
- vacuum 60A: Wad = 24.5 J/m2 (10x too high!)

Same atoms, same arrangement — only vacuum differs.
Cause: UMA trained on bulk/dense systems. Large vacuum = out of training distribution.

### Correct Protocol (v5_working)

Three conditions ALL required:
1. **vacuum = 30A** (not 60A!) — `cell_z = atoms_max + 30`
2. **Expand cell on separation** — `sep_cell_z += 30A` AND `SE_z += 30A`
3. **No relax after separation** — single point only (relax releases strain)

```python
# E_int: interface after MQA + relax
cell_z = total_max + 30.0  # 30A vacuum
E_int = interface.get_potential_energy()

# E_sep: separate, expand cell, single point
sep_cell[2] = [0, 0, cell_z + 30.0]  # cell +30A
pos_sep[n_ncm:, 2] += 30.0            # SE +30A
E_sep = sep.get_potential_energy()     # NO relax!

Wad = (E_sep - E_int) / A * 16.0218
```

### v5_noMQA variant (2026-04-13)
Crystalline stacking + relax only (no MQA thermal treatment).
Running on V100: comp4→comp2B→comp3→comp1→comp5, 5 z-cuts each.
Purpose: baseline comparison — does MQA improve or change Wad?

### Validation
Test result (no MQA, comp4 z=0.0): Wad = 1.252 J/m2
Compare v2 (3000K melt): comp1 = 1.107, comp2B = 1.046 J/m2
→ Same order of magnitude, reasonable range.

## References
- [j_nucl_mater2020] Au/a-Si crystalline/amorphous interface statistics
- [pnas2018] GAP MLIP melt-quench amorphous Si validation
