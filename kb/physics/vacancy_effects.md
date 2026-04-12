# Vacancy Effects in Argyrodite

## Structural Disorder vs Compositional Vacancy

### Structural Disorder (Li₆PS₅Cl)
- 48h sites: 48 positions, 24 Li atoms (occupancy = 0.5)
- All Li present (stoichiometric)
- Empty sites = Li hopping destinations (ionic conduction)
- Time-average: all sites visited equally (occ ≈ 0.5)
- **Not a vacancy** — it's positional disorder

### Compositional Vacancy (Li₅.₄PS₄.₄Cl₁.₆)
- 48h sites: 48 positions, 27 Li atoms (5 f.u.)
- Li truly missing — 3 fewer per 5 f.u. compared to Li₆
- 21 permanently empty sites (vs 24 for Li₆)
- Extra 3 empty sites = permanent vacancies
- Rietveld: site-dependent occ = 0.3-0.5 (non-uniform!)
  - Br⁻-adjacent sites: lower occ (weaker binding)
  - PS₄-remote sites: lower occ
- **Real vacancy** — permanent bond density reduction

## Impact on Mechanical Properties

### Bulk Modulus
- Li₆ average: ~26 GPa
- Li₅.₄ average: ~21 GPa
- ΔB0 = -5 GPa (19% reduction)

Three combined effects:
1. **S²⁻ reduction** (5.0 → 4.4/f.u.): Strong Li-S²⁻ bonds (q₁q₂=2) replaced by weak Li-Cl⁻ (q₁q₂=1) → B0↓
2. **Li vacancy**: Empty coordination sites → lower bond density → B0↓
3. **Extra Cl**: Li-Cl weaker than Li-S → B0↓

### Surface Energy
- Li₆: γ_SE ≈ 1.2 J/m²
- Li₅.₄: γ_SE ≈ 0.5 J/m²
- Factor 2.4× reduction! Vacancy drastically weakens surface bonds.

### Adhesion Energy
- Vacancy creates under-coordinated surface Li = "chemical anchor"
- These Li pull NCM O²⁻ across interface → new cross-interface bonds
- Result: Wad can be HIGHER for Li₅.₄ despite lower intrinsic strength
- z-cut variance: large for Li₅.₄ (depends on which vacancies face NCM)

### Elastic Constants
- Li ordering sensitivity amplified by vacancy
- comp5 (Li₅.₄, high Br): ΔC44 = 12.7 GPa between basins
- comp1 (Li₆, no Br): basin effects expected smaller

## Experimental Consequence: Grain Boundary Pinning

In polycrystalline pellets:
- Vacancy pins grain boundary (GB) sliding
- GB sliding = main deformation mechanism in polycrystals
- Pinning → higher macroscopic modulus
- Result: Li₅.₄ pellet STIFFER than Li₆ pellet (intrinsic trend reversed!)

This is the origin of the intrinsic-extrinsic discrepancy:
- Single crystal (DFT): Li₆ > Li₅.₄ (intrinsic bonding)
- Polycrystal (experiment): Li₅.₄ > Li₆ (GB pinning)

## Quantitative Decomposition

| Effect | ΔB0 (GPa) | Mechanism |
|--------|-----------|-----------|
| Compositional vacancy | ~5.4 | Bond density reduction |
| Li ordering (same comp) | ~1.3 | Local bonding variation |
| Br substitution (within family) | ~0.4-0.9 | Ionic radius / polarizability |
| GB pinning (extrinsic) | Reverses trend | Microstructure, not bonding |
