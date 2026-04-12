# Br Substitution Effects on Argyrodite Mechanical Properties

## Summary
Replacing Cl⁻ with Br⁻ in argyrodite solid electrolytes systematically weakens
mechanical properties. This document explains the physical mechanisms.

## Five Mechanisms

### 1. Ionic Radius Effect
- Br⁻ (1.96 Å) > Cl⁻ (1.81 Å) [Shannon 1976]
- Larger ion → longer Li-X bonds → weaker ionic interaction
- Li-Br: 2.71 Å vs Li-Cl: 2.51 Å (Δ = 0.20 Å) — direct DFT evidence

### 2. Charge Density Effect
- Br⁻ has same formal charge (-1) but distributed over larger volume
- Lower charge density → weaker Coulomb interaction
- Bader analysis confirms: |q(Cl)| > |q(Br)| in all compositions

### 3. Polarizability Effect
- Br⁻ is more polarizable than Cl⁻
- More polarizable → electronic cloud deforms under stress → less resistance
- Manifests as lower C44 (shear) more than C11 (compression)

### 4. Lattice Expansion
- Br substitution expands the lattice (larger ion)
- Larger unit cell → more space between PS₄ units → weaker inter-tetrahedral coupling
- B0 decreases with Br content within each family

### 5. Indirect: S²⁻ Replacement (in Li₅.₄ family)
- Li₆PS₅Cl → Li₅.₄PS₄.₄Cl₁.₆: not just halogen but S²⁻ also decreases
- Li-S²⁻ electrostatic force ∝ q₁×q₂ = 1×2 = 2 (strong)
- Li-Cl⁻ electrostatic force ∝ q₁×q₂ = 1×1 = 1 (half as strong)
- S²⁻ → Cl⁻ replacement halves the local ionic bond strength

## Quantitative Trends

### Bulk Modulus (B0)
- Li₆: comp1(26.2) → comp2(25.8) = -0.4 GPa per 50% Br
- Li₅.₄: Model C(21.7) > comp3(20.8) = comp4(20.8): Br effect ~0.9 GPa
- Cross-family: Li₆(26) → Li₅.₄(21) = -5 GPa (vacancy + S reduction dominant)

### Young's Modulus (E, 600K snapshot)
- 29.1 → 28.6 → 27.3 → 26.4 → 25.8 GPa
- Total reduction: 3.3 GPa (11%) from comp1 to comp5
- Perfectly monotonic with Br content

### Band Gap
- 2.28 → 2.04 → 2.10 → 1.96 → 2.06 eV
- Slight reduction, Br 4p contributes to VBM
- PS₄ framework (covalent) unchanged → electronic structure mostly preserved

## Intrinsic vs Extrinsic

**Within same family** (same vacancy level):
- Br↑ → E↓ in BOTH calculation AND experiment ✅
- Intrinsic (bonding) effect dominates

**Across families** (different vacancy level):
- Calculation: Li₆(29) > Li₅.₄(26) → vacancy weakens ✅
- Experiment: Li₅.₄ > Li₆ → REVERSED!
- Cause: Grain boundary (GB) pinning by vacancy
  → Vacancy pins GB sliding → polycrystal appears stiffer
  → Extrinsic (microstructure) effect dominates cross-family comparison

## Practical Implications for Coating

| Aspect | Br↑ Effect | Coating Implication |
|--------|-----------|---------------------|
| E (intrinsic) | ↓ 12% | ✅ More compliant SE |
| Wad | ↑ (vacancy) | ✅ Better adhesion |
| Variability | ↑ (Li ordering) | ⚠️ Less reproducible |
| Band gap | Slight ↓ | ≈ Neutral |

**Optimal: comp4** — balanced E reduction + adhesion gain + basin stability.
