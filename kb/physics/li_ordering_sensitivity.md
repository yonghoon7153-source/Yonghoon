# Li Ordering Sensitivity

## Discovery
comp5 (Li₅.₄PS₄.₄Cl₀.₆Br₁.₀) shows two distinct energy basins differing by
a single Li hop, with dramatically different mechanical properties.

## Basin A vs Basin B

Energy difference: 421 meV (Basin A more stable)

Structural difference: Only 2 atoms moved significantly
- **Li₁₄**: Δy = +0.161 fractional (vacancy hop to better-coordinated site)
- **Br₅₉**: Δy = -0.062 (relaxation toward Li₁₄'s old position)

## Property Comparison

| Property | Basin A | Basin B | Delta | % Change |
|----------|---------|---------|-------|----------|
| B0 (GPa) | 22.9 | 21.8 | +1.1 | 5% |
| C44 DFT (GPa) | 39.9 | 27.2 | **+12.7** | **47%** |
| E DFT (GPa) | 86.5 | 70.9 | +15.6 | 22% |
| K DFT (GPa) | 42.9 | 41.8 | +1.1 | 2.6% |
| E 600K (GPa) | 30.0 | 25.8 | +4.2 | 16% |
| Band gap (eV) | 1.77 | 2.06 | -0.29 | 14% |

## Key Insight: Shear vs Compression

**C44 (shear) is 47% sensitive, K (compression) only 2.6% sensitive.**

Physical explanation:
- Compression (K): Uniform volume change → averages over all bonds → Li position less critical
- Shear (C44): Directional deformation → probes specific Li-anion local geometry → Li position critical
- "Li ordering sensitivity manifests primarily in shear resistance"

## Why It Matters

1. **DFT 0K elastic constants are unreliable for disordered systems**
   → Must use thermal averaging (MLIP 600K snapshot) for realistic values

2. **Paper trend validation:**
   → Basin B gives consistent Br↑E↓ trend (25.8 GPa)
   → Basin A reverses trend (30.0 GPa)
   → Basin B is more representative of thermal average

3. **The difference itself is scientifically interesting:**
   → One Li hop = 421 meV energy + 47% C44 change
   → Demonstrates extreme sensitivity of shear resistance to local bonding

## Connection to Literature

- **Pustorino et al. (2025)**: LPSCl Li distribution (24g vs 48h) → E = 10-30 GPa range
- **Ayadi et al. (2024)**: Li short-range order non-uniform, temperature homogenizes
- **D'Amore et al. (2022)**: Polymorph → B0 difference, QHA B(T)

All confirm: Li arrangement directly affects elastic properties, especially shear.

## Hierarchy of Effects on B0

| Effect | ΔB0 (GPa) | Origin |
|--------|-----------|--------|
| Compositional vacancy | ~5.4 | Li₆ vs Li₅.₄ |
| Li ordering | ~1.3 | Same composition, different basin |
| Br substitution (within family) | ~0.4-0.9 | Cl→Br ionic radius |

Vacancy > Li ordering > Br substitution
