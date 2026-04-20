# Literature Coverage Database

Purpose: calibrate `plastic_coverage.py` k_spread against physical coverage values
reported in ASSB literature. Single source of truth for process conditions,
coverage metrics, and comparability notes.

## Files
- `coverage_db.json` — machine-readable entries (feed into calibration scripts)
- `README.md` — this file, human summary

## Our baseline
- Process: single-compression pellet, ~300 MPa
- DEM raw coverage: 85% (uncapped, non-physical)
- DEM plastic-capped coverage: 22% (geometric cap a² ≤ min(R1,R2)²)
- Metric: 3D surface-area coverage

## Entries

### 1. Lee 2024 Nature Communications (Tier 1, user-supplied text + SI)
| Process        | Coverage (2D) | Matches ours? |
|----------------|---------------|---------------|
| Hand-mixed     | 30.6%         | YES (single compression) |
| Wet-process    | 33.3%         | PARTIAL (slurry route)   |
| Dry-process    | 67.2%         | NO (high-shear mixing)   |

**Key insight from SI:**
- Coverage metric is **2D SEM perimeter**, not 3D surface
- Pressure: 7 TON on 13 mm = **~520 MPa** (not 300 MPa)
- AM radius range: 2-10 μm (representative R = 10 μm)
- SE radius range: 0.5-2 μm
- DEM validation uses Hertz-Mindlin with E_NCM=140 GPa, E_LPSCl=24 GPa
- Cell: Li-In | LPSCl (150 mg, 1.4 TON) | NMC532+LPSCl+CNT (80:18.5:1.5)

**Takeaway for k_spread:**
- Naive: target 30.6% vs DEM 22% → k_spread ≈ 1.18
- Corrected: 2D-to-3D factor ~1.3-1.5x → true 3D target ≈ 0.40 → k_spread ≈ 1.35
- DO NOT use dry-process 67.2% as calibration target (different process)

### 2. Bielefeld 2019 JPCC (Tier 1, user-supplied full text)
**Kind: computational microstructural modeling (NOT experiment).**

| Property                        | Value |
|---------------------------------|-------|
| Percolation exponent β          | 0.41 (3D site-percolation — validates our scaling law) |
| Threshold formula               | pc[vol%] = 7.83·ln(d/μm) + 36.67 |
| Optimal 5% porosity             | 80/20 wt NCM622:LPS |
| Optimal 10% porosity            | 82/18 wt |
| Optimal 20% porosity            | 86/14 wt |
| Operating band @ 20% porosity   | 69-79 vol% AM |
| AM shape                        | spheres (no overlap) |
| SE shape                        | convex polyhedra (with overlap, absorbed by AM at merge) |
| Metric                          | A_spec,a (m²/m³), not % coverage |

**Correction to earlier claim:** This paper does **NOT** give τ vs composition
or σ_eff vs composition curves directly. It gives percolation metrics. Still
valuable: validates our scaling law exponents (β, CN^1.5 dependence) and gives
optimal composition anchors that match our dataset.

**Coverage equivalent:** can be back-derived as
  coverage ≈ A_spec,a / (g_AM^V · 6 / d_AM)
but requires reading plot values (not tabulated).

**AM-SE interface coverage (digitized Fig 7, 8, 10):**
| Porosity | Peak A_spec,a | Peak composition | Coverage |
|----------|---------------|------------------|----------|
| 5%       | 5.8e5 m²/m³   | 65/35 vol%       | **78%**  |
| 10%      | 5.0e5 m²/m³   | 65/35 vol%       | **71%**  |
| 20%      | 3.3e5 m²/m³   | 70/30 vol%       | **49%**  |

Thickness (Fig 10, 20-140 μm) does NOT change peak coverage.

**Porosity sweep at 70/30 vol% composition (Fig 9, 5μm AM):**
| Porosity | AM total vol% | A_spec,a | Coverage |
|----------|---------------|----------|----------|
| 43%      | 40            | 0.01e5   | 0.2%     |
| 30%      | 49            | 0.50e5   | 8.5%     |
| 25%      | 52.5          | 2.00e5   | 31.7%    |
| 20%      | 56            | 3.30e5   | 49.1%    |
| 15%      | 59.5          | 4.30e5   | 60.2%    |
| 10%      | 63            | 5.20e5   | 68.8%    |
| 5%       | 66.5          | 5.90e5   | 73.9%    |
| 3%       | 67.9          | 6.00e5   | 73.6%    |

Regimes: ionic+electronic limited above φ≈34%; electronic-only limited
21-34%; well-connected below 21%. Monotonic coverage increase with
densification.

**k_spread calibration anchor (20% porosity, 5μm):**
- Our DEM plastic-capped: 22%
- Bielefeld 20%-porosity peak: 49%
- Ratio: 2.2× under-estimate → **k_spread ≈ 1.49**
- Cross-check with Lee 2024 hand-mixed (30.6% 2D → ~39% 3D): same order
- **Final k_spread range: 1.30-1.50**

### 3. Hlushkou 2018 JPS (Tier 1, user-supplied main + SI)
**Kind: experimental 3D FIB-SEM reconstruction + random-walk simulation.**

| Property                         | Value               |
|----------------------------------|---------------------|
| Pressure (electrochemistry)      | **276 MPa** (≈ our 300 MPa) |
| Chemistry                        | LCO/LPSI (partial NMC/LPSCl match) |
| AM                               | 5 μm LCO, LiNbO3-coated, irregular shape |
| SE                               | 0.67(0.75Li2S-0.25P2S5)-0.33LiI, σ=0.7 mS/cm |
| Initial mixture (vol)            | 38% LCO / 62% SE / 0% void |
| Reconstructed (vol)              | 33.1% LCO / 53.7% SE / **13.2% void** |
| τ_cond (EIS)                     | 1.6 ± 0.1           |
| τ_diff (FIB-SEM)                 | 1.74                |
| τ_Bruggeman at ε=0.537           | 1.34 (under-predicts) |
| Void-free τ_diff                 | 1.27                |
| D_eff / D_electrolyte (actual)   | 0.574               |
| D_eff / D_electrolyte (void-free)| 0.786               |

**Critical insight:** 276 MPa pressing ≠ full densification — **13.2% voids
remain**. Our DEM should target ~13% porosity, not ~20%, as our-match anchor.

**Coverage back-estimate (isotropic assumption):**
coverage_AM-SE ≈ V_SE/(V_SE+V_void) = 0.537/0.669 = **80.3%**
(upper bound — true value lower if voids cling to AM surface).

**Bruggeman failure:** τ=1.74 vs τ_B=1.34 → voids redistribute SE into
tortuous paths, not just reduce volume. Our C_blend(τ) term should allow
deviation from ε^(-0.5) at high void content.

## k_spread summary so far
| Source             | Porosity | Coverage | Method                    |
|--------------------|----------|----------|---------------------------|
| Lee 2024 (hand)    | ~unknown | 30.6% 2D | SEM perimeter (2D→3D ~39%) |
| Bielefeld 2019     | 20%      | 49%      | Computational GeoDict      |
| Bielefeld 2019     | 13%      | ~65%     | interpolated from Fig 9    |
| Hlushkou 2018      | 13.2%    | ~80% (UB)| FIB-SEM, isotropic bound   |
| **Our DEM plastic**| ~13-20%? | 22%      | DEM+δ/R cap                |

At 13% porosity: Bielefeld 65% vs Hlushkou 80% (UB) — real pressed pellet
matches or exceeds computational estimate. Our 22% DEM is **3-3.6×
under-estimate** at this porosity.

**Revised k_spread range: 1.55-1.90** (up from earlier 1.30-1.50 at fixed
20% porosity assumption).

## Pending entries (priority order)
- Tier 1: Strauss 2018 ACS Energy Lett (NCM622/argyrodite size vs inactive fraction)
- Tier 2: Minnmann/Neumann (NCM/LPSCl 3D FIB-SEM — chemistry-exact match)
- Tier 2: Shi 2020 JMCA (porosity-coverage correlation)
- Tier 3: Minnmann 2021 JECS (rate vs coverage)
- Tier 4: Review papers for cross-checks

## Open issue: 2D-to-3D conversion
For randomly sectioned spheres touching SE:
- Perimeter coverage (2D SEM) ≤ surface coverage (3D)
- Geometric factor depends on spatial isotropy of SE contacts
- Literature: Underwood stereology, E(3D) ≈ (4/π) · E(2D perimeter) ≈ 1.27x
- → Lee 2024 hand-mixed 30.6% (2D) ≈ 0.39 (3D) as first estimate

## How to use
```python
import json
db = json.load(open("docs/literature_coverage/coverage_db.json"))
for entry in db["entries"]:
    if entry["comparability_with_ours"]["process_match"] == "YES":
        print(entry["id"], entry["coverage_values_pct"])
```
