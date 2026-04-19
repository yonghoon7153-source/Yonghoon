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

## Pending entries (priority order)
- Tier 1: **ACS AMI 2018 (Dewald et al.)** — X-CT 84% coverage, 3D metric (direct k_spread target)
- Tier 2: Strauss 2020 ACS AMI (mechanical stability vs coverage)
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
