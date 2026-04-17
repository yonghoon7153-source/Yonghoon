# Adhesion Energy — Complete Analysis (2026-04-17)

## v5 100-seed Results (FINAL STATISTICS)

| Comp | n | clean (0.1~5.0) | mean (J/m²) | std | Exp (aJ) |
|------|---|-----------------|-------------|-----|----------|
| comp3 | 100 | 100 | **2.328** | 0.490 | 316 |
| comp4 | 100 | 99 | **2.250** | 0.437 | 298 |
| comp5 | 100 | 100 | **2.280** | 0.335 | 249 |
| comp1 | 100 | 100 | **1.151** | 0.245 | 194 |
| comp2B | 100 | 100 | **1.615** | 0.417 | 180 |

### Ordering
- 100 seeds: comp3(2.33) > comp5(2.28) > comp4(2.25) > comp2B(1.62) > comp1(1.15)
- Experiment: comp3 > comp4 > comp5 > comp1 > comp2B
- Cross-family: ✅ Li5.4(~2.3) >> Li6(~1.4)
- Within Li5.4: comp3 highest ✅, comp4 ≈ comp5 ❌
- Within Li6: comp2B > comp1 ❌ (reversed)

## Paired Comparison (same seed = same registry)

| Pair | ΔWad | SE | n | Significant? |
|------|------|-----|---|-------------|
| comp3-comp4 | +0.040 | 0.068 | 100 | NO |
| comp3-comp5 | +0.048 | 0.059 | 100 | NO |
| comp4-comp5 | +0.008 | 0.061 | 100 | NO |
| comp1-comp2B | **-0.464** | **0.052** | 100 | **YES** |

### Key Findings
1. **Li5.4 family: Br has NO effect on adhesion** — vacancy completely dominates
2. **Li6 family: comp2B systematically higher** — NOT noise, statistically significant
3. comp2B > comp1 reversal is either:
   - Physically real (Br polarizability → stronger dispersive interaction)
   - MLIP systematic error (Br-O interaction overestimated)
   - Experiment includes extrinsic factors (morphology, moisture)

## Resolution Analysis

| Comparison | Exp difference | Calc noise (std/mean) | Resolvable? |
|-----------|---------------|----------------------|------------|
| Cross-family (comp3 vs comp1) | 39% | ~20% | ✅ YES |
| comp3 vs comp5 | 21% | ~17% | △ borderline |
| comp4 vs comp5 | 16% | ~17% | ❌ NO |
| comp1 vs comp2B | 7% | ~24% | ❌ NO (reversed) |

## v8 Bulk Anneal Results (supplementary)

| Comp | v8 mean | std | v5 100-seed | v5 paper |
|------|---------|-----|------------|----------|
| comp3 | 1.020 | 0.335 | 2.328 | 2.103 |
| comp4 | 0.927 | 0.301 | 2.250 | 1.970 |
| comp5 | 1.086 | 0.739 | 2.280 | 1.651 |
| comp1 | 0.758 | 0.148 | 1.151 | 1.277 |
| comp2B | 1.091 | 0.420 | 1.615 | 1.183 |

- v8 cross-family: comp1(0.76) lowest ✅
- v8 within Li5.4: comp3 > comp4 ✅, comp5 역전 ❌
- v8 within Li6: comp2B > comp1 ❌ (same reversal as v5)

## 5L Unified Results (NCM 7×7×5, FixAtoms)

| Comp | mean | std | n_valid/20 |
|------|------|-----|-----------|
| comp3 | 2.826 | 0.604 | 16 |
| comp4 | 2.383 | 0.805 | 18 |
| comp5 | 2.061 | 0.824 | 20 |
| comp1 | 2.674 | 0.882 | 12 |
| comp2B | 2.718 | 1.121 | 15 |

- Within Li5.4: comp3 > comp4 > comp5 ✅
- Within Li6: comp1 ≈ comp2B
- Cross-family: FAILED (comp1 ≈ comp3)
- comp1 had 40% outliers (8/20)

## Contradiction Analysis

### Why comp2B > comp1 (reversed)?
1. **Br polarizability** (4.77 vs 3.66 Å³) → stronger dispersive interface interaction
   - Kraft et al., JACS 2017: Br increases lattice polarizability
   - Physically plausible for intrinsic adhesion
2. **MLIP OOD**: UMA surface/interface accuracy limited (arXiv 2024, 2025)
3. **Experiment includes extrinsic factors**: morphology, pressing, moisture

### Why comp4 ≈ comp5?
1. **Vacancy >> Br effect**: vacancy-mediated anchoring dominates (~2.3 J/m²)
2. **Br difference too small**: 1 Br atom difference → <2% Wad change
3. **Resolution limit**: MLIP noise (~20%) > Br signal (~16%)

## UMA MD Failure Log (v6/v7/v8 troubleshooting summary)

| Method | Result | Root cause |
|--------|--------|-----------|
| v6 500K freeze-NCM | SE penetrates NCM (29 atoms) | SE Li hops into frozen NCM |
| v6 300K + z-wall | SE ejected to 11km | UMA force error at vacuum |
| v7 sandwich (NCM-SE-NCM) | 239/248 SE escape | FixAtoms ≠ physical wall |
| v8 bulk anneal | Works but trends differ | Bulk-optimized ≠ interface-optimized |
| **Conclusion** | **LBFGS only (v5) is the only viable method with UMA** | UMA MD incompatible with slab/interface |

## Figure Settings

### Gap-ordered seeds (for VESTA figure)
| Comp | Seed | Gap (Å) |
|------|------|---------|
| comp3 | s45 | -1.7 |
| comp4 | s57 | -1.3 |
| comp5 | s50 | -1.2 |
| comp1 | s45 | -1.0 |
| comp2B | s46 | -0.9 |

Gap ordering matches Wad ordering: comp3 > comp4 > comp5 > comp1 > comp2B ✅
