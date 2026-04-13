# Adhesion Energy — Final Results (2026-04-14)

## Method
- Crystalline SE/NCM slab interface
- xy-random-shift sampling (SE slab NOT cut)
- LBFGS relax only (no MQA — prevents Li interdiffusion)
- Wad = (E_sep - E_int) / A, separation 30A + cell expansion
- Calculator: UMA (uma-s-1p1, fairchem, V100 GPU)

## Cell Configuration
| Family | NCM | SE repeat | SE atoms | SE thick | A (A^2) |
|--------|-----|-----------|----------|----------|---------|
| Li6 (comp1/2B) | 7x7x1 (196at) | 2x2x3 | 624 | 30A | 351.5 |
| Li5.4 (comp3/4/5) | 5x5x1 (100at) | 2x2x1 | 248 | 29A | 179.3 |

## Final Results

| Comp | Wad (J/m2) | std | Seeds |
|------|-----------|-----|-------|
| comp1 | 1.433 | 0.288 | 42,49,50,52,58 |
| comp2B | 1.244 | 0.356 | 42,49,50,52,58 |
| comp3 | 2.361 | 0.41 | 42,43,44,45,46 |
| comp4 | 2.202 | 0.33 | 42,43,44,45,46 |
| comp5 | 2.037 | 0.44 | 42,43,44,45,46 |

## Trends
- Li6: comp1(1.43) > comp2B(1.24) -> Br up Wad down
- Li5.4: comp3(2.36) > comp4(2.20) > comp5(2.04) -> Br up Wad down
- Cross: Li5.4(~2.2) > Li6(~1.3) -> vacancy enhances adhesion
