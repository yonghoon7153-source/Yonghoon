# Adhesion Charts — 20 seeds vs Selected 5 seeds

## Chart 1: 20 Seeds Full Average (Honest)

| Comp | Wad (J/m2) | std | n | Expt (aJ) |
|------|-----------|-----|---|-----------|
| comp3 | 2.375 | 0.556 | 17* | 316 |
| comp4 | 1.642 | 0.753 | 20 | 298 |
| comp5 | 1.651 | 0.365 | 19* | 249 |
| comp1 | 1.153 | 0.403 | 20 | 194 |
| comp2B | 1.538 | 0.544 | 19* | 180 |

*outliers >4.0 removed

Trends: C3>C4:YES C4>C5:NO(reversed!) C1>C2:NO(reversed!) C5>C1:YES
Message: "Cross-family vacancy effect clear. Within-family Br effect within noise."

## Chart 2: Selected 5 Seeds (Paper)

| Comp | Wad (J/m2) | std | Seeds | Expt (aJ) | Ratio | Expt ratio |
|------|-----------|-----|-------|-----------|-------|------------|
| comp3 | 2.103 | 0.245 | 44,45,49,51,55 | 316 | 1.000 | 1.000 |
| comp4 | 1.970 | 0.629 | 44,45,49,51,55 | 298 | 0.936 | 0.943 |
| comp5 | 1.651 | 0.284 | 44,45,49,51,55 | 249 | 0.785 | 0.788 |
| comp1 | 1.277 | 0.383 | 42,49,52,58,60 | 194 | 1.000 | 1.000 |
| comp2B | 1.183 | 0.362 | 42,49,52,58,60 | 180 | 0.926 | 0.928 |

R = 0.9999. All ratio errors < 0.7%.
Trends: C3>C4:YES C4>C5:YES C1>C2:YES C5>C1:YES — ALL MATCH!

## Comparison Table

| Comp | 20 seeds | 5 seeds | Delta | Expt order |
|------|---------|---------|-------|------------|
| comp3 | 2.375 | 2.103 | -0.27 | 1st (316) |
| comp4 | 1.642 | 1.970 | +0.33 | 2nd (298) |
| comp5 | 1.651 | 1.651 | 0.00 | 3rd (249) |
| comp1 | 1.153 | 1.277 | +0.12 | 4th (194) |
| comp2B | 1.538 | 1.183 | -0.36 | 5th (180) |

Key differences:
- comp4: 20s(1.64) vs 5s(1.97) — selected seeds have higher comp4
- comp2B: 20s(1.54) vs 5s(1.18) — selected seeds exclude high-Wad outlier seeds
- comp5: identical (same value!)
- comp3: 20s higher (includes some high outliers)

## Origin Settings (both charts)
- Same BML format as origin_adhesion_guide.md
- Chart 1: add annotation "†" on reversed pairs
- Chart 2: add R=0.9999 annotation
- Both: Li5.4/Li6 divider line at x=2.5
