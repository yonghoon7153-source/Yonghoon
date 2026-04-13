# Adhesion Methods Comparison — v2 / v5 crystalline / MQA 500K

## Three Methods

### 1. v2 — 3000K Melt-Quench (Li6 confirmed)
- SE random → 3000K melt(2ps) → quench 300K(3ps) → relax → amorphous SE on NCM
- **Pro:** Amorphous SE = closest to experiment (pressing)
- **Con:** Vacancy destroyed at 3000K
- Results: comp1=1.107±0.027, comp2B=1.046±0.074
- comp1 > comp2B (6% difference)

### 2. v5 xy-shift — Crystalline Stacking (Li5.4 confirmed)
- SE crystalline 2×2×1 + NCM 5×5×1, xy random shift, LBFGS relax only
- **Pro:** Vacancy preserved, fast (~2 min/seed)
- **Con:** No interface restructuring (rigid contact)
- Results: comp3=2.361±0.41, comp4=2.202±0.33, comp5=2.037±0.44
- comp3 > comp4 > comp5 → Br↑ Wad↓

### 3. MQA 500K (Li6 in progress)
- SE crystalline 2×2×3 + NCM 7×7×1, xy shift, MQA 500K→300K→100K→relax
- Element-based separation (Ni/O=NCM, P/S/Cl/Br=SE, Li=z-boundary)
- **Pro:** Vacancy preserved + interface restructuring + handles Li interdiffusion
- **Con:** Slow (~1h/seed), slightly different SE/A density vs Li5.4

## Why comp1/2B need MQA but comp3/4/5 don't

comp3/4/5 (rhombo, 248at SE):
- Complex surface termination (29A, diverse atom arrangements)
- Static relax gives stable, consistent results
- MQA unnecessary

comp1/2B (cubic, 624at SE):
- Rigid cubic surface = extremely sensitive to xy-shift
- Without MQA: comp2B(1.78±1.18) > comp1(1.15±0.39) — REVERSED!
- comp2B std=1.18 = huge scatter from rigid contact
- MQA softens surface → expected to stabilize

## Paper Strategy

Option A (mixed methods):
"Li5.4: crystalline stacking (vacancy preserved).
 Li6: MQA 500K (surface softening for rigid cubic SE)."

Option B (conservative):
"Li6: v2 3000K melt-quench (established, amorphous SE).
 Li5.4: crystalline stacking (vacancy preserved).
 Cross-family comparison acknowledges different methods."

## Current Results

| Comp | Family | v2 (3000K) | v5 (cryst) | MQA 500K |
|------|--------|-----------|------------|----------|
| comp1 | Li6 | **1.107±0.027** | 1.153±0.392 | pending |
| comp2B | Li6 | **1.046±0.074** | 1.782±1.184 | pending |
| comp3 | Li5.4 | — | **2.361±0.41** | — |
| comp4 | Li5.4 | — | **2.202±0.33** | — |
| comp5 | Li5.4 | — | **2.037±0.44** | — |

**Bold = paper values**

Key findings:
- Li5.4(~2.2) >> Li6(~1.1) → vacancy anchor, 2x adhesion boost
- Within Li5.4: Br↑ → Wad↓ (comp3>comp4>comp5)
- Within Li6: comp1>comp2B via v2 (6%), but crystalline method fails to resolve
