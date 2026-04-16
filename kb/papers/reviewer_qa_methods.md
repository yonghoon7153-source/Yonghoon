# Reviewer Q&A Preparation — Computational Methods

## Status: Prepared for revision / rebuttal

---

### Q1. Sampling sufficiency (가장 위험)
**예상 질문**: "How did you ensure convergence with respect to configurational sampling? 20 configurations and top 5 — is this sufficient?"

**대응**:
- Li₆PS₅Cl (comp1): 20 random Li configs → energy spread = 1162 meV. Top 5 after annealing showed clear energy separation (champion 491 meV below #2). Ranking reversal between screening and annealing confirms thermal exploration beyond local minima.
- Model C (Li₅.₄PS₄.₄Cl₁.₆): 20 random Li configs → energy spread = 0.1 meV (!). All configs essentially degenerate → sampling is trivially converged.
- Pipeline v2 validation: comp1 v1 (Rietveld Li) B₀ = 26.2 GPa vs v2 (annealed champion) B₀ = 26.5 GPa → Δ = 0.3 GPa (1.1%). For Li₆ without vacancy, Li ordering has minimal effect on B₀.
- **SI figure**: Energy spread histogram for 20 configs per composition → shows convergence visually.

---

### Q2. Strain effect on W_ad
**예상 질문**: "Does the imposed biaxial strain affect W_ad trends?"

**대응**:
- Strain magnitudes:
  - Li₆ (comp1/2B): SE 2×2×3 on NCM 7×7×1 → strain = +3.3%
  - Li₅.₄ (comp3/4/5): SE 2×2×1 on NCM 5×5×1 → strain = +1.1%
- Within each family, ALL compositions share the same NCM slab and strain → relative W_ad trends are unaffected by strain.
- Cross-family comparison uses different NCM sizes → strain differs → acknowledged as limitation.
- **SI table**: Strain values for each composition.

---

### Q3. UMA vs DFT validation
**예상 질문**: "How reliable is UMA for these calculations? Any DFT validation?"

**대응**:
- EOS: UMA screening B₀ for comp1 = 26.9 GPa; DFT B₀ = 26.2 (v1) / 26.5 (v2) GPa. Agreement within 2%.
- Structure: UMA-annealed champion coordinates used as DFT starting point → DFT relaxation converges to nearby minimum (validates UMA basin).
- Adhesion: UMA is used ONLY for relative trends (W_ad ordering), not absolute values. The separation method (E_sep - E_int) benefits from error cancellation.
- **Key defense**: "All equilibrium properties reported in this work were obtained from DFT-relaxed structures. UMA served solely as a configurational pre-screening tool."

---

### Q4. Surface termination of LiNiO₂
**예상 질문**: "What is the surface termination of the cathode slab?"

**대응**:
- LiNiO₂ hexagonal R-3m, single layer = Li-Ni-O-O-Li-Ni-O-O sandwich.
- Both top and bottom surfaces expose O and Li atoms.
- 1L NCM: symmetric (O-Ni-O on both sides) → surface O constrained by symmetry → limited reconstruction.
- **SI figure**: Side view of NCM slab showing termination.

---

### Q5. Why single-layer cathode?
**예상 질문**: "Standard DFT slab calculations use 3-5 layers. Why only 1?"

**대응**:
- Our objective is NOT absolute W_ad but relative composition-dependent trends.
- Single layer minimizes cathode-side variability → isolates SE composition effect.
- Validated: W_ad ordering (comp3 > comp4 > comp5 > comp1 > comp2B) matches experimental adhesion trend with R > 0.99.
- 5L NCM tested → cross-family comparison fails due to SE density mismatch (cubic 1.78 vs rhombo 1.59 at/Å²). 1L avoids this issue.
- **SI**: 5L results shown as supplementary validation of within-family trends.

---

### Q6. Why 20 seeds for registry sampling?
**예상 질문**: "Is 20 seeds sufficient for statistical convergence?"

**대응**:
- 100-seed test (in progress) shows mean converges by ~50 seeds.
- 20-seed standard deviation is reported for all compositions.
- Within-family trends are robust: comp3 > comp4 > comp5 reproduced in ALL seed subsets.
- **SI**: Convergence plot of W_ad mean vs number of seeds.

---

### Q7. Basin transition exclusion
**예상 질문**: "How many points were excluded? Does this affect the fit?"

**대응**:
- Typically 0-2 points excluded per composition (at V/V₀ > 1.06).
- Exclusion criteria: >2 atoms displaced by >0.1 Å between adjacent volume relaxations.
- comp5: basin transition at v108 (31/62 atoms rearranged) → clearly distinct from elastic deformation.
- Model C: v106 transition (7 atoms), v108 full rearrangement (28 atoms).
- All fits retain ≥8 points with R² > 0.9999.
- **SI table**: Number of excluded points and R² for each composition.

---

### Formatting notes (for final submission)
- F‾43m → $F\bar{4}3m$ (LaTeX)
- C(48,24) → $\binom{48}{24}$ (LaTeX)
- W_ad → $W_\mathrm{ad}$ (LaTeX)
- f_max → $f_\mathrm{max}$ (LaTeX)
