# Adhesion Troubleshooting Log

## Issue History (chronological)

### 1. v2 3000K Melt — Vacancy Destruction
- **Problem**: 3000K melts SE → amorphous → vacancy structure destroyed
- **Evidence**: v2 comp3(Li5.4) Wad≈1.0 = same as comp1(Li6) → no cross-family difference
- **Fix**: Use crystalline SE (v5 method)

### 2. v3/v4 Lattice Mismatch — PBC Artifacts
- **Problem**: SE rhombo 5×1×1 (7.1Å) vs NCM 3×3 (8.6Å) = ±20% mismatch
- **Evidence**: Wad = 0.006~8.6 (extreme scatter), seed44/45 always fail
- **Fix**: SE 2×2×1 repeat → strain < 3%

### 3. z-shift — Slab Cutting
- **Problem**: z-shift cuts crystalline SE → dangling bonds → nonphysical surfaces
- **Evidence**: comp4 Wad = 0.006~2.719, comp2B z=0.4 Wad=9.252
- **Fix**: xy-shift only (PBC-safe, preserves slab integrity)

### 4. UMA Vacuum Sensitivity
- **Problem**: UMA gives nonphysical energies with large vacuum (>40Å)
- **Evidence**: vacuum=30Å → Wad=1.25, vacuum=60Å → Wad=24.5 (same atoms!)
- **Fix**: Always use vacuum=30Å. UMA trained on bulk/dense systems.

### 5. Isolated Slab Method — UMA Incompatible
- **Problem**: E_ncm + E_se - E_int requires isolated slab energies
- **Evidence**: NCM 100at in 70Å cell → E=+248 eV (positive! nonsense)
- **Fix**: Use separation method: Wad = (E_sep - E_int) / A

### 6. Separation Without Cell Expansion — PBC Wrap
- **Problem**: SE moved +30Å but cell_z unchanged → SE wraps below NCM via PBC
- **Evidence**: Wad=40 J/m² (nonphysical)
- **Fix**: cell_z += 30Å when separating

### 7. Separation + Relax — Strain Release
- **Problem**: Relaxing after separation releases SE strain from NCM cell matching
- **Evidence**: E_sep < E_int → Wad negative (-1.9 J/m²)
- **Fix**: Single point only after separation (no relax)

### 8. New Calculator Required
- **Problem**: Reusing same UMA calculator after modifying atoms → shape mismatch error
- **Evidence**: "could not broadcast input array from shape (248,) into shape (348,)"
- **Fix**: Always create new_calc() for each energy evaluation

### 9. MQA 800K — Li Interdiffusion
- **Problem**: 800K MD causes Li to cross NCM-SE boundary (58/248 atoms)
- **Evidence**: NCM z=-18~+20, SE z=-9~+66 → completely mixed
- **Fix**: No MQA. Static relax only.

### 10. MQA 500K — Still Interdiffusion
- **Problem**: 500K still enough for Li hopping (Ea=0.2 eV << kT=0.043 eV)
- **Evidence**: z_boundary shifts from 17→10.3Å, element-based sep also fails
- **Fix**: No MD at any temperature for adhesion

### 11. SE Thickness Mismatch (2×2×1 vs 2×2×3)
- **Problem**: comp1/2B cubic SE 2×2×1 = 10Å (too thin!) → NCM penetration
- **Evidence**: Gap=-2.3Å, Wad=7.7 J/m²
- **Fix**: comp1/2B SE 2×2×3 = 30Å, comp3/4/5 SE 2×2×1 = 29Å (matched)

### 12. NCM 2L — Surface O Asymmetry
- **Problem**: 2L NCM top O moves freely toward SE, bottom O unconstrained
- **Evidence**: 2L Wad=4~11 J/m² (vs 1L ~1.2), Gap=-3.1Å
- **Physical reason**: Layered oxide vdW between layers → surface reconstruction
- **Fix**: 1L (symmetric constraint) or 5L+FixAtoms (bottom 3L fixed)

### 13. NCM 2L — No Bottom Vacuum
- **Problem**: Initially thought PBC was the issue
- **Analysis**: Bottom vacuum was actually OK (24.9Å for 1L, 44.3Å for 2L)
- **Real cause**: 2L free O has more freedom → over-reconstruction

### 14. NCM 5L + FixAtoms — Correct Approach
- **Solution**: Bottom 3 layers FIXED, top 2 layers FREE
- **Evidence**: Surface O moves only 3-5Å (vs 10Å+ without fix)
- **Standard**: Literature DFT slab method (3-5 layers, bottom fixed)

### 15. Cross-Family Cell Size Mismatch
- **Problem**: cubic (A=351Å²) vs rhombo (A=179Å²) → 2× area difference
- **Evidence**: 5L with different NCM → Li6(2.8) >> Li5.4(1.0) = reversed!
- **Root cause**: SE/A density 1.78 vs 1.38 (29% diff) + NCM size effect
- **Fix attempt**: NCM 7×7 unified (all comps same A=351)
  - comp1/2B: SE 2×2×3 (624at), strain +0.2%
  - comp3/4/5: SE 3×3×1 (558at), strain -5.7%
  - SE/A: 1.78 vs 1.59 (12% diff, improved from 29%)

### 16. Cubic → Rhombo Conversion — Impossible
- **Problem**: comp1 (4 f.u. cubic) cannot be exactly converted to 5 f.u. rhombo
- **Reason**: det(M) = 1.25 = non-integer → no integer supercell transformation
- **Attempts**: pymatgen SpacegroupAnalyzer (P1, symmetry broken), direct M matrix
- **Conclusion**: Different crystal systems = fundamentally different periodicity

### 17. NCM Standalone Relax — Surface Collapse
- **Problem**: Relaxing NCM slab alone → surface O escapes (no constraint from SE)
- **Evidence**: NCM thick=37Å (cell=28Å), O wraps via PBC
- **Fix**: Use pristine NCM (no standalone relax), interface relax only

### 18. v6 — NCM-frozen 500K SE Anneal — NCM Penetration
- **Problem**: Fix NCM, anneal SE at 500K → SE Li hops into frozen NCM lattice
- **Evidence**: 29/248 SE atoms below NCM_zmax (20 Li + 6 S + 3 P)
- **Result**: Wad = 45~67 J/m² (nonphysical, ripped-out atoms on separation)
- **Fix attempt**: Lower to 300K + z-wall (push back atoms below NCM_zmax+0.5Å)

### 19. v6 300K + z-wall — SE Vacuum Ejection
- **Problem**: z-wall prevents downward penetration, but SE atoms fly UPWARD into vacuum
- **Evidence**: SE Li z = 12~**11,080 Å**, S z up to 2,144 Å (cell_z = 75 Å!)
- **Root cause**: UMA gives nonphysical forces at vacuum boundary during MD
- LBFGS controls step size → OK; MD integrates force errors → cascade ejection
- **Conclusion**: z-wall on bottom doesn't fix top ejection. Temperature irrelevant.

### 20. v7 Sandwich (NCM-SE-NCM) — SE Escapes Through Frozen NCM
- **Problem**: NCM_top cap should prevent upward escape, but SE atoms pass through
- **Evidence**: 239/248 SE atoms escaped. SE z = -9,724~12,215 Å
- **Root cause**: FixAtoms ≠ physical wall. Fixed atoms don't repel. UMA force errors
  so large that SE atoms fly through frozen NCM lattice in a few MD steps.
- **Conclusion**: UMA MD fundamentally incompatible with slab/interface geometry.
  Force errors accumulate over timesteps → catastrophic for ANY temperature/constraint.

### 21. UMA MD — Fundamental Limitation (CONCLUSION)
- **Summary**: ALL attempts to use MD at NCM/SE interface fail with UMA:
  | Method | T(K) | Constraint | Result |
  |--------|------|------------|--------|
  | MQA (v2) | 800 | None | Li interdiffusion (58 atoms) |
  | MQA (v5) | 500 | None | Li crosses boundary |
  | v6 freeze-NCM | 500 | FixAtoms NCM | SE penetrates NCM (29 atoms) |
  | v6 + z-wall | 300 | FixAtoms + z-wall | SE ejected to 11 km |
  | v7 sandwich | 500 | NCM cap both sides | SE escapes through NCM (239 atoms) |
- **Root cause**: UMA trained on bulk/dense periodic systems. Interface + vacuum = 
  out-of-distribution → nonphysical forces. LBFGS OK (controlled steps), MD fails 
  (force integration accumulates errors).
- **Final fix**: LBFGS only (v5). No MD at interface with UMA.
- **Alternative (v8)**: Anneal SE in BULK PBC (UMA works) → stack → LBFGS only at interface.

## Current Best Protocol (2026-04-15, updated)

```
NCM: 7×7×5 pristine, bottom 3L FixAtoms
  → surface relax with FixAtoms (fmax=0.005)
  → z offset +15Å (bottom vacuum)
SE: crystalline slab from DFT V0
  → comp1/2B: 2×2×3 (624at, 30Å)
  → comp3/4/5: 3×3×1 (558at, 29Å) [for NCM 7×7 unified]
Interface: NCM + SE, gap=2.5Å, FixAtoms on NCM bottom 3L
  → LBFGS relax (fmax=0.01, steps=200)
Separation: SE +30Å, cell_z +30Å, new_calc(), single point
Wad = (E_sep - E_int) / A × 16.0218
xy-shift sampling: 10-20 seeds per composition
Calculator: UMA (uma-s-1p1), vacuum=30Å always
```
