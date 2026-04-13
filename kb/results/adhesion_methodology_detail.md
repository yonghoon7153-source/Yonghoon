# Adhesion Calculation — Detailed Methodology

## 7-3. Crystalline Slab xy-shift Method

### Step-by-Step Protocol

**Step 1: NCM slab construction**
- LiNiO2 R-3m, a=2.878A, c=14.19A
- NCM 7x7x1 (196 atoms) for Li6, NCM 5x5x1 (100 atoms) for Li5.4
- LBFGS relax (fmax=0.01, 100 steps)

**Step 2: SE slab preparation**
- Load DFT-relaxed crystalline SE from {comp}_slab.xyz
- Orthogonalize: a3 xy components removed (a3=[0,0,z] only)
  - Why: rhombo cell has a3=[-3.558,-2.054,29.047] -> 2x2x1 repeat skews slab
  - Fix: a3=[0,0,29.047] -> clean z-only stacking
- Repeat: 2x2x3 for Li6 (624at, 30A), 2x2x1 for Li5.4 (248at, 29A)
  - Why 2x2x3 for Li6: cubic cell a=10A -> 2x2x1=10A (too thin!) -> NCM penetration
  - 2x2x3=30A matches Li5.4's 29A

**Step 3: xy random shift**
```python
np.random.seed(seed)
dx = np.random.uniform(0, 1)
dy = np.random.uniform(0, 1)
frac_se[:,0] = (frac_se[:,0] + dx) % 1.0
frac_se[:,1] = (frac_se[:,1] + dy) % 1.0
# z is NOT touched! SE slab continuity preserved!
```
- Why xy-shift (not z-shift): z-shift cuts crystalline slab -> dangling bonds
  -> Wad = 0.006~9.252 (nonphysical). xy is PBC-periodic -> shift = registry change only.
- 5 seeds per composition -> 5 different NCM-SE contact patterns

**Step 4: Lattice matching (strain SE to NCM)**
```python
se_cell[0] = ncm_cell[0]
se_cell[1] = ncm_cell[1]
se.set_cell(se_cell, scale_atoms=True)
```
- Li6: strain = +0.2% (SE 2x2 a=20.12A -> NCM 7x7 a=20.15A)
- Li5.4: strain = +1.1% (SE 2x2 a=14.22A -> NCM 5x5 a=14.39A)
- SE strained to match NCM, not vice versa

**Step 5: Stacking (gap = 2.5A)**
- NCM z repositioned: z_min = 0
- SE z repositioned: z_min = NCM_z_max + 2.5A
- Cell z = SE_z_max + 30A (vacuum)
- Why 30A vacuum: UMA gives nonphysical energies with 60A vacuum!
  (tested: Wad=1.252 at 30A vs Wad=24.5 at 60A, same atoms)

**Step 6: LBFGS relax (NO MQA!)**
- LBFGS(fmax=0.01, max_steps=200)
- Why no MQA: 500K causes Li interdiffusion across NCM-SE boundary
  - z_boundary shifts from 17A to 10.3A
  - 58/248 atoms migrate -> structure destroyed
  - Element-based separation (Ni/O=NCM, P/S/Cl=SE, Li=z) also fails
  - Li ownership becomes ambiguous at any temperature > 300K
- Relax-only preserves crystalline structure + vacancy

**Step 7: E_int calculation**
- E_int = interface.get_potential_energy()

**Step 8: Separation + E_sep (single point)**
```python
sep_cell[2] = [0, 0, cell_z + 30.0]  # cell EXPANDS by 30A
pos_sep[n_ncm:, 2] += 30.0           # SE moves UP by 30A
E_sep = sep.get_potential_energy()    # single point, NO relax!
```
- Why cell expansion: without it, PBC wraps SE below NCM
- Why no relax after separation: relax releases SE strain -> Wad overestimated
  (tested: Wad=10 J/m2 with relax vs ~1-2 J/m2 without)
- Wad = (E_sep - E_int) / A * 16.0218 (eV/A2 -> J/m2)

### Key Discoveries During Development

**UMA vacuum sensitivity:**
- UMA (graph neural network MLIP) trained on bulk/dense systems
- Large vacuum (60A) = out of training distribution -> nonphysical energies
- Always use vacuum = 30A for slab calculations with UMA

**z-shift vs xy-shift:**
- z-shift cuts crystalline SE slab at arbitrary planes
  -> dangling bonds, nonphysical surfaces
  -> Wad ranges 0.006~9.252 (meaningless)
- xy-shift preserves SE slab continuity
  -> only NCM-SE registry changes
  -> physical variation from contact pattern

**3000K melt destroys vacancy effect:**
- v2 comp3(Li5.4) with 3000K melt: Wad ≈ 1.0 = same as comp1(Li6)!
- 3000K melts SE -> amorphous -> vacancy structure erased
- Crystalline method essential for Li5.4 vacancy preservation

## 7-4. Results

### Final Values (all crystalline xy-shift)

| Comp | Formula | Wad (J/m2) | std | Family |
|------|---------|-----------|-----|--------|
| comp1 | Li6PS5Cl | **1.433** | 0.288 | Li6 |
| comp2B | Li6PS5Cl0.5Br0.5 | **1.244** | 0.356 | Li6 |
| comp3 | Li5.4PS4.4Cl1.0Br0.6 | **2.361** | 0.41 | Li5.4 |
| comp4 | Li5.4PS4.4Cl0.8Br0.8 | **2.202** | 0.33 | Li5.4 |
| comp5 | Li5.4PS4.4Cl0.6Br1.0 | **2.037** | 0.44 | Li5.4 |

### Trends

**Br effect on adhesion (within family):**
- Li6: comp1(1.43) > comp2B(1.24) -> Br increases -> Wad decreases
- Li5.4: comp3(2.36) > comp4(2.20) > comp5(2.04) -> monotonic decrease with Br

**Vacancy effect on adhesion (cross-family):**
- Li5.4(~2.2 J/m2) >> Li6(~1.3 J/m2) -> vacancy nearly doubles adhesion
- Conservative estimate: Li5.4 has lower SE/A density (1.38 vs 1.78 at/A2)
  yet higher Wad -> vacancy effect even stronger than raw ratio

**Physical mechanism:**
- Li6 (no vacancy): surface Li fully coordinated -> no driving force for NCM bonding
- Li5.4 (vacancy): under-coordinated Li near vacancy = "chemical anchor"
  -> pulls NCM O2- across interface -> new cross-interface ionic bonds
  -> Wad enhanced by ~1.0 J/m2

**Experimental consistency:**
- Experiment: Li5.4(316 aJ) >> Li6(194 aJ) at r=10nm -> same trend!
- Br effect: within each family, Br increases -> Wad decreases in both calc and expt

### 3000K Melt Validation
v2(3000K melt) comp3: Wad ≈ 1.0 J/m2 = same as comp1(1.1)
-> Proves: amorphous SE erases vacancy effect
-> The ~1.2 J/m2 difference (v5 - v2 for Li5.4) = vacancy contribution quantified
