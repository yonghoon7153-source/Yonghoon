# EOS Fitting — Birch-Murnaghan Equation of State

## What It Does
Uniformly compress/expand the lattice, compute E(V) curve, fit 3rd-order
Birch-Murnaghan EOS to extract B0 (bulk modulus), B0' (pressure derivative), V0.

## Theory

Eulerian finite strain: f = (1/2)[(V0/V)^(2/3) - 1]

3rd order BM EOS:
```
E(V) = E0 + (9*V0*B0/16) * {[(V0/V)^(2/3)-1]^2
       x [6 + B0'*((V0/V)^(2/3)-1) - 4*(V0/V)^(2/3)]}
```

Parameter meanings:
- **E0**: Equilibrium energy (parabola minimum)
- **V0**: Equilibrium volume (minimum position)
- **B0**: Bulk modulus = parabola curvature = compression resistance (GPa)
- **B0'**: dB/dP = parabola asymmetry (B0'=4 is symmetric; >4 stiffens under pressure)

## Protocol

1. Start from MLIP-relaxed or DFT-relaxed champion structure
2. Scale volume uniformly: 96%, 97%, ..., 106% (skip 107-108% if basin transition)
3. At each volume: DFT relax (cell fixed, ions free) — `calculation='relax'`, `cell_dofree='none'`
4. Collect E(V) data points
5. Fit using scipy `curve_fit` with BM3 function
6. Quality: R^2 > 0.9999 required

## Basin Transition Warning

Volume expansion >+6% (v107, v108) can trigger basin transitions:
- comp5 v108: 31/62 atoms rearranged (full reconstruction)
- Model C v106: 7/62 atoms (local Li hop), v108: 28/62 atoms
- **Always check**: compare relaxed coordinates across volumes
- **Exclude** points where displacement analysis shows >2 atoms moving significantly

## Typical Values for Argyrodites

| Comp | B0 (GPa) | B0' | Notes |
|------|----------|-----|-------|
| comp1 | 26.2 | 4.17 | Li6, Cl only |
| comp2 | 25.8 | 4.22 | Li6, Cl+Br |
| comp3 | 20.8 | 6.79 | Li5.4, Cl-rich |
| comp4 | 20.8 | 6.33 | Li5.4, equal Cl/Br |
| comp5 | 22.9 | 5.21 | Li5.4, Br-rich (Basin A) |
| Model C | 21.7 | 4.31 | Li5.4, Cl only |

## Key References
- [birch1947] Birch, Phys. Rev. 71, 809 (1947) — Original BM EOS
- [minerals2019] Minerals 9, 745 (2019) — BM EOS tutorial
- [deng2016] Deng et al., J. Electrochem. Soc. 163, A67 (2016) — LPSCl DFT elastic
