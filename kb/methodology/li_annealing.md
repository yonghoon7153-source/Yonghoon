# Li Annealing — Thermal Li Sublattice Re-optimization

## Purpose
Re-optimize Li positions after halogen substitution without disturbing the
PS₄ framework or halogen placement. This is the key improvement of Pipeline v2 over v1.

## Why It's Needed

The two-stage disorder problem:
1. **Halogen disorder**: 4a/4c S/Cl/Br placement → 56-70 configs → enumerable
2. **Li disorder**: 48h site half-occupancy → C(48,24)=12 billion → NOT enumerable

Pipeline v1 uses Rietveld Li positions (optimized for original composition, not Br-substituted).
Pipeline v2 adds Li annealing to find the optimal Li arrangement for each halogen config.

## Protocol

1. From halogen screening, take top 5 Li configurations
2. MLIP MD at **500K** for 50-100 ps (Langevin thermostat, dt=2fs)
3. After MD, MLIP relax (FIRE optimizer)
4. Compare energies → select champion

## Temperature Selection: 500K

| Species | Ea (eV) | kT at 500K | Behavior |
|---------|---------|------------|----------|
| Li⁺ | ~0.2 | 0.043 | Active hopping ✅ |
| S²⁻ | ~3.5 (P-S) | 0.043 | Vibration only ✅ |
| Cl⁻/Br⁻ | ~1.0 | 0.043 | Frozen ✅ |

- **<500K**: Li hopping too slow, stuck in local minima
- **500K**: Optimal — Li freely explores, framework preserved
- **800K+**: Cl⁻/Br⁻ start hopping → halogen enumeration invalidated ⚠️
- **1500K+**: PS₄ framework can break ⚠️

## Key Finding: Ranking Reversal

comp1 (Li₆PS₅Cl) results:

| Config | Screening Energy (eV) | Screening Rank | Annealing Energy (eV) | Annealing Rank |
|--------|----------------------|----------------|----------------------|----------------|
| #0 | -217.468 | 1 | -217.533 | **1** |
| #1 | -217.283 | 2 | -216.953 | **5** ↓ |
| #8 | -217.200 | 3 | -217.030 | **3** |
| #15 | -217.175 | 4 | -217.042 | **2** ↑ |
| #9 | -217.131 | 5 | -217.013 | **4** |

**Interpretation:** 0K relaxation (screening) finds local minima. Thermal annealing
explores deeper basins. A config with high screening energy may sit near a deep basin
that 0K relaxation cannot reach.

## Energy Spread

Li₆PS₅Cl: 1162 meV spread across 20 random Li configurations.
This means Li arrangement can change the total energy by >1 eV — comparable to
the energy difference between compositions!

## Impact on Downstream Properties

| Property | Pipeline v1 | Pipeline v2 | Delta |
|----------|-------------|-------------|-------|
| comp1 B0 (MLIP) | — | 26.9 GPa | — |
| comp1 B0 (DFT v1) | 26.2 | TBD | Expected small (Li6 = no vacancy) |
| comp5 Basin problem | Severe | Expected reduced | — |

For Li₆ compositions (no vacancy), the effect is small (~0.7 GPa).
For Li₅.₄ compositions (with vacancy), the effect is expected to be larger
because vacancy creates more distinct local environments for Li.

## References
- [pustorino2025] Systematic Li structure sampling in LPSCl
- [ayadi2024] AIMD Li short-range order, temperature homogenization
