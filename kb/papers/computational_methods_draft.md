# Computational Methods — Paper Draft Section

## DFT Methodology

First-principles calculations were performed using Quantum ESPRESSO (QE) [1] within the generalized gradient approximation (GGA) parameterized by Perdew–Burke–Ernzerhof (PBE) [2]. Ultrasoft pseudopotentials (USPP) from the SSSP efficiency library [3] were employed with a kinetic energy cutoff of 60 Ry and a charge density cutoff of 480 Ry. Brillouin zone integration was performed using Γ-centered Monkhorst–Pack grids of 3×3×3 for cubic cells (52 atoms) and 2×2×1 for rhombohedral cells (62 atoms).

## Structure Generation

The argyrodite Li₆PS₅Cl (F‾43m) and halogen-rich Li₅.₄PS₄.₄Cl₁.₆₋ₓBrₓ structures were constructed following a two-stage disorder treatment. First, the S²⁻/halogen (Cl⁻, Br⁻) distribution over the 4a and 4c Wyckoff sites was enumerated using pymatgen [4], generating 56–70 symmetry-inequivalent configurations per composition. These were screened via machine-learning interatomic potential (MLIP) relaxation to identify the lowest-energy halogen arrangement. Second, the Li⁺ sublattice disorder on the 48h sites—which gives rise to C(48,24) ≈ 10⁹ possible arrangements for Li₆ compositions—was addressed by randomly sampling 20 configurations and ranking them by MLIP energy. The five lowest-energy candidates were then subjected to 500 K molecular dynamics (MD) annealing for 50 ps using the MACE universal potential [5], followed by quenching to 0 K. This thermal treatment allows Li⁺ ions (activation energy ~0.2 eV) to hop freely while preserving the PS₄³⁻ framework (P–S bond energy ~3.5 eV) and halogen positions [6]. The lowest-energy annealed structure was selected as the champion configuration for each composition.

## Equation of State

Equilibrium volumes and bulk moduli were determined by fitting the energy–volume relation to the third-order Birch–Murnaghan equation of state [7]. For each composition, the champion structure was uniformly scaled to 11 volumes spanning V/V₀ = 0.96–1.06. At each volume, ionic positions were relaxed with fixed cell parameters (force convergence < 10⁻⁴ Ry/Bohr, SCF convergence < 10⁻⁸ Ry). Volume points exhibiting basin transitions—identified by displacement analysis of relaxed coordinates across adjacent volumes—were excluded from the fit.

## Elastic Constants

Finite-temperature elastic constants were evaluated using a 600 K snapshot method with the MACE potential. A 2×2×1 supercell was equilibrated at 600 K via Langevin MD (timestep 2 fs, 30 ps), from which five equally spaced snapshots were extracted, quenched to 0 K, and relaxed. Elastic tensors were computed for each snapshot using the finite-strain approach (±0.5% strain, six independent deformations), and polycrystalline moduli were obtained via the Voigt–Reuss–Hill averaging scheme [8]. This approach naturally incorporates the thermal Li⁺ disorder present under experimental conditions, resolving the systematic overestimation of shear constants (C₄₄) inherent in ordered 0 K calculations [9,10].

## Interfacial Adhesion

The work of adhesion (Wad) between argyrodite SE and layered cathode (LiNiO₂) was evaluated using the Universal Materials Accelerator (UMA) potential [11]. Crystalline SE slabs were constructed from the DFT-relaxed equilibrium structures: 2×2×3 supercells for Li₆ compositions (624 atoms) and 2×2×1 for Li₅.₄ compositions (248 atoms). These were stacked onto a single-layer LiNiO₂ slab (5×5×1 or 7×7×1) with a 2.5 Å initial gap after applying lateral strain to match the cathode lattice. To generate statistically meaningful sampling, the SE slab was rigidly shifted along the xy plane by random fractional displacements for 20 independent seeds per composition, preserving the slab integrity while varying the interfacial contact registry [12]. Each configuration was relaxed using the L-BFGS optimizer (fmax < 0.01 eV/Å). The adhesion energy was computed via the separation method:

Wad = (E_sep − E_int) / A

where E_int is the energy of the relaxed interface, E_sep is the single-point energy after displacing the SE slab by 30 Å along z (with corresponding cell expansion), and A is the interfacial area. A single NCM layer was intentionally employed to minimize the cathode surface contribution, thereby isolating the SE composition effect on interfacial adhesion across five compositions under consistent conditions.

## References

[1] P. Giannozzi et al., "QUANTUM ESPRESSO: a modular and open-source software project for quantum simulations of materials," J. Phys.: Condens. Matter 21, 395502 (2009).
[2] J. P. Perdew, K. Burke, M. Ernzerhof, "Generalized Gradient Approximation Made Simple," Phys. Rev. Lett. 77, 3865 (1996).
[3] G. Prandini et al., "Precision and efficiency in solid-state pseudopotential calculations," npj Comput. Mater. 4, 72 (2018).
[4] S. P. Ong et al., "Python Materials Genomics (pymatgen)," Comput. Mater. Sci. 68, 314 (2013).
[5] I. Batatia et al., "MACE: Higher Order Equivariant Message Passing Neural Networks," NeurIPS (2022).
[6] P. Adeli et al., "Boosting Solid-State Diffusivity and Conductivity in Lithium Superionic Argyrodites by Halide Substitution," Angew. Chem. Int. Ed. 58, 8681 (2019).
[7] F. Birch, "Finite Elastic Strain of Cubic Crystals," Phys. Rev. 71, 809 (1947).
[8] R. Hill, "The Elastic Behaviour of a Crystalline Aggregate," Proc. Phys. Soc. A 65, 349 (1952).
[9] Z. Deng et al., "Elastic Properties of Alkali Superionic Conductor Electrolytes from First Principles Calculations," J. Electrochem. Soc. 163, A67 (2016).
[10] M. Torii et al., "First-Principles Investigation of Mechanical Properties and Anisotropy of Argyrodite Li₆PS₅Cl Crystal Electrolytes," J. Phys. Chem. C 129, 17882 (2025).
[11] Meta FAIR, "Universal Materials Accelerator (UMA)," fairchem (2024).
[12] A. Sakuda et al., "Sulfide Solid Electrolyte with Favorable Mechanical Property for All-Solid-State Lithium Battery," Sci. Rep. 3, 2261 (2013).

## Actual K-grid Used (for reference)

| Comp | Cell | Atoms | K-grid (SCF/relax) | Pipeline |
|------|------|-------|-------------------|----------|
| comp1 | cubic | 52 | 3×3×3 | v1, v2 |
| comp2 | cubic | 52 | 3×3×3 | v1 |
| comp3 | rhombo | 62 | 2×2×1 | v1 |
| comp4 | rhombo | 62 | 2×2×1 | v1 |
| comp5 | rhombo | 62 | 2×2×1 (v1), 3×3×3 (basin) | v1 |
| Model C | rhombo | 62 | 2×2×1 (v1), 6×6×3 (v2) | v1, v2 |
