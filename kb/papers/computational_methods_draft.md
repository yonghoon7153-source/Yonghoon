# Computational Methods — FINAL VERSION (2026-04-16)

## Status: CONFIRMED

First-principles calculations were performed using Quantum ESPRESSO (QE) [1] within the generalized gradient approximation (GGA) parameterized by Perdew–Burke–Ernzerhof (PBE) [2]. Ultrasoft pseudopotentials from the SSSP efficiency library [3] were employed with a plane-wave kinetic energy cutoff of 60 Ry and a charge density cutoff of 480 Ry. Brillouin zone sampling used Γ-centered Monkhorst–Pack grids of 3×3×3 for cubic cells (52 atoms) and 2×2×1 for rhombohedral cells (62 atoms).

**Structure generation and disorder treatment.** Five argyrodite compositions were investigated: Li₆PS₅Cl₁₋ₓBrₓ (x = 0, 0.5; F‾43m, cubic) and halogen-rich Li₅.₄PS₄.₄Cl₁.₆₋ᵧBrᵧ (y = 0.6, 0.8, 1.0; rhombohedral). Structural disorder was treated in two sequential stages. In the first stage, S²⁻/halogen (Cl⁻, Br⁻) distributions over the 4a and 4c Wyckoff sites were enumerated using pymatgen [4], yielding 56–70 symmetry-inequivalent configurations per composition. Each was relaxed with the UMA (Universal Model for Atoms) potential [5], and the lowest-energy halogen arrangement was retained. In the second stage, Li⁺ sublattice disorder on the 48h sites was addressed. Given the intractable combinatorial space (C(48,24) ≈ 10⁹ for Li₆ compositions), 20 configurations were randomly sampled and ranked by UMA energy. The five lowest-energy candidates were subjected to 500 K molecular dynamics (MD) annealing for 50 ps, then quenched to 0 K. This protocol enables efficient sampling of Li⁺ configurations, as Li-ion migration barriers (~0.2 eV) allow for active diffusion at 500 K, while the PS₄³⁻ framework and halogen sublattice remain intact due to their significantly higher bonding energies [6]. The lowest-energy quenched structure was adopted as the representative configuration for each composition. We note that UMA was employed as a surrogate model for configurational sampling and interfacial calculations, while all equilibrium structural properties were obtained from DFT-relaxed structures.

**Equation of state.** Equilibrium volumes and bulk moduli were determined by fitting energy–volume data to the third-order Birch–Murnaghan equation of state [7]. For each representative structure, the unit cell was isotropically scaled to 11 volumes spanning V/V₀ = 0.96–1.06. At each volume, ionic positions were relaxed at fixed cell parameters (force convergence < 10⁻⁴ Ry/Bohr; SCF convergence < 10⁻⁸ Ry). Volume points exhibiting discontinuities associated with basin transitions—identified via displacement analysis of relaxed coordinates between adjacent volumes—were excluded from the fit.

**Work of adhesion.** The work of adhesion (W_ad) between the argyrodite solid electrolyte (SE) and a LiNiO₂ cathode was evaluated using the UMA potential. SE slabs were constructed from DFT-relaxed equilibrium structures: 2×2×3 supercells for Li₆ compositions (624 atoms) and 2×2×1 for Li₅.₄ compositions (248 atoms). Each SE slab was placed onto a single-layer LiNiO₂ slab (5×5×1 or 7×7×1) with an initial gap of 2.5 Å. To construct coherent interfaces, the SE slab was biaxially strained to match the in-plane lattice parameters of the cathode. Interfacial registry was sampled by rigidly translating the SE slab along the xy plane using 20 random fractional displacements per composition [8]; each configuration was relaxed with the L-BFGS optimizer (f_max < 0.01 eV/Å). The work of adhesion was computed as:

W_ad = (E_sep − E_int) / A

where E_int is the energy of the relaxed interface, E_sep is the energy of the same structure with the SE slab displaced by 30 Å along z (with corresponding cell expansion) without further relaxation, and A is the interfacial area. A single cathode layer was employed as a uniform reference surface to isolate the effect of SE composition on interfacial adhesion across all five compositions.

## References

[1] P. Giannozzi et al., "QUANTUM ESPRESSO: a modular and open-source software project for quantum simulations of materials," J. Phys.: Condens. Matter 21, 395502 (2009).
[2] J. P. Perdew, K. Burke, M. Ernzerhof, "Generalized Gradient Approximation Made Simple," Phys. Rev. Lett. 77, 3865 (1996).
[3] G. Prandini et al., "Precision and efficiency in solid-state pseudopotential calculations," npj Comput. Mater. 4, 72 (2018).
[4] S. P. Ong et al., "Python Materials Genomics (pymatgen)," Comput. Mater. Sci. 68, 314 (2013).
[5] B. M. Wood et al., "UMA: A Family of Universal Models for Atoms," arXiv:2506.23971 (2025).
[6] P. Adeli et al., "Boosting Solid-State Diffusivity and Conductivity in Lithium Superionic Argyrodites by Halide Substitution," Angew. Chem. Int. Ed. 58, 8681 (2019).
[7] F. Birch, "Finite Elastic Strain of Cubic Crystals," Phys. Rev. 71, 809 (1947).
[8] A. Sakuda et al., "Sulfide Solid Electrolyte with Favorable Mechanical Property for All-Solid-State Lithium Battery," Sci. Rep. 3, 2261 (2013).
