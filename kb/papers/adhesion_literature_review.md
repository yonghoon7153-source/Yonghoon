# Adhesion Interface Modeling — Literature Summary

## Reference Papers for Cathode-SE Interface Adhesion

### 1. Xiao et al., Nature Reviews Materials 5, 105 (2020)
"Understanding interface stability in solid-state batteries"
- **Review**: Comprehensive overview of interface modeling in SSBs
- **Wad**: Work of adhesion as descriptor for heterogeneous interfacial stability
- **Method**: Ab initio (DFT) + thermodynamic analysis
- **Key point**: Computational models successfully predict SE stability; Wad quantifies adhesion
- DOI: 10.1038/s41578-019-0157-5

### 2. Schwietert et al., ACS Appl. Mater. Interfaces 15, 39, 45892 (2023)
"Electrolyte Coatings for High Adhesion Interfaces in Solid-State Batteries from First Principles"
- **Method**: Adhesion parameter from single-material slab cleavage energies (no direct interface calc!)
- **Slab**: Individual material slabs, surface energy calculation
- **Database**: 19,481 Li compounds screened, 945 slab terminations
- **Formula**: Adapted contact angle equation with cleavage energies
- **Code**: VASP, PBE, convergence issues for some large slabs
- DOI: 10.1021/acsami.3c04452
- PMC: PMC10520915

### 3. Haruyama et al., Chem. Mater. 26, 4248 (2014)
"Space-Charge Layer Effect at Interface between Oxide Cathode and Sulfide Electrolyte"
- **Materials**: LiCoO2 / β-Li3PS4
- **Slab**: 3-4 atomic layers per material
- **Vacuum**: ~15 Å
- **Method**: DFT+U (U=3.3 eV for Co 3d)
- **Fixed atoms**: Bottom layers fixed to simulate bulk
- **Key finding**: Space-charge layer forms at interface, Li depletion in SE side
- DOI: 10.1021/cm5016959

### 4. Lian et al., ACS Appl. Energy Mater. 3, 10, 9889 (2020)
"Comparative Study on Sulfide and Oxide Electrolyte Interfaces with Cathodes"
- **Materials**: LCO/LPS, LCO/Li3PO4, LCO/LLZO
- **Slab**: 3-5 atomic layers, systematic interface matching
- **Formula**: Wad = (E_slab1 + E_slab2 - E_interface) / A
- **Method**: DFT+U, direct interface calculation
- **Key finding**: All Wad positive (stable), unlike Li-metal/SE interfaces
- DOI: 10.1021/acsaem.0c02033

### 5. Meng et al., JACS 144, 43, 19972 (2022)
"Thermodynamics and Kinetics of the Cathode-Electrolyte Interface in All-Solid-State Li-S Batteries"
- **Materials**: Li2S cathode / SE interface
- **Method**: MTP (Moment Tensor Potential) MLIP
- **Slab**: >1000 atoms, large-scale interface
- **MD**: NPT, 600K, 5 ns, active learning approach
- **Key finding**: Most stable surface → lowest Li diffusion barrier at interface
- DOI: 10.1021/jacs.2c07482

### 6. Kim et al., J. Phys. Chem. C 126, 36, 15326 (2022)
"Interfacial Stability of Layered LiNixMnyCo1-x-yO2 Cathodes with Sulfide SE"
- **Materials**: NMC111/333/532/622/811 with Li6PS5Cl, Li3PS4
- **Method**: DFT+U
- **Slab**: Layered oxide cathode slabs, systematic study across compositions
- **Key finding**: Higher Ni content → less stable interface
- DOI: 10.1021/acs.jpcc.2c05336

### 7. Seymour et al., ChemSusChem 16, e202202215 (2023)
"Understanding and Engineering Interfacial Adhesion in Solid-State Batteries with Metallic Anodes"
- **Focus**: Li metal anode / SE adhesion (not cathode, but methodology relevant)
- **Method**: DFT slab models
- **Slab**: Multiple layers, bottom fixed
- **Key insight**: Adhesion depends on interface registry and surface termination
- DOI: 10.1002/cssc.202202215

### 8. Richards et al., Chem. Mater. 28, 266 (2016)
"Interface Stability in Solid-State Batteries"
- **Mo group** pioneering work on interface stability
- **Method**: Thermodynamic analysis + DFT
- **Approach**: Phase diagram-based interface stability prediction
- **Key**: Mutual decomposition energy as interface stability descriptor
- Note: Not slab-based but important reference for thermodynamic approach

---

## Standard DFT Slab Methodology Summary

### Typical Parameters in Literature
| Parameter | Range | Most Common |
|-----------|-------|-------------|
| Cathode layers | 3-6 | 4-5 |
| SE layers | 3-6 | 4-5 |
| Bottom fixed | 1-3 layers | 2 layers |
| Top free | 1-3 layers | 2 layers |
| Vacuum | 10-20 Å | 15 Å |
| DFT functional | PBE, PBE+U | PBE+U for oxides |
| Code | VASP, QE | VASP dominant |

### Wad Formulas Used
1. **Direct interface**: Wad = (E_slab1 + E_slab2 - E_interface) / A
   - Standard DFT approach
   - Requires three separate calculations
   - E_slab1, E_slab2 = isolated slab energies (same cell, relaxed)

2. **Separation method**: Wad = (E_separated - E_interface) / A
   - Used when isolated slab energies unreliable (e.g., MLIP)
   - Single interface calculation + rigid separation
   - Our method (UMA compatibility)

3. **Cleavage energy**: Wad ≈ f(γ_1, γ_2)
   - No interface calculation needed
   - Uses surface energies of individual materials
   - Schwietert et al. 2023 approach

### FixAtoms Standard
```
Typical 5-layer slab:
  Layer 1 (bottom): FIXED
  Layer 2:          FIXED
  Layer 3:          FREE or FIXED
  Layer 4:          FREE
  Layer 5 (top):    FREE

ASE implementation:
  from ase.constraints import FixAtoms
  constraint = FixAtoms(indices=[i for i,p in enumerate(atoms.positions) if p[2] < z_cutoff])
  atoms.set_constraint(constraint)

VASP implementation:
  Selective dynamics = T
  Bottom atoms: F F F
  Top atoms: T T T
```

### Key Physical Insights from Literature
1. **Layered oxide cathodes**: Surface O tends to reconstruct → need fixed bottom layers
2. **Sulfide SE**: PS4 tetrahedra are rigid → less surface reconstruction
3. **Interface Li**: Fast diffusion (Ea~0.2 eV) → any MD causes Li mixing
4. **Registry dependence**: Wad varies with lateral alignment (our xy-shift addresses this)
5. **Composition dependence**: Higher Ni → less stable interface (Kim 2022)
