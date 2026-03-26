# LIGGGHTS DEM Analysis Skill (Standard Mode)

## Overview

Comprehensive DEM simulation analysis toolkit for LIGGGHTS output files. Designed for solid-state battery composite cathode research with **AM (Active Material) + SE (Solid Electrolyte)** systems.

## Trigger Keywords

- LIGGGHTS, DEM simulation, atom dump, contact dump
- NCM, SE, solid-state battery, composite cathode
- coordination number, percolation, force chain
- Von Mises stress, fabric tensor, porosity
- DEM validity

## Contact Types (3 Types)

| Type | Description |
|------|-------------|
| AM-AM | Active material contacts (force chain backbone) |
| AM-SE | Interface contacts (electrochemical activity) |
| SE-SE | Ionic conduction network |

## Analysis Pipeline

```bash
# Step 1: Parse LIGGGHTS files
python scripts/parse_liggghts.py atom_*.liggghts contact_*.liggghts -o ./results

# Step 2: Basic Contact Analysis
python scripts/analyze_contacts.py results/atoms.csv results/contacts.csv \
    -o ./results -t "1:AM,2:SE" -s 1000

# Step 3: Basic Figures (Fig 1-4)
python scripts/generate_figures.py ./results -o ./figures -s 1000

# Step 4: Advanced Analysis + Figures (Fig 5-9)
python scripts/advanced_analysis.py results/atoms_analyzed.csv results/contacts_analyzed.csv \
    -o ./results -s 1000
python scripts/generate_advanced_figures.py ./results -o ./figures -s 1000

# Step 5: Thick Electrode Analysis (Fig 10-17)
python scripts/thick_electrode_analysis.py results/atoms_analyzed.csv results/contacts_analyzed.csv \
    -o ./results -n 10 -s 1000
python scripts/deep_electrode_analysis.py results/atoms_analyzed.csv results/contacts_analyzed.csv \
    -o ./results -s 1000
```

## Output Figures

### Basic Analysis (4 figures)

| Figure | Content | Key Metrics |
|--------|---------|-------------|
| fig1_comprehensive | Particle distribution, contact types, Z-profile | Type counts, contact breakdown |
| fig2_contact_network | Force vs Area, Coord vs Z, Force direction heatmap, Overlap vs Contact Area | R², Hertzian fit |
| fig3_ionic_conductivity | SE network analysis for ion transport | SE-SE coordination, percolation |
| fig4_mechanical | Force distribution and contact mechanics | Normal/tangential forces |

### Advanced Analysis (5 figures)

| Figure | Content | Key Metrics |
|--------|---------|-------------|
| fig5_stress_analysis | Stress tensor analysis | Von Mises, hydrostatic pressure |
| fig6_force_chains | Force chain visualization | Z-alignment, chain fraction |
| fig7_fabric_tensor | Structural anisotropy | Eigenvalues, anisotropy index |
| fig8_porosity_network | Packing, connectivity, percolation | NCM 99.9% vs SE 0% |
| fig9_torque_rolling | Contact dynamics | Rolling/sliding ratio, friction |

### Thick Electrode Analysis (3 figures)

| Figure | Content |
|--------|---------|
| thick_electrode_fig1 | Height profile overview |
| thick_electrode_fig2 | Top vs Bottom comparison |
| thick_electrode_fig3 | SE network by height |

### Deep Analysis (5 figures)

| Figure | Content |
|--------|---------|
| deep_fig1_am_utilization | AM surface coverage, dead AM |
| deep_fig2_contact_quality | Overlap analysis, Hertzian validation |
| deep_fig3_se_bottlenecks | SE network weak points, isolated SE |
| deep_fig4_boundary_heterogeneity | Wall effects, center-edge radial profiles |
| deep_fig5_force_network | Force transmission by type |

## Key Metrics & Thresholds

### Ionic Conductivity

| Metric | Threshold | Interpretation |
|--------|-----------|----------------|
| SE-SE Coordination | > 2.4 | Percolation threshold for 3D |
| SE Percolation | > 70% | Good ion pathway connectivity |
| NCM-SE Interface | > 80% | Active material utilization |

### Mechanical Properties

| Metric | Good | Warning |
|--------|------|---------|
| Overlap ratio | < 3% | > 5% |
| Hertzian correlation | > 0.9 | < 0.8 |
| Rolling fraction | 20-40% | > 80% |

### Thick Electrode

| Metric | Typical | Concern |
|--------|---------|---------|
| Top-Bottom packing diff | < 5% | > 10% |
| SE bottleneck fraction | < 20% | > 30% |
| AM surface coverage | 5-15% | < 2% |

## DEM Validity Assessment

### Highly Valid (Topology preserved after SE deformation)

- Contact network topology
- Coordination number distribution
- Force chain pathways
- Fabric tensor (structural anisotropy)
- Layer-wise uniformity (CV)

### Reference Only (Absolute values limited)

- Porosity (~35%): sphere packing geometric limit. Use CV for uniformity comparison.
- Packing density: does not reflect SE deformation. Use for relative layer comparison.

### Not Predictive

- Absolute electrode density (SE deformation not modeled)
- Direct ionic conductivity (depends on actual porosity)

### Porosity Gap

DEM porosity (~35%) vs Real electrode (<5%). Random Close Packing limit is ~36%. Real SE fills voids through plastic deformation.

### NCM vs SE Network

| Network | Percolation | Interpretation |
|---------|-------------|----------------|
| NCM | 99.9% | Electron pathway fully connected |
| SE | ~0% | Requires plastic deformation for ionic conduction |

## Scale Factor

- Default: 1000x (simulation uses mm, real is um)
- AM radius: 4 mm (sim) -> 4 um (real)
- SE radius: 0.8 mm (sim) -> 0.8 um (real)

## Dependencies

```
pandas, numpy, matplotlib, networkx
```

## Post-Processing Checklist

- [ ] Move all PNG files to figures folder
- [ ] Generate DEM_Analysis_Summary.md
- [ ] Copy results to outputs directory

## Summary MD Required Sections

1. System Overview (particle count, contact count, electrode height)
2. Contact Distribution (SE-SE, AM-SE, AM-AM ratios)
3. Coordination Numbers (AM, SE avg/min/max)
4. Ionic Conductivity Metrics (SE percolation, coordination)
5. Force Statistics (avg/max force)
6. Contact Quality (overlap, Hertzian correlation)
7. Key Findings
8. Generated Figures list
