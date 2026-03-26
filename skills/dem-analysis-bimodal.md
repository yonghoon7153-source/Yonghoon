# Bimodal DEM Analysis Skill

## Overview

Specialized DEM analysis for **bimodal active material particle systems** in solid-state battery composite cathode research. Handles systems with:
- **AM_P** (Primary/Large AM): Larger NCM particles
- **AM_S** (Secondary/Small AM): Smaller NCM particles
- **SE** (Solid Electrolyte): Fine SE particles filling gaps

## Trigger Keywords

- bimodal, two AM types, AM_P, AM_S
- primary AM, secondary AM, large AM, small AM
- size distribution, particle size ratio
- NCM bimodal, cathode bimodal

## Contact Types (6 Types)

| Type | Description |
|------|-------------|
| AM_P-AM_P | Large-Large AM contact |
| AM_P-AM_S | Large-Small AM contact (bimodal key!) |
| AM_P-SE | Large AM - SE interface |
| AM_S-AM_S | Small-Small AM contact |
| AM_S-SE | Small AM - SE interface |
| SE-SE | SE network for ionic conduction |

## Analysis Pipeline

```bash
# Step 1: Parse LIGGGHTS files
python scripts/parse_liggghts.py atom_*.liggghts contact_*.liggghts -o ./results

# Step 2: Bimodal Contact Analysis
python scripts/analyze_contacts_bimodal.py results/atoms.csv results/contacts.csv \
    -o ./results -t "1:AM_P,2:AM_S,3:SE" -s 1000

# Step 3: Basic Figures (Fig 1-4)
python scripts/generate_figures_bimodal.py ./results -o ./figures -s 1000

# Step 4: Advanced Analysis + Figures (Fig 5-9)
python scripts/advanced_analysis_bimodal.py results/atoms_analyzed.csv results/contacts_analyzed.csv \
    -o ./results -s 1000
python scripts/generate_advanced_figures_bimodal.py ./results -o ./figures -s 1000

# Step 5: Bimodal-Specific Analysis (Fig 10-11)
python scripts/bimodal_specific_analysis.py results/atoms_analyzed.csv results/contacts_analyzed.csv \
    -o ./results -s 1000

# Step 6: DEM Validity & SE Virtual Expansion (Fig S1)
python scripts/generate_figS1_bimodal.py --results ./results --output ./figures --scale 1000
```

## Output Figures

### Basic Analysis (4 figures)

| Figure | Content |
|--------|---------|
| fig1_comprehensive | AM_P/AM_S/SE distribution, 6 contact types |
| fig2_contact_network | Contact mechanics by all 6 types |
| fig3_ionic_conductivity | SE network analysis |
| fig4_mechanical | Force distribution |

### Advanced Analysis (5 figures)

| Figure | Content |
|--------|---------|
| fig5_stress_analysis | Stress tensor |
| fig6_force_chains | Force chain visualization |
| fig7_fabric_tensor | Structural anisotropy |
| fig8_porosity_network | Packing, connectivity |
| fig9_torque_rolling | Contact dynamics |

### Bimodal-Specific (2 figures)

| Figure | Content |
|--------|---------|
| fig10_bimodal_am_comparison | AM_P vs AM_S comparison |
| fig11_bimodal_packing | Bimodal packing efficiency |

### Supplementary (1 figure)

| Figure | Content |
|--------|---------|
| FigS1_DEM_validity | DEM validity & SE expansion |

## Performance Optimizations

| Operation | Before | After | Speedup |
|-----------|--------|-------|---------|
| fig11 SE shell analysis | O(AM x SE) nested loops | KDTree range query | 100-1000x |
| fig3 layer distribution | iterrows() + lookup | Vectorized mapping | 50-100x |
| Network analysis | iterrows() | Bulk edge addition | 10-50x |

## Scale Factor

- Default: 1000x (simulation uses mm, real is um)
- AM_P radius: 4-5 mm (sim) -> 4-5 um (real)
- AM_S radius: 2-2.5 mm (sim) -> 2-2.5 um (real)
- SE radius: 0.5-0.8 mm (sim) -> 0.5-0.8 um (real)

## Dependencies

```
pandas, numpy, matplotlib, scipy (for KDTree), networkx
```

## Performance Tips

1. **Large SE datasets (>50,000 particles)**: KDTree-based analysis handles this efficiently (~30 seconds for fig11)
2. **Memory**: KDTree requires ~8 bytes per coordinate per particle (~2.4 MB for 100k particles)
