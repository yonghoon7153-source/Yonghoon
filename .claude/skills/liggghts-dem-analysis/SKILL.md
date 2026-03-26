---
name: liggghts-dem-analysis
description: Comprehensive DEM (Discrete Element Method) simulation analysis for LIGGGHTS output files. Use when user provides LIGGGHTS dump files (atom_.liggghts, contact_.liggghts) or input scripts for solid-state battery composite cathode analysis. Supports NCM-SE composite analysis including contact classification, coordination number calculation, force distribution, ionic conductivity proxy metrics, stress tensor analysis, force chain identification, fabric tensor calculation, porosity analysis, network connectivity metrics, and DEM validity assessment. Generates publication-quality main figures. Triggers on keywords like LIGGGHTS, DEM simulation, atom dump, contact dump, NCM, SE, solid-state battery, composite cathode, coordination number, percolation, Von Mises stress, force chain, fabric tensor, porosity, DEM validity, 후막전극, thick electrode.
---

# LIGGGHTS DEM Analysis Skill

## MANDATORY WORKFLOW (반드시 준수)

이 skill을 사용할 때는 아래 전체 워크플로우를 항상 완료해야 합니다.

### Pre-Processing
- [ ] Contact 파일이 여러 part로 분할된 경우 → 병합 (part1 전체 + part2,3... 데이터만 추가)
- [ ] Atom/Contact 파일을 작업 디렉토리로 복사

### Analysis Pipeline (순차 실행 필수)

```bash
# Step 1: Parse
python parse_liggghts.py atom_*.liggghts contact_*.liggghts -o ./results

# Step 2: Basic Analysis
python analyze_contacts.py results/atoms.csv results/contacts.csv -o ./results -t "1:AM,2:SE" -s 1000

# Step 3: Basic Figures (Fig 1-4)
python generate_figures.py ./results -o ./figures -s 1000

# Step 4: Advanced Analysis + Figures (Fig 5-9)
python advanced_analysis.py results/atoms_analyzed.csv results/contacts_analyzed.csv -o ./results -s 1000
python generate_advanced_figures.py ./results -o ./figures -s 1000

# Step 5: Thick Electrode Analysis (Fig 10-17)
python thick_electrode_analysis.py results/atoms_analyzed.csv results/contacts_analyzed.csv -o ./results -n 10 -s 1000
python deep_electrode_analysis.py results/atoms_analyzed.csv results/contacts_analyzed.csv -o ./results -s 1000
```

### Post-Processing (필수)
- [ ] 모든 PNG 파일을 figures 폴더로 이동
- [ ] 분석 요약 마크다운 파일 생성 (DEM_Analysis_Summary.md)
- [ ] 모든 결과 파일을 outputs 디렉토리로 복사
- [ ] present_files로 사용자에게 제공

### Output Checklist

| Category | Files | Count |
|----------|-------|-------|
| Basic Figures | fig1~fig4 | 4 |
| Advanced Figures | fig5~fig9 | 5 |
| Thick Electrode | thick_electrode_fig1~3 | 3 |
| Deep Analysis | deep_fig1~5 | 5 |
| Summary | DEM_Analysis_Summary.md | 1 |
| **Total** | | **18 files** |

### Summary MD 필수 포함 항목

- System Overview (입자 수, 접촉 수, 전극 높이)
- Contact Distribution (SE-SE, AM-SE, AM-AM 비율)
- Coordination Numbers (AM, SE 평균/최소/최대)
- Ionic Conductivity Metrics (SE percolation, coordination)
- Force Statistics (평균/최대 힘)
- Contact Quality (overlap, Hertzian correlation)
- Key Findings (긍정적 / 모니터링 필요)
- Generated Figures 목록

---

## Overview

Comprehensive DEM (Discrete Element Method) simulation analysis toolkit for LIGGGHTS output files, specifically designed for solid-state battery composite cathode research. Supports both single multiparticle compression (multisphere) and thick film electrode (후막전극) analysis.

## Trigger Keywords

- LIGGGHTS, DEM simulation, atom dump, contact dump
- NCM, SE, solid-state battery, composite cathode
- 후막전극, thick electrode, thick film
- multisphere, single particle compression
- coordination number, percolation, force chain
- Von Mises stress, fabric tensor, porosity
- DEM validity, NCM network, SE network, ionic pathway

## Analysis Modes

### 1. Multisphere Mode (Single Particle Compression)
For analyzing compression of individual multiparticle aggregates.
- Typical scale: ~2,000-10,000 particles
- Focus: Particle-level stress, deformation, internal structure

### 2. Thick Electrode Mode (후막전극)
For analyzing full electrode layer compression.
- Typical scale: ~200,000+ particles
- Focus: Height gradients, top-bottom differences, ionic pathways
- Electrode height: typically 100-300 μm

## Quick Start

### Step 1: Parse LIGGGHTS Files
```bash
python scripts/parse_liggghts.py atom_*.liggghts contact_*.liggghts -o ./results
```

### Step 2: Basic Contact Analysis
```bash
python scripts/analyze_contacts.py results/atoms.csv results/contacts.csv -o ./results -t "1:AM,2:SE" -s 1000
```

### Step 3: Generate Basic Figures (Fig 1-4)
```bash
python scripts/generate_figures.py ./results -o ./figures
```

### Step 4: Advanced Analysis (Fig 5-9)
```bash
python scripts/advanced_analysis.py results/atoms_analyzed.csv results/contacts_analyzed.csv -o ./results -s 1000
python scripts/generate_advanced_figures.py ./results -o ./figures -s 1000
```

### Step 5: Thick Electrode Analysis (후막전극 전용)
```bash
python scripts/thick_electrode_analysis.py results/atoms_analyzed.csv results/contacts_analyzed.csv -o ./results -n 10 -s 1000
python scripts/deep_electrode_analysis.py results/atoms_analyzed.csv results/contacts_analyzed.csv -o ./results -s 1000
```

## Output Figures

### Basic Analysis (4 figures)

| Figure | Content | Key Metrics |
|--------|---------|-------------|
| fig1_comprehensive | Overview: particle distribution, contact types, Z-profile | Type counts, contact breakdown |
| fig2_contact_network | Contact mechanics analysis: (a) Force vs Area, (b) Coord vs Z with trends, (c) Force direction heatmap, (d) Overlap vs Contact Area with Hertzian fit | R², Avg δ, Avg A |
| fig3_ionic_conductivity | SE network analysis for ion transport | SE-SE coordination, percolation |
| fig4_mechanical | Force distribution and contact mechanics | Normal/tangential forces |

#### Figure 2 Panel Details

| Panel | Content | Key Metrics |
|-------|---------|-------------|
| (a) | Normal Force vs Contact Area by contact type (AM-AM, AM-SE, SE-SE) | Force-area relationship |
| (b) | Coordination Number vs Z Position with trend lines | Height uniformity |
| (c) | Force Direction Distribution (2D heatmap: Azimuthal θ vs Polar φ) | Isotropy check |
| (d) | AM-SE Interface: Overlap δ vs Contact Area with linear fit | R², Hertzian validation |

### Advanced Analysis (5 figures)

| Figure | Content | Key Metrics |
|--------|---------|-------------|
| fig5_stress_analysis | Stress tensor analysis | Von Mises, hydrostatic pressure |
| fig6_force_chains | Force chain visualization | Z-alignment, chain fraction |
| fig7_fabric_tensor | Structural anisotropy | Eigenvalues, anisotropy index |
| fig8_porosity_network | Packing, connectivity, NCM vs SE percolation | NCM 99.9% vs SE 0% |
| fig9_torque_rolling | Contact dynamics | Rolling/sliding ratio, friction |

#### Figure 8 Detailed Analysis

Figure 8 is a critical figure showing the contrast between electron (NCM) and ion (SE) transport pathways.

| Panel | Content | Key Finding |
|-------|---------|-------------|
| (a) | Porosity Profile | 0.347 avg, CV 5.5% (uniform) |
| (b) | Packing Density | 0.653 avg (near RCP limit) |
| (c) | Network Statistics | NCM 99.9% vs SE 0% percolation |
| (d) | Coordination Profile | NCM ~28, SE ~3.5 |

**Key Message**: NCM electron network is fully connected (99.9% percolation), but SE ion network requires plastic deformation (0% percolation in rigid state).

### Thick Electrode Analysis (3 figures)

| Figure | Content | Key Metrics |
|--------|---------|-------------|
| thick_electrode_fig1 | Height profile overview | Layer-wise distribution, coordination |
| thick_electrode_fig2 | Top vs Bottom comparison | Stress/coord/force gradients |
| thick_electrode_fig3 | SE network by height | SE-SE coordination, AM-SE interface |

### Deep Analysis (5 figures)

| Figure | Content | Key Metrics |
|--------|---------|-------------|
| deep_fig1_am_utilization | AM surface coverage | Dead AM, SE contacts per AM |
| deep_fig2_contact_quality | Overlap analysis | Hertzian validation, excessive overlap |
| deep_fig3_se_bottlenecks | SE network weak points | Bottleneck fraction, isolated SE |
| deep_fig4_boundary_heterogeneity | Wall effects, center-edge | Radial profiles |
| deep_fig5_force_network | Force transmission | Force by type, friction ratio |

## Key Metrics Reference

### Ionic Conductivity

| Metric | Threshold | Interpretation |
|--------|-----------|----------------|
| SE-SE Coordination | > 2.4 | Percolation threshold for 3D |
| SE Percolation | > 70% | Good ion pathway connectivity |
| NCM-SE Interface | > 80% | Active material utilization |

### Mechanical Properties

| Metric | Good | Warning |
|--------|------|---------|
| Overlap ratio | < 3% | > 5% excessive |
| Hertzian correlation | > 0.9 | < 0.8 model issue |
| Rolling fraction | 20-40% | > 80% not jammed |

### Thick Electrode Specific

| Metric | Typical | Concern |
|--------|---------|---------|
| Top-Bottom packing diff | < 5% | > 10% non-uniform |
| SE bottleneck fraction | < 20% | > 30% ion pathway issue |
| AM surface coverage | 5-15% | < 2% underutilized |

## Scale Factor

- Default: 1000x (simulation uses mm, real is μm)
- AM radius: 4 mm (sim) → 4 μm (real)
- SE radius: 0.8 mm (sim) → 0.8 μm (real)
- Forces: multiply by scale_factor

## File Structure

```
liggghts-dem-analysis/
├── SKILL.md
├── scripts/
│   ├── parse_liggghts.py          # Parse atom/contact dumps
│   ├── analyze_contacts.py         # Basic contact analysis
│   ├── generate_figures.py         # Fig 1-4
│   ├── advanced_analysis.py        # Stress, force chain, fabric
│   ├── generate_advanced_figures.py # Fig 5-9
│   ├── thick_electrode_analysis.py # Height gradient analysis
│   └── deep_electrode_analysis.py  # AM util, bottleneck, quality
└── references/
    ├── parameters.md               # Material properties
    ├── scale_factors.md            # Unit conversion
    ├── advanced_metrics.md         # Theory reference
    └── thick_electrode_guide.md    # Interpretation guide
```

## Example Results: Thick Electrode (214 μm)

### System Overview
- Particles: 217,130 (AM: 3,870, SE: 213,260)
- Contacts: 570,398
- Target pressure: 0.30 MPa

### Contact Distribution
- SE-SE: 60.2% (ionic network)
- AM-SE: 38.8% (interface)
- AM-AM: 1.0% (force chains)

### Key Findings

| Parameter | Bottom | Middle | Top |
|-----------|--------|--------|-----|
| Packing | 70.5% | 70.3% | 62.9% |
| AM Coord | 60.1 | 60.5 | 58.9 |
| SE-SE Coord | 3.22 | 3.20 | 3.30 |

### AM Utilization
- Active AM: 100% (all AM have SE contact)
- Surface coverage: 2.27% (low but typical for size ratio)
- SE contacts per AM: 57.2 average

### SE Network Quality
- Bottleneck fraction: 30.1% (below percolation threshold)
- Isolated SE: 1.1%
- Overall percolation: 98.8%

### Contact Quality
- Mean overlap: 1.32% (good)
- Hertzian correlation: 0.9998 (excellent)
- Excessive overlap: 3.41%

### Force Network
- AM-AM: 1,631 N (force chain backbone)
- AM-SE: 94 N (interface load)
- SE-SE: 9.6 N (network)

## Interpretation Guidelines

### Good Electrode Indicators
- SE-SE coordination > 2.4 throughout height
- Packing uniformity > 90%
- AM surface coverage > 5%
- Bottleneck fraction < 20%
- Hertzian correlation > 0.95

### Warning Signs
- Top layer packing drop > 10%
- SE bottleneck > 30%
- Excessive overlap > 5%
- Isolated SE > 5%
- Non-equilibrated state (high KE)

### Critical Issues
- SE-SE coordination < 2.4 (no percolation)
- Dead AM > 10% (inactive material)
- Overlap ratio > 10% (model breakdown)

---

## DEM Validity Assessment

### Overview

Rigid-sphere DEM has inherent limitations for solid-state battery electrodes where SE plastic deformation is critical. This section provides methods to assess DEM data validity.

### What DEM Can and Cannot Predict

#### Highly Valid (Topology preserved after SE deformation)

| Data | Validity | Reason |
|------|----------|--------|
| Contact network topology | Valid | Who connects to whom is preserved |
| Coordination number distribution | Valid | Contact counts maintained |
| Force chain pathways | Valid | Load transfer routes unchanged |
| Fabric tensor | Valid | Structural anisotropy independent of deformation |
| Layer-wise uniformity (CV) | Valid | Relative comparison valid |

#### Reference Only (Absolute values limited)

| Data | Limitation | Alternative Use |
|------|-----------|----------------|
| Porosity (~35%) | Sphere packing geometric limit | Use CV for uniformity |
| Packing density | Does not reflect SE deformation | Relative layer comparison |

#### Not Predictive

| Data | Reason |
|------|--------|
| Absolute electrode density | SE deformation not modeled |
| Direct ionic conductivity | Depends on actual porosity |

### Porosity Limitation

```
DEM Porosity: ~35% (sphere packing limit)
Real Electrode: <5% (SE plastic deformation)

Why the gap?
- Random Close Packing (RCP) limit: ~36%
- DEM uses rigid spheres that cannot deform
- Real SE fills voids through plastic deformation
- Overlap volume correction provides only ~0.2%p reduction
```

### NCM vs SE Network Analysis

Critical for understanding electron vs ion transport pathways:

```python
# NCM network analysis (electron conduction)
ncm_contacts = df[df['contact_type'] == 'NCM-NCM']
# Build graph, find components, check percolation

# SE network analysis (ion conduction)
se_contacts = df[df['contact_type'] == 'SE-SE']
# Build graph, find components, check percolation
```

#### Typical Results (Figure 8c)

| Network | Components | Largest (%) | Percolation | Meaning |
|---------|-----------|-------------|-------------|---------|
| NCM | 4 | 99.9% | 99.9% | Electron pathway OK |
| SE | 10,736 | 1.4% | 0.0% | Requires deformation |

#### Important Notes:
- Remove floating atoms (particles with no contacts) before analysis
- NCM percolation should be checked within the NCM-existing Z range
- Use contact-based filtering: `df_atom_filtered = df_atom[df_atom['id'].isin(contact_ids)]`

### Coordination Number Profile (Figure 8d)

| Particle | Avg Coordination | Range | Interpretation |
|----------|-----------------|-------|----------------|
| NCM | ~28 | 25-30 | High connectivity, stable network |
| SE | ~3.5 | 3.3-3.7 | Distributed around NCM surfaces |
| SE-SE only | ~2.0 | - | Linear chains, insufficient branching |

### SE Network Fragmentation Cause
- SE-SE avg degree ~2.0 (linear chains, not networks)
- ~35% SE particles isolated or endpoints
- NCM particles (larger) physically block SE-SE connections
- SE distributed around NCM surfaces
- 10,736 fragmented components vs NCM's 4 components

### Porosity Analysis Notes
- Remove floating atoms before calculating porosity
- Calculate porosity within actual particle distribution Z range
- Typical values: ~0.35 (near RCP limit ~0.36)
- CV < 15% indicates uniform compression

---

## Key Messages for Paper

| Figure | Core Message | Validity |
|--------|-------------|----------|
| Fig 1-2 | Uniform particle/composition distribution | High |
| Fig 3 | NCM-SE interface contacts | Very High |
| Fig 4 | NCM bears load, SE assists | High |
| Fig 5 | Stress concentration locations | High |
| Fig 6 | Z-direction force chains | Very High |
| Fig 7 | Isotropic structure (fabric tensor) | Very High |
| Fig 8 | Uniformity (CV 5.5%), connectivity | High |

### Selling Points
- NCM-SE contacts: 80,000+ (sufficient interface)
- NCM coordination: ~28 (stable electron network)
- Overall connectivity: 99.9% (nearly all particles connected)
- Force chains: Z-direction percolation (load transfer)
- Fabric tensor: Isotropic (uniform structure)
- Layer-wise CV: 5.5% (uniform compression)
