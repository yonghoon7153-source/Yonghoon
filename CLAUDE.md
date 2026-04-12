# Argyrodite Mechanical Properties Research — Knowledge Base

## Project Overview
PhD research (BML Lab, Hanyang University, 안용훈) on halogen-substituted argyrodite
solid electrolytes. Goal: predict mechanical properties (B0, E, Cij, Wad) via
multi-scale DFT/MLIP pipeline and understand how Cl/Br substitution affects them.

Paper: "Beyond Electrochemistry: Tailoring Mechanical Properties of Halogen-Substituted Argyrodites"

## Repository Structure
```
db/                          # Structured data (JSON) — machine-queryable
  compositions/              # Per-composition master files (comp1.json ... modelc.json)
  properties/                # Cross-composition property tables (eos.json, elastic.json, ...)
  literature/                # Reference database (refs.json)
  inputs/                    # DFT/MLIP input templates
    qe_templates/            # Quantum ESPRESSO input templates
    mlip_templates/          # MACE/UMA input templates
    adhesion_templates/      # SE/NCM adhesion calculation templates
  pipelines/                 # Pipeline definitions (v1.json, v2.json)
kb/                          # Knowledge base (Markdown) — human-readable
  methodology/               # Pipeline, methods, protocols
  physics/                   # Physical concepts, equations, interpretations
  results/                   # Per-composition and cross-composition analysis
  papers/                    # Paper drafts, notes, imported literature summaries
tools/                       # Python CLI tools
  kb_search.py               # Search across db/ and kb/ with evidence citations
  kb_index.py                # Auto-index new data, rebuild search index
  kb_export.py               # Export tables/figures for paper writing
tests/                       # Tests for tools
```

## Compositions
| ID     | Formula                      | Atoms | Cell    | f.u. | Family |
|--------|------------------------------|-------|---------|------|--------|
| comp1  | Li6PS5Cl                     | 52    | cubic   | 4    | Li6    |
| comp2  | Li6PS5Cl0.5Br0.5             | 52    | cubic   | 4    | Li6    |
| comp3  | Li5.4PS4.4Cl1.0Br0.6         | 62    | rhombo  | 5    | Li5.4  |
| comp4  | Li5.4PS4.4Cl0.8Br0.8         | 62    | rhombo  | 5    | Li5.4  |
| comp5  | Li5.4PS4.4Cl0.6Br1.0         | 62    | rhombo  | 5    | Li5.4  |
| modelc | Li5.4PS4.4Cl1.6              | 62    | rhombo  | 5    | Li5.4  |

## Key Results (Quick Reference)
- B0: comp1(26.2) > comp2(25.8) > comp5(22.9) > modelc(21.7) > comp3(20.8) = comp4(20.8) GPa
- E (600K snap): comp1(29.1) > comp2(28.6) > comp3(27.3) > comp4(26.4) > comp5(25.8) GPa
- Trend: Br increase -> B0/E decrease (intrinsic)
- Li ordering sensitivity: comp5 DC44 = 12.7 GPa (47%), DE = 15.6 GPa (22%)
- Band gap: 2.0-2.3 eV range, Br slightly reduces gap

## Pipeline
Two versions:
- **v1 (current paper):** Halogen enumerate -> MLIP screen -> DFT relax -> DFT EOS -> Post-processing
- **v2 (recommended):** + Li screening (20 random) + 500K MLIP MD annealing -> champion selection

## How to Use This Repo
1. **Query data:** `python tools/kb_search.py "comp5 bulk modulus"` or `python tools/kb_search.py "Br effect on elastic"`
2. **Add new data:** Edit JSON in db/, then `python tools/kb_index.py --rebuild`
3. **Export for paper:** `python tools/kb_export.py --table eos` or `--table elastic`
4. **Browse knowledge:** Read markdown files in kb/

## Conventions
- Energies in eV (or Ry where QE native)
- Pressures/moduli in GPa
- Distances in Angstrom
- Charges in |e|
- Surface/adhesion energies in J/m^2
- comp IDs: comp1, comp2, comp3, comp4, comp5, modelc
- Basin labels: A (most stable), B, C...
- Pipeline versions: v1, v2

## DFT Settings (QE)
- Functional: PBE
- Ecutwfc: 60 Ry, Ecutrho: 480 Ry
- K-grid: 6x6x3 (rhombo), 6x6x6 (cubic)
- Conv_thr: 1e-8 (relax), 1e-10 (tight SCF)
- Smearing: mv, degauss: 0.01 Ry
- PP: SSSP efficiency or PBE ONCV

## MLIP Settings
- Model: MACE-MP-0 (screening), UMA-s-1p2 (adhesion)
- MD: Langevin thermostat, dt=2fs
- Annealing: 500K, 50-100ps
- Elastic: 600K snapshot x5, quench, FIRE relax, finite-strain Cij
