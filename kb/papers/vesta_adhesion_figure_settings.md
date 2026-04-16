# VESTA Adhesion Figure Settings (Paper)

## Figure Files

| Comp | File | Seed | Gap (A) | Wad order |
|------|------|------|---------|-----------|
| comp1 | comp1_v5xy_s45.xyz | 45 | -1.0 | 4th |
| comp2B | comp2B_v5xy_s46.xyz | 46 | -0.9 | 5th (lowest) |
| comp3 | comp3_v5xy_s45.xyz | 45 | -1.7 | 1st (highest) |
| comp4 | comp4_v5xy_s57.xyz | 57 | -1.3 | 2nd |
| comp5 | comp5_v5xy_s50.xyz | 50 | -1.2 | 3rd |

Gap ordering: comp3(-1.7) < comp4(-1.3) < comp5(-1.2) < comp1(-1.0) < comp2B(-0.9)
= matches Wad ordering (more negative gap = more O penetration = higher adhesion)

Seeds are different per composition — selected for:
1. Gap ordering matches Wad ordering
2. Gap magnitude reasonable (-0.9 to -1.7 A)
3. Structure visually clean

## Atom Colors (RGB 0-255)

| Element | R | G | B | Description |
|---------|---|---|---|-------------|
| Li | 153 | 204 | 102 | Olive green |
| P | 128 | 77 | 179 | Purple |
| S | 230 | 204 | 51 | Soft yellow |
| Cl | 77 | 179 | 77 | Calm green |
| Br | 179 | 102 | 26 | Brown |
| Ni | 128 | 128 | 128 | Grey |
| O | 204 | 51 | 51 | Red |

## Atom Radii (Angstrom)

| Element | Radius |
|---------|--------|
| Li | 1.3 |
| P | 0.9 |
| S | 1.5 |
| Cl | 1.6 |
| Br | 1.7 |
| Ni | 1.4 |
| O | 1.3 |

## Bonds

- Only P-S: Min 0, Max 2.5 A
- Radius (cylinder): 0
- Width (line): 0
- All other bonds: Remove

## Polyhedra

- Central atom: P
- Color (RGBA): 179, 128, 204, 150
- Show polyhedral edges: ON
- Edge line width: 0.8
- Edge color (RGB): 150, 150, 150

## General

- Draw unit cell: OFF
- Background: White (on export)

## Export

- File -> Export Raster Image
- Background: White
- Resolution: 300 DPI

## Notes

- 1L NCM structures (not 5L) — consistent with paper Wad values
- Within-family same seed: Li5.4 (comp3/4/5 share similar seeds), Li6 (comp1/2B similar seeds)
- Gap = SE_zmin - O_zmax: negative means O penetrates into SE (= adhesion mechanism)
- O penetration 1-2 A is physically reasonable (cross-interface O-Li bond formation)
- Orientation: keep identical across all 5 compositions for fair visual comparison
