# VESTA Settings for Adhesion Interface Figure

## 1. File Load
- Open comp1_v5xy_s52.xyz (or other comp)
- VESTA will read as P1 symmetry (OK for interface structures)

## 2. Z-Crop (Boundary Settings)
**Edit → Boundary...**

| Axis | Min | Max |
|------|-----|-----|
| x (a) | 0.000 | 1.000 |
| y (b) | 0.000 | 1.000 |
| z (c) | 0.150 | 0.350 |

→ z = 0.15×80 = 12A to 0.35×80 = 28A → shows ~16A around interface
→ Adjust z range per composition:
  - comp1/2B: z(c) = 0.15 ~ 0.35 (NCM top + SE layer 1)
  - comp3/4/5: different cell_z, adjust accordingly

## 3. Atom Sizes (Properties → Atoms...)
**Objects → Properties → Atoms**

| Element | Radius (A) | Color | Role |
|---------|-----------|-------|------|
| Li | 0.40 | Light gray #C0C0C0 | Small, subtle |
| Ni | 0.55 | Silver #A0A0A0 | NCM metal |
| O | 0.45 | Red #FF3333 | NCM anion |
| P | 0.50 | Purple #9933CC | SE center |
| S | 0.55 | Yellow #FFCC00 | SE anion |
| Cl | 0.55 | Green #33CC33 | SE halogen |
| Br | 0.60 | Dark red #CC3300 | SE halogen (comp2-5) |

Set Style: Ball (not Ball and Stick for cleaner look)

## 4. Bond Settings (Properties → Bonds...)
**Objects → Properties → Bonds**

Turn OFF all bonds for cleaner interface view:
- Select all bond pairs → Delete
- Or: set all max distances to 0

If you want PS4 tetrahedra only:
| Pair | Min (A) | Max (A) | Style |
|------|---------|---------|-------|
| P-S | 0.0 | 2.3 | Stick, width 0.10 |
| All others | — | — | Delete |

## 5. Polyhedra (Optional)
**Objects → Properties → Polyhedra**

For PS4 tetrahedra visualization:
| Center | Vertex | Max dist |
|--------|--------|----------|
| P | S | 2.3 |

Color: Light purple, opacity 0.3

## 6. Background & Style
**Edit → Preference → Background**
- Background: White (#FFFFFF)
- Or: Light gray gradient for publication

**Style → Atoms**
- Ambient: 0.20
- Diffuse: 0.60
- Specular: 0.15
- Shininess: 10

## 7. Orientation
**Best view for interface:**
- Rotate to look along [110] or [100] direction
- Interface plane (xy) should be horizontal
- z-axis vertical → NCM at bottom, SE at top

**Camera:**
- Orthographic projection (not perspective)
- Edit → Preference → Display → Projection: Orthographic

## 8. Labels & Annotations
- Turn OFF atom labels (too cluttered with 820 atoms)
- Add manual text annotations in post-processing (PowerPoint/Illustrator):
  - "NCM" arrow pointing to bottom layer
  - "SE" arrow pointing to top layer
  - "Interface" bracket at z=17-19A region
  - Dashed line at NCM-SE boundary

## 9. Export
**File → Export Raster Image...**
- Format: TIFF or PNG
- Resolution: 300 dpi (for journal) or 600 dpi
- Size: 8 cm width (single column) or 17 cm (double column)
- Transparent background: optional

## 10. Multi-Panel Figure
For 5 compositions side-by-side:
1. Export each comp at SAME z-crop, SAME orientation, SAME atom sizes
2. Combine in PowerPoint/Illustrator
3. Label: comp1, comp2B, comp3, comp4, comp5
4. Add Wad values below each panel
5. Divider line between Li6 and Li5.4 groups

## Quick Checklist
- [ ] Z-crop: only interface region (~16A)
- [ ] Bonds: OFF (or PS4 only)
- [ ] Atom sizes: consistent across all 5 comps
- [ ] Background: white
- [ ] Projection: orthographic
- [ ] Same orientation for all 5
- [ ] Export: 300+ dpi TIFF
- [ ] Labels added in post-processing
