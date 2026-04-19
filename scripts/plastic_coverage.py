"""
Plastic coverage estimation from DEM contact dumps.

Physics: In DEM with hooke/hysteresis + reduced E_SE (porosity-calibrated),
the resulting particle overlap δ/R* already mimics real plastic-deformed state.
At each AM-SE contact we classify the regime by overlap_ratio = δ/R*:
  δ/R* < 0.003    → elastic (Hertzian point contact, a² = R*·δ)
  0.003-0.01      → elastic-plastic transition
  δ/R* > 0.01     → fully plastic (film formation, larger a)
Integrate film area per AM particle → plastic coverage fraction.

LIGGGHTS contact dump column layout (26 cols, compute cpl with
  pos id force force_normal force_tangential torque contactArea delta contactPoint):
  1-3   : pos1 (x1,y1,z1)
  4-6   : pos2 (x2,y2,z2)
  7-8   : id1, id2
  9     : periodic / ghost flag (ignore)
  10-12 : force (total)
  13-15 : force_normal
  16-18 : force_tangential
  19-21 : torque
  22    : contactArea
  23    : delta            ← key input
  24-26 : contactPoint
"""

from __future__ import annotations
import numpy as np
import os, sys, glob, argparse, json
from collections import defaultdict


# ---------- Physical constants (from literature + lab assumption) ----------
E_REAL_SE     = 24.0e9    # Pa, LPSCl Young's modulus (LAB VALUE: 24 GPa)
E_REAL_AM     = 140.0e9   # Pa, NCM Young's modulus
POISSON_SE    = 0.30
POISSON_AM    = 0.25
SIGMA_Y_SE    = 0.30e9    # Pa, LPSCl yield stress (H/2.8, H≈0.85 GPa)
H_REAL_SE     = 0.85e9    # Pa, LPSCl hardness (Tabor H ≈ 2.8 σ_y)

# Reduced modulus E* for AM-SE contact (dominated by softer SE)
#   1/E* = (1-ν₁²)/E₁ + (1-ν₂²)/E₂
_inv_Estar = (1 - POISSON_AM**2) / E_REAL_AM + (1 - POISSON_SE**2) / E_REAL_SE
E_STAR_AM_SE = 1.0 / _inv_Estar     # ≈ 22.4 GPa

# Plastic regime thresholds on δ/R* derived from Hertzian + Tabor
#   P_max = (2E*/π) · √(δ/R*)         [Hertzian peak contact pressure]
#   Yield onset:      P_max = 1.6 σ_y  → δ/R* = (0.8π σ_y/E*)²
#   Fully plastic:    P_mean = 2.8 σ_y → δ/R* = (2.1π σ_y/E*)²
_ratio = SIGMA_Y_SE / E_STAR_AM_SE   # ≈ 0.0134
DR_YIELD_ONSET    = (0.8 * np.pi * _ratio) ** 2   # ≈ 0.0011  (0.11%)
DR_FULLY_PLASTIC  = (2.1 * np.pi * _ratio) ** 2   # ≈ 0.0078  (0.78%)

# SE = solid electrolyte atom type in the DEM setup
SE_ATOM_TYPE = 3  # thin6/9 = 3 types (1 AM_P, 2 AM_S, 3 SE);  particulate12 = 2 types → override via CLI


# =============================================================
#   Parsers
# =============================================================
def parse_atom_dump(path: str) -> dict[int, dict]:
    """Parse LIGGGHTS atom dump. Returns {atom_id: {type, r, pos(np.ndarray)}}."""
    atoms: dict[int, dict] = {}
    with open(path) as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        if line.startswith("ITEM: ATOMS"):
            cols = line.strip().split()[2:]
            idx = {c: j for j, c in enumerate(cols)}
            for dl in lines[i + 1:]:
                v = dl.strip().split()
                if len(v) < len(cols):
                    break
                aid = int(v[idx["id"]])
                atoms[aid] = {
                    "type": int(v[idx["type"]]),
                    "r":    float(v[idx["radius"]]),
                    "pos":  np.array([float(v[idx["x"]]),
                                      float(v[idx["y"]]),
                                      float(v[idx["z"]])]),
                }
            break
    return atoms


def parse_contact_dump(path: str) -> list[dict]:
    """Parse LIGGGHTS contact dump (pair/gran/local). Returns list of contact dicts."""
    contacts: list[dict] = []
    with open(path) as f:
        lines = f.readlines()
    entries_start = None
    for i, line in enumerate(lines):
        if line.startswith("ITEM: ENTRIES"):
            entries_start = i + 1
            break
    if entries_start is None:
        return contacts
    for dl in lines[entries_start:]:
        v = dl.strip().split()
        if len(v) < 26:
            continue
        try:
            contacts.append({
                "id1":          int(float(v[6])),
                "id2":          int(float(v[7])),
                "force_normal": np.array([float(v[12]), float(v[13]), float(v[14])]),
                "contactArea":  float(v[21]),
                "delta":        float(v[22]),
            })
        except ValueError:
            continue
    return contacts


# =============================================================
#   Plastic coverage computation
# =============================================================
def film_area_from_overlap(delta: float, R_star: float,
                            R_min: float = None,
                            ligg_area: float = None,
                            mode: str = "capped") -> tuple[float, str]:
    """Return (contact film area in m², regime label).
    mode:
      'hertzian'  — pure elastic Hertzian (π R* δ). Underestimates plastic.
      'liggghts'  — use LIGGGHTS-reported contactArea directly. DEM-native.
      'capped'    — physics model with geometric cap a² ≤ R_min² (recommended).
    """
    dr = delta / R_star if R_star > 0 else 0.0
    if dr <= 0:
        return 0.0, "none"

    elastic_area = np.pi * R_star * delta   # Hertzian point contact (small overlap)

    if mode == "hertzian":
        regime = ("elastic" if dr < DR_YIELD_ONSET else
                  "transition" if dr < DR_FULLY_PLASTIC else "plastic")
        return elastic_area, regime

    if mode == "liggghts":
        regime = ("elastic" if dr < DR_YIELD_ONSET else
                  "transition" if dr < DR_FULLY_PLASTIC else "plastic")
        if ligg_area is not None and ligg_area > 0:
            return ligg_area, regime
        return elastic_area, regime

    # Default: 'capped' physics model
    # Geometric ceiling: film radius² can't exceed smaller particle's projected area
    cap_a2 = (R_min * R_min) if R_min else (R_star * R_star)

    if dr < DR_YIELD_ONSET:
        return elastic_area, "elastic"

    if dr < DR_FULLY_PLASTIC:
        # Smooth transition: interpolate elastic → capped plastic
        f = (dr - DR_YIELD_ONSET) / (DR_FULLY_PLASTIC - DR_YIELD_ONSET)
        plastic_a2 = min(R_star * R_star * dr / DR_FULLY_PLASTIC, cap_a2)
        return (1 - f) * elastic_area + f * np.pi * plastic_a2, "transition"

    # Fully plastic: area grows linearly with dr BEYOND threshold, but capped
    scale = dr / DR_FULLY_PLASTIC
    plastic_a2 = min(R_star * R_star * scale, cap_a2)
    return np.pi * plastic_a2, "plastic"


def compute_coverage(atom_path: str, contact_path: str,
                     se_type: int = SE_ATOM_TYPE,
                     mode: str = "capped",
                     dump_contacts: bool = False) -> dict:
    """Compute elastic + plastic coverage per AM particle for a single snapshot.
    If dump_contacts=True, includes per-contact list (for network solver input)."""
    atoms    = parse_atom_dump(atom_path)
    contacts = parse_contact_dump(contact_path)

    if not atoms or not contacts:
        return {"error": "empty atom or contact dump", "atom_path": atom_path}

    am_surface: dict[int, float] = {}
    for aid, a in atoms.items():
        if a["type"] != se_type:
            am_surface[aid] = 4.0 * np.pi * a["r"] ** 2  # full sphere surface

    elastic_sum  = defaultdict(float)
    plastic_sum  = defaultdict(float)
    regime_count = defaultdict(int)
    delta_stats  = []
    per_contact: list[dict] = []   # only populated if dump_contacts

    for c in contacts:
        a1, a2 = atoms.get(c["id1"]), atoms.get(c["id2"])
        if a1 is None or a2 is None:
            continue
        t1, t2 = a1["type"], a2["type"]
        is_se1, is_se2 = (t1 == se_type), (t2 == se_type)
        if is_se1 == is_se2:
            continue

        am_atom = a2 if is_se1 else a1
        am_id   = c["id2"] if is_se1 else c["id1"]
        se_atom = a1 if is_se1 else a2
        se_id   = c["id1"] if is_se1 else c["id2"]

        R_star = (am_atom["r"] * se_atom["r"]) / (am_atom["r"] + se_atom["r"])
        R_min  = min(am_atom["r"], se_atom["r"])
        delta  = c["delta"]
        if delta <= 0 or R_star <= 0:
            continue

        dr = delta / R_star
        delta_stats.append(dr)

        elastic_area = np.pi * R_star * delta  # Hertzian baseline
        plastic_area, regime = film_area_from_overlap(
            delta, R_star, R_min=R_min,
            ligg_area=c.get("contactArea"), mode=mode)

        elastic_sum[am_id] += elastic_area
        plastic_sum[am_id] += plastic_area
        regime_count[regime] += 1

        if dump_contacts:
            per_contact.append({
                "am_id": am_id, "se_id": se_id,
                "R_am": am_atom["r"], "R_se": se_atom["r"],
                "R_star": R_star, "R_min": R_min,
                "delta": delta, "delta_over_R": dr,
                "regime": regime,
                "elastic_area": elastic_area,
                "ligg_area":    c.get("contactArea", 0.0),
                "plastic_area": plastic_area,
                "amp":          plastic_area / elastic_area if elastic_area > 0 else 0.0,
            })

    elastic_cov = []
    plastic_cov = []
    for aid, surf in am_surface.items():
        if aid in elastic_sum:
            elastic_cov.append(min(elastic_sum[aid] / surf, 1.0))
            plastic_cov.append(min(plastic_sum[aid] / surf, 1.0))

    out = {
        "mode":               mode,
        "n_am":               len(am_surface),
        "n_am_with_contact":  len(elastic_cov),
        "n_contacts_am_se":   sum(regime_count.values()),
        "regime_counts":      dict(regime_count),
        "delta_over_R_mean":  float(np.mean(delta_stats))   if delta_stats else 0.0,
        "delta_over_R_med":   float(np.median(delta_stats)) if delta_stats else 0.0,
        "delta_over_R_max":   float(np.max(delta_stats))    if delta_stats else 0.0,
        "elastic_cov_mean":   float(np.mean(elastic_cov))   if elastic_cov else 0.0,
        "elastic_cov_med":    float(np.median(elastic_cov)) if elastic_cov else 0.0,
        "plastic_cov_mean":   float(np.mean(plastic_cov))   if plastic_cov else 0.0,
        "plastic_cov_med":    float(np.median(plastic_cov)) if plastic_cov else 0.0,
        "cov_amplification":  (float(np.mean(plastic_cov) / np.mean(elastic_cov))
                               if elastic_cov and np.mean(elastic_cov) > 0 else 0.0),
    }
    if dump_contacts:
        out["contacts"] = per_contact
    return out


# =============================================================
#   CLI (step-by-step verification)
# =============================================================
def _pick_latest(dirpath: str, pattern: str) -> str | None:
    files = sorted(glob.glob(os.path.join(dirpath, pattern)),
                   key=lambda p: int(''.join(c for c in os.path.basename(p)
                                              if c.isdigit()) or '0'))
    return files[-1] if files else None


def main_inspect(contact_path: str, n_show: int = 5) -> None:
    """Step 1: verify contact dump column mapping. Print first n contacts + stats."""
    contacts = parse_contact_dump(contact_path)
    print(f"=== Contact dump inspection: {contact_path} ===")
    print(f"Total contacts parsed: {len(contacts)}")
    if not contacts:
        print("  (nothing parsed — check column layout!)")
        return
    print(f"\nFirst {n_show} entries:")
    for i, c in enumerate(contacts[:n_show]):
        F_n = float(np.linalg.norm(c["force_normal"]))
        print(f"  [{i}] id1={c['id1']:5d}  id2={c['id2']:5d}  "
              f"|F_n|={F_n:.4e}  area={c['contactArea']:.4e}  delta={c['delta']:.4e}")
    deltas = np.array([c["delta"] for c in contacts if c["delta"] > 0])
    print(f"\nDelta statistics (N = {len(deltas)}):")
    print(f"  min = {deltas.min():.3e}")
    print(f"  med = {np.median(deltas):.3e}")
    print(f"  max = {deltas.max():.3e}")


def main_case(case_dir: str, se_type: int = SE_ATOM_TYPE) -> None:
    """Step 2: compute plastic coverage for one case (latest snapshot)."""
    atom_f    = _pick_latest(case_dir, "atom_*.liggghts")
    contact_f = _pick_latest(case_dir, "contact_*.liggghts")
    if not atom_f or not contact_f:
        print(f"!! No atom/contact dump found in {case_dir}")
        return
    print(f"=== Plastic coverage for {case_dir} ===")
    print(f"  atom    : {os.path.basename(atom_f)}")
    print(f"  contact : {os.path.basename(contact_f)}")
    print(f"  se_type : {se_type}")
    res = compute_coverage(atom_f, contact_f, se_type=se_type)
    print(json.dumps(res, indent=2))


def main_batch(root: str, pattern: str = "post_*",
               se_type: int = SE_ATOM_TYPE,
               csv_out: str = "plastic_coverage.csv") -> None:
    """Step 3: batch over all cases matching pattern under root, write CSV."""
    dirs = sorted(glob.glob(os.path.join(root, pattern)))
    print(f"=== Batch plastic coverage: {len(dirs)} directories ===")
    results = []
    for d in dirs:
        if not os.path.isdir(d):
            continue
        atom_f    = _pick_latest(d, "atom_*.liggghts")
        contact_f = _pick_latest(d, "contact_*.liggghts")
        if not atom_f or not contact_f:
            print(f"  SKIP {os.path.basename(d)} — missing dump")
            continue
        try:
            res = compute_coverage(atom_f, contact_f, se_type=se_type)
        except Exception as e:
            print(f"  FAIL {os.path.basename(d)}  {e}")
            continue
        res["case"] = os.path.basename(d)
        results.append(res)
        print(f"  OK  {res['case']:30s}  "
              f"elastic={res['elastic_cov_mean']:.3f}  "
              f"plastic={res['plastic_cov_mean']:.3f}  "
              f"amp={res['cov_amplification']:.2f}x")

    # Write CSV
    if results:
        import csv
        keys = ["case", "n_am", "n_am_with_contact", "n_contacts_am_se",
                "delta_over_R_mean", "delta_over_R_med", "delta_over_R_max",
                "elastic_cov_mean", "elastic_cov_med",
                "plastic_cov_mean", "plastic_cov_med",
                "cov_amplification"]
        with open(csv_out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=keys)
            w.writeheader()
            for r in results:
                w.writerow({k: r.get(k, "") for k in keys})
        print(f"\nWrote {csv_out} ({len(results)} rows)")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="DEM plastic coverage estimator")
    sub = ap.add_subparsers(dest="cmd")

    p1 = sub.add_parser("inspect", help="verify contact dump parsing")
    p1.add_argument("contact_file")
    p1.add_argument("--n-show", type=int, default=5)

    p2 = sub.add_parser("case", help="coverage for single case directory")
    p2.add_argument("case_dir")
    p2.add_argument("--se-type", type=int, default=SE_ATOM_TYPE)

    p3 = sub.add_parser("batch", help="coverage for all cases under a root")
    p3.add_argument("root")
    p3.add_argument("--pattern", default="post_*")
    p3.add_argument("--se-type", type=int, default=SE_ATOM_TYPE)
    p3.add_argument("--csv-out", default="plastic_coverage.csv")

    args = ap.parse_args()
    if args.cmd == "inspect":
        main_inspect(args.contact_file, args.n_show)
    elif args.cmd == "case":
        main_case(args.case_dir, args.se_type)
    elif args.cmd == "batch":
        main_batch(args.root, args.pattern, args.se_type, args.csv_out)
    else:
        ap.print_help()
