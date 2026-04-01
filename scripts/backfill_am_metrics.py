"""
Backfill AM-AM CN and φ_AM into existing full_metrics.json files.
Also run electronic conductivity regression analysis.
"""
import json, os, csv, sys
import numpy as np
from collections import defaultdict
from scipy import stats
from scipy.optimize import curve_fit

WEBAPP = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'webapp')


def calc_am_am_stats(atoms_path, contacts_path, type_map_str, scale=1000):
    """Calculate AM-AM CN and contact stats from raw CSV."""
    # Parse type_map
    type_map = {}
    if type_map_str:
        for part in type_map_str.split(','):
            if ':' in part:
                k, v = part.split(':', 1)
                type_map[int(k)] = v
    am_types = {t for t, name in type_map.items() if name.startswith('AM')}
    se_types = {t for t, name in type_map.items() if name.startswith('SE')}

    # Load atoms
    atoms = {}
    with open(atoms_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            aid = int(row['id'])
            atoms[aid] = {
                'type': int(row['type']),
                'radius': float(row['radius']),
            }

    # Load contacts — extract ALL fields
    contacts = []
    with open(contacts_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            fn_x = float(row.get('fn_x', 0))
            fn_y = float(row.get('fn_y', 0))
            fn_z = float(row.get('fn_z', 0))
            p1_x = float(row.get('p1_x', 0))
            p1_y = float(row.get('p1_y', 0))
            p1_z = float(row.get('p1_z', 0))
            p2_x = float(row.get('p2_x', 0))
            p2_y = float(row.get('p2_y', 0))
            p2_z = float(row.get('p2_z', 0))
            contacts.append({
                'id1': int(row['id1']),
                'id2': int(row['id2']),
                'contact_area': float(row.get('contact_area', 0)),
                'delta': float(row.get('delta', 0)),
                'fn': np.sqrt(fn_x**2 + fn_y**2 + fn_z**2),
                'hop_dist': np.sqrt((p1_x-p2_x)**2 + (p1_y-p2_y)**2 + (p1_z-p2_z)**2),
            })

    # AM-AM contact stats — extract everything
    cn = defaultdict(int)
    am_am_areas_sim = []
    am_am_deltas = []
    am_am_forces = []
    am_am_hops = []
    am_am_contact_radii = []

    for c in contacts:
        if c['id1'] in atoms and c['id2'] in atoms:
            if atoms[c['id1']]['type'] in am_types and atoms[c['id2']]['type'] in am_types:
                cn[c['id1']] += 1
                cn[c['id2']] += 1
                am_am_areas_sim.append(c['contact_area'])
                am_am_deltas.append(c['delta'])
                am_am_forces.append(c['fn'])
                am_am_hops.append(c['hop_dist'])
                # Contact radius from area: A = π*a² → a = √(A/π)
                if c['contact_area'] > 0:
                    am_am_contact_radii.append(np.sqrt(c['contact_area'] / np.pi))

    am_ids = [aid for aid, a in atoms.items() if a['type'] in am_types]
    se_ids = [aid for aid, a in atoms.items() if a['type'] in se_types]

    if not am_ids:
        return None

    values = np.array([cn.get(aid, 0) for aid in am_ids])

    # Scale: sim → real (μm, μm², μN)
    am_am_areas = [a * scale * scale for a in am_am_areas_sim]  # μm²
    am_am_deltas_real = [d * scale for d in am_am_deltas]  # μm
    am_am_forces_real = [f * scale * scale for f in am_am_forces]  # μN (F_real = F_sim * scale²)
    am_am_hops_real = [h * scale for h in am_am_hops]  # μm
    am_am_radii_real = [r * scale for r in am_am_contact_radii]  # μm

    # Contact pressure: P = F / A (MPa) — need consistent units
    am_am_pressures = []
    for i in range(len(am_am_areas)):
        if am_am_areas[i] > 0:
            # F in μN, A in μm² → P = F/A in μN/μm² = MPa
            am_am_pressures.append(am_am_forces_real[i] / am_am_areas[i])

    # φ_AM calculation (volume-based)
    v_am = sum(4/3 * np.pi * atoms[aid]['radius']**3 for aid in am_ids)
    v_se = sum(4/3 * np.pi * atoms[aid]['radius']**3 for aid in se_ids)
    v_total = v_am + v_se  # solid volume only

    # Bottleneck: minimum contact area per AM particle's contacts
    am_min_area = {}
    am_contact_areas = defaultdict(list)
    for c_idx, c in enumerate(contacts):
        if c['id1'] in atoms and c['id2'] in atoms:
            if atoms[c['id1']]['type'] in am_types and atoms[c['id2']]['type'] in am_types:
                a_real = am_am_areas[len(am_contact_areas[c['id1']]) + len(am_contact_areas[c['id2']])] if len(am_am_areas) > 0 else 0
    # Simpler approach: just compute min/max/percentiles of area distribution
    am_am_areas_arr = np.array(am_am_areas) if am_am_areas else np.array([0])

    return {
        'am_am_cn': float(np.mean(values)),
        'am_am_cn_std': float(np.std(values)),
        'am_am_n_contacts': int(np.sum(values)) // 2,
        # Area (μm²)
        'am_am_mean_area': float(np.mean(am_am_areas)) if am_am_areas else 0,
        'am_am_median_area': float(np.median(am_am_areas)) if am_am_areas else 0,
        'am_am_total_area': float(np.sum(am_am_areas)) if am_am_areas else 0,
        'am_am_min_area': float(np.min(am_am_areas_arr)) if len(am_am_areas) > 0 else 0,
        'am_am_max_area': float(np.max(am_am_areas_arr)) if len(am_am_areas) > 0 else 0,
        # Contact radius (μm)
        'am_am_mean_contact_radius': float(np.mean(am_am_radii_real)) if am_am_radii_real else 0,
        # Overlap / penetration (μm)
        'am_am_mean_delta': float(np.mean(am_am_deltas_real)) if am_am_deltas_real else 0,
        'am_am_max_delta': float(np.max(am_am_deltas_real)) if am_am_deltas_real else 0,
        # Normal force (μN)
        'am_am_mean_force': float(np.mean(am_am_forces_real)) if am_am_forces_real else 0,
        'am_am_max_force': float(np.max(am_am_forces_real)) if am_am_forces_real else 0,
        # Contact pressure (MPa)
        'am_am_mean_pressure': float(np.mean(am_am_pressures)) if am_am_pressures else 0,
        'am_am_max_pressure': float(np.max(am_am_pressures)) if am_am_pressures else 0,
        # Hop distance (μm)
        'am_am_mean_hop': float(np.mean(am_am_hops_real)) if am_am_hops_real else 0,
        # Derived
        'am_am_area_cv': float(np.std(am_am_areas) / np.mean(am_am_areas)) if am_am_areas and np.mean(am_am_areas) > 0 else 0,
        'phi_am_solid': float(v_am / v_total) if v_total > 0 else 0,
    }


def backfill_all():
    """Add AM-AM metrics to all existing full_metrics.json files."""
    count = 0
    for base in [os.path.join(WEBAPP, 'results'), os.path.join(WEBAPP, 'archive')]:
        if not os.path.isdir(base):
            continue
        for root, dirs, files in os.walk(base):
            if 'full_metrics.json' in files and 'atoms.csv' in files and 'contacts.csv' in files:
                met_path = os.path.join(root, 'full_metrics.json')
                atoms_path = os.path.join(root, 'atoms.csv')
                contacts_path = os.path.join(root, 'contacts.csv')

                with open(met_path) as f:
                    met = json.load(f)

                # Determine type_map
                type_map_str = '1:AM_S,2:SE'  # default
                # Try meta.json
                meta_path = os.path.join(root, 'meta.json')
                if not os.path.exists(meta_path):
                    # Dashboard case: meta in uploads
                    case_id = os.path.basename(root)
                    meta_path = os.path.join(WEBAPP, 'uploads', case_id, 'meta.json')
                if os.path.exists(meta_path):
                    with open(meta_path) as f:
                        meta = json.load(f)
                    type_map_str = meta.get('type_map', type_map_str)

                # Determine scale from meta
                scale = 1000
                if os.path.exists(meta_path):
                    with open(meta_path) as f:
                        meta2 = json.load(f)
                    scale = meta2.get('scale', 1000)

                stats = calc_am_am_stats(atoms_path, contacts_path, type_map_str, scale=scale)
                if stats:
                    for k, v in stats.items():
                        met[k] = v
                    # Also ensure phi_am is set
                    if 'phi_se' in met and 'porosity' in met:
                        met['phi_am'] = 1.0 - met['phi_se'] - met['porosity'] / 100.0

                    # AM percolation fraction from am_percolation_sets.json
                    am_perc_path = os.path.join(root, 'am_percolation_sets.json')
                    if os.path.exists(am_perc_path):
                        with open(am_perc_path) as f:
                            am_perc = json.load(f)
                        met['am_percolation_pct'] = am_perc.get('percolation_pct', 0)
                        met['am_n_components'] = am_perc.get('n_components', 0)

                    with open(met_path, 'w') as f:
                        json.dump(met, f, indent=2, default=str)
                    count += 1
                    name = os.path.basename(root)
                    am_perc_pct = met.get('am_percolation_pct', 0)
                    print(f"  {name}: CN={stats['am_am_cn']:.2f}, area={stats['am_am_mean_area']:.3f}μm², "
                          f"a={stats['am_am_mean_contact_radius']:.3f}μm, δ={stats['am_am_mean_delta']:.4f}μm, "
                          f"F={stats['am_am_mean_force']:.2f}μN, hop={stats['am_am_mean_hop']:.2f}μm")

    print(f"\nBackfilled {count} cases")
    return count


def electronic_regression():
    """Find best regression model for electronic conductivity."""
    rows = []

    for base in [os.path.join(WEBAPP, 'results'), os.path.join(WEBAPP, 'archive')]:
        if not os.path.isdir(base):
            continue
        for root, dirs, files in os.walk(base):
            if 'full_metrics.json' not in files:
                continue
            met_path = os.path.join(root, 'full_metrics.json')
            with open(met_path) as f:
                m = json.load(f)

            sigma_el = m.get('electronic_sigma_full_mScm', 0)
            if not sigma_el or sigma_el <= 0:
                continue

            phi_am = m.get('phi_am', 0)
            am_cn = m.get('am_am_cn', 0)
            am_area = m.get('am_am_mean_area', 0)
            phi_se = m.get('phi_se', 0)
            porosity = m.get('porosity', 0)
            se_cn = m.get('se_se_cn', 0)
            thickness = m.get('thickness_um', 0)
            am_perc = m.get('am_percolation_pct', 0)

            if phi_am <= 0 or am_cn <= 0:
                continue

            rows.append({
                'name': os.path.basename(root),
                'sigma_el': sigma_el,
                'phi_am': phi_am,
                'am_cn': am_cn,
                'am_area': am_area,
                'am_contact_radius': m.get('am_am_mean_contact_radius', 0),
                'am_delta': m.get('am_am_mean_delta', 0),
                'am_force': m.get('am_am_mean_force', 0),
                'am_pressure': m.get('am_am_mean_pressure', 0),
                'am_hop': m.get('am_am_mean_hop', 0),
                'am_area_cv': m.get('am_am_area_cv', 0),
                'am_total_area': m.get('am_am_total_area', 0),
                'am_n_contacts': m.get('am_am_n_contacts', 0),
                'phi_se': phi_se,
                'porosity': porosity,
                'se_cn': se_cn,
                'thickness': thickness,
                'am_perc': am_perc / 100.0 if am_perc > 0 else 0,
            })

    if len(rows) < 5:
        print(f"Only {len(rows)} cases with electronic data + AM metrics. Need at least 5.")
        return

    print(f"\n{'='*70}")
    print(f"ELECTRONIC CONDUCTIVITY REGRESSION ({len(rows)} cases)")
    print(f"{'='*70}")

    sigma_el = np.array([r['sigma_el'] for r in rows])
    phi_am = np.array([r['phi_am'] for r in rows])
    am_cn = np.array([r['am_cn'] for r in rows])
    am_area = np.array([r['am_area'] for r in rows])
    thickness = np.array([r['thickness'] for r in rows])

    SIGMA_AM = 50.0  # 0.05 S/cm = 50 mS/cm
    am_perc = np.array([r['am_perc'] for r in rows])

    ss_tot = np.sum((np.log(sigma_el) - np.mean(np.log(sigma_el)))**2)

    # 1. Bruggeman baseline: σ_el = σ_AM × φ_AM^n
    log_ratio = np.log(sigma_el / SIGMA_AM)
    log_phi = np.log(phi_am)
    s, i, r_val, _, _ = stats.linregress(log_phi, log_ratio)
    r2 = r_val**2
    print(f"\n1. Bruggeman: σ_el = σ_AM × φ_AM^{s:.2f}")
    print(f"   R² = {r2:.4f}, n_eff = {s:.2f}")

    # 2. φ_AM^a × CN^b
    X = np.column_stack([np.log(phi_am), np.log(am_cn), np.ones(len(rows))])
    b, _, _, _ = np.linalg.lstsq(X, np.log(sigma_el), rcond=None)
    pred = X @ b
    r2 = 1 - np.sum((np.log(sigma_el) - pred)**2) / ss_tot
    print(f"\n2. σ_el = exp({b[2]:.3f}) × φ_AM^{b[0]:.2f} × CN_AM^{b[1]:.2f}")
    print(f"   R² = {r2:.4f}")

    # 3. φ_AM^a × CN^b × Area^c
    valid_area = am_area > 0
    if np.sum(valid_area) > 5:
        X3 = np.column_stack([np.log(phi_am[valid_area]),
                              np.log(am_cn[valid_area]),
                              np.log(am_area[valid_area]),
                              np.ones(np.sum(valid_area))])
        y3 = np.log(sigma_el[valid_area])
        ss3 = np.sum((y3 - np.mean(y3))**2)
        b3, _, _, _ = np.linalg.lstsq(X3, y3, rcond=None)
        pred3 = X3 @ b3
        r2_3 = 1 - np.sum((y3 - pred3)**2) / ss3
        print(f"\n3. σ_el = exp({b3[3]:.3f}) × φ_AM^{b3[0]:.2f} × CN_AM^{b3[1]:.2f} × A_AM^{b3[2]:.2f}")
        print(f"   R² = {r2_3:.4f}")
    else:
        print("\n3. Skipped (no valid AM-AM area data)")

    # 4. + thickness (check for boundary effect)
    if np.all(thickness > 0):
        X4 = np.column_stack([np.log(phi_am), np.log(am_cn), np.log(thickness), np.ones(len(rows))])
        b4, _, _, _ = np.linalg.lstsq(X4, np.log(sigma_el), rcond=None)
        pred4 = X4 @ b4
        r2_4 = 1 - np.sum((np.log(sigma_el) - pred4)**2) / ss_tot
        print(f"\n4. σ_el = exp({b4[3]:.3f}) × φ_AM^{b4[0]:.2f} × CN_AM^{b4[1]:.2f} × T^{b4[2]:.2f}")
        print(f"   R² = {r2_4:.4f}  ⚠ T dependency = boundary effect?")

    # 5. φ_AM^a × CN^b × f_perc_AM^c (with AM percolation)
    valid_perc = am_perc > 0
    if np.sum(valid_perc) > 5:
        X5 = np.column_stack([np.log(phi_am[valid_perc]),
                              np.log(am_cn[valid_perc]),
                              np.log(am_perc[valid_perc]),
                              np.ones(np.sum(valid_perc))])
        y5 = np.log(sigma_el[valid_perc])
        ss5 = np.sum((y5 - np.mean(y5))**2)
        b5, _, _, _ = np.linalg.lstsq(X5, y5, rcond=None)
        pred5 = X5 @ b5
        r2_5 = 1 - np.sum((y5 - pred5)**2) / ss5
        print(f"\n5. σ_el = exp({b5[3]:.3f}) × φ_AM^{b5[0]:.2f} × CN_AM^{b5[1]:.2f} × f_perc_AM^{b5[2]:.2f}")
        print(f"   R² = {r2_5:.4f}")
    else:
        print(f"\n5. Skipped (only {np.sum(valid_perc)} cases with AM percolation data)")

    # 6. φ_AM^a × CN^b × Area^c × f_perc_AM^d (full model)
    valid_full = valid_area & valid_perc
    if np.sum(valid_full) > 5:
        X6 = np.column_stack([np.log(phi_am[valid_full]),
                              np.log(am_cn[valid_full]),
                              np.log(am_area[valid_full]),
                              np.log(am_perc[valid_full]),
                              np.ones(np.sum(valid_full))])
        y6 = np.log(sigma_el[valid_full])
        ss6 = np.sum((y6 - np.mean(y6))**2)
        b6, _, _, _ = np.linalg.lstsq(X6, y6, rcond=None)
        pred6 = X6 @ b6
        r2_6 = 1 - np.sum((y6 - pred6)**2) / ss6
        print(f"\n6. σ_el = exp({b6[4]:.3f}) × φ_AM^{b6[0]:.2f} × CN_AM^{b6[1]:.2f} × A_AM^{b6[2]:.2f} × f_perc^{b6[3]:.2f}")
        print(f"   R² = {r2_6:.4f}")

    # 7. Ionic analog: C × σ_AM × φ_AM × CN^a × √Area
    if np.sum(valid_area) > 5:
        try:
            def model7(X, C, a):
                phi, cn_val, area_val = X
                return np.log(C * SIGMA_AM * phi * cn_val**a * np.sqrt(area_val))

            X7 = (phi_am[valid_area], am_cn[valid_area], am_area[valid_area])
            y7 = np.log(sigma_el[valid_area])
            popt, _ = curve_fit(model7, X7, y7, p0=[0.01, 2.0], maxfev=10000)
            pred7 = model7(X7, *popt)
            ss7 = np.sum((y7 - np.mean(y7))**2)
            r2_7 = 1 - np.sum((y7 - pred7)**2) / ss7
            print(f"\n7. σ_el = {popt[0]:.4f} × σ_AM × φ_AM × CN_AM^{popt[1]:.2f} × √A_AM")
            print(f"   R² = {r2_7:.4f}  (ionic analog: 2 free params)")
        except Exception as e:
            print(f"\n7. curve_fit failed: {e}")

    # 8. Exclude thin (T < 50μm) cases — check if T dependency disappears
    thick_mask = thickness > 50
    if np.sum(thick_mask) > 10:
        X8 = np.column_stack([np.log(phi_am[thick_mask]),
                              np.log(am_cn[thick_mask]),
                              np.ones(np.sum(thick_mask))])
        y8 = np.log(sigma_el[thick_mask])
        ss8 = np.sum((y8 - np.mean(y8))**2)
        b8, _, _, _ = np.linalg.lstsq(X8, y8, rcond=None)
        pred8 = X8 @ b8
        r2_8 = 1 - np.sum((y8 - pred8)**2) / ss8
        print(f"\n8. [THICK ONLY T>50μm, n={np.sum(thick_mask)}]")
        print(f"   σ_el = exp({b8[2]:.3f}) × φ_AM^{b8[0]:.2f} × CN_AM^{b8[1]:.2f}")
        print(f"   R² = {r2_8:.4f}")

        if np.sum(thick_mask & valid_area) > 5:
            mask_ta = thick_mask & valid_area
            X8b = np.column_stack([np.log(phi_am[mask_ta]),
                                   np.log(am_cn[mask_ta]),
                                   np.log(am_area[mask_ta]),
                                   np.ones(np.sum(mask_ta))])
            y8b = np.log(sigma_el[mask_ta])
            ss8b = np.sum((y8b - np.mean(y8b))**2)
            b8b, _, _, _ = np.linalg.lstsq(X8b, y8b, rcond=None)
            pred8b = X8b @ b8b
            r2_8b = 1 - np.sum((y8b - pred8b)**2) / ss8b
            print(f"   + Area: exp({b8b[3]:.3f}) × φ_AM^{b8b[0]:.2f} × CN_AM^{b8b[1]:.2f} × A_AM^{b8b[2]:.2f}")
            print(f"   R² = {r2_8b:.4f}")

    # Per-case table
    print(f"\n{'='*70}")
    print("PER-CASE TABLE")
    print(f"{'='*70}")
    print(f"{'Name':25s} {'φ_AM':>6s} {'CN':>5s} {'Area':>7s} {'a(μm)':>6s} {'δ(μm)':>7s} {'F(μN)':>7s} {'hop':>6s} {'σ_el':>7s} {'T':>5s}")
    print("-" * 90)
    for r in sorted(rows, key=lambda x: x['sigma_el']):
        print(f"  {r['name'][:23]:23s} {r['phi_am']:6.3f} {r['am_cn']:5.2f} {r['am_area']:7.3f} "
              f"{r['am_contact_radius']:6.3f} {r['am_delta']:7.4f} {r['am_force']:7.2f} "
              f"{r['am_hop']:6.2f} {r['sigma_el']:7.3f} {r['thickness']:5.0f}")

    # Summary
    print(f"\n{'='*70}")
    print("RANKING")
    print(f"{'='*70}")
    print("Model with T: R² high but physically suspect (intensive property).")
    print("Model without T: need Area + AM percolation for good fit.")
    print("Recommendation: Model 6 (full) or Model 3 (area) if percolation data unavailable.")


def electronic_regression_physics():
    """Test physically-motivated fixed-exponent models for electronic conductivity.

    Only uses THICK cases (T > 50 um) and deduplicates by sigma_el value.
    """
    # ---- data loading (same pattern as electronic_regression) ----
    rows = []
    for base in [os.path.join(WEBAPP, 'results'), os.path.join(WEBAPP, 'archive')]:
        if not os.path.isdir(base):
            continue
        for root, dirs, files in os.walk(base):
            if 'full_metrics.json' not in files:
                continue
            met_path = os.path.join(root, 'full_metrics.json')
            with open(met_path) as f:
                m = json.load(f)

            sigma_el = m.get('electronic_sigma_full_mScm', 0)
            if not sigma_el or sigma_el <= 0:
                continue

            phi_am = m.get('phi_am', 0)
            am_cn = m.get('am_am_cn', 0)
            am_area = m.get('am_am_mean_area', 0)
            am_radius = m.get('am_am_mean_contact_radius', 0)
            am_delta = m.get('am_am_mean_delta', 0)
            am_force = m.get('am_am_mean_force', 0)
            am_hop = m.get('am_am_mean_hop', 0)
            am_total_area = m.get('am_am_total_area', 0)
            am_n_contacts = m.get('am_am_n_contacts', 0)
            thickness = m.get('thickness_um', 0)

            if phi_am <= 0 or am_cn <= 0:
                continue

            rows.append({
                'name': os.path.basename(root),
                'sigma_el': sigma_el,
                'phi_am': phi_am,
                'am_cn': am_cn,
                'am_area': am_area,
                'am_radius': am_radius,
                'am_delta': am_delta,
                'am_force': am_force,
                'am_hop': am_hop,
                'am_total_area': am_total_area,
                'am_n_contacts': am_n_contacts,
                'thickness': thickness,
            })

    # ---- filter: thick only (T > 50 um) ----
    rows = [r for r in rows if r['thickness'] > 50]

    # ---- deduplicate by sigma_el (rounded to 3 decimals) ----
    seen = set()
    unique_rows = []
    for r in rows:
        key = round(r['sigma_el'], 3)
        if key not in seen:
            seen.add(key)
            unique_rows.append(r)
    rows = unique_rows

    if len(rows) < 3:
        print(f"Only {len(rows)} thick unique cases. Need at least 3.")
        return

    # ---- extract arrays ----
    sigma_el = np.array([r['sigma_el'] for r in rows])
    phi_am = np.array([r['phi_am'] for r in rows])
    cn = np.array([r['am_cn'] for r in rows])
    am_area = np.array([r['am_area'] for r in rows])
    am_radius = np.array([r['am_radius'] for r in rows])
    am_delta = np.array([r['am_delta'] for r in rows])
    am_hop = np.array([r['am_hop'] for r in rows])
    am_total_area = np.array([r['am_total_area'] for r in rows])
    am_n_contacts = np.array([r['am_n_contacts'] for r in rows])

    SIGMA_AM = 50.0  # 0.05 S/cm = 50 mS/cm

    log_sigma = np.log(sigma_el)
    ss_tot = np.sum((log_sigma - np.mean(log_sigma))**2)

    print(f"\n{'='*70}")
    print(f"PHYSICS-MOTIVATED FIXED-EXPONENT MODELS  (THICK T>50um, n={len(rows)} unique)")
    print(f"{'='*70}")

    # ---- helper: fit single constant C via log-linear, return (C, R2) ----
    def fit_C(log_rhs):
        """Given log(sigma_el) = log(C) + log_rhs, solve for C and R2.

        log_rhs is an array of the same length as sigma_el.
        """
        # log(C) = mean(log_sigma - log_rhs)
        log_C = np.mean(log_sigma - log_rhs)
        pred = log_C + log_rhs
        ss_res = np.sum((log_sigma - pred)**2)
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
        return np.exp(log_C), r2

    # ---- helper: multi-parameter free fit, return (coeffs, R2) ----
    def fit_free(X_cols):
        """log(sigma_el) = b0 + b1*x1 + b2*x2 + ...  via OLS."""
        X = np.column_stack(list(X_cols) + [np.ones(len(rows))])
        b, _, _, _ = np.linalg.lstsq(X, log_sigma, rcond=None)
        pred = X @ b
        ss_res = np.sum((log_sigma - pred)**2)
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
        return b, r2

    # ---- fixed-exponent models (A-J) ----
    results = []

    # Check which rows have valid area > 0 for area-based models
    area_ok = am_area > 0

    # Model A: C * phi^1.5 * CN^2
    log_rhs_A = 1.5 * np.log(phi_am) + 2.0 * np.log(cn)
    C_A, r2_A = fit_C(log_rhs_A)
    results.append(('A', f'C * phi_AM^1.5 * CN^2  (Bruggeman n=1.5)',
                     C_A, r2_A, 1))

    # Model B: C * phi^(5/3) * CN^2
    log_rhs_B = (5.0/3.0) * np.log(phi_am) + 2.0 * np.log(cn)
    C_B, r2_B = fit_C(log_rhs_B)
    results.append(('B', f'C * phi_AM^(5/3) * CN^2  (empirical n=5/3)',
                     C_B, r2_B, 1))

    # Model C: C * sigma_AM * phi^1.5 * CN^2
    log_rhs_C = np.log(SIGMA_AM) + 1.5 * np.log(phi_am) + 2.0 * np.log(cn)
    C_C, r2_C = fit_C(log_rhs_C)
    results.append(('C', f'C * sigma_AM * phi_AM^1.5 * CN^2  (with sigma_AM)',
                     C_C, r2_C, 1))

    # Model D: C * sigma_AM * phi * CN^2
    log_rhs_D = np.log(SIGMA_AM) + 1.0 * np.log(phi_am) + 2.0 * np.log(cn)
    C_D, r2_D = fit_C(log_rhs_D)
    results.append(('D', f'C * sigma_AM * phi_AM * CN^2  (linear phi)',
                     C_D, r2_D, 1))

    # Model E: C * phi^1.5 * CN * sqrt(A_AM)  — only if area available
    if np.all(area_ok):
        log_rhs_E = 1.5 * np.log(phi_am) + 1.0 * np.log(cn) + 0.5 * np.log(am_area)
        C_E, r2_E = fit_C(log_rhs_E)
        results.append(('E', f'C * phi_AM^1.5 * CN * sqrt(A_AM)  (Kirkpatrick)',
                         C_E, r2_E, 1))
    else:
        results.append(('E', f'C * phi_AM^1.5 * CN * sqrt(A_AM)  (Kirkpatrick)',
                         np.nan, np.nan, 1))

    # Model F: C * phi^1.5 * CN^2 * sqrt(A_AM)
    if np.all(area_ok):
        log_rhs_F = 1.5 * np.log(phi_am) + 2.0 * np.log(cn) + 0.5 * np.log(am_area)
        C_F, r2_F = fit_C(log_rhs_F)
        results.append(('F', f'C * phi_AM^1.5 * CN^2 * sqrt(A_AM)  (+ constriction)',
                         C_F, r2_F, 1))
    else:
        results.append(('F', f'C * phi_AM^1.5 * CN^2 * sqrt(A_AM)  (+ constriction)',
                         np.nan, np.nan, 1))

    # Model G: C * (phi - 0.15)^1.5 * CN^2  — percolation threshold 0.15
    phi_shifted_G = phi_am - 0.15
    if np.all(phi_shifted_G > 0):
        log_rhs_G = 1.5 * np.log(phi_shifted_G) + 2.0 * np.log(cn)
        C_G, r2_G = fit_C(log_rhs_G)
        results.append(('G', f'C * (phi_AM-0.15)^1.5 * CN^2  (perc threshold)',
                         C_G, r2_G, 1))
    else:
        results.append(('G', f'C * (phi_AM-0.15)^1.5 * CN^2  (perc threshold)',
                         np.nan, np.nan, 1))

    # Model H: C * (phi - 0.20)^1.5 * CN^2  — percolation threshold 0.20
    phi_shifted_H = phi_am - 0.20
    if np.all(phi_shifted_H > 0):
        log_rhs_H = 1.5 * np.log(phi_shifted_H) + 2.0 * np.log(cn)
        C_H, r2_H = fit_C(log_rhs_H)
        results.append(('H', f'C * (phi_AM-0.20)^1.5 * CN^2  (perc threshold)',
                         C_H, r2_H, 1))
    else:
        results.append(('H', f'C * (phi_AM-0.20)^1.5 * CN^2  (perc threshold)',
                         np.nan, np.nan, 1))

    # Model I: C * phi^2 * CN^2
    log_rhs_I = 2.0 * np.log(phi_am) + 2.0 * np.log(cn)
    C_I, r2_I = fit_C(log_rhs_I)
    results.append(('I', f'C * phi_AM^2 * CN^2  (n=2)',
                     C_I, r2_I, 1))

    # Model J: C * phi^1.5 * CN^1.5
    log_rhs_J = 1.5 * np.log(phi_am) + 1.5 * np.log(cn)
    C_J, r2_J = fit_C(log_rhs_J)
    results.append(('J', f'C * phi_AM^1.5 * CN^1.5  (both 1.5)',
                     C_J, r2_J, 1))

    # ---- Models with contact radius (constriction R = ρ/(2a)) ----
    radius_ok = am_radius > 0

    # Model K1: C * phi^1.5 * CN² * a  (constriction: G ∝ a)
    if np.all(radius_ok):
        log_rhs_K1 = 1.5*np.log(phi_am) + 2*np.log(cn) + np.log(am_radius)
        C_K1, r2_K1 = fit_C(log_rhs_K1)
        results.append(('K1', f'C * phi^1.5 * CN² * a  (G∝a, Holm)',
                         C_K1, r2_K1, 1))

    # Model K2: C * phi^1.5 * CN² * √A  (same as F but explicit)
    # Already Model F

    # Model K3: C * phi^1.5 * CN² * δ  (overlap ∝ contact quality)
    delta_ok = am_delta > 0
    if np.all(delta_ok):
        log_rhs_K3 = 1.5*np.log(phi_am) + 2*np.log(cn) + np.log(am_delta)
        C_K3, r2_K3 = fit_C(log_rhs_K3)
        results.append(('K3', f'C * phi^1.5 * CN² * δ  (overlap)',
                         C_K3, r2_K3, 1))

    # Model K4: C * phi^1.5 * CN² / hop  (shorter hop = better)
    hop_ok = am_hop > 0
    if np.all(hop_ok):
        log_rhs_K4 = 1.5*np.log(phi_am) + 2*np.log(cn) - np.log(am_hop)
        C_K4, r2_K4 = fit_C(log_rhs_K4)
        results.append(('K4', f'C * phi^1.5 * CN² / hop  (shorter=better)',
                         C_K4, r2_K4, 1))

    # Model K5: C * phi^1.5 * CN² * N_contacts  (total edges in network)
    ncon_ok = am_n_contacts > 0
    if np.all(ncon_ok):
        log_rhs_K5 = 1.5*np.log(phi_am) + 2*np.log(cn) + np.log(am_n_contacts)
        C_K5, r2_K5 = fit_C(log_rhs_K5)
        results.append(('K5', f'C * phi^1.5 * CN² * N_contacts',
                         C_K5, r2_K5, 1))

    # Model K6: C * phi^1.5 * CN² * A_total  (total AM-AM interface)
    tarea_ok = am_total_area > 0
    if np.all(tarea_ok):
        log_rhs_K6 = 1.5*np.log(phi_am) + 2*np.log(cn) + np.log(am_total_area)
        C_K6, r2_K6 = fit_C(log_rhs_K6)
        results.append(('K6', f'C * phi^1.5 * CN² * A_total',
                         C_K6, r2_K6, 1))

    # ---- free-fit models ----

    # Model M1: C * phi^a * CN^b  (free a, b)
    b_K, r2_K = fit_free([np.log(phi_am), np.log(cn)])
    results.append(('M1', f'C * phi_AM^{b_K[0]:.2f} * CN^{b_K[1]:.2f}  (free a,b)',
                     np.exp(b_K[2]), r2_K, 3))

    # Model M2: C * phi^a * CN^b * A_AM^c  (free a, b, c)
    if np.all(area_ok):
        b_L, r2_L = fit_free([np.log(phi_am), np.log(cn), np.log(am_area)])
        results.append(('M2', f'C * phi^{b_L[0]:.2f} * CN^{b_L[1]:.2f} * A^{b_L[2]:.2f}  (free)',
                         np.exp(b_L[3]), r2_L, 4))

    # Model M3: C * phi^a * CN^b * a_contact^c  (free, with contact radius)
    if np.all(radius_ok):
        b_M3, r2_M3 = fit_free([np.log(phi_am), np.log(cn), np.log(am_radius)])
        results.append(('M3', f'C * phi^{b_M3[0]:.2f} * CN^{b_M3[1]:.2f} * a^{b_M3[2]:.2f}  (free+radius)',
                         np.exp(b_M3[3]), r2_M3, 4))

    # Model M4: C * phi^a * CN^b * δ^c  (free, with overlap)
    if np.all(delta_ok):
        b_M4, r2_M4 = fit_free([np.log(phi_am), np.log(cn), np.log(am_delta)])
        results.append(('M4', f'C * phi^{b_M4[0]:.2f} * CN^{b_M4[1]:.2f} * δ^{b_M4[2]:.2f}  (free+delta)',
                         np.exp(b_M4[3]), r2_M4, 4))

    # Model M5: C * phi^a * CN^b * hop^c  (free, with hop distance)
    if np.all(hop_ok):
        b_M5, r2_M5 = fit_free([np.log(phi_am), np.log(cn), np.log(am_hop)])
        results.append(('M5', f'C * phi^{b_M5[0]:.2f} * CN^{b_M5[1]:.2f} * hop^{b_M5[2]:.2f}  (free+hop)',
                         np.exp(b_M5[3]), r2_M5, 4))

    # Model M6: ALL — phi^a * CN^b * A^c * a^d * hop^e  (kitchen sink)
    if np.all(area_ok & radius_ok & hop_ok):
        b_M6, r2_M6 = fit_free([np.log(phi_am), np.log(cn), np.log(am_area),
                                 np.log(am_radius), np.log(am_hop)])
        results.append(('M6', f'phi^{b_M6[0]:.2f}*CN^{b_M6[1]:.2f}*A^{b_M6[2]:.2f}*a^{b_M6[3]:.2f}*hop^{b_M6[4]:.2f} (ALL)',
                         np.exp(b_M6[5]), r2_M6, 6))

    # ---- print individual results ----
    for tag, desc, C_val, r2_val, nparams in results:
        if np.isnan(r2_val):
            print(f"\n  Model {tag}: {desc}")
            print(f"    SKIPPED (insufficient data)")
        else:
            print(f"\n  Model {tag}: {desc}")
            print(f"    C = {C_val:.6e},  R^2 = {r2_val:.4f}  ({nparams} free param{'s' if nparams > 1 else ''})")

    # ---- comparison table sorted by R2 ----
    print(f"\n{'='*70}")
    print("COMPARISON TABLE — sorted by R^2 (descending)")
    print(f"{'='*70}")
    print(f"  {'Model':6s} {'#Params':>7s} {'C':>12s} {'R^2':>8s}   Description")
    print(f"  {'-'*6} {'-'*7} {'-'*12} {'-'*8}   {'-'*40}")

    ranked = sorted(results, key=lambda x: x[3] if not np.isnan(x[3]) else -999, reverse=True)
    for tag, desc, C_val, r2_val, nparams in ranked:
        if np.isnan(r2_val):
            print(f"  {tag:6s} {nparams:7d} {'N/A':>12s} {'N/A':>8s}   {desc}")
        else:
            print(f"  {tag:6s} {nparams:7d} {C_val:12.6e} {r2_val:8.4f}   {desc}")

    print(f"\n  Note: Models A-J have 1 free parameter (C only).")
    print(f"  Model K has 3 free params (C, a, b).  Model L has 4 (C, a, b, c).")
    print(f"  Best 'beautiful' formula = highest R^2 among 1-param models.\n")


if __name__ == '__main__':
    print("Step 1: Backfilling AM-AM metrics...")
    backfill_all()
    print("\nStep 2: Electronic conductivity regression...")
    electronic_regression()
    print("\nStep 3: Physics-motivated fixed-exponent models...")
    electronic_regression_physics()
