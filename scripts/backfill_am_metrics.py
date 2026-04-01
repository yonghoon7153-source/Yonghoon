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


def calc_am_am_stats(atoms_path, contacts_path, type_map_str):
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

    # Load contacts
    contacts = []
    with open(contacts_path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            contacts.append({
                'id1': int(row['id1']),
                'id2': int(row['id2']),
                'contact_area': float(row.get('contact_area', 0)),
            })

    # AM-AM CN
    cn = defaultdict(int)
    am_am_areas = []
    for c in contacts:
        if c['id1'] in atoms and c['id2'] in atoms:
            if atoms[c['id1']]['type'] in am_types and atoms[c['id2']]['type'] in am_types:
                cn[c['id1']] += 1
                cn[c['id2']] += 1
                am_am_areas.append(c['contact_area'])

    am_ids = [aid for aid, a in atoms.items() if a['type'] in am_types]
    se_ids = [aid for aid, a in atoms.items() if a['type'] in se_types]

    if not am_ids:
        return None

    values = np.array([cn.get(aid, 0) for aid in am_ids])

    # φ_AM calculation (volume-based)
    v_am = sum(4/3 * np.pi * atoms[aid]['radius']**3 for aid in am_ids)
    v_se = sum(4/3 * np.pi * atoms[aid]['radius']**3 for aid in se_ids)
    v_total = v_am + v_se  # solid volume only

    return {
        'am_am_cn': float(np.mean(values)),
        'am_am_cn_std': float(np.std(values)),
        'am_am_n_contacts': int(np.sum(values)) // 2,
        'am_am_mean_area': float(np.mean(am_am_areas)) if am_am_areas else 0,
        'am_am_total_area': float(np.sum(am_am_areas)) if am_am_areas else 0,
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

                stats = calc_am_am_stats(atoms_path, contacts_path, type_map_str)
                if stats:
                    for k, v in stats.items():
                        met[k] = v
                    # Also ensure phi_am is set
                    if 'phi_se' in met and 'porosity' in met:
                        met['phi_am'] = 1.0 - met['phi_se'] - met['porosity'] / 100.0

                    with open(met_path, 'w') as f:
                        json.dump(met, f, indent=2, default=str)
                    count += 1
                    name = os.path.basename(root)
                    print(f"  {name}: AM-AM CN={stats['am_am_cn']:.2f}, φ_AM_solid={stats['phi_am_solid']:.3f}, "
                          f"contacts={stats['am_am_n_contacts']}, area={stats['am_am_mean_area']:.2f}")

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

            if phi_am <= 0 or am_cn <= 0:
                continue

            rows.append({
                'name': os.path.basename(root),
                'sigma_el': sigma_el,
                'phi_am': phi_am,
                'am_cn': am_cn,
                'am_area': am_area,
                'phi_se': phi_se,
                'porosity': porosity,
                'se_cn': se_cn,
                'thickness': thickness,
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

    # Bruggeman baseline: σ_el_brug = σ_AM × φ_AM^n
    log_ratio = np.log(sigma_el / SIGMA_AM)
    log_phi = np.log(phi_am)

    s, i, r, _, _ = stats.linregress(log_phi, log_ratio)
    r2 = r**2
    print(f"\n1. Bruggeman: σ_el = σ_AM × φ_AM^{s:.2f}")
    print(f"   R² = {r2:.4f}, n_eff = {s:.2f}")

    # 2. φ_AM^a × CN^b
    ss_tot = np.sum((np.log(sigma_el) - np.mean(np.log(sigma_el)))**2)

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

    # 4. Add thickness
    if np.all(thickness > 0):
        X4 = np.column_stack([np.log(phi_am), np.log(am_cn), np.log(thickness), np.ones(len(rows))])
        b4, _, _, _ = np.linalg.lstsq(X4, np.log(sigma_el), rcond=None)
        pred4 = X4 @ b4
        r2_4 = 1 - np.sum((np.log(sigma_el) - pred4)**2) / ss_tot
        print(f"\n4. σ_el = exp({b4[3]:.3f}) × φ_AM^{b4[0]:.2f} × CN_AM^{b4[1]:.2f} × T^{b4[2]:.2f}")
        print(f"   R² = {r2_4:.4f}")

    # 5. Analog to ionic multi-scale: σ_AM × φ_AM × CN_AM^a × Area^b
    try:
        def model5(X, C, a, b_exp):
            phi, cn_val, area_val = X
            return np.log(C * SIGMA_AM * phi * cn_val**a * area_val**b_exp)

        valid5 = am_area > 0
        X5 = (phi_am[valid5], am_cn[valid5], am_area[valid5])
        y5 = np.log(sigma_el[valid5])
        popt, _ = curve_fit(model5, X5, y5, p0=[0.01, 2.0, 0.5], maxfev=10000)
        pred5 = model5(X5, *popt)
        ss5 = np.sum((y5 - np.mean(y5))**2)
        r2_5 = 1 - np.sum((y5 - pred5)**2) / ss5
        print(f"\n5. σ_el = {popt[0]:.4f} × σ_AM × φ_AM × CN_AM^{popt[1]:.2f} × A_AM^{popt[2]:.2f}")
        print(f"   R² = {r2_5:.4f}")
    except Exception as e:
        print(f"\n5. curve_fit failed: {e}")

    # Per-case table
    print(f"\n{'='*70}")
    print("PER-CASE TABLE")
    print(f"{'='*70}")
    print(f"{'Name':25s} {'φ_AM':>6s} {'CN_AM':>6s} {'Area':>8s} {'σ_el':>8s} {'T(μm)':>6s}")
    print("-" * 65)
    for r in sorted(rows, key=lambda x: x['sigma_el']):
        print(f"  {r['name'][:23]:23s} {r['phi_am']:6.3f} {r['am_cn']:6.2f} {r['am_area']:8.2f} {r['sigma_el']:8.3f} {r['thickness']:6.0f}")

    # Save best model info
    print(f"\n{'='*70}")
    print("BEST MODEL SUMMARY")
    print(f"{'='*70}")
    print("Use the model with highest R² and physical interpretability.")
    print("Ionic analog: σ_el = σ_AM × φ_AM × (contact quality correction)")


if __name__ == '__main__':
    print("Step 1: Backfilling AM-AM metrics...")
    backfill_all()
    print("\nStep 2: Electronic conductivity regression...")
    electronic_regression()
