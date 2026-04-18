"""
R_true fitting analysis using network solver results.
Run locally after run_all_network.py completes.
"""
import json, os, sys, numpy as np
from scipy import stats
from scipy.optimize import curve_fit

RESULTS_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'webapp', 'results')

def main():
    # Load network results
    all_path = os.path.join(RESULTS_DIR, 'network_conductivity_all.json')
    if not os.path.exists(all_path):
        print("ERROR: Run run_all_network.py first!")
        sys.exit(1)

    with open(all_path) as f:
        net_data = json.load(f)

    # Deduplicate
    seen = {}
    for d in net_data:
        name = d.get('name', '')
        if name not in seen or d.get('sigma_full') is not None:
            seen[name] = d
    unique = list(seen.values())
    print(f"Unique cases: {len(unique)}")

    # Match with full_metrics.json
    results = []
    for nd in unique:
        case_id = nd.get('case_id', '')
        sigma_full = nd.get('sigma_full')
        sigma_bulk_net = nd.get('sigma_bulk_net')
        R_brug = nd.get('R_brug_over_full')
        bulk_frac = nd.get('bulk_resistance_fraction')

        if sigma_full is None or R_brug is None:
            continue

        if case_id.startswith('archive:'):
            base = os.path.join(os.path.dirname(RESULTS_DIR), 'archive')
            metrics_path = os.path.join(base, case_id[8:], 'full_metrics.json')
        else:
            metrics_path = os.path.join(RESULTS_DIR, case_id, 'full_metrics.json')

        if not os.path.exists(metrics_path):
            continue

        with open(metrics_path) as f:
            m = json.load(f)

        gb_d = m.get('gb_density_mean', 0)
        T = m.get('thickness_um', 0)
        tau = m.get('tortuosity_recommended', m.get('tortuosity_mean', 0))
        phi_se = m.get('phi_se', 0)
        f_perc = m.get('percolation_pct', 0) / 100
        cn = m.get('se_se_cn', 0)
        hop_area = m.get('path_hop_area_mean', 0)
        bottleneck = m.get('path_hop_area_min_mean', 0)
        g_path = m.get('path_conductance_mean', 0)

        if gb_d <= 0 or T <= 0 or tau <= 0:
            continue

        results.append({
            'name': nd.get('name', ''),
            'sigma_full': sigma_full, 'sigma_full_mScm': nd.get('sigma_full_mScm', 0),
            'sigma_bulk_net': sigma_bulk_net,
            'R_brug': R_brug, 'bulk_frac': bulk_frac,
            'gb_d': gb_d, 'T': T, 'tau': tau,
            'phi_se': phi_se, 'f_perc': f_perc,
            'cn': cn, 'hop_area': hop_area,
            'bottleneck': bottleneck, 'g_path': g_path,
        })

    n = len(results)
    print(f"Matched cases with metrics: {n}")

    R_true = np.array([r['R_brug'] for r in results])
    gb_d = np.array([r['gb_d'] for r in results])
    T = np.array([r['T'] for r in results])
    tau = np.array([r['tau'] for r in results])
    phi_se = np.array([r['phi_se'] for r in results])
    cn = np.array([r['cn'] for r in results])
    hop_area = np.array([r['hop_area'] for r in results])
    sigma_full_mScm = np.array([r['sigma_full_mScm'] for r in results])

    logR = np.log(R_true)
    ss_tot = np.sum((logR - np.mean(logR))**2)
    ss_tot_R = np.sum((R_true - np.mean(R_true))**2)

    print(f"\nR_true range: {R_true.min():.2f} ~ {R_true.max():.2f}")
    print(f"σ_full range: {sigma_full_mScm.min():.4f} ~ {sigma_full_mScm.max():.4f} mS/cm")

    print("\n" + "="*70)
    print(f"MODEL FITTING (n={n})")
    print("="*70)

    model_results = []

    # A: Power Law (GB_d²×T)
    logx = np.log(gb_d**2 * T)
    s, i, r, _, _ = stats.linregress(logx, logR)
    r2 = r**2
    print(f"\nA: log(R) = {s:.4f}·log(GB_d²×T) + {i:.4f}  |  R²={r2:.4f}")
    model_results.append(('A: BLM (GB_d²×T)^α', r2, 2, f'α={s:.4f}'))

    # B: Additive R = 1 + C×(GB_d²×T)^α
    try:
        def mB(x, C, a): return 1 + C * x**a
        p, _ = curve_fit(mB, gb_d**2*T, R_true, p0=[0.1, 0.3], maxfev=10000)
        pred = mB(gb_d**2*T, *p)
        r2 = 1 - np.sum((R_true-pred)**2)/ss_tot_R
        print(f"B: R = 1 + {p[0]:.4f}×(GB_d²×T)^{p[1]:.4f}  |  R²={r2:.4f}")
        model_results.append(('B: 1+C×(GB_d²×T)^α', r2, 2, f'C={p[0]:.4f}, α={p[1]:.4f}'))
    except Exception as e:
        print(f"B: FAILED ({e})")

    # C: Multi-variable
    valid = (cn > 0) & (hop_area > 0)
    if valid.sum() >= 5:
        X = np.column_stack([np.log(gb_d[valid]), np.log(T[valid]),
                             np.log(cn[valid]), np.log(hop_area[valid]), np.ones(valid.sum())])
        b, _, _, _ = np.linalg.lstsq(X, logR[valid], rcond=None)
        pred = X @ b
        r2 = 1 - np.sum((logR[valid]-pred)**2)/np.sum((logR[valid]-np.mean(logR[valid]))**2)
        print(f"C: GB_d:{b[0]:.3f} T:{b[1]:.3f} CN:{b[2]:.3f} hop:{b[3]:.3f} c:{b[4]:.3f}  |  R²={r2:.4f}")
        model_results.append(('C: Multi-variable (5p)', r2, 5, f'GB_d:{b[0]:.2f},T:{b[1]:.2f},CN:{b[2]:.2f},hop:{b[3]:.2f}'))

    # D: GB_d only (no T)
    s, i, r, _, _ = stats.linregress(np.log(gb_d), logR)
    r2 = r**2
    print(f"D: log(R) = {s:.4f}·log(GB_d) + {i:.4f}  |  R²={r2:.4f}")
    model_results.append(('D: GB_d only', r2, 2, f'a={s:.4f}'))

    # E: R = 1 + β×GB_d^a (additive, no T)
    try:
        def mE(x, beta, a): return 1 + beta * x**a
        p, _ = curve_fit(mE, gb_d, R_true, p0=[1.0, 2.0], maxfev=10000)
        pred = mE(gb_d, *p)
        r2 = 1 - np.sum((R_true-pred)**2)/ss_tot_R
        print(f"E: R = 1 + {p[0]:.4f}×GB_d^{p[1]:.4f}  |  R²={r2:.4f}")
        model_results.append(('E: 1+β×GB_d^a (no T)', r2, 2, f'β={p[0]:.4f}, a={p[1]:.4f}'))
    except Exception as e:
        print(f"E: FAILED ({e})")

    # F: Bruggeman exponent
    sigma_ratio = sigma_full_mScm / 3.0  # σ_grain
    valid_f = (phi_se > 0) & (sigma_ratio > 0)
    s, i, r, _, _ = stats.linregress(np.log(phi_se[valid_f]), np.log(sigma_ratio[valid_f]))
    r2 = r**2
    print(f"F: σ/σ_bulk = exp({i:.3f})×φ_SE^{s:.3f}  |  R²={r2:.4f}  (n_eff≈{s:.1f}, 문헌≈3)")
    model_results.append(('F: Bruggeman φ_SE^n', r2, 2, f'n_eff={s:.2f}'))

    # G: GB_d, T 독립
    X = np.column_stack([np.log(gb_d), np.log(T), np.ones(n)])
    b, _, _, _ = np.linalg.lstsq(X, logR, rcond=None)
    pred = X @ b
    r2 = 1 - np.sum((logR-pred)**2)/ss_tot
    print(f"G: GB_d^{b[0]:.3f}×T^{b[1]:.3f}×exp({b[2]:.3f})  |  R²={r2:.4f}")
    model_results.append(('G: GB_d^a×T^b (M6)', r2, 3, f'a={b[0]:.3f},b={b[1]:.3f}'))

    # H: τ_eff² 역산
    tau_eff2 = 3.0 * phi_se / sigma_full_mScm  # σ_grain
    print(f"\nτ_eff² (역산): {tau_eff2.min():.2f} ~ {tau_eff2.max():.2f} (문헌: 4~5)")
    print(f"  τ_geo²: {(tau**2).min():.2f} ~ {(tau**2).max():.2f}")
    print(f"  R_contact = τ_eff²/τ_geo²: {(tau_eff2/(tau**2)).min():.2f} ~ {(tau_eff2/(tau**2)).max():.2f}")

    # Summary
    print("\n" + "="*70)
    print("RANKING")
    print("="*70)
    for rank, (name, r2, params, detail) in enumerate(sorted(model_results, key=lambda x: -x[1]), 1):
        print(f"  {rank}. {name:45s} R²={r2:.4f}  ({params}p)  {detail}")

    # Per-case table
    print("\n" + "="*70)
    print("PER-CASE TABLE")
    print("="*70)
    print(f"{'Name':25s} {'GB_d':>6s} {'T':>6s} {'τ':>5s} {'φ_SE':>6s} {'R_brug':>7s} {'σ_full':>8s} {'τ_eff²':>7s}")
    print("-"*75)
    for r in sorted(results, key=lambda x: x['sigma_full_mScm']):
        te2 = 3.0 * r['phi_se'] / r['sigma_full_mScm'] if r['sigma_full_mScm'] > 0 else 0
        print(f"  {r['name']:23s} {r['gb_d']:6.2f} {r['T']:6.0f} {r['tau']:5.2f} {r['phi_se']:6.3f} {r['R_brug']:7.2f} {r['sigma_full_mScm']:8.4f} {te2:7.2f}")


if __name__ == '__main__':
    main()
