"""
Ultimate scaling law finder.
Uses ALL DEM microstructural descriptors to find σ_eff relationships.

Approach:
1. Load ALL variables from full_metrics.json + network_conductivity.json
2. Correlation matrix → identify key drivers
3. Decompose: σ_eff = σ_brug × (1/R_contact)
4. Fit R_contact with every possible combination
5. Fit σ_eff directly with all variables
6. Find the physically meaningful minimum model
"""
import json, os, sys, numpy as np, warnings
from scipy import stats
from scipy.optimize import curve_fit
from itertools import combinations

warnings.filterwarnings('ignore')

WEBAPP = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'webapp')


def load_all_data():
    """Load and merge network results + full metrics for all cases."""
    with open(os.path.join(WEBAPP, 'results', 'network_conductivity_all.json')) as f:
        net_data = json.load(f)

    # Deduplicate
    seen = {}
    for d in net_data:
        name = d.get('name', '')
        if name not in seen:
            seen[name] = d
    unique = list(seen.values())

    rows = []
    for nd in unique:
        case_id = nd.get('case_id', '')
        if nd.get('sigma_full') is None:
            continue

        if case_id.startswith('archive:'):
            mp = os.path.join(WEBAPP, 'archive', case_id[8:], 'full_metrics.json')
        else:
            mp = os.path.join(WEBAPP, 'results', case_id, 'full_metrics.json')
        if not os.path.exists(mp):
            continue

        with open(mp) as f:
            m = json.load(f)

        gb_d = m.get('gb_density_mean', 0)
        T = m.get('thickness_um', 0)
        tau = m.get('tortuosity_recommended', m.get('tortuosity_mean', 0))
        if gb_d <= 0 or T <= 0 or tau <= 0:
            continue

        rows.append({
            'name': nd['name'],
            # Network solver
            'sigma_full': nd.get('sigma_full_mScm', nd['sigma_full'] * 1.3),
            'sigma_bulk_net': nd.get('sigma_bulk_net_mScm', 0),
            'R_brug': nd['R_brug_over_full'],
            'bulk_frac': nd['bulk_resistance_fraction'],
            # Bruggeman
            'phi_se': m.get('phi_se', 0),
            'f_perc': m.get('percolation_pct', 0) / 100,
            'tau': tau,
            'porosity': m.get('porosity', 0),
            'T': T,
            # Contact quality
            'gb_d': gb_d,
            'hop_area': m.get('path_hop_area_mean', 0),
            'bottleneck': m.get('path_hop_area_min_mean', 0),
            'g_path': m.get('path_conductance_mean', 0),
            # Connectivity
            'cn': m.get('se_se_cn', 0),
            'n_clusters': m.get('n_components', 0),
            # Interface
            'se_se_total': m.get('area_SE_SE_total', 0),
            # Stress (relative only)
            'stress_cv': m.get('stress_cv', 0),
            # Derived
            'tau2': tau**2,
            'gb_d2': gb_d**2,
            'gb_d2T': gb_d**2 * T,
            'constr_ratio': (1 - nd['bulk_resistance_fraction']) / nd['bulk_resistance_fraction']
                            if nd['bulk_resistance_fraction'] > 0 else 0,
            'sigma_brug': 1.3 * m.get('phi_se', 0) * m.get('percolation_pct', 0) / 100 / tau**2,
        })
    return rows


def correlation_analysis(rows):
    """Find what correlates with σ_full and R_contact."""
    targets = ['sigma_full', 'R_brug']
    features = ['phi_se', 'f_perc', 'tau', 'T', 'gb_d', 'cn', 'hop_area',
                'bottleneck', 'bulk_frac', 'porosity', 'stress_cv',
                'gb_d2T', 'constr_ratio', 'se_se_total', 'tau2']

    print("\n" + "="*70)
    print("CORRELATION ANALYSIS")
    print("="*70)

    for target in targets:
        print(f"\n--- Correlations with {target} ---")
        y = np.array([r[target] for r in rows])
        corrs = []
        for feat in features:
            x = np.array([r[feat] for r in rows])
            if np.std(x) == 0:
                continue
            valid = (x > 0) & (y > 0) & np.isfinite(x) & np.isfinite(y)
            if valid.sum() < 5:
                continue
            # Pearson on log-log
            try:
                r_log = np.corrcoef(np.log(x[valid]), np.log(y[valid]))[0, 1]
                r_lin = np.corrcoef(x[valid], y[valid])[0, 1]
                corrs.append((feat, r_lin, r_log))
            except:
                pass

        corrs.sort(key=lambda x: abs(x[2]), reverse=True)
        print(f"  {'Variable':20s} {'r(linear)':>10s} {'r(log-log)':>10s}")
        for feat, r_lin, r_log in corrs:
            marker = " ★" if abs(r_log) > 0.7 else ""
            print(f"  {feat:20s} {r_lin:10.3f} {r_log:10.3f}{marker}")


def fit_R_contact(rows):
    """Fit R_contact = σ_bulk_net / σ_full with various models."""
    n = len(rows)
    R = np.array([r['R_brug'] for r in rows])
    gb_d = np.array([r['gb_d'] for r in rows])
    T = np.array([r['T'] for r in rows])
    cn = np.array([r['cn'] for r in rows])
    hop = np.array([r['hop_area'] for r in rows])
    bn = np.array([r['bottleneck'] for r in rows])
    bf = np.array([r['bulk_frac'] for r in rows])
    cr = np.array([r['constr_ratio'] for r in rows])
    phi = np.array([r['phi_se'] for r in rows])
    tau = np.array([r['tau'] for r in rows])

    ss_tot = np.sum((R - np.mean(R))**2)
    ss_tot_log = np.sum((np.log(R) - np.mean(np.log(R)))**2)

    print("\n" + "="*70)
    print(f"R_CONTACT FITTING (n={n})")
    print("="*70)

    results = []

    # 1. R = a×constr_ratio + b (constriction/bulk ratio)
    s, i, r, _, _ = stats.linregress(cr, R)
    print(f"\n1. R = {s:.3f}×(constr/bulk) + {i:.3f}  |  R²={r**2:.4f}")
    results.append(('constr/bulk ratio (linear)', r**2, 2))

    # 2. log(R) = a×log(constr_ratio) + b
    s, i, r, _, _ = stats.linregress(np.log(cr), np.log(R))
    print(f"2. R ∝ (constr/bulk)^{s:.3f}  |  R²={r**2:.4f}")
    results.append(('constr/bulk ratio (power)', r**2, 2))

    # 3. R = f(hop_area)
    valid = hop > 0
    s, i, r, _, _ = stats.linregress(np.log(hop[valid]), np.log(R[valid]))
    print(f"3. R ∝ hop_area^{s:.3f}  |  R²={r**2:.4f}")
    results.append(('hop_area (power)', r**2, 2))

    # 4. R = f(CN)
    valid = cn > 0
    s, i, r, _, _ = stats.linregress(np.log(cn[valid]), np.log(R[valid]))
    print(f"4. R ∝ CN^{s:.3f}  |  R²={r**2:.4f}")
    results.append(('CN (power)', r**2, 2))

    # 5. R = f(hop_area × CN)
    valid = (hop > 0) & (cn > 0)
    x = hop[valid] * cn[valid]
    s, i, r, _, _ = stats.linregress(np.log(x), np.log(R[valid]))
    print(f"5. R ∝ (hop×CN)^{s:.3f}  |  R²={r**2:.4f}")
    results.append(('hop×CN (power)', r**2, 2))

    # 6. R = 1 + β/(hop_area^a)
    try:
        def m6(h, beta, a): return 1 + beta / h**a
        p, _ = curve_fit(m6, hop[valid], R[valid], p0=[0.5, 0.5], maxfev=10000)
        pred = m6(hop[valid], *p)
        r2 = 1 - np.sum((R[valid]-pred)**2)/np.sum((R[valid]-np.mean(R[valid]))**2)
        print(f"6. R = 1 + {p[0]:.4f}/hop^{p[1]:.4f}  |  R²={r2:.4f}")
        results.append(('1 + β/hop^a', r2, 2))
    except:
        print("6. FAILED")

    # 7. R = 1 + β/(hop×CN)^a
    try:
        x = hop[valid] * cn[valid]
        def m7(x, beta, a): return 1 + beta / x**a
        p, _ = curve_fit(m7, x, R[valid], p0=[1.0, 0.3], maxfev=10000)
        pred = m7(x, *p)
        r2 = 1 - np.sum((R[valid]-pred)**2)/np.sum((R[valid]-np.mean(R[valid]))**2)
        print(f"7. R = 1 + {p[0]:.4f}/(hop×CN)^{p[1]:.4f}  |  R²={r2:.4f}")
        results.append(('1 + β/(hop×CN)^a', r2, 2))
    except:
        print("7. FAILED")

    # 8. Multi: log(R) = a×log(hop) + b×log(CN) + c
    valid = (hop > 0) & (cn > 0)
    X = np.column_stack([np.log(hop[valid]), np.log(cn[valid]), np.ones(valid.sum())])
    b, _, _, _ = np.linalg.lstsq(X, np.log(R[valid]), rcond=None)
    pred = X @ b
    r2 = 1 - np.sum((np.log(R[valid])-pred)**2)/np.sum((np.log(R[valid])-np.mean(np.log(R[valid])))**2)
    print(f"8. R ∝ hop^{b[0]:.3f} × CN^{b[1]:.3f}  |  R²={r2:.4f}")
    results.append(('hop^a × CN^b', r2, 3))

    # 9. Multi: log(R) = a×log(hop) + b×log(CN) + c×log(GB_d) + d
    valid = (hop > 0) & (cn > 0)
    X = np.column_stack([np.log(hop[valid]), np.log(cn[valid]),
                         np.log(gb_d[valid]), np.ones(valid.sum())])
    b, _, _, _ = np.linalg.lstsq(X, np.log(R[valid]), rcond=None)
    pred = X @ b
    r2 = 1 - np.sum((np.log(R[valid])-pred)**2)/np.sum((np.log(R[valid])-np.mean(np.log(R[valid])))**2)
    print(f"9. R ∝ hop^{b[0]:.3f} × CN^{b[1]:.3f} × GB_d^{b[2]:.3f}  |  R²={r2:.4f}")
    results.append(('hop^a × CN^b × GB_d^c', r2, 4))

    # 10. R = f(bulk_frac) — direct from network decomposition
    s, i, r, _, _ = stats.linregress(bf, R)
    print(f"10. R = {s:.2f}×bulk_frac + {i:.2f}  |  R²={r**2:.4f}")
    results.append(('bulk_frac (linear)', r**2, 2))

    # 11. R ∝ bulk_frac^a
    s, i, r, _, _ = stats.linregress(np.log(bf), np.log(R))
    print(f"11. R ∝ bulk_frac^{s:.3f}  |  R²={r**2:.4f}")
    results.append(('bulk_frac (power)', r**2, 2))

    # Summary
    print(f"\n{'--- R_contact Ranking ---':^60}")
    for rank, (name, r2, p) in enumerate(sorted(results, key=lambda x: -x[1]), 1):
        print(f"  {rank:2d}. {name:35s} R²={r2:.4f} ({p}p)")

    return results


def fit_sigma_eff(rows):
    """Fit σ_eff directly with all variables."""
    n = len(rows)
    sf = np.array([r['sigma_full'] for r in rows])
    sb = np.array([r['sigma_brug'] for r in rows])
    phi = np.array([r['phi_se'] for r in rows])
    fp = np.array([r['f_perc'] for r in rows])
    tau = np.array([r['tau'] for r in rows])
    T = np.array([r['T'] for r in rows])
    gb_d = np.array([r['gb_d'] for r in rows])
    cn = np.array([r['cn'] for r in rows])
    hop = np.array([r['hop_area'] for r in rows])
    bn = np.array([r['bottleneck'] for r in rows])
    bf = np.array([r['bulk_frac'] for r in rows])

    log_sf = np.log(sf)
    ss_tot = np.sum((log_sf - np.mean(log_sf))**2)

    print("\n" + "="*70)
    print(f"σ_eff DIRECT FITTING (n={n})")
    print("="*70)

    results = []

    # 1. σ = c × φ^n (pure Bruggeman)
    s, i, r, _, _ = stats.linregress(np.log(phi), log_sf)
    print(f"\n1. σ ∝ φ_SE^{s:.2f}  |  R²={r**2:.4f}")
    results.append(('φ_SE^n', r**2, 2))

    # 2. σ = σ_brug / R → σ = σ_brug × f(contact)
    # σ = c × φ × f_perc / τ² × hop^a × CN^b
    valid = (hop > 0) & (cn > 0)
    X = np.column_stack([np.log(phi[valid]), np.log(fp[valid]),
                         np.log(tau[valid]), np.log(hop[valid]),
                         np.log(cn[valid]), np.ones(valid.sum())])
    b, _, _, _ = np.linalg.lstsq(X, log_sf[valid], rcond=None)
    pred = X @ b
    r2 = 1 - np.sum((log_sf[valid]-pred)**2)/np.sum((log_sf[valid]-np.mean(log_sf[valid]))**2)
    print(f"2. σ ∝ φ^{b[0]:.2f} × f_perc^{b[1]:.2f} × τ^{b[2]:.2f} × hop^{b[3]:.2f} × CN^{b[4]:.2f}")
    print(f"   R²={r2:.4f}")
    results.append(('φ×f_perc×τ×hop×CN', r2, 6))

    # 3. σ = σ_brug × (1 + β/hop^a)^-1  (Bruggeman + contact correction)
    # log(σ) = log(σ_brug) - log(R_contact)
    # where R_contact = f(hop, CN)
    log_sb = np.log(sb)
    residual = log_sf - log_sb  # This is -log(R_contact)
    # fit residual with hop and CN
    valid = (hop > 0) & (cn > 0) & np.isfinite(log_sb)
    X = np.column_stack([np.log(hop[valid]), np.log(cn[valid]), np.ones(valid.sum())])
    b, _, _, _ = np.linalg.lstsq(X, residual[valid], rcond=None)
    pred_full = log_sb[valid] + X @ b
    r2 = 1 - np.sum((log_sf[valid]-pred_full)**2)/np.sum((log_sf[valid]-np.mean(log_sf[valid]))**2)
    print(f"\n3. σ = σ_brug × hop^{b[0]:.3f} × CN^{b[1]:.3f} × exp({b[2]:.3f})")
    print(f"   = σ_bulk × φ_SE × f_perc / τ² × hop^{b[0]:.3f} × CN^{b[1]:.3f} × {np.exp(b[2]):.4f}")
    print(f"   R²={r2:.4f}")
    results.append(('σ_brug × hop^a × CN^b', r2, 4))

    # 4. Add GB_d
    X = np.column_stack([np.log(hop[valid]), np.log(cn[valid]),
                         np.log(gb_d[valid]), np.ones(valid.sum())])
    b, _, _, _ = np.linalg.lstsq(X, residual[valid], rcond=None)
    pred_full = log_sb[valid] + X @ b
    r2 = 1 - np.sum((log_sf[valid]-pred_full)**2)/np.sum((log_sf[valid]-np.mean(log_sf[valid]))**2)
    print(f"\n4. σ = σ_brug × hop^{b[0]:.3f} × CN^{b[1]:.3f} × GB_d^{b[2]:.3f}")
    print(f"   R²={r2:.4f}")
    results.append(('σ_brug × hop^a × CN^b × GB_d^c', r2, 5))

    # 5. Add bottleneck
    valid2 = valid & (bn > 0)
    X = np.column_stack([np.log(hop[valid2]), np.log(cn[valid2]),
                         np.log(bn[valid2]), np.ones(valid2.sum())])
    b, _, _, _ = np.linalg.lstsq(X, residual[valid2], rcond=None)
    pred_full = log_sb[valid2] + X @ b
    r2 = 1 - np.sum((log_sf[valid2]-pred_full)**2)/np.sum((log_sf[valid2]-np.mean(log_sf[valid2]))**2)
    print(f"\n5. σ = σ_brug × hop^{b[0]:.3f} × CN^{b[1]:.3f} × bottleneck^{b[2]:.3f}")
    print(f"   R²={r2:.4f}")
    results.append(('σ_brug × hop^a × CN^b × BN^c', r2, 5))

    # 6. Kitchen sink: all variables
    valid3 = valid2 & (T > 0)
    X = np.column_stack([np.log(hop[valid3]), np.log(cn[valid3]),
                         np.log(gb_d[valid3]), np.log(bn[valid3]),
                         np.log(T[valid3]), np.ones(valid3.sum())])
    b, _, _, _ = np.linalg.lstsq(X, residual[valid3], rcond=None)
    pred_full = log_sb[valid3] + X @ b
    r2 = 1 - np.sum((log_sf[valid3]-pred_full)**2)/np.sum((log_sf[valid3]-np.mean(log_sf[valid3]))**2)
    print(f"\n6. σ = σ_brug × hop^{b[0]:.3f} × CN^{b[1]:.3f} × GB_d^{b[2]:.3f} × BN^{b[3]:.3f} × T^{b[4]:.3f}")
    print(f"   R²={r2:.4f}")
    results.append(('Kitchen sink (6p)', r2, 7))

    # 7. Simplest meaningful: σ = σ_brug / (1 + β×something)
    # From decomposition: R_contact ≈ 1 + constr/bulk
    # constr/bulk ∝ 1/(hop_area × something)
    # Try: σ = σ_brug × hop^a (simplest contact correction)
    valid = hop > 0
    s, i, r, _, _ = stats.linregress(np.log(hop[valid]), residual[valid])
    r2_full = 1 - np.sum((log_sf[valid] - (log_sb[valid] + s*np.log(hop[valid]) + i))**2) / \
                   np.sum((log_sf[valid] - np.mean(log_sf[valid]))**2)
    print(f"\n7. σ = σ_brug × {np.exp(i):.4f} × hop^{s:.3f}  (simplest)")
    print(f"   = σ_bulk × φ_SE × f_perc / τ² × {np.exp(i):.4f} × hop_area^{s:.3f}")
    print(f"   R²={r2_full:.4f}")
    results.append(('σ_brug × hop^a (simplest)', r2_full, 3))

    # Summary
    print(f"\n{'--- σ_eff Ranking ---':^60}")
    for rank, (name, r2, p) in enumerate(sorted(results, key=lambda x: -x[1]), 1):
        star = " ★" if r2 > 0.9 else ""
        print(f"  {rank:2d}. {name:40s} R²={r2:.4f} ({p}p){star}")


def print_final_recommendation(rows):
    """Print the recommended final model."""
    print("\n" + "="*70)
    print("FINAL RECOMMENDATION")
    print("="*70)
    print("""
The ultimate model decomposes σ_eff into physically meaningful terms:

  σ_eff = σ_bulk × (φ_SE × f_perc / τ²) × (C × hop_area^a × CN^b)
          ├── Bruggeman ──────────┤   ├── Contact correction ──┤

  Bruggeman term:  captures geometry (tortuosity, volume fraction)
  Contact term:    captures inter-particle contact quality

  Where:
    φ_SE    → determined by AM:SE mass ratio
    f_perc  → determined by SE connectivity (SE size, composition)
    τ       → determined by packing geometry (P:S ratio, compaction)
    hop_area → determined by SE size + compaction pressure
    CN      → determined by SE size + composition + packing

  All variables directly computable from DEM output.
  σ_eff in mS/cm, directly comparable to EIS experiments.
""")


def main():
    rows = load_all_data()
    print(f"Loaded {len(rows)} unique cases")

    # 1. Correlation analysis
    correlation_analysis(rows)

    # 2. R_contact fitting
    fit_R_contact(rows)

    # 3. σ_eff direct fitting
    fit_sigma_eff(rows)

    # 4. Recommendation
    print_final_recommendation(rows)


if __name__ == '__main__':
    main()
