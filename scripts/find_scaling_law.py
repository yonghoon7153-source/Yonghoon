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
            'sigma_full': nd.get('sigma_full_mScm', nd['sigma_full'] * 3.0),
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
            'sigma_brug': 3.0 * m.get('phi_se', 0) * m.get('percolation_pct', 0) / 100 / tau**2,
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

    # ── SE size analysis: why does C vary with SE size? ──
    print(f"\n{'='*70}")
    print(f"SE SIZE ANALYSIS")
    print(f"{'='*70}")

    # Compute d_SE from hop_area (proxy: a_contact ∝ √A_hop, a ∝ d_SE via Hertz)
    # Or use GB_d as inverse proxy for d_SE
    # Test: does adding d_SE-related terms improve the model?

    # 8. σ_brug × hop^a × CN^b × GB_d^c with various fixed exponent combos
    valid = (hop > 0) & (cn > 0)
    log_sb = np.log(sb)
    residual = log_sf - log_sb

    # 8a. Free fit all 3 (baseline from model 4)
    X = np.column_stack([np.log(hop[valid]), np.log(cn[valid]),
                         np.log(gb_d[valid]), np.ones(valid.sum())])
    b, _, _, _ = np.linalg.lstsq(X, residual[valid], rcond=None)
    pred = log_sb[valid] + X @ b
    r2 = 1 - np.sum((log_sf[valid]-pred)**2)/np.sum((log_sf[valid]-np.mean(log_sf[valid]))**2)
    print(f"\n8a. Free: hop^{b[0]:.3f} × CN^{b[1]:.3f} × GB_d^{b[2]:.3f} × {np.exp(b[3]):.4f}")
    print(f"    R²={r2:.4f}")

    # 8b. Fixed (0.5, 2, 4/3) + C
    log_rhs = 0.5*np.log(hop[valid]) + 2*np.log(cn[valid]) + 4/3*np.log(gb_d[valid])
    ln_C = np.mean(residual[valid] - log_rhs)
    pred = log_sb[valid] + ln_C + log_rhs
    r2_fixed = 1 - np.sum((log_sf[valid]-pred)**2)/np.sum((log_sf[valid]-np.mean(log_sf[valid]))**2)
    print(f"8b. Fixed(0.5,2,4/3): C={np.exp(ln_C):.4f}, R²={r2_fixed:.4f}")

    # 8c. Test various GB_d exponents with hop=0.5, CN=2 fixed
    print(f"\n  GB_d exponent sweep (hop=0.5, CN=2 fixed):")
    for c_test in [0.5, 0.75, 1.0, 1.2, 1.24, 1.33, 1.5, 2.0]:
        log_rhs_t = 0.5*np.log(hop[valid]) + 2*np.log(cn[valid]) + c_test*np.log(gb_d[valid])
        ln_C_t = np.mean(residual[valid] - log_rhs_t)
        pred_t = log_sb[valid] + ln_C_t + log_rhs_t
        r2_t = 1 - np.sum((log_sf[valid]-pred_t)**2)/np.sum((log_sf[valid]-np.mean(log_sf[valid]))**2)
        print(f"    GB_d^{c_test:.2f}: C={np.exp(ln_C_t):.4f}, R²={r2_t:.4f}")

    # 8d. Test hop exponent sweep with CN=2, GB_d=4/3 fixed
    print(f"\n  hop exponent sweep (CN=2, GB_d=4/3 fixed):")
    for a_test in [0.3, 0.4, 0.5, 0.55, 0.58, 0.6, 0.7, 0.8, 1.0]:
        log_rhs_t = a_test*np.log(hop[valid]) + 2*np.log(cn[valid]) + 4/3*np.log(gb_d[valid])
        ln_C_t = np.mean(residual[valid] - log_rhs_t)
        pred_t = log_sb[valid] + ln_C_t + log_rhs_t
        r2_t = 1 - np.sum((log_sf[valid]-pred_t)**2)/np.sum((log_sf[valid]-np.mean(log_sf[valid]))**2)
        print(f"    hop^{a_test:.2f}: C={np.exp(ln_C_t):.4f}, R²={r2_t:.4f}")

    # 9. Try normalized combinations that cancel SE size
    # A_hop × GB_d² ∝ d_SE × (1/d_SE)² = 1/d_SE  (not cancel!)
    # A_hop × GB_d ∝ const (cancel!)
    # √A_hop × GB_d ∝ √d_SE × (1/d_SE) = 1/√d_SE
    # A_hop^(2/3) × GB_d ∝ d_SE^(2/3) × (1/d_SE) = d_SE^(-1/3)
    print(f"\n  SE-size-normalized combination tests:")

    # 9a. (A_hop × GB_d)^a × CN^b  (A_hop×GB_d should cancel d_SE)
    x_combo = hop[valid] * gb_d[valid]  # should be ~SE-size-independent
    X9 = np.column_stack([np.log(x_combo), np.log(cn[valid]), np.ones(valid.sum())])
    b9, _, _, _ = np.linalg.lstsq(X9, residual[valid], rcond=None)
    pred9 = log_sb[valid] + X9 @ b9
    r2_9 = 1 - np.sum((log_sf[valid]-pred9)**2)/np.sum((log_sf[valid]-np.mean(log_sf[valid]))**2)
    print(f"  9a. (A_hop×GB_d)^{b9[0]:.3f} × CN^{b9[1]:.3f}: R²={r2_9:.4f}")

    # 9b. (A_hop × GB_d²)^a × CN^b
    x_combo2 = hop[valid] * gb_d[valid]**2
    X9b = np.column_stack([np.log(x_combo2), np.log(cn[valid]), np.ones(valid.sum())])
    b9b, _, _, _ = np.linalg.lstsq(X9b, residual[valid], rcond=None)
    pred9b = log_sb[valid] + X9b @ b9b
    r2_9b = 1 - np.sum((log_sf[valid]-pred9b)**2)/np.sum((log_sf[valid]-np.mean(log_sf[valid]))**2)
    print(f"  9b. (A_hop×GB_d²)^{b9b[0]:.3f} × CN^{b9b[1]:.3f}: R²={r2_9b:.4f}")

    # 9c. (√A_hop × GB_d)^a × CN^b  (= a_contact × GB_d, physically: constriction × density)
    x_combo3 = np.sqrt(hop[valid]) * gb_d[valid]
    X9c = np.column_stack([np.log(x_combo3), np.log(cn[valid]), np.ones(valid.sum())])
    b9c, _, _, _ = np.linalg.lstsq(X9c, residual[valid], rcond=None)
    pred9c = log_sb[valid] + X9c @ b9c
    r2_9c = 1 - np.sum((log_sf[valid]-pred9c)**2)/np.sum((log_sf[valid]-np.mean(log_sf[valid]))**2)
    print(f"  9c. (√A_hop×GB_d)^{b9c[0]:.3f} × CN^{b9c[1]:.3f}: R²={r2_9c:.4f}")

    # 9d. bottleneck instead of A_hop
    valid2 = valid & (bn > 0)
    if valid2.sum() > 5:
        X9d = np.column_stack([np.log(bn[valid2]), np.log(cn[valid2]),
                               np.log(gb_d[valid2]), np.ones(valid2.sum())])
        b9d, _, _, _ = np.linalg.lstsq(X9d, (log_sf - log_sb)[valid2], rcond=None)
        pred9d = log_sb[valid2] + X9d @ b9d
        r2_9d = 1 - np.sum(((log_sf-log_sb)[valid2]-X9d@b9d)**2)/np.sum(((log_sf)[valid2]-np.mean((log_sf)[valid2]))**2)
        print(f"  9d. BN^{b9d[0]:.3f} × CN^{b9d[1]:.3f} × GB_d^{b9d[2]:.3f}: R²={r2_9d:.4f}  (bottleneck)")

    # 10. Per-SE-size R² (how good is the model within each SE size?)
    print(f"\n  Per-SE-size R² (fixed 0.5, 2, 4/3):")
    # Group by GB_d: SE 0.5μm → GB_d > 1.0, SE 1.0μm → 0.6~0.8, SE 1.5μm → GB_d < 0.6
    for se_label, gb_lo, gb_hi in [("SE 0.5μm", 1.0, 3.0), ("SE 1.0μm", 0.6, 1.0), ("SE 1.5μm", 0.0, 0.6)]:
        mask = valid & (gb_d >= gb_lo) & (gb_d < gb_hi)
        if mask.sum() < 3:
            continue
        log_rhs_m = 0.5*np.log(hop[mask]) + 2*np.log(cn[mask]) + 4/3*np.log(gb_d[mask])
        ln_C_m = np.mean(residual[mask] - log_rhs_m)
        pred_m = log_sb[mask] + ln_C_m + log_rhs_m
        ss_res_m = np.sum((log_sf[mask] - pred_m)**2)
        ss_tot_m = np.sum((log_sf[mask] - np.mean(log_sf[mask]))**2)
        r2_m = 1 - ss_res_m / ss_tot_m if ss_tot_m > 0 else 0
        print(f"    {se_label} (n={mask.sum()}, GB_d={gb_d[mask].min():.2f}~{gb_d[mask].max():.2f}): C={np.exp(ln_C_m):.4f}, R²={r2_m:.4f}")

    # ── NEW BEAUTIFUL FORMULA: (A_hop × GB_d²)^(3/5) × CN² ──
    print(f"\n{'='*70}")
    print(f"NEW FORMULA: σ = σ_brug × C × (A_hop × GB_d²)^(3/5) × CN²")
    print(f"{'='*70}")

    valid = (hop > 0) & (cn > 0)
    combo = hop[valid] * gb_d[valid]**2  # combined variable

    # Fixed 3/5 + C only (1 free param)
    log_rhs_new = 3/5 * np.log(combo) + 2 * np.log(cn[valid])
    ln_C_new = np.mean(residual[valid] - log_rhs_new)
    pred_new = log_sb[valid] + ln_C_new + log_rhs_new
    r2_new = 1 - np.sum((log_sf[valid]-pred_new)**2)/np.sum((log_sf[valid]-np.mean(log_sf[valid]))**2)
    print(f"\n  Fixed (3/5, 2): C={np.exp(ln_C_new):.4f}, R²={r2_new:.4f}")

    # Free fit for comparison
    X_new = np.column_stack([np.log(combo), np.log(cn[valid]), np.ones(valid.sum())])
    b_new, _, _, _ = np.linalg.lstsq(X_new, residual[valid], rcond=None)
    pred_new_free = log_sb[valid] + X_new @ b_new
    r2_new_free = 1 - np.sum((log_sf[valid]-pred_new_free)**2)/np.sum((log_sf[valid]-np.mean(log_sf[valid]))**2)
    print(f"  Free fit: (A_hop×GB_d²)^{b_new[0]:.3f} × CN^{b_new[1]:.3f}: C={np.exp(b_new[2]):.4f}, R²={r2_new_free:.4f}")

    # Exponent sweep for (A_hop × GB_d²)
    print(f"\n  (A_hop×GB_d²) exponent sweep (CN=2 fixed):")
    for e_test in [0.4, 0.5, 3/5, 0.65, 0.7, 0.8]:
        log_rhs_t = e_test * np.log(combo) + 2 * np.log(cn[valid])
        ln_C_t = np.mean(residual[valid] - log_rhs_t)
        pred_t = log_sb[valid] + ln_C_t + log_rhs_t
        r2_t = 1 - np.sum((log_sf[valid]-pred_t)**2)/np.sum((log_sf[valid]-np.mean(log_sf[valid]))**2)
        label = " ← 3/5" if abs(e_test - 3/5) < 0.001 else ""
        print(f"    ^{e_test:.3f}: C={np.exp(ln_C_t):.4f}, R²={r2_t:.4f}{label}")

    # Per-SE-size with new formula
    print(f"\n  Per-SE-size R² (NEW: (A_hop×GB_d²)^(3/5) × CN²):")
    for se_label, gb_lo, gb_hi in [("SE 0.5μm", 1.0, 3.0), ("SE 1.0μm", 0.6, 1.0), ("SE 1.5μm", 0.0, 0.6)]:
        mask = valid & (gb_d >= gb_lo) & (gb_d < gb_hi)
        if mask.sum() < 3:
            continue
        combo_m = hop[mask] * gb_d[mask]**2
        log_rhs_m = 3/5 * np.log(combo_m) + 2 * np.log(cn[mask])
        ln_C_m = np.mean((log_sf - log_sb)[mask] - log_rhs_m)
        pred_m = log_sb[mask] + ln_C_m + log_rhs_m
        ss_res_m = np.sum((log_sf[mask] - pred_m)**2)
        ss_tot_m = np.sum((log_sf[mask] - np.mean(log_sf[mask]))**2)
        r2_m = 1 - ss_res_m / ss_tot_m if ss_tot_m > 0 else 0
        print(f"    {se_label} (n={mask.sum()}): C={np.exp(ln_C_m):.4f}, R²={r2_m:.4f}")

    # Per-case accuracy
    s_pred_new = np.exp(pred_new)
    s_actual_new = np.exp(log_sf[valid])
    errors_new = np.abs(s_pred_new - s_actual_new) / s_actual_new * 100
    within_20 = np.sum(errors_new < 20)
    print(f"\n  Mean |error|: {np.mean(errors_new):.1f}%")
    print(f"  Within 20%: {within_20}/{len(errors_new)}")
    print(f"  Max |error|: {np.max(errors_new):.1f}%")

    # Compare old vs new
    print(f"\n  {'─'*50}")
    print(f"  OLD: √A_hop × CN² × GB_d^(4/3), R²={r2_fixed:.4f}")
    print(f"  NEW: (A_hop×GB_d²)^(3/5) × CN², R²={r2_new:.4f}")
    print(f"  ΔR² = {r2_new - r2_fixed:+.4f}")

    # ── UPGRADE ATTEMPTS ──
    print(f"\n{'='*70}")
    print(f"UPGRADE ATTEMPTS (beyond R²=0.923)")
    print(f"{'='*70}")

    valid2 = valid & (bn > 0)
    residual2 = (log_sf - log_sb)[valid2]

    # U1: (A_hop×GB_d²)^(3/5) × CN² × BN^d  (add bottleneck)
    combo2 = hop[valid2] * gb_d[valid2]**2
    X_u1 = np.column_stack([3/5*np.log(combo2), 2*np.log(cn[valid2]),
                            np.log(bn[valid2]), np.ones(valid2.sum())])
    b_u1, _, _, _ = np.linalg.lstsq(X_u1, residual2, rcond=None)
    pred_u1 = log_sb[valid2] + X_u1 @ b_u1
    r2_u1 = 1 - np.sum((log_sf[valid2]-pred_u1)**2)/np.sum((log_sf[valid2]-np.mean(log_sf[valid2]))**2)
    print(f"\n  U1: (A_hop×GB_d²)^(3/5) × CN² × BN^{b_u1[2]:.3f}")
    print(f"      R²={r2_u1:.4f} (+BN, but 3/5 and 2 fixed)")

    # U1b: (A_hop×GB_d²)^a × CN^b × BN^c (all free)
    X_u1b = np.column_stack([np.log(combo2), np.log(cn[valid2]),
                             np.log(bn[valid2]), np.ones(valid2.sum())])
    b_u1b, _, _, _ = np.linalg.lstsq(X_u1b, residual2, rcond=None)
    pred_u1b = log_sb[valid2] + X_u1b @ b_u1b
    r2_u1b = 1 - np.sum((log_sf[valid2]-pred_u1b)**2)/np.sum((log_sf[valid2]-np.mean(log_sf[valid2]))**2)
    print(f"  U1b: (A_hop×GB_d²)^{b_u1b[0]:.3f} × CN^{b_u1b[1]:.3f} × BN^{b_u1b[2]:.3f}")
    print(f"       R²={r2_u1b:.4f} (all free)")

    # U2: BN instead of A_hop: (BN × GB_d²)^a × CN^b
    combo_bn = bn[valid2] * gb_d[valid2]**2
    X_u2 = np.column_stack([np.log(combo_bn), np.log(cn[valid2]), np.ones(valid2.sum())])
    b_u2, _, _, _ = np.linalg.lstsq(X_u2, residual2, rcond=None)
    pred_u2 = log_sb[valid2] + X_u2 @ b_u2
    r2_u2 = 1 - np.sum((log_sf[valid2]-pred_u2)**2)/np.sum((log_sf[valid2]-np.mean(log_sf[valid2]))**2)
    print(f"\n  U2: (BN×GB_d²)^{b_u2[0]:.3f} × CN^{b_u2[1]:.3f}")
    print(f"      R²={r2_u2:.4f} (bottleneck replaces A_hop)")
    # Fixed 3/5 test
    log_rhs_u2f = 3/5*np.log(combo_bn) + 2*np.log(cn[valid2])
    ln_C_u2f = np.mean(residual2 - log_rhs_u2f)
    pred_u2f = log_sb[valid2] + ln_C_u2f + log_rhs_u2f
    r2_u2f = 1 - np.sum((log_sf[valid2]-pred_u2f)**2)/np.sum((log_sf[valid2]-np.mean(log_sf[valid2]))**2)
    print(f"  U2b: (BN×GB_d²)^(3/5) × CN² [fixed]: C={np.exp(ln_C_u2f):.4f}, R²={r2_u2f:.4f}")

    # U3: Geometric mean of A_hop and BN: (√(A_hop×BN) × GB_d²)^a × CN^b
    combo_geo = np.sqrt(hop[valid2] * bn[valid2]) * gb_d[valid2]**2
    X_u3 = np.column_stack([np.log(combo_geo), np.log(cn[valid2]), np.ones(valid2.sum())])
    b_u3, _, _, _ = np.linalg.lstsq(X_u3, residual2, rcond=None)
    pred_u3 = log_sb[valid2] + X_u3 @ b_u3
    r2_u3 = 1 - np.sum((log_sf[valid2]-pred_u3)**2)/np.sum((log_sf[valid2]-np.mean(log_sf[valid2]))**2)
    print(f"\n  U3: (√(A_hop×BN)×GB_d²)^{b_u3[0]:.3f} × CN^{b_u3[1]:.3f}")
    print(f"      R²={r2_u3:.4f} (geometric mean of A_hop & BN)")

    # U4: G_path (path conductance) directly — most physically direct
    g_path = np.array([r['g_path'] for r in rows])
    valid3 = valid & (g_path > 0)
    if valid3.sum() > 5:
        X_u4 = np.column_stack([np.log(g_path[valid3]), np.log(cn[valid3]),
                                np.log(gb_d[valid3]), np.ones(valid3.sum())])
        b_u4, _, _, _ = np.linalg.lstsq(X_u4, (log_sf-log_sb)[valid3], rcond=None)
        pred_u4 = log_sb[valid3] + X_u4 @ b_u4
        r2_u4 = 1 - np.sum((log_sf[valid3]-pred_u4)**2)/np.sum((log_sf[valid3]-np.mean(log_sf[valid3]))**2)
        print(f"\n  U4: G_path^{b_u4[0]:.3f} × CN^{b_u4[1]:.3f} × GB_d^{b_u4[2]:.3f}")
        print(f"      R²={r2_u4:.4f} (G_path = harmonic mean conductance)")

        # U4b: G_path only + CN
        X_u4b = np.column_stack([np.log(g_path[valid3]), np.log(cn[valid3]), np.ones(valid3.sum())])
        b_u4b, _, _, _ = np.linalg.lstsq(X_u4b, (log_sf-log_sb)[valid3], rcond=None)
        pred_u4b = log_sb[valid3] + X_u4b @ b_u4b
        r2_u4b = 1 - np.sum((log_sf[valid3]-pred_u4b)**2)/np.sum((log_sf[valid3]-np.mean(log_sf[valid3]))**2)
        print(f"  U4b: G_path^{b_u4b[0]:.3f} × CN^{b_u4b[1]:.3f}: R²={r2_u4b:.4f}")

    # U5: se_se_total (total SE-SE contact area) — network-level metric
    se_total = np.array([r['se_se_total'] for r in rows])
    valid4 = valid & (se_total > 0)
    if valid4.sum() > 5:
        X_u5 = np.column_stack([np.log(se_total[valid4]), np.log(cn[valid4]), np.ones(valid4.sum())])
        b_u5, _, _, _ = np.linalg.lstsq(X_u5, (log_sf-log_sb)[valid4], rcond=None)
        pred_u5 = log_sb[valid4] + X_u5 @ b_u5
        r2_u5 = 1 - np.sum((log_sf[valid4]-pred_u5)**2)/np.sum((log_sf[valid4]-np.mean(log_sf[valid4]))**2)
        print(f"\n  U5: SE_total^{b_u5[0]:.3f} × CN^{b_u5[1]:.3f}: R²={r2_u5:.4f} (total SE-SE area)")

    # U6: (A_hop × GB_d²)^(3/5) × CN² with f_perc separated from σ_brug
    # σ = σ_grain × φ_SE / τ² × C × f_perc^d × (A_hop×GB_d²)^(3/5) × CN²
    sigma_brug_no_fperc = np.array([SIGMA_BULK * phi[i] / tau[i]**2 if tau[i] > 0 else 0 for i in range(n)])
    valid5 = valid & (fp > 0) & (sigma_brug_no_fperc > 0)
    if valid5.sum() > 5:
        combo5 = hop[valid5] * gb_d[valid5]**2
        log_sb5 = np.log(sigma_brug_no_fperc[valid5])
        X_u6 = np.column_stack([np.log(fp[valid5]), 3/5*np.log(combo5),
                                2*np.log(cn[valid5]), np.ones(valid5.sum())])
        b_u6, _, _, _ = np.linalg.lstsq(X_u6, (log_sf - log_sb5)[valid5], rcond=None)
        pred_u6 = log_sb5[valid5] + X_u6 @ b_u6
        r2_u6 = 1 - np.sum((log_sf[valid5]-pred_u6)**2)/np.sum((log_sf[valid5]-np.mean(log_sf[valid5]))**2)
        print(f"\n  U6: f_perc^{b_u6[0]:.3f} × (A_hop×GB_d²)^(3/5) × CN²")
        print(f"      R²={r2_u6:.4f} (f_perc separated from σ_brug)")

    # Summary
    print(f"\n{'--- σ_eff Ranking ---':^60}")
    results.append(('NEW: σ_brug×(A_hop×GB_d²)^(3/5)×CN²', r2_new, 2))
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
