"""
Symbolic Regression: Auto-discover mathematical formulas
========================================================
DEM input parameters → Electrode properties

Uses PySR (genetic programming) to evolve mathematical formulas
that best fit the data. Like our manual scaling law search but
testing MILLIONS of combinations automatically.

Usage:
  python3 scripts/symbolic_regression.py              # Find formulas for all targets
  python3 scripts/symbolic_regression.py --target sigma_ion  # Specific target only
"""

import json, os, sys, warnings
import numpy as np
from pathlib import Path

warnings.filterwarnings('ignore')

WEBAPP = Path(__file__).parent.parent / 'webapp'


def load_data():
    """Load DEM input → output data (same as ml_predictor v2)."""
    rows = []
    for base in [WEBAPP / 'results', WEBAPP / 'archive']:
        if not base.is_dir():
            continue
        for met_path in base.rglob('full_metrics.json'):
            case_dir = met_path.parent
            ip_path = case_dir / 'input_params.json'
            try:
                with open(met_path) as f:
                    m = json.load(f)
                ip = {}
                if ip_path.exists():
                    with open(ip_path) as f:
                        ip = json.load(f)
            except:
                continue

            if not m.get('phi_se') or not m.get('porosity'):
                continue

            scale = ip.get('scale', 1000)
            r_se = ip.get('r_SE', 0)
            d_se = r_se / scale * 1e6 * 2 if r_se > 0 else 0
            r_am_p = ip.get('r_AM_P', 0)
            d_am_p = r_am_p / scale * 1e6 * 2 if r_am_p > 0 else 0
            r_am_s = ip.get('r_AM_S', 0)
            d_am_s = r_am_s / scale * 1e6 * 2 if r_am_s > 0 else 0

            am_se_str = ip.get('am_se_ratio', '')
            am_pct = 0
            if ':' in str(am_se_str):
                try:
                    am_pct = float(str(am_se_str).split(':')[0])
                except:
                    pass

            ps = m.get('ps_ratio', '')
            ps_frac = 0.5
            if isinstance(ps, str) and ':' in ps:
                try:
                    parts = ps.split(':')
                    p, s = float(parts[0]), float(parts[1])
                    ps_frac = p / (p + s) if (p + s) > 0 else 0.5
                except:
                    pass

            box_x = ip.get('box_x', 0.05)
            rve = box_x / scale * 1e6 if scale > 0 else 50

            mf = ip.get('mass_fractions', [0, 0, 0])
            mf_se = mf[2] if len(mf) >= 3 else 0

            if d_se <= 0:
                gb_d = m.get('gb_density_mean', 0)
                if gb_d > 0:
                    d_se = 1.0 / gb_d
            if am_pct <= 0:
                phi_am = m.get('phi_am', 0)
                phi_se = m.get('phi_se', 0)
                if phi_am > 0 and phi_se > 0:
                    am_pct = phi_am / (phi_am + phi_se) * 100

            if d_se <= 0 or am_pct <= 0:
                continue

            tau = m.get('tortuosity_mean', m.get('tortuosity_recommended', 0))
            sigma_ion = m.get('sigma_full_mScm', 0)

            # Preprocessing
            if tau > 10 or m.get('porosity', 0) > 30:
                continue
            if 0 < sigma_ion < 0.01:
                continue

            rows.append({
                'name': case_dir.name,
                'd_am_p': d_am_p, 'd_am_s': d_am_s, 'd_se': d_se,
                'am_pct': am_pct, 'se_pct': 100 - am_pct,
                'ps_frac': ps_frac, 'rve': rve, 'mf_se': mf_se,
                'porosity': m.get('porosity', 0),
                'phi_se': m.get('phi_se', 0),
                'phi_am': m.get('phi_am', 0),
                'thickness': m.get('thickness_um', 0),
                'tau': tau,
                'cn': m.get('se_se_cn', 0),
                'gb_d': m.get('gb_density_mean', 0),
                'g_path': m.get('path_conductance_mean', 0),
                'hop_area': m.get('path_hop_area_mean', 0),
                'f_perc': m.get('percolation_pct', 0),
                'sigma_ion': sigma_ion,
                'sigma_el': m.get('electronic_sigma_full_mScm', 0),
                'sigma_th': m.get('thermal_sigma_full_mScm', 0),
                'sigma_ratio': m.get('sigma_ratio', 0),
            })

    # Deduplicate
    seen = set()
    unique = []
    for r in rows:
        key = f"{r['phi_se']:.4f}_{r.get('thickness',0):.1f}_{r['tau']:.3f}"
        if key not in seen:
            seen.add(key)
            unique.append(r)
    return unique


def run_symbolic_regression(rows, target_name='sigma_ion'):
    """Run PySR symbolic regression."""
    try:
        from pysr import PySRRegressor
    except ImportError:
        print("PySR not installed. Trying gplearn fallback...")
        return run_gplearn(rows, target_name)

    # Features: DEM inputs
    feature_names = ['d_se', 'am_pct', 'ps_frac', 'd_am_p', 'd_am_s', 'rve']
    X = np.array([[r[f] for f in feature_names] for r in rows])
    y = np.array([r[target_name] for r in rows])

    # Filter valid
    valid = y > 0
    X, y = X[valid], y[valid]
    print(f"\n  Target: {target_name}, n={len(y)}")
    print(f"  Range: {y.min():.4f} ~ {y.max():.4f}")

    model = PySRRegressor(
        niterations=100,
        binary_operators=["+", "-", "*", "/", "^"],
        unary_operators=["sqrt", "log", "exp"],
        populations=30,
        population_size=50,
        maxsize=25,
        parsimony=0.01,  # prefer simpler formulas
        constraints={"^": (-1, 1)},  # limit exponent complexity
        extra_sympy_mappings={"inv": lambda x: 1/x},
        loss="loss(prediction, target) = (log(prediction) - log(target))^2",  # log-space loss
        warm_start=False,
        verbosity=1,
    )

    print(f"\n  Running PySR (genetic programming)...")
    print(f"  This may take 5-30 minutes...\n")
    model.fit(X, y, variable_names=feature_names)

    print(f"\n{'='*60}")
    print(f"DISCOVERED FORMULAS for {target_name}")
    print(f"{'='*60}")
    print(model)

    # Best formula
    best = model.get_best()
    print(f"\n  BEST: {best['equation']}")
    print(f"  R² = {1 - best['loss']:.4f}")

    return model


def run_gplearn(rows, target_name='sigma_ion'):
    """Fallback: gplearn symbolic regression."""
    try:
        from gplearn.genetic import SymbolicRegressor
    except ImportError:
        print("Neither PySR nor gplearn installed. Using manual search...")
        return run_manual_search(rows, target_name)

    feature_names = ['d_se', 'am_pct', 'ps_frac', 'd_am_p', 'd_am_s', 'rve']
    X = np.array([[r[f] for f in feature_names] for r in rows])
    y = np.array([r[target_name] for r in rows])

    valid = y > 0
    X, y = X[valid], y[valid]

    # Use log(y) for better fitting
    y_log = np.log(y)

    print(f"\n  Target: {target_name} (log space), n={len(y)}")

    sr = SymbolicRegressor(
        population_size=5000,
        generations=50,
        tournament_size=20,
        stopping_criteria=0.01,
        p_crossover=0.7,
        p_subtree_mutation=0.1,
        p_hoist_mutation=0.05,
        p_point_mutation=0.1,
        max_samples=1.0,
        parsimony_coefficient=0.01,
        function_set=['add', 'sub', 'mul', 'div', 'sqrt', 'log'],
        verbose=1,
        random_state=42,
        n_jobs=-1,
    )

    print(f"  Running gplearn (genetic programming)...")
    sr.fit(X, y_log)

    # Evaluate
    y_pred_log = sr.predict(X)
    y_pred = np.exp(y_pred_log)
    ss_res = np.sum((np.log(y) - y_pred_log) ** 2)
    ss_tot = np.sum((np.log(y) - np.mean(np.log(y))) ** 2)
    r2 = 1 - ss_res / ss_tot

    print(f"\n{'='*60}")
    print(f"DISCOVERED FORMULA for {target_name}")
    print(f"{'='*60}")
    print(f"  log({target_name}) = {sr._program}")
    print(f"  R² (log space) = {r2:.4f}")
    print(f"\n  Variables: {feature_names}")

    return sr


def run_manual_search(rows, target_name='sigma_ion'):
    """Manual exhaustive formula search (no external ML library needed)."""
    print(f"\n{'='*60}")
    print(f"EXHAUSTIVE FORMULA SEARCH for {target_name}")
    print(f"(no PySR/gplearn — using brute-force power law search)")
    print(f"{'='*60}")

    feature_names = ['d_se', 'am_pct', 'ps_frac', 'd_am_p', 'd_am_s', 'rve']
    X = np.array([[r[f] for f in feature_names] for r in rows])
    y = np.array([r[target_name] for r in rows])

    valid = y > 0
    for i in range(X.shape[1]):
        valid &= X[:, i] > 0
    X, y = X[valid], y[valid]
    n = len(y)
    print(f"  n={n}, features: {feature_names}")

    if n < 5:
        print("  Not enough data!")
        return None

    log_y = np.log(y)
    ss_tot = np.sum((log_y - np.mean(log_y)) ** 2)
    log_X = np.log(X)

    # ── 1. Single variable power law ──
    print(f"\n  --- Single variable: {target_name} = C × x^a ---")
    single_results = []
    for i, f in enumerate(feature_names):
        A = np.column_stack([log_X[:, i], np.ones(n)])
        b, _, _, _ = np.linalg.lstsq(A, log_y, rcond=None)
        pred = A @ b
        r2 = 1 - np.sum((log_y - pred) ** 2) / ss_tot
        single_results.append((f, b[0], np.exp(b[1]), r2))
    single_results.sort(key=lambda x: -x[3])
    for f, exp, C, r2 in single_results:
        print(f"    {target_name} = {C:.4f} × {f}^{exp:.3f}  R²={r2:.4f}")

    # ── 2. Two variable ──
    print(f"\n  --- Two variables: {target_name} = C × x1^a × x2^b ---")
    two_results = []
    for i in range(len(feature_names)):
        for j in range(i + 1, len(feature_names)):
            A = np.column_stack([log_X[:, i], log_X[:, j], np.ones(n)])
            b, _, _, _ = np.linalg.lstsq(A, log_y, rcond=None)
            pred = A @ b
            r2 = 1 - np.sum((log_y - pred) ** 2) / ss_tot
            two_results.append((feature_names[i], feature_names[j], b[0], b[1], np.exp(b[2]), r2))
    two_results.sort(key=lambda x: -x[5])
    for f1, f2, e1, e2, C, r2 in two_results[:10]:
        print(f"    {C:.4f} × {f1}^{e1:.2f} × {f2}^{e2:.2f}  R²={r2:.4f}")

    # ── 3. Three variable ──
    print(f"\n  --- Three variables: {target_name} = C × x1^a × x2^b × x3^c ---")
    three_results = []
    for i in range(len(feature_names)):
        for j in range(i + 1, len(feature_names)):
            for k in range(j + 1, len(feature_names)):
                A = np.column_stack([log_X[:, i], log_X[:, j], log_X[:, k], np.ones(n)])
                b, _, _, _ = np.linalg.lstsq(A, log_y, rcond=None)
                pred = A @ b
                r2 = 1 - np.sum((log_y - pred) ** 2) / ss_tot
                three_results.append((feature_names[i], feature_names[j], feature_names[k],
                                     b[0], b[1], b[2], np.exp(b[3]), r2))
    three_results.sort(key=lambda x: -x[7])
    for f1, f2, f3, e1, e2, e3, C, r2 in three_results[:10]:
        print(f"    {C:.4f} × {f1}^{e1:.2f} × {f2}^{e2:.2f} × {f3}^{e3:.2f}  R²={r2:.4f}")

    # ── 4. All variables ──
    print(f"\n  --- All variables ---")
    A = np.column_stack([log_X, np.ones(n)])
    b, _, _, _ = np.linalg.lstsq(A, log_y, rcond=None)
    pred = A @ b
    r2 = 1 - np.sum((log_y - pred) ** 2) / ss_tot
    parts = ' × '.join([f"{feature_names[i]}^{b[i]:.2f}" for i in range(len(feature_names))])
    print(f"    {np.exp(b[-1]):.4f} × {parts}  R²={r2:.4f}")

    # ── 5. Derived combinations ──
    print(f"\n  --- Smart combinations ---")
    combo_results = []

    # d_SE / d_AM ratio
    d_ratio = X[:, 0] / np.maximum(X[:, 3], X[:, 4])  # d_SE / max(d_AM_P, d_AM_S)
    valid_c = d_ratio > 0
    if valid_c.sum() > 5:
        for i, f in enumerate(feature_names):
            if f in ('d_am_p', 'd_am_s', 'd_se'):
                continue
            A = np.column_stack([np.log(d_ratio[valid_c]), log_X[valid_c, i], np.ones(valid_c.sum())])
            b_c, _, _, _ = np.linalg.lstsq(A, log_y[valid_c], rcond=None)
            pred_c = A @ b_c
            r2_c = 1 - np.sum((log_y[valid_c] - pred_c) ** 2) / np.sum((log_y[valid_c] - np.mean(log_y[valid_c])) ** 2)
            combo_results.append((f"(d_SE/d_AM)^{b_c[0]:.2f} × {f}^{b_c[1]:.2f}", r2_c))

    # SE volume estimation: mf_se / (mf_se + mf_am × ρ_SE/ρ_AM)
    # AM:SE mass → volume fraction proxy
    se_vol_proxy = (100 - X[:, 1]) / 100  # SE mass fraction proxy
    for i, f in enumerate(feature_names):
        if f == 'am_pct':
            continue
        A = np.column_stack([np.log(se_vol_proxy), log_X[:, i], np.ones(n)])
        b_c, _, _, _ = np.linalg.lstsq(A, log_y, rcond=None)
        pred_c = A @ b_c
        r2_c = 1 - np.sum((log_y - pred_c) ** 2) / ss_tot
        combo_results.append((f"SE_frac^{b_c[0]:.2f} × {f}^{b_c[1]:.2f}", r2_c))

    combo_results.sort(key=lambda x: -x[1])
    for formula, r2 in combo_results[:10]:
        print(f"    {formula}  R²={r2:.4f}")

    return three_results[0] if three_results else None


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Symbolic Regression for DEM')
    parser.add_argument('--target', default='all', help='Target variable (default: all key targets)')
    args = parser.parse_args()

    print("=" * 60)
    print("SYMBOLIC REGRESSION: DEM Input → Output Formula Discovery")
    print("=" * 60)

    rows = load_data()
    print(f"Loaded {len(rows)} cases (after preprocessing)")

    targets = ['sigma_ion', 'sigma_el', 'sigma_th', 'porosity', 'phi_se', 'tau', 'cn']
    if args.target != 'all':
        targets = [args.target]

    for target in targets:
        vals = [r[target] for r in rows if r[target] > 0]
        if len(vals) < 10:
            print(f"\n  {target}: not enough data ({len(vals)} cases)")
            continue
        run_manual_search(rows, target)  # Always do manual first (no deps)

        # Try PySR/gplearn if available
        try:
            from pysr import PySRRegressor
            print(f"\n  [PySR available — running advanced search...]")
            run_symbolic_regression(rows, target)
        except ImportError:
            try:
                from gplearn.genetic import SymbolicRegressor
                print(f"\n  [gplearn available — running genetic programming...]")
                run_gplearn(rows, target)
            except ImportError:
                pass


if __name__ == '__main__':
    main()
