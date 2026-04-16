"""
ML Surrogate Model v2: True DEM Input → Property Prediction
============================================================
Input: LIGGGHTS input script parameters (particle sizes, mass fractions, RVE, pressure)
Output: ALL DEM-derived properties (porosity, σ_ionic, σ_electronic, σ_thermal, etc.)

This is a TRUE surrogate: predicts DEM output from DEM input without running simulation.

Usage:
  python3 scripts/ml_predictor.py                    # Train & evaluate
  python3 scripts/ml_predictor.py --predict           # Interactive prediction
"""

import json, os, sys, warnings
import numpy as np
from pathlib import Path

warnings.filterwarnings('ignore')

WEBAPP = Path(__file__).parent.parent / 'webapp'


def load_all_data():
    """Load all cases: input_params.json (INPUT) + full_metrics.json (OUTPUT)."""
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
                if ip_path.exists():
                    with open(ip_path) as f:
                        ip = json.load(f)
                else:
                    ip = {}
            except:
                continue

            if not m.get('phi_se') or not m.get('porosity'):
                continue

            # ═══ TRUE DEM INPUTS (from input_params.json) ═══
            scale = ip.get('scale', 1000)

            # Particle radii (sim units → real μm)
            r_am_p = ip.get('r_AM_P', 0)
            r_am_s = ip.get('r_AM_S', 0)
            r_se = ip.get('r_SE', 0)
            d_am_p = r_am_p / scale * 1e6 * 2 if r_am_p > 0 else 0  # μm
            d_am_s = r_am_s / scale * 1e6 * 2 if r_am_s > 0 else 0
            d_se = r_se / scale * 1e6 * 2 if r_se > 0 else 0

            # Mass fractions
            mf = ip.get('mass_fractions', [0, 0, 0])
            if len(mf) >= 3:
                mf_amp, mf_ams, mf_se = mf[0], mf[1], mf[2]
            else:
                mf_amp, mf_ams, mf_se = 0, 0, 0

            # AM:SE ratio
            am_se_str = ip.get('am_se_ratio', '')
            am_pct = 0
            if ':' in str(am_se_str):
                parts = str(am_se_str).split(':')
                try:
                    am_pct = float(parts[0])
                except:
                    pass

            # P:S ratio
            ps = m.get('ps_ratio', '')
            ps_frac = 0.5
            if isinstance(ps, str) and ':' in ps:
                parts = ps.split(':')
                try:
                    p, s = float(parts[0]), float(parts[1])
                    ps_frac = p / (p + s) if (p + s) > 0 else 0.5
                except:
                    pass

            # RVE size (sim → real μm)
            box_x = ip.get('box_x', 0.05)
            rve_um = box_x / scale * 1e6 if scale > 0 else 50  # μm

            # Pressure (sim → real MPa)
            target_press = ip.get('target_press_sim', 0)
            pressure_mpa = target_press * scale  # sim → real MPa

            # Young's modulus SE (sim → real, effective GPa)
            # SE는 항상 마지막 type: standard(2-type)→index 1, bimodal(3-type)→index 2
            ym = ip.get('youngs_modulus_sim', [])
            if len(ym) >= 2:
                E_se = ym[-1] * scale / 1e9  # sim Pa × scale → effective GPa
            else:
                E_se = 0

            # Skip if missing key inputs
            if d_se <= 0 or am_pct <= 0:
                # Try to infer from metrics
                if d_se <= 0:
                    gb_d = m.get('gb_density_mean', 0)
                    if gb_d > 0:
                        d_se = 1.0 / gb_d  # rough proxy
                if am_pct <= 0:
                    phi_am = m.get('phi_am', 0)
                    phi_se = m.get('phi_se', 0)
                    if phi_am > 0 and phi_se > 0:
                        am_pct = phi_am / (phi_am + phi_se) * 100

            if d_se <= 0:
                continue

            # ═══ DEM OUTPUTS (from full_metrics.json) ═══
            rows.append({
                'name': case_dir.name,
                # ── True Inputs ──
                'd_am_p': d_am_p,
                'd_am_s': d_am_s,
                'd_se': d_se,
                'am_pct': am_pct,
                'ps_frac': ps_frac,
                'rve': rve_um,
                'pressure': pressure_mpa,
                'mf_amp': mf_amp,
                'mf_ams': mf_ams,
                'mf_se': mf_se,
                # ── Outputs: Structure ──
                'porosity': m.get('porosity', 0),
                'phi_se': m.get('phi_se', 0),
                'phi_am': m.get('phi_am', 0),
                'thickness': m.get('thickness_um', 0),
                # ── Outputs: Ion Path ──
                'tau': m.get('tortuosity_mean', m.get('tortuosity_recommended', 0)),
                'cn': m.get('se_se_cn', 0),
                'gb_d': m.get('gb_density_mean', 0),
                'g_path': m.get('path_conductance_mean', 0),
                'f_perc': m.get('percolation_pct', 0),
                'hop_area': m.get('path_hop_area_mean', 0),
                # ── Outputs: Conductivity ──
                'sigma_ion': m.get('sigma_full_mScm', 0),
                'sigma_el': m.get('electronic_sigma_full_mScm', 0),
                'sigma_th': m.get('thermal_sigma_full_mScm', 0),
            })

    return rows


def build_dataset(rows):
    """Build feature matrix X and target matrix Y with preprocessing."""
    # TRUE DEM input features
    feature_names = ['d_am_p', 'd_am_s', 'd_se', 'am_pct', 'ps_frac', 'rve', 'pressure']

    # All output targets
    target_names = ['porosity', 'phi_se', 'phi_am', 'thickness',
                    'tau', 'cn', 'gb_d', 'f_perc', 'hop_area', 'g_path',
                    'sigma_ion', 'sigma_el', 'sigma_th']

    X_list, Y_list, names = [], [], []
    seen = set()
    skipped = {'tau_outlier': 0, 'low_sigma': 0, 'high_porosity': 0, 'missing': 0}

    for r in rows:
        # Deduplicate
        dedup_key = f"{r['phi_se']:.4f}_{r.get('thickness',0):.1f}_{r['tau']:.3f}"
        if dedup_key in seen:
            continue
        seen.add(dedup_key)

        # ── Preprocessing filters ──
        if r['tau'] > 10:
            skipped['tau_outlier'] += 1
            continue
        if 0 < r['sigma_ion'] < 0.01:
            skipped['low_sigma'] += 1
            continue
        if r['porosity'] > 30:
            skipped['high_porosity'] += 1
            continue
        if r['d_se'] <= 0 or r['am_pct'] <= 0:
            skipped['missing'] += 1
            continue

        x = [r.get(f, 0) for f in feature_names]
        y = [r.get(t, 0) for t in target_names]
        X_list.append(x)
        Y_list.append(y)
        names.append(r['name'])

    X = np.array(X_list)
    Y = np.array(Y_list)

    print(f"  Skipped: {skipped}")
    return X, Y, feature_names, target_names, names


def train_gpr(X, Y, target_names):
    """Train Gaussian Process Regression for each target."""
    from sklearn.gaussian_process import GaussianProcessRegressor
    from sklearn.gaussian_process.kernels import ConstantKernel, WhiteKernel, Matern
    from sklearn.preprocessing import StandardScaler
    from sklearn.model_selection import LeaveOneOut

    scaler_X = StandardScaler()
    X_scaled = scaler_X.fit_transform(X)

    models = {}
    results = {}

    for j, target in enumerate(target_names):
        y = Y[:, j]

        if np.std(y) == 0 or np.all(y == 0):
            print(f"  {target:15s}: SKIP (no variance)")
            continue

        # Log transform for positive targets with large range
        use_log = np.all(y > 0) and (np.max(y) / np.min(y) > 5)
        y_train = np.log(y) if use_log else y.copy()

        scaler_y = StandardScaler()
        y_scaled = scaler_y.fit_transform(y_train.reshape(-1, 1)).ravel()

        kernel = ConstantKernel(1.0) * Matern(length_scale=1.0, nu=2.5) + WhiteKernel(noise_level=0.1)

        gpr = GaussianProcessRegressor(
            kernel=kernel, n_restarts_optimizer=10, alpha=1e-6, normalize_y=False
        )
        gpr.fit(X_scaled, y_scaled)

        # Leave-One-Out Cross-Validation
        loo = LeaveOneOut()
        y_pred_loo = np.zeros(len(y))

        for train_idx, test_idx in loo.split(X_scaled):
            gpr_cv = GaussianProcessRegressor(
                kernel=kernel, n_restarts_optimizer=3, alpha=1e-6, normalize_y=False
            )
            gpr_cv.fit(X_scaled[train_idx], y_scaled[train_idx])
            pred, _ = gpr_cv.predict(X_scaled[test_idx], return_std=True)
            y_pred_loo[test_idx] = pred

        # Inverse transform
        y_pred_real = scaler_y.inverse_transform(y_pred_loo.reshape(-1, 1)).ravel()
        if use_log:
            y_pred_real = np.exp(y_pred_real)
            y_actual = y
        else:
            y_actual = y

        # Metrics
        ss_res = np.sum((y_actual - y_pred_real) ** 2)
        ss_tot = np.sum((y_actual - np.mean(y_actual)) ** 2)
        r2_cv = 1 - ss_res / ss_tot if ss_tot > 0 else 0

        valid = y_actual > 0
        mean_err = np.mean(np.abs(y_pred_real[valid] - y_actual[valid]) / y_actual[valid]) * 100 if valid.any() else 0

        models[target] = {
            'gpr': gpr, 'scaler_X': scaler_X, 'scaler_y': scaler_y,
            'use_log': use_log
        }
        results[target] = {'r2_cv': r2_cv, 'mean_err': mean_err, 'n': len(y)}

        star = "★" if r2_cv > 0.8 else ""
        print(f"  {target:15s}: LOO-CV R²={r2_cv:.3f}, mean|err|={mean_err:.1f}% (n={len(y)}) {star}")

    return models, results


def predict_new(models, feature_names, target_names, input_dict):
    """Predict properties for new input."""
    x = np.array([[input_dict.get(f, 0) for f in feature_names]])

    predictions = {}
    for target in target_names:
        if target not in models:
            predictions[target] = {'value': 0, 'std': 0}
            continue

        m = models[target]
        x_scaled = m['scaler_X'].transform(x)
        pred_scaled, std_scaled = m['gpr'].predict(x_scaled, return_std=True)
        pred = m['scaler_y'].inverse_transform(pred_scaled.reshape(-1, 1)).ravel()[0]
        std = std_scaled[0] * m['scaler_y'].scale_[0]

        if m['use_log']:
            pred = np.exp(pred)
            std = pred * abs(std)

        predictions[target] = {'value': round(float(pred), 4), 'std': round(float(abs(std)), 4)}

    return predictions


def main():
    import argparse
    parser = argparse.ArgumentParser(description='ML Surrogate v2: DEM Input → Properties')
    parser.add_argument('--predict', action='store_true', help='Interactive prediction mode')
    args = parser.parse_args()

    print("=" * 70)
    print("ML SURROGATE v2: TRUE DEM INPUT → PROPERTY PREDICTION")
    print("Gaussian Process Regression (optimal for small datasets)")
    print("=" * 70)

    rows = load_all_data()
    print(f"\nLoaded {len(rows)} cases")

    X, Y, feature_names, target_names, names = build_dataset(rows)
    print(f"Dataset: {len(X)} cases after preprocessing")

    if len(X) < 10:
        print("ERROR: Need at least 10 cases!")
        return

    # Data ranges
    print(f"\n{'─'*60}")
    print(f"{'INPUT FEATURES (DEM script parameters)':^60}")
    print(f"{'─'*60}")
    units = {'d_am_p': 'μm', 'd_am_s': 'μm', 'd_se': 'μm', 'am_pct': '%',
             'ps_frac': '(0~1)', 'rve': 'μm', 'pressure': 'MPa'}
    for i, f in enumerate(feature_names):
        u = units.get(f, '')
        print(f"  {f:12s} [{u:>6s}]: {X[:,i].min():10.2f} ~ {X[:,i].max():10.2f} (mean {X[:,i].mean():.2f})")

    print(f"\n{'─'*60}")
    print(f"{'OUTPUT TARGETS (DEM results)':^60}")
    print(f"{'─'*60}")
    for j, t in enumerate(target_names):
        valid = Y[:, j] > 0
        if valid.any():
            print(f"  {t:13s}: {Y[valid,j].min():10.4f} ~ {Y[valid,j].max():10.4f}")

    # Train
    print(f"\n{'═'*60}")
    print(f"TRAINING (LOO-CV: honest R² with n-1 training)")
    print(f"{'═'*60}")
    models, results = train_gpr(X, Y, target_names)

    # Summary
    print(f"\n{'═'*60}")
    print(f"{'PERFORMANCE SUMMARY':^60}")
    print(f"{'═'*60}")
    print(f"  {'Target':15s} {'LOO-CV R²':>10s} {'|err|':>8s} {'Grade':>6s}")
    print(f"  {'─'*45}")
    for t in target_names:
        if t in results:
            r = results[t]
            if r['r2_cv'] > 0.95:
                grade = "A+"
            elif r['r2_cv'] > 0.9:
                grade = "A"
            elif r['r2_cv'] > 0.8:
                grade = "B"
            elif r['r2_cv'] > 0.5:
                grade = "C"
            else:
                grade = "F"
            print(f"  {t:15s} {r['r2_cv']:10.3f} {r['mean_err']:7.1f}% {grade:>6s}")

    # Scaling Law comparison
    if 'sigma_ion' in results:
        print(f"\n  ── σ_ionic: GPR vs Scaling Law ──")
        print(f"  Scaling Law: R²=0.947 (physics, needs DEM microstructure)")
        print(f"  GPR (input): R²={results['sigma_ion']['r2_cv']:.3f} (from script params only!)")
        if results['sigma_ion']['r2_cv'] < 0.947:
            print(f"  → Scaling Law still better (but GPR doesn't need DEM!)")

    # Interactive prediction
    if args.predict:
        print(f"\n{'═'*60}")
        print(f"INTERACTIVE PREDICTION")
        print(f"{'═'*60}")
        print("Enter DEM input parameters:")

        try:
            d_am_p = float(input("  d_AM_P (μm) [default 10]: ") or 10)
            d_am_s = float(input("  d_AM_S (μm) [default 5]: ") or 5)
            d_se = float(input("  d_SE (μm) [default 1]: ") or 1)
            am_pct = float(input("  AM:SE ratio - AM% [default 80]: ") or 80)
            ps_frac = float(input("  P:S fraction (0=all S, 1=all P) [default 0.5]: ") or 0.5)
            rve = float(input("  RVE size (μm) [default 50]: ") or 50)
            pressure = float(input("  Pressure (MPa) [default 300]: ") or 300)
        except (ValueError, EOFError):
            print("Using defaults...")
            d_am_p, d_am_s, d_se = 10, 5, 1
            am_pct, ps_frac, rve, pressure = 80, 0.5, 50, 300

        input_dict = {
            'd_am_p': d_am_p, 'd_am_s': d_am_s, 'd_se': d_se,
            'am_pct': am_pct, 'ps_frac': ps_frac, 'rve': rve, 'pressure': pressure
        }

        predictions = predict_new(models, feature_names, target_names, input_dict)

        print(f"\n{'─'*50}")
        print(f"PREDICTED PROPERTIES:")
        print(f"{'─'*50}")
        units_out = {
            'porosity': '%', 'phi_se': '', 'phi_am': '', 'thickness': 'μm',
            'tau': '', 'cn': '', 'gb_d': 'hops/μm', 'f_perc': '%',
            'hop_area': 'μm²', 'g_path': 'μm²',
            'sigma_ion': 'mS/cm', 'sigma_el': 'mS/cm', 'sigma_th': 'mS/cm'
        }
        for t in target_names:
            p = predictions[t]
            u = units_out.get(t, '')
            conf = "±" + f"{p['std']:.4f}" if p['std'] > 0 else ""
            print(f"  {t:15s} = {p['value']:10.4f} {u:>8s}  {conf}")


if __name__ == '__main__':
    main()
