"""
ML Surrogate Model for DEM Composite Cathode Properties
========================================================
Input: electrode design parameters (AM:SE ratio, P:S, SE size, thickness, pressure)
Output: predicted properties (σ_ionic, σ_electronic, σ_thermal, porosity, τ, CN, etc.)

Uses Gaussian Process Regression (GPR) — optimal for small datasets (n < 100).
Provides uncertainty estimates for each prediction.

Usage:
  python3 scripts/ml_predictor.py                    # Train & evaluate
  python3 scripts/ml_predictor.py --predict 80 5 0.5 150 300  # Predict for AM:SE=80:20, P:S=5:5, SE=0.5μm, T=150μm, P=300MPa
"""

import json, os, sys, warnings
import numpy as np
from pathlib import Path

warnings.filterwarnings('ignore')

WEBAPP = Path(__file__).parent.parent / 'webapp'


def load_all_data():
    """Load all cases from results + archive."""
    rows = []
    for base in [WEBAPP / 'results', WEBAPP / 'archive']:
        if not base.is_dir():
            continue
        for met_path in base.rglob('full_metrics.json'):
            try:
                with open(met_path) as f:
                    m = json.load(f)
            except:
                continue

            # Skip if no basic data
            if not m.get('phi_se') or not m.get('porosity'):
                continue

            # Extract input features
            # AM:SE ratio → from phi_se and phi_am
            phi_se = m.get('phi_se', 0)
            phi_am = m.get('phi_am', 0)
            if phi_se <= 0:
                continue

            am_se_mass = 0
            ps = m.get('ps_ratio', '')

            # P:S fraction (0=all S, 1=all P)
            ps_frac = 0.5  # default
            if isinstance(ps, str):
                if ':' in ps:
                    parts = ps.split(':')
                    try:
                        p, s = float(parts[0]), float(parts[1])
                        ps_frac = p / (p + s) if (p + s) > 0 else 0.5
                    except:
                        pass
                elif ps in ('P only', '10:0'):
                    ps_frac = 1.0
                elif ps in ('S only', '0:10'):
                    ps_frac = 0.0

            # SE size from input_params or GB_d proxy
            d_se = 0
            ip_path = met_path.parent / 'input_params.json'
            if ip_path.exists():
                try:
                    with open(ip_path) as f:
                        ip = json.load(f)
                    scale = ip.get('scale', 1000)
                    for k in ['r_SE']:
                        v = ip.get(k, 0)
                        if v > 0:
                            # v is in sim units (m, already scaled by DEM)
                            # Real radius (μm) = v / scale * 1e6, diameter = *2
                            d_se = v / scale * 1e6 * 2
                            break
                except:
                    pass
            if d_se <= 0:
                # Proxy from GB_d
                gb_d = m.get('gb_density_mean', 0)
                if gb_d > 0:
                    d_se = 1.0 / gb_d  # rough proxy

            thickness = m.get('thickness_um', 0)
            porosity = m.get('porosity', 0)

            # Output metrics
            sigma_ion = m.get('sigma_full_mScm', 0)
            sigma_el = m.get('electronic_sigma_full_mScm', 0)
            sigma_th = m.get('thermal_sigma_full_mScm', 0)
            tau = m.get('tortuosity_mean', m.get('tortuosity_recommended', 0))
            cn = m.get('se_se_cn', 0)
            gb_d = m.get('gb_density_mean', 0)
            g_path = m.get('path_conductance_mean', 0)
            f_perc = m.get('percolation_pct', 0)
            hop_area = m.get('path_hop_area_mean', 0)

            rows.append({
                'name': met_path.parent.name,
                # Input features
                'phi_se': phi_se,
                'phi_am': phi_am,
                'ps_frac': ps_frac,
                'd_se': d_se,
                'thickness': thickness,
                'porosity': porosity,
                # Output targets (from network solver)
                'sigma_ion': sigma_ion,
                'sigma_el': sigma_el,
                'sigma_th': sigma_th,
                # Intermediate outputs (from DEM analysis)
                'tau': tau,
                'cn': cn,
                'gb_d': gb_d,
                'g_path': g_path,
                'f_perc': f_perc,
                'hop_area': hop_area,
            })

    return rows


def build_dataset(rows):
    """Build feature matrix X and target matrix Y."""
    # Input features
    feature_names = ['phi_se', 'phi_am', 'ps_frac', 'd_se', 'thickness']
    # Output targets
    target_names = ['porosity', 'tau', 'cn', 'gb_d', 'f_perc', 'hop_area', 'g_path',
                    'sigma_ion', 'sigma_el', 'sigma_th']

    X_list, Y_list, names = [], [], []
    seen = set()
    skipped = {'tau_outlier': 0, 'low_sigma': 0, 'high_porosity': 0}
    for r in rows:
        # Deduplicate by rounding key metrics
        dedup_key = f"{r['phi_se']:.4f}_{r['thickness']:.1f}_{r['tau']:.3f}"
        if dedup_key in seen:
            continue
        seen.add(dedup_key)

        # ── Preprocessing filters ──
        # 1. tau outlier (> 10 = non-physical, DEM artifact)
        if r['tau'] > 10:
            skipped['tau_outlier'] += 1
            continue
        # 2. σ_ionic < 0.01 mS/cm (near percolation threshold)
        if 0 < r['sigma_ion'] < 0.01:
            skipped['low_sigma'] += 1
            continue
        # 3. porosity > 30% (electrode not properly formed)
        if r['porosity'] > 30:
            skipped['high_porosity'] += 1
            continue

        # Require valid inputs
        if r['phi_se'] <= 0 or r['thickness'] <= 0 or r['d_se'] <= 0:
            continue
        # Require at least some outputs
        if r['tau'] <= 0 and r['sigma_ion'] <= 0:
            continue

        x = [r[f] for f in feature_names]
        y = [r[t] for t in target_names]
        X_list.append(x)
        Y_list.append(y)
        names.append(r['name'])

    X = np.array(X_list)
    Y = np.array(Y_list)
    return X, Y, feature_names, target_names, names


def train_gpr(X, Y, target_names):
    """Train Gaussian Process Regression for each target."""
    from sklearn.gaussian_process import GaussianProcessRegressor
    from sklearn.gaussian_process.kernels import RBF, ConstantKernel, WhiteKernel, Matern
    from sklearn.preprocessing import StandardScaler
    from sklearn.model_selection import LeaveOneOut

    # Normalize inputs
    scaler_X = StandardScaler()
    X_scaled = scaler_X.fit_transform(X)

    models = {}
    results = {}

    for j, target in enumerate(target_names):
        y = Y[:, j]

        # Skip if all zeros or constant
        if np.std(y) == 0 or np.all(y == 0):
            print(f"  {target:15s}: SKIP (no variance)")
            continue

        # Use log transform for positive targets with large range
        use_log = np.all(y > 0) and (np.max(y) / np.min(y) > 5)
        if use_log:
            y_train = np.log(y)
        else:
            y_train = y.copy()

        # Normalize target
        scaler_y = StandardScaler()
        y_scaled = scaler_y.fit_transform(y_train.reshape(-1, 1)).ravel()

        # GPR kernel: Matern(ν=2.5) + WhiteKernel (noise)
        kernel = ConstantKernel(1.0) * Matern(length_scale=1.0, nu=2.5) + WhiteKernel(noise_level=0.1)

        gpr = GaussianProcessRegressor(
            kernel=kernel, n_restarts_optimizer=10, alpha=1e-6, normalize_y=False
        )
        gpr.fit(X_scaled, y_scaled)

        # Leave-One-Out Cross-Validation
        loo = LeaveOneOut()
        y_pred_loo = np.zeros(len(y))
        y_std_loo = np.zeros(len(y))

        for train_idx, test_idx in loo.split(X_scaled):
            gpr_cv = GaussianProcessRegressor(
                kernel=kernel, n_restarts_optimizer=5, alpha=1e-6, normalize_y=False
            )
            gpr_cv.fit(X_scaled[train_idx], y_scaled[train_idx])
            pred, std = gpr_cv.predict(X_scaled[test_idx], return_std=True)
            y_pred_loo[test_idx] = pred
            y_std_loo[test_idx] = std

        # Inverse transform predictions
        y_pred_real = scaler_y.inverse_transform(y_pred_loo.reshape(-1, 1)).ravel()
        if use_log:
            y_pred_real = np.exp(y_pred_real)
            y_actual = y
        else:
            y_actual = y

        # R² (LOO-CV)
        ss_res = np.sum((y_actual - y_pred_real) ** 2)
        ss_tot = np.sum((y_actual - np.mean(y_actual)) ** 2)
        r2_cv = 1 - ss_res / ss_tot if ss_tot > 0 else 0

        # Mean error
        valid = y_actual > 0
        if valid.any():
            mean_err = np.mean(np.abs(y_pred_real[valid] - y_actual[valid]) / y_actual[valid]) * 100
        else:
            mean_err = 0

        models[target] = {
            'gpr': gpr, 'scaler_X': scaler_X, 'scaler_y': scaler_y,
            'use_log': use_log, 'kernel': str(gpr.kernel_)
        }
        results[target] = {'r2_cv': r2_cv, 'mean_err': mean_err, 'n': len(y)}

        print(f"  {target:15s}: LOO-CV R²={r2_cv:.3f}, mean|err|={mean_err:.1f}% (n={len(y)})")

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
        std = std_scaled[0] * m['scaler_y'].scale_[0]  # approximate

        if m['use_log']:
            pred = np.exp(pred)
            std = pred * std  # approximate in original scale

        predictions[target] = {'value': round(pred, 4), 'std': round(abs(std), 4)}

    return predictions


def main():
    import argparse
    parser = argparse.ArgumentParser(description='ML Surrogate Model for DEM Cathode')
    parser.add_argument('--predict', nargs=5, type=float, metavar=('AM_PCT', 'P_RATIO', 'D_SE', 'T', 'PRESSURE'),
                        help='Predict: AM%% P-ratio(0-10) d_SE(μm) T(μm) Pressure(MPa)')
    args = parser.parse_args()

    print("=" * 70)
    print("ML SURROGATE MODEL FOR DEM COMPOSITE CATHODE")
    print("=" * 70)

    # Load data
    rows = load_all_data()
    print(f"\nLoaded {len(rows)} cases")

    X, Y, feature_names, target_names, names = build_dataset(rows)
    print(f"Valid dataset: {len(X)} cases (after preprocessing)")

    if len(X) < 10:
        print("ERROR: Need at least 10 cases for ML. Run more DEM simulations!")
        return

    # Show data range
    print(f"\n{'Feature':15s} {'min':>10s} {'max':>10s} {'mean':>10s}")
    print("-" * 50)
    for i, f in enumerate(feature_names):
        print(f"  {f:13s} {X[:,i].min():10.3f} {X[:,i].max():10.3f} {X[:,i].mean():10.3f}")

    print(f"\n{'Target':15s} {'min':>10s} {'max':>10s} {'mean':>10s}")
    print("-" * 50)
    for j, t in enumerate(target_names):
        valid = Y[:, j] > 0
        if valid.any():
            print(f"  {t:13s} {Y[valid,j].min():10.4f} {Y[valid,j].max():10.4f} {Y[valid,j].mean():10.4f}")

    # Train GPR
    print(f"\n--- Training Gaussian Process Regression (LOO-CV) ---")
    models, results = train_gpr(X, Y, target_names)

    # Summary
    print(f"\n{'=' * 50}")
    print(f"MODEL PERFORMANCE SUMMARY")
    print(f"{'=' * 50}")
    print(f"{'Target':15s} {'LOO-CV R²':>10s} {'Mean|err|':>10s}")
    print("-" * 40)
    for t in target_names:
        if t in results:
            r = results[t]
            star = " ★" if r['r2_cv'] > 0.8 else ""
            print(f"  {t:13s} {r['r2_cv']:10.3f} {r['mean_err']:9.1f}%{star}")

    # Predict new case
    if args.predict:
        am_pct, p_ratio, d_se, thickness, pressure = args.predict
        se_pct = 100 - am_pct
        phi_se_est = se_pct / 100 * 0.65  # rough estimate (porosity ~20%, SE density fraction)
        phi_am_est = am_pct / 100 * 0.65
        ps_frac = p_ratio / 10

        input_dict = {
            'phi_se': phi_se_est,
            'phi_am': phi_am_est,
            'ps_frac': ps_frac,
            'd_se': d_se,
            'thickness': thickness,
        }

        print(f"\n{'=' * 50}")
        print(f"PREDICTION: AM:SE={am_pct:.0f}:{se_pct:.0f}, P:S={p_ratio:.0f}:{10-p_ratio:.0f}, "
              f"d_SE={d_se}μm, T={thickness}μm")
        print(f"{'=' * 50}")

        predictions = predict_new(models, feature_names, target_names, input_dict)
        print(f"{'Property':20s} {'Predicted':>12s} {'±Uncertainty':>12s}")
        print("-" * 50)
        for t in target_names:
            p = predictions[t]
            print(f"  {t:18s} {p['value']:12.4f} ±{p['std']:10.4f}")

    # Comparison: GPR vs Scaling Law
    if 'sigma_ion' in results:
        print(f"\n--- GPR vs Scaling Law (ionic) ---")
        print(f"  Scaling Law: R²=0.947 (physics-based, 1 param)")
        print(f"  GPR LOO-CV:  R²={results['sigma_ion']['r2_cv']:.3f} (data-driven, kernel params)")
        if results['sigma_ion']['r2_cv'] > 0.947:
            print(f"  → GPR wins by {results['sigma_ion']['r2_cv'] - 0.947:.3f}")
        else:
            print(f"  → Scaling Law wins by {0.947 - results['sigma_ion']['r2_cv']:.3f}")
            print(f"  → Physics-based model is BETTER with this small dataset!")


if __name__ == '__main__':
    main()
