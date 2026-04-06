"""
2-Stage Hybrid Predictor: DEM Input → σ Prediction
===================================================
Stage 1: GPR (DEM input → microstructure)
Stage 2: Scaling Law (microstructure → conductivity)

Best of both worlds:
- GPR handles input→microstructure mapping (R²>0.95)
- Physics handles microstructure→σ (R²=0.947)

Usage:
  python3 scripts/hybrid_predictor.py                    # Train & evaluate
  python3 scripts/hybrid_predictor.py --predict           # Interactive
  python3 scripts/hybrid_predictor.py --sweep             # Parameter sweep
"""

import json, os, sys, warnings
import numpy as np
from pathlib import Path

warnings.filterwarnings('ignore')
WEBAPP = Path(__file__).parent.parent / 'webapp'

# ═══ Physical Constants ═══
SIGMA_GRAIN = 3.0       # mS/cm (SE grain interior)
SIGMA_AM = 50.0         # mS/cm (0.05 S/cm, NCM811)
K_AM = 4.0e-2           # W/(cm·K)
K_SE = 0.7e-2           # W/(cm·K)


def load_data():
    """Load all cases with DEM input + output."""
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

            if not m.get('phi_se'):
                continue

            scale = ip.get('scale', 1000)
            r_se = ip.get('r_SE', 0)
            d_se = r_se / scale * 1e6 * 2 if r_se > 0 and scale > 0 else 0
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

            # Infer missing
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

            # Loading (mAh/cm²) — infer from folder name or thickness
            loading = 0
            folder_path = str(met_path.parent)
            if '1mAh' in folder_path or 'thin' in folder_path.lower():
                loading = 1
            elif '8mAh' in folder_path:
                loading = 8
            elif '6mAh' in folder_path or 'real' in folder_path.lower() or 'Real' in folder_path:
                loading = 6
            elif 'particulate' in folder_path.lower():
                # Estimate from thickness: loading ≈ T(μm) × 0.05 (rough)
                t = m.get('thickness_um', 0)
                if t > 0:
                    loading = round(t * 0.05, 1)  # rough estimate
            if loading <= 0 and m.get('thickness_um', 0) > 0:
                loading = round(m['thickness_um'] * 0.05, 1)

            # d_am: use whichever is available
            d_am = max(d_am_p, d_am_s)

            rows.append({
                'name': case_dir.name,
                # Input
                'd_se': d_se, 'd_am': d_am, 'd_am_p': d_am_p, 'd_am_s': d_am_s,
                'am_pct': am_pct, 'ps_frac': ps_frac, 'rve': rve,
                'loading': loading,
                # Derived features
                'd_ratio': d_se / d_am if d_am > 0 else 0.1,  # SE/AM size ratio
                'am_loading': am_pct * loading / 100,  # AM mass per area proxy
                'se_density_proxy': (100 - am_pct) / max(d_se, 0.1),  # SE%/d_SE = packing density
                'layer_count': loading * 30 / max(d_se, 0.1),  # ~T/d_SE = number of SE layers
                # Microstructure (Stage 1 targets)
                'phi_se': m.get('phi_se', 0),
                'phi_am': m.get('phi_am', 0),
                'porosity': m.get('porosity', 0),
                'tau': tau,
                'cn': m.get('se_se_cn', 0),
                'gb_d': m.get('gb_density_mean', 0),
                'g_path': m.get('path_conductance_mean', 0),
                'hop_area': m.get('path_hop_area_mean', 0),
                'f_perc': m.get('percolation_pct', 0),
                'thickness': m.get('thickness_um', 0),
                'am_cn': m.get('am_am_cn', 0),
                # Conductivity (Stage 2 targets / ground truth)
                'sigma_ion': sigma_ion,
                'sigma_el': m.get('electronic_sigma_full_mScm', 0),
                'sigma_th': m.get('thermal_sigma_full_mScm', 0),
                # σ_brug = φ_SE × f_perc/100 / τ² (combines φ, f_perc, τ into one)
                'sigma_brug': m.get('sigma_ratio', 0),  # = φ×f_perc/τ² (ratio to σ_grain)
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


# ═══ STAGE 2: Physics-based Scaling Laws ═══

def scaling_law_ionic(phi_se, f_perc, tau, g_path, gb_d, cn, C=0.073):
    """σ_ionic = σ_brug × C × (G_path × GB_d²)^(1/4) × CN²"""
    if tau <= 0 or phi_se <= 0 or g_path <= 0 or gb_d <= 0 or cn <= 0:
        return 0
    sigma_brug = SIGMA_GRAIN * phi_se * (f_perc / 100) / tau ** 2
    return sigma_brug * C * (g_path * gb_d ** 2) ** 0.25 * cn ** 2

def scaling_law_electronic(phi_am, am_cn, thickness, d_am, C=0.015):
    """σ_el = C × σ_AM × φ_AM^(3/2) × CN_AM² × exp(π/(T/d))"""
    if phi_am <= 0 or am_cn <= 0 or d_am <= 0 or thickness <= 0:
        return 0
    ratio = thickness / d_am
    if ratio <= 0:
        return 0
    return C * SIGMA_AM * phi_am ** 1.5 * am_cn ** 2 * np.exp(np.pi / ratio)

def scaling_law_thermal(sigma_ion, phi_am, cn_se, C=286):
    """σ_th = C × σ_ion^(3/4) × φ_AM² / CN_SE"""
    if sigma_ion <= 0 or phi_am <= 0 or cn_se <= 0:
        return 0
    return C * sigma_ion ** 0.75 * phi_am ** 2 / cn_se


# ═══ STAGE 1: GPR Training ═══

def train_stage1(rows):
    """Train GPR for DEM input → microstructure."""
    from sklearn.gaussian_process import GaussianProcessRegressor
    from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel
    from sklearn.preprocessing import StandardScaler
    from sklearn.model_selection import LeaveOneOut

    input_features = ['d_se', 'd_am', 'am_pct', 'ps_frac', 'rve', 'loading',
                      'd_ratio', 'am_loading', 'se_density_proxy', 'layer_count']
    micro_targets = ['phi_se', 'phi_am', 'sigma_brug', 'tau', 'cn', 'gb_d', 'g_path', 'hop_area',
                     'f_perc', 'thickness', 'porosity', 'am_cn',
                     'sigma_ion']  # direct σ prediction for ensemble

    # Build arrays
    X = np.array([[r[f] for f in input_features] for r in rows])
    Y = {t: np.array([r[t] for r in rows]) for t in micro_targets}

    scaler_X = StandardScaler()
    X_scaled = scaler_X.fit_transform(X)

    models = {}
    print(f"\n  {'Target':15s} {'Train R²':>10s} {'LOO-CV R²':>10s} {'|err|':>8s}")
    print(f"  {'─'*48}")

    for target in micro_targets:
        y = Y[target]
        if np.std(y) == 0:
            continue

        use_log = np.all(y > 0) and (np.max(y) / np.min(y) > 3)
        y_train = np.log(y) if use_log else y.copy()

        scaler_y = StandardScaler()
        y_scaled = scaler_y.fit_transform(y_train.reshape(-1, 1)).ravel()

        kernel = ConstantKernel(1.0) * Matern(length_scale=1.0, nu=2.5) + WhiteKernel(0.1)
        gpr = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=5, alpha=1e-6)
        gpr.fit(X_scaled, y_scaled)

        # Train R²
        pred_train = gpr.predict(X_scaled)
        pred_train_real = scaler_y.inverse_transform(pred_train.reshape(-1, 1)).ravel()
        if use_log:
            pred_train_real = np.exp(pred_train_real)
        ss_res_t = np.sum((y - pred_train_real) ** 2)
        ss_tot = np.sum((y - np.mean(y)) ** 2)
        r2_train = 1 - ss_res_t / ss_tot if ss_tot > 0 else 0

        # LOO-CV
        loo = LeaveOneOut()
        y_pred_loo = np.zeros(len(y))
        for train_idx, test_idx in loo.split(X_scaled):
            gpr_cv = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=2, alpha=1e-6)
            gpr_cv.fit(X_scaled[train_idx], y_scaled[train_idx])
            y_pred_loo[test_idx] = gpr_cv.predict(X_scaled[test_idx])

        y_loo_real = scaler_y.inverse_transform(y_pred_loo.reshape(-1, 1)).ravel()
        if use_log:
            y_loo_real = np.exp(y_loo_real)

        ss_res_cv = np.sum((y - y_loo_real) ** 2)
        r2_cv = 1 - ss_res_cv / ss_tot if ss_tot > 0 else 0
        valid = y > 0
        mean_err = np.mean(np.abs(y_loo_real[valid] - y[valid]) / y[valid]) * 100 if valid.any() else 0

        star = "★" if r2_cv > 0.8 else ""
        print(f"  {target:15s} {r2_train:10.3f} {r2_cv:10.3f} {mean_err:7.1f}% {star}")

        models[target] = {
            'gpr': gpr, 'scaler_X': scaler_X, 'scaler_y': scaler_y,
            'use_log': use_log, 'r2_cv': r2_cv
        }

    return models, input_features, micro_targets


def predict_microstructure(models, input_features, input_dict):
    """Stage 1: Predict microstructure from DEM input."""
    x = np.array([[input_dict.get(f, 0) for f in input_features]])

    micro = {}
    for target, m in models.items():
        x_scaled = m['scaler_X'].transform(x)
        pred_scaled, std_scaled = m['gpr'].predict(x_scaled, return_std=True)
        pred = m['scaler_y'].inverse_transform(pred_scaled.reshape(-1, 1)).ravel()[0]
        std = std_scaled[0] * m['scaler_y'].scale_[0]
        if m['use_log']:
            pred = np.exp(pred)
            std = pred * abs(std)
        micro[target] = {'value': float(pred), 'std': float(abs(std))}

    return micro


def evaluate_hybrid(rows, models, input_features):
    """Evaluate full pipeline: input → GPR → Scaling Law → σ."""
    print(f"\n{'═'*60}")
    print(f"HYBRID EVALUATION (Stage1:GPR + Stage2:ScalingLaw)")
    print(f"{'═'*60}")

    sigma_actual = np.array([r['sigma_ion'] for r in rows])

    # Fit C from true microstructure data
    log_rhs = []
    valid_fit = []
    for i, r in enumerate(rows):
        if r['g_path'] > 0 and r['gb_d'] > 0 and r['cn'] > 0 and r['tau'] > 0 and r['phi_se'] > 0 and r['f_perc'] > 0 and sigma_actual[i] > 0.01:
            sb = SIGMA_GRAIN * r['phi_se'] * (r['f_perc']/100) / r['tau']**2
            rhs = np.log(sb) + 0.25*np.log(r['g_path']*r['gb_d']**2) + 2*np.log(r['cn'])
            log_rhs.append(rhs)
            valid_fit.append(i)
    if valid_fit:
        C_fit = np.exp(np.mean(np.log(sigma_actual[valid_fit]) - np.array(log_rhs)))
    else:
        C_fit = 0.073

    sigma_hybrid = []
    sigma_scaling_only = []

    for r in rows:
        # Stage 1: GPR predict microstructure
        input_dict = {f: r[f] for f in input_features}
        micro = predict_microstructure(models, input_features, input_dict)

        # Stage 2: Scaling law with GPR-predicted microstructure
        # Use σ_brug directly if predicted (bypasses τ problem!)
        if 'sigma_brug' in micro and micro['sigma_brug']['value'] > 0:
            sigma_brug_pred = micro['sigma_brug']['value'] * SIGMA_GRAIN
            sigma_h = sigma_brug_pred * C_fit * \
                      (micro['g_path']['value'] * micro['gb_d']['value']**2)**0.25 * \
                      micro['cn']['value']**2
        else:
            sigma_h = scaling_law_ionic(
                micro['phi_se']['value'], micro['f_perc']['value'],
                micro['tau']['value'], micro['g_path']['value'],
                micro['gb_d']['value'], micro['cn']['value'], C=C_fit
            )
        sigma_hybrid.append(sigma_h)

        # Scaling law with TRUE microstructure (for comparison)
        sigma_s = scaling_law_ionic(
            r['phi_se'], r['f_perc'], r['tau'],
            r['g_path'], r['gb_d'], r['cn'], C=C_fit
        )
        sigma_scaling_only.append(sigma_s)

    sigma_hybrid = np.array(sigma_hybrid)
    sigma_scaling_only = np.array(sigma_scaling_only)

    # R² calculations
    valid = (sigma_actual > 0.01) & (sigma_hybrid > 0) & (sigma_scaling_only > 0)

    def calc_r2(pred, actual, mask):
        log_p, log_a = np.log(pred[mask]), np.log(actual[mask])
        ss_res = np.sum((log_a - log_p) ** 2)
        ss_tot = np.sum((log_a - np.mean(log_a)) ** 2)
        return 1 - ss_res / ss_tot if ss_tot > 0 else 0

    def calc_err(pred, actual, mask):
        return np.mean(np.abs(pred[mask] - actual[mask]) / actual[mask]) * 100

    r2_hybrid = calc_r2(sigma_hybrid, sigma_actual, valid)
    r2_scaling = calc_r2(sigma_scaling_only, sigma_actual, valid)
    err_hybrid = calc_err(sigma_hybrid, sigma_actual, valid)
    err_scaling = calc_err(sigma_scaling_only, sigma_actual, valid)

    # GPR direct σ prediction
    sigma_gpr_direct = []
    for r in rows:
        input_dict = {f: r[f] for f in input_features}
        micro = predict_microstructure(models, input_features, input_dict)
        if 'sigma_ion' in micro and micro['sigma_ion']['value'] > 0:
            sigma_gpr_direct.append(micro['sigma_ion']['value'])
        else:
            sigma_gpr_direct.append(0)
    sigma_gpr_direct = np.array(sigma_gpr_direct)
    valid_gpr = valid & (sigma_gpr_direct > 0)
    r2_gpr = calc_r2(sigma_gpr_direct, sigma_actual, valid_gpr) if valid_gpr.any() else 0
    err_gpr = calc_err(sigma_gpr_direct, sigma_actual, valid_gpr) if valid_gpr.any() else 0

    # Ensemble: weighted average of Hybrid and GPR direct
    w_hybrid = 0.7  # physics-based gets more weight
    sigma_ensemble = w_hybrid * sigma_hybrid + (1 - w_hybrid) * sigma_gpr_direct
    valid_ens = valid & (sigma_ensemble > 0)
    r2_ensemble = calc_r2(sigma_ensemble, sigma_actual, valid_ens) if valid_ens.any() else 0
    err_ensemble = calc_err(sigma_ensemble, sigma_actual, valid_ens) if valid_ens.any() else 0

    # Optimal ensemble weight search
    best_w, best_r2_w = 0.5, 0
    for w in np.arange(0.0, 1.05, 0.05):
        s_ens = w * sigma_hybrid + (1 - w) * sigma_gpr_direct
        v_ens = valid & (s_ens > 0)
        if v_ens.any():
            r2_w = calc_r2(s_ens, sigma_actual, v_ens)
            if r2_w > best_r2_w:
                best_w, best_r2_w = w, r2_w

    sigma_best_ens = best_w * sigma_hybrid + (1 - best_w) * sigma_gpr_direct
    valid_best = valid & (sigma_best_ens > 0)
    err_best = calc_err(sigma_best_ens, sigma_actual, valid_best) if valid_best.any() else 0

    print(f"\n  σ_ionic prediction (n={valid.sum()}):")
    print(f"  {'Method':40s} {'R²':>8s} {'|err|':>8s}")
    print(f"  {'─'*60}")
    print(f"  {'Scaling Law (true microstructure)':40s} {r2_scaling:8.3f} {err_scaling:7.1f}%")
    print(f"  {'Hybrid (GPR micro + Scaling Law)':40s} {r2_hybrid:8.3f} {err_hybrid:7.1f}%")
    print(f"  {'GPR direct (input → σ)':40s} {r2_gpr:8.3f} {err_gpr:7.1f}%")
    print(f"  {'Ensemble (70% Hybrid + 30% GPR)':40s} {r2_ensemble:8.3f} {err_ensemble:7.1f}%")
    print(f"  {'Ensemble (optimal w={best_w:.2f})':40s} {best_r2_w:8.3f} {err_best:7.1f}%")

    if best_r2_w > 0.8:
        print(f"\n  ★ Best R²={best_r2_w:.3f} (w_hybrid={best_w:.2f}) — DEM 없이 σ_ionic 예측!")

    return best_r2_w, r2_scaling


def interactive_predict(models, input_features):
    """Interactive prediction mode."""
    print(f"\n{'═'*60}")
    print(f"INTERACTIVE PREDICTION")
    print(f"{'═'*60}")
    print("Enter electrode design parameters:\n")

    try:
        d_se = float(input("  SE diameter (μm) [0.5/1.0/1.5]: ") or 1.0)
        d_am = float(input("  AM diameter (μm) [5/10]: ") or 5.0)
        am_pct = float(input("  AM:SE - AM% [75/80/85]: ") or 80)
        ps_frac = float(input("  P:S fraction (0=S only, 1=P only) [0.5]: ") or 0.5)
        rve = float(input("  RVE size (μm) [50]: ") or 50)
        loading = float(input("  Loading (mAh/cm²) [1/6/8]: ") or 6)
    except (ValueError, EOFError):
        d_se, d_am, am_pct, ps_frac, rve, loading = 1.0, 5.0, 80, 0.5, 50, 6

    input_dict = {'d_se': d_se, 'd_am': d_am, 'am_pct': am_pct,
                  'ps_frac': ps_frac, 'rve': rve, 'loading': loading,
                  'd_ratio': d_se / d_am if d_am > 0 else 0.1,
                  'am_loading': am_pct * loading / 100,
                  'se_density_proxy': (100 - am_pct) / max(d_se, 0.1),
                  'layer_count': loading * 30 / max(d_se, 0.1)}

    print(f"\n{'─'*50}")
    print(f"INPUT: d_SE={d_se}μm, d_AM={d_am}μm, AM:SE={am_pct}:{100-am_pct}, P:S={ps_frac*10:.0f}:{(1-ps_frac)*10:.0f}, RVE={rve}μm")
    print(f"{'─'*50}")

    # Stage 1
    print(f"\n  ── Stage 1: GPR → Microstructure ──")
    micro = predict_microstructure(models, input_features, input_dict)
    for t in ['phi_se', 'phi_am', 'porosity', 'tau', 'cn', 'gb_d', 'g_path', 'hop_area', 'f_perc', 'thickness']:
        if t in micro:
            v = micro[t]
            print(f"    {t:15s} = {v['value']:10.4f} (±{v['std']:.4f})")

    # Stage 2
    print(f"\n  ── Stage 2: Scaling Law → Conductivity ──")

    sigma_ion = scaling_law_ionic(
        micro['phi_se']['value'], micro['f_perc']['value'],
        micro['tau']['value'], micro['g_path']['value'],
        micro['gb_d']['value'], micro['cn']['value']
    )
    print(f"    σ_ionic    = {sigma_ion:.4f} mS/cm")

    if micro.get('am_cn') and micro['am_cn']['value'] > 0:
        sigma_el = scaling_law_electronic(
            micro['phi_am']['value'], micro['am_cn']['value'],
            micro['thickness']['value'], d_am
        )
        print(f"    σ_electronic = {sigma_el:.2f} mS/cm")
    else:
        sigma_el = 0
        print(f"    σ_electronic = N/A (AM CN not predicted)")

    sigma_th = scaling_law_thermal(sigma_ion, micro['phi_am']['value'], micro['cn']['value'])
    print(f"    σ_thermal  = {sigma_th:.3f} mS/cm equiv")

    print(f"\n  ── Summary ──")
    print(f"    Rate-limiting: {'ionic' if sigma_ion < sigma_el else 'electronic'}")
    if sigma_ion > 0:
        print(f"    σ_el/σ_ion = {sigma_el/sigma_ion:.0f}× (electronic is {sigma_el/sigma_ion:.0f}× better)")


def parameter_sweep(models, input_features):
    """Sweep key parameters to find optimal design."""
    print(f"\n{'═'*60}")
    print(f"PARAMETER SWEEP: Finding optimal electrode design")
    print(f"{'═'*60}")

    best_sigma = 0
    best_params = {}

    print(f"\n  Sweeping d_SE × AM% × P:S...")
    print(f"  {'d_SE':>5s} {'AM%':>5s} {'P:S':>5s} {'φ_SE':>6s} {'τ':>6s} {'CN':>6s} {'σ_ion':>8s} {'σ_el':>8s}")
    print(f"  {'─'*55}")

    for d_se in [0.5, 1.0, 1.5, 3.0]:
        for am_pct in [62, 75, 80, 85]:
            for ps_frac in [0.0, 0.3, 0.5, 0.7, 1.0]:
                d_am_sw = 5.0
                input_dict = {'d_se': d_se, 'd_am': d_am_sw, 'am_pct': am_pct,
                             'ps_frac': ps_frac, 'rve': 50, 'loading': 6,
                             'd_ratio': d_se / d_am_sw,
                             'am_loading': am_pct * 6 / 100,
                             'se_density_proxy': (100 - am_pct) / max(d_se, 0.1),
                             'layer_count': 6 * 30 / max(d_se, 0.1)}
                micro = predict_microstructure(models, input_features, input_dict)

                sigma_ion = scaling_law_ionic(
                    micro['phi_se']['value'], micro['f_perc']['value'],
                    micro['tau']['value'], micro['g_path']['value'],
                    micro['gb_d']['value'], micro['cn']['value']
                )

                if sigma_ion > best_sigma:
                    best_sigma = sigma_ion
                    best_params = {'d_se': d_se, 'am_pct': am_pct, 'ps_frac': ps_frac}

                if sigma_ion > 0.1:  # only show promising
                    print(f"  {d_se:5.1f} {am_pct:5.0f} {ps_frac:5.1f} "
                          f"{micro['phi_se']['value']:6.3f} {micro['tau']['value']:6.2f} "
                          f"{micro['cn']['value']:6.2f} {sigma_ion:8.4f}")

    print(f"\n  ★ OPTIMAL: d_SE={best_params.get('d_se')}μm, "
          f"AM%={best_params.get('am_pct')}, P:S={best_params.get('ps_frac', 0.5)*10:.0f}:{(1-best_params.get('ps_frac', 0.5))*10:.0f}")
    print(f"    σ_ionic = {best_sigma:.4f} mS/cm")


def main():
    import argparse
    parser = argparse.ArgumentParser(description='2-Stage Hybrid Predictor')
    parser.add_argument('--predict', action='store_true', help='Interactive prediction')
    parser.add_argument('--sweep', action='store_true', help='Parameter sweep')
    args = parser.parse_args()

    print("═" * 60)
    print("2-STAGE HYBRID PREDICTOR")
    print("Stage 1: GPR (input → microstructure)")
    print("Stage 2: Scaling Law (microstructure → σ)")
    print("═" * 60)

    rows = load_data()
    print(f"\nData: {len(rows)} cases")

    if len(rows) < 10:
        print("Need more data!")
        return

    # Train Stage 1
    print(f"\n{'─'*60}")
    print(f"STAGE 1: Training GPR (DEM input → microstructure)")
    print(f"{'─'*60}")
    models, input_features, micro_targets = train_stage1(rows)

    # Evaluate hybrid
    evaluate_hybrid(rows, models, input_features)

    if args.predict:
        interactive_predict(models, input_features)

    if args.sweep:
        parameter_sweep(models, input_features)

    if not args.predict and not args.sweep:
        print(f"\n  Try: --predict (interactive) or --sweep (find optimal design)")


if __name__ == '__main__':
    main()
