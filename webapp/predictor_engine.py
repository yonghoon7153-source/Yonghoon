"""
ML Predictor Engine: 2-Stage Hybrid (GPR + Physics Scaling Laws)
Stage 1: GPR predicts microstructure from DEM design inputs
Stage 2: Physics scaling laws compute conductivity from microstructure
"""
import json
import os
import warnings
import numpy as np
from pathlib import Path

warnings.filterwarnings('ignore')

# ═══ Physical Constants ═══
SIGMA_GRAIN = 3.0       # mS/cm (SE grain interior)
SIGMA_AM = 50.0         # mS/cm (NCM811)

INPUT_FEATURES = [
    'd_se', 'd_am', 'am_pct', 'ps_frac', 'rve', 'loading',
    'd_ratio', 'am_loading', 'se_density_proxy', 'layer_count'
]

MICRO_TARGETS = [
    'phi_se', 'phi_am', 'sigma_brug', 'tau', 'cn', 'gb_d', 'g_path',
    'hop_area', 'f_perc', 'thickness', 'porosity', 'am_cn', 'sigma_ion'
]

# Cached models (module-level)
_cached_models = None
_training_data_count = 0


def derive_features(d_se, d_am, am_pct, ps_frac, rve, loading):
    """Compute derived features from basic inputs."""
    return {
        'd_se': d_se,
        'd_am': d_am,
        'am_pct': am_pct,
        'ps_frac': ps_frac,
        'rve': rve,
        'loading': loading,
        'd_ratio': d_se / d_am if d_am > 0 else 0.1,
        'am_loading': am_pct * loading / 100,
        'se_density_proxy': (100 - am_pct) / max(d_se, 0.1),
        'layer_count': loading * 30 / max(d_se, 0.1),
    }


def load_training_data(results_folder, archive_folder):
    """Load all cases with DEM input + output from results and archive."""
    rows = []
    for base in [Path(results_folder), Path(archive_folder)]:
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
            except Exception:
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
                except Exception:
                    pass

            ps = m.get('ps_ratio', '')
            ps_frac = 0.5
            if isinstance(ps, str) and ':' in ps:
                try:
                    parts = ps.split(':')
                    p, s = float(parts[0]), float(parts[1])
                    ps_frac = p / (p + s) if (p + s) > 0 else 0.5
                except Exception:
                    pass

            box_x = ip.get('box_x', 0.05)
            rve = box_x / scale * 1e6 if scale > 0 else 50

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

            if tau > 10 or m.get('porosity', 0) > 30:
                continue
            if 0 < sigma_ion < 0.01:
                continue

            # Loading
            loading = 0
            folder_path = str(met_path.parent)
            if '1mAh' in folder_path or 'thin' in folder_path.lower():
                loading = 1
            elif '8mAh' in folder_path:
                loading = 8
            elif '6mAh' in folder_path or 'real' in folder_path.lower() or 'Real' in folder_path:
                loading = 6
            elif 'particulate' in folder_path.lower():
                t = m.get('thickness_um', 0)
                if t > 0:
                    loading = round(t * 0.05, 1)
            if loading <= 0 and m.get('thickness_um', 0) > 0:
                loading = round(m['thickness_um'] * 0.05, 1)

            d_am = max(d_am_p, d_am_s)

            rows.append({
                'name': case_dir.name,
                'd_se': d_se, 'd_am': d_am, 'am_pct': am_pct, 'ps_frac': ps_frac,
                'rve': rve, 'loading': loading,
                'd_ratio': d_se / d_am if d_am > 0 else 0.1,
                'am_loading': am_pct * loading / 100,
                'se_density_proxy': (100 - am_pct) / max(d_se, 0.1),
                'layer_count': loading * 30 / max(d_se, 0.1),
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
                'sigma_ion': sigma_ion,
                'sigma_el': m.get('electronic_sigma_full_mScm', 0),
                'sigma_th': m.get('thermal_sigma_full_mScm', 0),
                'sigma_brug': m.get('sigma_ratio', 0),
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


def train_models(results_folder, archive_folder):
    """Train GPR models and cache them."""
    global _cached_models, _training_data_count
    from sklearn.gaussian_process import GaussianProcessRegressor
    from sklearn.gaussian_process.kernels import ConstantKernel, Matern, WhiteKernel
    from sklearn.preprocessing import StandardScaler

    rows = load_training_data(results_folder, archive_folder)
    _training_data_count = len(rows)

    if len(rows) < 3:
        _cached_models = None
        return {'success': False, 'message': f'Not enough training data ({len(rows)} cases, need >= 3)', 'count': len(rows)}

    X = np.array([[r[f] for f in INPUT_FEATURES] for r in rows])
    Y = {t: np.array([r[t] for r in rows]) for t in MICRO_TARGETS}

    scaler_X = StandardScaler()
    X_scaled = scaler_X.fit_transform(X)

    models = {}
    scores = {}

    for target in MICRO_TARGETS:
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

        # Train R2
        pred_train = gpr.predict(X_scaled)
        pred_real = scaler_y.inverse_transform(pred_train.reshape(-1, 1)).ravel()
        if use_log:
            pred_real = np.exp(pred_real)
        ss_res = np.sum((y - pred_real) ** 2)
        ss_tot = np.sum((y - np.mean(y)) ** 2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0

        models[target] = {
            'gpr': gpr, 'scaler_X': scaler_X, 'scaler_y': scaler_y,
            'use_log': use_log, 'r2': r2
        }
        scores[target] = round(r2, 3)

    # Fit C constant from data
    C_ionic = 0.073
    log_rhs = []
    log_act = []
    for r in rows:
        if r['g_path'] > 0 and r['gb_d'] > 0 and r['cn'] > 0 and r['sigma_ion'] > 0.01 and r['sigma_brug'] > 0:
            sb = r['sigma_brug'] * SIGMA_GRAIN
            rhs = np.log(sb) + 0.25 * np.log(r['g_path'] * r['gb_d'] ** 2) + 2 * np.log(r['cn'])
            log_rhs.append(rhs)
            log_act.append(np.log(r['sigma_ion']))
    if log_rhs:
        C_ionic = float(np.exp(np.mean(np.array(log_act) - np.array(log_rhs))))

    _cached_models = {
        'models': models,
        'C_ionic': C_ionic,
        'count': len(rows),
        'scores': scores,
    }

    return {'success': True, 'count': len(rows), 'scores': scores, 'C_ionic': round(C_ionic, 4)}


def get_data_count(results_folder, archive_folder):
    """Get count of available training cases without training."""
    global _training_data_count
    if _training_data_count > 0:
        return _training_data_count
    rows = load_training_data(results_folder, archive_folder)
    _training_data_count = len(rows)
    return _training_data_count


def predict(d_se, d_am, am_pct, ps_frac, loading, rve):
    """Run prediction using cached models."""
    if _cached_models is None:
        return {'error': 'Models not trained. Click "Train Models" first.'}

    models = _cached_models['models']
    C_ionic = _cached_models['C_ionic']

    input_dict = derive_features(d_se, d_am, am_pct, ps_frac, rve, loading)

    # Stage 1: GPR predict microstructure
    micro = {}
    for target, m in models.items():
        x = np.array([[input_dict.get(f, 0) for f in INPUT_FEATURES]])
        x_scaled = m['scaler_X'].transform(x)
        pred_scaled, std_scaled = m['gpr'].predict(x_scaled, return_std=True)
        pred = m['scaler_y'].inverse_transform(pred_scaled.reshape(-1, 1)).ravel()[0]
        std = std_scaled[0] * m['scaler_y'].scale_[0]
        if m['use_log']:
            pred = np.exp(pred)
            std = pred * abs(std)
        micro[target] = {'value': float(pred), 'std': float(abs(std))}

    # Stage 2: Physics scaling laws
    phi_se = micro.get('phi_se', {}).get('value', 0)
    phi_am = micro.get('phi_am', {}).get('value', 0)
    f_perc = micro.get('f_perc', {}).get('value', 100)
    tau = micro.get('tau', {}).get('value', 1)
    cn = micro.get('cn', {}).get('value', 0)
    gb_d = micro.get('gb_d', {}).get('value', 0)
    g_path = micro.get('g_path', {}).get('value', 0)
    am_cn = micro.get('am_cn', {}).get('value', 0)
    thickness = micro.get('thickness', {}).get('value', 0)
    porosity = micro.get('porosity', {}).get('value', 0)
    hop_area = micro.get('hop_area', {}).get('value', 0)
    sigma_brug_ratio = micro.get('sigma_brug', {}).get('value', 0)

    # Ionic conductivity
    sigma_brug = SIGMA_GRAIN * sigma_brug_ratio
    sigma_ionic = 0
    if sigma_brug > 0 and g_path > 0 and gb_d > 0 and cn > 0:
        sigma_ionic = sigma_brug * C_ionic * (g_path * gb_d ** 2) ** 0.25 * cn ** 2

    # Direct GPR prediction for ensemble
    sigma_gpr = micro.get('sigma_ion', {}).get('value', 0)

    # Ensemble (0.5 hybrid + 0.5 gpr by default)
    if sigma_ionic > 0 and sigma_gpr > 0:
        sigma_ionic_final = 0.5 * sigma_ionic + 0.5 * sigma_gpr
    elif sigma_ionic > 0:
        sigma_ionic_final = sigma_ionic
    else:
        sigma_ionic_final = sigma_gpr

    # Electronic conductivity
    sigma_electronic = 0
    if phi_am > 0 and am_cn > 0 and d_am > 0 and thickness > 0:
        ratio = thickness / d_am
        if ratio > 0:
            sigma_electronic = 0.015 * SIGMA_AM * phi_am ** 1.5 * am_cn ** 2 * np.exp(np.pi / ratio)

    # Thermal conductivity
    sigma_thermal = 0
    if sigma_ionic_final > 0 and phi_am > 0 and cn > 0:
        sigma_thermal = 286 * sigma_ionic_final ** 0.75 * phi_am ** 2 / cn

    # Uncertainty from GPR std
    sigma_ion_std = micro.get('sigma_ion', {}).get('std', 0)

    # Rate limiting
    rate_limiting = 'ionic'
    if sigma_ionic_final > 0 and sigma_electronic > 0:
        if sigma_electronic < sigma_ionic_final:
            rate_limiting = 'electronic'

    return {
        'microstructure': {
            'phi_se': _fmt(phi_se, micro.get('phi_se', {})),
            'phi_am': _fmt(phi_am, micro.get('phi_am', {})),
            'porosity': _fmt(porosity, micro.get('porosity', {})),
            'tau': _fmt(tau, micro.get('tau', {})),
            'cn': _fmt(cn, micro.get('cn', {})),
            'gb_d': _fmt(gb_d, micro.get('gb_d', {})),
            'g_path': _fmt(g_path, micro.get('g_path', {})),
            'hop_area': _fmt(hop_area, micro.get('hop_area', {})),
            'f_perc': _fmt(f_perc, micro.get('f_perc', {})),
            'thickness': _fmt(thickness, micro.get('thickness', {})),
            'am_cn': _fmt(am_cn, micro.get('am_cn', {})),
            'sigma_brug_ratio': _fmt(sigma_brug_ratio, micro.get('sigma_brug', {})),
        },
        'conductivity': {
            'sigma_ionic': round(sigma_ionic_final, 4),
            'sigma_ionic_hybrid': round(sigma_ionic, 4),
            'sigma_ionic_gpr': round(sigma_gpr, 4),
            'sigma_ionic_std': round(sigma_ion_std, 4),
            'sigma_electronic': round(sigma_electronic, 4),
            'sigma_thermal': round(sigma_thermal, 4),
        },
        'rate_limiting': rate_limiting,
        'input': input_dict,
    }


def _fmt(val, micro_entry):
    """Format a microstructure value with uncertainty."""
    return {
        'value': round(val, 4),
        'std': round(micro_entry.get('std', 0), 4),
    }


def sweep_optimal(top_n=5):
    """Sweep parameter space to find optimal design."""
    if _cached_models is None:
        return {'error': 'Models not trained.'}

    results = []
    for d_se in np.arange(0.3, 3.1, 0.3):
        for am_pct in range(60, 91, 5):
            for ps_frac in np.arange(0.0, 1.1, 0.2):
                for loading in np.arange(1, 11, 2):
                    pred = predict(d_se, 5.0, am_pct, ps_frac, loading, 50)
                    if 'error' in pred:
                        continue
                    sig = pred['conductivity']['sigma_ionic']
                    if sig > 0:
                        results.append({
                            'd_se': round(float(d_se), 1),
                            'd_am': 5.0,
                            'am_pct': int(am_pct),
                            'ps_frac': round(float(ps_frac), 1),
                            'loading': round(float(loading), 1),
                            'rve': 50,
                            'sigma_ionic': round(sig, 4),
                            'sigma_electronic': pred['conductivity']['sigma_electronic'],
                            'sigma_thermal': pred['conductivity']['sigma_thermal'],
                            'rate_limiting': pred['rate_limiting'],
                        })

    results.sort(key=lambda x: x['sigma_ionic'], reverse=True)
    return results[:top_n]
