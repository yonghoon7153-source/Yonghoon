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
    'hop_area', 'f_perc', 'thickness', 'porosity', 'am_cn', 'sigma_ion',
    'coverage'  # AM-SE contact fraction for Newman utilization
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
                # Coverage: fraction of AM surface in contact with SE
                'coverage': max(m.get('coverage_AM_P_mean', 0), m.get('coverage_AM_S_mean', 0), m.get('coverage_AM_mean', 0)) / 100 if max(m.get('coverage_AM_P_mean', 0), m.get('coverage_AM_S_mean', 0), m.get('coverage_AM_mean', 0)) > 0 else 0.20,
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

    # ── Utilization: Newman Electrode Model (Corrected) ──
    # Two versions: A) Newman simplified, B) PyBaMM (if available)
    #
    # Newman: util = tanh(ν)/ν where ν = T/L_c
    # L_c = characteristic length = √(κ_eff / (a_s × j₀ × F / (R×T_K)))
    # Ref: Newman & Tobias (1962), Fuller, Doyle, Newman (1994)
    #
    # Key correction: j₀ is FIXED (exchange current density, material property)
    # C-rate affects applied current i_app, which changes overpotential η
    # At high η (high C-rate), effective reaction rate increases → ν changes

    F_const = 96485      # C/mol
    R_const = 8.314      # J/(mol·K)
    T_kelvin = 298       # K
    r_AM_cm = d_am * 1e-4 / 2 if d_am > 0 else 2.5e-4  # μm → cm
    T_cm = thickness * 1e-4 if thickness > 0 else 0.01   # μm → cm
    kappa_eff = sigma_ionic_final * 1e-3 if sigma_ionic_final > 0 else 1e-6  # mS/cm → S/cm

    # Specific interfacial area: a_s = 3 × φ_AM / r_AM (spherical particles)
    # But in ASSB, not all AM surface contacts SE — apply coverage fraction
    # Use GPR-predicted coverage (NOT fixed 20%)
    a_s_geometric = 3 * phi_am / r_AM_cm if r_AM_cm > 0 and phi_am > 0 else 1e4  # cm⁻¹
    coverage_pred = micro.get('coverage', {}).get('value', 0.20) if isinstance(micro.get('coverage'), dict) else 0.20
    coverage_fraction = max(0.05, min(0.50, coverage_pred))  # clamp to reasonable range
    a_s = a_s_geometric * coverage_fraction  # effective interfacial area

    # Exchange current density: FIXED material property (NOT C-rate dependent!)
    # NCM811/LPSCl interface: j₀ ≈ 0.01 mA/cm² (10× lower than liquid LIB)
    # ASSB solid-solid interface has much higher charge-transfer resistance
    # Ref: Sakuda (2017), Kato (2016)
    j0 = 0.01e-3  # A/cm² (0.01 mA/cm², ASSB NCM/LPSCl interface)

    # Theoretical capacity for C-rate calculation
    # NCM811: 200 mAh/g, ρ_AM = 4.8 g/cm³
    Q_theo = 200  # mAh/g
    rho_AM = 4.8  # g/cm³
    # Volumetric capacity = Q × ρ × φ_AM
    Q_vol = Q_theo * rho_AM * phi_am * 1e-3  # Ah/cm³

    c_rates = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0]
    utilizations = {}
    for c_rate in c_rates:
        # Applied current density (A/cm²) = C-rate × Q_vol × T
        i_app = c_rate * Q_vol * T_cm  # A/cm²

        # Newman dimensionless parameter
        # ν² = a_s × j₀ × F × T² / (κ_eff × R × T_K)
        # This gives characteristic penetration depth vs electrode thickness
        if kappa_eff > 0 and a_s > 0 and j0 > 0:
            nu_sq = a_s * j0 * F_const * T_cm**2 / (kappa_eff * R_const * T_kelvin)
            nu = np.sqrt(nu_sq) if nu_sq > 0 else 0

            # Base utilization from Newman
            util_newman = np.tanh(nu) / nu if nu > 0.01 else 1.0

            # Additional C-rate penalty: ohmic drop limits
            # ΔV_ohmic = i_app × T / (2 × κ_eff)
            delta_V = i_app * T_cm / (2 * kappa_eff) if kappa_eff > 0 else 10
            V_cutoff = 0.3  # V (typical cutoff overpotential)
            ohmic_factor = min(1.0, V_cutoff / delta_V) if delta_V > 0 else 1.0

            util = util_newman * ohmic_factor
        else:
            util = 1.0

        utilizations[f'{c_rate}C'] = round(min(1.0, max(0.0, util)), 3)

    util_1C = utilizations.get('1.0C', 1.0)
    performance_score = sigma_ionic_final * util_1C

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
        'utilization': utilizations,
        'performance_score': round(performance_score, 4),
        'rate_limiting': rate_limiting,
        'input': input_dict,
    }


def _fmt(val, micro_entry):
    """Format a microstructure value with uncertainty."""
    return {
        'value': round(val, 4),
        'std': round(micro_entry.get('std', 0), 4),
    }


def sweep_optimal(top_n=5, fixed_params=None, sweep_keys=None, defaults=None):
    """Sweep unchecked parameters, keep checked ones fixed."""
    if _cached_models is None:
        return {'error': 'Models not trained.'}

    if fixed_params is None:
        fixed_params = {}
    if defaults is None:
        defaults = {}

    # Sweep ranges for each parameter
    sweep_ranges = {
        'd_se': np.arange(0.3, 3.1, 0.1),
        'd_am_p': np.arange(3, 13, 1),
        'd_am_s': np.arange(1, 9, 0.5),
        'am_pct': np.arange(60, 91, 2),
        'ps_frac': np.arange(0.0, 1.05, 0.1),
        'loading': np.arange(1, 11, 1),
        'rve': np.array([30, 50, 100]),
    }

    if sweep_keys is None:
        sweep_keys = ['d_se', 'am_pct', 'ps_frac', 'loading']

    # Build combinations for sweep params only
    import itertools
    sweep_vals = [sweep_ranges.get(k, np.array([defaults.get(k, 1.0)])) for k in sweep_keys]
    combos = list(itertools.product(*sweep_vals))

    # Base params from fixed + defaults
    base = {
        'd_se': 1.0, 'd_am_p': 10, 'd_am_s': 4, 'am_pct': 80,
        'ps_frac': 0.5, 'loading': 6, 'rve': 50,
    }
    base.update(defaults)
    base.update(fixed_params)

    results = []
    for combo in combos:
        params = dict(base)
        for i, k in enumerate(sweep_keys):
            params[k] = float(combo[i])
        d_am = max(params.get('d_am_p', 5), params.get('d_am_s', 4))
        pred = predict(params['d_se'], d_am, params['am_pct'],
                      params['ps_frac'], params['loading'], params['rve'])
        if 'error' in pred:
            continue
        sig = pred['conductivity']['sigma_ionic']
        if sig > 0:
            r = {
                'd_se': round(float(params['d_se']), 1),
                'am_pct': int(params['am_pct']),
                'ps_frac': round(float(params['ps_frac']), 1),
                'loading': round(float(params['loading']), 1),
                'sigma_ionic': round(sig, 4),
                'sigma_electronic': pred['conductivity']['sigma_electronic'],
                'sigma_thermal': pred['conductivity']['sigma_thermal'],
                'rate_limiting': pred['rate_limiting'],
            }
            # Add sweep params to result
            for k in sweep_keys:
                if k not in r:
                    r[k] = round(float(params[k]), 1)
            results.append(r)

    results.sort(key=lambda x: x['sigma_ionic'], reverse=True)
    return results[:top_n]
