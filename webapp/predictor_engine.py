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

# Core targets (always included)
CORE_TARGETS = [
    'phi_se', 'phi_am', 'sigma_brug', 'tau', 'cn', 'gb_d', 'g_path',
    'hop_area', 'f_perc', 'thickness', 'porosity', 'am_cn', 'sigma_ion',
    'coverage'
]
# Additional targets auto-discovered from full_metrics.json
# Will be populated dynamically during data loading
MICRO_TARGETS = list(CORE_TARGETS)  # starts with core, expands in train_models()

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
            # Auto-extract ALL numeric fields from full_metrics.json
            # Skip input features and already-mapped fields
            skip_keys = {'phi_se', 'phi_am', 'porosity', 'se_se_cn', 'tortuosity_mean',
                        'tortuosity_recommended', 'gb_density_mean', 'path_conductance_mean',
                        'path_hop_area_mean', 'percolation_pct', 'thickness_um', 'am_am_cn',
                        'sigma_full_mScm', 'electronic_sigma_full_mScm', 'thermal_sigma_full_mScm',
                        'sigma_ratio', 'sigma_full', 'coverage_AM_P_mean', 'coverage_AM_S_mean',
                        'coverage_AM_mean', 'ps_ratio'}
            for mk, mv in m.items():
                if mk not in skip_keys and isinstance(mv, (int, float)) and not isinstance(mv, bool):
                    safe_key = f'fm_{mk}'  # prefix to avoid collision
                    rows[-1][safe_key] = float(mv)
            rows[-1]['_all_fm_keys'] = [f'fm_{k}' for k, v in m.items()
                                         if k not in skip_keys and isinstance(v, (int, float)) and not isinstance(v, bool)]

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

    # Auto-discover all numeric targets from full_metrics
    global MICRO_TARGETS
    all_fm_keys = set()
    for r in rows:
        all_fm_keys.update(r.get('_all_fm_keys', []))
    # Filter: need at least 80% of cases to have this key, and non-zero variance
    extra_targets = []
    for fk in sorted(all_fm_keys):
        vals = [r.get(fk, 0) for r in rows]
        non_zero = sum(1 for v in vals if v != 0)
        if non_zero >= len(rows) * 0.5 and np.std(vals) > 0:
            extra_targets.append(fk)
    MICRO_TARGETS = list(CORE_TARGETS) + extra_targets
    print(f"  Targets: {len(CORE_TARGETS)} core + {len(extra_targets)} auto-discovered = {len(MICRO_TARGETS)} total")

    X = np.array([[r[f] for f in INPUT_FEATURES] for r in rows])
    Y = {t: np.array([r.get(t, 0) for r in rows]) for t in MICRO_TARGETS}

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


def predict(d_se, d_am, am_pct, ps_frac, loading, rve, temperature=298, additive='none'):
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

    # ── Temperature correction (Arrhenius) ──
    # σ(T) = σ(298) × exp(-Ea/R × (1/T - 1/298))
    # Ea_SE ≈ 0.30 eV (LPSCl argyrodite, Kraft 2017)
    # Ea_AM ≈ 0.50 eV (NCM811 electronic, rough)
    Ea_SE_eV = 0.30
    Ea_AM_eV = 0.50
    k_B = 8.617e-5  # eV/K
    if temperature != 298 and temperature > 0:
        arrhenius_SE = np.exp(-Ea_SE_eV / k_B * (1/temperature - 1/298))
        arrhenius_AM = np.exp(-Ea_AM_eV / k_B * (1/temperature - 1/298))
        sigma_ionic_final *= arrhenius_SE
    else:
        arrhenius_SE = 1.0
        arrhenius_AM = 1.0

    # Electronic conductivity
    sigma_electronic = 0
    if phi_am > 0 and am_cn > 0 and d_am > 0 and thickness > 0:
        ratio = thickness / d_am
        if ratio > 0:
            sigma_electronic = 0.015 * SIGMA_AM * phi_am ** 1.5 * am_cn ** 2 * np.exp(np.pi / ratio)
            sigma_electronic *= arrhenius_AM  # temperature correction

    # ── Conductive Additive / Binder Effect ──
    # Ref: Bielefeld 2023, Minnmann 2021, Kang 2024, Rosner 2026
    # KEY: C65 percolation threshold ~4wt%! Below that, NO electronic network!
    if additive == 'vgcf':
        # VGCF 1wt%: fiber morphology, poor percolation even at 10vol%
        # σ_el ≈ 0.4 mS/cm (Bielefeld 2023) — NOT percolating, just local enhancement
        # Ionic: minimal impact (fiber doesn't block SE paths much)
        sigma_electronic = max(sigma_electronic, 0.4)
        sigma_ionic_final *= 0.99  # ~1% ionic reduction (Nat. Commun. 2025: "minor influence")
    elif additive == 'c65':
        # C65 1wt%: BELOW percolation threshold (~4wt%!)
        # At 1wt%: sub-percolation, no e- network. σ_el barely improves.
        # At ≥4wt%: σ_el ≈ 70 mS/cm (Bielefeld 2023)
        # 1wt% → maybe 0.5~1 mS/cm boost (isolated clusters, not network)
        sigma_electronic = max(sigma_electronic, sigma_electronic + 1.0)
        sigma_ionic_final *= 0.97  # ~3% ionic reduction (SE contact disruption)
    elif additive == 'c65_4wt':
        # C65 4wt%: AT percolation threshold — full electronic network!
        # σ_el ≈ 70 mS/cm (Bielefeld 2023, directly measured)
        # But significant ionic penalty + SE decomposition at C/SE interface
        sigma_electronic = max(sigma_electronic, 70)
        sigma_ionic_final *= 0.85  # ~15% ionic reduction (pathway blocking + SE decomposition)
    elif additive == 'ptfe':
        # PTFE 0.5wt%: dry process binder (fibrillization)
        # Insulator but at 0.1-0.5wt% minimal content
        # Ref: Rosner 2026, Mun 2025 — 209.7 mAh/g, 97.4% retention
        # Effect: slight σ_el reduction (AM surface coating), σ_ion nearly unchanged
        sigma_electronic *= 0.7  # ~30% e- reduction from surface coating
        sigma_ionic_final *= 0.99  # ~1% ionic (dry process, no solvent damage)

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
        'temperature': temperature,
        'arrhenius_factor_SE': round(arrhenius_SE, 4),
        'input': input_dict,
    }


def _fmt(val, micro_entry):
    """Format a microstructure value with uncertainty."""
    return {
        'value': round(val, 4),
        'std': round(micro_entry.get('std', 0), 4),
    }


def generate_heatmap(x_param, y_param, target, fixed_params, n_points=20):
    """Generate 2D heatmap data for two parameters vs a target."""
    if _cached_models is None:
        return {'error': 'Models not trained.'}

    sweep_ranges = {
        'd_se': (0.3, 3.0),
        'd_am': (3, 12),
        'am_pct': (60, 90),
        'ps_frac': (0.0, 1.0),
        'loading': (1, 10),
        'rve': (30, 100),
    }

    base = {'d_se': 1.0, 'd_am': 5.0, 'am_pct': 80, 'ps_frac': 0.5, 'loading': 6, 'rve': 50}
    base.update(fixed_params or {})

    x_range = sweep_ranges.get(x_param, (0.5, 5.0))
    y_range = sweep_ranges.get(y_param, (0.5, 5.0))
    x_values = np.linspace(x_range[0], x_range[1], n_points).tolist()
    y_values = np.linspace(y_range[0], y_range[1], n_points).tolist()

    z_matrix = []
    for yi in y_values:
        row = []
        for xi in x_values:
            params = dict(base)
            params[x_param] = xi
            params[y_param] = yi
            pred = predict(
                d_se=params['d_se'], d_am=params['d_am'],
                am_pct=params['am_pct'], ps_frac=params['ps_frac'],
                loading=params['loading'], rve=params['rve'],
            )
            if 'error' in pred:
                row.append(0)
            elif target == 'sigma_ionic':
                row.append(pred['conductivity']['sigma_ionic'])
            elif target == 'sigma_electronic':
                row.append(pred['conductivity']['sigma_electronic'])
            elif target == 'energy':
                # Quick energy estimate
                phi_am_v = pred['microstructure'].get('phi_am', {}).get('value', 0.5)
                thick_v = pred['microstructure'].get('thickness', {}).get('value', 100)
                util_1c = pred.get('utilization', {}).get('1.0C', 1.0)
                row.append(200 * 4.8 * phi_am_v * thick_v * 1e-4 * 3.7 * util_1c / 1000)
            elif target == 'performance':
                row.append(pred.get('performance_score', 0))
            else:
                row.append(pred['conductivity'].get(target, 0))
            row[-1] = round(row[-1], 6)
        z_matrix.append(row)

    labels = {
        'd_se': 'd_SE (um)', 'd_am': 'd_AM (um)', 'am_pct': 'AM%',
        'ps_frac': 'P:S fraction', 'loading': 'Loading (mAh/cm2)', 'rve': 'RVE (um)',
        'sigma_ionic': 'sigma_ionic (mS/cm)', 'sigma_electronic': 'sigma_elec (mS/cm)',
        'energy': 'Energy (mWh/cm2)', 'performance': 'Perf Score',
    }

    return {
        'x_values': [round(v, 3) for v in x_values],
        'y_values': [round(v, 3) for v in y_values],
        'z_matrix': z_matrix,
        'x_label': labels.get(x_param, x_param),
        'y_label': labels.get(y_param, y_param),
        'z_label': labels.get(target, target),
    }


def sensitivity_analysis(base_params):
    """Vary each param +/-30% and measure sigma change."""
    if _cached_models is None:
        return {'error': 'Models not trained.'}

    d_am = base_params.get('d_am', max(base_params.get('d_am_p', 5), base_params.get('d_am_s', 4)))
    base_pred = predict(
        d_se=base_params['d_se'], d_am=d_am,
        am_pct=base_params['am_pct'], ps_frac=base_params['ps_frac'],
        loading=base_params['loading'], rve=base_params['rve'],
    )
    if 'error' in base_pred:
        return base_pred
    base_sigma = base_pred['conductivity']['sigma_ionic']

    params_to_vary = ['d_se', 'am_pct', 'ps_frac', 'loading']
    results = []
    for p in params_to_vary:
        val = base_params.get(p, 1.0)
        deltas = []
        for factor in [0.7, 1.3]:
            new_val = val * factor
            params_copy = dict(base_params)
            params_copy[p] = new_val
            d_am_c = max(params_copy.get('d_am_p', d_am), params_copy.get('d_am_s', d_am))
            pred = predict(
                d_se=params_copy.get('d_se', 1.0), d_am=d_am_c,
                am_pct=params_copy.get('am_pct', 80), ps_frac=params_copy.get('ps_frac', 0.5),
                loading=params_copy.get('loading', 6), rve=params_copy.get('rve', 50),
            )
            if 'error' not in pred:
                deltas.append(pred['conductivity']['sigma_ionic'] - base_sigma)
        avg_delta = np.mean([abs(d) for d in deltas]) if deltas else 0
        results.append({
            'param': p,
            'delta_sigma': round(float(avg_delta), 6),
            'delta_minus': round(float(deltas[0]), 6) if len(deltas) > 0 else 0,
            'delta_plus': round(float(deltas[1]), 6) if len(deltas) > 1 else 0,
        })

    results.sort(key=lambda x: x['delta_sigma'], reverse=True)
    return {'base_sigma': round(base_sigma, 6), 'sensitivities': results}


def compute_pareto(fixed_params=None, n_points=8):
    """Sweep combinations, compute sigma_ionic and energy_density, return Pareto front."""
    if _cached_models is None:
        return {'error': 'Models not trained.'}

    import itertools
    base = {'d_se': 1.0, 'd_am': 5.0, 'am_pct': 80, 'ps_frac': 0.5, 'loading': 6, 'rve': 50}
    base.update(fixed_params or {})

    d_se_vals = np.linspace(0.3, 3.0, n_points)
    am_pct_vals = np.linspace(60, 90, n_points)
    loading_vals = np.linspace(1, 10, n_points)

    all_points = []
    for d_se, am_pct, loading in itertools.product(d_se_vals, am_pct_vals, loading_vals):
        pred = predict(d_se=d_se, d_am=base['d_am'], am_pct=am_pct,
                       ps_frac=base['ps_frac'], loading=loading, rve=base['rve'])
        if 'error' in pred:
            continue
        sigma = pred['conductivity']['sigma_ionic']
        phi_am_v = pred['microstructure'].get('phi_am', {}).get('value', 0.5)
        thick_v = pred['microstructure'].get('thickness', {}).get('value', 100)
        util_1c = pred.get('utilization', {}).get('1.0C', 1.0)
        energy = 200 * 4.8 * phi_am_v * thick_v * 1e-4 * 3.7 * util_1c / 1000
        if sigma > 0 and energy > 0:
            all_points.append({
                'd_se': round(float(d_se), 2),
                'am_pct': round(float(am_pct), 1),
                'loading': round(float(loading), 1),
                'sigma_ionic': round(float(sigma), 4),
                'energy': round(float(energy), 4),
            })

    # Find Pareto front: no other point dominates in BOTH sigma AND energy
    pareto = []
    for p in all_points:
        dominated = False
        for q in all_points:
            if q['sigma_ionic'] >= p['sigma_ionic'] and q['energy'] >= p['energy'] and \
               (q['sigma_ionic'] > p['sigma_ionic'] or q['energy'] > p['energy']):
                dominated = True
                break
        if not dominated:
            pareto.append(p)

    pareto.sort(key=lambda x: x['sigma_ionic'])
    return {'pareto': pareto, 'all_count': len(all_points)}


def suggest_next(results_folder, archive_folder, top_n=3):
    """Suggest next DEM conditions using GPR uncertainty + predicted value (Expected Improvement)."""
    if _cached_models is None:
        return {'error': 'Models not trained.'}

    models = _cached_models['models']
    if 'sigma_ion' not in models:
        return {'error': 'sigma_ion model not available.'}

    m = models['sigma_ion']

    # Generate candidate points
    import itertools
    d_se_vals = np.linspace(0.3, 3.0, 10)
    am_pct_vals = np.linspace(60, 90, 7)
    ps_frac_vals = np.array([0.0, 0.3, 0.5, 0.7, 1.0])
    loading_vals = np.array([1, 3, 6, 8, 10])

    best_sigma = 0
    rows = load_training_data(results_folder, archive_folder)
    for r in rows:
        if r.get('sigma_ion', 0) > best_sigma:
            best_sigma = r['sigma_ion']

    candidates = []
    for d_se, am_pct, ps_frac, loading in itertools.product(d_se_vals, am_pct_vals, ps_frac_vals, loading_vals):
        d_am = 5.0
        rve = 50.0
        inp = derive_features(d_se, d_am, am_pct, ps_frac, rve, loading)
        x = np.array([[inp.get(f, 0) for f in INPUT_FEATURES]])
        x_scaled = m['scaler_X'].transform(x)
        pred_scaled, std_scaled = m['gpr'].predict(x_scaled, return_std=True)
        pred = m['scaler_y'].inverse_transform(pred_scaled.reshape(-1, 1)).ravel()[0]
        std_raw = std_scaled[0] * m['scaler_y'].scale_[0]
        if m['use_log']:
            pred = np.exp(pred)
            std_raw = pred * abs(std_raw)

        # Expected Improvement: EI = (pred - best) * Phi(z) + std * phi(z)
        # Simplified: score = pred + 1.5 * std (upper confidence bound)
        score = float(pred) + 1.5 * float(abs(std_raw))

        candidates.append({
            'd_se': round(float(d_se), 2),
            'am_pct': round(float(am_pct), 1),
            'ps_frac': round(float(ps_frac), 1),
            'loading': round(float(loading), 1),
            'predicted_sigma': round(float(pred), 4),
            'uncertainty': round(float(abs(std_raw)), 4),
            'acquisition_score': round(score, 4),
        })

    candidates.sort(key=lambda x: x['acquisition_score'], reverse=True)
    return {'suggestions': candidates[:top_n], 'best_known_sigma': round(best_sigma, 4)}


def export_training_csv(results_folder, archive_folder):
    """Export training data as CSV string."""
    rows = load_training_data(results_folder, archive_folder)
    if not rows:
        return ''

    # Collect all keys except internal ones
    skip = {'_all_fm_keys'}
    keys = []
    for r in rows:
        for k in r:
            if k not in skip and k not in keys:
                keys.append(k)

    import csv
    import io
    output = io.StringIO()
    writer = csv.DictWriter(output, fieldnames=keys, extrasaction='ignore')
    writer.writeheader()
    for r in rows:
        clean = {k: v for k, v in r.items() if k not in skip}
        writer.writerow(clean)
    return output.getvalue()


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
            # Calculate energy density & performance score
            util_1c = pred.get('utilization', {}).get('1.0C', 1.0)
            perf_score = sig * util_1c
            thickness_pred = pred['microstructure'].get('thickness', {})
            T_pred = thickness_pred.get('value', 100) if isinstance(thickness_pred, dict) else 100
            if T_pred <= 0:
                T_pred = params.get('loading', 6) * 20  # rough estimate: 1mAh≈20μm

            phi_am_pred = pred['microstructure'].get('phi_am', {})
            phi_am_val = phi_am_pred.get('value', 0.5) if isinstance(phi_am_pred, dict) else 0.5
            if phi_am_val <= 0:
                phi_am_val = (100 - params['am_pct']) / 100 * 0.65  # rough

            # Energy density: Q(mAh/cm²) = specific_cap × ρ_AM × φ_AM × T(cm)
            Q_per_cm2 = 200 * 4.8 * phi_am_val * T_pred * 1e-4  # mAh/cm²
            E_density = Q_per_cm2 * 3.7 * util_1c  # μWh/cm² → mWh/cm² (keep as mWh)

            # σ_electronic from prediction
            sigma_el_val = pred['conductivity'].get('sigma_electronic', 0)
            if sigma_el_val is None or sigma_el_val == 0:
                sigma_el_val = 0

            r = {
                'd_se': round(float(params['d_se']), 1),
                'am_pct': int(params['am_pct']),
                'ps_frac': round(float(params['ps_frac']), 1),
                'loading': round(float(params['loading']), 1),
                'sigma_ionic': round(sig, 4),
                'sigma_electronic': round(sigma_el_val, 2),
                'sigma_thermal': pred['conductivity'].get('sigma_thermal', 0),
                'util_1C': round(util_1c, 3),
                'perf_score': round(perf_score, 4),
                'energy_mWh_cm2': round(E_density, 2),
                'rate_limiting': pred['rate_limiting'],
            }
            # Add sweep params to result
            for k in sweep_keys:
                if k not in r:
                    r[k] = round(float(params[k]), 1)
            results.append(r)

    # Sort by performance score (σ × util@1C) — realistic optimum
    results.sort(key=lambda x: x.get('perf_score', x['sigma_ionic']), reverse=True)
    return results[:top_n]
