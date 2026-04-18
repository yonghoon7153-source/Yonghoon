"""
Thermal conductivity regression analysis.
Thermal uses ALL contacts (AM-AM, AM-SE, SE-SE) with material-specific k.
k_AM = 4.0e-2 W/(cm·K), k_SE = 0.7e-2 W/(cm·K), k_AM-SE = harmonic mean.
"""
import json, os, sys
import numpy as np
from scipy import stats
from scipy.optimize import curve_fit

WEBAPP = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'webapp')

K_AM = 4.0e-2    # W/(cm·K)
K_SE = 0.7e-2    # W/(cm·K)
K_RATIO = K_AM / K_SE  # ≈ 5.7
K_HARMONIC = 2 * K_AM * K_SE / (K_AM + K_SE)  # AM-SE harmonic mean


def load_thermal_data():
    rows = []
    for base in [os.path.join(WEBAPP, 'results'), os.path.join(WEBAPP, 'archive')]:
        if not os.path.isdir(base):
            continue
        for root, dirs, files in os.walk(base):
            if 'full_metrics.json' not in files:
                continue
            met_path = os.path.join(root, 'full_metrics.json')
            with open(met_path) as f:
                m = json.load(f)

            sigma_th = m.get('thermal_sigma_full_mScm', 0)
            if not sigma_th or sigma_th <= 0:
                continue

            phi_se = m.get('phi_se', 0)
            phi_am = m.get('phi_am', 0)
            porosity = m.get('porosity', 0)
            se_cn = m.get('se_se_cn', 0)
            am_cn = m.get('am_am_cn', 0)
            am_se_cn = m.get('am_se_cn_mean', 0)
            thickness = m.get('thickness_um', 0)

            # Contact areas
            area_se_se = m.get('area_SE_SE_total', 0)
            area_am_am = m.get('am_am_total_area', 0)
            area_am_se = m.get('area_AM_SE_total', m.get('area_AM_P_SE_total', 0) + m.get('area_AM_S_SE_total', 0))

            # Number of contacts
            n_se_se = m.get('area_SE_SE_n', 0)
            n_am_am = m.get('am_am_n_contacts', 0)
            n_am_se = m.get('area_AM_SE_n', m.get('area_AM_P_SE_n', 0) + m.get('area_AM_S_SE_n', 0))

            # Ionic metrics for reference
            sigma_ion = m.get('sigma_full_mScm', 0)
            tau = m.get('tortuosity_recommended', m.get('tortuosity_mean', 0))
            f_perc = m.get('percolation_pct', 0) / 100

            # AM particle size
            d_am = 0
            for key in ['r_AM_P', 'r_AM_S', 'r_AM']:
                r = m.get(key, 0)
                if r and r > d_am:
                    d_am = r
            d_am *= 2

            if phi_se <= 0 or se_cn <= 0:
                continue

            rows.append({
                'name': os.path.basename(root),
                'sigma_th': sigma_th,
                'phi_se': phi_se,
                'phi_am': phi_am,
                'porosity': porosity,
                'se_cn': se_cn,
                'am_cn': am_cn,
                'am_se_cn': am_se_cn,
                'thickness': thickness,
                'd_am': d_am,
                'area_se_se': area_se_se,
                'area_am_am': area_am_am,
                'area_am_se': area_am_se,
                'n_se_se': n_se_se,
                'n_am_am': n_am_am,
                'n_am_se': n_am_se,
                'sigma_ion': sigma_ion,
                'tau': tau,
                'f_perc': f_perc,
                # Derived
                'phi_total': phi_se + phi_am,
                'cn_total': se_cn + am_cn,  # rough total CN
                'k_eff_mix': phi_se * K_SE + phi_am * K_AM,  # linear mixing rule
                'k_eff_harm': K_SE * K_AM / (phi_am * K_SE + phi_se * K_AM) if (phi_am * K_SE + phi_se * K_AM) > 0 else 0,
                # Lichtenecker log mixing
                'k_eff_log': np.exp(phi_se * np.log(K_SE) + phi_am * np.log(K_AM)) if phi_se > 0 and phi_am > 0 else K_SE,
                # Thermal-weighted CN: k_weight per contact type
                'cn_thermal_w': (n_se_se * 1.0 + n_am_am * K_RATIO + n_am_se * 2*K_RATIO/(1+K_RATIO)),
                # Total contacts
                'n_total': n_se_se + n_am_am + n_am_se,
                # AM-SE bridge fraction
                'amse_frac': n_am_se / max(n_se_se + n_am_am + n_am_se, 1),
                # Area ratios
                'area_total': area_se_se + area_am_am + area_am_se,
                'amse_area_frac': area_am_se / max(area_se_se + area_am_am + area_am_se, 1),
            })
    return rows


def thermal_regression():
    rows_all = load_thermal_data()

    # Deduplicate
    seen = set()
    rows = []
    for r in rows_all:
        key = round(r['sigma_th'], 3)
        if key not in seen:
            seen.add(key)
            rows.append(r)

    if len(rows) < 5:
        print(f"Only {len(rows)} unique thermal cases. Need >= 5.")
        return

    print(f"\n{'='*70}")
    print(f"THERMAL CONDUCTIVITY REGRESSION ({len(rows)} unique cases)")
    print(f"{'='*70}")

    sigma_th = np.array([r['sigma_th'] for r in rows])
    phi_se = np.array([r['phi_se'] for r in rows])
    phi_am = np.array([r['phi_am'] for r in rows])
    phi_total = np.array([r['phi_total'] for r in rows])
    se_cn = np.array([r['se_cn'] for r in rows])
    am_cn = np.array([r['am_cn'] for r in rows])
    am_se_cn = np.array([r['am_se_cn'] for r in rows])
    thickness = np.array([r['thickness'] for r in rows])
    d_am = np.array([r['d_am'] for r in rows])
    tau = np.array([r['tau'] for r in rows])
    area_am_se = np.array([r['area_am_se'] for r in rows])
    k_mix = np.array([r['k_eff_mix'] for r in rows])
    n = len(rows)

    log_sigma = np.log(sigma_th)
    ss_tot = np.sum((log_sigma - np.mean(log_sigma))**2)

    cn_thermal_w = np.array([r['cn_thermal_w'] for r in rows])
    n_total = np.array([r['n_total'] for r in rows])
    amse_frac = np.array([r['amse_frac'] for r in rows])
    area_total = np.array([r['area_total'] for r in rows])
    amse_area_frac = np.array([r['amse_area_frac'] for r in rows])
    k_log = np.array([r['k_eff_log'] for r in rows])
    porosity_arr = np.array([r['porosity'] for r in rows])
    solid_frac = 1 - porosity_arr / 100

    K_SE_mScm = K_SE * 1e3  # W/(cm·K) → mW/(cm·K) = mS/cm equiv

    # ── Correlation analysis ──
    print("\n--- Correlation with log(σ_th) ---")
    features = {
        'φ_total': phi_total, 'τ': tau, 'φ_SE': phi_se,
        'A_AM-SE': area_am_se, 'T': thickness,
        'CN_SE': se_cn, 'φ_AM': phi_am, 'k_mix': k_mix,
        'CN_AM': am_cn, 'CN_AM-SE': am_se_cn,
        'CN_th_w': cn_thermal_w, 'N_total': n_total,
        'AMSE_frac': amse_frac, 'A_total': area_total,
        'AMSE_A%': amse_area_frac, 'k_log': k_log,
        'porosity': porosity_arr, 'solid': solid_frac,
    }
    corrs = []
    for name, vals in features.items():
        valid = (vals > 0) & np.isfinite(vals)
        if np.sum(valid) < 5:
            continue
        try:
            r_val = np.corrcoef(np.log(vals[valid]), log_sigma[valid])[0, 1]
            corrs.append((name, r_val))
        except:
            pass
    corrs.sort(key=lambda x: -abs(x[1]))
    for name, r_val in corrs:
        print(f"  {name:12s}: r = {r_val:+.3f}")

    # ── Model fitting ──
    results = []

    def fit_log(X_terms, label):
        """Fit log(σ) = const + Σ b_i × log(x_i), return R²."""
        X = np.column_stack(X_terms + [np.ones(n)])
        b, _, _, _ = np.linalg.lstsq(X, log_sigma, rcond=None)
        pred = X @ b
        r2 = 1 - np.sum((log_sigma - pred)**2) / ss_tot
        return b, r2

    def fit_C_fixed(log_rhs):
        """Fit log(σ) = log(C) + log_rhs, return C and R²."""
        log_C = np.mean(log_sigma - log_rhs)
        pred = log_C + log_rhs
        r2 = 1 - np.sum((log_sigma - pred)**2) / ss_tot
        return np.exp(log_C), r2

    # T1: Bruggeman with k_SE — σ_th = k_SE × φ_total^1.5
    log_rhs = 1.5 * np.log(phi_total) + np.log(K_SE_mScm)
    C, r2 = fit_C_fixed(log_rhs - np.log(K_SE_mScm))
    print(f"\n  T1: σ_th = C × φ_total^1.5,  C={C:.4f}, R²={r2:.4f}")
    results.append(('T1', f'C × φ_total^1.5', C, r2, 1))

    # T2: k_mix × φ_total^1.5 (linear mixing of k_AM, k_SE)
    valid_mix = k_mix > 0
    if np.all(valid_mix):
        log_rhs = 1.5 * np.log(phi_total) + np.log(k_mix * 1e3)
        C, r2 = fit_C_fixed(log_rhs - np.log(k_mix * 1e3))
        print(f"  T2: σ_th = C × k_mix × φ_total^1.5,  C={C:.4f}, R²={r2:.4f}")
        results.append(('T2', f'C × k_mix × φ_total^1.5', C, r2, 1))

    # T3: φ_SE^a × φ_AM^b (free, both phases)
    b3, r2_3 = fit_log([np.log(phi_se), np.log(phi_am)], 'T3')
    print(f"  T3: σ_th = C × φ_SE^{b3[0]:.2f} × φ_AM^{b3[1]:.2f},  R²={r2_3:.4f}")
    results.append(('T3', f'φ_SE^{b3[0]:.2f} × φ_AM^{b3[1]:.2f}', np.exp(b3[2]), r2_3, 3))

    # T4: φ_total^1.5 × CN_SE^2 (ionic analog)
    log_rhs = 1.5 * np.log(phi_total) + 2 * np.log(se_cn)
    C, r2 = fit_C_fixed(log_rhs)
    print(f"  T4: σ_th = C × φ_total^1.5 × CN_SE²,  C={C:.4f}, R²={r2:.4f}")
    results.append(('T4', f'C × φ_total^1.5 × CN_SE²', C, r2, 1))

    # T5: φ_SE^1.5 × CN_SE^2 + φ_AM contribution
    b5, r2_5 = fit_log([np.log(phi_se), np.log(se_cn), np.log(phi_am)], 'T5')
    print(f"  T5: σ_th = C × φ_SE^{b5[0]:.2f} × CN_SE^{b5[1]:.2f} × φ_AM^{b5[2]:.2f},  R²={r2_5:.4f}")
    results.append(('T5', f'φ_SE^{b5[0]:.2f} × CN_SE^{b5[1]:.2f} × φ_AM^{b5[2]:.2f}', np.exp(b5[3]), r2_5, 4))

    # T6: Include AM-SE CN (interfacial heat transfer)
    valid_amse = am_se_cn > 0
    if np.sum(valid_amse) > 10:
        mask = valid_amse
        b6, _, _, _ = np.linalg.lstsq(
            np.column_stack([np.log(phi_se[mask]), np.log(se_cn[mask]),
                           np.log(am_se_cn[mask]), np.ones(np.sum(mask))]),
            log_sigma[mask], rcond=None)
        pred6 = np.column_stack([np.log(phi_se[mask]), np.log(se_cn[mask]),
                                np.log(am_se_cn[mask]), np.ones(np.sum(mask))]) @ b6
        ss6 = np.sum((log_sigma[mask] - np.mean(log_sigma[mask]))**2)
        r2_6 = 1 - np.sum((log_sigma[mask] - pred6)**2) / ss6
        print(f"  T6: σ_th = C × φ_SE^{b6[0]:.2f} × CN_SE^{b6[1]:.2f} × CN_AM-SE^{b6[2]:.2f},  R²={r2_6:.4f}")
        results.append(('T6', f'φ_SE^{b6[0]:.2f} × CN_SE^{b6[1]:.2f} × CN_AM-SE^{b6[2]:.2f}', np.exp(b6[3]), r2_6, 4))

    # T7: ALL CNs — φ_total^1.5 × CN_SE^a × CN_AM^b
    valid_am = am_cn > 0
    if np.sum(valid_am) > 10:
        mask = valid_am
        b7, _, _, _ = np.linalg.lstsq(
            np.column_stack([np.log(phi_total[mask]), np.log(se_cn[mask]),
                           np.log(am_cn[mask]), np.ones(np.sum(mask))]),
            log_sigma[mask], rcond=None)
        pred7 = np.column_stack([np.log(phi_total[mask]), np.log(se_cn[mask]),
                                np.log(am_cn[mask]), np.ones(np.sum(mask))]) @ b7
        ss7 = np.sum((log_sigma[mask] - np.mean(log_sigma[mask]))**2)
        r2_7 = 1 - np.sum((log_sigma[mask] - pred7)**2) / ss7
        print(f"  T7: σ_th = C × φ_total^{b7[0]:.2f} × CN_SE^{b7[1]:.2f} × CN_AM^{b7[2]:.2f},  R²={r2_7:.4f}")
        results.append(('T7', f'φ_total^{b7[0]:.2f} × CN_SE^{b7[1]:.2f} × CN_AM^{b7[2]:.2f}', np.exp(b7[3]), r2_7, 4))

    # T8: φ_total^1.5 × (φ_SE × CN_SE² + k_ratio × φ_AM × CN_AM²) — weighted sum
    # This can't be log-linearized, use curve_fit
    try:
        def model_t8(X, C):
            phi_s, phi_a, cn_s, cn_a = X
            return np.log(C * (phi_s * cn_s**2 + K_RATIO * phi_a * cn_a**2))

        X8 = (phi_se, phi_am, se_cn, am_cn)
        popt8, _ = curve_fit(model_t8, X8, log_sigma, p0=[0.1], maxfev=10000)
        pred8 = model_t8(X8, *popt8)
        r2_8 = 1 - np.sum((log_sigma - pred8)**2) / ss_tot
        print(f"  T8: σ_th = C × (φ_SE×CN_SE² + {K_RATIO:.1f}×φ_AM×CN_AM²),  C={popt8[0]:.4f}, R²={r2_8:.4f}")
        results.append(('T8', f'C × (φ_SE×CN_SE² + k_r×φ_AM×CN_AM²)', popt8[0], r2_8, 1))
    except Exception as e:
        print(f"  T8: FAILED ({e})")

    # T9: Fixed beautiful formula — k_SE × C × φ_total^1.5 × CN_SE² × (1 + (k_AM/k_SE - 1) × φ_AM)
    try:
        def model_t9(X, C):
            phi_t, cn_s, phi_a = X
            k_correction = 1 + (K_RATIO - 1) * phi_a  # k increases with AM content
            return np.log(C * K_SE_mScm * phi_t**1.5 * cn_s**2 * k_correction)

        X9 = (phi_total, se_cn, phi_am)
        popt9, _ = curve_fit(model_t9, X9, log_sigma, p0=[0.02], maxfev=10000)
        pred9 = model_t9(X9, *popt9)
        r2_9 = 1 - np.sum((log_sigma - pred9)**2) / ss_tot
        k_corr_range = f"{1+(K_RATIO-1)*phi_am.min():.2f}~{1+(K_RATIO-1)*phi_am.max():.2f}"
        print(f"  T9: σ_th = {popt9[0]:.4f} × k_SE × φ^1.5 × CN_SE² × (1+{K_RATIO-1:.1f}×φ_AM),  R²={r2_9:.4f}")
        print(f"      k_correction range: {k_corr_range}")
        results.append(('T9', f'C×k_SE×φ^1.5×CN_SE²×(1+k_r×φ_AM)', popt9[0], r2_9, 1))
    except Exception as e:
        print(f"  T9: FAILED ({e})")

    # T10: With T/d correction (unified like electronic)
    valid_d = d_am > 0
    if np.sum(valid_d) > 10:
        T_over_d = thickness / d_am
        try:
            def model_t10(X, C, beta):
                phi_t, cn_s, phi_a, xi = X
                k_correction = 1 + (K_RATIO - 1) * phi_a
                return np.log(C * K_SE_mScm * phi_t**1.5 * cn_s**2 * k_correction) + beta/xi

            X10 = (phi_total[valid_d], se_cn[valid_d], phi_am[valid_d], T_over_d[valid_d])
            popt10, _ = curve_fit(model_t10, X10, log_sigma[valid_d], p0=[0.02, 1.0], maxfev=10000)
            pred10 = model_t10(X10, *popt10)
            ss10 = np.sum((log_sigma[valid_d] - np.mean(log_sigma[valid_d]))**2)
            r2_10 = 1 - np.sum((log_sigma[valid_d] - pred10)**2) / ss10
            print(f"  T10: T9 + exp({popt10[1]:.2f}/(T/d)),  C={popt10[0]:.4f}, R²={r2_10:.4f}")
            results.append(('T10', f'T9 + exp(β/(T/d))', popt10[0], r2_10, 2))
        except Exception as e:
            print(f"  T10: FAILED ({e})")

    # ── τ-based models (Bruggeman structure) ──
    valid_tau = tau > 0
    if np.sum(valid_tau) > 10:
        mask = valid_tau
        tau_m = tau[mask]
        log_sigma_m = log_sigma[mask]
        ss_m = np.sum((log_sigma_m - np.mean(log_sigma_m))**2)

        # T11: k_mix × φ_total / τ² (classic Bruggeman for 2-phase)
        k_mix_m = k_mix[mask] * 1e3  # mS/cm equiv
        log_rhs = np.log(k_mix_m) + np.log(phi_total[mask]) - 2*np.log(tau_m)
        log_C = np.mean(log_sigma_m - log_rhs)
        pred = log_C + log_rhs
        r2 = 1 - np.sum((log_sigma_m - pred)**2) / ss_m
        print(f"\n  T11: σ_th = {np.exp(log_C):.4f} × k_mix × φ_total / τ²,  R²={r2:.4f}")
        results.append(('T11', f'C × k_mix × φ_total / τ²', np.exp(log_C), r2, 1))

        # T12: k_mix × φ_total^1.5 / τ²
        log_rhs = np.log(k_mix_m) + 1.5*np.log(phi_total[mask]) - 2*np.log(tau_m)
        log_C = np.mean(log_sigma_m - log_rhs)
        pred = log_C + log_rhs
        r2 = 1 - np.sum((log_sigma_m - pred)**2) / ss_m
        print(f"  T12: σ_th = {np.exp(log_C):.4f} × k_mix × φ^1.5 / τ²,  R²={r2:.4f}")
        results.append(('T12', f'C × k_mix × φ^1.5 / τ²', np.exp(log_C), r2, 1))

        # T13: φ_total / τ² (pure geometry, no k)
        log_rhs = np.log(phi_total[mask]) - 2*np.log(tau_m)
        log_C = np.mean(log_sigma_m - log_rhs)
        pred = log_C + log_rhs
        r2 = 1 - np.sum((log_sigma_m - pred)**2) / ss_m
        print(f"  T13: σ_th = {np.exp(log_C):.4f} × φ_total / τ²,  R²={r2:.4f}")
        results.append(('T13', f'C × φ_total / τ²', np.exp(log_C), r2, 1))

        # T14: φ_total^a / τ^b (free)
        b14, _, _, _ = np.linalg.lstsq(
            np.column_stack([np.log(phi_total[mask]), np.log(tau_m), np.ones(np.sum(mask))]),
            log_sigma_m, rcond=None)
        pred14 = np.column_stack([np.log(phi_total[mask]), np.log(tau_m), np.ones(np.sum(mask))]) @ b14
        r2_14 = 1 - np.sum((log_sigma_m - pred14)**2) / ss_m
        print(f"  T14: σ_th = C × φ_total^{b14[0]:.2f} / τ^{-b14[1]:.2f},  R²={r2_14:.4f}")
        results.append(('T14', f'φ_total^{b14[0]:.2f} × τ^{b14[1]:.2f}', np.exp(b14[2]), r2_14, 3))

        # T15: φ_SE^a × φ_AM^b / τ^c (free, both phases + τ)
        b15, _, _, _ = np.linalg.lstsq(
            np.column_stack([np.log(phi_se[mask]), np.log(phi_am[mask]), np.log(tau_m), np.ones(np.sum(mask))]),
            log_sigma_m, rcond=None)
        pred15 = np.column_stack([np.log(phi_se[mask]), np.log(phi_am[mask]), np.log(tau_m), np.ones(np.sum(mask))]) @ b15
        r2_15 = 1 - np.sum((log_sigma_m - pred15)**2) / ss_m
        print(f"  T15: σ_th = C × φ_SE^{b15[0]:.2f} × φ_AM^{b15[1]:.2f} × τ^{b15[2]:.2f},  R²={r2_15:.4f}")
        results.append(('T15', f'φ_SE^{b15[0]:.2f} × φ_AM^{b15[1]:.2f} × τ^{b15[2]:.2f}', np.exp(b15[3]), r2_15, 4))

        # T16: k_mix × φ_total / τ² × (1 + α × A_AM-SE) — interface term
        valid_16 = mask & (area_am_se > 0)
        if np.sum(valid_16) > 10:
            try:
                def model_t16(X, C, alpha):
                    k_m, phi_t, tau_v, a_amse = X
                    return np.log(C * k_m * phi_t / tau_v**2 * (1 + alpha * a_amse))

                X16 = (k_mix[valid_16]*1e3, phi_total[valid_16], tau[valid_16], area_am_se[valid_16])
                popt16, _ = curve_fit(model_t16, X16, log_sigma[valid_16], p0=[1.0, 0.001], maxfev=10000)
                pred16 = model_t16(X16, *popt16)
                ss16 = np.sum((log_sigma[valid_16] - np.mean(log_sigma[valid_16]))**2)
                r2_16 = 1 - np.sum((log_sigma[valid_16] - pred16)**2) / ss16
                print(f"  T16: σ_th = {popt16[0]:.4f} × k_mix × φ/τ² × (1+{popt16[1]:.5f}×A_AM-SE),  R²={r2_16:.4f}")
                results.append(('T16', f'k_mix×φ/τ²×(1+α×A_AMSE)', popt16[0], r2_16, 2))
            except Exception as e:
                print(f"  T16: FAILED ({e})")

        # T17: Beautiful candidate — C × (φ_SE × k_SE + φ_AM × k_AM) × φ_total^0.5 / τ²
        # This is: k_eff_mix × Bruggeman(φ^0.5 / τ²)
        log_rhs = np.log(k_mix_m) + 0.5*np.log(phi_total[mask]) - 2*np.log(tau_m)
        log_C = np.mean(log_sigma_m - log_rhs)
        pred = log_C + log_rhs
        r2 = 1 - np.sum((log_sigma_m - pred)**2) / ss_m
        print(f"  T17: σ_th = {np.exp(log_C):.4f} × k_mix × φ^0.5 / τ²,  R²={r2:.4f}")
        results.append(('T17', f'C × k_mix × φ^0.5 / τ²', np.exp(log_C), r2, 1))

    # ── Contact-based thermal models ──

    # T18: Lichtenecker log-mixing × φ^1.5 / τ²
    valid_klog = k_log > 0
    if np.sum(valid_klog & valid_tau) > 10:
        mask = valid_klog & valid_tau
        log_rhs = np.log(k_log[mask]*1e3) + 1.5*np.log(phi_total[mask]) - 2*np.log(tau[mask])
        log_C = np.mean(log_sigma[mask] - log_rhs)
        pred = log_C + log_rhs
        ss_m = np.sum((log_sigma[mask] - np.mean(log_sigma[mask]))**2)
        r2 = 1 - np.sum((log_sigma[mask] - pred)**2) / ss_m
        print(f"\n  T18: σ_th = {np.exp(log_C):.4f} × k_Lichtenecker × φ^1.5 / τ²,  R²={r2:.4f}")
        results.append(('T18', f'C × k_Licht × φ^1.5 / τ²', np.exp(log_C), r2, 1))

    # T19: φ_SE × φ_AM × τ (interaction term: SE-AM bridge)
    if np.sum(valid_tau) > 10:
        mask = valid_tau
        b19, _, _, _ = np.linalg.lstsq(
            np.column_stack([np.log(phi_se[mask]*phi_am[mask]), np.log(tau[mask]), np.ones(np.sum(mask))]),
            log_sigma[mask], rcond=None)
        pred19 = np.column_stack([np.log(phi_se[mask]*phi_am[mask]), np.log(tau[mask]), np.ones(np.sum(mask))]) @ b19
        ss_m = np.sum((log_sigma[mask] - np.mean(log_sigma[mask]))**2)
        r2_19 = 1 - np.sum((log_sigma[mask] - pred19)**2) / ss_m
        print(f"  T19: σ_th = C × (φ_SE×φ_AM)^{b19[0]:.2f} / τ^{-b19[1]:.2f},  R²={r2_19:.4f}")
        results.append(('T19', f'(φ_SE×φ_AM)^{b19[0]:.2f} × τ^{b19[1]:.2f}', np.exp(b19[2]), r2_19, 3))

    # T20: Weighted CN thermal — CN_th_w^a (all contacts, k-weighted)
    valid_cnw = cn_thermal_w > 0
    if np.sum(valid_cnw & valid_tau) > 10:
        mask = valid_cnw & valid_tau
        b20, _, _, _ = np.linalg.lstsq(
            np.column_stack([np.log(phi_total[mask]), np.log(cn_thermal_w[mask]), np.log(tau[mask]), np.ones(np.sum(mask))]),
            log_sigma[mask], rcond=None)
        pred20 = np.column_stack([np.log(phi_total[mask]), np.log(cn_thermal_w[mask]), np.log(tau[mask]), np.ones(np.sum(mask))]) @ b20
        ss_m = np.sum((log_sigma[mask] - np.mean(log_sigma[mask]))**2)
        r2_20 = 1 - np.sum((log_sigma[mask] - pred20)**2) / ss_m
        print(f"  T20: σ_th = C × φ^{b20[0]:.2f} × CN_th_w^{b20[1]:.2f} / τ^{-b20[2]:.2f},  R²={r2_20:.4f}")
        results.append(('T20', f'φ^{b20[0]:.2f} × CN_th_w^{b20[1]:.2f} × τ^{b20[2]:.2f}', np.exp(b20[3]), r2_20, 4))

    # T21: N_total (total contacts) + φ + τ
    valid_nt = n_total > 0
    if np.sum(valid_nt & valid_tau) > 10:
        mask = valid_nt & valid_tau
        b21, _, _, _ = np.linalg.lstsq(
            np.column_stack([np.log(phi_total[mask]), np.log(n_total[mask]), np.log(tau[mask]), np.ones(np.sum(mask))]),
            log_sigma[mask], rcond=None)
        pred21 = np.column_stack([np.log(phi_total[mask]), np.log(n_total[mask]), np.log(tau[mask]), np.ones(np.sum(mask))]) @ b21
        ss_m = np.sum((log_sigma[mask] - np.mean(log_sigma[mask]))**2)
        r2_21 = 1 - np.sum((log_sigma[mask] - pred21)**2) / ss_m
        print(f"  T21: σ_th = C × φ^{b21[0]:.2f} × N_total^{b21[1]:.2f} / τ^{-b21[2]:.2f},  R²={r2_21:.4f}")
        results.append(('T21', f'φ^{b21[0]:.2f} × N_total^{b21[1]:.2f} × τ^{b21[2]:.2f}', np.exp(b21[3]), r2_21, 4))

    # T22: A_total (total contact area) + φ + τ
    valid_at = area_total > 0
    if np.sum(valid_at & valid_tau) > 10:
        mask = valid_at & valid_tau
        b22, _, _, _ = np.linalg.lstsq(
            np.column_stack([np.log(phi_total[mask]), np.log(area_total[mask]), np.log(tau[mask]), np.ones(np.sum(mask))]),
            log_sigma[mask], rcond=None)
        pred22 = np.column_stack([np.log(phi_total[mask]), np.log(area_total[mask]), np.log(tau[mask]), np.ones(np.sum(mask))]) @ b22
        ss_m = np.sum((log_sigma[mask] - np.mean(log_sigma[mask]))**2)
        r2_22 = 1 - np.sum((log_sigma[mask] - pred22)**2) / ss_m
        print(f"  T22: σ_th = C × φ^{b22[0]:.2f} × A_total^{b22[1]:.2f} / τ^{-b22[2]:.2f},  R²={r2_22:.4f}")
        results.append(('T22', f'φ^{b22[0]:.2f} × A_total^{b22[1]:.2f} × τ^{b22[2]:.2f}', np.exp(b22[3]), r2_22, 4))

    # T23: KITCHEN SINK — φ_SE, φ_AM, τ, CN_SE, AM-SE_area
    valid_all = valid_tau & (area_am_se > 0)
    if np.sum(valid_all) > 10:
        mask = valid_all
        b23, _, _, _ = np.linalg.lstsq(
            np.column_stack([np.log(phi_se[mask]), np.log(phi_am[mask]), np.log(tau[mask]),
                           np.log(se_cn[mask]), np.log(area_am_se[mask]), np.ones(np.sum(mask))]),
            log_sigma[mask], rcond=None)
        pred23 = np.column_stack([np.log(phi_se[mask]), np.log(phi_am[mask]), np.log(tau[mask]),
                                np.log(se_cn[mask]), np.log(area_am_se[mask]), np.ones(np.sum(mask))]) @ b23
        ss_m = np.sum((log_sigma[mask] - np.mean(log_sigma[mask]))**2)
        r2_23 = 1 - np.sum((log_sigma[mask] - pred23)**2) / ss_m
        print(f"  T23: φ_SE^{b23[0]:.2f} × φ_AM^{b23[1]:.2f} × τ^{b23[2]:.2f} × CN_SE^{b23[3]:.2f} × A_AMSE^{b23[4]:.2f},  R²={r2_23:.4f}")
        results.append(('T23', f'φ_SE^{b23[0]:.1f}×φ_AM^{b23[1]:.1f}×τ^{b23[2]:.1f}×CN^{b23[3]:.1f}×A^{b23[4]:.1f}', np.exp(b23[5]), r2_23, 6))

    # T24: solid_frac (1-porosity) based — simpler than φ_SE × φ_AM
    if np.sum(valid_tau) > 10:
        mask = valid_tau
        b24, _, _, _ = np.linalg.lstsq(
            np.column_stack([np.log(solid_frac[mask]), np.log(tau[mask]), np.ones(np.sum(mask))]),
            log_sigma[mask], rcond=None)
        pred24 = np.column_stack([np.log(solid_frac[mask]), np.log(tau[mask]), np.ones(np.sum(mask))]) @ b24
        ss_m = np.sum((log_sigma[mask] - np.mean(log_sigma[mask]))**2)
        r2_24 = 1 - np.sum((log_sigma[mask] - pred24)**2) / ss_m
        print(f"  T24: σ_th = C × (1-ε)^{b24[0]:.2f} / τ^{-b24[1]:.2f},  R²={r2_24:.4f}")
        results.append(('T24', f'(1-ε)^{b24[0]:.2f} × τ^{b24[1]:.2f}', np.exp(b24[2]), r2_24, 3))

    # ── Physics-based 2-phase models ──
    sigma_ion = np.array([r['sigma_ion'] for r in rows])
    f_perc = np.array([r['f_perc'] for r in rows])
    sigma_brug_se = np.array([K_SE_mScm * r['phi_se'] * r['f_perc'] / r['tau']**2
                              if r['tau'] > 0 else 0 for r in rows])

    print(f"\n  --- Physics-based 2-phase models ---")

    # T25: σ_th vs σ_ion correlation (same geometry hypothesis)
    valid_ion = sigma_ion > 0
    if np.sum(valid_ion) > 10:
        mask = valid_ion
        r_corr = np.corrcoef(np.log(sigma_ion[mask]), log_sigma[mask])[0, 1]
        b25, _, _, _ = np.linalg.lstsq(
            np.column_stack([np.log(sigma_ion[mask]), np.ones(np.sum(mask))]),
            log_sigma[mask], rcond=None)
        pred25 = np.column_stack([np.log(sigma_ion[mask]), np.ones(np.sum(mask))]) @ b25
        ss25 = np.sum((log_sigma[mask] - np.mean(log_sigma[mask]))**2)
        r2_25 = 1 - np.sum((log_sigma[mask] - pred25)**2) / ss25
        print(f"  T25: σ_th = C × σ_ion^{b25[0]:.2f}  (corr={r_corr:.3f}),  R²={r2_25:.4f}")
        results.append(('T25', f'C × σ_ion^{b25[0]:.2f}', np.exp(b25[1]), r2_25, 2))

    # T26: SE backbone model — k_SE × φ_SE^1.5 / τ² × (1 + (k_r-1) × φ_AM²)
    try:
        def model_t26(X, C, gamma):
            phi_s, phi_a, tau_v = X
            se_backbone = phi_s**1.5 / tau_v**2
            am_enhance = 1 + (K_RATIO - 1) * phi_a**gamma
            return np.log(C * K_SE_mScm * se_backbone * am_enhance)

        mask26 = tau > 0
        X26 = (phi_se[mask26], phi_am[mask26], tau[mask26])
        popt26, _ = curve_fit(model_t26, X26, log_sigma[mask26], p0=[1.0, 2.0], maxfev=10000)
        pred26 = model_t26(X26, *popt26)
        ss26 = np.sum((log_sigma[mask26] - np.mean(log_sigma[mask26]))**2)
        r2_26 = 1 - np.sum((log_sigma[mask26] - pred26)**2) / ss26
        print(f"  T26: σ_th = {popt26[0]:.4f} × k_SE × φ_SE^1.5/τ² × (1+{K_RATIO-1:.1f}×φ_AM^{popt26[1]:.2f}),  R²={r2_26:.4f}")
        results.append(('T26', f'k_SE×φ_SE^1.5/τ²×(1+k_r×φ_AM^γ)', popt26[0], r2_26, 2))
    except Exception as e:
        print(f"  T26: FAILED ({e})")

    # T27: Dual backbone — k_SE×φ_SE^1.5/τ² + k_AM×φ_AM^1.5  (series + parallel)
    try:
        def model_t27(X, C1, C2):
            phi_s, phi_a, tau_v = X
            se_path = C1 * K_SE_mScm * phi_s**1.5 / tau_v**2
            am_path = C2 * K_AM * 1e3 * phi_a**1.5
            return np.log(se_path + am_path)

        mask27 = tau > 0
        X27 = (phi_se[mask27], phi_am[mask27], tau[mask27])
        popt27, _ = curve_fit(model_t27, X27, log_sigma[mask27], p0=[1.0, 0.1], maxfev=10000)
        pred27 = model_t27(X27, *popt27)
        ss27 = np.sum((log_sigma[mask27] - np.mean(log_sigma[mask27]))**2)
        r2_27 = 1 - np.sum((log_sigma[mask27] - pred27)**2) / ss27
        print(f"  T27: σ_th = {popt27[0]:.4f}×k_SE×φ_SE^1.5/τ² + {popt27[1]:.4f}×k_AM×φ_AM^1.5,  R²={r2_27:.4f}")
        results.append(('T27', f'C1×k_SE×φ_SE^1.5/τ² + C2×k_AM×φ_AM^1.5', popt27[0], r2_27, 2))
    except Exception as e:
        print(f"  T27: FAILED ({e})")

    # T28: σ_brug_SE (ionic Bruggeman but with k_SE) as base
    valid_brug = sigma_brug_se > 0
    if np.sum(valid_brug) > 10:
        mask = valid_brug
        # T28a: σ_th = C × σ_brug_SE
        log_C = np.mean(log_sigma[mask] - np.log(sigma_brug_se[mask]))
        pred = log_C + np.log(sigma_brug_se[mask])
        ss28 = np.sum((log_sigma[mask] - np.mean(log_sigma[mask]))**2)
        r2_28 = 1 - np.sum((log_sigma[mask] - pred)**2) / ss28
        print(f"  T28: σ_th = {np.exp(log_C):.4f} × σ_brug_SE,  R²={r2_28:.4f}")
        results.append(('T28', f'C × σ_brug_SE(=k_SE×φ_SE×f_perc/τ²)', np.exp(log_C), r2_28, 1))

        # T28b: + AM correction
        try:
            def model_t28b(X, C, alpha):
                brug_se, phi_a = X
                return np.log(C * brug_se * (1 + alpha * phi_a))

            X28b = (sigma_brug_se[mask], phi_am[mask])
            popt28b, _ = curve_fit(model_t28b, X28b, log_sigma[mask], p0=[5.0, 5.0], maxfev=10000)
            pred28b = model_t28b(X28b, *popt28b)
            r2_28b = 1 - np.sum((log_sigma[mask] - pred28b)**2) / ss28
            print(f"  T28b: σ_th = {popt28b[0]:.4f} × σ_brug_SE × (1+{popt28b[1]:.2f}×φ_AM),  R²={r2_28b:.4f}")
            results.append(('T28b', f'C×σ_brug_SE×(1+α×φ_AM)', popt28b[0], r2_28b, 2))
        except Exception as e:
            print(f"  T28b: FAILED ({e})")

    # T29: σ_th_brug = (φ_SE×k_SE + φ_AM×k_AM) × φ_total^0.5 × f_perc / τ²
    # This is Bruggeman with WEIGHTED k for 2-phase
    try:
        valid_29 = (tau > 0) & (f_perc > 0)
        if np.sum(valid_29) > 10:
            mask = valid_29
            sigma_brug_2phase = k_mix[mask]*1e3 * phi_total[mask]**0.5 * f_perc[mask] / tau[mask]**2
            log_C = np.mean(log_sigma[mask] - np.log(sigma_brug_2phase))
            pred = log_C + np.log(sigma_brug_2phase)
            ss29 = np.sum((log_sigma[mask] - np.mean(log_sigma[mask]))**2)
            r2_29 = 1 - np.sum((log_sigma[mask] - pred)**2) / ss29
            print(f"  T29: σ_th = {np.exp(log_C):.4f} × k_mix × φ^0.5 × f_perc / τ²,  R²={r2_29:.4f}")
            results.append(('T29', f'C × k_mix × φ^0.5 × f_perc / τ²', np.exp(log_C), r2_29, 1))
    except Exception as e:
        print(f"  T29: FAILED ({e})")

    # T30: Probability-weighted k per hop
    # P(SE-SE) ∝ φ_SE², P(AM-AM) ∝ φ_AM², P(AM-SE) ∝ 2φ_SE×φ_AM
    # k_hop = φ_SE²×k_SE + φ_AM²×k_AM + 2×φ_SE×φ_AM×k_harm
    K_HARM = 2 * K_AM * K_SE / (K_AM + K_SE)
    k_hop = phi_se**2 * K_SE + phi_am**2 * K_AM + 2 * phi_se * phi_am * K_HARM
    k_hop_mScm = k_hop * 1e3

    valid_30 = (tau > 0) & (k_hop > 0)
    if np.sum(valid_30) > 10:
        mask = valid_30
        # T30: k_hop × φ_total / τ²
        log_rhs = np.log(k_hop_mScm[mask]) + np.log(phi_total[mask]) - 2*np.log(tau[mask])
        log_C = np.mean(log_sigma[mask] - log_rhs)
        pred = log_C + log_rhs
        ss30 = np.sum((log_sigma[mask] - np.mean(log_sigma[mask]))**2)
        r2_30 = 1 - np.sum((log_sigma[mask] - pred)**2) / ss30
        print(f"  T30: σ_th = {np.exp(log_C):.4f} × k_hop(prob) × φ/τ²,  R²={r2_30:.4f}")
        results.append(('T30', f'C × k_hop(prob-weighted) × φ/τ²', np.exp(log_C), r2_30, 1))

        # T30b: k_hop × φ^1.5 / τ²
        log_rhs = np.log(k_hop_mScm[mask]) + 1.5*np.log(phi_total[mask]) - 2*np.log(tau[mask])
        log_C = np.mean(log_sigma[mask] - log_rhs)
        pred = log_C + log_rhs
        r2_30b = 1 - np.sum((log_sigma[mask] - pred)**2) / ss30
        print(f"  T30b: σ_th = {np.exp(log_C):.4f} × k_hop × φ^1.5/τ²,  R²={r2_30b:.4f}")
        results.append(('T30b', f'C × k_hop × φ^1.5/τ²', np.exp(log_C), r2_30b, 1))

        # T30c: k_hop × φ^a / τ^b (free a, b)
        b30c, _, _, _ = np.linalg.lstsq(
            np.column_stack([np.log(k_hop_mScm[mask]), np.log(phi_total[mask]), np.log(tau[mask]), np.ones(np.sum(mask))]),
            log_sigma[mask], rcond=None)
        pred30c = np.column_stack([np.log(k_hop_mScm[mask]), np.log(phi_total[mask]), np.log(tau[mask]), np.ones(np.sum(mask))]) @ b30c
        r2_30c = 1 - np.sum((log_sigma[mask] - pred30c)**2) / ss30
        print(f"  T30c: σ_th = C × k_hop^{b30c[0]:.2f} × φ^{b30c[1]:.2f} / τ^{-b30c[2]:.2f},  R²={r2_30c:.4f}")
        results.append(('T30c', f'k_hop^{b30c[0]:.2f}×φ^{b30c[1]:.2f}×τ^{b30c[2]:.2f}', np.exp(b30c[3]), r2_30c, 4))

    # T31: THICK ONLY (T > 50μm)
    thick = thickness > 50
    valid_31 = thick & (tau > 0)
    if np.sum(valid_31) > 10:
        mask = valid_31
        log_sigma_m = log_sigma[mask]
        ss31 = np.sum((log_sigma_m - np.mean(log_sigma_m))**2)
        print(f"\n  --- THICK ONLY (T>50μm, n={np.sum(mask)}) ---")

        # T31a: φ_SE × φ_AM × τ
        b31, _, _, _ = np.linalg.lstsq(
            np.column_stack([np.log(phi_se[mask]), np.log(phi_am[mask]), np.log(tau[mask]), np.ones(np.sum(mask))]),
            log_sigma_m, rcond=None)
        pred31 = np.column_stack([np.log(phi_se[mask]), np.log(phi_am[mask]), np.log(tau[mask]), np.ones(np.sum(mask))]) @ b31
        r2_31 = 1 - np.sum((log_sigma_m - pred31)**2) / ss31
        print(f"  T31: σ_th = C × φ_SE^{b31[0]:.2f} × φ_AM^{b31[1]:.2f} / τ^{-b31[2]:.2f},  R²={r2_31:.4f}")
        results.append(('T31', f'[THICK] φ_SE^{b31[0]:.2f}×φ_AM^{b31[1]:.2f}×τ^{b31[2]:.2f}', np.exp(b31[3]), r2_31, 4))

        # T31b: k_hop × φ / τ²
        if np.sum(mask & (k_hop > 0)) > 10:
            mask2 = mask & (k_hop > 0)
            log_rhs = np.log(k_hop_mScm[mask2]) + np.log(phi_total[mask2]) - 2*np.log(tau[mask2])
            log_C = np.mean(log_sigma[mask2] - log_rhs)
            pred = log_C + log_rhs
            ss31b = np.sum((log_sigma[mask2] - np.mean(log_sigma[mask2]))**2)
            r2_31b = 1 - np.sum((log_sigma[mask2] - pred)**2) / ss31b
            print(f"  T31b: σ_th = {np.exp(log_C):.4f} × k_hop × φ/τ² [THICK],  R²={r2_31b:.4f}")
            results.append(('T31b', f'[THICK] C × k_hop × φ/τ²', np.exp(log_C), r2_31b, 1))

    # ══════════════════════════════════════════════════════════════════════
    # NEW PHYSICS MODELS (T32~T40): Fixed exponents + thin electrode correction
    # ══════════════════════════════════════════════════════════════════════
    print(f"\n  {'='*60}")
    print(f"  NEW PHYSICS MODELS (fixed exponents + thin correction)")
    print(f"  {'='*60}")

    # T32: Fixed exponents from T31 insight: φ_SE¹ × φ_AM² / √τ
    # T31 gave φ_SE^1.10 × φ_AM^1.82 / τ^0.51 → round to 1, 2, 0.5
    if np.sum(valid_tau) > 10:
        mask = valid_tau
        log_rhs = 1.0 * np.log(phi_se[mask]) + 2.0 * np.log(phi_am[mask]) - 0.5 * np.log(tau[mask])
        log_C = np.mean(log_sigma[mask] - log_rhs)
        pred = log_C + log_rhs
        ss_m = np.sum((log_sigma[mask] - np.mean(log_sigma[mask]))**2)
        r2_32 = 1 - np.sum((log_sigma[mask] - pred)**2) / ss_m
        print(f"\n  T32: σ_th = {np.exp(log_C):.4f} × φ_SE × φ_AM² / √τ  (0 free + C)")
        print(f"       R²={r2_32:.4f}")
        results.append(('T32', f'C × φ_SE × φ_AM² / √τ (fixed exp)', np.exp(log_C), r2_32, 1))

    # T32b: φ_SE^1.5 × φ_AM^1.5 / √τ (symmetric Bruggeman exponent)
    if np.sum(valid_tau) > 10:
        mask = valid_tau
        log_rhs = 1.5 * np.log(phi_se[mask]) + 1.5 * np.log(phi_am[mask]) - 0.5 * np.log(tau[mask])
        log_C = np.mean(log_sigma[mask] - log_rhs)
        pred = log_C + log_rhs
        ss_m = np.sum((log_sigma[mask] - np.mean(log_sigma[mask]))**2)
        r2 = 1 - np.sum((log_sigma[mask] - pred)**2) / ss_m
        print(f"  T32b: σ_th = {np.exp(log_C):.4f} × φ_SE^1.5 × φ_AM^1.5 / √τ,  R²={r2:.4f}")
        results.append(('T32b', f'C × φ_SE^1.5 × φ_AM^1.5 / √τ', np.exp(log_C), r2, 1))

    # T32c: φ_SE × φ_AM² / τ (τ^1 instead of τ^0.5)
    if np.sum(valid_tau) > 10:
        mask = valid_tau
        log_rhs = 1.0 * np.log(phi_se[mask]) + 2.0 * np.log(phi_am[mask]) - 1.0 * np.log(tau[mask])
        log_C = np.mean(log_sigma[mask] - log_rhs)
        pred = log_C + log_rhs
        ss_m = np.sum((log_sigma[mask] - np.mean(log_sigma[mask]))**2)
        r2 = 1 - np.sum((log_sigma[mask] - pred)**2) / ss_m
        print(f"  T32c: σ_th = {np.exp(log_C):.4f} × φ_SE × φ_AM² / τ,  R²={r2:.4f}")
        results.append(('T32c', f'C × φ_SE × φ_AM² / τ', np.exp(log_C), r2, 1))

    # T33: φ_SE × φ_AM² / √τ × exp(β/(T/d))  — thin electrode correction
    valid_33 = valid_tau & (d_am > 0) & (thickness > 0)
    if np.sum(valid_33) > 10:
        mask = valid_33
        T_over_d = thickness[mask] / d_am[mask]
        try:
            def model_t33(X, C, beta):
                phi_s, phi_a, tau_v, xi = X
                return np.log(C * phi_s * phi_a**2 / np.sqrt(tau_v)) + beta / xi

            X33 = (phi_se[mask], phi_am[mask], tau[mask], T_over_d)
            popt33, _ = curve_fit(model_t33, X33, log_sigma[mask], p0=[50.0, -1.0], maxfev=10000)
            pred33 = model_t33(X33, *popt33)
            ss33 = np.sum((log_sigma[mask] - np.mean(log_sigma[mask]))**2)
            r2_33 = 1 - np.sum((log_sigma[mask] - pred33)**2) / ss33
            print(f"\n  T33: σ_th = {popt33[0]:.4f} × φ_SE × φ_AM² / √τ × exp({popt33[1]:.3f}/(T/d))")
            print(f"       R²={r2_33:.4f}  (T/d range: {T_over_d.min():.1f}~{T_over_d.max():.1f})")
            results.append(('T33', f'φ_SE×φ_AM²/√τ×exp(β/(T/d))', popt33[0], r2_33, 2))
        except Exception as e:
            print(f"  T33: FAILED ({e})")

    # T33b: with π fixed as β (like electronic formula)
    if np.sum(valid_33) > 10:
        mask = valid_33
        T_over_d = thickness[mask] / d_am[mask]
        log_rhs = np.log(phi_se[mask]) + 2*np.log(phi_am[mask]) - 0.5*np.log(tau[mask]) + np.pi / T_over_d
        log_C = np.mean(log_sigma[mask] - log_rhs)
        pred = log_C + log_rhs
        ss33 = np.sum((log_sigma[mask] - np.mean(log_sigma[mask]))**2)
        r2 = 1 - np.sum((log_sigma[mask] - pred)**2) / ss33
        print(f"  T33b: σ_th = {np.exp(log_C):.6f} × φ_SE × φ_AM² / √τ × exp(π/(T/d)),  R²={r2:.4f}")
        results.append(('T33b', f'φ_SE×φ_AM²/√τ×exp(π/(T/d)) (β=π)', np.exp(log_C), r2, 1))

    # T34: σ_th = C × σ_ion^a × k_ratio^b — thermal ∝ ionic geometry
    # More general than T25 — include k_AM/k_SE effect
    valid_34 = (sigma_ion > 0) & valid_tau
    if np.sum(valid_34) > 10:
        mask = valid_34
        # T34a: σ_th = C × σ_ion^a (free a)
        b34, _, _, _ = np.linalg.lstsq(
            np.column_stack([np.log(sigma_ion[mask]), np.ones(np.sum(mask))]),
            log_sigma[mask], rcond=None)
        pred34 = np.column_stack([np.log(sigma_ion[mask]), np.ones(np.sum(mask))]) @ b34
        ss34 = np.sum((log_sigma[mask] - np.mean(log_sigma[mask]))**2)
        r2_34 = 1 - np.sum((log_sigma[mask] - pred34)**2) / ss34
        print(f"\n  T34: σ_th = {np.exp(b34[1]):.4f} × σ_ion^{b34[0]:.3f},  R²={r2_34:.4f}")
        results.append(('T34', f'C × σ_ion^{b34[0]:.2f}', np.exp(b34[1]), r2_34, 2))

        # T34b: σ_th = C × σ_ion^a × φ_AM^b (add AM contribution)
        b34b, _, _, _ = np.linalg.lstsq(
            np.column_stack([np.log(sigma_ion[mask]), np.log(phi_am[mask]), np.ones(np.sum(mask))]),
            log_sigma[mask], rcond=None)
        pred34b = np.column_stack([np.log(sigma_ion[mask]), np.log(phi_am[mask]), np.ones(np.sum(mask))]) @ b34b
        r2_34b = 1 - np.sum((log_sigma[mask] - pred34b)**2) / ss34
        print(f"  T34b: σ_th = C × σ_ion^{b34b[0]:.3f} × φ_AM^{b34b[1]:.3f},  R²={r2_34b:.4f}")
        results.append(('T34b', f'C × σ_ion^{b34b[0]:.2f} × φ_AM^{b34b[1]:.2f}', np.exp(b34b[2]), r2_34b, 3))

    # T35: k_hop × φ_SE / τ — SE backbone with probability-weighted k
    if np.sum(valid_30) > 10:
        mask = valid_30
        log_rhs = np.log(k_hop_mScm[mask]) + np.log(phi_se[mask]) - np.log(tau[mask])
        log_C = np.mean(log_sigma[mask] - log_rhs)
        pred = log_C + log_rhs
        ss35 = np.sum((log_sigma[mask] - np.mean(log_sigma[mask]))**2)
        r2_35 = 1 - np.sum((log_sigma[mask] - pred)**2) / ss35
        print(f"\n  T35: σ_th = {np.exp(log_C):.4f} × k_hop × φ_SE / τ,  R²={r2_35:.4f}")
        results.append(('T35', f'C × k_hop × φ_SE / τ', np.exp(log_C), r2_35, 1))

    # T36: CN_SE^2 term from ionic, adapted for thermal
    # σ_th = C × φ_SE^1.5 × CN_SE^2 × φ_AM^b
    if np.sum(valid_tau) > 10:
        mask = valid_tau
        b36, _, _, _ = np.linalg.lstsq(
            np.column_stack([np.log(phi_se[mask]), np.log(se_cn[mask]), np.log(phi_am[mask]),
                           np.log(tau[mask]), np.ones(np.sum(mask))]),
            log_sigma[mask], rcond=None)
        pred36 = np.column_stack([np.log(phi_se[mask]), np.log(se_cn[mask]), np.log(phi_am[mask]),
                                np.log(tau[mask]), np.ones(np.sum(mask))]) @ b36
        ss36 = np.sum((log_sigma[mask] - np.mean(log_sigma[mask]))**2)
        r2_36 = 1 - np.sum((log_sigma[mask] - pred36)**2) / ss36
        print(f"  T36: σ_th = C × φ_SE^{b36[0]:.2f} × CN_SE^{b36[1]:.2f} × φ_AM^{b36[2]:.2f} / τ^{-b36[3]:.2f}")
        print(f"       R²={r2_36:.4f}")
        results.append(('T36', f'φ_SE^{b36[0]:.1f}×CN^{b36[1]:.1f}×φ_AM^{b36[2]:.1f}×τ^{b36[3]:.1f}', np.exp(b36[4]), r2_36, 5))

    # T37: Hybrid — ionic analog: C × k_SE × φ_SE × f_perc / τ² × (1 + k_r × φ_AM^a)
    # This says thermal = ionic backbone × AM enhancement
    valid_37 = valid_tau & (f_perc > 0)
    if np.sum(valid_37) > 10:
        mask = valid_37
        try:
            def model_t37(X, C, gamma):
                phi_s, fp, tau_v, phi_a = X
                se_backbone = C * K_SE_mScm * phi_s * fp / tau_v**2
                am_enhance = 1 + (K_RATIO - 1) * phi_a**gamma
                return np.log(se_backbone * am_enhance)

            X37 = (phi_se[mask], f_perc[mask], tau[mask], phi_am[mask])
            popt37, _ = curve_fit(model_t37, X37, log_sigma[mask], p0=[1.0, 1.0], maxfev=10000)
            pred37 = model_t37(X37, *popt37)
            ss37 = np.sum((log_sigma[mask] - np.mean(log_sigma[mask]))**2)
            r2_37 = 1 - np.sum((log_sigma[mask] - pred37)**2) / ss37
            print(f"\n  T37: σ_th = {popt37[0]:.4f} × k_SE × φ_SE × f_perc / τ² × (1 + {K_RATIO-1:.1f} × φ_AM^{popt37[1]:.2f})")
            print(f"       R²={r2_37:.4f}  (ionic backbone × AM enhancement)")
            results.append(('T37', f'k_SE×φ_SE×f_perc/τ²×(1+k_r×φ_AM^γ)', popt37[0], r2_37, 2))
        except Exception as e:
            print(f"  T37: FAILED ({e})")

    # T38: Simple power law — φ_total^a / τ^b with CN_SE^c
    # Trying to see if CN_SE helps
    if np.sum(valid_tau) > 10:
        mask = valid_tau
        b38, _, _, _ = np.linalg.lstsq(
            np.column_stack([np.log(phi_total[mask]), np.log(tau[mask]),
                           np.log(se_cn[mask]), np.ones(np.sum(mask))]),
            log_sigma[mask], rcond=None)
        pred38 = np.column_stack([np.log(phi_total[mask]), np.log(tau[mask]),
                                np.log(se_cn[mask]), np.ones(np.sum(mask))]) @ b38
        ss38 = np.sum((log_sigma[mask] - np.mean(log_sigma[mask]))**2)
        r2_38 = 1 - np.sum((log_sigma[mask] - pred38)**2) / ss38
        print(f"  T38: σ_th = C × φ_total^{b38[0]:.2f} / τ^{-b38[1]:.2f} × CN_SE^{b38[2]:.2f},  R²={r2_38:.4f}")
        results.append(('T38', f'φ_total^{b38[0]:.1f}×τ^{b38[1]:.1f}×CN^{b38[2]:.1f}', np.exp(b38[3]), r2_38, 4))

    # T39: (φ_SE × φ_AM)^a × CN_SE^b / τ^c (interaction × connectivity)
    if np.sum(valid_tau) > 10:
        mask = valid_tau
        b39, _, _, _ = np.linalg.lstsq(
            np.column_stack([np.log(phi_se[mask] * phi_am[mask]), np.log(se_cn[mask]),
                           np.log(tau[mask]), np.ones(np.sum(mask))]),
            log_sigma[mask], rcond=None)
        pred39 = np.column_stack([np.log(phi_se[mask] * phi_am[mask]), np.log(se_cn[mask]),
                                np.log(tau[mask]), np.ones(np.sum(mask))]) @ b39
        ss39 = np.sum((log_sigma[mask] - np.mean(log_sigma[mask]))**2)
        r2_39 = 1 - np.sum((log_sigma[mask] - pred39)**2) / ss39
        print(f"  T39: σ_th = C × (φ_SE×φ_AM)^{b39[0]:.2f} × CN_SE^{b39[1]:.2f} / τ^{-b39[2]:.2f},  R²={r2_39:.4f}")
        results.append(('T39', f'(φ_SE×φ_AM)^{b39[0]:.1f}×CN^{b39[1]:.1f}×τ^{b39[2]:.1f}', np.exp(b39[3]), r2_39, 4))

    # T40: UNIFIED thermal = electronic + ionic structure
    # σ_th = C × k_hop × φ_total^1.5 × CN_SE^a / τ^b × exp(β/(T/d))
    valid_40 = valid_tau & (d_am > 0) & (thickness > 0) & (k_hop > 0)
    if np.sum(valid_40) > 10:
        mask = valid_40
        T_over_d = thickness[mask] / d_am[mask]
        try:
            def model_t40(X, C, a, b, beta):
                k_h, phi_t, cn_s, tau_v, xi = X
                return np.log(C * k_h * phi_t**1.5 * cn_s**a / tau_v**b) + beta / xi

            X40 = (k_hop_mScm[mask], phi_total[mask], se_cn[mask], tau[mask], T_over_d)
            popt40, _ = curve_fit(model_t40, X40, log_sigma[mask],
                                  p0=[0.1, 0.5, 1.0, -1.0], maxfev=10000)
            pred40 = model_t40(X40, *popt40)
            ss40 = np.sum((log_sigma[mask] - np.mean(log_sigma[mask]))**2)
            r2_40 = 1 - np.sum((log_sigma[mask] - pred40)**2) / ss40
            print(f"\n  T40: σ_th = {popt40[0]:.4f} × k_hop × φ^1.5 × CN_SE^{popt40[1]:.2f} / τ^{popt40[2]:.2f} × exp({popt40[3]:.3f}/(T/d))")
            print(f"       R²={r2_40:.4f}  (UNIFIED: k_hop + geometry + thin correction)")
            results.append(('T40', f'k_hop×φ^1.5×CN^{popt40[1]:.1f}/τ^{popt40[2]:.1f}×exp(β/(T/d))', popt40[0], r2_40, 4))
        except Exception as e:
            print(f"  T40: FAILED ({e})")

    # T41: THICK ONLY versions of best models
    if np.sum(valid_31) > 10:
        mask = valid_31
        log_sigma_m = log_sigma[mask]
        ss_thick = np.sum((log_sigma_m - np.mean(log_sigma_m))**2)
        print(f"\n  --- THICK ONLY new models (T>50μm, n={np.sum(mask)}) ---")

        # T41a: φ_SE × φ_AM² / √τ [THICK]
        log_rhs = np.log(phi_se[mask]) + 2*np.log(phi_am[mask]) - 0.5*np.log(tau[mask])
        log_C = np.mean(log_sigma_m - log_rhs)
        pred = log_C + log_rhs
        r2 = 1 - np.sum((log_sigma_m - pred)**2) / ss_thick
        print(f"  T41a: σ_th = {np.exp(log_C):.4f} × φ_SE × φ_AM² / √τ [THICK],  R²={r2:.4f}")
        results.append(('T41a', f'[THICK] C × φ_SE × φ_AM² / √τ', np.exp(log_C), r2, 1))

        # T41b: σ_ion^a × φ_AM^b [THICK]
        valid_41b = mask & (sigma_ion > 0)
        if np.sum(valid_41b) > 10:
            m2 = valid_41b
            b41b, _, _, _ = np.linalg.lstsq(
                np.column_stack([np.log(sigma_ion[m2]), np.log(phi_am[m2]), np.ones(np.sum(m2))]),
                log_sigma[m2], rcond=None)
            pred41b = np.column_stack([np.log(sigma_ion[m2]), np.log(phi_am[m2]), np.ones(np.sum(m2))]) @ b41b
            ss41b = np.sum((log_sigma[m2] - np.mean(log_sigma[m2]))**2)
            r2_41b = 1 - np.sum((log_sigma[m2] - pred41b)**2) / ss41b
            print(f"  T41b: σ_th = C × σ_ion^{b41b[0]:.3f} × φ_AM^{b41b[1]:.3f} [THICK],  R²={r2_41b:.4f}")
            results.append(('T41b', f'[THICK] σ_ion^{b41b[0]:.2f}×φ_AM^{b41b[1]:.2f}', np.exp(b41b[2]), r2_41b, 3))

        # T41c: CN_SE^a × φ_SE^b × φ_AM^c [THICK, no τ]
        b41c, _, _, _ = np.linalg.lstsq(
            np.column_stack([np.log(se_cn[mask]), np.log(phi_se[mask]), np.log(phi_am[mask]), np.ones(np.sum(mask))]),
            log_sigma_m, rcond=None)
        pred41c = np.column_stack([np.log(se_cn[mask]), np.log(phi_se[mask]), np.log(phi_am[mask]), np.ones(np.sum(mask))]) @ b41c
        r2_41c = 1 - np.sum((log_sigma_m - pred41c)**2) / ss_thick
        print(f"  T41c: σ_th = C × CN_SE^{b41c[0]:.2f} × φ_SE^{b41c[1]:.2f} × φ_AM^{b41c[2]:.2f} [THICK],  R²={r2_41c:.4f}")
        results.append(('T41c', f'[THICK] CN^{b41c[0]:.1f}×φ_SE^{b41c[1]:.1f}×φ_AM^{b41c[2]:.1f}', np.exp(b41c[3]), r2_41c, 4))

    # ══════════════════════════════════════════════════════════════════════
    # DEEPER PHYSICS: weighted contact area, re-parameterization, residuals
    # ══════════════════════════════════════════════════════════════════════
    print(f"\n  {'='*60}")
    print(f"  DEEPER PHYSICS MODELS (T42~T50)")
    print(f"  {'='*60}")

    area_se_se = np.array([r['area_se_se'] for r in rows])
    area_am_am = np.array([r['area_am_am'] for r in rows])

    # T42: Weighted contact area model (physically correct!)
    # σ_th ∝ (A_SE_SE × k_SE + A_AM_AM × k_AM + A_AM_SE × k_harm) / thickness
    K_HARM_VAL = 2 * K_AM * K_SE / (K_AM + K_SE)
    A_weighted = (area_se_se * K_SE + area_am_am * K_AM + area_am_se * K_HARM_VAL) * 1e3
    valid_aw = (A_weighted > 0) & (thickness > 0)
    if np.sum(valid_aw) > 10:
        mask = valid_aw
        # T42a: σ_th = C × A_weighted / T
        log_rhs = np.log(A_weighted[mask]) - np.log(thickness[mask])
        log_C = np.mean(log_sigma[mask] - log_rhs)
        pred = log_C + log_rhs
        ss42 = np.sum((log_sigma[mask] - np.mean(log_sigma[mask]))**2)
        r2_42 = 1 - np.sum((log_sigma[mask] - pred)**2) / ss42
        print(f"\n  T42: σ_th = {np.exp(log_C):.6f} × A_weighted / T,  R²={r2_42:.4f}")
        print(f"       A_weighted = A_SE×k_SE + A_AM×k_AM + A_AMSE×k_harm")
        results.append(('T42', f'C × A_weighted / T', np.exp(log_C), r2_42, 1))

        # T42b: free exponent
        b42b, _, _, _ = np.linalg.lstsq(
            np.column_stack([np.log(A_weighted[mask]), np.log(thickness[mask]), np.ones(np.sum(mask))]),
            log_sigma[mask], rcond=None)
        pred42b = np.column_stack([np.log(A_weighted[mask]), np.log(thickness[mask]), np.ones(np.sum(mask))]) @ b42b
        r2_42b = 1 - np.sum((log_sigma[mask] - pred42b)**2) / ss42
        print(f"  T42b: σ_th = C × A_weighted^{b42b[0]:.2f} × T^{b42b[1]:.2f},  R²={r2_42b:.4f}")
        results.append(('T42b', f'A_weighted^{b42b[0]:.1f}×T^{b42b[1]:.1f}', np.exp(b42b[2]), r2_42b, 3))

        # T42c: A_weighted only (no T dependence — testing thickness-independence)
        log_rhs_c = np.log(A_weighted[mask])
        log_C_c = np.mean(log_sigma[mask] - log_rhs_c)
        pred_c = log_C_c + log_rhs_c
        r2_42c = 1 - np.sum((log_sigma[mask] - pred_c)**2) / ss42
        print(f"  T42c: σ_th = {np.exp(log_C_c):.6f} × A_weighted (no T!),  R²={r2_42c:.4f}")
        results.append(('T42c', f'C × A_weighted (no T)', np.exp(log_C_c), r2_42c, 1))

    # T43: Re-parameterize with solid fraction + AM ratio
    # solid = φ_SE + φ_AM, x_AM = φ_AM / solid
    x_AM = phi_am / phi_total  # AM fraction among solids
    valid_43 = (x_AM > 0) & (x_AM < 1) & valid_tau
    if np.sum(valid_43) > 10:
        mask = valid_43
        ss43 = np.sum((log_sigma[mask] - np.mean(log_sigma[mask]))**2)

        # T43a: solid^a × x_AM^b / τ^c
        b43, _, _, _ = np.linalg.lstsq(
            np.column_stack([np.log(phi_total[mask]), np.log(x_AM[mask]), np.log(tau[mask]), np.ones(np.sum(mask))]),
            log_sigma[mask], rcond=None)
        pred43 = np.column_stack([np.log(phi_total[mask]), np.log(x_AM[mask]), np.log(tau[mask]), np.ones(np.sum(mask))]) @ b43
        r2_43 = 1 - np.sum((log_sigma[mask] - pred43)**2) / ss43
        print(f"\n  T43: σ_th = C × solid^{b43[0]:.2f} × x_AM^{b43[1]:.2f} / τ^{-b43[2]:.2f},  R²={r2_43:.4f}")
        print(f"       (solid=φ_SE+φ_AM, x_AM=φ_AM/solid)")
        results.append(('T43', f'solid^{b43[0]:.1f}×x_AM^{b43[1]:.1f}×τ^{b43[2]:.1f}', np.exp(b43[3]), r2_43, 4))

        # T43b: solid^a × (1-x_AM)^b / τ^c  (= SE fraction effect)
        b43b, _, _, _ = np.linalg.lstsq(
            np.column_stack([np.log(phi_total[mask]), np.log(1-x_AM[mask]), np.log(tau[mask]), np.ones(np.sum(mask))]),
            log_sigma[mask], rcond=None)
        pred43b = np.column_stack([np.log(phi_total[mask]), np.log(1-x_AM[mask]), np.log(tau[mask]), np.ones(np.sum(mask))]) @ b43b
        r2_43b = 1 - np.sum((log_sigma[mask] - pred43b)**2) / ss43
        print(f"  T43b: σ_th = C × solid^{b43b[0]:.2f} × x_SE^{b43b[1]:.2f} / τ^{-b43b[2]:.2f},  R²={r2_43b:.4f}")
        results.append(('T43b', f'solid^{b43b[0]:.1f}×x_SE^{b43b[1]:.1f}×τ^{b43b[2]:.1f}', np.exp(b43b[3]), r2_43b, 4))

    # T44: σ_th = C × σ_ion^0.5 × φ_AM² (clean hybrid)
    valid_44 = sigma_ion > 0
    if np.sum(valid_44) > 10:
        mask = valid_44
        ss44 = np.sum((log_sigma[mask] - np.mean(log_sigma[mask]))**2)
        log_rhs = 0.5 * np.log(sigma_ion[mask]) + 2.0 * np.log(phi_am[mask])
        log_C = np.mean(log_sigma[mask] - log_rhs)
        pred = log_C + log_rhs
        r2_44 = 1 - np.sum((log_sigma[mask] - pred)**2) / ss44
        print(f"\n  T44: σ_th = {np.exp(log_C):.4f} × σ_ion^0.5 × φ_AM²  (fixed exponents),  R²={r2_44:.4f}")
        results.append(('T44', f'C × σ_ion^0.5 × φ_AM² (fixed)', np.exp(log_C), r2_44, 1))

        # T44b: σ_ion^0.5 × φ_AM^1.5
        log_rhs_b = 0.5 * np.log(sigma_ion[mask]) + 1.5 * np.log(phi_am[mask])
        log_C_b = np.mean(log_sigma[mask] - log_rhs_b)
        pred_b = log_C_b + log_rhs_b
        r2_44b = 1 - np.sum((log_sigma[mask] - pred_b)**2) / ss44
        print(f"  T44b: σ_th = {np.exp(log_C_b):.4f} × σ_ion^0.5 × φ_AM^1.5,  R²={r2_44b:.4f}")
        results.append(('T44b', f'C × σ_ion^0.5 × φ_AM^1.5', np.exp(log_C_b), r2_44b, 1))

        # T44c: σ_ion^a × φ_AM^b × CN_SE^c (add connectivity)
        valid_44c = mask & (se_cn > 0)
        if np.sum(valid_44c) > 10:
            m2 = valid_44c
            b44c, _, _, _ = np.linalg.lstsq(
                np.column_stack([np.log(sigma_ion[m2]), np.log(phi_am[m2]),
                               np.log(se_cn[m2]), np.ones(np.sum(m2))]),
                log_sigma[m2], rcond=None)
            pred44c = np.column_stack([np.log(sigma_ion[m2]), np.log(phi_am[m2]),
                                     np.log(se_cn[m2]), np.ones(np.sum(m2))]) @ b44c
            ss44c = np.sum((log_sigma[m2] - np.mean(log_sigma[m2]))**2)
            r2_44c = 1 - np.sum((log_sigma[m2] - pred44c)**2) / ss44c
            print(f"  T44c: σ_th = C × σ_ion^{b44c[0]:.3f} × φ_AM^{b44c[1]:.3f} × CN_SE^{b44c[2]:.3f},  R²={r2_44c:.4f}")
            results.append(('T44c', f'σ_ion^{b44c[0]:.2f}×φ_AM^{b44c[1]:.2f}×CN^{b44c[2]:.2f}', np.exp(b44c[3]), r2_44c, 4))

    # ── T44c DEEP DIVE: fixed exponents + per-case accuracy ──
    if np.sum(valid_44) > 10 and np.sum(valid_44 & (se_cn > 0)) > 10:
        mask = valid_44 & (se_cn > 0)
        ss_dd = np.sum((log_sigma[mask] - np.mean(log_sigma[mask]))**2)

        # T44d: σ_ion^(3/4) × φ_AM² / CN_SE (fixed beautiful exponents)
        log_rhs_d = 0.75 * np.log(sigma_ion[mask]) + 2.0 * np.log(phi_am[mask]) - 1.0 * np.log(se_cn[mask])
        log_C_d = np.mean(log_sigma[mask] - log_rhs_d)
        pred_d = log_C_d + log_rhs_d
        r2_d = 1 - np.sum((log_sigma[mask] - pred_d)**2) / ss_dd
        print(f"\n  T44d: σ_th = {np.exp(log_C_d):.4f} × σ_ion^(3/4) × φ_AM² / CN_SE  (FIXED),  R²={r2_d:.4f}")
        results.append(('T44d', f'C × σ_ion^3/4 × φ_AM² / CN_SE (fixed)', np.exp(log_C_d), r2_d, 1))

        # T44e: σ_ion^(1/2) × φ_AM² / CN_SE
        log_rhs_e = 0.5 * np.log(sigma_ion[mask]) + 2.0 * np.log(phi_am[mask]) - 1.0 * np.log(se_cn[mask])
        log_C_e = np.mean(log_sigma[mask] - log_rhs_e)
        pred_e = log_C_e + log_rhs_e
        r2_e = 1 - np.sum((log_sigma[mask] - pred_e)**2) / ss_dd
        print(f"  T44e: σ_th = {np.exp(log_C_e):.4f} × σ_ion^(1/2) × φ_AM² / CN_SE,  R²={r2_e:.4f}")
        results.append(('T44e', f'C × σ_ion^1/2 × φ_AM² / CN_SE', np.exp(log_C_e), r2_e, 1))

        # T44f: σ_ion × φ_AM² / CN_SE (σ_ion^1)
        log_rhs_f = 1.0 * np.log(sigma_ion[mask]) + 2.0 * np.log(phi_am[mask]) - 1.0 * np.log(se_cn[mask])
        log_C_f = np.mean(log_sigma[mask] - log_rhs_f)
        pred_f = log_C_f + log_rhs_f
        r2_f = 1 - np.sum((log_sigma[mask] - pred_f)**2) / ss_dd
        print(f"  T44f: σ_th = {np.exp(log_C_f):.4f} × σ_ion × φ_AM² / CN_SE,  R²={r2_f:.4f}")
        results.append(('T44f', f'C × σ_ion × φ_AM² / CN_SE', np.exp(log_C_f), r2_f, 1))

        # T44g: σ_ion^(3/4) × φ_AM^(3/2) / CN_SE
        log_rhs_g = 0.75 * np.log(sigma_ion[mask]) + 1.5 * np.log(phi_am[mask]) - 1.0 * np.log(se_cn[mask])
        log_C_g = np.mean(log_sigma[mask] - log_rhs_g)
        pred_g = log_C_g + log_rhs_g
        r2_g = 1 - np.sum((log_sigma[mask] - pred_g)**2) / ss_dd
        print(f"  T44g: σ_th = {np.exp(log_C_g):.4f} × σ_ion^(3/4) × φ_AM^(3/2) / CN_SE,  R²={r2_g:.4f}")
        results.append(('T44g', f'C × σ_ion^3/4 × φ_AM^3/2 / CN_SE', np.exp(log_C_g), r2_g, 1))

        # T44h: σ_ion^a × φ_AM^b / CN_SE^c + f_perc (residual showed f_perc r=-0.27)
        valid_44h = mask & (f_perc > 0)
        if np.sum(valid_44h) > 10:
            m_h = valid_44h
            b44h, _, _, _ = np.linalg.lstsq(
                np.column_stack([np.log(sigma_ion[m_h]), np.log(phi_am[m_h]),
                               np.log(se_cn[m_h]), np.log(f_perc[m_h]), np.ones(np.sum(m_h))]),
                log_sigma[m_h], rcond=None)
            pred44h = np.column_stack([np.log(sigma_ion[m_h]), np.log(phi_am[m_h]),
                                     np.log(se_cn[m_h]), np.log(f_perc[m_h]), np.ones(np.sum(m_h))]) @ b44h
            ss44h = np.sum((log_sigma[m_h] - np.mean(log_sigma[m_h]))**2)
            r2_44h = 1 - np.sum((log_sigma[m_h] - pred44h)**2) / ss44h
            print(f"  T44h: σ_ion^{b44h[0]:.3f} × φ_AM^{b44h[1]:.3f} × CN^{b44h[2]:.3f} × f_perc^{b44h[3]:.3f},  R²={r2_44h:.4f}")
            results.append(('T44h', f'σ_ion^{b44h[0]:.2f}×φ_AM^{b44h[1]:.2f}×CN^{b44h[2]:.2f}×f_perc^{b44h[3]:.2f}', np.exp(b44h[4]), r2_44h, 5))

        # T44i: Expand σ_ion = σ_SE × φ_SE × f_perc / τ² → use components directly
        # σ_th = C × φ_SE^a × f_perc^b / τ^c × φ_AM^d / CN_SE^e
        valid_44i = mask & (f_perc > 0) & valid_tau
        if np.sum(valid_44i) > 10:
            m_i = valid_44i
            b44i, _, _, _ = np.linalg.lstsq(
                np.column_stack([np.log(phi_se[m_i]), np.log(f_perc[m_i]), np.log(tau[m_i]),
                               np.log(phi_am[m_i]), np.log(se_cn[m_i]), np.ones(np.sum(m_i))]),
                log_sigma[m_i], rcond=None)
            pred44i = np.column_stack([np.log(phi_se[m_i]), np.log(f_perc[m_i]), np.log(tau[m_i]),
                                     np.log(phi_am[m_i]), np.log(se_cn[m_i]), np.ones(np.sum(m_i))]) @ b44i
            ss44i = np.sum((log_sigma[m_i] - np.mean(log_sigma[m_i]))**2)
            r2_44i = 1 - np.sum((log_sigma[m_i] - pred44i)**2) / ss44i
            print(f"  T44i: φ_SE^{b44i[0]:.2f} × f_perc^{b44i[1]:.2f} / τ^{-b44i[2]:.2f} × φ_AM^{b44i[3]:.2f} / CN_SE^{-b44i[4]:.2f}")
            print(f"        R²={r2_44i:.4f}  (σ_ion components expanded)")
            results.append(('T44i', f'φ_SE^{b44i[0]:.1f}×f_perc^{b44i[1]:.1f}×τ^{b44i[2]:.1f}×φ_AM^{b44i[3]:.1f}×CN^{b44i[4]:.1f}', np.exp(b44i[5]), r2_44i, 6))

        # ── PER-CASE ACCURACY TABLE for T44c (champion) ──
        print(f"\n  {'='*70}")
        print(f"  T44c PER-CASE ACCURACY (σ_ion^0.76 × φ_AM^2 / CN_SE^1.17)")
        print(f"  {'='*70}")
        # Recompute T44c prediction
        b_champ = np.linalg.lstsq(
            np.column_stack([np.log(sigma_ion[mask]), np.log(phi_am[mask]),
                           np.log(se_cn[mask]), np.ones(np.sum(mask))]),
            log_sigma[mask], rcond=None)[0]
        pred_champ = np.column_stack([np.log(sigma_ion[mask]), np.log(phi_am[mask]),
                                    np.log(se_cn[mask]), np.ones(np.sum(mask))]) @ b_champ
        sigma_pred = np.exp(pred_champ)
        sigma_actual = sigma_th[mask]
        errors = (sigma_pred - sigma_actual) / sigma_actual * 100

        print(f"  {'Name':25s} {'σ_actual':>8s} {'σ_pred':>8s} {'error%':>8s} {'σ_ion':>7s} {'φ_AM':>6s} {'CN_SE':>6s}")
        print(f"  {'-'*68}")
        mask_idx = np.where(mask)[0]
        for j, idx in enumerate(mask_idx):
            print(f"  {rows[idx]['name'][:23]:23s} {sigma_actual[j]:8.3f} {sigma_pred[j]:8.3f} {errors[j]:+8.1f}% {sigma_ion[idx]:7.4f} {phi_am[idx]:6.3f} {se_cn[idx]:6.2f}")

        print(f"\n  Mean |error|: {np.mean(np.abs(errors)):.1f}%")
        print(f"  Max  |error|: {np.max(np.abs(errors)):.1f}%")
        print(f"  Cases within 10%: {np.sum(np.abs(errors) < 10)}/{len(errors)}")
        print(f"  Cases within 20%: {np.sum(np.abs(errors) < 20)}/{len(errors)}")

        # ── LOOCV for T44c ──
        print(f"\n  --- Leave-One-Out Cross-Validation for T44c ---")
        loocv_errors = []
        for j in range(np.sum(mask)):
            # Remove point j
            train_mask = np.ones(np.sum(mask), dtype=bool)
            train_mask[j] = False
            X_train = np.column_stack([np.log(sigma_ion[mask])[train_mask],
                                      np.log(phi_am[mask])[train_mask],
                                      np.log(se_cn[mask])[train_mask],
                                      np.ones(np.sum(train_mask))])
            y_train = log_sigma[mask][train_mask]
            b_cv, _, _, _ = np.linalg.lstsq(X_train, y_train, rcond=None)
            X_test = np.array([np.log(sigma_ion[mask])[j], np.log(phi_am[mask])[j],
                              np.log(se_cn[mask])[j], 1.0])
            pred_cv = X_test @ b_cv
            err_cv = (np.exp(pred_cv) - sigma_th[mask][j]) / sigma_th[mask][j] * 100
            loocv_errors.append(err_cv)
        loocv_errors = np.array(loocv_errors)
        ss_loocv = np.sum((sigma_th[mask] - np.exp(np.array([
            np.array([np.log(sigma_ion[mask])[j], np.log(phi_am[mask])[j],
                      np.log(se_cn[mask])[j], 1.0]) @
            np.linalg.lstsq(
                np.column_stack([np.log(sigma_ion[mask])[np.arange(np.sum(mask))!=j],
                                np.log(phi_am[mask])[np.arange(np.sum(mask))!=j],
                                np.log(se_cn[mask])[np.arange(np.sum(mask))!=j],
                                np.ones(np.sum(mask)-1)]),
                log_sigma[mask][np.arange(np.sum(mask))!=j], rcond=None)[0]
            for j in range(np.sum(mask))
        ])))**2)
        r2_loocv = 1 - ss_loocv / np.sum((sigma_th[mask] - np.mean(sigma_th[mask]))**2)
        print(f"  LOOCV R² = {r2_loocv:.4f} (train R² = {r2_44c:.4f})")
        print(f"  LOOCV Mean |error|: {np.mean(np.abs(loocv_errors)):.1f}%")
        print(f"  LOOCV Max  |error|: {np.max(np.abs(loocv_errors)):.1f}%")

    # ══════════════════════════════════════════════════════════════════════
    # T44 REFINEMENTS: push R² above 0.90
    # ══════════════════════════════════════════════════════════════════════
    print(f"\n  {'='*60}")
    print(f"  T44 REFINEMENTS — pushing R² above 0.90")
    print(f"  {'='*60}")

    if np.sum(valid_44) > 10 and np.sum(valid_44 & (se_cn > 0)) > 10:
        mask = valid_44 & (se_cn > 0)

        # R1: OUTLIER ANALYSIS — remove extreme low σ_th cases
        # Identify cases where σ_ion < 0.03 (SE backbone barely exists)
        low_ion = sigma_ion[mask] < 0.03
        n_low = np.sum(low_ion)
        mask_clean = mask.copy()
        # Remove the 5 lowest σ_ion cases
        ion_sorted = np.argsort(sigma_ion)
        low_5 = set(ion_sorted[:5])
        mask_r1 = np.array([i not in low_5 for i in range(n)]) & mask
        if np.sum(mask_r1) > 10:
            b_r1, _, _, _ = np.linalg.lstsq(
                np.column_stack([np.log(sigma_ion[mask_r1]), np.log(phi_am[mask_r1]),
                               np.log(se_cn[mask_r1]), np.ones(np.sum(mask_r1))]),
                log_sigma[mask_r1], rcond=None)
            pred_r1 = np.column_stack([np.log(sigma_ion[mask_r1]), np.log(phi_am[mask_r1]),
                                     np.log(se_cn[mask_r1]), np.ones(np.sum(mask_r1))]) @ b_r1
            ss_r1 = np.sum((log_sigma[mask_r1] - np.mean(log_sigma[mask_r1]))**2)
            r2_r1 = 1 - np.sum((log_sigma[mask_r1] - pred_r1)**2) / ss_r1
            print(f"\n  R1: T44c without 5 lowest σ_ion cases (n={np.sum(mask_r1)})")
            print(f"      σ_ion^{b_r1[0]:.3f} × φ_AM^{b_r1[1]:.3f} × CN^{b_r1[2]:.3f},  R²={r2_r1:.4f}")
            # Fixed exponents on clean data
            log_rhs_r1 = 0.75*np.log(sigma_ion[mask_r1]) + 2.0*np.log(phi_am[mask_r1]) - 1.0*np.log(se_cn[mask_r1])
            log_C_r1 = np.mean(log_sigma[mask_r1] - log_rhs_r1)
            pred_r1f = log_C_r1 + log_rhs_r1
            r2_r1f = 1 - np.sum((log_sigma[mask_r1] - pred_r1f)**2) / ss_r1
            print(f"  R1b: Fixed 3/4, 2, -1 on clean data:  C={np.exp(log_C_r1):.1f}, R²={r2_r1f:.4f}")
            results.append(('R1b', f'[CLEAN n={np.sum(mask_r1)}] σ_ion^3/4×φ_AM²/CN_SE', np.exp(log_C_r1), r2_r1f, 1))

        # R2: Analytic Bruggeman instead of σ_ion
        # σ_brug_analytic = φ_SE × f_perc / τ² (no network solver needed!)
        valid_r2 = mask & (f_perc > 0) & (tau > 0)
        if np.sum(valid_r2) > 10:
            m_r2 = valid_r2
            sigma_brug = phi_se[m_r2] * f_perc[m_r2] / tau[m_r2]**2  # normalized, no σ_SE
            b_r2, _, _, _ = np.linalg.lstsq(
                np.column_stack([np.log(sigma_brug), np.log(phi_am[m_r2]),
                               np.log(se_cn[m_r2]), np.ones(np.sum(m_r2))]),
                log_sigma[m_r2], rcond=None)
            pred_r2 = np.column_stack([np.log(sigma_brug), np.log(phi_am[m_r2]),
                                     np.log(se_cn[m_r2]), np.ones(np.sum(m_r2))]) @ b_r2
            ss_r2 = np.sum((log_sigma[m_r2] - np.mean(log_sigma[m_r2]))**2)
            r2_r2 = 1 - np.sum((log_sigma[m_r2] - pred_r2)**2) / ss_r2
            print(f"\n  R2: σ_brug^{b_r2[0]:.3f} × φ_AM^{b_r2[1]:.3f} × CN^{b_r2[2]:.3f}  (analytic Bruggeman)")
            print(f"      R²={r2_r2:.4f}  (σ_brug = φ_SE × f_perc / τ², no network solver)")
            results.append(('R2', f'σ_brug^{b_r2[0]:.2f}×φ_AM^{b_r2[1]:.2f}×CN^{b_r2[2]:.2f} (analytic)', 1.0, r2_r2, 4))

            # R2b: fixed exponents
            log_rhs_r2 = 0.75*np.log(sigma_brug) + 2.0*np.log(phi_am[m_r2]) - 1.0*np.log(se_cn[m_r2])
            log_C_r2 = np.mean(log_sigma[m_r2] - log_rhs_r2)
            pred_r2f = log_C_r2 + log_rhs_r2
            r2_r2f = 1 - np.sum((log_sigma[m_r2] - pred_r2f)**2) / ss_r2
            print(f"  R2b: Fixed (3/4, 2, -1): C={np.exp(log_C_r2):.1f}, R²={r2_r2f:.4f}")
            results.append(('R2b', f'C×σ_brug^3/4×φ_AM²/CN_SE (analytic, fixed)', np.exp(log_C_r2), r2_r2f, 1))

        # R3: AM-SE CN as additional variable
        valid_r3 = mask & (am_se_cn > 0)
        if np.sum(valid_r3) > 10:
            m_r3 = valid_r3
            b_r3, _, _, _ = np.linalg.lstsq(
                np.column_stack([np.log(sigma_ion[m_r3]), np.log(phi_am[m_r3]),
                               np.log(se_cn[m_r3]), np.log(am_se_cn[m_r3]), np.ones(np.sum(m_r3))]),
                log_sigma[m_r3], rcond=None)
            pred_r3 = np.column_stack([np.log(sigma_ion[m_r3]), np.log(phi_am[m_r3]),
                                     np.log(se_cn[m_r3]), np.log(am_se_cn[m_r3]), np.ones(np.sum(m_r3))]) @ b_r3
            ss_r3 = np.sum((log_sigma[m_r3] - np.mean(log_sigma[m_r3]))**2)
            r2_r3 = 1 - np.sum((log_sigma[m_r3] - pred_r3)**2) / ss_r3
            print(f"\n  R3: σ_ion^{b_r3[0]:.3f} × φ_AM^{b_r3[1]:.3f} × CN_SE^{b_r3[2]:.3f} × CN_AMSE^{b_r3[3]:.3f}")
            print(f"      R²={r2_r3:.4f}  (+AM-SE CN)")
            results.append(('R3', f'σ_ion^{b_r3[0]:.2f}×φ_AM^{b_r3[1]:.2f}×CN_SE^{b_r3[2]:.2f}×CN_AMSE^{b_r3[3]:.2f}', np.exp(b_r3[4]), r2_r3, 5))

        # R4: AM-AM CN instead of or in addition to CN_SE
        valid_r4 = mask & (am_cn > 0)
        if np.sum(valid_r4) > 10:
            m_r4 = valid_r4
            b_r4, _, _, _ = np.linalg.lstsq(
                np.column_stack([np.log(sigma_ion[m_r4]), np.log(phi_am[m_r4]),
                               np.log(se_cn[m_r4]), np.log(am_cn[m_r4]), np.ones(np.sum(m_r4))]),
                log_sigma[m_r4], rcond=None)
            pred_r4 = np.column_stack([np.log(sigma_ion[m_r4]), np.log(phi_am[m_r4]),
                                     np.log(se_cn[m_r4]), np.log(am_cn[m_r4]), np.ones(np.sum(m_r4))]) @ b_r4
            ss_r4 = np.sum((log_sigma[m_r4] - np.mean(log_sigma[m_r4]))**2)
            r2_r4 = 1 - np.sum((log_sigma[m_r4] - pred_r4)**2) / ss_r4
            print(f"  R4: σ_ion^{b_r4[0]:.3f} × φ_AM^{b_r4[1]:.3f} × CN_SE^{b_r4[2]:.3f} × CN_AM^{b_r4[3]:.3f}")
            print(f"      R²={r2_r4:.4f}  (+AM-AM CN)")
            results.append(('R4', f'σ_ion^{b_r4[0]:.2f}×φ_AM^{b_r4[1]:.2f}×CN_SE^{b_r4[2]:.2f}×CN_AM^{b_r4[3]:.2f}', np.exp(b_r4[4]), r2_r4, 5))

        # R5: Replace CN_SE with CN_ratio = CN_SE / CN_AM
        valid_r5 = mask & (am_cn > 0)
        if np.sum(valid_r5) > 10:
            m_r5 = valid_r5
            cn_ratio = se_cn[m_r5] / am_cn[m_r5]
            b_r5, _, _, _ = np.linalg.lstsq(
                np.column_stack([np.log(sigma_ion[m_r5]), np.log(phi_am[m_r5]),
                               np.log(cn_ratio), np.ones(np.sum(m_r5))]),
                log_sigma[m_r5], rcond=None)
            pred_r5 = np.column_stack([np.log(sigma_ion[m_r5]), np.log(phi_am[m_r5]),
                                     np.log(cn_ratio), np.ones(np.sum(m_r5))]) @ b_r5
            ss_r5 = np.sum((log_sigma[m_r5] - np.mean(log_sigma[m_r5]))**2)
            r2_r5 = 1 - np.sum((log_sigma[m_r5] - pred_r5)**2) / ss_r5
            print(f"  R5: σ_ion^{b_r5[0]:.3f} × φ_AM^{b_r5[1]:.3f} × (CN_SE/CN_AM)^{b_r5[2]:.3f}")
            print(f"      R²={r2_r5:.4f}  (CN ratio)")
            results.append(('R5', f'σ_ion^{b_r5[0]:.2f}×φ_AM^{b_r5[1]:.2f}×CN_ratio^{b_r5[2]:.2f}', np.exp(b_r5[3]), r2_r5, 4))

        # R6: Nonlinear — σ_th = C × σ_ion^a × φ_AM^b / (1 + d×CN_SE)
        try:
            def model_r6(X, C, a, b, d):
                s_ion, phi_a, cn_s = X
                return np.log(C * s_ion**a * phi_a**b / (1 + d * cn_s))

            X_r6 = (sigma_ion[mask], phi_am[mask], se_cn[mask])
            popt_r6, _ = curve_fit(model_r6, X_r6, log_sigma[mask],
                                   p0=[300, 0.75, 2.0, 0.1], maxfev=20000)
            pred_r6 = model_r6(X_r6, *popt_r6)
            ss_r6 = np.sum((log_sigma[mask] - np.mean(log_sigma[mask]))**2)
            r2_r6 = 1 - np.sum((log_sigma[mask] - pred_r6)**2) / ss_r6
            print(f"\n  R6: σ_th = {popt_r6[0]:.1f} × σ_ion^{popt_r6[1]:.3f} × φ_AM^{popt_r6[2]:.3f} / (1+{popt_r6[3]:.3f}×CN)")
            print(f"      R²={r2_r6:.4f}  (nonlinear CN term)")
            results.append(('R6', f'σ_ion^{popt_r6[1]:.2f}×φ_AM^{popt_r6[2]:.2f}/(1+d×CN)', popt_r6[0], r2_r6, 4))
        except Exception as e:
            print(f"  R6: FAILED ({e})")

        # R7: Contact area weighted by material: A_AMSE × k_harm + A_SESE × k_SE
        valid_r7 = mask & (area_am_se > 0) & (area_se_se > 0)
        if np.sum(valid_r7) > 10:
            m_r7 = valid_r7
            A_thermal = area_am_se[m_r7] * K_HARM_VAL + area_se_se[m_r7] * K_SE
            b_r7, _, _, _ = np.linalg.lstsq(
                np.column_stack([np.log(sigma_ion[m_r7]), np.log(phi_am[m_r7]),
                               np.log(A_thermal), np.ones(np.sum(m_r7))]),
                log_sigma[m_r7], rcond=None)
            pred_r7 = np.column_stack([np.log(sigma_ion[m_r7]), np.log(phi_am[m_r7]),
                                     np.log(A_thermal), np.ones(np.sum(m_r7))]) @ b_r7
            ss_r7 = np.sum((log_sigma[m_r7] - np.mean(log_sigma[m_r7]))**2)
            r2_r7 = 1 - np.sum((log_sigma[m_r7] - pred_r7)**2) / ss_r7
            print(f"\n  R7: σ_ion^{b_r7[0]:.3f} × φ_AM^{b_r7[1]:.3f} × A_thermal^{b_r7[2]:.3f}")
            print(f"      R²={r2_r7:.4f}  (A_thermal = A_AMSE×k_harm + A_SESE×k_SE)")
            results.append(('R7', f'σ_ion^{b_r7[0]:.2f}×φ_AM^{b_r7[1]:.2f}×A_th^{b_r7[2]:.2f}', np.exp(b_r7[3]), r2_r7, 4))

        # R8: LOOCV for T44d (fixed, 1p) — should be better than T44c LOOCV
        print(f"\n  --- LOOCV for T44d (fixed 3/4, 2, -1, 1p) ---")
        loocv_errors_d = []
        mask_idx = np.where(mask)[0]
        for j in range(np.sum(mask)):
            train_mask = np.ones(np.sum(mask), dtype=bool)
            train_mask[j] = False
            log_rhs_train = 0.75*np.log(sigma_ion[mask])[train_mask] + 2.0*np.log(phi_am[mask])[train_mask] - 1.0*np.log(se_cn[mask])[train_mask]
            log_C_cv = np.mean(log_sigma[mask][train_mask] - log_rhs_train)
            log_rhs_test = 0.75*np.log(sigma_ion[mask])[j] + 2.0*np.log(phi_am[mask])[j] - 1.0*np.log(se_cn[mask])[j]
            pred_cv = np.exp(log_C_cv + log_rhs_test)
            err = (pred_cv - sigma_th[mask][j]) / sigma_th[mask][j] * 100
            loocv_errors_d.append(err)
        loocv_errors_d = np.array(loocv_errors_d)
        # LOOCV R² for T44d
        pred_loocv_d = np.array([
            np.exp(np.mean(log_sigma[mask][np.arange(np.sum(mask))!=j] -
                          (0.75*np.log(sigma_ion[mask])[np.arange(np.sum(mask))!=j] +
                           2.0*np.log(phi_am[mask])[np.arange(np.sum(mask))!=j] -
                           1.0*np.log(se_cn[mask])[np.arange(np.sum(mask))!=j])) +
                  0.75*np.log(sigma_ion[mask])[j] + 2.0*np.log(phi_am[mask])[j] - 1.0*np.log(se_cn[mask])[j])
            for j in range(np.sum(mask))
        ])
        ss_loocv_d = np.sum((sigma_th[mask] - pred_loocv_d)**2)
        r2_loocv_d = 1 - ss_loocv_d / np.sum((sigma_th[mask] - np.mean(sigma_th[mask]))**2)
        print(f"  T44d LOOCV R² = {r2_loocv_d:.4f} (train R² = 0.8961)")
        print(f"  T44d LOOCV Mean |error|: {np.mean(np.abs(loocv_errors_d)):.1f}%")
        print(f"  T44d LOOCV Max  |error|: {np.max(np.abs(loocv_errors_d)):.1f}%")
        print(f"  T44d LOOCV within 10%: {np.sum(np.abs(loocv_errors_d)<10)}/{len(loocv_errors_d)}")
        print(f"  T44d LOOCV within 20%: {np.sum(np.abs(loocv_errors_d)<20)}/{len(loocv_errors_d)}")

    # T45: Additive model — C1 × φ_SE^a + C2 × φ_AM^b (can't log-linearize)
    try:
        def model_t45(X, C1, a, C2, b):
            phi_s, phi_a = X
            return np.log(C1 * phi_s**a + C2 * phi_a**b)

        popt45, _ = curve_fit(model_t45, (phi_se, phi_am), log_sigma,
                              p0=[10, 1.5, 5, 1.5], maxfev=20000)
        pred45 = model_t45((phi_se, phi_am), *popt45)
        r2_45 = 1 - np.sum((log_sigma - pred45)**2) / ss_tot
        print(f"\n  T45: σ_th = {popt45[0]:.4f}×φ_SE^{popt45[1]:.2f} + {popt45[2]:.4f}×φ_AM^{popt45[3]:.2f}")
        print(f"       R²={r2_45:.4f}  (additive, 4p)")
        results.append(('T45', f'{popt45[0]:.1f}×φ_SE^{popt45[1]:.1f}+{popt45[2]:.1f}×φ_AM^{popt45[3]:.1f}', 1.0, r2_45, 4))
    except Exception as e:
        print(f"  T45: FAILED ({e})")

    # T46: φ_SE × φ_AM (no τ, simplest interaction)
    log_rhs_46 = np.log(phi_se * phi_am)
    log_C_46 = np.mean(log_sigma - log_rhs_46)
    pred_46 = log_C_46 + log_rhs_46
    r2_46 = 1 - np.sum((log_sigma - pred_46)**2) / ss_tot
    print(f"\n  T46: σ_th = {np.exp(log_C_46):.4f} × φ_SE × φ_AM,  R²={r2_46:.4f}")
    results.append(('T46', f'C × φ_SE × φ_AM', np.exp(log_C_46), r2_46, 1))

    # T47: porosity-based (simple!) — (1-ε)^a
    if np.any(solid_frac > 0):
        b47, _, _, _ = np.linalg.lstsq(
            np.column_stack([np.log(solid_frac), np.ones(n)]),
            log_sigma, rcond=None)
        pred47 = np.column_stack([np.log(solid_frac), np.ones(n)]) @ b47
        r2_47 = 1 - np.sum((log_sigma - pred47)**2) / ss_tot
        print(f"  T47: σ_th = C × (1-ε)^{b47[0]:.2f},  R²={r2_47:.4f}")
        results.append(('T47', f'C × (1-ε)^{b47[0]:.1f}', np.exp(b47[1]), r2_47, 2))

    # T48: PER-CASE RESIDUAL from T15 (best full model)
    # Identify what's left unexplained
    if np.sum(valid_tau) > 10:
        mask = valid_tau
        b_t15 = np.linalg.lstsq(
            np.column_stack([np.log(phi_se[mask]), np.log(phi_am[mask]), np.log(tau[mask]), np.ones(np.sum(mask))]),
            log_sigma[mask], rcond=None)[0]
        pred_t15 = np.column_stack([np.log(phi_se[mask]), np.log(phi_am[mask]), np.log(tau[mask]), np.ones(np.sum(mask))]) @ b_t15
        residuals = log_sigma[mask] - pred_t15

        print(f"\n  T48: Per-case residuals from T15 (φ_SE^{b_t15[0]:.1f}×φ_AM^{b_t15[1]:.1f}×τ^{b_t15[2]:.1f})")
        print(f"       Residual std = {np.std(residuals):.4f} ({np.exp(np.std(residuals))-1:.1%} relative)")
        # Correlate residuals with remaining vars
        resid_corrs = []
        resid_features = {
            'T': thickness[mask], 'CN_SE': se_cn[mask], 'CN_AM': am_cn[mask],
            'A_AMSE': area_am_se[mask], 'A_total': area_total[mask],
            'f_perc': f_perc[mask], 'σ_ion': sigma_ion[mask],
        }
        for name, vals in resid_features.items():
            valid = (vals > 0) & np.isfinite(vals)
            if np.sum(valid) > 5:
                r_val = np.corrcoef(vals[valid], residuals[valid])[0, 1]
                resid_corrs.append((name, r_val))
        resid_corrs.sort(key=lambda x: -abs(x[1]))
        for name, r_val in resid_corrs:
            print(f"       Residual corr with {name:8s}: r={r_val:+.3f}")

    # T49: SE × AM product with porosity
    if np.sum(valid_tau) > 10:
        mask = valid_tau
        ss49 = np.sum((log_sigma[mask] - np.mean(log_sigma[mask]))**2)
        # (φ_SE × φ_AM)^a × (1-ε)^b / τ^c
        b49, _, _, _ = np.linalg.lstsq(
            np.column_stack([np.log(phi_se[mask]*phi_am[mask]), np.log(solid_frac[mask]),
                           np.log(tau[mask]), np.ones(np.sum(mask))]),
            log_sigma[mask], rcond=None)
        pred49 = np.column_stack([np.log(phi_se[mask]*phi_am[mask]), np.log(solid_frac[mask]),
                                np.log(tau[mask]), np.ones(np.sum(mask))]) @ b49
        r2_49 = 1 - np.sum((log_sigma[mask] - pred49)**2) / ss49
        print(f"\n  T49: σ_th = C × (φ_SE×φ_AM)^{b49[0]:.2f} × (1-ε)^{b49[1]:.2f} / τ^{-b49[2]:.2f},  R²={r2_49:.4f}")
        results.append(('T49', f'(φ_SE×φ_AM)^{b49[0]:.1f}×(1-ε)^{b49[1]:.1f}×τ^{b49[2]:.1f}', np.exp(b49[3]), r2_49, 4))

    # T50: SIMPLE BEAUTIFUL — φ_SE^1.5 × φ_AM^1.5 (product = geometric mean)
    # No τ, no CN — just volume fractions
    log_rhs_50 = 1.5 * np.log(phi_se) + 1.5 * np.log(phi_am)
    log_C_50 = np.mean(log_sigma - log_rhs_50)
    pred_50 = log_C_50 + log_rhs_50
    r2_50 = 1 - np.sum((log_sigma - pred_50)**2) / ss_tot
    print(f"\n  T50: σ_th = {np.exp(log_C_50):.4f} × φ_SE^1.5 × φ_AM^1.5  (0 free + C),  R²={r2_50:.4f}")
    print(f"       = {np.exp(log_C_50):.4f} × (φ_SE × φ_AM)^1.5")
    results.append(('T50', f'C × (φ_SE×φ_AM)^1.5 (Bruggeman product)', np.exp(log_C_50), r2_50, 1))

    # ══════════════════════════════════════════════════════════════════════
    # LITERATURE-BASED MODELS (Glover 2010, Lichtenecker, 3-phase Bruggeman)
    # ══════════════════════════════════════════════════════════════════════
    print(f"\n  {'='*60}")
    print(f"  LITERATURE-BASED MODELS")
    print(f"  {'='*60}")

    # L1: Additive Bruggeman (Glover, m=1.5) — ZERO free params
    k_pred_L1 = (K_AM * phi_am**1.5 + K_SE * phi_se**1.5) * 1e3  # mS/cm equiv
    log_pred_L1 = np.log(k_pred_L1)
    r2_L1 = 1 - np.sum((log_sigma - log_pred_L1)**2) / ss_tot
    print(f"\n  L1: k_eff = k_AM×φ_AM^1.5 + k_SE×φ_SE^1.5  (0 free params)")
    print(f"      R² = {r2_L1:.4f},  predicted range: {k_pred_L1.min():.2f}~{k_pred_L1.max():.2f}")
    results.append(('L1', f'k_AM×φ_AM^1.5 + k_SE×φ_SE^1.5 (Glover)', 1.0, r2_L1, 0))

    # L1b: + scaling constant C
    log_C_L1b = np.mean(log_sigma - log_pred_L1)
    pred_L1b = log_C_L1b + log_pred_L1
    r2_L1b = 1 - np.sum((log_sigma - pred_L1b)**2) / ss_tot
    print(f"  L1b: C × (k_AM×φ_AM^1.5 + k_SE×φ_SE^1.5),  C={np.exp(log_C_L1b):.4f}, R²={r2_L1b:.4f}")
    results.append(('L1b', f'C×(k_AM×φ_AM^1.5+k_SE×φ_SE^1.5)', np.exp(log_C_L1b), r2_L1b, 1))

    # L2: Generalized Archie — k_AM×φ_AM^m_AM + k_SE×φ_SE^m_SE (2 free params)
    try:
        def model_L2(X, m_am, m_se):
            phi_a, phi_s = X
            return np.log(K_AM * 1e3 * phi_a**m_am + K_SE * 1e3 * phi_s**m_se)

        popt_L2, _ = curve_fit(model_L2, (phi_am, phi_se), log_sigma, p0=[1.5, 1.5], maxfev=10000)
        pred_L2 = model_L2((phi_am, phi_se), *popt_L2)
        r2_L2 = 1 - np.sum((log_sigma - pred_L2)**2) / ss_tot
        print(f"\n  L2: k_eff = k_AM×φ_AM^{popt_L2[0]:.2f} + k_SE×φ_SE^{popt_L2[1]:.2f}  (Generalized Archie)")
        print(f"      R² = {r2_L2:.4f}  (2 free params: m_AM, m_SE)")
        results.append(('L2', f'k_AM×φ_AM^{popt_L2[0]:.1f}+k_SE×φ_SE^{popt_L2[1]:.1f} (Archie)', 1.0, r2_L2, 2))
    except Exception as e:
        print(f"  L2: FAILED ({e})")

    # L3: Lichtenecker — k_AM^φ_AM × k_SE^φ_SE × (1-ε)^α
    try:
        k_licht_base = K_AM**phi_am * K_SE**phi_se * 1e3  # mS/cm
        def model_L3(X, alpha):
            k_base, solid = X
            return np.log(k_base * solid**alpha)

        popt_L3, _ = curve_fit(model_L3, (k_licht_base, solid_frac), log_sigma, p0=[2.0], maxfev=10000)
        pred_L3 = model_L3((k_licht_base, solid_frac), *popt_L3)
        r2_L3 = 1 - np.sum((log_sigma - pred_L3)**2) / ss_tot
        print(f"\n  L3: k_eff = k_AM^φ_AM × k_SE^φ_SE × (1-ε)^{popt_L3[0]:.2f}  (Lichtenecker)")
        print(f"      R² = {r2_L3:.4f}  (1 free param: α)")
        results.append(('L3', f'k_AM^φ_AM × k_SE^φ_SE × (1-ε)^{popt_L3[0]:.1f}', 1.0, r2_L3, 1))
    except Exception as e:
        print(f"  L3: FAILED ({e})")

    # L4: 3-phase Bruggeman (numerical, 0 free params)
    try:
        from scipy.optimize import brentq
        k_brug_3phase = []
        for i in range(n):
            pa, ps, pp = phi_am[i], phi_se[i], 1 - phi_am[i] - phi_se[i]
            if pp < 0: pp = 0
            def brug_eq(k):
                t1 = pa * (K_AM - k) / (K_AM + 2*k) if (K_AM + 2*k) > 0 else 0
                t2 = ps * (K_SE - k) / (K_SE + 2*k) if (K_SE + 2*k) > 0 else 0
                t3 = -pp / 2
                return t1 + t2 + t3
            try:
                k_sol = brentq(brug_eq, 1e-10, K_AM)
                k_brug_3phase.append(k_sol * 1e3)
            except:
                k_brug_3phase.append(0)
        k_brug_3p = np.array(k_brug_3phase)
        valid_b3 = k_brug_3p > 0
        if np.sum(valid_b3) > 10:
            r2_B3_raw = 1 - np.sum((log_sigma[valid_b3] - np.log(k_brug_3p[valid_b3]))**2) / np.sum((log_sigma[valid_b3] - np.mean(log_sigma[valid_b3]))**2)
            print(f"\n  L4: 3-phase Bruggeman (numerical, 0 free params)")
            print(f"      R² = {r2_B3_raw:.4f},  range: {k_brug_3p[valid_b3].min():.2f}~{k_brug_3p[valid_b3].max():.2f}")
            results.append(('L4', f'3-phase Bruggeman (0 params)', 1.0, r2_B3_raw, 0))

            # L4b: + scaling C
            log_C_B3 = np.mean(log_sigma[valid_b3] - np.log(k_brug_3p[valid_b3]))
            pred_B3 = log_C_B3 + np.log(k_brug_3p[valid_b3])
            r2_B3 = 1 - np.sum((log_sigma[valid_b3] - pred_B3)**2) / np.sum((log_sigma[valid_b3] - np.mean(log_sigma[valid_b3]))**2)
            print(f"  L4b: C × Bruggeman_3phase,  C={np.exp(log_C_B3):.4f}, R²={r2_B3:.4f}")
            results.append(('L4b', f'C × Bruggeman_3phase', np.exp(log_C_B3), r2_B3, 1))
    except Exception as e:
        print(f"  L4: FAILED ({e})")

    # L5: Maxwell porosity correction — k_mix × (1-ε)/(1+ε/2)
    k_maxwell = k_mix * 1e3 * (1 - porosity_arr/100) / (1 + porosity_arr/200)
    valid_mx = k_maxwell > 0
    if np.sum(valid_mx) > 10:
        log_C_mx = np.mean(log_sigma[valid_mx] - np.log(k_maxwell[valid_mx]))
        pred_mx = log_C_mx + np.log(k_maxwell[valid_mx])
        ss_mx = np.sum((log_sigma[valid_mx] - np.mean(log_sigma[valid_mx]))**2)
        r2_mx = 1 - np.sum((log_sigma[valid_mx] - pred_mx)**2) / ss_mx
        print(f"\n  L5: C × k_mix × (1-ε)/(1+ε/2)  (Maxwell porosity),  C={np.exp(log_C_mx):.4f}, R²={r2_mx:.4f}")
        results.append(('L5', f'C × k_mix × (1-ε)/(1+ε/2)', np.exp(log_C_mx), r2_mx, 1))

    # L6: Generalized Archie + porosity correction — best candidate
    try:
        def model_L6(X, m_am, m_se, alpha):
            phi_a, phi_s, solid = X
            k_archie = K_AM * 1e3 * phi_a**m_am + K_SE * 1e3 * phi_s**m_se
            return np.log(k_archie * solid**alpha)

        popt_L6, _ = curve_fit(model_L6, (phi_am, phi_se, solid_frac), log_sigma,
                               p0=[1.5, 1.5, 1.0], maxfev=10000)
        pred_L6 = model_L6((phi_am, phi_se, solid_frac), *popt_L6)
        r2_L6 = 1 - np.sum((log_sigma - pred_L6)**2) / ss_tot
        print(f"\n  L6: (k_AM×φ_AM^{popt_L6[0]:.2f} + k_SE×φ_SE^{popt_L6[1]:.2f}) × (1-ε)^{popt_L6[2]:.2f}")
        print(f"      R² = {r2_L6:.4f}  (3 free params: m_AM, m_SE, α)")
        results.append(('L6', f'Archie+porosity (3p)', 1.0, r2_L6, 3))
    except Exception as e:
        print(f"  L6: FAILED ({e})")

    # L7: Additive Bruggeman + τ — k_AM×φ_AM^1.5 + k_SE×φ_SE^1.5, then /τ^b
    if np.sum(valid_tau) > 10:
        mask = valid_tau
        try:
            def model_L7(X, C, b):
                phi_a, phi_s, tau_v = X
                k_base = K_AM * 1e3 * phi_a**1.5 + K_SE * 1e3 * phi_s**1.5
                return np.log(C * k_base / tau_v**b)

            popt_L7, _ = curve_fit(model_L7, (phi_am[mask], phi_se[mask], tau[mask]),
                                   log_sigma[mask], p0=[1.0, 1.0], maxfev=10000)
            pred_L7 = model_L7((phi_am[mask], phi_se[mask], tau[mask]), *popt_L7)
            ss_L7 = np.sum((log_sigma[mask] - np.mean(log_sigma[mask]))**2)
            r2_L7 = 1 - np.sum((log_sigma[mask] - pred_L7)**2) / ss_L7
            print(f"\n  L7: C × (k_AM×φ_AM^1.5 + k_SE×φ_SE^1.5) / τ^{popt_L7[1]:.2f}")
            print(f"      C={popt_L7[0]:.4f}, R²={r2_L7:.4f}  (2 free params)")
            results.append(('L7', f'C×Archie_1.5/τ^{popt_L7[1]:.1f}', popt_L7[0], r2_L7, 2))
        except Exception as e:
            print(f"  L7: FAILED ({e})")

    # L8: ULTIMATE — Generalized Archie + τ + porosity
    if np.sum(valid_tau) > 10:
        mask = valid_tau
        try:
            def model_L8(X, m_am, m_se, b, alpha):
                phi_a, phi_s, tau_v, solid = X
                k_base = K_AM * 1e3 * phi_a**m_am + K_SE * 1e3 * phi_s**m_se
                return np.log(k_base * solid**alpha / tau_v**b)

            popt_L8, _ = curve_fit(model_L8, (phi_am[mask], phi_se[mask], tau[mask], solid_frac[mask]),
                                   log_sigma[mask], p0=[1.5, 1.5, 1.0, 1.0], maxfev=10000)
            pred_L8 = model_L8((phi_am[mask], phi_se[mask], tau[mask], solid_frac[mask]), *popt_L8)
            ss_L8 = np.sum((log_sigma[mask] - np.mean(log_sigma[mask]))**2)
            r2_L8 = 1 - np.sum((log_sigma[mask] - pred_L8)**2) / ss_L8
            print(f"\n  L8: (k_AM×φ_AM^{popt_L8[0]:.2f} + k_SE×φ_SE^{popt_L8[1]:.2f}) × (1-ε)^{popt_L8[3]:.2f} / τ^{popt_L8[2]:.2f}")
            print(f"      R² = {r2_L8:.4f}  (4 free params)")
            results.append(('L8', f'Archie+ε+τ (4p)', 1.0, r2_L8, 4))
        except Exception as e:
            print(f"  L8: FAILED ({e})")

    # ── Ranking ──
    print(f"\n{'='*70}")
    print("RANKING — sorted by R²")
    print(f"{'='*70}")
    for tag, desc, C, r2, np_ in sorted(results, key=lambda x: -x[3]):
        print(f"  {tag:5s} R²={r2:.4f} ({np_}p)  C={C:.6f}  {desc}")

    # ── Per-case table ──
    print(f"\n{'='*70}")
    print("PER-CASE TABLE")
    print(f"{'='*70}")
    print(f"{'Name':25s} {'φ_SE':>6s} {'φ_AM':>6s} {'CN_SE':>6s} {'CN_AM':>6s} {'σ_th':>7s} {'T':>5s}")
    print("-" * 70)
    for r in sorted(rows, key=lambda x: x['sigma_th']):
        print(f"  {r['name'][:23]:23s} {r['phi_se']:6.3f} {r['phi_am']:6.3f} {r['se_cn']:6.2f} {r['am_cn']:6.2f} {r['sigma_th']:7.3f} {r['thickness']:5.0f}")


if __name__ == '__main__':
    thermal_regression()
