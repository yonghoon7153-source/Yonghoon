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
