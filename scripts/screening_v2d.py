"""
v2.0d DEEP DIVE: φ_SE × φ_AM 결합 + CN 부호 역전 탐구
=====================================================
1. φ_SE^7 × φ_AM^5.5를 하나의 결합변수로 → 지수 감소
2. CN 음수의 물리적 근거 검증
3. Regularized fit (LOOCV) — overfitting 검출
4. Parsimonious model search (최소 변수)
"""
import json, os, sys, numpy as np, warnings
from pathlib import Path
from itertools import product, combinations

warnings.filterwarnings('ignore')
WEBAPP = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'webapp')
SIGMA_GRAIN = 3.0


def load_data():
    rows = []
    for base in [Path(WEBAPP) / 'results', Path(WEBAPP) / 'archive']:
        if not base.is_dir():
            continue
        for met_path in base.rglob('full_metrics.json'):
            try:
                with open(met_path) as f:
                    m = json.load(f)
            except:
                continue
            sigma_net = m.get('sigma_full_mScm', 0)
            if not sigma_net or sigma_net < 0.001:
                continue
            phi_se = m.get('phi_se', 0)
            phi_am = m.get('phi_am', 0)
            tau = m.get('tortuosity_recommended', m.get('tortuosity_mean', 0))
            f_perc = m.get('percolation_pct', 0) / 100
            cn = m.get('se_se_cn', 0)
            gb_d = m.get('gb_density_mean', 0)
            g_path = m.get('path_conductance_mean', 0)
            hop_area = m.get('path_hop_area_mean', 0)
            T = m.get('thickness_um', 0)
            porosity = m.get('porosity', 0)
            if tau <= 0 or phi_se <= 0 or cn <= 0 or gb_d <= 0 or T <= 0:
                continue
            rows.append({
                'name': met_path.parent.name,
                'sigma_net': sigma_net,
                'phi_se': phi_se, 'phi_am': max(phi_am, 0.01),
                'tau': tau, 'f_perc': max(f_perc, 0.5), 'cn': cn,
                'gb_d': gb_d, 'g_path': max(g_path, 1e-6),
                'hop_area': max(hop_area, 1e-6), 'T': T,
                'porosity': porosity,
                # Additional params
                'bottleneck': max(m.get('path_hop_area_min_mean', 1e-6), 1e-6),
                'am_cn': max(m.get('am_am_cn', 0.01), 0.01),
                'am_se_cn': max(m.get('am_se_cn_mean', m.get('area_AM전체_SE_n', 0.01)), 0.01),
                'coverage': max(m.get('coverage_AM_P_mean', m.get('coverage_AM_S_mean', m.get('coverage_AM_mean', 1))), 0.1),
                'se_se_area': max(m.get('area_SE_SE_mean', 1e-6), 1e-6),
                'se_se_total': max(m.get('area_SE_SE_total', 1e-6), 1e-6),
                'ionic_active': max(m.get('ionic_active_pct', 100), 1) / 100,
                'n_components': max(m.get('n_components', 1), 1),
                'largest_pct': max(m.get('largest_pct', 50), 1),
                'stress_cv': max(m.get('stress_cv', 1), 0.1),
                'bulk_frac': m.get('bulk_resistance_fraction', 0.25),
            })
    seen = set()
    unique = []
    for r in rows:
        key = f"{r['phi_se']:.4f}_{r['T']:.1f}_{r['tau']:.3f}"
        if key not in seen:
            seen.add(key)
            unique.append(r)
    return unique


def r2_log(actual, predicted):
    la, lp = np.log(actual), np.log(predicted)
    ss_res = np.sum((la - lp)**2)
    ss_tot = np.sum((la - np.mean(la))**2)
    return 1 - ss_res / ss_tot if ss_tot > 0 else -999


def fit_C(actual, rhs):
    valid = (rhs > 0) & np.isfinite(rhs)
    if valid.sum() < 3:
        return None
    return float(np.exp(np.mean(np.log(actual[valid] / rhs[valid]))))


def loocv_r2(log_sigma, X):
    """Leave-one-out cross-validation R² in log space."""
    n = len(log_sigma)
    errors = np.zeros(n)
    for i in range(n):
        mask = np.ones(n, bool)
        mask[i] = False
        try:
            coefs, _, _, _ = np.linalg.lstsq(X[mask], log_sigma[mask], rcond=None)
            errors[i] = log_sigma[i] - X[i] @ coefs
        except:
            errors[i] = 10  # large error
    ss_res = np.sum(errors**2)
    ss_tot = np.sum((log_sigma - np.mean(log_sigma))**2)
    return 1 - ss_res / ss_tot if ss_tot > 0 else -999


def regime_r2(rows, actual, predicted):
    tau = np.array([r['tau'] for r in rows])
    out = {}
    for label, mask in [('ALL', np.ones(len(rows), bool)),
                        ('τ≤1.5', tau <= 1.5),
                        ('1.5<τ≤2.5', (tau > 1.5) & (tau <= 2.5)),
                        ('τ>2.5', tau > 2.5)]:
        n = mask.sum()
        if n < 2:
            continue
        r2 = r2_log(actual[mask], predicted[mask])
        err = np.mean(np.abs(actual[mask] - predicted[mask]) / actual[mask]) * 100
        out[label] = (n, r2, err)
    return out


def main():
    rows = load_data()
    n = len(rows)
    print(f"Loaded {n} cases\n")

    sigma_net = np.array([r['sigma_net'] for r in rows])
    log_sigma = np.log(sigma_net)
    phi_se = np.array([r['phi_se'] for r in rows])
    phi_am = np.array([r['phi_am'] for r in rows])
    tau = np.array([r['tau'] for r in rows])
    cn = np.array([r['cn'] for r in rows])
    gb_d = np.array([r['gb_d'] for r in rows])
    g_path = np.array([r['g_path'] for r in rows])
    T = np.array([r['T'] for r in rows])
    f_perc = np.array([r['f_perc'] for r in rows])

    # ═══════════════════════════════════════════════
    # PART 1: φ_SE × φ_AM 결합변수 탐구
    # φ_SE^7 × φ_AM^5.5 = (φ_SE × φ_AM)^5.5 × φ_SE^1.5?
    # ═══════════════════════════════════════════════
    print("=" * 80)
    print("PART 1: φ_SE × φ_AM 결합변수")
    print("=" * 80)

    # Correlation check
    print(f"\nCorrelation matrix (log space):")
    vars_check = {'φ_SE': np.log(phi_se), 'φ_AM': np.log(phi_am), 'τ': np.log(tau),
                  'CN': np.log(cn), 'GB_d': np.log(gb_d), 'G_path': np.log(g_path)}
    names_check = list(vars_check.keys())
    mat = np.column_stack(list(vars_check.values()))
    corr = np.corrcoef(mat.T)
    # Print correlation with σ
    print(f"{'':12s}", end='')
    for nm in names_check:
        print(f"{nm:>8s}", end='')
    print(f"{'σ_net':>8s}")
    for i, nm in enumerate(names_check):
        print(f"{nm:12s}", end='')
        for j in range(len(names_check)):
            print(f"{corr[i,j]:8.2f}", end='')
        r_sigma = np.corrcoef(mat[:, i], log_sigma)[0, 1]
        print(f"{r_sigma:8.2f}")

    # Combined variables
    print(f"\n--- Combined variable test ---")
    combos = {
        'φ_SE×φ_AM': np.log(phi_se * phi_am),
        'φ_SE²×φ_AM': np.log(phi_se**2 * phi_am),
        'φ_SE×φ_AM²': np.log(phi_se * phi_am**2),
        'φ_SE/φ_AM': np.log(phi_se / phi_am),
        'φ_SE/(1-φ_SE)': np.log(phi_se / (1 - phi_se)),
        '(φ_SE×f_perc)²': np.log((phi_se * f_perc)**2),
        'φ_SE²×φ_AM/τ': np.log(phi_se**2 * phi_am / tau),
        'φ_SE×CN/τ': np.log(phi_se * cn / tau),
    }

    for cname, cvar in combos.items():
        r = np.corrcoef(cvar, log_sigma)[0, 1]
        print(f"  corr({cname:20s}, σ) = {r:.3f}")

    # ═══════════════════════════════════════════════
    # PART 2: CN 부호 역전 검증
    # φ_SE가 높으면 CN도 높음 → φ_SE^7에서 CN 이미 과다반영
    # → CN^(-0.7)이 보정하는 것일 수 있음
    # ═══════════════════════════════════════════════
    print(f"\n{'='*80}")
    print("PART 2: CN 부호 역전 검증")
    print("=" * 80)

    # Partial correlation: CN effect AFTER controlling φ_SE
    from numpy.linalg import lstsq
    # Residualize CN and σ against φ_SE
    X_phi = np.column_stack([np.log(phi_se), np.ones(n)])
    coef_cn, _, _, _ = lstsq(X_phi, np.log(cn), rcond=None)
    coef_sig, _, _, _ = lstsq(X_phi, log_sigma, rcond=None)
    resid_cn = np.log(cn) - X_phi @ coef_cn
    resid_sig = log_sigma - X_phi @ coef_sig
    partial_corr = np.corrcoef(resid_cn, resid_sig)[0, 1]
    print(f"\n  Partial corr(CN, σ | φ_SE) = {partial_corr:.3f}")
    print(f"  Raw corr(CN, σ) = {np.corrcoef(np.log(cn), log_sigma)[0,1]:.3f}")
    print(f"  Raw corr(CN, φ_SE) = {np.corrcoef(np.log(cn), np.log(phi_se))[0,1]:.3f}")

    if partial_corr < 0:
        print("  → CN has NEGATIVE partial effect on σ after controlling φ_SE")
        print("  → CN 증가는 φ_SE 증가를 동반. φ_SE 통제 후 CN↑ → σ↓")
        print("  → 물리: 같은 φ_SE에서 CN↑ = 입자 더 작음 → constriction↑ → σ↓")
    else:
        print("  → CN has POSITIVE partial effect on σ after controlling φ_SE")

    # Same for G_path
    coef_gp, _, _, _ = lstsq(X_phi, np.log(g_path), rcond=None)
    resid_gp = np.log(g_path) - X_phi @ coef_gp
    partial_gp = np.corrcoef(resid_gp, resid_sig)[0, 1]
    print(f"\n  Partial corr(G_path, σ | φ_SE) = {partial_gp:.3f}")

    # ═══════════════════════════════════════════════
    # PART 3: Parsimonious models (최소 변수)
    # 2변수로 R²>0.93 가능한가?
    # ═══════════════════════════════════════════════
    print(f"\n{'='*80}")
    print("PART 3: Parsimonious models (2-3 변수, LOOCV)")
    print("=" * 80)

    # ALL available parameters
    bottleneck = np.array([r['bottleneck'] for r in rows])
    am_cn = np.array([r['am_cn'] for r in rows])
    am_se_cn = np.array([r['am_se_cn'] for r in rows])
    coverage = np.array([r['coverage'] for r in rows])
    se_se_area = np.array([r['se_se_area'] for r in rows])
    ionic_active = np.array([r['ionic_active'] for r in rows])
    stress_cv = np.array([r['stress_cv'] for r in rows])

    var_pool = {
        # Basic
        'φ_SE': np.log(phi_se), 'φ_AM': np.log(phi_am),
        'τ': np.log(tau), 'CN': np.log(cn),
        'GB_d': np.log(gb_d), 'G_path': np.log(g_path),
        'f_perc': np.log(f_perc), 'T': np.log(T),
        'hop_area': np.log(np.array([r['hop_area'] for r in rows])),
        'BN': np.log(bottleneck),
        # AM-related
        'AM_CN': np.log(am_cn),
        'AM_SE_CN': np.log(am_se_cn),
        'coverage': np.log(coverage),
        # SE area
        'SE_area': np.log(se_se_area),
        'ionic_act': np.log(ionic_active),
        'stress': np.log(stress_cv),
        # Combined
        'φ_SE×φ_AM': np.log(phi_se * phi_am),
        'φ_SE²×φ_AM': np.log(phi_se**2 * phi_am),
        'CN/τ': np.log(cn / tau),
        'CN²/τ': np.log(cn**2 / tau),
        'G×d²': np.log(g_path * gb_d**2),
        'φ_SE×CN/τ': np.log(phi_se * cn / tau),
        'φ_SE×f/τ²': np.log(phi_se * f_perc / tau**2),  # = σ_brug/σ_grain
        'BN×d²': np.log(bottleneck * gb_d**2),
        'hop×d²': np.log(np.array([r['hop_area'] for r in rows]) * gb_d**2),
        'τ×GB_d': np.log(tau * gb_d),
        'φ_SE/τ': np.log(phi_se / tau),
        'CN×φ_SE': np.log(cn * phi_se),
    }

    results_23 = []
    for nvars in [2, 3]:
        for vnames in combinations(var_pool.keys(), nvars):
            X = np.column_stack([var_pool[v] for v in vnames] + [np.ones(n)])
            try:
                coefs, _, _, _ = lstsq(X, log_sigma, rcond=None)
                pred = np.exp(X @ coefs)
                r2 = r2_log(sigma_net, pred)
                r2_cv = loocv_r2(log_sigma, X)
                if r2 > 0.90:
                    regime = regime_r2(rows, sigma_net, pred)
                    tau_gt = regime.get('τ>2.5', (0, -99, 999))
                    results_23.append({
                        'nvars': nvars,
                        'vnames': vnames,
                        'coefs': coefs,
                        'r2': r2, 'r2_cv': r2_cv,
                        'tau_gt_r2': tau_gt[1],
                        'tau_gt_err': tau_gt[2] if tau_gt[0] > 0 else 999,
                        'pred': pred,
                        'C': np.exp(coefs[-1]),
                    })
            except:
                pass

    results_23.sort(key=lambda x: -x['r2_cv'])

    print(f"\nTop 2-3 variable models (sorted by LOOCV R²):")
    print(f"{'Rank':>4s} {'Vars':>4s} {'R²':>7s} {'LOOCV':>7s} {'τ>2.5':>7s} Formula")
    print("-" * 90)
    count = 0
    for r in results_23[:30]:
        count += 1
        terms = ' × '.join(f'{v}^{r["coefs"][i]:.2f}' for i, v in enumerate(r['vnames']))
        flag = '★' if r['tau_gt_r2'] > 0.5 else '⚠' if r['tau_gt_r2'] > 0 else ' '
        print(f"#{count:3d} {r['nvars']:4d} {r['r2']:7.4f} {r['r2_cv']:7.4f} {r['tau_gt_r2']:7.3f} {flag} {r['C']:.2f} × {terms}")

    # ═══════════════════════════════════════════════
    # PART 4: Fixed exponent candidates from best free fits
    # Round best free-fit exponents to clean fractions
    # ═══════════════════════════════════════════════
    print(f"\n{'='*80}")
    print("PART 4: Fixed exponent rounding (best → clean fractions)")
    print("=" * 80)

    # Take top 5 from PART 3 and try rounding
    for r in results_23[:5]:
        print(f"\n  Free fit: {' × '.join(f'{v}^{r['coefs'][i]:.2f}' for i, v in enumerate(r['vnames']))}")
        print(f"    R²={r['r2']:.4f}, LOOCV={r['r2_cv']:.4f}")

        # Try rounding each exponent to nearest 0.5
        rounded_options = []
        for rounding in [0.25, 0.5, 1.0]:
            rounded_exp = [round(r['coefs'][i] / rounding) * rounding for i in range(len(r['vnames']))]
            rhs = np.ones(n)
            for i, v in enumerate(r['vnames']):
                rhs *= np.exp(var_pool[v]) ** rounded_exp[i]
            C = fit_C(sigma_net, rhs)
            if C and C > 0:
                pred = C * rhs
                r2_fixed = r2_log(sigma_net, pred)
                regime = regime_r2(rows, sigma_net, pred)
                tau_gt = regime.get('τ>2.5', (0, -99, 999))
                terms = ' × '.join(f'{v}^{rounded_exp[i]}' for i, v in enumerate(r['vnames']))
                gap = r['r2'] - r2_fixed
                print(f"    Round({rounding}): {C:.2f} × {terms}  R²={r2_fixed:.4f} (gap={gap:.4f}) τ>2.5={tau_gt[1]:.3f}")

    # ═══════════════════════════════════════════════
    # PART 5: THE ULTIMATE — σ_brug based with φ_AM
    # Since σ_brug = σ_grain × φ_SE × f_perc / τ²
    # Test: σ_brug^a × φ_AM^b × (contact)^c
    # ═══════════════════════════════════════════════
    print(f"\n{'='*80}")
    print("PART 5: σ_brug^a × φ_AM^b × (contact)^c")
    print("=" * 80)

    sigma_brug = SIGMA_GRAIN * phi_se * f_perc / tau**2
    log_brug = np.log(sigma_brug)

    contact_pool = {
        'CN': np.log(cn), 'GB_d': np.log(gb_d),
        'G_path': np.log(g_path), 'CN/τ': np.log(cn/tau),
        'G×d²': np.log(g_path * gb_d**2),
    }

    results_brug = []
    # σ_brug^a × φ_AM^b × contact^c (free a, b, c)
    for cname, cvar in contact_pool.items():
        X = np.column_stack([log_brug, np.log(phi_am), cvar, np.ones(n)])
        try:
            coefs, _, _, _ = np.linalg.lstsq(X, log_sigma, rcond=None)
            pred = np.exp(X @ coefs)
            r2 = r2_log(sigma_net, pred)
            r2_cv = loocv_r2(log_sigma, X)
            regime = regime_r2(rows, sigma_net, pred)
            tau_gt = regime.get('τ>2.5', (0, -99, 999))
            results_brug.append({
                'formula': f'σ_brug^{coefs[0]:.2f} × φ_AM^{coefs[1]:.2f} × {cname}^{coefs[2]:.2f}',
                'r2': r2, 'r2_cv': r2_cv,
                'tau_gt_r2': tau_gt[1],
                'coefs': coefs, 'cname': cname, 'pred': pred,
            })
        except:
            pass

    # Also: σ_brug^a × φ_AM^b (no contact)
    X_simple = np.column_stack([log_brug, np.log(phi_am), np.ones(n)])
    coefs_s, _, _, _ = np.linalg.lstsq(X_simple, log_sigma, rcond=None)
    pred_s = np.exp(X_simple @ coefs_s)
    r2_s = r2_log(sigma_net, pred_s)
    r2_cv_s = loocv_r2(log_sigma, X_simple)
    regime_s = regime_r2(rows, sigma_net, pred_s)
    tau_gt_s = regime_s.get('τ>2.5', (0, -99, 999))
    print(f"\n  σ_brug^{coefs_s[0]:.2f} × φ_AM^{coefs_s[1]:.2f} × {np.exp(coefs_s[2]):.4f}")
    print(f"    R²={r2_s:.4f}, LOOCV={r2_cv_s:.4f}, τ>2.5 R²={tau_gt_s[1]:.3f}")
    for label, (nn, r2v, errv) in regime_s.items():
        print(f"    {label:12s}: n={nn:2d}  R²={r2v:.3f}  |err|={errv:.1f}%")

    results_brug.sort(key=lambda x: -x['r2_cv'])
    print(f"\n  σ_brug^a × φ_AM^b × contact^c (sorted by LOOCV):")
    for r in results_brug[:10]:
        flag = '★' if r['tau_gt_r2'] > 0.5 else '⚠'
        print(f"    {flag} R²={r['r2']:.4f} LOOCV={r['r2_cv']:.4f} τ>2.5={r['tau_gt_r2']:.3f}  {r['formula']}")

    # ═══════════════════════════════════════════════
    # PART 6: σ_brug^a × φ_AM^b fixed exponents
    # ═══════════════════════════════════════════════
    print(f"\n{'='*80}")
    print("PART 6: σ_brug^a × φ_AM^b (FIXED exponents, 1 free param)")
    print("=" * 80)

    best_fixed = []
    for a in np.arange(0.5, 3.1, 0.25):
        for b in np.arange(-2, 5.1, 0.25):
            rhs = sigma_brug**a * phi_am**b
            C = fit_C(sigma_net, rhs)
            if C is None or C <= 0:
                continue
            pred = C * rhs
            r2 = r2_log(sigma_net, pred)
            if r2 > 0.85:
                regime = regime_r2(rows, sigma_net, pred)
                tau_gt = regime.get('τ>2.5', (0, -99, 999))
                best_fixed.append({
                    'a': a, 'b': b, 'C': C, 'r2': r2,
                    'tau_gt_r2': tau_gt[1], 'pred': pred,
                })

    best_fixed.sort(key=lambda x: -x['r2'])
    print(f"\nTop σ_brug^a × φ_AM^b (1 free param):")
    for i, r in enumerate(best_fixed[:20]):
        flag = '★' if r['tau_gt_r2'] > 0.5 else '⚠' if r['tau_gt_r2'] > 0 else ' '
        print(f"  #{i+1:2d} {flag} R²={r['r2']:.4f}  σ_brug^{r['a']:.2f} × φ_AM^{r['b']:.2f} × {r['C']:.4f}   τ>2.5 R²={r['tau_gt_r2']:.3f}")

    # Best overall with contact
    print(f"\n  + contact correction:")
    for cname, cvar in contact_pool.items():
        for r in best_fixed[:3]:
            for c_exp in [-1, -0.5, -0.25, 0.25, 0.5, 1]:
                rhs2 = sigma_brug**r['a'] * phi_am**r['b'] * np.exp(cvar)**c_exp
                C2 = fit_C(sigma_net, rhs2)
                if C2 is None or C2 <= 0:
                    continue
                pred2 = C2 * rhs2
                r2_2 = r2_log(sigma_net, pred2)
                regime2 = regime_r2(rows, sigma_net, pred2)
                tau_gt2 = regime2.get('τ>2.5', (0, -99, 999))
                if r2_2 > r['r2'] + 0.005:
                    flag = '★' if tau_gt2[1] > 0.5 else '⚠'
                    print(f"    {flag} R²={r2_2:.4f} σ_brug^{r['a']:.2f}×φ_AM^{r['b']:.2f}×{cname}^{c_exp}×{C2:.4f}  τ>2.5={tau_gt2[1]:.3f}")


if __name__ == '__main__':
    main()
