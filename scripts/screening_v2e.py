"""
v2.0e PHYSICS-DRIVEN Screening
===============================
v4 후보의 물리적 근거를 찾고, 더 정교한 공식 탐색.

핵심 질문:
1. φ_SE × φ_AM이 왜 되는가? → 물리적 분해
2. τ² 대신 뭘 써야 하는가? → softened Bruggeman
3. σ_brug를 살릴 수 있는가? → τ^n (n<2) + φ_AM 보정
4. (φ_SE×φ_AM)^7의 지수 7을 분해할 수 있는가?
"""
import json, os, sys, numpy as np, warnings
from pathlib import Path
from itertools import product

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
            coverage = max(m.get('coverage_AM_P_mean', m.get('coverage_AM_S_mean', m.get('coverage_AM_mean', 1))), 0.1)
            se_area = max(m.get('area_SE_SE_mean', 1e-6), 1e-6)
            if tau <= 0 or phi_se <= 0 or cn <= 0 or gb_d <= 0 or T <= 0:
                continue
            rows.append({
                'name': met_path.parent.name, 'sigma_net': sigma_net,
                'phi_se': phi_se, 'phi_am': max(phi_am, 0.01),
                'tau': tau, 'f_perc': max(f_perc, 0.5), 'cn': cn,
                'gb_d': gb_d, 'g_path': max(g_path, 1e-6),
                'hop_area': max(hop_area, 1e-6), 'T': T,
                'porosity': porosity, 'coverage': coverage,
                'se_area': se_area,
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
    n = len(log_sigma)
    errors = np.zeros(n)
    for i in range(n):
        mask = np.ones(n, bool); mask[i] = False
        try:
            coefs, _, _, _ = np.linalg.lstsq(X[mask], log_sigma[mask], rcond=None)
            errors[i] = log_sigma[i] - X[i] @ coefs
        except:
            errors[i] = 10
    ss_res = np.sum(errors**2)
    ss_tot = np.sum((log_sigma - np.mean(log_sigma))**2)
    return 1 - ss_res / ss_tot if ss_tot > 0 else -999


def regime(rows, actual, predicted):
    tau = np.array([r['tau'] for r in rows])
    out = {}
    for label, mask in [('ALL', np.ones(len(rows), bool)),
                        ('thick', tau <= 1.5), ('mid', (tau > 1.5) & (tau <= 2.5)),
                        ('thin', tau > 2.5)]:
        n = mask.sum()
        if n < 2: continue
        r2 = r2_log(actual[mask], predicted[mask])
        err = np.mean(np.abs(actual[mask] - predicted[mask]) / actual[mask]) * 100
        out[label] = (n, r2, err)
    return out


def test_formula(name, sigma_net, rhs, rows):
    """Test a formula, report R², LOOCV (approx), per-regime."""
    C = fit_C(sigma_net, rhs)
    if C is None or C <= 0:
        return None
    pred = C * rhs
    r2 = r2_log(sigma_net, pred)
    reg = regime(rows, sigma_net, pred)
    thin_r2 = reg.get('thin', (0, -99, 0))[1]
    return {'name': name, 'r2': r2, 'C': C, 'pred': pred, 'reg': reg, 'thin_r2': thin_r2}


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
    hop_area = np.array([r['hop_area'] for r in rows])
    T = np.array([r['T'] for r in rows])
    f_perc = np.array([r['f_perc'] for r in rows])
    porosity = np.array([r['porosity'] for r in rows])
    coverage = np.array([r['coverage'] for r in rows])
    se_area = np.array([r['se_area'] for r in rows])

    results = []

    # ═══════════════════════════════════════
    # PART 1: SOFTENED BRUGGEMAN
    # σ = σ_grain × φ_SE × f_perc / τ^n × C × correction
    # Key: n < 2 to avoid τ² over-penalty
    # ═══════════════════════════════════════
    print("=" * 80)
    print("PART 1: SOFTENED BRUGGEMAN — τ^n with n = 0.5 ~ 3")
    print("σ = σ_grain × φ_SE × f_perc / τ^n × C × φ_AM^a × (contact)^b")
    print("=" * 80)

    for tau_n in [0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 2.5, 3]:
        sigma_soft = SIGMA_GRAIN * phi_se * f_perc / tau**tau_n
        # Pure softened Bruggeman
        r = test_formula(f'σ_grain×φ×f/τ^{tau_n}', sigma_net, sigma_soft, rows)
        if r and r['r2'] > 0.5:
            results.append(r)
        # + φ_AM^a
        for am_exp in [1, 2, 3, 4]:
            rhs = sigma_soft * phi_am**am_exp
            r = test_formula(f'σ_soft(τ^{tau_n}) × φ_AM^{am_exp}', sigma_net, rhs, rows)
            if r and r['r2'] > 0.8:
                results.append(r)
        # + φ_AM^a × (G_path×GB_d²)^b
        for am_exp in [1, 2, 3]:
            for gb_exp in [0.25, 0.5]:
                rhs = sigma_soft * phi_am**am_exp * (g_path * gb_d**2)**gb_exp
                r = test_formula(f'σ_soft(τ^{tau_n}) × φ_AM^{am_exp} × (G×d²)^{gb_exp}', sigma_net, rhs, rows)
                if r and r['r2'] > 0.85:
                    results.append(r)
        # + φ_AM^a × CN^b
        for am_exp in [1, 2, 3]:
            for cn_exp in [1, 2]:
                rhs = sigma_soft * phi_am**am_exp * cn**cn_exp
                r = test_formula(f'σ_soft(τ^{tau_n}) × φ_AM^{am_exp} × CN^{cn_exp}', sigma_net, rhs, rows)
                if r and r['r2'] > 0.85:
                    results.append(r)

    # ═══════════════════════════════════════
    # PART 2: φ_SE × φ_AM DECOMPOSITION
    # Is it really (φ_SE×φ_AM)^n, or φ_SE^a × φ_AM^b?
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("PART 2: φ_SE^a × φ_AM^b 분해 (a ≠ b)")
    print("=" * 80)

    for a in np.arange(2, 8.5, 0.5):
        for b in np.arange(0, 7.5, 0.5):
            rhs = phi_se**a * phi_am**b
            r = test_formula(f'φ_SE^{a:.1f} × φ_AM^{b:.1f}', sigma_net, rhs, rows)
            if r and r['r2'] > 0.88:
                results.append(r)

    # Best φ_SE^a × φ_AM^b + contact correction
    best_ab = sorted([r for r in results if r['name'].startswith('φ_SE^')], key=lambda x: -x['r2'])[:5]
    for base_r in best_ab:
        # Extract a, b from name
        parts = base_r['name'].split('×')
        a_val = float(parts[0].split('^')[1])
        b_val = float(parts[1].strip().split('^')[1])
        base_rhs = phi_se**a_val * phi_am**b_val

        for cname, cvar in [('G_path', g_path), ('CN/τ', cn/tau), ('hop_area', hop_area),
                            ('SE_area', se_area), ('G×d²', g_path*gb_d**2)]:
            for c_exp in [-0.5, -0.25, 0.25, 0.5]:
                rhs = base_rhs * cvar**c_exp
                name = f'φ_SE^{a_val:.1f}×φ_AM^{b_val:.1f}×{cname}^{c_exp}'
                r = test_formula(name, sigma_net, rhs, rows)
                if r and r['r2'] > base_r['r2']:
                    results.append(r)

    # ═══════════════════════════════════════
    # PART 3: PHYSICAL DECOMPOSITION of (φ×φ)^7
    # 7 = Bruggeman(1.5) + contact(2) + network(2) + finite_size(1.5)?
    # Test: φ_SE^1.5 × (φ_SE×CN)^2 × (φ_AM)^b
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("PART 3: PHYSICAL DECOMPOSITION — 지수 7의 기원")
    print("σ = C × φ_SE^n_brug × (φ_SE×CN²)^n_net × φ_AM^n_comp × G_path^n_con")
    print("=" * 80)

    # Decomposition: σ = C × φ^1.5 × CN^2 × φ_AM^b × correction
    # This recovers Bruggeman (φ^1.5) + current champion (CN²) + new (φ_AM)
    for brug_exp in [1, 1.5, 2]:
        for cn_exp in [0, 1, 2]:
            for am_exp in [0, 1, 2, 3, 4]:
                for corr_name, corr in [('1', np.ones(n)), ('G_path^0.25', g_path**0.25),
                                         ('G_path^-0.25', g_path**(-0.25)),
                                         ('(G×d²)^0.25', (g_path*gb_d**2)**0.25),
                                         ('SE_area^-0.25', se_area**(-0.25)),
                                         ('hop^0.25', hop_area**0.25)]:
                    rhs = phi_se**brug_exp * cn**cn_exp * phi_am**am_exp * corr
                    total_exp = brug_exp + cn_exp + am_exp
                    name = f'φ^{brug_exp}×CN^{cn_exp}×φ_AM^{am_exp}×{corr_name}'
                    r = test_formula(name, sigma_net, rhs, rows)
                    if r and r['r2'] > 0.88:
                        results.append(r)

    # ═══════════════════════════════════════
    # PART 4: PERCOLATION-WEIGHTED formulas
    # (φ_SE × f_perc) instead of φ_SE
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("PART 4: PERCOLATION-WEIGHTED — (φ_SE × f_perc)")
    print("=" * 80)

    phi_eff = phi_se * f_perc  # effective connected SE fraction

    for a in [3, 4, 5, 6, 7]:
        for b in [0, 1, 2, 3, 4]:
            rhs = phi_eff**a * phi_am**b
            r = test_formula(f'(φ×f)^{a} × φ_AM^{b}', sigma_net, rhs, rows)
            if r and r['r2'] > 0.88:
                results.append(r)
            # + contact
            for c_exp in [-0.25, 0.25]:
                for cname, cvar in [('G_path', g_path), ('SE_area', se_area)]:
                    rhs2 = phi_eff**a * phi_am**b * cvar**c_exp
                    r2 = test_formula(f'(φf)^{a}×φ_AM^{b}×{cname}^{c_exp}', sigma_net, rhs2, rows)
                    if r2 and r2['r2'] > 0.9:
                        results.append(r2)

    # ═══════════════════════════════════════
    # PART 5: MIXED — σ_brug(softened) × (φ_SE×φ_AM)^k
    # Best of both worlds?
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("PART 5: HYBRID — σ_soft × (φ_SE×φ_AM)^k")
    print("=" * 80)

    phi_prod = phi_se * phi_am
    for tau_n in [0.5, 1, 1.5, 2]:
        sigma_soft = SIGMA_GRAIN * phi_se * f_perc / tau**tau_n
        for k in [1, 2, 3, 4, 5]:
            rhs = sigma_soft * phi_prod**k
            r = test_formula(f'σ_soft(τ^{tau_n}) × (φ×φ_AM)^{k}', sigma_net, rhs, rows)
            if r and r['r2'] > 0.85:
                results.append(r)
            # + contact
            for c_exp in [-0.25, 0.25]:
                rhs2 = sigma_soft * phi_prod**k * g_path**c_exp
                r2 = test_formula(f'σ_soft(τ^{tau_n})×(φφ)^{k}×G^{c_exp}', sigma_net, rhs2, rows)
                if r2 and r2['r2'] > 0.9:
                    results.append(r2)

    # ═══════════════════════════════════════
    # RESULTS
    # ═══════════════════════════════════════
    results.sort(key=lambda x: -x['r2'])

    # Remove duplicates
    seen = set()
    unique_results = []
    for r in results:
        if r['name'] not in seen:
            seen.add(r['name'])
            unique_results.append(r)

    print(f"\n{'='*80}")
    print(f"TOP 40 PHYSICS-MOTIVATED FORMULAS (n={n})")
    print("=" * 80)

    for i, r in enumerate(unique_results[:40]):
        reg = r['reg']
        thin = reg.get('thin', (0, -99, 999))
        thick = reg.get('thick', (0, -99, 999))
        flag = '★★' if thin[1] > 0.8 else '★' if thin[1] > 0.5 else '⚠' if thin[1] > 0 else '  '
        print(f"\n#{i+1:2d} {flag} R²={r['r2']:.4f} C={r['C']:.4f}  thin={thin[1]:.3f}")
        print(f"   {r['name']}")
        for label, (nn, r2v, errv) in reg.items():
            f2 = '✓' if r2v > 0.85 else '⚠' if r2v > 0.5 else '✗'
            print(f"   {f2} {label:6s}: n={nn:2d} R²={r2v:.3f} |err|={errv:.1f}%")

    # ═══════════════════════════════════════
    # SPECIAL: σ_brug를 살리는 최고 모델
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("BRUG-COMPATIBLE: σ_brug를 포함하는 최고 모델")
    print("=" * 80)

    brug_results = [r for r in unique_results if 'σ_soft' in r['name'] or 'σ_grain' in r['name']]
    for i, r in enumerate(brug_results[:15]):
        reg = r['reg']
        thin = reg.get('thin', (0, -99, 999))
        print(f"  #{i+1} R²={r['r2']:.4f} thin={thin[1]:.3f}  {r['name']}")

    # ═══════════════════════════════════════
    # KEY COMPARISON TABLE
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("KEY COMPARISON")
    print("=" * 80)

    # v3 champion
    rhs_v3 = SIGMA_GRAIN * phi_se * f_perc / tau**2 * (g_path * gb_d**2)**0.25 * cn**2
    C_v3 = fit_C(sigma_net, rhs_v3)
    pred_v3 = C_v3 * rhs_v3
    reg_v3 = regime(rows, sigma_net, pred_v3)

    # v4 simple: (φ_SE×φ_AM)^7 × G_path^(-0.2)
    rhs_v4s = (phi_se * phi_am)**7 * g_path**(-0.2)
    C_v4s = fit_C(sigma_net, rhs_v4s)
    pred_v4s = C_v4s * rhs_v4s
    reg_v4s = regime(rows, sigma_net, pred_v4s)

    # Best from this screening
    best = unique_results[0] if unique_results else None

    print(f"\n{'Model':<45s} {'ALL':>7s} {'thick':>7s} {'mid':>7s} {'thin':>7s}")
    print("-" * 80)
    for label, reg_data in [('v3: σ_brug×C×(G×d²)^¼×CN²', reg_v3),
                            ('v4s: (φ_SE×φ_AM)^7 × G^(-0.2)', reg_v4s)]:
        vals = [f"{reg_data.get(k, (0,-99,0))[1]:.3f}" for k in ['ALL', 'thick', 'mid', 'thin']]
        print(f"  {label:<43s} {vals[0]:>7s} {vals[1]:>7s} {vals[2]:>7s} {vals[3]:>7s}")
    if best:
        vals = [f"{best['reg'].get(k, (0,-99,0))[1]:.3f}" for k in ['ALL', 'thick', 'mid', 'thin']]
        print(f"  BEST: {best['name'][:37]:<37s} {vals[0]:>7s} {vals[1]:>7s} {vals[2]:>7s} {vals[3]:>7s}")


if __name__ == '__main__':
    main()
