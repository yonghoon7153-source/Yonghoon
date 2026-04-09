"""
v2.0c ULTIMATE Screening: τ>2.5 정복을 위한 극한 탐색
====================================================
τ² penalty가 부족하다면? τ³, τ⁴, exp(τ) 등 비표준 correction 테스트.
thin100 15개 추가를 대비하여 방향성 탐색.
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
            if tau <= 0 or phi_se <= 0 or cn <= 0 or gb_d <= 0 or T <= 0 or g_path <= 0:
                continue
            rows.append({
                'name': met_path.parent.name,
                'sigma_net': sigma_net,
                'phi_se': phi_se, 'phi_am': max(phi_am, 0.01),
                'tau': tau, 'f_perc': f_perc, 'cn': cn,
                'gb_d': gb_d, 'g_path': g_path,
                'hop_area': max(hop_area, 1e-6), 'T': T,
                'porosity': porosity,
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
    print(f"Loaded {n} cases")

    sigma_net = np.array([r['sigma_net'] for r in rows])
    phi_se = np.array([r['phi_se'] for r in rows])
    phi_am = np.array([r['phi_am'] for r in rows])
    tau = np.array([r['tau'] for r in rows])
    cn = np.array([r['cn'] for r in rows])
    gb_d = np.array([r['gb_d'] for r in rows])
    g_path = np.array([r['g_path'] for r in rows])
    T = np.array([r['T'] for r in rows])
    f_perc = np.array([r['f_perc'] for r in rows])

    results = []

    # ═══════════════════════════════════════════════
    # PART 1: τ exponent beyond quadratic
    # σ_brug uses τ². What if real penalty is τ^n?
    # Test: σ_grain × φ_SE × f_perc / τ^n × correction
    # ═══════════════════════════════════════════════
    print("\n=== PART 1: Generalized τ exponent ===")
    print("σ = C × σ_grain × φ_SE × f_perc / τ^n × (G×d²)^a × CN^b × φ_AM^c")
    print("-" * 70)

    for tau_exp in [2, 2.5, 3, 3.5, 4, 5, 6]:
        sigma_gen = SIGMA_GRAIN * phi_se * f_perc / tau**tau_exp
        for a in [0, 1/4, 1/2]:
            for b in [0, 1, 2]:
                for c in [0, 1, 2, 3]:
                    rhs = sigma_gen * (g_path * gb_d**2)**a * cn**b * phi_am**c
                    C = fit_C(sigma_net, rhs)
                    if C is None or C <= 0:
                        continue
                    pred = C * rhs
                    r2 = r2_log(sigma_net, pred)
                    if r2 > 0.85:
                        results.append({
                            'part': 1, 'r2': r2,
                            'formula': f'σ_grain×φ×f/τ^{tau_exp} × C×(G×d²)^{a:.2f}×CN^{b}×φ_AM^{c}',
                            'C': C, 'pred': pred,
                            'tau_exp': tau_exp,
                        })

    # ═══════════════════════════════════════════════
    # PART 2: exp(-τ) type corrections
    # Maybe τ penalty is exponential, not power law
    # ═══════════════════════════════════════════════
    print("\n=== PART 2: Exponential τ corrections ===")
    print("σ = C × φ_SE^a × φ_AM^b × CN^c × exp(-d×τ)")
    print("-" * 70)

    for d_coef in [0.5, 1, 1.5, 2, 2.5, 3]:
        exp_tau = np.exp(-d_coef * tau)
        for a in [2, 3, 3.5, 4]:
            for b in [0, 1, 2, 3]:
                for c in [0, 0.5, 1]:
                    rhs = phi_se**a * phi_am**b * cn**c * exp_tau
                    C = fit_C(sigma_net, rhs)
                    if C is None or C <= 0:
                        continue
                    pred = C * rhs
                    r2 = r2_log(sigma_net, pred)
                    if r2 > 0.85:
                        results.append({
                            'part': 2, 'r2': r2,
                            'formula': f'C×φ_SE^{a}×φ_AM^{b}×CN^{c}×exp(-{d_coef}×τ)',
                            'C': C, 'pred': pred,
                        })

    # ═══════════════════════════════════════════════
    # PART 3: CN/τ^n as unified variable
    # CN/τ worked well. What about CN/τ^2, CN²/τ³?
    # ═══════════════════════════════════════════════
    print("\n=== PART 3: CN/τ^n generalized ===")
    print("-" * 70)

    for tn in [1, 1.5, 2, 2.5, 3]:
        cn_tau_n = cn / tau**tn
        for a in [3, 3.5, 4, 4.5, 5]:
            for b in [0, 1, 2, 3, 4]:
                for c in [0.5, 0.75, 1, 1.5, 2]:
                    rhs = phi_se**a * phi_am**b * cn_tau_n**c
                    C = fit_C(sigma_net, rhs)
                    if C is None or C <= 0:
                        continue
                    pred = C * rhs
                    r2 = r2_log(sigma_net, pred)
                    if r2 > 0.9:
                        results.append({
                            'part': 3, 'r2': r2,
                            'formula': f'C×φ_SE^{a}×φ_AM^{b}×(CN/τ^{tn})^{c}',
                            'C': C, 'pred': pred,
                        })

    # ═══════════════════════════════════════════════
    # PART 4: (φ_SE × f_perc)^a / τ^b — generalized Bruggeman
    # Standard: (φ×f)^1 / τ^2. What if exponents differ?
    # ═══════════════════════════════════════════════
    print("\n=== PART 4: Generalized Bruggeman ===")
    print("σ = C × (φ_SE×f_perc)^a / τ^b × (contact)^c × φ_AM^d")
    print("-" * 70)

    phi_f = phi_se * f_perc
    contact_vars = {
        'G×d²': g_path * gb_d**2,
        'CN': cn,
        'CN/τ': cn / tau,
        'g_path': g_path,
    }

    for cname, cvar in contact_vars.items():
        for a in [1, 1.5, 2, 2.5, 3]:
            for b in [2, 2.5, 3, 3.5, 4, 5]:
                for c in [0, 0.25, 0.5, 0.75, 1]:
                    for d in [0, 1, 2, 3]:
                        rhs = SIGMA_GRAIN * phi_f**a / tau**b * cvar**c * phi_am**d
                        C = fit_C(sigma_net, rhs)
                        if C is None or C <= 0:
                            continue
                        pred = C * rhs
                        r2 = r2_log(sigma_net, pred)
                        if r2 > 0.9:
                            results.append({
                                'part': 4, 'r2': r2,
                                'formula': f'σ_g×C×(φf)^{a}/τ^{b}×{cname}^{c}×φ_AM^{d}',
                                'C': C, 'pred': pred,
                            })

    # ═══════════════════════════════════════════════
    # PART 5: Free exponent fit (log-linear regression)
    # Let data decide ALL exponents
    # ═══════════════════════════════════════════════
    print("\n=== PART 5: Full free fit ===")

    log_sigma = np.log(sigma_net)
    var_pool = {
        'φ_SE': np.log(phi_se), 'φ_AM': np.log(phi_am),
        'τ': np.log(tau), 'CN': np.log(cn),
        'GB_d': np.log(gb_d), 'G_path': np.log(g_path),
        'f_perc': np.log(np.clip(f_perc, 0.5, 1)),
        'T': np.log(T),
    }

    from itertools import combinations
    for nvars in [3, 4, 5]:
        for vnames in combinations(var_pool.keys(), nvars):
            X = np.column_stack([var_pool[v] for v in vnames] + [np.ones(n)])
            try:
                coefs, _, _, _ = np.linalg.lstsq(X, log_sigma, rcond=None)
                pred = np.exp(X @ coefs)
                r2 = r2_log(sigma_net, pred)
                if r2 > 0.92:
                    terms = ' × '.join(f'{v}^{coefs[i]:.2f}' for i, v in enumerate(vnames))
                    results.append({
                        'part': 5, 'r2': r2,
                        'formula': f'{np.exp(coefs[-1]):.4f} × {terms}',
                        'C': np.exp(coefs[-1]), 'pred': pred,
                        'n_free': nvars + 1,
                    })
            except:
                pass

    # ═══════════════════════════════════════════════
    # RESULTS
    # ═══════════════════════════════════════════════
    results.sort(key=lambda x: -x['r2'])

    print("\n" + "=" * 90)
    print(f"TOP 30 UNIVERSAL FORMULAS (n={n})")
    print("=" * 90)

    seen = set()
    count = 0
    for r in results:
        key = r['formula'][:60]
        if key in seen:
            continue
        seen.add(key)
        count += 1
        if count > 30:
            break

        regime = regime_r2(rows, sigma_net, r['pred'])
        tau_gt = regime.get('τ>2.5', (0, -99, 999))
        tau_mid = regime.get('1.5<τ≤2.5', (0, -99, 999))
        tau_lo = regime.get('τ≤1.5', (0, -99, 999))
        fp = r.get('n_free', 1)

        # Quality flags
        universal = '★' if tau_gt[1] > 0.5 else '⚠' if tau_gt[1] > 0 else '✗'
        balanced = '★' if min(tau_lo[1], tau_mid[1]) > 0.85 else ' '

        print(f"\n#{count} {universal}{balanced} R²={r['r2']:.4f}  (Part {r['part']}, {fp}p, C={r['C']:.4f})")
        print(f"  {r['formula']}")
        for label, (nn, r2, err) in regime.items():
            flag = '✓' if r2 > 0.85 else '⚠' if r2 > 0.5 else '✗'
            print(f"  {flag} {label:12s}: n={nn:2d}  R²={r2:.3f}  |err|={err:.1f}%")

    # ═══════════════════════════════════════════════
    # SPECIAL: τ>2.5 최적 모델 (thin100 특화)
    # ═══════════════════════════════════════════════
    print("\n" + "=" * 90)
    print("THIN100 SPECIALISTS: τ>2.5에서 R² 최고인 모델 TOP 10")
    print("=" * 90)

    for r in results:
        regime = regime_r2(rows, sigma_net, r['pred'])
        tau_gt = regime.get('τ>2.5', (0, -99, 999))
        r['tau_gt_r2'] = tau_gt[1]

    thin_sorted = sorted(results, key=lambda x: -x.get('tau_gt_r2', -99))
    seen2 = set()
    count2 = 0
    for r in thin_sorted:
        key = r['formula'][:60]
        if key in seen2:
            continue
        seen2.add(key)
        count2 += 1
        if count2 > 10:
            break

        regime = regime_r2(rows, sigma_net, r['pred'])
        print(f"\n#{count2}  τ>2.5 R²={r.get('tau_gt_r2', -99):.3f}  |  ALL R²={r['r2']:.4f}  (Part {r['part']})")
        print(f"  {r['formula']}")
        for label, (nn, r2, err) in regime.items():
            print(f"    {label:12s}: n={nn:2d}  R²={r2:.3f}  |err|={err:.1f}%")


if __name__ == '__main__':
    main()
