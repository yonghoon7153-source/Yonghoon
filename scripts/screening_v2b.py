"""
v2.0b Focused Screening: φ_AM 추가 + Fixed Exponent 테스트
==========================================================
Goal: 현재 champion에 φ_AM을 추가하여 thin100까지 커버하는 최소 모델 찾기
- Overfitting 방지: fixed exponents (정수/반정수)로 테스트
- 공정 비교: 모든 R²를 log space에서 계산
- Per-regime breakdown: τ≤1.5 / 1.5<τ≤2.5 / τ>2.5
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
            if tau <= 0 or phi_se <= 0 or cn <= 0 or gb_d <= 0 or T <= 0 or g_path <= 0:
                continue
            sigma_brug = SIGMA_GRAIN * phi_se * f_perc / tau**2
            rows.append({
                'name': met_path.parent.name,
                'sigma_net': sigma_net, 'sigma_brug': sigma_brug,
                'phi_se': phi_se, 'phi_am': max(phi_am, 0.01),
                'tau': tau, 'f_perc': f_perc, 'cn': cn,
                'gb_d': gb_d, 'g_path': g_path,
                'hop_area': max(hop_area, 1e-6), 'T': T,
            })
    # Dedup
    seen = set()
    unique = []
    for r in rows:
        key = f"{r['phi_se']:.4f}_{r['T']:.1f}_{r['tau']:.3f}"
        if key not in seen:
            seen.add(key)
            unique.append(r)
    return unique


def r2_log(actual, predicted):
    """R² in log space (fair for data spanning orders of magnitude)."""
    la, lp = np.log(actual), np.log(predicted)
    ss_res = np.sum((la - lp)**2)
    ss_tot = np.sum((la - np.mean(la))**2)
    return 1 - ss_res / ss_tot if ss_tot > 0 else 0


def r2_lin(actual, predicted):
    """R² in linear space."""
    ss_res = np.sum((actual - predicted)**2)
    ss_tot = np.sum((actual - np.mean(actual))**2)
    return 1 - ss_res / ss_tot if ss_tot > 0 else 0


def fit_C_log(actual, rhs):
    """Fit C in log space: C = exp(mean(log(actual/rhs)))."""
    return float(np.exp(np.mean(np.log(actual / rhs))))


def per_regime(rows, actual, predicted):
    """Per-regime R² and error."""
    tau = np.array([r['tau'] for r in rows])
    regimes = [
        ('ALL', np.ones(len(rows), dtype=bool)),
        ('τ≤1.5', tau <= 1.5),
        ('1.5<τ≤2.5', (tau > 1.5) & (tau <= 2.5)),
        ('τ>2.5', tau > 2.5),
    ]
    results = {}
    for label, mask in regimes:
        n = mask.sum()
        if n < 2:
            continue
        r2l = r2_log(actual[mask], predicted[mask])
        r2n = r2_lin(actual[mask], predicted[mask])
        err = np.mean(np.abs(actual[mask] - predicted[mask]) / actual[mask]) * 100
        results[label] = {'n': n, 'r2_log': r2l, 'r2_lin': r2n, 'err': err}
    return results


def main():
    rows = load_data()
    n = len(rows)
    print(f"Loaded {n} cases (τ: {min(r['tau'] for r in rows):.2f}~{max(r['tau'] for r in rows):.2f})")
    print()

    sigma_net = np.array([r['sigma_net'] for r in rows])
    sigma_brug = np.array([r['sigma_brug'] for r in rows])
    phi_se = np.array([r['phi_se'] for r in rows])
    phi_am = np.array([r['phi_am'] for r in rows])
    tau = np.array([r['tau'] for r in rows])
    cn = np.array([r['cn'] for r in rows])
    gb_d = np.array([r['gb_d'] for r in rows])
    g_path = np.array([r['g_path'] for r in rows])
    hop_area = np.array([r['hop_area'] for r in rows])
    T = np.array([r['T'] for r in rows])
    f_perc = np.array([r['f_perc'] for r in rows])

    results = []

    # ═══════════════════════════════════════════════════════════
    # STRATEGY A: Current champion + φ_AM correction
    # σ = σ_brug × C × (G_path × GB_d²)^a × CN^b × φ_AM^c
    # ═══════════════════════════════════════════════════════════
    print("Strategy A: σ_brug × C × (G_path×GB_d²)^a × CN^b × φ_AM^c")
    print("-" * 70)

    exponent_grid_a = [0, 1/4, 1/3, 1/2, 2/3, 3/4, 1]
    exponent_grid_b = [0, 1, 3/2, 2, 5/2, 3]
    exponent_grid_c = [-2, -3/2, -1, -1/2, 0, 1/2, 1, 3/2, 2]

    for a, b, c in product(exponent_grid_a, exponent_grid_b, exponent_grid_c):
        if a == 0 and b == 0 and c == 0:
            continue
        rhs = sigma_brug * (g_path * gb_d**2)**a * cn**b * phi_am**c
        if np.any(rhs <= 0) or np.any(~np.isfinite(rhs)):
            continue
        C = fit_C_log(sigma_net, rhs)
        pred = C * rhs
        r2l = r2_log(sigma_net, pred)
        r2n = r2_lin(sigma_net, pred)
        results.append({
            'strategy': 'A',
            'formula': f'σ_brug × C × (G×d²)^{a:.2f} × CN^{b:.1f} × φ_AM^{c:.1f}',
            'short': f'a={a:.2f},b={b:.1f},c={c:.1f}',
            'r2_log': r2l, 'r2_lin': r2n, 'C': C,
            'n_free': 1,  # only C is free
            'pred': pred,
            'exponents': (a, b, c),
        })

    # ═══════════════════════════════════════════════════════════
    # STRATEGY B: φ_SE × φ_AM direct (no σ_brug)
    # σ = C × φ_SE^a × φ_AM^b × (contact_var)^c
    # ═══════════════════════════════════════════════════════════
    print("Strategy B: C × φ_SE^a × φ_AM^b × (contact)^c")
    print("-" * 70)

    contact_vars = {
        'CN': cn, 'GB_d': gb_d, 'G_path': g_path,
        'G×d²': g_path * gb_d**2, 'CN/τ': cn / tau,
        'τ×GB_d': tau * gb_d,
    }

    for cname, cvar in contact_vars.items():
        for a in [3/2, 2, 5/2, 3, 7/2, 4]:
            for b in [0, 1, 3/2, 2, 5/2, 3, 4]:
                for c in [-1, -1/2, 0, 1/4, 1/2, 3/4, 1]:
                    rhs = phi_se**a * phi_am**b * cvar**c
                    if np.any(rhs <= 0) or np.any(~np.isfinite(rhs)):
                        continue
                    C = fit_C_log(sigma_net, rhs)
                    pred = C * rhs
                    r2l = r2_log(sigma_net, pred)
                    results.append({
                        'strategy': 'B',
                        'formula': f'C × φ_SE^{a:.1f} × φ_AM^{b:.1f} × {cname}^{c:.2f}',
                        'short': f'φ_SE^{a:.1f}×φ_AM^{b:.1f}×{cname}^{c:.2f}',
                        'r2_log': r2l, 'r2_lin': r2_lin(sigma_net, pred), 'C': C,
                        'n_free': 1,
                        'pred': pred,
                    })

    # ═══════════════════════════════════════════════════════════
    # STRATEGY C: σ_brug × C × τ^a × (contact)^b × φ_AM^c
    # Extra τ correction beyond τ² in σ_brug
    # ═══════════════════════════════════════════════════════════
    print("Strategy C: σ_brug × C × τ^a × (contact)^b × φ_AM^c (extra τ correction)")
    print("-" * 70)

    for cname, cvar in {'G×d²': g_path * gb_d**2, 'CN': cn, 'CN/τ': cn/tau}.items():
        for a in [-3, -2, -3/2, -1, -1/2, 0, 1/2, 1]:
            for b in [0, 1/4, 1/2, 3/4, 1]:
                for c in [-2, -3/2, -1, -1/2, 0, 1/2, 1, 2]:
                    rhs = sigma_brug * tau**a * cvar**b * phi_am**c
                    if np.any(rhs <= 0) or np.any(~np.isfinite(rhs)):
                        continue
                    C = fit_C_log(sigma_net, rhs)
                    pred = C * rhs
                    r2l = r2_log(sigma_net, pred)
                    results.append({
                        'strategy': 'C',
                        'formula': f'σ_brug × C × τ^{a:.1f} × {cname}^{b:.2f} × φ_AM^{c:.1f}',
                        'short': f'τ^{a:.1f}×{cname}^{b:.2f}×φ_AM^{c:.1f}',
                        'r2_log': r2l, 'r2_lin': r2_lin(sigma_net, pred), 'C': C,
                        'n_free': 1,
                        'pred': pred,
                    })

    # ═══════════════════════════════════════════════════════════
    # STRATEGY D: Free fit (2 free exponents + C)
    # σ = σ_brug × C × (G_path × GB_d²)^a × CN^b
    # with a, b as free parameters (not fixed)
    # ═══════════════════════════════════════════════════════════
    print("Strategy D: Free-fit exponents (baseline)")
    print("-" * 70)

    # Free fit: log(σ/σ_brug) = a*log(G×d²) + b*log(CN) + log(C)
    log_ratio = np.log(sigma_net / sigma_brug)
    X_free = np.column_stack([np.log(g_path * gb_d**2), np.log(cn), np.ones(n)])
    coefs, _, _, _ = np.linalg.lstsq(X_free, log_ratio, rcond=None)
    pred_free = sigma_brug * np.exp(X_free @ coefs)
    r2_free = r2_log(sigma_net, pred_free)
    results.append({
        'strategy': 'D_free_2var',
        'formula': f'σ_brug × {np.exp(coefs[2]):.4f} × (G×d²)^{coefs[0]:.3f} × CN^{coefs[1]:.3f}',
        'short': f'FREE: (G×d²)^{coefs[0]:.3f}×CN^{coefs[1]:.3f}',
        'r2_log': r2_free, 'r2_lin': r2_lin(sigma_net, pred_free),
        'C': np.exp(coefs[2]), 'n_free': 3, 'pred': pred_free,
    })

    # Free fit with φ_AM: log(σ/σ_brug) = a*log(G×d²) + b*log(CN) + c*log(φ_AM) + log(C)
    X_free2 = np.column_stack([np.log(g_path * gb_d**2), np.log(cn), np.log(phi_am), np.ones(n)])
    coefs2, _, _, _ = np.linalg.lstsq(X_free2, log_ratio, rcond=None)
    pred_free2 = sigma_brug * np.exp(X_free2 @ coefs2)
    r2_free2 = r2_log(sigma_net, pred_free2)
    results.append({
        'strategy': 'D_free_3var',
        'formula': f'σ_brug × {np.exp(coefs2[3]):.4f} × (G×d²)^{coefs2[0]:.3f} × CN^{coefs2[1]:.3f} × φ_AM^{coefs2[2]:.3f}',
        'short': f'FREE: (G×d²)^{coefs2[0]:.3f}×CN^{coefs2[1]:.3f}×φ_AM^{coefs2[2]:.3f}',
        'r2_log': r2_free2, 'r2_lin': r2_lin(sigma_net, pred_free2),
        'C': np.exp(coefs2[3]), 'n_free': 4, 'pred': pred_free2,
    })

    # Free fit: φ_SE^a × φ_AM^b × CN^c (no σ_brug)
    X_free3 = np.column_stack([np.log(phi_se), np.log(phi_am), np.log(cn), np.ones(n)])
    coefs3, _, _, _ = np.linalg.lstsq(X_free3, np.log(sigma_net), rcond=None)
    pred_free3 = np.exp(X_free3 @ coefs3)
    r2_free3 = r2_log(sigma_net, pred_free3)
    results.append({
        'strategy': 'D_free_direct',
        'formula': f'{np.exp(coefs3[3]):.4f} × φ_SE^{coefs3[0]:.3f} × φ_AM^{coefs3[1]:.3f} × CN^{coefs3[2]:.3f}',
        'short': f'FREE DIRECT: φ_SE^{coefs3[0]:.2f}×φ_AM^{coefs3[1]:.2f}×CN^{coefs3[2]:.2f}',
        'r2_log': r2_free3, 'r2_lin': r2_lin(sigma_net, pred_free3),
        'C': np.exp(coefs3[3]), 'n_free': 4, 'pred': pred_free3,
    })

    # ═══════════════════════════════════════════════════════════
    # RESULTS
    # ═══════════════════════════════════════════════════════════
    results.sort(key=lambda x: -x['r2_log'])

    # Current champion baseline
    rhs_champ = sigma_brug * (g_path * gb_d**2)**0.25 * cn**2
    C_champ = fit_C_log(sigma_net, rhs_champ)
    pred_champ = C_champ * rhs_champ

    print("\n" + "=" * 80)
    print(f"CURRENT CHAMPION: σ_brug × C × (G×d²)^(1/4) × CN²")
    print(f"  C={C_champ:.4f}")
    regime = per_regime(rows, sigma_net, pred_champ)
    for label, r in regime.items():
        print(f"  {label:12s}: n={r['n']:2d}  R²_log={r['r2_log']:.3f}  R²_lin={r['r2_lin']:.3f}  |err|={r['err']:.1f}%")

    print("\n" + "=" * 80)
    print(f"TOP 30 CANDIDATES (n={n}, log-space R²)")
    print("=" * 80)

    seen = set()
    count = 0
    for r in results:
        if r['r2_log'] < 0.7:
            break
        key = r.get('short', r['formula'])
        if key in seen:
            continue
        seen.add(key)
        count += 1
        if count > 30:
            break

        print(f"\n#{count}  R²_log={r['r2_log']:.4f}  R²_lin={r['r2_lin']:.4f}  ({r['strategy']}, {r['n_free']}p, C={r['C']:.4f})")
        print(f"  {r['formula']}")

        # Per-regime
        regime = per_regime(rows, sigma_net, r['pred'])
        for label, rr in regime.items():
            flag = '⚠' if rr['r2_log'] < 0.5 else '✓' if rr['r2_log'] > 0.85 else ' '
            print(f"  {flag} {label:12s}: n={rr['n']:2d}  R²_log={rr['r2_log']:.3f}  |err|={rr['err']:.1f}%")

    # ═══════════════════════════════════════════════════════════
    # BEST FIXED-EXPONENT (1 free param) COMPARISON
    # ═══════════════════════════════════════════════════════════
    print("\n" + "=" * 80)
    print("BEST FIXED-EXPONENT MODELS (1 free param = C only)")
    print("=" * 80)

    fixed_only = [r for r in results if r['n_free'] == 1]
    fixed_only.sort(key=lambda x: -x['r2_log'])
    seen2 = set()
    for i, r in enumerate(fixed_only[:15]):
        key = r.get('short', r['formula'])
        if key in seen2:
            continue
        seen2.add(key)
        regime = per_regime(rows, sigma_net, r['pred'])
        tau_gt25 = regime.get('τ>2.5', {})
        tau_ok = tau_gt25.get('r2_log', -99) if tau_gt25 else -99
        print(f"\n  {r['formula']}")
        print(f"    R²_log={r['r2_log']:.4f}  C={r['C']:.4f}  τ>2.5: R²={tau_ok:.3f}")


if __name__ == '__main__':
    main()
