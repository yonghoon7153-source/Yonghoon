"""
Thin 개형 최적화: intra-group Spearman correlation
각 formula의 그룹 내 순위 상관을 평가
"""
import json, os, numpy as np, warnings
from pathlib import Path
from scipy.stats import spearmanr
warnings.filterwarnings('ignore')
WEBAPP = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'webapp')
SAM = 50.0

def load_data():
    rows = []
    for base in [Path(WEBAPP)/'results', Path(WEBAPP)/'archive']:
        if not base.is_dir(): continue
        for mp in base.rglob('full_metrics.json'):
            try:
                with open(mp) as f: m = json.load(f)
            except: continue
            sel = m.get('electronic_sigma_full_mScm', 0)
            if not sel or sel < 0.001: continue
            pa = max(m.get('phi_am', 0), 0.01)
            cn = max(m.get('am_am_cn', 0.01), 0.01)
            T = m.get('thickness_um', 0)
            tau = max(m.get('tortuosity_recommended', m.get('tortuosity_mean', 0)), 0.1)
            cov = max(m.get('coverage_AM_P_mean', m.get('coverage_AM_S_mean', m.get('coverage_AM_mean', 20))), 0.1) / 100
            r_am = max(m.get('r_AM_P', 0), m.get('r_AM_S', 0))
            d_am = r_am * 2 if r_am > 0.1 else 5.0
            por = max(m.get('porosity', 0), 0.1)
            ep = max(m.get('electronic_percolating_fraction', 0), 0.01)
            ea = max(m.get('electronic_active_fraction', 0), 0.01)
            ps = m.get('phi_se', 0)
            A = max(m.get('am_am_mean_area', 0), 0.01)
            d = max(m.get('am_am_mean_delta', 0), 0.001)
            hop = max(m.get('am_am_mean_hop', 0), 0.1)
            am_se = m.get('am_se_ratio', '')
            if pa <= 0 or cn <= 0 or T <= 0 or d_am <= 0: continue
            xi = T/d_am
            if xi >= 10: continue
            rows.append({
                'sel': sel, 'pa': pa, 'cn': cn, 'tau': tau, 'cov': cov,
                'xi': xi, 'por': por, 'ep': ep, 'ea': ea, 'ps': max(ps, 0.01),
                'A': A, 'd': d, 'hop': hop, 'group': am_se if am_se else f"T{T:.0f}",
                'name': mp.parent.name
            })
    seen = set(); u = []
    for r in rows:
        k = f"{r['pa']:.4f}_{r['xi']:.1f}"
        if k not in seen: seen.add(k); u.append(r)
    return u

def main():
    rows = load_data()
    if not rows: print("데이터 없음!"); return
    V = {k: np.array([r[k] for r in rows]) for k in rows[0] if k not in ('name','group')}
    names = [r['name'] for r in rows]; groups = [r['group'] for r in rows]
    n = len(rows)
    print(f"thin n={n}")

    # Group indices
    unique_groups = sorted(set(groups))
    group_masks = {}
    for g in unique_groups:
        mask = np.array([groups[i] == g for i in range(n)])
        if mask.sum() >= 3:
            group_masks[g] = mask

    print(f"Groups with ≥3 cases: {list(group_masks.keys())} ({sum(m.sum() for m in group_masks.values())} cases)")

    # Formula candidates
    def make_formulas():
        cn = V['cn']; pa = V['pa']; por = V['por']; cov = V['cov']
        tau = V['tau']; ep = V['ep']; ea = V['ea']; xi = V['xi']
        A = V['A']; d = V['d']; hop = V['hop']; ps = V['ps']
        return [
            # Simple single variables
            ('CN', cn),
            ('por', por),
            ('cov', cov),
            ('τ', tau),
            ('φ', pa),
            ('ep', ep),
            ('ea', ea),
            ('A', A),
            ('δ', d),
            ('hop', hop),
            ('ep×ea', ep*ea),
            ('CN×A', cn*A),
            ('√(CN×A)', np.sqrt(cn*A)),
            ('CN/φ', cn/pa),
            # 2-variable combos
            ('CN×por', cn*por),
            ('CN×τ', cn*tau),
            ('CN×cov', cn*cov),
            ('CN×ep', cn*ep),
            ('CN×ea', cn*ea),
            ('por×τ', por*tau),
            ('por×cov', por*cov),
            ('φ×CN', pa*cn),
            ('φ²×CN', pa**2*cn),
            ('CN²', cn**2),
            # 3-variable combos
            ('CN×por×cov', cn*por*cov),
            ('CN×por×τ', cn*por*tau),
            ('CN×por×ep', cn*por*ep),
            ('CN×τ×cov', cn*tau*cov),
            ('CN×ep×ea', cn*ep*ea),
            ('φ⁴×CN^1.5×cov×√τ', pa**4*cn**1.5*cov*tau**0.5),
            ('φ²×CN×cov', pa**2*cn*cov),
            # With A/δ
            ('CN×√A×por', cn*np.sqrt(A)*por),
            ('√(CN×A)×por', np.sqrt(cn*A)*por),
            ('CN×δ', cn*d),
            ('√(CN×A)×τ', np.sqrt(cn*A)*tau),
            # el_perc combinations
            ('ep^2×CN', ep**2*cn),
            ('(ep×ea)^0.75×CN', (ep*ea)**0.75*cn),
            ('ep×CN×por', ep*cn*por),
            ('ep×ea×CN×por', ep*ea*cn*por),
            # Inverse φ
            ('CN×por/φ', cn*por/pa),
            ('CN×por/φ²', cn*por/pa**2),
            ('CN/φ×por', cn/pa*por),
            # With T/d
            ('CN/√ξ', cn/xi**0.5),
            ('CN×por/√ξ', cn*por/xi**0.5),
            ('CN×τ/√ξ', cn*tau/xi**0.5),
            ('ep×CN/√ξ', ep*cn/xi**0.5),
        ]

    formulas = make_formulas()

    # Evaluate each formula: average intra-group Spearman ρ
    print(f"\n{'='*80}")
    print("Intra-group Spearman ρ (개형 정확도)")
    print("="*80)
    results = []
    for label, vals in formulas:
        rhos = []
        for g, mask in group_masks.items():
            sg = V['sel'][mask]
            fg = vals[mask]
            if len(sg) >= 3 and np.std(fg) > 1e-10:
                rho, _ = spearmanr(fg, sg)
                rhos.append(rho)
        if len(rhos) >= 1:
            mean_rho = np.mean(rhos)
            min_rho = min(rhos)
            consistent = all(r > 0 for r in rhos) or all(r < 0 for r in rhos)
            results.append((mean_rho, min_rho, label, rhos, consistent))

    results.sort(key=lambda x: -x[0])
    print(f"\n  {'Formula':35s} {'mean ρ':>7s} {'min ρ':>7s} {'consistent':>10s} {'per-group':>20s}")
    print("  "+"-"*85)
    for mean_rho, min_rho, label, rhos, consistent in results[:30]:
        flag = '★★' if consistent and mean_rho > 0.7 else '★' if mean_rho > 0.5 else ' '
        rho_str = ', '.join(f'{r:+.2f}' for r in rhos)
        print(f"  {flag}{label:34s} {mean_rho:+.3f}  {min_rho:+.3f}  {'YES' if consistent else 'no':>10s}  [{rho_str}]")

    # Also show R² for comparison
    print(f"\n  R² 비교 (global thin fit):")
    for mean_rho, min_rho, label, rhos, consistent in results[:10]:
        idx = [i for i, (l, v) in enumerate(formulas) if l == label][0]
        vals = formulas[idx][1]
        rhs = SAM * vals / V['xi']**0.5
        from functools import reduce
        v = (rhs > 0) & np.isfinite(rhs)
        if v.sum() < 3: continue
        C = float(np.exp(np.mean(np.log(V['sel'][v] / rhs[v]))))
        pred = C * rhs
        la = np.log(V['sel']); lp = np.log(pred)
        r2 = 1 - np.sum((la-lp)**2)/np.sum((la-np.mean(la))**2)
        print(f"    {label:35s} ρ={mean_rho:+.3f} R²={r2:.4f}")

if __name__ == '__main__':
    main()
