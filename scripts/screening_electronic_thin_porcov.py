"""
Thin: por^a × cov^b / (φ_SE^c × ξ^d) 개형 최적화
Spearman ρ + R² 동시 평가
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
            cn = max(m.get('am_am_cn', 0.01), 0.01)
            ep = max(m.get('electronic_percolating_fraction', 0), 0.01)
            T = m.get('thickness_um', 0)
            r_am = max(m.get('r_AM_P', 0), m.get('r_AM_S', 0))
            d_am = r_am * 2 if r_am > 0.1 else 5.0
            pa = max(m.get('phi_am', 0), 0.01)
            por = max(m.get('porosity', 0), 0.1)
            cov = max(m.get('coverage_AM_P_mean', m.get('coverage_AM_S_mean', m.get('coverage_AM_mean', 20))), 0.1) / 100
            ps = max(m.get('phi_se', 0), 0.01)
            am_se = m.get('am_se_ratio', '')
            if cn <= 0 or T <= 0 or d_am <= 0: continue
            xi = T/d_am
            if xi >= 10 or ep < 0.65: continue
            rows.append({'sel':sel,'por':por,'cov':cov,'ps':ps,'xi':xi,'pa':pa,'cn':cn,'ep':ep,
                         'group': am_se if am_se else f"T{T:.0f}"})
    seen = set(); u = []
    for r in rows:
        k = f"{r['pa']:.4f}_{r['xi']:.1f}"
        if k not in seen: seen.add(k); u.append(r)
    return u

def main():
    rows = load_data()
    if not rows: print("데이터 없음!"); return
    sel = np.array([r['sel'] for r in rows])
    por = np.array([r['por'] for r in rows])
    cov = np.array([r['cov'] for r in rows])
    ps = np.array([r['ps'] for r in rows])
    xi = np.array([r['xi'] for r in rows])
    pa = np.array([r['pa'] for r in rows])
    groups = [r['group'] for r in rows]
    n = len(rows)
    print(f"n={n}")

    # Group masks
    gm = {}
    for g in sorted(set(groups)):
        mask = np.array([groups[i] == g for i in range(n)])
        if mask.sum() >= 3: gm[g] = mask
    print(f"Groups: {list(gm.keys())}")

    def r2l(a, p):
        la, lp = np.log(a), np.log(p)
        return 1 - np.sum((la-lp)**2) / np.sum((la-np.mean(la))**2)
    def fitC(s, r):
        v = (r>0)&np.isfinite(r)
        return float(np.exp(np.mean(np.log(s[v]/r[v])))) if v.sum()>=3 else None
    def avg_spearman(vals):
        rhos = []
        for g, mask in gm.items():
            sg = sel[mask]; fg = vals[mask]
            if len(sg) >= 3 and np.std(fg) > 1e-10:
                rho, _ = spearmanr(fg, sg)
                rhos.append(rho)
        return np.mean(rhos) if rhos else -999, min(rhos) if rhos else -999, rhos

    # ═══════════════════════════════════════════
    # por^a × cov^b / (φ_SE^c × ξ^d)
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("por^a × cov^b / (φ_SE^c × ξ^d) — 개형+R² 동시")
    print("="*80)
    results = []
    for a in np.arange(0.1, 1.5, 0.1):
        for b in np.arange(0.1, 1.5, 0.1):
            for c in [0.5, 0.75, 1, 1.25, 1.5]:
                for d in [0.25, 0.5]:
                    rhs = SAM * por**a * cov**b / (ps**c * xi**d)
                    C = fitC(sel, rhs)
                    if C is None: continue
                    pred = C * rhs
                    r2 = r2l(sel, pred)
                    mean_rho, min_rho, rhos = avg_spearman(rhs)
                    if mean_rho > 0.3 and r2 > -5:
                        errs = (pred - sel) / sel * 100
                        results.append({
                            'r2': r2, 'rho': mean_rho, 'min_rho': min_rho,
                            'range': f"{errs.min():+.0f}~{errs.max():+.0f}%",
                            'rng': max(errs) - min(errs),
                            'label': f"por^{a:.1f}×cov^{b:.1f}/φ_SE^{c}×ξ^{d}",
                            'rhos': rhos
                        })

    # Sort by combined score: ρ × (1 + R²)
    results.sort(key=lambda x: -(x['rho'] * max(0.01, 1 + x['r2'])))
    print(f"\n  Top 20 (ρ×(1+R²) 순서 — 개형+level 균형):")
    print(f"  {'Formula':45s} {'ρ':>6s} {'R²':>7s} {'range':>12s} {'per-group ρ':>20s}")
    print("  "+"-"*95)
    for i, r in enumerate(results[:20]):
        rho_str = ', '.join(f'{rh:+.2f}' for rh in r['rhos'])
        flag = '★' if r['rho'] > 0.6 and r['r2'] > 0 else '●' if r['rho'] > 0.5 else ' '
        print(f"  {flag}{r['label']:44s} {r['rho']:+.3f} {r['r2']:7.3f} {r['range']:>12s} [{rho_str}]")

    # Also sort by R² only
    results.sort(key=lambda x: -x['r2'])
    print(f"\n  Top 10 (R² 순서):")
    for i, r in enumerate(results[:10]):
        rho_str = ', '.join(f'{rh:+.2f}' for rh in r['rhos'])
        print(f"  #{i+1} R²={r['r2']:7.3f} ρ={r['rho']:+.3f} {r['range']:>12s}  {r['label']}  [{rho_str}]")

    # Sort by range (개형 우선)
    results.sort(key=lambda x: x['rng'])
    print(f"\n  Top 10 (range 최소 — 개형 최우선):")
    for i, r in enumerate(results[:10]):
        rho_str = ', '.join(f'{rh:+.2f}' for rh in r['rhos'])
        print(f"  #{i+1} range={r['range']:>12s} R²={r['r2']:7.3f} ρ={r['rho']:+.3f}  {r['label']}  [{rho_str}]")

if __name__ == '__main__':
    main()
