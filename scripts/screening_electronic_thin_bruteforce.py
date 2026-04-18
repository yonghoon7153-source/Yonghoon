"""
Thin brute force: 모든 조합, 넓은 지수, φ_AM 포함
개형 우선 → signed error range 최소화
"""
import json, os, numpy as np, warnings
from pathlib import Path
from itertools import combinations
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
            cov = max(m.get('coverage_AM_P_mean', m.get('coverage_AM_S_mean', m.get('coverage_AM_mean', 20))), 0.1) / 100
            r_am = max(m.get('r_AM_P', 0), m.get('r_AM_S', 0))
            d_am = r_am * 2 if r_am > 0.1 else 5.0
            por = m.get('porosity', 0)
            ep = m.get('electronic_percolating_fraction', 0)
            ps = m.get('phi_se', 0)
            tau = max(m.get('tortuosity_recommended', m.get('tortuosity_mean', 0)), 0.1)
            if pa <= 0 or cn <= 0 or T <= 0: continue
            rows.append({
                'sel': sel, 'pa': pa, 'cn': cn, 'xi': T/d_am, 'por': max(por, 0.1),
                'ep': max(ep, 0.01), 'cov': cov, 'ps': max(ps, 0.01), 'tau': tau,
                'name': mp.parent.name
            })
    seen = set(); u = []
    for r in rows:
        k = f"{r['pa']:.4f}_{r['xi']:.1f}"
        if k not in seen: seen.add(k); u.append(r)
    return u

def r2l(a, p):
    la, lp = np.log(a), np.log(p); return 1-np.sum((la-lp)**2)/np.sum((la-np.mean(la))**2)
def fitC(s, r):
    v = (r>0)&np.isfinite(r); return float(np.exp(np.mean(np.log(s[v]/r[v])))) if v.sum()>=3 else None

def main():
    rows = load_data()
    if not rows: print("데이터 없음!"); return
    V = {k: np.array([r[k] for r in rows]) for k in rows[0] if k != 'name'}
    names = [r['name'] for r in rows]
    thin = V['xi'] < 10; sn = V['sel'][thin]
    T = {k: V[k][thin] for k in V}
    ok = T['ep'] >= 0.65; sn_ok = sn[ok]; n_ok = ok.sum()
    print(f"thin={thin.sum()}, perc={n_ok}")

    # ═══════════════════════════════════════════
    # 1. Free regression with φ_AM (3~5 vars)
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print(f"1. Free regression (ep≥0.65, n={n_ok})")
    print("="*80)
    var_pool = {
        'ep': np.log(T['ep'][ok]), 'φ': np.log(T['pa'][ok]),
        'CN': np.log(T['cn'][ok]), 'por': np.log(T['por'][ok]),
        'cov': np.log(T['cov'][ok]), 'ξ': np.log(T['xi'][ok]),
        'τ': np.log(T['tau'][ok]), 'φ_SE': np.log(T['ps'][ok]),
    }
    log_sn = np.log(sn_ok)
    results = []
    for nv in [3, 4, 5]:
        for vn in combinations(var_pool.keys(), nv):
            X = np.column_stack([var_pool[v] for v in vn] + [np.ones(n_ok)])
            try:
                coefs, _, _, _ = np.linalg.lstsq(X, log_sn, rcond=None)
                pred = np.exp(X @ coefs)
                r2 = r2l(sn_ok, pred)
                if r2 > 0.88:
                    errs = (pred - sn_ok) / sn_ok * 100
                    err_range = max(errs) - min(errs)
                    results.append({'r2': r2, 'coefs': coefs, 'vnames': vn,
                                    'range': f"{errs.min():+.0f}~{errs.max():+.0f}%",
                                    'err_range': err_range})
            except: pass
    # Sort by signed error range (개형 우선!)
    results.sort(key=lambda x: x['err_range'])
    print(f"\n  Top 20 (sorted by error RANGE — 개형 우선):")
    for i, r in enumerate(results[:20]):
        terms = ' × '.join(f"{v}^{r['coefs'][j]:.2f}" for j, v in enumerate(r['vnames']))
        print(f"  #{i+1} range={r['range']:>12s} R²={r['r2']:.4f}  ({len(r['vnames'])}v) {terms}")

    # Also sort by R²
    results.sort(key=lambda x: -x['r2'])
    print(f"\n  Top 10 (sorted by R²):")
    for i, r in enumerate(results[:10]):
        terms = ' × '.join(f"{v}^{r['coefs'][j]:.2f}" for j, v in enumerate(r['vnames']))
        print(f"  #{i+1} R²={r['r2']:.4f} {r['range']:>12s}  ({len(r['vnames'])}v) {terms}")

    # ═══════════════════════════════════════════
    # 2. φ_AM 포함 grid search
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("2. φ_AM 포함 grid search (개형 우선)")
    print("="*80)
    c2 = []
    for pa_e in [-2, -1, 0, 1, 2, 3]:
        for cn_e in [0.5, 1, 1.5]:
            for por_e in [0, 0.5, 1, 1.5, 2]:
                for cov_e in [0, 0.5, 1]:
                    for td_e in [0.25, 0.5]:
                        rhs = SAM * T['ep'][ok]**0.5 * T['pa'][ok]**pa_e * T['cn'][ok]**cn_e * T['por'][ok]**por_e * T['cov'][ok]**cov_e / T['xi'][ok]**td_e
                        C = fitC(sn_ok, rhs)
                        if C is None or C <= 0: continue
                        pred = C * rhs
                        r2 = r2l(sn_ok, pred)
                        if r2 > 0.85:
                            errs = (pred - sn_ok) / sn_ok * 100
                            err_range = max(errs) - min(errs)
                            ps = f'φ^{pa_e}×' if pa_e else ''
                            cs = f'×cov^{cov_e}' if cov_e else ''
                            pors = f'×por^{por_e}' if por_e else ''
                            c2.append({'r2':r2, 'err_range':err_range,
                                'range':f'{errs.min():+.0f}~{errs.max():+.0f}%',
                                'label':f'√ep×{ps}CN^{cn_e}{pors}{cs}/(T/d)^{td_e}'})
    # Sort by error range
    c2.sort(key=lambda x: x['err_range'])
    print(f"\n  Top 20 (개형 우선 — error range 최소):")
    for i, r in enumerate(c2[:20]):
        print(f"  #{i+1} range={r['range']:>12s} R²={r['r2']:.4f}  {r['label']}")
    c2.sort(key=lambda x: -x['r2'])
    print(f"\n  Top 10 (R² 우선):")
    for i, r in enumerate(c2[:10]):
        print(f"  #{i+1} R²={r['r2']:.4f} {r['range']:>12s}  {r['label']}")

    # ═══════════════════════════════════════════
    # 3. τ 포함 (thick에서는 √τ가 효과적)
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("3. τ 포함")
    print("="*80)
    c3 = []
    for cn_e in [0.5, 1, 1.5]:
        for por_e in [0.5, 1, 1.5]:
            for tau_e in [0.25, 0.5, 0.75, 1]:
                for cov_e in [0, 0.5]:
                    for td_e in [0.25, 0.5]:
                        rhs = SAM * T['ep'][ok]**0.5 * T['cn'][ok]**cn_e * T['por'][ok]**por_e * T['tau'][ok]**tau_e * T['cov'][ok]**cov_e / T['xi'][ok]**td_e
                        C = fitC(sn_ok, rhs)
                        if C is None: continue
                        pred = C * rhs; r2 = r2l(sn_ok, pred)
                        if r2 > 0.87:
                            errs = (pred - sn_ok) / sn_ok * 100
                            err_range = max(errs) - min(errs)
                            cs = f'×cov^{cov_e}' if cov_e else ''
                            c3.append({'r2':r2, 'err_range':err_range,
                                'range':f'{errs.min():+.0f}~{errs.max():+.0f}%',
                                'label':f'√ep×CN^{cn_e}×por^{por_e}×τ^{tau_e}{cs}/(T/d)^{td_e}'})
    c3.sort(key=lambda x: x['err_range'])
    print(f"\n  Top 15 (개형 우선):")
    for i, r in enumerate(c3[:15]):
        print(f"  #{i+1} range={r['range']:>12s} R²={r['r2']:.4f}  {r['label']}")

if __name__ == '__main__':
    main()
