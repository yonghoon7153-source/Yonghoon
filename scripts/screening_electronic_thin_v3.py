"""
Thin v3: f_p×CN×por×√cov/√ξ 기반 더 파기
- el_perc 지수 변형 (^½, ^1, ^2, (ep-0.5)^a)
- 5th variable 추가
- Free regression with el_perc
- 개형 우선 (signed error range 최소화)
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
            tau = max(m.get('tortuosity_recommended', m.get('tortuosity_mean', 0)), 0.1)
            cov = max(m.get('coverage_AM_P_mean', m.get('coverage_AM_S_mean', m.get('coverage_AM_mean', 20))), 0.1) / 100
            r_am = max(m.get('r_AM_P', 0), m.get('r_AM_S', 0))
            d_am = r_am * 2 if r_am > 0.1 else 5.0
            por = m.get('porosity', 0)
            ep = m.get('electronic_percolating_fraction', 0)
            ps = m.get('phi_se', 0)
            d = m.get('am_am_mean_delta', 0)
            A = m.get('am_am_mean_area', 0)
            hop = max(m.get('am_am_mean_hop', 0), 0.1)
            cn_std = m.get('am_am_cn_std', 0)
            if pa <= 0 or cn <= 0 or T <= 0: continue
            rows.append({
                'sel': sel, 'pa': pa, 'cn': cn, 'tau': tau, 'cov': cov,
                'xi': T/d_am, 'por': max(por, 0.1), 'ep': max(ep, 0.01),
                'ps': max(ps, 0.01), 'hop': hop,
                'd2a': d**2/A if A > 0 and d > 0 else 0.001,
                'cn_std': max(cn_std, 0.01),
                'name': mp.parent.name
            })
    seen = set(); u = []
    for r in rows:
        k = f"{r['pa']:.4f}_{r['xi']:.1f}"
        if k not in seen: seen.add(k); u.append(r)
    return u

def r2l(a, p):
    la, lp = np.log(a), np.log(p); return 1-np.sum((la-lp)**2)/np.sum((la-np.mean(la))**2)

def main():
    rows = load_data()
    if not rows: print("데이터 없음!"); return
    V = {k: np.array([r[k] for r in rows]) for k in rows[0] if k != 'name'}
    names = [r['name'] for r in rows]
    n = len(rows); thin = V['xi'] < 10; sn = V['sel'][thin]
    T = {k: V[k][thin] for k in V}
    n_thin = thin.sum()
    ok = T['ep'] >= 0.85; sn_ok = sn[ok]; n_ok = ok.sum()
    tn_names = [names[i] for i in np.where(thin)[0]]
    print(f"n={n}, thin={n_thin}, perc(≥0.85)={n_ok}")

    def fitC(s, r):
        v = (r>0)&np.isfinite(r)
        return float(np.exp(np.mean(np.log(s[v]/r[v])))) if v.sum()>=3 else None
    def eval_formula(label, rhs_all):
        rhs_ok = rhs_all[ok]
        C = fitC(sn_ok, rhs_ok)
        if C is None: return
        pred_ok = C * rhs_ok; pred_all = C * rhs_all
        r2_ok = r2l(sn_ok, pred_ok)
        errs_ok = (pred_ok - sn_ok) / sn_ok * 100
        err_range = f"{errs_ok.min():+.0f}~{errs_ok.max():+.0f}%"
        r2_all = r2l(sn, pred_all)
        return {'r2_ok': r2_ok, 'r2_all': r2_all, 'range': err_range,
                'rhs': rhs_all, 'C': C, 'label': label}

    # ═══════════════════════════════════════════
    # 1. el_perc 지수 변형
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("1. el_perc 지수 + (ep-threshold) 형태")
    print("="*80)
    base = T['cn'] * T['por'] * T['cov']**0.5 / T['xi']**0.5
    for label, ep_factor in [
        ('ep^0.5', T['ep']**0.5), ('ep^1', T['ep']),
        ('ep^1.5', T['ep']**1.5), ('ep^2', T['ep']**2),
        ('ep^3', T['ep']**3),
        ('(ep-0.5)^1', np.clip(T['ep']-0.5,0.001,None)),
        ('(ep-0.5)^2', np.clip(T['ep']-0.5,0.001,None)**2),
        ('(ep-0.6)^1', np.clip(T['ep']-0.6,0.001,None)),
        ('(ep-0.6)^2', np.clip(T['ep']-0.6,0.001,None)**2),
        ('(ep-0.7)^1', np.clip(T['ep']-0.7,0.001,None)),
        ('no ep', np.ones(n_thin)),
    ]:
        rhs = SAM * ep_factor * base
        r = eval_formula(f'{label}×CN×por×√cov/√ξ', rhs)
        if r:
            flag = '★' if r['r2_ok'] > 0.90 else '●' if r['r2_ok'] > 0.88 else ' '
            print(f"  {flag}{r['label']:45s} R²(p)={r['r2_ok']:.4f} all={r['r2_all']:.3f} {r['range']}")

    # ═══════════════════════════════════════════
    # 2. 5th variable 추가 (base = ep×CN×por×√cov/√ξ)
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("2. + 5th variable")
    print("="*80)
    base_ep = SAM * T['ep'] * T['cn'] * T['por'] * T['cov']**0.5 / T['xi']**0.5
    for xn, xa in [('×√τ', T['tau']**0.5), ('×τ', T['tau']),
                    ('×φ', T['pa']), ('×φ^-1', T['pa']**-1),
                    ('×√(d2a)', T['d2a']**0.5), ('×(d2a)^¼', T['d2a']**0.25),
                    ('×φ_SE', T['ps']), ('×√φ_SE', T['ps']**0.5),
                    ('×por^0.5', T['por']**0.5),
                    ('×cov^0.5', T['cov']**0.5),
                    ('×hop^¼', T['hop']**0.25), ('×hop^-¼', T['hop']**-0.25),
                    ('×CN_std^-½', T['cn_std']**-0.5),
                    ('×CN^¼', T['cn']**0.25),
                    ('/por^0.5', T['por']**-0.5)]:
        rhs = base_ep * xa
        r = eval_formula(f'ep×CN×por×√cov/√ξ{xn}', rhs)
        if r and r['r2_ok'] > 0.88:
            flag = '★' if r['r2_ok'] > 0.92 else '●' if r['r2_ok'] > 0.90 else ' '
            print(f"  {flag}{r['label']:50s} R²(p)={r['r2_ok']:.4f} all={r['r2_all']:.3f} {r['range']}")

    # ═══════════════════════════════════════════
    # 3. Free regression (ep 포함, perc only)
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print(f"3. Free regression (perc only, n={n_ok})")
    print("="*80)
    var_pool = {
        'ep': np.log(T['ep'][ok]), 'CN': np.log(T['cn'][ok]),
        'por': np.log(T['por'][ok]), 'cov': np.log(T['cov'][ok]),
        'ξ': np.log(T['xi'][ok]), 'τ': np.log(T['tau'][ok]),
        'φ': np.log(T['pa'][ok]), 'φ_SE': np.log(T['ps'][ok]),
        'd2a': np.log(T['d2a'][ok]), 'hop': np.log(T['hop'][ok]),
    }
    log_sn = np.log(sn_ok)
    results = []
    for nv in [3, 4, 5]:
        for vn in combinations(var_pool.keys(), nv):
            if 'ep' not in vn: continue  # ep 필수
            X = np.column_stack([var_pool[v] for v in vn] + [np.ones(n_ok)])
            try:
                coefs, _, _, _ = np.linalg.lstsq(X, log_sn, rcond=None)
                pred = np.exp(X @ coefs)
                r2 = r2l(sn_ok, pred)
                if r2 > 0.90:
                    errs = (pred - sn_ok) / sn_ok * 100
                    results.append({'r2': r2, 'coefs': coefs, 'vnames': vn,
                                    'range': f"{errs.min():+.0f}~{errs.max():+.0f}%"})
            except: pass
    results.sort(key=lambda x: -x['r2'])
    print(f"  Top 20:")
    for i, r in enumerate(results[:20]):
        terms = ' × '.join(f"{v}^{r['coefs'][j]:.2f}" for j, v in enumerate(r['vnames']))
        flag = '★' if r['r2'] > 0.94 else '●' if r['r2'] > 0.92 else ' '
        print(f"  {flag}#{i+1} R²={r['r2']:.4f} {r['range']:>12s}  ({len(r['vnames'])}v) {terms}")

    # ═══════════════════════════════════════════
    # 4. Best per-case
    # ═══════════════════════════════════════════
    if results:
        best = results[0]
        X = np.column_stack([var_pool[v] for v in best['vnames']] + [np.ones(n_ok)])
        pred = np.exp(X @ best['coefs'])
        errs = (pred - sn_ok) / sn_ok * 100
        ok_names = [tn_names[i] for i in np.where(ok)[0]]
        print(f"\n{'='*80}")
        print(f"4. Best free fit per-case")
        print("="*80)
        for j in np.argsort(-np.abs(errs)):
            print(f"  {ok_names[j]:35s} σ={sn_ok[j]:6.2f} pred={pred[j]:6.2f} err={errs[j]:+5.0f}% "
                  f"ep={T['ep'][ok][j]:.3f} CN={T['cn'][ok][j]:.2f}")

if __name__ == '__main__':
    main()
