"""
Thin 개형 분석: 그룹 내 P:S 변화에 따른 변수 추적
각 그룹에서 σ_el과 같은 방향으로 움직이는 변수가 뭔지
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
            n_contacts = m.get('am_am_n_contacts', 0)
            if pa <= 0 or cn <= 0 or T <= 0: continue
            # Group key: AM:SE ratio + thickness
            am_se = m.get('am_se_ratio', '')
            group = am_se if am_se else f"T{T:.0f}"
            rows.append({
                'sel': sel, 'pa': pa, 'cn': cn, 'tau': tau, 'cov': cov,
                'xi': T/d_am, 'por': max(por, 0.1), 'ep': max(ep, 0.01),
                'ps': max(ps, 0.01), 'd2a': d**2/A if A > 0 and d > 0 else 0.001,
                'hop': hop, 'cn_std': max(cn_std, 0.01),
                'n_c': max(n_contacts, 1), 'd': max(d, 0.001), 'A': max(A, 0.01),
                'cn_x_ep': cn * max(ep, 0.01),
                'cn_x_por': cn * max(por, 0.1),
                'group': group, 'name': mp.parent.name
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
    V = {k: np.array([r[k] for r in rows]) for k in rows[0] if k not in ('name','group')}
    names = [r['name'] for r in rows]; groups = [r['group'] for r in rows]
    thin = V['xi'] < 10; sn = V['sel'][thin]
    T = {k: V[k][thin] for k in V}
    tn_names = [names[i] for i in np.where(thin)[0]]
    tn_groups = [groups[i] for i in np.where(thin)[0]]
    n_thin = thin.sum()
    print(f"thin={n_thin}")

    # ═══════════════════════════════════════════
    # 1. 그룹별 변수 vs σ_el 상관 (intra-group)
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("1. Intra-group correlation (각 그룹 내에서)")
    print("="*80)

    vars_check = ['cn','por','cov','ep','ps','tau','d2a','hop','cn_std',
                   'pa','d','A','n_c','cn_x_ep','cn_x_por']
    # Pool intra-group correlations
    from collections import defaultdict
    corr_pool = defaultdict(list)

    unique_groups = sorted(set(tn_groups))
    for g in unique_groups:
        mask = np.array([tn_groups[i] == g for i in range(n_thin)])
        if mask.sum() < 3: continue
        sg = sn[mask]
        print(f"\n  Group: {g} (n={mask.sum()})")
        print(f"  σ_el: {sg.min():.1f} ~ {sg.max():.1f}")
        for vn in vars_check:
            arr = T[vn][mask]
            if arr.std() < 1e-10: continue
            # Rank correlation (Spearman)
            from scipy.stats import spearmanr
            rho, _ = spearmanr(arr, sg)
            if abs(rho) > 0.5:
                flag = '★' if abs(rho) > 0.8 else '●'
                print(f"    {flag}{vn:12s}: ρ={rho:+.3f}  range={arr.min():.4f}~{arr.max():.4f}")
                corr_pool[vn].append(rho)

    # Aggregate
    print(f"\n  ── Aggregate (모든 그룹에서 일관된 방향) ──")
    for vn in vars_check:
        if vn not in corr_pool: continue
        rhos = corr_pool[vn]
        if len(rhos) < 2: continue
        mean_rho = np.mean(rhos)
        consistent = all(r > 0 for r in rhos) or all(r < 0 for r in rhos)
        flag = '★★' if consistent and abs(mean_rho) > 0.6 else '★' if consistent else ' '
        print(f"    {flag}{vn:12s}: mean ρ={mean_rho:+.3f} ({len(rhos)} groups) {'CONSISTENT' if consistent else 'mixed'}")

    # ═══════════════════════════════════════════
    # 2. Consistent 변수로 formula
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("2. Intra-group consistent 변수로 thin formula")
    print("="*80)
    ok = T['ep'] >= 0.65; sn_ok = sn[ok]; n_ok = ok.sum()

    def fitC(s, r):
        v = (r>0)&np.isfinite(r)
        return float(np.exp(np.mean(np.log(s[v]/r[v])))) if v.sum()>=3 else None
    def loocv(s, r):
        n=len(s); la=np.log(s); lr=np.log(r); errs=[]
        for i in range(n):
            m=np.ones(n,bool); m[i]=False
            C=float(np.exp(np.mean(la[m]-lr[m])))
            errs.append((la[i]-np.log(C*r[i]))**2)
        return 1-np.sum(errs)/np.sum((la-np.mean(la))**2)

    # Based on what's consistent: try combinations
    c2 = []
    for a_cn in [0.5, 0.75, 1, 1.25, 1.5]:
        for vn1, va1, ves1 in [('por', T['por'][ok], [0.5, 1, 1.5, 2]),
                                ('ep', T['ep'][ok], [0.5, 1, 2]),
                                ('cn_x_ep', T['cn_x_ep'][ok], [0.5, 1]),
                                ('cn_x_por', T['cn_x_por'][ok], [0.5, 1])]:
            for v1 in ves1:
                for vn2, va2, ves2 in [('cov', T['cov'][ok], [0, 0.5, 1]),
                                        ('d2a', T['d2a'][ok], [0, 0.25]),
                                        ('hop', T['hop'][ok], [0, 0.25, 0.5]),
                                        ('tau', T['tau'][ok], [0, 0.5])]:
                    for v2 in ves2:
                        for td in [0.25, 0.5]:
                            rhs = SAM * T['cn'][ok]**a_cn * va1**v1 * va2**v2 / T['xi'][ok]**td
                            C = fitC(sn_ok, rhs)
                            if C is None: continue
                            r2 = r2l(sn_ok, C*rhs)
                            if r2 > 0.90:
                                s2 = f'×{vn2}^{v2}' if v2 else ''
                                c2.append({'r2':r2, 'rhs':rhs, 'C':C,
                                    'label':f'CN^{a_cn}×{vn1}^{v1}{s2}/(T/d)^{td}'})
    c2.sort(key=lambda x: -x['r2'])
    for i, r in enumerate(c2[:25]):
        cv = loocv(sn_ok, r['rhs']) if i < 10 else 0
        cs = f" LOOCV={cv:.4f}" if cv else ""
        pred = r['C'] * r['rhs']
        errs = (pred - sn_ok) / sn_ok * 100
        print(f"  #{i+1} R²={r['r2']:.4f}{cs} {errs.min():+.0f}~{errs.max():+.0f}%  {r['label']}")

if __name__ == '__main__':
    main()
