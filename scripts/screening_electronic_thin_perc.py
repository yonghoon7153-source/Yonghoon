"""
Thin: percolation 안 되는 케이스 제외 후 fitting
"""
import json, os, numpy as np, warnings
from pathlib import Path
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
            el_perc = m.get('electronic_percolating_fraction', 0)
            el_act = m.get('electronic_active_fraction', 0)
            ps = m.get('phi_se', 0)
            A = m.get('am_am_mean_area', 0)
            d = m.get('am_am_mean_delta', 0)
            hop = max(m.get('am_am_mean_hop', 0), 0.1)
            if pa <= 0 or cn <= 0 or T <= 0 or d <= 0 or A <= 0: continue
            rows.append({
                'sel': sel, 'pa': pa, 'cn': cn, 'tau': tau, 'cov': cov,
                'xi': T/d_am, 'por': max(por, 0.1), 'el_perc': el_perc,
                'el_act': el_act, 'ps': max(ps, 0.01),
                'd': d, 'A': A, 'hop': hop, 'd2a': d**2/A,
                'name': mp.parent.name
            })
    seen = set(); u = []
    for r in rows:
        k = f"{r['pa']:.4f}_{r['xi']:.1f}"
        if k not in seen: seen.add(k); u.append(r)
    return u

def r2l(a, p):
    la, lp = np.log(a), np.log(p); return 1-np.sum((la-lp)**2)/np.sum((la-np.mean(la))**2)
def fitC(a, r):
    v = (r>0)&np.isfinite(r); return float(np.exp(np.mean(np.log(a[v]/r[v])))) if v.sum()>=3 else None
def loocv_C(sn, rhs):
    n=len(sn); la=np.log(sn); lr=np.log(rhs); errs=[]
    for i in range(n):
        m=np.ones(n,bool); m[i]=False
        C_loo=float(np.exp(np.mean(la[m]-lr[m])))
        errs.append((la[i]-np.log(C_loo*rhs[i]))**2)
    return 1-np.sum(errs)/np.sum((la-np.mean(la))**2)

def main():
    rows = load_data()
    if not rows: print("데이터 없음!"); return
    V = {k: np.array([r[k] for r in rows]) for k in rows[0] if k != 'name'}
    names = [r['name'] for r in rows]
    n = len(rows); thin = V['xi'] < 10

    # ═══════════════════════════════════════════
    # 0. Thin el_perc 분포
    # ═══════════════════════════════════════════
    print(f"n={n}, thin={thin.sum()}")
    print(f"\n{'='*80}")
    print("0. Thin el_perc 분포")
    print("="*80)
    tn_idx = np.where(thin)[0]
    for i in np.argsort(V['el_perc'][thin]):
        j = tn_idx[i]
        flag = '✗' if V['el_perc'][j] < 0.90 else ' '
        print(f"  {flag}{names[j]:35s} el_perc={V['el_perc'][j]:.3f} el_act={V['el_act'][j]:.3f} "
              f"σ={V['sel'][j]:.2f} CN={V['cn'][j]:.2f} T/d={V['xi'][j]:.1f}")

    # ═══════════════════════════════════════════
    # 1. el_perc threshold scan
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("1. el_perc threshold scan")
    print("="*80)
    formulas = [
        ('CN×por^2.5×cov^1.5/√ξ [topo]', lambda T: SAM*T['cn']*T['por']**2.5*T['cov']**1.5/T['xi']**0.5),
        ('CN×por×φ_SE×√(d2a)/√ξ [hybrid]', lambda T: SAM*T['cn']*T['por']*T['ps']*np.sqrt(T['d2a'])/T['xi']**0.5),
        ('CN×√(d2a×hop/ξ) [current]', lambda T: SAM*T['cn']*np.sqrt(T['d2a']*T['hop']/T['xi'])),
        ('CN×(d2a)^¾/√ξ [clean]', lambda T: SAM*T['cn']*T['d2a']**0.75/T['xi']**0.5),
    ]

    for thresh in [0.0, 0.80, 0.85, 0.90, 0.95]:
        mask = thin & (V['el_perc'] >= thresh)
        nt = mask.sum()
        if nt < 5: continue
        sn = V['sel'][mask]
        T = {k: V[k][mask] for k in V}
        print(f"\n  ── el_perc ≥ {thresh} (n={nt}) ──")
        for label, fn in formulas:
            rhs = fn(T)
            C = fitC(sn, rhs)
            if C is None: continue
            r2 = r2l(sn, C*rhs); cv = loocv_C(sn, rhs)
            err = np.mean(np.abs(sn - C*rhs) / sn * 100)
            w20 = np.sum(np.abs(sn - C*rhs) / sn < 0.2)
            errs_signed = (C*rhs - sn) / sn * 100
            print(f"    {label:40s} R²={r2:.4f} LOOCV={cv:.4f} |err|={err:.0f}% ({w20}/{nt}) range={errs_signed.min():+.0f}~{errs_signed.max():+.0f}%")

    # ═══════════════════════════════════════════
    # 2. Best threshold에서 per-case
    # ═══════════════════════════════════════════
    for thresh in [0.90, 0.85]:
        mask = thin & (V['el_perc'] >= thresh)
        nt = mask.sum()
        if nt < 5: continue
        sn = V['sel'][mask]
        T = {k: V[k][mask] for k in V}
        tn_names = [names[i] for i in np.where(mask)[0]]

        print(f"\n{'='*80}")
        print(f"2. Per-case (el_perc ≥ {thresh}, n={nt})")
        print("="*80)
        for label, fn in formulas[:2]:
            rhs = fn(T); C = fitC(sn, rhs); pred = C * rhs
            r2 = r2l(sn, pred)
            errs = (pred - sn) / sn * 100
            print(f"\n  {label}  R²={r2:.4f}")
            for j in np.argsort(-np.abs(errs)):
                flag = '!' if abs(errs[j]) > 20 else ' '
                print(f"  {flag}{tn_names[j]:35s} σ={sn[j]:6.2f} pred={pred[j]:6.2f} err={errs[j]:+5.0f}% "
                      f"ep={T['el_perc'][j]:.3f} CN={T['cn'][j]:.2f} ξ={T['xi'][j]:.1f}")

    # ═══════════════════════════════════════════
    # 3. 전체 formula grid (best threshold)
    # ═══════════════════════════════════════════
    best_thresh = 0.90
    mask = thin & (V['el_perc'] >= best_thresh)
    nt = mask.sum()
    if nt >= 5:
        sn = V['sel'][mask]
        T = {k: V[k][mask] for k in V}
        print(f"\n{'='*80}")
        print(f"3. Formula grid (el_perc ≥ {best_thresh}, n={nt})")
        print("="*80)
        grid = []
        for a in [0.75, 1, 1.25]:
            for b in [1, 1.5, 2, 2.5]:
                for c in [0.5, 1, 1.5]:
                    for d in [0.25, 0.5]:
                        rhs = SAM * T['cn']**a * T['por']**b * T['cov']**c / T['xi']**d
                        C = fitC(sn, rhs)
                        if C is None: continue
                        r2 = r2l(sn, C*rhs)
                        if r2 > 0.90:
                            grid.append({'r2':r2, 'rhs':rhs,
                                'label':f'CN^{a}×por^{b}×cov^{c}/(T/d)^{d}'})
        # + d2a
        for a in [0.75, 1]:
            for b in [1, 1.5, 2]:
                for c in [0.5, 1]:
                    for d2 in [0.125, 0.25, 0.5]:
                        for e in [0.25, 0.5]:
                            rhs = SAM * T['cn']**a * T['por']**b * T['cov']**c * T['d2a']**d2 / T['xi']**e
                            C = fitC(sn, rhs)
                            if C is None: continue
                            r2 = r2l(sn, C*rhs)
                            if r2 > 0.92:
                                grid.append({'r2':r2, 'rhs':rhs,
                                    'label':f'CN^{a}×por^{b}×cov^{c}×(δ²/A)^{d2}/(T/d)^{e}'})
        grid.sort(key=lambda x: -x['r2'])
        for i, r in enumerate(grid[:20]):
            cv = loocv_C(sn, r['rhs']) if i < 10 else 0
            cs = f" LOOCV={cv:.4f}" if cv else ""
            print(f"  #{i+1} R²={r['r2']:.4f}{cs}  {r['label']}")

if __name__ == '__main__':
    main()
