"""
Thin: cn_std 기반 (intra-group ρ=0.958!)
cov 제거, cn_std 추가
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
            cov = max(m.get('coverage_AM_P_mean', m.get('coverage_AM_S_mean', m.get('coverage_AM_mean', 20))), 0.1) / 100
            r_am = max(m.get('r_AM_P', 0), m.get('r_AM_S', 0))
            d_am = r_am * 2 if r_am > 0.1 else 5.0
            por = m.get('porosity', 0)
            ep = m.get('electronic_percolating_fraction', 0)
            cn_std = m.get('am_am_cn_std', 0)
            d = m.get('am_am_mean_delta', 0)
            A = m.get('am_am_mean_area', 0)
            if pa <= 0 or cn <= 0 or T <= 0: continue
            rows.append({
                'sel': sel, 'cn': cn, 'xi': T/d_am, 'por': max(por, 0.1),
                'ep': max(ep, 0.01), 'cov': cov,
                'cn_std': max(float(str(cn_std)), 0.01) if cn_std else 0.01,
                'd': max(d, 0.001), 'A': max(A, 0.01),
                'd2a': d**2/A if A > 0 and d > 0 else 0.001,
                'name': mp.parent.name
            })
    seen = set(); u = []
    for r in rows:
        k = f"{r['cn']:.3f}_{r['xi']:.1f}"
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
    n_thin = thin.sum()
    print(f"thin={n_thin}, perc(≥0.65)={n_ok}")
    print(f"cn_std range: {T['cn_std'].min():.2f}~{T['cn_std'].max():.2f}")

    # ═══════════════════════════════════════════
    # 1. CN × por × cn_std 조합 (cov 없이!)
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("1. CN × por × cn_std (cov 없이!)")
    print("="*80)
    c1 = []
    for a in [0.75, 1, 1.25, 1.5]:
        for b in [0.5, 1, 1.5, 2]:
            for c in [0.25, 0.5, 0.75, 1]:
                for d in [0.25, 0.5]:
                    rhs = SAM * T['cn'][ok]**a * T['por'][ok]**b * T['cn_std'][ok]**c / T['xi'][ok]**d
                    C = fitC(sn_ok, rhs)
                    if C is None: continue
                    r2 = r2l(sn_ok, C*rhs)
                    if r2 > 0.88:
                        errs = (C*rhs - sn_ok) / sn_ok * 100
                        c1.append({'r2':r2, 'label':f'CN^{a}×por^{b}×cn_std^{c}/(T/d)^{d}',
                                   'range':f'{errs.min():+.0f}~{errs.max():+.0f}%'})
    c1.sort(key=lambda x: -x['r2'])
    for i, r in enumerate(c1[:15]):
        print(f"  #{i+1} R²={r['r2']:.4f} {r['range']:>12s}  {r['label']}")

    # ═══════════════════════════════════════════
    # 2. + ep (soft gate)
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("2. √ep × CN × por × cn_std (ep gate)")
    print("="*80)
    c2 = []
    for ep_e in [0.5, 1]:
        for a in [0.75, 1, 1.25]:
            for b in [0.5, 1, 1.5]:
                for c in [0.25, 0.5, 0.75, 1]:
                    for d in [0.25, 0.5]:
                        rhs = SAM * T['ep'][ok]**ep_e * T['cn'][ok]**a * T['por'][ok]**b * T['cn_std'][ok]**c / T['xi'][ok]**d
                        C = fitC(sn_ok, rhs)
                        if C is None: continue
                        r2 = r2l(sn_ok, C*rhs)
                        if r2 > 0.90:
                            errs = (C*rhs - sn_ok) / sn_ok * 100
                            eps = f'√ep×' if ep_e == 0.5 else 'ep×'
                            c2.append({'r2':r2, 'range':f'{errs.min():+.0f}~{errs.max():+.0f}%',
                                'label':f'{eps}CN^{a}×por^{b}×cn_std^{c}/(T/d)^{d}'})
    c2.sort(key=lambda x: -x['r2'])
    for i, r in enumerate(c2[:15]):
        print(f"  #{i+1} R²={r['r2']:.4f} {r['range']:>12s}  {r['label']}")

    # ═══════════════════════════════════════════
    # 3. + δ (consistent ρ=0.769)
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("3. + δ or d2a (soft)")
    print("="*80)
    c3 = []
    for a in [1, 1.25]:
        for b in [0.5, 1]:
            for c in [0.25, 0.5]:
                for d_e in [0, 0.25, 0.5]:
                    for td in [0.25, 0.5]:
                        rhs = SAM * T['ep'][ok]**0.5 * T['cn'][ok]**a * T['por'][ok]**b * T['cn_std'][ok]**c * T['d2a'][ok]**d_e / T['xi'][ok]**td
                        C = fitC(sn_ok, rhs)
                        if C is None: continue
                        r2 = r2l(sn_ok, C*rhs)
                        if r2 > 0.90:
                            errs = (C*rhs - sn_ok) / sn_ok * 100
                            ds = f'×(δ²/A)^{d_e}' if d_e else ''
                            c3.append({'r2':r2, 'range':f'{errs.min():+.0f}~{errs.max():+.0f}%',
                                'label':f'√ep×CN^{a}×por^{b}×cn_std^{c}{ds}/(T/d)^{td}'})
    c3.sort(key=lambda x: -x['r2'])
    for i, r in enumerate(c3[:15]):
        print(f"  #{i+1} R²={r['r2']:.4f} {r['range']:>12s}  {r['label']}")

    # ═══════════════════════════════════════════
    # 4. 비교: cov vs cn_std
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("4. 비교: cov vs cn_std")
    print("="*80)
    comparisons = [
        ('√ep×CN×por×√cov/√ξ [현재]', SAM*T['ep'][ok]**0.5*T['cn'][ok]*T['por'][ok]*T['cov'][ok]**0.5/T['xi'][ok]**0.5),
        ('√ep×CN×por×√cn_std/√ξ [NEW]', SAM*T['ep'][ok]**0.5*T['cn'][ok]*T['por'][ok]*T['cn_std'][ok]**0.5/T['xi'][ok]**0.5),
        ('√ep×CN×por×cn_std/√ξ', SAM*T['ep'][ok]**0.5*T['cn'][ok]*T['por'][ok]*T['cn_std'][ok]/T['xi'][ok]**0.5),
        ('√ep×CN×por×cn_std^¼/√ξ', SAM*T['ep'][ok]**0.5*T['cn'][ok]*T['por'][ok]*T['cn_std'][ok]**0.25/T['xi'][ok]**0.5),
        ('√ep×CN×por/√ξ [no 4th]', SAM*T['ep'][ok]**0.5*T['cn'][ok]*T['por'][ok]/T['xi'][ok]**0.5),
    ]
    tn_names = [names[i] for i in np.where(thin)[0]]
    ok_names = [tn_names[i] for i in np.where(ok)[0]]
    for label, rhs in comparisons:
        C = fitC(sn_ok, rhs)
        if C is None: continue
        pred = C * rhs
        r2 = r2l(sn_ok, pred)
        errs = (pred - sn_ok) / sn_ok * 100
        print(f"\n  {label}  R²={r2:.4f} range={errs.min():+.0f}~{errs.max():+.0f}%")
        for j in np.argsort(-np.abs(errs))[:5]:
            print(f"    {ok_names[j]:35s} σ={sn_ok[j]:6.2f} pred={pred[j]:6.2f} err={errs[j]:+5.0f}%")

if __name__ == '__main__':
    main()
