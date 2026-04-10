"""
Thin v2: 개형(trend) 안 맞는 문제 해결
=====================================
문제: CN×√(δ²×hop/(A×ξ))가 진폭을 과장
원인: δ²/A 변동이 크면 √ 안에서 증폭

접근:
1. δ²/A를 log로 soft화: ln(1+δ²/A) 등
2. 지수를 더 낮추기 (¾→½→¼)
3. 다른 base: por, cov 기반 + 접촉 보정
4. φ 추가해서 안정화
5. 그룹별 C 문제인지 vs 식 문제인지 구분
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
            el_perc = max(m.get('electronic_percolating_fraction', 0), 0.01)
            ps = m.get('phi_se', 0)
            A = m.get('am_am_mean_area', 0)
            d = m.get('am_am_mean_delta', 0)
            hop = max(m.get('am_am_mean_hop', 0), 0.1)
            if pa <= 0 or cn <= 0 or T <= 0 or d <= 0 or A <= 0: continue
            rows.append({
                'sel': sel, 'pa': pa, 'cn': cn, 'tau': tau, 'cov': cov,
                'xi': T/d_am, 'por': max(por, 0.1), 'el_perc': el_perc,
                'ps': max(ps, 0.01), 'd': d, 'A': A, 'hop': hop,
                'd2a': d**2/A, 'name': mp.parent.name
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
    n = len(rows); thin = V['xi'] < 10; sn = V['sel'][thin]
    T = {k: V[k][thin] for k in V}
    n_thin = thin.sum()
    print(f"n={n}, thin={n_thin}")
    print(f"δ²/A range: {T['d2a'].min():.5f}~{T['d2a'].max():.5f} (ratio={T['d2a'].max()/T['d2a'].min():.1f}x)")
    print(f"hop range: {T['hop'].min():.1f}~{T['hop'].max():.1f} (ratio={T['hop'].max()/T['hop'].min():.1f}x)")
    print(f"CN range: {T['cn'].min():.2f}~{T['cn'].max():.2f}")

    # ═══════════════════════════════════════════
    # 0. 현재 식 문제 진단
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("0. 현재 thin 식 진단")
    print("="*80)
    rhs0 = SAM * T['cn'] * np.sqrt(T['d2a'] * T['hop'] / T['xi'])
    C0 = fitC(sn, rhs0)
    pred0 = C0 * rhs0
    err0 = (pred0 - sn) / sn * 100  # signed error
    print(f"  Current: CN×√(δ²×hop/(A×ξ)), C={C0:.4f}, R²={r2l(sn,pred0):.4f}")
    print(f"\n  Case-by-case (signed error — +는 과대, -는 과소):")
    tn_names = [names[i] for i in np.where(thin)[0]]
    for j in np.argsort(-np.abs(err0)):
        print(f"    {tn_names[j]:35s} σ={sn[j]:6.2f} pred={pred0[j]:6.2f} err={err0[j]:+5.0f}% "
              f"d2a={T['d2a'][j]:.5f} hop={T['hop'][j]:5.1f} CN={T['cn'][j]:.2f} ξ={T['xi'][j]:.1f}")

    # ═══════════════════════════════════════════
    # 1. 접촉 역학 없는 식 (topology only)
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("1. Topology only (접촉 역학 없이)")
    print("="*80)
    topo = [
        ('CN×por^2.5×cov^1.5/√ξ', T['cn']*T['por']**2.5*T['cov']**1.5/T['xi']**0.5),
        ('CN×por^2×cov/√ξ', T['cn']*T['por']**2*T['cov']/T['xi']**0.5),
        ('CN×por^2×φ_SE^1.5/√ξ', T['cn']*T['por']**2*T['ps']**1.5/T['xi']**0.5),
        ('CN×por^3×cov/√ξ', T['cn']*T['por']**3*T['cov']/T['xi']**0.5),
        ('CN×por^2×cov^1.5×√τ/√ξ', T['cn']*T['por']**2*T['cov']**1.5*T['tau']**0.5/T['xi']**0.5),
        ('CN^1.25×por^2×cov^1.5/√ξ', T['cn']**1.25*T['por']**2*T['cov']**1.5/T['xi']**0.5),
        ('CN×por^2.5×cov/√ξ', T['cn']*T['por']**2.5*T['cov']/T['xi']**0.5),
        ('√[CN²×por⁶×φ_SE²×cov²/ξ]', np.sqrt(T['cn']**2*T['por']**6*T['ps']**2*T['cov']**2/T['xi'])),
    ]
    print(f"  {'Formula':45s} {'R²':>7s} {'LOOCV':>7s} {'|err|':>6s}")
    print("  "+"-"*65)
    for label, rhs_raw in sorted(topo, key=lambda x: -r2l(sn,fitC(sn,SAM*x[1])*(SAM*x[1])) if fitC(sn,SAM*x[1]) else -999):
        rhs = SAM * rhs_raw; C = fitC(sn, rhs)
        if C is None: continue
        r2 = r2l(sn, C*rhs); cv = loocv_C(sn, rhs)
        err = np.mean(np.abs(sn - C*rhs) / sn * 100)
        w20 = np.sum(np.abs(sn - C*rhs) / sn < 0.2)
        flag = '●' if r2 > 0.89 else ' '
        print(f"  {flag}{label:44s} {r2:7.4f} {cv:7.4f} {err:5.0f}% ({w20}/{n_thin})")
        # Show trend comparison
        if r2 > 0.89:
            pred = C * rhs
            errs = (pred - sn) / sn * 100
            print(f"    signed err range: {errs.min():+.0f}% ~ {errs.max():+.0f}%")

    # ═══════════════════════════════════════════
    # 2. 접촉 역학 soft화 (지수 낮추기)
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("2. 접촉 역학 soft (지수 ¼로)")
    print("="*80)
    soft = [
        ('CN×(d2a)^¼/√ξ', T['cn']*T['d2a']**0.25/T['xi']**0.5),
        ('CN×(d2a)^¼×hop^¼/√ξ', T['cn']*T['d2a']**0.25*T['hop']**0.25/T['xi']**0.5),
        ('CN×(d2a)^¼×por/√ξ', T['cn']*T['d2a']**0.25*T['por']/T['xi']**0.5),
        ('CN×(d2a)^¼×cov/√ξ', T['cn']*T['d2a']**0.25*T['cov']/T['xi']**0.5),
        ('CN×(d2a)^¼×por×cov/√ξ', T['cn']*T['d2a']**0.25*T['por']*T['cov']/T['xi']**0.5),
        ('CN×(d2a)^½×por/√ξ', T['cn']*T['d2a']**0.5*T['por']/T['xi']**0.5),
        ('CN×(d2a)^¼×por^1.5×cov/√ξ', T['cn']*T['d2a']**0.25*T['por']**1.5*T['cov']/T['xi']**0.5),
        ('CN×√(d2a)×por^1.5/√ξ', T['cn']*T['d2a']**0.5*T['por']**1.5/T['xi']**0.5),
    ]
    print(f"  {'Formula':45s} {'R²':>7s} {'LOOCV':>7s} {'|err|':>6s}")
    print("  "+"-"*65)
    for label, rhs_raw in sorted(soft, key=lambda x: -r2l(sn,fitC(sn,SAM*x[1])*(SAM*x[1])) if fitC(sn,SAM*x[1]) else -999):
        rhs = SAM * rhs_raw; C = fitC(sn, rhs)
        if C is None: continue
        r2 = r2l(sn, C*rhs); cv = loocv_C(sn, rhs)
        err = np.mean(np.abs(sn - C*rhs) / sn * 100)
        w20 = np.sum(np.abs(sn - C*rhs) / sn < 0.2)
        flag = '★' if r2 > 0.92 else '●' if r2 > 0.90 else ' '
        print(f"  {flag}{label:44s} {r2:7.4f} {cv:7.4f} {err:5.0f}% ({w20}/{n_thin})")

    # ═══════════════════════════════════════════
    # 3. Hybrid: topology + soft contact
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("3. Hybrid: topology base + contact 보정")
    print("="*80)
    hybrid = []
    # base = CN × por^a × cov^b / √ξ, + (d2a)^c
    for a in [1.5, 2, 2.5]:
        for b in [0.5, 1, 1.5]:
            for c in [0, 0.125, 0.25, 0.375, 0.5]:
                rhs = SAM * T['cn'] * T['por']**a * T['cov']**b * T['d2a']**c / T['xi']**0.5
                C = fitC(sn, rhs)
                if C is None: continue
                r2 = r2l(sn, C*rhs)
                if r2 > 0.90:
                    cs = f'×(δ²/A)^{c}' if c else ''
                    hybrid.append({'r2':r2, 'rhs':rhs,
                        'label':f'CN×por^{a}×cov^{b}{cs}/√ξ'})
    # + φ_SE
    for a in [1, 1.5, 2]:
        for b in [0.5, 1, 1.5]:
            for c in [0, 0.25, 0.5]:
                rhs = SAM * T['cn'] * T['por']**a * T['ps']**b * T['d2a']**c / T['xi']**0.5
                C = fitC(sn, rhs)
                if C is None: continue
                r2 = r2l(sn, C*rhs)
                if r2 > 0.90:
                    cs = f'×(δ²/A)^{c}' if c else ''
                    hybrid.append({'r2':r2, 'rhs':rhs,
                        'label':f'CN×por^{a}×φ_SE^{b}{cs}/√ξ'})
    hybrid.sort(key=lambda x: -x['r2'])
    for i, r in enumerate(hybrid[:20]):
        cv = loocv_C(sn, r['rhs']) if i < 10 else 0
        cs = f" LOOCV={cv:.4f}" if cv else ""
        err = np.mean(np.abs(sn - fitC(sn,r['rhs'])*r['rhs']) / sn * 100)
        w20 = np.sum(np.abs(sn - fitC(sn,r['rhs'])*r['rhs']) / sn < 0.2)
        print(f"  #{i+1} R²={r['r2']:.4f}{cs} |err|={err:.0f}% ({w20}/{n_thin})  {r['label']}")

    # ═══════════════════════════════════════════
    # 4. 그룹별 개형 비교 (best 3)
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("4. 개형 비교: best 식들의 per-case signed error")
    print("="*80)
    candidates = [
        ('CN×√(d2a×hop/ξ) [현재]', SAM*T['cn']*np.sqrt(T['d2a']*T['hop']/T['xi'])),
        ('CN×por^2.5×cov^1.5/√ξ [topo]', SAM*T['cn']*T['por']**2.5*T['cov']**1.5/T['xi']**0.5),
    ]
    if hybrid:
        candidates.append((f'Hybrid: {hybrid[0]["label"]}', hybrid[0]['rhs']))
    for label, rhs in candidates:
        C = fitC(sn, rhs); pred = C * rhs
        r2 = r2l(sn, pred)
        errs = (pred - sn) / sn * 100
        print(f"\n  {label}  R²={r2:.4f}")
        print(f"  {'Case':35s} {'σ':>6s} {'pred':>6s} {'err':>6s} {'ξ':>4s}")
        for j in np.argsort(-np.abs(errs)):
            flag = '!' if abs(errs[j]) > 20 else ' '
            print(f"  {flag}{tn_names[j]:34s} {sn[j]:6.2f} {pred[j]:6.2f} {errs[j]:+5.0f}% {T['xi'][j]:4.1f}")

if __name__ == '__main__':
    main()
