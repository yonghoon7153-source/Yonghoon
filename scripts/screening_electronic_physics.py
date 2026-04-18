"""
Electronic thin: 물리적 의미로 변수 그룹핑
==========================================
Group A: AM network factor (CN × por)
  → CN: AM-AM 접촉 수, por: AM이 움직일 수 있는 공간
Group B: SE interaction factor (φ_SE × cov)
  → φ_SE: SE가 AM을 구조화, cov: AM-SE 계면 품질
Group C: Geometry (T/d)
  → finite-size percolation
"""
import json,os,numpy as np,warnings
from pathlib import Path
warnings.filterwarnings('ignore')
WEBAPP=os.path.join(os.path.dirname(os.path.dirname(__file__)),'webapp')
SAM=50.0

def load_data():
    rows=[]
    for base in [Path(WEBAPP)/'results',Path(WEBAPP)/'archive']:
        if not base.is_dir(): continue
        for mp in base.rglob('full_metrics.json'):
            try:
                with open(mp) as f: m=json.load(f)
            except: continue
            sel=m.get('electronic_sigma_full_mScm',0)
            if not sel or sel<0.001: continue
            pa=max(m.get('phi_am',0),0.01);ps=m.get('phi_se',0)
            am_cn=max(m.get('am_am_cn',0.01),0.01)
            T=m.get('thickness_um',0)
            cov=max(m.get('coverage_AM_P_mean',m.get('coverage_AM_S_mean',m.get('coverage_AM_mean',20))),0.1)/100
            r_am=max(m.get('r_AM_P',0),m.get('r_AM_S',0))
            d_am=r_am*2 if r_am>0.1 else 5.0
            por=m.get('porosity',0)
            if pa<=0 or am_cn<=0 or T<=0: continue
            rows.append({'sel':sel,'pa':pa,'ps':ps,'am_cn':am_cn,
                'ratio':T/d_am,'cov':cov,'por':max(por,0.1)})
    seen=set();u=[]
    for r in rows:
        k=f"{r['pa']:.4f}_{r['ratio']:.1f}"
        if k not in seen: seen.add(k);u.append(r)
    return u

def r2l(a,p):
    la,lp=np.log(a),np.log(p);return 1-np.sum((la-lp)**2)/np.sum((la-np.mean(la))**2)
def fitC(a,r):
    v=(r>0)&np.isfinite(r);return float(np.exp(np.mean(np.log(a[v]/r[v])))) if v.sum()>=3 else None
def loocv_C(sn,rhs):
    n=len(sn);la=np.log(sn);lr=np.log(rhs);errs=[]
    for i in range(n):
        m=np.ones(n,bool);m[i]=False
        C_loo=float(np.exp(np.mean(la[m]-lr[m])))
        errs.append((la[i]-np.log(C_loo*rhs[i]))**2)
    return 1-np.sum(errs)/np.sum((la-np.mean(la))**2)

def main():
    rows=load_data()
    sel=np.array([r['sel'] for r in rows])
    pa=np.array([r['pa'] for r in rows]);ps=np.array([r['ps'] for r in rows])
    am_cn=np.array([r['am_cn'] for r in rows])
    ratio=np.array([r['ratio'] for r in rows])
    cov=np.array([r['cov'] for r in rows]);por=np.array([r['por'] for r in rows])
    thin=ratio<10
    sn=sel[thin];cn_n=am_cn[thin];rn=ratio[thin];porn=por[thin];cvn=cov[thin];psn=ps[thin]
    n_thin=thin.sum()
    print(f"Thin n={n_thin}\n")

    # ═══════════════════════════════════════
    # PHYSICS-BASED GROUPING
    # ═══════════════════════════════════════
    print("="*80)
    print("Physics groups:")
    print("  A = AM network: CN^a × por^b")
    print("  B = SE interaction: φ_SE^c × cov^d")
    print("  C = geometry: (T/d)^e")
    print("  σ = C × σ_AM × A × B × C")
    print("="*80)

    results=[]
    for a in [0.5, 0.75, 1, 1.25, 1.5, 2]:
        for b in [1, 1.5, 2, 2.5, 3]:
            for c in [0.5, 0.75, 1, 1.25, 1.5]:
                for d in [0.5, 0.75, 1, 1.25]:
                    for e in [-1, -0.75, -0.5, -0.25]:
                        A = cn_n**a * porn**b
                        B = psn**c * cvn**d
                        G = rn**e
                        rhs = SAM * A * B * G
                        C_fit = fitC(sn, rhs)
                        if C_fit is None: continue
                        pred = C_fit * rhs; r2 = r2l(sn, pred)
                        if r2 > 0.88:
                            results.append({
                                'a':a,'b':b,'c':c,'d':d,'e':e,
                                'r2':r2,'C':C_fit,'rhs':rhs,
                                'label': f'(CN^{a}×por^{b}) × (φ_SE^{c}×cov^{d}) / (T/d)^{-e}'
                            })

    results.sort(key=lambda x:-x['r2'])

    print(f"\nTop 20:")
    for i,r in enumerate(results[:20]):
        cv = loocv_C(sn,r['rhs']) if i<3 else 0
        cv_str = f" LOOCV={cv:.4f}" if cv else ""
        print(f"  #{i+1} R²={r['r2']:.4f}{cv_str}")
        print(f"     {r['label']}")

    # ═══════════════════════════════════════
    # COMBINED VARIABLES
    # (CN × por)^a × (φ_SE × cov)^b / (T/d)^c
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("COMBINED: (CN×por)^a × (φ_SE×cov)^b / (T/d)^c")
    print("="*80)

    cn_por = cn_n * porn
    se_cov = psn * cvn

    results2=[]
    for a in [1, 1.25, 1.5, 1.75, 2, 2.5, 3]:
        for b in [0.5, 0.75, 1, 1.25, 1.5, 2]:
            for c in [0.25, 0.5, 0.75, 1]:
                rhs = SAM * cn_por**a * se_cov**b / rn**c
                C_fit = fitC(sn, rhs)
                if C_fit is None: continue
                pred = C_fit*rhs; r2 = r2l(sn, pred)
                if r2 > 0.85:
                    results2.append({
                        'a':a,'b':b,'c':c,'r2':r2,'C':C_fit,'rhs':rhs,
                        'label': f'(CN×por)^{a} × (φ_SE×cov)^{b} / (T/d)^{c}'
                    })

    results2.sort(key=lambda x:-x['r2'])
    print(f"\nTop 10:")
    for i,r in enumerate(results2[:10]):
        cv = loocv_C(sn,r['rhs']) if i<3 else 0
        cv_str = f" LOOCV={cv:.4f}" if cv else ""
        print(f"  #{i+1} R²={r['r2']:.4f}{cv_str}  {r['label']}")

    # ═══════════════════════════════════════
    # √ and ⁴√ elegant forms
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("ELEGANT FORMS")
    print("="*80)

    elegant = [
        ('√[(CN×por)⁴ × (φ_SE×cov)² / (T/d)]',
         np.sqrt(cn_por**4 * se_cov**2 / rn)),
        ('√[(CN×por)³ × (φ_SE×cov)² / (T/d)]',
         np.sqrt(cn_por**3 * se_cov**2 / rn)),
        ('√[(CN×por)⁴ × (φ_SE×cov) / (T/d)]',
         np.sqrt(cn_por**4 * se_cov / rn)),
        ('⁴√[(CN×por)⁸ × (φ_SE×cov)⁴ / (T/d)²]',
         (cn_por**8 * se_cov**4 / rn**2)**0.25),
        ('⁴√[(CN×por)⁶ × (φ_SE×cov)⁵ / (T/d)²]',
         (cn_por**6 * se_cov**5 / rn**2)**0.25),
        ('(CN×por)² × (φ_SE×cov) / √(T/d)',
         cn_por**2 * se_cov / rn**0.5),
        ('(CN×por)^(3/2) × (φ_SE×cov) / √(T/d)',
         cn_por**1.5 * se_cov / rn**0.5),
        ('(CN×por)² × √(φ_SE×cov) / √(T/d)',
         cn_por**2 * se_cov**0.5 / rn**0.5),
        ('(CN×por)^(3/2) × √(φ_SE×cov/(T/d))',
         cn_por**1.5 * np.sqrt(se_cov/rn)),
        ('(CN×por)² × √(φ_SE×cov/(T/d))',
         cn_por**2 * np.sqrt(se_cov/rn)),
    ]

    print(f"\n  {'Formula':55s} {'R²':>7s} {'LOOCV':>7s}")
    print("  "+"-"*70)
    for label, rhs in sorted(elegant, key=lambda x: -r2l(sn,fitC(sn,x[1])*x[1]) if fitC(sn,x[1]) else -999):
        C_fit = fitC(sn, rhs)
        if C_fit is None: continue
        r2 = r2l(sn, C_fit*rhs)
        cv = loocv_C(sn, rhs) if r2 > 0.89 else 0
        cv_str = f"{cv:7.4f}" if cv else "      -"
        flag = '★' if r2 > 0.91 else '●' if r2 > 0.89 else ' '
        print(f"  {flag}{label:54s} {r2:7.4f} {cv_str}")

    # ═══════════════════════════════════════
    # BEST comparison
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("PHYSICAL MEANING")
    print(f"{'='*80}")
    print("""
  σ_el(thin) = C × σ_AM × [AM network]^a × [SE interaction]^b / [geometry]^c

  AM network = CN × porosity
    CN: AM 입자가 서로 접촉하는 수 (connectivity)
    porosity: AM이 연결되기 위한 자유 공간 (void)
    → 접촉이 많고 공간이 있어야 AM network 형성

  SE interaction = φ_SE × coverage
    φ_SE: SE가 AM 사이에 존재 → AM 구조 안정화
    coverage: AM-SE 계면 품질 → stress transfer → contact quality
    → SE가 적절히 있어야 AM이 잘 배열됨

  Geometry = T/d (AM 입자 층수)
    T/d < 10에서 finite-size effect
    → 층수 부족 → spanning cluster 어려움
    """)

if __name__=='__main__':
    main()
