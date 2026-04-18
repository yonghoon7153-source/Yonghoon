"""
Thin 식 예쁘게 만들기
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
            T=m.get('thickness_um',0);tau=m.get('tortuosity_recommended',m.get('tortuosity_mean',0))
            cov=max(m.get('coverage_AM_P_mean',m.get('coverage_AM_S_mean',m.get('coverage_AM_mean',20))),0.1)/100
            r_am=max(m.get('r_AM_P',0),m.get('r_AM_S',0))
            d_am=r_am*2 if r_am>0.1 else 5.0
            por=m.get('porosity',0)
            if pa<=0 or am_cn<=0 or T<=0 or tau<=0: continue
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

def main():
    rows=load_data()
    sel=np.array([r['sel'] for r in rows])
    pa=np.array([r['pa'] for r in rows]);ps=np.array([r['ps'] for r in rows])
    am_cn=np.array([r['am_cn'] for r in rows])
    ratio=np.array([r['ratio'] for r in rows])
    cov=np.array([r['cov'] for r in rows]);por=np.array([r['por'] for r in rows])
    thin=ratio<10
    sn=sel[thin];cn_n=am_cn[thin];rn=ratio[thin];porn=por[thin];cvn=cov[thin];psn=ps[thin]

    print(f"Thin n={thin.sum()}\n")
    print("깔끔한 식 후보들:")
    print(f"{'Formula':60s} {'R²':>7s}")
    print("-"*70)

    candidates = [
        # Current best (ugly)
        ('CN^(5/4) × φ_SE^(5/4) × cov^(5/4) × por³ / √(T/d)',
         cn_n**1.25 * psn**1.25 * cvn**1.25 * porn**3 / rn**0.5),

        # All ^1 (cleanest)
        ('CN × φ_SE × cov × por³ / √(T/d)',
         cn_n * psn * cvn * porn**3 / rn**0.5),
        ('CN × φ_SE × cov × por² / √(T/d)',
         cn_n * psn * cvn * porn**2 / rn**0.5),
        ('CN × φ_SE × √cov × por³ / √(T/d)',
         cn_n * psn * cvn**0.5 * porn**3 / rn**0.5),

        # √ form: √[CN² × φ_SE² × cov² × por⁶ / (T/d)]
        ('√[CN²×φ_SE²×cov²×por⁶/(T/d)]',
         np.sqrt(cn_n**2 * psn**2 * cvn**2 * porn**6 / rn)),

        # Mixed
        ('CN × φ_SE × cov × por^(5/2) / √(T/d)',
         cn_n * psn * cvn * porn**2.5 / rn**0.5),
        ('CN^(3/2) × φ_SE × cov × por^(5/2) / √(T/d)',
         cn_n**1.5 * psn * cvn * porn**2.5 / rn**0.5),
        ('CN × φ_SE × cov × por³ / (T/d)',
         cn_n * psn * cvn * porn**3 / rn),

        # (CN×φ_SE×cov) as combined
        ('(CN×φ_SE×cov) × por³ / √(T/d)',
         cn_n * psn * cvn * porn**3 / rn**0.5),

        # ⁴√ form
        ('⁴√[CN⁵×φ_SE⁵×cov⁵×por¹²/(T/d)²]',
         (cn_n**5 * psn**5 * cvn**5 * porn**12 / rn**2)**0.25),

        # ³√ form
        ('³√[CN⁴×φ_SE⁴×cov⁴×por⁹/(T/d)^(3/2)]',
         (cn_n**4 * psn**4 * cvn**4 * porn**9 / rn**1.5)**(1/3)),

        # φ_SE and cov combined differently
        ('CN × (φ_SE×cov)^(3/2) × por³ / √(T/d)',
         cn_n * (psn*cvn)**1.5 * porn**3 / rn**0.5),
        ('CN × √(φ_SE×cov) × por³ / √(T/d)',
         cn_n * np.sqrt(psn*cvn) * porn**3 / rn**0.5),

        # por² instead of por³
        ('CN^(3/2) × φ_SE × cov × por² / √(T/d)',
         cn_n**1.5 * psn * cvn * porn**2 / rn**0.5),

        # Elegant: √[CN × φ_SE × cov × por³ / (T/d)]
        ('√[CN²×φ_SE×cov×por⁶/(T/d)]',
         np.sqrt(cn_n**2 * psn * cvn * porn**6 / rn)),

        # Very clean: CN × por³ × √(φ_SE × cov / (T/d))
        ('CN × por³ × √(φ_SE×cov/(T/d))',
         cn_n * porn**3 * np.sqrt(psn * cvn / rn)),
    ]

    for label, rhs in sorted(candidates, key=lambda x: -r2l(sn, fitC(sn,x[1])*x[1]) if fitC(sn,x[1]) else -999):
        C = fitC(sn, rhs)
        if C is None: continue
        r2 = r2l(sn, C*rhs)
        flag = '★' if r2 > 0.91 else '●' if r2 > 0.90 else ' '
        print(f"  {flag}{label:58s} {r2:7.4f}")

if __name__ == '__main__':
    main()
