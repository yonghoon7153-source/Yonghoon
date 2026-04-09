"""
Electronic: exp 계수 + φ_AM 지수 정밀 스윕
============================================
Free fit: φ_AM^3.45 × CN^1.72 × exp(0.70×π/(T/d)) × φ_SE^0.81 × cov^0.18
→ exp 계수가 1이 아니라 0.7!
→ φ_AM이 1.5가 아니라 3.5!
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
            am_cn=max(m.get('am_am_cn',0.01),0.01);cn=m.get('se_se_cn',0)
            T=m.get('thickness_um',0);tau=m.get('tortuosity_recommended',m.get('tortuosity_mean',0))
            cov=max(m.get('coverage_AM_P_mean',m.get('coverage_AM_S_mean',m.get('coverage_AM_mean',20))),0.1)/100
            r_am=max(m.get('r_AM_P',0),m.get('r_AM_S',0))
            d_am=r_am*2 if r_am>0.1 else 5.0
            if pa<=0 or am_cn<=0 or T<=0: continue
            rows.append({'sel':sel,'pa':pa,'ps':ps,'am_cn':am_cn,'cn':cn,
                'T':T,'d_am':d_am,'tau':tau,'cov':cov,'ratio':T/d_am,'name':mp.parent.name})
    seen=set();u=[]
    for r in rows:
        k=f"{r['pa']:.4f}_{r['T']:.1f}"
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
    rows=load_data();n=len(rows);print(f"n={n}\n")
    sel=np.array([r['sel'] for r in rows])
    pa=np.array([r['pa'] for r in rows]);ps=np.array([r['ps'] for r in rows])
    am_cn=np.array([r['am_cn'] for r in rows]);cn=np.array([r['cn'] for r in rows])
    tau=np.array([r['tau'] for r in rows]);cov=np.array([r['cov'] for r in rows])
    ratio=np.array([r['ratio'] for r in rows])

    results=[]

    # ═══════════════════════════════════════
    # SWEEP: φ_AM^a × CN^b × exp(k×π/(T/d)) × φ_SE^c × cov^d
    # a: 1~4, b: 1~2.5, k: 0.5~1.2, c: 0~1, d: 0~0.5
    # ═══════════════════════════════════════
    print("FULL SWEEP: φ_AM^a × CN^b × exp(k×π/(T/d)) × φ_SE^c × cov^d")
    print("="*90)

    for a in np.arange(1.0, 4.1, 0.5):
        for b in np.arange(1.0, 2.6, 0.5):
            for k in np.arange(0.4, 1.3, 0.1):
                for c in [0, 0.25, 0.5, 0.75, 1.0]:
                    for d in [0, 0.25, 0.5]:
                        rhs = SAM * pa**a * am_cn**b * np.exp(k*np.pi/ratio) * ps**c * cov**d
                        C = fitC(sel, rhs)
                        if C is None or C <= 0: continue
                        pred = C * rhs; r2 = r2l(sel, pred)
                        if r2 > 0.90:
                            results.append({
                                'a':a,'b':b,'k':k,'c':c,'d':d,
                                'r2':r2,'C':C,'rhs':rhs,
                                'name':f'φ^{a}×CN^{b}×exp({k:.1f}π/r)×φ_SE^{c}×cov^{d}'
                            })

    results.sort(key=lambda x: -x['r2'])

    print(f"\nTop 20:")
    seen=set(); cnt=0
    for r in results:
        if r['name'] in seen: continue
        seen.add(r['name']); cnt+=1
        if cnt > 20: break
        cv = loocv_C(sel, r['rhs']) if cnt <= 5 else 0
        cv_str = f" LOOCV={cv:.4f}" if cv else ""
        print(f"  #{cnt:2d} R²={r['r2']:.4f}{cv_str}  {r['name']}  C={r['C']:.4f}")

    # ═══════════════════════════════════════
    # Also: WITHOUT exp, high φ_AM exponent
    # ═══════════════════════════════════════
    print(f"\n{'='*90}")
    print("WITHOUT exp (pure power law, high φ_AM)")
    print("="*90)

    results2 = []
    for a in np.arange(2.0, 5.1, 0.5):
        for b in np.arange(0.5, 2.6, 0.5):
            for c in [0, 0.5, 1.0]:
                for d in [0, 0.25, 0.5]:
                    for te in [0, 0.5, 1.0]:  # τ exponent
                        rhs = SAM * pa**a * am_cn**b * ps**c * cov**d * tau**te
                        C = fitC(sel, rhs)
                        if C is None or C <= 0: continue
                        pred = C*rhs; r2 = r2l(sel, pred)
                        if r2 > 0.88:
                            results2.append({
                                'r2':r2,'C':C,'rhs':rhs,
                                'name':f'NO_EXP: φ^{a}×CN^{b}×φ_SE^{c}×cov^{d}×τ^{te}'
                            })

    results2.sort(key=lambda x: -x['r2'])
    for i, r in enumerate(results2[:10]):
        print(f"  #{i+1} R²={r['r2']:.4f}  {r['name']}")

    # Best overall comparison
    print(f"\n{'='*90}")
    print("COMPARISON")
    print("="*90)
    cur = SAM * pa**1.5 * am_cn**2 * np.exp(np.pi/ratio)
    C_cur = fitC(sel, cur); r2_cur = r2l(sel, C_cur*cur)
    print(f"  Current:  φ^1.5×CN²×exp(π/r)            R²={r2_cur:.4f}")
    if results:
        b = results[0]
        print(f"  Best:     {b['name']:40s} R²={b['r2']:.4f}")
    if results2:
        b2 = results2[0]
        print(f"  Best(no exp): {b2['name']:36s} R²={b2['r2']:.4f}")

if __name__ == '__main__':
    main()
