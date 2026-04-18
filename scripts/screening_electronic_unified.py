"""
Electronic UNIFIED: thick top formula + exp for thin
====================================================
Thick best: φ_AM^4.4 × CN^1.4 × τ^0.78 × cov^1.14  R²=0.97
Thin에는 exp(π/(T/d)) 추가

통합: φ_AM^a × CN^b × τ^c × cov^d × exp(k×π/(T/d))
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
            if pa<=0 or am_cn<=0 or T<=0 or tau<=0: continue
            rows.append({'sel':sel,'pa':pa,'ps':ps,'am_cn':am_cn,'tau':tau,
                'cov':cov,'ratio':T/d_am,'name':mp.parent.name})
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
    rows=load_data();n=len(rows);print(f"n={n} (σ>0 only)\n")
    sel=np.array([r['sel'] for r in rows])
    pa=np.array([r['pa'] for r in rows]);ps=np.array([r['ps'] for r in rows])
    am_cn=np.array([r['am_cn'] for r in rows])
    tau=np.array([r['tau'] for r in rows]);cov=np.array([r['cov'] for r in rows])
    ratio=np.array([r['ratio'] for r in rows])
    names=[r['name'] for r in rows]

    thick=ratio>=10; thin=ratio<10

    # ═══════════════════════════════════════
    # v1 baseline
    # ═══════════════════════════════════════
    rhs_v1=SAM*pa**1.5*am_cn**2*np.exp(np.pi/ratio)
    C_v1=fitC(sel,rhs_v1);pred_v1=C_v1*rhs_v1
    print(f"v1: φ^1.5×CN²×exp(π/r)")
    print(f"  ALL={r2l(sel,pred_v1):.4f}  thick={r2l(sel[thick],pred_v1[thick]):.4f}  thin={r2l(sel[thin],pred_v1[thin]):.4f}\n")

    # ═══════════════════════════════════════
    # UNIFIED SWEEP: φ_AM^a × CN^b × τ^c × cov^d × exp(k×π/(T/d))
    # ═══════════════════════════════════════
    print("="*90)
    print("UNIFIED: φ_AM^a × CN^b × τ^c × cov^d × exp(k×π/(T/d))")
    print("="*90)

    results=[]
    for a in np.arange(1.5, 5.1, 0.5):
        for b in np.arange(1.0, 2.6, 0.5):
            for c in np.arange(-0.5, 1.6, 0.25):  # τ
                for d in np.arange(0, 1.6, 0.25):  # cov
                    for k in [0, 0.5, 1.0]:  # exp coefficient
                        if k > 0:
                            rhs=SAM*pa**a*am_cn**b*tau**c*cov**d*np.exp(k*np.pi/ratio)
                        else:
                            rhs=SAM*pa**a*am_cn**b*tau**c*cov**d
                        C=fitC(sel,rhs)
                        if C is None or C<=0: continue
                        pred=C*rhs
                        r2_all=r2l(sel,pred)
                        r2_thick=r2l(sel[thick],pred[thick]) if thick.sum()>=3 else -99
                        r2_thin=r2l(sel[thin],pred[thin]) if thin.sum()>=3 else -99
                        avg=(r2_thick+r2_thin)/2 if r2_thick>-50 and r2_thin>-50 else -99

                        if r2_all>0.88:
                            results.append({
                                'a':a,'b':b,'c':c,'d':d,'k':k,
                                'r2':r2_all,'thick':r2_thick,'thin':r2_thin,'avg':avg,
                                'C':C,'rhs':rhs,
                                'name':f'φ^{a}×CN^{b}×τ^{c}×cov^{d}×exp({k}π/r)' if k>0 else f'φ^{a}×CN^{b}×τ^{c}×cov^{d}'
                            })

    # Sort by avg(thick, thin)
    results.sort(key=lambda x:-x['avg'])

    print(f"\nTop 20 by avg(thick,thin):")
    seen=set();cnt=0
    for r in results:
        if r['name'] in seen: continue
        seen.add(r['name']);cnt+=1
        if cnt>20: break
        cv=loocv_C(sel,r['rhs']) if cnt<=3 else 0
        cv_str=f" LOOCV={cv:.4f}" if cv else ""
        flag='★' if r['avg']>0.93 else ' '
        print(f"  #{cnt:2d}{flag} avg={r['avg']:.3f} ALL={r['r2']:.4f} thick={r['thick']:.3f} thin={r['thin']:.3f}{cv_str}")
        print(f"      {r['name']}")

    # Top by ALL R²
    results.sort(key=lambda x:-x['r2'])
    print(f"\nTop 10 by ALL R²:")
    seen2=set();cnt2=0
    for r in results:
        if r['name'] in seen2: continue
        seen2.add(r['name']);cnt2+=1
        if cnt2>10: break
        print(f"  #{cnt2:2d} ALL={r['r2']:.4f} thick={r['thick']:.3f} thin={r['thin']:.3f}  {r['name']}")

    # ═══════════════════════════════════════
    # COMPARISON TABLE
    # ═══════════════════════════════════════
    print(f"\n{'='*90}")
    print("COMPARISON")
    print(f"{'='*90}")

    comparisons = [
        ('v1: φ^1.5×CN²×exp', SAM*pa**1.5*am_cn**2*np.exp(np.pi/ratio)),
    ]
    if results:
        best_avg = sorted(results, key=lambda x:-x['avg'])[0]
        comparisons.append((f'BEST avg: {best_avg["name"][:40]}', best_avg['rhs']))
        best_all = results[0]
        comparisons.append((f'BEST ALL: {best_all["name"][:40]}', best_all['rhs']))

    print(f"\n  {'Formula':50s} {'ALL':>7s} {'thick':>7s} {'thin':>7s} {'avg':>7s}")
    print("  "+"-"*80)
    for label, rhs in comparisons:
        C=fitC(sel,rhs);pred=C*rhs
        r2a=r2l(sel,pred);r2t=r2l(sel[thick],pred[thick]);r2n=r2l(sel[thin],pred[thin])
        avg=(r2t+r2n)/2
        print(f"  {label:50s} {r2a:7.4f} {r2t:7.3f} {r2n:7.3f} {avg:7.3f}")

if __name__=='__main__':
    main()
