"""
Coverage definition test: max vs mean vs weighted
"""
import json,os,numpy as np,warnings
from pathlib import Path
warnings.filterwarnings('ignore')
WEBAPP=os.path.join(os.path.dirname(os.path.dirname(__file__)),'webapp')
SG=3.0

def load_data():
    rows=[]
    for base in [Path(WEBAPP)/'results',Path(WEBAPP)/'archive']:
        if not base.is_dir(): continue
        for mp in base.rglob('full_metrics.json'):
            try:
                with open(mp) as f: m=json.load(f)
            except: continue
            sn=m.get('sigma_full_mScm',0)
            if not sn or sn<0.001: continue
            ps=m.get('phi_se',0)
            tau=m.get('tortuosity_recommended',m.get('tortuosity_mean',0))
            cn=m.get('se_se_cn',0)
            if tau<=0 or ps<=0 or cn<=0: continue

            cov_p=m.get('coverage_AM_P_mean',0)
            cov_s=m.get('coverage_AM_S_mean',0)
            cov_all=m.get('coverage_AM_mean',0)
            phi_am_p=m.get('phi_am_P',0)
            phi_am_s=m.get('phi_am_S',0)
            phi_am=m.get('phi_am',0)

            # max (current)
            cov_max=max(cov_p, cov_s, cov_all, 0.1)/100

            # mean
            vals=[v for v in [cov_p, cov_s, cov_all] if v>0]
            cov_mean=(sum(vals)/len(vals))/100 if vals else 0.01

            # weighted by phi_am
            if phi_am_p>0 and phi_am_s>0 and cov_p>0 and cov_s>0:
                total=phi_am_p+phi_am_s
                cov_weighted=(phi_am_p/total*cov_p + phi_am_s/total*cov_s)/100
            elif cov_all>0:
                cov_weighted=cov_all/100
            else:
                cov_weighted=cov_max

            rows.append({'sn':sn,'ps':ps,'tau':tau,'cn':cn,
                'cov_max':cov_max,'cov_mean':cov_mean,'cov_weighted':cov_weighted,
                'cov_p':cov_p,'cov_s':cov_s,'cov_all':cov_all,
                'name':mp.parent.name})
    seen=set();u=[]
    for r in rows:
        k=f"{r['ps']:.4f}_{r['tau']:.3f}"
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
    rows=load_data();n=len(rows)
    print(f"n={n}\n")

    sn=np.array([r['sn'] for r in rows])
    ps=np.array([r['ps'] for r in rows])
    tau=np.array([r['tau'] for r in rows])
    cn=np.array([r['cn'] for r in rows])
    phi_ex=np.clip(ps-0.18,0.001,None)

    thick=np.array([r['tau']<=1.5 for r in rows])
    thin=np.array([r['tau']>2.5 for r in rows])

    print("="*80)
    print("FORM X with different coverage definitions")
    print("σ = C × σ_grain × (φ-φc)^¾ × CN × √cov / √τ")
    print("="*80)

    for label, cov_key in [('max(P,S,all)', 'cov_max'),
                            ('mean(P,S,all)', 'cov_mean'),
                            ('weighted(φ_P×P + φ_S×S)', 'cov_weighted')]:
        cov=np.array([r[cov_key] for r in rows])
        rhs=SG*phi_ex**0.75*cn*np.sqrt(cov)/np.sqrt(tau)
        C=fitC(sn,rhs)
        if C is None: print(f"  {label}: FAILED"); continue
        pred=C*rhs
        r2=r2l(sn,pred)
        cv=loocv_C(sn,rhs)
        r2t=r2l(sn[thick],pred[thick]) if thick.sum()>=2 else -99
        r2n=r2l(sn[thin],pred[thin]) if thin.sum()>=2 else -99
        err=np.mean(np.abs(sn-pred)/sn*100)
        w20=np.sum(np.abs(sn-pred)/sn<0.2)

        print(f"\n  {label}:")
        print(f"    R²={r2:.4f}  LOOCV={cv:.4f}  C={C:.4f}")
        print(f"    thick={r2t:.3f}  thin={r2n:.3f}")
        print(f"    |err|={err:.1f}%  within20%={w20}/{n}")

    # Coverage stats
    print(f"\n{'='*80}")
    print("Coverage values comparison")
    print("="*80)
    print(f"  {'Case':>35s} {'cov_P':>7s} {'cov_S':>7s} {'cov_all':>7s} {'max':>7s} {'mean':>7s} {'wgt':>7s}")
    for r in sorted(rows, key=lambda x: -x['sn'])[:15]:
        print(f"  {r['name']:>35s} {r['cov_p']:7.1f} {r['cov_s']:7.1f} {r['cov_all']:7.1f} {r['cov_max']*100:7.1f} {r['cov_mean']*100:7.1f} {r['cov_weighted']*100:7.1f}")

if __name__=='__main__':
    main()
