"""
Champion Formula Comparison — 데이터 추가할 때마다 실행
=====================================================
v3, v4, FORM A, FORM X 전부 비교.
새 데이터 넣고 이것만 돌리면 어느 공식이 살아남는지 바로 확인.
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
            ps=m.get('phi_se',0);pa=max(m.get('phi_am',0),0.01)
            t=m.get('tortuosity_recommended',m.get('tortuosity_mean',0))
            fp=max(m.get('percolation_pct',0)/100,0.5)
            cn=m.get('se_se_cn',0);gd=m.get('gb_density_mean',0)
            gp=max(m.get('path_conductance_mean',0),1e-6)
            T=m.get('thickness_um',0)
            cov=max(m.get('coverage_AM_P_mean',m.get('coverage_AM_S_mean',m.get('coverage_AM_mean',1))),0.1)/100
            if t<=0 or ps<=0 or cn<=0 or gd<=0 or T<=0: continue
            rows.append({'sn':sn,'ps':ps,'pa':pa,'tau':t,'fp':fp,'cn':cn,'gd':gd,'gp':gp,'T':T,'cov':cov,'name':mp.parent.name})
    seen=set();u=[]
    for r in rows:
        k=f"{r['ps']:.4f}_{r['T']:.1f}_{r['tau']:.3f}"
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
    sn=np.array([r['sn'] for r in rows])
    ps=np.array([r['ps'] for r in rows]);pa=np.array([r['pa'] for r in rows])
    tau=np.array([r['tau'] for r in rows]);cn=np.array([r['cn'] for r in rows])
    gd=np.array([r['gd'] for r in rows]);gp=np.array([r['gp'] for r in rows])
    fp=np.array([r['fp'] for r in rows]);cov=np.array([r['cov'] for r in rows])
    thick=np.array([r['tau']<=1.5 for r in rows])
    mid=np.array([(r['tau']>1.5)&(r['tau']<=2.5) for r in rows])
    thin=np.array([r['tau']>2.5 for r in rows])
    sigma_brug=SG*ps*fp/tau**2
    phi_ex=np.clip(ps-0.18,0.001,None)

    champions = {
        'v3: σ_brug×(Gd²)^¼×CN²':
            sigma_brug*(gp*gd**2)**0.25*cn**2,
        'v4: σ_brug×√(1-φc/φ)×τ^1.5/f×CN^1.5':
            sigma_brug*np.clip(1-0.18/ps,0.001,None)**0.5*tau**1.5/fp*cn**1.5,
        'FORM A: √[φ(φ-φc)×CN²×cov/τ]':
            SG*np.sqrt(ps*phi_ex*cn**2*cov/tau),
        '★ FORM X: (φ-φc)^¾×CN×√cov/√τ':
            SG*phi_ex**0.75*cn*np.sqrt(cov)/np.sqrt(tau),
    }

    print(f"{'='*95}")
    print(f"CHAMPION COMPARISON — n={n} cases")
    print(f"τ range: {min(r['tau'] for r in rows):.2f}~{max(r['tau'] for r in rows):.2f}")
    print(f"thick(τ≤1.5): {thick.sum()}, mid(1.5<τ≤2.5): {mid.sum()}, thin(τ>2.5): {thin.sum()}")
    print(f"{'='*95}")

    print(f"\n{'Model':<45s} {'ALL':>6s} {'LOOCV':>6s} {'thick':>6s} {'mid':>6s} {'thin':>6s} {'avg3':>6s} {'|err|':>6s} {'C':>8s}")
    print("-"*95)

    for label,rhs in champions.items():
        C=fitC(sn,rhs)
        if C is None: print(f"  {label}: FAILED"); continue
        pred=C*rhs; r2a=r2l(sn,pred); cv=loocv_C(sn,rhs)
        r2t=r2l(sn[thick],pred[thick]) if thick.sum()>=2 else -99
        r2m=r2l(sn[mid],pred[mid]) if mid.sum()>=2 else -99
        r2n=r2l(sn[thin],pred[thin]) if thin.sum()>=2 else -99
        avg3=(r2t+r2m+r2n)/3 if min(r2t,r2m,r2n)>-50 else -99
        err=np.mean(np.abs(sn-pred)/sn*100)
        w20=np.sum(np.abs(sn-pred)/sn<0.2)
        print(f"  {label:<43s} {r2a:6.4f} {cv:6.4f} {r2t:6.3f} {r2m:6.3f} {r2n:6.3f} {avg3:6.3f} {err:5.1f}% {C:8.4f}")

    print(f"\n  Within 20% comparison:")
    for label,rhs in champions.items():
        C=fitC(sn,rhs)
        if C is None: continue
        pred=C*rhs; w20=np.sum(np.abs(sn-pred)/sn<0.2)
        print(f"    {label:<43s}: {w20}/{n} ({w20/n*100:.0f}%)")

    print(f"""
  ┌──────────────────────────────────────────────────────┐
  │ ★ CURRENT CHAMPION: FORM X                          │
  │                                                      │
  │ σ = C×σ_grain×(φ_SE-φ_c)^(3/4)×CN×√cov/√τ         │
  │   = C×σ_grain×⁴√[(φ-φc)³×CN⁴×cov²/τ²]             │
  │                                                      │
  │ C≈0.123, φ_c=0.18, σ_grain=3.0 mS/cm               │
  │ 1 free parameter                                     │
  └──────────────────────────────────────────────────────┘
    """)

if __name__=='__main__':
    main()
