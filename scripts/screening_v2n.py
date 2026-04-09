"""
v2.0n FORM A BROKEN → NEW CHAMPION 탐색
========================================
(φ-φc)^1.5이 φ×(φ-φc)보다 좋다!
Free fit: φ_SE 지수 ≈ 0 → φ_SE 불필요, (φ-φc)만 필요

New candidates:
FORM X: σ = C × σ_grain × (φ-φc)^(3/4) × CN × √cov / √τ
FORM Y: σ = C × σ_grain × ⁴√[(φ-φc)³ × CN⁴ × cov² / τ²]

+ α optimization on (φ-φc)^1.5 base
+ φ_c re-optimization
+ free fit convergence test
+ stress test 전부 재실행
"""
import json,os,numpy as np,warnings
from pathlib import Path
from scipy.optimize import minimize_scalar
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
    rows=load_data();n=len(rows);print(f"n={n}\n")
    sn=np.array([r['sn'] for r in rows])
    ps=np.array([r['ps'] for r in rows]);pa=np.array([r['pa'] for r in rows])
    tau=np.array([r['tau'] for r in rows]);cn=np.array([r['cn'] for r in rows])
    gd=np.array([r['gd'] for r in rows]);gp=np.array([r['gp'] for r in rows])
    fp=np.array([r['fp'] for r in rows]);cov=np.array([r['cov'] for r in rows])
    T=np.array([r['T'] for r in rows]);names=[r['name'] for r in rows]
    thick=np.array([r['tau']<=1.5 for r in rows])
    mid=np.array([(r['tau']>1.5)&(r['tau']<=2.5) for r in rows])
    thin=np.array([r['tau']>2.5 for r in rows])

    def full_test(label, rhs):
        C=fitC(sn,rhs)
        if C is None: print(f"  {label}: FAILED"); return None
        pred=C*rhs; r2a=r2l(sn,pred)
        cv=loocv_C(sn,rhs)
        r2t=r2l(sn[thick],pred[thick]) if thick.sum()>=2 else -99
        r2m=r2l(sn[mid],pred[mid]) if mid.sum()>=2 else -99
        r2n=r2l(sn[thin],pred[thin]) if thin.sum()>=2 else -99
        avg3=(r2t+r2m+r2n)/3
        err=np.mean(np.abs(sn-pred)/sn*100)
        w20=np.sum(np.abs(sn-pred)/sn<0.2)
        print(f"  {label}")
        print(f"    R²={r2a:.4f} LOOCV={cv:.4f} C={C:.4f}")
        print(f"    thick={r2t:.3f} mid={r2m:.3f} thin={r2n:.3f} avg={avg3:.3f}")
        print(f"    |err|={err:.1f}% within20%={w20}/{n}")
        return {'r2':r2a,'cv':cv,'thick':r2t,'mid':r2m,'thin':r2n,'avg3':avg3,'C':C,'pred':pred,'rhs':rhs}

    # ═══════════════════════════════════════
    # PART 1: (φ-φc)^p base — p sweep with coverage
    # ═══════════════════════════════════════
    print("="*90)
    print("PART 1: (φ-φc)^p × CN^a × cov^b / τ^c — systematic sweep")
    print("="*90)

    results = []
    for phic in [0.16, 0.17, 0.18, 0.19, 0.20]:
        phi_ex = np.clip(ps - phic, 0.001, None)
        for p in [0.5, 0.75, 1.0, 1.25, 1.5]:  # (φ-φc) exponent
            for a in [0.75, 1.0, 1.25, 1.5]:  # CN
                for b in [0, 0.25, 0.5, 0.75]:  # cov
                    for c in [0.25, 0.5, 0.75, 1.0]:  # τ
                        rhs = SG * phi_ex**p * cn**a * cov**b / tau**c
                        C=fitC(sn,rhs)
                        if C is None or C<=0: continue
                        pred=C*rhs; r2a=r2l(sn,pred)
                        if r2a<0.94: continue
                        r2t=r2l(sn[thick],pred[thick]) if thick.sum()>=2 else -99
                        r2m=r2l(sn[mid],pred[mid]) if mid.sum()>=2 else -99
                        r2n=r2l(sn[thin],pred[thin]) if thin.sum()>=2 else -99
                        avg3=(r2t+r2m+r2n)/3
                        results.append({
                            'phic':phic,'p':p,'a':a,'b':b,'c':c,
                            'r2':r2a,'thick':r2t,'mid':r2m,'thin':r2n,'avg3':avg3,'C':C,
                            'rhs':rhs,
                            'name':f'(φ-{phic})^{p}×CN^{a}×cov^{b}/τ^{c}'})

    results.sort(key=lambda x: -x['avg3'])
    print(f"\nTop 15 by avg(thick,mid,thin):")
    seen=set();cnt=0
    for r in results:
        if r['name'] in seen: continue
        seen.add(r['name']); cnt+=1
        if cnt>15: break
        flag='★' if r['avg3']>0.94 else '●'
        print(f"  #{cnt:2d}{flag} avg={r['avg3']:.3f} ALL={r['r2']:.4f} thick={r['thick']:.3f} mid={r['mid']:.3f} thin={r['thin']:.3f}")
        print(f"      {r['name']}")

    # ═══════════════════════════════════════
    # PART 2: HEAD-TO-HEAD — FORM A vs NEW candidates
    # ═══════════════════════════════════════
    print(f"\n{'='*90}")
    print("PART 2: HEAD-TO-HEAD COMPARISON")
    print("="*90)

    phi_ex18 = np.clip(ps-0.18, 0.001, None)
    phi_ex19 = np.clip(ps-0.19, 0.001, None)

    # FORM A: √[φ(φ-φc)×CN²×cov/τ]
    print(f"\nFORM A: √[φ(φ-φc)×CN²×cov/τ]")
    rA = full_test("FORM A", SG*np.sqrt(ps*phi_ex18*cn**2*cov/tau))

    # FORM X: (φ-φc)^(3/4) × CN × √cov / √τ
    print(f"\nFORM X: (φ-φc)^(3/4) × CN × √cov / √τ")
    rX = full_test("FORM X", SG*phi_ex18**0.75*cn*np.sqrt(cov)/np.sqrt(tau))

    # FORM Y: ⁴√[(φ-φc)³×CN⁴×cov²/τ²]
    print(f"\nFORM Y: ⁴√[(φ-φc)³×CN⁴×cov²/τ²]  (= FORM X)")
    rY = full_test("FORM Y", SG*(phi_ex18**3*cn**4*cov**2/tau**2)**0.25)

    # FORM Z: (φ-φc)^1 × CN × √cov / √τ  (simpler)
    print(f"\nFORM Z: (φ-φc) × CN × √cov / √τ")
    rZ = full_test("FORM Z", SG*phi_ex18*cn*np.sqrt(cov)/np.sqrt(tau))

    # FORM W: (φ-φc)^(3/4) × CN^(3/2) × cov^(3/8) / τ^(1/2) (from v2k best)
    print(f"\nFORM W: (φ-φc)^(3/4) × CN^(3/2) × cov^(3/8) / √τ")
    rW = full_test("FORM W", SG*phi_ex18**0.75*cn**1.5*cov**0.375/np.sqrt(tau))

    # FORM V: (φ-φc)^1 × CN^(5/4) × √cov / √τ
    print(f"\nFORM V: (φ-φc) × CN^(5/4) × √cov / √τ")
    rV = full_test("FORM V", SG*phi_ex18*cn**1.25*np.sqrt(cov)/np.sqrt(tau))

    # φ_c=0.19 versions
    print(f"\n--- φ_c=0.19 versions ---")
    print(f"\nFORM X19: (φ-0.19)^(3/4) × CN × √cov / √τ")
    rX19 = full_test("FORM X19", SG*phi_ex19**0.75*cn*np.sqrt(cov)/np.sqrt(tau))

    # ═══════════════════════════════════════
    # PART 3: α optimization on (φ-φc)^1.5 base
    # ═══════════════════════════════════════
    print(f"\n{'='*90}")
    print("PART 3: [(φ-φc)^1.5 × CN² × cov / τ]^α — α sweep")
    print("="*90)

    inner_new = phi_ex18**1.5 * cn**2 * cov / tau
    for alpha in np.arange(0.35, 0.76, 0.025):
        rhs=SG*inner_new**alpha; C=fitC(sn,rhs)
        if C is None: continue
        pred=C*rhs; r2=r2l(sn,pred)
        r2t=r2l(sn[thick],pred[thick]) if thick.sum()>=2 else -99
        r2n=r2l(sn[thin],pred[thin]) if thin.sum()>=2 else -99
        flag='★' if abs(alpha-0.5)<0.01 else ' '
        print(f"  {flag}α={alpha:.3f}: R²={r2:.4f} thick={r2t:.3f} thin={r2n:.3f}")

    def neg_r2(alpha):
        rhs=SG*inner_new**alpha;C=fitC(sn,rhs)
        if C is None: return 0
        return -r2l(sn,C*rhs)
    res=minimize_scalar(neg_r2,bounds=(0.3,0.8),method='bounded')
    print(f"\n  Optimal α = {res.x:.4f} (R²={-res.fun:.4f})")

    # ═══════════════════════════════════════
    # PART 4: Free fit with (φ-φc) only (no φ_SE)
    # ═══════════════════════════════════════
    print(f"\n{'='*90}")
    print("PART 4: Free fit — (φ-φc), CN, cov, τ")
    print("="*90)

    X = np.column_stack([np.log(phi_ex18), np.log(cn), np.log(cov), np.log(tau), np.ones(n)])
    coefs,_,_,_ = np.linalg.lstsq(X, np.log(sn), rcond=None)
    pred_free=np.exp(X@coefs); r2_free=r2l(sn,pred_free)
    print(f"  Free fit: R²={r2_free:.4f}")
    print(f"    (φ-φc): {coefs[0]:.3f}")
    print(f"    CN:     {coefs[1]:.3f}")
    print(f"    cov:    {coefs[2]:.3f}")
    print(f"    τ:      {coefs[3]:.3f}")
    print(f"    C:      {np.exp(coefs[4]):.4f}")

    # Nearest clean fractions
    ideal_X = np.array([0.75, 1.0, 0.5, -0.5])  # FORM X
    ideal_Z = np.array([1.0, 1.0, 0.5, -0.5])   # FORM Z
    print(f"\n    Distance to FORM X (3/4,1,1/2,-1/2): {np.sqrt(np.sum((coefs[:4]-ideal_X)**2)):.3f}")
    print(f"    Distance to FORM Z (1,1,1/2,-1/2):   {np.sqrt(np.sum((coefs[:4]-ideal_Z)**2)):.3f}")

    # ═══════════════════════════════════════
    # PART 5: Leave-SE-size-out on best candidate
    # ═══════════════════════════════════════
    print(f"\n{'='*90}")
    print("PART 5: Leave-SE-size-out (FORM X)")
    print("="*90)

    groups = {
        'SE~0.5μm': gd>1.0, 'SE~1.0μm': (gd>0.6)&(gd<=1.0), 'SE~1.5μm': gd<=0.6
    }
    for gname,gmask in groups.items():
        if gmask.sum()<3 or (~gmask).sum()<3: continue
        train=~gmask
        rhs_tr=SG*np.clip(ps[train]-0.18,0.001,None)**0.75*cn[train]*np.sqrt(cov[train])/np.sqrt(tau[train])
        C_tr=fitC(sn[train],rhs_tr)
        rhs_te=SG*np.clip(ps[gmask]-0.18,0.001,None)**0.75*cn[gmask]*np.sqrt(cov[gmask])/np.sqrt(tau[gmask])
        pred_te=C_tr*rhs_te
        r2_te=r2l(sn[gmask],pred_te)
        err_te=np.mean(np.abs(sn[gmask]-pred_te)/sn[gmask]*100)
        print(f"  {gname}: C={C_tr:.4f} R²_test={r2_te:.4f} |err|={err_te:.1f}%")

    # ═══════════════════════════════════════
    # PART 6: Per-case error for FORM X
    # ═══════════════════════════════════════
    if rX:
        print(f"\n{'='*90}")
        print("PART 6: Per-case error (FORM X)")
        print("="*90)
        per_err=np.abs(sn-rX['pred'])/sn*100
        print(f"  Mean={np.mean(per_err):.1f}% Median={np.median(per_err):.1f}%")
        print(f"  Within 20%: {np.sum(per_err<20)}/{n}")
        print(f"  Within 30%: {np.sum(per_err<30)}/{n}")
        print(f"  Within 50%: {np.sum(per_err<50)}/{n}")
        print(f"\n  Worst 5:")
        for idx in np.argsort(per_err)[-5:][::-1]:
            reg='thick' if tau[idx]<=1.5 else ('mid' if tau[idx]<=2.5 else 'THIN')
            print(f"    {names[idx]:35s} err={per_err[idx]:.1f}% τ={tau[idx]:.2f} φ={ps[idx]:.3f} [{reg}]")

    # ═══════════════════════════════════════
    # FINAL TABLE
    # ═══════════════════════════════════════
    print(f"\n{'='*90}")
    print("EVOLUTION TABLE: v3 → v4 → FORM A → FORM X")
    print("="*90)

    sigma_brug=SG*ps*fp/tau**2
    forms = {
        'v3: σ_brug×(Gd²)^¼×CN²': sigma_brug*(gp*gd**2)**0.25*cn**2,
        'v4: σ_brug×√(1-φc/φ)×τ^1.5/f×CN^1.5': sigma_brug*np.clip(1-0.18/ps,0.001,None)**0.5*tau**1.5/fp*cn**1.5,
        'FORM A: √[φ(φ-φc)×CN²×cov/τ]': SG*np.sqrt(ps*phi_ex18*cn**2*cov/tau),
        'FORM X: (φ-φc)^¾×CN×√cov/√τ': SG*phi_ex18**0.75*cn*np.sqrt(cov)/np.sqrt(tau),
    }

    print(f"\n  {'Model':<45s} {'ALL':>6s} {'LOOCV':>6s} {'thick':>6s} {'mid':>6s} {'thin':>6s} {'avg3':>6s}")
    print("  "+"-"*95)
    for label,rhs in forms.items():
        C=fitC(sn,rhs)
        if C is None: continue
        pred=C*rhs;r2a=r2l(sn,pred);cv=loocv_C(sn,rhs)
        r2t=r2l(sn[thick],pred[thick]) if thick.sum()>=2 else -99
        r2m=r2l(sn[mid],pred[mid]) if mid.sum()>=2 else -99
        r2n=r2l(sn[thin],pred[thin]) if thin.sum()>=2 else -99
        avg3=(r2t+r2m+r2n)/3
        print(f"  {label:<45s} {r2a:6.4f} {cv:6.4f} {r2t:6.3f} {r2m:6.3f} {r2n:6.3f} {avg3:6.3f}")

    print(f"""
  ┌──────────────────────────────────────────────────────────────┐
  │  FORM X: σ = C × σ_grain × (φ_SE - φ_c)^(3/4) × CN        │
  │              × √coverage / √τ                               │
  │                                                              │
  │  = C × σ_grain × ⁴√[ (φ_SE - φ_c)³ × CN⁴ × cov² / τ² ]  │
  │                                                              │
  │  Physical meaning:                                           │
  │  (φ-φc)^(3/4): 3D percolation (t=3/2 → ^3/4 after √)      │
  │  CN: network connectivity (linear)                           │
  │  √cov: AM-SE interface quality                               │
  │  1/√τ: softened tortuosity                                   │
  │                                                              │
  │  Elegance: ⁴√[integer powers inside]                        │
  │  (φ-φc)³, CN⁴, cov², τ⁻² — all integers!                  │
  └──────────────────────────────────────────────────────────────┘
    """)

if __name__ == '__main__':
    main()
