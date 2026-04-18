"""
Electronic 2-regime: 깔끔한 지수로 고정
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
            gd=m.get('gb_density_mean',0)
            r_am=max(m.get('r_AM_P',0),m.get('r_AM_S',0))
            d_am=r_am*2 if r_am>0.1 else 5.0
            por=m.get('porosity',0)
            if pa<=0 or am_cn<=0 or T<=0 or tau<=0: continue
            rows.append({'sel':sel,'pa':pa,'ps':ps,'am_cn':am_cn,'cn':cn,
                'T':T,'d_am':d_am,'tau':tau,'cov':cov,'gd':gd,
                'ratio':T/d_am,'por':max(por,0.1),'name':mp.parent.name})
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
    am_cn=np.array([r['am_cn'] for r in rows]);cn=np.array([r['cn'] for r in rows])
    tau=np.array([r['tau'] for r in rows]);cov=np.array([r['cov'] for r in rows])
    gd=np.array([r['gd'] for r in rows])
    ratio=np.array([r['ratio'] for r in rows])
    por=np.array([r['por'] for r in rows])

    thick=ratio>=10; thin=ratio<10
    n=len(rows)
    print(f"n={n}: thick={thick.sum()}, thin={thin.sum()}\n")

    # ═══════════════════════════════════════
    # THICK: φ_AM^a × CN^b × τ^c × cov^d
    # ═══════════════════════════════════════
    print("="*80)
    print("THICK (T/d ≥ 10): clean exponents")
    print("="*80)

    results_t=[]
    st=sel[thick]; pt=pa[thick]; ct=am_cn[thick]; tt=tau[thick]; cvt=cov[thick]

    for a in [3, 3.5, 4, 4.5, 5]:
        for b in [1, 1.5, 2]:
            for c in [0, 0.5, 0.75, 1]:
                for d in [0, 0.5, 0.75, 1]:
                    rhs=SAM*pt**a*ct**b*tt**c*cvt**d
                    C=fitC(st,rhs)
                    if C is None: continue
                    pred=C*rhs;r2=r2l(st,pred)
                    if r2>0.95:
                        results_t.append({'a':a,'b':b,'c':c,'d':d,'r2':r2,'C':C,
                            'name':f'φ^{a}×CN^{b}×τ^{c}×cov^{d}'})

    results_t.sort(key=lambda x:-x['r2'])
    print(f"\n  {'Formula':40s} {'R²':>7s} {'C':>8s}")
    print("  "+"-"*60)
    for r in results_t[:15]:
        print(f"  {r['name']:40s} {r['r2']:7.4f} {r['C']:8.4f}")

    # ═══════════════════════════════════════
    # THIN: AM_CN^a × (T/d)^b × var^c × var^d
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("THIN (T/d < 10): clean exponents")
    print("="*80)

    results_n=[]
    sn=sel[thin]; pn=pa[thin]; psn=ps[thin]; cn_n=am_cn[thin]
    rn=ratio[thin]; cvn=cov[thin]; tn=tau[thin]; porn=por[thin]; gdn=gd[thin]

    # CN^a × (T/d)^b × third^c
    thirds = [('φ_SE',psn),('por',porn),('cov',cvn),('τ',tn),('φ_AM',pn),('GB_d',gdn)]
    for a in [1, 1.5, 2]:
        for b in [-1, -0.75, -0.5, -0.25]:
            for tname,tvar in thirds:
                for c in [-1, -0.5, 0, 0.5, 1, 1.5, 2]:
                    rhs=SAM*cn_n**a*rn**b*tvar**c
                    C=fitC(sn,rhs)
                    if C is None: continue
                    pred=C*rhs;r2=r2l(sn,pred)
                    if r2>0.83:
                        results_n.append({'r2':r2,'C':C,
                            'name':f'CN^{a}×(T/d)^{b}×{tname}^{c}'})
                    # + fourth variable
                    for fname,fvar in thirds:
                        if fname==tname: continue
                        for d in [-0.5, 0.5, 1]:
                            rhs2=rhs*fvar**d
                            C2=fitC(sn,rhs2)
                            if C2 is None: continue
                            pred2=C2*rhs2;r2_2=r2l(sn,pred2)
                            if r2_2>0.87:
                                results_n.append({'r2':r2_2,'C':C2,
                                    'name':f'CN^{a}×(T/d)^{b}×{tname}^{c}×{fname}^{d}'})

    results_n.sort(key=lambda x:-x['r2'])
    seen=set()
    print(f"\n  {'Formula':50s} {'R²':>7s} {'C':>8s}")
    print("  "+"-"*70)
    cnt=0
    for r in results_n:
        if r['name'] in seen: continue
        seen.add(r['name']);cnt+=1
        if cnt>15: break
        print(f"  {r['name']:50s} {r['r2']:7.4f} {r['C']:8.4f}")

    # ═══════════════════════════════════════
    # FINAL: best clean for each regime
    # ═══════════════════════════════════════
    if results_t and results_n:
        bt=results_t[0]; bn=results_n[0]
        print(f"\n{'='*80}")
        print("ELECTRONIC 2-REGIME CHAMPION")
        print(f"{'='*80}")

        # LOOCV
        rhs_t=SAM*pt**bt['a']*ct**bt['b']*tt**bt['c']*cvt**bt['d']
        cv_t=loocv_C(st,rhs_t)

        print(f"""
  ┌─────────────────────────────────────────────────────────────┐
  │  T/d ≥ 10 (thick): R²={bt['r2']:.4f}, LOOCV={cv_t:.4f}                │
  │  σ_el = {bt['C']:.4f} × σ_AM × {bt['name']:30s}       │
  │                                                             │
  │  T/d < 10 (thin): R²={bn['r2']:.4f}                                │
  │  σ_el = {bn['C']:.4f} × σ_AM × {bn['name']:30s}          │
  │                                                             │
  │  물리:                                                      │
  │  Thick = AM topology 지배 (φ_AM, CN, τ, cov)               │
  │  Thin = finite-size 지배 (CN, T/d, por/φ_SE)               │
  └─────────────────────────────────────────────────────────────┘
        """)

if __name__=='__main__':
    main()
