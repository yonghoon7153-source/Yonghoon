"""
Electronic ULTIMATE: τ추가 + percolation + el_active + 극한 미세조정
=================================================================
partial corr(τ | φ_AM,CN) = +0.678 — τ가 빠져있다!
"""
import json,os,numpy as np,warnings
from pathlib import Path
from itertools import combinations
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
            gd=m.get('gb_density_mean',0);gp=max(m.get('path_conductance_mean',0),1e-6)
            r_am=max(m.get('r_AM_P',0),m.get('r_AM_S',0))
            d_am=r_am*2 if r_am>0.1 else 5.0
            el_act=max(m.get('electronic_active_fraction',0),0.01)
            el_perc=max(m.get('electronic_percolating_fraction',0),0.01)
            if pa<=0 or am_cn<=0 or T<=0 or tau<=0: continue
            rows.append({'sel':sel,'pa':pa,'ps':ps,'am_cn':am_cn,'cn':cn,
                'T':T,'d_am':d_am,'tau':tau,'cov':cov,'gd':gd,'gp':gp,
                'ratio':T/d_am,'el_act':el_act,'el_perc':el_perc,'name':mp.parent.name})
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
    gd=np.array([r['gd'] for r in rows]);gp=np.array([r['gp'] for r in rows])
    el_act=np.array([r['el_act'] for r in rows])
    el_perc=np.array([r['el_perc'] for r in rows])
    log_sel=np.log(sel)

    results=[]
    def test(name,rhs):
        C=fitC(sel,rhs)
        if C is None or C<=0: return
        pred=C*rhs;r2=r2l(sel,pred)
        if r2<0.90: return
        results.append({'name':name,'r2':r2,'C':C,'rhs':rhs})

    # Current best baseline
    rhs_cur=SAM*pa**2.5*am_cn**2*np.exp(np.pi/ratio)*ps**0.5*cov**0.25
    C_cur=fitC(sel,rhs_cur);r2_cur=r2l(sel,C_cur*rhs_cur)
    print(f"Current best: R²={r2_cur:.4f}\n")

    # ═══════════════════════════════════════
    # PART 1: +τ (partial corr=+0.678!)
    # τ↑ → electronic↑ (SE tortuous → more AM space)
    # ═══════════════════════════════════════
    print("PART 1: + τ^e to current best")
    print("="*80)
    for e in np.arange(-1.0, 2.1, 0.25):
        rhs=rhs_cur*tau**e
        C=fitC(sel,rhs)
        if C is None: continue
        pred=C*rhs;r2=r2l(sel,pred)
        delta=r2-r2_cur
        if abs(delta)>0.001:
            flag='★' if delta>0 else ' '
            print(f"  {flag}τ^{e:+.2f}: R²={r2:.4f} (Δ={delta:+.4f})")
        if r2>0.90:
            results.append({'name':f'best+τ^{e:.2f}','r2':r2,'C':C,'rhs':rhs})

    # ═══════════════════════════════════════
    # PART 2: Full sweep with τ
    # φ_AM^a × CN^b × exp(π/(T/d)) × φ_SE^c × cov^d × τ^e
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("PART 2: Full sweep including τ")
    print("="*80)
    for a in [2.0, 2.5, 3.0, 3.5]:
        for b in [1.5, 2.0]:
            for c in [0, 0.5]:
                for d in [0, 0.25]:
                    for e in [0, 0.25, 0.5, 0.75, 1.0]:
                        rhs=SAM*pa**a*am_cn**b*np.exp(np.pi/ratio)*ps**c*cov**d*tau**e
                        test(f'φ^{a}×CN^{b}×exp×φ_SE^{c}×cov^{d}×τ^{e}',rhs)

    # ═══════════════════════════════════════
    # PART 3: AM percolation threshold + τ
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("PART 3: (φ_AM-φc) + τ")
    print("="*80)
    for phic in [0.20, 0.25, 0.30]:
        pa_ex=np.clip(pa-phic,0.001,None)
        for a in [1, 1.5, 2, 2.5]:
            for b in [1.5, 2]:
                for e in [0, 0.5, 1]:
                    rhs=SAM*pa_ex**a*am_cn**b*np.exp(np.pi/ratio)*tau**e
                    test(f'(φ-{phic})^{a}×CN^{b}×exp×τ^{e}',rhs)
                    for c in [0.25, 0.5]:
                        rhs2=rhs*cov**c
                        test(f'(φ-{phic})^{a}×CN^{b}×exp×τ^{e}×cov^{c}',rhs2)

    # ═══════════════════════════════════════
    # PART 4: el_active × CN^b (direct)
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("PART 4: el_active based")
    print("="*80)
    for a in [0.5, 1, 1.5, 2]:
        for b in [1, 1.5, 2]:
            for c in [0, 0.5]:
                for e in [0, 0.5, 1]:
                    rhs=SAM*el_act**a*am_cn**b*cov**c*tau**e
                    test(f'el_act^{a}×CN^{b}×cov^{c}×τ^{e}',rhs)
                    rhs2=rhs*np.exp(np.pi/ratio)
                    test(f'el_act^{a}×CN^{b}×cov^{c}×τ^{e}×exp',rhs2)

    # ═══════════════════════════════════════
    # PART 5: Ionic FORM X parallel
    # σ_el = C × σ_AM × ⁴√[(φ_AM-φc)³ × CN⁴ × cov² × τ²]
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("PART 5: FORM X parallel for electronic")
    print("="*80)
    for phic in [0, 0.20, 0.25]:
        pa_ex=np.clip(pa-phic,0.001,None) if phic>0 else pa
        # ⁴√[(φ-φc)^a × CN^b × cov^c × τ^d] × exp^e
        for a in [2, 3, 4, 5]:
            for b in [3, 4, 5]:
                for c in [0, 1, 2]:
                    for d in [0, 1, 2]:
                        inner=pa_ex**a*am_cn**b*cov**c*tau**d
                        rhs=SAM*inner**0.25
                        test(f'⁴√[(φ-{phic})^{a}×CN^{b}×cov^{c}×τ^{d}]',rhs)
                        rhs2=rhs*np.exp(np.pi/ratio)
                        test(f'⁴√[(φ-{phic})^{a}×CN^{b}×cov^{c}×τ^{d}]×exp',rhs2)
                        rhs3=rhs*np.exp(0.5*np.pi/ratio)
                        test(f'⁴√[...]×exp(0.5π/r)',rhs3)

    # ═══════════════════════════════════════
    # PART 6: GB_d, G_path (SE network → AM structure)
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("PART 6: SE network variables")
    print("="*80)
    base=SAM*pa**2.5*am_cn**2*np.exp(np.pi/ratio)
    for vn,vv in [('GB_d',gd),('G_path',gp),('SE_CN',cn),('1/GB_d',1/gd)]:
        for e in [-0.5,-0.25,0.25,0.5]:
            rhs=base*vv**e
            test(f'base×{vn}^{e}',rhs)
            rhs2=rhs*tau**0.5
            test(f'base×{vn}^{e}×τ^0.5',rhs2)

    # ═══════════════════════════════════════
    # RESULTS
    # ═══════════════════════════════════════
    results.sort(key=lambda x:-x['r2'])

    print(f"\n{'='*90}")
    print(f"TOP 30 (n={n})")
    print(f"{'='*90}")
    seen=set();cnt=0
    for r in results:
        if r['name'] in seen: continue
        seen.add(r['name']);cnt+=1
        if cnt>30: break
        cv=loocv_C(sel,r['rhs']) if cnt<=5 else 0
        cv_str=f" LOOCV={cv:.4f}" if cv else ""
        delta=r['r2']-r2_cur
        print(f"  #{cnt:2d} R²={r['r2']:.4f} ({delta:+.4f}){cv_str}  {r['name']}  C={r['C']:.4f}")

    # Free fit with τ included
    print(f"\n{'='*90}")
    print("FREE FIT with τ")
    print(f"{'='*90}")
    X=np.column_stack([np.log(pa),np.log(am_cn),np.pi/ratio,np.log(ps),np.log(cov),np.log(tau),np.ones(n)])
    coefs,_,_,_=np.linalg.lstsq(X,log_sel,rcond=None)
    pred=np.exp(X@coefs);r2=r2l(sel,pred)
    print(f"  φ_AM:{coefs[0]:.2f} CN:{coefs[1]:.2f} exp_k:{coefs[2]:.2f} φ_SE:{coefs[3]:.2f} cov:{coefs[4]:.2f} τ:{coefs[5]:.2f}")
    print(f"  R²={r2:.4f}")

    # EVOLUTION TABLE
    print(f"\n{'='*90}")
    print("EVOLUTION")
    print(f"{'='*90}")
    formulas={
        'v1: φ^1.5×CN²×exp': SAM*pa**1.5*am_cn**2*np.exp(np.pi/ratio),
        'v2: +φ_SE+cov': SAM*pa**2.5*am_cn**2*np.exp(np.pi/ratio)*ps**0.5*cov**0.25,
    }
    if results:
        formulas[f'v3: {results[0]["name"][:35]}']=results[0]['rhs']

    for label,rhs in formulas.items():
        C=fitC(sel,rhs);pred=C*rhs;r2=r2l(sel,pred);cv=loocv_C(sel,rhs)
        print(f"  {label:45s} R²={r2:.4f} LOOCV={cv:.4f}")

if __name__=='__main__':
    main()
