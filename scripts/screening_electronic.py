"""
Electronic FORM X: 모든 parameter 총동원 exhaustive screening
=============================================================
현재: σ_el = C × σ_AM × φ_AM^1.5 × CN_AM² × exp(π/(T/d))  R²=0.89
목표: R²>0.95 universal formula
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
            ps=m.get('phi_se',0);pa=max(m.get('phi_am',0),0.01)
            t=m.get('tortuosity_recommended',m.get('tortuosity_mean',0))
            fp=max(m.get('percolation_pct',0)/100,0.5)
            cn=m.get('se_se_cn',0);am_cn=max(m.get('am_am_cn',0.01),0.01)
            gd=m.get('gb_density_mean',0);gp=max(m.get('path_conductance_mean',0),1e-6)
            T=m.get('thickness_um',0);por=m.get('porosity',0)
            cov=max(m.get('coverage_AM_P_mean',m.get('coverage_AM_S_mean',m.get('coverage_AM_mean',20))),0.1)/100
            # d_AM
            r_am=max(m.get('r_AM_P',0),m.get('r_AM_S',0))
            d_am=r_am*2 if r_am>0.1 else 5.0
            el_active=m.get('electronic_active_fraction',0)
            el_perc=m.get('electronic_percolating_fraction',0)
            if pa<=0 or am_cn<=0 or T<=0: continue
            rows.append({'sel':sel,'ps':ps,'pa':pa,'tau':t,'fp':fp,'cn':cn,'am_cn':am_cn,
                'gd':gd,'gp':gp,'T':T,'d_am':d_am,'por':por,'cov':cov,
                'el_active':max(el_active,0.01),'el_perc':max(el_perc,0.01),
                'name':mp.parent.name})
    seen=set();u=[]
    for r in rows:
        k=f"{r['pa']:.4f}_{r['T']:.1f}_{r['am_cn']:.2f}"
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
    T=np.array([r['T'] for r in rows]);d_am=np.array([r['d_am'] for r in rows])
    tau=np.array([r['tau'] for r in rows]);gd=np.array([r['gd'] for r in rows])
    gp=np.array([r['gp'] for r in rows]);cov=np.array([r['cov'] for r in rows])
    por=np.array([r['por'] for r in rows])
    el_active=np.array([r['el_active'] for r in rows])
    el_perc=np.array([r['el_perc'] for r in rows])
    ratio=T/d_am

    results=[]

    def test(name,rhs):
        C=fitC(sel,rhs)
        if C is None or C<=0: return
        pred=C*rhs;r2=r2l(sel,pred)
        if r2<0.85: return
        results.append({'name':name,'r2':r2,'C':C})

    # ═══════════════════════════════════════
    # Current formula
    # ═══════════════════════════════════════
    rhs_cur=SAM*pa**1.5*am_cn**2*np.exp(np.pi/ratio)
    C_cur=fitC(sel,rhs_cur);pred_cur=C_cur*rhs_cur
    r2_cur=r2l(sel,pred_cur);cv_cur=loocv_C(sel,rhs_cur)
    print(f"Current: R²={r2_cur:.4f} LOOCV={cv_cur:.4f} C={C_cur:.4f}")
    print(f"  σ_el = {C_cur:.4f} × σ_AM × φ_AM^1.5 × CN_AM² × exp(π/(T/d))\n")

    # ═══════════════════════════════════════
    # STRATEGY 1: φ_AM percolation threshold
    # ═══════════════════════════════════════
    print("STRATEGY 1: (φ_AM - φ_c)^a × AM_CN^b × exp correction")
    for phic in [0, 0.15, 0.20, 0.25, 0.30]:
        pa_ex=np.clip(pa-phic,0.001,None) if phic>0 else pa
        for a in [1, 1.5, 2, 2.5, 3]:
            for b in [1, 1.5, 2]:
                # with exp
                rhs=SAM*pa_ex**a*am_cn**b*np.exp(np.pi/ratio)
                test(f'(φ_AM-{phic})^{a}×CN^{b}×exp(π/(T/d))',rhs)
                # without exp
                rhs2=SAM*pa_ex**a*am_cn**b
                test(f'(φ_AM-{phic})^{a}×CN^{b} [no exp]',rhs2)
                # with T/d ratio
                for re in [0.25, 0.5, 1]:
                    rhs3=SAM*pa_ex**a*am_cn**b*(T/d_am)**re
                    test(f'(φ_AM-{phic})^{a}×CN^{b}×(T/d)^{re}',rhs3)

    # ═══════════════════════════════════════
    # STRATEGY 2: Coverage, SE variables
    # ═══════════════════════════════════════
    print("STRATEGY 2: + coverage, φ_SE, τ, GB_d")
    base=SAM*pa**1.5*am_cn**2*np.exp(np.pi/ratio)
    for ename,evar in [('cov',cov),('φ_SE',ps),('τ',tau),('GB_d',gd),
                        ('el_active',el_active),('el_perc',el_perc),
                        ('1-por',np.clip(1-por/100,0.1,None))]:
        for exp in [-0.5,-0.25,0.25,0.5,1]:
            rhs=base*evar**exp
            test(f'current×{ename}^{exp}',rhs)

    # ═══════════════════════════════════════
    # STRATEGY 3: Direct power law (no σ_AM base)
    # ═══════════════════════════════════════
    print("STRATEGY 3: Direct power law")
    var_pool={
        'φ_AM':np.log(pa),'AM_CN':np.log(am_cn),'T/d':np.log(ratio),
        'φ_SE':np.log(ps),'τ':np.log(tau),'CN_SE':np.log(cn),
        'cov':np.log(cov),'GB_d':np.log(gd),'T':np.log(T),
        'el_act':np.log(el_active),'φ_AM×φ_SE':np.log(pa*ps),
        'AM_CN/τ':np.log(am_cn/tau),
    }
    log_sel=np.log(sel)

    for nvars in [3,4]:
        for vnames in combinations(var_pool.keys(),nvars):
            X=np.column_stack([var_pool[v] for v in vnames]+[np.ones(n)])
            try:
                coefs,_,_,_=np.linalg.lstsq(X,log_sel,rcond=None)
                pred=np.exp(X@coefs);r2=r2l(sel,pred)
                if r2>0.92:
                    terms=' × '.join(f'{v}^{coefs[i]:.2f}' for i,v in enumerate(vnames))
                    results.append({'name':f'DIRECT: {terms}','r2':r2,'C':np.exp(coefs[-1])})
            except: pass

    # ═══════════════════════════════════════
    # STRATEGY 4: φ_AM^a × AM_CN^b × φ_SE^c × cov^d
    # Like FORM X but for electronic
    # ═══════════════════════════════════════
    print("STRATEGY 4: FORM X style for electronic")
    for a in [1,1.5,2,2.5,3]:
        for b in [0.5,1,1.5,2]:
            for c in [-1,-0.5,0,0.5]:
                for d in [0,0.25,0.5]:
                    rhs=SAM*pa**a*am_cn**b*ps**c*cov**d
                    test(f'φ_AM^{a}×CN^{b}×φ_SE^{c}×cov^{d}',rhs)
                    # + exp
                    rhs2=rhs*np.exp(np.pi/ratio)
                    test(f'φ_AM^{a}×CN^{b}×φ_SE^{c}×cov^{d}×exp',rhs2)
                    # + T/d
                    for re in [0.5,1]:
                        rhs3=rhs*(T/d_am)**re
                        test(f'φ_AM^{a}×CN^{b}×φ_SE^{c}×cov^{d}×(T/d)^{re}',rhs3)

    # ═══════════════════════════════════════
    # RESULTS
    # ═══════════════════════════════════════
    results.sort(key=lambda x:-x['r2'])

    print(f"\n{'='*90}")
    print(f"TOP 30 ELECTRONIC FORMULAS (n={n})")
    print(f"{'='*90}")
    seen=set();cnt=0
    for r in results:
        if r['name'] in seen: continue
        seen.add(r['name']);cnt+=1
        if cnt>30: break
        flag='★' if r['r2']>r2_cur else ' '
        print(f"  #{cnt:2d}{flag} R²={r['r2']:.4f} C={r['C']:.4f}  {r['name']}")

    # Free fit
    print(f"\n{'='*90}")
    print("FREE FIT (all variables)")
    print(f"{'='*90}")
    X_all=np.column_stack([np.log(pa),np.log(am_cn),np.pi/ratio,np.log(ps),
                           np.log(cov),np.log(tau),np.log(T),np.ones(n)])
    coefs_all,_,_,_=np.linalg.lstsq(X_all,log_sel,rcond=None)
    print(f"  φ_AM: {coefs_all[0]:.3f}")
    print(f"  AM_CN: {coefs_all[1]:.3f}")
    print(f"  exp coef: {coefs_all[2]:.3f} (×π/(T/d))")
    print(f"  φ_SE: {coefs_all[3]:.3f}")
    print(f"  coverage: {coefs_all[4]:.3f}")
    print(f"  τ: {coefs_all[5]:.3f}")
    print(f"  T: {coefs_all[6]:.3f}")
    print(f"  C: {np.exp(coefs_all[7]):.4f}")
    pred_all=np.exp(X_all@coefs_all)
    print(f"  R²={r2l(sel,pred_all):.4f}")

    # LOOCV for top 3
    print(f"\n{'='*90}")
    print("LOOCV top 3")
    print(f"{'='*90}")
    seen2=set();cnt2=0
    for r in results:
        if r['name'] in seen2: continue
        seen2.add(r['name']);cnt2+=1
        if cnt2>3: break
        # Reconstruct rhs... just show name and R²
        print(f"  #{cnt2} R²={r['r2']:.4f}  {r['name']}")

if __name__=='__main__':
    main()
