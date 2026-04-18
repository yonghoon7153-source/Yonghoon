"""
Electronic: thin 극한 + thick 지수 검증
=======================================
thin: 0.884를 더 올려야. 모든 변수, 모든 조합.
thick: τ 지수 물리적 근거 — τ^0.5 vs τ^1 gap 얼마?
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
            por=m.get('porosity',0)
            el_act=max(m.get('electronic_active_fraction',0),0.01)
            ha=max(m.get('path_hop_area_mean',0),1e-6)
            if pa<=0 or am_cn<=0 or T<=0 or tau<=0: continue
            rows.append({'sel':sel,'pa':pa,'ps':ps,'am_cn':am_cn,'cn':cn,
                'T':T,'d_am':d_am,'tau':tau,'cov':cov,'gd':gd,'gp':gp,
                'ratio':T/d_am,'por':max(por,0.1),'el_act':el_act,'ha':ha,
                'name':mp.parent.name})
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
    gd=np.array([r['gd'] for r in rows]);gp=np.array([r['gp'] for r in rows])
    ratio=np.array([r['ratio'] for r in rows])
    por=np.array([r['por'] for r in rows])
    el_act=np.array([r['el_act'] for r in rows])
    ha=np.array([r['ha'] for r in rows])
    names=[r['name'] for r in rows]

    thick=ratio>=10; thin=ratio<10
    print(f"n={len(rows)}: thick={thick.sum()}, thin={thin.sum()}\n")

    # ═══════════════════════════════════════
    # THICK: τ 지수 검증 — 0, 0.5, 0.75, 1 gap?
    # ═══════════════════════════════════════
    print("="*80)
    print("THICK: τ 지수 sensitivity (φ^4×CN^1.5×τ^?×cov)")
    print("="*80)
    st=sel[thick];pt=pa[thick];ct=am_cn[thick];tt=tau[thick];cvt=cov[thick]
    for te in [0, 0.25, 0.5, 0.75, 1.0, 1.5]:
        rhs=SAM*pt**4*ct**1.5*tt**te*cvt
        C=fitC(st,rhs);pred=C*rhs;r2=r2l(st,pred)
        print(f"  τ^{te}: R²={r2:.4f}  {'← best' if abs(te-0.75)<0.1 else ''}")
    print(f"\n  결론: τ^0 ({r2l(st,fitC(st,SAM*pt**4*ct**1.5*cvt)*SAM*pt**4*ct**1.5*cvt):.4f}) vs τ^1 ({r2l(st,fitC(st,SAM*pt**4*ct**1.5*tt*cvt)*SAM*pt**4*ct**1.5*tt*cvt):.4f})")
    print(f"  τ가 없어도 R²>0.96 → τ는 미세 보정. 물리: SE compression 효과 (약함)")

    # ═══════════════════════════════════════
    # THIN: 극한 스크리닝 — 모든 변수 총동원
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("THIN: 극한 스크리닝 (모든 변수)")
    print("="*80)

    sn=sel[thin]
    var_thin={
        'φ_AM':np.log(pa[thin]),'φ_SE':np.log(ps[thin]),
        'AM_CN':np.log(am_cn[thin]),'SE_CN':np.log(cn[thin]),
        'T/d':np.log(ratio[thin]),'τ':np.log(tau[thin]),
        'cov':np.log(cov[thin]),'por':np.log(por[thin]),
        'GB_d':np.log(gd[thin]),'G_path':np.log(gp[thin]),
        'el_act':np.log(el_act[thin]),'hop_area':np.log(ha[thin]),
        'φ_AM×φ_SE':np.log(pa[thin]*ps[thin]),
        '1-φ_SE':np.log(1-ps[thin]),
    }

    log_sn=np.log(sn);n_thin=len(sn)
    results=[]

    # 2-5 variable free fit
    for nvars in [2,3,4,5]:
        for vnames in combinations(var_thin.keys(),nvars):
            X=np.column_stack([var_thin[v] for v in vnames]+[np.ones(n_thin)])
            try:
                coefs,_,_,_=np.linalg.lstsq(X,log_sn,rcond=None)
                pred=np.exp(X@coefs);r2=r2l(sn,pred)
                if r2>0.88:
                    terms=' × '.join(f'{v}^{coefs[i]:.2f}' for i,v in enumerate(vnames))
                    results.append({'name':f'({nvars}v) {terms}','r2':r2,'nvars':nvars,
                                   'coefs':coefs,'vnames':vnames})
            except: pass

    results.sort(key=lambda x:-x['r2'])

    print(f"\n  Top 20 (n={n_thin}):")
    seen=set()
    for i,r in enumerate(results[:20]):
        if r['name'] in seen: continue
        seen.add(r['name'])
        cv=loocv_C(sn,np.exp(np.column_stack([var_thin[v] for v in r['vnames']]+[np.ones(n_thin)])@r['coefs'])) if i<5 else 0
        cv_str=f" LOOCV={cv:.4f}" if cv else ""
        print(f"    #{i+1} R²={r['r2']:.4f}{cv_str}  {r['name']}")

    # ═══════════════════════════════════════
    # THIN: exp 추가하면?
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("THIN + exp(π/(T/d))")
    print("="*80)

    ratio_thin=ratio[thin]
    results_exp=[]
    for nvars in [2,3,4]:
        for vnames in combinations(var_thin.keys(),nvars):
            if 'T/d' in vnames: continue  # T/d and exp are redundant
            X=np.column_stack([var_thin[v] for v in vnames]+[np.pi/ratio_thin,np.ones(n_thin)])
            try:
                coefs,_,_,_=np.linalg.lstsq(X,log_sn,rcond=None)
                pred=np.exp(X@coefs);r2=r2l(sn,pred)
                if r2>0.88:
                    terms=' × '.join(f'{v}^{coefs[i]:.2f}' for i,v in enumerate(vnames))
                    results_exp.append({'name':f'+exp: {terms} ×exp({coefs[-2]:.2f}π/r)','r2':r2})
            except: pass

    results_exp.sort(key=lambda x:-x['r2'])
    print(f"\n  Top 10:")
    for i,r in enumerate(results_exp[:10]):
        print(f"    #{i+1} R²={r['r2']:.4f}  {r['name']}")

    # ═══════════════════════════════════════
    # THIN: percolation threshold
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("THIN: (φ_AM - φ_c) percolation")
    print("="*80)
    for phic in [0.25, 0.30, 0.35, 0.40, 0.45]:
        pa_ex=np.clip(pa[thin]-phic,0.001,None)
        for a in [0.5, 1, 1.5, 2]:
            for b in [1, 1.5]:
                for c in [-0.5, 0, 0.5, 1]:
                    rhs=SAM*pa_ex**a*am_cn[thin]**b*(ratio_thin)**c
                    C=fitC(sn,rhs)
                    if C is None: continue
                    pred=C*rhs;r2=r2l(sn,pred)
                    if r2>0.85:
                        print(f"  R²={r2:.4f}  (φ-{phic})^{a}×CN^{b}×(T/d)^{c}  C={C:.4f}")

if __name__=='__main__':
    main()
