"""
Electronic thin: 획기적 접근
============================
1. σ_ionic을 변수로 (SE network 품질 → AM 구조에 영향)
2. el_active, el_perc 직접 사용
3. σ_ionic × AM_CN × ... (ionic-electronic coupling)
4. 완전 자유 탐색 (σ_ionic 포함)
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
            sion=m.get('sigma_full_mScm',0)
            pa=max(m.get('phi_am',0),0.01);ps=m.get('phi_se',0)
            am_cn=max(m.get('am_am_cn',0.01),0.01);cn=m.get('se_se_cn',0)
            T=m.get('thickness_um',0);tau=m.get('tortuosity_recommended',m.get('tortuosity_mean',0))
            cov=max(m.get('coverage_AM_P_mean',m.get('coverage_AM_S_mean',m.get('coverage_AM_mean',20))),0.1)/100
            gd=m.get('gb_density_mean',0);gp=max(m.get('path_conductance_mean',0),1e-6)
            r_am=max(m.get('r_AM_P',0),m.get('r_AM_S',0))
            d_am=r_am*2 if r_am>0.1 else 5.0
            por=m.get('porosity',0)
            el_act=max(m.get('electronic_active_fraction',0),0.01)
            el_perc=max(m.get('electronic_percolating_fraction',0),0.01)
            bulk_frac=m.get('bulk_resistance_fraction',0.25)
            if pa<=0 or am_cn<=0 or T<=0: continue
            rows.append({'sel':sel,'sion':max(sion,0.001),'pa':pa,'ps':ps,'am_cn':am_cn,'cn':cn,
                'T':T,'d_am':d_am,'tau':max(tau,0.1),'cov':cov,'gd':gd,'gp':gp,
                'ratio':T/d_am,'por':max(por,0.1),'el_act':el_act,'el_perc':el_perc,
                'bulk_frac':max(bulk_frac,0.01),'name':mp.parent.name})
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
    sion=np.array([r['sion'] for r in rows])
    pa=np.array([r['pa'] for r in rows]);ps=np.array([r['ps'] for r in rows])
    am_cn=np.array([r['am_cn'] for r in rows]);cn=np.array([r['cn'] for r in rows])
    tau=np.array([r['tau'] for r in rows]);cov=np.array([r['cov'] for r in rows])
    gd=np.array([r['gd'] for r in rows]);gp=np.array([r['gp'] for r in rows])
    ratio=np.array([r['ratio'] for r in rows])
    por=np.array([r['por'] for r in rows])
    el_act=np.array([r['el_act'] for r in rows])
    el_perc=np.array([r['el_perc'] for r in rows])
    names=[r['name'] for r in rows]

    thick=ratio>=10; thin=ratio<10
    print(f"n={len(rows)}: thick={thick.sum()}, thin={thin.sum()}\n")

    # ═══════════════════════════════════════
    # THIN with σ_ionic
    # ═══════════════════════════════════════
    sn=sel[thin]; sion_n=sion[thin]; cn_n=am_cn[thin]; rn=ratio[thin]
    porn=por[thin]; cvn=cov[thin]; psn=ps[thin]; pan=pa[thin]
    taun=tau[thin]; gdn=gd[thin]; gpn=gp[thin]
    el_act_n=el_act[thin]; el_perc_n=el_perc[thin]
    n_thin=thin.sum()

    print("="*80)
    print("THIN: σ_ionic 포함 스크리닝")
    print("="*80)

    var_pool={
        'σ_ion':np.log(sion_n),'φ_AM':np.log(pan),'φ_SE':np.log(psn),
        'AM_CN':np.log(cn_n),'SE_CN':np.log(cn[thin]),
        'T/d':np.log(rn),'τ':np.log(taun),'cov':np.log(cvn),
        'por':np.log(porn),'GB_d':np.log(gdn),'G_path':np.log(gpn),
        'el_act':np.log(el_act_n),'el_perc':np.log(el_perc_n),
    }

    log_sn=np.log(sn)
    results=[]
    for nvars in [2,3,4]:
        for vnames in combinations(var_pool.keys(),nvars):
            X=np.column_stack([var_pool[v] for v in vnames]+[np.ones(n_thin)])
            try:
                coefs,_,_,_=np.linalg.lstsq(X,log_sn,rcond=None)
                pred=np.exp(X@coefs);r2=r2l(sn,pred)
                if r2>0.90:
                    terms=' × '.join(f'{v}^{coefs[i]:.2f}' for i,v in enumerate(vnames))
                    results.append({'name':f'({nvars}v) {terms}','r2':r2,'nvars':nvars,
                                   'coefs':coefs,'vnames':vnames})
            except: pass

    results.sort(key=lambda x:-x['r2'])
    print(f"\n  Top 20 (n={n_thin}):")
    for i,r in enumerate(results[:20]):
        cv=loocv_C(sn,np.exp(np.column_stack([var_pool[v] for v in r['vnames']]+[np.ones(n_thin)])@r['coefs'])) if i<5 else 0
        cv_str=f" LOOCV={cv:.4f}" if cv else ""
        has_sion='★' if 'σ_ion' in r['name'] else ' '
        has_el='●' if 'el_' in r['name'] else ' '
        print(f"    #{i+1}{has_sion}{has_el} R²={r['r2']:.4f}{cv_str}  {r['name']}")

    # ═══════════════════════════════════════
    # THICK with σ_ionic
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("THICK: σ_ionic 포함")
    print("="*80)

    st=sel[thick]; sion_t=sion[thick]; cn_t=am_cn[thick]
    cvt=cov[thick]; pt=pa[thick]; taut=tau[thick]

    var_thick={
        'σ_ion':np.log(sion_t),'φ_AM':np.log(pt),
        'AM_CN':np.log(cn_t),'τ':np.log(taut),'cov':np.log(cvt),
    }

    results_t=[]
    for nvars in [2,3,4]:
        for vnames in combinations(var_thick.keys(),nvars):
            X=np.column_stack([var_thick[v] for v in vnames]+[np.ones(thick.sum())])
            try:
                coefs,_,_,_=np.linalg.lstsq(X,np.log(st),rcond=None)
                pred=np.exp(X@coefs);r2=r2l(st,pred)
                if r2>0.95:
                    terms=' × '.join(f'{v}^{coefs[i]:.2f}' for i,v in enumerate(vnames))
                    results_t.append({'name':f'({nvars}v) {terms}','r2':r2})
            except: pass

    results_t.sort(key=lambda x:-x['r2'])
    print(f"\n  Top 10:")
    for i,r in enumerate(results_t[:10]):
        has_sion='★' if 'σ_ion' in r['name'] else ' '
        print(f"    #{i+1}{has_sion} R²={r['r2']:.4f}  {r['name']}")

    # ═══════════════════════════════════════
    # σ_ionic 직접 사용 fixed exponents
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("σ_ionic × AM_CN^a × por^b (thin, fixed)")
    print("="*80)
    for a in [0.5, 1, 1.5, 2]:
        for b in [0, 1, 2, 3]:
            for c in [-0.5, 0, 0.5]:
                rhs=SAM*sion_n**a*cn_n**1*porn**b*rn**c
                C=fitC(sn,rhs)
                if C is None: continue
                r2=r2l(sn,C*rhs)
                if r2>0.90:
                    print(f"  R²={r2:.4f}  σ_ion^{a}×CN×por^{b}×(T/d)^{c}")

if __name__=='__main__':
    main()
