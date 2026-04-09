"""
Electronic DEEP DIVE: 식 자체를 의심한다
=========================================
현재 가정: σ_el = f(φ_AM, CN_AM) × exp(π/(T/d))
의심:
- exp가 정말 필요한가? 아니면 다른 형태?
- SE가 AM에 미치는 구조적 영향은?
- φ_AM × φ_SE 같은 2상 경쟁?
- percolation threshold on φ_AM?
- electronic active fraction이 직접 변수?
"""
import json,os,numpy as np,warnings
from pathlib import Path
from itertools import combinations
from scipy.optimize import minimize_scalar
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
            cn=m.get('se_se_cn',0);am_cn=max(m.get('am_am_cn',0.01),0.01)
            gd=m.get('gb_density_mean',0);gp=max(m.get('path_conductance_mean',0),1e-6)
            T=m.get('thickness_um',0);por=m.get('porosity',0)
            cov=max(m.get('coverage_AM_P_mean',m.get('coverage_AM_S_mean',m.get('coverage_AM_mean',20))),0.1)/100
            r_am=max(m.get('r_AM_P',0),m.get('r_AM_S',0))
            d_am=r_am*2 if r_am>0.1 else 5.0
            el_active=m.get('electronic_active_fraction',0)
            el_perc=m.get('electronic_percolating_fraction',0)
            ha=max(m.get('path_hop_area_mean',0),1e-6)
            sa=max(m.get('area_SE_SE_mean',1e-6),1e-6)
            if pa<=0 or am_cn<=0 or T<=0: continue
            rows.append({'sel':sel,'ps':ps,'pa':pa,'tau':t,'cn':cn,'am_cn':am_cn,
                'gd':gd,'gp':gp,'T':T,'d_am':d_am,'por':por,'cov':cov,
                'el_active':max(el_active,0.01),'el_perc':max(el_perc,0.01),
                'ha':ha,'sa':sa,'ratio':T/d_am,'name':mp.parent.name})
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
    por=np.array([r['por'] for r in rows]);ratio=T/d_am
    el_active=np.array([r['el_active'] for r in rows])
    el_perc=np.array([r['el_perc'] for r in rows])
    log_sel=np.log(sel)

    # Correlation matrix
    print("="*80)
    print("CORRELATION with σ_el (log space)")
    print("="*80)
    vars_corr = {
        'φ_AM':np.log(pa),'φ_SE':np.log(ps),'AM_CN':np.log(am_cn),
        'SE_CN':np.log(cn),'τ':np.log(tau),'T':np.log(T),'T/d':np.log(ratio),
        'cov':np.log(cov),'GB_d':np.log(gd),'por':np.log(np.clip(por,0.1,None)),
        'el_act':np.log(el_active),'el_perc':np.log(el_perc),
        'φ_AM×φ_SE':np.log(pa*ps),'AM_CN/τ':np.log(am_cn/tau),
        'exp(π/r)':np.pi/ratio,
    }
    for vn,vv in sorted(vars_corr.items(), key=lambda x:-abs(np.corrcoef(x[1],log_sel)[0,1])):
        r=np.corrcoef(vv,log_sel)[0,1]
        print(f"  corr({vn:15s}, σ_el) = {r:+.3f}")

    results=[]
    def test(name,rhs):
        C=fitC(sel,rhs)
        if C is None or C<=0: return
        pred=C*rhs;r2=r2l(sel,pred)
        if r2<0.88: return
        results.append({'name':name,'r2':r2,'C':C,'rhs':rhs})

    # ═══════════════════════════════════════
    # PART 1: exp 없이 가능한가?
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("PART 1: exp(π/(T/d)) 없이")
    print("="*80)
    for a in [1.5,2,2.5,3,3.5]:
        for b in [1,1.5,2]:
            for c in [0,0.5,1]:
                for d in [0,0.25,0.5]:
                    rhs=SAM*pa**a*am_cn**b*ps**c*cov**d
                    test(f'NO_EXP: φ_AM^{a}×CN^{b}×φ_SE^{c}×cov^{d}',rhs)

    # ═══════════════════════════════════════
    # PART 2: exp를 다른 형태로
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("PART 2: exp 대신 (T/d)^k")
    print("="*80)
    for a in [1.5,2,2.5,3]:
        for b in [1.5,2]:
            for k in [0.25,0.5,0.75,1,1.5]:
                rhs=SAM*pa**a*am_cn**b*ratio**k
                test(f'RATIO: φ_AM^{a}×CN^{b}×(T/d)^{k}',rhs)
                for c in [0.25,0.5]:
                    rhs2=rhs*cov**c
                    test(f'RATIO: φ_AM^{a}×CN^{b}×(T/d)^{k}×cov^{c}',rhs2)
                for c in [0.25,0.5]:
                    rhs3=rhs*ps**c
                    test(f'RATIO: φ_AM^{a}×CN^{b}×(T/d)^{k}×φ_SE^{c}',rhs3)

    # ═══════════════════════════════════════
    # PART 3: φ_AM percolation threshold
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("PART 3: (φ_AM - φ_c)^a")
    print("="*80)
    for phic in np.arange(0.15,0.40,0.05):
        pa_ex=np.clip(pa-phic,0.001,None)
        for a in [1,1.5,2,2.5,3]:
            for b in [1,1.5,2]:
                rhs=SAM*pa_ex**a*am_cn**b*np.exp(np.pi/ratio)
                test(f'PERC: (φ_AM-{phic:.2f})^{a}×CN^{b}×exp',rhs)
                rhs2=SAM*pa_ex**a*am_cn**b
                test(f'PERC: (φ_AM-{phic:.2f})^{a}×CN^{b} [no exp]',rhs2)
                for k in [0.5,1]:
                    rhs3=SAM*pa_ex**a*am_cn**b*ratio**k
                    test(f'PERC: (φ_AM-{phic:.2f})^{a}×CN^{b}×(T/d)^{k}',rhs3)

    # ═══════════════════════════════════════
    # PART 4: IONIC FORM X 미러
    # σ_el = C × σ_AM × (φ_AM-φc)^a × AM_CN^b × √cov / √τ_AM?
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("PART 4: FORM X mirror for electronic")
    print("="*80)
    for phic in [0,0.20,0.25,0.30]:
        pa_ex=np.clip(pa-phic,0.001,None) if phic>0 else pa
        for a in [0.75,1,1.5,2]:
            for b in [0.5,1,1.5]:
                for c in [0,0.25,0.5]:
                    rhs=SAM*pa_ex**a*am_cn**b*cov**c
                    test(f'MIRROR: (φ_AM-{phic})^{a}×CN^{b}×cov^{c}',rhs)
                    # +φ_SE
                    for se in [0.25,0.5]:
                        rhs2=rhs*ps**se
                        test(f'MIRROR: (φ_AM-{phic})^{a}×CN^{b}×cov^{c}×φ_SE^{se}',rhs2)

    # ═══════════════════════════════════════
    # PART 5: el_active / el_perc 직접 사용
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("PART 5: el_active, el_perc 직접")
    print("="*80)
    for a in [0.5,1,1.5,2]:
        for b in [0.5,1,1.5]:
            rhs=SAM*el_active**a*am_cn**b
            test(f'ACTIVE: el_act^{a}×CN^{b}',rhs)
            rhs2=SAM*el_perc**a*am_cn**b
            test(f'PERC_F: el_perc^{a}×CN^{b}',rhs2)
            for c in [0.25,0.5]:
                rhs3=SAM*el_active**a*am_cn**b*cov**c
                test(f'ACTIVE: el_act^{a}×CN^{b}×cov^{c}',rhs3)

    # ═══════════════════════════════════════
    # PART 6: 2상 경쟁 — φ_AM × φ_SE
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("PART 6: 2상 경쟁")
    print("="*80)
    for a in [2,3,4,5]:
        for b in [0,1,2,3]:
            rhs=SAM*(pa**a)*(ps**b)*am_cn
            test(f'2PHASE: φ_AM^{a}×φ_SE^{b}×CN',rhs)
            rhs2=rhs*np.exp(np.pi/ratio)
            test(f'2PHASE: φ_AM^{a}×φ_SE^{b}×CN×exp',rhs2)

    # ═══════════════════════════════════════
    # PART 7: Free fit — 모든 변수
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("PART 7: Free fit combinations (3-5 vars)")
    print("="*80)
    var_pool={
        'φ_AM':np.log(pa),'AM_CN':np.log(am_cn),'φ_SE':np.log(ps),
        'cov':np.log(cov),'T/d':np.log(ratio),'τ':np.log(tau),
        'el_act':np.log(el_active),'GB_d':np.log(gd),
        'φ_AM²':2*np.log(pa),'AM_CN²':2*np.log(am_cn),
        'φ_AM×φ_SE':np.log(pa*ps),'φ_AM-φ_SE':np.log(pa)-np.log(ps),
    }
    for nvars in [3,4,5]:
        for vnames in combinations(var_pool.keys(),nvars):
            X=np.column_stack([var_pool[v] for v in vnames]+[np.ones(n)])
            try:
                coefs,_,_,_=np.linalg.lstsq(X,log_sel,rcond=None)
                pred=np.exp(X@coefs);r2=r2l(sel,pred)
                if r2>0.93:
                    terms=' × '.join(f'{v}^{coefs[i]:.2f}' for i,v in enumerate(vnames))
                    results.append({'name':f'FREE({nvars}v): {terms}','r2':r2,'C':np.exp(coefs[-1]),'rhs':np.exp(X@coefs)/np.exp(coefs[-1])})
            except: pass

    # ═══════════════════════════════════════
    # RESULTS
    # ═══════════════════════════════════════
    results.sort(key=lambda x:-x['r2'])

    print(f"\n{'='*90}")
    print(f"TOP 30 (n={n})")
    print(f"{'='*90}")
    cur_r2=0.8998
    seen=set();cnt=0
    for r in results:
        if r['name'] in seen: continue
        seen.add(r['name']);cnt+=1
        if cnt>30: break
        delta=r['r2']-cur_r2
        flag='★★' if r['r2']>0.93 else '★' if r['r2']>cur_r2 else ' '
        print(f"  #{cnt:2d}{flag} R²={r['r2']:.4f} ({delta:+.4f}) C={r['C']:.4f}  {r['name']}")

    # LOOCV for top 5
    print(f"\n{'='*90}")
    print("LOOCV top 5")
    print(f"{'='*90}")
    seen2=set();cnt2=0
    for r in results:
        if r['name'] in seen2 or 'rhs' not in r: continue
        seen2.add(r['name']);cnt2+=1
        if cnt2>5: break
        cv=loocv_C(sel,r['rhs'])
        print(f"  #{cnt2} R²={r['r2']:.4f} LOOCV={cv:.4f} gap={r['r2']-cv:.4f}  {r['name']}")

    # Partial correlations
    print(f"\n{'='*90}")
    print("PARTIAL CORRELATIONS (after controlling φ_AM + AM_CN)")
    print(f"{'='*90}")
    from numpy.linalg import lstsq
    X_base=np.column_stack([np.log(pa),np.log(am_cn),np.ones(n)])
    coef_base,_,_,_=lstsq(X_base,log_sel,rcond=None)
    resid_sel=log_sel-X_base@coef_base
    for vn,vv in vars_corr.items():
        if vn in ('φ_AM','AM_CN'): continue
        coef_v,_,_,_=lstsq(X_base,vv,rcond=None)
        resid_v=vv-X_base@coef_v
        pc=np.corrcoef(resid_v,resid_sel)[0,1]
        if abs(pc)>0.15:
            print(f"  partial({vn:15s} | φ_AM,CN) = {pc:+.3f}")

if __name__=='__main__':
    main()
