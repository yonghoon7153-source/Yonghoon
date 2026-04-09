"""
Electronic 2-Regime: thick/thin 각각 최적 + 두께효과 제거
========================================================
T/d ≥ 10: thick regime (exp 불필요, 순수 topology)
T/d < 10: thin regime (finite-size 지배)

각 regime에서 T/d 의존성을 제거하고 핵심 변수만 추출
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
            if pa<=0 or am_cn<=0 or T<=0 or tau<=0: continue
            rows.append({'sel':sel,'pa':pa,'ps':ps,'am_cn':am_cn,'cn':cn,
                'T':T,'d_am':d_am,'tau':tau,'cov':cov,'gd':gd,'gp':gp,
                'ratio':T/d_am,'por':por,'name':mp.parent.name})
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

def sweep(sel, var_pool, label, threshold=0.80):
    n=len(sel); log_sel=np.log(sel); results=[]
    # Fixed exponent grid
    var_names=list(var_pool.keys())
    for nvars in [2,3,4]:
        for vnames in combinations(var_names, nvars):
            X=np.column_stack([var_pool[v] for v in vnames]+[np.ones(n)])
            try:
                coefs,_,_,_=np.linalg.lstsq(X,log_sel,rcond=None)
                pred=np.exp(X@coefs);r2=r2l(sel,pred)
                if r2>threshold:
                    terms=' × '.join(f'{v}^{coefs[i]:.2f}' for i,v in enumerate(vnames))
                    results.append({'name':terms,'r2':r2,'C':np.exp(coefs[-1]),'nvars':nvars,'coefs':coefs,'vnames':vnames})
            except: pass
    results.sort(key=lambda x:-x['r2'])
    print(f"\n  Top 10 free fit ({label}, n={n}):")
    for i,r in enumerate(results[:10]):
        cv=loocv_C(sel,np.exp(np.column_stack([var_pool[v] for v in r['vnames']]+[np.ones(n)])@r['coefs'])) if i<3 else 0
        cv_str=f" LOOCV={cv:.4f}" if cv else ""
        print(f"    #{i+1} R²={r['r2']:.4f}{cv_str} ({r['nvars']}v)  {r['name']}")
    return results

def main():
    rows=load_data();n=len(rows)
    sel=np.array([r['sel'] for r in rows])
    pa=np.array([r['pa'] for r in rows]);ps=np.array([r['ps'] for r in rows])
    am_cn=np.array([r['am_cn'] for r in rows]);cn=np.array([r['cn'] for r in rows])
    tau=np.array([r['tau'] for r in rows]);cov=np.array([r['cov'] for r in rows])
    gd=np.array([r['gd'] for r in rows]);gp=np.array([r['gp'] for r in rows])
    ratio=np.array([r['ratio'] for r in rows])
    por=np.array([r['por'] for r in rows])
    names=[r['name'] for r in rows]

    thick=ratio>=10; thin=ratio<10
    print(f"n={n}: thick={thick.sum()}, thin={thin.sum()}\n")

    # ═══════════════════════════════════════
    # THICK (T/d ≥ 10): 두께효과 무시, 순수 topology
    # T, d_AM, T/d, exp 전부 제외
    # ═══════════════════════════════════════
    print("="*80)
    print("THICK REGIME (T/d ≥ 10): 두께효과 제거")
    print("변수: φ_AM, AM_CN, τ, cov, φ_SE, GB_d, porosity")
    print("="*80)

    var_thick = {
        'φ_AM': np.log(pa[thick]), 'AM_CN': np.log(am_cn[thick]),
        'τ': np.log(tau[thick]), 'cov': np.log(cov[thick]),
        'φ_SE': np.log(ps[thick]), 'GB_d': np.log(gd[thick]),
        'por': np.log(np.clip(por[thick],0.1,None)),
        'φ_AM×φ_SE': np.log(pa[thick]*ps[thick]),
        'SE_CN': np.log(cn[thick]),
    }
    res_thick = sweep(sel[thick], var_thick, "thick")

    # ═══════════════════════════════════════
    # THIN (T/d < 10): T/d 포함, finite-size
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("THIN REGIME (T/d < 10): finite-size 포함")
    print("변수: φ_AM, AM_CN, T/d, cov, φ_SE, τ, porosity")
    print("="*80)

    var_thin = {
        'φ_AM': np.log(pa[thin]), 'AM_CN': np.log(am_cn[thin]),
        'T/d': np.log(ratio[thin]), 'cov': np.log(cov[thin]),
        'φ_SE': np.log(ps[thin]), 'τ': np.log(tau[thin]),
        'por': np.log(np.clip(por[thin],0.1,None)),
        'GB_d': np.log(gd[thin]),
        'SE_CN': np.log(cn[thin]),
    }
    res_thin = sweep(sel[thin], var_thin, "thin", threshold=0.75)

    # ═══════════════════════════════════════
    # 2-REGIME 결합 테스트
    # ═══════════════════════════════════════
    if res_thick and res_thin:
        print(f"\n{'='*80}")
        print("2-REGIME COMBINED")
        print("="*80)

        # Best thick + best thin
        bt=res_thick[0]; bn=res_thin[0]

        # Predict all using regime-appropriate formula
        pred_combined = np.zeros(n)
        for i in range(n):
            if thick[i]:
                # Use thick formula
                X_i = np.array([np.log(eval(f"rows[{i}]['{k.split('^')[0]}']" if '×' not in k else "1")) for k in bt['vnames']])
                # Simpler: just reconstruct
                pass
            # Too complex to reconstruct generically, just report R² per regime

        print(f"\n  Thick best: R²={bt['r2']:.4f}  {bt['name']}")
        print(f"  Thin best:  R²={bn['r2']:.4f}  {bn['name']}")

        # Physical interpretation
        print(f"\n  === PHYSICAL INTERPRETATION ===")
        print(f"\n  THICK (T/d ≥ 10):")
        print(f"    σ_el = C × σ_AM × {bt['name']}")
        print(f"    → AM topology 지배: φ_AM^~4 × CN^~1.5")
        print(f"    → exp 불필요 (두께 효과 없음)")
        print(f"    → 논리: 충분한 층수에서 AM percolation은 안정적. 접촉 품질과 구조가 결정.")

        print(f"\n  THIN (T/d < 10):")
        print(f"    σ_el = C × σ_AM × {bn['name']}")
        print(f"    → finite-size 지배: T/d가 핵심")
        print(f"    → 논리: 적은 층수에서 AM spanning cluster 형성이 T/d에 민감.")

        print(f"\n  통합 표현:")
        print(f"    σ_el = σ_AM × {{ C_thick × f_thick(φ,CN,τ,cov)  if T/d ≥ 10")
        print(f"                   {{ C_thin × f_thin(φ,CN,T/d,...)  if T/d < 10")

    # ═══════════════════════════════════════
    # Per-case for best thick/thin
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("PER-CASE")
    print(f"{'='*80}")

    # Thick per-case with v1
    rhs_v1_t = SAM * pa[thick]**1.5 * am_cn[thick]**2 * np.exp(np.pi/ratio[thick])
    C_v1_t = fitC(sel[thick], rhs_v1_t)
    pred_v1_t = C_v1_t * rhs_v1_t
    print(f"\n  Thick cases (v1 R²={r2l(sel[thick],pred_v1_t):.4f}):")
    for i, idx in enumerate(np.where(thick)[0]):
        err = abs(sel[idx]-pred_v1_t[i])/sel[idx]*100
        if err > 30:
            print(f"    {names[idx]:35s} σ={sel[idx]:.2f} pred={pred_v1_t[i]:.2f} err={err:.0f}%")

    # Thin per-case
    rhs_v1_n = SAM * pa[thin]**1.5 * am_cn[thin]**2 * np.exp(np.pi/ratio[thin])
    C_v1_n = fitC(sel[thin], rhs_v1_n)
    pred_v1_n = C_v1_n * rhs_v1_n
    print(f"\n  Thin cases (v1 R²={r2l(sel[thin],pred_v1_n):.4f}):")
    names_thin = [names[i] for i in np.where(thin)[0]]
    for i, idx in enumerate(np.where(thin)[0]):
        err = abs(sel[idx]-pred_v1_n[i])/sel[idx]*100
        print(f"    {names[idx]:35s} σ={sel[idx]:.2f} pred={pred_v1_n[i]:.2f} err={err:.0f}% T/d={ratio[idx]:.1f} φ={pa[idx]:.3f}")

if __name__=='__main__':
    main()
