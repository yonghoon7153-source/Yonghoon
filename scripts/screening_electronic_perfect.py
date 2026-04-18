"""
Electronic: thick + thin 둘 다 무결점까지
==========================================
"""
import json,os,numpy as np,warnings
from pathlib import Path
from itertools import product
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
    names=[r['name'] for r in rows]

    thick=ratio>=10; thin=ratio<10
    n=len(rows)
    print(f"n={n}: thick={thick.sum()}, thin={thin.sum()}\n")

    # ═══════════════════════════════════════
    # THICK: 0.125 단위 미세조정
    # ═══════════════════════════════════════
    print("="*80)
    print("THICK: 0.125 단위 미세조정")
    print("="*80)

    st=sel[thick];pt=pa[thick];ct=am_cn[thick];tt=tau[thick];cvt=cov[thick]
    results_t=[]
    step=0.125
    for a in np.arange(3.5,5.1,step):
        for b in np.arange(1.0,2.1,step):
            for c in np.arange(0,1.1,step):
                for d in np.arange(0.5,1.5,step):
                    rhs=SAM*pt**a*ct**b*tt**c*cvt**d
                    C=fitC(st,rhs)
                    if C is None: continue
                    pred=C*rhs;r2=r2l(st,pred)
                    if r2>0.965:
                        results_t.append({'a':a,'b':b,'c':c,'d':d,'r2':r2,'C':C,'rhs':rhs})

    results_t.sort(key=lambda x:-x['r2'])
    print(f"\n  Top 10:")
    for i,r in enumerate(results_t[:10]):
        cv=loocv_C(st,r['rhs']) if i<3 else 0
        cv_str=f" LOOCV={cv:.4f}" if cv else ""
        print(f"    #{i+1} R²={r['r2']:.4f}{cv_str}  φ^{r['a']}×CN^{r['b']}×τ^{r['c']}×cov^{r['d']}  C={r['C']:.4f}")

    # Ablation: 각 항 제거
    if results_t:
        best_t=results_t[0]
        print(f"\n  Ablation (from best):")
        base_r2=best_t['r2']
        for label,rhs_abl in [
            ('Remove τ',SAM*pt**best_t['a']*ct**best_t['b']*cvt**best_t['d']),
            ('Remove cov',SAM*pt**best_t['a']*ct**best_t['b']*tt**best_t['c']),
            ('Remove CN',SAM*pt**best_t['a']*tt**best_t['c']*cvt**best_t['d']),
            ('Remove φ_AM',SAM*ct**best_t['b']*tt**best_t['c']*cvt**best_t['d']),
        ]:
            C=fitC(st,rhs_abl);r2=r2l(st,C*rhs_abl) if C else -99
            print(f"    {label:15s}: R²={r2:.4f} (Δ={r2-base_r2:+.4f})")

    # ═══════════════════════════════════════
    # THIN: 0.25 단위 미세조정 + 5변수
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("THIN: 0.25 단위 미세조정")
    print("="*80)

    sn=sel[thin];pn=pa[thin];psn=ps[thin];cn_n=am_cn[thin]
    rn=ratio[thin];cvn=cov[thin];tn=tau[thin];porn=por[thin];gdn=gd[thin]

    results_n=[]
    # Base: CN^a × (T/d)^b × por^c × cov^d × φ_SE^e
    for a in np.arange(0.75,2.1,0.25):
        for b in np.arange(-1.0,-0.1,0.25):
            for c in np.arange(1.0,3.1,0.25):
                for d in np.arange(0,1.5,0.25):
                    for e in np.arange(0,2.1,0.25):
                        rhs=SAM*cn_n**a*rn**b*porn**c*cvn**d*psn**e
                        C=fitC(sn,rhs)
                        if C is None: continue
                        pred=C*rhs;r2=r2l(sn,pred)
                        if r2>0.88:
                            results_n.append({'a':a,'b':b,'c':c,'d':d,'e':e,'r2':r2,'C':C,'rhs':rhs})

    results_n.sort(key=lambda x:-x['r2'])
    print(f"\n  Top 15:")
    for i,r in enumerate(results_n[:15]):
        cv=loocv_C(sn,r['rhs']) if i<3 else 0
        cv_str=f" LOOCV={cv:.4f}" if cv else ""
        print(f"    #{i+1} R²={r['r2']:.4f}{cv_str}  CN^{r['a']}×(T/d)^{r['b']}×por^{r['c']}×cov^{r['d']}×φ_SE^{r['e']}  C={r['C']:.6f}")

    # Ablation thin
    if results_n:
        best_n=results_n[0]
        print(f"\n  Ablation (from best):")
        base_r2n=best_n['r2']
        ablations=[
            ('Remove φ_SE',SAM*cn_n**best_n['a']*rn**best_n['b']*porn**best_n['c']*cvn**best_n['d']),
            ('Remove cov',SAM*cn_n**best_n['a']*rn**best_n['b']*porn**best_n['c']*psn**best_n['e']),
            ('Remove por',SAM*cn_n**best_n['a']*rn**best_n['b']*cvn**best_n['d']*psn**best_n['e']),
            ('Remove T/d',SAM*cn_n**best_n['a']*porn**best_n['c']*cvn**best_n['d']*psn**best_n['e']),
        ]
        for label,rhs_abl in ablations:
            C=fitC(sn,rhs_abl);r2=r2l(sn,C*rhs_abl) if C else -99
            print(f"    {label:15s}: R²={r2:.4f} (Δ={r2-base_r2n:+.4f})")

    # Per-case thin
    if results_n:
        best_n=results_n[0]
        C_n=best_n['C'];pred_n=C_n*best_n['rhs']
        err_n=np.abs(sn-pred_n)/sn*100
        print(f"\n  Per-case thin (best):")
        names_thin=[names[i] for i in np.where(thin)[0]]
        for i in np.argsort(-err_n):
            print(f"    {names_thin[i]:35s} σ={sn[i]:.2f} pred={pred_n[i]:.2f} err={err_n[i]:.0f}% T/d={rn[i]:.1f}")

    # ═══════════════════════════════════════
    # FINAL CHAMPION
    # ═══════════════════════════════════════
    if results_t and results_n:
        bt=results_t[0];bn=results_n[0]
        cv_t=loocv_C(st,bt['rhs']);cv_n=loocv_C(sn,bn['rhs'])
        err_t=np.mean(np.abs(st-bt['C']*bt['rhs'])/st*100)
        err_n=np.mean(np.abs(sn-bn['C']*bn['rhs'])/sn*100)
        w20_t=np.sum(np.abs(st-bt['C']*bt['rhs'])/st<0.2)
        w20_n=np.sum(np.abs(sn-bn['C']*bn['rhs'])/sn<0.2)

        print(f"\n{'='*80}")
        print("ELECTRONIC 2-REGIME CHAMPION")
        print(f"{'='*80}")
        print(f"""
  T/d ≥ 10 (THICK, n={thick.sum()}):
    σ_el = {bt['C']:.4f} × σ_AM × φ_AM^{bt['a']} × CN^{bt['b']} × τ^{bt['c']} × cov^{bt['d']}
    R²={bt['r2']:.4f}, LOOCV={cv_t:.4f}, |err|={err_t:.1f}%, within20%={w20_t}/{thick.sum()}

  T/d < 10 (THIN, n={thin.sum()}):
    σ_el = {bn['C']:.6f} × σ_AM × CN^{bn['a']} × (T/d)^{bn['b']} × por^{bn['c']} × cov^{bn['d']} × φ_SE^{bn['e']}
    R²={bn['r2']:.4f}, LOOCV={cv_n:.4f}, |err|={err_n:.1f}%, within20%={w20_n}/{thin.sum()}
        """)

if __name__=='__main__':
    main()
