"""
Electronic & Thermal Scaling Law Stress Test
=============================================
Same methodology as ionic FORM X development.
Test: exponent optimization, percolation threshold, LOOCV, per-regime.
"""
import json,os,numpy as np,warnings
from pathlib import Path
from scipy.optimize import minimize_scalar
warnings.filterwarnings('ignore')
WEBAPP=os.path.join(os.path.dirname(os.path.dirname(__file__)),'webapp')
SG=3.0; SAM=50.0

def load_data():
    rows=[]
    for base in [Path(WEBAPP)/'results',Path(WEBAPP)/'archive']:
        if not base.is_dir(): continue
        for mp in base.rglob('full_metrics.json'):
            try:
                with open(mp) as f: m=json.load(f)
            except: continue
            sn=m.get('sigma_full_mScm',0)
            sel=m.get('electronic_sigma_full_mScm',0)
            sth=m.get('thermal_sigma_full_mScm',0)
            ps=m.get('phi_se',0);pa=max(m.get('phi_am',0),0.01)
            t=m.get('tortuosity_recommended',m.get('tortuosity_mean',0))
            cn=m.get('se_se_cn',0);am_cn=max(m.get('am_am_cn',0.01),0.01)
            T=m.get('thickness_um',0)
            # d_AM from stored r_AM values
            r_am=max(m.get('r_AM_P',0),m.get('r_AM_S',0))
            d_am=r_am*2 if r_am>0.1 else 5.0
            if t<=0 or ps<=0 or T<=0: continue
            rows.append({'sn':sn,'sel':sel,'sth':sth,'ps':ps,'pa':pa,'tau':t,
                'cn':cn,'am_cn':am_cn,'T':T,'d_am':d_am,'name':mp.parent.name})
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
    print(f"Loaded {n} cases\n")

    # ═══════════════════════════════════════
    # ELECTRONIC
    # ═══════════════════════════════════════
    el_rows=[r for r in rows if r['sel'] and r['sel']>0 and r['pa']>0 and r['am_cn']>0 and r['d_am']>0]
    print(f"{'='*80}")
    print(f"ELECTRONIC: {len(el_rows)} cases with σ_el data")
    print(f"{'='*80}")

    if len(el_rows)>=3:
        sel=np.array([r['sel'] for r in el_rows])
        pa=np.array([r['pa'] for r in el_rows])
        am_cn=np.array([r['am_cn'] for r in el_rows])
        T=np.array([r['T'] for r in el_rows])
        d_am=np.array([r['d_am'] for r in el_rows])
        tau_el=np.array([r['tau'] for r in el_rows])
        ratio=T/d_am

        # Current formula: σ_el = C × σ_AM × φ_AM^1.5 × CN_AM² × exp(π/(T/d))
        print("\n--- Current Formula ---")
        rhs_cur=SAM*pa**1.5*am_cn**2*np.exp(np.pi/ratio)
        C_cur=fitC(sel,rhs_cur);pred_cur=C_cur*rhs_cur
        r2_cur=r2l(sel,pred_cur);cv_cur=loocv_C(sel,rhs_cur)
        print(f"  σ_el = {C_cur:.4f} × σ_AM × φ_AM^1.5 × CN_AM² × exp(π/(T/d))")
        print(f"  R²={r2_cur:.4f}, LOOCV={cv_cur:.4f}, C={C_cur:.4f}")

        # Free fit
        print("\n--- Free fit ---")
        X=np.column_stack([np.log(pa),np.log(am_cn),np.pi/ratio,np.ones(len(el_rows))])
        coefs,_,_,_=np.linalg.lstsq(X,np.log(sel),rcond=None)
        pred_free=np.exp(X@coefs);r2_free=r2l(sel,pred_free)
        print(f"  φ_AM^{coefs[0]:.2f} × CN_AM^{coefs[1]:.2f} × exp({coefs[2]:.2f}×π/(T/d))")
        print(f"  R²={r2_free:.4f}")

        # Alternatives
        print("\n--- Exponent sweep ---")
        for pa_e in [1, 1.5, 2]:
            for cn_e in [1, 1.5, 2]:
                for use_exp in [True, False]:
                    if use_exp:
                        rhs=SAM*pa**pa_e*am_cn**cn_e*np.exp(np.pi/ratio)
                        label=f'φ_AM^{pa_e}×CN^{cn_e}×exp(π/(T/d))'
                    else:
                        rhs=SAM*pa**pa_e*am_cn**cn_e
                        label=f'φ_AM^{pa_e}×CN^{cn_e} (no exp)'
                    C=fitC(sel,rhs)
                    if C is None: continue
                    pred=C*rhs;r2=r2l(sel,pred)
                    flag='★' if r2>r2_cur else ' '
                    print(f"  {flag}R²={r2:.4f} {label}")

        # Percolation threshold on φ_AM?
        print("\n--- φ_AM percolation threshold ---")
        for phic in [0.15, 0.20, 0.25, 0.30]:
            pa_ex=np.clip(pa-phic,0.001,None)
            rhs=SAM*pa_ex**1.5*am_cn**2*np.exp(np.pi/ratio)
            C=fitC(sel,rhs)
            if C is None: continue
            pred=C*rhs;r2=r2l(sel,pred)
            print(f"  φ_c={phic}: R²={r2:.4f}")

        # Per-case error
        print(f"\n--- Per-case (current formula) ---")
        errs=np.abs(sel-pred_cur)/sel*100
        print(f"  Mean={np.mean(errs):.1f}%, Within 30%: {np.sum(errs<30)}/{len(el_rows)}")
    else:
        print("  Insufficient electronic data")

    # ═══════════════════════════════════════
    # THERMAL
    # ═══════════════════════════════════════
    th_rows=[r for r in rows if r['sth'] and r['sth']>0 and r['sn'] and r['sn']>0.01 and r['cn']>0 and r['pa']>0]
    print(f"\n{'='*80}")
    print(f"THERMAL: {len(th_rows)} cases with σ_th data")
    print(f"{'='*80}")

    if len(th_rows)>=3:
        sth=np.array([r['sth'] for r in th_rows])
        sion=np.array([r['sn'] for r in th_rows])
        pa_th=np.array([r['pa'] for r in th_rows])
        cn_th=np.array([r['cn'] for r in th_rows])
        tau_th=np.array([r['tau'] for r in th_rows])

        # Current: σ_th = C × σ_ion^0.75 × φ_AM² / CN_SE
        print("\n--- Current Formula ---")
        rhs_cur=sion**0.75*pa_th**2/cn_th
        C_cur=fitC(sth,rhs_cur);pred_cur=C_cur*rhs_cur
        r2_cur=r2l(sth,pred_cur);cv_cur=loocv_C(sth,rhs_cur)
        print(f"  σ_th = {C_cur:.1f} × σ_ion^0.75 × φ_AM² / CN_SE")
        print(f"  R²={r2_cur:.4f}, LOOCV={cv_cur:.4f}, C={C_cur:.1f}")

        # Free fit
        print("\n--- Free fit ---")
        X=np.column_stack([np.log(sion),np.log(pa_th),np.log(cn_th),np.ones(len(th_rows))])
        coefs,_,_,_=np.linalg.lstsq(X,np.log(sth),rcond=None)
        pred_free=np.exp(X@coefs);r2_free=r2l(sth,pred_free)
        print(f"  σ_ion^{coefs[0]:.2f} × φ_AM^{coefs[1]:.2f} × CN^{coefs[2]:.2f}")
        print(f"  R²={r2_free:.4f}, C={np.exp(coefs[3]):.1f}")

        # Alternatives
        print("\n--- Exponent sweep ---")
        for si_e in [0.5, 0.75, 1.0]:
            for pa_e in [1, 1.5, 2, 2.5]:
                for cn_e in [-2, -1.5, -1, -0.5, 0]:
                    rhs=sion**si_e*pa_th**pa_e*cn_th**cn_e
                    C=fitC(sth,rhs)
                    if C is None: continue
                    pred=C*rhs;r2=r2l(sth,pred)
                    if r2>r2_cur-0.01:
                        flag='★' if r2>r2_cur else ' '
                        print(f"  {flag}R²={r2:.4f} σ_ion^{si_e}×φ_AM^{pa_e}×CN^{cn_e}")

        # CN sign verification
        print("\n--- CN sign: does it flip for thermal? ---")
        from numpy.linalg import lstsq
        X_phi=np.column_stack([np.log(sion),np.log(pa_th),np.ones(len(th_rows))])
        coef_cn,_,_,_=lstsq(X_phi,np.log(cn_th),rcond=None)
        coef_sig,_,_,_=lstsq(X_phi,np.log(sth),rcond=None)
        resid_cn=np.log(cn_th)-X_phi@coef_cn
        resid_sig=np.log(sth)-X_phi@coef_sig
        partial=np.corrcoef(resid_cn,resid_sig)[0,1]
        print(f"  Partial corr(CN, σ_th | σ_ion, φ_AM) = {partial:.3f}")
        print(f"  {'NEGATIVE → CN hurts thermal (as expected)' if partial<0 else 'POSITIVE → CN helps thermal'}")

        # Per-case
        errs=np.abs(sth-pred_cur)/sth*100
        print(f"\n--- Per-case (current formula) ---")
        print(f"  Mean={np.mean(errs):.1f}%, Within 30%: {np.sum(errs<30)}/{len(th_rows)}")

    print(f"\n{'='*80}")
    print("SUMMARY")
    print(f"{'='*80}")
    print("Run this after adding new data to verify Electronic/Thermal formulas still hold.")

if __name__=='__main__':
    main()
