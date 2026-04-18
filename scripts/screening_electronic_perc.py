"""
Electronic: Finite-Size Percolation Threshold
==============================================
φ_c(T/d) = φ_c∞ + a × (T/d)^(-1/ν)
σ_el = C × σ_AM × (φ_AM - φ_c(T/d))^p × CN^q × ...

When φ_AM < φ_c(T/d): σ → 0 (percolation 실패)
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
            # Include σ=0 cases too!
            pa=max(m.get('phi_am',0),0.01);ps=m.get('phi_se',0)
            am_cn=max(m.get('am_am_cn',0.01),0.01);cn=m.get('se_se_cn',0)
            T=m.get('thickness_um',0);tau=m.get('tortuosity_recommended',m.get('tortuosity_mean',0))
            cov=max(m.get('coverage_AM_P_mean',m.get('coverage_AM_S_mean',m.get('coverage_AM_mean',20))),0.1)/100
            r_am=max(m.get('r_AM_P',0),m.get('r_AM_S',0))
            d_am=r_am*2 if r_am>0.1 else 5.0
            el_act=m.get('electronic_active_fraction',0) or 0
            el_perc=m.get('electronic_percolating_fraction',0) or 0
            if pa<=0 or T<=0: continue
            rows.append({'sel':max(sel,0),'pa':pa,'ps':ps,'am_cn':am_cn,'cn':cn,
                'T':T,'d_am':d_am,'tau':tau,'cov':cov,'ratio':T/d_am,
                'el_act':el_act,'el_perc':el_perc,'name':mp.parent.name})
    seen=set();u=[]
    for r in rows:
        k=f"{r['pa']:.4f}_{r['T']:.1f}"
        if k not in seen: seen.add(k);u.append(r)
    return u

def r2l(a,p):
    # For cases with σ=0, use linear R²
    ss_res=np.sum((a-p)**2); ss_tot=np.sum((a-np.mean(a))**2)
    return 1-ss_res/ss_tot if ss_tot>0 else -999
def r2l_log(a,p):
    mask=(a>0)&(p>0)
    if mask.sum()<3: return -999
    la,lp=np.log(a[mask]),np.log(p[mask])
    return 1-np.sum((la-lp)**2)/np.sum((la-np.mean(la))**2)
def fitC(a,r):
    mask=(a>0)&(r>0)&np.isfinite(r)
    if mask.sum()<3: return None
    return float(np.exp(np.mean(np.log(a[mask]/r[mask]))))

def main():
    rows=load_data();n=len(rows)
    sel=np.array([r['sel'] for r in rows])
    pa=np.array([r['pa'] for r in rows]);ps=np.array([r['ps'] for r in rows])
    am_cn=np.array([r['am_cn'] for r in rows])
    ratio=np.array([r['ratio'] for r in rows])
    cov=np.array([r['cov'] for r in rows])
    tau=np.array([r['tau'] for r in rows])
    el_perc=np.array([r['el_perc'] for r in rows])
    names=[r['name'] for r in rows]

    n_zero=np.sum(sel==0);n_pos=np.sum(sel>0)
    print(f"n={n} (σ>0: {n_pos}, σ=0: {n_zero})")
    print(f"φ_AM range: {pa.min():.3f}~{pa.max():.3f}")
    print(f"T/d range: {ratio.min():.1f}~{ratio.max():.1f}")
    print(f"el_perc=0 cases: {np.sum(el_perc==0)}\n")

    # Show which cases have σ=0
    print("Cases with σ_el=0 or el_perc=0:")
    for i in range(n):
        if sel[i]==0 or el_perc[i]==0:
            print(f"  {names[i]:35s} σ={sel[i]:.4f} φ_AM={pa[i]:.3f} T/d={ratio[i]:.1f} el_perc={el_perc[i]:.2f}")

    # ═══════════════════════════════════════
    # PART 1: φ_c(T/d) = φ_c∞ + a/(T/d)^b
    # Sweep φ_c∞, a, b
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("PART 1: φ_c(T/d) = φ_c∞ + a/(T/d)^b")
    print("="*80)

    results=[]
    for phic_inf in np.arange(0.15, 0.40, 0.05):
        for a_coef in np.arange(0.0, 3.1, 0.5):
            for b_exp in [0.5, 1.0, 1.5]:  # 1/ν ≈ 1.14 for 3D
                phic = phic_inf + a_coef / ratio**b_exp
                phi_excess = pa - phic

                # Cases where φ_AM < φ_c → predict σ=0
                pred_zero = phi_excess <= 0
                actual_zero = sel == 0

                # Accuracy of zero prediction
                correct_zero = np.sum(pred_zero & actual_zero)
                false_zero = np.sum(pred_zero & ~actual_zero)  # predicted 0 but actually >0
                missed_zero = np.sum(~pred_zero & actual_zero)  # predicted >0 but actually 0

                # For positive cases, fit scaling law
                pos_mask = (~pred_zero) & (sel > 0)
                if pos_mask.sum() < 5: continue

                phi_ex_pos = np.clip(phi_excess[pos_mask], 0.001, None)
                for p in [1, 1.5, 2]:
                    for q in [1, 1.5, 2]:
                        rhs = SAM * phi_ex_pos**p * am_cn[pos_mask]**q
                        C = fitC(sel[pos_mask], rhs)
                        if C is None or C <= 0: continue

                        # Full prediction (including zeros)
                        pred_full = np.zeros(n)
                        for i in range(n):
                            if phi_excess[i] > 0:
                                pred_full[i] = C * SAM * max(phi_excess[i],0.001)**p * am_cn[i]**q
                            # else: pred_full[i] = 0 (already)

                        r2_lin = r2l(sel, pred_full)
                        r2_log = r2l_log(sel, pred_full)

                        if r2_lin > 0.80:
                            results.append({
                                'phic_inf':phic_inf,'a':a_coef,'b':b_exp,'p':p,'q':q,
                                'r2_lin':r2_lin,'r2_log':r2_log,'C':C,
                                'correct_zero':correct_zero,'false_zero':false_zero,'missed_zero':missed_zero,
                                'name':f'φc={phic_inf}+{a_coef}/(T/d)^{b_exp} | (φ-φc)^{p}×CN^{q}',
                                'pred':pred_full,
                            })

    results.sort(key=lambda x:-x['r2_lin'])

    print(f"\nTop 15:")
    for i,r in enumerate(results[:15]):
        print(f"  #{i+1} R²_lin={r['r2_lin']:.4f} R²_log={r['r2_log']:.4f} C={r['C']:.4f}")
        print(f"     {r['name']}")
        print(f"     zero: correct={r['correct_zero']}, false={r['false_zero']}, missed={r['missed_zero']}")

    # ═══════════════════════════════════════
    # PART 2: Same but with exp(k×π/(T/d)) added
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("PART 2: + exp(k×π/(T/d))")
    print("="*80)

    results2=[]
    for phic_inf in [0.20, 0.25, 0.30]:
        for a_coef in [0, 0.5, 1.0, 1.5]:
            for b_exp in [1.0]:
                phic = phic_inf + a_coef / ratio**b_exp
                phi_excess = pa - phic
                pos_mask = (phi_excess > 0) & (sel > 0)
                if pos_mask.sum() < 5: continue
                phi_ex_pos = np.clip(phi_excess[pos_mask], 0.001, None)

                for p in [1, 1.5, 2]:
                    for q in [1.5, 2]:
                        for k in [0.5, 1.0]:
                            rhs = SAM * phi_ex_pos**p * am_cn[pos_mask]**q * np.exp(k*np.pi/ratio[pos_mask])
                            C = fitC(sel[pos_mask], rhs)
                            if C is None or C<=0: continue

                            pred_full = np.zeros(n)
                            for i in range(n):
                                if phi_excess[i] > 0:
                                    pred_full[i] = C*SAM*max(phi_excess[i],0.001)**p*am_cn[i]**q*np.exp(k*np.pi/ratio[i])

                            r2_lin = r2l(sel, pred_full)
                            r2_log = r2l_log(sel, pred_full)
                            pred_zero = phi_excess<=0; actual_zero = sel==0
                            cz=np.sum(pred_zero&actual_zero); fz=np.sum(pred_zero&~actual_zero)

                            if r2_lin > 0.80:
                                results2.append({
                                    'r2_lin':r2_lin,'r2_log':r2_log,'C':C,
                                    'correct_zero':cz,'false_zero':fz,
                                    'name':f'φc={phic_inf}+{a_coef}/(T/d) | (φ-φc)^{p}×CN^{q}×exp({k}π/r)',
                                    'pred':pred_full,
                                })

    results2.sort(key=lambda x:-x['r2_lin'])
    print(f"\nTop 10:")
    for i,r in enumerate(results2[:10]):
        print(f"  #{i+1} R²_lin={r['r2_lin']:.4f} R²_log={r['r2_log']:.4f}  zero:ok={r['correct_zero']},false={r['false_zero']}  {r['name']}")

    # ═══════════════════════════════════════
    # PART 3: Per-case error for best
    # ═══════════════════════════════════════
    best = (results + results2)
    best.sort(key=lambda x:-x['r2_lin'])
    if best:
        b = best[0]
        print(f"\n{'='*80}")
        print(f"BEST: {b['name']}")
        print(f"R²_lin={b['r2_lin']:.4f}, R²_log={b['r2_log']:.4f}")
        print(f"{'='*80}")
        pred = b['pred']
        print(f"\n{'Case':>35s} {'σ_act':>8s} {'σ_pred':>8s} {'err%':>7s}")
        for i in np.argsort(-np.abs(sel-pred)):
            err = abs(sel[i]-pred[i])/max(sel[i],0.001)*100 if sel[i]>0 else (0 if pred[i]==0 else 999)
            flag = '✓' if (sel[i]==0 and pred[i]==0) else ('✗' if (sel[i]==0)!=(pred[i]==0) else ' ')
            print(f"  {flag}{names[i]:34s} {sel[i]:8.2f} {pred[i]:8.2f} {err:6.1f}%")

    # ═══════════════════════════════════════
    # COMPARISON
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("COMPARISON: v1 vs finite-size percolation")
    print(f"{'='*80}")

    # v1 (no percolation check)
    rhs_v1 = SAM * pa**1.5 * am_cn**2 * np.exp(np.pi/ratio)
    C_v1 = fitC(sel, rhs_v1); pred_v1 = np.where(sel>0, C_v1*rhs_v1, C_v1*rhs_v1)
    # v1 always predicts >0, even for σ=0 cases
    r2_v1_lin = r2l(sel, pred_v1)
    r2_v1_log = r2l_log(sel, pred_v1)

    print(f"  v1 (no perc check):  R²_lin={r2_v1_lin:.4f} R²_log={r2_v1_log:.4f}")
    if best:
        print(f"  Best (perc check):   R²_lin={best[0]['r2_lin']:.4f} R²_log={best[0]['r2_log']:.4f}")

if __name__=='__main__':
    main()
