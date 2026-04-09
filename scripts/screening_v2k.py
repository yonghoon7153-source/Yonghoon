"""
v2.0k NO LIMITS: coverage 확정 + 2nd 추가변수 + outlier 공략
=============================================================
"""
import json, os, numpy as np, warnings
from pathlib import Path
warnings.filterwarnings('ignore')
WEBAPP = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'webapp')
SG = 3.0

def load_data():
    rows = []
    for base in [Path(WEBAPP)/'results', Path(WEBAPP)/'archive']:
        if not base.is_dir(): continue
        for mp in base.rglob('full_metrics.json'):
            try:
                with open(mp) as f: m = json.load(f)
            except: continue
            sn = m.get('sigma_full_mScm', 0)
            if not sn or sn < 0.001: continue
            ps=m.get('phi_se',0); pa=max(m.get('phi_am',0),0.01)
            t=m.get('tortuosity_recommended',m.get('tortuosity_mean',0))
            fp=max(m.get('percolation_pct',0)/100,0.5)
            cn=m.get('se_se_cn',0); gd=m.get('gb_density_mean',0)
            gp=max(m.get('path_conductance_mean',0),1e-6)
            ha=max(m.get('path_hop_area_mean',0),1e-6)
            bn=max(m.get('path_hop_area_min_mean',0),1e-6)
            T=m.get('thickness_um',0)
            sa=max(m.get('area_SE_SE_mean',1e-6),1e-6)
            cov=max(m.get('coverage_AM_P_mean',m.get('coverage_AM_S_mean',m.get('coverage_AM_mean',1))),0.1)
            stcv=max(m.get('stress_cv',0.1),0.1)
            por=m.get('porosity',0)
            am_cn=max(m.get('am_am_cn',0.01),0.01)
            if t<=0 or ps<=0 or cn<=0 or gd<=0 or T<=0: continue
            rows.append({'sn':sn,'ps':ps,'pa':pa,'tau':t,'fp':fp,'cn':cn,'gd':gd,
                'gp':gp,'ha':ha,'bn':bn,'T':T,'sa':sa,'cov':cov,'stcv':stcv,
                'por':por,'am_cn':am_cn,'name':mp.parent.name})
    seen=set(); u=[]
    for r in rows:
        k=f"{r['ps']:.4f}_{r['T']:.1f}_{r['tau']:.3f}"
        if k not in seen: seen.add(k); u.append(r)
    return u

def r2l(a,p):
    la,lp=np.log(a),np.log(p); return 1-np.sum((la-lp)**2)/np.sum((la-np.mean(la))**2)
def fitC(a,r):
    v=(r>0)&np.isfinite(r); return float(np.exp(np.mean(np.log(a[v]/r[v])))) if v.sum()>=3 else None
def loocv_C(sn, rhs):
    n=len(sn); la=np.log(sn); lr=np.log(rhs); errs=[]
    for i in range(n):
        m=np.ones(n,bool); m[i]=False
        C_loo=float(np.exp(np.mean(la[m]-lr[m])))
        errs.append((la[i]-np.log(C_loo*rhs[i]))**2)
    return 1-np.sum(errs)/np.sum((la-np.mean(la))**2)

def main():
    rows=load_data(); n=len(rows); print(f"n={n}\n")
    sn=np.array([r['sn'] for r in rows])
    ps=np.array([r['ps'] for r in rows]); pa=np.array([r['pa'] for r in rows])
    tau=np.array([r['tau'] for r in rows]); cn=np.array([r['cn'] for r in rows])
    gd=np.array([r['gd'] for r in rows]); gp=np.array([r['gp'] for r in rows])
    ha=np.array([r['ha'] for r in rows]); bn=np.array([r['bn'] for r in rows])
    T=np.array([r['T'] for r in rows]); fp=np.array([r['fp'] for r in rows])
    sa=np.array([r['sa'] for r in rows]); cov=np.array([r['cov'] for r in rows])
    stcv=np.array([r['stcv'] for r in rows]); por=np.array([r['por'] for r in rows])
    am_cn=np.array([r['am_cn'] for r in rows]); names_arr=[r['name'] for r in rows]
    thick=np.array([r['tau']<=1.5 for r in rows])
    mid=np.array([(r['tau']>1.5)&(r['tau']<=2.5) for r in rows])
    thin=np.array([r['tau']>2.5 for r in rows])
    sigma_brug = SG * ps * fp / tau**2

    # ═══════════════════════════════════════
    # PART 1: v4+coverage 위에 2nd 변수 추가
    # Base: σ_brug×√(1-0.18/φ)×τ^1.5/f×CN^1.25×cov^0.375
    # ═══════════════════════════════════════
    print("PART 1: v4+coverage + 2nd extra variable")
    print("="*95)

    R_perc = np.clip(1-0.18/ps, 0.001, None)**0.5
    base_rhs = sigma_brug * R_perc * tau**1.5 / fp * cn**1.25 * (cov/100)**0.375
    C_base = fitC(sn, base_rhs); pred_base = C_base*base_rhs
    r2_base = r2l(sn, pred_base)
    print(f"  Base R²={r2_base:.4f}\n")

    extras = {
        # Contact quality
        'G_path': gp, 'GB_d': gd, 'G×d²': gp*gd**2,
        'hop_area': ha, 'bottleneck': bn, 'SE_area': sa,
        'hop×d²': ha*gd**2, 'BN×d²': bn*gd**2,
        # Network
        'AM_CN': am_cn, 'CN/τ': cn/tau, 'CN²/τ': cn**2/tau,
        # Geometry
        'T': T, 'φ_AM': pa, 'φ_SE×φ_AM': ps*pa,
        'porosity': np.clip(por,0.1,None),
        'φ_SE/τ': ps/tau, '1-por': np.clip(1-por/100,0.1,None),
        # Stress
        'stress_CV': stcv,
        # Derived
        'τ×GB_d': tau*gd, 'φ_SE²': ps**2,
        'CN×φ_SE': cn*ps, 'cov×CN': (cov/100)*cn,
    }
    for ename, evar in extras.items():
        for exp in [-0.5, -0.25, -0.125, 0.125, 0.25, 0.5]:
            rhs2 = base_rhs * evar**exp
            C2 = fitC(sn, rhs2)
            if C2 is None or C2<=0: continue
            pred2 = C2*rhs2
            r2_2 = r2l(sn, pred2)
            delta = r2_2 - r2_base
            if abs(delta) > 0.001:
                r2t = r2l(sn[thick],pred2[thick]) if thick.sum()>=2 else -99
                r2m = r2l(sn[mid],pred2[mid]) if mid.sum()>=2 else -99
                r2n = r2l(sn[thin],pred2[thin]) if thin.sum()>=2 else -99
                sign = '+' if delta>0 else '-'
                print(f"  {sign} {ename:>12s}^{exp:+.3f}: ΔR²={delta:+.4f} ALL={r2_2:.4f} thick={r2t:.3f} mid={r2m:.3f} thin={r2n:.3f}")

    # ═══════════════════════════════════════
    # PART 2: Outlier 공략 — thin_8 variants 무엇이 다른가?
    # ═══════════════════════════════════════
    print(f"\n{'='*95}")
    print("PART 2: Outlier Analysis — worst 5 cases 특성")
    print("="*95)

    per_err = np.abs(sn - pred_base) / sn * 100
    worst5 = np.argsort(per_err)[-5:][::-1]
    best5 = np.argsort(per_err)[:5]

    print(f"\n  {'':35s} {'σ_act':>7s} {'σ_pred':>7s} {'err%':>6s} {'τ':>5s} {'φ_SE':>5s} {'CN':>4s} {'cov%':>5s} {'GB_d':>5s} {'G_p':>8s}")
    print(f"  WORST 5:")
    for idx in worst5:
        print(f"  {names_arr[idx]:35s} {sn[idx]:7.4f} {pred_base[idx]:7.4f} {per_err[idx]:5.1f}% {tau[idx]:5.2f} {ps[idx]:5.3f} {cn[idx]:4.1f} {cov[idx]:5.1f} {gd[idx]:5.2f} {gp[idx]:8.6f}")
    print(f"  BEST 5:")
    for idx in best5:
        print(f"  {names_arr[idx]:35s} {sn[idx]:7.4f} {pred_base[idx]:7.4f} {per_err[idx]:5.1f}% {tau[idx]:5.2f} {ps[idx]:5.3f} {cn[idx]:4.1f} {cov[idx]:5.1f} {gd[idx]:5.2f} {gp[idx]:8.6f}")

    # What's different about outliers?
    print(f"\n  Outlier vs Average comparison:")
    for label, var, arr in [('τ',tau,tau), ('φ_SE',ps,ps), ('CN',cn,cn),
                            ('cov%',cov,cov), ('GB_d',gd,gd), ('G_path',gp,gp),
                            ('φ_AM',pa,pa), ('f_perc',fp,fp)]:
        avg_all = np.mean(arr)
        avg_out = np.mean(arr[worst5])
        ratio = avg_out/avg_all if avg_all>0 else 0
        print(f"    {label:8s}: all_avg={avg_all:.3f}  outlier_avg={avg_out:.3f}  ratio={ratio:.2f}")

    # ═══════════════════════════════════════
    # PART 3: ULTIMATE 미세조정 (0.0625 단위)
    # v4+coverage 위에서 모든 지수 극한 탐색
    # ═══════════════════════════════════════
    print(f"\n{'='*95}")
    print("PART 3: ULTIMATE 미세조정 (v4+coverage)")
    print("="*95)

    results = []
    step = 0.125
    for tb in np.arange(1.0, 2.01, step):
        for cne in np.arange(1.0, 2.01, step):
            for cova in np.arange(0.125, 0.626, step):
                for pe in [0.25, 0.5, 0.75]:  # percolation exponent
                    rhs = sigma_brug * np.clip(1-0.18/ps,0.001,None)**pe * tau**tb / fp * cn**cne * (cov/100)**cova
                    C = fitC(sn, rhs)
                    if C is None or C<=0: continue
                    pred = C*rhs
                    r2a = r2l(sn, pred)
                    if r2a < 0.94: continue
                    r2t = r2l(sn[thick],pred[thick]) if thick.sum()>=2 else -99
                    r2m = r2l(sn[mid],pred[mid]) if mid.sum()>=2 else -99
                    r2n = r2l(sn[thin],pred[thin]) if thin.sum()>=2 else -99
                    avg3 = (r2t+r2m+r2n)/3
                    results.append({
                        'pe':pe,'tb':tb,'cne':cne,'cova':cova,
                        'r2':r2a,'thick':r2t,'mid':r2m,'thin':r2n,'avg3':avg3,
                        'C':C,'rhs':rhs})

    results.sort(key=lambda x: -x['avg3'])

    print(f"\n  Top 10 by avg(thick,mid,thin):")
    for i,r in enumerate(results[:10]):
        eff_tau = r['tb']-2
        cv = loocv_C(sn, r['rhs']) if i<3 else 0
        print(f"\n  #{i+1} avg={r['avg3']:.3f} ALL={r['r2']:.4f} C={r['C']:.4f} {'LOOCV='+f'{cv:.4f}' if cv else ''}")
        print(f"     perc^{r['pe']} × τ^{r['tb']} × CN^{r['cne']} × cov^{r['cova']} [eff τ^{eff_tau:.3f}]")
        print(f"     thick={r['thick']:.3f} mid={r['mid']:.3f} thin={r['thin']:.3f}")

    # ═══════════════════════════════════════
    # PART 4: CLEAN FRACTIONS 테스트
    # Best를 깔끔한 분수로 고정
    # ═══════════════════════════════════════
    print(f"\n{'='*95}")
    print("PART 4: Clean fraction candidates")
    print("="*95)

    clean = [
        ('√perc × τ^(3/2) × CN^(5/4) × cov^(3/8)',   0.5, 1.5, 1.25, 0.375),
        ('√perc × τ^(3/2) × CN^(5/4) × ⁴√cov',       0.5, 1.5, 1.25, 0.25),
        ('√perc × τ^(3/2) × CN^(3/2) × ⁴√cov',       0.5, 1.5, 1.5,  0.25),
        ('√perc × τ^(5/4) × CN^(5/4) × ⁴√cov',       0.5, 1.25, 1.25, 0.25),
        ('√perc × τ^(5/4) × CN^(5/4) × cov^(3/8)',    0.5, 1.25, 1.25, 0.375),
        ('perc^(3/4) × τ^(3/2) × CN^(5/4) × ⁴√cov',  0.75, 1.5, 1.25, 0.25),
        ('√perc × τ^(3/2) × CN × ⁴√cov',              0.5, 1.5, 1.0,  0.25),
        ('√perc × τ^(3/2) × CN^(3/2) × cov^(3/8)',    0.5, 1.5, 1.5,  0.375),
        # Without coverage (v4 champion for comparison)
        ('√perc × τ^(3/2) × CN^(3/2) [NO COV]',       0.5, 1.5, 1.5,  0),
    ]

    print(f"\n  {'Formula':55s} {'ALL':>6s} {'LOOCV':>6s} {'thick':>6s} {'mid':>6s} {'thin':>6s} {'avg3':>6s}")
    print("  "+"-"*100)
    for label, pe, tb, cne, cova in clean:
        rhs = sigma_brug * np.clip(1-0.18/ps,0.001,None)**pe * tau**tb / fp * cn**cne
        if cova > 0:
            rhs = rhs * (cov/100)**cova
        C = fitC(sn, rhs)
        if C is None: continue
        pred = C*rhs
        r2a = r2l(sn, pred)
        cv = loocv_C(sn, rhs)
        r2t = r2l(sn[thick],pred[thick]) if thick.sum()>=2 else -99
        r2m = r2l(sn[mid],pred[mid]) if mid.sum()>=2 else -99
        r2n = r2l(sn[thin],pred[thin]) if thin.sum()>=2 else -99
        avg3 = (r2t+r2m+r2n)/3
        flag = '★' if avg3>0.93 else ' '
        print(f"  {flag}{label:54s} {r2a:6.4f} {cv:6.4f} {r2t:6.3f} {r2m:6.3f} {r2n:6.3f} {avg3:6.3f}")

    # ═══════════════════════════════════════
    # FINAL: 최종 v4++ 공식
    # ═══════════════════════════════════════
    if results:
        best = results[0]
        R = np.clip(1-0.18/ps,0.001,None)**best['pe']
        rhs_final = sigma_brug * R * tau**best['tb'] / fp * cn**best['cne'] * (cov/100)**best['cova']
        C_final = fitC(sn, rhs_final)
        cv_final = loocv_C(sn, rhs_final)
        pred_final = C_final * rhs_final
        eff_tau = best['tb'] - 2

        print(f"\n{'='*95}")
        print(f"ABSOLUTE BEST v4++")
        print(f"{'='*95}")
        print(f"""
  σ = σ_brug × C × (1-0.18/φ)^{best['pe']} × τ^{best['tb']} / f_perc × CN^{best['cne']} × coverage^{best['cova']}

  Effective:
  σ = C × σ_grain × φ_SE × (1-0.18/φ)^{best['pe']} × CN^{best['cne']} × cov^{best['cova']} / τ^{-eff_tau}

  C = {C_final:.4f}
  R² = {best['r2']:.4f}, LOOCV = {cv_final:.4f}
  thick = {best['thick']:.3f}, mid = {best['mid']:.3f}, thin = {best['thin']:.3f}
  avg = {best['avg3']:.3f}

  Mean |err| = {np.mean(np.abs(sn-pred_final)/sn*100):.1f}%
  Median |err| = {np.median(np.abs(sn-pred_final)/sn*100):.1f}%
  Within 20% = {np.sum(np.abs(sn-pred_final)/sn<0.2)}/{n}
  Within 30% = {np.sum(np.abs(sn-pred_final)/sn<0.3)}/{n}
        """)


if __name__ == '__main__':
    main()
