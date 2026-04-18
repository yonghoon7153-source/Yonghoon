"""
v2.0i ABSOLUTE FINAL: 무결점 스크리닝
=====================================
v4_brug champion을 base로:
σ = σ_brug × C × √(1-φ_c/φ) × τ^1.5 / f_perc × CN^1.5

이 위에 모든 parameter를 하나씩 올려서 R² 개선되는지 확인.
+ 지수 미세 조정 (0.25 단위)
+ 모든 φ_c 값
+ LOOCV 전부
+ Per-case error 분석 (어느 케이스가 outlier인지)
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
            am_cn=max(m.get('am_am_cn',0.01),0.01)
            por=m.get('porosity',0)
            stcv=max(m.get('stress_cv',0.1),0.1)
            if t<=0 or ps<=0 or cn<=0 or gd<=0 or T<=0: continue
            rows.append({'sn':sn,'ps':ps,'pa':pa,'tau':t,'fp':fp,'cn':cn,'gd':gd,
                'gp':gp,'ha':ha,'bn':bn,'T':T,'sa':sa,'cov':cov,'am_cn':am_cn,
                'por':por,'stcv':stcv,'name':mp.parent.name})
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
    """LOOCV for 1-param (C only) model."""
    n=len(sn); la=np.log(sn); lr=np.log(rhs)
    errs=[]
    for i in range(n):
        m=np.ones(n,bool); m[i]=False
        C_loo=float(np.exp(np.mean(la[m]-lr[m])))
        errs.append((la[i]-np.log(C_loo*rhs[i]))**2)
    ss_r=np.sum(errs); ss_t=np.sum((la-np.mean(la))**2)
    return 1-ss_r/ss_t if ss_t>0 else -999

def main():
    rows = load_data()
    n = len(rows); print(f"n={n}\n")
    sn=np.array([r['sn'] for r in rows])
    ps=np.array([r['ps'] for r in rows]); pa=np.array([r['pa'] for r in rows])
    tau=np.array([r['tau'] for r in rows]); cn=np.array([r['cn'] for r in rows])
    gd=np.array([r['gd'] for r in rows]); gp=np.array([r['gp'] for r in rows])
    ha=np.array([r['ha'] for r in rows]); bn=np.array([r['bn'] for r in rows])
    T=np.array([r['T'] for r in rows]); fp=np.array([r['fp'] for r in rows])
    sa=np.array([r['sa'] for r in rows]); cov=np.array([r['cov'] for r in rows])
    am_cn=np.array([r['am_cn'] for r in rows]); por=np.array([r['por'] for r in rows])
    stcv=np.array([r['stcv'] for r in rows])

    thick=np.array([r['tau']<=1.5 for r in rows])
    mid=np.array([(r['tau']>1.5)&(r['tau']<=2.5) for r in rows])
    thin=np.array([r['tau']>2.5 for r in rows])
    names_arr = [r['name'] for r in rows]

    sigma_brug = SG * ps * fp / tau**2
    gp_gd2 = gp * gd**2

    # ═══════════════════════════════════════
    # PHASE 1: v4_brug 미세 조정
    # 현재 best: σ_brug × (1-0.18/φ)^0.5 × τ^1.5 / f × CN^1.5
    # 지수를 0.25 단위로 미세 조정
    # ═══════════════════════════════════════
    print("PHASE 1: v4_brug 미세 조정 (0.25 단위)")
    print("="*95)

    results = []
    for phi_c in [0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20]:
        R_perc_raw = np.clip(1-phi_c/ps, 0.001, None)
        for p in [0.25, 0.5, 0.75, 1.0]:
            R_perc = R_perc_raw**p
            for tb in [1.0, 1.25, 1.5, 1.75, 2.0]:  # τ^b in correction
                for fc in [0, 0.5, 1.0]:  # f_perc^c
                    for cne in [1.0, 1.25, 1.5, 1.75, 2.0]:  # CN
                        for ame in [0, 0.25, 0.5]:  # φ_AM (was 0 in champion)
                            rhs = sigma_brug * R_perc * tau**tb / fp**fc * cn**cne * pa**ame
                            C = fitC(sn, rhs)
                            if C is None or C<=0: continue
                            pred = C*rhs
                            r2a = r2l(sn, pred)
                            if r2a < 0.93: continue
                            r2t = r2l(sn[thick],pred[thick]) if thick.sum()>=2 else -99
                            r2m = r2l(sn[mid],pred[mid]) if mid.sum()>=2 else -99
                            r2n = r2l(sn[thin],pred[thin]) if thin.sum()>=2 else -99
                            avg3 = (r2t+r2m+r2n)/3 if min(r2t,r2m,r2n)>-50 else -99
                            results.append({
                                'name': f'σ_brug×(1-{phi_c}/φ)^{p}×τ^{tb}/f^{fc}×CN^{cne}×φ_AM^{ame}',
                                'r2':r2a,'C':C,'thick':r2t,'mid':r2m,'thin':r2n,'avg3':avg3,
                                'phi_c':phi_c,'p':p,'tb':tb,'fc':fc,'cne':cne,'ame':ame,'pred':pred,'rhs':rhs})

    results.sort(key=lambda x: -x['avg3'])

    print(f"\nTop 15 by avg(thick,mid,thin):")
    seen=set(); cnt=0
    for r in results:
        if r['name'] in seen: continue
        seen.add(r['name']); cnt+=1
        if cnt>15: break
        print(f"  #{cnt:2d} avg={r['avg3']:.3f} ALL={r['r2']:.4f} thick={r['thick']:.3f} mid={r['mid']:.3f} thin={r['thin']:.3f}")
        print(f"      {r['name']}")

    # ═══════════════════════════════════════
    # PHASE 2: champion 위에 추가 변수 1개씩 올리기
    # ═══════════════════════════════════════
    print(f"\n{'='*95}")
    print("PHASE 2: v4_brug champion + 추가 변수 1개")
    print("="*95)

    # Base: best from PHASE 1
    best = results[0] if results else None
    if best:
        base_rhs = best['rhs']
        base_r2 = best['r2']
        print(f"Base: {best['name']}, R²={base_r2:.4f}\n")

        extras = {
            'G_path': gp, 'GB_d': gd, '(G×d²)': gp_gd2,
            'hop_area': ha, 'bottleneck': bn, 'SE_area': sa,
            'coverage': cov, 'AM_CN': am_cn, 'T': T,
            'porosity': np.clip(por, 0.1, None), 'stress_CV': stcv,
            '1/G_path': 1/gp, '1/GB_d': 1/gd,
        }
        for ename, evar in extras.items():
            for exp in [-0.5, -0.25, 0.25, 0.5]:
                rhs2 = base_rhs * evar**exp
                C2 = fitC(sn, rhs2)
                if C2 is None or C2<=0: continue
                pred2 = C2*rhs2
                r2_2 = r2l(sn, pred2)
                delta = r2_2 - base_r2
                if abs(delta) > 0.001:
                    r2t2 = r2l(sn[thick],pred2[thick]) if thick.sum()>=2 else -99
                    r2n2 = r2l(sn[thin],pred2[thin]) if thin.sum()>=2 else -99
                    sign = '+' if delta>0 else '-'
                    print(f"  {sign} {ename}^{exp:+.2f}: ΔR²={delta:+.4f} ALL={r2_2:.4f} thick={r2t2:.3f} thin={r2n2:.3f}")

    # ═══════════════════════════════════════
    # PHASE 3: LOOCV + Per-case error
    # ═══════════════════════════════════════
    print(f"\n{'='*95}")
    print("PHASE 3: LOOCV + Per-case error analysis")
    print("="*95)

    top5 = []
    seen3=set()
    for r in results:
        if r['name'] in seen3: continue
        seen3.add(r['name'])
        top5.append(r)
        if len(top5)>=5: break

    for i, r in enumerate(top5):
        cv = loocv_C(sn, r['rhs'])
        pred = r['C'] * r['rhs']
        per_case_err = np.abs(sn - pred) / sn * 100
        worst_idx = np.argsort(per_case_err)[-5:]  # top 5 worst

        print(f"\n  #{i+1} LOOCV={cv:.4f} (train={r['r2']:.4f}, gap={r['r2']-cv:.4f})")
        print(f"      {r['name']}")
        print(f"      Worst 5 cases:")
        for wi in worst_idx[::-1]:
            regime = 'thick' if tau[wi]<=1.5 else ('mid' if tau[wi]<=2.5 else 'THIN')
            print(f"        {names_arr[wi]:30s} err={per_case_err[wi]:.1f}% τ={tau[wi]:.2f} φ={ps[wi]:.3f} [{regime}]")

    # ═══════════════════════════════════════
    # PHASE 4: 최종 후보 3개 head-to-head
    # ═══════════════════════════════════════
    print(f"\n{'='*95}")
    print("PHASE 4: FINAL HEAD-TO-HEAD")
    print("="*95)

    candidates = {}

    # v3
    rhs_v3 = sigma_brug * (gp_gd2)**0.25 * cn**2
    C_v3 = fitC(sn, rhs_v3); p_v3 = C_v3*rhs_v3
    cv_v3 = loocv_C(sn, rhs_v3)
    candidates['v3: σ_brug×(Gd²)^¼×CN²'] = (p_v3, rhs_v3, C_v3, cv_v3)

    # v4 direct
    phi_ex = np.clip(ps-0.18, 0.001, None)
    rhs_v4d = phi_ex * pa * cn**2
    C_v4d = fitC(sn, rhs_v4d); p_v4d = C_v4d*rhs_v4d
    cv_v4d = loocv_C(sn, rhs_v4d)
    candidates['v4_direct: (φ-0.18)×φ_AM×CN²'] = (p_v4d, rhs_v4d, C_v4d, cv_v4d)

    # v4_brug champion (from PHASE 1)
    if best:
        cv_best = loocv_C(sn, best['rhs'])
        candidates[f'v4_brug: {best["name"][:45]}'] = (best['pred'], best['rhs'], best['C'], cv_best)

    print(f"\n  {'Model':<55s} {'ALL':>6s} {'LOOCV':>6s} {'thick':>6s} {'mid':>6s} {'thin':>6s} {'C':>8s}")
    print("  "+"-"*100)
    for label, (pred, rhs, C, cv) in candidates.items():
        r2a = r2l(sn, pred)
        r2t = r2l(sn[thick],pred[thick]) if thick.sum()>=2 else -99
        r2m = r2l(sn[mid],pred[mid]) if mid.sum()>=2 else -99
        r2n = r2l(sn[thin],pred[thin]) if thin.sum()>=2 else -99
        print(f"  {label:<55s} {r2a:6.3f} {cv:6.3f} {r2t:6.3f} {r2m:6.3f} {r2n:6.3f} {C:8.4f}")

    # Effective formula decomposition for winner
    if best:
        print(f"\n  ★ v4_brug EFFECTIVE FORMULA:")
        eff_tau = best['tb'] - 2
        eff_fp = 1 - best['fc']
        print(f"    σ = {best['C']:.4f} × σ_grain × φ_SE^{1} × (1-{best['phi_c']}/φ)^{best['p']} × CN^{best['cne']} / τ^{-eff_tau}")
        if eff_fp != 1:
            print(f"    × f_perc^{eff_fp}")
        print(f"\n    = {best['C']:.4f} × {SG} × φ_SE × √(1-{best['phi_c']}/φ_SE) × CN^{best['cne']} / τ^{-eff_tau}")
        print(f"    φ_c = {best['phi_c']}, C = {best['C']:.4f}")
        print(f"\n    Physical interpretation:")
        print(f"    - σ_grain = 3.0 mS/cm (LPSCl grain interior)")
        print(f"    - φ_SE: SE volume fraction")
        print(f"    - (1-φ_c/φ)^{best['p']}: percolation proximity factor")
        print(f"    - CN^{best['cne']}: network redundancy (parallel paths)")
        print(f"    - τ^{eff_tau}: softened tortuosity (Bruggeman τ^-2 → τ^{eff_tau})")


if __name__ == '__main__':
    main()
