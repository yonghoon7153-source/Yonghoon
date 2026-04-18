"""
v2.0j BEYOND CHAMPION: coverage + outlier 해결 + 극한
======================================================
v4 champion 위에:
1. coverage 추가 (유일한 ΔR²>0 변수)
2. outlier 2개 (thin_8 variants) 원인 분석
3. φ_c 정밀 탐색 (0.01 단위)
4. 지수 0.125 단위 극한 미세조정
5. 2-term correction: √(1-φ_c/φ) × coverage^a 결합
6. 전체 parameter 1개씩 ablation study
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
            T=m.get('thickness_um',0)
            cov=max(m.get('coverage_AM_P_mean',m.get('coverage_AM_S_mean',m.get('coverage_AM_mean',1))),0.1)
            stcv=max(m.get('stress_cv',0.1),0.1)
            por=m.get('porosity',0)
            if t<=0 or ps<=0 or cn<=0 or gd<=0 or T<=0: continue
            rows.append({'sn':sn,'ps':ps,'pa':pa,'tau':t,'fp':fp,'cn':cn,'gd':gd,
                'gp':gp,'T':T,'cov':cov,'stcv':stcv,'por':por,'name':mp.parent.name})
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
    T=np.array([r['T'] for r in rows]); fp=np.array([r['fp'] for r in rows])
    cov=np.array([r['cov'] for r in rows])
    stcv=np.array([r['stcv'] for r in rows])
    por=np.array([r['por'] for r in rows])
    names_arr=[r['name'] for r in rows]

    thick=np.array([r['tau']<=1.5 for r in rows])
    mid=np.array([(r['tau']>1.5)&(r['tau']<=2.5) for r in rows])
    thin=np.array([r['tau']>2.5 for r in rows])
    sigma_brug = SG * ps * fp / tau**2

    results = []

    # ═══════════════════════════════════════
    # PART 1: v4 champion + coverage
    # ═══════════════════════════════════════
    print("PART 1: v4 champion + coverage^a")
    print("="*95)

    for phi_c in np.arange(0.14, 0.21, 0.01):
        R_perc = np.clip(1-phi_c/ps, 0.001, None)**0.5
        for tb in [1.25, 1.375, 1.5]:
            for cne in [1.25, 1.375, 1.5]:
                for cova in [0, 0.125, 0.25, 0.375, 0.5]:
                    rhs = sigma_brug * R_perc * tau**tb / fp * cn**cne * cov**(cova/100)  # coverage in %
                    # coverage is in %, need to handle
                    rhs2 = sigma_brug * R_perc * tau**tb / fp * cn**cne * (cov/100)**cova
                    C = fitC(sn, rhs2)
                    if C is None or C<=0: continue
                    pred = C*rhs2
                    r2a = r2l(sn, pred)
                    if r2a < 0.93: continue
                    r2t = r2l(sn[thick],pred[thick]) if thick.sum()>=2 else -99
                    r2m = r2l(sn[mid],pred[mid]) if mid.sum()>=2 else -99
                    r2n = r2l(sn[thin],pred[thin]) if thin.sum()>=2 else -99
                    avg3 = (r2t+r2m+r2n)/3
                    results.append({
                        'name': f'v4+cov: φ_c={phi_c:.2f} τ^{tb} CN^{cne} cov^{cova}',
                        'r2':r2a,'C':C,'thick':r2t,'mid':r2m,'thin':r2n,'avg3':avg3,
                        'phi_c':phi_c,'tb':tb,'cne':cne,'cova':cova,'pred':pred,'rhs':rhs2})

    results.sort(key=lambda x: -x['avg3'])
    print(f"\nTop 20 by avg(thick,mid,thin):")
    seen=set(); cnt=0
    for r in results:
        if r['name'] in seen: continue
        seen.add(r['name']); cnt+=1
        if cnt>20: break
        print(f"  #{cnt:2d} avg={r['avg3']:.3f} ALL={r['r2']:.4f} thick={r['thick']:.3f} mid={r['mid']:.3f} thin={r['thin']:.3f} cov^{r['cova']}")
        print(f"      {r['name']}")

    # ═══════════════════════════════════════
    # PART 2: φ_c 정밀 탐색 (0.005 단위)
    # ═══════════════════════════════════════
    print(f"\n{'='*95}")
    print("PART 2: φ_c 정밀 탐색")
    print("="*95)

    phic_results = []
    for phi_c in np.arange(0.10, 0.25, 0.005):
        R_perc = np.clip(1-phi_c/ps, 0.001, None)**0.5
        # Use champion exponents: τ^1.5, CN^1.5
        rhs = sigma_brug * R_perc * tau**1.5 / fp * cn**1.5
        C = fitC(sn, rhs)
        if C is None or C<=0: continue
        pred = C*rhs
        r2a = r2l(sn, pred)
        r2t = r2l(sn[thick],pred[thick]) if thick.sum()>=2 else -99
        r2n = r2l(sn[thin],pred[thin]) if thin.sum()>=2 else -99
        phic_results.append({'phi_c':phi_c,'r2':r2a,'thick':r2t,'thin':r2n})

    print(f"  {'φ_c':>6s} {'ALL':>7s} {'thick':>7s} {'thin':>7s}")
    for r in phic_results:
        flag = '★' if r['r2']>0.94 and r['thin']>0.85 else ' '
        print(f"  {flag}{r['phi_c']:5.3f} {r['r2']:7.4f} {r['thick']:7.3f} {r['thin']:7.3f}")

    # ═══════════════════════════════════════
    # PART 3: Ablation study (각 항 제거 시 R² 변화)
    # ═══════════════════════════════════════
    print(f"\n{'='*95}")
    print("PART 3: Ablation Study — 각 항 제거 시 ΔR²")
    print("="*95)

    # Full champion
    R_perc_full = np.clip(1-0.18/ps, 0.001, None)**0.5
    rhs_full = sigma_brug * R_perc_full * tau**1.5 / fp * cn**1.5
    C_full = fitC(sn, rhs_full); pred_full = C_full*rhs_full
    r2_full = r2l(sn, pred_full)
    print(f"  Full: σ_brug × √(1-0.18/φ) × τ^1.5/f × CN^1.5 → R²={r2_full:.4f}")

    ablations = {
        'Remove √(1-φ_c/φ)': sigma_brug * tau**1.5 / fp * cn**1.5,
        'Remove τ^1.5 (keep σ_brug τ^-2)': sigma_brug * R_perc_full / fp * cn**1.5,
        'Remove /f_perc': sigma_brug * R_perc_full * tau**1.5 * cn**1.5,
        'Remove CN^1.5': sigma_brug * R_perc_full * tau**1.5 / fp,
        'Remove σ_brug entirely': R_perc_full * tau**1.5 / fp * cn**1.5,
        'Only σ_brug (no correction)': sigma_brug,
        'Only (φ-0.18)×φ_AM×CN² (v4_direct)': np.clip(ps-0.18,0.001,None)*pa*cn**2,
    }

    for label, rhs_abl in ablations.items():
        C_abl = fitC(sn, rhs_abl)
        if C_abl is None or C_abl<=0:
            print(f"  {label}: FAILED")
            continue
        pred_abl = C_abl*rhs_abl
        r2_abl = r2l(sn, pred_abl)
        delta = r2_abl - r2_full
        r2t = r2l(sn[thick],pred_abl[thick]) if thick.sum()>=2 else -99
        r2n = r2l(sn[thin],pred_abl[thin]) if thin.sum()>=2 else -99
        print(f"  {label:40s}: R²={r2_abl:.4f} (Δ={delta:+.4f}) thick={r2t:.3f} thin={r2n:.3f}")

    # ═══════════════════════════════════════
    # PART 4: LOOCV + outlier deep dive
    # ═══════════════════════════════════════
    print(f"\n{'='*95}")
    print("PART 4: Per-case error (전체)")
    print("="*95)

    pred_champ = C_full * rhs_full
    per_err = np.abs(sn - pred_champ) / sn * 100
    cv = loocv_C(sn, rhs_full)
    print(f"  LOOCV = {cv:.4f} (train={r2_full:.4f})")
    print(f"  Mean |err| = {np.mean(per_err):.1f}%, Median = {np.median(per_err):.1f}%")
    print(f"  Within 20%: {np.sum(per_err<20)}/{n} ({np.sum(per_err<20)/n*100:.0f}%)")
    print(f"  Within 30%: {np.sum(per_err<30)}/{n} ({np.sum(per_err<30)/n*100:.0f}%)")
    print(f"  Within 50%: {np.sum(per_err<50)}/{n} ({np.sum(per_err<50)/n*100:.0f}%)")

    print(f"\n  All cases sorted by error:")
    print(f"  {'Case':>35s} {'σ_act':>8s} {'σ_pred':>8s} {'err%':>7s} {'τ':>5s} {'φ_SE':>6s} {'CN':>5s} {'regime'}")
    sorted_idx = np.argsort(per_err)[::-1]
    for idx in sorted_idx:
        regime = 'thick' if tau[idx]<=1.5 else ('mid' if tau[idx]<=2.5 else 'THIN')
        flag = '>>>' if per_err[idx]>50 else '> ' if per_err[idx]>30 else '  '
        print(f"  {flag}{names_arr[idx]:33s} {sn[idx]:8.4f} {pred_champ[idx]:8.4f} {per_err[idx]:6.1f}% {tau[idx]:5.2f} {ps[idx]:6.3f} {cn[idx]:5.1f} {regime}")

    # ═══════════════════════════════════════
    # PART 5: σ_brug 형태로 최종 정리
    # ═══════════════════════════════════════
    print(f"\n{'='*95}")
    print("PART 5: 최종 공식 정리")
    print("="*95)

    print(f"""
  ┌─────────────────────────────────────────────────────────┐
  │  v4 CHAMPION (Universal Ionic Conductivity Formula)     │
  │                                                         │
  │  σ_ion = σ_brug × C × √(1 - φ_c/φ_SE) × τ^(3/2)     │
  │          / f_perc × CN^(3/2)                            │
  │                                                         │
  │  Effective form:                                        │
  │  σ_ion = C × σ_grain × φ_SE × √(1-φ_c/φ_SE)          │
  │          × CN^(3/2) / √τ                               │
  │                                                         │
  │  C = {C_full:.4f}, φ_c = 0.18                          │
  │  σ_grain = 3.0 mS/cm (LPSCl grain interior)            │
  │                                                         │
  │  R² = {r2_full:.4f} (ALL, n={n})                        │
  │  LOOCV = {cv:.4f}                                       │
  │  thick = {r2l(sn[thick],pred_champ[thick]):.3f}         │
  │  thin  = {r2l(sn[thin],pred_champ[thin]):.3f}           │
  │  1 free parameter (C only)                              │
  └─────────────────────────────────────────────────────────┘

  Physical decomposition (from σ_brug):
  ┌────────────────────────────────────────┐
  │ σ_brug = σ_grain × φ × f_perc / τ²    │  ← Bruggeman base
  │ × √(1-φ_c/φ)                          │  ← Percolation threshold
  │ × τ^(3/2) / f_perc                    │  ← τ² over-penalty cancel
  │ × CN^(3/2)                            │  ← Network connectivity
  │ × C                                    │  ← Constriction geometry
  └────────────────────────────────────────┘

  After cancellation:
  - τ^(3/2) × τ^(-2) = τ^(-1/2)  ← softened from τ^(-2)
  - f_perc^(-1) × f_perc = f_perc^0  ← cancelled
  - φ × (1-φ_c/φ)^(1/2) = φ × √(φ-φ_c)/√φ = √φ × √(φ-φ_c)
    ← φ^(1/2) × (φ-φ_c)^(1/2)

  So: σ = C × σ_grain × √φ_SE × √(φ_SE-φ_c) × CN^(3/2) / √τ
  """)


if __name__ == '__main__':
    main()
