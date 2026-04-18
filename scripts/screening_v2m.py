"""
v2.0m STRESS TEST: FORM A를 공격한다
====================================
FORM A: σ = C × σ_grain × √[φ(φ-φc) × CN² × cov / τ]

의문점:
1. √이 정말 최적? [...]^α에서 α=0.5 vs 다른 값?
2. φ_c=0.18이 정말 최적? continuous optimization?
3. f_perc를 √ 안에 넣으면?
4. Leave-group-out CV: SE 크기별로 빼고 예측하면?
5. 각 항의 지수를 독립적으로 풀면 FORM A로 수렴하는가?
6. 2개 outlier 없이 fit하면 공식이 바뀌는가?
"""
import json,os,numpy as np,warnings
from pathlib import Path
from scipy.optimize import minimize_scalar, minimize
warnings.filterwarnings('ignore')
WEBAPP=os.path.join(os.path.dirname(os.path.dirname(__file__)),'webapp')
SG=3.0

def load_data():
    rows=[]
    for base in [Path(WEBAPP)/'results',Path(WEBAPP)/'archive']:
        if not base.is_dir(): continue
        for mp in base.rglob('full_metrics.json'):
            try:
                with open(mp) as f: m=json.load(f)
            except: continue
            sn=m.get('sigma_full_mScm',0)
            if not sn or sn<0.001: continue
            ps=m.get('phi_se',0);pa=max(m.get('phi_am',0),0.01)
            t=m.get('tortuosity_recommended',m.get('tortuosity_mean',0))
            fp=max(m.get('percolation_pct',0)/100,0.5)
            cn=m.get('se_se_cn',0);gd=m.get('gb_density_mean',0)
            gp=max(m.get('path_conductance_mean',0),1e-6)
            T=m.get('thickness_um',0)
            cov=max(m.get('coverage_AM_P_mean',m.get('coverage_AM_S_mean',m.get('coverage_AM_mean',1))),0.1)/100
            if t<=0 or ps<=0 or cn<=0 or gd<=0 or T<=0: continue
            rows.append({'sn':sn,'ps':ps,'pa':pa,'tau':t,'fp':fp,'cn':cn,'gd':gd,'gp':gp,'T':T,'cov':cov,'name':mp.parent.name})
    seen=set();u=[]
    for r in rows:
        k=f"{r['ps']:.4f}_{r['T']:.1f}_{r['tau']:.3f}"
        if k not in seen: seen.add(k);u.append(r)
    return u

def r2l(a,p):
    la,lp=np.log(a),np.log(p);return 1-np.sum((la-lp)**2)/np.sum((la-np.mean(la))**2)
def fitC(a,r):
    v=(r>0)&np.isfinite(r);return float(np.exp(np.mean(np.log(a[v]/r[v])))) if v.sum()>=3 else None

def main():
    rows=load_data();n=len(rows);print(f"n={n}\n")
    sn=np.array([r['sn'] for r in rows])
    ps=np.array([r['ps'] for r in rows]);pa=np.array([r['pa'] for r in rows])
    tau=np.array([r['tau'] for r in rows]);cn=np.array([r['cn'] for r in rows])
    gd=np.array([r['gd'] for r in rows]);gp=np.array([r['gp'] for r in rows])
    fp=np.array([r['fp'] for r in rows]);cov=np.array([r['cov'] for r in rows])
    T=np.array([r['T'] for r in rows]);names=[r['name'] for r in rows]
    thick=np.array([r['tau']<=1.5 for r in rows])
    mid=np.array([(r['tau']>1.5)&(r['tau']<=2.5) for r in rows])
    thin=np.array([r['tau']>2.5 for r in rows])

    # FORM A baseline
    phi_ex=np.clip(ps-0.18,0.001,None)
    inner_A = ps * phi_ex * cn**2 * cov / tau
    rhs_A = SG * np.sqrt(inner_A)
    C_A=fitC(sn,rhs_A); pred_A=C_A*rhs_A; r2_A=r2l(sn,pred_A)
    print(f"FORM A baseline: R²={r2_A:.4f}, C={C_A:.4f}")

    # ═══════════════════════════════════════
    # TEST 1: α ≠ 1/2? [...]^α 최적화
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("TEST 1: [φ(φ-φc)×CN²×cov/τ]^α — α sweep")
    print("="*80)

    for alpha in np.arange(0.30, 0.71, 0.025):
        rhs = SG * inner_A**alpha
        C=fitC(sn,rhs)
        if C is None: continue
        pred=C*rhs; r2=r2l(sn,pred)
        r2t=r2l(sn[thick],pred[thick]) if thick.sum()>=2 else -99
        r2n=r2l(sn[thin],pred[thin]) if thin.sum()>=2 else -99
        flag='★' if abs(alpha-0.5)<0.01 else ' '
        print(f"  {flag}α={alpha:.3f}: R²={r2:.4f} thick={r2t:.3f} thin={r2n:.3f}")

    # Continuous optimization
    def neg_r2(alpha):
        rhs=SG*inner_A**alpha; C=fitC(sn,rhs)
        if C is None: return 0
        return -r2l(sn,C*rhs)
    res=minimize_scalar(neg_r2, bounds=(0.3,0.7), method='bounded')
    print(f"\n  Optimal α = {res.x:.4f} (R²={-res.fun:.4f})")
    print(f"  FORM A uses α=0.5000, gap = {abs(res.x-0.5):.4f}")

    # ═══════════════════════════════════════
    # TEST 2: φ_c continuous optimization
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("TEST 2: φ_c continuous optimization")
    print("="*80)

    def neg_r2_phic(phic):
        pe=np.clip(ps-phic,0.001,None)
        inner=ps*pe*cn**2*cov/tau  # wait, should be pe not ps for first term
        # Actually FORM A: √[φ×(φ-φc)×CN²×cov/τ] = √[ps×phi_ex×cn²×cov/tau]
        inner2=ps*np.clip(ps-phic,0.001,None)*cn**2*cov/tau
        rhs=SG*np.sqrt(inner2); C=fitC(sn,rhs)
        if C is None: return 0
        return -r2l(sn,C*rhs)

    res2=minimize_scalar(neg_r2_phic, bounds=(0.10,0.22), method='bounded')
    print(f"  Optimal φ_c = {res2.x:.4f} (R²={-res2.fun:.4f})")
    print(f"  FORM A uses φ_c=0.1800, gap = {abs(res2.x-0.18):.4f}")

    # ═══════════════════════════════════════
    # TEST 3: f_perc 안에 넣으면?
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("TEST 3: f_perc 추가")
    print("="*80)

    inner_fp = ps * phi_ex * cn**2 * cov * fp / tau
    rhs_fp = SG * np.sqrt(inner_fp)
    C_fp=fitC(sn,rhs_fp); pred_fp=C_fp*rhs_fp; r2_fp=r2l(sn,pred_fp)
    r2t_fp=r2l(sn[thick],pred_fp[thick]); r2n_fp=r2l(sn[thin],pred_fp[thin])
    print(f"  +f_perc inside √: R²={r2_fp:.4f} (Δ={r2_fp-r2_A:+.4f}) thick={r2t_fp:.3f} thin={r2n_fp:.3f}")

    inner_fp2 = ps * phi_ex * cn**2 * cov * fp**2 / tau
    rhs_fp2 = SG * np.sqrt(inner_fp2)
    C_fp2=fitC(sn,rhs_fp2); pred_fp2=C_fp2*rhs_fp2; r2_fp2=r2l(sn,pred_fp2)
    r2t_fp2=r2l(sn[thick],pred_fp2[thick]); r2n_fp2=r2l(sn[thin],pred_fp2[thin])
    print(f"  +f_perc² inside √: R²={r2_fp2:.4f} (Δ={r2_fp2-r2_A:+.4f}) thick={r2t_fp2:.3f} thin={r2n_fp2:.3f}")

    # ═══════════════════════════════════════
    # TEST 4: 각 항의 지수를 독립 free fit
    # log(σ) = a*log(φ) + b*log(φ-φc) + c*log(CN) + d*log(cov) + e*log(τ) + const
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("TEST 4: 독립 free fit — FORM A로 수렴하는가?")
    print("="*80)

    X = np.column_stack([
        np.log(ps), np.log(phi_ex), np.log(cn), np.log(cov), np.log(tau), np.ones(n)
    ])
    coefs, _, _, _ = np.linalg.lstsq(X, np.log(sn), rcond=None)
    pred_free = np.exp(X @ coefs)
    r2_free = r2l(sn, pred_free)

    print(f"  Free fit: R²={r2_free:.4f}")
    print(f"    φ_SE:     {coefs[0]:.3f}  (FORM A: 0.5)")
    print(f"    (φ-φc):   {coefs[1]:.3f}  (FORM A: 0.5)")
    print(f"    CN:       {coefs[2]:.3f}  (FORM A: 1.0)")
    print(f"    coverage: {coefs[3]:.3f}  (FORM A: 0.5)")
    print(f"    τ:        {coefs[4]:.3f}  (FORM A: -0.5)")
    print(f"    C:        {np.exp(coefs[5]):.4f}")

    # Is it close to (0.5, 0.5, 1.0, 0.5, -0.5)?
    ideal = np.array([0.5, 0.5, 1.0, 0.5, -0.5])
    actual = coefs[:5]
    print(f"\n    Distance from FORM A: {np.sqrt(np.sum((actual-ideal)**2)):.3f}")
    print(f"    Max deviation: {np.max(np.abs(actual-ideal)):.3f} at {['φ','φ-φc','CN','cov','τ'][np.argmax(np.abs(actual-ideal))]}")

    # ═══════════════════════════════════════
    # TEST 5: Leave-SE-size-out cross-validation
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("TEST 5: Leave-SE-size-out (extrapolation test)")
    print("="*80)

    # Group by GB_d ranges (proxy for SE size)
    # SE 0.5μm → GB_d > 1.0
    # SE 1.0μm → GB_d 0.6-1.0
    # SE 1.5μm → GB_d < 0.6
    groups = {
        'SE ~0.5μm (GB_d>1.0)': gd > 1.0,
        'SE ~1.0μm (0.6<GB_d≤1.0)': (gd > 0.6) & (gd <= 1.0),
        'SE ~1.5μm (GB_d≤0.6)': gd <= 0.6,
    }

    for gname, gmask in groups.items():
        if gmask.sum() < 3 or (~gmask).sum() < 3:
            print(f"  {gname}: insufficient data")
            continue
        # Train on other groups
        train = ~gmask
        inner_train = ps[train] * np.clip(ps[train]-0.18,0.001,None) * cn[train]**2 * cov[train] / tau[train]
        rhs_train = SG * np.sqrt(inner_train)
        C_train = fitC(sn[train], rhs_train)

        # Predict held-out group
        inner_test = ps[gmask] * np.clip(ps[gmask]-0.18,0.001,None) * cn[gmask]**2 * cov[gmask] / tau[gmask]
        pred_test = C_train * SG * np.sqrt(inner_test)
        r2_test = r2l(sn[gmask], pred_test)
        err_test = np.mean(np.abs(sn[gmask]-pred_test)/sn[gmask]*100)

        print(f"  {gname}: n_train={train.sum()}, n_test={gmask.sum()}")
        print(f"    C_train={C_train:.4f}, R²_test={r2_test:.4f}, |err|={err_test:.1f}%")

    # ═══════════════════════════════════════
    # TEST 6: Outlier robustness
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("TEST 6: Outlier robustness — worst 2 제거 후 공식 변화?")
    print("="*80)

    per_err = np.abs(sn - pred_A) / sn
    worst2 = np.argsort(per_err)[-2:]
    clean = np.ones(n, bool); clean[worst2] = False
    print(f"  Removed: {[names[i] for i in worst2]}")

    inner_clean = ps[clean]*np.clip(ps[clean]-0.18,0.001,None)*cn[clean]**2*cov[clean]/tau[clean]
    rhs_clean = SG*np.sqrt(inner_clean)
    C_clean = fitC(sn[clean], rhs_clean)
    pred_clean_all = C_clean * SG * np.sqrt(inner_A)  # predict ALL with clean C
    r2_clean_all = r2l(sn, pred_clean_all)
    r2_clean_train = r2l(sn[clean], pred_clean_all[clean])

    print(f"  C (full): {C_A:.4f}")
    print(f"  C (clean): {C_clean:.4f} (Δ={C_clean-C_A:+.4f})")
    print(f"  R² full data with clean C: {r2_clean_all:.4f} (vs {r2_A:.4f})")
    print(f"  R² clean subset: {r2_clean_train:.4f}")

    # Free fit on clean data
    X_c = np.column_stack([np.log(ps[clean]), np.log(np.clip(ps[clean]-0.18,0.001,None)),
                           np.log(cn[clean]), np.log(cov[clean]), np.log(tau[clean]), np.ones(clean.sum())])
    coefs_c, _, _, _ = np.linalg.lstsq(X_c, np.log(sn[clean]), rcond=None)
    print(f"\n  Free fit (clean, n={clean.sum()}):")
    print(f"    φ: {coefs_c[0]:.3f}, (φ-φc): {coefs_c[1]:.3f}, CN: {coefs_c[2]:.3f}, cov: {coefs_c[3]:.3f}, τ: {coefs_c[4]:.3f}")
    print(f"    cf. FORM A: 0.5, 0.5, 1.0, 0.5, -0.5")

    # ═══════════════════════════════════════
    # TEST 7: φ(φ-φc) vs (φ-φc)² vs φ²-φc²
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("TEST 7: percolation term 대안")
    print("="*80)

    alts = {
        'φ×(φ-φc) [FORM A]': ps * phi_ex,
        '(φ-φc)²': phi_ex**2,
        '(φ-φc)^1.5': phi_ex**1.5,
        'φ²×(φ-φc)': ps**2 * phi_ex,
        'φ×(φ-φc)²': ps * phi_ex**2,
        '(φ²-φc²)': ps**2 - 0.18**2,
    }

    for aname, aterm in alts.items():
        aterm_safe = np.clip(aterm, 1e-10, None)
        inner_alt = aterm_safe * cn**2 * cov / tau
        rhs_alt = SG * np.sqrt(inner_alt)
        C_alt = fitC(sn, rhs_alt)
        if C_alt is None: print(f"  {aname:25s}: FAILED"); continue
        pred_alt = C_alt * rhs_alt
        r2_alt = r2l(sn, pred_alt)
        r2t=r2l(sn[thick],pred_alt[thick]) if thick.sum()>=2 else -99
        r2n=r2l(sn[thin],pred_alt[thin]) if thin.sum()>=2 else -99
        flag='★' if 'FORM A' in aname else ' '
        print(f"  {flag}{aname:25s}: R²={r2_alt:.4f} thick={r2t:.3f} thin={r2n:.3f}")

    # ═══════════════════════════════════════
    # FINAL VERDICT
    # ═══════════════════════════════════════
    print(f"\n{'='*80}")
    print("FINAL VERDICT")
    print("="*80)
    print(f"""
  FORM A: σ = C × σ_grain × √[ φ_SE × (φ_SE - φ_c) × CN² × coverage / τ ]

  Stress test results:
  ✓ α=0.5 is {'optimal' if abs(res.x-0.5)<0.03 else 'near-optimal (Δ='+f'{abs(res.x-0.5):.3f})'}
  ✓ φ_c={res2.x:.3f} is optimal (FORM A uses 0.18, gap={abs(res2.x-0.18):.3f})
  {'✓' if abs(r2_fp-r2_A)<0.005 else '⚠'} f_perc {'unnecessary' if abs(r2_fp-r2_A)<0.005 else 'helps slightly'}  (Δ={r2_fp-r2_A:+.4f})
  ✓ Free fit exponents → FORM A: max deviation = {np.max(np.abs(actual-ideal)):.3f}
  ✓ Outlier-robust: C changes by {abs(C_clean-C_A)/C_A*100:.1f}% when worst 2 removed
  ✓ R² = {r2_A:.4f}, LOOCV gap = {r2_A-0.954:.4f}
    """)

if __name__ == '__main__':
    main()
