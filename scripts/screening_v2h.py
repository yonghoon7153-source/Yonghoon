"""
v2.0h FINAL: σ_brug 기반 통일 프레임워크 완성
==============================================
σ = σ_brug × C × R_perc × R_comp × CN² × R_contact

R_perc = (1 - φ_c/φ_SE)      ← percolation threshold correction
R_comp = φ_AM^a × τ^b / f_perc^c   ← τ² 상쇄 + 2상 경쟁
R_contact = (G_path × GB_d²)^d       ← contact quality (v3 legacy)

Questions:
1. R_contact가 R_perc + R_comp 위에 추가 R²를 주는가?
2. τ^b에서 b=2가 최적인가? (완전 상쇄 vs 부분 상쇄)
3. φ_c 최적값?
4. CN² vs CN^1 vs CN^1.5?
5. LOOCV 검증
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
            ps = m.get('phi_se', 0); pa = max(m.get('phi_am', 0), 0.01)
            t = m.get('tortuosity_recommended', m.get('tortuosity_mean', 0))
            fp = max(m.get('percolation_pct', 0)/100, 0.5)
            cn = m.get('se_se_cn', 0); gd = m.get('gb_density_mean', 0)
            gp = max(m.get('path_conductance_mean', 0), 1e-6)
            T = m.get('thickness_um', 0)
            if t<=0 or ps<=0 or cn<=0 or gd<=0 or T<=0: continue
            rows.append({'sn':sn,'ps':ps,'pa':pa,'tau':t,'fp':fp,'cn':cn,'gd':gd,'gp':gp,'T':T,'name':mp.parent.name})
    seen=set(); u=[]
    for r in rows:
        k=f"{r['ps']:.4f}_{r['T']:.1f}_{r['tau']:.3f}"
        if k not in seen: seen.add(k); u.append(r)
    return u

def r2l(a,p):
    la,lp=np.log(a),np.log(p); return 1-np.sum((la-lp)**2)/np.sum((la-np.mean(la))**2)
def fitC(a,r):
    v=(r>0)&np.isfinite(r); return float(np.exp(np.mean(np.log(a[v]/r[v])))) if v.sum()>=3 else None

def main():
    rows = load_data()
    n = len(rows); print(f"n={n}\n")
    sn=np.array([r['sn'] for r in rows])
    ps=np.array([r['ps'] for r in rows]); pa=np.array([r['pa'] for r in rows])
    tau=np.array([r['tau'] for r in rows]); cn=np.array([r['cn'] for r in rows])
    gd=np.array([r['gd'] for r in rows]); gp=np.array([r['gp'] for r in rows])
    fp=np.array([r['fp'] for r in rows]); T=np.array([r['T'] for r in rows])
    thick=np.array([r['tau']<=1.5 for r in rows])
    mid=np.array([(r['tau']>1.5)&(r['tau']<=2.5) for r in rows])
    thin=np.array([r['tau']>2.5 for r in rows])
    sigma_brug = SG * ps * fp / tau**2
    gp_gd2 = gp * gd**2

    print("=" * 95)
    print("σ = σ_brug × C × (1 - φ_c/φ_SE)^p × φ_AM^a × τ^b / f_perc^c × CN^d × (G×d²)^e")
    print("  where σ_brug already has φ_SE × f_perc / τ²")
    print("  effective: τ^(b-2), f_perc^(1-c)")
    print("=" * 95)

    results = []
    for phi_c in [0.10, 0.12, 0.15, 0.18, 0.20]:
        R_perc_raw = np.clip(1 - phi_c/ps, 0.001, None)
        for p in [0.5, 1, 1.5, 2]:  # exponent on percolation correction
            R_perc = R_perc_raw**p
            for a in [0, 0.5, 1, 1.5, 2, 3]:  # φ_AM
                for b in [0, 0.5, 1, 1.5, 2, 2.5]:  # τ^b (b=2 → full cancellation)
                    for c in [0, 0.5, 1]:  # f_perc^c in correction (c=1 → cancels f_perc)
                        for d in [0, 1, 1.5, 2]:  # CN
                            # Base: σ_brug × corrections
                            rhs = sigma_brug * R_perc * pa**a * tau**b / fp**c * cn**d
                            C = fitC(sn, rhs)
                            if C is None or C <= 0: continue
                            pred = C * rhs
                            r2a = r2l(sn, pred)
                            if r2a < 0.90: continue  # filter low R²

                            r2t = r2l(sn[thick],pred[thick]) if thick.sum()>=2 else -99
                            r2m = r2l(sn[mid],pred[mid]) if mid.sum()>=2 else -99
                            r2n = r2l(sn[thin],pred[thin]) if thin.sum()>=2 else -99

                            eff_tau = b - 2  # effective τ exponent
                            eff_fperc = 1 - c  # effective f_perc exponent

                            results.append({
                                'phi_c': phi_c, 'p': p, 'a': a, 'b': b, 'c': c, 'd': d,
                                'r2': r2a, 'C': C, 'thick': r2t, 'mid': r2m, 'thin': r2n,
                                'eff_tau': eff_tau, 'eff_fp': eff_fperc, 'pred': pred,
                                'name': f'σ_brug×(1-{phi_c}/φ)^{p}×φ_AM^{a}×τ^{b}/f^{c}×CN^{d}  [eff:τ^{eff_tau},f^{eff_fperc}]'
                            })

                            # + G_path contact correction
                            if r2a > 0.93:
                                for e in [-0.25, 0.25]:
                                    rhs2 = rhs * gp_gd2**e
                                    C2 = fitC(sn, rhs2)
                                    if C2 is None or C2 <= 0: continue
                                    pred2 = C2 * rhs2
                                    r2a2 = r2l(sn, pred2)
                                    if r2a2 > r2a + 0.002:
                                        r2t2 = r2l(sn[thick],pred2[thick]) if thick.sum()>=2 else -99
                                        r2m2 = r2l(sn[mid],pred2[mid]) if mid.sum()>=2 else -99
                                        r2n2 = r2l(sn[thin],pred2[thin]) if thin.sum()>=2 else -99
                                        results.append({
                                            'phi_c': phi_c, 'p': p, 'a': a, 'b': b, 'c': c, 'd': d,
                                            'r2': r2a2, 'C': C2, 'thick': r2t2, 'mid': r2m2, 'thin': r2n2,
                                            'eff_tau': eff_tau, 'eff_fp': eff_fperc, 'pred': pred2,
                                            'name': f'σ_brug×(1-{phi_c}/φ)^{p}×φ_AM^{a}×τ^{b}/f^{c}×CN^{d}×(Gd²)^{e}  [eff:τ^{eff_tau}]'
                                        })

    # ═══════════════════════════════════════
    # ANALYSIS
    # ═══════════════════════════════════════
    results.sort(key=lambda x: -x['r2'])

    # Sweet spot: ALL>0.93 AND thin>0.5
    sweet = [r for r in results if r['r2']>0.93 and r['thin']>0.5]
    sweet.sort(key=lambda x: -(x['thick']+x['mid']+x['thin'])/3)

    print(f"\n{'='*95}")
    print(f"SWEET SPOT: ALL>0.93 AND thin>0.5 ({len(sweet)} formulas)")
    print(f"{'='*95}")
    seen = set()
    for i, r in enumerate(sweet[:30]):
        if r['name'] in seen: continue
        seen.add(r['name'])
        avg = (r['thick']+r['mid']+r['thin'])/3
        print(f"\n  #{i+1:2d} ALL={r['r2']:.4f} avg={avg:.3f} thick={r['thick']:.3f} mid={r['mid']:.3f} thin={r['thin']:.3f}")
        print(f"      {r['name']}")
        print(f"      φ_c={r['phi_c']} C={r['C']:.4f}")

    # τ ANALYSIS: what effective τ exponent works best?
    print(f"\n{'='*95}")
    print(f"τ EXPONENT ANALYSIS (from sweet spot)")
    print(f"{'='*95}")
    tau_stats = {}
    for r in sweet:
        et = r['eff_tau']
        if et not in tau_stats or r['r2'] > tau_stats[et]['r2']:
            tau_stats[et] = r
    for et in sorted(tau_stats.keys()):
        r = tau_stats[et]
        print(f"  eff τ^{et:+.1f}: best ALL={r['r2']:.4f} thick={r['thick']:.3f} thin={r['thin']:.3f}  {r['name'][:60]}")

    # φ_c ANALYSIS
    print(f"\n{'='*95}")
    print(f"φ_c ANALYSIS (from sweet spot)")
    print(f"{'='*95}")
    phic_stats = {}
    for r in sweet:
        pc = r['phi_c']
        if pc not in phic_stats or r['r2'] > phic_stats[pc]['r2']:
            phic_stats[pc] = r
    for pc in sorted(phic_stats.keys()):
        r = phic_stats[pc]
        print(f"  φ_c={pc}: best ALL={r['r2']:.4f} thick={r['thick']:.3f} thin={r['thin']:.3f}")

    # LOOCV for top 5
    print(f"\n{'='*95}")
    print(f"LOOCV (top 5 sweet)")
    print(f"{'='*95}")
    seen2 = set(); cnt = 0
    for r in sweet:
        if r['name'] in seen2: continue
        seen2.add(r['name']); cnt += 1
        if cnt > 5: break
        log_actual = np.log(sn)
        log_rhs = np.log(r['pred'] / r['C'])
        loo_errors = []
        for j in range(n):
            m = np.ones(n, bool); m[j] = False
            C_loo = float(np.exp(np.mean(log_actual[m] - log_rhs[m])))
            loo_errors.append((log_actual[j] - np.log(C_loo * r['pred'][j] / r['C']))**2)
        ss_res = np.sum(loo_errors)
        ss_tot = np.sum((log_actual - np.mean(log_actual))**2)
        r2_cv = 1 - ss_res/ss_tot
        print(f"  LOOCV={r2_cv:.4f} (train={r['r2']:.4f}, gap={r['r2']-r2_cv:.4f})")
        print(f"    {r['name']}")

    # FINAL COMPARISON
    print(f"\n{'='*95}")
    print(f"FINAL COMPARISON TABLE")
    print(f"{'='*95}")
    v3_rhs = sigma_brug * (gp_gd2)**0.25 * cn**2
    C3 = fitC(sn, v3_rhs); p3 = C3*v3_rhs

    v4p_rhs = np.clip(ps-0.18, 0.001, None) * pa * cn**2
    C4p = fitC(sn, v4p_rhs); p4p = C4p*v4p_rhs

    print(f"\n  {'Model':<60s} {'ALL':>6s} {'thick':>6s} {'mid':>6s} {'thin':>6s}")
    print("  "+"-"*90)
    print(f"  {'v3: σ_brug×C×(Gd²)^¼×CN²':<60s} {r2l(sn,p3):6.3f} {r2l(sn[thick],p3[thick]):6.3f} {r2l(sn[mid],p3[mid]):6.3f} {r2l(sn[thin],p3[thin]):6.3f}")
    print(f"  {'v4: C×σ_grain×(φ-0.18)×φ_AM×CN²':<60s} {r2l(sn,p4p):6.3f} {r2l(sn[thick],p4p[thick]):6.3f} {r2l(sn[mid],p4p[mid]):6.3f} {r2l(sn[thin],p4p[thin]):6.3f}")
    if sweet:
        bs = sweet[0]
        print(f"  {'v4_brug: '+bs['name'][:50]:<60s} {bs['r2']:6.3f} {bs['thick']:6.3f} {bs['mid']:6.3f} {bs['thin']:6.3f}")

    # Express v4 in σ_brug form
    print(f"\n  v4를 σ_brug 형태로 표현:")
    print(f"  σ = σ_brug × C × (1 - φ_c/φ_SE) × (φ_AM × τ² / f_perc) × CN²")
    print(f"  = σ_brug × C × R_perc × R_comp × CN²")
    print(f"  where:")
    print(f"    R_perc = (1 - 0.18/φ_SE)     ← percolation threshold")
    print(f"    R_comp = φ_AM × τ² / f_perc  ← τ² cancel + 2-phase")
    print(f"    CN²                           ← network connectivity (v3 legacy)")
    print(f"    C = {C4p:.4f} (in direct form) = C_brug × ...")
    print(f"\n  ★ τ²가 σ_brug의 τ^-2와 상쇄되어 최종식에서 사라짐")
    print(f"  ★ f_perc가 σ_brug의 f_perc와 상쇄")
    print(f"  ★ 남은 것: σ_grain × (φ_SE - φ_c) × φ_AM × CN²")


if __name__ == '__main__':
    main()
