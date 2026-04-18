"""
v2.0f BALANCED SEARCH: thick + thin 동시 최적화
================================================
최적화 기준: √(R²_thick × R²_thin) — 둘 다 좋아야 높음
"""
import json, os, numpy as np, warnings
from pathlib import Path
from itertools import product

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
            ps, pa = m.get('phi_se', 0), max(m.get('phi_am', 0), 0.01)
            t = m.get('tortuosity_recommended', m.get('tortuosity_mean', 0))
            fp = max(m.get('percolation_pct', 0)/100, 0.5)
            cn = m.get('se_se_cn', 0)
            gd = m.get('gb_density_mean', 0)
            gp = max(m.get('path_conductance_mean', 0), 1e-6)
            ha = max(m.get('path_hop_area_mean', 0), 1e-6)
            T = m.get('thickness_um', 0)
            sa = max(m.get('area_SE_SE_mean', 1e-6), 1e-6)
            bn = max(m.get('path_hop_area_min_mean', 1e-6), 1e-6)
            if t<=0 or ps<=0 or cn<=0 or gd<=0 or T<=0: continue
            rows.append({'name': mp.parent.name, 'sn': sn,
                'ps': ps, 'pa': pa, 'tau': t, 'fp': fp, 'cn': cn,
                'gd': gd, 'gp': gp, 'ha': ha, 'T': T, 'sa': sa, 'bn': bn})
    seen = set()
    u = []
    for r in rows:
        k = f"{r['ps']:.4f}_{r['T']:.1f}_{r['tau']:.3f}"
        if k not in seen: seen.add(k); u.append(r)
    return u


def r2l(a, p):
    la, lp = np.log(a), np.log(p)
    return 1 - np.sum((la-lp)**2) / np.sum((la-np.mean(la))**2)


def fitC(a, r):
    v = (r>0)&np.isfinite(r)
    return float(np.exp(np.mean(np.log(a[v]/r[v])))) if v.sum()>=3 else None


def main():
    rows = load_data()
    n = len(rows)
    print(f"n={n}\n")

    sn = np.array([r['sn'] for r in rows])
    ps = np.array([r['ps'] for r in rows])
    pa = np.array([r['pa'] for r in rows])
    tau = np.array([r['tau'] for r in rows])
    cn = np.array([r['cn'] for r in rows])
    gd = np.array([r['gd'] for r in rows])
    gp = np.array([r['gp'] for r in rows])
    ha = np.array([r['ha'] for r in rows])
    T = np.array([r['T'] for r in rows])
    fp = np.array([r['fp'] for r in rows])
    sa = np.array([r['sa'] for r in rows])
    bn = np.array([r['bn'] for r in rows])

    thick = np.array([r['tau'] <= 1.5 for r in rows])
    mid = np.array([(r['tau'] > 1.5) & (r['tau'] <= 2.5) for r in rows])
    thin = np.array([r['tau'] > 2.5 for r in rows])

    results = []

    # Derived variables
    ps_pa = ps * pa           # φ_SE × φ_AM
    cn_tau = cn / tau          # CN/τ
    ps_tau = ps / tau          # φ_SE/τ
    ps_f = ps * fp             # φ_SE × f_perc
    gp_gd2 = gp * gd**2       # G_path × GB_d²

    def test(name, rhs, nfree=1):
        C = fitC(sn, rhs)
        if C is None or C<=0: return
        pred = C * rhs
        r2_all = r2l(sn, pred)
        r2_t = r2l(sn[thick], pred[thick]) if thick.sum()>=2 else -99
        r2_m = r2l(sn[mid], pred[mid]) if mid.sum()>=2 else -99
        r2_n = r2l(sn[thin], pred[thin]) if thin.sum()>=2 else -99
        # Balanced score: geometric mean of thick, mid, thin
        scores = [r2_t, r2_m, r2_n]
        valid_scores = [s for s in scores if s > 0]
        balanced = np.prod(valid_scores)**(1/len(valid_scores)) if valid_scores else 0
        err_all = np.mean(np.abs(sn-pred)/sn)*100
        results.append({'name': name, 'r2': r2_all, 'C': C, 'nf': nfree,
            'thick': r2_t, 'mid': r2_m, 'thin': r2_n, 'bal': balanced, 'err': err_all})

    # ═══════════════════════════════════════
    # IDEA 1: (φ_SE/τ)^a × φ_AM^b × CN^c
    # "Quality-adjusted Bruggeman"
    # ═══════════════════════════════════════
    print("IDEA 1: (φ_SE/τ)^a × φ_AM^b × CN^c")
    for a in np.arange(1, 5.1, 0.5):
        for b in np.arange(0, 5.1, 0.5):
            for c in np.arange(-1, 3.1, 0.5):
                rhs = ps_tau**a * pa**b * cn**c
                test(f'(φ/τ)^{a:.1f}×φ_AM^{b:.1f}×CN^{c:.1f}', rhs)

    # ═══════════════════════════════════════
    # IDEA 2: (φ_SE×φ_AM)^a × (CN/τ)^b × G_path^c
    # ═══════════════════════════════════════
    print("IDEA 2: (φ×φ_AM)^a × (CN/τ)^b × contact^c")
    for a in np.arange(3, 8.1, 0.5):
        for b in np.arange(0, 3.1, 0.25):
            for c_name, c_var in [('1', np.ones(n)), ('G^-¼', gp**(-0.25)),
                                   ('G^¼', gp**0.25), ('SA^-¼', sa**(-0.25)),
                                   ('(Gd²)^-¼', gp_gd2**(-0.25)), ('BN^-¼', bn**(-0.25))]:
                rhs = ps_pa**a * cn_tau**b * c_var
                test(f'(φφ)^{a:.1f}×(CN/τ)^{b:.2f}×{c_name}', rhs)

    # ═══════════════════════════════════════
    # IDEA 3: φ_SE^a × φ_AM^b / τ^c × CN^d
    # Generalized with independent τ, CN exponents
    # ═══════════════════════════════════════
    print("IDEA 3: φ_SE^a × φ_AM^b / τ^c × CN^d")
    for a in np.arange(2, 7.1, 0.5):
        for b in np.arange(0, 5.1, 0.5):
            for c in np.arange(0, 3.1, 0.5):
                for d in np.arange(-1, 3.1, 0.5):
                    rhs = ps**a * pa**b / tau**c * cn**d
                    test(f'φ^{a:.1f}×φ_AM^{b:.1f}/τ^{c:.1f}×CN^{d:.1f}', rhs)

    # ═══════════════════════════════════════
    # IDEA 4: σ_brug_soft × φ_AM^a × CN^b
    # σ_brug but with τ^1 or τ^1.5 instead of τ^2
    # ═══════════════════════════════════════
    print("IDEA 4: σ_soft(τ^n) × φ_AM^a × CN^b × contact^c")
    for tn in [0.5, 0.75, 1, 1.25, 1.5]:
        s_soft = SG * ps * fp / tau**tn
        for a in np.arange(0, 5.1, 0.5):
            for b in np.arange(0, 3.1, 0.5):
                rhs = s_soft * pa**a * cn**b
                test(f'σ_soft(τ^{tn})×φ_AM^{a:.1f}×CN^{b:.1f}', rhs)
                # + G_path correction
                for gc in [-0.25, 0.25]:
                    rhs2 = rhs * gp**gc
                    test(f'σ_soft(τ^{tn})×φ_AM^{a:.1f}×CN^{b:.1f}×G^{gc}', rhs2)

    # ═══════════════════════════════════════
    # IDEA 5: PERCOLATION THEORY
    # σ ∝ (φ - φ_c)^t where t≈2 in 3D
    # φ_c from data? φ_c ≈ 0.15 (minimum φ_SE in data)
    # ═══════════════════════════════════════
    print("IDEA 5: (φ_SE - φ_c)^t × φ_AM^a × CN^b")
    for phi_c in [0.10, 0.12, 0.15, 0.18, 0.20]:
        phi_excess = np.clip(ps - phi_c, 0.01, None)
        for t in [1.5, 2, 2.5, 3, 4, 5, 6]:
            for a in [0, 1, 2, 3]:
                for b in [0, 1, 2]:
                    rhs = phi_excess**t * pa**a * cn**b
                    test(f'(φ-{phi_c})^{t}×φ_AM^{a}×CN^{b}', rhs)

    # ═══════════════════════════════════════
    # IDEA 6b: v3 + τ^k × φ_AM^m 보정 (USER INSIGHT!)
    # σ_brug는 그대로, 별도 τ 보정으로 τ² 과잉 상쇄
    # ═══════════════════════════════════════
    print("IDEA 6b: v3_champion × τ^k × φ_AM^m (σ_brug 유지, τ 보정 추가)")
    sigma_brug = SG * ps * fp / tau**2
    v3_rhs = sigma_brug * (gp * gd**2)**0.25 * cn**2  # v3 without C

    for k in np.arange(0, 2.1, 0.25):
        for m in np.arange(0, 5.1, 0.5):
            rhs = v3_rhs * tau**k * pa**m
            test(f'v3×τ^{k:.2f}×φ_AM^{m:.1f} (eff τ^-{2-k:.2f})', rhs)
            # + additional contact correction
            for gc in [-0.25, 0.25]:
                rhs2 = rhs * gp**gc
                test(f'v3×τ^{k:.2f}×φ_AM^{m:.1f}×G^{gc} (eff τ^-{2-k:.2f})', rhs2)

    # Also: v3 × (CN/τ)^k instead of CN² × extra_τ
    # This replaces CN² with (CN/τ)^k — unified variable
    print("IDEA 6c: σ_brug × C × (G×d²)^¼ × (CN/τ)^k × φ_AM^m")
    for k in np.arange(0.5, 3.1, 0.25):
        for m in np.arange(0, 5.1, 0.5):
            rhs = sigma_brug * (gp * gd**2)**0.25 * cn_tau**k * pa**m
            test(f'σ_brug×(G×d²)^¼×(CN/τ)^{k:.2f}×φ_AM^{m:.1f}', rhs)

    # ═══════════════════════════════════════
    # IDEA 6: SIGMOID-like τ correction
    # Instead of /τ^n, use 1/(1+τ^n) or exp(-τ^k)
    # These saturate for large τ instead of diverging
    # ═══════════════════════════════════════
    print("IDEA 6: Saturating τ corrections")
    for k in [0.5, 1, 1.5, 2]:
        tau_corr = 1 / (1 + (tau/2)**k)  # sigmoid at τ=2
        for a in np.arange(3, 7.1, 0.5):
            for b in np.arange(0, 5.1, 1):
                rhs = ps**a * pa**b * tau_corr
                test(f'φ^{a:.1f}×φ_AM^{b}×1/(1+(τ/2)^{k})', rhs)
                # + CN
                for c in [0.5, 1, 2]:
                    rhs2 = rhs * cn**c
                    test(f'φ^{a:.1f}×φ_AM^{b}×CN^{c}/(1+(τ/2)^{k})', rhs2)

    # ═══════════════════════════════════════
    # RESULTS — sorted by BALANCED score
    # ═══════════════════════════════════════
    results.sort(key=lambda x: -x['bal'])

    print(f"\n{'='*90}")
    print(f"TOP 30 BY BALANCED SCORE √(thick×mid×thin) — n={n}")
    print(f"{'='*90}")
    seen = set()
    count = 0
    for r in results:
        if r['name'] in seen: continue
        seen.add(r['name'])
        count += 1
        if count > 30: break
        flag = '★' if r['bal'] > 0.85 else ' '
        print(f"\n#{count:2d}{flag} BAL={r['bal']:.3f} ALL={r['r2']:.3f} C={r['C']:.4f} |err|={r['err']:.0f}%")
        print(f"   {r['name']}")
        print(f"   thick={r['thick']:.3f}  mid={r['mid']:.3f}  thin={r['thin']:.3f}")

    # TOP by ALL R²
    results.sort(key=lambda x: -x['r2'])
    print(f"\n{'='*90}")
    print(f"TOP 15 BY ALL R²")
    print(f"{'='*90}")
    seen2 = set()
    count2 = 0
    for r in results:
        if r['name'] in seen2: continue
        seen2.add(r['name'])
        count2 += 1
        if count2 > 15: break
        print(f"  #{count2:2d} ALL={r['r2']:.4f} BAL={r['bal']:.3f} thick={r['thick']:.3f} thin={r['thin']:.3f}  {r['name']}")

    # ═══════════════════════════════════════
    # SWEET SPOT: R²_all > 0.90 AND thin > 0.5
    # ═══════════════════════════════════════
    sweet = [r for r in results if r['r2'] > 0.90 and r['thin'] > 0.5]
    sweet.sort(key=lambda x: -x['bal'])
    print(f"\n{'='*90}")
    print(f"SWEET SPOT: ALL > 0.90 AND thin > 0.5 ({len(sweet)} formulas)")
    print(f"{'='*90}")
    seen3 = set()
    count3 = 0
    for r in sweet:
        if r['name'] in seen3: continue
        seen3.add(r['name'])
        count3 += 1
        if count3 > 20: break
        print(f"  #{count3:2d} BAL={r['bal']:.3f} ALL={r['r2']:.3f} thick={r['thick']:.3f} mid={r['mid']:.3f} thin={r['thin']:.3f}  {r['name']}")


if __name__ == '__main__':
    main()
