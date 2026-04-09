"""
v2.0g UNIFIED FRAMEWORK: Modified Bruggeman + Percolation + Contact
===================================================================
σ = σ_grain × (φ_SE - φ_c)^n_brug × f_perc^a / τ^b × C × (contact)^c × CN^d × φ_AM^e

Bruggeman을 버리지 않고, percolation threshold를 추가하여 확장.
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
            ps = m.get('phi_se', 0); pa = max(m.get('phi_am', 0), 0.01)
            t = m.get('tortuosity_recommended', m.get('tortuosity_mean', 0))
            fp = max(m.get('percolation_pct', 0)/100, 0.5)
            cn = m.get('se_se_cn', 0); gd = m.get('gb_density_mean', 0)
            gp = max(m.get('path_conductance_mean', 0), 1e-6)
            T = m.get('thickness_um', 0)
            sa = max(m.get('area_SE_SE_mean', 1e-6), 1e-6)
            if t<=0 or ps<=0 or cn<=0 or gd<=0 or T<=0: continue
            rows.append({'sn': sn, 'ps': ps, 'pa': pa, 'tau': t, 'fp': fp,
                'cn': cn, 'gd': gd, 'gp': gp, 'T': T, 'sa': sa, 'name': mp.parent.name})
    seen = set(); u = []
    for r in rows:
        k = f"{r['ps']:.4f}_{r['T']:.1f}_{r['tau']:.3f}"
        if k not in seen: seen.add(k); u.append(r)
    return u


def r2l(a, p):
    la, lp = np.log(a), np.log(p)
    return 1 - np.sum((la-lp)**2)/np.sum((la-np.mean(la))**2)

def fitC(a, r):
    v = (r>0)&np.isfinite(r)
    return float(np.exp(np.mean(np.log(a[v]/r[v])))) if v.sum()>=3 else None

def loocv(log_s, X):
    n = len(log_s); errs = np.zeros(n)
    for i in range(n):
        m = np.ones(n, bool); m[i] = False
        try:
            c, _, _, _ = np.linalg.lstsq(X[m], log_s[m], rcond=None)
            errs[i] = log_s[i] - X[i]@c
        except: errs[i] = 10
    ss_r = np.sum(errs**2); ss_t = np.sum((log_s-np.mean(log_s))**2)
    return 1 - ss_r/ss_t if ss_t>0 else -999


def main():
    rows = load_data()
    n = len(rows)
    print(f"n={n}\n")

    sn = np.array([r['sn'] for r in rows])
    ls = np.log(sn)
    ps = np.array([r['ps'] for r in rows])
    pa = np.array([r['pa'] for r in rows])
    tau = np.array([r['tau'] for r in rows])
    cn = np.array([r['cn'] for r in rows])
    gd = np.array([r['gd'] for r in rows])
    gp = np.array([r['gp'] for r in rows])
    fp = np.array([r['fp'] for r in rows])
    sa = np.array([r['sa'] for r in rows])

    thick = np.array([r['tau']<=1.5 for r in rows])
    mid = np.array([(r['tau']>1.5)&(r['tau']<=2.5) for r in rows])
    thin = np.array([r['tau']>2.5 for r in rows])

    gp_gd2 = gp * gd**2

    results = []

    def test(name, rhs, nf=1):
        C = fitC(sn, rhs)
        if C is None or C<=0: return
        pred = C*rhs
        r2a = r2l(sn, pred)
        r2t = r2l(sn[thick], pred[thick]) if thick.sum()>=2 else -99
        r2m = r2l(sn[mid], pred[mid]) if mid.sum()>=2 else -99
        r2n = r2l(sn[thin], pred[thin]) if thin.sum()>=2 else -99
        vs = [s for s in [r2t,r2m,r2n] if s>0]
        bal = np.prod(vs)**(1/len(vs)) if vs else 0
        err = np.mean(np.abs(sn-pred)/sn)*100
        results.append({'name':name,'r2':r2a,'C':C,'nf':nf,
            'thick':r2t,'mid':r2m,'thin':r2n,'bal':bal,'err':err,'pred':pred})

    # ═══════════════════════════════════════
    # FRAMEWORK: σ_grain × (φ-φ_c)^n × f_perc^a / τ^b × correction
    # ═══════════════════════════════════════
    for phi_c in [0.12, 0.15, 0.18, 0.20]:
        phi_excess = np.clip(ps - phi_c, 0.001, None)

        # ── LEVEL 1: Pure modified Bruggeman ──
        for nb in [1, 1.5, 2]:  # Bruggeman exponent on (φ-φ_c)
            for fa in [0, 0.5, 1]:  # f_perc exponent
                for tb in [0, 0.5, 1, 1.5, 2]:  # τ exponent
                    base = SG * phi_excess**nb * fp**fa / tau**tb
                    test(f'σ_g×(φ-{phi_c})^{nb}/τ^{tb}×f^{fa}', base)

                    # ── LEVEL 2: + φ_AM ──
                    for am in [1, 2, 3]:
                        rhs = base * pa**am
                        test(f'σ_g×(φ-{phi_c})^{nb}/τ^{tb}×f^{fa}×φ_AM^{am}', rhs)

                        # ── LEVEL 3: + CN ──
                        for cne in [0.5, 1, 2]:
                            rhs2 = rhs * cn**cne
                            test(f'σ_g×(φ-{phi_c})^{nb}/τ^{tb}×f^{fa}×φ_AM^{am}×CN^{cne}', rhs2)

                            # ── LEVEL 4: + contact correction ──
                            for gc in [0.25, -0.25]:
                                rhs3 = rhs2 * gp_gd2**gc
                                test(f'σ_g×(φ-{phi_c})^{nb}/τ^{tb}×f^{fa}×φ_AM^{am}×CN^{cne}×(Gd²)^{gc}', rhs3)
                                rhs3b = rhs2 * gp**gc
                                test(f'σ_g×(φ-{phi_c})^{nb}/τ^{tb}×f^{fa}×φ_AM^{am}×CN^{cne}×G^{gc}', rhs3b)

    # ═══════════════════════════════════════
    # RESULTS
    # ═══════════════════════════════════════
    results.sort(key=lambda x: -x['r2'])

    # Deduplicate
    seen = set(); unique = []
    for r in results:
        if r['name'] not in seen: seen.add(r['name']); unique.append(r)

    # SWEET SPOT
    sweet = [r for r in unique if r['r2']>0.90 and r['thin']>0.5]
    sweet.sort(key=lambda x: -x['bal'])

    print(f"{'='*90}")
    print(f"SWEET SPOT: ALL>0.90 AND thin>0.5 ({len(sweet)} formulas)")
    print(f"{'='*90}")
    for i, r in enumerate(sweet[:25]):
        print(f"  #{i+1:2d} BAL={r['bal']:.3f} ALL={r['r2']:.3f} thick={r['thick']:.3f} mid={r['mid']:.3f} thin={r['thin']:.3f} C={r['C']:.3f}")
        print(f"      {r['name']}")

    # TOP by ALL R²
    print(f"\n{'='*90}")
    print(f"TOP 20 by ALL R²")
    print(f"{'='*90}")
    for i, r in enumerate(unique[:20]):
        flag = '★' if r['thin']>0.5 else '⚠' if r['thin']>0 else '✗'
        print(f"  #{i+1:2d}{flag} ALL={r['r2']:.4f} thick={r['thick']:.3f} mid={r['mid']:.3f} thin={r['thin']:.3f}  {r['name']}")

    # LOOCV for top 5 sweet spot
    print(f"\n{'='*90}")
    print(f"LOOCV VALIDATION (top 5 sweet spot)")
    print(f"{'='*90}")
    for i, r in enumerate(sweet[:5]):
        # Rebuild X matrix for LOOCV
        # Extract formula components from name
        pred = r['pred']
        log_pred = np.log(pred / r['C'])  # = log(rhs)
        X = np.column_stack([log_pred, np.ones(n)])  # log(σ) = 1*log(rhs) + log(C)
        # Actually just compute LOOCV directly
        log_actual = np.log(sn)
        log_rhs = np.log(pred / r['C'])
        # C is the only free param → LOOCV for C
        loo_errors = []
        for j in range(n):
            m = np.ones(n, bool); m[j] = False
            C_loo = float(np.exp(np.mean(log_actual[m] - log_rhs[m])))
            pred_loo = C_loo * (pred[j] / r['C'])
            loo_errors.append((np.log(sn[j]) - np.log(pred_loo))**2)
        ss_res = np.sum(loo_errors)
        ss_tot = np.sum((log_actual - np.mean(log_actual))**2)
        r2_cv = 1 - ss_res/ss_tot
        print(f"  #{i+1} LOOCV={r2_cv:.4f} (train={r['r2']:.4f}, gap={r['r2']-r2_cv:.4f})")
        print(f"     {r['name']}")

    # ═══════════════════════════════════════
    # COMPARISON TABLE
    # ═══════════════════════════════════════
    print(f"\n{'='*90}")
    print(f"FINAL COMPARISON")
    print(f"{'='*90}")

    # v3
    sigma_brug = SG * ps * fp / tau**2
    v3_rhs = sigma_brug * (gp*gd**2)**0.25 * cn**2
    C3 = fitC(sn, v3_rhs); p3 = C3*v3_rhs

    # v4 percolation
    phi_ex = np.clip(ps-0.18, 0.001, None)
    v4_rhs = phi_ex**1.5 * pa**2 * cn
    C4 = fitC(sn, v4_rhs); p4 = C4*v4_rhs

    # Best sweet spot
    if sweet:
        bs = sweet[0]

    print(f"\n{'Model':<55s} {'ALL':>6s} {'thick':>6s} {'mid':>6s} {'thin':>6s}")
    print("-"*85)
    print(f"  {'v3: σ_brug×C×(Gd²)^¼×CN²':<53s} {r2l(sn,p3):6.3f} {r2l(sn[thick],p3[thick]):6.3f} {r2l(sn[mid],p3[mid]):6.3f} {r2l(sn[thin],p3[thin]):6.3f}")
    print(f"  {'v4_perc: (φ-0.18)^1.5×φ_AM²×CN':<53s} {r2l(sn,p4):6.3f} {r2l(sn[thick],p4[thick]):6.3f} {r2l(sn[mid],p4[mid]):6.3f} {r2l(sn[thin],p4[thin]):6.3f}")
    if sweet:
        print(f"  {'v4_unified: '+bs['name'][:40]:<53s} {bs['r2']:6.3f} {bs['thick']:6.3f} {bs['mid']:6.3f} {bs['thin']:6.3f}")


if __name__ == '__main__':
    main()
