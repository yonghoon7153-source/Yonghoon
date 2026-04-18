"""
Electronic Universal: T/d 없이 thick+thin 하나의 식
===================================================
T/d가 φ, CN, por, cov의 proxy일 수 있다.
thin이면 자동으로 φ↑, por↑, CN 변화 → T/d 불필요?
"""
import json, os, numpy as np, warnings
from pathlib import Path
from itertools import combinations
warnings.filterwarnings('ignore')
WEBAPP = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'webapp')
SAM = 50.0

def load_data():
    rows = []
    for base in [Path(WEBAPP)/'results', Path(WEBAPP)/'archive']:
        if not base.is_dir(): continue
        for mp in base.rglob('full_metrics.json'):
            try:
                with open(mp) as f: m = json.load(f)
            except: continue
            sel = m.get('electronic_sigma_full_mScm', 0)
            if not sel or sel < 0.001: continue
            pa = max(m.get('phi_am', 0), 0.01)
            ps = m.get('phi_se', 0)
            cn = max(m.get('am_am_cn', 0.01), 0.01)
            cn_se = m.get('se_se_cn', 0)
            T = m.get('thickness_um', 0)
            tau = m.get('tortuosity_recommended', m.get('tortuosity_mean', 0))
            cov = max(m.get('coverage_AM_P_mean',
                          m.get('coverage_AM_S_mean',
                               m.get('coverage_AM_mean', 20))), 0.1) / 100
            r_am = max(m.get('r_AM_P', 0), m.get('r_AM_S', 0))
            d_am = r_am * 2 if r_am > 0.1 else 5.0
            por = m.get('porosity', 0)
            el_perc = max(m.get('electronic_percolating_fraction', 0), 0.01)
            sion = m.get('sigma_full_mScm', 0)
            am_area = m.get('am_am_mean_area', 0)
            am_cr = m.get('am_am_mean_contact_radius', 0)
            am_delta = m.get('am_am_mean_delta', 0)
            am_hop = m.get('am_am_mean_hop', 0)
            am_pres = m.get('am_am_mean_pressure', 0)
            if pa <= 0 or cn <= 0 or T <= 0: continue
            rows.append({
                'sel': sel, 'pa': pa, 'ps': max(ps, 0.01), 'cn': cn, 'cn_se': max(cn_se, 0.01),
                'tau': max(tau, 0.1), 'cov': cov, 'ratio': T / d_am,
                'por': max(por, 0.1), 'el_perc': el_perc,
                'sion': max(sion, 0.001),
                'am_area': max(am_area, 0.01), 'am_cr': max(am_cr, 0.01),
                'am_delta': max(am_delta, 0.001), 'am_hop': max(am_hop, 0.1),
                'am_pres': max(am_pres, 0.01),
                'name': mp.parent.name
            })
    seen = set(); u = []
    for r in rows:
        k = f"{r['pa']:.4f}_{r['ratio']:.1f}"
        if k not in seen: seen.add(k); u.append(r)
    return u

def r2l(a, p):
    la, lp = np.log(a), np.log(p)
    return 1 - np.sum((la - lp)**2) / np.sum((la - np.mean(la))**2)
def fitC(a, r):
    v = (r > 0) & np.isfinite(r)
    return float(np.exp(np.mean(np.log(a[v] / r[v])))) if v.sum() >= 3 else None
def loocv_C(sn, rhs):
    n = len(sn); la = np.log(sn); lr = np.log(rhs); errs = []
    for i in range(n):
        m = np.ones(n, bool); m[i] = False
        C_loo = float(np.exp(np.mean(la[m] - lr[m])))
        errs.append((la[i] - np.log(C_loo * rhs[i]))**2)
    return 1 - np.sum(errs) / np.sum((la - np.mean(la))**2)

def main():
    rows = load_data()
    if not rows:
        print("데이터 없음!"); return

    sel = np.array([r['sel'] for r in rows])
    pa = np.array([r['pa'] for r in rows])
    ps = np.array([r['ps'] for r in rows])
    cn = np.array([r['cn'] for r in rows])
    cn_se = np.array([r['cn_se'] for r in rows])
    tau = np.array([r['tau'] for r in rows])
    cov = np.array([r['cov'] for r in rows])
    ratio = np.array([r['ratio'] for r in rows])
    por = np.array([r['por'] for r in rows])
    el_perc = np.array([r['el_perc'] for r in rows])
    sion = np.array([r['sion'] for r in rows])
    am_area = np.array([r['am_area'] for r in rows])
    am_cr = np.array([r['am_cr'] for r in rows])
    am_delta = np.array([r['am_delta'] for r in rows])
    am_hop = np.array([r['am_hop'] for r in rows])
    am_pres = np.array([r['am_pres'] for r in rows])
    n = len(rows)

    thick = ratio >= 10; thin = ratio < 10
    print(f"n={n}: thick={thick.sum()}, thin={thin.sum()}")

    # ═══════════════════════════════════════════
    # 1. T/d 없는 power law grid search
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("1. T/d 없는 Universal power law")
    print("="*80)
    c1 = []
    for a in [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5]:
        for b in [1, 1.5, 2, 2.5]:
            for vn, va, ves in [
                ('cov', cov, [0, 0.5, 0.75, 1, 1.5]),
                ('τ', tau, [0, 0.25, 0.5, -0.5]),
                ('por', por, [0, 0.5, 1, -0.5, -1]),
            ]:
                for ve in ves:
                    rhs = SAM * pa**a * cn**b * va**ve
                    C = fitC(sel, rhs)
                    if C is None: continue
                    r2 = r2l(sel, C * rhs)
                    if r2 > 0.82:
                        vs = f'×{vn}^{ve}' if ve else ''
                        c1.append({'r2': r2, 'rhs': rhs,
                                   'label': f'φ^{a}×CN^{b}{vs}'})
    c1.sort(key=lambda x: -x['r2'])
    for i, r in enumerate(c1[:25]):
        cv = loocv_C(sel, r['rhs']) if i < 10 else 0
        cs = f" LOOCV={cv:.4f}" if cv else ""
        print(f"  #{i+1} R²={r['r2']:.4f}{cs}  {r['label']}")

    # ═══════════════════════════════════════════
    # 2. 2변수 추가 (T/d 없이)
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("2. φ^a × CN^b × X × Y (T/d 없이, 4변수)")
    print("="*80)
    c2 = []
    for a in [2, 3, 4]:
        for b in [1.5, 2]:
            for xn, xa, xes in [('cov', cov, [0.5, 1]), ('τ', tau, [0.5, -0.5]),
                                 ('por', por, [-0.5, 0.5, 1]), ('φ_SE', ps, [0.5, 1])]:
                for xe in xes:
                    for yn, ya, yes in [('cov', cov, [0.5, 1]), ('τ', tau, [0.5, -0.5]),
                                         ('por', por, [-0.5, 0.5]), ('el_perc', el_perc, [0.5, 1])]:
                        if xn == yn: continue
                        for ye in yes:
                            rhs = SAM * pa**a * cn**b * xa**xe * ya**ye
                            C = fitC(sel, rhs)
                            if C is None: continue
                            r2 = r2l(sel, C * rhs)
                            if r2 > 0.85:
                                c2.append({'r2': r2, 'rhs': rhs,
                                    'label': f'φ^{a}×CN^{b}×{xn}^{xe}×{yn}^{ye}'})
    c2.sort(key=lambda x: -x['r2'])
    for i, r in enumerate(c2[:20]):
        cv = loocv_C(sel, r['rhs']) if i < 10 else 0
        cs = f" LOOCV={cv:.4f}" if cv else ""
        print(f"  #{i+1} R²={r['r2']:.4f}{cs}  {r['label']}")

    # ═══════════════════════════════════════════
    # 3. Free regression (T/d 제외)
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("3. Free regression (T/d 제외)")
    print("="*80)
    var_pool = {
        'φ_AM': np.log(pa), 'CN': np.log(cn), 'τ': np.log(tau),
        'cov': np.log(cov), 'por': np.log(por), 'φ_SE': np.log(ps),
        'el_perc': np.log(el_perc), 'σ_ion': np.log(sion),
        'A_c': np.log(am_area), 'δ': np.log(am_delta), 'hop': np.log(am_hop),
    }
    log_sel = np.log(sel)
    c3 = []
    for nv in [3, 4, 5]:
        for vn in combinations(var_pool.keys(), nv):
            X = np.column_stack([var_pool[v] for v in vn] + [np.ones(n)])
            try:
                coefs, _, _, _ = np.linalg.lstsq(X, log_sel, rcond=None)
                pred = np.exp(X @ coefs)
                r2 = r2l(sel, pred)
                if r2 > 0.85:
                    c3.append({'r2': r2, 'pred': pred, 'coefs': coefs, 'vnames': vn})
            except: pass
    c3.sort(key=lambda x: -x['r2'])
    print(f"  Top 30 (n={n}):")
    for i, r in enumerate(c3[:30]):
        cv = loocv_C(sel, r['pred']) if i < 10 else 0
        cs = f" LOOCV={cv:.4f}" if cv else ""
        terms = ' × '.join(f"{v}^{r['coefs'][j]:.2f}" for j, v in enumerate(r['vnames']))
        clean = all(abs(r['coefs'][j] - round(r['coefs'][j]*4)/4) < 0.08 for j in range(len(r['vnames'])))
        flag = '★' if clean else ' '
        print(f"    #{i+1}{flag} R²={r['r2']:.4f}{cs}  ({len(r['vnames'])}v) {terms}")

    # ═══════════════════════════════════════════
    # 4. Free regression (T/d 포함, 비교용)
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("4. Free regression (T/d 포함, 비교)")
    print("="*80)
    var_pool_td = dict(var_pool)
    var_pool_td['T/d'] = np.log(ratio)
    c4 = []
    for nv in [3, 4, 5]:
        for vn in combinations(var_pool_td.keys(), nv):
            X = np.column_stack([var_pool_td[v] for v in vn] + [np.ones(n)])
            try:
                coefs, _, _, _ = np.linalg.lstsq(X, log_sel, rcond=None)
                pred = np.exp(X @ coefs)
                r2 = r2l(sel, pred)
                if r2 > 0.87:
                    c4.append({'r2': r2, 'pred': pred, 'coefs': coefs, 'vnames': vn})
            except: pass
    c4.sort(key=lambda x: -x['r2'])
    print(f"  Top 20 (n={n}):")
    for i, r in enumerate(c4[:20]):
        cv = loocv_C(sel, r['pred']) if i < 10 else 0
        cs = f" LOOCV={cv:.4f}" if cv else ""
        terms = ' × '.join(f"{v}^{r['coefs'][j]:.2f}" for j, v in enumerate(r['vnames']))
        has_td = '★' if 'T/d' not in r['vnames'] else ' '
        print(f"    #{i+1}{has_td} R²={r['r2']:.4f}{cs}  ({len(r['vnames'])}v) {terms}")

    # ═══════════════════════════════════════════
    # 5. FORM X 스타일 (⁴√, √, T/d 없이)
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("5. FORM X 스타일 (T/d 없이)")
    print("="*80)
    formx = [
        ('φ⁴×CN^(3/2)×cov', pa**4 * cn**1.5 * cov),
        ('φ⁴×CN^(3/2)×cov×√τ', pa**4 * cn**1.5 * cov * tau**0.5),
        ('φ⁴×CN^(3/2)×√cov', pa**4 * cn**1.5 * cov**0.5),
        ('φ²×CN²×√cov', pa**2 * cn**2 * cov**0.5),
        ('φ²×CN²×cov', pa**2 * cn**2 * cov),
        ('φ³×CN²×cov', pa**3 * cn**2 * cov),
        ('φ⁴×CN²×cov', pa**4 * cn**2 * cov),
        ('φ²×CN²', pa**2 * cn**2),
        ('⁴√[φ¹⁶×CN⁶×cov⁴]', (pa**16 * cn**6 * cov**4)**0.25),
        ('⁴√[φ¹⁶×CN⁶×cov⁴×τ²]', (pa**16 * cn**6 * cov**4 * tau**2)**0.25),
        ('⁴√[φ⁸×CN⁸×cov²]', (pa**8 * cn**8 * cov**2)**0.25),
        ('⁴√[φ¹²×CN⁸×cov⁴]', (pa**12 * cn**8 * cov**4)**0.25),
        ('√[φ⁴×CN²×cov]', np.sqrt(pa**4 * cn**2 * cov)),
        ('√[φ⁸×CN³×cov²]', np.sqrt(pa**8 * cn**3 * cov**2)),
        ('φ²×CN^(3/2)×cov', pa**2 * cn**1.5 * cov),
        ('φ³×CN^(3/2)×√cov', pa**3 * cn**1.5 * cov**0.5),
        ('φ³×CN^(3/2)×cov×√τ', pa**3 * cn**1.5 * cov * tau**0.5),
        ('φ⁴×CN×cov', pa**4 * cn * cov),
        ('φ²×CN²×cov/√τ', pa**2 * cn**2 * cov / tau**0.5),
        ('φ³×CN²×√cov', pa**3 * cn**2 * cov**0.5),
        ('φ²×CN²×por^(-0.5)', pa**2 * cn**2 * por**(-0.5)),
        ('φ⁴×CN^(3/2)×cov/por^0.5', pa**4 * cn**1.5 * cov / por**0.5),
    ]
    print(f"  {'Formula':45s} {'R²':>7s} {'LOOCV':>7s}")
    print("  " + "-" * 60)
    for label, rhs_raw in sorted(formx, key=lambda x: -r2l(sel, fitC(sel, SAM*x[1])*(SAM*x[1])) if fitC(sel, SAM*x[1]) else -999):
        rhs = SAM * rhs_raw
        C = fitC(sel, rhs)
        if C is None: continue
        r2 = r2l(sel, C * rhs)
        cv = loocv_C(sel, rhs) if r2 > 0.82 else 0
        cs = f"{cv:7.4f}" if cv else "      -"
        flag = '★' if r2 > 0.90 else '●' if r2 > 0.85 else ' '
        print(f"  {flag}{label:44s} {r2:7.4f} {cs}  C={C:.6f}")

    # ═══════════════════════════════════════════
    # 6. Per-regime breakdown (best T/d-free formula)
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("6. Best formula: thick vs thin 분리 R²")
    print("="*80)
    best_formulas = [
        ('φ⁴×CN^(3/2)×cov×√τ', SAM * pa**4 * cn**1.5 * cov * tau**0.5),
        ('φ⁴×CN^(3/2)×cov', SAM * pa**4 * cn**1.5 * cov),
        ('φ²×CN²×√cov', SAM * pa**2 * cn**2 * cov**0.5),
        ('φ²×CN²×cov', SAM * pa**2 * cn**2 * cov),
        ('φ³×CN²×cov', SAM * pa**3 * cn**2 * cov),
    ]
    for label, rhs in best_formulas:
        C = fitC(sel, rhs)
        if C is None: continue
        pred = C * rhs
        r2_all = r2l(sel, pred)
        # Thick/thin breakdown
        if thick.sum() >= 3:
            r2_tk = r2l(sel[thick], pred[thick])
            err_tk = np.mean(np.abs(sel[thick] - pred[thick]) / sel[thick] * 100)
        else:
            r2_tk = 0; err_tk = 0
        if thin.sum() >= 3:
            r2_tn = r2l(sel[thin], pred[thin])
            err_tn = np.mean(np.abs(sel[thin] - pred[thin]) / sel[thin] * 100)
        else:
            r2_tn = 0; err_tn = 0
        print(f"  {label:35s} ALL R²={r2_all:.4f} | thick R²={r2_tk:.4f} err={err_tk:.0f}% | thin R²={r2_tn:.4f} err={err_tn:.0f}%")

if __name__ == '__main__':
    main()
