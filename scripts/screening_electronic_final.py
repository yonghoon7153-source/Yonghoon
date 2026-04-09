"""
Electronic FINAL: 깔끔한 식 탐색
================================
발견: F_n/A_contact = P_contact (접촉 압력)이 thin의 핵심 변수
후막은 topology(φ, CN)만으로 R²=0.97

접근:
  A. Contact pressure P 기반 (thin)
  B. δ (penetration depth) 기반 (thin)
  C. 접촉 역학 없이 (기존 변수만, thin)
  D. Universal 1-regime (exp(π/(T/d)) + contact)
  E. FORM X 스타일 (⁴√, √)
  F. Thick 최종 확인
"""
import json, os, numpy as np, warnings
from pathlib import Path
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
            am_cn = max(m.get('am_am_cn', 0.01), 0.01)
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
            am_force = m.get('am_am_mean_force', 0)
            am_pres = m.get('am_am_mean_pressure', 0)
            am_hop = m.get('am_am_mean_hop', 0)
            # AM-AM path metrics
            am_gd = m.get('am_gb_density_mean', 0)
            am_gc = m.get('am_path_conductance_mean', 0)
            am_bn = m.get('am_path_bottleneck_mean', 0)
            # Ionic path metrics (SE-SE, proxy for AM structure)
            ion_gd = m.get('gb_density_mean', 0)
            ion_gc = m.get('path_conductance_mean', 0)
            if pa <= 0 or am_cn <= 0 or T <= 0: continue
            P_contact = am_force / am_area if am_area > 0 and am_force > 0 else max(am_pres, 0.01)
            rows.append({
                'sel': sel, 'pa': pa, 'ps': ps, 'am_cn': am_cn,
                'T': T, 'd_am': d_am, 'tau': max(tau, 0.1), 'cov': cov,
                'ratio': T / d_am, 'por': max(por, 0.1),
                'el_perc': el_perc, 'sion': max(sion, 0.001),
                'am_area': max(am_area, 0.01), 'am_cr': max(am_cr, 0.01),
                'am_delta': max(am_delta, 0.001), 'am_force': max(am_force, 0.01),
                'P_contact': max(P_contact, 0.01), 'am_hop': max(am_hop, 0.1),
                'am_gd': max(am_gd, 0.01), 'am_gc': max(am_gc, 0.001),
                'am_bn': max(am_bn, 0.001),
                'ion_gd': max(ion_gd, 0.01), 'ion_gc': max(ion_gc, 1e-6),
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
    am_cn = np.array([r['am_cn'] for r in rows])
    tau = np.array([r['tau'] for r in rows])
    cov = np.array([r['cov'] for r in rows])
    ratio = np.array([r['ratio'] for r in rows])
    por = np.array([r['por'] for r in rows])
    el_perc = np.array([r['el_perc'] for r in rows])
    P_c = np.array([r['P_contact'] for r in rows])
    am_delta = np.array([r['am_delta'] for r in rows])
    am_cr = np.array([r['am_cr'] for r in rows])
    am_hop = np.array([r['am_hop'] for r in rows])
    am_area = np.array([r['am_area'] for r in rows])
    # AM-AM path metrics
    am_gd = np.array([r['am_gd'] for r in rows])
    am_gc = np.array([r['am_gc'] for r in rows])
    am_bn = np.array([r['am_bn'] for r in rows])
    ion_gd = np.array([r['ion_gd'] for r in rows])
    ion_gc = np.array([r['ion_gc'] for r in rows])

    thick = ratio >= 10; thin = ratio < 10
    n = len(rows)
    print(f"n={n}: thick={thick.sum()}, thin={thin.sum()}")
    print(f"P_contact: {P_c.min():.0f}~{P_c.max():.0f} MPa")
    print(f"δ: {am_delta.min():.4f}~{am_delta.max():.4f} µm")
    print(f"a_c: {am_cr.min():.3f}~{am_cr.max():.3f} µm")

    sn = sel[thin]; pn = pa[thin]; psn = ps[thin]
    cn_n = am_cn[thin]; rn = ratio[thin]; porn = por[thin]
    cvn = cov[thin]; tn = tau[thin]; eln = el_perc[thin]
    P_n = P_c[thin]; dn = am_delta[thin]; cr_n = am_cr[thin]
    hop_n = am_hop[thin]; area_n = am_area[thin]
    gd_n = am_gd[thin]; gc_n = am_gc[thin]; bn_n = am_bn[thin]
    igc_n = ion_gc[thin]
    n_thin = thin.sum()
    print(f"AM Gc: {gc_n.min():.4f}~{gc_n.max():.4f}, BN: {bn_n.min():.4f}~{bn_n.max():.4f}")

    # A0. AM path conductance (NEW!)
    print(f"\n{'='*80}")
    print(f"A0. THIN: AM path conductance Gc, bottleneck BN (n={n_thin})")
    print("="*80)
    cA0 = []
    for var_name, var_arr in [('Gc', gc_n), ('BN', bn_n), ('√Gc', gc_n**0.5), ('√BN', bn_n**0.5)]:
        for v_e in [0.25, 0.5, 0.75, 1, 1.5, 2]:
            for cn_e in [0, 0.5, 1]:
                for td_e in [-0.5, -0.25, 0]:
                    rhs = SAM * var_arr**v_e * cn_n**cn_e * rn**td_e
                    C = fitC(sn, rhs)
                    if C is None: continue
                    r2 = r2l(sn, C * rhs)
                    if r2 > 0.85:
                        cn_s = f'×CN^{cn_e}' if cn_e else ''
                        td_s = f'×(T/d)^{td_e}' if td_e else ''
                        cv = loocv_C(sn, rhs)
                        cA0.append({'r2': r2, 'rhs': rhs,
                                    'label': f'{var_name}^{v_e}{cn_s}{td_s}', 'cv': cv})
    # Also try Gc × P combinations
    for gc_e in [0.25, 0.5, 1]:
        for p_e in [0.5, 1, 2]:
            for td_e in [-0.5, -0.25, 0]:
                rhs = SAM * gc_n**gc_e * P_n**p_e * rn**td_e
                C = fitC(sn, rhs)
                if C is None: continue
                r2 = r2l(sn, C * rhs)
                if r2 > 0.85:
                    td_s = f'×(T/d)^{td_e}' if td_e else ''
                    cv = loocv_C(sn, rhs)
                    cA0.append({'r2': r2, 'rhs': rhs,
                                'label': f'Gc^{gc_e}×P^{p_e}{td_s}', 'cv': cv})
    cA0.sort(key=lambda x: -x['r2'])
    for i, r in enumerate(cA0[:20]):
        cs = f" LOOCV={r['cv']:.4f}" if r.get('cv') else ""
        print(f"  #{i+1} R²={r['r2']:.4f}{cs}  {r['label']}")

    # A. Contact pressure
    print(f"\n{'='*80}")
    print(f"A. THIN: Contact pressure (n={n_thin})")
    print("="*80)
    cA = []
    for p_e in [0.5, 0.75, 1, 1.5, 2, 3, 4]:
        for cn_e in [0, 0.5, 0.75, 1, 1.5]:
            for td_e in [-0.5, -0.25, 0]:
                for xn, xa in [('', np.ones(n_thin)), ('×el_perc', eln),
                               ('×√cov', cvn**0.5), ('×cov', cvn), ('×√por', porn**0.5)]:
                    rhs = SAM * P_n**p_e * cn_n**cn_e * rn**td_e * xa
                    C = fitC(sn, rhs)
                    if C is None: continue
                    r2 = r2l(sn, C * rhs)
                    if r2 > 0.88:
                        cn_s = f'×CN^{cn_e}' if cn_e else ''
                        td_s = f'×(T/d)^{td_e}' if td_e else ''
                        cA.append({'r2': r2, 'rhs': rhs, 'label': f'P^{p_e}{cn_s}{td_s}{xn}'})
    cA.sort(key=lambda x: -x['r2'])
    for i, r in enumerate(cA[:20]):
        cv = loocv_C(sn, r['rhs']) if i < 10 else 0
        cs = f" LOOCV={cv:.4f}" if cv else ""
        print(f"  #{i+1} R²={r['r2']:.4f}{cs}  {r['label']}")

    # B. δ (penetration)
    print(f"\n{'='*80}")
    print(f"B. THIN: δ penetration (n={n_thin})")
    print("="*80)
    cB = []
    for d_e in [0.25, 0.5, 0.75, 1, 1.5, 2]:
        for cn_e in [0.5, 0.75, 1, 1.5]:
            for td_e in [-0.5, -0.25, 0]:
                for xn, xa in [('', np.ones(n_thin)), ('×el_perc', eln), ('×cov', cvn)]:
                    rhs = SAM * dn**d_e * cn_n**cn_e * rn**td_e * xa
                    C = fitC(sn, rhs)
                    if C is None: continue
                    r2 = r2l(sn, C * rhs)
                    if r2 > 0.85:
                        cB.append({'r2': r2, 'rhs': rhs,
                                   'label': f'δ^{d_e}×CN^{cn_e}' + (f'×(T/d)^{td_e}' if td_e else '') + xn})
    cB.sort(key=lambda x: -x['r2'])
    for i, r in enumerate(cB[:15]):
        cv = loocv_C(sn, r['rhs']) if i < 5 else 0
        cs = f" LOOCV={cv:.4f}" if cv else ""
        print(f"  #{i+1} R²={r['r2']:.4f}{cs}  {r['label']}")

    # C. No contact mechanics
    print(f"\n{'='*80}")
    print(f"C. THIN: 접촉 역학 없이 (n={n_thin})")
    print("="*80)
    cC = []
    for cn_e in [0.75, 1, 1.25, 1.5]:
        for por_e in [1, 1.5, 2, 2.5, 3]:
            for td_e in [-0.75, -0.5, -0.25, 0]:
                for vn, va, ves in [('cov', cvn, [0, 0.5, 0.75, 1, 1.5]),
                                     ('φ_SE', psn, [0, 0.5, 1, 1.5, 2]),
                                     ('el_perc', eln, [0, 0.5, 1])]:
                    for ve in ves:
                        rhs = SAM * cn_n**cn_e * porn**por_e * rn**td_e * va**ve
                        C = fitC(sn, rhs)
                        if C is None: continue
                        r2 = r2l(sn, C * rhs)
                        if r2 > 0.88:
                            cC.append({'r2': r2, 'rhs': rhs,
                                       'label': f'CN^{cn_e}×por^{por_e}' +
                                       (f'×(T/d)^{td_e}' if td_e else '') +
                                       (f'×{vn}^{ve}' if ve else '')})
    cC.sort(key=lambda x: -x['r2'])
    for i, r in enumerate(cC[:15]):
        cv = loocv_C(sn, r['rhs']) if i < 5 else 0
        cs = f" LOOCV={cv:.4f}" if cv else ""
        print(f"  #{i+1} R²={r['r2']:.4f}{cs}  {r['label']}")

    # D. Universal
    print(f"\n{'='*80}")
    print(f"D. UNIVERSAL (n={n})")
    print("="*80)
    exp_td = np.exp(np.pi / ratio)
    cD = []
    for a in [1, 1.5, 2, 3, 4]:
        for b in [1, 1.5, 2]:
            for cn, ca in [('', np.ones(n)), ('×cov', cov), ('×√cov', cov**0.5),
                           ('×P^0.5', P_c**0.5), ('×δ^0.5', am_delta**0.5),
                           ('×√a_c', am_cr**0.5), ('×√(δ/hop)', np.sqrt(am_delta/am_hop)),
                           ('×Gc^0.5', am_gc**0.5), ('×Gc', am_gc),
                           ('×BN^0.5', am_bn**0.5)]:
                for fn, fa in [('×exp(π/ξ)', exp_td),
                               ('×(T/d)^-0.5', ratio**-0.5),
                               ('', np.ones(n))]:
                    rhs = SAM * pa**a * am_cn**b * fa * ca
                    C = fitC(sel, rhs)
                    if C is None: continue
                    r2 = r2l(sel, C * rhs)
                    if r2 > 0.86:
                        cv = loocv_C(sel, rhs)
                        cD.append({'r2': r2, 'rhs': rhs, 'cv': cv,
                                   'label': f'φ^{a}×CN^{b}{fn}{cn}'})
    cD.sort(key=lambda x: -x['r2'])
    for i, r in enumerate(cD[:20]):
        cs = f" LOOCV={r['cv']:.4f}" if r.get('cv') else ""
        print(f"  #{i+1} R²={r['r2']:.4f}{cs}  {r['label']}")

    # E. FORM X style
    print(f"\n{'='*80}")
    print(f"E. FORM X 스타일 — thin (n={n_thin})")
    print("="*80)
    formx = [
        ('⁴√[P⁴×CN⁴/(T/d)]', (P_n**4 * cn_n**4 / rn)**0.25),
        ('⁴√[P⁴×CN⁴×cov²/(T/d)]', (P_n**4 * cn_n**4 * cvn**2 / rn)**0.25),
        ('⁴√[P³×CN⁴/(T/d)]', (P_n**3 * cn_n**4 / rn)**0.25),
        ('⁴√[P²×CN⁴×por⁴/(T/d)]', (P_n**2 * cn_n**4 * porn**4 / rn)**0.25),
        ('⁴√[δ²×CN⁴/(T/d)]', (dn**2 * cn_n**4 / rn)**0.25),
        ('⁴√[δ⁴×CN⁴/(T/d)²]', (dn**4 * cn_n**4 / rn**2)**0.25),
        ('⁴√[a_c⁴×CN⁴/(T/d)]', (cr_n**4 * cn_n**4 / rn)**0.25),
        ('⁴√[CN⁴×por⁸/(T/d)²]', (cn_n**4 * porn**8 / rn**2)**0.25),
        ('⁴√[CN⁴×por¹²/(T/d)²]', (cn_n**4 * porn**12 / rn**2)**0.25),
        ('⁴√[CN⁴×(φ_SE×cov)⁴/(T/d)]', (cn_n**4 * (psn*cvn)**4 / rn)**0.25),
        ('⁴√[el_p⁴×CN⁴/(T/d)]', (eln**4 * cn_n**4 / rn)**0.25),
        ('P×CN/(T/d)^¼', P_n * cn_n / rn**0.25),
        ('P^¾×CN×√cov/(T/d)^¼', P_n**0.75 * cn_n * cvn**0.5 / rn**0.25),
        ('√[P²×CN²/(T/d)]', np.sqrt(P_n**2 * cn_n**2 / rn)),
        ('√[P²×CN²×cov/(T/d)]', np.sqrt(P_n**2 * cn_n**2 * cvn / rn)),
        ('√[δ×CN²×por²/(T/d)]', np.sqrt(dn * cn_n**2 * porn**2 / rn)),
        ('√[CN²×por⁶/(T/d)]', np.sqrt(cn_n**2 * porn**6 / rn)),
        ('√[CN²×por⁶×φ_SE²×cov²/(T/d)]', np.sqrt(cn_n**2 * porn**6 * psn**2 * cvn**2 / rn)),
        ('δ^½×CN×(T/d)^(-¼)', dn**0.5 * cn_n * rn**(-0.25)),
        ('a_c×CN/(T/d)^¼', cr_n * cn_n / rn**0.25),
        # AM path conductance
        ('Gc×CN/(T/d)^½', gc_n * cn_n / rn**0.5),
        ('√[Gc²×CN²/(T/d)]', np.sqrt(gc_n**2 * cn_n**2 / rn)),
        ('Gc^½×CN/(T/d)^¼', gc_n**0.5 * cn_n / rn**0.25),
        ('⁴√[Gc⁴×CN⁴/(T/d)]', (gc_n**4 * cn_n**4 / rn)**0.25),
        ('Gc^¾×CN/(T/d)^¼', gc_n**0.75 * cn_n / rn**0.25),
        ('BN^½×CN/(T/d)^¼', bn_n**0.5 * cn_n / rn**0.25),
        ('√[BN×CN²/(T/d)]', np.sqrt(bn_n * cn_n**2 / rn)),
        ('Gc×P/(T/d)^½', gc_n * P_n / rn**0.5),
        ('√[Gc×P²/(T/d)]', np.sqrt(gc_n * P_n**2 / rn)),
    ]
    print(f"  {'Formula':50s} {'R²':>7s} {'LOOCV':>7s}")
    print("  " + "-" * 66)
    for label, rhs_raw in sorted(formx, key=lambda x: -r2l(sn, fitC(sn, SAM*x[1])*(SAM*x[1])) if fitC(sn, SAM*x[1]) else -999):
        rhs = SAM * rhs_raw
        C = fitC(sn, rhs)
        if C is None: continue
        r2 = r2l(sn, C * rhs)
        cv = loocv_C(sn, rhs) if r2 > 0.85 else 0
        cs = f"{cv:7.4f}" if cv else "      -"
        f = '★' if r2 > 0.93 else '●' if r2 > 0.90 else ' '
        print(f"  {f}{label:49s} {r2:7.4f} {cs}")

    # F. Thick final
    print(f"\n{'='*80}")
    print(f"F. THICK (n={thick.sum()})")
    print("="*80)
    st = sel[thick]; pt = pa[thick]; ct = am_cn[thick]
    tt = tau[thick]; cvt = cov[thick]; P_t = P_c[thick]
    dt = am_delta[thick]; crt = am_cr[thick]
    tf = [
        ('φ⁴×CN^(3/2)×cov', pt**4 * ct**1.5 * cvt),
        ('φ⁴×CN^(3/2)×cov×√τ', pt**4 * ct**1.5 * cvt * tt**0.5),
        ('φ⁴×CN^(3/2)×cov×√P', pt**4 * ct**1.5 * cvt * P_t**0.5),
        ('φ⁴×CN^(3/2)×cov×√δ', pt**4 * ct**1.5 * cvt * dt**0.5),
        ('φ⁴×CN^(3/2)×cov×√a_c', pt**4 * ct**1.5 * cvt * crt**0.5),
        ('⁴√[φ¹⁶×CN⁶×cov⁴×τ²]', (pt**16 * ct**6 * cvt**4 * tt**2)**0.25),
    ]
    print(f"  {'Formula':40s} {'R²':>7s} {'LOOCV':>7s}")
    print("  " + "-" * 56)
    for label, rhs_raw in sorted(tf, key=lambda x: -r2l(st, fitC(st, SAM*x[1])*(SAM*x[1])) if fitC(st, SAM*x[1]) else -999):
        rhs = SAM * rhs_raw
        C = fitC(st, rhs)
        if C is None: continue
        r2 = r2l(st, C * rhs)
        cv = loocv_C(st, rhs)
        f = '★' if r2 > 0.96 else '●' if r2 > 0.95 else ' '
        print(f"  {f}{label:39s} {r2:7.4f} {cv:7.4f}  C={C:.4f}")

if __name__ == '__main__':
    main()
