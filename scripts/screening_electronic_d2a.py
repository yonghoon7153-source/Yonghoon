"""
δ²/A 심층: thin + thick + universal 전부
=========================================
thin: CN × (T/d)^(-½) × (δ²/A)^¾ × hop = R²=0.944
이 변수가 thick에서도 되는지? universal 가능한지?
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
            cn = max(m.get('am_am_cn', 0.01), 0.01)
            T = m.get('thickness_um', 0)
            tau = m.get('tortuosity_recommended', m.get('tortuosity_mean', 0))
            cov = max(m.get('coverage_AM_P_mean', m.get('coverage_AM_S_mean', m.get('coverage_AM_mean', 20))), 0.1) / 100
            r_am = max(m.get('r_AM_P', 0), m.get('r_AM_S', 0))
            d_am = r_am * 2 if r_am > 0.1 else 5.0
            por = m.get('porosity', 0)
            el_perc = max(m.get('electronic_percolating_fraction', 0), 0.01)
            am_area = m.get('am_am_mean_area', 0)
            am_delta = m.get('am_am_mean_delta', 0)
            am_hop = m.get('am_am_mean_hop', 0)
            am_cn_std = m.get('am_am_cn_std', 0)
            ps = m.get('phi_se', 0)
            if pa <= 0 or cn <= 0 or T <= 0 or am_delta <= 0 or am_area <= 0: continue
            d2a = am_delta**2 / am_area
            rows.append({
                'sel': sel, 'pa': pa, 'cn': cn, 'tau': max(tau, 0.1),
                'cov': cov, 'xi': T/d_am, 'por': max(por, 0.1),
                'el_perc': el_perc, 'ps': max(ps, 0.01),
                'd2a': d2a, 'hop': max(am_hop, 0.1),
                'cn_std': max(am_cn_std, 0.01),
                'd': am_delta, 'A': am_area,
                'name': mp.parent.name
            })
    seen = set(); u = []
    for r in rows:
        k = f"{r['pa']:.4f}_{r['xi']:.1f}"
        if k not in seen: seen.add(k); u.append(r)
    return u

def r2l(a, p):
    la, lp = np.log(a), np.log(p); return 1-np.sum((la-lp)**2)/np.sum((la-np.mean(la))**2)
def fitC(a, r):
    v = (r>0)&np.isfinite(r); return float(np.exp(np.mean(np.log(a[v]/r[v])))) if v.sum()>=3 else None
def loocv_C(sn, rhs):
    n=len(sn); la=np.log(sn); lr=np.log(rhs); errs=[]
    for i in range(n):
        m=np.ones(n,bool); m[i]=False
        C_loo=float(np.exp(np.mean(la[m]-lr[m])))
        errs.append((la[i]-np.log(C_loo*rhs[i]))**2)
    return 1-np.sum(errs)/np.sum((la-np.mean(la))**2)

def scan(label, sel, candidates):
    """Print sorted results"""
    candidates.sort(key=lambda x: -x['r2'])
    for i, r in enumerate(candidates[:20]):
        cv = loocv_C(sel, r['rhs']) if i < 10 else 0
        cs = f" LOOCV={cv:.4f}" if cv else ""
        print(f"  #{i+1} R²={r['r2']:.4f}{cs}  {r['label']}")

def main():
    rows = load_data()
    if not rows: print("데이터 없음!"); return
    V = {k: np.array([r[k] for r in rows]) for k in rows[0].keys() if k != 'name'}
    n = len(rows)
    thick = V['xi'] >= 10; thin = V['xi'] < 10
    print(f"n={n}: thick={thick.sum()}, thin={thin.sum()}")
    print(f"d2a(thin): {V['d2a'][thin].min():.6f}~{V['d2a'][thin].max():.6f}")
    print(f"d2a(thick): {V['d2a'][thick].min():.6f}~{V['d2a'][thick].max():.6f}")
    print(f"d2a correlation: thin r={np.corrcoef(np.log(V['d2a'][thin]),np.log(V['sel'][thin]))[0,1]:.3f}, "
          f"thick r={np.corrcoef(np.log(V['d2a'][thick]),np.log(V['sel'][thick]))[0,1]:.3f}, "
          f"all r={np.corrcoef(np.log(V['d2a']),np.log(V['sel']))[0,1]:.3f}")

    # ═══════════════════════════════════════════
    # A. THIN: d2a fine-tuning
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print(f"A. THIN: d2a fine-tuning (n={thin.sum()})")
    print("="*80)
    sn = V['sel'][thin]
    T = {k: V[k][thin] for k in V}

    # A1. CN × (T/d)^b × d2a^c × hop^d
    print(f"\n  [A1] CN^a × (T/d)^b × (δ²/A)^c × hop^d")
    c = []
    for a in [0.5, 0.75, 1, 1.25]:
        for b in [-0.75, -0.5, -0.25]:
            for cc in [0.5, 0.625, 0.75, 0.875, 1]:
                for dd in [0, 0.5, 0.75, 1, 1.25]:
                    rhs = SAM * T['cn']**a * T['xi']**b * T['d2a']**cc * T['hop']**dd
                    C = fitC(sn, rhs)
                    if C is None: continue
                    r2 = r2l(sn, C*rhs)
                    if r2 > 0.92:
                        hs = f'×hop^{dd}' if dd else ''
                        c.append({'r2':r2, 'rhs':rhs,
                                  'label':f'CN^{a}×(T/d)^{b}×(δ²/A)^{cc}{hs}'})
    scan("A1", sn, c)

    # A2. + extra variable
    print(f"\n  [A2] Best + extra")
    if c:
        best = c[0]
        for xn, xa in [('×φ', T['pa']), ('×φ²', T['pa']**2), ('×φ³', T['pa']**3),
                        ('×√τ', T['tau']**0.5), ('×cov', T['cov']), ('×√cov', T['cov']**0.5),
                        ('×el_p^½', T['el_perc']**0.5), ('×el_p', T['el_perc']),
                        ('×por^½', T['por']**0.5), ('×CN_std^-½', T['cn_std']**-0.5),
                        ('×φ_SE', T['ps']), ('×φ^-1', T['pa']**-1)]:
            rhs = best['rhs'] * xa
            C = fitC(sn, rhs)
            if C is None: continue
            r2 = r2l(sn, C*rhs)
            cv = loocv_C(sn, rhs)
            flag = '★' if r2 > best['r2'] else ' '
            print(f"    {flag}{best['label']}{xn}: R²={r2:.4f} LOOCV={cv:.4f}")

    # A3. FORM X style
    print(f"\n  [A3] FORM X thin")
    formx_t = [
        ('CN×(δ²/A)^¾/√(T/d)×hop', T['cn'] * T['d2a']**0.75 / T['xi']**0.5 * T['hop']),
        ('CN×(δ²/A)^¾/√(T/d)', T['cn'] * T['d2a']**0.75 / T['xi']**0.5),
        ('⁴√[CN⁴×(δ²/A)³/(T/d)²]×hop', (T['cn']**4 * T['d2a']**3 / T['xi']**2)**0.25 * T['hop']),
        ('⁴√[CN⁴×(δ²/A)³×hop⁴/(T/d)²]', (T['cn']**4 * T['d2a']**3 * T['hop']**4 / T['xi']**2)**0.25),
        ('√[CN²×(δ²/A)^(3/2)×hop²/(T/d)]', np.sqrt(T['cn']**2 * T['d2a']**1.5 * T['hop']**2 / T['xi'])),
        ('CN×(δ²/A)×hop/√(T/d)', T['cn'] * T['d2a'] * T['hop'] / T['xi']**0.5),
        ('CN×√(δ²/A)×hop/(T/d)^¼', T['cn'] * T['d2a']**0.5 * T['hop'] / T['xi']**0.25),
        ('CN×(δ²/A)^¾×hop^¾/√(T/d)', T['cn'] * T['d2a']**0.75 * T['hop']**0.75 / T['xi']**0.5),
        # without hop
        ('CN×(δ²/A)/√(T/d)', T['cn'] * T['d2a'] / T['xi']**0.5),
        ('CN×(δ²/A)^½/(T/d)^¼', T['cn'] * T['d2a']**0.5 / T['xi']**0.25),
    ]
    print(f"  {'Formula':55s} {'R²':>7s} {'LOOCV':>7s}")
    print("  " + "-" * 70)
    for label, rhs_raw in sorted(formx_t, key=lambda x: -r2l(sn, fitC(sn,SAM*x[1])*(SAM*x[1])) if fitC(sn,SAM*x[1]) else -999):
        rhs = SAM * rhs_raw
        C = fitC(sn, rhs)
        if C is None: continue
        r2 = r2l(sn, C*rhs)
        cv = loocv_C(sn, rhs) if r2 > 0.88 else 0
        cs = f"{cv:7.4f}" if cv else "      -"
        flag = '★' if r2 > 0.93 else '●' if r2 > 0.90 else ' '
        print(f"  {flag}{label:54s} {r2:7.4f} {cs}")

    # ═══════════════════════════════════════════
    # B. THICK: d2a가 도움되는가?
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print(f"B. THICK: d2a 효과 (n={thick.sum()})")
    print("="*80)
    st = V['sel'][thick]
    K = {k: V[k][thick] for k in V}

    thick_forms = [
        ('φ⁴×CN^(3/2)×cov×√τ [기존]', K['pa']**4 * K['cn']**1.5 * K['cov'] * K['tau']**0.5),
        ('φ⁴×CN^(3/2)×cov', K['pa']**4 * K['cn']**1.5 * K['cov']),
        # + d2a
        ('φ⁴×CN^(3/2)×cov×√τ×(δ²/A)^¼', K['pa']**4 * K['cn']**1.5 * K['cov'] * K['tau']**0.5 * K['d2a']**0.25),
        ('φ⁴×CN^(3/2)×cov×√τ×√(δ²/A)', K['pa']**4 * K['cn']**1.5 * K['cov'] * K['tau']**0.5 * K['d2a']**0.5),
        ('φ⁴×CN^(3/2)×(δ²/A)^½×√τ', K['pa']**4 * K['cn']**1.5 * K['d2a']**0.5 * K['tau']**0.5),
        ('φ⁴×(δ²/A)^¾×√τ', K['pa']**4 * K['d2a']**0.75 * K['tau']**0.5),
        ('φ⁴×CN×(δ²/A)^¾×√τ', K['pa']**4 * K['cn'] * K['d2a']**0.75 * K['tau']**0.5),
        # d2a replaces cov?
        ('φ⁴×CN^(3/2)×(δ²/A)^½', K['pa']**4 * K['cn']**1.5 * K['d2a']**0.5),
        ('φ⁴×CN^(3/2)×√(δ²/A)×cov', K['pa']**4 * K['cn']**1.5 * K['d2a']**0.5 * K['cov']),
    ]
    print(f"  {'Formula':50s} {'R²':>7s} {'LOOCV':>7s}")
    print("  " + "-" * 66)
    for label, rhs_raw in sorted(thick_forms, key=lambda x: -r2l(st, fitC(st,SAM*x[1])*(SAM*x[1])) if fitC(st,SAM*x[1]) else -999):
        rhs = SAM * rhs_raw
        C = fitC(st, rhs)
        if C is None: continue
        r2 = r2l(st, C*rhs)
        cv = loocv_C(st, rhs)
        flag = '★' if r2 > 0.96 else '●' if r2 > 0.94 else ' '
        print(f"  {flag}{label:49s} {r2:7.4f} {cv:7.4f}")

    # ═══════════════════════════════════════════
    # C. UNIVERSAL: d2a + φ + T/d
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print(f"C. UNIVERSAL: d2a 기반 (n={n})")
    print("="*80)
    cu = []
    for a in [2, 3, 3.5, 4]:
        for b in [-0.5, -0.25, 0]:
            for cc in [0.5, 0.75, 1]:
                for dd in [0, 0.5, 0.75, 1]:
                    for xn, xa, xes in [('CN', V['cn'], [0, 0.5, 1]),
                                         ('τ', V['tau'], [0, 0.25, 0.5]),
                                         ('cov', V['cov'], [0, 0.5, 1])]:
                        for xe in xes:
                            rhs = SAM * V['pa']**a * V['xi']**b * V['d2a']**cc * V['hop']**dd * xa**xe
                            C = fitC(V['sel'], rhs)
                            if C is None: continue
                            r2 = r2l(V['sel'], C*rhs)
                            if r2 > 0.89:
                                hs = f'×hop^{dd}' if dd else ''
                                ts = f'×(T/d)^{b}' if b else ''
                                xs = f'×{xn}^{xe}' if xe else ''
                                cu.append({'r2':r2, 'rhs':rhs,
                                    'label':f'φ^{a}{ts}×(δ²/A)^{cc}{hs}{xs}'})
    cu.sort(key=lambda x: -x['r2'])
    for i, r in enumerate(cu[:20]):
        cv = loocv_C(V['sel'], r['rhs']) if i < 10 else 0
        cs = f" LOOCV={cv:.4f}" if cv else ""
        pred = fitC(V['sel'], r['rhs']) * r['rhs']
        r2_tk = r2l(V['sel'][thick], pred[thick]) if thick.sum()>=3 else 0
        r2_tn = r2l(V['sel'][thin], pred[thin]) if thin.sum()>=3 else 0
        print(f"  #{i+1} R²={r['r2']:.4f}{cs} tk={r2_tk:.3f} tn={r2_tn:.3f}  {r['label']}")

    # ═══════════════════════════════════════════
    # D. FORM X style universal with d2a
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("D. FORM X universal with d2a")
    print("="*80)
    formx_u = [
        ('φ⁴×(δ²/A)^¾×√τ/(T/d)^¼', V['pa']**4 * V['d2a']**0.75 * V['tau']**0.5 / V['xi']**0.25),
        ('φ⁴×CN×(δ²/A)^¾/(T/d)^¼', V['pa']**4 * V['cn'] * V['d2a']**0.75 / V['xi']**0.25),
        ('φ⁴×(δ²/A)×√τ/(T/d)^¼', V['pa']**4 * V['d2a'] * V['tau']**0.5 / V['xi']**0.25),
        ('φ³×CN×(δ²/A)^¾×hop/(T/d)^½', V['pa']**3 * V['cn'] * V['d2a']**0.75 * V['hop'] / V['xi']**0.5),
        ('φ⁴×(δ²/A)^¾×hop^½/(T/d)^¼', V['pa']**4 * V['d2a']**0.75 * V['hop']**0.5 / V['xi']**0.25),
        ('⁴√[φ¹⁶×(δ²/A)³×τ²/(T/d)]', (V['pa']**16 * V['d2a']**3 * V['tau']**2 / V['xi'])**0.25),
        ('⁴√[φ¹⁶×(δ²/A)³×hop⁴/(T/d)²]', (V['pa']**16 * V['d2a']**3 * V['hop']**4 / V['xi']**2)**0.25),
        ('√[φ⁸×(δ²/A)^(3/2)×τ/(T/d)^½]', np.sqrt(V['pa']**8 * V['d2a']**1.5 * V['tau'] / V['xi']**0.5)),
        ('φ⁴×√(δ²/A)×cov×√τ', V['pa']**4 * V['d2a']**0.5 * V['cov'] * V['tau']**0.5),
        ('φ⁴×CN^½×(δ²/A)^¾×√τ/(T/d)^¼', V['pa']**4 * V['cn']**0.5 * V['d2a']**0.75 * V['tau']**0.5 / V['xi']**0.25),
    ]
    print(f"  {'Formula':55s} {'R²':>7s} {'LOOCV':>7s} {'tk':>6s} {'tn':>6s}")
    print("  " + "-" * 80)
    for label, rhs_raw in sorted(formx_u, key=lambda x: -r2l(V['sel'], fitC(V['sel'],SAM*x[1])*(SAM*x[1])) if fitC(V['sel'],SAM*x[1]) else -999):
        rhs = SAM * rhs_raw
        C = fitC(V['sel'], rhs)
        if C is None: continue
        pred = C * rhs
        r2 = r2l(V['sel'], pred)
        cv = loocv_C(V['sel'], rhs)
        r2_tk = r2l(V['sel'][thick], pred[thick]) if thick.sum()>=3 else 0
        r2_tn = r2l(V['sel'][thin], pred[thin]) if thin.sum()>=3 else 0
        flag = '★' if min(r2_tk,r2_tn) > 0.85 else '●' if r2 > 0.88 else ' '
        print(f"  {flag}{label:54s} {r2:7.4f} {cv:7.4f} {r2_tk:6.3f} {r2_tn:6.3f}")

    # ═══════════════════════════════════════════
    # E. GRAND CHAMPION TABLE
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("E. GRAND CHAMPION")
    print("="*80)
    print(f"\n  2-REGIME:")
    # Thick best
    rhs_tk = SAM * V['pa'][thick]**4 * V['cn'][thick]**1.5 * V['cov'][thick] * V['tau'][thick]**0.5
    C_tk = fitC(V['sel'][thick], rhs_tk)
    r2_tk = r2l(V['sel'][thick], C_tk*rhs_tk)
    cv_tk = loocv_C(V['sel'][thick], rhs_tk)
    print(f"    THICK: φ⁴×CN^(3/2)×cov×√τ  R²={r2_tk:.4f} LOOCV={cv_tk:.4f} C={C_tk:.4f}")
    # Thin best from A1
    if c:
        best_tn = c[0]
        C_tn = fitC(sn, best_tn['rhs'])
        r2_tn = best_tn['r2']
        cv_tn = loocv_C(sn, best_tn['rhs'])
        print(f"    THIN:  {best_tn['label']}  R²={r2_tn:.4f} LOOCV={cv_tn:.4f}")
    # Universal best
    if cu:
        best_u = cu[0]
        C_u = fitC(V['sel'], best_u['rhs'])
        pred_u = C_u * best_u['rhs']
        r2_u = best_u['r2']
        cv_u = loocv_C(V['sel'], best_u['rhs'])
        r2u_tk = r2l(V['sel'][thick], pred_u[thick])
        r2u_tn = r2l(V['sel'][thin], pred_u[thin])
        print(f"\n  UNIVERSAL:")
        print(f"    {best_u['label']}  R²={r2_u:.4f} LOOCV={cv_u:.4f} tk={r2u_tk:.3f} tn={r2u_tn:.3f}")

if __name__ == '__main__':
    main()
