"""
Electronic Mix: T/d + 접촉역학 짬뽕
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
            cn = max(m.get('am_am_cn', 0.01), 0.01)
            T = m.get('thickness_um', 0)
            tau = m.get('tortuosity_recommended', m.get('tortuosity_mean', 0))
            cov = max(m.get('coverage_AM_P_mean', m.get('coverage_AM_S_mean', m.get('coverage_AM_mean', 20))), 0.1) / 100
            r_am = max(m.get('r_AM_P', 0), m.get('r_AM_S', 0))
            d_am = r_am * 2 if r_am > 0.1 else 5.0
            por = m.get('porosity', 0)
            am_area = m.get('am_am_mean_area', 0)
            am_delta = m.get('am_am_mean_delta', 0)
            if pa <= 0 or cn <= 0 or T <= 0 or am_delta <= 0 or am_area <= 0: continue
            rows.append({
                'sel': sel, 'pa': pa, 'cn': cn, 'tau': max(tau, 0.1),
                'cov': cov, 'ratio': T/d_am, 'por': max(por, 0.1),
                'A': am_area, 'delta': am_delta, 'name': mp.parent.name
            })
    seen = set(); u = []
    for r in rows:
        k = f"{r['pa']:.4f}_{r['ratio']:.1f}"
        if k not in seen: seen.add(k); u.append(r)
    return u

def r2l(a, p):
    la, lp = np.log(a), np.log(p)
    return 1 - np.sum((la-lp)**2) / np.sum((la-np.mean(la))**2)
def fitC(a, r):
    v = (r>0) & np.isfinite(r)
    return float(np.exp(np.mean(np.log(a[v]/r[v])))) if v.sum()>=3 else None
def loocv_C(sn, rhs):
    n=len(sn); la=np.log(sn); lr=np.log(rhs); errs=[]
    for i in range(n):
        m=np.ones(n,bool); m[i]=False
        C_loo=float(np.exp(np.mean(la[m]-lr[m])))
        errs.append((la[i]-np.log(C_loo*rhs[i]))**2)
    return 1-np.sum(errs)/np.sum((la-np.mean(la))**2)

def main():
    rows = load_data()
    if not rows: print("데이터 없음!"); return
    sel = np.array([r['sel'] for r in rows])
    pa = np.array([r['pa'] for r in rows])
    cn = np.array([r['cn'] for r in rows])
    tau = np.array([r['tau'] for r in rows])
    cov = np.array([r['cov'] for r in rows])
    xi = np.array([r['ratio'] for r in rows])  # T/d
    por = np.array([r['por'] for r in rows])
    A = np.array([r['A'] for r in rows])
    d = np.array([r['delta'] for r in rows])
    n = len(rows)
    thick = xi >= 10; thin = xi < 10
    print(f"n={n}: thick={thick.sum()}, thin={thin.sum()}")

    # ═══════════════════════════════════════════
    # 1. φ^a × (T/d)^b × δ^c / A^e (4변수)
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("1. φ^a × (T/d)^b × δ^c / A^e")
    print("="*80)
    c1 = []
    for a in np.arange(2.5, 5.5, 0.25):
        for b in np.arange(-0.5, 0.25, 0.25):
            for c in np.arange(1, 3.5, 0.25):
                for e in np.arange(0.25, 2, 0.25):
                    rhs = SAM * pa**a * xi**b * d**c / A**e
                    C = fitC(sel, rhs)
                    if C is None: continue
                    r2 = r2l(sel, C*rhs)
                    if r2 > 0.89:
                        c1.append({'r2':r2, 'rhs':rhs,
                                   'label':f'φ^{a}×(T/d)^{b}×δ^{c}/A^{e}'})
    c1.sort(key=lambda x: -x['r2'])
    for i, r in enumerate(c1[:20]):
        cv = loocv_C(sel, r['rhs']) if i < 10 else 0
        cs = f" LOOCV={cv:.4f}" if cv else ""
        print(f"  #{i+1} R²={r['r2']:.4f}{cs}  {r['label']}")

    # ═══════════════════════════════════════════
    # 2. + τ, CN, cov (5변수)
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("2. φ^a × (T/d)^b × δ^c / A^e × X")
    print("="*80)
    c2 = []
    for a in [3.5, 3.75, 4]:
        for b in [-0.25, -0.5]:
            for c in [2, 2.25, 2.5]:
                for e in [0.75, 1]:
                    base = SAM * pa**a * xi**b * d**c / A**e
                    for xn, xa in [('×√τ', tau**0.5), ('×τ^¾', tau**0.75),
                                    ('×CN^½', cn**0.5), ('×CN', cn),
                                    ('×√cov', cov**0.5), ('×cov', cov),
                                    ('×CN^-½', cn**-0.5)]:
                        rhs = base * xa
                        C = fitC(sel, rhs)
                        if C is None: continue
                        r2 = r2l(sel, C*rhs)
                        if r2 > 0.90:
                            c2.append({'r2':r2, 'rhs':rhs,
                                       'label':f'φ^{a}×(T/d)^{b}×δ^{c}/A^{e}{xn}'})
    c2.sort(key=lambda x: -x['r2'])
    for i, r in enumerate(c2[:20]):
        cv = loocv_C(sel, r['rhs']) if i < 10 else 0
        cs = f" LOOCV={cv:.4f}" if cv else ""
        print(f"  #{i+1} R²={r['r2']:.4f}{cs}  {r['label']}")

    # ═══════════════════════════════════════════
    # 3. exp(π/(T/d)) × 접촉역학
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("3. exp(π/(T/d)) × 접촉역학")
    print("="*80)
    exp_td = np.exp(np.pi / xi)
    c3 = []
    for a in [2, 2.5, 3, 3.5, 4]:
        for c in [1, 1.5, 2, 2.5]:
            for e in [0, 0.5, 0.75, 1]:
                for xn, xa in [('', np.ones(n)), ('×√τ', tau**0.5),
                                ('×CN', cn), ('×√cov', cov**0.5)]:
                    rhs = SAM * pa**a * exp_td * d**c / np.clip(A**e, 1e-10, None) * xa
                    C = fitC(sel, rhs)
                    if C is None: continue
                    r2 = r2l(sel, C*rhs)
                    if r2 > 0.90:
                        a_s = f'/A^{e}' if e else ''
                        c3.append({'r2':r2, 'rhs':rhs,
                                   'label':f'φ^{a}×exp(π/ξ)×δ^{c}{a_s}{xn}'})
    c3.sort(key=lambda x: -x['r2'])
    for i, r in enumerate(c3[:20]):
        cv = loocv_C(sel, r['rhs']) if i < 10 else 0
        cs = f" LOOCV={cv:.4f}" if cv else ""
        print(f"  #{i+1} R²={r['r2']:.4f}{cs}  {r['label']}")

    # ═══════════════════════════════════════════
    # 4. FORM X style 깔끔한 식
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("4. FORM X style (T/d + δ + A)")
    print("="*80)
    formx = [
        ('φ⁴×δ²/(A×√(T/d))', pa**4 * d**2 / (A * xi**0.5)),
        ('φ⁴×δ²×√τ/(A×√(T/d))', pa**4 * d**2 * tau**0.5 / (A * xi**0.5)),
        ('φ⁴×δ²/(A^¾×√(T/d))', pa**4 * d**2 / (A**0.75 * xi**0.5)),
        ('φ⁴×δ²×√τ/(A^¾×√(T/d))', pa**4 * d**2 * tau**0.5 / (A**0.75 * xi**0.5)),
        ('⁴√[φ¹⁶×δ⁸/(A⁴×(T/d)²)]', (pa**16 * d**8 / (A**4 * xi**2))**0.25),
        ('⁴√[φ¹⁶×δ⁸×τ²/(A⁴×(T/d)²)]', (pa**16 * d**8 * tau**2 / (A**4 * xi**2))**0.25),
        ('⁴√[φ¹⁶×δ⁸/(A³×(T/d)²)]', (pa**16 * d**8 / (A**3 * xi**2))**0.25),
        ('⁴√[φ¹⁶×δ⁸×τ²/(A³×(T/d)²)]', (pa**16 * d**8 * tau**2 / (A**3 * xi**2))**0.25),
        ('√[φ⁸×δ⁴/(A²×(T/d))]', np.sqrt(pa**8 * d**4 / (A**2 * xi))),
        ('√[φ⁸×δ⁴×τ/(A²×(T/d))]', np.sqrt(pa**8 * d**4 * tau / (A**2 * xi))),
        ('√[φ⁸×δ⁴/(A^(3/2)×(T/d))]', np.sqrt(pa**8 * d**4 / (A**1.5 * xi))),
        ('φ⁴×exp(π/ξ)×δ²/A', pa**4 * exp_td * d**2 / A),
        ('φ⁴×exp(π/ξ)×δ/√A', pa**4 * exp_td * d / A**0.5),
        ('φ⁴×exp(π/ξ)×√τ×δ²/A', pa**4 * exp_td * tau**0.5 * d**2 / A),
        ('φ²×CN²×√cov×exp(π/ξ) [기존]', pa**2 * cn**2 * cov**0.5 * exp_td),
        ('φ⁴×CN^(3/2)×cov×√τ [thick only]', pa**4 * cn**1.5 * cov * tau**0.5),
        # (T/d)^-0.25
        ('φ⁴×δ²/(A^¾×(T/d)^¼)', pa**4 * d**2 / (A**0.75 * xi**0.25)),
        ('φ⁴×√τ×δ²/(A^¾×(T/d)^¼)', pa**4 * tau**0.5 * d**2 / (A**0.75 * xi**0.25)),
    ]
    print(f"  {'Formula':50s} {'R²':>7s} {'LOOCV':>7s}")
    print("  " + "-" * 66)
    for label, rhs_raw in sorted(formx, key=lambda x: -r2l(sel, fitC(sel,SAM*x[1])*(SAM*x[1])) if fitC(sel,SAM*x[1]) else -999):
        rhs = SAM * rhs_raw
        C = fitC(sel, rhs)
        if C is None: continue
        r2 = r2l(sel, C*rhs)
        cv = loocv_C(sel, rhs) if r2 > 0.85 else 0
        cs = f"{cv:7.4f}" if cv else "      -"
        flag = '★' if r2 > 0.92 else '●' if r2 > 0.89 else ' '
        print(f"  {flag}{label:49s} {r2:7.4f} {cs}  C={C:.4f}")

    # ═══════════════════════════════════════════
    # 5. thick/thin breakdown (top 5)
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("5. thick vs thin")
    print("="*80)
    all_cands = (c1[:3] if c1 else []) + (c2[:3] if c2 else []) + (c3[:3] if c3 else [])
    # Add FORM X manually
    for label, rhs_raw in [
        ('φ⁴×δ²/(A^¾×√(T/d))', pa**4 * d**2 / (A**0.75 * xi**0.5)),
        ('φ⁴×√τ×δ²/(A^¾×(T/d)^¼)', pa**4 * tau**0.5 * d**2 / (A**0.75 * xi**0.25)),
        ('φ⁴×exp(π/ξ)×δ²/A', pa**4 * exp_td * d**2 / A),
    ]:
        rhs = SAM * rhs_raw
        C = fitC(sel, rhs)
        if C: all_cands.append({'r2': r2l(sel, C*rhs), 'rhs': rhs, 'label': label})

    for r in sorted(all_cands, key=lambda x: -x['r2'])[:10]:
        C = fitC(sel, r['rhs']); pred = C * r['rhs']
        r2_all = r2l(sel, pred)
        r2_tk = r2l(sel[thick], pred[thick]) if thick.sum()>=3 else 0
        r2_tn = r2l(sel[thin], pred[thin]) if thin.sum()>=3 else 0
        err_tk = np.mean(np.abs(sel[thick]-pred[thick])/sel[thick]*100)
        err_tn = np.mean(np.abs(sel[thin]-pred[thin])/sel[thin]*100)
        print(f"  {r['label']}")
        print(f"    ALL={r2_all:.4f} | thick={r2_tk:.4f} err={err_tk:.0f}% | thin={r2_tn:.4f} err={err_tn:.0f}%")

if __name__ == '__main__':
    main()
