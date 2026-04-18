"""
Electronic τ-based Universal: δ와 A 분리, τ가 T/d proxy
=======================================================
핵심: δ/R로 합치면 R²=0.79, δ와 A 따로면 R²=0.91
τ ↔ T/d: r=-0.72 → τ가 thin/thick를 자연스럽게 분리
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
            am_hop = m.get('am_am_mean_hop', 0)
            if pa <= 0 or cn <= 0 or T <= 0 or am_delta <= 0 or am_area <= 0: continue
            rows.append({
                'sel': sel, 'pa': pa, 'cn': cn, 'tau': max(tau, 0.1),
                'cov': cov, 'ratio': T/d_am, 'por': max(por, 0.1),
                'A': am_area, 'delta': am_delta, 'hop': max(am_hop, 0.1),
                'name': mp.parent.name
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
    ratio = np.array([r['ratio'] for r in rows])
    por = np.array([r['por'] for r in rows])
    A = np.array([r['A'] for r in rows])
    delta = np.array([r['delta'] for r in rows])
    hop = np.array([r['hop'] for r in rows])
    n = len(rows)
    thick = ratio >= 10; thin = ratio < 10
    print(f"n={n}: thick={thick.sum()}, thin={thin.sum()}")

    # ═══════════════════════════════════════════
    # 1. φ^a × τ^b × δ^c / A^d (핵심 4변수)
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("1. φ^a × τ^b × δ^c / A^d")
    print("="*80)
    c1 = []
    for a in np.arange(3, 5.5, 0.25):
        for b in np.arange(0.25, 1.25, 0.25):
            for c in np.arange(1, 3.5, 0.25):
                for d in np.arange(0.25, 2, 0.25):
                    rhs = SAM * pa**a * tau**b * delta**c / A**d
                    C = fitC(sel, rhs)
                    if C is None: continue
                    r2 = r2l(sel, C*rhs)
                    if r2 > 0.87:
                        c1.append({'r2':r2, 'rhs':rhs, 'a':a, 'b':b, 'c':c, 'd':d,
                                   'label':f'φ^{a}×τ^{b}×δ^{c}/A^{d}'})
    c1.sort(key=lambda x: -x['r2'])
    for i, r in enumerate(c1[:25]):
        cv = loocv_C(sel, r['rhs']) if i < 10 else 0
        cs = f" LOOCV={cv:.4f}" if cv else ""
        print(f"  #{i+1} R²={r['r2']:.4f}{cs}  {r['label']}")

    # ═══════════════════════════════════════════
    # 2. + CN (5변수)
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("2. φ^a × τ^b × δ^c / A^d × CN^e")
    print("="*80)
    c2 = []
    if c1:
        # Top 10 from c1에 CN 추가
        for base in c1[:10]:
            for e in [-0.5, 0, 0.5, 1]:
                rhs = base['rhs'] * cn**e
                C = fitC(sel, rhs)
                if C is None: continue
                r2 = r2l(sel, C*rhs)
                if r2 > 0.88:
                    cn_s = f'×CN^{e}' if e else ''
                    c2.append({'r2':r2, 'rhs':rhs,
                               'label':f"{base['label']}{cn_s}"})
    c2.sort(key=lambda x: -x['r2'])
    for i, r in enumerate(c2[:15]):
        cv = loocv_C(sel, r['rhs']) if i < 8 else 0
        cs = f" LOOCV={cv:.4f}" if cv else ""
        print(f"  #{i+1} R²={r['r2']:.4f}{cs}  {r['label']}")

    # ═══════════════════════════════════════════
    # 3. + cov, hop (6변수)
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("3. + cov, hop 추가")
    print("="*80)
    if c1:
        best = c1[0]
        for xn, xa in [('×cov', cov), ('×√cov', cov**0.5),
                        ('×hop^0.5', hop**0.5), ('×hop^-0.5', hop**-0.5),
                        ('×por^-0.5', por**-0.5), ('×por^0.5', por**0.5),
                        ('×CN', cn), ('×CN^0.5', cn**0.5),
                        ('×CN^-0.5', cn**-0.5)]:
            rhs = best['rhs'] * xa
            C = fitC(sel, rhs)
            if C is None: continue
            r2 = r2l(sel, C*rhs)
            cv = loocv_C(sel, rhs)
            flag = '★' if r2 > best['r2'] else ' '
            print(f"  {flag}{best['label']}{xn}: R²={r2:.4f} LOOCV={cv:.4f}")

    # ═══════════════════════════════════════════
    # 4. FORM X style (⁴√, √)
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("4. FORM X style (δ, A 분리)")
    print("="*80)
    formx = [
        ('⁴√[φ¹⁶×τ²×δ⁸/A⁴]', (pa**16 * tau**2 * delta**8 / A**4)**0.25),
        ('⁴√[φ²⁰×τ⁴×δ⁸/A⁴]', (pa**20 * tau**4 * delta**8 / A**4)**0.25),
        ('⁴√[φ¹⁶×τ²×δ⁴/A²]', (pa**16 * tau**2 * delta**4 / A**2)**0.25),
        ('⁴√[φ¹⁶×τ⁴×δ⁸/A⁴]', (pa**16 * tau**4 * delta**8 / A**4)**0.25),
        ('φ⁴×√τ×δ²/A', pa**4 * tau**0.5 * delta**2 / A),
        ('φ⁴×√τ×δ^(3/2)/A^(3/4)', pa**4 * tau**0.5 * delta**1.5 / A**0.75),
        ('φ⁴×τ^(3/4)×δ²/A', pa**4 * tau**0.75 * delta**2 / A),
        ('φ⁴×√τ×δ/√A', pa**4 * tau**0.5 * delta / A**0.5),
        ('φ⁴×√τ×√δ/⁴√A', pa**4 * tau**0.5 * delta**0.5 / A**0.25),
        ('√[φ⁸×τ×δ⁴/A²]', np.sqrt(pa**8 * tau * delta**4 / A**2)),
        ('√[φ⁸×τ×δ²/A]', np.sqrt(pa**8 * tau * delta**2 / A)),
        ('√[φ⁸×τ×δ³/A^(3/2)]', np.sqrt(pa**8 * tau * delta**3 / A**1.5)),
        ('φ⁴×√τ×(δ/√A)²', pa**4 * tau**0.5 * (delta / A**0.5)**2),
        ('φ⁴×√τ×(δ²/A)^(3/4)', pa**4 * tau**0.5 * (delta**2/A)**0.75),
        ('φ⁴×√τ×(δ²/A)^(1/2)', pa**4 * tau**0.5 * (delta**2/A)**0.5),
        # with CN
        ('φ⁴×√τ×δ²/(A×CN^0.5)', pa**4 * tau**0.5 * delta**2 / (A * cn**0.5)),
        ('φ⁴×√τ×δ²×CN^0.5/A', pa**4 * tau**0.5 * delta**2 * cn**0.5 / A),
        # with cov
        ('φ⁴×√τ×δ²×cov/A', pa**4 * tau**0.5 * delta**2 * cov / A),
        ('φ⁴×√τ×δ²×√cov/A', pa**4 * tau**0.5 * delta**2 * cov**0.5 / A),
    ]
    print(f"  {'Formula':50s} {'R²':>7s} {'LOOCV':>7s}")
    print("  " + "-" * 66)
    for label, rhs_raw in sorted(formx, key=lambda x: -r2l(sel, fitC(sel,SAM*x[1])*(SAM*x[1])) if fitC(sel,SAM*x[1]) else -999):
        rhs = SAM * rhs_raw
        C = fitC(sel, rhs)
        if C is None: continue
        r2 = r2l(sel, C*rhs)
        cv = loocv_C(sel, rhs) if r2 > 0.82 else 0
        cs = f"{cv:7.4f}" if cv else "      -"
        flag = '★' if r2 > 0.88 else '●' if r2 > 0.85 else ' '
        print(f"  {flag}{label:49s} {r2:7.4f} {cs}  C={C:.4f}")

    # ═══════════════════════════════════════════
    # 5. thick/thin breakdown
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("5. thick vs thin 분리")
    print("="*80)
    checks = c1[:5] if c1 else []
    for r in checks:
        rhs = r['rhs']; C = fitC(sel, rhs); pred = C * rhs
        r2_all = r2l(sel, pred)
        r2_tk = r2l(sel[thick], pred[thick]) if thick.sum()>=3 else 0
        r2_tn = r2l(sel[thin], pred[thin]) if thin.sum()>=3 else 0
        err_tk = np.mean(np.abs(sel[thick]-pred[thick])/sel[thick]*100)
        err_tn = np.mean(np.abs(sel[thin]-pred[thin])/sel[thin]*100)
        print(f"  {r['label']}")
        print(f"    ALL={r2_all:.4f} | thick={r2_tk:.4f} err={err_tk:.0f}% | thin={r2_tn:.4f} err={err_tn:.0f}%")

    # ═══════════════════════════════════════════
    # 6. Per-case (best formula)
    # ═══════════════════════════════════════════
    if c1:
        best = c1[0]
        C = fitC(sel, best['rhs']); pred = C * best['rhs']
        err = np.abs(sel - pred) / sel * 100
        print(f"\n{'='*80}")
        print(f"6. Per-case: {best['label']} (R²={best['r2']:.4f})")
        print("="*80)
        print(f"  {'Name':35s} {'σ':>6s} {'pred':>6s} {'err%':>5s} {'T/d':>5s} {'δ':>6s} {'A':>6s} {'τ':>5s}")
        for i in np.argsort(-err):
            r = rows[i]
            regime = 'TK' if ratio[i]>=10 else 'tn'
            print(f"  {r['name']:35s} {sel[i]:6.2f} {pred[i]:6.2f} {err[i]:5.0f}% {ratio[i]:5.1f} {delta[i]:6.4f} {A[i]:6.3f} {tau[i]:5.2f} {regime}")

if __name__ == '__main__':
    main()
