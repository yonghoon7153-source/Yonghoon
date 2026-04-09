"""
Electronic Hertz: δ/R 물리적 의미 + universal 심층 탐구
======================================================
Hertz contact theory:
  F = (4/3) E* √R δ^(3/2)
  A = π R δ
  a = √(Rδ)
  P_mean = F/A = (4/3π) E* √(δ/R)

따라서:
  δ/R = 상대 침투 깊이 = (3πP/(4E*))²
  a/R = √(δ/R) = contact radius ratio
  δ²/A = δ²/(πRδ) = δ/(πR) ∝ δ/R

모든 contact mechanics 변수가 결국 δ/R 하나로 귀결!

물리적 의미:
  σ_el ∝ φ^a × (δ/R)^b × CN^c × √τ
  = φ^a × (상대 침투)^b × (연결 수)^c × (경로 직진성)^(1/2)

  φ: AM이 얼마나 있는가 (volume)
  δ/R: AM끼리 얼마나 세게 눌리는가 (contact quality)  
  CN: AM끼리 몇 개 접촉하는가 (connectivity)
  τ: 경로가 얼마나 직진하는가 (path efficiency)
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
            T = m.get('thickness_um', 0)
            tau = m.get('tortuosity_recommended', m.get('tortuosity_mean', 0))
            cov = max(m.get('coverage_AM_P_mean', m.get('coverage_AM_S_mean', m.get('coverage_AM_mean', 20))), 0.1) / 100
            r_am = max(m.get('r_AM_P', 0), m.get('r_AM_S', 0))
            d_am = r_am * 2 if r_am > 0.1 else 5.0
            por = m.get('porosity', 0)
            el_perc = max(m.get('electronic_percolating_fraction', 0), 0.01)
            am_delta = m.get('am_am_mean_delta', 0)
            am_area = m.get('am_am_mean_area', 0)
            am_cr = m.get('am_am_mean_contact_radius', 0)
            am_hop = m.get('am_am_mean_hop', 0)
            if pa <= 0 or cn <= 0 or T <= 0 or am_delta <= 0: continue
            R = d_am / 2
            dR = am_delta / R
            aR = am_cr / R if am_cr > 0 else np.sqrt(dR)
            d2A = am_delta**2 / am_area if am_area > 0 else dR / np.pi
            rows.append({
                'sel': sel, 'pa': pa, 'ps': max(ps, 0.01), 'cn': cn,
                'tau': max(tau, 0.1), 'cov': cov, 'ratio': T/d_am,
                'por': max(por, 0.1), 'el_perc': el_perc,
                'dR': dR, 'd2A': d2A, 'aR': aR,
                'delta': am_delta, 'A': max(am_area, 0.01),
                'hop': max(am_hop, 0.1), 'R': R,
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
    dR = np.array([r['dR'] for r in rows])
    d2A = np.array([r['d2A'] for r in rows])
    aR = np.array([r['aR'] for r in rows])
    n = len(rows)
    thick = ratio >= 10; thin = ratio < 10
    print(f"n={n}: thick={thick.sum()}, thin={thin.sum()}")
    print(f"δ/R: {dR.min():.5f}~{dR.max():.5f}")
    print(f"a/R: {aR.min():.4f}~{aR.max():.4f}")

    # ═══════════════════════════════════════════
    # 0. Correlation (δ/R vs σ_el, T/d)
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("0. Correlation matrix (log space)")
    print("="*80)
    log_s = np.log(sel)
    vars_check = [('φ', pa), ('CN', cn), ('τ', tau), ('cov', cov),
                  ('T/d', ratio), ('por', por), ('δ/R', dR), ('a/R', aR)]
    print(f"  {'':10s}", end='')
    for name, _ in vars_check:
        print(f" {name:>7s}", end='')
    print(f" {'σ_el':>7s}")
    for name1, arr1 in vars_check:
        print(f"  {name1:10s}", end='')
        for name2, arr2 in vars_check:
            r = np.corrcoef(np.log(arr1), np.log(arr2))[0,1]
            print(f" {r:+.3f} ", end='')
        r = np.corrcoef(np.log(arr1), log_s)[0,1]
        print(f" {r:+.3f}")

    # ═══════════════════════════════════════════
    # 1. φ^a × (δ/R)^b — 2변수만
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("1. φ^a × (δ/R)^b (2변수)")
    print("="*80)
    for a in np.arange(2, 6, 0.5):
        for b in np.arange(0.5, 4, 0.25):
            rhs = SAM * pa**a * dR**b
            C = fitC(sel, rhs)
            if C is None: continue
            r2 = r2l(sel, C*rhs)
            if r2 > 0.82:
                cv = loocv_C(sel, rhs)
                print(f"  φ^{a}×(δ/R)^{b}: R²={r2:.4f} LOOCV={cv:.4f}")

    # ═══════════════════════════════════════════
    # 2. φ^a × (δ/R)^b × X^c — 3변수
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("2. φ^a × (δ/R)^b × X^c (3변수)")
    print("="*80)
    c2 = []
    for a in [3, 3.5, 4, 4.5]:
        for b in [0.5, 0.75, 1, 1.5, 2]:
            for xn, xa, xes in [('CN', cn, [0.5, 1, 1.5]),
                                 ('√τ', tau**0.5, [1]),
                                 ('τ', tau, [0.5, 0.75]),
                                 ('cov', cov, [0.5, 1]),
                                 ('el_perc', el_perc, [0.5, 1])]:
                for xe in xes:
                    rhs = SAM * pa**a * dR**b * xa**xe
                    C = fitC(sel, rhs)
                    if C is None: continue
                    r2 = r2l(sel, C*rhs)
                    if r2 > 0.84:
                        c2.append({'r2': r2, 'rhs': rhs,
                                   'label': f'φ^{a}×(δ/R)^{b}×{xn}^{xe}'})
    c2.sort(key=lambda x: -x['r2'])
    for i, r in enumerate(c2[:20]):
        cv = loocv_C(sel, r['rhs']) if i < 10 else 0
        cs = f" LOOCV={cv:.4f}" if cv else ""
        print(f"  #{i+1} R²={r['r2']:.4f}{cs}  {r['label']}")

    # ═══════════════════════════════════════════
    # 3. φ^a × (δ/R)^b × CN^c × τ^d — 4변수
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("3. φ^a × (δ/R)^b × CN^c × τ^d (4변수, 깔끔)")
    print("="*80)
    c3 = []
    for a in [3, 3.5, 4, 4.5, 5]:
        for b in [0.5, 0.75, 1, 1.25, 1.5, 2]:
            for c in [0, 0.5, 1]:
                for d in [0, 0.25, 0.5, 0.75]:
                    rhs = SAM * pa**a * dR**b * cn**c * tau**d
                    C = fitC(sel, rhs)
                    if C is None: continue
                    r2 = r2l(sel, C*rhs)
                    if r2 > 0.85:
                        cn_s = f'×CN^{c}' if c else ''
                        tau_s = f'×τ^{d}' if d else ''
                        c3.append({'r2': r2, 'rhs': rhs,
                                   'label': f'φ^{a}×(δ/R)^{b}{cn_s}{tau_s}'})
    c3.sort(key=lambda x: -x['r2'])
    for i, r in enumerate(c3[:20]):
        cv = loocv_C(sel, r['rhs']) if i < 10 else 0
        cs = f" LOOCV={cv:.4f}" if cv else ""
        print(f"  #{i+1} R²={r['r2']:.4f}{cs}  {r['label']}")

    # ═══════════════════════════════════════════
    # 4. FORM X style: ⁴√[φ^a × (δ/R)^b × ...]
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("4. ⁴√ FORM X style")
    print("="*80)
    formx = [
        ('⁴√[φ¹⁶×(δ/R)⁴×CN⁴×τ²]', (pa**16 * dR**4 * cn**4 * tau**2)**0.25),
        ('⁴√[φ¹⁶×(δ/R)⁴×τ²×cov²]', (pa**16 * dR**4 * tau**2 * cov**2)**0.25),
        ('⁴√[φ¹⁶×(δ/R)⁸×τ²]', (pa**16 * dR**8 * tau**2)**0.25),
        ('⁴√[φ²⁰×(δ/R)⁴×τ²]', (pa**20 * dR**4 * tau**2)**0.25),
        ('⁴√[φ¹⁶×(δ/R)⁴×CN²×τ²]', (pa**16 * dR**4 * cn**2 * tau**2)**0.25),
        ('φ⁴×(δ/R)×√τ', pa**4 * dR * tau**0.5),
        ('φ⁴×(δ/R)×CN×√τ', pa**4 * dR * cn * tau**0.5),
        ('φ⁴×√(δ/R)×CN×√τ', pa**4 * dR**0.5 * cn * tau**0.5),
        ('φ⁴×(δ/R)^(3/4)×√τ', pa**4 * dR**0.75 * tau**0.5),
        ('φ⁴×(δ/R)^(3/4)×CN^(1/2)×√τ', pa**4 * dR**0.75 * cn**0.5 * tau**0.5),
        ('φ⁴×(δ/R)×cov×√τ', pa**4 * dR * cov * tau**0.5),
        ('(φ-0.2)³×(δ/R)×CN×√τ', np.clip(pa-0.2,0.001,None)**3 * dR * cn * tau**0.5),
        ('√[φ⁸×(δ/R)²×τ]', np.sqrt(pa**8 * dR**2 * tau)),
        ('√[φ⁸×(δ/R)²×CN²×τ]', np.sqrt(pa**8 * dR**2 * cn**2 * tau)),
        ('√[φ⁸×(δ/R)×CN²×τ]', np.sqrt(pa**8 * dR * cn**2 * tau)),
    ]
    print(f"  {'Formula':50s} {'R²':>7s} {'LOOCV':>7s}")
    print("  " + "-" * 66)
    for label, rhs_raw in sorted(formx, key=lambda x: -r2l(sel, fitC(sel, SAM*x[1])*(SAM*x[1])) if fitC(sel, SAM*x[1]) else -999):
        rhs = SAM * rhs_raw
        C = fitC(sel, rhs)
        if C is None: continue
        r2 = r2l(sel, C*rhs)
        cv = loocv_C(sel, rhs) if r2 > 0.82 else 0
        cs = f"{cv:7.4f}" if cv else "      -"
        flag = '★' if r2 > 0.88 else '●' if r2 > 0.84 else ' '
        print(f"  {flag}{label:49s} {r2:7.4f} {cs}  C={C:.6f}")

    # ═══════════════════════════════════════════
    # 5. Per-regime & per-case
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("5. Best → thick/thin/per-case")
    print("="*80)
    # Pick top from c3
    if c3:
        best = c3[0]
        C = fitC(sel, best['rhs']); pred = C * best['rhs']
        r2_tk = r2l(sel[thick], pred[thick]) if thick.sum()>=3 else 0
        r2_tn = r2l(sel[thin], pred[thin]) if thin.sum()>=3 else 0
        err_tk = np.mean(np.abs(sel[thick]-pred[thick])/sel[thick]*100)
        err_tn = np.mean(np.abs(sel[thin]-pred[thin])/sel[thin]*100)
        print(f"\n  BEST: {best['label']}  (R²={best['r2']:.4f})")
        print(f"  thick: R²={r2_tk:.4f}, |err|={err_tk:.0f}%")
        print(f"  thin:  R²={r2_tn:.4f}, |err|={err_tn:.0f}%")
        print(f"\n  Per-case (worst 10):")
        err = np.abs(sel - pred) / sel * 100
        for i in np.argsort(-err)[:10]:
            r = rows[i]
            regime = 'THICK' if ratio[i] >= 10 else 'thin'
            print(f"    {r['name']:35s} σ={sel[i]:.2f} pred={pred[i]:.2f} err={err[i]:.0f}% "
                  f"T/d={ratio[i]:.1f} δ/R={dR[i]:.4f} {regime}")

    # ═══════════════════════════════════════════
    # PHYSICAL MEANING
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("PHYSICAL MEANING")
    print("="*80)
    print("""
  σ_el = C × σ_AM × φ^a × (δ/R)^b × CN^c × τ^d

  각 항의 물리:
  ─────────────────────────────────────────────────
  φ_AM^a (a≈4):  AM volume fraction
    → percolation theory: σ ∝ (φ-φc)^t, t≈2 for 3D
    → φ^4 ≈ (φ-φc)^2 × φ^2 (percolation + dilution)

  (δ/R)^b (b≈1~2):  상대 침투 깊이
    → Hertz: δ/R = (3πP/(4E*))² → pressing 강도 반영
    → 접촉 면적 A = πRδ, 접촉 반경 a = √(Rδ)
    → constriction resistance R_c = 1/(2σa) ∝ 1/√(Rδ) = 1/(R√(δ/R))
    → 높은 δ/R → 넓은 접촉 → 낮은 저항
    → thin electrode: pressing force 불균일 → δ/R 변동 큼

  CN^c (c≈0~1):  AM-AM coordination number
    → 전기적 path 수 (parallel resistors)
    → CN=0이면 percolation 불가

  τ^d (d≈0.5):  tortuosity
    → ionic τ가 electronic에도 영향? SE packing이 AM 배열을 결정
    → SE τ 높으면 → AM도 불규칙 → electronic path 비효율

  T/d 없는 이유:
    → T/d ∝ 입자 층수 → finite-size effect
    → 하지만 δ/R이 이미 thin/thick 차이를 흡수
    → thin: pressing 불균일 → δ/R 작음 → σ 낮음 (자동 반영!)
    → thick: pressing 균일 → δ/R 큼 → σ 높음
    """)

if __name__ == '__main__':
    main()
