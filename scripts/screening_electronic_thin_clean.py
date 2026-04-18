"""
Thin 깔끔 식: δ²/A 변형 + 대안 조합
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
            tau = max(m.get('tortuosity_recommended', m.get('tortuosity_mean', 0)), 0.1)
            cov = max(m.get('coverage_AM_P_mean', m.get('coverage_AM_S_mean', m.get('coverage_AM_mean', 20))), 0.1) / 100
            r_am = max(m.get('r_AM_P', 0), m.get('r_AM_S', 0))
            d_am = r_am * 2 if r_am > 0.1 else 5.0
            por = m.get('porosity', 0)
            el_perc = max(m.get('electronic_percolating_fraction', 0), 0.01)
            A = m.get('am_am_mean_area', 0)
            d = m.get('am_am_mean_delta', 0)
            hop = max(m.get('am_am_mean_hop', 0), 0.1)
            a_c = max(m.get('am_am_mean_contact_radius', 0), 0.01)
            F = max(m.get('am_am_mean_force', 0), 0.01)
            if pa <= 0 or cn <= 0 or T <= 0 or d <= 0 or A <= 0: continue
            R = d_am / 2
            rows.append({
                'sel': sel, 'pa': pa, 'cn': cn, 'tau': tau, 'cov': cov,
                'xi': T/d_am, 'por': max(por, 0.1), 'el_perc': el_perc,
                'd': d, 'A': A, 'hop': hop, 'a_c': a_c, 'F': F, 'R': R,
                # Derived ratios
                'd2A': d**2/A, 'dA': d/A, 'd3A2': d**3/A**2,
                'dR': d/R, 'aR': a_c/R, 'dsqA': d/np.sqrt(A),
                'Fd': F*d, 'FA': F/A,  # pressure
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

def main():
    rows = load_data()
    if not rows: print("데이터 없음!"); return
    V = {k: np.array([r[k] for r in rows]) for k in rows[0] if k != 'name'}
    n = len(rows)
    thin = V['xi'] < 10; sn = V['sel'][thin]; n_thin = thin.sum()
    T = {k: V[k][thin] for k in V}
    print(f"n={n}, thin={n_thin}")

    # ═══════════════════════════════════════════
    # 1. δ²/A vs 대안 contact quality ratios
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("1. Contact quality ratio 비교 (CN × Q^a / √(T/d))")
    print("="*80)
    ratios = [
        ('δ²/A', T['d2A']),
        ('δ/A', T['dA']),
        ('δ³/A²', T['d3A2']),
        ('δ/√A', T['dsqA']),
        ('δ/R', T['dR']),
        ('a/R', T['aR']),
        ('δ', T['d']),
        ('A', T['A']),
        ('a_c', T['a_c']),
        ('F×δ', T['Fd']),
        ('F/A (P)', T['FA']),
        ('hop', T['hop']),
        ('δ×hop', T['d']*T['hop']),
        ('δ²×hop/A', T['d2A']*T['hop']),
        ('a_c×δ', T['a_c']*T['d']),
    ]
    print(f"  {'Ratio Q':20s} {'Q^½':>7s} {'Q^¾':>7s} {'Q^1':>7s} {'Q^(3/2)':>7s} {'best_a':>7s} {'best_R²':>7s}")
    print("  " + "-" * 80)
    for name, arr in ratios:
        best_r2 = -999; best_a = 0
        results = {}
        for a in [0.5, 0.75, 1, 1.5]:
            rhs = SAM * T['cn'] * arr**a / T['xi']**0.5
            C = fitC(sn, rhs)
            if C is None: results[a] = 0; continue
            r2 = r2l(sn, C*rhs); results[a] = r2
            if r2 > best_r2: best_r2 = r2; best_a = a
        print(f"  {name:20s} {results.get(0.5,0):7.4f} {results.get(0.75,0):7.4f} "
              f"{results.get(1,0):7.4f} {results.get(1.5,0):7.4f} {best_a:7.2f} {best_r2:7.4f}")

    # ═══════════════════════════════════════════
    # 2. Top ratio + hop 조합
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("2. Top ratio × hop^b / (T/d)^c × CN^d")
    print("="*80)
    c2 = []
    for qn, qa in [('δ²/A', T['d2A']), ('δ/√A', T['dsqA']), ('δ×hop', T['d']*T['hop']),
                     ('a_c×δ', T['a_c']*T['d']), ('δ²×hop/A', T['d2A']*T['hop'])]:
        for a in [0.5, 0.625, 0.75, 1]:
            for b in [0, 0.25, 0.5]:
                for c in [0.25, 0.375, 0.5]:
                    for d in [0.75, 0.875, 1, 1.125]:
                        rhs = SAM * T['cn']**d * qa**a * T['hop']**b / T['xi']**c
                        C = fitC(sn, rhs)
                        if C is None: continue
                        r2 = r2l(sn, C*rhs)
                        if r2 > 0.92:
                            hs = f'×hop^{b}' if b else ''
                            c2.append({'r2':r2, 'rhs':rhs,
                                'label':f'CN^{d}×({qn})^{a}{hs}/(T/d)^{c}'})
    c2.sort(key=lambda x: -x['r2'])
    for i, r in enumerate(c2[:25]):
        cv = loocv_C(sn, r['rhs']) if i < 10 else 0
        cs = f" LOOCV={cv:.4f}" if cv else ""
        print(f"  #{i+1} R²={r['r2']:.4f}{cs}  {r['label']}")

    # ═══════════════════════════════════════════
    # 3. FORM X style 깔끔한 식 모음
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("3. 깔끔한 식 모음 (정수/반정수 지수만)")
    print("="*80)
    clean = [
        # δ²/A base
        ('CN×(δ²/A)^¾/√ξ', T['cn']*T['d2A']**0.75/T['xi']**0.5),
        ('CN×(δ²/A)^¾×√hop/√ξ', T['cn']*T['d2A']**0.75*T['hop']**0.5/T['xi']**0.5),
        ('CN×√(δ²/A)/√ξ', T['cn']*T['d2A']**0.5/T['xi']**0.5),
        ('CN×√(δ²/A)×√hop/√ξ', T['cn']*T['d2A']**0.5*T['hop']**0.5/T['xi']**0.5),
        ('CN×(δ²/A)/√ξ', T['cn']*T['d2A']/T['xi']**0.5),
        # ⁴√
        ('⁴√[CN⁴×(δ²/A)³/ξ²]', (T['cn']**4*T['d2A']**3/T['xi']**2)**0.25),
        ('⁴√[CN⁴×(δ²/A)³×hop²/ξ²]', (T['cn']**4*T['d2A']**3*T['hop']**2/T['xi']**2)**0.25),
        ('⁴√[CN⁴×(δ²/A)²×hop²/ξ²]', (T['cn']**4*T['d2A']**2*T['hop']**2/T['xi']**2)**0.25),
        # √
        ('√[CN²×(δ²/A)^(3/2)/ξ]', np.sqrt(T['cn']**2*T['d2A']**1.5/T['xi'])),
        ('√[CN²×(δ²/A)×hop/ξ]', np.sqrt(T['cn']**2*T['d2A']*T['hop']/T['xi'])),
        ('√[CN²×(δ²/A)^(3/2)×hop/ξ]', np.sqrt(T['cn']**2*T['d2A']**1.5*T['hop']/T['xi'])),
        # δ/√A base (alternative)
        ('CN×(δ/√A)^(3/2)/√ξ', T['cn']*(T['dsqA'])**1.5/T['xi']**0.5),
        ('CN×(δ/√A)/√ξ', T['cn']*T['dsqA']/T['xi']**0.5),
        ('CN×(δ/√A)×√hop/√ξ', T['cn']*T['dsqA']*T['hop']**0.5/T['xi']**0.5),
        # a_c×δ base
        ('CN×(aδ)^¾/√ξ', T['cn']*(T['a_c']*T['d'])**0.75/T['xi']**0.5),
        ('CN×√(aδ)×√hop/√ξ', T['cn']*np.sqrt(T['a_c']*T['d'])*T['hop']**0.5/T['xi']**0.5),
        # δ only (simpler)
        ('CN×δ/√ξ', T['cn']*T['d']/T['xi']**0.5),
        ('CN×δ×√hop/√ξ', T['cn']*T['d']*T['hop']**0.5/T['xi']**0.5),
        ('CN×√δ×hop/√ξ', T['cn']*T['d']**0.5*T['hop']/T['xi']**0.5),
        # hop×δ combined
        ('CN×(δ×hop)^½/√ξ', T['cn']*np.sqrt(T['d']*T['hop'])/T['xi']**0.5),
        ('CN×(δ×hop)^¾/√ξ', T['cn']*(T['d']*T['hop'])**0.75/T['xi']**0.5),
        ('CN×(δ²×hop/A)^½/√ξ', T['cn']*np.sqrt(T['d2A']*T['hop'])/T['xi']**0.5),
        ('CN×(δ²×hop/A)^¾/√ξ', T['cn']*(T['d2A']*T['hop'])**0.75/T['xi']**0.5),
    ]
    print(f"  {'Formula':45s} {'R²':>7s} {'LOOCV':>7s} {'|err|':>6s}")
    print("  " + "-" * 65)
    for label, rhs_raw in sorted(clean, key=lambda x: -r2l(sn, fitC(sn,SAM*x[1])*(SAM*x[1])) if fitC(sn,SAM*x[1]) else -999):
        rhs = SAM * rhs_raw
        C = fitC(sn, rhs)
        if C is None: continue
        r2 = r2l(sn, C*rhs); cv = loocv_C(sn, rhs)
        err = np.mean(np.abs(sn - C*rhs) / sn * 100)
        w20 = np.sum(np.abs(sn - C*rhs) / sn < 0.2)
        flag = '★' if r2 > 0.92 else '●' if r2 > 0.90 else ' '
        print(f"  {flag}{label:44s} {r2:7.4f} {cv:7.4f} {err:5.0f}% ({w20}/{n_thin})")

    # ═══════════════════════════════════════════
    # 4. 물리적 의미 정리
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("4. 물리적 의미")
    print("="*80)
    print("""
  δ²/A = (침투깊이)² / (접촉면적)
    - Hertz: δ²/A = δ²/(πRδ) = δ/(πR) ∝ 상대침투
    - DEM: A와 δ가 독립 → δ²/A가 δ/R보다 정보량 많음
    - 물리: 깊이 파고들수록 & 면적 작을수록 → 접촉 집중도

  δ/√A = 침투/√면적
    - √A ∝ a (접촉반경), so δ/√A ∝ δ/a
    - 물리: 접촉반경 대비 침투깊이 = 접촉 종횡비

  a_c×δ = 접촉반경 × 침투깊이
    - Hertz: a²=Rδ → a×δ = √(Rδ)×δ = δ^(3/2)×√R
    - 물리: 접촉의 크기 scale

  δ×hop = 침투 × 간격
    - 물리: 접촉 강도 × network 밀도
    """)

if __name__ == '__main__':
    main()
