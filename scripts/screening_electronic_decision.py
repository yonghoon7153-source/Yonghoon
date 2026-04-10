"""
Electronic DECISION: thick + thin 최종 결정
============================================
0.125 단위 미세조정, ablation, per-case, LOOCV
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
                'cov': cov, 'xi': T/d_am, 'por': max(por, 0.1),
                'd2a': am_delta**2/am_area, 'hop': max(am_hop, 0.1),
                'd': am_delta, 'A': am_area, 'name': mp.parent.name
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
    n = len(rows); names = [r['name'] for r in rows]
    thick = V['xi'] >= 10; thin = V['xi'] < 10
    st = V['sel'][thick]; sn = V['sel'][thin]
    print(f"n={n}: thick={thick.sum()}, thin={thin.sum()}")

    # ╔═══════════════════════════════════════════╗
    # ║  THICK: φ^a × CN^b × cov^c × τ^d        ║
    # ╚═══════════════════════════════════════════╝
    print(f"\n{'='*80}")
    print(f"THICK 미세조정 (n={thick.sum()}, step=0.125)")
    print("="*80)
    T = {k: V[k][thick] for k in V}
    ct = []
    for a in np.arange(3.5, 5.0, 0.125):
        for b in np.arange(1.0, 2.125, 0.125):
            for c in np.arange(0.5, 1.5, 0.125):
                for d in np.arange(0, 1.0, 0.125):
                    rhs = SAM * T['pa']**a * T['cn']**b * T['cov']**c * T['tau']**d
                    C = fitC(st, rhs)
                    if C is None: continue
                    r2 = r2l(st, C*rhs)
                    if r2 > 0.965:
                        ct.append({'r2':r2,'rhs':rhs,'a':a,'b':b,'c':c,'d':d,
                                   'label':f'φ^{a}×CN^{b}×cov^{c}×τ^{d}'})
    ct.sort(key=lambda x: -x['r2'])
    print(f"\n  Top 15:")
    for i, r in enumerate(ct[:15]):
        cv = loocv_C(st, r['rhs']) if i < 10 else 0
        cs = f" LOOCV={cv:.4f}" if cv else ""
        print(f"  #{i+1} R²={r['r2']:.4f}{cs}  {r['label']}  C={fitC(st,r['rhs']):.4f}")

    # Thick ablation
    if ct:
        best = ct[0]
        print(f"\n  Ablation (from #{1}: {best['label']}):")
        base_r2 = best['r2']
        ablations = [
            ('Remove τ', SAM*T['pa']**best['a']*T['cn']**best['b']*T['cov']**best['c']),
            ('Remove cov', SAM*T['pa']**best['a']*T['cn']**best['b']*T['tau']**best['d']),
            ('Remove CN', SAM*T['pa']**best['a']*T['cov']**best['c']*T['tau']**best['d']),
            ('Remove φ', SAM*T['cn']**best['b']*T['cov']**best['c']*T['tau']**best['d']),
        ]
        for label, rhs in ablations:
            C = fitC(st, rhs); r2 = r2l(st, C*rhs) if C else -99
            print(f"    {label:15s}: R²={r2:.4f} (Δ={r2-base_r2:+.4f})")

    # Thick 깔끔 후보
    print(f"\n  깔끔한 지수 후보:")
    clean_thick = [
        ('φ⁴×CN^(3/2)×cov×√τ', T['pa']**4*T['cn']**1.5*T['cov']*T['tau']**0.5),
        ('φ⁴×CN^(3/2)×cov', T['pa']**4*T['cn']**1.5*T['cov']),
        ('φ⁴×CN^(3/2)×cov^¾×τ^½', T['pa']**4*T['cn']**1.5*T['cov']**0.75*T['tau']**0.5),
        ('φ⁴×CN^(3/2)×cov^¾×τ^¾', T['pa']**4*T['cn']**1.5*T['cov']**0.75*T['tau']**0.75),
        ('φ⁴×CN^(5/4)×cov×τ^½', T['pa']**4*T['cn']**1.25*T['cov']*T['tau']**0.5),
        ('φ⁴×CN^(5/4)×cov^¾×τ^½', T['pa']**4*T['cn']**1.25*T['cov']**0.75*T['tau']**0.5),
        ('φ⁴×CN^(7/4)×cov×τ^½', T['pa']**4*T['cn']**1.75*T['cov']*T['tau']**0.5),
        ('⁴√[φ¹⁶×CN⁶×cov⁴×τ²]', (T['pa']**16*T['cn']**6*T['cov']**4*T['tau']**2)**0.25),
        ('⁴√[φ¹⁶×CN⁶×cov³×τ²]', (T['pa']**16*T['cn']**6*T['cov']**3*T['tau']**2)**0.25),
        ('⁴√[φ¹⁶×CN⁶×cov³×τ³]', (T['pa']**16*T['cn']**6*T['cov']**3*T['tau']**3)**0.25),
        ('⁴√[φ¹⁶×CN⁵×cov⁴×τ²]', (T['pa']**16*T['cn']**5*T['cov']**4*T['tau']**2)**0.25),
    ]
    for label, rhs_raw in sorted(clean_thick, key=lambda x: -r2l(st,fitC(st,SAM*x[1])*(SAM*x[1])) if fitC(st,SAM*x[1]) else -999):
        rhs = SAM * rhs_raw
        C = fitC(st, rhs); r2 = r2l(st, C*rhs); cv = loocv_C(st, rhs)
        flag = '★' if r2 > 0.968 else '●' if r2 > 0.96 else ' '
        print(f"  {flag}{label:40s} R²={r2:.4f} LOOCV={cv:.4f} C={C:.4f}")

    # Thick per-case
    if ct:
        best = ct[0]; C = fitC(st, best['rhs']); pred = C * best['rhs']
        err = np.abs(st - pred) / st * 100
        print(f"\n  Per-case (best):")
        tk_idx = np.where(thick)[0]
        for j in np.argsort(-err):
            i = tk_idx[j]
            print(f"    {names[i]:35s} σ={st[j]:.2f} pred={pred[j]:.2f} err={err[j]:.0f}% "
                  f"φ={V['pa'][i]:.3f} CN={V['cn'][i]:.2f} cov={V['cov'][i]:.3f} τ={V['tau'][i]:.2f}")

    # ╔═══════════════════════════════════════════╗
    # ║  THIN: CN^a × (δ²/A)^b × hop^c / (T/d)^d ║
    # ╚═══════════════════════════════════════════╝
    print(f"\n\n{'='*80}")
    print(f"THIN 미세조정 (n={thin.sum()}, step=0.125)")
    print("="*80)
    S = {k: V[k][thin] for k in V}
    cn2 = []
    for a in np.arange(0.5, 1.625, 0.125):
        for b in np.arange(0.5, 1.125, 0.125):
            for c in np.arange(0, 1.25, 0.125):
                for d in np.arange(0.25, 0.75, 0.125):
                    rhs = SAM * S['cn']**a * S['d2a']**b * S['hop']**c / S['xi']**d
                    C = fitC(sn, rhs)
                    if C is None: continue
                    r2 = r2l(sn, C*rhs)
                    if r2 > 0.90:
                        hs = f'×hop^{c}' if c else ''
                        cn2.append({'r2':r2,'rhs':rhs,'a':a,'b':b,'c':c,'d':d,
                                    'label':f'CN^{a}×(δ²/A)^{b}{hs}/(T/d)^{d}'})
    cn2.sort(key=lambda x: -x['r2'])
    print(f"\n  Top 20:")
    for i, r in enumerate(cn2[:20]):
        cv = loocv_C(sn, r['rhs']) if i < 10 else 0
        cs = f" LOOCV={cv:.4f}" if cv else ""
        print(f"  #{i+1} R²={r['r2']:.4f}{cs}  {r['label']}  C={fitC(sn,r['rhs']):.6f}")

    # Thin ablation
    if cn2:
        best = cn2[0]
        print(f"\n  Ablation (from #{1}: {best['label']}):")
        base_r2 = best['r2']
        parts = [
            ('Remove hop', SAM*S['cn']**best['a']*S['d2a']**best['b']/S['xi']**best['d']),
            ('Remove (T/d)', SAM*S['cn']**best['a']*S['d2a']**best['b']*S['hop']**best['c']),
            ('Remove (δ²/A)', SAM*S['cn']**best['a']*S['hop']**best['c']/S['xi']**best['d']),
            ('Remove CN', SAM*S['d2a']**best['b']*S['hop']**best['c']/S['xi']**best['d']),
        ]
        for label, rhs in parts:
            C = fitC(sn, rhs); r2 = r2l(sn, C*rhs) if C else -99
            print(f"    {label:15s}: R²={r2:.4f} (Δ={r2-base_r2:+.4f})")

    # Thin 깔끔 후보
    print(f"\n  깔끔한 지수 후보:")
    clean_thin = [
        ('CN×(δ²/A)^¾/√(T/d)', S['cn']*S['d2a']**0.75/S['xi']**0.5),
        ('CN×(δ²/A)^½/√(T/d)', S['cn']*S['d2a']**0.5/S['xi']**0.5),
        ('CN×(δ²/A)/√(T/d)', S['cn']*S['d2a']/S['xi']**0.5),
        ('CN×(δ²/A)^¾×√hop/√(T/d)', S['cn']*S['d2a']**0.75*S['hop']**0.5/S['xi']**0.5),
        ('CN×(δ²/A)^½×√hop/√(T/d)', S['cn']*S['d2a']**0.5*S['hop']**0.5/S['xi']**0.5),
        ('CN×(δ²/A)^(5/8)×√hop/√(T/d)', S['cn']*S['d2a']**0.625*S['hop']**0.5/S['xi']**0.5),
        ('⁴√[CN⁴×(δ²/A)³/(T/d)²]', (S['cn']**4*S['d2a']**3/S['xi']**2)**0.25),
        ('⁴√[CN⁴×(δ²/A)³×hop²/(T/d)²]', (S['cn']**4*S['d2a']**3*S['hop']**2/S['xi']**2)**0.25),
        ('√[CN²×(δ²/A)^(3/2)/T/d]', np.sqrt(S['cn']**2*S['d2a']**1.5/S['xi'])),
        ('√[CN²×(δ²/A)×hop/(T/d)]', np.sqrt(S['cn']**2*S['d2a']*S['hop']/S['xi'])),
        # without hop
        ('CN^¾×(δ²/A)^¾/√(T/d)', S['cn']**0.75*S['d2a']**0.75/S['xi']**0.5),
        ('CN^(5/4)×(δ²/A)^¾/√(T/d)', S['cn']**1.25*S['d2a']**0.75/S['xi']**0.5),
        # without T/d
        ('CN×(δ²/A)^¾×√hop', S['cn']*S['d2a']**0.75*S['hop']**0.5),
        ('CN×(δ²/A)^¾', S['cn']*S['d2a']**0.75),
    ]
    for label, rhs_raw in sorted(clean_thin, key=lambda x: -r2l(sn,fitC(sn,SAM*x[1])*(SAM*x[1])) if fitC(sn,SAM*x[1]) else -999):
        rhs = SAM * rhs_raw
        C = fitC(sn, rhs); r2 = r2l(sn, C*rhs); cv = loocv_C(sn, rhs)
        flag = '★' if r2 > 0.92 else '●' if r2 > 0.90 else ' '
        print(f"  {flag}{label:45s} R²={r2:.4f} LOOCV={cv:.4f} C={C:.6f}")

    # Thin per-case
    if cn2:
        best = cn2[0]; C = fitC(sn, best['rhs']); pred = C * best['rhs']
        err = np.abs(sn - pred) / sn * 100
        w20 = np.sum(err < 20)
        print(f"\n  Per-case (best, {w20}/{thin.sum()} within 20%):")
        tn_idx = np.where(thin)[0]
        for j in np.argsort(-err):
            i = tn_idx[j]
            print(f"    {names[i]:35s} σ={sn[j]:.2f} pred={pred[j]:.2f} err={err[j]:.0f}% "
                  f"T/d={V['xi'][i]:.1f} CN={V['cn'][i]:.2f} δ²/A={V['d2a'][i]:.5f} hop={V['hop'][i]:.1f}")

    # ╔═══════════════════════════════════════════╗
    # ║  FINAL DECISION                           ║
    # ╚═══════════════════════════════════════════╝
    print(f"\n\n{'='*80}")
    print("FINAL DECISION")
    print("="*80)
    # Thick
    rhs_tk = SAM * V['pa'][thick]**4 * V['cn'][thick]**1.5 * V['cov'][thick] * V['tau'][thick]**0.5
    C_tk = fitC(st, rhs_tk); cv_tk = loocv_C(st, rhs_tk)
    err_tk = np.mean(np.abs(st - C_tk*rhs_tk) / st * 100)
    w20_tk = np.sum(np.abs(st - C_tk*rhs_tk) / st < 0.2)
    # Thin (clean)
    rhs_tn = SAM * V['cn'][thin] * V['d2a'][thin]**0.75 / V['xi'][thin]**0.5
    C_tn = fitC(sn, rhs_tn); cv_tn = loocv_C(sn, rhs_tn)
    err_tn = np.mean(np.abs(sn - C_tn*rhs_tn) / sn * 100)
    w20_tn = np.sum(np.abs(sn - C_tn*rhs_tn) / sn < 0.2)
    # Thin (with hop)
    if cn2:
        best_tn = cn2[0]
        C_tn2 = fitC(sn, best_tn['rhs']); cv_tn2 = loocv_C(sn, best_tn['rhs'])
        err_tn2 = np.mean(np.abs(sn - C_tn2*best_tn['rhs']) / sn * 100)
        w20_tn2 = np.sum(np.abs(sn - C_tn2*best_tn['rhs']) / sn < 0.2)

    print(f"""
  ┌────────────────────────────────────────────────────────────────┐
  │  THICK (T/d ≥ 10, n={thick.sum()})                               │
  │                                                                │
  │  σ_el = {C_tk:.4f} × σ_AM × φ⁴ × CN^(3/2) × cov × √τ         │
  │  = {C_tk:.4f} × σ_AM × ⁴√[φ¹⁶ × CN⁶ × cov⁴ × τ²]            │
  │                                                                │
  │  R² = {r2l(st,C_tk*rhs_tk):.4f}  LOOCV = {cv_tk:.4f}                          │
  │  |err| = {err_tk:.1f}%  within 20%: {w20_tk}/{thick.sum()}                          │
  │                                                                │
  │  물리: φ⁴ = percolation+dilution                               │
  │       CN^(3/2) = network connectivity                         │
  │       cov = AM-SE interface → AM 구조 안정화                   │
  │       √τ = SE packing → AM 배열 영향                          │
  └────────────────────────────────────────────────────────────────┘

  ┌────────────────────────────────────────────────────────────────┐
  │  THIN (T/d < 10, n={thin.sum()})                                 │
  │                                                                │
  │  OPTION A (깔끔):                                              │
  │  σ_el = {C_tn:.4f} × σ_AM × CN × (δ²/A)^¾ / √(T/d)          │
  │  R² = {r2l(sn,C_tn*rhs_tn):.4f}  LOOCV = {cv_tn:.4f}  |err| = {err_tn:.0f}%  {w20_tn}/{thin.sum()}          │
  │                                                                │""")
    if cn2:
        print(f"""  │  OPTION B (정밀):                                              │
  │  σ_el = {C_tn2:.6f} × σ_AM × {best_tn['label'][:45]:45s}│
  │  R² = {best_tn['r2']:.4f}  LOOCV = {cv_tn2:.4f}  |err| = {err_tn2:.0f}%  {w20_tn2}/{thin.sum()}          │
  │                                                                │""")
    print(f"""  │  물리: CN = AM-AM 접촉 수                                    │
  │       (δ²/A)^¾ = contact quality (침투²/면적)                 │
  │       1/√(T/d) = finite-size correction                       │
  │       hop = AM 입자 간 거리 (network spacing)                 │
  └────────────────────────────────────────────────────────────────┘
    """)

if __name__ == '__main__':
    main()
