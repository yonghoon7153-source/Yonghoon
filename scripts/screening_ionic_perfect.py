"""
Ionic FORM X 무결점: R²=0.96+ 달성
===================================
현재 FORM X: σ = C × σ_grain × (φ-φc)^¾ × CN × √cov / √τ
R²≈0.94 (n=47)

개선 시도:
  1. φ_c 최적화 (0.10~0.25)
  2. 지수 미세조정 (α=0.5~1.0 for φ-φc, etc)
  3. Path conductance, bottleneck 추가
  4. GB density 추가
  5. 새 변수 조합
  6. cov 대신 bottleneck?
"""
import json, os, numpy as np, warnings
from pathlib import Path
from itertools import combinations
warnings.filterwarnings('ignore')
WEBAPP = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'webapp')
SGRAIN = 3.0  # σ_grain mS/cm

def load_data():
    rows = []
    for base in [Path(WEBAPP)/'results', Path(WEBAPP)/'archive']:
        if not base.is_dir(): continue
        for mp in base.rglob('full_metrics.json'):
            try:
                with open(mp) as f: m = json.load(f)
            except: continue
            sion = m.get('sigma_full_mScm', 0)
            if not sion or sion < 0.01: continue
            ps = m.get('phi_se', 0)
            cn = m.get('se_se_cn', 0)
            tau = m.get('tortuosity_recommended', m.get('tortuosity_mean', 0))
            cov = max(m.get('coverage_AM_P_mean',
                          m.get('coverage_AM_S_mean',
                               m.get('coverage_AM_mean', 20))), 0.1) / 100
            gd = m.get('gb_density_mean', 0)
            gp = max(m.get('path_conductance_mean', 0), 1e-6)
            bn = max(m.get('path_hop_area_min_mean', 0), 1e-4)  # bottleneck
            ha = max(m.get('path_hop_area_mean', 0), 1e-4)
            perc = m.get('percolation_pct', 0) / 100
            r_am = max(m.get('r_AM_P', 0), m.get('r_AM_S', 0))
            d_am = r_am * 2 if r_am > 0.1 else 5.0
            T = m.get('thickness_um', 0)
            por = m.get('porosity', 0)
            pa = max(m.get('phi_am', 0), 0.01)
            if ps <= 0 or cn <= 0 or tau <= 0 or T <= 0: continue
            rows.append({
                'sion': sion, 'ps': ps, 'cn': cn, 'tau': max(tau, 0.1),
                'cov': cov, 'gd': max(gd, 0.01), 'gp': gp, 'bn': bn, 'ha': ha,
                'perc': max(perc, 0.01), 'ratio': T / d_am,
                'por': max(por, 0.1), 'pa': pa,
                'name': mp.parent.name
            })
    seen = set(); u = []
    for r in rows:
        k = f"{r['ps']:.4f}_{r['ratio']:.1f}"
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

    sion = np.array([r['sion'] for r in rows])
    ps = np.array([r['ps'] for r in rows])
    cn = np.array([r['cn'] for r in rows])
    tau = np.array([r['tau'] for r in rows])
    cov = np.array([r['cov'] for r in rows])
    gd = np.array([r['gd'] for r in rows])
    gp = np.array([r['gp'] for r in rows])
    bn = np.array([r['bn'] for r in rows])
    ha = np.array([r['ha'] for r in rows])
    perc = np.array([r['perc'] for r in rows])
    ratio = np.array([r['ratio'] for r in rows])
    por = np.array([r['por'] for r in rows])
    pa = np.array([r['pa'] for r in rows])
    n = len(rows)
    print(f"n={n}")
    print(f"BN (bottleneck): {bn.min():.4f}~{bn.max():.4f} µm²")
    print(f"Gc (path cond): {gp.min():.6f}~{gp.max():.6f} µm²")
    print(f"Gd: {gd.min():.2f}~{gd.max():.2f}")

    # ═══════════════════════════════════════════
    # 1. FORM X baseline + φ_c scan
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("1. FORM X φ_c scan")
    print("="*80)
    for phi_c in np.arange(0.05, 0.30, 0.01):
        phi_ex = np.clip(ps - phi_c, 0.001, None)
        rhs = SGRAIN * phi_ex**0.75 * cn * cov**0.5 / tau**0.5
        C = fitC(sion, rhs)
        if C is None: continue
        r2 = r2l(sion, C * rhs)
        cv = loocv_C(sion, rhs)
        flag = '★' if r2 > 0.95 else '●' if r2 > 0.93 else ' '
        print(f"  {flag}φ_c={phi_c:.2f}: R²={r2:.4f} LOOCV={cv:.4f} C={C:.4f}")

    # ═══════════════════════════════════════════
    # 2. 지수 미세조정
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("2. 지수 미세조정 (φ-φc)^a × CN^b × cov^c / τ^d")
    print("="*80)
    results = []
    for phi_c in [0.15, 0.16, 0.17, 0.18, 0.19, 0.20]:
        phi_ex = np.clip(ps - phi_c, 0.001, None)
        for a in [0.5, 0.625, 0.75, 0.875, 1.0]:
            for b in [0.75, 1.0, 1.25]:
                for c in [0.25, 0.5, 0.75]:
                    for d in [0.25, 0.5, 0.75]:
                        rhs = SGRAIN * phi_ex**a * cn**b * cov**c / tau**d
                        C = fitC(sion, rhs)
                        if C is None: continue
                        r2 = r2l(sion, C * rhs)
                        if r2 > 0.93:
                            results.append({'r2': r2, 'rhs': rhs, 'phi_c': phi_c,
                                'label': f'φc={phi_c} (φ-φc)^{a}×CN^{b}×cov^{c}/τ^{d}'})
    results.sort(key=lambda x: -x['r2'])
    for i, r in enumerate(results[:20]):
        cv = loocv_C(sion, r['rhs']) if i < 10 else 0
        cs = f" LOOCV={cv:.4f}" if cv else ""
        print(f"  #{i+1} R²={r['r2']:.4f}{cs}  {r['label']}")

    # ═══════════════════════════════════════════
    # 3. Path conductance / bottleneck 추가
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("3. FORM X + path metrics")
    print("="*80)
    phi_ex18 = np.clip(ps - 0.18, 0.001, None)
    for extra_name, extra_arr in [
        ('×Gc^0.25', gp**0.25), ('×Gc^0.5', gp**0.5), ('×Gc^0.75', gp**0.75),
        ('×BN^0.25', bn**0.25), ('×BN^0.5', bn**0.5),
        ('×Gd^-0.25', gd**-0.25), ('×Gd^-0.5', gd**-0.5),
        ('×HA^0.25', ha**0.25), ('×HA^0.5', ha**0.5),
        ('×perc^0.5', perc**0.5), ('×perc', perc),
    ]:
        rhs = SGRAIN * phi_ex18**0.75 * cn * cov**0.5 / tau**0.5 * extra_arr
        C = fitC(sion, rhs)
        if C is None: continue
        r2 = r2l(sion, C * rhs)
        cv = loocv_C(sion, rhs)
        flag = '★' if r2 > 0.95 else '●' if r2 > 0.93 else ' '
        print(f"  {flag}FORM X {extra_name}: R²={r2:.4f} LOOCV={cv:.4f}")

    # ═══════════════════════════════════════════
    # 4. cov를 bottleneck으로 대체?
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("4. cov → bottleneck/path_cond 대체")
    print("="*80)
    for phi_c in [0.15, 0.18, 0.20]:
        phi_ex = np.clip(ps - phi_c, 0.001, None)
        replacements = [
            ('(φ-φc)^¾×CN×√BN/√τ', phi_ex**0.75 * cn * bn**0.5 / tau**0.5),
            ('(φ-φc)^¾×CN×√Gc/√τ', phi_ex**0.75 * cn * gp**0.5 / tau**0.5),
            ('(φ-φc)^¾×CN×√HA/√τ', phi_ex**0.75 * cn * ha**0.5 / tau**0.5),
            ('(φ-φc)^¾×CN×BN^¼/√τ', phi_ex**0.75 * cn * bn**0.25 / tau**0.5),
            ('(φ-φc)^¾×CN×Gc^¼/√τ', phi_ex**0.75 * cn * gp**0.25 / tau**0.5),
            ('(φ-φc)^¾×CN×√cov×√BN/√τ', phi_ex**0.75 * cn * cov**0.5 * bn**0.5 / tau**0.5),
            ('(φ-φc)^¾×CN×√(cov×BN)/√τ', phi_ex**0.75 * cn * np.sqrt(cov*bn) / tau**0.5),
            ('(φ-φc)^¾×CN/√τ (no cov)', phi_ex**0.75 * cn / tau**0.5),
        ]
        for label, rhs_raw in replacements:
            rhs = SGRAIN * rhs_raw
            C = fitC(sion, rhs)
            if C is None: continue
            r2 = r2l(sion, C * rhs)
            cv = loocv_C(sion, rhs)
            flag = '★' if r2 > 0.95 else '●' if r2 > 0.93 else ' '
            print(f"  {flag}φc={phi_c} {label}: R²={r2:.4f} LOOCV={cv:.4f}")

    # ═══════════════════════════════════════════
    # 5. Free regression (all variables)
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("5. Free regression")
    print("="*80)
    phi_ex18 = np.clip(ps - 0.18, 0.001, None)
    var_pool = {
        '(φ-φc)': np.log(phi_ex18), 'CN': np.log(cn), 'τ': np.log(tau),
        'cov': np.log(cov), 'Gd': np.log(gd), 'Gc': np.log(gp),
        'BN': np.log(bn), 'HA': np.log(ha), 'perc': np.log(perc),
        'por': np.log(por),
    }
    log_sion = np.log(sion)
    results_free = []
    for nv in [3, 4, 5]:
        for vn in combinations(var_pool.keys(), nv):
            X = np.column_stack([var_pool[v] for v in vn] + [np.ones(n)])
            try:
                coefs, _, _, _ = np.linalg.lstsq(X, log_sion, rcond=None)
                pred = np.exp(X @ coefs)
                r2 = r2l(sion, pred)
                if r2 > 0.93:
                    results_free.append({'r2': r2, 'pred': pred, 'coefs': coefs, 'vnames': vn})
            except: pass
    results_free.sort(key=lambda x: -x['r2'])
    print(f"  Top 20 (n={n}):")
    for i, r in enumerate(results_free[:20]):
        cv = loocv_C(sion, r['pred']) if i < 10 else 0
        cs = f" LOOCV={cv:.4f}" if cv else ""
        terms = ' × '.join(f"{v}^{r['coefs'][j]:.2f}" for j, v in enumerate(r['vnames']))
        clean = all(abs(r['coefs'][j] - round(r['coefs'][j]*4)/4) < 0.08 for j in range(len(r['vnames'])))
        flag = '★' if clean else ' '
        print(f"    #{i+1}{flag} R²={r['r2']:.4f}{cs}  ({len(r['vnames'])}v) {terms}")

    # ═══════════════════════════════════════════
    # 6. FORM X 스타일 (⁴√)
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("6. ⁴√ 스타일")
    print("="*80)
    formx_cands = [
        ('⁴√[(φ-φc)³×CN⁴×cov²/τ²]', (phi_ex18**3 * cn**4 * cov**2 / tau**2)**0.25),
        ('⁴√[(φ-φc)³×CN⁴×BN²/τ²]', (phi_ex18**3 * cn**4 * bn**2 / tau**2)**0.25),
        ('⁴√[(φ-φc)³×CN⁴×Gc²/τ²]', (phi_ex18**3 * cn**4 * gp**2 / tau**2)**0.25),
        ('⁴√[(φ-φc)³×CN⁴×cov×BN/τ²]', (phi_ex18**3 * cn**4 * cov * bn / tau**2)**0.25),
        ('⁴√[(φ-φc)³×CN⁴×(cov×BN)²/τ²]', (phi_ex18**3 * cn**4 * (cov*bn)**2 / tau**2)**0.25),
        ('⁴√[(φ-φc)³×CN⁴×cov²/τ]', (phi_ex18**3 * cn**4 * cov**2 / tau)**0.25),
        ('⁴√[(φ-φc)⁴×CN⁴×cov²/τ²]', (phi_ex18**4 * cn**4 * cov**2 / tau**2)**0.25),
        ('⁴√[(φ-φc)²×CN⁴×cov²/τ²]', (phi_ex18**2 * cn**4 * cov**2 / tau**2)**0.25),
        ('⁴√[(φ-φc)³×CN⁴×HA²/τ²]', (phi_ex18**3 * cn**4 * ha**2 / tau**2)**0.25),
    ]
    print(f"  {'Formula':50s} {'R²':>7s} {'LOOCV':>7s}")
    print("  " + "-" * 66)
    for label, rhs_raw in sorted(formx_cands, key=lambda x: -r2l(sion, fitC(sion, SGRAIN*x[1])*(SGRAIN*x[1])) if fitC(sion, SGRAIN*x[1]) else -999):
        rhs = SGRAIN * rhs_raw
        C = fitC(sion, rhs)
        if C is None: continue
        r2 = r2l(sion, C * rhs)
        cv = loocv_C(sion, rhs)
        flag = '★' if r2 > 0.95 else '●' if r2 > 0.93 else ' '
        print(f"  {flag}{label:49s} {r2:7.4f} {cv:7.4f}  C={C:.4f}")

    # ═══════════════════════════════════════════
    # 7. Per-case error (current FORM X)
    # ═══════════════════════════════════════════
    rhs_formx = SGRAIN * phi_ex18**0.75 * cn * cov**0.5 / tau**0.5
    C_formx = fitC(sion, rhs_formx)
    pred_formx = C_formx * rhs_formx
    err = np.abs(sion - pred_formx) / sion * 100
    print(f"\n{'='*80}")
    print(f"7. FORM X per-case (φc=0.18, C={C_formx:.4f})")
    print("="*80)
    for i in np.argsort(-err)[:10]:
        r = rows[i]
        print(f"  {r['name']:35s} σ={sion[i]:.3f} pred={pred_formx[i]:.3f} "
              f"err={err[i]:.0f}% φ_SE={r['ps']:.3f} CN={r['cn']:.1f} τ={r['tau']:.2f} cov={r['cov']:.3f}")

if __name__ == '__main__':
    main()
