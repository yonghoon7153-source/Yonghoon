"""
Thin 집중 분석: 왜 R²=0.77이 한계인가?
residual 패턴 → 빠진 물리 찾기
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
            ps = m.get('phi_se', 0)
            sion = m.get('sigma_full_mScm', 0)
            am_area = m.get('am_am_mean_area', 0)
            am_delta = m.get('am_am_mean_delta', 0)
            am_hop = m.get('am_am_mean_hop', 0)
            am_pres = m.get('am_am_mean_pressure', 0)
            am_cn_std = m.get('am_am_cn_std', 0)
            if pa <= 0 or cn <= 0 or T <= 0 or am_delta <= 0 or am_area <= 0: continue
            rows.append({
                'sel': sel, 'pa': pa, 'cn': cn, 'tau': max(tau, 0.1),
                'cov': cov, 'ratio': T/d_am, 'por': max(por, 0.1),
                'el_perc': el_perc, 'ps': max(ps, 0.01), 'sion': max(sion, 0.001),
                'A': am_area, 'delta': am_delta, 'hop': max(am_hop, 0.1),
                'P': max(am_pres, 0.01), 'cn_std': max(am_cn_std, 0.01),
                'd_am': d_am, 'name': mp.parent.name
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
    xi = np.array([r['ratio'] for r in rows])
    por = np.array([r['por'] for r in rows])
    el_perc = np.array([r['el_perc'] for r in rows])
    ps = np.array([r['ps'] for r in rows])
    sion = np.array([r['sion'] for r in rows])
    A = np.array([r['A'] for r in rows])
    d = np.array([r['delta'] for r in rows])
    hop = np.array([r['hop'] for r in rows])
    P = np.array([r['P'] for r in rows])
    cn_std = np.array([r['cn_std'] for r in rows])
    d_am = np.array([r['d_am'] for r in rows])
    n = len(rows)
    thick = xi >= 10; thin = xi < 10
    print(f"n={n}: thick={thick.sum()}, thin={thin.sum()}")

    # ═══════════════════════════════════════════
    # 1. Base formula residual 분석 (thin만)
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("1. Residual 분석: φ^3.5×(T/d)^(-¼)×δ^(9/4)/A에서 thin residual")
    print("="*80)
    rhs_base = SAM * pa**3.5 * xi**(-0.25) * d**2.25 / A
    C_base = fitC(sel, rhs_base)
    pred_base = C_base * rhs_base
    resid = np.log(sel) - np.log(pred_base)  # log residual

    # Thin만
    sn = sel[thin]; rn = resid[thin]
    print(f"\n  Thin residual correlation (log(σ_actual/σ_pred) vs variables):")
    for name, arr in [('φ', pa), ('CN', cn), ('τ', tau), ('cov', cov), ('T/d', xi),
                       ('por', por), ('el_perc', el_perc), ('φ_SE', ps), ('σ_ion', sion),
                       ('δ', d), ('A', A), ('hop', hop), ('P', P), ('CN_std', cn_std),
                       ('d_AM', d_am), ('δ/hop', d/hop), ('A/d_AM²', A/d_am**2)]:
        corr = np.corrcoef(np.log(arr[thin]), rn)[0, 1]
        flag = '★' if abs(corr) > 0.4 else '●' if abs(corr) > 0.25 else ' '
        print(f"  {flag}{name:12s}: r={corr:+.3f}")

    # ═══════════════════════════════════════════
    # 2. Residual에서 발견된 변수로 보정
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("2. 보정항 추가: base × X^a (전체 n에서)")
    print("="*80)
    for xn, xa in [('CN', cn), ('cov', cov), ('por', por), ('el_perc', el_perc),
                    ('φ_SE', ps), ('σ_ion', sion), ('hop', hop), ('P', P),
                    ('CN_std', cn_std), ('δ/hop', d/hop), ('A/d²', A/d_am**2),
                    ('√τ', tau**0.5), ('τ', tau)]:
        for a in [-1, -0.75, -0.5, -0.25, 0.25, 0.5, 0.75, 1]:
            rhs = rhs_base * xa**a
            C = fitC(sel, rhs)
            if C is None: continue
            r2 = r2l(sel, C*rhs)
            if r2 > 0.905:
                cv = loocv_C(sel, rhs)
                # thick/thin
                pred = C * rhs
                r2_tk = r2l(sel[thick], pred[thick]) if thick.sum()>=3 else 0
                r2_tn = r2l(sel[thin], pred[thin]) if thin.sum()>=3 else 0
                print(f"  ★ ×{xn}^{a}: R²={r2:.4f} LOOCV={cv:.4f} | tk={r2_tk:.3f} tn={r2_tn:.3f}")

    # ═══════════════════════════════════════════
    # 3. Thin-only free regression on residual
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("3. Thin-only: residual 보정 free regression")
    print("="*80)
    var_thin = {
        'CN': np.log(cn[thin]), 'τ': np.log(tau[thin]), 'cov': np.log(cov[thin]),
        'por': np.log(por[thin]), 'el_perc': np.log(el_perc[thin]),
        'φ_SE': np.log(ps[thin]), 'σ_ion': np.log(sion[thin]),
        'hop': np.log(hop[thin]), 'P': np.log(P[thin]),
        'CN_std': np.log(cn_std[thin]), 'δ/hop': np.log(d[thin]/hop[thin]),
    }
    results = []
    for nv in [1, 2]:
        for vn in combinations(var_thin.keys(), nv):
            X = np.column_stack([var_thin[v] for v in vn] + [np.ones(thin.sum())])
            try:
                coefs, _, _, _ = np.linalg.lstsq(X, rn, rcond=None)
                pred_resid = X @ coefs
                # Apply correction to full dataset
                # For thin: use regression, for thick: no correction (resid=0)
                ss_resid = np.sum((rn - pred_resid)**2)
                ss_tot = np.sum((rn - np.mean(rn))**2)
                r2_resid = 1 - ss_resid / ss_tot if ss_tot > 0 else 0
                if r2_resid > 0.15:
                    terms = ' + '.join(f"{coefs[j]:.2f}×{v}" for j, v in enumerate(vn))
                    results.append({'r2': r2_resid, 'vnames': vn, 'coefs': coefs, 'terms': terms})
            except: pass
    results.sort(key=lambda x: -x['r2'])
    print(f"  Top 15 (thin residual R²):")
    for i, r in enumerate(results[:15]):
        print(f"    #{i+1} R²={r['r2']:.4f}  {r['terms']}")

    # ═══════════════════════════════════════════
    # 4. 전체 + 보정항 조합 (best base × thin correction)
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("4. 전체 fine-tuning: φ^a × (T/d)^b × δ^c / A^e × X^f × Y^g")
    print("="*80)
    c4 = []
    for a in [3, 3.25, 3.5, 3.75, 4]:
        for b in [-0.5, -0.25, 0]:
            for c in [2, 2.25, 2.5]:
                for e in [0.75, 1, 1.25]:
                    for xn, xa, xes in [('por', por, [-0.5, 0.5, 1]),
                                         ('cov', cov, [-0.5, 0.5, 1]),
                                         ('el_perc', el_perc, [0.5, 1]),
                                         ('φ_SE', ps, [0.5, 1, 1.5]),
                                         ('hop', hop, [-0.5, 0.5]),
                                         ('τ', tau, [0.25, 0.5])]:
                        for xe in xes:
                            rhs = SAM * pa**a * xi**b * d**c / A**e * xa**xe
                            C = fitC(sel, rhs)
                            if C is None: continue
                            r2 = r2l(sel, C*rhs)
                            if r2 > 0.90:
                                td_s = f'×(T/d)^{b}' if b else ''
                                c4.append({'r2':r2, 'rhs':rhs,
                                    'label':f'φ^{a}{td_s}×δ^{c}/A^{e}×{xn}^{xe}'})
    c4.sort(key=lambda x: -x['r2'])
    for i, r in enumerate(c4[:20]):
        cv = loocv_C(sel, r['rhs']) if i < 10 else 0
        cs = f" LOOCV={cv:.4f}" if cv else ""
        C = fitC(sel, r['rhs']); pred = C * r['rhs']
        r2_tk = r2l(sel[thick], pred[thick]) if thick.sum()>=3 else 0
        r2_tn = r2l(sel[thin], pred[thin]) if thin.sum()>=3 else 0
        print(f"  #{i+1} R²={r['r2']:.4f}{cs} tk={r2_tk:.3f} tn={r2_tn:.3f}  {r['label']}")

    # ═══════════════════════════════════════════
    # 5. 6변수 (base + 2 extras)
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("5. 6변수: φ^a × (T/d)^b × δ^c / A^e × X × Y")
    print("="*80)
    if c4:
        for base in c4[:5]:
            for yn, ya in [('×CN^0.5', cn**0.5), ('×CN^-0.5', cn**-0.5),
                           ('×√τ', tau**0.5), ('×por^-0.5', por**-0.5),
                           ('×el_perc^0.5', el_perc**0.5),
                           ('×φ_SE^0.5', ps**0.5)]:
                rhs = base['rhs'] * ya
                C = fitC(sel, rhs)
                if C is None: continue
                r2 = r2l(sel, C*rhs)
                if r2 > base['r2'] + 0.005:
                    pred = C * rhs
                    r2_tk = r2l(sel[thick], pred[thick]) if thick.sum()>=3 else 0
                    r2_tn = r2l(sel[thin], pred[thin]) if thin.sum()>=3 else 0
                    cv = loocv_C(sel, rhs)
                    print(f"  R²={r2:.4f} LOOCV={cv:.4f} tk={r2_tk:.3f} tn={r2_tn:.3f}  {base['label']}{yn}")

    # ═══════════════════════════════════════════
    # 6. Thin per-case table (best formula)
    # ═══════════════════════════════════════════
    best_all = c4[0] if c4 else None
    if best_all:
        print(f"\n{'='*80}")
        print(f"6. Thin per-case: {best_all['label']}")
        print("="*80)
        C = fitC(sel, best_all['rhs']); pred = C * best_all['rhs']
        err = np.abs(sel - pred) / sel * 100
        names = [r['name'] for r in rows]
        print(f"  {'Name':35s} {'σ':>6s} {'pred':>6s} {'err%':>5s} {'T/d':>5s} {'CN':>5s} {'por':>5s} {'cov':>5s} {'el_p':>5s}")
        for i in np.where(thin)[0]:
            j = np.argsort(-err[thin])
            # Just print thin cases sorted by error
        thin_idx = np.where(thin)[0]
        for i in thin_idx[np.argsort(-err[thin_idx])]:
            r = rows[i]
            print(f"  {r['name']:35s} {sel[i]:6.2f} {pred[i]:6.2f} {err[i]:5.0f}% {xi[i]:5.1f} {cn[i]:5.2f} {por[i]:5.1f} {cov[i]:5.3f} {el_perc[i]:5.3f}")

if __name__ == '__main__':
    main()
