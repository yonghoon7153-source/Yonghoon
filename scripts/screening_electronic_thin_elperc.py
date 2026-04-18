"""
Thin with el_perc: f_p × topology 최적화
el_perc ≥ 0.85에서 C fitting, 전체에서 R² 평가
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
            el_perc = m.get('electronic_percolating_fraction', 0)
            ps = m.get('phi_se', 0)
            A = m.get('am_am_mean_area', 0)
            d = m.get('am_am_mean_delta', 0)
            hop = max(m.get('am_am_mean_hop', 0), 0.1)
            if pa <= 0 or cn <= 0 or T <= 0: continue
            rows.append({
                'sel': sel, 'pa': pa, 'cn': cn, 'tau': tau, 'cov': cov,
                'xi': T/d_am, 'por': max(por, 0.1), 'ep': max(el_perc, 0.01),
                'ps': max(ps, 0.01), 'd2a': d**2/A if A > 0 and d > 0 else 0.001,
                'hop': hop, 'name': mp.parent.name
            })
    seen = set(); u = []
    for r in rows:
        k = f"{r['pa']:.4f}_{r['xi']:.1f}"
        if k not in seen: seen.add(k); u.append(r)
    return u

def r2l(a, p):
    la, lp = np.log(a), np.log(p); return 1-np.sum((la-lp)**2)/np.sum((la-np.mean(la))**2)
def fitC_filtered(sn, rhs, ep, thresh=0.85):
    """Fit C only on percolating cases"""
    ok = (ep >= thresh) & (rhs > 0) & np.isfinite(rhs)
    if ok.sum() < 3: return fitC_all(sn, rhs)
    return float(np.exp(np.mean(np.log(sn[ok] / rhs[ok]))))
def fitC_all(sn, rhs):
    v = (rhs > 0) & np.isfinite(rhs)
    return float(np.exp(np.mean(np.log(sn[v] / rhs[v])))) if v.sum() >= 3 else None
def loocv_filtered(sn, rhs, ep, thresh=0.85):
    ok = ep >= thresh
    if ok.sum() < 5: return 0
    sn_ok = sn[ok]; rhs_ok = rhs[ok]
    n = len(sn_ok); la = np.log(sn_ok); lr = np.log(rhs_ok); errs = []
    for i in range(n):
        m = np.ones(n, bool); m[i] = False
        C_loo = float(np.exp(np.mean(la[m] - lr[m])))
        errs.append((la[i] - np.log(C_loo * rhs_ok[i]))**2)
    return 1 - np.sum(errs) / np.sum((la - np.mean(la))**2)

def main():
    rows = load_data()
    if not rows: print("데이터 없음!"); return
    V = {k: np.array([r[k] for r in rows]) for k in rows[0] if k != 'name'}
    names = [r['name'] for r in rows]
    n = len(rows); thin = V['xi'] < 10; sn = V['sel'][thin]
    T = {k: V[k][thin] for k in V}
    n_thin = thin.sum()
    perc_ok = T['ep'] >= 0.85
    print(f"n={n}, thin={n_thin}, percolating(≥0.85)={perc_ok.sum()}")

    # ═══════════════════════════════════════════
    # 1. f_p × topology grid (C from perc only, R² on perc only)
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print(f"1. f_p × topology (C fit on ep≥0.85, R² on ep≥0.85)")
    print("="*80)
    sn_ok = sn[perc_ok]
    c1 = []
    for a in [0.75, 1, 1.25, 1.5]:
        for b in [0.5, 1, 1.5, 2, 2.5]:
            for c_e in [0, 0.5, 1, 1.5]:
                for d_e in [0.25, 0.5]:
                    rhs_all = SAM * T['ep'] * T['cn']**a * T['por']**b * T['cov']**c_e / T['xi']**d_e
                    rhs_ok = rhs_all[perc_ok]
                    C = fitC_all(sn_ok, rhs_ok)
                    if C is None: continue
                    r2 = r2l(sn_ok, C*rhs_ok)
                    if r2 > 0.90:
                        cs = f'×cov^{c_e}' if c_e else ''
                        c1.append({'r2':r2, 'rhs_all':rhs_all, 'rhs_ok':rhs_ok, 'C':C,
                            'label':f'f_p×CN^{a}×por^{b}{cs}/(T/d)^{d_e}'})
    c1.sort(key=lambda x: -x['r2'])
    for i, r in enumerate(c1[:20]):
        cv = loocv_filtered(sn, r['rhs_all'], T['ep']) if i < 10 else 0
        cs = f" LOOCV={cv:.4f}" if cv else ""
        # Also show R² on ALL thin (including low perc)
        pred_all = r['C'] * r['rhs_all']
        r2_all = r2l(sn, pred_all)
        print(f"  #{i+1} R²={r['r2']:.4f}{cs} (all={r2_all:.3f})  {r['label']}")

    # ═══════════════════════════════════════════
    # 2. 깔끔한 식 + per-case
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("2. 깔끔한 후보 (per-case signed error)")
    print("="*80)
    clean = [
        ('f_p×CN^(5/4)×por×√cov/√ξ', T['ep']*T['cn']**1.25*T['por']*T['cov']**0.5/T['xi']**0.5),
        ('f_p×CN×por×√cov/√ξ', T['ep']*T['cn']*T['por']*T['cov']**0.5/T['xi']**0.5),
        ('f_p×CN×por×cov/√ξ', T['ep']*T['cn']*T['por']*T['cov']/T['xi']**0.5),
        ('f_p×CN×por^1.5×√cov/√ξ', T['ep']*T['cn']*T['por']**1.5*T['cov']**0.5/T['xi']**0.5),
        ('f_p×CN×por^2×cov/√ξ', T['ep']*T['cn']*T['por']**2*T['cov']/T['xi']**0.5),
        ('f_p×CN^(5/4)×por^1.5×√cov/√ξ', T['ep']*T['cn']**1.25*T['por']**1.5*T['cov']**0.5/T['xi']**0.5),
        ('f_p×CN^(5/4)×por×cov/√ξ', T['ep']*T['cn']**1.25*T['por']*T['cov']/T['xi']**0.5),
        ('f_p×CN×por^2.5×cov^1.5/√ξ', T['ep']*T['cn']*T['por']**2.5*T['cov']**1.5/T['xi']**0.5),
        # √ forms
        ('f_p×√[CN²×por²×cov/ξ]', T['ep']*np.sqrt(T['cn']**2*T['por']**2*T['cov']/T['xi'])),
        ('f_p×⁴√[CN⁵×por⁴×cov²/ξ²]', T['ep']*(T['cn']**5*T['por']**4*T['cov']**2/T['xi']**2)**0.25),
    ]
    tn_names = [names[i] for i in np.where(thin)[0]]
    print(f"  {'Formula':45s} {'R²(p)':>7s} {'R²(a)':>7s} {'LOOCV':>7s} {'range':>15s}")
    print("  " + "-" * 80)
    for label, rhs_raw in sorted(clean, key=lambda x: -r2l(sn[perc_ok], fitC_all(sn[perc_ok],SAM*x[1][perc_ok])*(SAM*x[1][perc_ok])) if fitC_all(sn[perc_ok],SAM*x[1][perc_ok]) else -999):
        rhs = SAM * rhs_raw
        C = fitC_filtered(sn, rhs, T['ep'])
        if C is None: continue
        pred = C * rhs
        r2_p = r2l(sn[perc_ok], pred[perc_ok])
        r2_a = r2l(sn, pred)
        cv = loocv_filtered(sn, rhs, T['ep'])
        errs = (pred[perc_ok] - sn[perc_ok]) / sn[perc_ok] * 100
        flag = '★' if r2_p > 0.93 else '●' if r2_p > 0.90 else ' '
        print(f"  {flag}{label:44s} {r2_p:7.4f} {r2_a:7.4f} {cv:7.4f} {errs.min():+.0f}~{errs.max():+.0f}%")

    # ═══════════════════════════════════════════
    # 3. Best per-case (perc only)
    # ═══════════════════════════════════════════
    if c1:
        best = c1[0]
        pred = best['C'] * best['rhs_all']
        print(f"\n{'='*80}")
        print(f"3. Best per-case: {best['label']} R²={best['r2']:.4f}")
        print("="*80)
        errs = (pred - sn) / sn * 100
        for j in np.argsort(-np.abs(errs)):
            flag = '✗' if T['ep'][j] < 0.85 else ('!' if abs(errs[j]) > 20 else ' ')
            print(f"  {flag}{tn_names[j]:35s} σ={sn[j]:6.2f} pred={pred[j]:6.2f} err={errs[j]:+5.0f}% "
                  f"ep={T['ep'][j]:.3f} CN={T['cn'][j]:.2f} por={T['por'][j]:.1f} ξ={T['xi'][j]:.1f}")

if __name__ == '__main__':
    main()
