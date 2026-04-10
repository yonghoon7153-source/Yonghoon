"""
Thin 극한: R²>0.95 목표, 모든 변수 총동원
2-regime OK → thin만 집중
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
            am_cr = m.get('am_am_mean_contact_radius', 0)
            cn_se = m.get('se_se_cn', 0)
            if pa <= 0 or cn <= 0 or T <= 0 or am_delta <= 0 or am_area <= 0: continue
            R = d_am / 2
            rows.append({
                'sel': sel, 'pa': pa, 'cn': cn, 'tau': max(tau, 0.1),
                'cov': cov, 'xi': T/d_am, 'por': max(por, 0.1),
                'el_perc': el_perc, 'ps': max(ps, 0.01), 'sion': max(sion, 0.001),
                'A': am_area, 'd': am_delta, 'hop': max(am_hop, 0.1),
                'P': max(am_pres, 0.01), 'cn_std': max(am_cn_std, 0.01),
                'a_c': max(am_cr, 0.01), 'cn_se': max(cn_se, 0.01),
                'dR': am_delta/R, 'aR': max(am_cr,0.01)/R,
                'd2A': am_delta**2/am_area if am_area > 0 else 0.001,
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
    # All arrays
    keys = ['sel','pa','cn','tau','cov','xi','por','el_perc','ps','sion',
            'A','d','hop','P','cn_std','a_c','cn_se','dR','aR','d2A']
    V = {k: np.array([r[k] for r in rows]) for k in keys}
    n = len(rows)
    thin = V['xi'] < 10
    n_thin = thin.sum()
    print(f"n={n}, thin={n_thin}")

    # Thin arrays
    T = {k: V[k][thin] for k in keys}
    sn = T['sel']

    # ═══════════════════════════════════════════
    # 1. 모든 변수 상관관계 (thin only)
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print(f"1. Thin correlation with σ_el (n={n_thin})")
    print("="*80)
    log_sn = np.log(sn)
    for name in ['pa','cn','tau','cov','xi','por','el_perc','ps','sion',
                  'A','d','hop','P','cn_std','a_c','cn_se','dR','aR','d2A']:
        corr = np.corrcoef(np.log(T[name]), log_sn)[0,1]
        flag = '★' if abs(corr) > 0.5 else '●' if abs(corr) > 0.3 else ' '
        print(f"  {flag}{name:10s}: r={corr:+.3f}")

    # ═══════════════════════════════════════════
    # 2. Free regression (thin only, 20변수)
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("2. Free regression thin only")
    print("="*80)
    var_pool = {k: np.log(T[k]) for k in
                ['pa','cn','tau','cov','xi','por','el_perc','ps','sion',
                 'A','d','hop','P','cn_std','a_c','dR','d2A']}
    results = []
    for nv in [3, 4, 5]:
        for vn in combinations(var_pool.keys(), nv):
            X = np.column_stack([var_pool[v] for v in vn] + [np.ones(n_thin)])
            try:
                coefs, _, _, _ = np.linalg.lstsq(X, log_sn, rcond=None)
                pred = np.exp(X @ coefs)
                r2 = r2l(sn, pred)
                if r2 > 0.92:
                    results.append({'r2': r2, 'coefs': coefs, 'vnames': vn, 'pred': pred})
            except: pass
    results.sort(key=lambda x: -x['r2'])
    print(f"  Top 30 (n={n_thin}):")
    for i, r in enumerate(results[:30]):
        cv = loocv_C(sn, r['pred']) if i < 15 else 0
        cs = f" LOOCV={cv:.4f}" if cv else ""
        terms = ' × '.join(f"{v}^{r['coefs'][j]:.2f}" for j, v in enumerate(r['vnames']))
        print(f"    #{i+1} R²={r['r2']:.4f}{cs}  ({len(r['vnames'])}v) {terms}")

    # ═══════════════════════════════════════════
    # 3. P + δ/A 조합 (thin only)
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("3. P + δ/A 조합 (thin only, 깔끔 지수)")
    print("="*80)
    c3 = []
    for p_e in [0, 0.5, 1, 1.5, 2, 3]:
        for d_e in [0, 0.5, 1, 1.5, 2]:
            for a_e in [0, 0.25, 0.5, 0.75, 1]:
                for xn, xa, xes in [('CN', T['cn'], [0, 0.5, 1]),
                                     ('(T/d)', T['xi'], [-0.5, -0.25, 0]),
                                     ('el_p', T['el_perc'], [0, 0.5, 1]),
                                     ('CN_std', T['cn_std'], [0, -0.5])]:
                    for xe in xes:
                        if p_e == 0 and d_e == 0: continue
                        rhs = SAM * T['P']**p_e * T['d']**d_e / np.clip(T['A']**a_e, 1e-10, None) * xa**xe
                        C = fitC(sn, rhs)
                        if C is None: continue
                        r2 = r2l(sn, C*rhs)
                        if r2 > 0.92:
                            ps = f'P^{p_e}×' if p_e else ''
                            ds = f'δ^{d_e}' if d_e else ''
                            As = f'/A^{a_e}' if a_e else ''
                            xs = f'×{xn}^{xe}' if xe else ''
                            c3.append({'r2':r2, 'rhs':rhs,
                                       'label':f'{ps}{ds}{As}{xs}'})
    c3.sort(key=lambda x: -x['r2'])
    for i, r in enumerate(c3[:20]):
        cv = loocv_C(sn, r['rhs']) if i < 10 else 0
        cs = f" LOOCV={cv:.4f}" if cv else ""
        print(f"  #{i+1} R²={r['r2']:.4f}{cs}  {r['label']}")

    # ═══════════════════════════════════════════
    # 4. CN_std 기반 (thin only)
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("4. CN_std 기반 (접촉 균일성)")
    print("="*80)
    c4 = []
    for a in [0, 1, 2, 3]:
        for b in [0, 0.5, 1]:
            for c in [-1, -0.75, -0.5, -0.25, 0]:
                for d_e in [0, 0.5, 1, 1.5]:
                    for a_e in [0, 0.5, 0.75]:
                        for td in [-0.5, -0.25, 0]:
                            if a == 0 and b == 0 and c == 0: continue
                            rhs = SAM * T['pa']**a * T['cn']**b * T['cn_std']**c * T['d']**d_e / np.clip(T['A']**a_e, 1e-10, None) * T['xi']**td
                            C = fitC(sn, rhs)
                            if C is None: continue
                            r2 = r2l(sn, C*rhs)
                            if r2 > 0.90:
                                ps = f'φ^{a}×' if a else ''
                                cs2 = f'CN^{b}×' if b else ''
                                cs3 = f'CN_std^{c}×' if c else ''
                                ds = f'δ^{d_e}' if d_e else ''
                                As = f'/A^{a_e}' if a_e else ''
                                ts = f'×(T/d)^{td}' if td else ''
                                c4.append({'r2':r2, 'rhs':rhs,
                                           'label':f'{ps}{cs2}{cs3}{ds}{As}{ts}'})
    c4.sort(key=lambda x: -x['r2'])
    for i, r in enumerate(c4[:20]):
        cv = loocv_C(sn, r['rhs']) if i < 10 else 0
        cs = f" LOOCV={cv:.4f}" if cv else ""
        print(f"  #{i+1} R²={r['r2']:.4f}{cs}  {r['label']}")

    # ═══════════════════════════════════════════
    # 5. FORM X style (thin only)
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("5. FORM X style (thin only)")
    print("="*80)
    formx = [
        ('P³×CN/(T/d)^½', T['P']**3 * T['cn'] / T['xi']**0.5),
        ('P³×CN/(T/d)^½/CN_std^½', T['P']**3 * T['cn'] / T['xi']**0.5 / T['cn_std']**0.5),
        ('P²×δ×CN/(T/d)^½', T['P']**2 * T['d'] * T['cn'] / T['xi']**0.5),
        ('δ²×CN/(A×(T/d)^¼)', T['d']**2 * T['cn'] / (T['A'] * T['xi']**0.25)),
        ('δ²×el_p/(A×(T/d)^¼)', T['d']**2 * T['el_perc'] / (T['A'] * T['xi']**0.25)),
        ('δ²×CN/(A^¾×(T/d)^¼)', T['d']**2 * T['cn'] / (T['A']**0.75 * T['xi']**0.25)),
        ('P×δ×CN/(T/d)^½', T['P'] * T['d'] * T['cn'] / T['xi']**0.5),
        ('P^(3/2)×CN/(T/d)^½', T['P']**1.5 * T['cn'] / T['xi']**0.5),
        ('√P×δ×CN/(T/d)^¼', T['P']**0.5 * T['d'] * T['cn'] / T['xi']**0.25),
        ('δ^(3/2)×CN/(A^½×(T/d)^¼)', T['d']**1.5 * T['cn'] / (T['A']**0.5 * T['xi']**0.25)),
        ('P×CN/((T/d)^½×CN_std^½)', T['P'] * T['cn'] / (T['xi']**0.5 * T['cn_std']**0.5)),
        ('P²×CN/(T/d)^½', T['P']**2 * T['cn'] / T['xi']**0.5),
        ('⁴√[P⁴×CN⁴/(T/d)²]', (T['P']**4 * T['cn']**4 / T['xi']**2)**0.25),
        ('⁴√[P⁸×CN⁴/(T/d)²×CN_std²]', (T['P']**8 * T['cn']**4 / (T['xi']**2 * T['cn_std']**2))**0.25),
        ('⁴√[δ⁸×CN⁴/(A⁴×(T/d))]', (T['d']**8 * T['cn']**4 / (T['A']**4 * T['xi']))**0.25),
        ('por^2×CN×(T/d)^(-½)×cov', T['por']**2 * T['cn'] * T['xi']**(-0.5) * T['cov']),
        ('por^(5/2)×CN×(T/d)^(-½)×cov^(3/2)', T['por']**2.5 * T['cn'] * T['xi']**(-0.5) * T['cov']**1.5),
    ]
    print(f"  {'Formula':50s} {'R²':>7s} {'LOOCV':>7s}")
    print("  " + "-" * 66)
    for label, rhs_raw in sorted(formx, key=lambda x: -r2l(sn, fitC(sn,SAM*x[1])*(SAM*x[1])) if fitC(sn,SAM*x[1]) else -999):
        rhs = SAM * rhs_raw
        C = fitC(sn, rhs)
        if C is None: continue
        r2 = r2l(sn, C*rhs)
        cv = loocv_C(sn, rhs) if r2 > 0.88 else 0
        cs = f"{cv:7.4f}" if cv else "      -"
        flag = '★' if r2 > 0.93 else '●' if r2 > 0.90 else ' '
        print(f"  {flag}{label:49s} {r2:7.4f} {cs}  C={C:.6f}")

    # ═══════════════════════════════════════════
    # 6. Champion comparison
    # ═══════════════════════════════════════════
    print(f"\n{'='*80}")
    print("6. CHAMPION TABLE")
    print("="*80)
    champs = [
        ('P³×CN/(T/d)^½ [이전 best]', SAM * T['P']**3 * T['cn'] / T['xi']**0.5),
        ('por^(5/2)×CN×cov^(3/2)/(T/d)^½ [no contact]', SAM * T['por']**2.5 * T['cn'] * T['cov']**1.5 / T['xi']**0.5),
        ('φ²×CN²×√cov×exp(π/ξ) [universal]', SAM * T['pa']**2 * T['cn']**2 * T['cov']**0.5 * np.exp(np.pi/T['xi'])),
    ]
    # Add top from free regression
    if results:
        best_free = results[0]
        rhs_f = SAM * np.exp(np.column_stack([var_pool[v] for v in best_free['vnames']] + [np.ones(n_thin)]) @ best_free['coefs'])
        champs.append((f"Free fit best ({len(best_free['vnames'])}v)", rhs_f))
    # Add top from c3
    if c3: champs.append((f"P+δ/A best: {c3[0]['label']}", c3[0]['rhs']))
    if c4: champs.append((f"CN_std best: {c4[0]['label']}", c4[0]['rhs']))

    print(f"  {'Formula':55s} {'R²':>7s} {'LOOCV':>7s} {'|err|':>6s}")
    print("  " + "-" * 75)
    for label, rhs in sorted(champs, key=lambda x: -r2l(sn, fitC(sn,x[1])*x[1]) if fitC(sn,x[1]) else -999):
        C = fitC(sn, rhs)
        if C is None: continue
        pred = C * rhs
        r2 = r2l(sn, pred)
        cv = loocv_C(sn, rhs)
        err = np.mean(np.abs(sn - pred) / sn * 100)
        w20 = np.sum(np.abs(sn - pred) / sn < 0.2)
        flag = '★' if r2 > 0.93 else '●' if r2 > 0.90 else ' '
        print(f"  {flag}{label:54s} {r2:7.4f} {cv:7.4f} {err:5.0f}% ({w20}/{n_thin})")

if __name__ == '__main__':
    main()
