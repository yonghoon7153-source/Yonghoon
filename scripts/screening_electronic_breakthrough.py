"""
Electronic BREAKTHROUGH: FORM X 급 quantum jump를 목표
======================================================
Ionic FORM X 교훈:
  - Bruggeman(φ/τ²) → percolation (φ-φc)^α: 관점 전환이 핵심
  - τ² → √τ: 과도한 penalty 완화
  - coverage 추가: 계면 품질

Electronic에도 같은 전략:
  1. AM percolation threshold: (φ_AM - φ_c_AM)^α
  2. el_active_fraction, el_perc_fraction: network solver가 직접 계산한 값
  3. σ_ionic coupling: SE network이 AM 배열을 지배
  4. bulk_resistance_fraction: R_bulk/(R_bulk+R_constr) → 구조 정보
  5. Thick에서도 percolation 접근
  6. T/d 대신 (T/d - threshold)^α 시도
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
            am_cn = max(m.get('am_am_cn', 0.01), 0.01)
            cn = m.get('se_se_cn', 0)
            T = m.get('thickness_um', 0)
            tau = m.get('tortuosity_recommended', m.get('tortuosity_mean', 0))
            cov = max(m.get('coverage_AM_P_mean',
                          m.get('coverage_AM_S_mean',
                               m.get('coverage_AM_mean', 20))), 0.1) / 100
            r_am = max(m.get('r_AM_P', 0), m.get('r_AM_S', 0))
            d_am = r_am * 2 if r_am > 0.1 else 5.0
            por = m.get('porosity', 0)

            # Network solver electronic metrics
            el_act = max(m.get('electronic_active_fraction', 0), 0.01)
            el_perc = max(m.get('electronic_percolating_fraction', 0), 0.01)
            el_bulk_frac = m.get('electronic_bulk_frac', 0.25)

            # Ionic sigma (SE network quality → proxy for AM structure)
            sion = m.get('sigma_full_mScm', 0)

            # GB density, path conductance
            gd = m.get('gb_density_mean', 0)
            gp = max(m.get('path_conductance_mean', 0), 1e-6)

            # AM-AM contact mechanics (constriction resistance 직접 관련)
            am_area = max(m.get('am_am_mean_area', 0), 0.01)  # µm²
            am_cr = max(m.get('am_am_mean_contact_radius', 0), 0.01)  # µm
            am_delta = max(m.get('am_am_mean_delta', 0), 0.001)  # µm penetration
            am_force = max(m.get('am_am_mean_force', 0), 0.001)  # µN
            am_pres = max(m.get('am_am_mean_pressure', 0), 0.01)  # MPa
            am_hop = max(m.get('am_am_mean_hop', 0), 0.1)  # µm
            am_area_cv = m.get('am_am_area_cv', 0.5)
            am_perc_pct = max(m.get('am_percolation_pct', 0), 0.1)  # 0-100%

            if pa <= 0 or am_cn <= 0 or T <= 0: continue

            rows.append({
                'sel': sel, 'pa': pa, 'ps': ps, 'am_cn': am_cn, 'cn': cn,
                'T': T, 'd_am': d_am, 'tau': max(tau, 0.1), 'cov': cov,
                'ratio': T / d_am, 'por': max(por, 0.1),
                'el_act': el_act, 'el_perc': el_perc,
                'el_bulk_frac': max(el_bulk_frac, 0.01),
                'sion': max(sion, 0.001), 'gd': gd, 'gp': gp,
                'am_area': am_area, 'am_cr': am_cr, 'am_delta': am_delta,
                'am_force': am_force, 'am_pres': am_pres, 'am_hop': am_hop,
                'am_area_cv': max(am_area_cv, 0.01),
                'am_perc': am_perc_pct / 100,  # → 0~1 fraction
                'name': mp.parent.name
            })

    seen = set(); u = []
    for r in rows:
        k = f"{r['pa']:.4f}_{r['ratio']:.1f}"
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

def fit_multivar(log_y, var_dict, var_names, min_r2=0.90):
    """Multi-variable log-linear regression: log(y) = Σ ai*log(xi) + c"""
    n = len(log_y)
    X = np.column_stack([var_dict[v] for v in var_names] + [np.ones(n)])
    try:
        coefs, _, _, _ = np.linalg.lstsq(X, log_y, rcond=None)
        pred = X @ coefs
        ss_res = np.sum((log_y - pred)**2)
        ss_tot = np.sum((log_y - np.mean(log_y))**2)
        r2 = 1 - ss_res / ss_tot
        if r2 > min_r2:
            return {'r2': r2, 'coefs': coefs, 'vnames': var_names, 'pred': np.exp(pred)}
    except:
        pass
    return None


def main():
    rows = load_data()
    if not rows:
        print("데이터 없음! 서버에서 실행하세요.")
        return

    # Arrays
    sel = np.array([r['sel'] for r in rows])
    pa = np.array([r['pa'] for r in rows])
    ps = np.array([r['ps'] for r in rows])
    am_cn = np.array([r['am_cn'] for r in rows])
    cn = np.array([r['cn'] for r in rows])
    tau = np.array([r['tau'] for r in rows])
    cov = np.array([r['cov'] for r in rows])
    ratio = np.array([r['ratio'] for r in rows])
    por = np.array([r['por'] for r in rows])
    el_act = np.array([r['el_act'] for r in rows])
    el_perc = np.array([r['el_perc'] for r in rows])
    el_bulk_frac = np.array([r['el_bulk_frac'] for r in rows])
    sion = np.array([r['sion'] for r in rows])
    # AM-AM contact mechanics
    am_area = np.array([r['am_area'] for r in rows])
    am_cr = np.array([r['am_cr'] for r in rows])
    am_delta = np.array([r['am_delta'] for r in rows])
    am_force = np.array([r['am_force'] for r in rows])
    am_pres = np.array([r['am_pres'] for r in rows])
    am_hop = np.array([r['am_hop'] for r in rows])
    am_area_cv = np.array([r['am_area_cv'] for r in rows])
    am_perc = np.array([r['am_perc'] for r in rows])
    names = [r['name'] for r in rows]

    thick = ratio >= 10; thin = ratio < 10
    n = len(rows)
    print(f"n={n}: thick={thick.sum()}, thin={thin.sum()}")
    print(f"el_act range: {el_act.min():.3f} ~ {el_act.max():.3f}")
    print(f"el_perc range: {el_perc.min():.3f} ~ {el_perc.max():.3f}")
    print(f"σ_ion range: {sion.min():.4f} ~ {sion.max():.4f}")
    print(f"bulk_frac range: {el_bulk_frac.min():.3f} ~ {el_bulk_frac.max():.3f}")
    print(f"am_area range: {am_area.min():.3f} ~ {am_area.max():.3f} µm²")
    print(f"am_cr range: {am_cr.min():.3f} ~ {am_cr.max():.3f} µm")
    print(f"am_perc range: {am_perc.min():.3f} ~ {am_perc.max():.3f}")

    # ═══════════════════════════════════════════════════════════════
    # PART 1: THICK — percolation approach (φ_AM - φ_c)
    # ═══════════════════════════════════════════════════════════════
    print(f"\n{'='*80}")
    print("THICK: (φ_AM - φ_c) percolation approach")
    print("="*80)

    st = sel[thick]; pt = pa[thick]; ct = am_cn[thick]
    tt = tau[thick]; cvt = cov[thick]; port = por[thick]
    el_act_t = el_act[thick]; el_perc_t = el_perc[thick]
    sion_t = sion[thick]; bf_t = el_bulk_frac[thick]

    # 1a. φ_c scan for thick
    print("\n  [1a] φ_c scan: (φ_AM-φ_c)^a × CN^b × cov^c")
    for phi_c in [0.0, 0.10, 0.15, 0.18, 0.20, 0.25, 0.30]:
        phi_ex = np.clip(pt - phi_c, 0.001, None)
        for a in [2, 3, 3.5, 4, 4.5]:
            for b in [1, 1.5, 2]:
                for c in [0, 0.5, 1]:
                    rhs = SAM * phi_ex**a * ct**b * cvt**c
                    C = fitC(st, rhs)
                    if C is None: continue
                    r2 = r2l(st, C * rhs)
                    if r2 > 0.95:
                        cv = loocv_C(st, rhs)
                        print(f"    φ_c={phi_c:.2f} (φ-φc)^{a}×CN^{b}×cov^{c}: R²={r2:.4f} LOOCV={cv:.4f}")

    # 1b. Thick with τ, por, el_perc
    print("\n  [1b] Extended thick: + τ, por, el_perc, σ_ion, bulk_frac")
    var_thick = {
        'φ_AM': np.log(pt), 'CN': np.log(ct), 'τ': np.log(tt),
        'cov': np.log(cvt), 'por': np.log(port),
        'el_perc': np.log(el_perc_t), 'σ_ion': np.log(sion_t),
    }
    results_thick = []
    for nv in [3, 4, 5]:
        for vn in combinations(var_thick.keys(), nv):
            res = fit_multivar(np.log(st), var_thick, vn, min_r2=0.96)
            if res:
                results_thick.append(res)
    results_thick.sort(key=lambda x: -x['r2'])
    print(f"  Top 15 (n={thick.sum()}):")
    for i, r in enumerate(results_thick[:15]):
        cv = loocv_C(st, r['pred']) if i < 5 else 0
        cv_str = f" LOOCV={cv:.4f}" if cv else ""
        terms = ' × '.join(f"{v}^{r['coefs'][j]:.2f}" for j, v in enumerate(r['vnames']))
        print(f"    #{i+1} R²={r['r2']:.4f}{cv_str}  {terms}")

    # 1c. Thick clean exponents
    print("\n  [1c] Thick 깔끔한 지수:")
    for a, b, c, d in [
        (4, 1.5, 1, 0), (4, 1.5, 0.5, 0), (4, 1.5, 0.75, 0),
        (4, 1, 1, 0), (3.5, 1.5, 1, 0), (4.5, 1.5, 1, 0),
        (4, 1.5, 1, 0.5), (4, 1.5, 1, -0.5),  # with τ
        (3, 2, 1, 0), (3, 1.5, 1, 0),
    ]:
        rhs = SAM * pt**a * ct**b * cvt**c * tt**d
        C = fitC(st, rhs)
        if C is None: continue
        r2 = r2l(st, C * rhs)
        cv = loocv_C(st, rhs)
        t_str = f"×τ^{d}" if d != 0 else ""
        print(f"    φ^{a}×CN^{b}×cov^{c}{t_str}: R²={r2:.4f} LOOCV={cv:.4f} C={C:.6f}")

    # ═══════════════════════════════════════════════════════════════
    # PART 2: THIN — BREAKTHROUGH APPROACHES
    # ═══════════════════════════════════════════════════════════════
    print(f"\n{'='*80}")
    print("THIN: BREAKTHROUGH 접근들")
    print("="*80)

    sn = sel[thin]; pn = pa[thin]; psn = ps[thin]
    cn_n = am_cn[thin]; rn = ratio[thin]; porn = por[thin]
    cvn = cov[thin]; tn = tau[thin]; cn_se_n = cn[thin]
    el_act_n = el_act[thin]; el_perc_n = el_perc[thin]
    sion_n = sion[thin]; bf_n = el_bulk_frac[thin]
    # AM contact mechanics (thin)
    am_area_n = am_area[thin]; am_cr_n = am_cr[thin]
    am_delta_n = am_delta[thin]; am_force_n = am_force[thin]
    am_pres_n = am_pres[thin]; am_hop_n = am_hop[thin]
    am_perc_n = am_perc[thin]
    n_thin = thin.sum()
    names_thin = [names[i] for i in np.where(thin)[0]]

    # ───────────────────────────────────────────
    # 2a. APPROACH 1: el_perc 직접 사용 (network solver truth)
    # σ_el = C × σ_AM × el_perc^a × CN^b × ...
    # ───────────────────────────────────────────
    print(f"\n  [2a] el_perc 직접 사용 (n={n_thin})")
    for a in [0.5, 0.75, 1, 1.5, 2]:
        for b in [0.5, 0.75, 1, 1.5]:
            for c in [1, 1.5, 2, 2.5, 3]:
                for d in [0, 0.5, 0.75, 1]:
                    rhs = SAM * el_perc_n**a * cn_n**b * porn**c * cvn**d
                    C = fitC(sn, rhs)
                    if C is None: continue
                    r2 = r2l(sn, C * rhs)
                    if r2 > 0.88:
                        cv = loocv_C(sn, rhs)
                        print(f"    el_perc^{a}×CN^{b}×por^{c}×cov^{d}: R²={r2:.4f} LOOCV={cv:.4f}")

    # ───────────────────────────────────────────
    # 2b. APPROACH 2: σ_ionic coupling
    # ionic network quality → AM 배열 quality
    # ───────────────────────────────────────────
    print(f"\n  [2b] σ_ionic coupling")
    for a in [0.25, 0.5, 0.75, 1]:
        for b in [0.5, 0.75, 1, 1.5]:
            for c in [1, 1.5, 2, 2.5, 3]:
                for d in [-0.5, -0.25, 0]:
                    rhs = SAM * sion_n**a * cn_n**b * porn**c * rn**d
                    C = fitC(sn, rhs)
                    if C is None: continue
                    r2 = r2l(sn, C * rhs)
                    if r2 > 0.88:
                        cv = loocv_C(sn, rhs)
                        print(f"    σ_ion^{a}×CN^{b}×por^{c}×(T/d)^{d}: R²={r2:.4f} LOOCV={cv:.4f}")

    # ───────────────────────────────────────────
    # 2c. APPROACH 3: (φ_AM - φ_c) percolation for thin
    # ───────────────────────────────────────────
    print(f"\n  [2c] AM percolation threshold (thin)")
    for phi_c in [0.0, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35]:
        phi_ex = np.clip(pn - phi_c, 0.001, None)
        for a in [1, 1.5, 2, 2.5, 3, 4]:
            for b in [0.5, 1, 1.5]:
                for c in [1, 2, 3]:
                    for d in [-0.5, 0, 0.5]:
                        rhs = SAM * phi_ex**a * cn_n**b * porn**c * rn**d
                        C = fitC(sn, rhs)
                        if C is None: continue
                        r2 = r2l(sn, C * rhs)
                        if r2 > 0.88:
                            cv = loocv_C(sn, rhs)
                            print(f"    φc={phi_c:.2f} (φ-φc)^{a}×CN^{b}×por^{c}×(T/d)^{d}: R²={r2:.4f} LOOCV={cv:.4f}")

    # ───────────────────────────────────────────
    # 2d. APPROACH 4: Full free regression (all variables)
    # ───────────────────────────────────────────
    print(f"\n  [2d] Free regression (모든 변수 조합)")
    var_thin = {
        'φ_AM': np.log(pn), 'φ_SE': np.log(psn), 'CN': np.log(cn_n),
        'T/d': np.log(rn), 'por': np.log(porn), 'cov': np.log(cvn),
        'τ': np.log(tn), 'el_perc': np.log(el_perc_n),
        'σ_ion': np.log(sion_n),
        # AM-AM contact mechanics
        'A_contact': np.log(am_area_n), 'r_contact': np.log(am_cr_n),
        'δ_pen': np.log(am_delta_n), 'F_n': np.log(am_force_n),
        'P_contact': np.log(am_pres_n), 'hop': np.log(am_hop_n),
        'am_perc': np.log(am_perc_n),
    }

    results_thin = []
    for nv in [3, 4, 5]:
        for vn in combinations(var_thin.keys(), nv):
            res = fit_multivar(np.log(sn), var_thin, vn, min_r2=0.88)
            if res:
                results_thin.append(res)

    results_thin.sort(key=lambda x: -x['r2'])
    print(f"  Top 30 (n={n_thin}):")
    for i, r in enumerate(results_thin[:30]):
        cv = loocv_C(sn, r['pred']) if i < 10 else 0
        cv_str = f" LOOCV={cv:.4f}" if cv else ""
        terms = ' × '.join(f"{v}^{r['coefs'][j]:.2f}" for j, v in enumerate(r['vnames']))
        n_vars = len(r['vnames'])
        # 지수가 깔끔한지 체크
        clean = all(abs(r['coefs'][j] - round(r['coefs'][j]*4)/4) < 0.1
                     for j in range(n_vars))
        flag = '★' if clean else ' '
        print(f"    #{i+1}{flag} R²={r['r2']:.4f}{cv_str}  ({n_vars}v) {terms}")

    # ───────────────────────────────────────────
    # 2e. APPROACH 5: FORM X 스타일 (4th root approach)
    # σ = C × σ_AM × ⁴√[A^a × B^b × C^c / D^d]
    # ───────────────────────────────────────────
    print(f"\n  [2e] FORM X 스타일 (4th root)")
    for A_exp, A_name, A_arr in [
        (3, '(φ_AM-0.2)³', np.clip(pn-0.2, 0.001, None)**3),
        (4, 'CN⁴', cn_n**4),
        (3, 'CN³', cn_n**3),
        (4, 'por⁴', porn**4),
        (3, 'por³', porn**3),
        (8, 'por⁸', porn**8),
        (12, 'por¹²', porn**12),
    ]:
        for B_exp, B_name, B_arr in [
            (4, 'CN⁴', cn_n**4),
            (2, 'CN²', cn_n**2),
            (4, '(φ_SE×cov)⁴', (psn*cvn)**4),
            (2, '(φ_SE×cov)²', (psn*cvn)**2),
            (4, 'por⁴', porn**4),
            (8, 'por⁸', porn**8),
        ]:
            if A_name == B_name: continue
            for C_exp, C_name, C_arr in [
                (2, 'cov²', cvn**2),
                (1, 'cov', cvn),
                (2, '(T/d)²', rn**2),
                (2, '(φ_SE)²', psn**2),
                (1, '1', np.ones(n_thin)),
            ]:
                for D_exp, D_name, D_arr in [
                    (2, '(T/d)²', rn**2),
                    (1, '(T/d)', rn),
                    (1, '1', np.ones(n_thin)),
                ]:
                    inside = A_arr * B_arr * C_arr / np.clip(D_arr, 1e-10, None)
                    inside = np.clip(inside, 1e-30, None)
                    rhs = SAM * inside**0.25
                    C = fitC(sn, rhs)
                    if C is None: continue
                    r2 = r2l(sn, C * rhs)
                    if r2 > 0.88:
                        cv = loocv_C(sn, rhs)
                        label = f"⁴√[{A_name}×{B_name}×{C_name}/{D_name}]"
                        print(f"    R²={r2:.4f} LOOCV={cv:.4f}  {label}")

    # ───────────────────────────────────────────
    # 2f. APPROACH 6: √ form (square root)
    # σ = C × σ_AM × √[A × B / C]
    # ───────────────────────────────────────────
    print(f"\n  [2f] √ form")
    for A_name, A_arr in [
        ('CN²×por⁶', cn_n**2 * porn**6),
        ('CN²×por⁴', cn_n**2 * porn**4),
        ('CN³×por⁶', cn_n**3 * porn**6),
        ('CN⁴×por⁶', cn_n**4 * porn**6),
        ('CN²×por⁸', cn_n**2 * porn**8),
        ('(φ-0.2)²×CN²×por⁶', np.clip(pn-0.2,0.001,None)**2 * cn_n**2 * porn**6),
        ('el_perc²×CN²×por⁴', el_perc_n**2 * cn_n**2 * porn**4),
        ('σ_ion×CN×por⁴', sion_n * cn_n * porn**4),
    ]:
        for B_name, B_arr in [
            ('φ_SE²×cov²', psn**2 * cvn**2),
            ('φ_SE×cov', psn * cvn),
            ('cov²', cvn**2),
            ('cov', cvn),
            ('1', np.ones(n_thin)),
        ]:
            for D_name, D_arr in [
                ('(T/d)', rn), ('(T/d)²', rn**2), ('1', np.ones(n_thin)),
            ]:
                inside = A_arr * B_arr / np.clip(D_arr, 1e-10, None)
                rhs = SAM * np.sqrt(np.clip(inside, 1e-30, None))
                C = fitC(sn, rhs)
                if C is None: continue
                r2 = r2l(sn, C * rhs)
                if r2 > 0.88:
                    cv = loocv_C(sn, rhs)
                    label = f"√[{A_name}×{B_name}/{D_name}]"
                    print(f"    R²={r2:.4f} LOOCV={cv:.4f}  {label}")

    # ───────────────────────────────────────────
    # 2g. APPROACH 7: (T/d - tc) threshold percolation
    # T/d가 power law가 아니라 percolation이라면?
    # ───────────────────────────────────────────
    print(f"\n  [2g] T/d percolation: (T/d - tc)^α")
    for tc in [0, 1, 2, 3, 4, 5]:
        td_ex = np.clip(rn - tc, 0.01, None)
        for a in [0.5, 0.75, 1, 1.5, 2]:
            for b in [1, 1.5, 2, 3]:
                for c in [0.5, 1, 1.25]:
                    for d in [0.5, 1, 1.25]:
                        rhs = SAM * td_ex**a * cn_n**b * porn**c * psn**d
                        C = fitC(sn, rhs)
                        if C is None: continue
                        r2 = r2l(sn, C * rhs)
                        if r2 > 0.88:
                            cv = loocv_C(sn, rhs)
                            print(f"    tc={tc} (T/d-{tc})^{a}×CN^{b}×por^{c}×φ_SE^{d}: R²={r2:.4f} LOOCV={cv:.4f}")

    # ───────────────────────────────────────────
    # 2h. APPROACH 8: Dimensional analysis / Buckingham Pi
    # Π₁ = σ_el/(σ_AM × el_perc)
    # Π₂ = CN × por × (T/d)
    # ───────────────────────────────────────────
    print(f"\n  [2h] Dimensionless groups")
    # σ_el / (σ_AM × el_perc) should be ~ f(CN, por, T/d, cov)
    sigma_ratio_n = sn / (SAM * el_perc_n)
    var_pi = {
        'CN': np.log(cn_n), 'por': np.log(porn), 'T/d': np.log(rn),
        'cov': np.log(cvn), 'φ_SE': np.log(psn), 'τ': np.log(tn),
    }
    results_pi = []
    for nv in [2, 3, 4]:
        for vn in combinations(var_pi.keys(), nv):
            res = fit_multivar(np.log(sigma_ratio_n), var_pi, vn, min_r2=0.85)
            if res:
                results_pi.append(res)
    results_pi.sort(key=lambda x: -x['r2'])
    print(f"  σ/(σ_AM×el_perc) = f(vars) (n={n_thin}):")
    for i, r in enumerate(results_pi[:10]):
        terms = ' × '.join(f"{v}^{r['coefs'][j]:.2f}" for j, v in enumerate(r['vnames']))
        print(f"    #{i+1} R²={r['r2']:.4f}  {terms}")

    # ───────────────────────────────────────────
    # 2i. APPROACH 9: σ_el = σ_AM × el_perc × g(structure)
    # network solver가 이미 percolation을 계산했으니
    # 남은 것은 structure factor
    # ───────────────────────────────────────────
    print(f"\n  [2i] σ_AM × el_perc × structure_factor")
    for a in [1, 1.5, 2, 3]:
        for b in [0, 0.5, 1]:
            for c in [0, 0.5, 1]:
                rhs = SAM * el_perc_n * cn_n**a * porn**b * cvn**c
                C = fitC(sn, rhs)
                if C is None: continue
                r2 = r2l(sn, C * rhs)
                if r2 > 0.85:
                    cv = loocv_C(sn, rhs)
                    print(f"    el_perc×CN^{a}×por^{b}×cov^{c}: R²={r2:.4f} LOOCV={cv:.4f}")

    # ═══════════════════════════════════════════════════════════════
    # PART 3: UNIVERSAL (thick+thin in one formula?)
    # ═══════════════════════════════════════════════════════════════
    print(f"\n{'='*80}")
    print("UNIVERSAL: thick+thin 하나의 식?")
    print("="*80)

    var_all = {
        'φ_AM': np.log(pa), 'CN': np.log(am_cn), 'por': np.log(por),
        'cov': np.log(cov), 'T/d': np.log(ratio), 'τ': np.log(tau),
        'φ_SE': np.log(ps), 'el_perc': np.log(el_perc),
        'σ_ion': np.log(sion),
        'A_contact': np.log(am_area), 'r_contact': np.log(am_cr),
        'am_perc': np.log(am_perc),
    }
    results_uni = []
    for nv in [4, 5, 6]:
        for vn in combinations(var_all.keys(), nv):
            res = fit_multivar(np.log(sel), var_all, vn, min_r2=0.90)
            if res:
                results_uni.append(res)
    results_uni.sort(key=lambda x: -x['r2'])
    print(f"  Top 20 (n={n}):")
    for i, r in enumerate(results_uni[:20]):
        cv = loocv_C(sel, r['pred']) if i < 5 else 0
        cv_str = f" LOOCV={cv:.4f}" if cv else ""
        terms = ' × '.join(f"{v}^{r['coefs'][j]:.2f}" for j, v in enumerate(r['vnames']))
        print(f"    #{i+1} R²={r['r2']:.4f}{cv_str}  ({len(r['vnames'])}v) {terms}")

    # ═══════════════════════════════════════════════════════════════
    # PART 4: Per-case diagnostics
    # ═══════════════════════════════════════════════════════════════
    print(f"\n{'='*80}")
    print("DATA INSPECTION")
    print("="*80)
    print(f"\n  {'Name':35s} {'σ_el':>8s} {'φ_AM':>6s} {'CN':>5s} {'T/d':>5s} {'por':>5s} {'cov':>5s} {'φ_SE':>5s} {'el_p':>5s} {'σ_ion':>7s} {'A_c':>6s} {'r_c':>5s} {'regime':>6s}")
    print("  " + "-" * 125)
    for i in np.argsort(-sel):
        r = rows[i]
        regime = 'THICK' if ratio[i] >= 10 else 'thin'
        print(f"  {r['name']:35s} {r['sel']:8.3f} {r['pa']:6.3f} {r['am_cn']:5.2f} "
              f"{r['ratio']:5.1f} {r['por']:5.3f} {r['cov']:5.3f} {r['ps']:5.3f} "
              f"{r['el_perc']:5.3f} {r['sion']:7.4f} {r['am_area']:6.3f} {r['am_cr']:5.3f} {regime:>6s}")


if __name__ == '__main__':
    main()
