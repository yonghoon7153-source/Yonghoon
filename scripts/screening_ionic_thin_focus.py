#!/usr/bin/env python3
"""
Ionic Scaling Law — 박막 + 85:15 집중 개선
==========================================
현재 FORM X v3: CN^(2-¼lnCN) × (φ-φc)^(¾+¼lnτ) × cov^¼
R²=0.932, w20=47/54

문제점:
  1. 85:15 (φ_SE≈0.15 < φ_c=0.185) → (φ-φc) clip → 과소예측
  2. 박막(높은 τ) case에서 전반적 mismatch

시도:
  A. φ_c 최적화 (0.185 → lower)
  B. Soft percolation (hard clip → smooth transition)
  C. 추가 변수 (thickness, φ_AM, gb_density...)
  D. 2-regime (thick/thin ionic)
  E. Full quadratic with regularization
"""
import json, os, sys
import numpy as np
from pathlib import Path

SG = 3.0  # σ_grain mS/cm

def load_ionic_data():
    base = Path(__file__).resolve().parent.parent / 'webapp'
    rows = []
    for d in [base / 'results', base / 'archive']:
        if not d.is_dir():
            continue
        for mp in sorted(d.rglob('full_metrics.json')):
            try:
                m = json.loads(mp.read_text())
            except:
                continue
            sn = m.get('sigma_full_mScm', 0)
            if not sn or sn < 0.001:
                continue
            ps = m.get('phi_se', 0)
            pa = max(m.get('phi_am', 0), 0.01)
            tau = m.get('tortuosity_recommended', m.get('tortuosity_mean', 0))
            cn = m.get('se_se_cn', 0)
            gd = max(m.get('gb_density_mean', 0), 0.001)
            gp = max(m.get('path_conductance_mean', 0), 1e-6)
            T = m.get('thickness_um', 0)
            cov_vals = [v for v in [m.get('coverage_AM_P_mean', 0),
                                     m.get('coverage_AM_S_mean', 0),
                                     m.get('coverage_AM_mean', 0)] if v > 0]
            cov = (sum(cov_vals) / len(cov_vals) / 100) if cov_vals else 0.20
            por = max(m.get('porosity', 10), 0.1)
            fp = max(m.get('percolation_pct', 0) / 100, 0.01)
            r_se = max(m.get('r_SE', 0), 0.1)
            r_am = max(m.get('r_AM_P', 0), m.get('r_AM_S', 0), m.get('r_AM', 0), 0.1)
            if tau <= 0 or ps <= 0 or cn <= 0 or T <= 0:
                continue
            ratio = T / (r_am * 2) if r_am > 0.1 else 1.0
            rows.append({
                'sn': sn, 'ps': ps, 'pa': pa, 'tau': tau, 'cn': cn,
                'gd': gd, 'gp': gp, 'T': T, 'cov': cov, 'por': por,
                'fp': fp, 'r_se': r_se, 'r_am': r_am, 'ratio': ratio,
                'name': mp.parent.name
            })
    seen = set(); u = []
    for r in rows:
        k = f"{r['ps']:.4f}_{r['T']:.1f}_{r['tau']:.3f}"
        if k not in seen:
            seen.add(k); u.append(r)
    return u


def r2_log(actual, predicted):
    la, lp = np.log(actual), np.log(predicted)
    ss_res = np.sum((la - lp) ** 2)
    ss_tot = np.sum((la - np.mean(la)) ** 2)
    return 1 - ss_res / ss_tot if ss_tot > 1e-12 else -999


def fit_C_log(actual, rhs):
    v = (rhs > 0) & np.isfinite(rhs)
    if v.sum() < 3:
        return None
    return float(np.exp(np.mean(np.log(actual[v] / rhs[v]))))


def score(sn, rhs, label="", verbose=False):
    C = fit_C_log(sn, rhs)
    if C is None:
        return None
    pred = C * rhs
    r2 = r2_log(sn, pred)
    errs = np.abs(pred - sn) / sn * 100
    w20 = int(np.sum(errs < 20))
    n = len(sn)
    if verbose:
        print(f"  {label}")
        print(f"    R²={r2:.4f}  C={C:.4f}  |err|={np.mean(errs):.1f}%  w20={w20}/{n}")
    return {'r2': r2, 'C': C, 'w20': w20, 'err': np.mean(errs),
            'pred': pred, 'label': label}


def main():
    rows = load_ionic_data()
    n = len(rows)
    print(f"Total ionic data: n={n}")

    sn = np.array([r['sn'] for r in rows])
    ps = np.array([r['ps'] for r in rows])
    pa = np.array([r['pa'] for r in rows])
    tau = np.array([r['tau'] for r in rows])
    cn = np.array([r['cn'] for r in rows])
    cov = np.array([r['cov'] for r in rows])
    gd = np.array([r['gd'] for r in rows])
    gp = np.array([r['gp'] for r in rows])
    T = np.array([r['T'] for r in rows])
    por = np.array([r['por'] for r in rows])
    fp = np.array([r['fp'] for r in rows])
    ratio = np.array([r['ratio'] for r in rows])
    names = [r['name'] for r in rows]

    # Identify 85:15 cases
    is_85 = ps < 0.20  # φ_SE < 0.20 → likely 85:15
    print(f"85:15 cases (φ_SE < 0.20): {is_85.sum()}")
    print(f"Other cases: {(~is_85).sum()}")

    # Log variables
    lCN = np.log(cn)
    lTau = np.log(tau)
    lCov = np.log(cov)
    lGd = np.log(gd)
    lGp = np.log(gp)
    lT = np.log(T)
    lPor = np.log(por)
    lRatio = np.log(ratio)
    lFp = np.log(fp)
    ls = np.log(sn)

    # ═══════════════════════════════════════════════════════════
    print(f"\n{'='*70}")
    print("PHASE 1: φ_c 최적화 (현재 0.185)")
    print(f"{'='*70}")
    # ═══════════════════════════════════════════════════════════

    best_phic = []
    for phic in np.arange(0.05, 0.20, 0.005):
        phi_ex = np.clip(ps - phic, 0.001, None)
        lPhi = np.log(phi_ex)

        # Current v3: CN^2 × (φ-φc)^¾ × cov^¼ + interactions
        rhs = SG * phi_ex**0.75 * cn**2 * cov**0.25 * \
              np.exp(-0.25 * lCN**2 + 0.25 * lPhi * lTau)
        res = score(sn, rhs)
        if res:
            best_phic.append((res['r2'], res['w20'], phic, res['err']))

        # Original FORM X: (φ-φc)^¾ × CN × √cov / √τ
        rhs2 = SG * phi_ex**0.75 * cn * cov**0.5 / tau**0.5
        res2 = score(sn, rhs2)
        if res2:
            best_phic.append((res2['r2'], res2['w20'], phic, res2['err']))

    best_phic.sort(key=lambda x: (-x[1], -x[0]))
    print("Top by w20 then R²:")
    seen = set()
    for r2, w20, phic, err in best_phic[:20]:
        k = f"{phic:.3f}"
        if k in seen: continue
        seen.add(k)
        flag = " ★" if phic < 0.15 else ""
        print(f"  φ_c={phic:.3f}  R²={r2:.4f}  w20={w20}/{n}  |err|={err:.1f}%{flag}")

    # ═══════════════════════════════════════════════════════════
    print(f"\n{'='*70}")
    print("PHASE 2: Soft Percolation Functions")
    print(f"{'='*70}")
    # ═══════════════════════════════════════════════════════════

    soft_results = []

    for phic in [0.10, 0.12, 0.14, 0.15, 0.16, 0.185]:
        # A. Softplus: log(1 + exp(k*(φ-φc))) / k
        for k in [5, 10, 20, 50]:
            sp = np.log(1 + np.exp(k * (ps - phic))) / k
            sp = np.clip(sp, 0.001, None)
            for a in [0.5, 0.75, 1.0]:
                rhs = SG * sp**a * cn * cov**0.5 / tau**0.5
                res = score(sn, rhs)
                if res and res['r2'] > 0.90:
                    soft_results.append((res['r2'], res['w20'], res['err'],
                                         f"softplus(k={k},φc={phic})^{a} × CN × √cov/√τ"))

        # B. Exponential: φ × exp(-φc/φ)
        exp_phi = ps * np.exp(-phic / ps)
        for a in [0.5, 0.75, 1.0, 1.5]:
            rhs = SG * exp_phi**a * cn * cov**0.5 / tau**0.5
            res = score(sn, rhs)
            if res and res['r2'] > 0.90:
                soft_results.append((res['r2'], res['w20'], res['err'],
                                     f"[φ·exp(-{phic}/φ)]^{a} × CN × √cov/√τ"))

        # C. Sigmoid: φ^a / (φ^a + φc^a)
        for a in [2, 4, 6, 10]:
            sig = ps**a / (ps**a + phic**a)
            for b in [0.5, 0.75, 1.0]:
                rhs = SG * sig**b * cn * cov**0.5 / tau**0.5
                res = score(sn, rhs)
                if res and res['r2'] > 0.90:
                    soft_results.append((res['r2'], res['w20'], res['err'],
                                         f"sigmoid(a={a},φc={phic})^{b} × CN × √cov/√τ"))

        # D. Power smooth: (φ/φc)^n × clip(φ-φc)
        phi_ex = np.clip(ps - phic, 0.001, None)
        for pw in [0.5, 1.0, 2.0]:
            smooth = (ps / phic)**pw * phi_ex
            for a in [0.5, 0.75, 1.0]:
                rhs = SG * smooth**a * cn * cov**0.5 / tau**0.5
                res = score(sn, rhs)
                if res and res['r2'] > 0.90:
                    soft_results.append((res['r2'], res['w20'], res['err'],
                                         f"[(φ/{phic})^{pw}×(φ-{phic})]^{a} × CN × √cov/√τ"))

    soft_results.sort(key=lambda x: (-x[1], -x[0]))
    print("Top 20:")
    for r2, w20, err, label in soft_results[:20]:
        print(f"  R²={r2:.4f}  w20={w20}/{n}  |err|={err:.1f}%  {label}")

    # ═══════════════════════════════════════════════════════════
    print(f"\n{'='*70}")
    print("PHASE 3: Soft Percolation + Interaction Terms")
    print(f"{'='*70}")
    # ═══════════════════════════════════════════════════════════

    phase3 = []
    for phic in [0.10, 0.12, 0.14, 0.16, 0.185]:
        # Best soft functions from each type
        sp10 = np.log(np.clip(1 + np.exp(10 * (ps - phic)), 1e-10, None)) / 10
        sp10 = np.clip(sp10, 0.001, None)
        exp_phi = ps * np.exp(-phic / ps)
        hard = np.clip(ps - phic, 0.001, None)

        for phi_name, phi_val in [('hard', hard), ('softplus', sp10), ('exp', exp_phi)]:
            lP = np.log(phi_val)
            for cn_e in [1.0, 1.5, 2.0]:
                for phi_e in [0.5, 0.75, 1.0]:
                    for cov_e in [0, 0.25, 0.5]:
                        for tau_e in [0, -0.25, -0.5, 0.25]:
                            base = phi_e * lP + cn_e * lCN + cov_e * lCov + tau_e * lTau
                            # No interaction
                            rhs0 = SG * np.exp(base)
                            res0 = score(sn, rhs0)
                            if res0 and res0['r2'] > 0.92:
                                phase3.append((res0['r2'], res0['w20'], res0['err'],
                                    f"{phi_name}(φc={phic})^{phi_e}×CN^{cn_e}×cov^{cov_e}×τ^{tau_e}"))

                            # + CN² saturation
                            for cn2 in [-0.25, -0.5]:
                                rhs1 = SG * np.exp(base + cn2 * lCN**2)
                                res1 = score(sn, rhs1)
                                if res1 and res1['r2'] > 0.92:
                                    phase3.append((res1['r2'], res1['w20'], res1['err'],
                                        f"{phi_name}(φc={phic})^{phi_e}×CN^{cn_e}×cov^{cov_e}×τ^{tau_e}+CN²={cn2}"))

                            # + CN×τ
                            for cnt in [-0.25, -0.5]:
                                rhs2 = SG * np.exp(base + cnt * lCN * lTau)
                                res2 = score(sn, rhs2)
                                if res2 and res2['r2'] > 0.92:
                                    phase3.append((res2['r2'], res2['w20'], res2['err'],
                                        f"{phi_name}(φc={phic})^{phi_e}×CN^{cn_e}×cov^{cov_e}×τ^{tau_e}+CN×τ={cnt}"))

                            # + φ×τ
                            for pt in [0.25, 0.5]:
                                rhs3 = SG * np.exp(base + pt * lP * lTau)
                                res3 = score(sn, rhs3)
                                if res3 and res3['r2'] > 0.92:
                                    phase3.append((res3['r2'], res3['w20'], res3['err'],
                                        f"{phi_name}(φc={phic})^{phi_e}×CN^{cn_e}×cov^{cov_e}×τ^{tau_e}+φ×τ={pt}"))

                            # + CN² + φ×τ combo
                            for cn2 in [-0.25]:
                                for pt in [0.25]:
                                    rhs4 = SG * np.exp(base + cn2 * lCN**2 + pt * lP * lTau)
                                    res4 = score(sn, rhs4)
                                    if res4 and res4['r2'] > 0.92:
                                        phase3.append((res4['r2'], res4['w20'], res4['err'],
                                            f"{phi_name}(φc={phic})^{phi_e}×CN^{cn_e}×cov^{cov_e}×τ^{tau_e}+CN²={cn2}+φτ={pt}"))

                            # + CN×τ + φ×τ combo
                            for cnt in [-0.25, -0.5]:
                                for pt in [0.25]:
                                    rhs5 = SG * np.exp(base + cnt * lCN * lTau + pt * lP * lTau)
                                    res5 = score(sn, rhs5)
                                    if res5 and res5['r2'] > 0.92:
                                        phase3.append((res5['r2'], res5['w20'], res5['err'],
                                            f"{phi_name}(φc={phic})^{phi_e}×CN^{cn_e}×cov^{cov_e}×τ^{tau_e}+CNτ={cnt}+φτ={pt}"))

    phase3.sort(key=lambda x: (-x[1], -x[0]))
    print(f"Total candidates: {len(phase3)}")
    print("Top 30 by w20 then R²:")
    for r2, w20, err, label in phase3[:30]:
        print(f"  R²={r2:.4f}  w20={w20}/{n}  |err|={err:.1f}%  {label}")

    # ═══════════════════════════════════════════════════════════
    print(f"\n{'='*70}")
    print("PHASE 4: 추가 변수 탐색 (gd, gp, por, fp, ratio, T)")
    print(f"{'='*70}")
    # ═══════════════════════════════════════════════════════════

    # Best base from current: (φ-0.185)^0.75 × CN × √cov / √τ
    phi_ex = np.clip(ps - 0.185, 0.001, None)
    lPhi = np.log(phi_ex)
    base_rhs = 0.75 * lPhi + 1.0 * lCN + 0.5 * lCov - 0.5 * lTau

    extra_vars = {
        'gd': lGd, 'gp': lGp, 'por': lPor, 'fp': lFp,
        'ratio': lRatio, 'T': lT, 'φ_AM': np.log(pa),
        'φ_SE': np.log(ps), 'CN/τ': np.log(cn / tau),
        'CN×cov': np.log(cn * cov),
    }

    phase4 = []
    for vname, vlog in extra_vars.items():
        if not np.all(np.isfinite(vlog)):
            continue
        for e in [-1.0, -0.5, -0.25, 0.25, 0.5, 1.0]:
            rhs = SG * np.exp(base_rhs + e * vlog)
            res = score(sn, rhs)
            if res and res['r2'] > 0.92:
                phase4.append((res['r2'], res['w20'], res['err'],
                               f"FORMX + {vname}^{e}"))

    # 2-variable additions
    for v1, l1 in extra_vars.items():
        if not np.all(np.isfinite(l1)):
            continue
        for v2, l2 in extra_vars.items():
            if v2 <= v1 or not np.all(np.isfinite(l2)):
                continue
            for e1 in [-0.5, -0.25, 0.25, 0.5]:
                for e2 in [-0.5, -0.25, 0.25, 0.5]:
                    rhs = SG * np.exp(base_rhs + e1 * l1 + e2 * l2)
                    res = score(sn, rhs)
                    if res and res['r2'] > 0.93:
                        phase4.append((res['r2'], res['w20'], res['err'],
                                       f"FORMX + {v1}^{e1} + {v2}^{e2}"))

    phase4.sort(key=lambda x: (-x[1], -x[0]))
    print("Top 20:")
    for r2, w20, err, label in phase4[:20]:
        print(f"  R²={r2:.4f}  w20={w20}/{n}  |err|={err:.1f}%  {label}")

    # ═══════════════════════════════════════════════════════════
    print(f"\n{'='*70}")
    print("PHASE 5: Full Quadratic (OLS, all interactions)")
    print(f"{'='*70}")
    # ═══════════════════════════════════════════════════════════

    # 4 variables: CN, (φ-φc), cov, τ
    phi_ex = np.clip(ps - 0.185, 0.001, None)
    lP = np.log(phi_ex)
    y = ls - np.log(SG)

    # Design matrix
    X = np.column_stack([
        np.ones(n), lCN, lP, lCov, lTau,
        lCN**2, lP**2, lCov**2, lTau**2,
        lCN * lP, lCN * lCov, lCN * lTau, lP * lCov, lP * lTau, lCov * lTau
    ])
    labels_q = ['const', 'CN', 'φ', 'cov', 'τ',
                'CN²', 'φ²', 'cov²', 'τ²',
                'CN×φ', 'CN×cov', 'CN×τ', 'φ×cov', 'φ×τ', 'cov×τ']

    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    pred_q = X @ beta
    ss_res = np.sum((y - pred_q) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r2_q = 1 - ss_res / ss_tot

    pred_sn = SG * np.exp(pred_q)
    errs_q = np.abs(pred_sn - sn) / sn * 100
    w20_q = int(np.sum(errs_q < 20))

    print(f"Full quadratic OLS: R²={r2_q:.4f}  w20={w20_q}/{n}  |err|={np.mean(errs_q):.1f}%")
    print("Coefficients:")
    for l, b in zip(labels_q, beta):
        print(f"  {l:>10s} = {b:+.4f}{'  ★' if abs(b) > 0.3 else ''}")

    # LOOCV
    loocv_err = []
    for i in range(n):
        mask = np.ones(n, bool); mask[i] = False
        b_loo = np.linalg.lstsq(X[mask], y[mask], rcond=None)[0]
        p_loo = X[i] @ b_loo
        loocv_err.append((y[i] - p_loo) ** 2)
    loocv_r2 = 1 - np.sum(loocv_err) / ss_tot
    print(f"LOOCV R² = {loocv_r2:.4f}")

    # Per-case details
    print(f"\nPer-case (sorted by error):")
    idx_sorted = np.argsort(errs_q)[::-1]
    for i in idx_sorted[:15]:
        flag = "★" if errs_q[i] > 20 else " "
        print(f"  {flag} {names[i]:40s} pred={pred_sn[i]:.4f} act={sn[i]:.4f} err={errs_q[i]:.0f}% φ_SE={ps[i]:.3f} τ={tau[i]:.2f} CN={cn[i]:.1f}")

    # ═══════════════════════════════════════════════════════════
    print(f"\n{'='*70}")
    print("PHASE 5b: Full Quadratic with DIFFERENT φ_c values")
    print(f"{'='*70}")
    # ═══════════════════════════════════════════════════════════

    for phic in [0.05, 0.08, 0.10, 0.12, 0.14, 0.16, 0.185]:
        phi_ex_c = np.clip(ps - phic, 0.001, None)
        lP_c = np.log(phi_ex_c)
        X_c = np.column_stack([
            np.ones(n), lCN, lP_c, lCov, lTau,
            lCN**2, lP_c**2, lCov**2, lTau**2,
            lCN * lP_c, lCN * lCov, lCN * lTau, lP_c * lCov, lP_c * lTau, lCov * lTau
        ])
        beta_c = np.linalg.lstsq(X_c, y, rcond=None)[0]
        pred_c = X_c @ beta_c
        r2_c = 1 - np.sum((y - pred_c)**2) / ss_tot
        pred_sn_c = SG * np.exp(pred_c)
        errs_c = np.abs(pred_sn_c - sn) / sn * 100
        w20_c = int(np.sum(errs_c < 20))

        # LOOCV
        loo = []
        for i in range(n):
            m = np.ones(n, bool); m[i] = False
            b = np.linalg.lstsq(X_c[m], y[m], rcond=None)[0]
            loo.append((y[i] - X_c[i] @ b)**2)
        loo_r2 = 1 - np.sum(loo) / ss_tot
        print(f"  φ_c={phic:.3f}: R²={r2_c:.4f}  LOOCV={loo_r2:.4f}  w20={w20_c}/{n}  |err|={np.mean(errs_c):.1f}%")

    # ═══════════════════════════════════════════════════════════
    print(f"\n{'='*70}")
    print("PHASE 6: Soft Percolation × Full Interaction Sweep")
    print(f"{'='*70}")
    # ═══════════════════════════════════════════════════════════

    # Best soft percolation candidates combined with full interaction sweeps
    phase6 = []

    for phic in [0.08, 0.10, 0.12, 0.14, 0.16]:
        for soft_type in ['hard', 'softplus10', 'softplus20', 'exp']:
            if soft_type == 'hard':
                phi_val = np.clip(ps - phic, 0.001, None)
            elif soft_type == 'softplus10':
                phi_val = np.clip(np.log(1 + np.exp(10 * (ps - phic))) / 10, 0.001, None)
            elif soft_type == 'softplus20':
                phi_val = np.clip(np.log(1 + np.exp(20 * (ps - phic))) / 20, 0.001, None)
            elif soft_type == 'exp':
                phi_val = ps * np.exp(-phic / ps)

            lP = np.log(phi_val)

            for cn_e in [1.0, 1.5, 2.0]:
                for phi_e in [0.5, 0.75, 1.0]:
                    for cov_e in [0, 0.25, 0.5]:
                        for tau_e in [-0.5, -0.25, 0, 0.25]:
                            base = np.log(SG) + phi_e * lP + cn_e * lCN + cov_e * lCov + tau_e * lTau

                            # With all 2-way interactions
                            for cn2 in [0, -0.25]:
                                for cnt in [0, -0.25, -0.5]:
                                    for pt in [0, 0.25]:
                                        rhs_log = base + cn2 * lCN**2 + cnt * lCN * lTau + pt * lP * lTau
                                        rhs = np.exp(rhs_log)
                                        C = fit_C_log(sn, rhs)
                                        if C is None or C <= 0:
                                            continue
                                        pred = C * rhs
                                        r2 = r2_log(sn, pred)
                                        errs = np.abs(pred - sn) / sn * 100
                                        w20 = int(np.sum(errs < 20))
                                        if w20 >= 47:
                                            terms = f"{soft_type}(φc={phic})^{phi_e}×CN^{cn_e}×cov^{cov_e}×τ^{tau_e}"
                                            if cn2: terms += f"+CN²={cn2}"
                                            if cnt: terms += f"+CNτ={cnt}"
                                            if pt: terms += f"+φτ={pt}"
                                            phase6.append((r2, w20, np.mean(errs), terms))

    phase6.sort(key=lambda x: (-x[1], -x[0]))
    print(f"Candidates with w20≥47: {len(phase6)}")
    print("Top 30:")
    for r2, w20, err, label in phase6[:30]:
        print(f"  R²={r2:.4f}  w20={w20}/{n}  |err|={err:.1f}%  {label}")

    # ═══════════════════════════════════════════════════════════
    print(f"\n{'='*70}")
    print("PHASE 7: 85:15 per-case analysis (worst cases)")
    print(f"{'='*70}")
    # ═══════════════════════════════════════════════════════════

    # Use current FORM X for per-case analysis
    phi_ex = np.clip(ps - 0.185, 0.001, None)
    lPhi = np.log(phi_ex)
    rhs_v3 = SG * np.exp(0.75 * lPhi + 2.0 * lCN + 0.25 * lCov
                          - 0.25 * lCN**2 + 0.25 * lPhi * lTau)
    C_v3 = fit_C_log(sn, rhs_v3)
    pred_v3 = C_v3 * rhs_v3
    errs_v3 = np.abs(pred_v3 - sn) / sn * 100

    print(f"\nCurrent v3 worst cases:")
    idx_sorted = np.argsort(errs_v3)[::-1]
    for i in idx_sorted[:20]:
        flag = "★" if errs_v3[i] > 20 else " "
        print(f"  {flag} err={errs_v3[i]:5.0f}%  pred={pred_v3[i]:.4f}  act={sn[i]:.4f}  "
              f"φ_SE={ps[i]:.3f}  τ={tau[i]:.2f}  CN={cn[i]:.1f}  cov={cov[i]:.2f}  "
              f"gp={gp[i]:.4f}  {names[i]}")

    # Correlation of residuals with additional variables
    log_residual = np.log(sn / pred_v3)
    print(f"\nResidual correlation (log(actual/predicted) vs variables):")
    for vname, vvals in [('gd', gd), ('gp', gp), ('por', por), ('fp', fp),
                          ('ratio', ratio), ('T', T), ('φ_AM', pa), ('φ_SE', ps),
                          ('cov', cov), ('CN', cn), ('τ', tau)]:
        corr = np.corrcoef(np.log(vvals), log_residual)[0, 1]
        flag = " ★★" if abs(corr) > 0.3 else (" ★" if abs(corr) > 0.2 else "")
        print(f"  {vname:>10s}: r={corr:+.3f}{flag}")


if __name__ == '__main__':
    main()
