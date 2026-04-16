#!/usr/bin/env python3
"""
Electronic conductivity — systematic parameter sweep for thick & thin regimes.
Ionic FORM X급의 정확도를 목표로 모든 변수/지수 조합을 탐색.
"""
import os, json, sys
import numpy as np
from pathlib import Path
from itertools import product

def load_all_electronic():
    base = Path(__file__).resolve().parent.parent / 'webapp'
    rows = []
    for search_dir in [base / 'results', base / 'archive']:
        if not search_dir.is_dir():
            continue
        for mp in sorted(search_dir.rglob('full_metrics.json')):
            try:
                m = json.loads(mp.read_text())
            except:
                continue
            sel = m.get('electronic_sigma_full_mScm', 0)
            if not sel or sel < 0.001:
                continue
            pa = m.get('phi_am', 0)
            cn = m.get('am_am_cn', 0)
            tau = max(m.get('tortuosity_recommended', m.get('tortuosity_mean', 1)), 0.1)
            cov_vals = [v for v in [m.get('coverage_AM_P_mean', 0), m.get('coverage_AM_S_mean', 0), m.get('coverage_AM_mean', 0)] if v > 0]
            cov = (sum(cov_vals) / len(cov_vals) / 100) if cov_vals else 0.2
            delta = max(m.get('am_am_mean_delta', 0), 0.001)
            area = max(m.get('am_am_mean_area', 0), 0.001)
            hop = max(m.get('am_am_mean_hop', 0), 0.01)
            gc = max(m.get('am_path_conductance_mean', 0), 0.0001)
            bottleneck = max(m.get('am_path_bottleneck_mean', 0), 0.0001)
            ram = max(m.get('r_AM_P', 0), m.get('r_AM_S', 0), m.get('r_AM', 0))
            dam = ram * 2 if ram > 0.1 else 5.0
            T = m.get('thickness_um', 0)
            por = max(m.get('porosity', 10), 0.1)
            phi_se = max(m.get('phi_se', 0.2), 0.01)
            ep = max(m.get('electronic_percolating_fraction', 0), 0.01)
            ea = max(m.get('electronic_active_fraction', 0), 0.01)
            cn_std = max(float(m.get('am_am_cn_std', 0) or 0), 0.01)
            n_contacts = max(int(float(m.get('am_am_n_contacts', 0) or 0)), 1)
            area_cv = max(float(m.get('am_am_area_cv', 0) or 0), 0.01)
            contact_r = max(float(m.get('am_am_mean_contact_radius', 0) or 0), 0.001)
            force = max(float(m.get('am_am_mean_force', 0) or 0), 0.001)
            pressure = max(float(m.get('am_am_mean_pressure', 0) or 0), 0.001)
            se_cn = max(float(m.get('se_se_cn', 0) or 0), 0.01)

            if pa <= 0 or cn <= 0 or T <= 0 or dam <= 0:
                continue

            ratio = T / dam
            rows.append({
                'sigma': sel, 'phi_am': pa, 'cn': cn, 'tau': tau, 'cov': cov,
                'delta': delta, 'area': area, 'hop': hop, 'ratio': ratio,
                'por': por, 'phi_se': phi_se, 'ep': ep, 'ea': ea,
                'gc': gc, 'bottleneck': bottleneck, 'cn_std': cn_std,
                'n_contacts': n_contacts, 'area_cv': area_cv,
                'contact_r': contact_r, 'force': force, 'pressure': pressure,
                'se_cn': se_cn, 'T': T, 'dam': dam,
                'path': str(mp.parent),
            })

    # Deduplicate by (phi_am, ratio)
    seen = set()
    unique = []
    for r in rows:
        k = f"{r['phi_am']:.4f}_{r['ratio']:.1f}"
        if k not in seen:
            seen.add(k)
            unique.append(r)
    return unique


def fit_and_score(log_sigma, log_rhs):
    """Fit C (1 free param) and return R², C."""
    if len(log_sigma) < 3:
        return -999, 0
    lnC = np.mean(log_sigma - log_rhs)
    pred = lnC + log_rhs
    ss_res = np.sum((log_sigma - pred) ** 2)
    ss_tot = np.sum((log_sigma - np.mean(log_sigma)) ** 2)
    if ss_tot < 1e-12:
        return -999, 0
    r2 = 1 - ss_res / ss_tot
    return r2, np.exp(lnC)


def sweep_thick(data):
    """Sweep thick regime (ratio >= 10)."""
    thick = [r for r in data if r['ratio'] >= 10]
    print(f"\n{'='*60}")
    print(f"THICK REGIME SWEEP (T/d >= 10, n={len(thick)})")
    print(f"{'='*60}")

    if len(thick) < 5:
        print("  Not enough thick data!")
        return

    s = np.log(np.array([r['sigma'] for r in thick]))

    # Available variables
    vars_dict = {
        'phi_am': np.log(np.array([r['phi_am'] for r in thick])),
        'cn': np.log(np.array([r['cn'] for r in thick])),
        'tau': np.log(np.array([r['tau'] for r in thick])),
        'cov': np.log(np.array([r['cov'] for r in thick])),
        'area': np.log(np.array([r['area'] for r in thick])),
        'contact_r': np.log(np.array([r['contact_r'] for r in thick])),
        'delta': np.log(np.array([r['delta'] for r in thick])),
        'hop': np.log(np.array([r['hop'] for r in thick])),
        'ep': np.log(np.array([r['ep'] for r in thick])),
        'ea': np.log(np.array([r['ea'] for r in thick])),
        'gc': np.log(np.array([r['gc'] for r in thick])),
        'por': np.log(np.array([r['por'] for r in thick])),
        'phi_se': np.log(np.array([r['phi_se'] for r in thick])),
        'force': np.log(np.array([r['force'] for r in thick])),
        'pressure': np.log(np.array([r['pressure'] for r in thick])),
        'n_contacts': np.log(np.array([r['n_contacts'] for r in thick])),
    }

    # Phase 1: 2-variable sweep (φ_am^a × X^b)
    print("\n--- Phase 1: φ_am^a × X^b ---")
    exponents_a = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0]
    exponents_b = [-1.0, -0.5, -0.25, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0]

    results = []
    for var_name, var_log in vars_dict.items():
        if var_name == 'phi_am':
            continue
        for a, b in product(exponents_a, exponents_b):
            log_rhs = a * vars_dict['phi_am'] + b * var_log
            r2, C = fit_and_score(s, log_rhs)
            results.append((r2, f"φ^{a} × {var_name}^{b}", a, var_name, b, C))

    results.sort(reverse=True)
    print(f"  Top 20:")
    for r2, formula, *_ in results[:20]:
        print(f"    R²={r2:.4f}  {formula}")

    # Phase 2: 3-variable sweep (φ_am^a × X^b × Y^c) — top 5 × all Y
    print("\n--- Phase 2: φ_am^a × X^b × Y^c ---")
    top5_2var = results[:5]
    results3 = []
    for _, _, a, var1_name, b, _ in top5_2var:
        for var2_name, var2_log in vars_dict.items():
            if var2_name in ('phi_am', var1_name):
                continue
            for c in [-0.5, -0.25, 0.25, 0.5, 0.75, 1.0, 1.5]:
                log_rhs = a * vars_dict['phi_am'] + b * vars_dict[var1_name] + c * var2_log
                r2, C = fit_and_score(s, log_rhs)
                results3.append((r2, f"φ^{a} × {var1_name}^{b} × {var2_name}^{c}", C))

    results3.sort(reverse=True)
    print(f"  Top 20:")
    for r2, formula, C in results3[:20]:
        print(f"    R²={r2:.4f}  C={C:.4f}  {formula}")

    # Phase 3: 4-variable sweep (best 3-var + Z^d)
    print("\n--- Phase 3: Best 3-var + Z^d ---")
    # Parse best 3-var formula to extend
    results4 = []
    for _, formula3, _ in results3[:3]:
        # Extract variables from formula string
        parts = formula3.split(' × ')
        used_vars = set()
        for p in parts:
            vname = p.split('^')[0].replace('φ', 'phi_am')
            used_vars.add(vname)

        for var_name, var_log in vars_dict.items():
            if var_name in used_vars:
                continue
            for d in [-0.5, -0.25, 0.25, 0.5, 1.0]:
                # Recompute 3-var log_rhs + d * var_log
                # Use free fit for simplicity
                X = np.column_stack([vars_dict[v.split('^')[0].replace('φ', 'phi_am')]
                                      for v in parts] + [var_log])
                # OLS with free exponents for 4 vars
                from numpy.linalg import lstsq
                X_full = np.column_stack([X, np.ones(len(s))])
                try:
                    coefs, _, _, _ = lstsq(X_full, s, rcond=None)
                    pred = X_full @ coefs
                    ss_res = np.sum((s - pred) ** 2)
                    ss_tot = np.sum((s - np.mean(s)) ** 2)
                    r2 = 1 - ss_res / ss_tot
                    results4.append((r2, f"{formula3} + {var_name} [free fit {len(parts)+1} vars]"))
                except:
                    pass

    results4.sort(reverse=True)
    if results4:
        print(f"  Top 10:")
        for r2, formula in results4[:10]:
            print(f"    R²={r2:.4f}  {formula}")

    # Phase 4: Free exponent OLS on all variables
    print("\n--- Phase 4: Free OLS (all variables) ---")
    all_vars = list(vars_dict.keys())
    X = np.column_stack([vars_dict[v] for v in all_vars])
    X_full = np.column_stack([X, np.ones(len(s))])
    from numpy.linalg import lstsq
    coefs, _, _, _ = lstsq(X_full, s, rcond=None)
    pred = X_full @ coefs
    ss_res = np.sum((s - pred) ** 2)
    ss_tot = np.sum((s - np.mean(s)) ** 2)
    r2_free = 1 - ss_res / ss_tot
    print(f"  R²={r2_free:.4f} (free OLS, {len(all_vars)} variables)")
    print(f"  Exponents:")
    for v, c in zip(all_vars, coefs[:-1]):
        print(f"    {v:15s}: {c:+.3f}")
    print(f"    {'lnC':15s}: {coefs[-1]:+.3f} (C={np.exp(coefs[-1]):.4f})")


def sweep_thin(data):
    """Sweep thin regime (ratio < 10)."""
    thin = [r for r in data if r['ratio'] < 10]
    print(f"\n{'='*60}")
    print(f"THIN REGIME SWEEP (T/d < 10, n={len(thin)})")
    print(f"{'='*60}")

    if len(thin) < 5:
        print("  Not enough thin data!")
        return

    s = np.log(np.array([r['sigma'] for r in thin]))

    vars_dict = {
        'phi_am': np.log(np.array([r['phi_am'] for r in thin])),
        'cn': np.log(np.array([r['cn'] for r in thin])),
        'tau': np.log(np.array([r['tau'] for r in thin])),
        'cov': np.log(np.array([r['cov'] for r in thin])),
        'area': np.log(np.array([r['area'] for r in thin])),
        'contact_r': np.log(np.array([r['contact_r'] for r in thin])),
        'delta': np.log(np.array([r['delta'] for r in thin])),
        'hop': np.log(np.array([r['hop'] for r in thin])),
        'ep': np.log(np.array([r['ep'] for r in thin])),
        'ea': np.log(np.array([r['ea'] for r in thin])),
        'gc': np.log(np.array([r['gc'] for r in thin])),
        'por': np.log(np.array([r['por'] for r in thin])),
        'phi_se': np.log(np.array([r['phi_se'] for r in thin])),
        'force': np.log(np.array([r['force'] for r in thin])),
        'pressure': np.log(np.array([r['pressure'] for r in thin])),
        'ratio': np.log(np.array([r['ratio'] for r in thin])),
        'bottleneck': np.log(np.array([r['bottleneck'] for r in thin])),
        'n_contacts': np.log(np.array([r['n_contacts'] for r in thin])),
    }

    # Phase 1: 2-variable sweep
    print("\n--- Phase 1: X^a × Y^b ---")
    var_names = list(vars_dict.keys())
    exponents = [-1.0, -0.5, -0.25, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0]

    results = []
    for i, v1 in enumerate(var_names):
        for j, v2 in enumerate(var_names):
            if j <= i:
                continue
            for a, b in product(exponents, exponents):
                log_rhs = a * vars_dict[v1] + b * vars_dict[v2]
                r2, C = fit_and_score(s, log_rhs)
                results.append((r2, f"{v1}^{a} × {v2}^{b}", v1, a, v2, b, C))

    results.sort(reverse=True)
    print(f"  Top 30:")
    for r2, formula, *_ in results[:30]:
        print(f"    R²={r2:.4f}  {formula}")

    # Phase 2: 3-variable (top 5 + Z)
    print("\n--- Phase 2: Best 2-var + Z^c ---")
    results3 = []
    for _, _, v1, a, v2, b, _ in results[:5]:
        for v3 in var_names:
            if v3 in (v1, v2):
                continue
            for c in [-1.0, -0.5, -0.25, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0]:
                log_rhs = a * vars_dict[v1] + b * vars_dict[v2] + c * vars_dict[v3]
                r2, C = fit_and_score(s, log_rhs)
                results3.append((r2, f"{v1}^{a} × {v2}^{b} × {v3}^{c}", C))

    results3.sort(reverse=True)
    print(f"  Top 30:")
    for r2, formula, C in results3[:30]:
        print(f"    R²={r2:.4f}  C={C:.4f}  {formula}")

    # Phase 3: Free OLS
    print("\n--- Phase 3: Free OLS (all variables) ---")
    from numpy.linalg import lstsq
    all_vars = list(vars_dict.keys())
    X = np.column_stack([vars_dict[v] for v in all_vars])
    X_full = np.column_stack([X, np.ones(len(s))])
    coefs, _, _, _ = lstsq(X_full, s, rcond=None)
    pred = X_full @ coefs
    ss_res = np.sum((s - pred) ** 2)
    ss_tot = np.sum((s - np.mean(s)) ** 2)
    r2_free = 1 - ss_res / ss_tot
    print(f"  R²={r2_free:.4f} (free OLS, {len(all_vars)} variables)")
    print(f"  Exponents:")
    for v, c in zip(all_vars, coefs[:-1]):
        print(f"    {v:15s}: {c:+.3f}")
    print(f"    {'lnC':15s}: {coefs[-1]:+.3f} (C={np.exp(coefs[-1]):.4f})")

    # Spearman rank correlation for trend analysis
    print("\n--- Spearman Rank Correlation (σ vs each variable) ---")
    from scipy.stats import spearmanr
    sigma_arr = np.array([r['sigma'] for r in thin])
    for vname in sorted(vars_dict.keys()):
        vals = np.exp(vars_dict[vname])
        rho, pval = spearmanr(vals, sigma_arr)
        sig = '***' if pval < 0.001 else '**' if pval < 0.01 else '*' if pval < 0.05 else ''
        print(f"    {vname:15s}: ρ={rho:+.3f} p={pval:.4f} {sig}")


if __name__ == '__main__':
    data = load_all_electronic()
    print(f"Total electronic data: {len(data)} unique cases")
    print(f"  Thick (T/d >= 10): {sum(1 for r in data if r['ratio'] >= 10)}")
    print(f"  Thin  (T/d <  10): {sum(1 for r in data if r['ratio'] < 10)}")

    sweep_thick(data)
    sweep_thin(data)

    print("\n" + "=" * 60)
    print("SWEEP COMPLETE")
    print("=" * 60)
