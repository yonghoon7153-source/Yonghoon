"""
v2.0 Exhaustive Scaling Law Screening
=====================================
Goal: Find a formula that works for ALL cases including extreme thin (τ=3~5).

Approach:
1. Load all cases with network solver results
2. Test every 1/2/3-variable combination with σ_brug as base
3. Also test direct (non-σ_brug-based) formulas
4. Report top candidates with R², per-regime accuracy
"""
import json, os, sys, numpy as np, warnings
from scipy import stats
from itertools import combinations, product

warnings.filterwarnings('ignore')
WEBAPP = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'webapp')
SIGMA_GRAIN = 3.0  # mS/cm


def load_data():
    """Load all cases with network solver σ_ionic from results + archive."""
    from pathlib import Path
    rows = []
    for base in [Path(WEBAPP) / 'results', Path(WEBAPP) / 'archive']:
        if not base.is_dir():
            continue
        for met_path in base.rglob('full_metrics.json'):
            try:
                with open(met_path) as f:
                    m = json.load(f)
            except:
                continue

            sigma_net = m.get('sigma_full_mScm', 0)
            if not sigma_net or sigma_net < 0.001:
                continue

            phi_se = m.get('phi_se', 0)
            phi_am = m.get('phi_am', 0)
            tau = m.get('tortuosity_recommended', m.get('tortuosity_mean', 0))
            f_perc = m.get('percolation_pct', 0) / 100
            cn = m.get('se_se_cn', 0)
            gb_d = m.get('gb_density_mean', 0)
            g_path = m.get('path_conductance_mean', 0)
            hop_area = m.get('path_hop_area_mean', 0)
            bottleneck = m.get('path_hop_area_min_mean', 0)
            T = m.get('thickness_um', 0)
            porosity = m.get('porosity', 0)
            am_cn = m.get('am_am_cn', 0)
            bulk_frac = m.get('bulk_resistance_fraction', 0.25)

            if tau <= 0 or phi_se <= 0 or cn <= 0 or gb_d <= 0 or T <= 0:
                continue

            sigma_brug = SIGMA_GRAIN * phi_se * f_perc / tau**2

            rows.append({
                'name': met_path.parent.name,
                'sigma_net': sigma_net,
                'sigma_brug': sigma_brug,
                'phi_se': phi_se,
                'phi_am': phi_am,
                'tau': tau,
                'f_perc': f_perc,
                'cn': cn,
                'gb_d': gb_d,
                'g_path': g_path if g_path > 0 else 1e-6,
                'hop_area': hop_area if hop_area > 0 else 1e-6,
                'bottleneck': bottleneck if bottleneck > 0 else 1e-6,
                'T': T,
                'porosity': porosity,
                'am_cn': am_cn,
                'bulk_frac': bulk_frac,
            })

    # Dedup
    seen = set()
    unique = []
    for r in rows:
        key = f"{r['phi_se']:.4f}_{r['T']:.1f}_{r['tau']:.3f}"
        if key not in seen:
            seen.add(key)
            unique.append(r)
    return unique


def compute_r2(actual, predicted):
    ss_res = np.sum((actual - predicted)**2)
    ss_tot = np.sum((actual - np.mean(actual))**2)
    return 1 - ss_res / ss_tot if ss_tot > 0 else 0


def screening(rows):
    n = len(rows)
    sigma_net = np.array([r['sigma_net'] for r in rows])
    log_sigma = np.log(sigma_net)

    # Base variables (log space)
    variables = {
        'phi_se': np.log(np.array([r['phi_se'] for r in rows])),
        'tau': np.log(np.array([r['tau'] for r in rows])),
        'cn': np.log(np.array([r['cn'] for r in rows])),
        'gb_d': np.log(np.array([r['gb_d'] for r in rows])),
        'g_path': np.log(np.array([r['g_path'] for r in rows])),
        'hop_area': np.log(np.array([r['hop_area'] for r in rows])),
        'bottleneck': np.log(np.array([r['bottleneck'] for r in rows])),
        'T': np.log(np.array([r['T'] for r in rows])),
        'f_perc': np.log(np.array([max(r['f_perc'], 0.5) for r in rows])),
        'phi_am': np.log(np.array([max(r['phi_am'], 0.01) for r in rows])),
    }

    # Combined variables
    variables['g_path_gb2'] = variables['g_path'] + 2*variables['gb_d']  # G_path × GB_d²
    variables['hop_gb2'] = variables['hop_area'] + 2*variables['gb_d']  # A_hop × GB_d²
    variables['bn_gb2'] = variables['bottleneck'] + 2*variables['gb_d']  # BN × GB_d²
    variables['gb2_T'] = 2*variables['gb_d'] + variables['T']  # GB_d² × T
    variables['tau_gb'] = variables['tau'] + variables['gb_d']  # τ × GB_d
    variables['cn_tau'] = variables['cn'] - variables['tau']  # CN / τ
    variables['phi_tau2'] = variables['phi_se'] - 2*variables['tau']  # φ/τ² (σ_brug without f_perc)

    var_names = list(variables.keys())

    results = []

    # ─── Strategy 1: σ_brug × C × X^a × Y^b (2-var + σ_brug base) ───
    sigma_brug = np.array([r['sigma_brug'] for r in rows])
    log_brug = np.log(sigma_brug)
    log_ratio = log_sigma - log_brug  # log(σ_net / σ_brug)

    for v1, v2 in combinations(var_names, 2):
        X = np.column_stack([variables[v1], variables[v2], np.ones(n)])
        try:
            coefs, res, rank, sv = np.linalg.lstsq(X, log_ratio, rcond=None)
            pred = log_brug + X @ coefs
            r2 = compute_r2(log_sigma, pred)
            C = np.exp(coefs[2])
            results.append({
                'type': 'brug_based',
                'vars': f'{v1}^{coefs[0]:.2f} × {v2}^{coefs[1]:.2f}',
                'formula': f'σ_brug × {C:.4f} × {v1}^{coefs[0]:.2f} × {v2}^{coefs[1]:.2f}',
                'r2': r2,
                'n_params': 3,
                'coefs': coefs.tolist(),
                'v1': v1, 'v2': v2,
            })
        except:
            pass

    # ─── Strategy 2: σ_brug × C × X^a (1-var) ───
    for v1 in var_names:
        X = np.column_stack([variables[v1], np.ones(n)])
        try:
            coefs, res, rank, sv = np.linalg.lstsq(X, log_ratio, rcond=None)
            pred = log_brug + X @ coefs
            r2 = compute_r2(log_sigma, pred)
            C = np.exp(coefs[1])
            results.append({
                'type': 'brug_1var',
                'vars': f'{v1}^{coefs[0]:.2f}',
                'formula': f'σ_brug × {C:.4f} × {v1}^{coefs[0]:.2f}',
                'r2': r2,
                'n_params': 2,
                'coefs': coefs.tolist(),
                'v1': v1,
            })
        except:
            pass

    # ─── Strategy 3: Direct power law (no σ_brug base) ───
    for v1, v2, v3 in combinations(var_names, 3):
        X = np.column_stack([variables[v1], variables[v2], variables[v3], np.ones(n)])
        try:
            coefs, res, rank, sv = np.linalg.lstsq(X, log_sigma, rcond=None)
            pred = X @ coefs
            r2 = compute_r2(log_sigma, pred)
            results.append({
                'type': 'direct_3var',
                'vars': f'{v1}^{coefs[0]:.2f} × {v2}^{coefs[1]:.2f} × {v3}^{coefs[2]:.2f}',
                'formula': f'{np.exp(coefs[3]):.4f} × {v1}^{coefs[0]:.2f} × {v2}^{coefs[1]:.2f} × {v3}^{coefs[2]:.2f}',
                'r2': r2,
                'n_params': 4,
                'coefs': coefs.tolist(),
            })
        except:
            pass

    # ─── Strategy 4: σ_brug × C × X^a × Y^b × Z^c (3-var) ───
    top_2var = sorted([r for r in results if r['type'] == 'brug_based'], key=lambda x: -x['r2'])[:10]
    for v1, v2, v3 in combinations(var_names, 3):
        X = np.column_stack([variables[v1], variables[v2], variables[v3], np.ones(n)])
        try:
            coefs, res, rank, sv = np.linalg.lstsq(X, log_ratio, rcond=None)
            pred = log_brug + X @ coefs
            r2 = compute_r2(log_sigma, pred)
            C = np.exp(coefs[3])
            results.append({
                'type': 'brug_3var',
                'vars': f'{v1}^{coefs[0]:.2f} × {v2}^{coefs[1]:.2f} × {v3}^{coefs[2]:.2f}',
                'formula': f'σ_brug × {C:.4f} × {v1}^{coefs[0]:.2f} × {v2}^{coefs[1]:.2f} × {v3}^{coefs[2]:.2f}',
                'r2': r2,
                'n_params': 4,
                'coefs': coefs.tolist(),
            })
        except:
            pass

    results.sort(key=lambda x: -x['r2'])
    return results


def main():
    rows = load_data()
    print(f"Loaded {len(rows)} cases")
    print(f"  τ range: {min(r['tau'] for r in rows):.2f} ~ {max(r['tau'] for r in rows):.2f}")
    print(f"  σ range: {min(r['sigma_net'] for r in rows):.4f} ~ {max(r['sigma_net'] for r in rows):.4f} mS/cm")

    # Identify thin100 cases
    thin100 = [r for r in rows if r['tau'] > 3]
    thick = [r for r in rows if r['tau'] <= 3]
    print(f"  Thick (τ≤3): {len(thick)}, Thin extreme (τ>3): {len(thin100)}")
    print()

    results = screening(rows)

    # ─── Top results ───
    print("=" * 80)
    print(f"TOP 20 FORMULAS (n={len(rows)}, ALL data including thin100)")
    print("=" * 80)

    seen_formulas = set()
    count = 0
    for r in results:
        if r['r2'] < 0.5:
            break
        # Dedup similar formulas
        key = f"{r['type']}_{r['vars']}"
        if key in seen_formulas:
            continue
        seen_formulas.add(key)

        # Per-regime R²
        sigma_net_all = np.array([row['sigma_net'] for row in rows])
        sigma_brug_all = np.array([row['sigma_brug'] for row in rows])

        count += 1
        if count > 20:
            break
        print(f"\n#{count}  R²={r['r2']:.4f}  ({r['type']}, {r['n_params']}p)")
        print(f"  {r['formula']}")

    # ─── Champion comparison ───
    print("\n" + "=" * 80)
    print("CHAMPION COMPARISON")
    print("=" * 80)

    # Current champion: σ_brug × C × (G_path × GB_d²)^(1/4) × CN²
    sigma_net_all = np.array([r['sigma_net'] for r in rows])
    sigma_brug_all = np.array([r['sigma_brug'] for r in rows])
    g_path = np.array([r['g_path'] for r in rows])
    gb_d = np.array([r['gb_d'] for r in rows])
    cn = np.array([r['cn'] for r in rows])
    tau = np.array([r['tau'] for r in rows])

    # Fixed exponent champion
    rhs = sigma_brug_all * (g_path * gb_d**2)**0.25 * cn**2
    C = np.mean(sigma_net_all / rhs)
    pred = C * rhs
    r2_champ = compute_r2(sigma_net_all, pred)
    print(f"\nCurrent champion (fixed 1/4, 2): R²={r2_champ:.4f}, C={C:.4f}")

    # Per-regime
    for label, mask_fn in [('τ≤2', lambda r: r['tau'] <= 2),
                           ('2<τ≤3', lambda r: 2 < r['tau'] <= 3),
                           ('τ>3', lambda r: r['tau'] > 3),
                           ('thick (τ≤1.5)', lambda r: r['tau'] <= 1.5)]:
        idx = [i for i, r in enumerate(rows) if mask_fn(r)]
        if idx:
            r2_sub = compute_r2(sigma_net_all[idx], pred[idx])
            err = np.mean(np.abs(sigma_net_all[idx] - pred[idx]) / sigma_net_all[idx]) * 100
            print(f"  {label}: n={len(idx)}, R²={r2_sub:.3f}, |err|={err:.1f}%")

    print()


if __name__ == '__main__':
    main()
