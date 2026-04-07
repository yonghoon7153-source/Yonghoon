#!/usr/bin/env python3
"""
DEM Contact Analysis - Standard mode (AM + SE, 2 types)
Full pipeline: porosity, interface area, coverage, percolation, tortuosity, ionic active AM

Usage:
    python analyze_contacts.py results/atoms.csv results/contacts.csv -o ./results -t "1:AM,2:SE" -s 1000
"""
import argparse
import os
import sys
import json
import numpy as np
import pandas as pd
from collections import defaultdict

sys.path.insert(0, os.path.dirname(__file__))
from dem_analysis_core import run_full_analysis


def load_atoms_raw(csv_path):
    df = pd.read_csv(csv_path)
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    df['id'] = df['id'].astype(int)
    df['type'] = df['type'].astype(int)
    atoms = {}
    for _, row in df.iterrows():
        atom = {
            'type': int(row['type']),
            'x': row['x'], 'y': row['y'], 'z': row['z'],
            'radius': row['radius'],
        }
        # Add stress if available (sim units)
        if 'c_strs[1]' in row and not pd.isna(row['c_strs[1]']):
            vol = (4.0 / 3.0) * np.pi * row['radius']**3
            if vol > 0:
                atom['sigma_xx'] = row['c_strs[1]'] / vol
                atom['sigma_yy'] = row['c_strs[2]'] / vol
                atom['sigma_zz'] = row['c_strs[3]'] / vol
        atoms[int(row['id'])] = atom
    return atoms, df


def load_contacts_raw(csv_path):
    df = pd.read_csv(csv_path)
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors='coerce')
    df['id1'] = df['id1'].astype(int)
    df['id2'] = df['id2'].astype(int)
    contacts = []
    for _, row in df.iterrows():
        contacts.append({
            'id1': int(row['id1']), 'id2': int(row['id2']),
            'fn': np.sqrt(row['fn_x']**2 + row['fn_y']**2 + row['fn_z']**2),
            'fn_x': row.get('fn_x', 0), 'fn_y': row.get('fn_y', 0), 'fn_z': row.get('fn_z', 0),
            'ft': np.sqrt(row['ft_x']**2 + row['ft_y']**2 + row['ft_z']**2),
            'contact_area': row['contact_area'],
            'delta': row['delta'],
        })
    return contacts, df


def save_results(results, atoms_raw, contacts_raw, df_atom, df_contact,
                 type_map, scale, output_dir):
    area_conv = 1.0 / (scale ** 2) * 1e12

    # Read box dimensions from input_params.json (default 0.05)
    _box_x, _box_y = 0.05, 0.05
    _ip_path = os.path.join(output_dir, 'input_params.json')
    if os.path.exists(_ip_path):
        with open(_ip_path) as _f:
            _ip = json.load(_f)
        _box_x = _ip.get('box_x', 0.05)
        _box_y = _ip.get('box_y', 0.05)

    # Contact Summary
    rows = []
    total_n = 0
    total_area = 0
    am_types_count = sum(1 for v in type_map.values() if 'AM' in v)
    is_bimodal = am_types_count > 1

    for ct, v in results['interface'].items():
        # Standard 모드에서 AM전체-SE는 AM-SE와 동일하므로 생략
        if ct == 'AM전체-SE' and not is_bimodal:
            continue
        rows.append({
            '접촉유형': ct, '접촉수': v['n_contacts'],
            '접촉면적_mean(μm²)': round(v['mean_area'], 4),
            '접촉면적_total(μm²)': round(v['total_area'], 2),
        })
        if ct != 'AM전체-SE':
            total_n += v['n_contacts']
            total_area += v['total_area']
    rows.append({
        '접촉유형': 'All',
        '접촉수': total_n,
        '접촉면적_mean(μm²)': round(total_area / total_n, 4) if total_n > 0 else 0,
        '접촉면적_total(μm²)': round(total_area, 2),
    })
    pd.DataFrame(rows).to_csv(os.path.join(output_dir, 'contact_summary.csv'), index=False)

    # Atom Statistics (입자수, 반지름, 영률)
    # Load Young's modulus from input_params.json if available
    input_params = {}
    params_path = os.path.join(output_dir, 'input_params.json')
    if os.path.exists(params_path):
        with open(params_path) as f:
            input_params = json.load(f)

    e_sim_list = input_params.get('youngs_modulus_sim', [])
    # SE real E ≈ 24 GPa
    se_real_e = '24 GPa'

    rows = []
    type_keys = sorted(type_map.keys())
    for t in type_keys:
        name = type_map[t]
        sub = {aid: a for aid, a in atoms_raw.items() if a['type'] == t}
        if not sub: continue
        rs = [a['radius'] for a in sub.values()]
        r_real = round(np.mean(rs) * scale, 2)

        # Young's modulus: sim value × scale → real value
        e_idx = t - 1  # type 1 → index 0
        if e_idx < len(e_sim_list):
            e_sim = e_sim_list[e_idx]
            e_real = e_sim * scale  # × 1000
            if e_real >= 1e9:
                e_str = f"{e_real/1e9:.1f} GPa"
            else:
                e_str = f"{e_real/1e6:.1f} MPa"
            # SE: 유효영률 옆에 실제값 표기
            if name == 'SE':
                e_str = f"{e_str} ({se_real_e})"
        else:
            e_str = '-'

        rows.append({
            '입자유형': name,
            '입자수': len(sub),
            '반지름(μm)': r_real,
            '영률': e_str,
        })
    pd.DataFrame(rows).to_csv(os.path.join(output_dir, 'atom_statistics.csv'), index=False)

    # Coordination Summary
    coord = defaultdict(int)
    for c in contacts_raw:
        coord[c['id1']] += 1
        coord[c['id2']] += 1
    rows = []
    for t, name in type_map.items():
        ids = [aid for aid, a in atoms_raw.items() if a['type'] == t]
        vals = np.array([coord.get(aid, 0) for aid in ids])
        rows.append({
            '입자유형': name, '입자수': len(ids),
            '배위수_mean': round(float(np.mean(vals)), 1),
            '배위수_std': round(float(np.std(vals)), 1),
            '배위수_min': int(np.min(vals)),
            '배위수_max': int(np.max(vals)),
        })
    pd.DataFrame(rows).to_csv(os.path.join(output_dir, 'coordination_summary.csv'), index=False)

    # Network Summary (comprehensive)
    # 구조 → 계면 → 이온경로(연결성→경로효율→경로품질) → 활성도 → 응력
    perc = results['percolation']
    cn = results['se_se_cn']
    tau = results['tortuosity']
    ionic = results['ionic_active']
    eff_cond = results.get('effective_conductivity')
    rows = [
        # ── 구조 ──
        {'지표': '── 구조 ──', '값': ''},
        {'지표': 'Porosity(%)', '값': round(results['porosity'], 2)},
        {'지표': '전극두께(μm)', '값': round(results['thickness_um'], 2)},
        # ── 계면 ──
        {'지표': '── 계면 ──', '값': ''},
        {'지표': 'AM-SE Total(μm²)', '값': round(results['interface'].get('AM전체-SE', {}).get('total_area', 0), 2)},
        {'지표': 'SE-SE Total(μm²)', '값': round(results['interface'].get('SE-SE', {}).get('total_area', 0), 2)},
    ]
    for lbl, v in results['coverage'].items():
        rows.append({'지표': f'Coverage {lbl}(%)', '값': round(v['mean'], 1)})
    rows += [
        # ── 이온경로: 연결성 ──
        {'지표': '── 이온경로: 연결성 ──', '값': ''},
        {'지표': 'SE-SE CN mean', '값': round(cn['mean'], 2)},
        {'지표': 'SE Cluster 수', '값': f"{perc.get('n_large_components', '-')}(≥10) / {perc['n_components']}"},
        {'지표': 'SE Percolation(%)', '값': round(perc['percolation_pct'], 1)},
        {'지표': 'Top Reachable(%)', '값': round(perc['top_reachable_pct'], 1)},
        # ── 이온경로: 경로 효율 ──
        {'지표': '── 이온경로: 경로 효율 ──', '값': ''},
        {'지표': 'Tortuosity mean', '값': round(tau['mean'], 2) if tau['mean'] else 'N/A'},
        {'지표': 'Tortuosity median', '값': round(tau.get('median', 0), 2) if tau.get('median') else 'N/A'},
        {'지표': 'Tortuosity std', '값': round(tau['std'], 2) if tau['std'] else 'N/A'},
        {'지표': 'GB Density(hops/μm)', '값': '-'},  # placeholder, updated after cluster calc
        # ── 이온경로: 경로 품질 ──
        {'지표': '── 이온경로: 경로 품질 ──', '값': ''},
        {'지표': 'Path Hop Area mean(μm²)', '값': '-'},
        {'지표': 'Path Bottleneck(μm²)', '값': '-'},
        {'지표': 'Path Conductance(μm²)', '값': '-'},
        # ── 활성도 ──
        {'지표': '── 활성도 ──', '값': ''},
    ]
    am_risk = results.get('am_isolation_risk')
    if am_risk:
        rows.append({'지표': 'AM-SE CN mean', '값': round(am_risk['am_se_cn_mean'], 2)})
    rows.append({'지표': 'Ionic Active AM(%)', '값': round(ionic['active_pct'], 1)})
    if am_risk:
        rows.append({'지표': 'AM Vulnerable(%)', '값': round(am_risk['vulnerable_pct'], 1)})
    # ── 이온전도 ──
    if eff_cond:
        rows.append({'지표': '── 이온전도 ──', '값': ''})
        rows.append({'지표': 'SE Volume Fraction', '값': round(eff_cond['phi_se'], 3)})
        rows.append({'지표': 'σ_eff/σ_bulk', '값': round(eff_cond['sigma_ratio'], 4)})
    # ── Network Solver (자동 추가) ──
    # These values come from network_conductivity.py merge → full_metrics.json
    # Read from full_metrics if already computed
    met_path = os.path.join(output_dir, 'full_metrics.json')
    if os.path.exists(met_path):
        with open(met_path) as _mf:
            _met = json.load(_mf)
        if _met.get('sigma_full_mScm'):
            rows.append({'지표': '── Network Solver ──', '값': ''})
            rows.append({'지표': 'σ_ionic (mS/cm)', '값': round(_met['sigma_full_mScm'], 4)})
            if _met.get('R_brug_over_full'):
                rows.append({'지표': 'R_brug (과대추정 배수)', '값': f"{_met['R_brug_over_full']:.1f}×"})
            if _met.get('bulk_resistance_fraction'):
                rows.append({'지표': 'Constriction 비율(%)', '값': round((1-_met['bulk_resistance_fraction'])*100, 1)})
            if _met.get('electronic_sigma_full_mScm'):
                rows.append({'지표': 'σ_electronic (mS/cm)', '값': round(_met['electronic_sigma_full_mScm'], 2)})
            if _met.get('thermal_sigma_full_mScm'):
                rows.append({'지표': 'σ_thermal (mS/cm equiv)', '값': round(_met['thermal_sigma_full_mScm'], 3)})
    # ── 응력 ──
    stress = results.get('stress')
    if stress:
        rows.append({'지표': '── 응력 ──', '값': ''})
        rows.append({'지표': 'Stress CV(%)', '값': round(stress['vm_cv'], 1)})
        for tn in ['AM_P', 'AM_S', 'SE']:
            if tn in stress['type_stress']:
                rows.append({'지표': f'σ_{tn}/σ_mean', '값': round(stress['type_stress'][tn]['ratio'], 3)})
    pd.DataFrame(rows).to_csv(os.path.join(output_dir, 'network_summary.csv'), index=False)

    # Auto-detect P:S ratio from mass (count × volume × density)
    # NCM811 density = 4800 kg/m³ for both AM_P and AM_S
    am_density = 4800  # kg/m³ (NCM811)
    am_p_atoms = [a for a in atoms_raw.values() if type_map.get(a['type']) == 'AM_P']
    am_s_atoms = [a for a in atoms_raw.values() if type_map.get(a['type']) == 'AM_S']

    # Standard mode: type_map has 'AM' (not AM_P/AM_S)
    # Determine if it's P or S from radius (AM_P ≈ 6μm, AM_S ≈ 2μm in sim units)
    am_generic = [a for a in atoms_raw.values() if type_map.get(a['type']) == 'AM']
    if am_generic and not am_p_atoms and not am_s_atoms:
        avg_r = np.mean([a['radius'] for a in am_generic])
        # In sim units: AM_P radius ≈ 0.006, AM_S radius ≈ 0.002
        if avg_r > 0.004:  # large AM → P only
            am_p_atoms = am_generic
        else:  # small AM → S only
            am_s_atoms = am_generic

    am_p_mass = sum(am_density * 4/3 * np.pi * a['radius']**3 for a in am_p_atoms)
    am_s_mass = sum(am_density * 4/3 * np.pi * a['radius']**3 for a in am_s_atoms)

    if am_p_mass > 0 and am_s_mass > 0:
        total_am_mass = am_p_mass + am_s_mass
        p_frac = am_p_mass / total_am_mass
        # Round to nearest 10% step: 0.3 → 3:7, 0.5 → 5:5, 0.7 → 7:3
        p_pct = round(p_frac * 10)
        s_pct = 10 - p_pct
        ps_ratio = f"{p_pct}:{s_pct}"
    elif am_p_mass > 0 and am_s_mass == 0:
        ps_ratio = "10:0"
    elif am_s_mass > 0 and am_p_mass == 0:
        ps_ratio = "0:10"
    else:
        ps_ratio = ""

    # Auto-detect AM:SE mass ratio
    se_density = 2000  # kg/m³ (Li₆PS₅Cl)
    se_atoms = [a for a in atoms_raw.values() if type_map.get(a['type']) == 'SE']
    total_am_mass_all = am_p_mass + am_s_mass
    se_mass = sum(se_density * 4/3 * np.pi * a['radius']**3 for a in se_atoms)
    if total_am_mass_all > 0 and se_mass > 0:
        total_mass = total_am_mass_all + se_mass
        am_pct = round(total_am_mass_all / total_mass * 100 / 5) * 5  # nearest 5%
        se_pct = 100 - am_pct
        am_se_ratio = f"{am_pct}:{se_pct}"
    else:
        am_se_ratio = ""

    # Particle counts and radii per type
    particle_info = {}
    for t, name in type_map.items():
        sub = [a for a in atoms_raw.values() if a['type'] == t]
        if sub:
            particle_info[f'n_{name}'] = len(sub)
            particle_info[f'r_{name}'] = round(np.mean([a['radius'] for a in sub]) * scale, 2)

    # Full metrics JSON
    metrics = {
        'porosity': results['porosity'],
        'thickness_um': results['thickness_um'],
        'plate_z_source': results['plate_z_source'],
        'ps_ratio': ps_ratio,
        'am_se_ratio': am_se_ratio,
        'se_se_cn': cn['mean'],
        'percolation_pct': perc['percolation_pct'],
        'top_reachable_pct': perc['top_reachable_pct'],
        'n_components': perc['n_components'],
        'n_large_components': perc.get('n_large_components', 0),
        'tortuosity_mean': tau['mean'],
        'tortuosity_median': tau.get('median'),
        'tortuosity_std': tau['std'],
        'tortuosity_use_median': tau.get('use_median', False),
        'tortuosity_recommended': tau.get('recommended', tau['mean']),
        'ionic_active_pct': ionic['active_pct'],
    }

    # GB density and path conductance from cluster paths (exact per-hop contact area)
    perc_clusters_path = os.path.join(output_dir, 'se_clusters.json')
    # Will be computed later, add placeholder; actual values filled after cluster computation

    metrics.update(particle_info)
    for ct, v in results['interface'].items():
        safe = ct.replace('-', '_')
        metrics[f'area_{safe}_total'] = v['total_area']
        metrics[f'area_{safe}_n'] = v['n_contacts']
        metrics[f'area_{safe}_mean'] = v['mean_area']
    for lbl, v in results['coverage'].items():
        metrics[f'coverage_{lbl}_mean'] = v['mean']
        metrics[f'coverage_{lbl}_std'] = v['std']
    # Stress (relative)
    stress = results.get('stress')
    if stress:
        metrics['stress_cv'] = stress['vm_cv']
        for tn, sv in stress['type_stress'].items():
            metrics[f'stress_ratio_{tn}'] = sv['ratio']
        metrics['stress_z_layer_cv'] = stress['z_layer_cv']
    # New metrics
    force_dist = results.get('force_dist', {})
    for ct, v in force_dist.items():
        safe = ct.replace('-', '_')
        metrics[f'fn_{safe}_mean'] = v['mean']
        metrics[f'fn_{safe}_max'] = v['max']
    cp = results.get('contact_pressure', {})
    if cp.get('overall'):
        metrics['contact_pressure_mean'] = cp['overall']['mean']
        metrics['contact_pressure_max'] = cp['overall']['max']
    overlap = results.get('overlap_ratio')
    if overlap:
        metrics['overlap_mean'] = overlap['mean']
        metrics['overlap_max'] = overlap['max']
        metrics['overlap_pct_above_5'] = overlap['pct_above_5']
    am_risk = results.get('am_isolation_risk')
    if am_risk:
        metrics['am_vulnerable_pct'] = am_risk['vulnerable_pct']
        metrics['am_se_cn_mean'] = am_risk['am_se_cn_mean']
    eff_cond = results.get('effective_conductivity')
    if eff_cond:
        metrics['sigma_ratio'] = eff_cond['sigma_ratio']
        metrics['phi_se'] = eff_cond['phi_se']
        metrics['phi_am'] = 1.0 - eff_cond['phi_se'] - results['porosity'] / 100.0
    am_am_cn = results.get('am_am_cn')
    if am_am_cn and am_am_cn.get('n_am', 0) > 0:
        metrics['am_am_cn'] = am_am_cn['mean']
        metrics['am_am_cn_std'] = am_am_cn['std']
        metrics['am_am_n_contacts'] = am_am_cn['n_contacts']
        metrics['am_am_mean_area'] = am_am_cn['mean_area']
        metrics['am_am_total_area'] = am_am_cn['total_area']
    # Target pressure from input_params (sim → real MPa)
    if input_params.get('target_press_sim'):
        metrics['target_pressure_mpa'] = round(input_params['target_press_sim'] * scale, 1)

    # ── Data quality warnings ──
    warnings = []
    tau = results['tortuosity']
    perc_data = results['percolation']
    overlap = results.get('overlap_ratio')

    # 1. Tortuosity anomaly: mean > 5 or std/mean > 0.5
    if tau.get('mean') and tau['mean'] > 5:
        med_info = f" (median={tau.get('median', 0):.2f}, using median)" if tau.get('use_median') else ""
        warnings.append({'type': 'tortuosity_extreme', 'severity': 'critical',
                        'msg': f"Tortuosity mean={tau['mean']:.1f} (>5): RVE artifact{med_info}"})
    elif tau.get('use_median'):
        warnings.append({'type': 'tortuosity_unstable', 'severity': 'warning',
                        'msg': f"Tortuosity mean={tau['mean']:.2f} vs median={tau['median']:.2f}: high dispersion, using median"})
    elif tau.get('mean') and tau['mean'] > 3:
        warnings.append({'type': 'tortuosity_high', 'severity': 'warning',
                        'msg': f"Tortuosity mean={tau['mean']:.1f} (>3): unusually tortuous paths"})

    # 2. Low percolation
    if perc_data.get('percolation_pct', 100) < 85:
        warnings.append({'type': 'percolation_low', 'severity': 'critical',
                        'msg': f"Percolation={perc_data['percolation_pct']:.1f}% (<85%): SE network severely fragmented"})
    elif perc_data.get('percolation_pct', 100) < 95:
        warnings.append({'type': 'percolation_marginal', 'severity': 'warning',
                        'msg': f"Percolation={perc_data['percolation_pct']:.1f}% (<95%): SE network connectivity marginal"})

    # 3. High porosity (poor compaction)
    if results.get('porosity', 0) > 25:
        warnings.append({'type': 'porosity_high', 'severity': 'warning',
                        'msg': f"Porosity={results['porosity']:.1f}% (>25%): insufficient compaction"})

    # 5. Too few particles (RVE representativeness)
    total_particles = sum(particle_info.get(f'n_{name}', 0) for name in type_map.values())
    if total_particles < 500:
        warnings.append({'type': 'rve_small', 'severity': 'warning',
                        'msg': f"Total particles={total_particles} (<500): RVE may not be representative"})

    # 6. SE-SE CN near percolation threshold
    if cn.get('mean', 10) < 3.5:
        warnings.append({'type': 'cn_critical', 'severity': 'critical',
                        'msg': f"SE-SE CN={cn['mean']:.1f} (<3.5): near percolation threshold, network may collapse"})

    # 7. Electronic Active AM (from network solver, if available)
    met_path = os.path.join(output_dir, 'full_metrics.json')
    if os.path.exists(met_path):
        with open(met_path) as _mf:
            _met = json.load(_mf)
        el_active = _met.get('electronic_active_fraction')
        if el_active is not None:
            el_pct = el_active * 100
            if el_pct < 10:
                warnings.append({'type': 'electronic_dead', 'severity': 'critical',
                                'msg': f"Electronic Active AM={el_pct:.0f}% (<10%): 도전재 필수! AM-AM percolation 없음"})
            elif el_pct < 50:
                warnings.append({'type': 'electronic_low', 'severity': 'critical',
                                'msg': f"Electronic Active AM={el_pct:.0f}% (<50%): 대량 dead AM, 도전재 강력 권장"})
            elif el_pct < 80:
                warnings.append({'type': 'electronic_marginal', 'severity': 'warning',
                                'msg': f"Electronic Active AM={el_pct:.0f}% (<80%): 일부 dead AM, 도전재 권장"})

    if warnings:
        metrics['warnings'] = warnings
        metrics['warning_count'] = len(warnings)
        for w in warnings:
            print(f"  ⚠ [{w['severity']}] {w['msg']}")

    with open(os.path.join(output_dir, 'full_metrics.json'), 'w') as f:
        json.dump(metrics, f, indent=2, default=str)

    # Save percolation sets for 3D viewer
    perc = results['percolation']
    perc_sets = {
        'top_reachable': list(perc.get('top_reachable_se', set())),
        'bottom_se': list(perc.get('bottom_se', set())),
        'top_se': list(perc.get('top_se', set())),
    }
    with open(os.path.join(output_dir, 'percolation_sets.json'), 'w') as f:
        json.dump(perc_sets, f)

    # Save SE cluster data for 3D viewer (cluster membership + paths)
    if 'graph' in perc:
        import networkx as nx
        G = perc['graph']
        components = list(nx.connected_components(G))
        bottom_se = perc.get('bottom_se', set())
        top_se = perc.get('top_se', set())

        clusters = []
        for comp in components:
            if len(comp) < 2:
                continue
            comp_list = list(comp)
            has_bottom = len(comp & bottom_se) > 0
            has_top = len(comp & top_se) > 0
            percolating = has_bottom and has_top

            cluster_info = {
                'ids': comp_list,
                'size': len(comp),
                'percolating': percolating,
                'path': None,
            }

            # Find multiple shortest paths, sorted by closeness to τ mean
            if percolating:
                src_candidates = list(comp & bottom_se)
                tgt_candidates = list(comp & top_se)
                import random
                random.seed(42)
                random.shuffle(src_candidates)
                random.shuffle(tgt_candidates)

                # Build contact area lookup: (id1,id2) → area (sim units)
                contact_area_map = {}
                for c in contacts_raw:
                    i1, i2 = c['id1'], c['id2']
                    contact_area_map[(min(i1,i2), max(i1,i2))] = c['contact_area']

                area_conv = 1.0 / (scale ** 2) * 1e12  # sim m² → real μm²

                all_paths = []
                seen_pairs = set()
                for si in range(min(15, len(src_candidates))):
                    for ti in range(min(15, len(tgt_candidates))):
                        if len(all_paths) >= 30:
                            break
                        src = src_candidates[si]
                        tgt = tgt_candidates[ti]
                        pair = (src, tgt)
                        if pair in seen_pairs:
                            continue
                        seen_pairs.add(pair)
                        try:
                            path = nx.shortest_path(G, src, tgt, weight='distance')
                            path_len = 0
                            # Exact per-hop contact area and resistance
                            hop_areas = []  # real μm²
                            sum_inv_a = 0  # Σ(1/A_i) for resistance
                            for ki in range(len(path)-1):
                                a1, a2 = atoms_raw[path[ki]], atoms_raw[path[ki+1]]
                                dx = abs(a1['x'] - a2['x'])
                                dy = abs(a1['y'] - a2['y'])
                                dz = a1['z'] - a2['z']
                                dx = min(dx, _box_x - dx)
                                dy = min(dy, _box_y - dy)
                                path_len += np.sqrt(dx**2 + dy**2 + dz**2)
                                # Lookup actual contact area
                                pair_key = (min(path[ki], path[ki+1]), max(path[ki], path[ki+1]))
                                ca = contact_area_map.get(pair_key, 0)
                                ca_real = ca * area_conv  # → μm²
                                hop_areas.append(ca_real)
                                if ca_real > 0:
                                    sum_inv_a += 1.0 / ca_real
                            z_dist = abs(atoms_raw[tgt]['z'] - atoms_raw[src]['z'])
                            n_hop = len(path) - 1
                            if z_dist > 0:
                                gb_density = round(n_hop / (z_dist * scale), 3)  # hops/μm
                                path_conductance = round(1.0 / sum_inv_a, 6) if sum_inv_a > 0 else 0  # μm² (effective)
                                all_paths.append({
                                    'ids': path,
                                    'tortuosity': round(path_len / z_dist, 2),
                                    'path_length': round(path_len * scale, 1),
                                    'z_distance': round(z_dist * scale, 1),
                                    'n_hop': n_hop,
                                    'gb_density': gb_density,
                                    'path_conductance': path_conductance,
                                    'hop_area_mean': round(np.mean(hop_areas), 4) if hop_areas else 0,
                                    'hop_area_min': round(min(hop_areas), 4) if hop_areas else 0,
                                })
                        except nx.NetworkXNoPath:
                            pass
                    if len(all_paths) >= 30:
                        break
                if all_paths:
                    # Sort: best τ (1~10), mean τ (11~20), worst τ (21~30)
                    tau_mean = results['tortuosity'].get('mean', 0) or 0
                    all_paths.sort(key=lambda p: p['tortuosity'])  # ascending
                    best = all_paths[:10]  # lowest τ
                    worst = all_paths[-10:] if len(all_paths) > 20 else all_paths[10:]
                    # mean-close: sort middle by closeness to mean
                    middle = [p for p in all_paths[10:] if p not in worst] if len(all_paths) > 10 else []
                    middle.sort(key=lambda p: abs(p['tortuosity'] - tau_mean))
                    middle = middle[:10]
                    paths_list = best + middle + worst
                    # Tag each path
                    for i, p in enumerate(paths_list):
                        if i < 10:
                            p['category'] = 'best'
                        elif i < 20:
                            p['category'] = 'mean'
                        else:
                            p['category'] = 'worst'
                    cluster_info['path'] = paths_list[0]
                    cluster_info['paths'] = paths_list

            clusters.append(cluster_info)

        # Sort by size (largest first)
        clusters.sort(key=lambda c: c['size'], reverse=True)

        # Build SE id → cluster index map
        se_cluster_map = {}
        for i, cl in enumerate(clusters):
            for sid in cl['ids']:
                se_cluster_map[sid] = i

        with open(os.path.join(output_dir, 'se_clusters.json'), 'w') as f:
            json.dump({'clusters': clusters, 'se_cluster_map': se_cluster_map}, f)

        # Update metrics with GB density and path conductance from percolating cluster paths
        all_perc_paths = []
        for cl in clusters:
            if cl.get('paths'):
                all_perc_paths.extend(cl['paths'])
        if all_perc_paths:
            gb_densities = [p['gb_density'] for p in all_perc_paths if 'gb_density' in p]
            conductances = [p['path_conductance'] for p in all_perc_paths if p.get('path_conductance', 0) > 0]
            hop_areas = [p['hop_area_mean'] for p in all_perc_paths if 'hop_area_mean' in p]
            hop_area_mins = [p['hop_area_min'] for p in all_perc_paths if 'hop_area_min' in p]

            metrics_update = {}
            if gb_densities:
                metrics_update['gb_density_mean'] = round(float(np.mean(gb_densities)), 3)
            if conductances:
                metrics_update['path_conductance_mean'] = round(float(np.mean(conductances)), 6)
            if hop_areas:
                metrics_update['path_hop_area_mean'] = round(float(np.mean(hop_areas)), 4)
            if hop_area_mins:
                metrics_update['path_hop_area_min_mean'] = round(float(np.mean(hop_area_mins)), 4)

            # Re-read and update full_metrics.json
            metrics_path = os.path.join(output_dir, 'full_metrics.json')
            if os.path.exists(metrics_path):
                with open(metrics_path) as f:
                    existing = json.load(f)
                existing.update(metrics_update)
                with open(metrics_path, 'w') as f:
                    json.dump(existing, f, indent=2, default=str)

            # Update network_summary.csv placeholders with actual values
            ns_path = os.path.join(output_dir, 'network_summary.csv')
            if os.path.exists(ns_path):
                ns_df = pd.read_csv(ns_path)
                update_map = {}
                if gb_densities:
                    update_map['GB Density(hops/μm)'] = round(float(np.mean(gb_densities)), 3)
                if hop_areas:
                    update_map['Path Hop Area mean(μm²)'] = round(float(np.mean(hop_areas)), 4)
                if hop_area_mins:
                    update_map['Path Bottleneck(μm²)'] = round(float(np.mean(hop_area_mins)), 4)
                if conductances:
                    update_map['Path Conductance(μm²)'] = round(float(np.mean(conductances)), 6)
                for label, val in update_map.items():
                    mask = ns_df['지표'] == label
                    if mask.any():
                        ns_df.loc[mask, '값'] = str(val)
                ns_df.to_csv(ns_path, index=False)

    # Force chain data for 3D viewer
    force_chains = []
    fn_values = []
    for c in contacts_raw:
        fn = c.get('fn', 0) or np.sqrt(c.get('fn_x', 0)**2 + c.get('fn_y', 0)**2 + c.get('fn_z', 0)**2)
        fn_values.append(fn)
    fn_threshold = np.percentile(fn_values, 90) if fn_values else 0  # top 10%
    print(f"  Force chains: {len(fn_values)} contacts, threshold={fn_threshold:.6f}, max fn={max(fn_values) if fn_values else 0:.6f}")
    for c in contacts_raw:
        if c['id1'] in atoms_raw and c['id2'] in atoms_raw:
            fn = c.get('fn', 0) or np.sqrt(c.get('fn_x', 0)**2 + c.get('fn_y', 0)**2 + c.get('fn_z', 0)**2)
            if fn >= fn_threshold:
                a1, a2 = atoms_raw[c['id1']], atoms_raw[c['id2']]
                t1 = type_map.get(a1['type'], '?')
                t2 = type_map.get(a2['type'], '?')
                force_chains.append({
                    'p1': [a1['x'] * scale, a1['z'] * scale, a1['y'] * scale],
                    'p2': [a2['x'] * scale, a2['z'] * scale, a2['y'] * scale],
                    'fn': round(fn * 1e6 / scale**2, 3),  # μN
                    'type': '-'.join(sorted([t1, t2])),
                })
    print(f"  Force chains saved: {len(force_chains)} chains")
    with open(os.path.join(output_dir, 'force_chains.json'), 'w') as f:
        json.dump(force_chains, f)

    # Save tortuosity sample paths for 3D viewer (legacy)
    tau = results['tortuosity']
    if tau.get('mean') and 'graph' in perc:
        import networkx as nx
        G = perc['graph']
        reach_bottom = list(perc.get('bottom_se', set()) & perc.get('top_reachable_se', set()))
        reach_top = list(perc.get('top_se', set()) & perc.get('top_reachable_se', set()))
        paths_data = []
        for i in range(min(5, len(reach_bottom), len(reach_top))):
            src = reach_bottom[i % len(reach_bottom)]
            tgt = reach_top[i % len(reach_top)]
            try:
                path = nx.shortest_path(G, src, tgt, weight='distance')
                path_len = 0
                for ki in range(len(path)-1):
                    a1, a2 = atoms_raw[path[ki]], atoms_raw[path[ki+1]]
                    dx = abs(a1['x'] - a2['x'])
                    dy = abs(a1['y'] - a2['y'])
                    dz = a1['z'] - a2['z']
                    dx = min(dx, _box_x - dx)
                    dy = min(dy, _box_y - dy)
                    path_len += np.sqrt(dx**2 + dy**2 + dz**2)
                z_dist = abs(atoms_raw[tgt]['z'] - atoms_raw[src]['z'])
                if z_dist > 0:
                    paths_data.append({
                        'ids': path,
                        'tortuosity': round(path_len / z_dist, 2),
                        'path_length': round(path_len * scale, 1),
                        'z_distance': round(z_dist * scale, 1),
                    })
            except nx.NetworkXNoPath:
                pass
        with open(os.path.join(output_dir, 'tortuosity_paths.json'), 'w') as f:
            json.dump(paths_data, f)

    # Analyzed CSVs
    df_atom['type_name'] = df_atom['type'].map(type_map)
    coord_s = pd.Series(coord, name='coordination')
    df_atom = df_atom.merge(coord_s, left_on='id', right_index=True, how='left')
    df_atom['coordination'] = df_atom['coordination'].fillna(0).astype(int)
    df_atom.to_csv(os.path.join(output_dir, 'atoms_analyzed.csv'), index=False)

    id_to_type = {aid: type_map.get(a['type'], '?') for aid, a in atoms_raw.items()}
    df_contact['type1'] = df_contact['id1'].map(id_to_type)
    df_contact['type2'] = df_contact['id2'].map(id_to_type)
    df_contact['contact_type'] = df_contact.apply(
        lambda r: '-'.join(sorted([str(r['type1']), str(r['type2'])])), axis=1)
    df_contact.to_csv(os.path.join(output_dir, 'contacts_analyzed.csv'), index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('atoms_csv')
    parser.add_argument('contacts_csv')
    parser.add_argument('-o', '--output', default='./results')
    parser.add_argument('-t', '--type-map', default='1:AM,2:SE')
    parser.add_argument('-s', '--scale', type=float, default=1000)
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)
    type_map = {}
    for item in args.type_map.split(','):
        k, v = item.split(':')
        type_map[int(k)] = v.strip()

    print(f"Type map: {type_map}, Scale: {args.scale}")

    atoms_raw, df_atom = load_atoms_raw(args.atoms_csv)
    print(f"  {len(atoms_raw)} atoms")
    contacts_raw, df_contact = load_contacts_raw(args.contacts_csv)
    print(f"  {len(contacts_raw)} contacts")

    print("\n=== Full Analysis ===")
    results = run_full_analysis(atoms_raw, contacts_raw, type_map, args.scale, args.output)

    print("\nSaving...")
    save_results(results, atoms_raw, contacts_raw, df_atom, df_contact,
                 type_map, args.scale, args.output)
    print("Done!")


if __name__ == '__main__':
    main()
