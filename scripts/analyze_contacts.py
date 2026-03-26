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
        atoms[int(row['id'])] = {
            'type': int(row['type']),
            'x': row['x'], 'y': row['y'], 'z': row['z'],
            'radius': row['radius'],
        }
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
            'ft': np.sqrt(row['ft_x']**2 + row['ft_y']**2 + row['ft_z']**2),
            'contact_area': row['contact_area'],
            'delta': row['delta'],
        })
    return contacts, df


def save_results(results, atoms_raw, contacts_raw, df_atom, df_contact,
                 type_map, scale, output_dir):
    area_conv = 1.0 / (scale ** 2) * 1e12

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
    perc = results['percolation']
    cn = results['se_se_cn']
    tau = results['tortuosity']
    ionic = results['ionic_active']
    # 구조 → 계면 → 이온경로 → 활성도 순서
    rows = [
        # 1. 구조 (전극 기본 정보)
        {'지표': 'Porosity(%)', '값': round(results['porosity'], 2)},
        {'지표': '전극두께(μm)', '값': round(results['thickness_um'], 2)},
        # 2. 계면 (반응 면적)
        {'지표': 'AM-SE Total(μm²)', '값': round(results['interface'].get('AM전체-SE', {}).get('total_area', 0), 2)},
        {'지표': 'SE-SE Total(μm²)', '값': round(results['interface'].get('SE-SE', {}).get('total_area', 0), 2)},
    ]
    for lbl, v in results['coverage'].items():
        rows.append({'지표': f'Coverage {lbl}(%)', '값': round(v['mean'], 1)})
    rows += [
        # 3. 이온경로 (SE 네트워크)
        {'지표': 'SE-SE CN mean', '값': round(cn['mean'], 2)},
        {'지표': 'SE Cluster 수', '값': perc['n_components']},
        {'지표': 'SE Percolation(%)', '값': round(perc['percolation_pct'], 1)},
        {'지표': 'Top Reachable(%)', '값': round(perc['top_reachable_pct'], 1)},
        {'지표': 'Tortuosity mean', '값': round(tau['mean'], 2) if tau['mean'] else 'N/A'},
        {'지표': 'Tortuosity std', '값': round(tau['std'], 2) if tau['std'] else 'N/A'},
        # 4. 활성도 (최종 성능)
        {'지표': 'Ionic Active AM(%)', '값': round(ionic['active_pct'], 1)},
    ]
    pd.DataFrame(rows).to_csv(os.path.join(output_dir, 'network_summary.csv'), index=False)

    # Auto-detect P:S ratio from mass (count × volume × density)
    # NCM811 density = 4800 kg/m³ for both AM_P and AM_S
    am_density = 4800  # kg/m³ (NCM811)
    am_p_atoms = [a for a in atoms_raw.values() if type_map.get(a['type']) == 'AM_P']
    am_s_atoms = [a for a in atoms_raw.values() if type_map.get(a['type']) == 'AM_S']

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

    # Full metrics JSON
    metrics = {
        'porosity': results['porosity'],
        'thickness_um': results['thickness_um'],
        'plate_z_source': results['plate_z_source'],
        'ps_ratio': ps_ratio,
        'se_se_cn': cn['mean'],
        'percolation_pct': perc['percolation_pct'],
        'top_reachable_pct': perc['top_reachable_pct'],
        'n_components': perc['n_components'],
        'tortuosity_mean': tau['mean'],
        'tortuosity_std': tau['std'],
        'ionic_active_pct': ionic['active_pct'],
    }
    for ct, v in results['interface'].items():
        safe = ct.replace('-', '_')
        metrics[f'area_{safe}_total'] = v['total_area']
        metrics[f'area_{safe}_n'] = v['n_contacts']
        metrics[f'area_{safe}_mean'] = v['mean_area']
    for lbl, v in results['coverage'].items():
        metrics[f'coverage_{lbl}_mean'] = v['mean']
        metrics[f'coverage_{lbl}_std'] = v['std']
    with open(os.path.join(output_dir, 'full_metrics.json'), 'w') as f:
        json.dump(metrics, f, indent=2, default=str)

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
