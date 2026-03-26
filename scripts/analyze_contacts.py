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
    for ct, v in results['interface'].items():
        rows.append({
            '접촉유형': ct, '접촉수': v['n_contacts'],
            '접촉면적_mean(μm²)': round(v['mean_area'], 4),
            '접촉면적_total(μm²)': round(v['total_area'], 2),
        })
    pd.DataFrame(rows).to_csv(os.path.join(output_dir, 'contact_summary.csv'), index=False)

    # Atom Statistics
    rows = []
    for t, name in type_map.items():
        sub = {aid: a for aid, a in atoms_raw.items() if a['type'] == t}
        if not sub: continue
        zs = [a['z'] for a in sub.values()]
        rs = [a['radius'] for a in sub.values()]
        rows.append({
            '입자유형': name, '입자수': len(sub),
            '반지름(μm)': round(np.mean(rs) * scale, 2),
            'Z_min(μm)': round(min(zs) * scale, 1),
            'Z_max(μm)': round(max(zs) * scale, 1),
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

    # Force Summary
    force_by_type = defaultdict(list)
    for c in contacts_raw:
        if c['id1'] in atoms_raw and c['id2'] in atoms_raw:
            t1 = type_map.get(atoms_raw[c['id1']]['type'], '?')
            t2 = type_map.get(atoms_raw[c['id2']]['type'], '?')
            ct = '-'.join(sorted([t1, t2]))
            force_by_type[ct].append(c['fn'])
    rows = []
    for ct, forces in force_by_type.items():
        f = np.array(forces) * scale
        rows.append({
            '접촉유형': ct,
            '수직력_mean(μN)': round(float(np.mean(f)), 2),
            '수직력_median(μN)': round(float(np.median(f)), 2),
            '수직력_max(μN)': round(float(np.max(f)), 2),
        })
    pd.DataFrame(rows).to_csv(os.path.join(output_dir, 'force_summary.csv'), index=False)

    # Network Summary (comprehensive)
    perc = results['percolation']
    cn = results['se_se_cn']
    tau = results['tortuosity']
    ionic = results['ionic_active']
    rows = [
        {'지표': 'Porosity(%)', '값': round(results['porosity'], 2)},
        {'지표': '전극두께(μm)', '값': round(results['thickness_um'], 2)},
        {'지표': '두께기준', '값': results['plate_z_source']},
        {'지표': 'SE-SE CN mean', '값': round(cn['mean'], 2)},
        {'지표': 'SE Percolation(%)', '값': round(perc['percolation_pct'], 1)},
        {'지표': 'Top Reachable(%)', '값': round(perc['top_reachable_pct'], 1)},
        {'지표': '네트워크 수', '값': perc['n_components']},
        {'지표': 'Tortuosity mean', '값': round(tau['mean'], 2) if tau['mean'] else 'N/A'},
        {'지표': 'Tortuosity std', '값': round(tau['std'], 2) if tau['std'] else 'N/A'},
        {'지표': 'Ionic Active AM(%)', '값': round(ionic['active_pct'], 1)},
    ]
    for lbl, v in results['coverage'].items():
        rows.append({'지표': f'Coverage {lbl}(%)', '값': round(v['mean'], 1)})
        rows.append({'지표': f'Coverage {lbl} std', '값': round(v['std'], 1)})
    pd.DataFrame(rows).to_csv(os.path.join(output_dir, 'network_summary.csv'), index=False)

    # Full metrics JSON
    metrics = {
        'porosity': results['porosity'],
        'thickness_um': results['thickness_um'],
        'plate_z_source': results['plate_z_source'],
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
