#!/usr/bin/env python3
"""
Analyze contacts from parsed LIGGGHTS CSV files.
Classifies contacts, computes coordination numbers, force statistics, network metrics.

Usage:
    python analyze_contacts.py results/atoms.csv results/contacts.csv -o ./results -t "1:AM,2:SE" -s 1000
"""
import argparse
import os
import sys
import numpy as np
import pandas as pd
import json

def main():
    parser = argparse.ArgumentParser(description='Analyze DEM contacts')
    parser.add_argument('atoms_csv', help='Path to atoms.csv')
    parser.add_argument('contacts_csv', help='Path to contacts.csv')
    parser.add_argument('-o', '--output', default='./results')
    parser.add_argument('-t', '--type-map', default='1:AM,2:SE',
                        help='Type mapping, e.g. "1:AM,2:SE" or "1:AM_P,2:AM_S,3:SE"')
    parser.add_argument('-s', '--scale', type=float, default=1000,
                        help='Scale factor (sim mm -> real um)')
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)
    scale = args.scale

    # Parse type map
    type_map = {}
    for item in args.type_map.split(','):
        k, v = item.split(':')
        type_map[int(k)] = v.strip()
    print(f"Type map: {type_map}")

    # ─── Load Data ──────────────────────────────────────────────────────────
    print("Loading atoms...")
    df_atom = pd.read_csv(args.atoms_csv)
    # Ensure numeric
    for col in df_atom.columns:
        if col not in ['id', 'type']:
            df_atom[col] = pd.to_numeric(df_atom[col], errors='coerce')
    df_atom['id'] = df_atom['id'].astype(int)
    df_atom['type'] = df_atom['type'].astype(int)
    df_atom['type_name'] = df_atom['type'].map(type_map)
    print(f"  Atoms: {len(df_atom)}")
    for t, name in type_map.items():
        n = (df_atom['type'] == t).sum()
        print(f"    Type {t} ({name}): {n}")

    print("Loading contacts...")
    df_contact = pd.read_csv(args.contacts_csv)
    for col in df_contact.columns:
        df_contact[col] = pd.to_numeric(df_contact[col], errors='coerce')
    df_contact['id1'] = df_contact['id1'].astype(int)
    df_contact['id2'] = df_contact['id2'].astype(int)
    print(f"  Contacts: {len(df_contact)}")

    # ─── Scale Conversion ───────────────────────────────────────────────────
    # Positions: sim (mm-scale) -> real (um)
    pos_cols_atom = ['x', 'y', 'z']
    for c in pos_cols_atom:
        if c in df_atom.columns:
            df_atom[c] = df_atom[c] * scale
    if 'radius' in df_atom.columns:
        df_atom['radius'] = df_atom['radius'] * scale

    pos_cols_contact = ['p1_x', 'p1_y', 'p1_z', 'p2_x', 'p2_y', 'p2_z',
                        'cp_x', 'cp_y', 'cp_z']
    for c in pos_cols_contact:
        if c in df_contact.columns:
            df_contact[c] = df_contact[c] * scale

    if 'delta' in df_contact.columns:
        df_contact['delta'] = df_contact['delta'] * scale
    if 'contact_area' in df_contact.columns:
        df_contact['contact_area'] = df_contact['contact_area'] * (scale ** 2)

    # Force: scale factor
    force_cols = ['fx', 'fy', 'fz', 'fn_x', 'fn_y', 'fn_z', 'ft_x', 'ft_y', 'ft_z']
    for c in force_cols:
        if c in df_contact.columns:
            df_contact[c] = df_contact[c] * scale

    # Stress: already in simulation units, scale to MPa
    stress_cols = ['c_strs[1]', 'c_strs[2]', 'c_strs[3]']
    for c in stress_cols:
        if c in df_atom.columns:
            # LIGGGHTS stress/atom is in Pa*volume, need to divide by volume
            # For per-atom stress: σ = strs / volume, volume = 4/3 π r³
            pass  # We'll compute per-atom stress below

    # Compute per-atom stress (σ = stress_atom / volume)
    if 'radius' in df_atom.columns and all(c in df_atom.columns for c in stress_cols):
        vol = (4.0 / 3.0) * np.pi * (df_atom['radius'] / scale) ** 3  # volume in sim units
        for c in stress_cols:
            col_name = c.replace('c_strs[1]', 'sigma_xx').replace('c_strs[2]', 'sigma_yy').replace('c_strs[3]', 'sigma_zz')
            # stress/atom output is already stress*volume in LIGGGHTS
            # divide by volume to get stress, then convert units
            df_atom[col_name] = df_atom[c] / vol / 1e6  # Pa -> MPa

    # ─── Contact Classification ─────────────────────────────────────────────
    print("Classifying contacts...")
    id_to_type = df_atom.set_index('id')['type_name'].to_dict()

    df_contact['type1'] = df_contact['id1'].map(id_to_type)
    df_contact['type2'] = df_contact['id2'].map(id_to_type)

    # Standardize contact type (alphabetical order)
    def classify_contact(row):
        t1, t2 = sorted([str(row['type1']), str(row['type2'])])
        return f"{t1}-{t2}"

    df_contact['contact_type'] = df_contact.apply(classify_contact, axis=1)
    contact_types = df_contact['contact_type'].value_counts()
    print("  Contact type distribution:")
    for ct, n in contact_types.items():
        print(f"    {ct}: {n} ({100*n/len(df_contact):.1f}%)")

    # ─── Computed Columns ───────────────────────────────────────────────────
    # Normal force magnitude
    if all(c in df_contact.columns for c in ['fn_x', 'fn_y', 'fn_z']):
        df_contact['fn_mag'] = np.sqrt(
            df_contact['fn_x']**2 + df_contact['fn_y']**2 + df_contact['fn_z']**2)

    # Tangential force magnitude
    if all(c in df_contact.columns for c in ['ft_x', 'ft_y', 'ft_z']):
        df_contact['ft_mag'] = np.sqrt(
            df_contact['ft_x']**2 + df_contact['ft_y']**2 + df_contact['ft_z']**2)

    # Total force magnitude
    if all(c in df_contact.columns for c in ['fx', 'fy', 'fz']):
        df_contact['f_mag'] = np.sqrt(
            df_contact['fx']**2 + df_contact['fy']**2 + df_contact['fz']**2)

    # Ft/Fn ratio
    if 'fn_mag' in df_contact.columns and 'ft_mag' in df_contact.columns:
        df_contact['ft_fn_ratio'] = df_contact['ft_mag'] / df_contact['fn_mag'].replace(0, np.nan)

    # ─── Coordination Number ────────────────────────────────────────────────
    print("Computing coordination numbers...")
    contacts_per_atom = pd.concat([
        df_contact[['id1']].rename(columns={'id1': 'id'}),
        df_contact[['id2']].rename(columns={'id2': 'id'})
    ])
    coord_num = contacts_per_atom.groupby('id').size().rename('coordination')
    df_atom = df_atom.merge(coord_num, left_on='id', right_index=True, how='left')
    df_atom['coordination'] = df_atom['coordination'].fillna(0).astype(int)

    # Type-specific coordination (e.g., SE-SE contacts for SE particles)
    type_names = list(type_map.values())
    for tn in type_names:
        # Contacts where this type is involved
        mask1 = df_contact['type1'] == tn
        mask2 = df_contact['type2'] == tn
        same_contacts = df_contact[(mask1 & (df_contact['type2'] == tn)) |
                                    (mask2 & (df_contact['type1'] == tn))]
        ids = pd.concat([same_contacts['id1'], same_contacts['id2']])
        ids = ids[ids.isin(df_atom[df_atom['type_name'] == tn]['id'])]
        same_coord = ids.groupby(ids).size().rename(f'coord_{tn}_{tn}')
        df_atom = df_atom.merge(same_coord, left_on='id', right_index=True, how='left')
        df_atom[f'coord_{tn}_{tn}'] = df_atom[f'coord_{tn}_{tn}'].fillna(0).astype(int)

    coord_stats = df_atom.groupby('type_name')['coordination'].agg(['mean', 'std', 'min', 'max'])
    print("  Coordination by type:")
    print(coord_stats.to_string())

    # ─── Network Analysis ───────────────────────────────────────────────────
    print("Analyzing networks...")
    import networkx as nx

    network_results = {}
    for tn in type_names:
        # Build graph for same-type contacts
        type_ids = set(df_atom[df_atom['type_name'] == tn]['id'])
        same_contacts = df_contact[
            (df_contact['type1'] == tn) & (df_contact['type2'] == tn)
        ]
        G = nx.Graph()
        G.add_nodes_from(type_ids)
        edges = list(zip(same_contacts['id1'].values, same_contacts['id2'].values))
        G.add_edges_from(edges)

        components = list(nx.connected_components(G))
        largest = max(components, key=len) if components else set()
        n_total = len(type_ids)
        n_connected = sum(1 for n in type_ids if G.degree(n) > 0)
        n_largest = len(largest)

        # Percolation check (does largest component span Z range?)
        if n_largest > 0:
            z_vals = df_atom[df_atom['id'].isin(largest)]['z']
            z_range = z_vals.max() - z_vals.min()
            total_z = df_atom['z'].max() - df_atom['z'].min()
            percolation = z_range / total_z * 100 if total_z > 0 else 0
        else:
            percolation = 0

        network_results[tn] = {
            'total': n_total,
            'connected': n_connected,
            'connected_pct': 100 * n_connected / n_total if n_total > 0 else 0,
            'components': len(components),
            'largest_pct': 100 * n_largest / n_total if n_total > 0 else 0,
            'percolation_pct': percolation,
        }
        print(f"  {tn}: {n_connected}/{n_total} connected ({100*n_connected/n_total:.1f}%), "
              f"{len(components)} components, largest {100*n_largest/n_total:.1f}%, "
              f"percolation {percolation:.1f}%")

    # ─── Save Results ───────────────────────────────────────────────────────
    print("Saving results...")

    # Analyzed CSVs
    df_atom.to_csv(os.path.join(args.output, 'atoms_analyzed.csv'), index=False)
    df_contact.to_csv(os.path.join(args.output, 'contacts_analyzed.csv'), index=False)

    # Contact summary
    contact_summary = []
    for ct in df_contact['contact_type'].unique():
        subset = df_contact[df_contact['contact_type'] == ct]
        row = {
            '접촉유형': ct,
            '접촉수': len(subset),
            '비율(%)': f"{100*len(subset)/len(df_contact):.1f}",
        }
        if 'fn_mag' in subset.columns:
            row['수직력_mean(μN)'] = f"{subset['fn_mag'].mean():.2f}"
            row['수직력_max(μN)'] = f"{subset['fn_mag'].max():.2f}"
        if 'contact_area' in subset.columns:
            row['접촉면적_mean(μm²)'] = f"{subset['contact_area'].mean():.4f}"
            row['접촉면적_total(μm²)'] = f"{subset['contact_area'].sum():.2f}"
        if 'delta' in subset.columns:
            row['겹침량_mean(μm)'] = f"{subset['delta'].mean():.4f}"
        contact_summary.append(row)
    pd.DataFrame(contact_summary).to_csv(
        os.path.join(args.output, 'contact_summary.csv'), index=False)

    # Coordination summary
    coord_summary = []
    for tn in type_names:
        sub = df_atom[df_atom['type_name'] == tn]
        coord_summary.append({
            '입자유형': tn,
            '입자수': len(sub),
            '배위수_mean': f"{sub['coordination'].mean():.1f}",
            '배위수_std': f"{sub['coordination'].std():.1f}",
            '배위수_min': sub['coordination'].min(),
            '배위수_max': sub['coordination'].max(),
        })
    pd.DataFrame(coord_summary).to_csv(
        os.path.join(args.output, 'coordination_summary.csv'), index=False)

    # Network summary
    net_rows = []
    for tn, data in network_results.items():
        net_rows.append({
            '입자유형': tn,
            '전체수': data['total'],
            '연결된수': data['connected'],
            '연결비율(%)': f"{data['connected_pct']:.1f}",
            '네트워크수': data['components'],
            '최대네트워크(%)': f"{data['largest_pct']:.1f}",
            'Percolation(%)': f"{data['percolation_pct']:.1f}",
        })
    pd.DataFrame(net_rows).to_csv(
        os.path.join(args.output, 'network_summary.csv'), index=False)

    # Force summary
    force_summary = []
    for ct in df_contact['contact_type'].unique():
        subset = df_contact[df_contact['contact_type'] == ct]
        row = {'접촉유형': ct}
        if 'fn_mag' in subset.columns:
            row['수직력_mean(μN)'] = f"{subset['fn_mag'].mean():.2f}"
            row['수직력_median(μN)'] = f"{subset['fn_mag'].median():.2f}"
            row['수직력_max(μN)'] = f"{subset['fn_mag'].max():.2f}"
        if 'ft_mag' in subset.columns:
            row['전단력_mean(μN)'] = f"{subset['ft_mag'].mean():.2f}"
        if 'ft_fn_ratio' in subset.columns:
            valid = subset['ft_fn_ratio'].dropna()
            row['전단/수직_median'] = f"{valid.median():.3f}" if len(valid) > 0 else 'N/A'
        force_summary.append(row)
    pd.DataFrame(force_summary).to_csv(
        os.path.join(args.output, 'force_summary.csv'), index=False)

    # Atom statistics
    atom_stats = []
    for tn in type_names:
        sub = df_atom[df_atom['type_name'] == tn]
        row = {
            '입자유형': tn,
            '입자수': len(sub),
            '반지름(μm)': f"{sub['radius'].mean():.2f}" if 'radius' in sub.columns else 'N/A',
            'Z_min(μm)': f"{sub['z'].min():.1f}",
            'Z_max(μm)': f"{sub['z'].max():.1f}",
        }
        if 'sigma_zz' in sub.columns:
            row['σ_zz_mean(MPa)'] = f"{sub['sigma_zz'].mean():.1f}"
        atom_stats.append(row)
    pd.DataFrame(atom_stats).to_csv(
        os.path.join(args.output, 'atom_statistics.csv'), index=False)

    print(f"\nAll results saved to {args.output}")
    print("Done!")


if __name__ == '__main__':
    main()
