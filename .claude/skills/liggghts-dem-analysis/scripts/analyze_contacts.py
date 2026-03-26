#!/usr/bin/env python3
"""
Contact Analysis for LIGGGHTS DEM Simulations
Analyzes contact networks for NCM-SE composite cathode systems.
"""

import numpy as np
import pandas as pd
import argparse
import json
import os
from collections import defaultdict

def calculate_force_magnitudes(df_contact):
    """Calculate force magnitudes from components."""
    df = df_contact.copy()

    # Total force magnitude
    df['F_total'] = np.sqrt(df['fx']**2 + df['fy']**2 + df['fz']**2)

    # Normal force magnitude
    df['F_normal'] = np.sqrt(df['fn_x']**2 + df['fn_y']**2 + df['fn_z']**2)

    # Tangential force magnitude
    df['F_tangential'] = np.sqrt(df['ft_x']**2 + df['ft_y']**2 + df['ft_z']**2)

    # Torque magnitude
    df['torque_mag'] = np.sqrt(df['torque_x']**2 + df['torque_y']**2 + df['torque_z']**2)

    return df


def classify_contacts(df_contact, df_atom, type_names=None):
    """
    Classify contacts by particle type pairs.

    Args:
        df_contact: Contact dataframe
        df_atom: Atom dataframe with 'id' and 'type' columns
        type_names: Optional dict mapping type numbers to names (e.g., {1: 'NCM', 2: 'SE'})
    """
    if type_names is None:
        type_names = {1: 'NCM', 2: 'SE'}

    id_to_type = dict(zip(df_atom['id'], df_atom['type']))

    contact_types = []
    for _, row in df_contact.iterrows():
        type1 = id_to_type.get(row['id1'], -1)
        type2 = id_to_type.get(row['id2'], -1)

        # Sort types for consistent naming
        types_sorted = tuple(sorted([type1, type2]))

        name1 = type_names.get(types_sorted[0], f'Type{types_sorted[0]}')
        name2 = type_names.get(types_sorted[1], f'Type{types_sorted[1]}')

        if types_sorted[0] == -1 or types_sorted[1] == -1:
            contact_types.append('Unknown')
        else:
            contact_types.append(f'{name1}-{name2}')

    df_contact['contact_type'] = contact_types
    return df_contact


def calculate_coordination_number(df_contact, df_atom):
    """Calculate coordination number for each particle."""
    coord_count = defaultdict(int)

    for _, row in df_contact.iterrows():
        coord_count[row['id1']] += 1
        coord_count[row['id2']] += 1

    df_atom['coordination'] = df_atom['id'].map(lambda x: coord_count.get(x, 0))
    return df_atom


def calculate_type_specific_coordination(df_contact, df_atom, type_names=None):
    """
    Calculate coordination number for specific contact types.
    E.g., NCM-SE coordination for NCM particles, SE-SE coordination for SE particles.
    """
    if type_names is None:
        type_names = {1: 'NCM', 2: 'SE'}

    id_to_type = dict(zip(df_atom['id'], df_atom['type']))

    # Initialize coordination counters for each type pair
    coord_by_type = defaultdict(lambda: defaultdict(int))

    for _, row in df_contact.iterrows():
        type1 = id_to_type.get(row['id1'], -1)
        type2 = id_to_type.get(row['id2'], -1)

        if type1 != -1 and type2 != -1:
            # Count for particle 1
            name2 = type_names.get(type2, f'Type{type2}')
            coord_by_type[row['id1']][f'coord_{name2}'] += 1

            # Count for particle 2
            name1 = type_names.get(type1, f'Type{type1}')
            coord_by_type[row['id2']][f'coord_{name1}'] += 1

    # Add to dataframe
    for coord_col in ['coord_NCM', 'coord_SE']:
        df_atom[coord_col] = df_atom['id'].map(lambda x: coord_by_type[x].get(coord_col, 0))

    return df_atom


def calculate_statistics(df_contact, df_atom, type_names=None):
    """Calculate comprehensive statistics."""
    if type_names is None:
        type_names = {1: 'NCM', 2: 'SE'}

    stats = {}

    # Total counts
    stats['total_atoms'] = len(df_atom)
    stats['total_contacts'] = len(df_contact)

    # Particle type distribution
    type_counts = df_atom['type'].value_counts().to_dict()
    for ptype, count in type_counts.items():
        name = type_names.get(ptype, f'Type{ptype}')
        stats[f'n_{name}'] = count

    # Contact type distribution
    contact_counts = df_contact['contact_type'].value_counts().to_dict()
    for ctype, count in contact_counts.items():
        stats[f'{ctype}_contacts'] = count
        stats[f'{ctype}_fraction'] = count / stats['total_contacts'] * 100

    # Coordination number statistics by type
    for ptype in df_atom['type'].unique():
        name = type_names.get(ptype, f'Type{ptype}')
        type_atoms = df_atom[df_atom['type'] == ptype]
        stats[f'{name}_avg_coordination'] = type_atoms['coordination'].mean()
        stats[f'{name}_max_coordination'] = type_atoms['coordination'].max()
        stats[f'{name}_min_coordination'] = type_atoms['coordination'].min()

    # Force statistics
    stats['avg_F_normal'] = df_contact['F_normal'].mean()
    stats['std_F_normal'] = df_contact['F_normal'].std()
    stats['max_F_normal'] = df_contact['F_normal'].max()
    stats['avg_F_tangential'] = df_contact['F_tangential'].mean()

    # Force statistics by contact type
    for ctype in df_contact['contact_type'].unique():
        ct_data = df_contact[df_contact['contact_type'] == ctype]
        stats[f'{ctype}_avg_F_normal'] = ct_data['F_normal'].mean()
        stats[f'{ctype}_avg_contact_area'] = ct_data['contact_area'].mean()

    # Contact area statistics
    stats['avg_contact_area'] = df_contact['contact_area'].mean()
    stats['total_contact_area'] = df_contact['contact_area'].sum()

    # Overlap (delta) statistics
    if 'delta' in df_contact.columns:
        stats['avg_delta'] = df_contact['delta'].mean()
        stats['max_delta'] = df_contact['delta'].max()

    return stats


def calculate_ionic_conductivity_metrics(df_contact, df_atom, type_names=None):
    """
    Calculate metrics relevant to ionic conductivity estimation.
    """
    if type_names is None:
        type_names = {1: 'NCM', 2: 'SE'}

    metrics = {}

    # Get type IDs
    ncm_type = [k for k, v in type_names.items() if v == 'NCM'][0] if 'NCM' in type_names.values() else 1
    se_type = [k for k, v in type_names.items() if v == 'SE'][0] if 'SE' in type_names.values() else 2

    ncm_ids = set(df_atom[df_atom['type'] == ncm_type]['id'])
    se_ids = set(df_atom[df_atom['type'] == se_type]['id'])

    # NCM-SE contacts
    ncm_se_contacts = df_contact[df_contact['contact_type'] == 'NCM-SE']
    se_se_contacts = df_contact[df_contact['contact_type'] == 'SE-SE']

    # Total NCM-SE interface area
    metrics['total_NCM_SE_area'] = ncm_se_contacts['contact_area'].sum()
    metrics['avg_NCM_SE_area'] = ncm_se_contacts['contact_area'].mean() if len(ncm_se_contacts) > 0 else 0

    # NCM-SE contact area per NCM particle
    ncm_contact_area = defaultdict(float)
    for _, row in ncm_se_contacts.iterrows():
        if row['id1'] in ncm_ids:
            ncm_contact_area[row['id1']] += row['contact_area']
        if row['id2'] in ncm_ids:
            ncm_contact_area[row['id2']] += row['contact_area']

    # NCM particles with SE contact
    ncm_with_contact = len([pid for pid in ncm_ids if ncm_contact_area.get(pid, 0) > 0])
    metrics['NCM_with_SE_contact'] = ncm_with_contact
    metrics['NCM_with_SE_contact_fraction'] = ncm_with_contact / len(ncm_ids) * 100 if len(ncm_ids) > 0 else 0

    # SE-SE coordination (for percolation analysis)
    se_se_coord = defaultdict(int)
    for _, row in se_se_contacts.iterrows():
        if row['id1'] in se_ids:
            se_se_coord[row['id1']] += 1
        if row['id2'] in se_ids:
            se_se_coord[row['id2']] += 1

    coord_values = [se_se_coord.get(pid, 0) for pid in se_ids]
    metrics['SE_SE_avg_coordination'] = np.mean(coord_values) if coord_values else 0

    # Percolation analysis (threshold typically ~2.4 for 3D)
    percolation_threshold = 2.4
    se_connected = sum(1 for c in coord_values if c > 0)
    se_percolated = sum(1 for c in coord_values if c >= percolation_threshold)

    metrics['SE_connected'] = se_connected
    metrics['SE_connected_fraction'] = se_connected / len(se_ids) * 100 if len(se_ids) > 0 else 0
    metrics['SE_percolated'] = se_percolated
    metrics['SE_percolated_fraction'] = se_percolated / len(se_ids) * 100 if len(se_ids) > 0 else 0
    metrics['percolation_threshold'] = percolation_threshold

    return metrics


def main():
    parser = argparse.ArgumentParser(description='Analyze LIGGGHTS contact data')
    parser.add_argument('atom_file', help='Path to atoms CSV file')
    parser.add_argument('contact_file', help='Path to contacts CSV file')
    parser.add_argument('--output-dir', '-o', default='./results',
                        help='Output directory')
    parser.add_argument('--type-names', '-t', type=str, default='1:NCM,2:SE',
                        help='Type name mapping (e.g., "1:NCM,2:SE")')
    parser.add_argument('--scale-factor', '-s', type=float, default=1000,
                        help='Scale factor for force/stress (default: 1000)')
    args = parser.parse_args()

    # Parse type names
    type_names = {}
    for pair in args.type_names.split(','):
        k, v = pair.split(':')
        type_names[int(k)] = v

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    print("="*60)
    print("LIGGGHTS Contact Analysis")
    print("="*60)

    # Load data
    print(f"\n[1] Loading data...")
    df_atom = pd.read_csv(args.atom_file)
    df_contact = pd.read_csv(args.contact_file)
    print(f"    Atoms: {len(df_atom)}")
    print(f"    Contacts: {len(df_contact)}")

    # Calculate force magnitudes
    print(f"\n[2] Calculating force magnitudes...")
    df_contact = calculate_force_magnitudes(df_contact)

    # Classify contacts
    print(f"\n[3] Classifying contacts...")
    df_contact = classify_contacts(df_contact, df_atom, type_names)
    contact_dist = df_contact['contact_type'].value_counts()
    for ctype, count in contact_dist.items():
        print(f"    {ctype}: {count} ({100*count/len(df_contact):.1f}%)")

    # Calculate coordination numbers
    print(f"\n[4] Calculating coordination numbers...")
    df_atom = calculate_coordination_number(df_contact, df_atom)
    df_atom = calculate_type_specific_coordination(df_contact, df_atom, type_names)

    for ptype in df_atom['type'].unique():
        name = type_names.get(ptype, f'Type{ptype}')
        avg_coord = df_atom[df_atom['type'] == ptype]['coordination'].mean()
        print(f"    {name} avg coordination: {avg_coord:.2f}")

    # Calculate statistics
    print(f"\n[5] Computing statistics...")
    stats = calculate_statistics(df_contact, df_atom, type_names)

    # Calculate ionic conductivity metrics
    print(f"\n[6] Computing ionic conductivity metrics...")
    ionic_metrics = calculate_ionic_conductivity_metrics(df_contact, df_atom, type_names)

    print(f"    SE-SE avg coordination: {ionic_metrics['SE_SE_avg_coordination']:.2f}")
    print(f"    SE percolated: {ionic_metrics['SE_percolated_fraction']:.1f}%")
    print(f"    NCM with SE contact: {ionic_metrics['NCM_with_SE_contact_fraction']:.1f}%")

    # Save results
    print(f"\n[7] Saving results...")

    # Save processed dataframes
    df_atom.to_csv(os.path.join(args.output_dir, 'atoms_analyzed.csv'), index=False)
    df_contact.to_csv(os.path.join(args.output_dir, 'contacts_analyzed.csv'), index=False)

    # Save statistics
    all_stats = {
        'basic_stats': stats,
        'ionic_metrics': ionic_metrics,
        'scale_factor': args.scale_factor,
        'type_names': type_names
    }
    with open(os.path.join(args.output_dir, 'analysis_stats.json'), 'w') as f:
        json.dump(all_stats, f, indent=2, default=float)

    print(f"    Saved to {args.output_dir}/")

    # Print summary
    print("\n" + "="*60)
    print("ANALYSIS SUMMARY")
    print("="*60)

    print(f"\n--- Particle Statistics ---")
    for ptype in df_atom['type'].unique():
        name = type_names.get(ptype, f'Type{ptype}')
        count = len(df_atom[df_atom['type'] == ptype])
        print(f"{name}: {count} particles")

    print(f"\n--- Contact Statistics ---")
    print(f"Total contacts: {stats['total_contacts']}")
    for ctype in df_contact['contact_type'].unique():
        count = stats.get(f'{ctype}_contacts', 0)
        frac = stats.get(f'{ctype}_fraction', 0)
        print(f"{ctype}: {count} ({frac:.1f}%)")

    print(f"\n--- Force Statistics (×{args.scale_factor} for real values) ---")
    print(f"Avg normal force: {stats['avg_F_normal']*args.scale_factor:.4f} N")
    print(f"Max normal force: {stats['max_F_normal']*args.scale_factor:.4f} N")

    print(f"\n--- Ionic Conductivity Metrics ---")
    print(f"SE-SE coordination: {ionic_metrics['SE_SE_avg_coordination']:.2f} (threshold: {ionic_metrics['percolation_threshold']})")
    print(f"SE percolated: {ionic_metrics['SE_percolated']}/{len(df_atom[df_atom['type']==2])} ({ionic_metrics['SE_percolated_fraction']:.1f}%)")
    print(f"NCM-SE interface: {ionic_metrics['NCM_with_SE_contact_fraction']:.1f}% NCM active")

    return df_atom, df_contact, all_stats


if __name__ == "__main__":
    main()
