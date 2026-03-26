#!/usr/bin/env python3
"""
Advanced DEM Analysis for LIGGGHTS Simulations
Includes: Stress Tensor, Force Chain, Fabric Tensor, Porosity, Network Analysis
"""

import numpy as np
import pandas as pd
import argparse
import json
import os
from collections import defaultdict
from scipy.spatial import Voronoi, ConvexHull
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import shortest_path, connected_components
import warnings
warnings.filterwarnings('ignore')


# =============================================================================
# 1. STRESS TENSOR ANALYSIS
# =============================================================================

def analyze_stress_tensor(df_atom, scale_factor=1000):
    """
    Analyze full stress tensor: Von Mises, Hydrostatic, Anisotropy

    Stress columns: c_strs[1-6] = σ_xx, σ_yy, σ_zz, σ_xy, σ_xz, σ_yz
    """
    results = {}

    # Check for stress columns
    stress_cols = [f'c_strs[{i}]' for i in range(1, 7)]
    if not all(col in df_atom.columns for col in stress_cols):
        print("    Warning: Stress tensor columns not found")
        return None

    # Extract stress components (scaled)
    s_xx = df_atom['c_strs[1]'].values * scale_factor
    s_yy = df_atom['c_strs[2]'].values * scale_factor
    s_zz = df_atom['c_strs[3]'].values * scale_factor
    s_xy = df_atom['c_strs[4]'].values * scale_factor
    s_xz = df_atom['c_strs[5]'].values * scale_factor
    s_yz = df_atom['c_strs[6]'].values * scale_factor

    # Hydrostatic pressure (mean stress)
    hydrostatic = (s_xx + s_yy + s_zz) / 3
    df_atom['hydrostatic_pressure'] = hydrostatic

    # Von Mises stress
    von_mises = np.sqrt(0.5 * ((s_xx - s_yy)**2 + (s_yy - s_zz)**2 +
                                (s_zz - s_xx)**2 + 6*(s_xy**2 + s_xz**2 + s_yz**2)))
    df_atom['von_mises_stress'] = von_mises

    # Deviatoric stress components
    s_dev_xx = s_xx - hydrostatic
    s_dev_yy = s_yy - hydrostatic
    s_dev_zz = s_zz - hydrostatic

    # Stress anisotropy (ratio of axial to lateral)
    lateral_avg = (s_xx + s_yy) / 2
    anisotropy = np.where(np.abs(lateral_avg) > 1e-10,
                          s_zz / lateral_avg, 0)
    df_atom['stress_anisotropy'] = anisotropy

    # Stress invariants
    I1 = s_xx + s_yy + s_zz  # First invariant
    I2 = s_xx*s_yy + s_yy*s_zz + s_zz*s_xx - s_xy**2 - s_xz**2 - s_yz**2

    # Statistics by particle type
    for ptype in df_atom['type'].unique():
        mask = df_atom['type'] == ptype
        prefix = f'type{ptype}'

        results[f'{prefix}_von_mises_avg'] = float(von_mises[mask].mean())
        results[f'{prefix}_von_mises_max'] = float(von_mises[mask].max())
        results[f'{prefix}_von_mises_std'] = float(von_mises[mask].std())
        results[f'{prefix}_hydrostatic_avg'] = float(hydrostatic[mask].mean())
        results[f'{prefix}_anisotropy_avg'] = float(np.nanmean(anisotropy[mask]))

    # Overall statistics
    results['von_mises_avg'] = float(von_mises.mean())
    results['von_mises_max'] = float(von_mises.max())
    results['von_mises_std'] = float(von_mises.std())
    results['hydrostatic_avg'] = float(hydrostatic.mean())
    results['stress_uniformity'] = float(1 - von_mises.std() / (von_mises.mean() + 1e-10))

    # Stress heterogeneity index
    results['heterogeneity_index'] = float(von_mises.std() / von_mises.mean()) if von_mises.mean() > 0 else 0

    return results


# =============================================================================
# 2. FORCE CHAIN ANALYSIS
# =============================================================================

def analyze_force_chains(df_contact, df_atom, percentile=90):
    """
    Identify and analyze force chains (high-force contact networks)
    """
    results = {}

    if 'F_normal' not in df_contact.columns:
        df_contact['F_normal'] = np.sqrt(
            df_contact['fn_x']**2 + df_contact['fn_y']**2 + df_contact['fn_z']**2
        )

    # Force threshold (top percentile)
    force_threshold = np.percentile(df_contact['F_normal'], percentile)
    results['force_threshold'] = float(force_threshold)
    results['percentile_used'] = percentile

    # High-force contacts
    high_force = df_contact[df_contact['F_normal'] >= force_threshold]
    results['n_high_force_contacts'] = len(high_force)
    results['high_force_fraction'] = len(high_force) / len(df_contact) * 100

    # Force chain direction analysis
    if len(high_force) > 0:
        # Calculate force direction vectors
        fx = high_force['fn_x'].values
        fy = high_force['fn_y'].values
        fz = high_force['fn_z'].values
        f_mag = np.sqrt(fx**2 + fy**2 + fz**2)

        # Normalized directions
        nx = fx / (f_mag + 1e-10)
        ny = fy / (f_mag + 1e-10)
        nz = fz / (f_mag + 1e-10)

        # Z-direction preference (compression direction)
        z_alignment = np.abs(nz).mean()
        results['z_alignment'] = float(z_alignment)

        # Horizontal vs vertical ratio
        horizontal = np.sqrt(nx**2 + ny**2).mean()
        results['horizontal_component'] = float(horizontal)
        results['vertical_component'] = float(np.abs(nz).mean())

        # Force chain connectivity
        chain_particles = set(high_force['id1'].values) | set(high_force['id2'].values)
        results['n_chain_particles'] = len(chain_particles)
        results['chain_particle_fraction'] = len(chain_particles) / len(df_atom) * 100

    # Build force chain graph
    chain_contacts = high_force[['id1', 'id2', 'F_normal']].copy()
    results['force_chain_contacts'] = chain_contacts.to_dict('records')

    return results


# =============================================================================
# 3. FABRIC TENSOR ANALYSIS
# =============================================================================

def analyze_fabric_tensor(df_contact, df_atom):
    """
    Calculate fabric tensor to quantify structural anisotropy
    """
    results = {}

    # Contact branch vectors (center to center)
    x1, y1, z1 = df_contact['x1'].values, df_contact['y1'].values, df_contact['z1'].values
    x2, y2, z2 = df_contact['x2'].values, df_contact['y2'].values, df_contact['z2'].values

    # Branch vectors
    bx = x2 - x1
    by = y2 - y1
    bz = z2 - z1
    b_mag = np.sqrt(bx**2 + by**2 + bz**2)

    # Normalized branch vectors
    nx = bx / (b_mag + 1e-10)
    ny = by / (b_mag + 1e-10)
    nz = bz / (b_mag + 1e-10)

    # Fabric tensor components (symmetric 3x3)
    N = len(df_contact)
    F_xx = np.sum(nx * nx) / N
    F_yy = np.sum(ny * ny) / N
    F_zz = np.sum(nz * nz) / N
    F_xy = np.sum(nx * ny) / N
    F_xz = np.sum(nx * nz) / N
    F_yz = np.sum(ny * nz) / N

    fabric_tensor = np.array([
        [F_xx, F_xy, F_xz],
        [F_xy, F_yy, F_yz],
        [F_xz, F_yz, F_zz]
    ])

    # Eigenvalue analysis
    eigenvalues, eigenvectors = np.linalg.eigh(fabric_tensor)
    eigenvalues = np.sort(eigenvalues)[::-1]  # Descending order

    results['fabric_tensor'] = fabric_tensor.tolist()
    results['fabric_eigenvalues'] = eigenvalues.tolist()

    # Anisotropy parameters
    # Deviator anisotropy
    a_d = eigenvalues[0] - eigenvalues[2]
    results['deviatoric_anisotropy'] = float(a_d)

    # Intermediate anisotropy
    if eigenvalues[0] != eigenvalues[2]:
        a_i = (eigenvalues[0] - 2*eigenvalues[1] + eigenvalues[2]) / (eigenvalues[0] - eigenvalues[2])
    else:
        a_i = 0
    results['intermediate_anisotropy'] = float(a_i)

    # Overall anisotropy (0 = isotropic, 1 = fully anisotropic)
    trace = np.trace(fabric_tensor)
    dev_tensor = fabric_tensor - (trace/3) * np.eye(3)
    anisotropy = np.sqrt(1.5 * np.sum(dev_tensor**2)) / trace if trace > 0 else 0
    results['anisotropy_index'] = float(anisotropy)

    # Principal direction (major eigenvector)
    major_direction = eigenvectors[:, 0].tolist()
    results['principal_direction'] = major_direction

    # Fabric by contact type
    for ct in df_contact['contact_type'].unique():
        mask = df_contact['contact_type'] == ct
        ct_nx = nx[mask]
        ct_ny = ny[mask]
        ct_nz = nz[mask]
        n_ct = np.sum(mask)

        if n_ct > 0:
            ct_Fzz = np.sum(ct_nz * ct_nz) / n_ct
            ct_Fxx = np.sum(ct_nx * ct_nx) / n_ct
            results[f'{ct}_fabric_zz'] = float(ct_Fzz)
            results[f'{ct}_fabric_xx'] = float(ct_Fxx)
            results[f'{ct}_fabric_ratio'] = float(ct_Fzz / ct_Fxx) if ct_Fxx > 0 else 0

    return results


# =============================================================================
# 4. TORQUE AND ROLLING ANALYSIS
# =============================================================================

def analyze_torque_rolling(df_contact, scale_factor=1000):
    """
    Analyze torque distribution and rolling behavior
    """
    results = {}

    # Calculate torque magnitude if not present
    if 'torque_mag' not in df_contact.columns:
        df_contact['torque_mag'] = np.sqrt(
            df_contact['torque_x']**2 +
            df_contact['torque_y']**2 +
            df_contact['torque_z']**2
        ) * scale_factor

    # Torque statistics
    results['torque_avg'] = float(df_contact['torque_mag'].mean())
    results['torque_max'] = float(df_contact['torque_mag'].max())
    results['torque_std'] = float(df_contact['torque_mag'].std())

    # Rolling vs sliding indicator
    # High F_t/F_n with low torque = sliding
    # Low F_t/F_n with high torque = rolling

    F_t = df_contact['F_tangential'] if 'F_tangential' in df_contact.columns else \
          np.sqrt(df_contact['ft_x']**2 + df_contact['ft_y']**2 + df_contact['ft_z']**2)
    F_n = df_contact['F_normal'] if 'F_normal' in df_contact.columns else \
          np.sqrt(df_contact['fn_x']**2 + df_contact['fn_y']**2 + df_contact['fn_z']**2)

    friction_ratio = F_t / (F_n + 1e-10)

    results['friction_ratio_avg'] = float(friction_ratio.mean())
    results['friction_ratio_max'] = float(friction_ratio.max())

    # Classify contacts
    friction_threshold = 0.3  # Below this = rolling-dominated
    rolling_contacts = np.sum(friction_ratio < friction_threshold)
    sliding_contacts = np.sum(friction_ratio >= friction_threshold)

    results['rolling_contacts'] = int(rolling_contacts)
    results['sliding_contacts'] = int(sliding_contacts)
    results['rolling_fraction'] = float(rolling_contacts / len(df_contact) * 100)

    # By contact type
    for ct in df_contact['contact_type'].unique():
        mask = df_contact['contact_type'] == ct
        results[f'{ct}_friction_ratio'] = float(friction_ratio[mask].mean())
        results[f'{ct}_torque_avg'] = float(df_contact.loc[mask, 'torque_mag'].mean())

    return results


# =============================================================================
# 5. POROSITY AND PACKING ANALYSIS
# =============================================================================

def analyze_porosity_packing(df_atom, n_layers=10, box_bounds=None):
    """
    Calculate layer-wise porosity and packing density
    """
    results = {}

    # Determine box bounds
    if box_bounds is None:
        x_min, x_max = df_atom['x'].min(), df_atom['x'].max()
        y_min, y_max = df_atom['y'].min(), df_atom['y'].max()
        z_min, z_max = df_atom['z'].min(), df_atom['z'].max()
    else:
        x_min, x_max = box_bounds['x']
        y_min, y_max = box_bounds['y']
        z_min, z_max = box_bounds['z']

    # Layer analysis along Z
    z_bins = np.linspace(z_min, z_max, n_layers + 1)
    layer_thickness = z_bins[1] - z_bins[0]
    layer_volume = (x_max - x_min) * (y_max - y_min) * layer_thickness

    layer_porosity = []
    layer_packing = []
    layer_z = []

    for i in range(n_layers):
        z_low, z_high = z_bins[i], z_bins[i+1]
        layer_atoms = df_atom[(df_atom['z'] >= z_low) & (df_atom['z'] < z_high)]

        if len(layer_atoms) > 0:
            # Solid volume (sum of sphere volumes)
            radii = layer_atoms['radius'].values
            solid_volume = np.sum(4/3 * np.pi * radii**3)

            packing_fraction = solid_volume / layer_volume
            porosity = 1 - packing_fraction
        else:
            packing_fraction = 0
            porosity = 1

        layer_packing.append(float(packing_fraction))
        layer_porosity.append(float(porosity))
        layer_z.append(float((z_low + z_high) / 2))

    results['layer_z'] = layer_z
    results['layer_porosity'] = layer_porosity
    results['layer_packing'] = layer_packing

    # Overall statistics
    valid_packing = [p for p in layer_packing if p > 0]
    if valid_packing:
        results['avg_packing_fraction'] = float(np.mean(valid_packing))
        results['avg_porosity'] = float(1 - np.mean(valid_packing))
        results['packing_uniformity'] = float(1 - np.std(valid_packing) / np.mean(valid_packing))

    # Packing by particle type
    for ptype in df_atom['type'].unique():
        type_atoms = df_atom[df_atom['type'] == ptype]
        type_volume = np.sum(4/3 * np.pi * type_atoms['radius'].values**3)
        results[f'type{ptype}_solid_volume'] = float(type_volume)

    return results


# =============================================================================
# 6. NETWORK ANALYSIS
# =============================================================================

def analyze_contact_network(df_contact, df_atom, type_names=None):
    """
    Graph-based network analysis: connectivity, shortest paths, centrality
    """
    if type_names is None:
        type_names = {1: 'NCM', 2: 'SE'}

    results = {}

    # Build adjacency matrix
    all_ids = set(df_contact['id1'].values) | set(df_contact['id2'].values)
    id_to_idx = {pid: idx for idx, pid in enumerate(sorted(all_ids))}
    n_particles = len(all_ids)

    # Sparse adjacency matrix
    rows = [id_to_idx[id1] for id1 in df_contact['id1']]
    cols = [id_to_idx[id2] for id2 in df_contact['id2']]
    data = np.ones(len(df_contact))

    adj_matrix = csr_matrix((data, (rows, cols)), shape=(n_particles, n_particles))
    adj_matrix = adj_matrix + adj_matrix.T  # Symmetric

    # Connected components
    n_components, labels = connected_components(adj_matrix, directed=False)
    results['n_connected_components'] = int(n_components)

    # Largest component size
    component_sizes = np.bincount(labels)
    results['largest_component_size'] = int(component_sizes.max())
    results['largest_component_fraction'] = float(component_sizes.max() / n_particles * 100)

    # SE-only network analysis
    id_to_type = dict(zip(df_atom['id'], df_atom['type']))
    se_type = [k for k, v in type_names.items() if v == 'SE'][0] if 'SE' in type_names.values() else 2

    se_contacts = df_contact[df_contact['contact_type'] == 'SE-SE']
    if len(se_contacts) > 0:
        se_ids = set(se_contacts['id1'].values) | set(se_contacts['id2'].values)
        se_id_to_idx = {pid: idx for idx, pid in enumerate(sorted(se_ids))}
        n_se = len(se_ids)

        se_rows = [se_id_to_idx[id1] for id1 in se_contacts['id1']]
        se_cols = [se_id_to_idx[id2] for id2 in se_contacts['id2']]
        se_data = np.ones(len(se_contacts))

        se_adj = csr_matrix((se_data, (se_rows, se_cols)), shape=(n_se, n_se))
        se_adj = se_adj + se_adj.T

        # SE connectivity
        se_n_comp, se_labels = connected_components(se_adj, directed=False)
        results['SE_n_components'] = int(se_n_comp)

        se_comp_sizes = np.bincount(se_labels)
        results['SE_largest_component'] = int(se_comp_sizes.max())
        results['SE_percolation_fraction'] = float(se_comp_sizes.max() / n_se * 100)

        # Average shortest path (sample for large networks)
        if n_se <= 500:
            dist_matrix = shortest_path(se_adj, directed=False, unweighted=True)
            finite_dists = dist_matrix[np.isfinite(dist_matrix) & (dist_matrix > 0)]
            if len(finite_dists) > 0:
                results['SE_avg_path_length'] = float(finite_dists.mean())
                results['SE_max_path_length'] = float(finite_dists.max())
        else:
            # Sample-based estimation
            sample_size = min(100, n_se)
            sample_indices = np.random.choice(n_se, sample_size, replace=False)
            dist_sample = shortest_path(se_adj, directed=False, unweighted=True, indices=sample_indices)
            finite_dists = dist_sample[np.isfinite(dist_sample) & (dist_sample > 0)]
            if len(finite_dists) > 0:
                results['SE_avg_path_length_estimated'] = float(finite_dists.mean())

    # Degree distribution (coordination number from graph)
    degrees = np.array(adj_matrix.sum(axis=1)).flatten()
    results['avg_degree'] = float(degrees.mean())
    results['max_degree'] = int(degrees.max())
    results['degree_std'] = float(degrees.std())

    # Clustering coefficient (local)
    # For each node, fraction of neighbors that are connected
    clustering_coeffs = []
    for i in range(min(500, n_particles)):  # Sample for large networks
        neighbors = adj_matrix[i].nonzero()[1]
        k = len(neighbors)
        if k >= 2:
            # Count edges among neighbors
            neighbor_adj = adj_matrix[neighbors][:, neighbors]
            n_edges = neighbor_adj.sum() / 2
            max_edges = k * (k - 1) / 2
            clustering_coeffs.append(n_edges / max_edges)

    if clustering_coeffs:
        results['avg_clustering_coefficient'] = float(np.mean(clustering_coeffs))

    return results


# =============================================================================
# 7. CONTACT POINT DISTRIBUTION ANALYSIS
# =============================================================================

def analyze_contact_points(df_contact, n_bins=20):
    """
    Analyze spatial distribution of contact points
    """
    results = {}

    cp_x = df_contact['cp_x'].values
    cp_y = df_contact['cp_y'].values
    cp_z = df_contact['cp_z'].values

    # Z-profile of contact density
    z_bins = np.linspace(cp_z.min(), cp_z.max(), n_bins + 1)
    z_hist, _ = np.histogram(cp_z, bins=z_bins)
    z_centers = (z_bins[:-1] + z_bins[1:]) / 2

    results['contact_z_profile'] = {
        'z': z_centers.tolist(),
        'count': z_hist.tolist(),
        'density': (z_hist / z_hist.sum()).tolist()
    }

    # Contact point clustering by type
    for ct in df_contact['contact_type'].unique():
        mask = df_contact['contact_type'] == ct
        ct_z = cp_z[mask]
        results[f'{ct}_z_mean'] = float(ct_z.mean())
        results[f'{ct}_z_std'] = float(ct_z.std())

    # Radial distribution (from center)
    center_x = (cp_x.max() + cp_x.min()) / 2
    center_y = (cp_y.max() + cp_y.min()) / 2

    radial_dist = np.sqrt((cp_x - center_x)**2 + (cp_y - center_y)**2)
    results['radial_mean'] = float(radial_dist.mean())
    results['radial_std'] = float(radial_dist.std())

    return results


# =============================================================================
# MAIN FUNCTION
# =============================================================================

def run_advanced_analysis(atoms_csv, contacts_csv, output_dir, scale_factor=1000, type_names=None):
    """
    Run all advanced analyses
    """
    if type_names is None:
        type_names = {1: 'NCM', 2: 'SE'}

    print("="*60)
    print("Advanced DEM Analysis")
    print("="*60)

    # Load data
    print("\n[1] Loading data...")
    df_atom = pd.read_csv(atoms_csv)
    df_contact = pd.read_csv(contacts_csv)
    print(f"    Atoms: {len(df_atom)}, Contacts: {len(df_contact)}")

    # Ensure force magnitudes exist
    if 'F_normal' not in df_contact.columns:
        df_contact['F_normal'] = np.sqrt(
            df_contact['fn_x']**2 + df_contact['fn_y']**2 + df_contact['fn_z']**2
        )
    if 'F_tangential' not in df_contact.columns:
        df_contact['F_tangential'] = np.sqrt(
            df_contact['ft_x']**2 + df_contact['ft_y']**2 + df_contact['ft_z']**2
        )

    all_results = {'scale_factor': scale_factor}

    # 1. Stress Tensor Analysis
    print("\n[2] Analyzing stress tensor...")
    stress_results = analyze_stress_tensor(df_atom, scale_factor)
    if stress_results:
        all_results['stress_tensor'] = stress_results
        print(f"    Von Mises avg: {stress_results['von_mises_avg']:.4f}")
        print(f"    Heterogeneity index: {stress_results['heterogeneity_index']:.4f}")

    # 2. Force Chain Analysis
    print("\n[3] Analyzing force chains...")
    force_chain_results = analyze_force_chains(df_contact, df_atom, percentile=90)
    all_results['force_chains'] = {k: v for k, v in force_chain_results.items()
                                   if k != 'force_chain_contacts'}
    print(f"    High-force contacts (top 10%): {force_chain_results['n_high_force_contacts']}")
    print(f"    Z-alignment: {force_chain_results.get('z_alignment', 'N/A'):.4f}")

    # Save force chain contacts separately
    chain_df = pd.DataFrame(force_chain_results.get('force_chain_contacts', []))
    if len(chain_df) > 0:
        chain_df.to_csv(os.path.join(output_dir, 'force_chain_contacts.csv'), index=False)

    # 3. Fabric Tensor Analysis
    print("\n[4] Analyzing fabric tensor...")
    fabric_results = analyze_fabric_tensor(df_contact, df_atom)
    all_results['fabric_tensor'] = fabric_results
    print(f"    Anisotropy index: {fabric_results['anisotropy_index']:.4f}")
    print(f"    Principal direction: {fabric_results['principal_direction']}")

    # 4. Torque and Rolling Analysis
    print("\n[5] Analyzing torque and rolling...")
    torque_results = analyze_torque_rolling(df_contact, scale_factor)
    all_results['torque_rolling'] = torque_results
    print(f"    Rolling fraction: {torque_results['rolling_fraction']:.1f}%")
    print(f"    Avg friction ratio: {torque_results['friction_ratio_avg']:.4f}")

    # 5. Porosity and Packing Analysis
    print("\n[6] Analyzing porosity and packing...")
    porosity_results = analyze_porosity_packing(df_atom, n_layers=15)
    all_results['porosity_packing'] = porosity_results
    print(f"    Avg packing fraction: {porosity_results.get('avg_packing_fraction', 'N/A'):.4f}")
    print(f"    Packing uniformity: {porosity_results.get('packing_uniformity', 'N/A'):.4f}")

    # 6. Network Analysis
    print("\n[7] Analyzing contact network...")
    network_results = analyze_contact_network(df_contact, df_atom, type_names)
    all_results['network'] = network_results
    print(f"    Connected components: {network_results['n_connected_components']}")
    print(f"    SE percolation: {network_results.get('SE_percolation_fraction', 'N/A'):.1f}%")

    # 7. Contact Point Distribution
    print("\n[8] Analyzing contact point distribution...")
    cp_results = analyze_contact_points(df_contact)
    all_results['contact_points'] = cp_results

    # Save results
    print("\n[9] Saving results...")

    # Save updated atom dataframe with stress analysis
    df_atom.to_csv(os.path.join(output_dir, 'atoms_advanced.csv'), index=False)

    # Save all results
    output_path = os.path.join(output_dir, 'advanced_analysis.json')
    with open(output_path, 'w') as f:
        json.dump(all_results, f, indent=2, default=float)

    print(f"    Saved to {output_dir}/")
    print("="*60)

    return all_results, df_atom, df_contact


def main():
    parser = argparse.ArgumentParser(description='Advanced DEM Analysis')
    parser.add_argument('atoms_csv', help='Path to analyzed atoms CSV')
    parser.add_argument('contacts_csv', help='Path to analyzed contacts CSV')
    parser.add_argument('--output-dir', '-o', default='./results',
                        help='Output directory')
    parser.add_argument('--scale-factor', '-s', type=float, default=1000,
                        help='Scale factor for stress/force')
    parser.add_argument('--type-names', '-t', type=str, default='1:NCM,2:SE',
                        help='Type name mapping')
    args = parser.parse_args()

    # Parse type names
    type_names = {}
    for pair in args.type_names.split(','):
        k, v = pair.split(':')
        type_names[int(k)] = v

    os.makedirs(args.output_dir, exist_ok=True)

    run_advanced_analysis(
        args.atoms_csv,
        args.contacts_csv,
        args.output_dir,
        args.scale_factor,
        type_names
    )


if __name__ == "__main__":
    main()
