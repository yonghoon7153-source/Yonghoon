"""
DEM-Native Ionic Transport Framework
=====================================
Resistor network solver for effective ionic conductivity in ASSB composite cathodes.

Three decomposition runs:
  1. FULL: R_bulk + R_constriction per edge → σ_full (ground truth)
  2. BULK_ONLY: R_bulk per edge (R_constriction=0) → σ_bulk_net (≈Bruggeman?)
  3. CONSTRICTION_ONLY: R_constriction per edge (R_bulk=0) → σ_constr_net (≈current proxy?)

Comparison with Bruggeman and BLM+Constriction scaling validates/refines existing models.

References:
  - Nature Comm 2025: Resistor network for SSB composites (Münster)
  - Holm 1967: Maxwell constriction resistance
  - De Levie: Transmission line model for porous electrodes
"""

import numpy as np
import json
import os
import sys
from scipy import sparse
from scipy.sparse.linalg import spsolve


# LPSCl argyrodite grain interior conductivity (NOT pellet value)
SIGMA_BULK_DEFAULT = 3.0e-3  # S/cm (grain interior, ionic)

# NCM electronic conductivity (typical, SOC-dependent)
SIGMA_AM_ELECTRONIC = 0.05  # S/cm (50 mS/cm, discharged NCM)

# Thermal conductivity (W/(m·K) = W/(m·K) × 10⁻⁴ = W/(cm·K))
K_AM_THERMAL = 4.0e-2   # W/(cm·K) ≈ 4 W/(m·K), NCM
K_SE_THERMAL = 0.7e-2   # W/(cm·K) ≈ 0.7 W/(m·K), LPSCl (Ketter 2025)


def build_network(atoms_raw, contacts_raw, target_types, scale,
                  plate_z, box_x=0.05, box_y=0.05, boundary_factor=2.0,
                  mode='ionic', type_map=None, results_dir=None):
    """
    Build resistor network from DEM data.
    mode='ionic': SE-SE network only (uses percolation_sets.json for boundaries)
    mode='electronic': AM-AM network only
    mode='thermal': ALL contacts (AM-AM, AM-SE, SE-SE)
    Returns nodes, edges, bottom/top boundary sets.
    """
    if mode == 'thermal':
        target_ids = list(atoms_raw.keys())
    else:
        target_ids = [aid for aid, a in atoms_raw.items() if a['type'] in target_types]

    if not target_ids:
        return None

    # Boundary detection: use percolation_sets.json if available (ionic)
    bottom_ids = None
    top_ids = None
    if mode == 'ionic' and results_dir:
        perc_path = os.path.join(results_dir, 'percolation_sets.json')
        if os.path.exists(perc_path):
            with open(perc_path) as f:
                perc_sets = json.load(f)
            bottom_ids = {int(x) for x in perc_sets.get('bottom_se', [])} & set(target_ids)
            top_ids = {int(x) for x in perc_sets.get('top_reachable', [])} & set(target_ids)
            if bottom_ids and top_ids:
                print(f"  Boundaries from percolation_sets.json: {len(bottom_ids)} bottom, {len(top_ids)} top_reachable")

    # Fallback: z-coordinate based
    if not bottom_ids or not top_ids:
        r_ref = atoms_raw[target_ids[0]]['radius']
        z_bottom = 0.0 + r_ref * boundary_factor
        z_top = plate_z - r_ref * boundary_factor
        bottom_ids = {aid for aid in target_ids if atoms_raw[aid]['z'] <= z_bottom}
        top_ids = {aid for aid in target_ids if atoms_raw[aid]['z'] >= z_top}

    # Determine SE types for thermal mode
    se_type_set = set()
    if type_map:
        se_type_set = {k for k, v in type_map.items() if v == 'SE'}

    # Build contact area map
    contact_area_map = {}
    for c in contacts_raw:
        id1, id2 = c['id1'], c['id2']
        if id1 in atoms_raw and id2 in atoms_raw:
            if mode == 'thermal':
                # All contacts for thermal
                pair = (min(id1, id2), max(id1, id2))
                ca = c.get('contact_area', 0)
                if ca > 0:
                    contact_area_map[pair] = ca
            else:
                if atoms_raw[id1]['type'] in target_types and atoms_raw[id2]['type'] in target_types:
                    pair = (min(id1, id2), max(id1, id2))
                    ca = c.get('contact_area', 0)
                    if ca > 0:
                        contact_area_map[pair] = ca

    # Build edges with physical resistances
    # All distances in μm, areas in μm², resistivity in Ω·μm
    # ρ = 1/σ, σ = 1.3e-3 S/cm = 1.3e-7 S/μm → ρ = 7.69e6 Ω·μm
    # But we normalize: set ρ=1, then σ_eff comes out as ratio to σ_bulk

    edges = []
    for pair, ca_sim in contact_area_map.items():
        id1, id2 = pair
        a1, a2 = atoms_raw[id1], atoms_raw[id2]

        # Hop distance (μm) with periodic boundary
        dx = abs(a1['x'] - a2['x'])
        dy = abs(a1['y'] - a2['y'])
        dz = a1['z'] - a2['z']
        dx = min(dx, box_x - dx)
        dy = min(dy, box_y - dy)
        d_ij = np.sqrt(dx**2 + dy**2 + dz**2) * scale  # μm

        # Contact area (μm²) and contact radius (μm)
        A_contact = ca_sim * scale**2  # μm²
        a_contact = np.sqrt(A_contact / np.pi)  # μm

        # Particle radii (μm)
        r1 = a1['radius'] * scale
        r2 = a2['radius'] * scale

        # Thermal mode: material-specific conductivity weighting
        # k_AM ≈ 4.0 W/m·K, k_SE ≈ 0.7 W/m·K
        # AM-AM: weight=k_AM/k_SE, SE-SE: weight=1, AM-SE: harmonic mean
        if mode == 'thermal' and se_type_set:
            t1_is_se = a1['type'] in se_type_set
            t2_is_se = a2['type'] in se_type_set
            k_ratio = K_AM_THERMAL / K_SE_THERMAL  # ~5.7
            if not t1_is_se and not t2_is_se:
                # AM-AM: high thermal conductivity
                k_weight = k_ratio
            elif t1_is_se and t2_is_se:
                # SE-SE: baseline
                k_weight = 1.0
            else:
                # AM-SE: harmonic mean
                k_weight = 2 * k_ratio / (1 + k_ratio)
        else:
            k_weight = 1.0

        # Normalized resistances (ρ=1, scaled by k_weight for thermal):
        # R_bulk = d / (k_weight × π × r²)
        R_bulk_1 = (d_ij / 2) / (k_weight * np.pi * r1**2) if r1 > 0 else 0
        R_bulk_2 = (d_ij / 2) / (k_weight * np.pi * r2**2) if r2 > 0 else 0
        R_bulk = R_bulk_1 + R_bulk_2

        # R_constriction = 1 / (k_weight × 2a)  (Maxwell spreading)
        R_constriction = 1.0 / (k_weight * 2 * a_contact) if a_contact > 0 else 1e12

        edges.append({
            'id1': id1, 'id2': id2,
            'R_bulk': R_bulk,
            'R_constriction': R_constriction,
            'R_total': R_bulk + R_constriction,
            'd_ij': d_ij,
            'A_contact': A_contact,
        })

    return {
        'nodes': target_ids,
        'edges': edges,
        'bottom': bottom_ids,
        'top': top_ids,
        'plate_z': plate_z,
        'box_x': box_x,
        'box_y': box_y,
        'scale': scale,
    }


def solve_network(network_data, mode='full'):
    """
    Solve resistor network for effective conductance.

    mode: 'full' (R_bulk + R_constriction),
          'bulk_only' (R_bulk, R_constriction=0),
          'constriction_only' (R_constriction, R_bulk=0)

    Returns:
        G_eff: effective conductance (normalized, ρ=1)
        sigma_ratio: σ_eff / σ_bulk
    """
    nodes = network_data['nodes']
    edges = network_data['edges']
    bottom = network_data['bottom']
    top = network_data['top']
    scale = network_data['scale']
    plate_z = network_data['plate_z']
    box_x = network_data['box_x']
    box_y = network_data['box_y']

    if not bottom or not top or not edges:
        return None, None

    # Build networkx graph to find percolating component
    import networkx as nx
    G_nx = nx.Graph()
    for e in edges:
        G_nx.add_edge(e['id1'], e['id2'])

    # Find components that connect bottom to top
    perc_nodes = set()
    for comp in nx.connected_components(G_nx):
        has_bot = len(comp & bottom) > 0
        has_top = len(comp & top) > 0
        if has_bot and has_top:
            perc_nodes |= comp

    if not perc_nodes:
        print("  No percolating component found")
        return None, None

    # Filter to percolating nodes only
    perc_bottom = bottom & perc_nodes
    perc_top = top & perc_nodes
    perc_edges = [e for e in edges if e['id1'] in perc_nodes and e['id2'] in perc_nodes]

    print(f"  Percolating component: {len(perc_nodes)} nodes, {len(perc_edges)} edges")

    # Node index mapping (percolating only)
    all_ids = list(perc_nodes)
    id_to_idx = {nid: i for i, nid in enumerate(all_ids)}
    N = len(all_ids)

    # Virtual source (idx=N) and sink (idx=N+1)
    source_idx = N
    sink_idx = N + 1
    total_nodes = N + 2

    # Build conductance matrix (sparse)
    row, col, val = [], [], []

    def add_conductance(i, j, g):
        if g <= 0:
            return
        # Add g to (i,i), (j,j) and subtract from (i,j), (j,i)
        row.extend([i, j, i, j])
        col.extend([i, j, j, i])
        val.extend([g, g, -g, -g])

    for e in perc_edges:
        i = id_to_idx[e['id1']]
        j = id_to_idx[e['id2']]

        if mode == 'full':
            R = e['R_total']
        elif mode == 'bulk_only':
            R = e['R_bulk'] if e['R_bulk'] > 0 else 1e-12
        elif mode == 'constriction_only':
            R = e['R_constriction']
        else:
            R = e['R_total']

        if R > 0:
            g = 1.0 / R
            add_conductance(i, j, g)

    # Connect bottom SE to source with large conductance (low resistance)
    g_boundary = 1e6  # effectively zero resistance
    for bid in perc_bottom:
        add_conductance(id_to_idx[bid], source_idx, g_boundary)

    # Connect top SE to sink
    for tid in perc_top:
        add_conductance(id_to_idx[tid], sink_idx, g_boundary)

    # Build sparse Laplacian
    L = sparse.csr_matrix((val, (row, col)), shape=(total_nodes, total_nodes))

    # Right-hand side: inject current at source, extract at sink
    b = np.zeros(total_nodes)
    b[source_idx] = 1.0
    b[sink_idx] = -1.0

    # Ground one node to make system solvable (pin sink to V=0)
    # Modify: remove sink equation, set V_sink = 0
    # Instead, use pseudo-approach: fix V_sink = 0 by replacing its row
    L_mod = L.tolil()
    L_mod[sink_idx, :] = 0
    L_mod[sink_idx, sink_idx] = 1.0
    b[sink_idx] = 0.0
    L_mod = L_mod.tocsr()

    try:
        V = spsolve(L_mod, b)
    except Exception as e:
        print(f"  Network solve failed: {e}")
        return None, None

    V_source = V[source_idx]
    V_sink = V[sink_idx]  # = 0

    if V_source <= 0:
        return None, None

    # G_eff = I / ΔV = 1.0 / V_source  (since I=1, V_sink=0)
    G_eff = 1.0 / V_source

    # Convert to σ_eff / σ_bulk
    # G_eff is in normalized units (ρ=1)
    # σ_eff = G_eff × L / A where L = plate_z*scale (μm), A = box_x*box_y*scale² (μm²)
    T_um = plate_z * scale
    A_um2 = box_x * box_y * scale**2

    # σ_ratio = σ_eff / σ_bulk = G_eff × T / A  (dimensionless when ρ=1)
    sigma_ratio = G_eff * T_um / A_um2

    return G_eff, sigma_ratio


def run_decomposition(atoms_raw, contacts_raw, target_types, scale,
                      plate_z, box_x=0.05, box_y=0.05,
                      sigma_bulk=SIGMA_BULK_DEFAULT, results_dir=None,
                      type_map=None):
    """
    Run full decomposition analysis:
    1. FULL (bulk + constriction)
    2. BULK_ONLY (Bruggeman-like)
    3. CONSTRICTION_ONLY (current proxy-like)

    Also computes Bruggeman prediction for comparison.
    """
    print(f"  Building resistor network ({len(target_types)} target types)...")
    net = build_network(atoms_raw, contacts_raw, target_types, scale,
                        plate_z, box_x, box_y, results_dir=results_dir,
                        type_map=type_map)

    if net is None:
        print("  No network found")
        return None

    n_nodes = len(net['nodes'])
    n_edges = len(net['edges'])
    n_bottom = len(net['bottom'])
    n_top = len(net['top'])
    print(f"  Network: {n_nodes} nodes, {n_edges} edges, {n_bottom} bottom, {n_top} top")

    # Edge statistics
    R_bulks = [e['R_bulk'] for e in net['edges']]
    R_constrs = [e['R_constriction'] for e in net['edges']]
    R_totals = [e['R_total'] for e in net['edges']]

    bulk_frac = np.mean([rb/(rb+rc) for rb, rc in zip(R_bulks, R_constrs) if rb+rc > 0])
    print(f"  R_bulk fraction: {bulk_frac:.1%} (vs R_constriction: {1-bulk_frac:.1%})")

    # === Run 1: FULL ===
    print("  Solving FULL network (bulk + constriction)...")
    G_full, sigma_full = solve_network(net, mode='full')

    # === Run 2: BULK ONLY ===
    print("  Solving BULK_ONLY network...")
    G_bulk, sigma_bulk_net = solve_network(net, mode='bulk_only')

    # === Run 3: CONSTRICTION ONLY ===
    print("  Solving CONSTRICTION_ONLY network...")
    G_constr, sigma_constr_net = solve_network(net, mode='constriction_only')

    # === Bruggeman prediction ===
    V_se = sum(4/3 * np.pi * atoms_raw[aid]['radius']**3
               for aid in net['nodes'])
    V_box = box_x * box_y * plate_z
    phi_se = V_se / V_box if V_box > 0 else 0

    # Results
    results = {
        'n_nodes': n_nodes,
        'n_edges': n_edges,
        'n_bottom': n_bottom,
        'n_top': n_top,
        'phi_se': round(phi_se, 4),
        'bulk_resistance_fraction': round(bulk_frac, 4),
        'sigma_full': round(sigma_full, 8) if sigma_full else None,
        'sigma_bulk_net': round(sigma_bulk_net, 8) if sigma_bulk_net else None,
        'sigma_constr_net': round(sigma_constr_net, 8) if sigma_constr_net else None,
        'sigma_full_mScm': round(sigma_full * sigma_bulk * 1000, 6) if sigma_full else None,
        'sigma_bulk_net_mScm': round(sigma_bulk_net * sigma_bulk * 1000, 6) if sigma_bulk_net else None,
        'sigma_constr_net_mScm': round(sigma_constr_net * sigma_bulk * 1000, 6) if sigma_constr_net else None,
    }

    # Bruggeman comparison
    if sigma_bulk_net and sigma_full:
        results['R_brug_over_full'] = round(sigma_bulk_net / sigma_full, 4)
        results['R_brug_over_constr'] = round(sigma_bulk_net / sigma_constr_net, 4) if sigma_constr_net else None

    # Print summary
    print(f"\n  ═══ Decomposition Results ═══")
    print(f"  φ_SE = {phi_se:.4f}")
    print(f"  R_bulk fraction per edge: {bulk_frac:.1%}")
    print(f"")
    print(f"  {'Mode':<20s} {'σ/σ_bulk':>12s} {'σ (mS/cm)':>12s}")
    print(f"  {'─'*44}")
    if sigma_full:
        print(f"  {'FULL (truth)':20s} {sigma_full:12.6f} {sigma_full*sigma_bulk*1000:12.6f}")
    if sigma_bulk_net:
        print(f"  {'BULK_ONLY':20s} {sigma_bulk_net:12.6f} {sigma_bulk_net*sigma_bulk*1000:12.6f}")
    if sigma_constr_net:
        print(f"  {'CONSTRICTION_ONLY':20s} {sigma_constr_net:12.6f} {sigma_constr_net*sigma_bulk*1000:12.6f}")
    print(f"")
    if sigma_bulk_net and sigma_full:
        print(f"  Bruggeman overestimation (σ_bulk_net/σ_full): {sigma_bulk_net/sigma_full:.2f}×")
    if sigma_constr_net and sigma_full:
        print(f"  Constriction underestimation (σ_constr/σ_full): {sigma_constr_net/sigma_full:.4f}×")

    return results


if __name__ == '__main__':
    import argparse
    sys.path.insert(0, os.path.dirname(__file__))
    from analyze_contacts import load_atoms_raw, load_contacts_raw

    parser = argparse.ArgumentParser(description='DEM-Native Ionic Transport Solver')
    parser.add_argument('atoms_csv', help='atoms.csv path')
    parser.add_argument('contacts_csv', help='contacts.csv path')
    parser.add_argument('-o', '--output', required=True, help='Output directory')
    parser.add_argument('-t', '--type-map', default='1:AM_S,2:SE', help='Type map')
    parser.add_argument('-s', '--scale', type=int, default=1000, help='Scale factor')
    args = parser.parse_args()

    # Parse type map
    type_map = {}
    for pair in args.type_map.split(','):
        k, v = pair.split(':')
        type_map[int(k)] = v.strip()

    target_types = [k for k, v in type_map.items() if v == 'SE']
    am_types = [k for k, v in type_map.items() if 'AM' in v]

    # Load data
    atoms_raw, _ = load_atoms_raw(args.atoms_csv)
    contacts_raw, _ = load_contacts_raw(args.contacts_csv)
    print(f"Loaded {len(atoms_raw)} atoms, {len(contacts_raw)} contacts")

    # Get plate_z
    mesh_file = os.path.join(args.output, 'mesh_info.json')
    if os.path.exists(mesh_file):
        with open(mesh_file) as f:
            plate_z = json.load(f)['plate_z']
    else:
        plate_z = max(a['z'] + a['radius'] for a in atoms_raw.values())

    # Get box dimensions
    box_x, box_y = 0.05, 0.05
    ip_path = os.path.join(args.output, 'input_params.json')
    if os.path.exists(ip_path):
        with open(ip_path) as f:
            ip = json.load(f)
        box_x = ip.get('box_x', 0.05)
        box_y = ip.get('box_y', 0.05)

    print(f"box={box_x}×{box_y}, plate_z={plate_z:.6f}, scale={args.scale}")

    # === IONIC (SE-SE network) ===
    print("\n" + "="*50)
    print("IONIC CONDUCTIVITY (SE-SE network)")
    print("="*50)
    results = run_decomposition(atoms_raw, contacts_raw, target_types, args.scale,
                                plate_z, box_x, box_y, sigma_bulk=SIGMA_BULK_DEFAULT,
                                results_dir=args.output, type_map=type_map)

    # === ELECTRONIC (AM-AM network) ===
    results_el = None
    if am_types:
        print("\n" + "="*50)
        print("ELECTRONIC CONDUCTIVITY (AM-AM network)")
        print("="*50)
        try:
            results_el = run_decomposition(atoms_raw, contacts_raw, am_types, args.scale,
                                           plate_z, box_x, box_y, sigma_bulk=SIGMA_AM_ELECTRONIC,
                                           type_map=type_map)
        except Exception as e:
            print(f"  Electronic solver failed: {e}")

    # === THERMAL (ALL contacts) ===
    results_th = None
    try:
        print("\n" + "="*50)
        print("THERMAL CONDUCTIVITY (ALL contacts)")
        print("="*50)
        all_types = list(type_map.keys())
        results_th = run_decomposition(atoms_raw, contacts_raw, all_types, args.scale,
                                       plate_z, box_x, box_y, sigma_bulk=K_SE_THERMAL,
                                       type_map=type_map)
    except Exception as e:
        print(f"  Thermal solver failed: {e}")

    # Save results (ionic is primary)
    if results:
        if results_el:
            results['electronic_sigma_full'] = results_el.get('sigma_full')
            results['electronic_sigma_full_mScm'] = results_el.get('sigma_full_mScm')
            results['electronic_R_brug'] = results_el.get('R_brug_over_full')
            results['electronic_bulk_frac'] = results_el.get('bulk_resistance_fraction')
            results['electronic_n_nodes'] = results_el.get('n_nodes')
            results['electronic_n_edges'] = results_el.get('n_edges')
        if results_th:
            results['thermal_sigma_full'] = results_th.get('sigma_full')
            results['thermal_sigma_full_mScm'] = results_th.get('sigma_full_mScm')
            results['thermal_R_brug'] = results_th.get('R_brug_over_full')
            results['thermal_bulk_frac'] = results_th.get('bulk_resistance_fraction')

        out_path = os.path.join(args.output, 'network_conductivity.json')
        with open(out_path, 'w') as f:
            json.dump(results, f, indent=2)
        print(f"\nResults saved: {out_path}")
