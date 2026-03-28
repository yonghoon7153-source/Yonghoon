#!/usr/bin/env python3
"""
Core DEM analysis functions.
Used by both analyze_contacts.py and analyze_contacts_bimodal.py.

Scale convention:
  - Length: sim(mm) → real(μm), factor = SCALE = 1000
  - E: sim = real × 0.001
  - P: sim = real × 0.001
  - Area: sim(m²) → real(μm²) = sim / SCALE² × 1e12
  - Force: F_real(μN) = F_sim(N) × SCALE  (∵ F_real_N = F_sim/SCALE)
  - Stress: σ_real(MPa) = σ_sim(Pa) × SCALE / 1e6
"""
import numpy as np
import json
import os
from collections import defaultdict
import networkx as nx


# ─── Porosity & Thickness ──────────────────────────────────────────────────

def calc_porosity(atoms, plate_z, box_xy=0.05):
    """
    Porosity = 1 - V_solid / V_box
    V_box uses mesh plate_z (= actual electrode thickness).
    All in sim units.
    """
    V_solid = sum(4/3 * np.pi * a['radius']**3 for a in atoms.values())
    V_box = box_xy * box_xy * plate_z
    porosity = (1 - V_solid / V_box) * 100
    return porosity


def get_plate_z(results_dir, atoms, scale):
    """Get plate_z from mesh_info.json, or estimate from atoms."""
    mesh_file = os.path.join(results_dir, 'mesh_info.json')
    if os.path.exists(mesh_file):
        with open(mesh_file) as f:
            info = json.load(f)
        return info['plate_z'], 'mesh'
    # Fallback: estimate from atom positions
    z_max = max(a['z'] + a['radius'] for a in atoms.values())
    return z_max, 'estimated'


# ─── Interface Area ────────────────────────────────────────────────────────

def calc_interface_area(atoms, contacts, type_map, scale):
    """
    Contact area by type. Returns dict with real μm² values.
    contactArea (sim m²) → real μm² = sim / scale² × 1e12
    """
    area_conv = 1.0 / (scale ** 2) * 1e12  # sim m² → real μm²
    ca = defaultdict(float)
    cc = defaultdict(int)

    for c in contacts:
        if c['id1'] in atoms and c['id2'] in atoms:
            t1 = type_map.get(atoms[c['id1']]['type'], '?')
            t2 = type_map.get(atoms[c['id2']]['type'], '?')
            ct = '-'.join(sorted([t1, t2]))
            ca[ct] += c['contact_area']
            cc[ct] += 1

    result = {}
    for ct in ca:
        result[ct] = {
            'total_area': ca[ct] * area_conv,
            'n_contacts': cc[ct],
            'mean_area': (ca[ct] / cc[ct]) * area_conv if cc[ct] > 0 else 0,
        }

    # AM_total-SE
    am_se_total = sum(v['total_area'] for k, v in result.items()
                      if 'SE' in k and k != 'SE-SE')
    am_se_n = sum(v['n_contacts'] for k, v in result.items()
                  if 'SE' in k and k != 'SE-SE')
    result['AM전체-SE'] = {
        'total_area': am_se_total,
        'n_contacts': am_se_n,
        'mean_area': am_se_total / am_se_n if am_se_n > 0 else 0,
    }

    return result


# ─── Coverage ──────────────────────────────────────────────────────────────

def calc_coverage(atoms, contacts, type_map, scale):
    """
    Coverage = SE접촉면적 / (전체표면적 - AM-AM접촉면적) per AM particle.
    Returns per-type mean, std, min, max.
    """
    am_types = [k for k, v in type_map.items() if 'AM' in v]
    se_types = [k for k, v in type_map.items() if v == 'SE']

    am_surf = {}
    am_se_area = defaultdict(float)
    am_am_area = defaultdict(float)

    for aid, a in atoms.items():
        if a['type'] in am_types:
            am_surf[aid] = 4 * np.pi * a['radius']**2

    for c in contacts:
        if c['id1'] not in atoms or c['id2'] not in atoms:
            continue
        t1, t2 = atoms[c['id1']]['type'], atoms[c['id2']]['type']

        if t1 in am_types and t2 in se_types:
            am_se_area[c['id1']] += c['contact_area']
        elif t2 in am_types and t1 in se_types:
            am_se_area[c['id2']] += c['contact_area']

        if t1 in am_types and t2 in am_types:
            am_am_area[c['id1']] += c['contact_area']
            am_am_area[c['id2']] += c['contact_area']

    result = {}
    for lbl in set(type_map.values()):
        if 'AM' not in lbl:
            continue
        tids = [k for k, v in type_map.items() if v == lbl]
        covs = []
        for aid, a in atoms.items():
            if a['type'] in tids:
                free = am_surf.get(aid, 0) - am_am_area.get(aid, 0)
                se = am_se_area.get(aid, 0)
                covs.append(min(se / free * 100, 100) if free > 0 else 0)
        if covs:
            covs = np.array(covs)
            result[lbl] = {
                'mean': float(np.mean(covs)),
                'std': float(np.std(covs)),
                'min': float(np.min(covs)),
                'max': float(np.max(covs)),
                'n': len(covs),
            }
    return result


# ─── SE-SE Coordination Number ─────────────────────────────────────────────

def calc_se_se_cn(atoms, contacts, se_types):
    cn = defaultdict(int)
    for c in contacts:
        if c['id1'] in atoms and c['id2'] in atoms:
            if atoms[c['id1']]['type'] in se_types and atoms[c['id2']]['type'] in se_types:
                cn[c['id1']] += 1
                cn[c['id2']] += 1

    se_ids = [aid for aid, a in atoms.items() if a['type'] in se_types]
    values = np.array([cn.get(aid, 0) for aid in se_ids])

    return {
        'mean': float(np.mean(values)),
        'std': float(np.std(values)),
        'ge2_pct': float(np.sum(values >= 2) / len(values) * 100) if len(values) > 0 else 0,
    }


# ─── SE Percolation ────────────────────────────────────────────────────────

def calc_percolation(atoms, contacts, se_types, plate_z, boundary_factor=2.0):
    """
    SE-SE contact graph → percolation analysis.
    Boundary: bottom = z <= r_SE × boundary_factor, top = z >= plate_z - r_SE × boundary_factor
    """
    se_ids = [aid for aid, a in atoms.items() if a['type'] in se_types]
    if not se_ids:
        return {'se_count': 0, 'percolation_pct': 0, 'top_reachable_pct': 0,
                'n_components': 0, 'largest_pct': 0,
                'top_reachable_se': set(), 'bottom_se': set(), 'top_se': set(), 'graph': nx.Graph()}

    r_se = atoms[se_ids[0]]['radius']
    z_bottom = 0.0 + r_se * boundary_factor
    z_top = plate_z - r_se * boundary_factor

    G = nx.Graph()
    G.add_nodes_from(se_ids)
    box_x = 0.05  # periodic boundary box size
    box_y = 0.05
    for c in contacts:
        if c['id1'] in atoms and c['id2'] in atoms:
            if atoms[c['id1']]['type'] in se_types and atoms[c['id2']]['type'] in se_types:
                # Add edge with periodic-corrected distance as weight
                d = _periodic_dist(atoms[c['id1']], atoms[c['id2']], box_x, box_y)
                G.add_edge(c['id1'], c['id2'], distance=d)

    bottom_se = {aid for aid in se_ids if atoms[aid]['z'] <= z_bottom}
    top_se = {aid for aid in se_ids if atoms[aid]['z'] >= z_top}

    components = list(nx.connected_components(G))
    largest = max(len(c) for c in components) if components else 0
    n_large = sum(1 for c in components if len(c) >= 10)

    percolating_se = set()
    top_reachable_se = set()

    for comp in components:
        has_bottom = len(comp & bottom_se) > 0
        has_top = len(comp & top_se) > 0
        if has_bottom and has_top:
            percolating_se.update(comp)
        if has_top:
            top_reachable_se.update(comp)

    n = len(se_ids)
    return {
        'se_count': n,
        'percolation_pct': len(percolating_se) / n * 100,
        'top_reachable_pct': len(top_reachable_se) / n * 100,
        'n_components': len(components),
        'n_large_components': n_large,
        'largest_pct': largest / n * 100,
        'top_reachable_se': top_reachable_se,
        'bottom_se': bottom_se,
        'top_se': top_se,
        'graph': G,
    }


# ─── Tortuosity ────────────────────────────────────────────────────────────

def _periodic_dist(a1, a2, box_x=0.05, box_y=0.05):
    """Distance between two atoms with periodic boundary in x,y (NOT z)."""
    dx = abs(a1['x'] - a2['x'])
    dy = abs(a1['y'] - a2['y'])
    dz = a1['z'] - a2['z']
    # Minimum image convention for periodic x,y
    dx = min(dx, box_x - dx)
    dy = min(dy, box_y - dy)
    return np.sqrt(dx**2 + dy**2 + dz**2)


def calc_tortuosity(atoms, perc_result, n_samples=200, box_xy=0.05):
    """tortuosity = path length / z distance.
    Uses distance-weighted shortest path (physical shortest distance)
    with periodic boundary correction.
    Includes all top-reachable SE (not just percolating) as sources."""
    G = perc_result['graph']
    top_reachable_se = perc_result['top_reachable_se']
    bottom_se = perc_result['bottom_se']
    top_se = perc_result['top_se']

    if not top_reachable_se:
        return {'mean': None, 'std': None, 'n_samples': 0}

    reach_top = list(top_se & top_reachable_se)
    if not reach_top:
        return {'mean': None, 'std': None, 'n_samples': 0}

    # Sources: lowest-z SE in each top-reachable component (not just bottom_se)
    components = list(nx.connected_components(G))
    src_candidates = []
    for comp in components:
        comp_top_reach = comp & top_reachable_se
        if not comp_top_reach:
            continue
        # Pick the lowest-z SE in this component as source
        lowest = min(comp_top_reach, key=lambda sid: atoms[sid]['z'])
        src_candidates.append(lowest)

    if not src_candidates:
        return {'mean': None, 'std': None, 'n_samples': 0}

    taus = []
    for i in range(min(n_samples, max(len(src_candidates), len(reach_top)))):
        src = src_candidates[i % len(src_candidates)]
        tgt = reach_top[i % len(reach_top)]
        try:
            path = nx.shortest_path(G, src, tgt, weight='distance')
            path_len = sum(
                _periodic_dist(atoms[path[k]], atoms[path[k+1]], box_xy, box_xy)
                for k in range(len(path) - 1))
            z_dist = abs(atoms[tgt]['z'] - atoms[src]['z'])
            if z_dist > 0:
                taus.append(path_len / z_dist)
        except nx.NetworkXNoPath:
            pass

    return {
        'mean': float(np.mean(taus)) if taus else None,
        'std': float(np.std(taus)) if taus else None,
        'n_samples': len(taus),
    }


# ─── Ionic Active AM ──────────────────────────────────────────────────────

def calc_ionic_active_am(atoms, contacts, perc_result, se_types, am_types, type_map):
    """
    Ionic Active AM = AM connected to SE that reaches top (SE pellet).
    """
    top_reachable_se = perc_result['top_reachable_se']
    am_ids = [aid for aid, a in atoms.items() if a['type'] in am_types]

    am_se_contacts = {}
    for c in contacts:
        if c['id1'] in atoms and c['id2'] in atoms:
            t1, t2 = atoms[c['id1']]['type'], atoms[c['id2']]['type']
            if t1 in am_types and t2 in se_types:
                am_se_contacts.setdefault(c['id1'], set()).add(c['id2'])
            elif t2 in am_types and t1 in se_types:
                am_se_contacts.setdefault(c['id2'], set()).add(c['id1'])

    ionic_active = set()
    ionic_dead = set()
    no_se = set()

    for aid in am_ids:
        if aid in am_se_contacts:
            if len(am_se_contacts[aid] & top_reachable_se) > 0:
                ionic_active.add(aid)
            else:
                ionic_dead.add(aid)
        else:
            no_se.add(aid)

    n_am = len(am_ids)
    result = {
        'n_am': n_am,
        'active_pct': len(ionic_active) / n_am * 100 if n_am > 0 else 0,
        'dead_pct': len(ionic_dead) / n_am * 100 if n_am > 0 else 0,
        'no_se_pct': len(no_se) / n_am * 100 if n_am > 0 else 0,
    }

    for lbl in set(type_map.values()):
        if 'AM' not in lbl:
            continue
        tids = [k for k, v in type_map.items() if v == lbl]
        sub = {aid for aid in am_ids if atoms[aid]['type'] in tids}
        sub_active = sub & ionic_active
        result[f'{lbl}_active_pct'] = len(sub_active) / len(sub) * 100 if sub else 0

    return result


# ─── Contact Force Distribution ───────────────────────────────────────────

def calc_contact_force_distribution(atoms, contacts, type_map, scale):
    """Normal force distribution by contact type. Returns stats in real units (μN)."""
    # Force conversion: sim N → real μN = sim / scale² × 1e6
    # In DEM: F_sim = E_sim * L_sim² proportional quantities
    # Real force = F_sim / scale² (length²) × 1e6 (N→μN)
    force_conv = 1.0 / (scale ** 2) * 1e6

    forces_by_type = defaultdict(list)
    for c in contacts:
        if c['id1'] in atoms and c['id2'] in atoms:
            t1 = type_map.get(atoms[c['id1']]['type'], '?')
            t2 = type_map.get(atoms[c['id2']]['type'], '?')
            ct = '-'.join(sorted([t1, t2]))
            fn = np.sqrt(c.get('fn_x', 0)**2 + c.get('fn_y', 0)**2 + c.get('fn_z', 0)**2)
            forces_by_type[ct].append(fn * force_conv)

    result = {}
    for ct, forces in forces_by_type.items():
        f = np.array(forces)
        result[ct] = {
            'mean': float(np.mean(f)),
            'std': float(np.std(f)),
            'max': float(np.max(f)),
            'n': len(f),
        }
    return result


# ─── Contact Pressure ────────────────────────────────────────────────────

def calc_contact_pressure(atoms, contacts, type_map, scale):
    """Contact pressure = Fn / contact_area per contact. Returns stats in MPa."""
    # Pressure = Force / Area. In sim units pressure is already Pa-like.
    # Real pressure = F_real(N) / A_real(m²)
    # F_real = F_sim / scale², A_real = A_sim / scale²  → P_real = F_sim / A_sim (same ratio!)
    # Convert to MPa: / 1e6
    pressures_by_type = defaultdict(list)
    for c in contacts:
        if c['id1'] in atoms and c['id2'] in atoms:
            ca = c.get('contact_area', 0)
            if ca <= 0:
                continue
            t1 = type_map.get(atoms[c['id1']]['type'], '?')
            t2 = type_map.get(atoms[c['id2']]['type'], '?')
            ct = '-'.join(sorted([t1, t2]))
            fn = np.sqrt(c.get('fn_x', 0)**2 + c.get('fn_y', 0)**2 + c.get('fn_z', 0)**2)
            pressure = fn / ca / 1e6  # Pa → MPa
            pressures_by_type[ct].append(pressure)

    result = {}
    for ct, pressures in pressures_by_type.items():
        p = np.array(pressures)
        result[ct] = {
            'mean': float(np.mean(p)),
            'std': float(np.std(p)),
            'max': float(np.max(p)),
        }
    # Overall
    all_p = []
    for v in pressures_by_type.values():
        all_p.extend(v)
    if all_p:
        ap = np.array(all_p)
        result['overall'] = {
            'mean': float(np.mean(ap)),
            'std': float(np.std(ap)),
            'max': float(np.max(ap)),
        }
    return result


# ─── Overlap Ratio ────────────────────────────────────────────────────────

def calc_overlap_ratio(atoms, contacts):
    """δ/R ratio for each contact. High values (>5%) indicate DEM inaccuracy."""
    ratios = []
    for c in contacts:
        delta = abs(c.get('delta', 0))
        if delta <= 0:
            continue
        if c['id1'] in atoms and c['id2'] in atoms:
            r1 = atoms[c['id1']]['radius']
            r2 = atoms[c['id2']]['radius']
            r_eff = min(r1, r2)
            if r_eff > 0:
                ratios.append(delta / r_eff * 100)  # percentage

    if not ratios:
        return None
    r = np.array(ratios)
    return {
        'mean': float(np.mean(r)),
        'std': float(np.std(r)),
        'max': float(np.max(r)),
        'pct_above_5': float(np.sum(r > 5) / len(r) * 100),
    }


# ─── AM Isolation Risk ───────────────────────────────────────────────────

def calc_am_isolation_risk(atoms, contacts, type_map):
    """AM particles connected to only 1 SE = vulnerable (single point of failure)."""
    am_types = [k for k, v in type_map.items() if 'AM' in v]
    se_types = [k for k, v in type_map.items() if v == 'SE']

    am_se_count = defaultdict(int)
    for c in contacts:
        if c['id1'] in atoms and c['id2'] in atoms:
            t1 = atoms[c['id1']]['type']
            t2 = atoms[c['id2']]['type']
            if t1 in am_types and t2 in se_types:
                am_se_count[c['id1']] += 1
            elif t2 in am_types and t1 in se_types:
                am_se_count[c['id2']] += 1

    am_ids = [aid for aid, a in atoms.items() if a['type'] in am_types]
    n_am = len(am_ids)
    if n_am == 0:
        return None

    no_se = sum(1 for aid in am_ids if am_se_count.get(aid, 0) == 0)
    single_se = sum(1 for aid in am_ids if am_se_count.get(aid, 0) == 1)
    counts = [am_se_count.get(aid, 0) for aid in am_ids]

    return {
        'no_se_pct': float(no_se / n_am * 100),
        'single_se_pct': float(single_se / n_am * 100),
        'vulnerable_pct': float((no_se + single_se) / n_am * 100),
        'am_se_cn_mean': float(np.mean(counts)),
        'am_se_cn_std': float(np.std(counts)),
    }


# ─── Effective Ionic Conductivity ─────────────────────────────────────────

def calc_effective_conductivity(atoms, perc_result, porosity, tortuosity_result, type_map, plate_z, box_xy=0.05):
    """Estimate σ_eff / σ_bulk using Bruggeman-like relation.
    σ_eff/σ_bulk = ε_SE^α / τ²  where α ≈ 1 (connected fraction)."""
    se_types = [k for k, v in type_map.items() if v == 'SE']
    se_ids = [aid for aid, a in atoms.items() if a['type'] in se_types]

    if not se_ids:
        return None

    # SE volume fraction
    v_electrode = box_xy * box_xy * plate_z
    v_se = sum(4/3 * np.pi * atoms[aid]['radius']**3 for aid in se_ids)
    phi_se = v_se / v_electrode if v_electrode > 0 else 0

    tau = tortuosity_result.get('mean')
    if not tau or tau <= 0:
        return None

    # Connected SE fraction (percolating)
    perc_pct = perc_result.get('percolation_pct', 0) / 100

    # Effective conductivity ratio
    sigma_ratio = phi_se * perc_pct / (tau ** 2)

    return {
        'phi_se': float(phi_se),
        'tau': float(tau),
        'perc_fraction': float(perc_pct),
        'sigma_ratio': float(sigma_ratio),
    }


# ─── Von Mises Stress Analysis ─────────────────────────────────────────────

def calc_von_mises_stress(atoms_raw, type_map, scale, plate_z, n_layers=10):
    """
    Simplified Von Mises from diagonal stress components.
    σ_VM = √(σ_xx² + σ_yy² + σ_zz² - σ_xx·σ_yy - σ_yy·σ_zz - σ_xx·σ_zz)

    atoms_raw must have 'sigma_xx', 'sigma_yy', 'sigma_zz' keys (sim Pa).
    Returns: overall CV, type ratios, z-layer CV profile.
    Note: absolute values are not reliable (effective E used), only relative comparison.
    """
    # Check if stress data available
    sample = next(iter(atoms_raw.values()))
    if 'sigma_xx' not in sample:
        return None

    # Compute Von Mises for each atom (sim units, relative only)
    vm_data = {}
    for aid, a in atoms_raw.items():
        sxx = a.get('sigma_xx', 0)
        syy = a.get('sigma_yy', 0)
        szz = a.get('sigma_zz', 0)
        vm = np.sqrt(sxx**2 + syy**2 + szz**2 - sxx*syy - syy*szz - sxx*szz)
        vm_data[aid] = vm

    all_vm = np.array(list(vm_data.values()))
    vm_mean = float(np.mean(all_vm))
    vm_std = float(np.std(all_vm))
    vm_cv = (vm_std / vm_mean * 100) if vm_mean > 0 else 0

    # Type-specific mean (for ratio calculation)
    type_stress = {}
    for t_name in set(type_map.values()):
        t_ids = [aid for aid, a in atoms_raw.items() if type_map.get(a['type']) == t_name]
        if t_ids:
            t_vm = np.array([vm_data[aid] for aid in t_ids])
            type_stress[t_name] = {
                'mean': float(np.mean(t_vm)),
                'ratio': float(np.mean(t_vm) / vm_mean) if vm_mean > 0 else 0,
            }

    # Z-layer CV profile
    z_layer_cv = []
    z_values = np.array([a['z'] for a in atoms_raw.values()])
    z_min, z_max = 0.0, plate_z
    layer_edges = np.linspace(z_min, z_max, n_layers + 1)

    atom_ids = list(atoms_raw.keys())
    atom_z = np.array([atoms_raw[aid]['z'] for aid in atom_ids])
    atom_vm = np.array([vm_data[aid] for aid in atom_ids])

    for i in range(n_layers):
        mask = (atom_z >= layer_edges[i]) & (atom_z < layer_edges[i+1])
        if mask.sum() > 1:
            layer_vm = atom_vm[mask]
            layer_mean = np.mean(layer_vm)
            layer_cv = (np.std(layer_vm) / layer_mean * 100) if layer_mean > 0 else 0
            z_mid = (layer_edges[i] + layer_edges[i+1]) / 2 * scale  # to μm
            z_layer_cv.append({
                'z_mid_um': float(z_mid),
                'cv': float(layer_cv),
                'mean_normalized': float(layer_mean / vm_mean) if vm_mean > 0 else 0,
            })

    return {
        'vm_cv': vm_cv,
        'vm_mean': vm_mean,
        'type_stress': type_stress,
        'z_layer_cv': z_layer_cv,
    }


# ─── Full Analysis ─────────────────────────────────────────────────────────

def run_full_analysis(atoms_raw, contacts_raw, type_map, scale, results_dir, box_xy=0.05):
    """
    Run all analyses. atoms_raw/contacts_raw are dicts in SIM units.
    Returns comprehensive results dict.
    """
    se_types = [k for k, v in type_map.items() if v == 'SE']
    am_types = [k for k, v in type_map.items() if 'AM' in v]

    # Get plate_z
    plate_z, pz_source = get_plate_z(results_dir, atoms_raw, scale)
    thickness_um = plate_z * scale  # sim → μm

    print(f"  plate_z = {plate_z:.6f} ({pz_source}), thickness = {thickness_um:.1f} μm")

    # 1. Porosity
    porosity = calc_porosity(atoms_raw, plate_z, box_xy)
    print(f"  Porosity: {porosity:.2f}%")

    # 2. Interface Area
    iface = calc_interface_area(atoms_raw, contacts_raw, type_map, scale)
    for ct, v in iface.items():
        print(f"  {ct}: {v['n_contacts']} contacts, total {v['total_area']:.1f} μm²")

    # 3. Coverage
    cov = calc_coverage(atoms_raw, contacts_raw, type_map, scale)
    for lbl, v in cov.items():
        print(f"  Coverage {lbl}: {v['mean']:.1f}% ± {v['std']:.1f}%")

    # 4. SE-SE CN
    cn = calc_se_se_cn(atoms_raw, contacts_raw, se_types)
    print(f"  SE-SE CN: {cn['mean']:.2f} ± {cn['std']:.2f}")

    # 5. Percolation
    perc = calc_percolation(atoms_raw, contacts_raw, se_types, plate_z)
    print(f"  Percolation: {perc['percolation_pct']:.1f}%, Top Reachable: {perc['top_reachable_pct']:.1f}%")
    print(f"  Components: {perc['n_components']}, Largest: {perc['largest_pct']:.1f}%")

    # 6. Tortuosity
    tau = calc_tortuosity(atoms_raw, perc)
    tau_str = f"{tau['mean']:.2f} ± {tau['std']:.2f}" if tau['mean'] else "N/A"
    print(f"  Tortuosity: {tau_str} ({tau['n_samples']} samples)")

    # 7. Ionic Active AM
    ionic = calc_ionic_active_am(atoms_raw, contacts_raw, perc, se_types, am_types, type_map)
    print(f"  Ionic Active AM: {ionic['active_pct']:.1f}%")

    # 8. Von Mises Stress (relative)
    stress = calc_von_mises_stress(atoms_raw, type_map, scale, plate_z)
    if stress:
        print(f"  Stress CV: {stress['vm_cv']:.1f}%")
        for tn, sv in stress['type_stress'].items():
            print(f"    {tn}: σ/σ_mean = {sv['ratio']:.2f}")
    else:
        print("  Stress: N/A (no stress data)")

    # 9. Contact Force Distribution
    force_dist = calc_contact_force_distribution(atoms_raw, contacts_raw, type_map, scale)
    for ct, v in force_dist.items():
        print(f"  Force {ct}: {v['mean']:.3f} ± {v['std']:.3f} μN (max {v['max']:.3f})")

    # 10. Contact Pressure
    contact_pressure = calc_contact_pressure(atoms_raw, contacts_raw, type_map, scale)
    if 'overall' in contact_pressure:
        cp = contact_pressure['overall']
        print(f"  Contact Pressure: {cp['mean']:.1f} ± {cp['std']:.1f} MPa (max {cp['max']:.1f})")

    # 11. Overlap Ratio
    overlap = calc_overlap_ratio(atoms_raw, contacts_raw)
    if overlap:
        print(f"  Overlap δ/R: {overlap['mean']:.2f}% ± {overlap['std']:.2f}% (max {overlap['max']:.2f}%, >5%: {overlap['pct_above_5']:.1f}%)")

    # 12. AM Isolation Risk
    am_risk = calc_am_isolation_risk(atoms_raw, contacts_raw, type_map)
    if am_risk:
        print(f"  AM Risk: vulnerable {am_risk['vulnerable_pct']:.1f}% (no SE: {am_risk['no_se_pct']:.1f}%, single SE: {am_risk['single_se_pct']:.1f}%)")

    # 13. Effective Ionic Conductivity
    eff_cond = calc_effective_conductivity(atoms_raw, perc, porosity, tau, type_map, plate_z, box_xy)
    if eff_cond:
        print(f"  σ_eff/σ_bulk: {eff_cond['sigma_ratio']:.4f} (φ_SE={eff_cond['phi_se']:.3f}, τ={eff_cond['tau']:.2f})")

    return {
        'plate_z': plate_z,
        'plate_z_source': pz_source,
        'thickness_um': thickness_um,
        'porosity': porosity,
        'interface': iface,
        'coverage': cov,
        'se_se_cn': cn,
        'percolation': perc,
        'tortuosity': tau,
        'ionic_active': ionic,
        'stress': stress,
        'force_dist': force_dist,
        'contact_pressure': contact_pressure,
        'overlap_ratio': overlap,
        'am_isolation_risk': am_risk,
        'effective_conductivity': eff_cond,
    }
