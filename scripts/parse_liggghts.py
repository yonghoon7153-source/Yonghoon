#!/usr/bin/env python3
"""
Parse LIGGGHTS atom, contact, and mesh dump files into CSV.

Usage:
    python parse_liggghts.py atom_*.liggghts contact_*.liggghts -o ./results
    python parse_liggghts.py atom_*.liggghts contact_*.liggghts mesh_*.stl -o ./results
"""
import argparse
import os
import sys
import re
import json
import numpy as np


def parse_atom_file(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()
    headers, rows = None, []
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith('ITEM: ATOMS'):
            headers = line.replace('ITEM: ATOMS', '').strip().split()
            i += 1
            while i < len(lines) and not lines[i].strip().startswith('ITEM:'):
                vals = lines[i].strip().split()
                if len(vals) == len(headers):
                    rows.append(vals)
                i += 1
            continue
        i += 1
    return headers, rows


def parse_contact_file(filepath):
    col_map = {
        'c_cpl[1]': 'p1_x', 'c_cpl[2]': 'p1_y', 'c_cpl[3]': 'p1_z',
        'c_cpl[4]': 'p2_x', 'c_cpl[5]': 'p2_y', 'c_cpl[6]': 'p2_z',
        'c_cpl[7]': 'id1', 'c_cpl[8]': 'id2',
        'c_cpl[9]': 'periodic_flag',
        'c_cpl[10]': 'fx', 'c_cpl[11]': 'fy', 'c_cpl[12]': 'fz',
        'c_cpl[13]': 'fn_x', 'c_cpl[14]': 'fn_y', 'c_cpl[15]': 'fn_z',
        'c_cpl[16]': 'ft_x', 'c_cpl[17]': 'ft_y', 'c_cpl[18]': 'ft_z',
        'c_cpl[19]': 'torque_x', 'c_cpl[20]': 'torque_y', 'c_cpl[21]': 'torque_z',
        'c_cpl[22]': 'contact_area', 'c_cpl[23]': 'delta',
        'c_cpl[24]': 'cp_x', 'c_cpl[25]': 'cp_y', 'c_cpl[26]': 'cp_z',
    }
    with open(filepath, 'r') as f:
        lines = f.readlines()
    headers, rows = None, []
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith('ITEM: ENTRIES'):
            raw = line.replace('ITEM: ENTRIES', '').strip().split()
            headers = [col_map.get(h, h) for h in raw]
            i += 1
            while i < len(lines) and not lines[i].strip().startswith('ITEM:'):
                vals = lines[i].strip().split()
                if len(vals) == len(headers):
                    rows.append(vals)
                i += 1
            continue
        i += 1
    return headers, rows


def parse_mesh_stl(filepath):
    """Parse mesh STL -> plate z coordinate."""
    with open(filepath, 'r') as f:
        lines = f.readlines()
    zs = []
    for line in lines:
        if 'vertex' in line.lower():
            parts = line.split()
            if len(parts) >= 4:
                zs.append(float(parts[3]))
    return np.mean(zs) if zs else None


def find_last_file(file_list):
    if len(file_list) == 1:
        return file_list[0]
    def extract_num(fname):
        nums = re.findall(r'(\d+)', os.path.basename(fname))
        return int(nums[-1]) if nums else 0
    return sorted(file_list, key=extract_num)[-1]


def main():
    parser = argparse.ArgumentParser(description='Parse LIGGGHTS dump files to CSV')
    parser.add_argument('files', nargs='+', help='atom_*.liggghts, contact_*.liggghts, mesh_*.stl')
    parser.add_argument('-o', '--output', default='./results')
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    atom_files = [f for f in args.files if 'atom' in os.path.basename(f).lower() and f.endswith('.liggghts')]
    contact_files = [f for f in args.files if 'contact' in os.path.basename(f).lower() and f.endswith('.liggghts')]
    mesh_files = [f for f in args.files if f.lower().endswith('.stl')]

    if not atom_files:
        print("ERROR: No atom files found", file=sys.stderr); sys.exit(1)
    if not contact_files:
        print("ERROR: No contact files found", file=sys.stderr); sys.exit(1)

    atom_file = find_last_file(atom_files)
    contact_file = find_last_file(contact_files)

    print(f"Parsing atom file: {atom_file}")
    headers_a, rows_a = parse_atom_file(atom_file)
    if not headers_a or not rows_a:
        print("ERROR: Failed to parse atom file", file=sys.stderr); sys.exit(1)
    print(f"  -> {len(rows_a)} atoms")

    with open(os.path.join(args.output, 'atoms.csv'), 'w') as f:
        f.write(','.join(headers_a) + '\n')
        for row in rows_a:
            f.write(','.join(row) + '\n')

    print(f"Parsing contact file: {contact_file}")
    headers_c, rows_c = parse_contact_file(contact_file)
    if not headers_c or not rows_c:
        print("ERROR: Failed to parse contact file", file=sys.stderr); sys.exit(1)
    print(f"  -> {len(rows_c)} contacts")

    with open(os.path.join(args.output, 'contacts.csv'), 'w') as f:
        f.write(','.join(headers_c) + '\n')
        for row in rows_c:
            f.write(','.join(row) + '\n')

    # Mesh STL (optional)
    if mesh_files:
        mesh_file = find_last_file(mesh_files)
        print(f"Parsing mesh file: {mesh_file}")
        plate_z = parse_mesh_stl(mesh_file)
        if plate_z is not None:
            with open(os.path.join(args.output, 'mesh_info.json'), 'w') as f:
                json.dump({'plate_z': plate_z, 'source': os.path.basename(mesh_file)}, f, indent=2)
            print(f"  -> plate_z = {plate_z:.6f}")
    else:
        print("No mesh STL files (will estimate thickness from atoms)")

    print("Parse complete!")


if __name__ == '__main__':
    main()
