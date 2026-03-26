#!/usr/bin/env python3
"""
Parse LIGGGHTS atom and contact dump files into CSV.
Handles both standard (2-type) and bimodal (3-type) systems.

Usage:
    python parse_liggghts.py atom_*.liggghts contact_*.liggghts -o ./results
"""
import argparse
import os
import sys

def parse_atom_file(filepath):
    """Parse a LIGGGHTS atom dump file. Returns (header_list, rows_list)."""
    rows = []
    headers = None
    with open(filepath, 'r') as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line == 'ITEM: TIMESTEP':
            i += 2  # skip timestep value
        elif line.startswith('ITEM: NUMBER OF ATOMS'):
            i += 1
            n_atoms = int(lines[i].strip())
            i += 1
        elif line.startswith('ITEM: BOX BOUNDS'):
            i += 3  # skip 3 bound lines
        elif line.startswith('ITEM: ATOMS'):
            # Parse header
            parts = line.replace('ITEM: ATOMS', '').strip().split()
            headers = parts
            i += 1
            # Read atom data
            while i < len(lines) and not lines[i].strip().startswith('ITEM:'):
                vals = lines[i].strip().split()
                if len(vals) == len(headers):
                    rows.append(vals)
                i += 1
            continue
        i += 1

    return headers, rows


def parse_contact_file(filepath):
    """Parse a LIGGGHTS contact/local dump file. Returns (header_list, rows_list)."""
    rows = []
    headers = None
    with open(filepath, 'r') as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        line = lines[i].strip()
        if line == 'ITEM: TIMESTEP':
            i += 2
        elif line.startswith('ITEM: NUMBER OF ENTRIES'):
            i += 1
            i += 1
        elif line.startswith('ITEM: BOX BOUNDS'):
            i += 3
        elif line.startswith('ITEM: ENTRIES'):
            raw_headers = line.replace('ITEM: ENTRIES', '').strip().split()
            # Map cpl columns to meaningful names
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
            headers = [col_map.get(h, h) for h in raw_headers]
            i += 1
            while i < len(lines) and not lines[i].strip().startswith('ITEM:'):
                vals = lines[i].strip().split()
                if len(vals) == len(headers):
                    rows.append(vals)
                i += 1
            continue
        i += 1

    return headers, rows


def find_last_file(file_list):
    """From a list of dump files, return the one with the latest timestep."""
    if len(file_list) == 1:
        return file_list[0]
    # Sort by filename (typically contains timestep number)
    import re
    def extract_num(fname):
        nums = re.findall(r'(\d+)', os.path.basename(fname))
        return int(nums[-1]) if nums else 0
    return sorted(file_list, key=extract_num)[-1]


def main():
    parser = argparse.ArgumentParser(description='Parse LIGGGHTS dump files to CSV')
    parser.add_argument('files', nargs='+', help='atom_*.liggghts and contact_*.liggghts files')
    parser.add_argument('-o', '--output', default='./results', help='Output directory')
    parser.add_argument('--all-timesteps', action='store_true',
                        help='Parse all timesteps (default: last only)')
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    # Separate atom and contact files
    atom_files = [f for f in args.files if 'atom' in os.path.basename(f).lower()]
    contact_files = [f for f in args.files if 'contact' in os.path.basename(f).lower()]

    if not atom_files:
        print("ERROR: No atom files found", file=sys.stderr)
        sys.exit(1)
    if not contact_files:
        print("ERROR: No contact files found", file=sys.stderr)
        sys.exit(1)

    # Use last timestep file
    atom_file = find_last_file(atom_files)
    contact_file = find_last_file(contact_files)

    print(f"Parsing atom file: {atom_file}")
    headers_a, rows_a = parse_atom_file(atom_file)
    if not headers_a or not rows_a:
        print("ERROR: Failed to parse atom file", file=sys.stderr)
        sys.exit(1)
    print(f"  -> {len(rows_a)} atoms parsed")

    print(f"Parsing contact file: {contact_file}")
    headers_c, rows_c = parse_contact_file(contact_file)
    if not headers_c or not rows_c:
        print("ERROR: Failed to parse contact file", file=sys.stderr)
        sys.exit(1)
    print(f"  -> {len(rows_c)} contacts parsed")

    # Write CSVs
    atoms_csv = os.path.join(args.output, 'atoms.csv')
    with open(atoms_csv, 'w') as f:
        f.write(','.join(headers_a) + '\n')
        for row in rows_a:
            f.write(','.join(row) + '\n')
    print(f"  -> Saved: {atoms_csv}")

    contacts_csv = os.path.join(args.output, 'contacts.csv')
    with open(contacts_csv, 'w') as f:
        f.write(','.join(headers_c) + '\n')
        for row in rows_c:
            f.write(','.join(row) + '\n')
    print(f"  -> Saved: {contacts_csv}")

    print("Parse complete!")


if __name__ == '__main__':
    main()
