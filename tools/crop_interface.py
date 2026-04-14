#!/usr/bin/env python3
"""
Crop adhesion interface xyz files for VESTA visualization.
Extracts atoms within z-range around NCM-SE interface.

Usage:
    python tools/crop_interface.py comp1_v5xy_s52.xyz
    python tools/crop_interface.py comp1_v5xy_s52.xyz --zpad 8
    python tools/crop_interface.py *.xyz --zpad 10
"""

import sys
import os
import argparse


def read_xyz(filename):
    """Read extended xyz file, return atoms list and header."""
    with open(filename) as f:
        lines = f.readlines()

    n_atoms = int(lines[0].strip())
    header = lines[1].strip()
    atoms = []
    for line in lines[2:2 + n_atoms]:
        parts = line.split()
        if len(parts) >= 4:
            species = parts[0]
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
            rest = parts[4:] if len(parts) > 4 else []
            atoms.append((species, x, y, z, rest))

    return n_atoms, header, atoms


def find_interface(atoms, n_ncm):
    """Find NCM-SE interface z-position."""
    ncm_atoms = atoms[:n_ncm]
    se_atoms = atoms[n_ncm:]

    ncm_z = [a[3] for a in ncm_atoms]
    se_z = [a[3] for a in se_atoms]

    ncm_top = max(ncm_z)
    se_bottom = min(se_z)
    interface_z = (ncm_top + se_bottom) / 2

    return ncm_top, se_bottom, interface_z


def detect_n_ncm(atoms):
    """Detect NCM/SE boundary by finding where P/S/Cl start after Ni/O region."""
    ncm_elements = {'Li', 'Ni', 'O'}
    se_unique = {'P', 'S', 'Cl', 'Br'}

    # Find first P atom (SE starts here)
    for i, (sp, x, y, z, rest) in enumerate(atoms):
        if sp == 'P':
            return i

    # Fallback: count Ni atoms * 4 (each NCM unit has Li+Ni+2O = 4 atoms)
    n_ni = sum(1 for sp, *_ in atoms if sp == 'Ni')
    return n_ni * 4


def crop_and_save(filename, z_pad=8.0, n_ncm=None):
    """Crop atoms around interface and save new xyz."""
    n_atoms, header, atoms = read_xyz(filename)

    if n_ncm is None:
        n_ncm = detect_n_ncm(atoms)

    ncm_top, se_bottom, interface_z = find_interface(atoms, n_ncm)

    z_min = interface_z - z_pad
    z_max = interface_z + z_pad

    # Crop
    cropped = []
    for sp, x, y, z, rest in atoms:
        if z_min <= z <= z_max:
            cropped.append((sp, x, y, z, rest))

    # Count by type
    ncm_count = sum(1 for i, (sp, x, y, z, rest) in enumerate(atoms)
                     if i < n_ncm and z_min <= z <= z_max)
    se_count = len(cropped) - ncm_count

    # Output filename
    base = os.path.splitext(filename)[0]
    outfile = f"{base}_crop.xyz"

    with open(outfile, 'w') as f:
        f.write(f"{len(cropped)}\n")
        # Modify header to note cropping
        f.write(f'{header} comment="cropped z={z_min:.1f}~{z_max:.1f}"\n')
        for sp, x, y, z, rest in cropped:
            rest_str = '  '.join(rest) if rest else ''
            if rest_str:
                f.write(f"{sp}  {x:.8f}  {y:.8f}  {z:.8f}  {rest_str}\n")
            else:
                f.write(f"{sp}  {x:.8f}  {y:.8f}  {z:.8f}\n")

    print(f"=== {filename} ===")
    print(f"  Total: {n_atoms} atoms (NCM: {n_ncm}, SE: {n_atoms - n_ncm})")
    print(f"  NCM top: {ncm_top:.1f} A, SE bottom: {se_bottom:.1f} A")
    print(f"  Interface: {interface_z:.1f} A, Gap: {se_bottom - ncm_top:.1f} A")
    print(f"  Crop: z = {z_min:.1f} ~ {z_max:.1f} A ({2*z_pad:.0f} A window)")
    print(f"  Cropped: {len(cropped)} atoms (NCM: {ncm_count}, SE: {se_count})")
    print(f"  Saved: {outfile}")
    print()

    return outfile


def main():
    parser = argparse.ArgumentParser(description="Crop adhesion interface for VESTA")
    parser.add_argument("files", nargs="+", help="xyz files to crop")
    parser.add_argument("--zpad", type=float, default=8.0,
                        help="z padding around interface (A, default 8)")
    parser.add_argument("--ncm", type=int, default=None,
                        help="Number of NCM atoms (auto-detect if not set)")
    args = parser.parse_args()

    for f in args.files:
        if os.path.exists(f):
            crop_and_save(f, z_pad=args.zpad, n_ncm=args.ncm)
        else:
            print(f"File not found: {f}")


if __name__ == "__main__":
    main()
