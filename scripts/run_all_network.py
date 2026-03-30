"""
Batch runner: network conductivity solver for all cases.
Reads type_map from meta.json, runs network_conductivity.py on each.
"""
import os
import json
import sys
import subprocess
import time

WEBAPP_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'webapp')
RESULTS_DIR = os.path.join(WEBAPP_DIR, 'results')
ARCHIVE_DIR = os.path.join(WEBAPP_DIR, 'archive')
SCRIPT = os.path.join(os.path.dirname(__file__), 'network_conductivity.py')


def find_cases():
    """Find all cases with atoms.csv + contacts.csv + full_metrics.json."""
    cases = []

    # Dashboard cases
    if os.path.isdir(RESULTS_DIR):
        for case_id in os.listdir(RESULTS_DIR):
            case_dir = os.path.join(RESULTS_DIR, case_id)
            if not os.path.isdir(case_dir):
                continue
            atoms = os.path.join(case_dir, 'atoms.csv')
            contacts = os.path.join(case_dir, 'contacts.csv')
            metrics = os.path.join(case_dir, 'full_metrics.json')
            if os.path.exists(atoms) and os.path.exists(contacts) and os.path.exists(metrics):
                # Get type_map from upload meta.json
                upload_dir = os.path.join(WEBAPP_DIR, 'uploads', case_id)
                meta_file = os.path.join(upload_dir, 'meta.json')
                type_map = None
                if os.path.exists(meta_file):
                    with open(meta_file) as f:
                        meta = json.load(f)
                    type_map = meta.get('type_map', '')

                # Get case name
                name = case_id
                if os.path.exists(meta_file):
                    with open(meta_file) as f:
                        name = json.load(f).get('name', case_id)

                cases.append({
                    'id': case_id,
                    'name': name,
                    'dir': case_dir,
                    'atoms': atoms,
                    'contacts': contacts,
                    'type_map': type_map,
                })

    # Archive cases
    if os.path.isdir(ARCHIVE_DIR):
        for root, dirs, files in os.walk(ARCHIVE_DIR):
            if 'full_metrics.json' in files and 'atoms.csv' in files:
                contacts_csv = os.path.join(root, 'contacts.csv')
                if not os.path.exists(contacts_csv):
                    continue
                rel = os.path.relpath(root, ARCHIVE_DIR)
                meta_file = os.path.join(root, 'meta.json')
                type_map = None
                if os.path.exists(meta_file):
                    with open(meta_file) as f:
                        meta = json.load(f)
                    type_map = meta.get('type_map', '')
                name = os.path.basename(root)
                cases.append({
                    'id': f'archive:{rel}',
                    'name': name,
                    'dir': root,
                    'atoms': os.path.join(root, 'atoms.csv'),
                    'contacts': contacts_csv,
                    'type_map': type_map,
                })

    return cases


def infer_type_map(type_map_str):
    """Convert type_map string to CLI arg. e.g. '1:AM_P,2:AM_S,3:SE' """
    if not type_map_str:
        return '1:AM_S,2:SE'  # default
    return type_map_str


def main():
    cases = find_cases()
    print(f"Found {len(cases)} cases with atoms+contacts+metrics")
    print()

    results = []
    errors = []

    for i, case in enumerate(cases):
        type_map = infer_type_map(case['type_map'])
        print(f"[{i+1}/{len(cases)}] {case['name']} (type_map={type_map})")

        # Skip if already computed
        out_file = os.path.join(case['dir'], 'network_conductivity.json')
        if os.path.exists(out_file):
            with open(out_file) as f:
                existing = json.load(f)
            if existing.get('sigma_full') is not None:
                print(f"  → Already computed: σ_full={existing['sigma_full']:.6f}")
                existing['name'] = case['name']
                existing['case_id'] = case['id']
                results.append(existing)
                continue

        cmd = [
            sys.executable, SCRIPT,
            case['atoms'], case['contacts'],
            '-o', case['dir'],
            '-t', type_map,
            '-s', '1000'
        ]

        t0 = time.time()
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
            elapsed = time.time() - t0

            if proc.returncode == 0 and os.path.exists(out_file):
                with open(out_file) as f:
                    res = json.load(f)
                res['name'] = case['name']
                res['case_id'] = case['id']
                results.append(res)
                sigma = res.get('sigma_full_mScm', 'N/A')
                r_brug = res.get('R_brug_over_full', 'N/A')
                print(f"  → σ_full={sigma} mS/cm, R_brug={r_brug}× ({elapsed:.1f}s)")
            else:
                print(f"  → FAILED: {proc.stderr[-200:] if proc.stderr else 'unknown'}")
                errors.append(case['name'])
        except subprocess.TimeoutExpired:
            print(f"  → TIMEOUT (>300s)")
            errors.append(case['name'])

    # Summary
    print(f"\n{'='*60}")
    print(f"SUMMARY: {len(results)} succeeded, {len(errors)} failed")
    print(f"{'='*60}")

    if results:
        print(f"\n{'Name':30s} {'σ_full':>10s} {'σ_bulk_net':>10s} {'R_brug':>8s} {'bulk%':>6s}")
        print('-' * 70)
        for r in sorted(results, key=lambda x: x.get('sigma_full_mScm', 0) or 0):
            name = r.get('name', '?')[:28]
            sf = r.get('sigma_full_mScm')
            sb = r.get('sigma_bulk_net_mScm')
            rb = r.get('R_brug_over_full')
            bf = r.get('bulk_resistance_fraction')
            print(f"  {name:28s} {sf:10.5f} {sb:10.5f} {rb:8.2f} {bf:6.1%}" if sf else f"  {name:28s} {'N/A':>10s}")

    if errors:
        print(f"\nFailed cases: {errors}")

    # Save combined results
    combined_path = os.path.join(RESULTS_DIR, 'network_conductivity_all.json')
    with open(combined_path, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f"\nCombined results saved: {combined_path}")


if __name__ == '__main__':
    main()
