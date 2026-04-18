#!/usr/bin/env python3
"""
Backfill e_se_eff_gpa into existing full_metrics.json files.
Reads youngs_modulus_sim from input_params.json, computes E_SE effective GPa.
"""
import os, json, glob

SCALE = 1000  # DEM scale factor

def backfill():
    bases = [
        os.path.join(os.path.dirname(__file__), '..', 'webapp', 'results'),
        os.path.join(os.path.dirname(__file__), '..', 'webapp', 'archive'),
    ]
    updated = 0
    for base in bases:
        for root, dirs, files in os.walk(base):
            if 'full_metrics.json' not in files or 'input_params.json' not in files:
                continue
            ip_path = os.path.join(root, 'input_params.json')
            met_path = os.path.join(root, 'full_metrics.json')
            with open(ip_path) as f:
                ip = json.load(f)
            with open(met_path) as f:
                met = json.load(f)

            ym = ip.get('youngs_modulus_sim', [])
            if len(ym) < 2:
                continue
            e_se_eff_gpa = round(ym[-1] * SCALE / 1e9, 4)

            if met.get('e_se_eff_gpa') == e_se_eff_gpa:
                continue  # already correct

            met['e_se_eff_gpa'] = e_se_eff_gpa
            with open(met_path, 'w') as f:
                json.dump(met, f, indent=2, default=str)
            print(f"  Updated {root}: E_SE_eff = {e_se_eff_gpa} GPa")
            updated += 1

    print(f"\nDone. Updated {updated} cases.")

if __name__ == '__main__':
    backfill()
