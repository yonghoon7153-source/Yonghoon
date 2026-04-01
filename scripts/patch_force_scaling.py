"""
Patch force/pressure values in existing full_metrics.json files.

Old code: force_conv = 1/scale² × 1e6 = 1.0 (for scale=1000) → F_sim 그대로
Correct:  force_conv = 1e6/scale = 1000 → ×1000 차이

Old code: pressure = F_sim / A_sim / 1e6 → P_sim / 1e6
Correct:  pressure = F_sim / A_sim × scale / 1e6 → ×scale 차이

Fix: force values ×1000, pressure values ×1000
"""
import json, os

WEBAPP = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'webapp')
SCALE_FIX = 1000  # multiply by this to correct


def patch_all():
    count = 0
    for base in [os.path.join(WEBAPP, 'results'), os.path.join(WEBAPP, 'archive')]:
        if not os.path.isdir(base):
            continue
        for root, dirs, files in os.walk(base):
            if 'full_metrics.json' not in files:
                continue
            met_path = os.path.join(root, 'full_metrics.json')
            with open(met_path) as f:
                met = json.load(f)

            changed = False

            # Force keys: fn_*_mean, fn_*_max, fn_*_std
            for key in list(met.keys()):
                if key.startswith('fn_') and ('_mean' in key or '_max' in key or '_std' in key):
                    if met[key] and isinstance(met[key], (int, float)) and met[key] != 0:
                        met[key] = met[key] * SCALE_FIX
                        changed = True

            # Pressure keys: contact_pressure_mean, contact_pressure_max
            for key in ['contact_pressure_mean', 'contact_pressure_max']:
                if key in met and met[key] and isinstance(met[key], (int, float)) and met[key] != 0:
                    met[key] = met[key] * SCALE_FIX
                    changed = True

            if changed:
                with open(met_path, 'w') as f:
                    json.dump(met, f, indent=2, default=str)
                count += 1
                name = os.path.basename(root)
                fn_mean = met.get('fn_SE_SE_mean', met.get('fn_AM_SE_mean', '?'))
                p_mean = met.get('contact_pressure_mean', '?')
                print(f"  {name}: fn={fn_mean}, P={p_mean} MPa")

    print(f"\nPatched {count} cases (force ×{SCALE_FIX}, pressure ×{SCALE_FIX})")


if __name__ == '__main__':
    patch_all()
