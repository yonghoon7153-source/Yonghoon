"""Merge network_conductivity.json into full_metrics.json for each case."""
import json, os

WEBAPP = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'webapp')

def merge_all():
    count = 0
    for base in [os.path.join(WEBAPP, 'results'), os.path.join(WEBAPP, 'archive')]:
        if not os.path.isdir(base):
            continue
        for root, dirs, files in os.walk(base):
            if 'network_conductivity.json' in files and 'full_metrics.json' in files:
                net_path = os.path.join(root, 'network_conductivity.json')
                met_path = os.path.join(root, 'full_metrics.json')
                with open(net_path) as f:
                    net = json.load(f)
                with open(met_path) as f:
                    met = json.load(f)
                
                changed = False
                for key in ['sigma_full', 'sigma_full_mScm', 'sigma_bulk_net',
                           'sigma_bulk_net_mScm', 'R_brug_over_full',
                           'bulk_resistance_fraction',
                           'electronic_sigma_full_mScm', 'electronic_R_brug',
                           'thermal_sigma_full_mScm', 'thermal_R_brug']:
                    if key in net and net[key] is not None:
                        met[key] = net[key]
                        changed = True
                
                if changed:
                    with open(met_path, 'w') as f:
                        json.dump(met, f, indent=2, default=str)
                    count += 1
    
    print(f"Merged {count} cases")

if __name__ == '__main__':
    merge_all()
