#!/usr/bin/env python3
"""
Electronic THIN — advanced formula sweep with compound expressions.
Targets: τ(ρ=0.929), ratio(ρ=-0.756), por(ρ=0.720), gc(ρ=0.544), δ(ρ=0.454)
"""
import numpy as np, json, sys, os
sys.path.insert(0, os.path.dirname(__file__))
from screening_electronic_sweep import load_all_electronic

SIGMA_AM = 50.0
THRESHOLD = 8

def fit_r2(log_s, log_rhs):
    if len(log_s) < 4: return -999, 0
    lnC = np.mean(log_s - log_rhs)
    ss_res = np.sum((log_s - lnC - log_rhs)**2)
    ss_tot = np.sum((log_s - np.mean(log_s))**2)
    return (1 - ss_res/ss_tot if ss_tot > 1e-12 else -999), np.exp(lnC)

data = load_all_electronic()
thin = [r for r in data if r['ratio'] < THRESHOLD]
print(f"Thin data: n={len(thin)} (T/d < {THRESHOLD})")

s = np.log(np.array([r['sigma'] for r in thin]))
# Variables
cn = np.array([r['cn'] for r in thin])
tau = np.array([r['tau'] for r in thin])
delta = np.array([r['delta'] for r in thin])
ratio = np.array([r['ratio'] for r in thin])
por = np.array([r['por'] for r in thin])
phi_am = np.array([r['phi_am'] for r in thin])
phi_se = np.array([r['phi_se'] for r in thin])
cov = np.array([r['cov'] for r in thin])
area = np.array([r['area'] for r in thin])
contact_r = np.array([r['contact_r'] for r in thin])
gc = np.array([r['gc'] for r in thin])
bottleneck = np.array([r['bottleneck'] for r in thin])
hop = np.array([r['hop'] for r in thin])
force = np.array([r['force'] for r in thin])
ep = np.array([r['ep'] for r in thin])

# Compound variables
d2_over_a = delta**2 / area
cn_x_delta = cn * delta
tau_x_delta = tau * delta
cn_over_ratio = cn / ratio
gc_x_tau = gc * tau
por_x_cov = por/100 * cov
phi_am_minus = np.clip(phi_am - 0.15, 0.001, None)  # percolation excess

log_vars = {
    'cn': np.log(cn), 'tau': np.log(tau), 'delta': np.log(delta),
    'ratio': np.log(ratio), 'por': np.log(por), 'phi_am': np.log(phi_am),
    'phi_se': np.log(phi_se), 'cov': np.log(cov), 'area': np.log(area),
    'contact_r': np.log(contact_r), 'gc': np.log(gc),
    'bottleneck': np.log(bottleneck), 'hop': np.log(hop),
    'force': np.log(force), 'ep': np.log(ep),
    # Compound
    'd2/A': np.log(d2_over_a), 'cn×δ': np.log(cn_x_delta),
    'τ×δ': np.log(tau_x_delta), 'cn/ξ': np.log(cn_over_ratio),
    'gc×τ': np.log(gc_x_tau), 'por×cov': np.log(por_x_cov),
    '(φ-φc)': np.log(phi_am_minus),
}

print("\n" + "="*60)
print("Phase 1: IONIC-style formulas (compound expressions)")
print("="*60)

results = []
# Type A: (φ-φc)^a × X^b × Y^c
for a in [0.5, 0.75, 1.0, 1.5, 2.0]:
    lp = a * log_vars['(φ-φc)']
    for vn, vl in log_vars.items():
        if vn == '(φ-φc)': continue
        for b in [-1.0, -0.5, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0]:
            r2, C = fit_r2(s, lp + b * vl)
            results.append((r2, f"(φ-φc)^{a} × {vn}^{b}", C))
            # + 3rd variable
            for vn2, vl2 in log_vars.items():
                if vn2 in ('(φ-φc)', vn): continue
                for c in [-0.5, 0.25, 0.5, 1.0]:
                    r2_3, C3 = fit_r2(s, lp + b*vl + c*vl2)
                    results.append((r2_3, f"(φ-φc)^{a} × {vn}^{b} × {vn2}^{c}", C3))

results.sort(reverse=True)
print("Top 30:")
for r2, formula, C in results[:30]:
    print(f"  R²={r2:.4f}  C={C:.4f}  {formula}")

print("\n" + "="*60)
print("Phase 2: τ-based formulas (Spearman champion)")
print("="*60)

results2 = []
for a in [0.5, 0.75, 1.0, 1.5, 2.0]:
    lt = a * log_vars['tau']
    for vn, vl in log_vars.items():
        if vn == 'tau': continue
        for b in [-1.0, -0.5, -0.25, 0.25, 0.5, 0.75, 1.0]:
            for vn2, vl2 in log_vars.items():
                if vn2 in ('tau', vn): continue
                for c in [-0.5, -0.25, 0.25, 0.5, 1.0]:
                    r2, C = fit_r2(s, lt + b*vl + c*vl2)
                    results2.append((r2, f"τ^{a} × {vn}^{b} × {vn2}^{c}", C))

results2.sort(reverse=True)
print("Top 30:")
for r2, formula, C in results2[:30]:
    print(f"  R²={r2:.4f}  C={C:.4f}  {formula}")

print("\n" + "="*60)
print("Phase 3: FORM-X style (σ_grain × (φ-φc)^a × CN^b × f(contact) / f(geometry))")
print("="*60)

results3 = []
SGRAIN = 3.0
for a in [0.5, 0.75, 1.0, 1.5]:
    for b in [0.5, 0.75, 1.0, 1.5, 2.0]:
        base = a * log_vars['(φ-φc)'] + b * log_vars['cn']
        # Contact quality terms
        for contact_name, contact_log in [('√δ', 0.5*log_vars['delta']),
                                            ('√A', 0.5*log_vars['area']),
                                            ('√gc', 0.5*log_vars['gc']),
                                            ('√a_c', 0.5*log_vars['contact_r']),
                                            ('δ', log_vars['delta']),
                                            ('gc^0.25', 0.25*log_vars['gc'])]:
            # Geometry correction terms
            for geo_name, geo_log in [('1/√ξ', -0.5*log_vars['ratio']),
                                       ('1/ξ', -1.0*log_vars['ratio']),
                                       ('1/√por', -0.5*log_vars['por']),
                                       ('√τ', 0.5*log_vars['tau']),
                                       ('τ', log_vars['tau']),
                                       ('1/√τ', -0.5*log_vars['tau'])]:
                rhs = base + contact_log + geo_log
                r2, C = fit_r2(s, rhs)
                results3.append((r2, f"(φ-φc)^{a}×CN^{b}×{contact_name}×{geo_name}", C))

results3.sort(reverse=True)
print("Top 30:")
for r2, formula, C in results3[:30]:
    print(f"  R²={r2:.4f}  C={C:.4f}  {formula}")

print("\n" + "="*60)
print("Phase 4: Best candidates — detailed comparison")
print("="*60)

# Test the top formulas and show per-case predictions
candidates = [
    ("CN×√δ/√ξ (current)", lambda: np.log(SIGMA_AM) + log_vars['cn'] + 0.5*log_vars['delta'] - 0.5*log_vars['ratio']),
    ("τ×√δ×√gc", lambda: np.log(SIGMA_AM) + log_vars['tau'] + 0.5*log_vars['delta'] + 0.5*log_vars['gc']),
    ("CN×√δ/√ξ×√τ", lambda: np.log(SIGMA_AM) + log_vars['cn'] + 0.5*log_vars['delta'] - 0.5*log_vars['ratio'] + 0.5*log_vars['tau']),
    ("(φ-φc)^0.5×CN×√δ/√ξ", lambda: np.log(SIGMA_AM) + 0.5*log_vars['(φ-φc)'] + log_vars['cn'] + 0.5*log_vars['delta'] - 0.5*log_vars['ratio']),
]

# Add top results from each phase
if results:
    top1 = results[0]
    candidates.append((f"Phase1 best: {top1[1]}", None))
if results2:
    top2 = results2[0]
    candidates.append((f"Phase2 best: {top2[1]}", None))
if results3:
    top3 = results3[0]
    candidates.append((f"Phase3 best: {top3[1]}", None))

for name, fn in candidates:
    if fn is None:
        continue
    rhs = fn()
    r2, C = fit_r2(s, rhs)
    pred = np.exp(np.log(C) + rhs)
    actual = np.exp(s)
    errs = np.abs(pred - actual) / actual * 100
    print(f"\n  {name}: R²={r2:.4f}, C={C:.4f}, |err|={np.mean(errs):.1f}%")
    for i, r in enumerate(thin):
        print(f"    T/d={r['ratio']:.1f} σ_pred={pred[i]:.2f} σ_act={actual[i]:.2f} err={errs[i]:.0f}%")

print("\n" + "="*60)
print("DONE")
