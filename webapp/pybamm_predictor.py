"""
PyBaMM ASSB Full Electrochemical Model
=======================================
Doyle-Fuller-Newman (DFN) model customized for All-Solid-State Battery.

Key ASSB modifications:
- Solid electrolyte: constant conductivity (no concentration dependence)
- Transference number = 1 (single-ion conductor)
- No electrolyte diffusion limitation
- Interface kinetics: lower j₀ than liquid LIB

Outputs: voltage curves, capacity, utilization, rate capability
"""

import numpy as np
import json


def run_assb_simulation(sigma_ionic, phi_se, phi_am, thickness, d_am,
                         tau=1.5, porosity=15, c_rates=None, temperature=298):
    """
    Run full electrochemical simulation for ASSB composite cathode.

    Args:
        sigma_ionic: effective ionic conductivity (mS/cm) from DEM/scaling law
        phi_se: SE volume fraction
        phi_am: AM volume fraction
        thickness: electrode thickness (μm)
        d_am: AM particle diameter (μm)
        tau: tortuosity
        porosity: porosity (%)
        c_rates: list of C-rates
        temperature: K

    Returns:
        dict with voltage curves, capacity, utilization per C-rate
    """
    if c_rates is None:
        c_rates = [0.05, 0.1, 0.2, 0.5, 1.0, 2.0]

    try:
        import pybamm
        return _run_pybamm_dfn(pybamm, sigma_ionic, phi_se, phi_am, thickness,
                                d_am, tau, porosity, c_rates, temperature)
    except ImportError:
        print("PyBaMM not available. Using Newman analytical model.")
        return _newman_analytical(sigma_ionic, phi_am, thickness, d_am, c_rates)


def _run_pybamm_dfn(pybamm, sigma_ionic, phi_se, phi_am, thickness,
                     d_am, tau, porosity, c_rates, temperature):
    """Full DFN model with ASSB parameters."""

    T_m = thickness * 1e-6       # μm → m
    r_AM_m = d_am * 1e-6 / 2    # radius in m
    kappa_SI = sigma_ionic * 0.1  # mS/cm → S/m
    pore_fraction = porosity / 100

    results = {}

    for c_rate in c_rates:
        try:
            # Use DFN (Doyle-Fuller-Newman) — full porous electrode model
            model = pybamm.lithium_ion.DFN()
            param = model.default_parameter_values

            # ═══ POSITIVE ELECTRODE (ASSB Composite Cathode) ═══
            param["Positive electrode thickness [m]"] = T_m
            param["Positive particle radius [m]"] = r_AM_m
            param["Positive electrode porosity"] = max(0.01, pore_fraction)
            param["Positive electrode active material volume fraction"] = phi_am
            param["Positive electrode Bruggeman coefficient (electrolyte)"] = tau

            # ═══ ELECTROLYTE: ASSB Solid Electrolyte ═══
            # Constant conductivity (not concentration-dependent like liquid)
            param["Electrolyte conductivity [S.m-1]"] = kappa_SI

            # High diffusivity → no concentration polarization in solid SE
            param["Electrolyte diffusivity [m2.s-1]"] = 1e-9

            # Transference number = 1 for single-ion conductor (Li⁺ in LPSCl)
            # In liquid: t+ ≈ 0.3-0.4. In solid SE: t+ = 1
            param["Cation transference number"] = 1.0

            # Thermodynamic factor = 1 (ideal solution, no activity correction)
            param["Thermodynamic factor"] = 1.0

            # ═══ CELL GEOMETRY (ASSB) ═══
            param["Negative electrode thickness [m]"] = 20e-6  # thin Li metal
            param["Separator thickness [m]"] = 30e-6  # thin SE separator

            # ═══ C-RATE ═══
            param["Current function [A]"] = param["Nominal cell capacity [A.h]"] * c_rate

            # ═══ SOLVE ═══
            sim = pybamm.Simulation(model, parameter_values=param)
            t_end = min(3600 / c_rate * 1.1, 36000)  # slightly over to catch cutoff
            sol = sim.solve([0, t_end])

            # ═══ EXTRACT RESULTS ═══
            t = sol["Time [s]"].entries
            V = sol["Voltage [V]"].entries
            Q = sol["Discharge capacity [A.h]"].entries
            nominal = param["Nominal cell capacity [A.h]"]

            # Utilization = total discharged / nominal
            util = Q[-1] / nominal if nominal > 0 else 0

            # Find time indices for voltage curve sampling
            n_points = min(100, len(t))
            idx = np.linspace(0, len(t)-1, n_points, dtype=int)

            results[f'{c_rate}C'] = {
                'utilization': round(float(min(1.0, max(0.0, util))), 3),
                'capacity_mAh': round(float(Q[-1] * 1000), 2),
                'voltage_final': round(float(V[-1]), 3),
                'voltage_mean': round(float(np.mean(V)), 3),
                'energy_Wh': round(float(np.trapz(V, Q)), 4),
                'voltage_curve': {
                    'time_s': [round(float(t[i]), 1) for i in idx],
                    'voltage_V': [round(float(V[i]), 4) for i in idx],
                    'capacity_mAh': [round(float(Q[i]*1000), 2) for i in idx],
                },
                'method': 'PyBaMM_DFN',
            }

        except Exception as e:
            error_msg = str(e)[:100]
            print(f"  PyBaMM {c_rate}C failed: {error_msg}")
            # Newman fallback
            util = _newman_single(sigma_ionic, phi_am, thickness, d_am, c_rate)
            results[f'{c_rate}C'] = {
                'utilization': util,
                'method': f'Newman_fallback',
                'error': error_msg,
            }

    # ═══ RATE CAPABILITY SUMMARY ═══
    summary = {
        'c_rates': list(results.keys()),
        'utilizations': {k: v['utilization'] for k, v in results.items()},
        'has_voltage_curves': any('voltage_curve' in v for v in results.values()),
    }

    return {'per_crate': results, 'summary': summary}


def _newman_analytical(sigma_ionic, phi_am, thickness, d_am, c_rates):
    """Newman analytical model fallback."""
    results = {}
    for c_rate in c_rates:
        util = _newman_single(sigma_ionic, phi_am, thickness, d_am, c_rate)
        results[f'{c_rate}C'] = {
            'utilization': util,
            'method': 'Newman_analytical',
        }
    return {'per_crate': results, 'summary': {
        'c_rates': list(results.keys()),
        'utilizations': {k: v['utilization'] for k, v in results.items()},
    }}


def _newman_single(sigma_ionic, phi_am, thickness, d_am, c_rate):
    """Corrected Newman utilization (ASSB parameters)."""
    F = 96485; R = 8.314; T_K = 298

    r_AM_cm = d_am * 1e-4 / 2 if d_am > 0 else 2.5e-4
    T_cm = thickness * 1e-4 if thickness > 0 else 0.01
    kappa = sigma_ionic * 1e-3 if sigma_ionic > 0 else 1e-6

    a_s_geo = 3 * phi_am / r_AM_cm if r_AM_cm > 0 and phi_am > 0 else 1e4
    a_s = a_s_geo * 0.20  # coverage ~20%
    j0 = 0.01e-3  # A/cm², ASSB interface

    if kappa > 0 and a_s > 0:
        nu_sq = a_s * j0 * F * T_cm**2 / (kappa * R * T_K)
        nu = np.sqrt(nu_sq) if nu_sq > 0 else 0
        util_newman = np.tanh(nu) / nu if nu > 0.01 else 1.0
    else:
        util_newman = 1.0

    Q_vol = 200 * 4.8 * phi_am * 1e-3
    i_app = c_rate * Q_vol * T_cm
    delta_V = i_app * T_cm / (2 * kappa) if kappa > 0 else 10
    ohmic_factor = min(1.0, 0.3 / delta_V) if delta_V > 0 else 1.0

    return round(min(1.0, max(0.0, util_newman * ohmic_factor)), 3)


if __name__ == '__main__':
    print("=" * 60)
    print("PyBaMM ASSB Full Electrochemical Simulation")
    print("=" * 60)

    tests = [
        {'name': '후막 80:20 (8mAh)', 'sigma': 0.19, 'phi_se': 0.30, 'phi_am': 0.52,
         'T': 160, 'd_am': 10, 'tau': 1.22, 'porosity': 17},
        {'name': '박막 80:25 (1mAh)', 'sigma': 0.17, 'phi_se': 0.31, 'phi_am': 0.52,
         'T': 20, 'd_am': 10, 'tau': 1.87, 'porosity': 16},
        {'name': 'SE 1.5μm Real', 'sigma': 0.08, 'phi_se': 0.29, 'phi_am': 0.52,
         'T': 115, 'd_am': 10, 'tau': 1.49, 'porosity': 16},
        {'name': '후막 75:25 (8mAh)', 'sigma': 0.26, 'phi_se': 0.36, 'phi_am': 0.46,
         'T': 183, 'd_am': 10, 'tau': 1.20, 'porosity': 18},
    ]

    for tc in tests:
        print(f"\n{'─'*50}")
        print(f"{tc['name']}: σ_ion={tc['sigma']}mS/cm, T={tc['T']}μm, τ={tc['tau']}")
        print(f"{'─'*50}")

        result = run_assb_simulation(
            sigma_ionic=tc['sigma'], phi_se=tc['phi_se'], phi_am=tc['phi_am'],
            thickness=tc['T'], d_am=tc['d_am'], tau=tc['tau'], porosity=tc['porosity'],
            c_rates=[0.05, 0.1, 0.2, 0.5, 1.0, 2.0]
        )

        for cr, r in result['per_crate'].items():
            util_pct = r['utilization'] * 100
            method = r['method']
            extra = ''
            if 'voltage_mean' in r:
                extra = f"V_mean={r['voltage_mean']}V, E={r.get('energy_Wh', 0):.4f}Wh"
            elif 'error' in r:
                extra = r['error'][:50]
            print(f"  {cr:6s}: util={util_pct:5.1f}% [{method}] {extra}")

    print(f"\n{'='*60}")
    print("Done!")
