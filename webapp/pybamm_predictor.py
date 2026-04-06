"""
PyBaMM-based electrode utilization prediction.
Uses Doyle-Fuller-Newman (DFN) model with DEM-derived parameters.

Usage:
  from pybamm_predictor import predict_with_pybamm
  result = predict_with_pybamm(sigma_ionic=0.15, phi_se=0.3, phi_am=0.5,
                                thickness=150, d_am=5, c_rates=[0.1, 0.5, 1, 2])
"""

import numpy as np

def predict_with_pybamm(sigma_ionic, phi_se, phi_am, thickness, d_am,
                         tau=1.5, c_rates=None, temperature=298):
    """
    Run PyBaMM SPMe model with DEM-derived electrode parameters.

    Args:
        sigma_ionic: effective ionic conductivity (mS/cm)
        phi_se: SE volume fraction
        phi_am: AM volume fraction
        thickness: electrode thickness (μm)
        d_am: AM particle diameter (μm)
        tau: tortuosity
        c_rates: list of C-rates to evaluate
        temperature: K

    Returns:
        dict with utilization per C-rate, capacity, voltage curves
    """
    if c_rates is None:
        c_rates = [0.1, 0.2, 0.5, 1.0, 2.0]

    try:
        import pybamm
    except ImportError:
        return _newman_fallback(sigma_ionic, phi_am, thickness, d_am, tau, c_rates)

    # Convert units
    T_m = thickness * 1e-6       # μm → m
    r_AM_m = d_am * 1e-6 / 2    # μm → m, radius
    kappa_eff = sigma_ionic * 0.1  # mS/cm → S/m

    results = {}

    for c_rate in c_rates:
        try:
            # Use SPM (not SPMe) since we override electrolyte transport
            # SPMe electrolyte corrections are for liquid; ASSB has no
            # concentration gradients in solid electrolyte
            model = pybamm.lithium_ion.SPM()

            # Update parameters with DEM values
            param = model.default_parameter_values

            # ── Positive electrode (our ASSB composite cathode) ──
            param["Positive electrode thickness [m]"] = T_m
            param["Positive particle radius [m]"] = r_AM_m
            param["Positive electrode porosity"] = max(0.01, 1 - phi_am - phi_se)
            param["Positive electrode active material volume fraction"] = phi_am
            param["Positive electrode Bruggeman coefficient (electrolyte)"] = tau

            # ── Negative electrode: thin Li metal anode for ASSB ──
            param["Negative electrode thickness [m]"] = 20e-6  # 20 μm Li metal

            # ── Separator: thin SE separator layer ──
            param["Separator thickness [m]"] = 30e-6  # 30 μm SE separator

            # ── Electrolyte: override with ASSB solid electrolyte values ──
            # Conductivity: use our DEM-derived σ_ionic (NOT liquid LiPF6)
            kappa_val = sigma_ionic * 0.1  # mS/cm -> S/m
            param["Electrolyte conductivity [S.m-1]"] = kappa_val

            # Diffusivity: ASSB solid electrolyte has no concentration gradient
            # Use high value to effectively eliminate concentration polarization
            param["Electrolyte diffusivity [m2.s-1]"] = 1e-10  # effectively instant in solid

            # ── Set C-rate ──
            param["Current function [A]"] = param["Nominal cell capacity [A.h]"] * c_rate

            # ── Solve ──
            sim = pybamm.Simulation(model, parameter_values=param)
            sol = sim.solve([0, 3600 / c_rate])  # discharge time

            # ── Extract capacity utilization ──
            # Use discharge capacity relative to theoretical capacity
            discharge_capacity = sol["Discharge capacity [A.h]"].entries[-1]
            nominal = param["Nominal cell capacity [A.h]"]

            # Calculate theoretical Q per area for our electrode
            # NCM811: 200 mAh/g, rho_AM = 4.8 g/cm³
            Q_theo_mAh_g = 200
            rho_AM = 4.8e6  # g/m³
            electrode_area = param.get("Electrode width [m]", 1.0) * param.get("Electrode height [m]", 1.0)
            Q_per_area = Q_theo_mAh_g * 1e-3 * rho_AM * phi_am * T_m  # Ah/m²

            # Use nominal capacity ratio as utilization
            util = discharge_capacity / nominal if nominal > 0 else 0

            # Voltage curve
            t = sol["Time [s]"].entries
            V = sol["Voltage [V]"].entries

            results[f'{c_rate}C'] = {
                'utilization': round(float(min(1.0, max(0.0, util))), 3),
                'capacity_Ah': round(float(discharge_capacity), 4),
                'voltage_final': round(float(V[-1]), 3),
                'method': 'PyBaMM_SPM',
            }

        except Exception as e:
            # Fallback to Newman if PyBaMM fails for this C-rate
            # Each C-rate is independent — one failure doesn't kill all
            results[f'{c_rate}C'] = {
                'utilization': _newman_single(sigma_ionic, phi_am, thickness, d_am, c_rate),
                'method': f'Newman_fallback ({str(e)[:50]})',
            }

    return results


def _newman_fallback(sigma_ionic, phi_am, thickness, d_am, tau, c_rates):
    """Newman analytical model when PyBaMM not available."""
    results = {}
    for c_rate in c_rates:
        util = _newman_single(sigma_ionic, phi_am, thickness, d_am, c_rate)
        results[f'{c_rate}C'] = {
            'utilization': util,
            'method': 'Newman_analytical',
        }
    return results


def _newman_single(sigma_ionic, phi_am, thickness, d_am, c_rate):
    """Single C-rate Newman utilization calculation."""
    F = 96485
    R = 8.314
    T_K = 298

    r_AM_cm = d_am * 1e-4 / 2 if d_am > 0 else 2.5e-4
    T_cm = thickness * 1e-4 if thickness > 0 else 0.01
    kappa = sigma_ionic * 1e-3 if sigma_ionic > 0 else 1e-6  # mS/cm → S/cm

    a_s_geometric = 3 * phi_am / r_AM_cm if r_AM_cm > 0 and phi_am > 0 else 1e4
    coverage_fraction = 0.20  # from DEM: ~20% of AM surface contacts SE
    a_s = a_s_geometric * coverage_fraction  # effective interfacial area
    j0 = 0.01e-3  # A/cm² (0.01 mA/cm², ASSB NCM/LPSCl interface)

    # Newman: ν = T × √(a_s × j₀ × F / (κ × R × T))
    if kappa > 0 and a_s > 0:
        nu_sq = a_s * j0 * F * T_cm**2 / (kappa * R * T_K)
        nu = np.sqrt(nu_sq) if nu_sq > 0 else 0
        util_newman = np.tanh(nu) / nu if nu > 0.01 else 1.0
    else:
        util_newman = 1.0

    # Ohmic C-rate penalty
    Q_vol = 200 * 4.8 * phi_am * 1e-3  # Ah/cm³
    i_app = c_rate * Q_vol * T_cm
    delta_V = i_app * T_cm / (2 * kappa) if kappa > 0 else 10
    ohmic_factor = min(1.0, 0.3 / delta_V) if delta_V > 0 else 1.0

    return round(min(1.0, max(0.0, util_newman * ohmic_factor)), 3)


if __name__ == '__main__':
    print("=" * 50)
    print("PyBaMM Electrode Utilization Test")
    print("=" * 50)

    # Test with typical DEM values
    test_cases = [
        {'name': '후막 80:20', 'sigma': 0.19, 'phi_se': 0.30, 'phi_am': 0.52, 'T': 160, 'd_am': 10},
        {'name': '박막 80:25', 'sigma': 0.17, 'phi_se': 0.31, 'phi_am': 0.52, 'T': 20, 'd_am': 10},
        {'name': 'SE 1.5μm',  'sigma': 0.08, 'phi_se': 0.29, 'phi_am': 0.52, 'T': 115, 'd_am': 10},
    ]

    for tc in test_cases:
        print(f"\n{tc['name']}: σ={tc['sigma']}, T={tc['T']}μm")
        result = predict_with_pybamm(
            sigma_ionic=tc['sigma'], phi_se=tc['phi_se'], phi_am=tc['phi_am'],
            thickness=tc['T'], d_am=tc['d_am']
        )
        for cr, r in result.items():
            print(f"  {cr}: util={r['utilization']:.1%} [{r['method']}]")
