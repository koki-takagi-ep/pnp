#!/usr/bin/env python3
"""
Compare voltage-charge characteristics: Standard PB vs Bikerman model.
Shows how steric effects cause capacitance saturation at high voltages.
"""

import subprocess
import numpy as np
import matplotlib.pyplot as plt
import re
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent / 'styles'))
from plot_style import setup_plot_style, setup_axis_style, set_labels, FIGURE_SIZES

def run_solver(phi0_mV, model='standard', ion_size=0.7):
    """Run PNP solver and extract surface charge density."""
    cmd = [
        './build/pnp_solver',
        '--phi0', str(phi0_mV),
        '--phi-right', '0',
        '--closed-system',
        '--c0', '1.0',
        '--L', '50',
        '--N', '4001',  # Increased for better accuracy at high voltages
        '--dual-electrode',
        '--model', model
    ]
    if model == 'bikerman':
        cmd.extend(['--ion-size', str(ion_size)])

    result = subprocess.run(cmd, capture_output=True, text=True)
    output = result.stdout

    match = re.search(r'Left electrode:\s+([+-]?\d+\.?\d*(?:e[+-]?\d+)?)\s*μC/cm²', output)
    if match:
        return float(match.group(1))
    return None

def gouy_chapman_analytical(phi_mV, c0_M=1.0, eps_r=12, T=298.15):
    """Analytical Gouy-Chapman surface charge density.

    σ = √(8 ε₀ εᵣ kB T c₀ NA) × sinh(eψ₀ / 2kBT)

    where ψ₀ is the surface-to-bulk potential (= Δφ/2 for symmetric dual electrode).
    """
    # Use same constants as the solver (pnp_solver.hpp)
    eps_0 = 8.8541878128e-12
    kB = 1.380649e-23
    e = 1.602176634e-19
    NA = 6.02214076e23

    c0 = c0_M * 1000  # mol/m³
    psi0 = phi_mV * 1e-3 / 2  # Surface-to-bulk potential (half of applied)
    phi_T = kB * T / e  # Thermal voltage

    # Gouy-Chapman: sinh argument is eψ₀/(2kBT) = ψ₀/(2φT)
    sigma = np.sqrt(8 * eps_0 * eps_r * kB * T * c0 * NA) * np.sinh(psi0 / (2 * phi_T))
    return sigma * 1e2  # Convert to μC/cm²

def main():
    setup_plot_style()

    voltages = np.linspace(0, 300, 31)

    print("=" * 70)
    print("  Voltage-Charge Comparison: Standard PB vs Bikerman")
    print("=" * 70)

    # Standard PB
    print("\n[Standard Poisson-Boltzmann]")
    sigma_std = []
    for V in voltages:
        if V == 0:
            sigma = 0.0
        else:
            sigma = run_solver(V, model='standard')
        sigma_std.append(sigma if sigma else 0)
        print(f"  {V:6.0f} mV: {sigma_std[-1]:8.3f} μC/cm²")

    # Bikerman model
    print("\n[Bikerman Model (a = 0.7 nm)]")
    sigma_bik = []
    for V in voltages:
        if V == 0:
            sigma = 0.0
        else:
            sigma = run_solver(V, model='bikerman', ion_size=0.7)
        sigma_bik.append(sigma if sigma else 0)
        print(f"  {V:6.0f} mV: {sigma_bik[-1]:8.3f} μC/cm²")

    # Analytical Gouy-Chapman
    sigma_gc = [gouy_chapman_analytical(V) for V in voltages]

    # Convert to arrays
    V_arr = np.array(voltages)
    sigma_std = np.array(sigma_std)
    sigma_bik = np.array(sigma_bik)
    sigma_gc = np.array(sigma_gc)

    # Plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.4, 3.7))

    # Left: σ vs V
    ax1.plot(V_arr, sigma_gc, '--', color='gray', linewidth=1, label='Gouy-Chapman (analytical)')
    ax1.plot(V_arr, sigma_std, 'o-', color='C0', linewidth=1.5, markersize=4, label='Standard PB (numerical)')
    ax1.plot(V_arr, sigma_bik, 's-', color='C1', linewidth=1.5, markersize=4, label='Bikerman (a=0.7 nm)')

    set_labels(ax1, r'Applied voltage $\Delta\phi$ (mV)', r'Surface charge $\sigma$ ($\mu$C/cm$^2$)')
    setup_axis_style(ax1)
    ax1.legend(loc='upper left', frameon=False, fontsize=8)
    ax1.set_xlim(0, 310)
    ax1.set_ylim(0, None)

    # Right: Differential capacitance (dσ/dV)
    dV = V_arr[1] - V_arr[0]
    C_std = np.gradient(sigma_std, dV) * 1000  # μC/cm² per V = μF/cm²
    C_bik = np.gradient(sigma_bik, dV) * 1000
    C_gc = np.gradient(sigma_gc, dV) * 1000

    ax2.plot(V_arr[1:-1], C_gc[1:-1], '--', color='gray', linewidth=1, label='Gouy-Chapman')
    ax2.plot(V_arr[1:-1], C_std[1:-1], 'o-', color='C0', linewidth=1.5, markersize=4, label='Standard PB')
    ax2.plot(V_arr[1:-1], C_bik[1:-1], 's-', color='C1', linewidth=1.5, markersize=4, label='Bikerman')

    set_labels(ax2, r'Applied voltage $\Delta\phi$ (mV)', r'Diff. capacitance $C$ ($\mu$F/cm$^2$)')
    setup_axis_style(ax2)
    ax2.legend(loc='upper left', frameon=False, fontsize=8)
    ax2.set_xlim(0, 310)
    ax2.set_ylim(0, None)

    plt.tight_layout()

    results_dir = Path(__file__).parent.parent / 'results'
    plt.savefig(results_dir / 'voltage_charge_comparison.png', dpi=300, bbox_inches='tight')
    plt.savefig(results_dir / 'voltage_charge_comparison.svg', bbox_inches='tight')
    print(f"\nSaved: results/voltage_charge_comparison.png")

    # Save data
    data = np.column_stack([V_arr, sigma_std, sigma_bik, sigma_gc])
    np.savetxt(results_dir / 'voltage_charge_comparison.csv', data,
               header='Voltage_mV,Sigma_Standard,Sigma_Bikerman,Sigma_GC_analytical',
               delimiter=',', comments='')

if __name__ == '__main__':
    main()
