#!/usr/bin/env python3
"""
Plot surface charge density vs applied voltage.
Runs the PNP solver at multiple voltages and extracts charge from output.
"""

import subprocess
import numpy as np
import matplotlib.pyplot as plt
import re
from pathlib import Path
import sys

# Add styles directory to path
sys.path.insert(0, str(Path(__file__).parent.parent / 'styles'))
from plot_style import setup_plot_style, setup_axis_style, set_labels, FIGURE_SIZES

def run_solver(phi0_mV, phi_right_mV=0, c0=1.0, L=50, N=1001):
    """Run PNP solver and extract surface charge density."""
    cmd = [
        './build/pnp_solver',
        '--phi0', str(phi0_mV),
        '--phi-right', str(phi_right_mV),
        '--closed-system',
        '--c0', str(c0),
        '--L', str(L),
        '--N', str(N),
        '--dual-electrode'
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    output = result.stdout

    # Extract surface charge density from output
    # Looking for: "Left electrode:  +5.12994 μC/cm²"
    match = re.search(r'Left electrode:\s+([+-]?\d+\.?\d*(?:e[+-]?\d+)?)\s*μC/cm²', output)
    if match:
        sigma = float(match.group(1))
        return sigma
    else:
        print(f"Warning: Could not extract charge for phi0={phi0_mV} mV")
        print(output)
        return None

def main():
    setup_plot_style()

    # Voltage range: 0 to 200 mV (applied between electrodes)
    # For closed system with phi_right=0, the EDL potential is phi0/2
    voltages = np.linspace(0, 200, 21)  # 0 to 200 mV in 10 mV steps

    print("=" * 60)
    print("  Voltage-Charge Characteristic")
    print("=" * 60)
    print(f"{'Voltage [mV]':>15} {'σ [μC/cm²]':>15}")
    print("-" * 35)

    charges = []
    valid_voltages = []

    for V in voltages:
        if V == 0:
            # At zero voltage, charge is zero
            sigma = 0.0
        else:
            sigma = run_solver(V, phi_right_mV=0)

        if sigma is not None:
            charges.append(sigma)
            valid_voltages.append(V)
            print(f"{V:>15.1f} {sigma:>15.4f}")

    print("-" * 35)

    # Convert to numpy arrays
    V_array = np.array(valid_voltages)
    sigma_array = np.array(charges)

    # Create plot
    fig, ax = plt.subplots(figsize=FIGURE_SIZES['single'])

    ax.plot(V_array, sigma_array, 'o-', color='C0', linewidth=1.5, markersize=5)

    # Add reference line for linear capacitor (using C at low voltage)
    # C ≈ ε₀εᵣ/λ_D for Gouy-Chapman at low potential
    # At 1M, λ_D ≈ 0.119 nm, so C ≈ 89 μF/cm² (single EDL)
    # For series: C_total ≈ 44.5 μF/cm²
    # But we use numerical value from low-voltage slope
    if len(V_array) > 2:
        # Linear fit at low voltages
        low_idx = V_array < 50
        if np.sum(low_idx) > 2:
            slope, intercept = np.polyfit(V_array[low_idx], sigma_array[low_idx], 1)
            V_linear = np.array([0, V_array.max()])
            sigma_linear = slope * V_linear + intercept
            ax.plot(V_linear, sigma_linear, '--', color='gray', linewidth=1,
                    label=f'Linear (C = {slope*1000:.1f} μF/cm²)')
            ax.legend(loc='lower right', frameon=False)

    set_labels(ax, r'Applied voltage $\Delta\phi$ (mV)', r'Surface charge $\sigma$ ($\mu$C/cm$^2$)')
    setup_axis_style(ax)

    ax.set_xlim(0, V_array.max() * 1.05)
    ax.set_ylim(0, sigma_array.max() * 1.1)

    plt.tight_layout()

    # Save
    results_dir = Path(__file__).parent.parent / 'results'
    plt.savefig(results_dir / 'voltage_charge.png', dpi=300, bbox_inches='tight')
    plt.savefig(results_dir / 'voltage_charge.svg', bbox_inches='tight')
    print(f"\nSaved: results/voltage_charge.png")
    print(f"Saved: results/voltage_charge.svg")

    # Also save data
    data_file = results_dir / 'voltage_charge.csv'
    np.savetxt(data_file, np.column_stack([V_array, sigma_array]),
               header='Voltage_mV,Sigma_uC_cm2', delimiter=',', comments='')
    print(f"Saved: results/voltage_charge.csv")

if __name__ == '__main__':
    main()
