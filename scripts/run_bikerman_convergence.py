#!/usr/bin/env python3
"""
Grid convergence study for Bikerman model.
Compares numerical surface charge with analytical solution.
"""

import subprocess
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys
import re

sys.path.insert(0, str(Path(__file__).parent.parent / 'styles'))
from plot_style import setup_plot_style, setup_axis_style, set_labels, FIGURE_SIZES

# Physical constants (same as solver)
eps_0 = 8.8541878128e-12  # F/m
kB = 1.380649e-23         # J/K
e = 1.602176634e-19       # C
NA = 6.02214076e23        # 1/mol

def bikerman_analytical_sigma(phi0_mV, c0_M=1.0, eps_r=12, T=298.15, a_nm=0.7):
    """
    Analytical Bikerman surface charge density.

    From first integral of Bikerman equation:
    (dψ/dξ)² = (2/ν) ln(1 - ν + ν cosh(ψ))

    Surface charge:
    σ = sqrt(2ε kB T c0 NA) * sqrt((2/ν) ln(1 - ν + ν cosh(ψ0)))

    where:
    - ν = 2 a³ c0 NA (packing fraction)
    - ψ0 = e φ0 / (kB T) (dimensionless potential)

    For dual electrode (closed system), surface-to-bulk potential is φ0/2.
    """
    eps = eps_0 * eps_r
    c0 = c0_M * 1000  # mol/m³
    a = a_nm * 1e-9   # m

    # Packing fraction
    nu = 2 * a**3 * c0 * NA

    # Dimensionless surface potential (surface to bulk = half of applied)
    psi0 = (phi0_mV * 1e-3 / 2) / (kB * T / e)

    # Analytical surface charge from first integral
    # σ = sqrt(2ε kB T c0 NA) * sqrt((2/ν) ln(g(ψ0)))
    # where g(ψ) = 1 - ν + ν cosh(ψ)
    prefactor = np.sqrt(2 * eps * kB * T * c0 * NA)
    g_psi0 = 1 - nu + nu * np.cosh(psi0)
    log_term = (2.0 / nu) * np.log(g_psi0)
    sigma = prefactor * np.sqrt(log_term)

    # Convert C/m² to μC/cm²: 1 C/m² = 1e6 μC / 1e4 cm² = 100 μC/cm²
    return sigma * 100


def run_solver(phi0_mV, N, model='bikerman', ion_size=0.7, stretch=3.0):
    """Run PNP solver and extract surface charge density."""
    cmd = [
        './build/pnp_solver',
        '--phi0', str(phi0_mV),
        '--phi-right', '0',
        '--closed-system',
        '--dual-electrode',
        '--c0', '1.0',
        '--L', '50',
        '--N', str(N),
        '--stretch', str(stretch),
        '--model', model,
        '--ion-size', str(ion_size)
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    output = result.stdout

    # Extract surface charge density
    match = re.search(r'Left electrode:\s+([+-]?\d+\.?\d*(?:e[+-]?\d+)?)\s*μC/cm²', output)
    if match:
        return float(match.group(1))
    return None


def main():
    setup_plot_style()

    # Test parameters
    phi0_mV = 100  # Applied voltage
    c0_M = 1.0
    eps_r = 12
    a_nm = 0.7
    stretch = 3.0  # Use non-uniform grid

    # Grid sizes for convergence study
    N_values = [101, 201, 401, 801, 1601, 3201]

    # Analytical solution
    sigma_analytical = bikerman_analytical_sigma(phi0_mV, c0_M, eps_r, a_nm=a_nm)

    print("=" * 70)
    print("  Bikerman Model Grid Convergence Study")
    print("=" * 70)
    print(f"  Applied voltage: {phi0_mV} mV")
    print(f"  Bulk concentration: {c0_M} M")
    print(f"  Ion size: {a_nm} nm")
    print(f"  Grid stretching: {stretch}")
    print(f"  Analytical σ: {sigma_analytical:.6f} μC/cm²")
    print("=" * 70)
    print(f"{'N':>8} {'σ_num [μC/cm²]':>18} {'Error [%]':>12} {'Order':>8}")
    print("-" * 50)

    sigma_values = []
    errors = []

    for N in N_values:
        sigma_num = run_solver(phi0_mV, N, model='bikerman', ion_size=a_nm, stretch=stretch)
        if sigma_num is not None:
            error = (sigma_num - sigma_analytical) / sigma_analytical * 100
            sigma_values.append(sigma_num)
            errors.append(abs(error))

            # Calculate convergence order
            if len(errors) > 1:
                # Order = log2(error_prev / error_curr) when grid is halved
                order = np.log2(errors[-2] / errors[-1]) if errors[-1] > 1e-10 else float('nan')
                print(f"{N:>8} {sigma_num:>18.6f} {error:>+12.4f} {order:>8.2f}")
            else:
                print(f"{N:>8} {sigma_num:>18.6f} {error:>+12.4f} {'—':>8}")

    print("-" * 50)
    print(f"Analytical: {sigma_analytical:.6f} μC/cm²")

    # Calculate average convergence order (excluding first)
    if len(errors) > 2:
        orders = []
        for i in range(1, len(errors)):
            if errors[i] > 1e-10:
                orders.append(np.log2(errors[i-1] / errors[i]))
        if orders:
            avg_order = np.mean(orders)
            print(f"Average convergence order: {avg_order:.2f}")

    # Plot convergence
    fig, ax = plt.subplots(figsize=FIGURE_SIZES['single'])

    N_arr = np.array(N_values[:len(errors)])
    err_arr = np.array(errors)

    ax.loglog(N_arr, err_arr, 'o-', color='C0', linewidth=1.5, markersize=6, label='Bikerman')

    # Reference lines
    N_ref = np.array([N_arr[0], N_arr[-1]])
    # 2nd order reference
    err_ref_2 = err_arr[0] * (N_ref[0] / N_ref)**2
    ax.loglog(N_ref, err_ref_2, '--', color='gray', linewidth=1, label='2nd order')

    set_labels(ax, r'Grid points $N$', r'Relative error (\%)')
    setup_axis_style(ax)
    ax.legend(loc='upper right', frameon=False)
    ax.set_title(f'Bikerman convergence ($\\phi_0$ = {phi0_mV} mV)', fontsize=10)

    plt.tight_layout()

    results_dir = Path(__file__).parent.parent / 'results'
    plt.savefig(results_dir / 'bikerman_convergence.png', dpi=300, bbox_inches='tight')
    plt.savefig(results_dir / 'bikerman_convergence.svg', bbox_inches='tight')
    print(f"\nSaved: results/bikerman_convergence.png")

    # Save data
    data = np.column_stack([N_arr, np.array(sigma_values), err_arr])
    np.savetxt(results_dir / 'bikerman_convergence.csv', data,
               header='N,Sigma_numerical,Error_percent', delimiter=',', comments='')
    print(f"Saved: results/bikerman_convergence.csv")


if __name__ == '__main__':
    main()
