#!/usr/bin/env python3
"""
Comprehensive grid convergence study for Bikerman model.
1. Surface charge error vs analytical solution
2. L2 error of potential profile using Richardson extrapolation
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
    """Analytical Bikerman surface charge density."""
    eps = eps_0 * eps_r
    c0 = c0_M * 1000  # mol/m³
    a = a_nm * 1e-9   # m

    nu = 2 * a**3 * c0 * NA
    psi0 = (phi0_mV * 1e-3 / 2) / (kB * T / e)

    prefactor = np.sqrt(2 * eps * kB * T * c0 * NA)
    g_psi0 = 1 - nu + nu * np.cosh(psi0)
    log_term = (2.0 / nu) * np.log(g_psi0)
    sigma = prefactor * np.sqrt(log_term)

    return sigma * 100  # Convert to μC/cm²


def run_solver(phi0_mV, N, model='bikerman', ion_size=0.7, stretch=3.0, output_file=None):
    """Run PNP solver and return output file path."""
    if output_file is None:
        output_file = f'results/conv_bikerman_{N}.dat'

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
        '--ion-size', str(ion_size),
        '--output', output_file
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    # Extract surface charge
    match = re.search(r'Left electrode:\s+([+-]?\d+\.?\d*(?:e[+-]?\d+)?)\s*μC/cm²', result.stdout)
    sigma = float(match.group(1)) if match else None

    return output_file, sigma


def load_potential_profile(filename):
    """Load potential profile from solver output."""
    data = np.loadtxt(filename, comments='#')
    x = data[:, 0]      # nm
    phi = data[:, 2]    # mV
    return x, phi


def interpolate_to_grid(x_fine, phi_fine, x_coarse):
    """Interpolate fine solution to coarse grid points."""
    return np.interp(x_coarse, x_fine, phi_fine)


def compute_l2_error(phi_num, phi_ref, dx=None):
    """Compute L2 error."""
    diff = phi_num - phi_ref
    if dx is not None:
        # Weighted L2 norm
        l2 = np.sqrt(np.sum(diff**2 * dx) / np.sum(dx))
    else:
        # Simple RMS
        l2 = np.sqrt(np.mean(diff**2))
    return l2


def main():
    setup_plot_style()

    # Test parameters
    phi0_mV = 100
    c0_M = 1.0
    eps_r = 12
    a_nm = 0.7
    stretch = 3.0

    # Grid sizes
    N_values = [51, 101, 201, 401, 801, 1601, 3201, 6401]

    # Analytical surface charge
    sigma_analytical = bikerman_analytical_sigma(phi0_mV, c0_M, eps_r, a_nm=a_nm)

    print("=" * 80)
    print("  Bikerman Model Comprehensive Convergence Study")
    print("=" * 80)
    print(f"  Applied voltage: {phi0_mV} mV")
    print(f"  Ion size: {a_nm} nm")
    print(f"  Analytical σ: {sigma_analytical:.6f} μC/cm²")
    print("=" * 80)

    # Run all simulations
    results = {}
    for N in N_values:
        print(f"Running N = {N}...", end=" ", flush=True)
        output_file, sigma = run_solver(phi0_mV, N, ion_size=a_nm, stretch=stretch)
        x, phi = load_potential_profile(output_file)
        results[N] = {'x': x, 'phi': phi, 'sigma': sigma, 'file': output_file}
        print(f"σ = {sigma:.4f} μC/cm²")

    # Reference solution (finest grid)
    N_ref = N_values[-1]
    x_ref = results[N_ref]['x']
    phi_ref = results[N_ref]['phi']

    print("\n" + "=" * 80)
    print("  Results")
    print("=" * 80)
    print(f"{'N':>8} {'σ [μC/cm²]':>14} {'σ err [%]':>12} {'L2 err [mV]':>14} {'σ order':>10} {'L2 order':>10}")
    print("-" * 80)

    sigma_errors = []
    l2_errors = []

    for i, N in enumerate(N_values[:-1]):  # Exclude reference grid
        sigma_num = results[N]['sigma']
        sigma_err = abs(sigma_num - sigma_analytical) / sigma_analytical * 100
        sigma_errors.append(sigma_err)

        # Interpolate reference to coarse grid for L2 comparison
        x_coarse = results[N]['x']
        phi_coarse = results[N]['phi']
        phi_ref_interp = interpolate_to_grid(x_ref, phi_ref, x_coarse)
        l2_err = compute_l2_error(phi_coarse, phi_ref_interp)
        l2_errors.append(l2_err)

        # Convergence orders
        if i > 0:
            sigma_order = np.log2(sigma_errors[-2] / sigma_errors[-1]) if sigma_errors[-1] > 1e-10 else float('nan')
            l2_order = np.log2(l2_errors[-2] / l2_errors[-1]) if l2_errors[-1] > 1e-10 else float('nan')
            print(f"{N:>8} {sigma_num:>14.6f} {sigma_err:>12.4f} {l2_err:>14.6f} {sigma_order:>10.2f} {l2_order:>10.2f}")
        else:
            print(f"{N:>8} {sigma_num:>14.6f} {sigma_err:>12.4f} {l2_err:>14.6f} {'—':>10} {'—':>10}")

    print("-" * 80)
    print(f"Reference (N={N_ref}): σ = {results[N_ref]['sigma']:.6f} μC/cm²")

    # Average convergence orders
    if len(sigma_errors) > 2:
        sigma_orders = [np.log2(sigma_errors[i] / sigma_errors[i+1])
                        for i in range(len(sigma_errors)-1) if sigma_errors[i+1] > 1e-10]
        l2_orders = [np.log2(l2_errors[i] / l2_errors[i+1])
                     for i in range(len(l2_errors)-1) if l2_errors[i+1] > 1e-10]
        print(f"\nAverage σ convergence order: {np.mean(sigma_orders):.2f}")
        print(f"Average L2 convergence order: {np.mean(l2_orders):.2f}")

    # Plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.4, 3.7))

    N_arr = np.array(N_values[:-1])

    # Left: Surface charge error
    ax1.loglog(N_arr, sigma_errors, 'o-', color='C0', linewidth=1.5, markersize=6)
    # 2nd order reference
    N_ref_line = np.array([N_arr[0], N_arr[-1]])
    err_ref_2 = sigma_errors[0] * (N_ref_line[0] / N_ref_line)**2
    ax1.loglog(N_ref_line, err_ref_2, '--', color='gray', linewidth=1, label='2nd order')
    set_labels(ax1, r'Grid points $N$', r'Surface charge error (\%)')
    ax1.legend(loc='upper right', frameon=False, fontsize=8)
    ax1.set_title('(a) Surface charge vs analytical', fontsize=10)
    # Remove auto minor locator for log scale
    ax1.minorticks_off()

    # Right: L2 error (Richardson)
    ax2.loglog(N_arr, l2_errors, 's-', color='C1', linewidth=1.5, markersize=6)
    err_ref_2_l2 = l2_errors[0] * (N_ref_line[0] / N_ref_line)**2
    ax2.loglog(N_ref_line, err_ref_2_l2, '--', color='gray', linewidth=1, label='2nd order')
    set_labels(ax2, r'Grid points $N$', r'L2 error vs $N$=6401 (mV)')
    ax2.legend(loc='upper right', frameon=False, fontsize=8)
    ax2.set_title('(b) Potential profile (Richardson)', fontsize=10)
    ax2.minorticks_off()

    plt.tight_layout()

    results_dir = Path(__file__).parent.parent / 'results'
    plt.savefig(results_dir / 'bikerman_convergence_full.png', dpi=300, bbox_inches='tight')
    plt.savefig(results_dir / 'bikerman_convergence_full.svg', bbox_inches='tight')
    print(f"\nSaved: results/bikerman_convergence_full.png")

    # Save data
    data = np.column_stack([N_arr, sigma_errors, l2_errors])
    np.savetxt(results_dir / 'bikerman_convergence_full.csv', data,
               header='N,Sigma_error_percent,L2_error_mV', delimiter=',', comments='')
    print(f"Saved: results/bikerman_convergence_full.csv")


if __name__ == '__main__':
    main()
