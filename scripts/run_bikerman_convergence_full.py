#!/usr/bin/env python3
"""
Comprehensive grid convergence study for Bikerman model.
1. Surface charge error vs analytical solution
2. L2 error of potential profile vs analytical solution
"""

import subprocess
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys
import re
from scipy.integrate import solve_ivp

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


def bikerman_analytical_profile(x_nm, phi0_mV, c0_M=1.0, eps_r=12, T=298.15, a_nm=0.7, L_nm=50.0):
    """
    Analytical Bikerman potential profile by integrating the first integral.

    The first integral of Bikerman equation gives:
        |dφ/dx| = sqrt(2 kB T c0 NA / ε) * sqrt((2/ν) ln(1 - ν + ν cosh(ψ)))

    For dual-electrode (closed system) with left > right potential:
        - dφ/dx < 0 everywhere (potential decreases from left to right)
        - φ(0) = φ0/2 (left electrode, relative to bulk)
        - φ(L) = -φ0/2 (right electrode)
        - φ_bulk at center = 0

    We integrate from both electrodes towards the center.
    """
    eps = eps_0 * eps_r
    c0 = c0_M * 1000  # mol/m³
    a = a_nm * 1e-9   # m
    L = L_nm * 1e-9   # m
    phi_T = kB * T / e  # Thermal voltage

    # Packing fraction
    nu = 2 * a**3 * c0 * NA

    # Prefactor for |dφ/dx|
    prefactor = np.sqrt(2 * kB * T * c0 * NA / eps)

    def dphidx(x, phi):
        """RHS of dφ/dx equation. dφ/dx < 0 everywhere for this setup."""
        psi = phi / phi_T  # Dimensionless potential
        g = 1 - nu + nu * np.cosh(psi)
        # Avoid log of small numbers
        if g < 1e-10:
            g = 1e-10
        magnitude = prefactor * np.sqrt((2.0 / nu) * np.log(g))
        # dφ/dx is always negative (potential decreases from left to right)
        return -magnitude

    # Surface potential (relative to bulk)
    phi_surface = phi0_mV * 1e-3 / 2  # V

    # Integrate from left electrode (x=0) to center (x=L/2)
    # Initial condition: φ(0) = phi_surface > 0, φ decreases towards 0
    x_left = np.linspace(0, L/2, 20001)
    sol_left = solve_ivp(
        dphidx, [0, L/2], [phi_surface],
        t_eval=x_left, method='DOP853', rtol=1e-12, atol=1e-14
    )

    # Integrate from right electrode (x=L) towards center (x=L/2)
    # Initial condition: φ(L) = -phi_surface < 0
    # Since dφ/dx < 0 and we integrate backwards (x decreasing), φ increases towards 0
    x_right = np.linspace(L, L/2, 20001)
    sol_right = solve_ivp(
        dphidx, [L, L/2], [-phi_surface],
        t_eval=x_right, method='DOP853', rtol=1e-12, atol=1e-14
    )

    # Combine solutions
    x_full = np.concatenate([sol_left.t, sol_right.t[::-1][1:]])  # Avoid duplicate at center
    phi_full = np.concatenate([sol_left.y[0], sol_right.y[0][::-1][1:]])

    # Convert to user coordinates: add bulk potential (which is φ0/2 for symmetric case)
    # In user coords: left = φ0, right = 0, bulk = φ0/2
    phi_bulk_user = phi0_mV * 1e-3 / 2
    phi_user = phi_full + phi_bulk_user

    # Interpolate to requested x points
    x_full_nm = x_full * 1e9  # Convert to nm
    phi_interp = np.interp(x_nm, x_full_nm, phi_user * 1e3)  # Convert to mV

    return phi_interp


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

    # Grid sizes (exclude finest for L2 since we use Richardson extrapolation)
    N_values = [51, 101, 201, 401, 801, 1601, 3201]

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

    print("\n" + "=" * 80)
    print("  Results")
    print("  - Surface charge error: vs analytical solution (Bikerman)")
    print("  - L2 error: vs finest grid (Richardson extrapolation)")
    print("=" * 80)
    print(f"{'N':>8} {'σ [μC/cm²]':>14} {'σ err [%]':>12} {'L2 err [mV]':>14} {'σ order':>10} {'L2 order':>10}")
    print("-" * 80)

    sigma_errors = []
    l2_errors = []

    # Reference for L2: finest grid
    N_finest = N_values[-1]
    x_ref = results[N_finest]['x']
    phi_ref = results[N_finest]['phi']

    for i, N in enumerate(N_values[:-1]):  # Exclude finest grid
        sigma_num = results[N]['sigma']
        sigma_err = abs(sigma_num - sigma_analytical) / sigma_analytical * 100
        sigma_errors.append(sigma_err)

        # Compute L2 error vs finest grid (Richardson extrapolation)
        x_num = results[N]['x']
        phi_num = results[N]['phi']
        # Interpolate finest grid to current grid points
        phi_ref_interp = np.interp(x_num, x_ref, phi_ref)
        l2_err = compute_l2_error(phi_num, phi_ref_interp)
        l2_errors.append(l2_err)

        # Convergence orders
        if i > 0:
            sigma_order = np.log2(sigma_errors[-2] / sigma_errors[-1]) if sigma_errors[-1] > 1e-10 else float('nan')
            l2_order = np.log2(l2_errors[-2] / l2_errors[-1]) if l2_errors[-1] > 1e-10 else float('nan')
            print(f"{N:>8} {sigma_num:>14.6f} {sigma_err:>12.4f} {l2_err:>14.6f} {sigma_order:>10.2f} {l2_order:>10.2f}")
        else:
            print(f"{N:>8} {sigma_num:>14.6f} {sigma_err:>12.4f} {l2_err:>14.6f} {'—':>10} {'—':>10}")

    # Add finest grid row
    sigma_num = results[N_finest]['sigma']
    sigma_err = abs(sigma_num - sigma_analytical) / sigma_analytical * 100
    sigma_order = np.log2(sigma_errors[-1] / sigma_err) if sigma_err > 1e-10 else float('nan')
    print(f"{N_finest:>8} {sigma_num:>14.6f} {sigma_err:>12.4f} {'(reference)':>14} {sigma_order:>10.2f} {'—':>10}")
    sigma_errors.append(sigma_err)

    print("-" * 80)
    print(f"Analytical σ: {sigma_analytical:.6f} μC/cm²")

    # Average convergence orders (exclude first point)
    if len(sigma_errors) > 2:
        sigma_orders = [np.log2(sigma_errors[i] / sigma_errors[i+1])
                        for i in range(len(sigma_errors)-2) if sigma_errors[i+1] > 1e-10]
        l2_orders = [np.log2(l2_errors[i] / l2_errors[i+1])
                     for i in range(len(l2_errors)-1) if l2_errors[i+1] > 1e-10]
        print(f"\nAverage σ convergence order: {np.mean(sigma_orders):.2f}")
        print(f"Average L2 convergence order: {np.mean(l2_orders):.2f}")

    # Plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.4, 3.7))

    # For surface charge: use all N values
    N_arr_sigma = np.array(N_values)
    L = 50.0  # Domain length in nm
    dx_arr_sigma = L / (N_arr_sigma - 1)

    # For L2: exclude finest grid (it's the reference)
    N_arr_l2 = np.array(N_values[:-1])
    dx_arr_l2 = L / (N_arr_l2 - 1)

    # Left: Surface charge error vs analytical
    ax1.loglog(dx_arr_sigma, sigma_errors, 'o-', color='C0', linewidth=1.5, markersize=6)
    # 2nd order reference
    dx_ref_line = np.array([dx_arr_sigma[0], dx_arr_sigma[-1]])
    err_ref_2 = sigma_errors[-1] * (dx_ref_line / dx_ref_line[-1])**2
    ax1.loglog(dx_ref_line, err_ref_2, '--', color='gray', linewidth=1, label='2nd order')
    set_labels(ax1, r'Grid spacing $\Delta x$ (nm)', r'Surface charge error (\%)')
    ax1.legend(loc='lower right', frameon=False, fontsize=8)
    ax1.set_title('(a) vs analytical $\\sigma$', fontsize=10)
    ax1.minorticks_off()

    # Right: L2 error vs finest grid (Richardson)
    ax2.loglog(dx_arr_l2, l2_errors, 's-', color='C1', linewidth=1.5, markersize=6)
    dx_ref_line_l2 = np.array([dx_arr_l2[0], dx_arr_l2[-1]])
    err_ref_2_l2 = l2_errors[-1] * (dx_ref_line_l2 / dx_ref_line_l2[-1])**2
    ax2.loglog(dx_ref_line_l2, err_ref_2_l2, '--', color='gray', linewidth=1, label='2nd order')
    set_labels(ax2, r'Grid spacing $\Delta x$ (nm)', r'L2 error vs finest (mV)')
    ax2.legend(loc='lower right', frameon=False, fontsize=8)
    ax2.set_title('(b) Richardson extrapolation', fontsize=10)
    ax2.minorticks_off()

    plt.tight_layout()

    results_dir = Path(__file__).parent.parent / 'results'
    plt.savefig(results_dir / 'bikerman_convergence_full.png', dpi=300, bbox_inches='tight')
    plt.savefig(results_dir / 'bikerman_convergence_full.svg', bbox_inches='tight')
    print(f"\nSaved: results/bikerman_convergence_full.png")

    # Save data (pad l2_errors with nan for finest grid)
    l2_errors_padded = l2_errors + [np.nan]
    data = np.column_stack([N_arr_sigma, sigma_errors, l2_errors_padded])
    np.savetxt(results_dir / 'bikerman_convergence_full.csv', data,
               header='N,Sigma_error_percent,L2_error_vs_finest_mV', delimiter=',', comments='')
    print(f"Saved: results/bikerman_convergence_full.csv")


if __name__ == '__main__':
    main()
