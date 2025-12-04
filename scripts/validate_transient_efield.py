#!/usr/bin/env python3
"""
Validation script for E-field transient solver.

Benchmark 1: Convergence to steady-state
- Run transient simulation until t_final
- Compare final state with steady-state solver result
- Report L2 error and max error

Benchmark 2: Charge conservation
- Track total charge throughout simulation
- Verify Q_err stays bounded
"""

import subprocess
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

# Add styles to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'styles'))
from plot_style import setup_plot_style, setup_axis_style, set_labels, FIGURE_SIZES

def run_solver(args):
    """Run the PNP solver with given arguments."""
    cmd = ['./build/pnp_solver'] + args
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result

def load_results(filename):
    """Load results from .dat file."""
    # Count header lines (lines starting with #)
    with open(filename, 'r') as f:
        header_lines = 0
        for line in f:
            if line.startswith('#'):
                header_lines += 1
            else:
                break

    data = np.loadtxt(filename, skiprows=header_lines)
    # Column format:
    # 0: x [nm], 1: x/lambda_D, 2: phi [mV], 3: phi_norm,
    # 4: c+ [mol/m³], 5: c- [mol/m³], 6: c+/c0, 7: c-/c0, ...
    return {
        'x': data[:, 0],      # nm
        'phi': data[:, 2],    # mV
        'c_plus': data[:, 6], # normalized c+/c0
        'c_minus': data[:, 7] # normalized c-/c0
    }

def compute_errors(transient, steady):
    """Compute L2 and max errors between transient and steady state."""
    # Interpolate steady to transient grid if needed
    if len(transient['x']) != len(steady['x']):
        from scipy.interpolate import interp1d
        for key in ['phi', 'c_plus', 'c_minus']:
            f = interp1d(steady['x'], steady[key], kind='cubic', fill_value='extrapolate')
            steady[key] = f(transient['x'])
        steady['x'] = transient['x'].copy()

    errors = {}
    for key in ['phi', 'c_plus', 'c_minus']:
        diff = transient[key] - steady[key]
        dx = np.diff(transient['x'])
        dx = np.append(dx, dx[-1])  # Extend to same length

        # L2 error
        l2 = np.sqrt(np.sum(diff**2 * dx) / np.sum(dx))
        # Max error
        linf = np.max(np.abs(diff))
        # Relative L2
        ref_l2 = np.sqrt(np.sum(steady[key]**2 * dx) / np.sum(dx))
        rel_l2 = l2 / ref_l2 if ref_l2 > 1e-10 else l2

        errors[key] = {'L2': l2, 'Linf': linf, 'rel_L2': rel_l2}

    return errors

def main():
    setup_plot_style()

    # Parameters for validation
    phi0 = 100  # mV
    N = 101     # Grid points
    t_final = 10.0  # µs (EDL charging time ~ L²/D ~ 25 µs)

    print("=" * 60)
    print("  E-field Transient Solver Validation")
    print("=" * 60)
    print(f"\nParameters:")
    print(f"  phi0 = {phi0} mV")
    print(f"  N = {N} grid points")
    print(f"  t_final = {t_final} µs")

    # Create output directory
    os.makedirs('results/validation', exist_ok=True)

    # ========================================
    # Benchmark 1: Steady-state solver
    # ========================================
    print("\n" + "-" * 40)
    print("Running steady-state solver...")
    run_solver([
        f'--phi0', str(phi0),
        '--phi-right', '0',
        '--closed-system',
        f'--N', str(N),
        '--output', 'results/validation/steady_state.dat'
    ])
    steady = load_results('results/validation/steady_state.dat')
    print(f"  phi range: [{steady['phi'].min():.2f}, {steady['phi'].max():.2f}] mV")
    print(f"  c+/c0 max: {steady['c_plus'].max():.4f}")

    # ========================================
    # Benchmark 2: Transient solver
    # ========================================
    print("\n" + "-" * 40)
    print("Running transient solver (E-field formulation)...")
    result = run_solver([
        f'--phi0', str(phi0),
        '--phi-right', '0',
        '--closed-system',
        f'--N', str(N),
        '--efield',
        '--dt', '0.1',
        f'--t-final', str(t_final),
        '--output', 'results/validation/transient_final.dat'
    ])

    # Parse output for diagnostics
    # Format: "  t = 2.51 ns, iter = 32, c+/c0 = 1.04e+00, c-/c0 = 1.04e+00, Q_err = 2.34e-13"
    lines = result.stdout.split('\n')
    q_errors = []
    times = []
    c_max = []
    for line in lines:
        if 't = ' in line and 'Q_err' in line:
            parts = line.split(',')
            # parts[0]: "  t = X ns", parts[1]: " iter = N"
            # parts[2]: " c+/c0 = X", parts[3]: " c-/c0 = X", parts[4]: " Q_err = X"
            t = float(parts[0].split('=')[1].strip().replace(' ns', ''))
            c = float(parts[2].split('=')[1].strip())  # c+/c0
            q = float(parts[4].split('=')[1].strip())  # Q_err (was parts[3], should be parts[4])
            times.append(t)
            q_errors.append(q)
            c_max.append(c)

    transient = load_results('results/validation/transient_final.dat')
    print(f"  phi range: [{transient['phi'].min():.2f}, {transient['phi'].max():.2f}] mV")
    print(f"  c+/c0 max: {transient['c_plus'].max():.4f}")

    # ========================================
    # Compute errors
    # ========================================
    print("\n" + "-" * 40)
    print("Error Analysis (Transient vs Steady-State):")
    errors = compute_errors(transient, steady)

    print(f"\n  Potential φ:")
    print(f"    L2 error:     {errors['phi']['L2']:.4e} mV")
    print(f"    Max error:    {errors['phi']['Linf']:.4e} mV")
    print(f"    Relative L2:  {errors['phi']['rel_L2']:.4e}")

    print(f"\n  Cation c+/c0:")
    print(f"    L2 error:     {errors['c_plus']['L2']:.4e}")
    print(f"    Max error:    {errors['c_plus']['Linf']:.4e}")
    print(f"    Relative L2:  {errors['c_plus']['rel_L2']:.4e}")

    print(f"\n  Anion c-/c0:")
    print(f"    L2 error:     {errors['c_minus']['L2']:.4e}")
    print(f"    Max error:    {errors['c_minus']['Linf']:.4e}")
    print(f"    Relative L2:  {errors['c_minus']['rel_L2']:.4e}")

    # ========================================
    # Charge conservation analysis
    # ========================================
    print("\n" + "-" * 40)
    print("Charge Conservation:")
    if len(q_errors) > 0:
        print(f"    Initial Q_err: {q_errors[0]:.4e} C/m²")
        print(f"    Final Q_err:   {q_errors[-1]:.4e} C/m²")
        print(f"    Max Q_err:     {max(q_errors):.4e} C/m²")

    # ========================================
    # Create validation plots
    # ========================================
    fig, axes = plt.subplots(2, 2, figsize=(7.4, 7.4))

    # (a) Potential comparison
    ax = axes[0, 0]
    ax.plot(steady['x'], steady['phi'], 'k-', lw=1.5, label='Steady-state')
    ax.plot(transient['x'], transient['phi'], 'r--', lw=1, label='Transient (final)')
    ax.set_xlim([0, 5])
    set_labels(ax, r'$x$ (nm)', r'$\phi$ (mV)')
    ax.legend(fontsize=8)
    ax.set_title('(a) Potential profile', fontsize=10)
    setup_axis_style(ax)

    # (b) Concentration comparison
    ax = axes[0, 1]
    ax.plot(steady['x'], steady['c_plus'], 'b-', lw=1.5, label=r'$n_+$ steady')
    ax.plot(steady['x'], steady['c_minus'], 'r-', lw=1.5, label=r'$n_-$ steady')
    ax.plot(transient['x'], transient['c_plus'], 'b--', lw=1, label=r'$n_+$ transient')
    ax.plot(transient['x'], transient['c_minus'], 'r--', lw=1, label=r'$n_-$ transient')
    ax.set_xlim([0, 5])
    set_labels(ax, r'$x$ (nm)', r'$n/n_0$')
    ax.legend(fontsize=8, loc='upper right')
    ax.set_title('(b) Concentration profiles', fontsize=10)
    setup_axis_style(ax)

    # (c) Charge conservation error
    ax = axes[1, 0]
    if len(times) > 0:
        ax.semilogy(np.array(times)/1000, q_errors, 'k-', lw=1)
        set_labels(ax, r'$t$ (µs)', r'$|Q_{\mathrm{err}}|$ (C/m²)')
        ax.set_title('(c) Charge conservation error', fontsize=10)
    setup_axis_style(ax)

    # (d) Concentration evolution
    ax = axes[1, 1]
    if len(times) > 0:
        ax.plot(np.array(times)/1000, c_max, 'b-', lw=1)
        # Add steady-state reference
        ax.axhline(steady['c_plus'].max(), color='k', ls='--', lw=0.5, label='Steady-state')
        set_labels(ax, r'$t$ (µs)', r'$\max(n_+/n_0)$')
        ax.set_title('(d) Max concentration evolution', fontsize=10)
        ax.legend(fontsize=8)
    setup_axis_style(ax)

    plt.tight_layout()
    plt.savefig('results/validation/transient_validation.png', dpi=150, bbox_inches='tight')
    plt.savefig('results/validation/transient_validation.svg', bbox_inches='tight')
    print("\nPlots saved to results/validation/transient_validation.png")

    # ========================================
    # Summary
    # ========================================
    print("\n" + "=" * 60)
    print("  Validation Summary")
    print("=" * 60)

    # Check if validation passed
    phi_pass = errors['phi']['rel_L2'] < 0.01  # 1% relative error
    c_pass = errors['c_plus']['rel_L2'] < 0.01 and errors['c_minus']['rel_L2'] < 0.01
    q_pass = max(q_errors) < 1e-9 if len(q_errors) > 0 else False

    print(f"\n  [{'PASS' if phi_pass else 'FAIL'}] Potential convergence: rel_L2 = {errors['phi']['rel_L2']:.2e} (< 1%)")
    print(f"  [{'PASS' if c_pass else 'FAIL'}] Concentration convergence: rel_L2 = {max(errors['c_plus']['rel_L2'], errors['c_minus']['rel_L2']):.2e} (< 1%)")
    print(f"  [{'PASS' if q_pass else 'FAIL'}] Charge conservation: max_err = {max(q_errors):.2e} C/m² (< 1e-9)")

    all_pass = phi_pass and c_pass and q_pass
    print(f"\n  Overall: {'ALL TESTS PASSED' if all_pass else 'SOME TESTS FAILED'}")

    return 0 if all_pass else 1

if __name__ == '__main__':
    sys.exit(main())
