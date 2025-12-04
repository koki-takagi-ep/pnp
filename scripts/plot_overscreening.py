#!/usr/bin/env python3
"""
Plot overscreening/crowding results from the Modified Poisson-Fermi solver.
Reproduces figures similar to Bazant et al., PRL 106, 046102 (2011).

Usage:
    python3 scripts/plot_overscreening.py
    python3 scripts/plot_overscreening.py results/pnp_results.dat
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
import sys
import subprocess

# Add project root to path for styles import
script_dir = Path(__file__).parent
project_dir = script_dir.parent
sys.path.insert(0, str(project_dir))

from styles.plot_style import setup_plot_style, setup_axis_style, set_labels, FIGURE_SIZES


def load_data(filename):
    """Load results from the solver output file."""
    data = np.loadtxt(filename, comments='#')
    return {
        'x_nm': data[:, 0],
        'x_norm': data[:, 1],       # x / lambda_D
        'phi_mV': data[:, 2],
        'phi_norm': data[:, 3],     # e*phi / (kB*T)
        'c_plus': data[:, 4],
        'c_minus': data[:, 5],
        'c_plus_norm': data[:, 6],  # c+ / c0
        'c_minus_norm': data[:, 7], # c- / c0
        'rho': data[:, 8],          # charge density [C/m^3]
        'phi_gc_mV': data[:, 9],
    }


def run_solver(model, phi0_mV, delta_c=10, N=1001, L=50, ion_size=0.7,
               output_file="results/temp_mpf.dat"):
    """Run the PNP solver with specified parameters."""
    cmd = [
        "./build/pnp_solver",
        "--model", model,
        "--phi0", str(phi0_mV),
        "--N", str(N),
        "--L", str(L),
        "--ion-size", str(ion_size),
        "--output", output_file,
    ]
    if model in ["mpf", "modified-pf"]:
        cmd.extend(["--delta-c", str(delta_c)])

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error running solver: {result.stderr}")
        return None
    return output_file


def plot_charge_density_dimensionless(data, output_dir, title=None):
    """
    Plot dimensionless charge density vs dimensionless distance.
    Similar to Bazant 2011 Fig. 2(a).
    """
    setup_plot_style()
    fig, ax = plt.subplots(figsize=FIGURE_SIZES['single'])

    x = data['x_norm']  # x / lambda_D
    c_plus = data['c_plus_norm']
    c_minus = data['c_minus_norm']

    # Dimensionless charge density: rho_tilde = (n+ - n-) / (2*n0)
    rho_tilde = (c_plus - c_minus) / 2.0

    ax.plot(x, rho_tilde, 'b-', linewidth=1)
    ax.axhline(y=0, color='gray', linestyle=':', linewidth=0.5)

    set_labels(ax, r'$x / \lambda_D$', r'$\tilde{\rho}$')

    # Set x limit to show the near-electrode region
    x_max = min(20, x.max())
    ax.set_xlim([0, x_max])

    if title:
        ax.set_title(title, fontsize=10, fontweight='bold')

    setup_axis_style(ax)
    plt.tight_layout()
    plt.savefig(output_dir / 'charge_density_dimensionless.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("Saved: charge_density_dimensionless.png")

    return rho_tilde


def plot_potential_dimensionless(data, output_dir, title=None):
    """
    Plot dimensionless potential vs dimensionless distance.
    """
    setup_plot_style()
    fig, ax = plt.subplots(figsize=FIGURE_SIZES['single'])

    x = data['x_norm']  # x / lambda_D
    phi = data['phi_norm']  # e*phi / (kB*T)

    ax.plot(x, phi, 'b-', linewidth=1)
    ax.axhline(y=0, color='gray', linestyle=':', linewidth=0.5)

    set_labels(ax, r'$x / \lambda_D$', r'$\tilde{\phi}$')

    x_max = min(20, x.max())
    ax.set_xlim([0, x_max])

    if title:
        ax.set_title(title, fontsize=10, fontweight='bold')

    setup_axis_style(ax)
    plt.tight_layout()
    plt.savefig(output_dir / 'potential_dimensionless.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("Saved: potential_dimensionless.png")


def plot_comparison_delta_c(output_dir, phi0_mV=100, delta_c_values=[0.1, 1, 5, 10, 20]):
    """
    Compare charge density profiles for different delta_c values.
    Reproduces Bazant 2011 Fig. 2(a) style.
    """
    setup_plot_style()
    fig, ax = plt.subplots(figsize=FIGURE_SIZES['single'])

    colors = plt.cm.viridis(np.linspace(0, 0.9, len(delta_c_values)))

    for i, delta_c in enumerate(delta_c_values):
        output_file = f"results/mpf_delta_{delta_c:.1f}.dat"

        # Run solver
        run_solver("mpf", phi0_mV, delta_c=delta_c, output_file=output_file)

        # Load and plot
        if Path(output_file).exists():
            data = load_data(output_file)
            x = data['x_norm']
            c_plus = data['c_plus_norm']
            c_minus = data['c_minus_norm']
            rho_tilde = (c_plus - c_minus) / 2.0

            ax.plot(x, rho_tilde, color=colors[i], linewidth=1,
                    label=rf'$\delta_c = {delta_c}$')

    ax.axhline(y=0, color='gray', linestyle=':', linewidth=0.5)
    set_labels(ax, r'$x / \lambda_D$', r'$\tilde{\rho}$')

    ax.set_xlim([0, 10])
    ax.legend(loc='best', fontsize=8)

    setup_axis_style(ax)
    plt.tight_layout()
    plt.savefig(output_dir / 'overscreening_delta_c_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("Saved: overscreening_delta_c_comparison.png")


def plot_comparison_voltage(output_dir, delta_c=10, phi0_values=[50, 100, 200, 500, 1000]):
    """
    Compare charge density profiles for different surface potentials.
    Shows transition from overscreening to crowding.
    """
    setup_plot_style()
    fig, ax = plt.subplots(figsize=FIGURE_SIZES['single'])

    colors = plt.cm.plasma(np.linspace(0, 0.9, len(phi0_values)))

    for i, phi0 in enumerate(phi0_values):
        output_file = f"results/mpf_phi_{phi0}.dat"

        # Run solver
        run_solver("mpf", phi0, delta_c=delta_c, output_file=output_file)

        # Load and plot
        if Path(output_file).exists():
            data = load_data(output_file)
            x = data['x_norm']
            c_plus = data['c_plus_norm']
            c_minus = data['c_minus_norm']
            rho_tilde = (c_plus - c_minus) / 2.0

            ax.plot(x, rho_tilde, color=colors[i], linewidth=1,
                    label=rf'$\tilde{{V}} = {phi0/25.69:.1f}$')

    ax.axhline(y=0, color='gray', linestyle=':', linewidth=0.5)
    set_labels(ax, r'$x / \lambda_D$', r'$\tilde{\rho}$')

    ax.set_xlim([0, 15])
    ax.legend(loc='best', fontsize=8)

    setup_axis_style(ax)
    plt.tight_layout()
    plt.savefig(output_dir / 'overscreening_voltage_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("Saved: overscreening_voltage_comparison.png")


def plot_model_comparison(output_dir, phi0_mV=100, delta_c=10):
    """
    Compare Standard PB, Bikerman, and Modified PF models.
    """
    setup_plot_style()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.4, 3.7))

    models = [
        ("standard", "Standard PB", "b-"),
        ("bikerman", "Bikerman", "g--"),
        ("mpf", "Modified PF", "r-"),
    ]

    for model, label, style in models:
        output_file = f"results/model_{model}.dat"

        # Run solver
        if model == "mpf":
            run_solver(model, phi0_mV, delta_c=delta_c, output_file=output_file)
        else:
            run_solver(model, phi0_mV, output_file=output_file)

        if Path(output_file).exists():
            data = load_data(output_file)
            x = data['x_norm']
            phi = data['phi_norm']
            c_plus = data['c_plus_norm']
            c_minus = data['c_minus_norm']
            rho_tilde = (c_plus - c_minus) / 2.0

            ax1.plot(x, phi, style, linewidth=1, label=label)
            ax2.plot(x, rho_tilde, style, linewidth=1, label=label)

    # Potential plot
    ax1.axhline(y=0, color='gray', linestyle=':', linewidth=0.5)
    set_labels(ax1, r'$x / \lambda_D$', r'$\tilde{\phi}$')
    ax1.set_xlim([0, 10])
    ax1.legend(loc='best', fontsize=8)
    ax1.set_title('(a) Potential', fontsize=10, fontweight='bold')
    setup_axis_style(ax1)

    # Charge density plot
    ax2.axhline(y=0, color='gray', linestyle=':', linewidth=0.5)
    set_labels(ax2, r'$x / \lambda_D$', r'$\tilde{\rho}$')
    ax2.set_xlim([0, 10])
    ax2.legend(loc='best', fontsize=8)
    ax2.set_title('(b) Charge Density', fontsize=10, fontweight='bold')
    setup_axis_style(ax2)

    plt.tight_layout()
    plt.savefig(output_dir / 'model_comparison_mpf.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("Saved: model_comparison_mpf.png")


def detect_overscreening(data):
    """
    Detect overscreening by checking for sign changes in charge density.
    Returns the number of oscillations and their locations.
    """
    c_plus = data['c_plus_norm']
    c_minus = data['c_minus_norm']
    rho_tilde = (c_plus - c_minus) / 2.0
    x = data['x_norm']

    # Find sign changes
    sign_changes = []
    for i in range(1, len(rho_tilde)):
        if rho_tilde[i-1] * rho_tilde[i] < 0:
            # Linear interpolation for zero crossing
            x_cross = x[i-1] + (x[i] - x[i-1]) * (-rho_tilde[i-1]) / (rho_tilde[i] - rho_tilde[i-1])
            sign_changes.append(x_cross)

    return len(sign_changes), sign_changes


def main():
    results_dir = project_dir / 'results'
    results_dir.mkdir(exist_ok=True)

    # Single result file analysis
    if len(sys.argv) > 1:
        result_file = Path(sys.argv[1])
        if not result_file.exists():
            print(f"Error: Result file not found: {result_file}")
            sys.exit(1)

        print(f"Loading data from: {result_file}")
        data = load_data(result_file)

        # Check for overscreening
        n_oscillations, crossings = detect_overscreening(data)
        print(f"\nOverscreening Analysis:")
        print(f"  Number of charge density oscillations: {n_oscillations}")
        if crossings:
            print(f"  Zero-crossing positions (x/λ_D): {[f'{x:.2f}' for x in crossings]}")
        else:
            print("  No overscreening detected (monotonic charge density)")

        plot_charge_density_dimensionless(data, results_dir)
        plot_potential_dimensionless(data, results_dir)

    else:
        print("Running comprehensive Modified Poisson-Fermi analysis...")
        print("=" * 50)

        # Run comparison plots
        print("\n1. Comparing different delta_c values...")
        plot_comparison_delta_c(results_dir)

        print("\n2. Comparing different surface potentials...")
        plot_comparison_voltage(results_dir)

        print("\n3. Comparing PB, Bikerman, and MPF models...")
        plot_model_comparison(results_dir)

        # Also run a default case and analyze
        print("\n4. Analyzing default case (phi0=100mV, delta_c=10)...")
        default_file = "results/pnp_results.dat"
        run_solver("mpf", 100, delta_c=10, output_file=default_file)

        if Path(default_file).exists():
            data = load_data(default_file)
            n_oscillations, crossings = detect_overscreening(data)
            print(f"\n  Overscreening Analysis:")
            print(f"    Number of charge density oscillations: {n_oscillations}")
            if crossings:
                print(f"    Zero-crossing positions (x/λ_D): {[f'{x:.2f}' for x in crossings]}")
            else:
                print("    No overscreening detected")

            plot_charge_density_dimensionless(data, results_dir)
            plot_potential_dimensionless(data, results_dir)

    print(f"\nAll plots saved to: {results_dir}")


if __name__ == '__main__':
    main()
