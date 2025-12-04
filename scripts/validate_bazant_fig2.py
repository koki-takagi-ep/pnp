#!/usr/bin/env python3
"""
Validate Modified Poisson-Fermi solver against Bazant 2011 PRL Fig.2

This script compares our MPF solver results against digitized data from
Bazant et al., PRL 106, 046102 (2011) Figure 2.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
import sys
import subprocess

script_dir = Path(__file__).parent
project_dir = script_dir.parent
sys.path.insert(0, str(project_dir))

from styles.plot_style import setup_plot_style, setup_axis_style, set_labels, FIGURE_SIZES


def load_fig2a_data():
    """Load digitized data from Bazant 2011 Fig.2a using numpy"""
    csv_path = project_dir / 'material/Bazant-PRL-2011/Poisson-Fermi-eq/ver3/Fig2a.csv'

    # Read raw data, skipping the two header rows
    raw_data = np.genfromtxt(csv_path, delimiter=',', skip_header=2)

    data = {}
    # Column pairs: (X, Y) for each case
    # Columns: V10_dc10, V100_dc10, V10_dc0, V100_dc0, V1_dc10, V1_dc0
    col_names = ['V10_dc10', 'V100_dc10', 'V10_dc0', 'V100_dc0', 'V1_dc10', 'V1_dc0']

    for i, name in enumerate(col_names):
        x_col = raw_data[:, 2*i]
        y_col = raw_data[:, 2*i + 1]
        # Remove NaN entries
        mask = ~(np.isnan(x_col) | np.isnan(y_col))
        data[name] = {
            'x': x_col[mask],
            'y': y_col[mask]
        }
    return data


def run_solver(phi0_mV, delta_c, a_nm=1.0, output_file=None):
    """Run MPF solver with Bazant 2011 parameters"""
    if output_file is None:
        output_file = f"results/bazant_V{phi0_mV}_dc{delta_c}.dat"

    # Bazant 2011 parameters (from paper)
    # gamma = 0.5 (packing fraction for symmetric electrolyte)
    # This corresponds to n_ref = 1/(2*v) at maximum packing
    # For a = 1 nm, v = 0.83 * a^3 = 0.83 nm^3

    cmd = [
        "./build/pnp_solver",
        "--model", "mpf",
        "--phi0", str(phi0_mV),
        "--delta-c", str(delta_c),
        "--ion-size", str(a_nm),  # Ion diameter in nm
        "--L", "30",              # Domain length (30 nm = 30 ion diameters for a=1nm)
        "--N", "501",             # Grid points
        "--eps", "5",             # Permittivity (typical for ionic liquids)
        "--c0", "0.5",            # Concentration to get gamma ~ 0.5
        "--output", output_file,
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Solver error: {result.stderr}")
        return None
    return output_file


def load_solver_result(filename):
    """Load solver results"""
    data = np.loadtxt(filename, comments='#')
    return {
        'x_nm': data[:, 0],
        'x_norm': data[:, 1],  # x/lambda_D
        'phi_mV': data[:, 2],
        'phi_norm': data[:, 3],
        'c_plus': data[:, 4],
        'c_minus': data[:, 5],
        'c_plus_norm': data[:, 6],
        'c_minus_norm': data[:, 7],
        'rho': data[:, 8],
    }


def convert_to_bazant_units(solver_data, a_nm, lambda_D_nm):
    """Convert solver output to Bazant's dimensionless units

    Bazant uses:
    - x/a for position (ion diameter units)
    - ρ̃ = (n+ - n-)/n_ref for charge density
    """
    x_a = solver_data['x_nm'] / a_nm  # x/a

    # Dimensionless charge density
    # rho_tilde = (c+ - c-) / c0 = c_plus_norm - c_minus_norm
    rho_tilde = solver_data['c_plus_norm'] - solver_data['c_minus_norm']

    return x_a, rho_tilde


def plot_comparison():
    """Generate comparison plot against Bazant 2011 Fig.2a"""
    setup_plot_style()
    fig, ax = plt.subplots(figsize=FIGURE_SIZES['single'])

    # Load reference data
    ref_data = load_fig2a_data()

    # Run our solver for V=100, delta_c=10 case
    print("Running solver for V=100 mV, δ_c=10...")
    output_file = run_solver(100, 10)

    if output_file and Path(output_file).exists():
        solver_data = load_solver_result(output_file)

        # Get Debye length from file header
        with open(output_file, 'r') as f:
            for line in f:
                if 'Debye length' in line:
                    lambda_D_nm = float(line.split(':')[1].strip().split()[0])
                    break

        a_nm = 0.7  # Ion diameter
        x_a, rho_tilde = convert_to_bazant_units(solver_data, a_nm, lambda_D_nm)

        # Plot our result
        ax.plot(x_a, rho_tilde, 'b-', linewidth=1.5, label='This work')

        # Plot reference data
        ax.plot(ref_data['V100_dc10']['x'], ref_data['V100_dc10']['y'],
                'ko', markersize=3, label='Bazant 2011')

    ax.axhline(y=0, color='gray', linestyle=':', linewidth=0.5)
    ax.set_xlim([0, 5])
    set_labels(ax, r'$x/a$', r'$\tilde{\rho}$')
    ax.legend(loc='best', fontsize=8)
    ax.set_title(r'$V = 100$ mV, $\delta_c = 10$', fontsize=10)

    setup_axis_style(ax)
    plt.tight_layout()

    output_dir = project_dir / 'results'
    plt.savefig(output_dir / 'bazant_fig2_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("Saved: bazant_fig2_comparison.png")


def plot_multi_comparison():
    """Compare multiple cases from Bazant 2011 Fig.2"""
    setup_plot_style()
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.7))

    ref_data = load_fig2a_data()

    # Case 1: V=100 mV, delta_c=10
    print("Running solver for V=100 mV, δ_c=10...")
    output_file = run_solver(100, 10)
    if output_file and Path(output_file).exists():
        solver_data = load_solver_result(output_file)
        a_nm = 0.7
        x_a, rho_tilde = convert_to_bazant_units(solver_data, a_nm, 0.1)

        ax = axes[0]
        ax.plot(x_a, rho_tilde, 'b-', linewidth=1.5, label='This work')
        ax.plot(ref_data['V100_dc10']['x'], ref_data['V100_dc10']['y'],
                'ko', markersize=3, label='Bazant 2011')
        ax.axhline(y=0, color='gray', linestyle=':', linewidth=0.5)
        ax.set_xlim([0, 5])
        ax.set_ylim([-1, 2.5])
        set_labels(ax, r'$x/a$', r'$\tilde{\rho}$')
        ax.legend(loc='best', fontsize=8)
        ax.set_title(r'(a) $\tilde{V} = 4$, $\delta_c = 10$', fontsize=10)
        setup_axis_style(ax)

    # Case 2: V=100 mV, delta_c=0 (Bikerman limit)
    print("Running solver for V=100 mV, δ_c=0 (Bikerman)...")
    output_file2 = "results/bazant_V100_dc0.dat"
    cmd = [
        "./build/pnp_solver",
        "--model", "bikerman",
        "--phi0", "100",
        "--ion-size", "0.7",
        "--L", "30",
        "--N", "501",
        "--eps", "5",
        "--c0", "0.5",
        "--output", output_file2,
    ]
    subprocess.run(cmd, capture_output=True, text=True)

    if Path(output_file2).exists():
        solver_data = load_solver_result(output_file2)
        a_nm = 0.7
        x_a, rho_tilde = convert_to_bazant_units(solver_data, a_nm, 0.1)

        ax = axes[1]
        ax.plot(x_a, rho_tilde, 'r-', linewidth=1.5, label='This work')
        ax.plot(ref_data['V100_dc0']['x'], ref_data['V100_dc0']['y'],
                'ko', markersize=3, label='Bazant 2011')
        ax.axhline(y=0, color='gray', linestyle=':', linewidth=0.5)
        ax.set_xlim([0, 5])
        ax.set_ylim([-0.5, 2.5])
        set_labels(ax, r'$x/a$', r'$\tilde{\rho}$')
        ax.legend(loc='best', fontsize=8)
        ax.set_title(r'(b) $\tilde{V} = 4$, $\delta_c = 0$ (Bikerman)', fontsize=10)
        setup_axis_style(ax)

    plt.tight_layout()

    output_dir = project_dir / 'results'
    plt.savefig(output_dir / 'bazant_fig2_multi_comparison.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("Saved: bazant_fig2_multi_comparison.png")


if __name__ == '__main__':
    plot_comparison()
    plot_multi_comparison()
    print("\nValidation plots saved to results/")
