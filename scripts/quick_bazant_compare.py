#!/usr/bin/env python3
"""Quick comparison with Bazant 2011 Fig.2a reference data

Bazant 2011 PRL notation:
- Ṽ = eφ₀/(kBT) : dimensionless voltage (Ṽ=10 means φ₀≈257 mV at T=298K)
- δc = lc/λD : dimensionless correlation length
- γ = a³c₀NA : packing fraction (γ=0.5 in paper)
- ρ̃ = (c- - c+)/(2γc₀) : dimensionless charge density (note: anion - cation)
- At saturation (crowding limit): |ρ̃| = 1/γ = 2 for γ=0.5

Note: Bazant defines ρ̃ with opposite sign (c- - c+) so positive voltage
gives positive ρ̃ at the electrode surface (anion accumulation).
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
import sys

script_dir = Path(__file__).parent
project_dir = script_dir.parent
sys.path.insert(0, str(project_dir))

from styles.plot_style import setup_plot_style, setup_axis_style, set_labels, FIGURE_SIZES


def load_reference_data(case='V10_dc10'):
    """Load reference data from Bazant 2011 Fig.2a

    Available cases: V10_dc10, V100_dc10, V10_dc0, V100_dc0, V1_dc10, V1_dc0
    """
    csv_path = project_dir / 'material/Bazant-PRL-2011/reference/result/Fig2a.csv'
    raw_data = np.genfromtxt(csv_path, delimiter=',', skip_header=2)

    # Column mapping
    col_map = {
        'V10_dc10': (0, 1),
        'V100_dc10': (2, 3),
        'V10_dc0': (4, 5),
        'V100_dc0': (6, 7),
        'V1_dc10': (8, 9),
        'V1_dc0': (10, 11),
    }

    xi, yi = col_map[case]
    x = raw_data[:, xi]
    y = raw_data[:, yi]
    mask = ~(np.isnan(x) | np.isnan(y))
    return x[mask], y[mask]


def load_solver_result(filename):
    """Load our solver output"""
    data = np.loadtxt(filename, comments='#')
    return {
        'x_nm': data[:, 0],
        'c_plus_norm': data[:, 6],
        'c_minus_norm': data[:, 7],
    }


def plot_single_comparison():
    """Single panel comparison for V=10, delta_c=10"""
    setup_plot_style()
    fig, ax = plt.subplots(figsize=FIGURE_SIZES['single'])

    # Load reference (V=10, δc=10)
    ref_x, ref_y = load_reference_data('V10_dc10')

    # Load our result (gamma=0.5 version)
    result_file = project_dir / 'results/bazant_V10_dc10_gamma05.dat'
    if result_file.exists():
        data = load_solver_result(result_file)

        # Convert to Bazant units
        a_nm = 0.7  # Ion diameter
        gamma = 0.5  # Packing fraction (Bazant uses γ=0.5)

        x_a = data['x_nm'] / a_nm  # x/a

        # Bazant uses ρ̃ = (c- - c+)/(2γc₀), opposite sign convention
        rho_tilde = (data['c_minus_norm'] - data['c_plus_norm']) / (2 * gamma)

        # Plot our result
        ax.plot(x_a, rho_tilde, 'b-', linewidth=1.5, label='This work (MPF)')

    # Plot reference
    ax.plot(ref_x, ref_y, 'ko', markersize=3, label='Bazant 2011')

    ax.axhline(y=0, color='gray', linestyle=':', linewidth=0.5)
    ax.axhline(y=2, color='red', linestyle=':', linewidth=0.5, alpha=0.5)
    ax.axhline(y=-2, color='red', linestyle=':', linewidth=0.5, alpha=0.5)

    ax.set_xlim([0, 5])
    ax.set_ylim([-1.5, 2.5])
    set_labels(ax, r'$x/a$', r'$\tilde{\rho}$')
    ax.legend(loc='best', fontsize=8)
    ax.set_title(r'Bazant 2011 Fig.2a: $\tilde{V}=10$, $\delta_c=10$', fontsize=10)

    setup_axis_style(ax)
    plt.tight_layout()

    output_file = project_dir / 'results/bazant_V10_comparison.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_file}")


def plot_multi_comparison():
    """Two-panel comparison: MPF (delta_c=10) vs Bikerman (delta_c=0)"""
    setup_plot_style()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.4, 3.7))

    a_nm = 0.7  # Ion diameter
    gamma = 0.5  # Packing fraction

    # Panel (a): MPF with delta_c = 10
    ref_x, ref_y = load_reference_data('V10_dc10')
    result_file = project_dir / 'results/bazant_V10_dc10_gamma05.dat'
    if result_file.exists():
        data = load_solver_result(result_file)
        x_a = data['x_nm'] / a_nm
        rho_tilde = (data['c_minus_norm'] - data['c_plus_norm']) / (2 * gamma)
        ax1.plot(x_a, rho_tilde, 'b-', linewidth=1.5, label='This work')

    ax1.plot(ref_x, ref_y, 'ko', markersize=3, label='Bazant 2011')
    ax1.axhline(y=0, color='gray', linestyle=':', linewidth=0.5)
    ax1.axhline(y=2, color='red', linestyle=':', linewidth=0.5, alpha=0.3)
    ax1.set_xlim([0, 5])
    ax1.set_ylim([-1.0, 2.5])
    set_labels(ax1, r'$x/a$', r'$\tilde{\rho}$')
    ax1.legend(loc='best', fontsize=8)
    ax1.set_title(r'(a) $\tilde{V}=10$, $\delta_c=10$', fontsize=10)
    setup_axis_style(ax1)

    # Panel (b): Bikerman (delta_c = 0)
    ref_x, ref_y = load_reference_data('V10_dc0')
    result_file = project_dir / 'results/bazant_V10_dc0.dat'
    if result_file.exists():
        data = load_solver_result(result_file)
        x_a = data['x_nm'] / a_nm
        rho_tilde = (data['c_minus_norm'] - data['c_plus_norm']) / (2 * gamma)
        ax2.plot(x_a, rho_tilde, 'r-', linewidth=1.5, label='This work')

    ax2.plot(ref_x, ref_y, 'ko', markersize=3, label='Bazant 2011')
    ax2.axhline(y=0, color='gray', linestyle=':', linewidth=0.5)
    ax2.axhline(y=2, color='red', linestyle=':', linewidth=0.5, alpha=0.3)
    ax2.set_xlim([0, 5])
    ax2.set_ylim([-1.0, 2.5])
    set_labels(ax2, r'$x/a$', r'$\tilde{\rho}$')
    ax2.legend(loc='best', fontsize=8)
    ax2.set_title(r'(b) $\tilde{V}=10$, $\delta_c=0$ (Bikerman)', fontsize=10)
    setup_axis_style(ax2)

    plt.tight_layout()
    output_file = project_dir / 'results/bazant_fig2_multi_comparison.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_file}")


def main():
    plot_single_comparison()
    plot_multi_comparison()


if __name__ == '__main__':
    main()
