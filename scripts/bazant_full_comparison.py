#!/usr/bin/env python3
"""Full comparison with Bazant 2011 Fig.2a reference data

Compares Ṽ = 10 and Ṽ = 100 cases for both MPF (δc=10) and Bikerman (δc=0).
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

# Scaling factor: Bazant's λD is √2 larger than ours
DEBYE_LENGTH_SCALE = 1.0 / np.sqrt(2)


def load_reference_data(case='V10_dc10'):
    """Load reference data from Bazant 2011 Fig.2a

    Available cases: V10_dc10, V100_dc10, V10_dc0, V100_dc0, V1_dc10, V1_dc0
    """
    csv_path = project_dir / 'material/Bazant-PRL-2011/reference/result/Fig2a.csv'
    raw_data = np.genfromtxt(csv_path, delimiter=',', skip_header=2)

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


def plot_full_comparison():
    """Four-panel comparison: V=10 and V=100, MPF and Bikerman"""
    setup_plot_style()
    fig, axes = plt.subplots(2, 2, figsize=(7.4, 7.4))

    a_nm = 0.7  # Ion diameter
    gamma = 0.5  # Packing fraction

    # File mapping: (case, result_file, color, model_label)
    panels = [
        # Row 1: V=10
        (axes[0, 0], 'V10_dc10', 'results/bazant_mpf_dc10.dat', 'blue',
         r'(a) $\tilde{V}=10$, $\delta_c=10$ (MPF)', [-1.5, 2.5]),
        (axes[0, 1], 'V10_dc0', 'results/bazant_bikerman.dat', 'red',
         r'(b) $\tilde{V}=10$, $\delta_c=0$ (Bikerman)', [-0.5, 2.5]),
        # Row 2: V=100
        (axes[1, 0], 'V100_dc10', 'results/bazant_mpf_V100_dc10.dat', 'blue',
         r'(c) $\tilde{V}=100$, $\delta_c=10$ (MPF)', [-1.5, 2.5]),
        (axes[1, 1], 'V100_dc0', 'results/bazant_bikerman_V100.dat', 'red',
         r'(d) $\tilde{V}=100$, $\delta_c=0$ (Bikerman)', [-0.5, 2.5]),
    ]

    for ax, ref_case, result_path, color, title, ylim in panels:
        # Load reference
        ref_x, ref_y = load_reference_data(ref_case)

        # Load our result
        result_file = project_dir / result_path
        if result_file.exists():
            data = load_solver_result(result_file)
            x_a = data['x_nm'] / a_nm * DEBYE_LENGTH_SCALE
            rho_tilde = (data['c_minus_norm'] - data['c_plus_norm']) / (2 * gamma)
            ax.plot(x_a, rho_tilde, '-', color=color, linewidth=1.5, label='This work')
        else:
            print(f"Warning: {result_file} not found")

        ax.plot(ref_x, ref_y, 'ko', markersize=3, label='Bazant 2011')
        ax.axhline(y=0, color='gray', linestyle=':', linewidth=0.5)
        ax.axhline(y=2, color='gray', linestyle=':', linewidth=0.5, alpha=0.3)
        ax.set_xlim([0, 5])
        ax.set_ylim(ylim)
        set_labels(ax, r'$x/a$', r'$\tilde{\rho}$')
        ax.legend(loc='best', fontsize=8)
        ax.set_title(title, fontsize=10)
        setup_axis_style(ax)

    plt.tight_layout()
    output_file = project_dir / 'results/bazant_full_comparison.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_file}")


def main():
    plot_full_comparison()


if __name__ == '__main__':
    main()
