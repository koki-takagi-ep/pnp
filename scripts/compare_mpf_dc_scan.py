#!/usr/bin/env python3
"""Compare MPF results with different δc values against Bazant 2011"""

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

# Scaling factor: Bazant's λD is √2 larger than ours
DEBYE_LENGTH_SCALE = 1.0 / np.sqrt(2)


def load_reference_data():
    """Load Bazant 2011 Fig.2a V=10, δc=10 data"""
    csv_path = project_dir / 'material/Bazant-PRL-2011/reference/result/Fig2a.csv'
    raw_data = np.genfromtxt(csv_path, delimiter=',', skip_header=2)
    x = raw_data[:, 0]
    y = raw_data[:, 1]
    mask = ~(np.isnan(x) | np.isnan(y))
    return x[mask], y[mask]


def load_solver_result(filename):
    """Load solver output"""
    data = np.loadtxt(filename, comments='#')
    return {
        'x_nm': data[:, 0],
        'c_plus_norm': data[:, 6],
        'c_minus_norm': data[:, 7],
    }


def run_mpf_solver(delta_c, output_file):
    """Run MPF solver with given δc"""
    cmd = [
        str(project_dir / 'build/pnp_solver'),
        '--model', 'mpf',
        '--phi0', '257',
        '--delta-c', str(delta_c),
        '--ion-size', '0.7',
        '--L', '10',
        '--N', '401',
        '--eps', '10',
        '--c0', '2.42',
        '--output', str(output_file)
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.returncode == 0


def main():
    setup_plot_style()
    fig, ax = plt.subplots(figsize=FIGURE_SIZES['single'])

    # Load reference
    ref_x, ref_y = load_reference_data()
    ax.plot(ref_x, ref_y, 'ko', markersize=4, label='Bazant 2011', zorder=10)

    a_nm = 0.7
    gamma = 0.5

    # Test different δc values
    dc_values = [10.0, 14.14, 7.07]
    colors = ['blue', 'red', 'green']
    labels = [r'$\delta_c = 10$', r'$\delta_c = 14.14$', r'$\delta_c = 7.07$']

    for dc, color, label in zip(dc_values, colors, labels):
        output_file = project_dir / f'results/mpf_dc{dc:.2f}.dat'

        if run_mpf_solver(dc, output_file):
            data = load_solver_result(output_file)
            x_a = data['x_nm'] / a_nm * DEBYE_LENGTH_SCALE
            rho_tilde = (data['c_minus_norm'] - data['c_plus_norm']) / (2 * gamma)
            ax.plot(x_a, rho_tilde, '-', color=color, linewidth=1.5, label=label)

    ax.axhline(y=0, color='gray', linestyle=':', linewidth=0.5)
    ax.axhline(y=2, color='gray', linestyle=':', linewidth=0.5, alpha=0.3)

    ax.set_xlim([0, 5])
    ax.set_ylim([-1.5, 2.5])
    set_labels(ax, r'$x/a$', r'$\tilde{\rho}$')
    ax.legend(loc='best', fontsize=8)
    ax.set_title(r'MPF: $\tilde{V}=10$, varying $\delta_c$', fontsize=10)

    setup_axis_style(ax)
    plt.tight_layout()

    output_file = project_dir / 'results/bazant_mpf_compare_dc.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_file}")


if __name__ == '__main__':
    main()
