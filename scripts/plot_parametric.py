#!/usr/bin/env python3
"""
Plot parametric study results from the 1D PNP solver.

This script generates comparison plots for different:
1. Surface potentials
2. Bulk concentrations
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
from pathlib import Path
import sys

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
        'x_norm': data[:, 1],
        'phi_mV': data[:, 2],
        'phi_norm': data[:, 3],
        'c_plus': data[:, 4],
        'c_minus': data[:, 5],
        'c_plus_norm': data[:, 6],
        'c_minus_norm': data[:, 7],
        'rho': data[:, 8],
        'phi_gc_mV': data[:, 9],
    }


def plot_potential_comparison(results_dir):
    """Compare potential profiles for different surface potentials."""
    setup_plot_style()
    potentials = [25, 50, 100, 200]
    colors = plt.cm.viridis(np.linspace(0.2, 0.9, len(potentials)))

    fig, ax = plt.subplots(figsize=FIGURE_SIZES['single'])

    for phi0, color in zip(potentials, colors):
        filename = results_dir / f'pnp_phi{phi0}mV.dat'
        if not filename.exists():
            print(f"Warning: {filename} not found, skipping...")
            continue

        data = load_data(filename)

        # Absolute potential
        ax.plot(data['x_norm'], data['phi_mV'], color=color,
                label=f'{phi0} mV', linewidth=1)

    set_labels(ax, r'$x / \lambda_D$', r'$\phi$ (mV)')
    ax.legend(loc='best', fontsize=8)
    ax.set_xlim([0, 10])
    setup_axis_style(ax)

    plt.tight_layout()
    plt.savefig(results_dir / 'parametric_potential.png', dpi=300, bbox_inches='tight')
    plt.savefig(results_dir / 'parametric_potential.svg', bbox_inches='tight')
    plt.close()
    print("Saved: parametric_potential.png/svg")


def plot_concentration_comparison(results_dir):
    """Compare concentration profiles for different surface potentials (non-normalized)."""
    setup_plot_style()
    potentials = [25, 50, 100, 200]
    colors = plt.cm.viridis(np.linspace(0.2, 0.9, len(potentials)))

    # Cation concentration plot
    fig, ax = plt.subplots(figsize=FIGURE_SIZES['single'])

    for phi0, color in zip(potentials, colors):
        filename = results_dir / f'pnp_phi{phi0}mV.dat'
        if not filename.exists():
            continue

        data = load_data(filename)
        ax.plot(data['x_norm'], data['c_plus'], color=color,
                label=f'{phi0} mV', linewidth=1)

    set_labels(ax, r'$x / \lambda_D$', r'$n_+$ (mol/L)')
    ax.legend(loc='best', fontsize=8)
    ax.set_xlim([0, 10])
    ax.set_yscale('log')
    setup_axis_style(ax)

    plt.tight_layout()
    plt.savefig(results_dir / 'parametric_cation.png', dpi=300, bbox_inches='tight')
    plt.savefig(results_dir / 'parametric_cation.svg', bbox_inches='tight')
    plt.close()
    print("Saved: parametric_cation.png/svg")

    # Anion concentration plot
    fig, ax = plt.subplots(figsize=FIGURE_SIZES['single'])

    for phi0, color in zip(potentials, colors):
        filename = results_dir / f'pnp_phi{phi0}mV.dat'
        if not filename.exists():
            continue

        data = load_data(filename)
        ax.plot(data['x_norm'], data['c_minus'], color=color,
                label=f'{phi0} mV', linewidth=1)

    set_labels(ax, r'$x / \lambda_D$', r'$n_-$ (mol/L)')
    ax.legend(loc='best', fontsize=8)
    ax.set_xlim([0, 10])
    ax.set_yscale('log')
    setup_axis_style(ax)

    plt.tight_layout()
    plt.savefig(results_dir / 'parametric_anion.png', dpi=300, bbox_inches='tight')
    plt.savefig(results_dir / 'parametric_anion.svg', bbox_inches='tight')
    plt.close()
    print("Saved: parametric_anion.png/svg")


def plot_bulk_concentration_effect(results_dir):
    """Compare results for different bulk concentrations."""
    setup_plot_style()
    concentrations = [0.1, 1.0, 5.0]
    labels = ['0.1 M', '1.0 M', '5.0 M']
    colors = ['blue', 'green', 'red']

    fig, ax = plt.subplots(figsize=FIGURE_SIZES['single'])

    for c0, label, color in zip(concentrations, labels, colors):
        filename = results_dir / f'pnp_c0_{c0}M.dat'
        if not filename.exists():
            print(f"Warning: {filename} not found, skipping...")
            continue

        data = load_data(filename)

        # Potential vs distance in nm
        ax.plot(data['x_nm'], data['phi_mV'], color=color,
                label=label, linewidth=1)

    set_labels(ax, r'$x$ (nm)', r'$\phi$ (mV)')
    ax.legend(loc='best', fontsize=8)
    ax.set_xlim([0, 50])
    setup_axis_style(ax)

    plt.tight_layout()
    plt.savefig(results_dir / 'parametric_bulk_concentration.png', dpi=300, bbox_inches='tight')
    plt.savefig(results_dir / 'parametric_bulk_concentration.svg', bbox_inches='tight')
    plt.close()
    print("Saved: parametric_bulk_concentration.png/svg")


def plot_debye_screening(results_dir):
    """Illustrate Debye screening effect."""
    setup_plot_style()
    concentrations = [0.1, 1.0, 5.0]
    labels = ['0.1 M', '1.0 M', '5.0 M']
    colors = ['blue', 'green', 'red']

    fig, ax = plt.subplots(figsize=FIGURE_SIZES['single'])

    for c0, label, color in zip(concentrations, labels, colors):
        filename = results_dir / f'pnp_c0_{c0}M.dat'
        if not filename.exists():
            continue

        data = load_data(filename)

        # Plot normalized potential on semilog
        mask = data['phi_mV'] > 0.1  # Avoid log(0) issues
        ax.semilogy(data['x_nm'][mask], data['phi_mV'][mask], color=color,
                    label=label, linewidth=1)

    set_labels(ax, r'$x$ (nm)', r'$\phi$ (mV)')
    ax.legend(loc='best', fontsize=8)
    ax.set_xlim([0, 50])
    setup_axis_style(ax)

    plt.tight_layout()
    plt.savefig(results_dir / 'debye_screening.png', dpi=300, bbox_inches='tight')
    plt.savefig(results_dir / 'debye_screening.svg', bbox_inches='tight')
    plt.close()
    print("Saved: debye_screening.png/svg")


def main():
    # Setup paths
    results_dir = project_dir / 'results'

    if not results_dir.exists():
        print(f"Error: Results directory not found: {results_dir}")
        print("Run 'make parametric' first.")
        sys.exit(1)

    print("Generating parametric study plots...")
    plot_potential_comparison(results_dir)
    plot_concentration_comparison(results_dir)
    plot_bulk_concentration_effect(results_dir)
    plot_debye_screening(results_dir)

    print(f"\nAll parametric plots saved to: {results_dir}")


if __name__ == '__main__':
    main()
