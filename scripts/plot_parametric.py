#!/usr/bin/env python3
"""
Plot parametric study results from the 1D PNP solver.

This script generates comparison plots for different:
1. Surface potentials
2. Bulk concentrations
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys

# Configure matplotlib for publication-quality figures
plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 12,
    'axes.labelsize': 14,
    'axes.titlesize': 14,
    'legend.fontsize': 10,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'lines.linewidth': 2,
    'lines.markersize': 6,
    'figure.figsize': (8, 6),
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.grid': True,
    'grid.alpha': 0.3,
    'axes.axisbelow': True,
})


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
    potentials = [25, 50, 100, 200]
    colors = plt.cm.viridis(np.linspace(0.2, 0.9, len(potentials)))

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    for phi0, color in zip(potentials, colors):
        filename = results_dir / f'pnp_phi{phi0}mV.dat'
        if not filename.exists():
            print(f"Warning: {filename} not found, skipping...")
            continue

        data = load_data(filename)

        # Absolute potential
        axes[0].plot(data['x_norm'], data['phi_mV'], color=color,
                     label=f'{phi0} mV', linewidth=2)

        # Normalized potential
        axes[1].plot(data['x_norm'], data['phi_mV'] / phi0, color=color,
                     label=f'{phi0} mV', linewidth=2)

    axes[0].set_xlabel(r'$x / \lambda_D$')
    axes[0].set_ylabel(r'$\phi$ [mV]')
    axes[0].set_title('(a) Electric Potential')
    axes[0].legend(title='Surface potential')
    axes[0].set_xlim([0, 10])

    axes[1].set_xlabel(r'$x / \lambda_D$')
    axes[1].set_ylabel(r'$\phi / \phi_0$')
    axes[1].set_title('(b) Normalized Potential')
    axes[1].legend(title='Surface potential')
    axes[1].set_xlim([0, 10])

    plt.tight_layout()
    plt.savefig(results_dir / 'parametric_potential.png')
    plt.savefig(results_dir / 'parametric_potential.pdf')
    plt.close()
    print("Saved: parametric_potential.png/pdf")


def plot_concentration_comparison(results_dir):
    """Compare concentration profiles for different surface potentials."""
    potentials = [25, 50, 100, 200]
    colors = plt.cm.viridis(np.linspace(0.2, 0.9, len(potentials)))

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    for phi0, color in zip(potentials, colors):
        filename = results_dir / f'pnp_phi{phi0}mV.dat'
        if not filename.exists():
            continue

        data = load_data(filename)

        # Cation concentration
        axes[0].plot(data['x_norm'], data['c_plus_norm'], color=color,
                     label=f'{phi0} mV', linewidth=2)

        # Anion concentration
        axes[1].plot(data['x_norm'], data['c_minus_norm'], color=color,
                     label=f'{phi0} mV', linewidth=2)

    for ax in axes:
        ax.axhline(y=1.0, color='gray', linestyle=':', linewidth=1)
        ax.set_xlabel(r'$x / \lambda_D$')
        ax.set_xlim([0, 10])
        ax.set_yscale('log')

    axes[0].set_ylabel(r'$c_+ / c_0$')
    axes[0].set_title('(a) Cation Concentration')
    axes[0].legend(title='Surface potential')

    axes[1].set_ylabel(r'$c_- / c_0$')
    axes[1].set_title('(b) Anion Concentration')
    axes[1].legend(title='Surface potential')

    plt.tight_layout()
    plt.savefig(results_dir / 'parametric_concentration.png')
    plt.savefig(results_dir / 'parametric_concentration.pdf')
    plt.close()
    print("Saved: parametric_concentration.png/pdf")


def plot_bulk_concentration_effect(results_dir):
    """Compare results for different bulk concentrations."""
    concentrations = [0.1, 1.0, 5.0]
    labels = ['0.1 M', '1.0 M', '5.0 M']
    colors = ['blue', 'green', 'red']

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    for c0, label, color in zip(concentrations, labels, colors):
        filename = results_dir / f'pnp_c0_{c0}M.dat'
        if not filename.exists():
            print(f"Warning: {filename} not found, skipping...")
            continue

        data = load_data(filename)

        # Potential vs distance in nm
        axes[0].plot(data['x_nm'], data['phi_mV'], color=color,
                     label=label, linewidth=2)

        # Potential vs normalized distance
        axes[1].plot(data['x_norm'], data['phi_mV'], color=color,
                     label=label, linewidth=2)

    axes[0].set_xlabel(r'$x$ [nm]')
    axes[0].set_ylabel(r'$\phi$ [mV]')
    axes[0].set_title('(a) Potential vs. Distance')
    axes[0].legend(title='Concentration')
    axes[0].set_xlim([0, 50])

    axes[1].set_xlabel(r'$x / \lambda_D$')
    axes[1].set_ylabel(r'$\phi$ [mV]')
    axes[1].set_title('(b) Potential vs. Normalized Distance')
    axes[1].legend(title='Concentration')
    axes[1].set_xlim([0, 10])

    plt.tight_layout()
    plt.savefig(results_dir / 'parametric_bulk_concentration.png')
    plt.savefig(results_dir / 'parametric_bulk_concentration.pdf')
    plt.close()
    print("Saved: parametric_bulk_concentration.png/pdf")


def plot_debye_screening(results_dir):
    """Illustrate Debye screening effect."""
    concentrations = [0.1, 1.0, 5.0]
    labels = ['0.1 M', '1.0 M', '5.0 M']
    colors = ['blue', 'green', 'red']

    fig, ax = plt.subplots(figsize=(8, 6))

    for c0, label, color in zip(concentrations, labels, colors):
        filename = results_dir / f'pnp_c0_{c0}M.dat'
        if not filename.exists():
            continue

        data = load_data(filename)

        # Plot normalized potential on semilog
        mask = data['phi_mV'] > 0.1  # Avoid log(0) issues
        ax.semilogy(data['x_nm'][mask], data['phi_mV'][mask], color=color,
                    label=label, linewidth=2)

    ax.set_xlabel(r'$x$ [nm]')
    ax.set_ylabel(r'$\phi$ [mV]')
    ax.set_title('Debye Screening Effect')
    ax.legend(title='Concentration')
    ax.set_xlim([0, 50])

    plt.tight_layout()
    plt.savefig(results_dir / 'debye_screening.png')
    plt.savefig(results_dir / 'debye_screening.pdf')
    plt.close()
    print("Saved: debye_screening.png/pdf")


def main():
    # Setup paths
    script_dir = Path(__file__).parent
    project_dir = script_dir.parent
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
