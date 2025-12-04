#!/usr/bin/env python3
"""
Plot results from the 1D PNP solver.
Uses standard plot style from styles/plot_style.py
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
from matplotlib.ticker import AutoMinorLocator


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


def plot_potential(data, output_dir):
    """Plot electric potential profile with GC comparison."""
    setup_plot_style()
    fig, ax = plt.subplots(figsize=FIGURE_SIZES['single'])

    ax.plot(data['x_nm'], data['phi_mV'], 'b-', label='PNP numerical', linewidth=1)
    ax.plot(data['x_nm'], data['phi_gc_mV'], 'r--', label='Gouy-Chapman', linewidth=1)

    set_labels(ax, r'$x$ (nm)', r'$\phi$ (mV)')
    ax.legend(loc='best', fontsize=8)
    ax.set_xlim([0, min(50, data['x_nm'].max())])

    setup_axis_style(ax)
    plt.tight_layout()
    plt.savefig(output_dir / 'potential_profile.png', dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / 'potential_profile.svg', bbox_inches='tight')
    plt.close()
    print("Saved: potential_profile.png/svg")


def plot_concentrations(data, output_dir, c0_mM=1.0):
    """Plot ion concentration profiles (non-normalized, in mol/L)."""
    setup_plot_style()
    fig, ax = plt.subplots(figsize=FIGURE_SIZES['single'])

    # Convert from normalized to actual concentration (mol/L = M)
    # c_plus and c_minus are already in mol/L from solver
    ax.plot(data['x_nm'], data['c_plus'], 'r-', label=r'$n_+$ (cation)', linewidth=1)
    ax.plot(data['x_nm'], data['c_minus'], 'b-', label=r'$n_-$ (anion)', linewidth=1)

    set_labels(ax, r'$x$ (nm)', r'$n$ (mol/L)')
    ax.legend(loc='best', fontsize=8)
    ax.set_xlim([0, min(50, data['x_nm'].max())])
    ax.set_yscale('log')

    # Apply axis style without minor locator for log scale
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_tick_params(which='major', direction='in', length=6, width=1)
    ax.yaxis.set_tick_params(which='major', direction='in', length=6, width=1)
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))

    plt.tight_layout()
    plt.savefig(output_dir / 'concentration_profiles.png', dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / 'concentration_profiles.svg', bbox_inches='tight')
    plt.close()
    print("Saved: concentration_profiles.png/svg")


def plot_charge_density(data, output_dir):
    """Plot charge density profile."""
    setup_plot_style()
    fig, ax = plt.subplots(figsize=FIGURE_SIZES['single'])

    ax.plot(data['x_nm'], data['rho'] / 1e6, 'g-', linewidth=1)
    ax.axhline(y=0, color='gray', linestyle=':', linewidth=0.5)

    set_labels(ax, r'$x$ (nm)', r'$\rho$ (MC/m$^3$)')
    ax.set_xlim([0, min(50, data['x_nm'].max())])

    setup_axis_style(ax)
    plt.tight_layout()
    plt.savefig(output_dir / 'charge_density.png', dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / 'charge_density.svg', bbox_inches='tight')
    plt.close()
    print("Saved: charge_density.png/svg")


def plot_combined(data, output_dir):
    """Create a combined 2x2 plot."""
    setup_plot_style()
    fig, axes = plt.subplots(2, 2, figsize=FIGURE_SIZES['2x2'])

    # (a) Potential
    ax = axes[0, 0]
    ax.plot(data['x_nm'], data['phi_mV'], 'b-', label='PNP', linewidth=1)
    ax.plot(data['x_nm'], data['phi_gc_mV'], 'r--', label='G-C', linewidth=1)
    set_labels(ax, r'$x$ (nm)', r'$\phi$ (mV)')
    ax.set_title('(a) Electric Potential', fontsize=10, fontweight='bold')
    ax.legend(loc='best', fontsize=8)
    ax.set_xlim([0, min(50, data['x_nm'].max())])
    setup_axis_style(ax)

    # (b) Normalized potential
    ax = axes[0, 1]
    ax.plot(data['x_nm'], data['phi_norm'], 'b-', linewidth=1)
    set_labels(ax, r'$x$ (nm)', r'$e\phi / k_B T$')
    ax.set_title('(b) Normalized Potential', fontsize=10, fontweight='bold')
    ax.set_xlim([0, min(50, data['x_nm'].max())])
    setup_axis_style(ax)

    # (c) Concentrations (non-normalized)
    ax = axes[1, 0]
    ax.plot(data['x_nm'], data['c_plus'], 'r-', label=r'$n_+$', linewidth=1)
    ax.plot(data['x_nm'], data['c_minus'], 'b-', label=r'$n_-$', linewidth=1)
    set_labels(ax, r'$x$ (nm)', r'$n$ (mol/L)')
    ax.set_title('(c) Ion Concentrations', fontsize=10, fontweight='bold')
    ax.legend(loc='best', fontsize=8)
    ax.set_xlim([0, min(50, data['x_nm'].max())])
    ax.set_yscale('log')
    # Axis style without minor locator for log y-scale
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_tick_params(which='major', direction='in', length=6, width=1)
    ax.yaxis.set_tick_params(which='major', direction='in', length=6, width=1)
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))

    # (d) Charge density
    ax = axes[1, 1]
    ax.plot(data['x_nm'], data['rho'] / 1e6, 'g-', linewidth=1)
    ax.axhline(y=0, color='gray', linestyle=':', linewidth=0.5)
    set_labels(ax, r'$x$ (nm)', r'$\rho$ (MC/m$^3$)')
    ax.set_title('(d) Space Charge Density', fontsize=10, fontweight='bold')
    ax.set_xlim([0, min(50, data['x_nm'].max())])
    setup_axis_style(ax)

    plt.tight_layout()
    plt.savefig(output_dir / 'combined_results.png', dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / 'combined_results.svg', bbox_inches='tight')
    plt.close()
    print("Saved: combined_results.png/svg")


def compute_error_metrics(data):
    """Compute L2 and L-infinity error metrics."""
    phi_num = data['phi_mV']
    phi_gc = data['phi_gc_mV']

    error = phi_num - phi_gc

    # L2 error (RMS)
    L2_error = np.sqrt(np.mean(error**2))

    # Relative L2 error
    L2_ref = np.sqrt(np.mean(phi_gc**2))
    L2_rel = L2_error / L2_ref if L2_ref > 1e-10 else 0.0

    # L-infinity error (max)
    Linf_error = np.max(np.abs(error))

    return {
        'L2': L2_error,
        'L2_rel': L2_rel,
        'Linf': Linf_error
    }


def plot_error_analysis(data, output_dir):
    """Plot error analysis comparing PNP with Gouy-Chapman."""
    setup_plot_style()
    fig, ax = plt.subplots(figsize=FIGURE_SIZES['single'])

    error = np.abs(data['phi_mV'] - data['phi_gc_mV'])
    ax.plot(data['x_nm'], error, 'k-', linewidth=1)

    # Compute and display error metrics
    metrics = compute_error_metrics(data)
    textstr = (f"$L_2$ error: {metrics['L2']:.3f} mV\n"
               f"$L_2$ rel: {metrics['L2_rel']*100:.2f}%\n"
               f"$L_\\infty$ error: {metrics['Linf']:.3f} mV")
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    ax.text(0.95, 0.95, textstr, transform=ax.transAxes, fontsize=9,
            verticalalignment='top', horizontalalignment='right', bbox=props)

    set_labels(ax, r'$x$ (nm)', r'$|\phi_{PNP} - \phi_{GC}|$ (mV)')
    ax.set_xlim([0, min(50, data['x_nm'].max())])

    setup_axis_style(ax)
    plt.tight_layout()
    plt.savefig(output_dir / 'error_analysis.png', dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / 'error_analysis.svg', bbox_inches='tight')
    plt.close()
    print("Saved: error_analysis.png/svg")


def plot_grid_distribution(data, output_dir):
    """Plot grid point distribution."""
    setup_plot_style()
    fig, ax = plt.subplots(figsize=FIGURE_SIZES['single'])

    x_nm = data['x_nm']
    ax.plot(range(len(x_nm)), x_nm, 'b.-', linewidth=1, markersize=2)

    set_labels(ax, 'Grid index', r'$x$ (nm)')

    setup_axis_style(ax)
    plt.tight_layout()
    plt.savefig(output_dir / 'grid_distribution.png', dpi=300, bbox_inches='tight')
    plt.savefig(output_dir / 'grid_distribution.svg', bbox_inches='tight')
    plt.close()
    print("Saved: grid_distribution.png/svg")


def main():
    results_dir = project_dir / 'results'

    if len(sys.argv) > 1:
        result_file = Path(sys.argv[1])
    else:
        result_file = results_dir / 'pnp_results.dat'

    if not result_file.exists():
        print(f"Error: Result file not found: {result_file}")
        print("Run the solver first: make run")
        sys.exit(1)

    print(f"Loading data from: {result_file}")
    data = load_data(result_file)

    # Compute and print error metrics
    metrics = compute_error_metrics(data)
    print("\n========================================")
    print("  Error Analysis (vs. Gouy-Chapman)")
    print("========================================")
    print(f"  L2 error:          {metrics['L2']:.4f} mV")
    print(f"  Relative L2 error: {metrics['L2_rel']*100:.4f} %")
    print(f"  L-inf (max) error: {metrics['Linf']:.4f} mV")

    output_dir = results_dir
    output_dir.mkdir(exist_ok=True)

    print("\nGenerating plots...")
    plot_potential(data, output_dir)
    plot_concentrations(data, output_dir)
    plot_charge_density(data, output_dir)
    plot_combined(data, output_dir)
    plot_error_analysis(data, output_dir)
    plot_grid_distribution(data, output_dir)

    print(f"\nAll plots saved to: {output_dir}")


if __name__ == '__main__':
    main()
