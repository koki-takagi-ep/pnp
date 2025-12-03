#!/usr/bin/env python3
"""
Compare Standard PB vs Bikerman model concentration profiles.
Shows the crowding effect limiting ion concentration at the electrode surface.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent / 'styles'))
from plot_style import setup_plot_style, setup_axis_style, set_labels, FIGURE_SIZES

def load_results(filename):
    """Load solver output data."""
    data = np.loadtxt(filename, comments='#')
    return {
        'x': data[:, 0],           # nm
        'x_lambda': data[:, 1],    # x/λD
        'phi': data[:, 2],         # mV
        'c_plus': data[:, 6],      # c+/c0
        'c_minus': data[:, 7],     # c-/c0
        'rho': data[:, 8],         # C/m³
    }

def main():
    setup_plot_style()

    results_dir = Path(__file__).parent.parent / 'results'

    # Load both model results
    std = load_results(results_dir / 'standard_pb_results.dat')
    bik = load_results(results_dir / 'bikerman_results.dat')

    # Create figure with 2 panels
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.4, 3.7))

    # Left panel: Concentration profiles (zoom into left EDL)
    # Show first 2 nm only
    mask_std = std['x'] <= 2.0
    mask_bik = bik['x'] <= 2.0

    ax1.plot(std['x'][mask_std], std['c_minus'][mask_std], '-', color='C0',
             linewidth=1.5, label=r'Standard PB ($n_-$)')
    ax1.plot(bik['x'][mask_bik], bik['c_minus'][mask_bik], '--', color='C1',
             linewidth=1.5, label=r'Bikerman ($n_-$)')
    ax1.plot(std['x'][mask_std], std['c_plus'][mask_std], '-', color='C0',
             linewidth=1, alpha=0.5)
    ax1.plot(bik['x'][mask_bik], bik['c_plus'][mask_bik], '--', color='C1',
             linewidth=1, alpha=0.5)

    # Add horizontal line for Bikerman max concentration
    # max = 1/ν where ν = 2a³c₀NA, for a=0.7nm, c₀=1M: ν≈0.41, so max≈2.4
    nu = 0.413  # packing fraction from solver output
    c_max = 1.0 / nu
    ax1.axhline(y=c_max, color='gray', linestyle=':', linewidth=1,
                label=f'Max ($1/\\nu = {c_max:.1f}$)')

    set_labels(ax1, r'$x$ (nm)', r'$n/n_0$')
    setup_axis_style(ax1)
    ax1.legend(loc='upper right', frameon=False, fontsize=8)
    ax1.set_xlim(0, 2.0)
    ax1.set_ylim(0, 8)
    ax1.set_title('Left EDL', fontsize=10)

    # Right panel: Potential profiles (full domain)
    ax2.plot(std['x'], std['phi'], '-', color='C0', linewidth=1.5, label='Standard PB')
    ax2.plot(bik['x'], bik['phi'], '--', color='C1', linewidth=1.5, label='Bikerman')

    set_labels(ax2, r'$x$ (nm)', r'$\phi$ (mV)')
    setup_axis_style(ax2)
    ax2.legend(loc='center right', frameon=False, fontsize=8)
    ax2.set_xlim(0, 50)
    ax2.set_title('Potential profile', fontsize=10)

    plt.tight_layout()

    # Save
    plt.savefig(results_dir / 'model_comparison.png', dpi=300, bbox_inches='tight')
    plt.savefig(results_dir / 'model_comparison.svg', bbox_inches='tight')
    print(f"Saved: results/model_comparison.png")
    print(f"Saved: results/model_comparison.svg")

if __name__ == '__main__':
    main()
