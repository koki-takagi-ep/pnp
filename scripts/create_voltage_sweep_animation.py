#!/usr/bin/env python3
"""
Create GIF animation by sweeping surface voltage from 0 to target.
Uses steady-state solver at each voltage step for stability.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from pathlib import Path
import subprocess
import sys
import os

# Check if imageio is available
try:
    import imageio.v2 as imageio
    HAS_IMAGEIO = True
except ImportError:
    try:
        import imageio
        HAS_IMAGEIO = True
    except ImportError:
        HAS_IMAGEIO = False


def setup_axis_style(ax):
    """Setup axis style."""
    ax.tick_params(axis='both', which='major', direction='in', length=6, width=1, labelsize=9)
    ax.tick_params(axis='both', which='minor', direction='in', length=3, width=0.5)
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)


def load_data(filename):
    """Load results from the solver output file."""
    data = np.loadtxt(filename, comments='#')
    return {
        'x_nm': data[:, 0],
        'phi_mV': data[:, 2],
        'c_plus_norm': data[:, 6],
        'c_minus_norm': data[:, 7],
        'phi_gc_mV': data[:, 9],
    }


def create_frame(data, phi0, frame_idx, output_path, phi0_target, x_max=50):
    """Create a single frame for the animation."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # (a) Electric potential
    ax = axes[0]
    ax.plot(data['x_nm'], data['phi_mV'], 'b-', linewidth=2, label='PNP (numerical)')
    ax.plot(data['x_nm'], data['phi_gc_mV'], 'r--', linewidth=1.5, alpha=0.7, label='Gouy-Chapman')
    ax.set_xlabel(r'$x$ [nm]', fontsize=11, fontweight='bold')
    ax.set_ylabel(r'$\phi$ [mV]', fontsize=11, fontweight='bold')
    ax.set_title('Electric Potential', fontsize=11, fontweight='bold')
    ax.legend(loc='upper right', fontsize=9)
    ax.set_xlim([0, x_max])
    ax.set_ylim([0, phi0_target * 1.1])
    setup_axis_style(ax)

    # (b) Ion concentrations
    ax = axes[1]
    ax.plot(data['x_nm'], data['c_plus_norm'], 'r-', linewidth=2, label=r'$c_+/c_0$ (cation)')
    ax.plot(data['x_nm'], data['c_minus_norm'], 'b-', linewidth=2, label=r'$c_-/c_0$ (anion)')
    ax.axhline(y=1.0, color='gray', linestyle=':', linewidth=1)
    ax.set_xlabel(r'$x$ [nm]', fontsize=11, fontweight='bold')
    ax.set_ylabel(r'$c / c_0$', fontsize=11, fontweight='bold')
    ax.set_title('Ion Concentrations (EDL region)', fontsize=11, fontweight='bold')
    ax.legend(loc='upper right', fontsize=9)
    ax.set_xlim([0, min(20, x_max)])
    ax.set_yscale('log')
    ax.set_ylim([0.01, 100])
    setup_axis_style(ax)

    # Add voltage annotation
    fig.suptitle(f'Surface Potential: {phi0:.1f} mV', fontsize=14, fontweight='bold', y=0.98)

    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(output_path, dpi=100, bbox_inches='tight')
    plt.close()


def main():
    # Parameters
    phi0_target = 100  # Target surface potential [mV]
    c0 = 0.01  # Bulk concentration [mol/L]
    N = 501
    L = 50
    stretch = 3.0
    n_steps = 50  # Number of voltage steps

    if len(sys.argv) > 1:
        phi0_target = float(sys.argv[1])
    if len(sys.argv) > 2:
        n_steps = int(sys.argv[2])

    output_dir = Path('results/voltage_sweep')
    output_dir.mkdir(parents=True, exist_ok=True)

    frame_dir = output_dir / 'frames'
    frame_dir.mkdir(exist_ok=True)

    # Voltage values to sweep
    phi0_values = np.linspace(0, phi0_target, n_steps + 1)
    # Skip phi0=0 (trivial case)
    phi0_values = phi0_values[1:]

    print(f"Creating voltage sweep animation")
    print(f"  Target voltage: {phi0_target} mV")
    print(f"  Bulk concentration: {c0} mol/L")
    print(f"  Number of steps: {n_steps}")
    print()

    frame_files = []

    for idx, phi0 in enumerate(phi0_values):
        # Run steady-state solver
        result_file = output_dir / f'result_{idx:03d}.dat'
        cmd = [
            './build/pnp_solver',
            '--phi0', str(phi0),
            '--c0', str(c0),
            '--N', str(N),
            '--L', str(L),
            '--stretch', str(stretch),
            '--output', str(result_file)
        ]

        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"Error at phi0={phi0}: {result.stderr}")
            continue

        # Load results and create frame
        if result_file.exists():
            data = load_data(result_file)
            frame_path = frame_dir / f'frame_{idx:05d}.png'
            create_frame(data, phi0, idx, frame_path, phi0_target)
            frame_files.append(str(frame_path))

            if (idx + 1) % 10 == 0:
                print(f"  Completed {idx + 1}/{len(phi0_values)} (φ0 = {phi0:.1f} mV)")

    print(f"\nCreating GIF animation...")

    # Create GIF
    if HAS_IMAGEIO and frame_files:
        images = []
        for frame_file in frame_files:
            images.append(imageio.imread(frame_file))

        # Add pause at end
        for _ in range(20):
            images.append(images[-1])

        output_gif = Path('results/edl_voltage_sweep.gif')
        imageio.mimsave(str(output_gif), images, fps=10, loop=0)
        print(f"Animation saved to: {output_gif}")

        # Copy to standard location
        import shutil
        shutil.copy(output_gif, 'results/edl_evolution.gif')
        print(f"Also copied to: results/edl_evolution.gif")
    else:
        print("Error: imageio not available or no frames generated")

    print("\nDone!")


if __name__ == '__main__':
    main()
