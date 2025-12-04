#!/usr/bin/env python3
"""
Validate Modified Poisson-Fermi solver against Bazant 2011 PRL Fig.2

This script compares our MPF solver results against digitized data from
Bazant et al., PRL 106, 046102 (2011) Figure 2.

Key dimensionless parameters from Bazant 2011:
- Ṽ = eφ₀/(kBT) : dimensionless voltage (Ṽ=100 corresponds to φ₀≈2570 mV)
- δc = lc/λD : dimensionless correlation length
- γ = 2a³c₀NA : packing fraction (γ=0.5 in paper)
- ρ̃ = (n+ - n-)/(2γnmax) : dimensionless charge density, where nmax = 1/a³

At saturation (crowding limit): ρ̃ = ±1/γ = ±2 for γ=0.5
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

# Physical constants
kB = 1.380649e-23  # J/K
e = 1.602176634e-19  # C
NA = 6.02214076e23  # 1/mol
T = 298.15  # K
phi_T = kB * T / e  # Thermal voltage ≈ 25.7 mV


def load_fig2a_data():
    """Load digitized data from Bazant 2011 Fig.2a"""
    # Use the correct path to reference data
    csv_path = project_dir / 'material/Bazant-PRL-2011/reference/result/Fig2a.csv'

    if not csv_path.exists():
        # Fallback to old path
        csv_path = project_dir / 'material/Bazant-PRL-2011/Poisson-Fermi-eq/ver3/Fig2a.csv'

    # Read raw data, skipping the two header rows
    raw_data = np.genfromtxt(csv_path, delimiter=',', skip_header=2)

    data = {}
    # Column pairs: (X, Y) for each case
    # V10_deltac10, V100_deltac10, V10_deltac0, V100_deltac0, V1_deltac10, V1_deltac0
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


def run_mpf_solver(V_tilde, delta_c, gamma=0.5, output_file=None):
    """
    Run MPF solver with Bazant 2011 parameters

    Parameters:
    -----------
    V_tilde : float
        Dimensionless voltage eφ₀/(kBT)
    delta_c : float
        Dimensionless correlation length lc/λD
    gamma : float
        Packing fraction (default 0.5)
    """
    # Convert dimensionless voltage to mV
    phi0_mV = V_tilde * phi_T * 1000  # Convert V to mV

    if output_file is None:
        output_file = f"results/bazant_V{V_tilde}_dc{delta_c}.dat"

    # For γ = 0.5 and a = 0.7 nm:
    # γ = 2a³c₀NA → c₀ = γ/(2a³NA) = 0.5/(2×(0.7e-9)³×6.02e23) ≈ 1.21 M
    a_nm = 0.7  # Ion diameter in nm
    a_m = a_nm * 1e-9
    c0_M = gamma / (2 * a_m**3 * NA) / 1000  # Convert to mol/L

    # Debye length: λD = sqrt(εε₀kBT/(2e²c₀NA))
    # For ionic liquids, ε ≈ 10-15
    eps_r = 10.0

    # Domain length should be several times the expected EDL thickness
    # For high voltage, crowding layer ~ a, oscillation wavelength ~ 2π√(2δc) × λD
    L_nm = max(20, 10 * delta_c)  # nm

    cmd = [
        "./build/pnp_solver",
        "--model", "mpf",
        "--phi0", str(phi0_mV),
        "--delta-c", str(delta_c),
        "--ion-size", str(a_nm),
        "--L", str(L_nm),
        "--N", "2001",  # High resolution
        "--eps", str(eps_r),
        "--c0", str(c0_M),
        "--output", output_file,
    ]

    print(f"  Running: Ṽ={V_tilde}, δc={delta_c}, φ₀={phi0_mV:.1f} mV, c₀={c0_M:.3f} M")

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  Solver error: {result.stderr}")
        return None
    return output_file


def run_bikerman_solver(V_tilde, gamma=0.5, output_file=None):
    """Run Bikerman solver (δc = 0 limit)"""
    phi0_mV = V_tilde * phi_T * 1000

    if output_file is None:
        output_file = f"results/bazant_V{V_tilde}_dc0.dat"

    a_nm = 0.7
    a_m = a_nm * 1e-9
    c0_M = gamma / (2 * a_m**3 * NA) / 1000
    eps_r = 10.0
    L_nm = 20

    cmd = [
        "./build/pnp_solver",
        "--model", "bikerman",
        "--phi0", str(phi0_mV),
        "--ion-size", str(a_nm),
        "--L", str(L_nm),
        "--N", "2001",
        "--eps", str(eps_r),
        "--c0", str(c0_M),
        "--output", output_file,
    ]

    print(f"  Running Bikerman: Ṽ={V_tilde}, φ₀={phi0_mV:.1f} mV, c₀={c0_M:.3f} M")

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  Solver error: {result.stderr}")
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


def get_lambda_D(filename):
    """Extract Debye length from solver output header"""
    with open(filename, 'r') as f:
        for line in f:
            if 'Debye length' in line:
                return float(line.split(':')[1].strip().split()[0])
    return None


def convert_to_bazant_units(solver_data, a_nm, gamma=0.5):
    """
    Convert solver output to Bazant's dimensionless units

    Bazant uses:
    - x/a for position (ion diameter units)
    - ρ̃ = (c+ - c-)/(2γc₀) for charge density

    At saturation: c- = cmax = 1/(a³NA) → ρ̃ = -cmax/(2γc₀) = -1/γ = -2
    """
    x_a = solver_data['x_nm'] / a_nm  # x/a

    # Dimensionless charge density: ρ̃ = (c+ - c-)/(2γc₀)
    # Since c_plus_norm = c+/c₀ and c_minus_norm = c-/c₀:
    # ρ̃ = (c_plus_norm - c_minus_norm) / (2γ)
    rho_tilde = (solver_data['c_plus_norm'] - solver_data['c_minus_norm']) / (2 * gamma)

    return x_a, rho_tilde


def plot_fig2a_comparison():
    """
    Generate comprehensive comparison plot against Bazant 2011 Fig.2a

    Tests multiple cases:
    - Ṽ = 1, 10, 100
    - δc = 0 (Bikerman) and δc = 10 (MPF with overscreening)
    """
    setup_plot_style()
    fig, ax = plt.subplots(figsize=FIGURE_SIZES['single'])

    # Load reference data
    ref_data = load_fig2a_data()
    gamma = 0.5
    a_nm = 0.7

    # Define line styles matching Bazant's figure
    # δc = 10: solid lines, δc = 0: dashed lines
    cases = [
        ('V100_dc10', 100, 10, 'k-', 'Ṽ=100, δc=10'),
        ('V10_dc10', 10, 10, 'b-', 'Ṽ=10, δc=10'),
        ('V1_dc10', 1, 10, 'g-', 'Ṽ=1, δc=10'),
        ('V100_dc0', 100, 0, 'k--', 'Ṽ=100, δc=0'),
        ('V10_dc0', 10, 0, 'b--', 'Ṽ=10, δc=0'),
        ('V1_dc0', 1, 0, 'g--', 'Ṽ=1, δc=0'),
    ]

    print("Running MPF solver for Bazant Fig.2a comparison...")

    for ref_key, V_tilde, delta_c, style, label in cases:
        # Run solver
        if delta_c > 0:
            output_file = run_mpf_solver(V_tilde, delta_c, gamma)
        else:
            output_file = run_bikerman_solver(V_tilde, gamma)

        if output_file and Path(output_file).exists():
            solver_data = load_solver_result(output_file)
            x_a, rho_tilde = convert_to_bazant_units(solver_data, a_nm, gamma)

            # Plot our result
            ax.plot(x_a, rho_tilde, style, linewidth=1.5, label=f'{label} (this work)')

        # Plot reference data (small markers)
        if ref_key in ref_data:
            ax.plot(ref_data[ref_key]['x'], ref_data[ref_key]['y'],
                    'o', markersize=2, color='gray', alpha=0.5)

    ax.axhline(y=0, color='gray', linestyle=':', linewidth=0.5)
    ax.axhline(y=2, color='red', linestyle=':', linewidth=0.5, alpha=0.5)  # Saturation
    ax.axhline(y=-2, color='red', linestyle=':', linewidth=0.5, alpha=0.5)

    ax.set_xlim([0, 5])
    ax.set_ylim([-1.5, 2.5])
    set_labels(ax, r'$x/a$', r'$\tilde{\rho}$')
    ax.legend(loc='lower right', fontsize=6, ncol=2)
    ax.set_title('Comparison with Bazant 2011 PRL Fig.2a', fontsize=10)

    setup_axis_style(ax)
    plt.tight_layout()

    output_dir = project_dir / 'results'
    plt.savefig(output_dir / 'bazant_fig2a_validation.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("Saved: bazant_fig2a_validation.png")


def plot_single_case_comparison(V_tilde=100, delta_c=10):
    """
    Detailed comparison for a single case (Ṽ=100, δc=10)
    """
    setup_plot_style()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.4, 3.7))

    ref_data = load_fig2a_data()
    gamma = 0.5
    a_nm = 0.7

    print(f"\nRunning detailed comparison for Ṽ={V_tilde}, δc={delta_c}...")

    # MPF case
    output_file = run_mpf_solver(V_tilde, delta_c, gamma)

    if output_file and Path(output_file).exists():
        solver_data = load_solver_result(output_file)
        lambda_D = get_lambda_D(output_file)
        x_a, rho_tilde = convert_to_bazant_units(solver_data, a_nm, gamma)

        print(f"  λD = {lambda_D:.4f} nm, a/λD = {a_nm/lambda_D:.3f}")

        # Left panel: Charge density
        ax1.plot(x_a, rho_tilde, 'b-', linewidth=1.5, label='This work (MPF)')

        ref_key = f'V{V_tilde}_dc{delta_c}'
        if ref_key in ref_data:
            ax1.plot(ref_data[ref_key]['x'], ref_data[ref_key]['y'],
                    'ko', markersize=3, label='Bazant 2011')

        ax1.axhline(y=0, color='gray', linestyle=':', linewidth=0.5)
        ax1.axhline(y=2, color='red', linestyle=':', linewidth=0.5, alpha=0.5, label='Crowding limit')
        ax1.set_xlim([0, 5])
        ax1.set_ylim([-1.5, 2.5])
        set_labels(ax1, r'$x/a$', r'$\tilde{\rho}$')
        ax1.legend(loc='lower right', fontsize=8)
        ax1.set_title(f'(a) Charge density: Ṽ={V_tilde}, δc={delta_c}', fontsize=10)
        setup_axis_style(ax1)

        # Right panel: Ion concentrations
        c_plus_norm = solver_data['c_plus_norm']
        c_minus_norm = solver_data['c_minus_norm']
        c_total = c_plus_norm + c_minus_norm

        ax2.plot(x_a, c_minus_norm, 'b-', linewidth=1.5, label=r'$\tilde{c}_-$')
        ax2.plot(x_a, c_plus_norm, 'r--', linewidth=1.5, label=r'$\tilde{c}_+$')
        ax2.plot(x_a, c_total, 'g-.', linewidth=1.5, label=r'$\tilde{c}_+ + \tilde{c}_-$')

        # Maximum concentration line
        c_max_norm = 1 / gamma  # cmax/c0 = 1/(γ) when γ = 2a³c₀NA
        ax2.axhline(y=c_max_norm, color='gray', linestyle=':', linewidth=0.5)

        ax2.set_xlim([0, 5])
        ax2.set_ylim([0, 5])
        set_labels(ax2, r'$x/a$', r'$\tilde{c}$')
        ax2.legend(loc='upper right', fontsize=8)
        ax2.set_title(f'(b) Ion concentrations: Ṽ={V_tilde}', fontsize=10)
        setup_axis_style(ax2)

    plt.tight_layout()

    output_dir = project_dir / 'results'
    plt.savefig(output_dir / 'bazant_fig2_detailed.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("Saved: bazant_fig2_detailed.png")


if __name__ == '__main__':
    # Test single high-voltage case first
    plot_single_case_comparison(V_tilde=100, delta_c=10)

    # Full comparison (may take longer)
    plot_fig2a_comparison()

    print("\nValidation plots saved to results/")
