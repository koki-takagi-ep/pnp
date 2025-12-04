#!/usr/bin/env python3
"""
Plot dual-electrode (capacitor) simulation results
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

# EMI-BF4 ionic liquid physical parameters
lambda_D_pm = 119      # Debye length [pm]
d_EMI_pm = 760         # EMI+ cation diameter [pm] (~0.76 nm)
d_BF4_pm = 460         # BF4- anion diameter [pm] (~0.46 nm)

# Load data
data_file = project_dir / 'results' / 'dual_electrode.dat'
data = np.loadtxt(data_file, comments='#')

x_nm = data[:, 0]          # x [nm]
phi_mV = data[:, 2]        # phi [mV]
c_plus = data[:, 4]        # c+ [mol/m^3]
c_minus = data[:, 5]       # c- [mol/m^3]
c_plus_ratio = data[:, 6]  # c+/c0
c_minus_ratio = data[:, 7] # c-/c0
dx_nm = data[:, 10]        # dx [nm] (grid spacing)

# Extract actual electrode potentials from data
phi_left_mV = phi_mV[0]    # Anode (left electrode) potential
phi_right_mV = phi_mV[-1]  # Cathode (right electrode) potential

# Try to extract bulk potential from header comments
phi_bulk_mV = 0.0
try:
    with open(data_file, 'r') as f:
        for line in f:
            if line.startswith('#') and 'Bulk potential:' in line:
                phi_bulk_mV = float(line.split(':')[1].replace('mV', '').strip())
                break
except:
    phi_bulk_mV = phi_mV[len(phi_mV)//2]

print(f"Electrode potentials: Left={phi_left_mV:.1f} mV, Right={phi_right_mV:.1f} mV")
print(f"Bulk potential: {phi_bulk_mV:.1f} mV")

results_dir = project_dir / 'results'

# ============================================================================
# Figure 1: Capacitor Structure Overview (2x2 layout)
# ============================================================================
setup_plot_style()
fig, axes = plt.subplots(2, 2, figsize=(7, 7))

# (a) Full potential profile
ax1 = axes[0, 0]
ax1.plot(x_nm, phi_mV, 'b-', linewidth=1)
ax1.axhline(y=phi_bulk_mV, color='gray', linestyle='--', linewidth=0.5, label=f'Bulk ({phi_bulk_mV:.0f} mV)')
set_labels(ax1, r'$x$ (nm)', r'$\phi$ (mV)')
ax1.set_title('(a) Electric Potential', fontsize=10, fontweight='bold')
L_nm = x_nm[-1]  # Domain length from data
ax1.set_xlim([0, L_nm])
phi_min, phi_max = min(phi_mV.min(), phi_right_mV), max(phi_mV.max(), phi_left_mV)
phi_margin = (phi_max - phi_min) * 0.1
ax1.set_ylim([phi_min - phi_margin, phi_max + phi_margin])
ax1.legend(loc='best', fontsize=8)
ax1.annotate(f'Anode\n{phi_left_mV:+.0f} mV', xy=(0, phi_left_mV), xytext=(2, phi_left_mV-10), fontsize=8, color='red')
ax1.annotate(f'Cathode\n{phi_right_mV:+.0f} mV', xy=(L_nm, phi_right_mV), xytext=(L_nm*0.7, phi_right_mV+10), fontsize=8, color='blue')
setup_axis_style(ax1)

# (b) Full concentration profile
ax2 = axes[0, 1]
ax2.semilogy(x_nm, c_plus_ratio, 'r-', linewidth=1, label=r'$n_+/n_0$ (cation)')
ax2.semilogy(x_nm, c_minus_ratio, 'b-', linewidth=1, label=r'$n_-/n_0$ (anion)')
ax2.axhline(y=1, color='gray', linestyle='--', linewidth=0.5)
set_labels(ax2, r'$x$ (nm)', r'$n/n_0$')
ax2.set_title('(b) Ion Concentrations', fontsize=10, fontweight='bold')
ax2.set_xlim([0, L_nm])
ax2.set_ylim([1e-1, 1e1])
ax2.legend(loc='upper right', fontsize=8)
setup_axis_style(ax2)

# (c) Left electrode EDL (zoomed)
ax3 = axes[1, 0]
mask_left = x_nm <= 0.5
x_left_pm = x_nm[mask_left] * 1000  # pm

ax3.plot(x_left_pm, phi_mV[mask_left], 'k-', linewidth=1, label=r'$\phi$ (mV)')
ax3.axhline(y=phi_bulk_mV, color='gray', linestyle='--', linewidth=0.5)
set_labels(ax3, 'Distance from anode (pm)', r'$\phi$ (mV)')
ax3.set_title(f'(c) Left EDL (Anode)', fontsize=10, fontweight='bold')
ax3.set_xlim([0, 500])
left_edl_ylim = [min(phi_bulk_mV, phi_left_mV) - 5, max(phi_bulk_mV, phi_left_mV) + 5]
ax3.set_ylim(left_edl_ylim)

ax3.axvline(x=lambda_D_pm, color='purple', linestyle=':', linewidth=1, label=f'$\\lambda_D$ = {lambda_D_pm} pm')
setup_axis_style(ax3)

ax3b = ax3.twinx()
ax3b.plot(x_left_pm, c_plus_ratio[mask_left], 'r--', linewidth=1, marker='o', markersize=2, label=r'$n_+/n_0$')
ax3b.plot(x_left_pm, c_minus_ratio[mask_left], 'b--', linewidth=1, marker='s', markersize=2, label=r'$n_-/n_0$')
ax3b.set_ylabel(r'$n/n_0$', fontsize=10, fontweight='bold')
ax3b.set_ylim([0, 8])
ax3b.fill_between(x_left_pm, 0, c_minus_ratio[mask_left], alpha=0.15, color='blue')

lines1, labels1 = ax3.get_legend_handles_labels()
lines2, labels2 = ax3b.get_legend_handles_labels()
ax3.legend(lines1 + lines2, labels1 + labels2, loc='right', fontsize=7)

# (d) Right electrode EDL (zoomed, mirror image)
ax4 = axes[1, 1]
mask_right = x_nm >= (L_nm - 0.5)
x_right_pm = (L_nm - x_nm[mask_right]) * 1000

ax4.plot(x_right_pm[::-1], phi_mV[mask_right][::-1], 'k-', linewidth=1, label=r'$\phi$ (mV)')
ax4.axhline(y=phi_bulk_mV, color='gray', linestyle='--', linewidth=0.5)
set_labels(ax4, 'Distance from cathode (pm)', r'$\phi$ (mV)')
ax4.set_title(f'(d) Right EDL (Cathode)', fontsize=10, fontweight='bold')
ax4.set_xlim([500, 0])  # Inverted x-axis
right_edl_ylim = [min(phi_bulk_mV, phi_right_mV) - 5, max(phi_bulk_mV, phi_right_mV) + 5]
ax4.set_ylim(right_edl_ylim)

ax4.axvline(x=lambda_D_pm, color='purple', linestyle=':', linewidth=1, label=f'$\\lambda_D$ = {lambda_D_pm} pm')
setup_axis_style(ax4)

ax4b = ax4.twinx()
ax4b.plot(x_right_pm[::-1], c_plus_ratio[mask_right][::-1], 'r--', linewidth=1, marker='o', markersize=2, label=r'$n_+/n_0$')
ax4b.plot(x_right_pm[::-1], c_minus_ratio[mask_right][::-1], 'b--', linewidth=1, marker='s', markersize=2, label=r'$n_-/n_0$')
ax4b.set_ylabel(r'$n/n_0$', fontsize=10, fontweight='bold')
ax4b.set_ylim([0, 8])
ax4b.fill_between(x_right_pm[::-1], 0, c_plus_ratio[mask_right][::-1], alpha=0.15, color='red')

lines1, labels1 = ax4.get_legend_handles_labels()
lines2, labels2 = ax4b.get_legend_handles_labels()
ax4.legend(lines1 + lines2, labels1 + labels2, loc='center left', fontsize=7)

plt.tight_layout()
plt.savefig(results_dir / 'dual_electrode.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: results/dual_electrode.png")

# ============================================================================
# Figure 2: EDL Structure at Both Electrodes (Linear Scale)
# ============================================================================
setup_plot_style()
fig, axes = plt.subplots(1, 2, figsize=(7, 3.5))

# Left EDL
ax1 = axes[0]
mask_left = x_nm <= 0.5
x_left_pm = x_nm[mask_left] * 1000

ax1.plot(x_left_pm, c_plus_ratio[mask_left], 'r-', linewidth=1, marker='o', markersize=2, label=r'$n_+/n_0$ (cation)')
ax1.plot(x_left_pm, c_minus_ratio[mask_left], 'b-', linewidth=1, marker='s', markersize=2, label=r'$n_-/n_0$ (anion)')
ax1.axhline(y=1, color='gray', linestyle='--', linewidth=0.5)
set_labels(ax1, 'Distance from anode (pm)', r'$n/n_0$')
ax1.set_title(f'(a) Left EDL (Anode, $\\phi$ = {phi_left_mV:+.0f} mV)', fontsize=10, fontweight='bold')
ax1.set_xlim([0, 500])
ax1.set_ylim([0, 8])
ax1.fill_between(x_left_pm, 0, c_minus_ratio[mask_left], alpha=0.2, color='blue')
ax1.axvline(x=lambda_D_pm, color='purple', linestyle=':', linewidth=1, label=f'$\\lambda_D$ = {lambda_D_pm} pm')
ax1.legend(fontsize=7)
setup_axis_style(ax1)

# Right EDL (mirror image)
ax2 = axes[1]
mask_right = x_nm >= (L_nm - 0.5)
x_right_pm = (L_nm - x_nm[mask_right]) * 1000

ax2.plot(x_right_pm[::-1], c_plus_ratio[mask_right][::-1], 'r-', linewidth=1, marker='o', markersize=2, label=r'$n_+/n_0$ (cation)')
ax2.plot(x_right_pm[::-1], c_minus_ratio[mask_right][::-1], 'b-', linewidth=1, marker='s', markersize=2, label=r'$n_-/n_0$ (anion)')
ax2.axhline(y=1, color='gray', linestyle='--', linewidth=0.5)
set_labels(ax2, 'Distance from cathode (pm)', r'$n/n_0$')
ax2.set_title(f'(b) Right EDL (Cathode, $\\phi$ = {phi_right_mV:+.0f} mV)', fontsize=10, fontweight='bold')
ax2.set_xlim([500, 0])  # Inverted x-axis
ax2.set_ylim([0, 8])
ax2.fill_between(x_right_pm[::-1], 0, c_plus_ratio[mask_right][::-1], alpha=0.2, color='red')
ax2.axvline(x=lambda_D_pm, color='purple', linestyle=':', linewidth=1, label=f'$\\lambda_D$ = {lambda_D_pm} pm')
ax2.legend(loc='upper left', fontsize=7)
setup_axis_style(ax2)

plt.tight_layout()
plt.savefig(results_dir / 'dual_electrode_edl.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: results/dual_electrode_edl.png")

# ============================================================================
# Figure 3: Grid Distribution (2x2 layout)
# ============================================================================
setup_plot_style()
fig, axes = plt.subplots(2, 2, figsize=(7, 7))

# (a) Grid spacing (full domain, log scale)
ax1 = axes[0, 0]
ax1.semilogy(x_nm, dx_nm * 1000, 'b-', linewidth=1, marker='o', markersize=2)
set_labels(ax1, r'$x$ (nm)', r'$\Delta x$ (pm)')
ax1.set_title('(a) Grid Spacing (Full Domain)', fontsize=10, fontweight='bold')
ax1.set_xlim([0, L_nm])
ax1.axhline(y=119, color='purple', linestyle=':', linewidth=1, label=f'$\\lambda_D$ = 119 pm')
ax1.legend(fontsize=8)
setup_axis_style(ax1)

# (b) Grid point density
ax2 = axes[0, 1]
point_density = 1.0 / dx_nm
ax2.semilogy(x_nm, point_density, 'r-', linewidth=1)
set_labels(ax2, r'$x$ (nm)', 'Grid point density (points/nm)')
ax2.set_title('(b) Grid Point Density', fontsize=10, fontweight='bold')
ax2.set_xlim([0, L_nm])
setup_axis_style(ax2)

# (c) Left electrode zoom
ax3 = axes[1, 0]
mask_left = x_nm <= 1.0
x_left_pm = x_nm[mask_left] * 1000
dx_left_pm = dx_nm[mask_left] * 1000

ax3.plot(x_left_pm, dx_left_pm, 'b-', linewidth=1, marker='o', markersize=3)
set_labels(ax3, 'Distance from anode (pm)', r'$\Delta x$ (pm)')
ax3.set_title('(c) Grid Near Left Electrode', fontsize=10, fontweight='bold')
ax3.set_xlim([0, 1000])
ax3.axvline(x=119, color='purple', linestyle=':', linewidth=1, label=f'$\\lambda_D$ = 119 pm')
ax3.axhline(y=119, color='purple', linestyle=':', linewidth=1, alpha=0.5)
ax3.legend(fontsize=8)
setup_axis_style(ax3)

# (d) Right electrode zoom (mirror image)
ax4 = axes[1, 1]
mask_right = x_nm >= (L_nm - 1.0)
x_right_pm = (L_nm - x_nm[mask_right]) * 1000
dx_right_pm = dx_nm[mask_right] * 1000

ax4.plot(x_right_pm[::-1], dx_right_pm[::-1], 'b-', linewidth=1, marker='o', markersize=3)
set_labels(ax4, 'Distance from cathode (pm)', r'$\Delta x$ (pm)')
ax4.set_title('(d) Grid Near Right Electrode', fontsize=10, fontweight='bold')
ax4.set_xlim([1000, 0])  # Inverted x-axis
ax4.axvline(x=119, color='purple', linestyle=':', linewidth=1, label=f'$\\lambda_D$ = 119 pm')
ax4.axhline(y=119, color='purple', linestyle=':', linewidth=1, alpha=0.5)
ax4.legend(loc='upper left', fontsize=8)
setup_axis_style(ax4)

plt.tight_layout()
plt.savefig(results_dir / 'dual_electrode_grid.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: results/dual_electrode_grid.png")

# ============================================================================
# Figure 4: Grid Point Position Visualization
# ============================================================================
setup_plot_style()
fig, ax = plt.subplots(figsize=(7, 2))

for x in x_nm:
    ax.axvline(x=x, color='black', linewidth=0.3, alpha=0.3)

ax.axvline(x=0, color='red', linewidth=3, label=f'Anode ({phi_left_mV:+.0f} mV)')
ax.axvline(x=L_nm, color='blue', linewidth=3, label=f'Cathode ({phi_right_mV:+.0f} mV)')

lambda_D_nm = 0.119
edl_width = 5 * lambda_D_nm
ax.axvspan(0, edl_width, alpha=0.2, color='red', label=f'EDL region (5$\\lambda_D$ = {edl_width*1000:.0f} pm)')
ax.axvspan(L_nm - edl_width, L_nm, alpha=0.2, color='blue')

set_labels(ax, r'$x$ (nm)', '')
ax.set_yticks([])
ax.set_xlim([-1, L_nm + 1])
ax.legend(loc='upper center', ncol=3, fontsize=8)
setup_axis_style(ax)

plt.tight_layout()
plt.savefig(results_dir / 'dual_electrode_grid_points.png', dpi=300, bbox_inches='tight')
plt.close()
print("Saved: results/dual_electrode_grid_points.png")

# Print grid statistics
print("\n--- Grid Statistics ---")
print(f"Total grid points: {len(x_nm)}")
print(f"Min Δx: {dx_nm.min()*1000:.2f} pm (at electrodes)")
print(f"Max Δx: {dx_nm.max()*1000:.2f} pm (at bulk center)")
print(f"Δx ratio (max/min): {dx_nm.max()/dx_nm.min():.1f}")

n_left_edl = np.sum(x_nm <= edl_width)
n_right_edl = np.sum(x_nm >= L_nm - edl_width)
print(f"Points within 5λ_D of left electrode: {n_left_edl}")
print(f"Points within 5λ_D of right electrode: {n_right_edl}")
