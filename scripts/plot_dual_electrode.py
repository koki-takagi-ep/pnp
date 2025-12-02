#!/usr/bin/env python3
"""
Plot dual-electrode (capacitor) simulation results
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams['font.size'] = 11

# EMI-BF4 ionic liquid physical parameters
lambda_D_pm = 119      # Debye length [pm]
d_EMI_pm = 760         # EMI+ cation diameter [pm] (~0.76 nm)
d_BF4_pm = 460         # BF4- anion diameter [pm] (~0.46 nm)

# Load data
data = np.loadtxt('results/dual_electrode.dat', comments='#')

x_nm = data[:, 0]          # x [nm]
phi_mV = data[:, 2]        # phi [mV]
c_plus = data[:, 4]        # c+ [mol/m^3]
c_minus = data[:, 5]       # c- [mol/m^3]
c_plus_ratio = data[:, 6]  # c+/c0
c_minus_ratio = data[:, 7] # c-/c0
dx_nm = data[:, 10]        # dx [nm] (grid spacing)

# ============================================================================
# Figure 1: Capacitor Structure Overview (2x2 layout)
# ============================================================================
fig = plt.figure(figsize=(14, 10))

# Top left: Full potential profile
ax1 = fig.add_subplot(2, 2, 1)
ax1.plot(x_nm, phi_mV, 'b-', linewidth=2)
ax1.axhline(y=0, color='gray', linestyle='--', linewidth=0.5)
ax1.set_xlabel('x [nm]')
ax1.set_ylabel('φ [mV]')
ax1.set_title('(a) Electric Potential Profile (Full Domain)')
ax1.set_xlim([0, 100])
ax1.set_ylim([-60, 60])
ax1.grid(True, alpha=0.3)
ax1.annotate('Anode\n+50 mV', xy=(0, 50), xytext=(5, 40), fontsize=10, color='red')
ax1.annotate('Cathode\n-50 mV', xy=(100, -50), xytext=(80, -40), fontsize=10, color='blue')

# Top right: Full concentration profile
ax2 = fig.add_subplot(2, 2, 2)
ax2.semilogy(x_nm, c_plus_ratio, 'r-', linewidth=2, label='c$_+$/c$_0$ (cation)')
ax2.semilogy(x_nm, c_minus_ratio, 'b-', linewidth=2, label='c$_-$/c$_0$ (anion)')
ax2.axhline(y=1, color='gray', linestyle='--', linewidth=0.5, label='bulk')
ax2.set_xlabel('x [nm]')
ax2.set_ylabel('c/c$_0$ [-]')
ax2.set_title('(b) Ion Concentration Profiles (Full Domain)')
ax2.set_xlim([0, 100])
ax2.set_ylim([1e-1, 1e1])
ax2.grid(True, alpha=0.3)
ax2.legend(loc='upper right')

# Bottom left: Left electrode EDL (zoomed) - Linear scale
ax3 = fig.add_subplot(2, 2, 3)
mask_left = x_nm <= 0.5
x_left_pm = x_nm[mask_left] * 1000  # pm

ax3.plot(x_left_pm, c_plus_ratio[mask_left], 'r-', linewidth=2, marker='o', markersize=3, label='c$_+$/c$_0$ (EMI$^+$)')
ax3.plot(x_left_pm, c_minus_ratio[mask_left], 'b-', linewidth=2, marker='s', markersize=3, label='c$_-$/c$_0$ (BF$_4^-$)')
ax3.axhline(y=1, color='gray', linestyle='--', linewidth=0.5, label='bulk')
ax3.set_xlabel('Distance from anode [pm]')
ax3.set_ylabel('c/c$_0$ [-]')
ax3.set_title('(c) Left EDL (Anode, φ = +50 mV)\nBF$_4^-$ accumulation, EMI$^+$ depletion')
ax3.set_xlim([0, 500])
ax3.set_ylim([0, 8])
ax3.grid(True, alpha=0.3)
ax3.fill_between(x_left_pm, 0, c_minus_ratio[mask_left], alpha=0.2, color='blue')

# Add characteristic length markers
ax3.axvline(x=lambda_D_pm, color='purple', linestyle=':', linewidth=2, label=f'λ$_D$ = {lambda_D_pm} pm')
ax3.axvline(x=d_EMI_pm, color='red', linestyle='-.', linewidth=1.5, alpha=0.7, label=f'd$_{{EMI^+}}$ = {d_EMI_pm} pm')
ax3.axvline(x=d_BF4_pm, color='blue', linestyle='-.', linewidth=1.5, alpha=0.7, label=f'd$_{{BF_4^-}}$ = {d_BF4_pm} pm')
ax3.legend(loc='upper right', fontsize=9)

# Bottom right: Right electrode EDL (zoomed) - Linear scale
ax4 = fig.add_subplot(2, 2, 4)
mask_right = x_nm >= 99.5
x_right_pm = (100 - x_nm[mask_right]) * 1000  # Distance from cathode in pm

ax4.plot(x_right_pm[::-1], c_plus_ratio[mask_right][::-1], 'r-', linewidth=2, marker='o', markersize=3, label='c$_+$/c$_0$ (EMI$^+$)')
ax4.plot(x_right_pm[::-1], c_minus_ratio[mask_right][::-1], 'b-', linewidth=2, marker='s', markersize=3, label='c$_-$/c$_0$ (BF$_4^-$)')
ax4.axhline(y=1, color='gray', linestyle='--', linewidth=0.5, label='bulk')
ax4.set_xlabel('Distance from cathode [pm]')
ax4.set_ylabel('c/c$_0$ [-]')
ax4.set_title('(d) Right EDL (Cathode, φ = -50 mV)\nEMI$^+$ accumulation, BF$_4^-$ depletion')
ax4.set_xlim([0, 500])
ax4.set_ylim([0, 8])
ax4.grid(True, alpha=0.3)
ax4.fill_between(x_right_pm[::-1], 0, c_plus_ratio[mask_right][::-1], alpha=0.2, color='red')

# Add characteristic length markers
ax4.axvline(x=lambda_D_pm, color='purple', linestyle=':', linewidth=2, label=f'λ$_D$ = {lambda_D_pm} pm')
ax4.axvline(x=d_EMI_pm, color='red', linestyle='-.', linewidth=1.5, alpha=0.7, label=f'd$_{{EMI^+}}$ = {d_EMI_pm} pm')
ax4.axvline(x=d_BF4_pm, color='blue', linestyle='-.', linewidth=1.5, alpha=0.7, label=f'd$_{{BF_4^-}}$ = {d_BF4_pm} pm')
ax4.legend(loc='upper right', fontsize=9)

plt.suptitle('Figure 1: Dual-Electrode Model - EMI-BF$_4$ Ionic Liquid\n(φ$_L$ = +50 mV, φ$_R$ = -50 mV, c$_0$ = 1 M, λ$_D$ = 119 pm)',
             fontsize=14, y=0.98)
plt.tight_layout()
plt.savefig('results/dual_electrode.png', dpi=150, bbox_inches='tight')
plt.close()

print("Saved: results/dual_electrode.png")

# ============================================================================
# Figure 2: EDL Structure at Both Electrodes (Linear Scale)
# ============================================================================
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# Left EDL with linear scale
ax1 = axes[0]
mask_left = x_nm <= 0.5
x_left_pm = x_nm[mask_left] * 1000  # pm

ax1.plot(x_left_pm, c_plus_ratio[mask_left], 'r-', linewidth=2, marker='o', markersize=3, label='c$_+$/c$_0$ (cation)')
ax1.plot(x_left_pm, c_minus_ratio[mask_left], 'b-', linewidth=2, marker='s', markersize=3, label='c$_-$/c$_0$ (anion)')
ax1.axhline(y=1, color='gray', linestyle='--', linewidth=0.5, label='bulk')
ax1.set_xlabel('Distance from anode [pm]')
ax1.set_ylabel('c/c$_0$ [-]')
ax1.set_title('(a) Left EDL (Anode, φ = +50 mV)\nAnion accumulation, Cation depletion')
ax1.set_xlim([0, 500])
ax1.set_ylim([0, 8])
ax1.grid(True, alpha=0.3)
ax1.legend()
ax1.fill_between(x_left_pm, 0, c_minus_ratio[mask_left], alpha=0.2, color='blue')

# Add Debye length marker
lambda_D_pm = 119  # pm
ax1.axvline(x=lambda_D_pm, color='purple', linestyle=':', linewidth=2, label=f'λ$_D$ = {lambda_D_pm} pm')
ax1.legend()

# Right EDL with linear scale
ax2 = axes[1]
mask_right = x_nm >= 99.5
x_right_pm = (100 - x_nm[mask_right]) * 1000  # Distance from cathode in pm

ax2.plot(x_right_pm[::-1], c_plus_ratio[mask_right][::-1], 'r-', linewidth=2, marker='o', markersize=3, label='c$_+$/c$_0$ (cation)')
ax2.plot(x_right_pm[::-1], c_minus_ratio[mask_right][::-1], 'b-', linewidth=2, marker='s', markersize=3, label='c$_-$/c$_0$ (anion)')
ax2.axhline(y=1, color='gray', linestyle='--', linewidth=0.5, label='bulk')
ax2.set_xlabel('Distance from cathode [pm]')
ax2.set_ylabel('c/c$_0$ [-]')
ax2.set_title('(b) Right EDL (Cathode, φ = -50 mV)\nCation accumulation, Anion depletion')
ax2.set_xlim([0, 500])
ax2.set_ylim([0, 8])
ax2.grid(True, alpha=0.3)
ax2.legend()
ax2.fill_between(x_right_pm[::-1], 0, c_plus_ratio[mask_right][::-1], alpha=0.2, color='red')

# Add Debye length marker
ax2.axvline(x=lambda_D_pm, color='purple', linestyle=':', linewidth=2, label=f'λ$_D$ = {lambda_D_pm} pm')
ax2.legend()

plt.suptitle('Figure 2: EDL Structure at Both Electrodes (Linear Scale)', fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig('results/dual_electrode_edl.png', dpi=150, bbox_inches='tight')
plt.close()

print("Saved: results/dual_electrode_edl.png")

# ============================================================================
# Figure 3: Grid Distribution (2x2 layout)
# ============================================================================
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Top left: Grid spacing across full domain (log scale)
ax1 = axes[0, 0]
ax1.semilogy(x_nm, dx_nm * 1000, 'b-', linewidth=2, marker='o', markersize=2)
ax1.set_xlabel('x [nm]')
ax1.set_ylabel('Δx [pm]')
ax1.set_title('(a) Grid Spacing Distribution (Full Domain)')
ax1.set_xlim([0, 100])
ax1.grid(True, alpha=0.3)
ax1.axhline(y=119, color='purple', linestyle=':', linewidth=2, label=f'λ$_D$ = 119 pm')
ax1.legend()

# Top right: Grid point density (points per nm)
ax2 = axes[0, 1]
# Calculate point density as 1/dx
point_density = 1.0 / dx_nm  # points per nm
ax2.semilogy(x_nm, point_density, 'r-', linewidth=2)
ax2.set_xlabel('x [nm]')
ax2.set_ylabel('Grid point density [points/nm]')
ax2.set_title('(b) Grid Point Density (Full Domain)')
ax2.set_xlim([0, 100])
ax2.grid(True, alpha=0.3)

# Bottom left: Zoomed view near left electrode
ax3 = axes[1, 0]
mask_left = x_nm <= 1.0
x_left_pm = x_nm[mask_left] * 1000
dx_left_pm = dx_nm[mask_left] * 1000

ax3.plot(x_left_pm, dx_left_pm, 'b-', linewidth=2, marker='o', markersize=4)
ax3.set_xlabel('Distance from anode [pm]')
ax3.set_ylabel('Δx [pm]')
ax3.set_title('(c) Grid Spacing Near Left Electrode (Anode)')
ax3.set_xlim([0, 1000])
ax3.grid(True, alpha=0.3)
ax3.axvline(x=119, color='purple', linestyle=':', linewidth=2, label=f'λ$_D$ = 119 pm')
ax3.axhline(y=119, color='purple', linestyle=':', linewidth=2, alpha=0.5)
ax3.legend()

# Bottom right: Zoomed view near right electrode
ax4 = axes[1, 1]
mask_right = x_nm >= 99.0
x_right_pm = (100 - x_nm[mask_right]) * 1000  # Distance from cathode
dx_right_pm = dx_nm[mask_right] * 1000

ax4.plot(x_right_pm[::-1], dx_right_pm[::-1], 'b-', linewidth=2, marker='o', markersize=4)
ax4.set_xlabel('Distance from cathode [pm]')
ax4.set_ylabel('Δx [pm]')
ax4.set_title('(d) Grid Spacing Near Right Electrode (Cathode)')
ax4.set_xlim([0, 1000])
ax4.grid(True, alpha=0.3)
ax4.axvline(x=119, color='purple', linestyle=':', linewidth=2, label=f'λ$_D$ = 119 pm')
ax4.axhline(y=119, color='purple', linestyle=':', linewidth=2, alpha=0.5)
ax4.legend()

plt.suptitle('Figure 3: Grid Distribution for Dual-Electrode Model\n(Symmetric tanh-stretching, fine near both electrodes)',
             fontsize=14, y=0.98)
plt.tight_layout()
plt.savefig('results/dual_electrode_grid.png', dpi=150, bbox_inches='tight')
plt.close()

print("Saved: results/dual_electrode_grid.png")

# ============================================================================
# Figure 4: Grid Point Position Visualization
# ============================================================================
fig, ax = plt.subplots(figsize=(14, 4))

# Plot grid points as vertical lines
for x in x_nm:
    ax.axvline(x=x, color='blue', linewidth=0.5, alpha=0.5)

# Mark electrode positions
ax.axvline(x=0, color='red', linewidth=3, label='Anode (+50 mV)')
ax.axvline(x=100, color='blue', linewidth=3, label='Cathode (-50 mV)')

# Highlight EDL regions (5*lambda_D from each electrode)
lambda_D_nm = 0.119
edl_width = 5 * lambda_D_nm
ax.axvspan(0, edl_width, alpha=0.2, color='red', label=f'EDL region (5λ$_D$ = {edl_width*1000:.0f} pm)')
ax.axvspan(100 - edl_width, 100, alpha=0.2, color='blue')

ax.set_xlabel('x [nm]')
ax.set_ylabel('')
ax.set_yticks([])
ax.set_title(f'Figure 4: Grid Point Distribution (N = {len(x_nm)} points)')
ax.set_xlim([-1, 101])
ax.legend(loc='upper center', ncol=3)
ax.grid(True, alpha=0.3, axis='x')

plt.tight_layout()
plt.savefig('results/dual_electrode_grid_points.png', dpi=150, bbox_inches='tight')
plt.close()

print("Saved: results/dual_electrode_grid_points.png")

# Print grid statistics
print("\n--- Grid Statistics ---")
print(f"Total grid points: {len(x_nm)}")
print(f"Min Δx: {dx_nm.min()*1000:.2f} pm (at electrodes)")
print(f"Max Δx: {dx_nm.max()*1000:.2f} pm (at bulk center)")
print(f"Δx ratio (max/min): {dx_nm.max()/dx_nm.min():.1f}")

# Count points in EDL regions
n_left_edl = np.sum(x_nm <= edl_width)
n_right_edl = np.sum(x_nm >= 100 - edl_width)
print(f"Points within 5λ_D of left electrode: {n_left_edl}")
print(f"Points within 5λ_D of right electrode: {n_right_edl}")
