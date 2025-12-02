#!/usr/bin/env python3
"""
Plot dual-electrode (capacitor) simulation results
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams['font.size'] = 11

# Load data
data = np.loadtxt('results/dual_electrode.dat', comments='#')

x_nm = data[:, 0]          # x [nm]
phi_mV = data[:, 2]        # phi [mV]
c_plus = data[:, 4]        # c+ [mol/m^3]
c_minus = data[:, 5]       # c- [mol/m^3]
c_plus_ratio = data[:, 6]  # c+/c0
c_minus_ratio = data[:, 7] # c-/c0

# Create 2x2 figure with overview and zoomed views
fig = plt.figure(figsize=(14, 10))

# Top left: Full potential profile
ax1 = fig.add_subplot(2, 2, 1)
ax1.plot(x_nm, phi_mV, 'b-', linewidth=2)
ax1.axhline(y=0, color='gray', linestyle='--', linewidth=0.5)
ax1.set_xlabel('x [nm]')
ax1.set_ylabel('φ [mV]')
ax1.set_title('Electric Potential Profile (Full Domain)')
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
ax2.set_title('Ion Concentration Profiles (Full Domain)')
ax2.set_xlim([0, 100])
ax2.set_ylim([1e-1, 1e1])
ax2.grid(True, alpha=0.3)
ax2.legend(loc='upper right')

# Bottom left: Left electrode EDL (zoomed)
ax3 = fig.add_subplot(2, 2, 3)
mask_left = x_nm <= 1.0
ax3.plot(x_nm[mask_left] * 1000, phi_mV[mask_left], 'b-', linewidth=2)  # Convert to pm for visibility
ax3.set_xlabel('x [pm]')
ax3.set_ylabel('φ [mV]')
ax3.set_title('Left EDL Structure (Anode, zoomed)')
ax3.set_xlim([0, 1000])
ax3.grid(True, alpha=0.3)

# Add second y-axis for concentration
ax3b = ax3.twinx()
ax3b.semilogy(x_nm[mask_left] * 1000, c_plus_ratio[mask_left], 'r--', linewidth=2, label='c$_+$/c$_0$')
ax3b.semilogy(x_nm[mask_left] * 1000, c_minus_ratio[mask_left], 'g--', linewidth=2, label='c$_-$/c$_0$')
ax3b.set_ylabel('c/c$_0$ [-]', color='green')
ax3b.set_ylim([1e-1, 1e1])
ax3b.legend(loc='right')

# Bottom right: Right electrode EDL (zoomed)
ax4 = fig.add_subplot(2, 2, 4)
mask_right = x_nm >= 99.0
x_right_pm = (100 - x_nm[mask_right]) * 1000  # Distance from right electrode in pm
ax4.plot(x_right_pm[::-1], phi_mV[mask_right][::-1], 'b-', linewidth=2)
ax4.set_xlabel('Distance from cathode [pm]')
ax4.set_ylabel('φ [mV]')
ax4.set_title('Right EDL Structure (Cathode, zoomed)')
ax4.set_xlim([0, 1000])
ax4.grid(True, alpha=0.3)

# Add second y-axis for concentration
ax4b = ax4.twinx()
ax4b.semilogy(x_right_pm[::-1], c_plus_ratio[mask_right][::-1], 'r--', linewidth=2, label='c$_+$/c$_0$')
ax4b.semilogy(x_right_pm[::-1], c_minus_ratio[mask_right][::-1], 'g--', linewidth=2, label='c$_-$/c$_0$')
ax4b.set_ylabel('c/c$_0$ [-]', color='green')
ax4b.set_ylim([1e-1, 1e1])
ax4b.legend(loc='right')

plt.suptitle('Dual-Electrode Model: Capacitor Structure\n(φ$_L$ = +50 mV, φ$_R$ = -50 mV, c$_0$ = 1 M, λ$_D$ = 0.12 nm)',
             fontsize=14, y=0.98)
plt.tight_layout()
plt.savefig('results/dual_electrode.png', dpi=150, bbox_inches='tight')
plt.close()

print("Saved: results/dual_electrode.png")

# Create additional detailed EDL structure plot
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
ax1.set_title('Left EDL (Anode, φ = +50 mV)\nAnion accumulation, Cation depletion')
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
ax2.set_title('Right EDL (Cathode, φ = -50 mV)\nCation accumulation, Anion depletion')
ax2.set_xlim([0, 500])
ax2.set_ylim([0, 8])
ax2.grid(True, alpha=0.3)
ax2.legend()
ax2.fill_between(x_right_pm[::-1], 0, c_plus_ratio[mask_right][::-1], alpha=0.2, color='red')

# Add Debye length marker
ax2.axvline(x=lambda_D_pm, color='purple', linestyle=':', linewidth=2, label=f'λ$_D$ = {lambda_D_pm} pm')
ax2.legend()

plt.suptitle('EDL Structure at Both Electrodes (Linear Scale)', fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig('results/dual_electrode_edl.png', dpi=150, bbox_inches='tight')
plt.close()

print("Saved: results/dual_electrode_edl.png")
