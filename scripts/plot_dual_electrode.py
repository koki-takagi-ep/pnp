#!/usr/bin/env python3
"""
Plot dual-electrode (capacitor) simulation results
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams['font.size'] = 12

# Load data
data = np.loadtxt('results/dual_electrode.dat', comments='#')

x_nm = data[:, 0]          # x [nm]
phi_mV = data[:, 2]        # phi [mV]
c_plus = data[:, 4]        # c+ [mol/m^3]
c_minus = data[:, 5]       # c- [mol/m^3]
c_plus_ratio = data[:, 6]  # c+/c0
c_minus_ratio = data[:, 7] # c-/c0

# Create figure
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Left plot: Potential profile
ax1 = axes[0]
ax1.plot(x_nm, phi_mV, 'b-', linewidth=2, label='PNP numerical')
ax1.axhline(y=0, color='gray', linestyle='--', linewidth=0.5)
ax1.set_xlabel('x [nm]')
ax1.set_ylabel('φ [mV]')
ax1.set_title('Electric Potential Profile')
ax1.set_xlim([0, 100])
ax1.set_ylim([-60, 60])
ax1.grid(True, alpha=0.3)
ax1.legend()

# Add annotations
ax1.annotate('φ = +50 mV\n(anode)', xy=(0, 50), xytext=(10, 45),
            fontsize=10, ha='left')
ax1.annotate('φ = -50 mV\n(cathode)', xy=(100, -50), xytext=(70, -45),
            fontsize=10, ha='left')

# Right plot: Ion concentration profiles
ax2 = axes[1]
ax2.semilogy(x_nm, c_plus_ratio, 'r-', linewidth=2, label='c$_+$/c$_0$ (cation)')
ax2.semilogy(x_nm, c_minus_ratio, 'b-', linewidth=2, label='c$_-$/c$_0$ (anion)')
ax2.axhline(y=1, color='gray', linestyle='--', linewidth=0.5, label='bulk')
ax2.set_xlabel('x [nm]')
ax2.set_ylabel('c/c$_0$ [-]')
ax2.set_title('Ion Concentration Profiles')
ax2.set_xlim([0, 100])
ax2.set_ylim([1e-2, 1e2])
ax2.grid(True, alpha=0.3)
ax2.legend()

# Add annotations for EDL regions
ax2.axvspan(0, 1, alpha=0.2, color='red', label='_nolegend_')
ax2.axvspan(99, 100, alpha=0.2, color='blue', label='_nolegend_')
ax2.text(0.5, 50, 'EDL\n(anode)', fontsize=9, ha='center', va='top')
ax2.text(99.5, 50, 'EDL\n(cathode)', fontsize=9, ha='center', va='top')

plt.suptitle('Dual-Electrode Model: Capacitor Structure\n(φ$_L$ = +50 mV, φ$_R$ = -50 mV, c$_0$ = 1 M)',
             fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig('results/dual_electrode.png', dpi=150, bbox_inches='tight')
plt.close()

print("Saved: results/dual_electrode.png")

# Also create a zoomed view of EDL regions
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Left EDL (anode)
ax1 = axes[0]
mask_left = x_nm <= 2
ax1.semilogy(x_nm[mask_left], c_plus_ratio[mask_left], 'r-', linewidth=2, label='c$_+$/c$_0$')
ax1.semilogy(x_nm[mask_left], c_minus_ratio[mask_left], 'b-', linewidth=2, label='c$_-$/c$_0$')
ax1.axhline(y=1, color='gray', linestyle='--', linewidth=0.5)
ax1.set_xlabel('x [nm]')
ax1.set_ylabel('c/c$_0$ [-]')
ax1.set_title('Left EDL (Anode, φ = +50 mV)')
ax1.set_xlim([0, 2])
ax1.set_ylim([1e-1, 1e1])
ax1.grid(True, alpha=0.3)
ax1.legend()

# Right EDL (cathode)
ax2 = axes[1]
mask_right = x_nm >= 98
ax2.semilogy(x_nm[mask_right], c_plus_ratio[mask_right], 'r-', linewidth=2, label='c$_+$/c$_0$')
ax2.semilogy(x_nm[mask_right], c_minus_ratio[mask_right], 'b-', linewidth=2, label='c$_-$/c$_0$')
ax2.axhline(y=1, color='gray', linestyle='--', linewidth=0.5)
ax2.set_xlabel('x [nm]')
ax2.set_ylabel('c/c$_0$ [-]')
ax2.set_title('Right EDL (Cathode, φ = -50 mV)')
ax2.set_xlim([98, 100])
ax2.set_ylim([1e-1, 1e1])
ax2.grid(True, alpha=0.3)
ax2.legend()

plt.suptitle('EDL Structure at Both Electrodes', fontsize=14, y=1.02)
plt.tight_layout()
plt.savefig('results/dual_electrode_edl.png', dpi=150, bbox_inches='tight')
plt.close()

print("Saved: results/dual_electrode_edl.png")
