#!/usr/bin/env python3
"""
Grid convergence analysis for 1D PNP solver
Plots L2 error vs number of cells to evaluate numerical accuracy
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import os

# Read convergence data
data_file = 'results/convergence_data.csv'
if not os.path.exists(data_file):
    print(f"Error: {data_file} not found. Run 'bash scripts/run_convergence.sh' first.")
    exit(1)

# Load data (skip header)
data = np.genfromtxt(data_file, delimiter=',', skip_header=1, dtype=None, encoding='utf-8')
N = np.array([row[0] for row in data])
dx = np.array([float(row[1]) for row in data])
L2_error = np.array([float(row[2]) for row in data])

# Calculate convergence order
order = []
for i in range(1, len(N)):
    p = np.log(L2_error[i-1] / L2_error[i]) / np.log(N[i] / N[i-1])
    order.append(p)

print("=" * 60)
print("Grid Convergence Analysis - 1D Poisson-Boltzmann Solver")
print("=" * 60)
print(f"{'N':>6} {'dx [nm]':>10} {'L2 [mV]':>14} {'Order':>8}")
print("-" * 60)
for i in range(len(N)):
    if i == 0:
        print(f"{N[i]:>6} {dx[i]:>10.4f} {L2_error[i]:>14.6e} {'--':>8}")
    else:
        print(f"{N[i]:>6} {dx[i]:>10.4f} {L2_error[i]:>14.6e} {order[i-1]:>8.2f}")

print("-" * 60)
print(f"Average convergence order (last 3): {np.mean(order[-3:]):.2f}")
print(f"Theoretical order: 2.00 (2nd-order central difference)")
print("=" * 60)

# Create figure
fig, ax = plt.subplots(figsize=(8, 6))

# Plot L2 error vs dx (grid spacing)
ax.loglog(dx, L2_error, 'bo-', markersize=10, linewidth=2, label='Numerical L2 Error')

# Add reference slopes
dx_ref = np.array([0.05, 3.0])
# First order reference (slope = 1 on log-log, error ~ dx^1)
L2_1st = L2_error[0] * (dx_ref / dx[0]) ** 1
ax.loglog(dx_ref, L2_1st, 'k--', linewidth=1.5, alpha=0.6, label='1st order slope')
# Second order reference (slope = 2 on log-log, error ~ dx^2)
L2_2nd = L2_error[0] * (dx_ref / dx[0]) ** 2
ax.loglog(dx_ref, L2_2nd, 'k:', linewidth=2, alpha=0.8, label='2nd order slope')

ax.set_xlabel(r'$\Delta x$ [nm]', fontsize=12)
ax.set_ylabel('L2 Error [mV]', fontsize=12)
ax.set_title('Grid Convergence Analysis\n1D Poisson-Boltzmann Equation (2nd-order FDM)', fontsize=14)
ax.legend(fontsize=10, loc='lower right')
ax.grid(True, which='both', linestyle='-', alpha=0.3)
ax.set_xlim([0.04, 3.0])

# Add text annotation for average order
avg_order = np.mean(order[-3:])
ax.text(0.05, 0.95, f'Observed order: {avg_order:.2f}\nTheoretical: 2.00',
        transform=ax.transAxes, fontsize=11,
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5),
        verticalalignment='top')

plt.tight_layout()

# Save figure
os.makedirs('results', exist_ok=True)
plt.savefig('results/grid_convergence.png', dpi=150, bbox_inches='tight')
plt.savefig('results/grid_convergence.svg', bbox_inches='tight')
print(f"\nFigures saved to:")
print(f"  - results/grid_convergence.png")
print(f"  - results/grid_convergence.svg")
