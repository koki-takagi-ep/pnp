#!/usr/bin/env python3
"""
Grid convergence analysis for 1D PNP solver
Plots L2 error vs number of cells to evaluate numerical accuracy
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

# Read convergence data
data = pd.read_csv('results/convergence_data.csv')
N = data['N'].values
L2_error = data['L2_error_mV'].values

# Calculate grid spacing (assuming uniform grid over 100 nm domain)
dx = 100.0 / (N - 1)  # nm

# Calculate convergence order
order = []
for i in range(1, len(N)):
    p = np.log(L2_error[i-1] / L2_error[i]) / np.log(N[i] / N[i-1])
    order.append(p)

print("Grid Convergence Analysis")
print("=" * 50)
print(f"{'N':>6} {'dx [nm]':>10} {'L2 [mV]':>12} {'Order':>8}")
print("-" * 50)
for i in range(len(N)):
    if i == 0:
        print(f"{N[i]:>6} {dx[i]:>10.4f} {L2_error[i]:>12.6f} {'--':>8}")
    else:
        print(f"{N[i]:>6} {dx[i]:>10.4f} {L2_error[i]:>12.6f} {order[i-1]:>8.2f}")

print("-" * 50)
print(f"Average convergence order: {np.mean(order):.2f}")

# Create figure
fig, ax = plt.subplots(figsize=(8, 6))

# Plot L2 error vs N
ax.loglog(N, L2_error, 'bo-', markersize=8, linewidth=2, label='L2 Error')

# Add reference slopes
N_ref = np.array([50, 2000])
# First order reference
L2_1st = L2_error[0] * (N[0] / N_ref) ** 1
ax.loglog(N_ref, L2_1st, 'k--', linewidth=1.5, alpha=0.7, label='1st order')
# Second order reference
L2_2nd = L2_error[0] * (N[0] / N_ref) ** 2
ax.loglog(N_ref, L2_2nd, 'k:', linewidth=1.5, alpha=0.7, label='2nd order')

ax.set_xlabel('Number of Cells N', fontsize=12)
ax.set_ylabel('L2 Error [mV]', fontsize=12)
ax.set_title('Grid Convergence Analysis\n(1D Poisson-Boltzmann Solver)', fontsize=14)
ax.legend(fontsize=10)
ax.grid(True, which='both', linestyle='-', alpha=0.3)
ax.set_xlim([40, 2500])

# Add text annotation for average order
avg_order = np.mean(order)
ax.text(0.05, 0.05, f'Average order: {avg_order:.2f}',
        transform=ax.transAxes, fontsize=11,
        bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

plt.tight_layout()

# Save figure
os.makedirs('results', exist_ok=True)
plt.savefig('results/grid_convergence.png', dpi=150, bbox_inches='tight')
plt.savefig('results/grid_convergence.svg', bbox_inches='tight')
print(f"\nFigure saved to results/grid_convergence.png")

plt.show()
