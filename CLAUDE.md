# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build and Run Commands

```bash
make              # Build the solver
make run          # Build and run with default parameters
make clean        # Clean build artifacts

# Run with specific options
./build/pnp_solver --phi0 100 --phi-right 0 --closed-system --c0 1.0

# Python environment (use .venv)
source .venv/bin/activate
python3 scripts/plot_dual_electrode.py
```

## Architecture Overview

This is a 1D Poisson-Nernst-Planck (PNP) solver for simulating electric double layers (EDL) in ionic liquids.

### Core Components

**Solver (`src/pnp_solver.cpp`, `include/pnp_solver.hpp`)**:
- `PNPSolver1D` class: Main solver with Newton-Raphson for steady-state and multiple transient schemes
- `solve()`: Steady-state Poisson-Boltzmann with optional charge neutrality constraint (closed system)
- `solve_transient_*()`: Various transient solvers (Gummel, Slotboom, Shen-Xu schemes)

**Key Physical Concepts**:
- **Open system**: Right boundary at bulk concentration (Dirichlet c=c₀)
- **Closed system** (`--closed-system`): Zero-flux at both boundaries, total charge conserved
  - Bulk potential φ_bulk is self-consistently determined from charge neutrality
  - For symmetric 1:1 electrolyte: φ_bulk = (φ_L + φ_R)/2

**Coordinate System**:
- Internal calculation uses bulk-centered frame (φ_bulk = 0)
- Output converts to user-specified coordinates (e.g., 100/0 mV)
- Boltzmann distribution uses bulk-relative potential: c± = c₀·exp(∓(φ-φ_bulk)/φ_T)

### Key Member Variables (PNPSolver1D)

- `phi_[]`: Electric potential array [V] - stored in user coordinate system after solve()
- `phi_bulk_`: Self-consistently determined bulk potential [V] (non-zero for closed system)
- `c_plus_[]`, `c_minus_[]`: Ion concentrations [mol/m³]
- `lambda_D_`: Debye length [m]
- `phi_T_`: Thermal voltage kT/e ≈ 25.7 mV at 298K

### Grid System

Non-uniform grid with tanh-stretching near electrodes:
- `--dual-electrode`: Symmetric grid, fine mesh near both electrodes
- Single-electrode (default): Fine mesh only near left electrode (x=0)

## Project Rules

### README Updates

**When adding/updating figures or plots, always update README.md:**
- Add new figures in `results/` to the relevant README section
- Use sequential numbering (Figure 1, Figure 2, ...)
- Include descriptive captions
- Use `?v=N` query parameter to bust image cache when updating

### Pre-commit Checklist

1. Commit code changes
2. Verify README.md is up to date
3. Include figures and data files in commit
4. Push to branch

## File Structure

- `src/pnp_solver.cpp`: Solver implementation (~2700 lines)
- `src/main.cpp`: CLI entry point with argument parsing
- `include/pnp_solver.hpp`: Class definition and parameters struct
- `scripts/plot_*.py`: Visualization scripts
- `results/`: Output data (.dat) and figures (.png)
- `docs/`: Reference papers (PDF)
