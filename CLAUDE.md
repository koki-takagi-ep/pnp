# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build and Run Commands

```bash
make              # Build the solver
make run          # Build and run with default parameters
make clean        # Clean build artifacts

# Run with specific options
./build/pnp_solver --phi0 100 --phi-right 0 --closed-system --c0 1.0

# Transient solvers
./build/pnp_solver --gummel --dt 0.01 --t-final 0.2    # Gummel iteration
./build/pnp_solver --slotboom --dt 0.01 --t-final 0.2  # Slotboom transformation
./build/pnp_solver --shenxu --dt 0.01 --t-final 0.2    # Shen-Xu positivity-preserving

# Bikerman model (steric effects)
./build/pnp_solver --model bikerman --ion-size 0.7

# Python environment (use .venv)
source .venv/bin/activate
python3 scripts/plot_dual_electrode.py

# Convergence tests
bash scripts/run_convergence.sh           # Basic convergence study
bash scripts/run_comprehensive_convergence.sh  # Full parametric study
```

## Key CLI Options

| Option | Description | Default |
|--------|-------------|---------|
| `--phi0 <mV>` | Left electrode potential | 100 |
| `--phi-right <mV>` | Right electrode potential | 0 |
| `--c0 <mol/L>` | Bulk concentration | 1.0 |
| `--L <nm>` | Domain length | 50 |
| `--eps <value>` | Relative permittivity | 12 |
| `--closed-system` | Zero-flux BC at both ends | off |
| `--model <type>` | standard or bikerman | standard |
| `--N <points>` | Grid points | 1001 |
| `--stretch <factor>` | Grid stretching factor | 3.0 |

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

### Workflow: Commit and Push

**作業が完了したら必ず以下の手順を実行すること：**

1. README.md を更新（必要に応じて）
2. 変更をコミット
3. リモートにプッシュ

**これを徹底すること。作業後にプッシュを忘れないこと。**

### README Updates

**When adding/updating figures or plots, always update README.md:**
- Add new figures in `results/` to the relevant README section
- Use sequential numbering (Figure 1, Figure 2, ...)
- Include descriptive captions
- Use `?v=N` query parameter to bust image cache when updating

### Pre-commit Checklist

1. Update README.md if needed
2. Commit code changes
3. Include figures and data files in commit
4. **Push to remote branch（必須）**

## File Structure

- `src/pnp_solver.cpp`: Solver implementation (~2700 lines)
- `src/main.cpp`: CLI entry point with argument parsing
- `include/pnp_solver.hpp`: Class definition and parameters struct
- `scripts/plot_*.py`: Visualization scripts
- `styles/`: Plot style utilities (see below)
- `results/`: Output data (.dat) and figures (.png)
- `docs/`: Reference papers (PDF)

## Notation Convention

**Use number density n instead of concentration c in documentation and plots:**
- Axis labels: `$n/n_0$` (not `$c/c_0$`)
- Legend: `$n_+$`, `$n_-$` (not `$c_+$`, `$c_-$`)
- README equations use n [m⁻³], with note: n = c × NA × 10³

**Avoid underscores outside math mode in Markdown:**
- Use Unicode subscripts: λD, kB, εᵣ, NA, φT (not λ_D, k_B, ε_r)
- Or wrap in math mode: `$\lambda_D$`

## Plot Style Guide

Standard plot style is defined in `styles/plot_style.py`. Use this for all figures.

```python
from styles.plot_style import setup_plot_style, setup_axis_style, set_labels, FIGURE_SIZES

setup_plot_style()  # Set global matplotlib defaults

fig, ax = plt.subplots(figsize=FIGURE_SIZES['single'])  # 3.7 x 3.7 inches

# ... plot data ...

set_labels(ax, r'$x$ (units)', r'$y$ (units)')
setup_axis_style(ax)  # Ticks: direction='in', both sides, minor ticks
plt.tight_layout()
```

**Key style elements:**
- Figure size: 3.7 x 3.7 inches (square, single panel)
- Font: 10pt, bold for axis labels
- Ticks: direction='in', on both sides, major length=6, minor length=3
- Minor ticks: AutoMinorLocator(5) - do NOT use on log scales
- Colors: viridis colormap (`get_viridis_colors(n)`)
- Line width: 1
- No legend by default (add only if needed)
- All plots must be square format
