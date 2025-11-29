# 1D Poisson-Nernst-Planck Solver for Ionic Liquids

A C++ implementation of the Poisson-Nernst-Planck (PNP) equations for simulating the electric double layer (EDL) in ionic liquids.

## Features

- Newton-Raphson solver for the nonlinear Poisson-Boltzmann equation
- Non-uniform grid with clustering near the interface
- Gouy-Chapman analytical solution for validation
- **Bikerman model** for steric effects (ion size exclusion)
- **Transient solver** using implicit Scharfetter-Gummel scheme
- L2/L-infinity error analysis
- Python visualization scripts

## Theory

### Governing Equations

The Poisson-Nernst-Planck equations describe ion transport in electrolytes:

**Poisson Equation** (electrostatic potential):
```
∇²φ = -ρ/ε = -(e/ε)(z₊c₊ + z₋c₋)
```

**Nernst-Planck Equation** (ion flux):
```
J_i = -D_i(∇c_i + (z_i e c_i)/(k_B T) ∇φ)
```

**Continuity Equation**:
```
∂c_i/∂t = -∇·J_i
```

### 1D Formulation

For a 1D domain (x-direction only):

```
d²φ/dx² = -(e/ε)(z₊c₊ + z₋c₋)

∂c_±/∂t = D_± ∂/∂x (∂c_±/∂x ± (e c_±)/(k_B T) ∂φ/∂x)
```

### Steady-State Solution

At equilibrium with zero flux (J = 0), the ion concentrations follow Boltzmann distribution:

```
c_± = c₀ exp(∓eφ/(k_B T))
```

Substituting into Poisson's equation yields the **Poisson-Boltzmann equation**:

```
d²φ/dx² = (2 e N_A c₀/ε) sinh(eφ/(k_B T))
```

### Non-Dimensionalization

**Characteristic scales**:
- Length: Debye length λ_D = √(ε k_B T / (2 e² c₀ N_A))
- Potential: Thermal voltage φ_T = k_B T / e ≈ 25.7 mV at 298 K
- Concentration: Bulk concentration c₀

**Non-dimensional variables**:
- ξ = x / λ_D
- ψ = φ / φ_T = eφ / (k_B T)

**Non-dimensional Poisson-Boltzmann equation**:
```
d²ψ/dξ² = sinh(ψ)
```

### Gouy-Chapman Analytical Solution

For a 1:1 electrolyte with surface potential ψ₀, the analytical solution is:

```
tanh(ψ/4) = tanh(ψ₀/4) exp(-ξ)
```

Or equivalently:
```
φ(x) = (4 k_B T / e) arctanh[tanh(eφ₀/(4k_B T)) exp(-x/λ_D)]
```

### Bikerman Model (Modified Poisson-Boltzmann)

The standard Poisson-Boltzmann equation predicts unphysically high concentrations near charged surfaces. The Bikerman model introduces steric effects by accounting for finite ion size.

**Modified Boltzmann Distribution**:
```
c_± = c₀ exp(∓ψ) / g(ψ)
```

where the crowding function is:
```
g(ψ) = 1 - ν + ν cosh(ψ)
```

**Packing Fraction**:
```
ν = 2 a³ c₀ N_A
```

where `a` is the ion diameter (typically 0.5-1.0 nm for ionic liquids).

**Modified Poisson-Boltzmann Equation**:
```
d²ψ/dξ² = sinh(ψ) / g(ψ)
```

**Physical Interpretation**:
- At low potentials (|ψ| << 1): g(ψ) ≈ 1, reduces to standard PB
- At high potentials: g(ψ) limits the maximum concentration to ~1/ν × c₀
- Prevents unphysical crowding where c > 1/(a³ N_A)

**Reference**: Kilic, Bazant & Ajdari, *Phys. Rev. E* 75, 021502 (2007)

## Numerical Method

### Discretization

**Non-uniform Grid**:

The grid is stretched to cluster points near the interface (x = 0):
```
x_i = L × [1 - (1 - ξ_i)^β]
```
where ξ_i = i/(N-1) ∈ [0,1] and β > 1 is the stretching factor.

**Second Derivative (Non-uniform Grid)**:

For interior points:
```
d²φ/dx² ≈ [φ_{i+1} - φ_i)/Δx⁺ - (φ_i - φ_{i-1})/Δx⁻] / Δx_avg
```
where:
- Δx⁺ = x_{i+1} - x_i
- Δx⁻ = x_i - x_{i-1}
- Δx_avg = (Δx⁺ + Δx⁻) / 2

### Newton-Raphson Method

The nonlinear Poisson-Boltzmann equation is solved using Newton-Raphson iteration:

**Residual**:
```
F(φ) = d²φ/dx² - κ² φ_T sinh(φ/φ_T) = 0
```

**Jacobian**:
```
J = dF/dφ = d²/dx² - κ² cosh(φ/φ_T)
```

**Newton Update**:
```
J · δφ = -F
φ^{n+1} = φ^n + α · δφ
```
where α ∈ (0, 1] is an adaptive damping factor for stability.

### Algorithm (Steady-State)

1. Initialize with Gouy-Chapman analytical solution
2. Build Jacobian matrix and residual vector
3. Solve tridiagonal system using Thomas algorithm
4. Update solution with adaptive damping
5. Check convergence (relative change < tolerance)
6. Repeat until converged

### Transient Solver

The transient solver uses an implicit scheme for stability:

**Time Discretization**:
```
(c_i^{n+1} - c_i^n) / Δt = D ∂/∂x [∂c/∂x + z e c/(kT) ∂φ/∂x]
```

**Scharfetter-Gummel Scheme**:

For drift-diffusion problems, the Scharfetter-Gummel scheme provides stable discretization:
```
J_{i+1/2} = D/Δx [B(v Δx) c_{i+1} - B(-v Δx) c_i]
```

where `B(x) = x/(exp(x)-1)` is the Bernoulli function and `v = zeE/(kT)`.

**Algorithm at Each Time Step**:
1. Solve quasi-static Poisson equation for φ
2. Update c₊ using implicit Nernst-Planck with Scharfetter-Gummel
3. Update c₋ using implicit Nernst-Planck with Scharfetter-Gummel
4. Check for steady state (Δc/c₀ < tolerance)

**Characteristic Time Scale**:
```
τ_D = λ_D² / D ≈ 0.1-1 ns (for typical ionic liquids)
```

## Building and Running

### Requirements

- C++17 compatible compiler (g++ recommended)
- Python 3 with NumPy and Matplotlib (for visualization)

### Build

```bash
make
```

### Run

```bash
# Default parameters (1 M, 100 mV, steady-state)
make run

# Standard Poisson-Boltzmann
./build/pnp_solver --c0 0.1 --phi0 100 --output results/standard.dat

# Bikerman model with steric effects
./build/pnp_solver --c0 0.1 --phi0 100 --model bikerman --ion-size 0.7 --output results/bikerman.dat

# Transient simulation
./build/pnp_solver --c0 0.1 --phi0 100 --transient --dt 0.1 --t-final 1.0 --output results/transient.dat

# Available options:
#   --phi0 <value>      Surface potential in mV (default: 100)
#   --c0 <value>        Bulk concentration in mol/L (default: 1.0)
#   --eps <value>       Relative permittivity (default: 12)
#   --L <value>         Domain length in nm (default: 100)
#   --N <value>         Number of grid points (default: 1001)
#   --model <type>      Model: standard or bikerman (default: standard)
#   --ion-size <value>  Ion diameter in nm for Bikerman (default: 0.7)
#   --transient         Run transient simulation
#   --dt <value>        Time step in ns (default: 0.1)
#   --t-final <value>   Final time in microseconds (default: 1.0)
```

### Visualize

```bash
python3 scripts/plot_results.py
```

## Results

### Error Metrics

The solver computes the following error metrics against the Gouy-Chapman analytical solution:

**L2 Error (RMS)**:
```
L2 = √(Σ(φ_num - φ_GC)² / N)
```

**Relative L2 Error**:
```
L2_rel = L2 / √(Σφ_GC² / N)
```

**L-infinity Error (Maximum)**:
```
L∞ = max|φ_num - φ_GC|
```

### Validation with Gouy-Chapman Theory

The numerical solution is validated against the analytical Gouy-Chapman solution for dilute electrolytes.

**Test Case**: c₀ = 0.1 M, φ₀ = 100 mV, ε_r = 12

| Parameter | Value |
|-----------|-------|
| Debye length | 0.376 nm |
| Thermal voltage | 25.7 mV |
| Normalized potential | 3.89 |
| L2 error | ~1.5 mV |
| Relative L2 error | ~3.5% |
| L∞ (max) error | ~4 mV |
| Convergence | 4 iterations |

### Electric Double Layer Structure

![Combined Results](results/combined_results.png)

*Figure: (a) Electric potential profile comparing numerical PNP solution with Gouy-Chapman theory. (b) Normalized potential. (c) Ion concentration profiles showing cation depletion and anion accumulation near the positively charged surface. (d) Space charge density.*

### Key Observations

1. **Potential decay**: The potential decays exponentially from the surface with characteristic length λ_D
2. **Ion distribution**: Counterions (anions for positive surface) are enriched, co-ions (cations) are depleted
3. **Electroneutrality**: The bulk region (x >> λ_D) is electroneutral with c₊ = c₋ = c₀
4. **Boltzmann statistics**: Ion concentrations follow Boltzmann distribution in equilibrium

## File Structure

```
pnp/
├── include/
│   └── pnp_solver.hpp      # Header file
├── src/
│   ├── pnp_solver.cpp      # Solver implementation
│   └── main.cpp            # Main program
├── scripts/
│   ├── plot_results.py     # Visualization
│   └── plot_parametric.py  # Parametric study plots
├── results/                 # Output data and figures
├── Makefile
└── README.md
```

## Future Work

1. **Carnahan-Starling model**: More accurate equation of state for hard spheres

2. **2D/3D extension**: Extend to higher dimensions for more complex geometries

3. **Applied voltage**: Non-equilibrium steady states with current flow

4. **Asymmetric ions**: Different sizes for cations and anions

5. **Multiple species**: Extension to multi-component electrolytes

## References

1. Newman, J., & Thomas-Alyea, K. E. (2004). *Electrochemical Systems* (3rd ed.). Wiley.

2. Bazant, M. Z., Kilic, M. S., Storey, B. D., & Ajdari, A. (2009). Towards an understanding of induced-charge electrokinetics at large applied voltages in concentrated solutions. *Advances in Colloid and Interface Science*, 152(1-2), 48-88.

3. Kilic, M. S., Bazant, M. Z., & Ajdari, A. (2007). Steric effects in the dynamics of electrolytes at large applied voltages. *Physical Review E*, 75(2), 021502.

4. Bazant, M. Z., Storey, B. D., & Kornyshev, A. A. (2011). Double layer in ionic liquids: Overscreening versus crowding. *Physical Review Letters*, 106(4), 046102.

## License

MIT License
