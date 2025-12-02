/**
 * @file pnp_solver.cpp
 * @brief Implementation of 1D Poisson-Nernst-Planck solver
 *
 * Solves the Poisson-Boltzmann equation using Newton-Raphson iteration
 * on a non-uniform grid with clustering near the interface (x=0).
 *
 * Grid stretching uses: x = L * (1 - (1-ξ)^β) where ξ ∈ [0,1], β > 1
 * This gives finer spacing near x=0 and coarser spacing near x=L.
 */

#include "pnp_solver.hpp"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <limits>
#include <sstream>

namespace pnp {

PNPSolver1D::PNPSolver1D(const PNPParameters& params)
    : params_(params)
{
    beta_ = 1.0 / (PhysicalConstants::kB * params_.T);
    eps_ = params_.eps_r * PhysicalConstants::eps0;
    phi_T_ = PhysicalConstants::kB * params_.T / PhysicalConstants::e;

    // Compute Debye length: λ_D = sqrt(ε*kB*T / (2*e²*c0*NA))
    double c0_particles = params_.c0 * PhysicalConstants::NA;
    lambda_D_ = std::sqrt(eps_ * PhysicalConstants::kB * params_.T /
                          (2.0 * PhysicalConstants::e * PhysicalConstants::e * c0_particles));

    // Compute packing fraction for Bikerman model: ν = 2 * a³ * c0 * NA
    // a is ion diameter, factor of 2 for both + and - ions
    double a3 = params_.a * params_.a * params_.a;
    nu_ = 2.0 * a3 * params_.c0 * PhysicalConstants::NA;

    std::cout << "========================================\n";
    std::cout << "  1D PNP Solver Initialization\n";
    std::cout << "========================================\n";
    std::cout << "  Model: " << get_model_name() << "\n";
    std::cout << "  Debye length: " << lambda_D_ * 1e9 << " nm\n";
    std::cout << "  Thermal voltage: " << phi_T_ * 1000.0 << " mV\n";
    if (params_.model == ModelType::BIKERMAN) {
        std::cout << "  Ion diameter: " << params_.a * 1e9 << " nm\n";
        std::cout << "  Packing fraction: " << nu_ << "\n";
    }
    std::cout << "  Grid points: " << params_.N << "\n";
    std::cout << "  Grid stretching: " << params_.grid_stretch << "\n";

    initialize();
}

std::string PNPSolver1D::get_model_name() const {
    switch (params_.model) {
        case ModelType::STANDARD_PB: return "Standard Poisson-Boltzmann";
        case ModelType::BIKERMAN: return "Modified PB (Bikerman steric model)";
        default: return "Unknown";
    }
}

void PNPSolver1D::create_nonuniform_grid() {
    // Create non-uniform grid
    // For single electrode (open system): cluster near x=0 only
    // For dual electrode (closed system or phi_right != 0): cluster near both ends

    int n = params_.N;
    double L = params_.L;
    double beta = params_.grid_stretch;

    x_.resize(n);
    dx_.resize(n);

    // Check if we need symmetric grid (dual electrode model)
    bool symmetric_grid = params_.closed_system ||
                          (std::abs(params_.phi_right) > 1e-10);

    if (symmetric_grid && beta > 1.0) {
        // Symmetric grid: cluster near both x=0 and x=L
        // Using hyperbolic tangent stretching for symmetric clustering
        // x = L/2 * [1 + tanh(α*(2ξ-1)) / tanh(α)]
        // where α controls the stretching (larger = more clustering at ends)
        double alpha = std::log(beta + std::sqrt(beta * beta - 1));  // Convert beta to alpha
        if (alpha < 0.5) alpha = 0.5;  // Minimum stretching

        for (int i = 0; i < n; ++i) {
            double xi = static_cast<double>(i) / (n - 1);  // ξ ∈ [0, 1]
            // Symmetric stretching using tanh
            double eta = 2.0 * xi - 1.0;  // η ∈ [-1, 1]
            x_[i] = L * 0.5 * (1.0 + std::tanh(alpha * eta) / std::tanh(alpha));
        }
        std::cout << "  Grid type: symmetric (dual-electrode)\n";
    } else {
        // Original single-sided clustering near x=0
        // Using transformation: x = L * [1 - (1 - ξ)^β] where ξ ∈ [0, 1]
        // β > 1 gives finer spacing near x=0
        for (int i = 0; i < n; ++i) {
            double xi = static_cast<double>(i) / (n - 1);  // ξ ∈ [0, 1]
            x_[i] = L * (1.0 - std::pow(1.0 - xi, beta));
        }
        std::cout << "  Grid type: single-sided (left electrode)\n";
    }

    // Compute local grid spacing
    dx_[0] = x_[1] - x_[0];
    for (int i = 1; i < n - 1; ++i) {
        dx_[i] = 0.5 * (x_[i + 1] - x_[i - 1]);  // Central difference spacing
    }
    dx_[n - 1] = x_[n - 1] - x_[n - 2];

    // Print grid statistics
    double dx_min = *std::min_element(dx_.begin(), dx_.end());
    double dx_max = *std::max_element(dx_.begin(), dx_.end());

    std::cout << "  Grid spacing range: " << dx_min * 1e9 << " - " << dx_max * 1e9 << " nm\n";
    std::cout << "  Min spacing / Debye length: " << dx_min / lambda_D_ << "\n";

    // Count points within 5 Debye lengths of each boundary
    int points_in_left_edl = 0;
    int points_in_right_edl = 0;
    for (int i = 0; i < n; ++i) {
        if (x_[i] < 5.0 * lambda_D_) points_in_left_edl++;
        if (x_[i] > L - 5.0 * lambda_D_) points_in_right_edl++;
    }
    std::cout << "  Points within 5*λ_D of left: " << points_in_left_edl << "\n";
    if (symmetric_grid) {
        std::cout << "  Points within 5*λ_D of right: " << points_in_right_edl << "\n";
    }
}

void PNPSolver1D::initialize() {
    create_nonuniform_grid();

    int n = params_.N;
    phi_.resize(n);
    c_plus_.resize(n);
    c_minus_.resize(n);

    // Initialize with Gouy-Chapman analytical solution
    double gamma0 = std::tanh(params_.phi_left / (4.0 * phi_T_));

    for (int i = 0; i < n; ++i) {
        double gamma = gamma0 * std::exp(-x_[i] / lambda_D_);
        gamma = std::clamp(gamma, -0.9999, 0.9999);
        phi_[i] = 4.0 * phi_T_ * std::atanh(gamma);
    }

    // Enforce boundary conditions
    phi_[0] = params_.phi_left;
    phi_[n - 1] = params_.phi_right;

    update_concentrations_from_phi();
}

void PNPSolver1D::update_concentrations_from_phi() {
    for (int i = 0; i < params_.N; ++i) {
        double phi_norm = std::clamp(phi_[i] / phi_T_, -50.0, 50.0);

        if (params_.model == ModelType::BIKERMAN) {
            // Bikerman model: c± = c₀ * exp(∓ψ) / g(ψ)
            // where g(ψ) = 1 - ν + ν*cosh(ψ)
            double cosh_val = std::cosh(phi_norm);
            double g = 1.0 - nu_ + nu_ * cosh_val;
            g = std::max(g, 0.01);

            c_plus_[i] = params_.c0 * std::exp(-params_.z_plus * phi_norm) / g;
            c_minus_[i] = params_.c0 * std::exp(-params_.z_minus * phi_norm) / g;
        } else {
            // Standard Boltzmann distribution
            c_plus_[i] = params_.c0 * std::exp(-params_.z_plus * phi_norm);
            c_minus_[i] = params_.c0 * std::exp(-params_.z_minus * phi_norm);
        }
    }
}

double PNPSolver1D::bernoulli(double x) {
    if (std::abs(x) < 1e-6) {
        return 1.0 - 0.5 * x + x * x / 12.0;
    } else if (x > 50.0) {
        return x * std::exp(-x);
    } else if (x < -50.0) {
        return -x;
    } else {
        return x / (std::exp(x) - 1.0);
    }
}

std::vector<double> PNPSolver1D::solve_tridiagonal(
    const std::vector<double>& a,
    const std::vector<double>& b,
    const std::vector<double>& c,
    const std::vector<double>& d
) const {
    int n = b.size();
    std::vector<double> c_prime(n), d_prime(n), x(n);

    c_prime[0] = c[0] / b[0];
    d_prime[0] = d[0] / b[0];

    for (int i = 1; i < n; ++i) {
        double denom = b[i] - a[i] * c_prime[i - 1];
        if (std::abs(denom) < 1e-30) {
            denom = (denom >= 0) ? 1e-30 : -1e-30;
        }
        c_prime[i] = c[i] / denom;
        d_prime[i] = (d[i] - a[i] * d_prime[i - 1]) / denom;
    }

    x[n - 1] = d_prime[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        x[i] = d_prime[i] - c_prime[i] * x[i + 1];
    }

    return x;
}

void PNPSolver1D::solve_poisson() {
    // Placeholder - actual solving is in solve()
}

void PNPSolver1D::update_concentrations() {
    update_concentrations_from_phi();
}

double PNPSolver1D::compute_residual() const {
    double kappa2 = 1.0 / (lambda_D_ * lambda_D_);
    double max_residual = 0.0;

    for (int i = 1; i < params_.N - 1; ++i) {
        // Non-uniform grid second derivative
        double dx_plus = x_[i + 1] - x_[i];
        double dx_minus = x_[i] - x_[i - 1];
        double dx_avg = 0.5 * (dx_plus + dx_minus);

        double laplacian = (phi_[i + 1] - phi_[i]) / dx_plus - (phi_[i] - phi_[i - 1]) / dx_minus;
        laplacian /= dx_avg;

        double phi_norm = std::clamp(phi_[i] / phi_T_, -50.0, 50.0);
        double rhs = kappa2 * phi_T_ * std::sinh(phi_norm);

        double residual = std::abs(laplacian - rhs);
        max_residual = std::max(max_residual, residual);
    }

    return max_residual;
}

bool PNPSolver1D::solve() {
    std::cout << "\nStarting Newton-Raphson solver...\n";
    std::cout << "  Normalized surface potential: " << params_.phi_left / phi_T_ << "\n";

    double kappa2 = 1.0 / (lambda_D_ * lambda_D_);
    int n = params_.N;

    residual_history_.clear();

    for (int iter = 0; iter < params_.max_iter; ++iter) {
        std::vector<double> a(n, 0.0), b(n, 0.0), c(n, 0.0), rhs(n, 0.0);

        // Left boundary (Dirichlet)
        b[0] = 1.0;
        rhs[0] = 0.0;

        // Interior points with non-uniform grid
        for (int i = 1; i < n - 1; ++i) {
            double dx_plus = x_[i + 1] - x_[i];
            double dx_minus = x_[i] - x_[i - 1];
            double dx_avg = 0.5 * (dx_plus + dx_minus);

            double phi_norm = std::clamp(phi_[i] / phi_T_, -50.0, 50.0);
            double sinh_val = std::sinh(phi_norm);
            double cosh_val = std::cosh(phi_norm);

            // Laplacian coefficient for non-uniform grid
            double lap_minus = 1.0 / (dx_minus * dx_avg);  // coeff for φ[i-1]
            double lap_plus = 1.0 / (dx_plus * dx_avg);    // coeff for φ[i+1]
            double lap_center = -lap_minus - lap_plus;      // coeff for φ[i]

            // Laplacian value
            double laplacian = (phi_[i + 1] - phi_[i]) / dx_plus - (phi_[i] - phi_[i - 1]) / dx_minus;
            laplacian /= dx_avg;

            double source_term, jacobian_term;

            if (params_.model == ModelType::BIKERMAN) {
                // Bikerman model with steric effects
                // F(φ) = d²φ/dx² - κ²*φ_T*sinh(ψ) / g(ψ) = 0
                // where g(ψ) = 1 - ν + ν*cosh(ψ)
                double g = 1.0 - nu_ + nu_ * cosh_val;
                g = std::max(g, 0.01);  // Prevent division by zero

                source_term = kappa2 * phi_T_ * sinh_val / g;

                // Jacobian: dF/dφ = d²/dx² - κ² * (cosh(ψ) - ν*(cosh(ψ) - 1)) / g²
                double dg_dpsi = nu_ * sinh_val;
                double df_dpsi = (cosh_val * g - sinh_val * dg_dpsi) / (g * g);
                jacobian_term = kappa2 * df_dpsi;
            } else {
                // Standard Poisson-Boltzmann
                // F(φ) = d²φ/dx² - κ²*φ_T*sinh(ψ) = 0
                source_term = kappa2 * phi_T_ * sinh_val;
                jacobian_term = kappa2 * cosh_val;
            }

            // Jacobian matrix entries
            a[i] = lap_minus;
            c[i] = lap_plus;
            b[i] = lap_center - jacobian_term;

            // Residual: -F
            rhs[i] = -(laplacian - source_term);
        }

        // Right boundary (Dirichlet)
        b[n - 1] = 1.0;
        rhs[n - 1] = 0.0;

        // Solve for Newton correction
        std::vector<double> delta_phi = solve_tridiagonal(a, b, c, rhs);

        // Update with adaptive damping
        double max_delta = 0.0;
        for (int i = 1; i < n - 1; ++i) {
            max_delta = std::max(max_delta, std::abs(delta_phi[i]));
        }

        double damping = 1.0;
        if (max_delta > phi_T_) {
            damping = phi_T_ / max_delta;
        }

        for (int i = 1; i < n - 1; ++i) {
            phi_[i] += damping * delta_phi[i];
        }

        double residual = compute_residual();
        residual_history_.push_back(residual);

        double rel_change = max_delta / std::max(std::abs(params_.phi_left), phi_T_);

        if (iter < 10 || iter % 10 == 0) {
            std::cout << "  Iter " << iter << ": rel_change = " << rel_change
                      << ", residual = " << residual << "\n";
        }

        if (rel_change < params_.tol && iter > 0) {
            std::cout << "  Converged after " << iter + 1 << " iterations.\n";
            update_concentrations_from_phi();
            return true;
        }
    }

    std::cout << "  Warning: Did not converge.\n";
    update_concentrations_from_phi();
    return false;
}

void PNPSolver1D::solve_transient(double dt, double t_final) {
    /**
     * Transient PNP solver using implicit time stepping
     *
     * Equations:
     *   ∂c±/∂t = D± ∂/∂x [∂c±/∂x ± (e c±)/(kT) ∂φ/∂x]
     *   ∇²φ = -ρ/ε  (quasi-static Poisson)
     *
     * Initial conditions:
     *   - Uniform bulk concentration throughout domain
     *   - Linear potential profile from surface to bulk
     *
     * Algorithm at each time step:
     *   1. Solve Poisson equation for φ given current c±
     *   2. Update c± using implicit Nernst-Planck discretization
     *   3. Repeat until steady state or t_final reached
     */

    int n_steps = static_cast<int>(t_final / dt);
    int n = params_.N;
    double e = PhysicalConstants::e;
    double kT = PhysicalConstants::kB * params_.T;

    std::cout << "\n========================================\n";
    std::cout << "  Starting Transient Solver\n";
    std::cout << "========================================\n";
    std::cout << "  Time step: " << dt * 1e9 << " ns\n";
    std::cout << "  Total time: " << t_final * 1e6 << " µs\n";
    std::cout << "  Number of steps: " << n_steps << "\n";

    // Set initial conditions for transient analysis
    // Physical initial state: uniform concentration, linear potential
    std::cout << "  Initial conditions: uniform c0, linear φ profile\n";

    for (int i = 0; i < n; ++i) {
        c_plus_[i] = params_.c0;
        c_minus_[i] = params_.c0;
        double xi = x_[i] / params_.L;
        phi_[i] = params_.phi_left * (1.0 - xi) + params_.phi_right * xi;
    }
    phi_[0] = params_.phi_left;
    phi_[n - 1] = params_.phi_right;

    // Store previous concentrations for convergence check
    std::vector<double> c_plus_old(n), c_minus_old(n);

    for (int step = 0; step < n_steps; ++step) {
        // Store old values
        c_plus_old = c_plus_;
        c_minus_old = c_minus_;

        // Step 1: Solve Poisson equation (quasi-static assumption)
        // ∇²φ = -ρ/ε with BCs: φ(0) = φ_left, φ(L) = φ_right = 0
        // Direct solve: A * φ = b
        {
            std::vector<double> a(n, 0.0), b(n, 0.0), c(n, 0.0), rhs(n, 0.0);

            // Left boundary: φ(0) = φ_left (Dirichlet)
            b[0] = 1.0;
            rhs[0] = params_.phi_left;

            for (int i = 1; i < n - 1; ++i) {
                double dx_plus = x_[i + 1] - x_[i];
                double dx_minus = x_[i] - x_[i - 1];
                double dx_avg = 0.5 * (dx_plus + dx_minus);

                // Charge density from current concentrations
                double rho_i = PhysicalConstants::NA * e *
                               (params_.z_plus * c_plus_[i] + params_.z_minus * c_minus_[i]);

                // Laplacian coefficients: d²φ/dx² ≈ (φ_{i+1} - φ_i)/dx+ - (φ_i - φ_{i-1})/dx-) / dx_avg
                double lap_minus = 1.0 / (dx_minus * dx_avg);
                double lap_plus = 1.0 / (dx_plus * dx_avg);
                double lap_center = -lap_minus - lap_plus;

                // A * φ = b, where A is the Laplacian matrix and b = -ρ/ε
                a[i] = lap_minus;
                c[i] = lap_plus;
                b[i] = lap_center;
                rhs[i] = -rho_i / eps_;
            }

            // Right boundary: φ(L) = φ_right = 0 (Dirichlet)
            b[n - 1] = 1.0;
            rhs[n - 1] = params_.phi_right;

            phi_ = solve_tridiagonal(a, b, c, rhs);
        }

        // Step 2: Update concentrations using implicit Nernst-Planck
        // Using Scharfetter-Gummel discretization for stability
        for (int species = 0; species < 2; ++species) {
            double D = (species == 0) ? params_.D_plus : params_.D_minus;
            int z = (species == 0) ? params_.z_plus : params_.z_minus;
            std::vector<double>& conc = (species == 0) ? c_plus_ : c_minus_;
            const std::vector<double>& conc_old = (species == 0) ? c_plus_old : c_minus_old;

            std::vector<double> a(n, 0.0), b(n, 0.0), c_vec(n, 0.0), rhs(n, 0.0);

            // Left boundary (x=0): Zero flux (blocking electrode)
            double dx0 = x_[1] - x_[0];
            double v0 = -z * e / kT * (phi_[1] - phi_[0]) / dx0;
            double B0_p = bernoulli(v0 * dx0);
            double B0_m = bernoulli(-v0 * dx0);

            double alpha0 = D * dt / (dx0 * dx0);
            b[0] = 1.0 + alpha0 * B0_m;
            c_vec[0] = -alpha0 * B0_p;
            rhs[0] = conc_old[0];

            for (int i = 1; i < n - 1; ++i) {
                double dx_plus = x_[i + 1] - x_[i];
                double dx_minus = x_[i] - x_[i - 1];
                double dx_avg = 0.5 * (dx_plus + dx_minus);

                // Electric field and drift velocity
                double v_plus = -z * e / kT * (phi_[i + 1] - phi_[i]) / dx_plus;
                double v_minus = -z * e / kT * (phi_[i] - phi_[i - 1]) / dx_minus;

                double B_plus_p = bernoulli(v_plus * dx_plus);
                double B_plus_m = bernoulli(-v_plus * dx_plus);
                double B_minus_p = bernoulli(v_minus * dx_minus);
                double B_minus_m = bernoulli(-v_minus * dx_minus);

                double alpha_plus = D * dt / (dx_plus * dx_avg);
                double alpha_minus = D * dt / (dx_minus * dx_avg);

                a[i] = -alpha_minus * B_minus_p;
                c_vec[i] = -alpha_plus * B_plus_p;
                b[i] = 1.0 + alpha_minus * B_minus_m + alpha_plus * B_plus_m;
                rhs[i] = conc_old[i];
            }

            // Right boundary (x=L): Dirichlet (bulk concentration)
            b[n - 1] = 1.0;
            rhs[n - 1] = params_.c0;

            conc = solve_tridiagonal(a, b, c_vec, rhs);

            // Ensure positivity
            for (int i = 0; i < n; ++i) {
                conc[i] = std::max(conc[i], 1e-10 * params_.c0);
            }
        }

        // Print progress
        if (step % std::max(1, n_steps / 20) == 0 || step == n_steps - 1) {
            double max_change_plus = 0.0, max_change_minus = 0.0;
            for (int i = 0; i < n; ++i) {
                max_change_plus = std::max(max_change_plus,
                    std::abs(c_plus_[i] - c_plus_old[i]) / params_.c0);
                max_change_minus = std::max(max_change_minus,
                    std::abs(c_minus_[i] - c_minus_old[i]) / params_.c0);
            }

            double time_us = (step + 1) * dt * 1e6;
            std::cout << "  t = " << std::fixed << std::setprecision(3) << time_us
                      << " µs: Δc+/c0 = " << std::scientific << std::setprecision(2)
                      << max_change_plus << ", Δc-/c0 = " << max_change_minus << "\n";

            // Check for steady state
            if (max_change_plus < 1e-8 && max_change_minus < 1e-8) {
                std::cout << "  Steady state reached at step " << step << "\n";
                break;
            }
        }
    }

    std::cout << "  Transient simulation completed.\n";
}

void PNPSolver1D::solve_transient_with_snapshots(double dt, double t_final,
                                                  const std::string& output_dir,
                                                  int snapshot_interval) {
    /**
     * Transient PNP solver with snapshot output for animation
     *
     * Physical setup:
     *   - Initial condition: uniform bulk concentration, zero potential
     *   - At t=0, surface potential is applied (step function)
     *   - System evolves toward equilibrium EDL structure
     *
     * Numerical scheme:
     *   - Implicit time stepping for stability
     *   - Scharfetter-Gummel scheme for drift-diffusion flux
     *   - Quasi-static Poisson solver at each time step
     *
     * Boundary conditions:
     *   - Left (x=0): Fixed potential (Dirichlet), zero ion flux (blocking electrode)
     *   - Right (x=L): Zero potential, bulk concentration (Dirichlet)
     *
     * References:
     *   - Bazant et al., Phys. Rev. E 70, 021506 (2004) - EDL charging dynamics
     *   - Scharfetter & Gummel, IEEE Trans. ED 16, 64 (1969) - SG scheme
     */

    int n_steps = static_cast<int>(t_final / dt);
    int n = params_.N;
    double e = PhysicalConstants::e;
    double kT = PhysicalConstants::kB * params_.T;

    std::cout << "\n========================================\n";
    std::cout << "  Starting Transient Solver (with snapshots)\n";
    std::cout << "========================================\n";
    std::cout << "  Time step: " << dt * 1e9 << " ns\n";
    std::cout << "  Total time: " << t_final * 1e6 << " µs\n";
    std::cout << "  Number of steps: " << n_steps << "\n";
    std::cout << "  Snapshot interval: " << snapshot_interval << " steps\n";
    std::cout << "  Output directory: " << output_dir << "\n";

    // ============================================
    // Set initial conditions for transient analysis
    // ============================================
    // Physical initial state: uniform concentration, zero potential everywhere
    // Then at t=0+, surface potential is applied as step function

    std::cout << "\n  Setting initial conditions:\n";
    std::cout << "    - Uniform concentration: c = c0 = " << params_.c0 << " mol/m³\n";
    std::cout << "    - Initial potential: φ = 0 (then step to φ0 at surface)\n";

    // Initialize with uniform bulk concentration
    for (int i = 0; i < n; ++i) {
        c_plus_[i] = params_.c0;
        c_minus_[i] = params_.c0;
    }

    // Initialize potential: linear profile from surface to bulk
    // This provides a smooth initial condition that respects BCs
    for (int i = 0; i < n; ++i) {
        // Linear interpolation: φ(x) = φ_left * (1 - x/L) + φ_right * (x/L)
        double xi = x_[i] / params_.L;
        phi_[i] = params_.phi_left * (1.0 - xi) + params_.phi_right * xi;
    }

    // Enforce boundary conditions exactly
    phi_[0] = params_.phi_left;
    phi_[n - 1] = params_.phi_right;

    // Save time info file
    std::ofstream time_file(output_dir + "/time_info.dat");
    time_file << "# Snapshot time information\n";
    time_file << "# snapshot_index  time_ns\n";

    // Store previous concentrations for convergence check
    std::vector<double> c_plus_old(n), c_minus_old(n);

    int snapshot_count = 0;
    bool steady_state_reached = false;

    // Save initial state (t=0)
    {
        std::ostringstream ss;
        ss << output_dir << "/snapshot_" << std::setfill('0') << std::setw(5) << snapshot_count << ".dat";
        save_results(ss.str());
        time_file << snapshot_count << "\t" << 0.0 << "\n";
        snapshot_count++;
    }

    for (int step = 0; step < n_steps; ++step) {
        // Store old values
        c_plus_old = c_plus_;
        c_minus_old = c_minus_;

        // Step 1: Solve Poisson equation (quasi-static assumption)
        // ∇²φ = -ρ/ε with BCs: φ(0) = φ_left, φ(L) = φ_right = 0
        // Direct solve: A * φ = b
        {
            std::vector<double> a(n, 0.0), b(n, 0.0), c(n, 0.0), rhs(n, 0.0);

            // Left boundary: φ(0) = φ_left (Dirichlet)
            b[0] = 1.0;
            rhs[0] = params_.phi_left;

            for (int i = 1; i < n - 1; ++i) {
                double dx_plus = x_[i + 1] - x_[i];
                double dx_minus = x_[i] - x_[i - 1];
                double dx_avg = 0.5 * (dx_plus + dx_minus);

                // Charge density from current concentrations
                double rho_i = PhysicalConstants::NA * e *
                               (params_.z_plus * c_plus_[i] + params_.z_minus * c_minus_[i]);

                // Laplacian coefficients
                double lap_minus = 1.0 / (dx_minus * dx_avg);
                double lap_plus = 1.0 / (dx_plus * dx_avg);
                double lap_center = -lap_minus - lap_plus;

                // A * φ = b, where A is the Laplacian matrix and b = -ρ/ε
                a[i] = lap_minus;
                c[i] = lap_plus;
                b[i] = lap_center;
                rhs[i] = -rho_i / eps_;
            }

            // Right boundary: φ(L) = φ_right = 0 (Dirichlet)
            b[n - 1] = 1.0;
            rhs[n - 1] = params_.phi_right;

            phi_ = solve_tridiagonal(a, b, c, rhs);
        }

        // Step 2: Update concentrations using implicit Nernst-Planck
        // Using Scharfetter-Gummel discretization for stability
        for (int species = 0; species < 2; ++species) {
            double D = (species == 0) ? params_.D_plus : params_.D_minus;
            int z = (species == 0) ? params_.z_plus : params_.z_minus;
            std::vector<double>& conc = (species == 0) ? c_plus_ : c_minus_;
            const std::vector<double>& conc_old = (species == 0) ? c_plus_old : c_minus_old;

            std::vector<double> a(n, 0.0), b(n, 0.0), c_vec(n, 0.0), rhs(n, 0.0);

            // Left boundary (x=0): Zero flux (blocking electrode)
            // dc/dx + z*e/(kT) * c * dφ/dx = 0
            // Use ghost point approach: flux at i=0.5 is zero
            // This means c[0] adjusts to maintain zero flux
            double dx0 = x_[1] - x_[0];
            double v0 = -z * e / kT * (phi_[1] - phi_[0]) / dx0;  // drift velocity / D
            double B0_p = bernoulli(v0 * dx0);
            double B0_m = bernoulli(-v0 * dx0);

            // Zero flux BC: J_{0.5} = 0 => B(v)*c[1] - B(-v)*c[0] = 0
            // For time stepping: (c - c_old)/dt = D/dx * (flux_in - flux_out)
            // At i=0: flux_in = 0 (wall), flux_out = J_{0.5}
            double alpha0 = D * dt / (dx0 * dx0);
            b[0] = 1.0 + alpha0 * B0_m;
            c_vec[0] = -alpha0 * B0_p;
            rhs[0] = conc_old[0];

            for (int i = 1; i < n - 1; ++i) {
                double dx_plus = x_[i + 1] - x_[i];
                double dx_minus = x_[i] - x_[i - 1];
                double dx_avg = 0.5 * (dx_plus + dx_minus);

                // Electric field and drift velocity
                double v_plus = -z * e / kT * (phi_[i + 1] - phi_[i]) / dx_plus;
                double v_minus = -z * e / kT * (phi_[i] - phi_[i - 1]) / dx_minus;

                double B_plus_p = bernoulli(v_plus * dx_plus);
                double B_plus_m = bernoulli(-v_plus * dx_plus);
                double B_minus_p = bernoulli(v_minus * dx_minus);
                double B_minus_m = bernoulli(-v_minus * dx_minus);

                // Flux: J_{i+1/2} = D/dx * [B(v)*c_{i+1} - B(-v)*c_i]
                double alpha_plus = D * dt / (dx_plus * dx_avg);
                double alpha_minus = D * dt / (dx_minus * dx_avg);

                a[i] = -alpha_minus * B_minus_p;
                c_vec[i] = -alpha_plus * B_plus_p;
                b[i] = 1.0 + alpha_minus * B_minus_m + alpha_plus * B_plus_m;
                rhs[i] = conc_old[i];
            }

            // Right boundary (x=L): Dirichlet (bulk concentration)
            b[n - 1] = 1.0;
            rhs[n - 1] = params_.c0;

            conc = solve_tridiagonal(a, b, c_vec, rhs);

            // Ensure positivity
            for (int i = 0; i < n; ++i) {
                conc[i] = std::max(conc[i], 1e-10 * params_.c0);
            }
        }

        // Save snapshot at specified intervals
        if ((step + 1) % snapshot_interval == 0) {
            std::ostringstream ss;
            ss << output_dir << "/snapshot_" << std::setfill('0') << std::setw(5) << snapshot_count << ".dat";
            save_results(ss.str());
            double time_ns = (step + 1) * dt * 1e9;
            time_file << snapshot_count << "\t" << time_ns << "\n";
            snapshot_count++;
        }

        // Check for steady state and print progress
        if (step % std::max(1, n_steps / 20) == 0 || step == n_steps - 1) {
            double max_change_plus = 0.0, max_change_minus = 0.0;
            for (int i = 0; i < n; ++i) {
                max_change_plus = std::max(max_change_plus,
                    std::abs(c_plus_[i] - c_plus_old[i]) / params_.c0);
                max_change_minus = std::max(max_change_minus,
                    std::abs(c_minus_[i] - c_minus_old[i]) / params_.c0);
            }

            double time_ns = (step + 1) * dt * 1e9;
            std::cout << "  t = " << std::fixed << std::setprecision(1) << time_ns
                      << " ns: Δc+/c0 = " << std::scientific << std::setprecision(2)
                      << max_change_plus << ", Δc-/c0 = " << max_change_minus << "\n";

            if (max_change_plus < 1e-8 && max_change_minus < 1e-8) {
                std::cout << "  Steady state reached at step " << step << "\n";
                steady_state_reached = true;

                // Save final steady state
                std::ostringstream ss;
                ss << output_dir << "/snapshot_" << std::setfill('0') << std::setw(5) << snapshot_count << ".dat";
                save_results(ss.str());
                double time_ns_final = (step + 1) * dt * 1e9;
                time_file << snapshot_count << "\t" << time_ns_final << "\n";
                snapshot_count++;
                break;
            }
        }
    }

    time_file.close();
    std::cout << "  Total snapshots saved: " << snapshot_count << "\n";
    std::cout << "  Transient simulation with snapshots completed.\n";
}

void PNPSolver1D::solve_transient_gummel(double dt, double t_final,
                                          const std::string& output_dir,
                                          int snapshot_interval) {
    /**
     * Transient PNP solver using Gummel iteration (fully implicit)
     *
     * At each time step, iterate between Poisson and NP equations until convergence.
     * This ensures proper coupling and unconditional stability.
     *
     * Algorithm:
     *   For each time step n -> n+1:
     *     1. Initialize: φ^(0) = φ^n, c^(0) = c^n
     *     2. Gummel loop (k = 0, 1, 2, ...):
     *        a. Solve Poisson: ∇²φ^(k+1) = -ρ(c^(k))/ε
     *        b. Solve NP for c+: (c+ - c+^n)/dt = ∇·(D∇c+ + D*z*e/kT * c+ * ∇φ^(k+1))
     *        c. Solve NP for c-: (c- - c-^n)/dt = ∇·(D∇c- + D*z*e/kT * c- * ∇φ^(k+1))
     *        d. Check convergence: ||c^(k+1) - c^(k)|| < tol
     *     3. Update: φ^(n+1) = φ^(k+1), c^(n+1) = c^(k+1)
     */

    int n_steps = static_cast<int>(t_final / dt);
    int n = params_.N;
    double e = PhysicalConstants::e;
    double kT = PhysicalConstants::kB * params_.T;

    const int max_gummel_iter = 200;
    const double gummel_tol = 1e-6;
    const double omega = 0.5;  // Under-relaxation factor for stability

    std::cout << "\n========================================\n";
    std::cout << "  Transient Solver (Gummel Iteration)\n";
    std::cout << "========================================\n";
    std::cout << "  Time step: " << dt * 1e9 << " ns\n";
    std::cout << "  Total time: " << t_final * 1e9 << " ns\n";
    std::cout << "  Number of steps: " << n_steps << "\n";
    std::cout << "  Snapshot interval: " << snapshot_interval << " steps\n";
    std::cout << "  Gummel tolerance: " << gummel_tol << "\n";

    // Initialize with uniform bulk concentration
    for (int i = 0; i < n; ++i) {
        c_plus_[i] = params_.c0;
        c_minus_[i] = params_.c0;
    }

    // Solve Poisson with initial uniform concentration to get starting potential
    // With ρ = 0, this gives linear profile
    {
        std::vector<double> a_p(n, 0.0), b_p(n, 0.0), c_p(n, 0.0), rhs_p(n, 0.0);
        b_p[0] = 1.0;
        rhs_p[0] = params_.phi_left;

        for (int i = 1; i < n - 1; ++i) {
            double dx_plus = x_[i + 1] - x_[i];
            double dx_minus = x_[i] - x_[i - 1];
            double dx_avg = 0.5 * (dx_plus + dx_minus);

            double rho_i = PhysicalConstants::NA * e *
                           (params_.z_plus * c_plus_[i] + params_.z_minus * c_minus_[i]);

            a_p[i] = 1.0 / (dx_minus * dx_avg);
            c_p[i] = 1.0 / (dx_plus * dx_avg);
            b_p[i] = -(a_p[i] + c_p[i]);
            rhs_p[i] = -rho_i / eps_;
        }

        b_p[n - 1] = 1.0;
        rhs_p[n - 1] = params_.phi_right;
        phi_ = solve_tridiagonal(a_p, b_p, c_p, rhs_p);
    }

    // Time info file
    std::ofstream time_file(output_dir + "/time_info.dat");
    time_file << "# snapshot_index  time_ns\n";

    // Store old time-step values
    std::vector<double> c_plus_n(n), c_minus_n(n), phi_n(n);
    // Store Gummel iteration values
    std::vector<double> c_plus_k(n), c_minus_k(n);

    int snapshot_count = 0;

    // Save initial state
    {
        std::ostringstream ss;
        ss << output_dir << "/snapshot_" << std::setfill('0') << std::setw(5) << snapshot_count << ".dat";
        save_results(ss.str());
        time_file << snapshot_count << "\t" << 0.0 << "\n";
        snapshot_count++;
    }

    for (int step = 0; step < n_steps; ++step) {
        // Store values at time n
        c_plus_n = c_plus_;
        c_minus_n = c_minus_;
        phi_n = phi_;

        // Gummel iteration
        int gummel_iter = 0;
        double gummel_error = 1.0;

        while (gummel_iter < max_gummel_iter && gummel_error > gummel_tol) {
            // Store current iteration values
            c_plus_k = c_plus_;
            c_minus_k = c_minus_;

            // Step 1: Solve Poisson equation with current concentrations
            {
                std::vector<double> a_p(n, 0.0), b_p(n, 0.0), c_p(n, 0.0), rhs_p(n, 0.0);

                b_p[0] = 1.0;
                rhs_p[0] = params_.phi_left;

                for (int i = 1; i < n - 1; ++i) {
                    double dx_plus = x_[i + 1] - x_[i];
                    double dx_minus = x_[i] - x_[i - 1];
                    double dx_avg = 0.5 * (dx_plus + dx_minus);

                    double rho_i = PhysicalConstants::NA * e *
                                   (params_.z_plus * c_plus_[i] + params_.z_minus * c_minus_[i]);

                    a_p[i] = 1.0 / (dx_minus * dx_avg);
                    c_p[i] = 1.0 / (dx_plus * dx_avg);
                    b_p[i] = -(a_p[i] + c_p[i]);
                    rhs_p[i] = -rho_i / eps_;
                }

                b_p[n - 1] = 1.0;
                rhs_p[n - 1] = params_.phi_right;
                phi_ = solve_tridiagonal(a_p, b_p, c_p, rhs_p);
            }

            // Step 2: Solve Nernst-Planck for each species (implicit in time)
            for (int species = 0; species < 2; ++species) {
                double D = (species == 0) ? params_.D_plus : params_.D_minus;
                int z = (species == 0) ? params_.z_plus : params_.z_minus;
                std::vector<double>& conc = (species == 0) ? c_plus_ : c_minus_;
                const std::vector<double>& conc_n = (species == 0) ? c_plus_n : c_minus_n;

                std::vector<double> a_c(n, 0.0), b_c(n, 0.0), c_c(n, 0.0), rhs_c(n, 0.0);

                // Left BC: zero flux (blocking electrode)
                // Use one-sided discretization for flux = 0
                double dx0 = x_[1] - x_[0];
                double v0 = -z * e / kT * (phi_[1] - phi_[0]) / dx0;
                double B0_p = bernoulli(v0 * dx0);
                double B0_m = bernoulli(-v0 * dx0);
                double alpha0 = D * dt / (dx0 * dx0);

                b_c[0] = 1.0 + alpha0 * B0_m;
                c_c[0] = -alpha0 * B0_p;
                rhs_c[0] = conc_n[0];

                // Interior points
                for (int i = 1; i < n - 1; ++i) {
                    double dx_plus = x_[i + 1] - x_[i];
                    double dx_minus = x_[i] - x_[i - 1];
                    double dx_avg = 0.5 * (dx_plus + dx_minus);

                    double v_plus = -z * e / kT * (phi_[i + 1] - phi_[i]) / dx_plus;
                    double v_minus = -z * e / kT * (phi_[i] - phi_[i - 1]) / dx_minus;

                    double B_plus_p = bernoulli(v_plus * dx_plus);
                    double B_plus_m = bernoulli(-v_plus * dx_plus);
                    double B_minus_p = bernoulli(v_minus * dx_minus);
                    double B_minus_m = bernoulli(-v_minus * dx_minus);

                    double alpha_plus = D * dt / (dx_plus * dx_avg);
                    double alpha_minus = D * dt / (dx_minus * dx_avg);

                    a_c[i] = -alpha_minus * B_minus_p;
                    c_c[i] = -alpha_plus * B_plus_p;
                    b_c[i] = 1.0 + alpha_minus * B_minus_m + alpha_plus * B_plus_m;
                    rhs_c[i] = conc_n[i];
                }

                // Right BC: depends on closed_system flag
                if (params_.closed_system) {
                    // Zero-flux BC (blocking electrode)
                    double dx_last = x_[n - 1] - x_[n - 2];
                    double v_last = -z * e / kT * (phi_[n - 1] - phi_[n - 2]) / dx_last;
                    double B_last_p = bernoulli(v_last * dx_last);
                    double B_last_m = bernoulli(-v_last * dx_last);
                    double alpha_last = D * dt / (dx_last * dx_last);

                    a_c[n - 1] = -alpha_last * B_last_p;
                    b_c[n - 1] = 1.0 + alpha_last * B_last_m;
                    c_c[n - 1] = 0.0;
                    rhs_c[n - 1] = conc_n[n - 1];
                } else {
                    // Dirichlet BC (bulk concentration)
                    b_c[n - 1] = 1.0;
                    rhs_c[n - 1] = params_.c0;
                }

                std::vector<double> conc_new = solve_tridiagonal(a_c, b_c, c_c, rhs_c);

                // Under-relaxation and positivity enforcement
                for (int i = 0; i < n; ++i) {
                    conc_new[i] = std::max(conc_new[i], 1e-12 * params_.c0);
                    conc[i] = omega * conc_new[i] + (1.0 - omega) * conc[i];
                }
            }

            // Compute Gummel convergence error
            gummel_error = 0.0;
            for (int i = 0; i < n; ++i) {
                double err_plus = std::abs(c_plus_[i] - c_plus_k[i]) / (params_.c0 + 1e-20);
                double err_minus = std::abs(c_minus_[i] - c_minus_k[i]) / (params_.c0 + 1e-20);
                gummel_error = std::max(gummel_error, std::max(err_plus, err_minus));
            }

            gummel_iter++;
        }

        // Save snapshot
        if ((step + 1) % snapshot_interval == 0) {
            std::ostringstream ss;
            ss << output_dir << "/snapshot_" << std::setfill('0') << std::setw(5) << snapshot_count << ".dat";
            save_results(ss.str());
            double time_ns = (step + 1) * dt * 1e9;
            time_file << snapshot_count << "\t" << time_ns << "\n";
            snapshot_count++;
        }

        // Progress output
        if (step % std::max(1, n_steps / 20) == 0 || step == n_steps - 1) {
            double time_ns = (step + 1) * dt * 1e9;
            double max_c_plus = *std::max_element(c_plus_.begin(), c_plus_.end());
            double min_c_plus = *std::min_element(c_plus_.begin(), c_plus_.end());
            std::cout << "  t = " << std::fixed << std::setprecision(2) << time_ns
                      << " ns, Gummel iter = " << gummel_iter
                      << ", c+/c0 range: [" << std::scientific << std::setprecision(2)
                      << min_c_plus / params_.c0 << ", " << max_c_plus / params_.c0 << "]\n";
        }
    }

    time_file.close();
    std::cout << "  Total snapshots: " << snapshot_count << "\n";
    std::cout << "  Transient (Gummel) completed.\n";
}

void PNPSolver1D::solve_transient_continuation(int n_steps,
                                                const std::string& output_dir) {
    /**
     * Pseudo-transient solver using continuation method
     *
     * Strategy: Start from zero potential and gradually increase to target.
     * At each step, solve steady-state PNP (which is stable).
     * This mimics the physical process of charging an EDL.
     *
     * The "time" axis represents the fraction of charging completed,
     * which is physically meaningful as the RC charging time.
     */

    int n = params_.N;
    double phi0_target = params_.phi_left;

    std::cout << "\n========================================\n";
    std::cout << "  Transient Solver (Continuation Method)\n";
    std::cout << "========================================\n";
    std::cout << "  Target potential: " << phi0_target * 1e3 << " mV\n";
    std::cout << "  Number of steps: " << n_steps << "\n";
    std::cout << "  Output directory: " << output_dir << "\n";

    // Time info file
    std::ofstream time_file(output_dir + "/time_info.dat");
    time_file << "# snapshot_index  time_ns  phi0_mV\n";

    // Estimate RC time constant: τ_RC = λ_D * L / D
    double tau_RC = lambda_D_ * params_.L / params_.D_plus;
    std::cout << "  Estimated RC time: " << tau_RC * 1e9 << " ns\n";

    // Initialize with uniform concentration
    for (int i = 0; i < n; ++i) {
        c_plus_[i] = params_.c0;
        c_minus_[i] = params_.c0;
        phi_[i] = 0.0;
    }

    int snapshot_count = 0;

    // Save initial state (phi0 = 0)
    params_.phi_left = 0.0;
    {
        std::ostringstream ss;
        ss << output_dir << "/snapshot_" << std::setfill('0') << std::setw(5) << snapshot_count << ".dat";
        save_results(ss.str());
        time_file << snapshot_count << "\t" << 0.0 << "\t" << 0.0 << "\n";
        snapshot_count++;
    }

    // Gradually increase surface potential
    for (int step = 1; step <= n_steps; ++step) {
        double frac = static_cast<double>(step) / n_steps;

        // Exponential charging profile: V(t) = V0 * (1 - exp(-t/τ))
        // Inverted: t/τ = -ln(1 - V/V0) = -ln(1 - frac)
        double t_over_tau = -std::log(1.0 - frac + 1e-10);
        double time_ns = t_over_tau * tau_RC * 1e9;

        // Set current surface potential
        params_.phi_left = phi0_target * frac;

        // Update boundary condition in solution
        phi_[0] = params_.phi_left;

        // Use Boltzmann distribution as initial guess
        update_concentrations_from_phi();

        // Solve steady-state PNP with current BCs
        // Uses Newton-Raphson which is stable
        solve();

        // Save snapshot
        std::ostringstream ss;
        ss << output_dir << "/snapshot_" << std::setfill('0') << std::setw(5) << snapshot_count << ".dat";
        save_results(ss.str());
        time_file << snapshot_count << "\t" << time_ns << "\t" << params_.phi_left * 1e3 << "\n";
        snapshot_count++;

        if (step % std::max(1, n_steps / 10) == 0) {
            std::cout << "  Step " << step << "/" << n_steps
                      << ", φ0 = " << params_.phi_left * 1e3 << " mV"
                      << ", t = " << time_ns << " ns\n";
        }
    }

    // Restore original target potential
    params_.phi_left = phi0_target;

    time_file.close();
    std::cout << "  Total snapshots: " << snapshot_count << "\n";
    std::cout << "  Continuation method completed.\n";
}

void PNPSolver1D::solve_transient_newton(double dt, double t_final,
                                          const std::string& output_dir,
                                          int snapshot_interval) {
    /**
     * Fully implicit Newton method for transient PNP equations.
     *
     * References:
     * - Scharfetter & Gummel (1969) - IEEE Trans. Electron Devices
     * - Bousquet et al. (2018) - SIAM J. Sci. Comput.
     *
     * Governing equations:
     * - Poisson: ∇²φ = -ρ/ε = -(e·NA/ε)(z₊c₊ + z₋c₋)
     * - NP (conservative form): ∂c/∂t = -∇·J where J is outward flux
     *   J₊ = -D₊(∇c₊ + (z₊e/kT)c₊∇φ)  [drift-diffusion]
     *
     * Scharfetter-Gummel flux (from i to i+1):
     *   J_{i+1/2} = (D/dx)[B(v)c_{i+1} - B(-v)c_i]
     *   where v = z·e·(φ_{i+1} - φ_i)/(kT)
     *
     * NP discretization (conservative form):
     *   (c_i^{n+1} - c_i^n)/dt = -(J_{i+1/2} - J_{i-1/2})/dx
     *   = (D/dx²)[B(-v_{i+1/2})c_i - B(v_{i+1/2})c_{i+1}]
     *   - (D/dx²)[B(-v_{i-1/2})c_{i-1} - B(v_{i-1/2})c_i]
     */

    int n_steps = static_cast<int>(t_final / dt);
    int n = params_.N;
    double e_charge = PhysicalConstants::e;
    double kT = PhysicalConstants::kB * params_.T;
    double NA = PhysicalConstants::NA;

    // Newton parameters
    const int max_newton_iter = 50;
    const double newton_tol = 1e-8;
    const double min_damping = 0.01;

    std::cout << "\n========================================\n";
    std::cout << "  Transient Solver (Fully Implicit Newton)\n";
    std::cout << "========================================\n";
    std::cout << "  Time step: " << dt * 1e9 << " ns\n";
    std::cout << "  Total time: " << t_final * 1e9 << " ns\n";
    std::cout << "  Number of steps: " << n_steps << "\n";
    std::cout << "  Snapshot interval: " << snapshot_interval << " steps\n";
    std::cout << "  Newton tolerance: " << newton_tol << "\n";

    // Use UNIFORM grid for stability (non-uniform causes issues with SG scheme)
    std::vector<double> x_uniform(n), dx_uniform(n);
    double L = params_.L;
    double dx = L / (n - 1);
    for (int i = 0; i < n; ++i) {
        x_uniform[i] = i * dx;
        dx_uniform[i] = dx;
    }

    // Initialize with uniform concentration and apply potential instantly
    std::vector<double> phi(n), c_plus(n), c_minus(n);
    std::vector<double> phi_old(n), c_plus_old(n), c_minus_old(n);

    for (int i = 0; i < n; ++i) {
        c_plus[i] = params_.c0;
        c_minus[i] = params_.c0;
        // Linear initial potential (ramps from phi_left to 0)
        phi[i] = params_.phi_left * (1.0 - x_uniform[i] / L);
    }

    // Time info file
    std::ofstream time_file(output_dir + "/time_info.dat");
    time_file << "# snapshot_index  time_ns\n";

    int snapshot_count = 0;

    // Helper: Bernoulli function B(x) = x / (exp(x) - 1)
    // Key properties: B(0) = 1, B(x) + B(-x) = x
    auto bernoulli = [](double x) -> double {
        if (std::abs(x) < 1e-6) {
            // Taylor expansion: B(x) ≈ 1 - x/2 + x²/12 - x⁴/720
            return 1.0 - x/2.0 + x*x/12.0 - x*x*x*x/720.0;
        }
        if (x > 700.0) return 0.0;  // Prevent overflow
        if (x < -700.0) return -x;  // B(x) ≈ -x for large negative x
        return x / (std::exp(x) - 1.0);
    };

    // Helper: Derivative of Bernoulli function dB/dx
    auto bernoulli_deriv = [&bernoulli](double x) -> double {
        if (std::abs(x) < 1e-6) {
            // Taylor expansion: dB/dx ≈ -1/2 + x/6 - x³/180
            return -0.5 + x/6.0 - x*x*x/180.0;
        }
        if (std::abs(x) > 700.0) return 0.0;  // Saturated region
        double ex = std::exp(x);
        double B = x / (ex - 1.0);
        return (1.0 - B - B * ex) / x;
    };

    // Save initial state
    {
        // Copy to member variables for save_results
        x_ = x_uniform;
        phi_ = phi;
        c_plus_ = c_plus;
        c_minus_ = c_minus;

        std::ostringstream ss;
        ss << output_dir << "/snapshot_" << std::setfill('0') << std::setw(5) << snapshot_count << ".dat";
        save_results(ss.str());
        time_file << snapshot_count << "\t" << 0.0 << "\n";
        snapshot_count++;
    }

    // System size: 3 unknowns per grid point (φ, c+, c-)
    // But boundary points are fixed, so we have effective unknowns
    // For simplicity, use full 3n system with boundary conditions embedded

    for (int step = 0; step < n_steps; ++step) {
        double time = (step + 1) * dt;

        // Store old values
        phi_old = phi;
        c_plus_old = c_plus;
        c_minus_old = c_minus;

        // Newton iteration
        int newton_iter = 0;
        double newton_error = 1.0;
        double damping = 1.0;

        while (newton_iter < max_newton_iter && newton_error > newton_tol) {
            // Build residual vector F and Jacobian J
            // Using block structure: for each node i, variables are [φ_i, c+_i, c-_i]

            // Residual vector (3n)
            std::vector<double> F(3*n, 0.0);

            // Jacobian: block tridiagonal with 3x3 blocks
            // Store as: diag_lower[i] (3x3), diag_main[i] (3x3), diag_upper[i] (3x3)
            std::vector<std::array<std::array<double, 3>, 3>> A_lower(n);  // sub-diagonal
            std::vector<std::array<std::array<double, 3>, 3>> A_main(n);   // main diagonal
            std::vector<std::array<std::array<double, 3>, 3>> A_upper(n);  // super-diagonal

            // Initialize to zero
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < 3; ++j) {
                    for (int k = 0; k < 3; ++k) {
                        A_lower[i][j][k] = 0.0;
                        A_main[i][j][k] = 0.0;
                        A_upper[i][j][k] = 0.0;
                    }
                }
            }

            // --- Left boundary (i=0): Dirichlet for φ, zero-flux for c ---
            //
            // Zero-flux BC: J_{1/2} = 0 at left boundary
            // This means no flux enters/exits from the left
            //
            // For stability, we use NP equation at i=0 with ghost cell approach:
            // ∂c₀/∂t = -∇·J = -(J_{1/2} - J_{-1/2})/dx
            // With J_{-1/2} = 0 (no flux from left), we have:
            // ∂c₀/∂t = -J_{1/2}/dx
            //
            // But J_{1/2} = 0 is enforced, so:
            // (c₀^{n+1} - c₀^n)/dt = 0 implies c₀ follows from interior dynamics

            // F[0]: φ_0 = φ_left (Dirichlet)
            F[0] = phi[0] - params_.phi_left;
            A_main[0][0][0] = 1.0;

            // F[1]: NP equation for c+ at i=0 with zero-flux BC (blocking electrode)
            //
            // Conservation: ∂c/∂t = -∇·J  (outward flux convention)
            // At left boundary with zero-flux: J_{-1/2} = 0
            // Discretization: (c₀ - c₀^n)/dt = -(J_{1/2} - J_{-1/2})/dx = -J_{1/2}/dx
            //
            // SG flux (rightward): J_{1/2} = (D/dx)[B(v)c₁ - B(-v)c₀]
            // where v = z·e·(φ₁ - φ₀)/(kT)
            {
                double D_plus = params_.D_plus;
                double z_plus = params_.z_plus;
                double alpha_plus = z_plus * e_charge / kT;

                // Right flux J_{1/2}: from cell 0 to cell 1
                double v_plus = alpha_plus * (phi[1] - phi[0]);
                double Bp = bernoulli(v_plus);
                double Bm = bernoulli(-v_plus);
                double dBp = bernoulli_deriv(v_plus);
                double dBm = bernoulli_deriv(-v_plus);
                double J_half = (D_plus/dx) * (Bp * c_plus[1] - Bm * c_plus[0]);

                // NP equation: (c - c_old)/dt = -J_{1/2}/dx
                // Residual: F = (c - c_old)/dt + J_{1/2}/dx = 0
                F[1] = (c_plus[0] - c_plus_old[0])/dt + J_half/dx;

                // Jacobian: ∂F/∂(unknowns)
                // ∂J_{1/2}/∂c₀ = -(D/dx)*B(-v)
                // ∂J_{1/2}/∂c₁ = (D/dx)*B(v)
                // ∂J_{1/2}/∂φ₀ = (D/dx)*α*[-dB(v)*c₁ - dB(-v)*c₀]
                // ∂J_{1/2}/∂φ₁ = (D/dx)*α*[dB(v)*c₁ + dB(-v)*c₀]

                // ∂F/∂c₀ = 1/dt + (1/dx)*(∂J/∂c₀) = 1/dt - (D/dx²)*B(-v)
                A_main[0][1][1] = 1.0/dt - (D_plus/(dx*dx)) * Bm;
                // ∂F/∂c₁ = (1/dx)*(∂J/∂c₁) = (D/dx²)*B(v)
                A_upper[0][1][1] = (D_plus/(dx*dx)) * Bp;
                // ∂F/∂φ₀ = (1/dx)*(∂J/∂φ₀)
                A_main[0][1][0] = (D_plus/(dx*dx)) * alpha_plus * (-dBp * c_plus[1] - dBm * c_plus[0]);
                // ∂F/∂φ₁ = (1/dx)*(∂J/∂φ₁)
                A_upper[0][1][0] = (D_plus/(dx*dx)) * alpha_plus * (dBp * c_plus[1] + dBm * c_plus[0]);
            }

            // F[2]: NP equation for c- at i=0 with zero-flux BC
            {
                double D_minus = params_.D_minus;
                double z_minus = params_.z_minus;
                double alpha_minus = z_minus * e_charge / kT;

                double v_minus = alpha_minus * (phi[1] - phi[0]);
                double Bp = bernoulli(v_minus);
                double Bm = bernoulli(-v_minus);
                double dBp = bernoulli_deriv(v_minus);
                double dBm = bernoulli_deriv(-v_minus);
                double J_half = (D_minus/dx) * (Bp * c_minus[1] - Bm * c_minus[0]);

                F[2] = (c_minus[0] - c_minus_old[0])/dt + J_half/dx;

                A_main[0][2][2] = 1.0/dt - (D_minus/(dx*dx)) * Bm;
                A_upper[0][2][2] = (D_minus/(dx*dx)) * Bp;
                A_main[0][2][0] = (D_minus/(dx*dx)) * alpha_minus * (-dBp * c_minus[1] - dBm * c_minus[0]);
                A_upper[0][2][0] = (D_minus/(dx*dx)) * alpha_minus * (dBp * c_minus[1] + dBm * c_minus[0]);
            }

            // --- Interior points (i=1 to n-2) ---
            for (int i = 1; i < n - 1; ++i) {
                int idx = 3 * i;
                double dx_m = x_uniform[i] - x_uniform[i-1];
                double dx_p = x_uniform[i+1] - x_uniform[i];
                double dx_avg = 0.5 * (dx_m + dx_p);

                // --- Poisson equation ---
                // F_φ = (φ_{i+1} - 2φ_i + φ_{i-1})/dx² + (e*NA/ε)(z+*c+ + z-*c-)
                double laplacian_coeff = 1.0 / (dx * dx);
                double rho_coeff = e_charge * NA / eps_;

                F[idx] = laplacian_coeff * (phi[i+1] - 2.0*phi[i] + phi[i-1])
                       + rho_coeff * (params_.z_plus * c_plus[i] + params_.z_minus * c_minus[i]);

                // ∂F_φ/∂φ_{i-1}
                A_lower[i][0][0] = laplacian_coeff;
                // ∂F_φ/∂φ_i
                A_main[i][0][0] = -2.0 * laplacian_coeff;
                // ∂F_φ/∂φ_{i+1}
                A_upper[i][0][0] = laplacian_coeff;
                // ∂F_φ/∂c+_i
                A_main[i][0][1] = rho_coeff * params_.z_plus;
                // ∂F_φ/∂c-_i
                A_main[i][0][2] = rho_coeff * params_.z_minus;

                // --- NP equation for c+ ---
                //
                // Conservation law: ∂c/∂t + ∇·J = 0
                // Discrete: (c - c_old)/dt + (J_{i+1/2} - J_{i-1/2})/dx = 0
                //
                // SG flux (positive = rightward flow):
                // J_{i+1/2} = (D/dx)[B(v)c_{i+1} - B(-v)c_i]
                //   where v = z·e·(φ_{i+1} - φ_i)/(kT)

                double D_plus = params_.D_plus;
                double z_plus = params_.z_plus;
                double alpha_plus = z_plus * e_charge / kT;

                // Right flux (i+1/2): from i to i+1
                double v_plus_r = alpha_plus * (phi[i+1] - phi[i]);
                double Bp_r = bernoulli(v_plus_r);
                double Bm_r = bernoulli(-v_plus_r);
                double dBp_r = bernoulli_deriv(v_plus_r);
                double dBm_r = bernoulli_deriv(-v_plus_r);
                double J_plus_r = (D_plus/dx) * (Bp_r * c_plus[i+1] - Bm_r * c_plus[i]);

                // Left flux (i-1/2): from i-1 to i
                double v_plus_l = alpha_plus * (phi[i] - phi[i-1]);
                double Bp_l = bernoulli(v_plus_l);
                double Bm_l = bernoulli(-v_plus_l);
                double dBp_l = bernoulli_deriv(v_plus_l);
                double dBm_l = bernoulli_deriv(-v_plus_l);
                double J_plus_l = (D_plus/dx) * (Bp_l * c_plus[i] - Bm_l * c_plus[i-1]);

                // Residual: F = (c - c_old)/dt + (J_r - J_l)/dx
                F[idx+1] = (c_plus[i] - c_plus_old[i])/dt + (J_plus_r - J_plus_l)/dx;

                // Jacobian for NP+
                // F = (c_i - c_i^n)/dt + (J_r - J_l)/dx
                //
                // J_r = (D/dx)[B(v_r)c_{i+1} - B(-v_r)c_i]
                // J_l = (D/dx)[B(v_l)c_i - B(-v_l)c_{i-1}]
                //
                // ∂F/∂c_i = 1/dt + (1/dx)[∂J_r/∂c_i - ∂J_l/∂c_i]
                //         = 1/dt + (D/dx²)[-B(-v_r) - B(v_l)]
                //         = 1/dt - (D/dx²)[B(-v_r) + B(v_l)]
                A_main[i][1][1] = 1.0/dt - (D_plus/(dx*dx)) * (Bm_r + Bp_l);

                // ∂F/∂c_{i-1} = (1/dx)[-∂J_l/∂c_{i-1}] = (D/dx²)*B(-v_l)
                A_lower[i][1][1] = (D_plus/(dx*dx)) * Bm_l;

                // ∂F/∂c_{i+1} = (1/dx)[∂J_r/∂c_{i+1}] = (D/dx²)*B(v_r)
                A_upper[i][1][1] = (D_plus/(dx*dx)) * Bp_r;

                // ∂J_r/∂φ_i = (D/dx)*α*(-dB(v_r)c_{i+1} - dB(-v_r)c_i)
                // ∂J_r/∂φ_{i+1} = (D/dx)*α*(dB(v_r)c_{i+1} + dB(-v_r)c_i)
                // ∂J_l/∂φ_{i-1} = (D/dx)*α*(-dB(v_l)c_i - dB(-v_l)c_{i-1})
                // ∂J_l/∂φ_i = (D/dx)*α*(dB(v_l)c_i + dB(-v_l)c_{i-1})
                //
                // ∂F/∂φ_i = (1/dx)[∂J_r/∂φ_i - ∂J_l/∂φ_i]
                A_main[i][1][0] = (D_plus/(dx*dx)) * alpha_plus *
                    ((-dBp_r * c_plus[i+1] - dBm_r * c_plus[i]) - (dBp_l * c_plus[i] + dBm_l * c_plus[i-1]));
                // ∂F/∂φ_{i-1} = -(1/dx)*∂J_l/∂φ_{i-1}
                A_lower[i][1][0] = -(D_plus/(dx*dx)) * alpha_plus *
                    (-dBp_l * c_plus[i] - dBm_l * c_plus[i-1]);
                // ∂F/∂φ_{i+1} = (1/dx)*∂J_r/∂φ_{i+1}
                A_upper[i][1][0] = (D_plus/(dx*dx)) * alpha_plus *
                    (dBp_r * c_plus[i+1] + dBm_r * c_plus[i]);

                // --- NP equation for c- ---
                double D_minus = params_.D_minus;
                double z_minus = params_.z_minus;
                double alpha_minus = z_minus * e_charge / kT;

                // Right flux (i+1/2)
                double v_minus_r = alpha_minus * (phi[i+1] - phi[i]);
                double Bp_mr = bernoulli(v_minus_r);
                double Bm_mr = bernoulli(-v_minus_r);
                double dBp_mr = bernoulli_deriv(v_minus_r);
                double dBm_mr = bernoulli_deriv(-v_minus_r);
                double J_minus_r = (D_minus/dx) * (Bp_mr * c_minus[i+1] - Bm_mr * c_minus[i]);

                // Left flux (i-1/2)
                double v_minus_l = alpha_minus * (phi[i] - phi[i-1]);
                double Bp_ml = bernoulli(v_minus_l);
                double Bm_ml = bernoulli(-v_minus_l);
                double dBp_ml = bernoulli_deriv(v_minus_l);
                double dBm_ml = bernoulli_deriv(-v_minus_l);
                double J_minus_l = (D_minus/dx) * (Bp_ml * c_minus[i] - Bm_ml * c_minus[i-1]);

                F[idx+2] = (c_minus[i] - c_minus_old[i])/dt + (J_minus_r - J_minus_l)/dx;

                // Jacobian for NP- (same structure as NP+)
                A_main[i][2][2] = 1.0/dt - (D_minus/(dx*dx)) * (Bm_mr + Bp_ml);
                A_lower[i][2][2] = (D_minus/(dx*dx)) * Bm_ml;
                A_upper[i][2][2] = (D_minus/(dx*dx)) * Bp_mr;

                A_main[i][2][0] = (D_minus/(dx*dx)) * alpha_minus *
                    ((-dBp_mr * c_minus[i+1] - dBm_mr * c_minus[i]) - (dBp_ml * c_minus[i] + dBm_ml * c_minus[i-1]));
                A_lower[i][2][0] = -(D_minus/(dx*dx)) * alpha_minus *
                    (-dBp_ml * c_minus[i] - dBm_ml * c_minus[i-1]);
                A_upper[i][2][0] = (D_minus/(dx*dx)) * alpha_minus *
                    (dBp_mr * c_minus[i+1] + dBm_mr * c_minus[i]);
            }

            // --- Right boundary (i=n-1): Dirichlet for all ---
            int idx_last = 3 * (n - 1);

            // φ_{n-1} = φ_right = 0
            F[idx_last] = phi[n-1] - params_.phi_right;
            A_main[n-1][0][0] = 1.0;

            // c+_{n-1} = c0
            F[idx_last+1] = c_plus[n-1] - params_.c0;
            A_main[n-1][1][1] = 1.0;

            // c-_{n-1} = c0
            F[idx_last+2] = c_minus[n-1] - params_.c0;
            A_main[n-1][2][2] = 1.0;

            // --- Solve block tridiagonal system J * δU = -F ---
            // Using block Thomas algorithm

            // Forward elimination
            std::vector<std::array<std::array<double, 3>, 3>> A_mod(n);
            std::vector<std::array<double, 3>> F_mod(n);

            // Copy initial values
            for (int j = 0; j < 3; ++j) {
                F_mod[0][j] = -F[j];
                for (int k = 0; k < 3; ++k) {
                    A_mod[0][j][k] = A_main[0][j][k];
                }
            }

            for (int i = 1; i < n; ++i) {
                // Compute L = A_lower[i] * A_mod[i-1]^{-1}
                // Then A_mod[i] = A_main[i] - L * A_upper[i-1]
                // F_mod[i] = -F[3i:3i+3] - L * F_mod[i-1]

                // Invert A_mod[i-1] (3x3 matrix)
                std::array<std::array<double, 3>, 3> A_inv;
                double det = A_mod[i-1][0][0] * (A_mod[i-1][1][1]*A_mod[i-1][2][2] - A_mod[i-1][1][2]*A_mod[i-1][2][1])
                           - A_mod[i-1][0][1] * (A_mod[i-1][1][0]*A_mod[i-1][2][2] - A_mod[i-1][1][2]*A_mod[i-1][2][0])
                           + A_mod[i-1][0][2] * (A_mod[i-1][1][0]*A_mod[i-1][2][1] - A_mod[i-1][1][1]*A_mod[i-1][2][0]);

                if (std::abs(det) < 1e-30) {
                    std::cerr << "Warning: Singular matrix at node " << i << ", det = " << det << "\n";
                    det = (det > 0) ? 1e-30 : -1e-30;
                }

                A_inv[0][0] = (A_mod[i-1][1][1]*A_mod[i-1][2][2] - A_mod[i-1][1][2]*A_mod[i-1][2][1]) / det;
                A_inv[0][1] = (A_mod[i-1][0][2]*A_mod[i-1][2][1] - A_mod[i-1][0][1]*A_mod[i-1][2][2]) / det;
                A_inv[0][2] = (A_mod[i-1][0][1]*A_mod[i-1][1][2] - A_mod[i-1][0][2]*A_mod[i-1][1][1]) / det;
                A_inv[1][0] = (A_mod[i-1][1][2]*A_mod[i-1][2][0] - A_mod[i-1][1][0]*A_mod[i-1][2][2]) / det;
                A_inv[1][1] = (A_mod[i-1][0][0]*A_mod[i-1][2][2] - A_mod[i-1][0][2]*A_mod[i-1][2][0]) / det;
                A_inv[1][2] = (A_mod[i-1][0][2]*A_mod[i-1][1][0] - A_mod[i-1][0][0]*A_mod[i-1][1][2]) / det;
                A_inv[2][0] = (A_mod[i-1][1][0]*A_mod[i-1][2][1] - A_mod[i-1][1][1]*A_mod[i-1][2][0]) / det;
                A_inv[2][1] = (A_mod[i-1][0][1]*A_mod[i-1][2][0] - A_mod[i-1][0][0]*A_mod[i-1][2][1]) / det;
                A_inv[2][2] = (A_mod[i-1][0][0]*A_mod[i-1][1][1] - A_mod[i-1][0][1]*A_mod[i-1][1][0]) / det;

                // L = A_lower[i] * A_inv
                std::array<std::array<double, 3>, 3> L;
                for (int j = 0; j < 3; ++j) {
                    for (int k = 0; k < 3; ++k) {
                        L[j][k] = 0.0;
                        for (int m = 0; m < 3; ++m) {
                            L[j][k] += A_lower[i][j][m] * A_inv[m][k];
                        }
                    }
                }

                // A_mod[i] = A_main[i] - L * A_upper[i-1]
                for (int j = 0; j < 3; ++j) {
                    for (int k = 0; k < 3; ++k) {
                        A_mod[i][j][k] = A_main[i][j][k];
                        for (int m = 0; m < 3; ++m) {
                            A_mod[i][j][k] -= L[j][m] * A_upper[i-1][m][k];
                        }
                    }
                }

                // F_mod[i] = -F[3i:] - L * F_mod[i-1]
                for (int j = 0; j < 3; ++j) {
                    F_mod[i][j] = -F[3*i + j];
                    for (int m = 0; m < 3; ++m) {
                        F_mod[i][j] -= L[j][m] * F_mod[i-1][m];
                    }
                }
            }

            // Back substitution
            std::vector<std::array<double, 3>> delta(n);

            // Solve A_mod[n-1] * delta[n-1] = F_mod[n-1]
            {
                int i = n - 1;
                std::array<std::array<double, 3>, 3> A_inv;
                double det = A_mod[i][0][0] * (A_mod[i][1][1]*A_mod[i][2][2] - A_mod[i][1][2]*A_mod[i][2][1])
                           - A_mod[i][0][1] * (A_mod[i][1][0]*A_mod[i][2][2] - A_mod[i][1][2]*A_mod[i][2][0])
                           + A_mod[i][0][2] * (A_mod[i][1][0]*A_mod[i][2][1] - A_mod[i][1][1]*A_mod[i][2][0]);

                if (std::abs(det) < 1e-30) {
                    det = (det > 0) ? 1e-30 : -1e-30;
                }

                A_inv[0][0] = (A_mod[i][1][1]*A_mod[i][2][2] - A_mod[i][1][2]*A_mod[i][2][1]) / det;
                A_inv[0][1] = (A_mod[i][0][2]*A_mod[i][2][1] - A_mod[i][0][1]*A_mod[i][2][2]) / det;
                A_inv[0][2] = (A_mod[i][0][1]*A_mod[i][1][2] - A_mod[i][0][2]*A_mod[i][1][1]) / det;
                A_inv[1][0] = (A_mod[i][1][2]*A_mod[i][2][0] - A_mod[i][1][0]*A_mod[i][2][2]) / det;
                A_inv[1][1] = (A_mod[i][0][0]*A_mod[i][2][2] - A_mod[i][0][2]*A_mod[i][2][0]) / det;
                A_inv[1][2] = (A_mod[i][0][2]*A_mod[i][1][0] - A_mod[i][0][0]*A_mod[i][1][2]) / det;
                A_inv[2][0] = (A_mod[i][1][0]*A_mod[i][2][1] - A_mod[i][1][1]*A_mod[i][2][0]) / det;
                A_inv[2][1] = (A_mod[i][0][1]*A_mod[i][2][0] - A_mod[i][0][0]*A_mod[i][2][1]) / det;
                A_inv[2][2] = (A_mod[i][0][0]*A_mod[i][1][1] - A_mod[i][0][1]*A_mod[i][1][0]) / det;

                for (int j = 0; j < 3; ++j) {
                    delta[i][j] = 0.0;
                    for (int k = 0; k < 3; ++k) {
                        delta[i][j] += A_inv[j][k] * F_mod[i][k];
                    }
                }
            }

            for (int i = n - 2; i >= 0; --i) {
                // delta[i] = A_mod[i]^{-1} * (F_mod[i] - A_upper[i] * delta[i+1])
                std::array<double, 3> rhs;
                for (int j = 0; j < 3; ++j) {
                    rhs[j] = F_mod[i][j];
                    for (int k = 0; k < 3; ++k) {
                        rhs[j] -= A_upper[i][j][k] * delta[i+1][k];
                    }
                }

                std::array<std::array<double, 3>, 3> A_inv;
                double det = A_mod[i][0][0] * (A_mod[i][1][1]*A_mod[i][2][2] - A_mod[i][1][2]*A_mod[i][2][1])
                           - A_mod[i][0][1] * (A_mod[i][1][0]*A_mod[i][2][2] - A_mod[i][1][2]*A_mod[i][2][0])
                           + A_mod[i][0][2] * (A_mod[i][1][0]*A_mod[i][2][1] - A_mod[i][1][1]*A_mod[i][2][0]);

                if (std::abs(det) < 1e-30) {
                    det = (det > 0) ? 1e-30 : -1e-30;
                }

                A_inv[0][0] = (A_mod[i][1][1]*A_mod[i][2][2] - A_mod[i][1][2]*A_mod[i][2][1]) / det;
                A_inv[0][1] = (A_mod[i][0][2]*A_mod[i][2][1] - A_mod[i][0][1]*A_mod[i][2][2]) / det;
                A_inv[0][2] = (A_mod[i][0][1]*A_mod[i][1][2] - A_mod[i][0][2]*A_mod[i][1][1]) / det;
                A_inv[1][0] = (A_mod[i][1][2]*A_mod[i][2][0] - A_mod[i][1][0]*A_mod[i][2][2]) / det;
                A_inv[1][1] = (A_mod[i][0][0]*A_mod[i][2][2] - A_mod[i][0][2]*A_mod[i][2][0]) / det;
                A_inv[1][2] = (A_mod[i][0][2]*A_mod[i][1][0] - A_mod[i][0][0]*A_mod[i][1][2]) / det;
                A_inv[2][0] = (A_mod[i][1][0]*A_mod[i][2][1] - A_mod[i][1][1]*A_mod[i][2][0]) / det;
                A_inv[2][1] = (A_mod[i][0][1]*A_mod[i][2][0] - A_mod[i][0][0]*A_mod[i][2][1]) / det;
                A_inv[2][2] = (A_mod[i][0][0]*A_mod[i][1][1] - A_mod[i][0][1]*A_mod[i][1][0]) / det;

                for (int j = 0; j < 3; ++j) {
                    delta[i][j] = 0.0;
                    for (int k = 0; k < 3; ++k) {
                        delta[i][j] += A_inv[j][k] * rhs[k];
                    }
                }
            }

            // Compute error and apply damping
            newton_error = 0.0;
            double max_delta_phi = 0.0, max_delta_c = 0.0;
            for (int i = 0; i < n; ++i) {
                max_delta_phi = std::max(max_delta_phi, std::abs(delta[i][0]));
                max_delta_c = std::max(max_delta_c, std::abs(delta[i][1]/params_.c0));
                max_delta_c = std::max(max_delta_c, std::abs(delta[i][2]/params_.c0));
            }
            newton_error = std::max(max_delta_phi/phi_T_, max_delta_c);

            // Adaptive damping to ensure positivity
            damping = 1.0;
            for (int i = 0; i < n; ++i) {
                if (c_plus[i] + damping * delta[i][1] < 0.1 * params_.c0 * 1e-6) {
                    double new_damp = 0.9 * c_plus[i] / std::abs(delta[i][1] + 1e-30);
                    damping = std::min(damping, std::max(new_damp, min_damping));
                }
                if (c_minus[i] + damping * delta[i][2] < 0.1 * params_.c0 * 1e-6) {
                    double new_damp = 0.9 * c_minus[i] / std::abs(delta[i][2] + 1e-30);
                    damping = std::min(damping, std::max(new_damp, min_damping));
                }
            }

            // Update solution
            for (int i = 0; i < n; ++i) {
                phi[i] += damping * delta[i][0];
                c_plus[i] += damping * delta[i][1];
                c_minus[i] += damping * delta[i][2];

                // Ensure strict positivity
                c_plus[i] = std::max(c_plus[i], params_.c0 * 1e-12);
                c_minus[i] = std::max(c_minus[i], params_.c0 * 1e-12);
            }

            newton_iter++;
        }

        // Save snapshot
        if ((step + 1) % snapshot_interval == 0) {
            x_ = x_uniform;
            phi_ = phi;
            c_plus_ = c_plus;
            c_minus_ = c_minus;

            std::ostringstream ss;
            ss << output_dir << "/snapshot_" << std::setfill('0') << std::setw(5) << snapshot_count << ".dat";
            save_results(ss.str());
            double time_ns = (step + 1) * dt * 1e9;
            time_file << snapshot_count << "\t" << time_ns << "\n";
            snapshot_count++;
        }

        // Progress
        if (step % std::max(1, n_steps / 20) == 0 || step == n_steps - 1) {
            double time_ns = (step + 1) * dt * 1e9;
            double c_plus_surface = c_plus[0] / params_.c0;
            double c_minus_surface = c_minus[0] / params_.c0;
            std::cout << "  t = " << std::fixed << std::setprecision(2) << time_ns
                      << " ns, Newton iter = " << newton_iter
                      << ", c+/c0 = " << std::scientific << std::setprecision(2) << c_plus_surface
                      << ", c-/c0 = " << c_minus_surface
                      << ", damping = " << std::fixed << std::setprecision(2) << damping << "\n";
        }
    }

    // Final update to member variables
    x_ = x_uniform;
    phi_ = phi;
    c_plus_ = c_plus;
    c_minus_ = c_minus;

    time_file.close();
    std::cout << "  Total snapshots: " << snapshot_count << "\n";
    std::cout << "  Fully implicit Newton method completed.\n";
}

void PNPSolver1D::solve_transient_slotboom(double dt, double t_final,
                                            const std::string& output_dir,
                                            int snapshot_interval) {
    /**
     * Transient PNP solver using Scharfetter-Gummel discretization with
     * Slotboom transformation for guaranteed positivity and stability.
     *
     * The key insight is that the Nernst-Planck flux can be written as:
     *   J± = -D± * exp(∓φ/φT) * ∇u±
     * where u± = c± * exp(±φ/φT) are Slotboom variables.
     *
     * Using harmonic mean approximation for the exponential factor at cell faces:
     *   exp(-φ/φT)|_{i+1/2} ≈ 2 / (exp(φᵢ/φT) + exp(φᵢ₊₁/φT))
     *
     * This leads to the Scharfetter-Gummel flux formula:
     *   J_{i+1/2} = (D/dx) * [B(v)*c_{i+1} - B(-v)*c_i]
     * where v = (φᵢ - φᵢ₊₁)/φT and B(x) = x/(exp(x)-1) is the Bernoulli function.
     *
     * The scheme is:
     * - Unconditionally stable (implicit Euler for time)
     * - Positivity preserving (M-matrix property)
     * - Mass conserving
     *
     * References:
     * - Scharfetter & Gummel (1969), IEEE Trans. Electron Devices
     * - Liu & Wang (2021), Numerische Mathematik
     */

    int n_steps = static_cast<int>(t_final / dt + 0.5);
    int n = params_.N;

    std::cout << "\n========================================\n";
    std::cout << "  Transient Solver (Scharfetter-Gummel)\n";
    std::cout << "========================================\n";
    std::cout << "  Time step: " << dt * 1e9 << " ns\n";
    std::cout << "  Total time: " << t_final * 1e9 << " ns\n";
    std::cout << "  Number of steps: " << n_steps << "\n";
    std::cout << "  Snapshot interval: " << snapshot_interval << " steps\n";

    // Use uniform grid
    double L = params_.L;
    double dx = L / (n - 1);
    std::vector<double> x_uniform(n);
    for (int i = 0; i < n; ++i) {
        x_uniform[i] = i * dx;
    }

    // Initialize
    std::vector<double> phi(n), c_plus(n), c_minus(n);

    for (int i = 0; i < n; ++i) {
        c_plus[i] = params_.c0;
        c_minus[i] = params_.c0;
        phi[i] = params_.phi_left * (1.0 - x_uniform[i] / L);
    }

    // Time info file
    std::ofstream time_file(output_dir + "/time_info.dat");
    time_file << "# snapshot_index  time_ns\n";

    int snapshot_count = 0;

    // Bernoulli function: B(x) = x / (exp(x) - 1)
    // For numerical stability, use Taylor series near x=0
    auto bernoulli_fn = [](double x) -> double {
        if (std::abs(x) < 1e-6) {
            // Taylor: B(x) ≈ 1 - x/2 + x²/12 - x⁴/720
            return 1.0 - x/2.0 + x*x/12.0;
        }
        return x / (std::exp(x) - 1.0);
    };

    // Helper: Solve tridiagonal system Ax = d
    auto solve_tridiag = [](const std::vector<double>& a,  // sub-diagonal
                            const std::vector<double>& b,  // main diagonal
                            const std::vector<double>& c,  // super-diagonal
                            std::vector<double>& d,        // rhs, overwritten with solution
                            int n) {
        std::vector<double> c_mod(n), d_mod(n);
        c_mod[0] = c[0] / b[0];
        d_mod[0] = d[0] / b[0];
        for (int i = 1; i < n; ++i) {
            double denom = b[i] - a[i] * c_mod[i-1];
            if (std::abs(denom) < 1e-30) {
                // Near-singular matrix, use regularization
                denom = (denom >= 0) ? 1e-30 : -1e-30;
            }
            if (i < n - 1) c_mod[i] = c[i] / denom;
            d_mod[i] = (d[i] - a[i] * d_mod[i-1]) / denom;
        }
        d[n-1] = d_mod[n-1];
        for (int i = n - 2; i >= 0; --i) {
            d[i] = d_mod[i] - c_mod[i] * d[i+1];
        }
    };

    // Solve Poisson equation: ∇²φ = -ρ/ε
    double e_charge = PhysicalConstants::e;
    double rho_coeff = e_charge * PhysicalConstants::NA / eps_;
    double dx2 = dx * dx;

    auto solve_poisson = [&]() {
        const int max_iter = 100;
        const double tol = 1e-10;

        for (int iter = 0; iter < max_iter; ++iter) {
            double max_residual = 0.0;

            std::vector<double> a(n, 0.0), b(n, 0.0), c(n, 0.0), rhs(n, 0.0);

            // Left boundary: φ = φ_left
            b[0] = 1.0;
            rhs[0] = params_.phi_left - phi[0];

            // Interior points
            for (int i = 1; i < n - 1; ++i) {
                double laplacian = (phi[i+1] - 2.0*phi[i] + phi[i-1]) / dx2;
                double rho = rho_coeff * (params_.z_plus * c_plus[i] + params_.z_minus * c_minus[i]);
                double residual = laplacian + rho;
                max_residual = std::max(max_residual, std::abs(residual));

                a[i] = 1.0 / dx2;
                b[i] = -2.0 / dx2;
                c[i] = 1.0 / dx2;
                rhs[i] = -residual;
            }

            // Right boundary: φ = 0
            b[n-1] = 1.0;
            rhs[n-1] = params_.phi_right - phi[n-1];

            if (max_residual < tol) break;

            solve_tridiag(a, b, c, rhs, n);
            for (int i = 0; i < n; ++i) {
                phi[i] += rhs[i];
            }
        }
    };

    // Solve NP equation using Scharfetter-Gummel discretization (implicit)
    //
    // The Nernst-Planck equation is:
    //   ∂c/∂t = -∇·J  where  J = -D(∇c + (ze/kT)c∇φ)
    //         = ∇·[D(∇c + (ze/kT)c∇φ)]
    //
    // For the flux from cell i to cell i+1 (in +x direction):
    //   J_{i+1/2} = (D/dx)[B(-v)c_{i+1} - B(v)c_i]
    // where v = z*(φ_{i+1} - φ_i)/φT
    //
    // This formula gives positive flux (from i to i+1) when c_i > c_{i+1} (diffusion)
    // and also accounts for drift due to electric field.
    //
    // Conservation: (c_i^{n+1} - c_i^n)/dt = -(J_{i+1/2} - J_{i-1/2})/dx
    //             = (1/dx)[J_{i-1/2} - J_{i+1/2}]

    auto solve_np_implicit = [&](std::vector<double>& c, int z, double D) {
        std::vector<double> a_coef(n, 0.0), b_coef(n, 0.0), c_coef(n, 0.0), rhs(n, 0.0);
        double r = D * dt / (dx * dx);

        // Compute SG coefficients at each face
        // v_{i+1/2} = z * (φ_{i+1} - φ_i) / φT  (potential difference in +x direction)
        std::vector<double> B_v(n-1), B_mv(n-1);  // B(v) and B(-v) at faces

        for (int i = 0; i < n - 1; ++i) {
            double v = z * (phi[i+1] - phi[i]) / phi_T_;
            B_v[i] = bernoulli_fn(v);
            B_mv[i] = bernoulli_fn(-v);
        }

        // The SG flux from i to i+1:
        // J_{i+1/2} = (D/dx)[B(-v_{i+1/2})*c_{i+1} - B(v_{i+1/2})*c_i]
        //
        // Mass conservation at cell i:
        // (c_i^{n+1} - c_i^n)/dt = (J_{i-1/2} - J_{i+1/2})/dx
        //
        // J_{i-1/2} = (D/dx)[B(-v_{i-1/2})*c_i - B(v_{i-1/2})*c_{i-1}]
        // J_{i+1/2} = (D/dx)[B(-v_{i+1/2})*c_{i+1} - B(v_{i+1/2})*c_i]
        //
        // (c_i^{n+1} - c_i^n)/dt = (D/dx²)[B(-v_{i-1/2})*c_i - B(v_{i-1/2})*c_{i-1}
        //                                 - B(-v_{i+1/2})*c_{i+1} + B(v_{i+1/2})*c_i]
        //
        // Collecting terms:
        // c_i^{n+1}[1] + c_{i-1}^{n+1}[-r*B(v_{i-1/2})]
        //             + c_i^{n+1}[r*(B(-v_{i-1/2}) + B(v_{i+1/2}))]
        //             + c_{i+1}^{n+1}[-r*B(-v_{i+1/2})] = c_i^n
        //
        // But wait, for implicit Euler we need:
        // c_i^{n+1} - dt*(RHS with c^{n+1}) = c_i^n
        //
        // Let's define: coefficient of c_{i-1} is alpha, c_i is beta, c_{i+1} is gamma
        // alpha = -r * B(v_{i-1/2})
        // beta = 1 + r * (B(-v_{i-1/2}) + B(v_{i+1/2}))
        // gamma = -r * B(-v_{i+1/2})

        // Left boundary (i=0): zero-flux => J_{-1/2} = 0
        // Conservation: (c_0^{n+1} - c_0^n)/dt = (J_{-1/2} - J_{1/2})/dx = -J_{1/2}/dx
        // c_0^{n+1} + (dt/dx)*J_{1/2}^{n+1} = c_0^n
        // c_0^{n+1} + r*[B(-v_{1/2})*c_1 - B(v_{1/2})*c_0] = c_0^n
        // c_0*(1 + r*B(v_{1/2})) - r*B(-v_{1/2})*c_1 = c_0^n
        // Wait, this doesn't match the M-matrix property. Let me reconsider.
        //
        // Actually, with zero-flux at left: J_{-1/2} = 0
        // The cell 0 balance is:
        // (c_0^{n+1} - c_0^n)/dt = (J_{-1/2} - J_{1/2})/dx = -J_{1/2}/dx
        //
        // J_{1/2} = (D/dx)[B(-v_{1/2})*c_1 - B(v_{1/2})*c_0]  (flux leaving cell 0)
        //
        // c_0^{n+1} - c_0^n = -(dt/dx)*J_{1/2}
        //                  = -(D*dt/dx²)[B(-v_{1/2})*c_1 - B(v_{1/2})*c_0]
        //                  = r*[B(v_{1/2})*c_0 - B(-v_{1/2})*c_1]
        //
        // c_0*(1 - r*B(v_{1/2})) + r*B(-v_{1/2})*c_1 = c_0^n  -- WRONG SIGN
        //
        // Hmm, the issue is that M-matrix requires positive diagonal and negative off-diagonals.
        // Let me reconsider the flux direction convention.
        //
        // Standard SG convention: J_{i+1/2} is the flux from LEFT to RIGHT
        // J_{i+1/2} = (D/dx) * [B(v)*c_i - B(-v)*c_{i+1}]  where v = z*(φ_{i+1}-φ_i)/φT
        //
        // Wait, I had it backwards. Let me use the standard convention:
        // J_{i→i+1} = (D/dx) * [B(v)*c_i - B(-v)*c_{i+1}]
        //
        // For cation (z=+1) with positive potential on left (φ_i > φ_{i+1}):
        // v = (φ_{i+1} - φ_i)/φT < 0
        // B(v) → large (exponential growth), B(-v) → small
        // So J > 0 means flow from i to i+1, i.e., cation repelled from high potential
        // Wait, that's wrong for cations which should be attracted to negative electrode.
        //
        // Let me think again. For cations at positive electrode (high φ):
        // - Electric field points away from electrode (E = -∇φ > 0 if φ decreases)
        // - Force on cation = z*e*E = +e*E, pointing away from electrode
        // - So cations are repelled from positive electrode
        // - This means c+ should decrease near positive electrode
        //
        // But we have positive potential at left! That means cations should be repelled,
        // not attracted. Let me re-check the expected physics...
        //
        // At t=0, we apply +100mV at left electrode.
        // For equilibrium (Boltzmann):
        //   c+ = c0 * exp(-z*e*φ/kT) = c0 * exp(-φ/φT) → decreased at positive electrode
        //   c- = c0 * exp(+φ/φT) → increased at positive electrode
        //
        // OK so c+ should decrease and c- should increase near the left electrode.
        // Let's check if SG discretization gives this behavior.
        //
        // For c+ (z=+1): drift is towards lower potential (opposite to E field)
        // v = +1 * (φ_{i+1} - φ_i)/φT
        // At left boundary, φ_0 > φ_1, so v < 0
        // B(v) is large, B(-v) is small
        // J_{1/2} = (D/dx)[B(v)*c_0 - B(-v)*c_1]
        //         ≈ (D/dx)*B(v)*c_0  (B(v) >> B(-v))
        // This is positive flux from 0 to 1, meaning cations flow away from left electrode
        // Good! That matches expected physics.
        //
        // Now for the matrix assembly. Conservation:
        // (c_i - c_i^n)/dt = -(J_{i+1/2} - J_{i-1/2})/dx  (flux OUT minus flux IN)
        //
        // J_{i+1/2} = (D/dx)[B(v_{i+1/2})*c_i - B(-v_{i+1/2})*c_{i+1}]  (leaving cell i)
        // J_{i-1/2} = (D/dx)[B(v_{i-1/2})*c_{i-1} - B(-v_{i-1/2})*c_i]  (leaving cell i-1 = entering cell i)
        //
        // Net flux out of cell i: J_{i+1/2} - J_{i-1/2}
        // = (D/dx)[B(v_{i+1/2})*c_i - B(-v_{i+1/2})*c_{i+1} - B(v_{i-1/2})*c_{i-1} + B(-v_{i-1/2})*c_i]
        //
        // (c_i - c_i^n)/dt = -(D/dx²)[B(v_{i+1/2})*c_i - B(-v_{i+1/2})*c_{i+1}
        //                            - B(v_{i-1/2})*c_{i-1} + B(-v_{i-1/2})*c_i]
        //
        // c_i + r*[B(v_{i+1/2}) + B(-v_{i-1/2})]*c_i - r*B(-v_{i+1/2})*c_{i+1} - r*B(v_{i-1/2})*c_{i-1} = c_i^n
        //
        // Rearranging: a*c_{i-1} + b*c_i + c*c_{i+1} = d
        // where:
        //   a = -r * B(v_{i-1/2})
        //   b = 1 + r * (B(v_{i+1/2}) + B(-v_{i-1/2}))
        //   c = -r * B(-v_{i+1/2})
        //   d = c_i^n
        //
        // For M-matrix: need b > 0 and a, c <= 0
        // B(x) > 0 for all x, so:
        //   a = -r * B(v_{i-1/2}) < 0 ✓
        //   c = -r * B(-v_{i+1/2}) < 0 ✓
        //   b = 1 + r*(positive) > 0 ✓
        // Great, this is an M-matrix!

        // Left boundary (i=0): zero-flux, J_{-1/2} = 0
        // (c_0 - c_0^n)/dt = -J_{1/2}/dx  (only outgoing flux)
        // c_0 + r*[B(v_{1/2})*c_0 - B(-v_{1/2})*c_1] = c_0^n
        // c_0*(1 + r*B(v_{1/2})) - r*B(-v_{1/2})*c_1 = c_0^n

        b_coef[0] = 1.0 + r * B_v[0];
        c_coef[0] = -r * B_mv[0];
        rhs[0] = c[0];

        // Interior points
        for (int i = 1; i < n - 1; ++i) {
            a_coef[i] = -r * B_v[i-1];
            b_coef[i] = 1.0 + r * (B_v[i] + B_mv[i-1]);
            c_coef[i] = -r * B_mv[i];
            rhs[i] = c[i];
        }

        // Right boundary: Dirichlet c = c0
        b_coef[n-1] = 1.0;
        rhs[n-1] = params_.c0;

        solve_tridiag(a_coef, b_coef, c_coef, rhs, n);

        // Ensure positivity
        for (int i = 0; i < n; ++i) {
            c[i] = std::max(rhs[i], params_.c0 * 1e-12);
        }
    };

    // Save initial state
    {
        x_ = x_uniform;
        phi_ = phi;
        c_plus_ = c_plus;
        c_minus_ = c_minus;
        std::ostringstream ss;
        ss << output_dir << "/snapshot_" << std::setfill('0') << std::setw(5) << snapshot_count << ".dat";
        save_results(ss.str());
        time_file << snapshot_count << "\t" << 0.0 << "\n";
        snapshot_count++;
    }

    // Time stepping with Gummel iteration for Poisson-NP coupling
    for (int step = 0; step < n_steps; ++step) {
        double time = (step + 1) * dt;

        // Gummel iteration: iterate between Poisson and NP until convergence
        const int max_gummel_iter = 50;
        const double gummel_tol = 1e-6;

        for (int gummel_iter = 0; gummel_iter < max_gummel_iter; ++gummel_iter) {
            std::vector<double> c_plus_prev = c_plus;
            std::vector<double> c_minus_prev = c_minus;

            // Step 1: Solve Poisson equation for φ
            solve_poisson();

            // Step 2: Solve NP for c+ (cation, z=+1)
            solve_np_implicit(c_plus, params_.z_plus, params_.D_plus);

            // Step 3: Solve NP for c- (anion, z=-1)
            solve_np_implicit(c_minus, params_.z_minus, params_.D_minus);

            // Check Gummel convergence
            double max_change = 0.0;
            for (int i = 0; i < n; ++i) {
                max_change = std::max(max_change,
                    std::abs(c_plus[i] - c_plus_prev[i]) / (c_plus_prev[i] + 1e-30));
                max_change = std::max(max_change,
                    std::abs(c_minus[i] - c_minus_prev[i]) / (c_minus_prev[i] + 1e-30));
            }

            if (max_change < gummel_tol) {
                break;
            }
        }  // end Gummel iteration

        // Check for NaN
        bool has_nan = false;
        for (int i = 0; i < n; ++i) {
            if (std::isnan(c_plus[i]) || std::isnan(c_minus[i]) || std::isnan(phi[i])) {
                has_nan = true;
                break;
            }
        }
        if (has_nan) {
            std::cout << "  ERROR: NaN detected at step " << step << "!\n";
            break;
        }

        // Save snapshot
        if ((step + 1) % snapshot_interval == 0 || step == n_steps - 1) {
            x_ = x_uniform;
            phi_ = phi;
            c_plus_ = c_plus;
            c_minus_ = c_minus;

            std::ostringstream ss;
            ss << output_dir << "/snapshot_" << std::setfill('0') << std::setw(5) << snapshot_count << ".dat";
            save_results(ss.str());
            time_file << snapshot_count << "\t" << time * 1e9 << "\n";
            snapshot_count++;
        }

        // Progress
        if (step % std::max(1, n_steps / 20) == 0 || step == n_steps - 1) {
            double time_ns = time * 1e9;
            double c_plus_surface = c_plus[0] / params_.c0;
            double c_minus_surface = c_minus[0] / params_.c0;
            std::cout << "  t = " << std::fixed << std::setprecision(2) << time_ns
                      << " ns, c+/c0 = " << std::scientific << std::setprecision(2) << c_plus_surface
                      << ", c-/c0 = " << c_minus_surface << "\n";
        }
    }

    // Final update
    x_ = x_uniform;
    phi_ = phi;
    c_plus_ = c_plus;
    c_minus_ = c_minus;

    time_file.close();
    std::cout << "  Total snapshots: " << snapshot_count << "\n";
    std::cout << "  Scharfetter-Gummel solver completed.\n";
}

void PNPSolver1D::solve_transient_shenxu(double dt, double t_final,
                                          const std::string& output_dir,
                                          int snapshot_interval) {
    /**
     * Shen-Xu Positivity Preserving Scheme for PNP equations
     *
     * Based on: Shen & Xu (2021), Numer. Math. 148:671-697
     * "Unconditionally positivity preserving and energy dissipative schemes
     *  for Poisson-Nernst-Planck equations"
     *
     * The key insight is to reformulate the scheme as:
     *   (c^{n+1} - c^n)/δt = ∇·(D c^n ∇(log c^{n+1} + z*φ^{n+1}/φT))
     *
     * Setting w = log(c), the equation becomes:
     *   (exp(w^{n+1}) - exp(w^n))/δt = ∇·(D exp(w^n) ∇(w^{n+1} + z*φ^{n+1}/φT))
     *
     * This can be linearized by defining the effective diffusion coefficient:
     *   M^n = D * c^n = D * exp(w^n)
     *
     * For implementation, we use the exponential fitting approach:
     * Let ψ = z*φ/φT be the dimensionless potential.
     *
     * The flux at face i+1/2:
     *   J_{i+1/2} = -M_{i+1/2} * (∇w + z*∇φ/φT)
     *             = -M_{i+1/2} * ∇(w + ψ)
     *
     * Using Scharfetter-Gummel for the w variable:
     *   This reduces to standard SG discretization for c!
     *
     * The paper's innovation is treating M^n = D*c^n as fixed coefficient,
     * making the system for log(c^{n+1}) LINEAR but still guaranteeing positivity.
     */

    int n_steps = static_cast<int>(t_final / dt + 0.5);
    int n = params_.N;

    std::cout << "\n========================================\n";
    std::cout << "  Transient Solver (Shen-Xu Scheme)\n";
    std::cout << "========================================\n";
    std::cout << "  Time step: " << dt * 1e9 << " ns\n";
    std::cout << "  Total time: " << t_final * 1e9 << " ns\n";
    std::cout << "  Number of steps: " << n_steps << "\n";
    std::cout << "  Snapshot interval: " << snapshot_interval << " steps\n";

    // Use uniform grid
    double L = params_.L;
    double dx = L / (n - 1);
    std::vector<double> x_uniform(n);
    for (int i = 0; i < n; ++i) {
        x_uniform[i] = i * dx;
    }

    // Initialize: start from equilibrium state
    // First compute steady-state solution as initial condition
    std::vector<double> phi(n), c_plus(n), c_minus(n);

    // Initialize with linear potential profile
    for (int i = 0; i < n; ++i) {
        double xi = x_uniform[i] / L;
        phi[i] = params_.phi_left * (1.0 - xi);
        // Boltzmann distribution as initial guess
        c_plus[i] = params_.c0 * std::exp(-params_.z_plus * phi[i] / phi_T_);
        c_minus[i] = params_.c0 * std::exp(-params_.z_minus * phi[i] / phi_T_);
        // Clamp to reasonable values
        c_plus[i] = std::clamp(c_plus[i], params_.c0 * 1e-10, params_.c0 * 1e10);
        c_minus[i] = std::clamp(c_minus[i], params_.c0 * 1e-10, params_.c0 * 1e10);
    }

    // Log concentrations: w± = log(c±)
    std::vector<double> w_plus(n), w_minus(n);
    for (int i = 0; i < n; ++i) {
        w_plus[i] = std::log(c_plus[i]);
        w_minus[i] = std::log(c_minus[i]);
    }

    // Time info file
    std::ofstream time_file(output_dir + "/time_info.dat");
    time_file << "# snapshot_index  time_ns\n";

    int snapshot_count = 0;

    // Helper: Solve tridiagonal system Ax = d
    auto solve_tridiag = [](const std::vector<double>& a,
                            const std::vector<double>& b,
                            const std::vector<double>& c,
                            std::vector<double>& d,
                            int n) {
        std::vector<double> c_mod(n), d_mod(n);
        c_mod[0] = c[0] / b[0];
        d_mod[0] = d[0] / b[0];
        for (int i = 1; i < n; ++i) {
            double denom = b[i] - a[i] * c_mod[i-1];
            if (std::abs(denom) < 1e-30) {
                denom = (denom >= 0) ? 1e-30 : -1e-30;
            }
            if (i < n - 1) c_mod[i] = c[i] / denom;
            d_mod[i] = (d[i] - a[i] * d_mod[i-1]) / denom;
        }
        d[n-1] = d_mod[n-1];
        for (int i = n - 2; i >= 0; --i) {
            d[i] = d_mod[i] - c_mod[i] * d[i+1];
        }
    };

    // Coefficients for Poisson equation
    double e_charge = PhysicalConstants::e;
    double rho_coeff = e_charge * PhysicalConstants::NA / eps_;
    double dx2 = dx * dx;

    // Solve Poisson equation: ∇²φ = -ρ/ε
    auto solve_poisson = [&]() {
        const int max_iter = 100;
        const double tol = 1e-10;

        for (int iter = 0; iter < max_iter; ++iter) {
            double max_residual = 0.0;

            std::vector<double> a(n, 0.0), b(n, 0.0), c(n, 0.0), rhs(n, 0.0);

            // Left boundary: φ = φ_left
            b[0] = 1.0;
            rhs[0] = params_.phi_left - phi[0];

            // Interior points
            for (int i = 1; i < n - 1; ++i) {
                double laplacian = (phi[i+1] - 2.0*phi[i] + phi[i-1]) / dx2;
                double rho = rho_coeff * (params_.z_plus * c_plus[i] + params_.z_minus * c_minus[i]);
                double residual = laplacian + rho;
                max_residual = std::max(max_residual, std::abs(residual));

                a[i] = 1.0 / dx2;
                b[i] = -2.0 / dx2;
                c[i] = 1.0 / dx2;
                rhs[i] = -residual;
            }

            // Right boundary: φ = 0
            b[n-1] = 1.0;
            rhs[n-1] = params_.phi_right - phi[n-1];

            if (max_residual < tol) break;

            solve_tridiag(a, b, c, rhs, n);
            for (int i = 0; i < n; ++i) {
                phi[i] += rhs[i];
            }
        }
    };

    /**
     * Solve NP equation using modified Shen-Xu scheme.
     *
     * The Shen-Xu scheme: (c^{n+1} - c^n)/δt = ∇·(D c^n ∇(log c^{n+1} + z*φ/φT))
     *
     * Let w = log(c), ψ = z*φ/φT. Then c = exp(w), and:
     *   exp(w^{n+1}) - exp(w^n) = δt * ∇·(D exp(w^n) ∇(w^{n+1} + ψ))
     *
     * Define M^n_{i+1/2} = D * (c^n_i + c^n_{i+1})/2 (arithmetic mean)
     *
     * The discretized equation at interior point i:
     *   exp(w_i^{n+1}) - exp(w_i^n) = (δt/dx²) * [M_{i+1/2}((w_{i+1}+ψ_{i+1}) - (w_i+ψ_i))
     *                                           - M_{i-1/2}((w_i+ψ_i) - (w_{i-1}+ψ_{i-1}))]
     *
     * This is a nonlinear equation for w^{n+1}. We solve it using Newton iteration.
     */
    auto solve_np_shenxu = [&](std::vector<double>& c, std::vector<double>& w,
                                const std::vector<double>& c_old, int z, double D) {
        const int max_newton = 50;
        const double newton_tol = 1e-8;

        // Mobility coefficient M_{i+1/2} = D * (c_old_i + c_old_{i+1})/2
        std::vector<double> M_face(n-1);
        for (int i = 0; i < n - 1; ++i) {
            M_face[i] = D * 0.5 * (c_old[i] + c_old[i+1]);
        }

        double r = dt / dx2;

        // Dimensionless potential: ψ = z * φ / φT
        std::vector<double> psi(n);
        for (int i = 0; i < n; ++i) {
            psi[i] = z * phi[i] / phi_T_;
        }

        // Safe exp function to prevent overflow
        auto safe_exp = [](double x) -> double {
            const double MAX_EXP = 50.0;
            if (x > MAX_EXP) return std::exp(MAX_EXP);
            if (x < -MAX_EXP) return std::exp(-MAX_EXP);
            return std::exp(x);
        };

        for (int newton_iter = 0; newton_iter < max_newton; ++newton_iter) {
            std::vector<double> F(n, 0.0);  // Residual
            std::vector<double> a(n, 0.0), b(n, 0.0), cc(n, 0.0);  // Jacobian

            // Left boundary (i=0): zero flux
            // At zero-flux boundary, the outward flux is zero:
            // M_{1/2} * ((w_1+ψ_1) - (w_0+ψ_0)) / dx = 0
            // This means w_0 + ψ_0 = w_1 + ψ_1 for zero flux
            // But for the mass balance, we use the original equation with J_{-1/2}=0:
            // exp(w_0^{n+1}) - exp(w_0^n) = (δt/dx²) * M_{1/2}((w_1+ψ_1) - (w_0+ψ_0))
            {
                double exp_w = safe_exp(w[0]);
                double exp_w_old = safe_exp(std::log(c_old[0]));
                double flux_right = M_face[0] * ((w[1] + psi[1]) - (w[0] + psi[0]));
                F[0] = exp_w - exp_w_old - r * flux_right;

                // Jacobian: ∂F/∂w_0 = exp(w_0) + r * M_{1/2}
                //           ∂F/∂w_1 = -r * M_{1/2}
                b[0] = exp_w + r * M_face[0];
                cc[0] = -r * M_face[0];
            }

            // Interior points
            for (int i = 1; i < n - 1; ++i) {
                double exp_w = safe_exp(w[i]);
                double exp_w_old = safe_exp(std::log(c_old[i]));
                double flux_right = M_face[i] * ((w[i+1] + psi[i+1]) - (w[i] + psi[i]));
                double flux_left = M_face[i-1] * ((w[i] + psi[i]) - (w[i-1] + psi[i-1]));
                F[i] = exp_w - exp_w_old - r * (flux_right - flux_left);

                // Jacobian
                a[i] = -r * M_face[i-1];
                b[i] = exp_w + r * (M_face[i-1] + M_face[i]);
                cc[i] = -r * M_face[i];
            }

            // Right boundary: Dirichlet c = c0
            // => w = log(c0)
            {
                double w_bc = std::log(params_.c0);
                F[n-1] = w[n-1] - w_bc;
                b[n-1] = 1.0;
            }

            // Check convergence
            double max_F = 0.0;
            for (int i = 0; i < n; ++i) {
                max_F = std::max(max_F, std::abs(F[i]));
            }
            if (max_F < newton_tol) break;

            // Solve J * δw = -F
            for (int i = 0; i < n; ++i) {
                F[i] = -F[i];
            }
            solve_tridiag(a, b, cc, F, n);

            // Damped update with line search
            double alpha = 1.0;
            for (int i = 0; i < n; ++i) {
                // Limit step size to prevent too large changes
                double max_step = 2.0;  // Maximum change in log(c) per iteration
                if (std::abs(F[i]) > max_step) {
                    F[i] = (F[i] > 0) ? max_step : -max_step;
                }
                w[i] += alpha * F[i];
            }
        }

        // Recover c from w: c = exp(w)
        for (int i = 0; i < n; ++i) {
            c[i] = safe_exp(w[i]);
            c[i] = std::max(c[i], params_.c0 * 1e-15);
        }
    };

    // Save initial state
    {
        x_ = x_uniform;
        phi_ = phi;
        c_plus_ = c_plus;
        c_minus_ = c_minus;
        std::ostringstream ss;
        ss << output_dir << "/snapshot_" << std::setfill('0') << std::setw(5) << snapshot_count << ".dat";
        save_results(ss.str());
        time_file << snapshot_count << "\t" << 0.0 << "\n";
        snapshot_count++;
    }

    // Time stepping
    for (int step = 0; step < n_steps; ++step) {
        double time = (step + 1) * dt;

        // Store old concentrations for Shen-Xu scheme
        std::vector<double> c_plus_old = c_plus;
        std::vector<double> c_minus_old = c_minus;

        // Gummel-like iteration for Poisson-NP coupling
        const int max_gummel_iter = 20;
        const double gummel_tol = 1e-6;

        for (int gummel_iter = 0; gummel_iter < max_gummel_iter; ++gummel_iter) {
            std::vector<double> c_plus_prev = c_plus;
            std::vector<double> c_minus_prev = c_minus;

            // Step 1: Solve Poisson equation for φ
            solve_poisson();

            // Step 2: Solve NP for c+ using Shen-Xu scheme
            // Note: c_old is from the previous TIME step, not Gummel iteration
            solve_np_shenxu(c_plus, w_plus, c_plus_old, params_.z_plus, params_.D_plus);

            // Step 3: Solve NP for c- using Shen-Xu scheme
            solve_np_shenxu(c_minus, w_minus, c_minus_old, params_.z_minus, params_.D_minus);

            // Check Gummel convergence
            double max_change = 0.0;
            for (int i = 0; i < n; ++i) {
                max_change = std::max(max_change,
                    std::abs(c_plus[i] - c_plus_prev[i]) / (c_plus_prev[i] + 1e-30));
                max_change = std::max(max_change,
                    std::abs(c_minus[i] - c_minus_prev[i]) / (c_minus_prev[i] + 1e-30));
            }

            if (max_change < gummel_tol) {
                break;
            }
        }

        // Check for NaN or negative concentrations
        bool has_issue = false;
        for (int i = 0; i < n; ++i) {
            if (std::isnan(c_plus[i]) || std::isnan(c_minus[i]) || std::isnan(phi[i]) ||
                c_plus[i] <= 0.0 || c_minus[i] <= 0.0) {
                has_issue = true;
                break;
            }
        }
        if (has_issue) {
            std::cout << "  ERROR: Invalid values detected at step " << step << "!\n";
            break;
        }

        // Save snapshot
        if ((step + 1) % snapshot_interval == 0 || step == n_steps - 1) {
            x_ = x_uniform;
            phi_ = phi;
            c_plus_ = c_plus;
            c_minus_ = c_minus;

            std::ostringstream ss;
            ss << output_dir << "/snapshot_" << std::setfill('0') << std::setw(5) << snapshot_count << ".dat";
            save_results(ss.str());
            time_file << snapshot_count << "\t" << time * 1e9 << "\n";
            snapshot_count++;
        }

        // Progress
        if (step % std::max(1, n_steps / 20) == 0 || step == n_steps - 1) {
            double time_ns = time * 1e9;
            double c_plus_surface = c_plus[0] / params_.c0;
            double c_minus_surface = c_minus[0] / params_.c0;
            std::cout << "  t = " << std::fixed << std::setprecision(2) << time_ns
                      << " ns, c+/c0 = " << std::scientific << std::setprecision(3) << c_plus_surface
                      << ", c-/c0 = " << c_minus_surface << "\n";
        }
    }

    // Final update
    x_ = x_uniform;
    phi_ = phi;
    c_plus_ = c_plus;
    c_minus_ = c_minus;

    time_file.close();
    std::cout << "  Total snapshots: " << snapshot_count << "\n";
    std::cout << "  Shen-Xu solver completed.\n";
}

double PNPSolver1D::get_debye_length() const {
    return lambda_D_;
}

std::vector<double> PNPSolver1D::get_rho() const {
    std::vector<double> rho(params_.N);
    double e = PhysicalConstants::e;
    double NA = PhysicalConstants::NA;

    for (int i = 0; i < params_.N; ++i) {
        rho[i] = e * NA * (params_.z_plus * c_plus_[i] + params_.z_minus * c_minus_[i]);
    }
    return rho;
}

std::vector<double> PNPSolver1D::gouy_chapman_solution(double phi0) const {
    std::vector<double> phi_gc(params_.N);
    double gamma0 = std::tanh(phi0 / (4.0 * phi_T_));

    for (int i = 0; i < params_.N; ++i) {
        double gamma = gamma0 * std::exp(-x_[i] / lambda_D_);
        gamma = std::clamp(gamma, -0.9999, 0.9999);
        phi_gc[i] = 4.0 * phi_T_ * std::atanh(gamma);
    }

    return phi_gc;
}

void PNPSolver1D::save_results(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }

    double e = PhysicalConstants::e;
    double NA = PhysicalConstants::NA;

    file << "# 1D PNP Solver Results\n";
    file << "# Debye length: " << lambda_D_ * 1e9 << " nm\n";
    file << "# Surface potential: " << params_.phi_left * 1000.0 << " mV\n";
    file << "# Bulk concentration: " << params_.c0 << " mol/m^3\n";
    file << "# Thermal voltage: " << phi_T_ * 1000.0 << " mV\n";
    file << "# Grid: non-uniform with stretching factor " << params_.grid_stretch << "\n";
    file << "#\n";
    file << "# Columns:\n";
    file << "# 1: x [nm]\n";
    file << "# 2: x/lambda_D [-]\n";
    file << "# 3: phi [mV]\n";
    file << "# 4: phi_normalized = e*phi/(kB*T) [-]\n";
    file << "# 5: c+ [mol/m^3]\n";
    file << "# 6: c- [mol/m^3]\n";
    file << "# 7: c+/c0 [-]\n";
    file << "# 8: c-/c0 [-]\n";
    file << "# 9: rho [C/m^3]\n";
    file << "# 10: phi_GC [mV] (Gouy-Chapman analytical)\n";

    file << std::scientific << std::setprecision(8);

    std::vector<double> phi_gc = gouy_chapman_solution(params_.phi_left);

    for (int i = 0; i < params_.N; ++i) {
        double x_nm = x_[i] * 1e9;
        double x_norm = x_[i] / lambda_D_;
        double phi_mV = phi_[i] * 1000.0;
        double phi_norm = phi_[i] / phi_T_;
        double rho = e * NA * (params_.z_plus * c_plus_[i] + params_.z_minus * c_minus_[i]);
        double phi_gc_mV = phi_gc[i] * 1000.0;

        file << x_nm << "\t"
             << x_norm << "\t"
             << phi_mV << "\t"
             << phi_norm << "\t"
             << c_plus_[i] << "\t"
             << c_minus_[i] << "\t"
             << c_plus_[i] / params_.c0 << "\t"
             << c_minus_[i] / params_.c0 << "\t"
             << rho << "\t"
             << phi_gc_mV << "\n";
    }

    file.close();
    std::cout << "Results saved to " << filename << "\n";
}

double PNPSolver1D::compute_L2_error() const {
    // L2 error = sqrt(mean((phi - phi_GC)^2))
    std::vector<double> phi_gc = gouy_chapman_solution(params_.phi_left);

    double sum_sq = 0.0;
    for (int i = 0; i < params_.N; ++i) {
        double diff = phi_[i] - phi_gc[i];
        sum_sq += diff * diff;
    }

    return std::sqrt(sum_sq / params_.N);
}

double PNPSolver1D::compute_relative_L2_error() const {
    // Relative L2 error = L2 / sqrt(mean(phi_GC^2))
    std::vector<double> phi_gc = gouy_chapman_solution(params_.phi_left);

    double sum_sq_error = 0.0;
    double sum_sq_ref = 0.0;

    for (int i = 0; i < params_.N; ++i) {
        double diff = phi_[i] - phi_gc[i];
        sum_sq_error += diff * diff;
        sum_sq_ref += phi_gc[i] * phi_gc[i];
    }

    if (sum_sq_ref < 1e-30) {
        return 0.0;
    }

    return std::sqrt(sum_sq_error / sum_sq_ref);
}

double PNPSolver1D::compute_Linf_error() const {
    // L-infinity error = max(|phi - phi_GC|)
    std::vector<double> phi_gc = gouy_chapman_solution(params_.phi_left);

    double max_error = 0.0;
    for (int i = 0; i < params_.N; ++i) {
        max_error = std::max(max_error, std::abs(phi_[i] - phi_gc[i]));
    }

    return max_error;
}

} // namespace pnp
