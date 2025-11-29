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

    std::cout << "========================================\n";
    std::cout << "  1D PNP Solver Initialization\n";
    std::cout << "========================================\n";
    std::cout << "  Debye length: " << lambda_D_ * 1e9 << " nm\n";
    std::cout << "  Thermal voltage: " << phi_T_ * 1000.0 << " mV\n";
    std::cout << "  Grid points: " << params_.N << "\n";
    std::cout << "  Grid stretching: " << params_.grid_stretch << "\n";

    initialize();
}

void PNPSolver1D::create_nonuniform_grid() {
    // Create non-uniform grid with clustering near x=0 (interface)
    // Using transformation: x = L * [1 - (1 - ξ)^β] where ξ ∈ [0, 1]
    // β > 1 gives finer spacing near x=0

    int n = params_.N;
    double L = params_.L;
    double beta = params_.grid_stretch;

    x_.resize(n);
    dx_.resize(n);

    for (int i = 0; i < n; ++i) {
        double xi = static_cast<double>(i) / (n - 1);  // ξ ∈ [0, 1]
        x_[i] = L * (1.0 - std::pow(1.0 - xi, beta));
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

    // Count points within 5 Debye lengths
    int points_in_edl = 0;
    for (int i = 0; i < n; ++i) {
        if (x_[i] < 5.0 * lambda_D_) points_in_edl++;
    }
    std::cout << "  Points within 5*λ_D: " << points_in_edl << "\n";
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
        c_plus_[i] = params_.c0 * std::exp(-params_.z_plus * phi_norm);
        c_minus_[i] = params_.c0 * std::exp(-params_.z_minus * phi_norm);
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
        // Build Jacobian and residual for Newton step on non-uniform grid
        // F(φ) = d²φ/dx² - κ²*φ_T*sinh(φ/φ_T) = 0

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

            // Jacobian entries for non-uniform grid
            // d²φ/dx² ≈ [φ(i+1) - φ(i)]/dx+ - [φ(i) - φ(i-1)]/dx- ] / dx_avg
            a[i] = 1.0 / (dx_minus * dx_avg);                    // coeff for φ[i-1]
            c[i] = 1.0 / (dx_plus * dx_avg);                     // coeff for φ[i+1]
            b[i] = -1.0 / (dx_minus * dx_avg) - 1.0 / (dx_plus * dx_avg) - kappa2 * cosh_val;  // coeff for φ[i]

            // Residual
            double laplacian = (phi_[i + 1] - phi_[i]) / dx_plus - (phi_[i] - phi_[i - 1]) / dx_minus;
            laplacian /= dx_avg;

            rhs[i] = -(laplacian - kappa2 * phi_T_ * sinh_val);
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
    int n_steps = static_cast<int>(t_final / dt);
    std::cout << "Starting transient solver: " << n_steps << " time steps\n";

    for (int step = 0; step < n_steps; ++step) {
        solve();
        if (step % 1000 == 0) {
            std::cout << "Time step " << step << " / " << n_steps << "\n";
        }
    }
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

} // namespace pnp
