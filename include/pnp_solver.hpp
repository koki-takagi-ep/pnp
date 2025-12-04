/**
 * @file pnp_solver.hpp
 * @brief 1D Poisson-Nernst-Planck equation solver for ionic liquids
 *
 * This solver computes the electric double layer (EDL) structure
 * by solving the coupled Poisson and Nernst-Planck equations.
 *
 * Features:
 * - Non-uniform grid with clustering near interface
 * - Newton-Raphson solver for Poisson-Boltzmann equation
 * - Gouy-Chapman analytical solution for validation
 */

#ifndef PNP_SOLVER_HPP
#define PNP_SOLVER_HPP

#include <vector>
#include <string>
#include <cmath>
#include <utility>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <algorithm>

namespace pnp {

/**
 * @brief Physical constants (SI units)
 */
struct PhysicalConstants {
    static constexpr double e = 1.602176634e-19;      // Elementary charge [C]
    static constexpr double kB = 1.380649e-23;        // Boltzmann constant [J/K]
    static constexpr double eps0 = 8.8541878128e-12;  // Vacuum permittivity [F/m]
    static constexpr double NA = 6.02214076e23;       // Avogadro number [1/mol]
};

/**
 * @brief Model type for the solver
 */
enum class ModelType {
    STANDARD_PB,     // Standard Poisson-Boltzmann
    BIKERMAN,        // Modified PB with steric effects (Bikerman model)
    MODIFIED_PF      // Modified Poisson-Fermi with electrostatic correlations (Bazant 2011)
};

/**
 * @brief Parameters for the PNP solver
 */
struct PNPParameters {
    // Domain parameters
    double L;           // Domain length [m]
    int N;              // Number of grid points

    // Physical parameters
    double T;           // Temperature [K]
    double eps_r;       // Relative permittivity
    double c0;          // Bulk concentration [mol/m^3]
    double D_plus;      // Diffusivity of cation [m^2/s]
    double D_minus;     // Diffusivity of anion [m^2/s]
    int z_plus;         // Valence of cation
    int z_minus;        // Valence of anion

    // Boundary conditions
    double phi_left;    // Potential at left boundary [V]
    double phi_right;   // Potential at right boundary [V]

    // Solver parameters
    double tol;         // Convergence tolerance
    int max_iter;       // Maximum iterations
    double omega;       // Relaxation factor

    // Grid stretching parameter (higher = more clustering near interface)
    double grid_stretch;

    // Boundary condition type
    bool closed_system; // If true, use zero-flux BC at both ends (capacitor)
                        // If false, use Dirichlet c=c0 at right (open system)

    // Model selection
    ModelType model;    // Standard PB or Modified (Bikerman)

    // Ion size parameters for Bikerman model
    double a;           // Ion diameter [m]

    // Correlation length for Modified Poisson-Fermi model (Bazant 2011)
    double l_c;         // Electrostatic correlation length [m]
    double delta_c;     // Dimensionless correlation length l_c / lambda_D

    // Default constructor with typical values for ionic liquids
    PNPParameters()
        : L(100e-9)       // 100 nm domain
        , N(501)          // Grid points
        , T(298.15)       // Room temperature
        , eps_r(12.0)     // Typical for ionic liquids
        , c0(1000.0)      // 1 M = 1000 mol/m^3
        , D_plus(1e-10)   // Typical for ionic liquids
        , D_minus(1e-10)
        , z_plus(1)
        , z_minus(-1)
        , phi_left(0.1)   // 100 mV surface potential
        , phi_right(0.0)
        , tol(1e-10)
        , max_iter(1000)
        , omega(0.3)
        , grid_stretch(3.0)  // Stretching factor for non-uniform grid
        , closed_system(false) // Default: open system with Dirichlet BC at right
        , model(ModelType::STANDARD_PB)
        , a(0.7e-9)       // ~0.7 nm typical ion diameter for ionic liquids
        , l_c(0.0)        // Correlation length (0 = no correlation, Bikerman limit)
        , delta_c(10.0)   // Default dimensionless correlation length (Bazant 2011)
    {}
};

/**
 * @brief 1D Poisson-Nernst-Planck solver
 *
 * Solves the Poisson-Boltzmann equation using Newton-Raphson method
 * on a non-uniform grid with clustering near the interface.
 */
class PNPSolver1D {
public:
    explicit PNPSolver1D(const PNPParameters& params);

    bool solve();
    void solve_transient(double dt, double t_final);
    void solve_transient_with_snapshots(double dt, double t_final,
                                         const std::string& output_dir,
                                         int snapshot_interval = 10);
    void solve_transient_gummel(double dt, double t_final,
                                const std::string& output_dir,
                                int snapshot_interval = 10);
    void solve_transient_continuation(int n_steps,
                                      const std::string& output_dir);
    void solve_transient_newton(double dt, double t_final,
                                const std::string& output_dir,
                                int snapshot_interval = 10);

    /**
     * @brief Transient solver using Slotboom transformation
     *
     * Uses Slotboom variables u± = c± * exp(±φ/φT) to transform the NP equations
     * into self-adjoint diffusion form, which is more stable and preserves positivity.
     *
     * References:
     * - Slotboom (1969), Electronics Letters
     * - Liu & Wang (2021), Numerische Mathematik
     */
    void solve_transient_slotboom(double dt, double t_final,
                                  const std::string& output_dir,
                                  int snapshot_interval = 10);

    /**
     * @brief Transient solver using Shen-Xu positivity preserving scheme
     *
     * Implements the unconditionally positivity preserving and energy dissipative
     * scheme from Shen & Xu (2021). The key idea is to treat c^n explicitly as
     * a coefficient and solve for log(c^{n+1}) implicitly.
     *
     * Scheme: (c^{n+1} - c^n)/δt = ∇·(D c^n ∇(log c^{n+1} + z*φ^{n+1}/φT))
     *
     * This leads to solving a strictly convex minimization problem at each step,
     * guaranteeing unique solution and positivity regardless of timestep size.
     *
     * References:
     * - Shen & Xu (2021), Numer. Math. 148:671-697
     */
    void solve_transient_shenxu(double dt, double t_final,
                                const std::string& output_dir,
                                int snapshot_interval = 10);

    /**
     * @brief Stable transient solver using electric field formulation
     *
     * Uses electric field E as primary variable instead of potential φ.
     * This approach, inspired by PoNPs (Toyoura & Ueno, Kyoto University),
     * treats Poisson equation as first-order: dE/dx = ρ/ε
     *
     * Features:
     * - Backward Euler implicit time integration
     * - Simple arithmetic mean flux (stable for all Peclet numbers)
     * - Newton-Raphson nonlinear solver
     * - Adaptive time stepping
     * - Explicit charge conservation monitoring
     *
     * @param dt_init Initial time step [s]
     * @param t_final Final simulation time [s]
     * @param output_dir Directory for output files
     * @param snapshot_interval Steps between snapshots
     */
    void solve_transient_efield(double dt_init, double t_final,
                                const std::string& output_dir,
                                int snapshot_interval = 10);

    const std::vector<double>& get_x() const { return x_; }
    const std::vector<double>& get_phi() const { return phi_; }
    const std::vector<double>& get_c_plus() const { return c_plus_; }
    const std::vector<double>& get_c_minus() const { return c_minus_; }

    std::vector<double> get_rho() const;
    double get_debye_length() const;
    std::vector<double> gouy_chapman_solution(double phi0) const;
    void save_results(const std::string& filename) const;
    const std::vector<double>& get_residual_history() const { return residual_history_; }

    /**
     * @brief Compute L2 error norm against Gouy-Chapman solution
     * @return L2 norm error in volts
     */
    double compute_L2_error() const;

    /**
     * @brief Compute relative L2 error norm
     * @return Relative L2 error (dimensionless)
     */
    double compute_relative_L2_error() const;

    /**
     * @brief Compute L-infinity (max) error
     * @return Maximum absolute error in volts
     */
    double compute_Linf_error() const;

    /**
     * @brief Get the packing fraction for Bikerman model
     */
    double get_packing_fraction() const { return nu_; }

    /**
     * @brief Get self-consistently determined bulk potential (for closed system)
     */
    double get_bulk_potential() const { return phi_bulk_; }

    /**
     * @brief Get model name as string
     */
    std::string get_model_name() const;

    /**
     * @brief Compute surface charge density at electrodes using Gauss's law
     * @return pair of (sigma_left, sigma_right) in C/m^2
     *
     * Uses σ = -ε₀εᵣ(dφ/dx) at electrode surfaces with 2nd-order FD
     * For closed system: sigma_left + sigma_right ≈ 0 (electroneutrality)
     */
    std::pair<double, double> compute_surface_charge() const;

    /**
     * @brief Compute differential capacitance of each EDL
     * @return pair of (C_left, C_right) in F/m^2
     *
     * C = |σ| / |Δφ_EDL| where Δφ_EDL is potential drop across EDL
     */
    std::pair<double, double> compute_capacitance() const;

    /**
     * @brief Compute total capacitance of the capacitor (two EDLs in series)
     * @return Total capacitance in F/m^2
     */
    double compute_total_capacitance() const;

private:
    PNPParameters params_;

    // Grid and solution vectors
    std::vector<double> x_;        // Spatial grid [m]
    std::vector<double> dx_;       // Local grid spacing [m] (for non-uniform grid)
    std::vector<double> phi_;      // Electric potential [V]
    std::vector<double> c_plus_;   // Cation concentration [mol/m^3]
    std::vector<double> c_minus_;  // Anion concentration [mol/m^3]

    // Convergence tracking
    std::vector<double> residual_history_;

    // Derived quantities
    double beta_;      // 1/(kB*T) [1/J]
    double eps_;       // Permittivity [F/m]
    double lambda_D_;  // Debye length [m]
    double phi_T_;     // Thermal voltage [V]
    double nu_;        // Packing fraction for Bikerman model
    double phi_bulk_;  // Self-consistently determined bulk potential [V] (for closed system)
    double l_c_;       // Correlation length [m] for Modified Poisson-Fermi

    void initialize();

    /**
     * @brief Solve Modified Poisson-Fermi equation (4th order BVP)
     *
     * Solves: (δ_c² ∇² - 1) ∇² φ̃ = -ρ̃(φ̃)
     * where ρ̃ = sinh(φ̃) / [1 + 2γ sinh²(φ̃/2)]
     *
     * The 4th order ODE is converted to a first-order system:
     *   y₁ = φ, y₂ = φ', y₃ = φ'', y₄ = φ'''
     *
     * Boundary conditions:
     *   y₁(0) = V (surface potential)
     *   y₄(0) = 0 (flat charge density at surface)
     *   y₁(∞) = 0 (bulk)
     *   y₂(∞) = 0 (bulk)
     *
     * Reference: Bazant, Storey, Kornyshev, PRL 106, 046102 (2011)
     */
    bool solve_modified_pf();
    void create_nonuniform_grid();
    void update_concentrations_from_phi();
    void solve_poisson();
    void update_concentrations();
    double compute_residual() const;

    std::vector<double> solve_tridiagonal(
        const std::vector<double>& a,
        const std::vector<double>& b,
        const std::vector<double>& c,
        const std::vector<double>& d
    ) const;

    static double bernoulli(double x);
};

} // namespace pnp

#endif // PNP_SOLVER_HPP
