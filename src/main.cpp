/**
 * @file main.cpp
 * @brief Main program for 1D Poisson-Nernst-Planck solver
 *
 * Simulates the electric double layer (EDL) in ionic liquids
 * using the PNP equations.
 */

#include "pnp_solver.hpp"
#include <iostream>
#include <cstdlib>
#include <string>
#include <filesystem>

void print_usage(const char* program_name) {
    std::cout << "Usage: " << program_name << " [options]\n"
              << "\nOptions:\n"
              << "  --phi0 <value>      Left electrode potential in mV (default: 100)\n"
              << "  --phi-right <value> Right electrode potential in mV (default: 0)\n"
              << "  --c0 <value>        Bulk concentration in mol/L (default: 1.0)\n"
              << "  --eps <value>       Relative permittivity (default: 12)\n"
              << "  --L <value>         Domain length in nm (default: 50)\n"
              << "  --N <value>         Number of grid points (default: 1001)\n"
              << "  --output <file>     Output filename (default: results/pnp_results.dat)\n"
              << "  --model <type>      Model type: standard, bikerman, or mpf (default: standard)\n"
              << "  --ion-size <value>  Ion diameter in nm for Bikerman/MPF model (default: 0.7)\n"
              << "  --delta-c <value>   Dimensionless correlation length for MPF model (default: 10)\n"
              << "  --stretch <value>   Grid stretching factor (default: 3.0, use 1.0 for uniform)\n"
              << "  --closed-system     Use zero-flux BC at both ends (capacitor model)\n"
              << "  --transient         Run transient simulation instead of steady-state\n"
              << "  --dt <value>        Time step in ns for transient (default: 0.1)\n"
              << "  --t-final <value>   Final time in µs for transient (default: 1.0)\n"
              << "  --animation         Save snapshots for GIF animation (transient mode)\n"
              << "  --snapshot-dir <dir> Directory for animation snapshots (default: results/snapshots)\n"
              << "  --snapshot-interval <n> Save snapshot every n steps (default: 10)\n"
              << "  --gummel            Use Gummel iteration for transient (more stable)\n"
              << "  --continuation <n>  Use continuation method with n steps (most stable)\n"
              << "  --newton            Use fully implicit Newton method for transient\n"
              << "  --slotboom          Use Slotboom transformation for transient (stable)\n"
              << "  --shenxu            Use Shen-Xu positivity preserving scheme (most stable)\n"
              << "  --efield            Use E-field formulation (inspired by PoNPs, stable)\n"
              << "  --help              Show this help message\n"
              << std::endl;
}

int main(int argc, char* argv[]) {
    // Default parameters
    double phi0_mV = 100.0;       // Left electrode potential [mV]
    double phi_right_mV = 0.0;    // Right electrode potential [mV]
    double c0_molL = 1.0;         // Bulk concentration [mol/L]
    double eps_r = 12.0;          // Relative permittivity
    double L_nm = 50.0;           // Domain length [nm]
    int N = 1001;                 // Grid points
    std::string output_file = "results/pnp_results.dat";
    std::string model_type = "standard";  // Model type: standard or bikerman
    double ion_size_nm = 0.7;     // Ion diameter for Bikerman model [nm]
    double delta_c = 10.0;        // Dimensionless correlation length for MPF model
    double grid_stretch = 3.0;    // Grid stretching factor
    bool closed_system = false;   // Use zero-flux BC at both ends
    bool run_transient = false;   // Run transient simulation
    double dt_ns = 0.1;           // Time step [ns]
    double t_final_us = 1.0;      // Final time [µs]
    bool save_animation = false;  // Save snapshots for animation
    std::string snapshot_dir = "results/snapshots";  // Snapshot output directory
    int snapshot_interval = 10;   // Save every n steps
    bool use_gummel = false;      // Use Gummel iteration for transient
    bool use_continuation = false; // Use continuation method
    int continuation_steps = 50;   // Number of continuation steps
    bool use_newton = false;       // Use fully implicit Newton method
    bool use_slotboom = false;     // Use Slotboom transformation method
    bool use_shenxu = false;       // Use Shen-Xu positivity preserving scheme
    bool use_efield = false;       // Use E-field formulation (PoNPs-inspired)

    // Parse command line arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help") {
            print_usage(argv[0]);
            return 0;
        } else if (arg == "--phi0" && i + 1 < argc) {
            phi0_mV = std::atof(argv[++i]);
        } else if (arg == "--phi-right" && i + 1 < argc) {
            phi_right_mV = std::atof(argv[++i]);
        } else if (arg == "--c0" && i + 1 < argc) {
            c0_molL = std::atof(argv[++i]);
        } else if (arg == "--eps" && i + 1 < argc) {
            eps_r = std::atof(argv[++i]);
        } else if (arg == "--L" && i + 1 < argc) {
            L_nm = std::atof(argv[++i]);
        } else if (arg == "--N" && i + 1 < argc) {
            N = std::atoi(argv[++i]);
        } else if (arg == "--output" && i + 1 < argc) {
            output_file = argv[++i];
        } else if (arg == "--model" && i + 1 < argc) {
            model_type = argv[++i];
        } else if (arg == "--ion-size" && i + 1 < argc) {
            ion_size_nm = std::atof(argv[++i]);
        } else if (arg == "--delta-c" && i + 1 < argc) {
            delta_c = std::atof(argv[++i]);
        } else if (arg == "--stretch" && i + 1 < argc) {
            grid_stretch = std::atof(argv[++i]);
        } else if (arg == "--closed-system") {
            closed_system = true;
        } else if (arg == "--transient") {
            run_transient = true;
        } else if (arg == "--dt" && i + 1 < argc) {
            dt_ns = std::atof(argv[++i]);
        } else if (arg == "--t-final" && i + 1 < argc) {
            t_final_us = std::atof(argv[++i]);
        } else if (arg == "--animation") {
            save_animation = true;
            run_transient = true;  // Animation requires transient mode
        } else if (arg == "--snapshot-dir" && i + 1 < argc) {
            snapshot_dir = argv[++i];
        } else if (arg == "--snapshot-interval" && i + 1 < argc) {
            snapshot_interval = std::atoi(argv[++i]);
        } else if (arg == "--gummel") {
            use_gummel = true;
            run_transient = true;
        } else if (arg == "--continuation" && i + 1 < argc) {
            use_continuation = true;
            continuation_steps = std::atoi(argv[++i]);
            run_transient = true;
        } else if (arg == "--newton") {
            use_newton = true;
            run_transient = true;
        } else if (arg == "--slotboom") {
            use_slotboom = true;
            run_transient = true;
        } else if (arg == "--shenxu") {
            use_shenxu = true;
            run_transient = true;
        } else if (arg == "--efield") {
            use_efield = true;
            run_transient = true;
        }
    }

    std::cout << "========================================\n";
    std::cout << "  1D Poisson-Nernst-Planck Solver\n";
    std::cout << "  for Ionic Liquid Electric Double Layer\n";
    std::cout << "========================================\n\n";

    // Setup parameters
    pnp::PNPParameters params;
    params.phi_left = phi0_mV * 1e-3;      // Convert mV to V
    params.phi_right = phi_right_mV * 1e-3; // Convert mV to V
    params.c0 = c0_molL * 1000.0;          // Convert mol/L to mol/m^3
    params.eps_r = eps_r;
    params.L = L_nm * 1e-9;                // Convert nm to m
    params.N = N;
    params.a = ion_size_nm * 1e-9;         // Convert nm to m
    params.grid_stretch = grid_stretch;
    params.closed_system = closed_system;

    // Set model type
    if (model_type == "bikerman") {
        params.model = pnp::ModelType::BIKERMAN;
    } else if (model_type == "mpf" || model_type == "modified-pf") {
        params.model = pnp::ModelType::MODIFIED_PF;
        params.delta_c = delta_c;
    } else {
        params.model = pnp::ModelType::STANDARD_PB;
    }

    std::cout << "Parameters:\n";
    std::cout << "  Left electrode potential: " << phi0_mV << " mV\n";
    std::cout << "  Right electrode potential: " << phi_right_mV << " mV\n";
    std::cout << "  Bulk concentration: " << c0_molL << " mol/L\n";
    std::cout << "  Relative permittivity: " << eps_r << "\n";
    std::cout << "  Domain length: " << L_nm << " nm\n";
    std::cout << "  Grid points: " << N << "\n";
    std::cout << "  Model: " << model_type << "\n";
    if (model_type == "bikerman") {
        std::cout << "  Ion size: " << ion_size_nm << " nm\n";
    } else if (model_type == "mpf" || model_type == "modified-pf") {
        std::cout << "  Ion size: " << ion_size_nm << " nm\n";
        std::cout << "  Correlation length (delta_c): " << delta_c << "\n";
    }
    if (closed_system) {
        std::cout << "  Boundary: closed system (zero-flux at both ends)\n";
    } else {
        std::cout << "  Boundary: open system (Dirichlet c=c0 at right)\n";
    }
    if (run_transient) {
        std::cout << "  Mode: transient\n";
        std::cout << "  Time step: " << dt_ns << " ns\n";
        std::cout << "  Final time: " << t_final_us << " µs\n";
        if (save_animation) {
            std::cout << "  Animation: enabled\n";
            std::cout << "  Snapshot directory: " << snapshot_dir << "\n";
            std::cout << "  Snapshot interval: " << snapshot_interval << " steps\n";
        }
    } else {
        std::cout << "  Mode: steady-state\n";
    }
    std::cout << std::endl;

    // Create solver
    pnp::PNPSolver1D solver(params);

    // Print Debye length
    double lambda_D = solver.get_debye_length();
    std::cout << "Computed Debye length: " << lambda_D * 1e9 << " nm\n";
    std::cout << "Domain / Debye length ratio: " << params.L / lambda_D << "\n\n";

    // Solve PNP equations
    bool converged = true;
    if (run_transient) {
        // Transient simulation
        double dt = dt_ns * 1e-9;           // Convert ns to s
        double t_final = t_final_us * 1e-6; // Convert µs to s

        // Create snapshot directory if needed
        std::filesystem::create_directories(snapshot_dir);

        if (use_efield) {
            // Use E-field formulation (inspired by PoNPs, stable)
            solver.solve_transient_efield(dt, t_final, snapshot_dir, snapshot_interval);
            std::cout << "\nSnapshots saved to: " << snapshot_dir << "\n";
            std::cout << "Use 'python3 scripts/create_animation.py' to generate GIF.\n";
        } else if (use_shenxu) {
            // Use Shen-Xu positivity preserving scheme (unconditionally stable)
            solver.solve_transient_shenxu(dt, t_final, snapshot_dir, snapshot_interval);
            std::cout << "\nSnapshots saved to: " << snapshot_dir << "\n";
            std::cout << "Use 'python3 scripts/create_animation.py' to generate GIF.\n";
        } else if (use_slotboom) {
            // Use Slotboom transformation (most stable for true transient)
            solver.solve_transient_slotboom(dt, t_final, snapshot_dir, snapshot_interval);
            std::cout << "\nSnapshots saved to: " << snapshot_dir << "\n";
            std::cout << "Use 'python3 scripts/create_animation.py' to generate GIF.\n";
        } else if (use_newton) {
            // Use fully implicit Newton method
            solver.solve_transient_newton(dt, t_final, snapshot_dir, snapshot_interval);
            std::cout << "\nSnapshots saved to: " << snapshot_dir << "\n";
            std::cout << "Use 'python3 scripts/create_animation.py' to generate GIF.\n";
        } else if (use_continuation) {
            // Use continuation method (most stable)
            solver.solve_transient_continuation(continuation_steps, snapshot_dir);
            std::cout << "\nSnapshots saved to: " << snapshot_dir << "\n";
            std::cout << "Use 'python3 scripts/create_animation.py' to generate GIF.\n";
        } else if (use_gummel) {
            // Use Gummel iteration (more stable)
            solver.solve_transient_gummel(dt, t_final, snapshot_dir, snapshot_interval);
            std::cout << "\nSnapshots saved to: " << snapshot_dir << "\n";
            std::cout << "Use 'python3 scripts/create_animation.py' to generate GIF.\n";
        } else if (save_animation) {
            solver.solve_transient_with_snapshots(dt, t_final, snapshot_dir, snapshot_interval);
            std::cout << "\nSnapshots saved to: " << snapshot_dir << "\n";
            std::cout << "Use 'python3 scripts/create_animation.py' to generate GIF.\n";
        } else {
            solver.solve_transient(dt, t_final);
        }
    } else {
        // Steady-state solution
        converged = solver.solve();

        if (converged) {
            std::cout << "\nSolution converged successfully!\n";
        } else {
            std::cout << "\nWarning: Solution may not have fully converged.\n";
        }
    }

    // Save results
    solver.save_results(output_file);

    // Print some statistics
    const auto& phi = solver.get_phi();
    const auto& c_plus = solver.get_c_plus();
    const auto& c_minus = solver.get_c_minus();

    std::cout << "\n========================================\n";
    std::cout << "Solution Summary:\n";
    std::cout << "========================================\n";

    // Get bulk potential for closed systems
    double phi_bulk = solver.get_bulk_potential();
    int n_mid = phi.size() / 2;

    std::cout << "  Left electrode (x=0):\n";
    std::cout << "    Potential: " << phi[0] * 1000.0 << " mV\n";
    std::cout << "    c+/c0: " << c_plus[0] / params.c0 << "\n";
    std::cout << "    c-/c0: " << c_minus[0] / params.c0 << "\n";
    std::cout << "  Bulk (x=L/2):\n";
    std::cout << "    Potential: " << phi[n_mid] * 1000.0 << " mV\n";
    std::cout << "    c+/c0: " << c_plus[n_mid] / params.c0 << "\n";
    std::cout << "    c-/c0: " << c_minus[n_mid] / params.c0 << "\n";
    std::cout << "  Right electrode (x=L):\n";
    std::cout << "    Potential: " << phi.back() * 1000.0 << " mV\n";
    std::cout << "    c+/c0: " << c_plus.back() / params.c0 << "\n";
    std::cout << "    c-/c0: " << c_minus.back() / params.c0 << "\n";

    if (params.closed_system) {
        std::cout << "  Bulk potential (self-consistent): " << phi_bulk * 1000.0 << " mV\n";
    }

    // Compute and display charge and capacitance
    auto [sigma_left, sigma_right] = solver.compute_surface_charge();
    auto [C_left, C_right] = solver.compute_capacitance();
    double C_total = solver.compute_total_capacitance();

    std::cout << "\nCharge Analysis:\n";
    std::cout << "  Surface charge density:\n";
    std::cout << "    Left electrode:  " << std::showpos << sigma_left * 1e2 << std::noshowpos << " μC/cm²\n";
    std::cout << "    Right electrode: " << std::showpos << sigma_right * 1e2 << std::noshowpos << " μC/cm²\n";
    std::cout << "  Electroneutrality check: Σσ = " << (sigma_left + sigma_right) * 1e2 << " μC/cm²\n";
    std::cout << "\nCapacitance:\n";
    std::cout << "  Single EDL:  " << C_left * 1e2 << " μF/cm² (left), " << C_right * 1e2 << " μF/cm² (right)\n";
    std::cout << "  Total (series): " << C_total * 1e2 << " μF/cm²\n";

    // Compute error metrics against Gouy-Chapman theory
    double L2_error = solver.compute_L2_error();
    double L2_rel_error = solver.compute_relative_L2_error();
    double Linf_error = solver.compute_Linf_error();

    std::cout << "\nError Analysis (vs. Gouy-Chapman theory):\n";
    std::cout << "  L2 error:          " << L2_error * 1000.0 << " mV\n";
    std::cout << "  Relative L2 error: " << L2_rel_error * 100.0 << " %\n";
    std::cout << "  L-inf (max) error: " << Linf_error * 1000.0 << " mV\n";

    std::cout << "\nResults saved to: " << output_file << "\n";
    std::cout << "Use 'python3 scripts/plot_results.py' to visualize.\n";

    return converged ? 0 : 1;
}
