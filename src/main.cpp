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
              << "  --phi0 <value>      Surface potential in mV (default: 100)\n"
              << "  --c0 <value>        Bulk concentration in mol/L (default: 1.0)\n"
              << "  --eps <value>       Relative permittivity (default: 12)\n"
              << "  --L <value>         Domain length in nm (default: 100)\n"
              << "  --N <value>         Number of grid points (default: 1001)\n"
              << "  --output <file>     Output filename (default: results/pnp_results.dat)\n"
              << "  --model <type>      Model type: standard or bikerman (default: standard)\n"
              << "  --ion-size <value>  Ion diameter in nm for Bikerman model (default: 0.7)\n"
              << "  --stretch <value>   Grid stretching factor (default: 3.0, use 1.0 for uniform)\n"
              << "  --transient         Run transient simulation instead of steady-state\n"
              << "  --dt <value>        Time step in ns for transient (default: 0.1)\n"
              << "  --t-final <value>   Final time in µs for transient (default: 1.0)\n"
              << "  --animation         Save snapshots for GIF animation (transient mode)\n"
              << "  --snapshot-dir <dir> Directory for animation snapshots (default: results/snapshots)\n"
              << "  --snapshot-interval <n> Save snapshot every n steps (default: 10)\n"
              << "  --help              Show this help message\n"
              << std::endl;
}

int main(int argc, char* argv[]) {
    // Default parameters
    double phi0_mV = 100.0;       // Surface potential [mV]
    double c0_molL = 1.0;         // Bulk concentration [mol/L]
    double eps_r = 12.0;          // Relative permittivity
    double L_nm = 100.0;          // Domain length [nm]
    int N = 1001;                 // Grid points
    std::string output_file = "results/pnp_results.dat";
    std::string model_type = "standard";  // Model type: standard or bikerman
    double ion_size_nm = 0.7;     // Ion diameter for Bikerman model [nm]
    double grid_stretch = 3.0;    // Grid stretching factor
    bool run_transient = false;   // Run transient simulation
    double dt_ns = 0.1;           // Time step [ns]
    double t_final_us = 1.0;      // Final time [µs]
    bool save_animation = false;  // Save snapshots for animation
    std::string snapshot_dir = "results/snapshots";  // Snapshot output directory
    int snapshot_interval = 10;   // Save every n steps

    // Parse command line arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--help") {
            print_usage(argv[0]);
            return 0;
        } else if (arg == "--phi0" && i + 1 < argc) {
            phi0_mV = std::atof(argv[++i]);
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
        } else if (arg == "--stretch" && i + 1 < argc) {
            grid_stretch = std::atof(argv[++i]);
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
        }
    }

    std::cout << "========================================\n";
    std::cout << "  1D Poisson-Nernst-Planck Solver\n";
    std::cout << "  for Ionic Liquid Electric Double Layer\n";
    std::cout << "========================================\n\n";

    // Setup parameters
    pnp::PNPParameters params;
    params.phi_left = phi0_mV * 1e-3;      // Convert mV to V
    params.c0 = c0_molL * 1000.0;          // Convert mol/L to mol/m^3
    params.eps_r = eps_r;
    params.L = L_nm * 1e-9;                // Convert nm to m
    params.N = N;
    params.a = ion_size_nm * 1e-9;         // Convert nm to m
    params.grid_stretch = grid_stretch;

    // Set model type
    if (model_type == "bikerman") {
        params.model = pnp::ModelType::BIKERMAN;
    } else {
        params.model = pnp::ModelType::STANDARD_PB;
    }

    std::cout << "Parameters:\n";
    std::cout << "  Surface potential: " << phi0_mV << " mV\n";
    std::cout << "  Bulk concentration: " << c0_molL << " mol/L\n";
    std::cout << "  Relative permittivity: " << eps_r << "\n";
    std::cout << "  Domain length: " << L_nm << " nm\n";
    std::cout << "  Grid points: " << N << "\n";
    std::cout << "  Model: " << model_type << "\n";
    if (model_type == "bikerman") {
        std::cout << "  Ion size: " << ion_size_nm << " nm\n";
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

        if (save_animation) {
            // Create snapshot directory
            std::filesystem::create_directories(snapshot_dir);
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
    std::cout << "  Potential at surface: " << phi[0] * 1000.0 << " mV\n";
    std::cout << "  Potential at bulk: " << phi.back() * 1000.0 << " mV\n";
    std::cout << "  c+ at surface: " << c_plus[0] << " mol/m^3 (ratio: "
              << c_plus[0] / params.c0 << ")\n";
    std::cout << "  c- at surface: " << c_minus[0] << " mol/m^3 (ratio: "
              << c_minus[0] / params.c0 << ")\n";
    std::cout << "  c+ at bulk: " << c_plus.back() << " mol/m^3\n";
    std::cout << "  c- at bulk: " << c_minus.back() << " mol/m^3\n";

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
