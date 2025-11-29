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

void print_usage(const char* program_name) {
    std::cout << "Usage: " << program_name << " [options]\n"
              << "\nOptions:\n"
              << "  --phi0 <value>      Surface potential in mV (default: 100)\n"
              << "  --c0 <value>        Bulk concentration in mol/L (default: 1.0)\n"
              << "  --eps <value>       Relative permittivity (default: 12)\n"
              << "  --L <value>         Domain length in nm (default: 100)\n"
              << "  --N <value>         Number of grid points (default: 1001)\n"
              << "  --output <file>     Output filename (default: results/pnp_results.dat)\n"
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

    std::cout << "Parameters:\n";
    std::cout << "  Surface potential: " << phi0_mV << " mV\n";
    std::cout << "  Bulk concentration: " << c0_molL << " mol/L\n";
    std::cout << "  Relative permittivity: " << eps_r << "\n";
    std::cout << "  Domain length: " << L_nm << " nm\n";
    std::cout << "  Grid points: " << N << "\n";
    std::cout << std::endl;

    // Create solver
    pnp::PNPSolver1D solver(params);

    // Print Debye length
    double lambda_D = solver.get_debye_length();
    std::cout << "Computed Debye length: " << lambda_D * 1e9 << " nm\n";
    std::cout << "Domain / Debye length ratio: " << params.L / lambda_D << "\n\n";

    // Solve steady-state PNP equations
    bool converged = solver.solve();

    if (converged) {
        std::cout << "\nSolution converged successfully!\n";
    } else {
        std::cout << "\nWarning: Solution may not have fully converged.\n";
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

    // Compare with Gouy-Chapman at surface
    auto phi_gc = solver.gouy_chapman_solution(params.phi_left);
    double max_error = 0.0;
    for (size_t i = 0; i < phi.size(); ++i) {
        max_error = std::max(max_error, std::abs(phi[i] - phi_gc[i]));
    }
    std::cout << "\nComparison with Gouy-Chapman theory:\n";
    std::cout << "  Maximum potential error: " << max_error * 1000.0 << " mV\n";

    std::cout << "\nResults saved to: " << output_file << "\n";
    std::cout << "Use 'python3 scripts/plot_results.py' to visualize.\n";

    return converged ? 0 : 1;
}
