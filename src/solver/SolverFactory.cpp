#include "solver/SolverFactory.hpp"

#include <algorithm>
#include <cctype>
#include <iostream>
#include <stdexcept>

#include "solver/AnalyticalSolver.hpp"
#include "solver/GodunovKolganRodionovSolver.hpp"
#include "solver/GodunovSolver.hpp"

namespace {
/**
 * @brief Validates that the solver and reconstruction combination is compatible
 * 
 * Compatibility rules:
 * - godunov: P0
 * - godunov-kolgan: P0, P1, ENO, WENO
 * - godunov-kolgan-rodionov: P1, ENO, WENO
 * - analytical: no reconstruction needed
 * 
 * @param solver Solver type (lowercase)
 * @param reconstruction Reconstruction scheme (lowercase)
 * @throws std::runtime_error if combination is incompatible
 */
void ValidateSolverReconstructionCompatibility(const std::string& solver,
                                               const std::string& reconstruction) {
    // Analytical solver doesn't use reconstruction
    if (solver == "analytical") {
        return;
    }
    
    // Godunov solver only supports P0
    if (solver == "godunov") {
        if (reconstruction.find("p0") == std::string::npos) {
            std::cerr << "\nError: Incompatible solver-reconstruction pair!\n";
            std::cerr << "  Solver: " << solver << "\n";
            std::cerr << "  Reconstruction: " << reconstruction << "\n\n";
            std::cerr << "The 'godunov' solver only supports P0 (piecewise constant) reconstruction.\n";
            std::cerr << "To see all compatible solver-reconstruction pairs, use:\n";
            std::cerr << "  ./cfd_numerical_solver --list-solvers\n\n";
            throw std::runtime_error("Incompatible solver-reconstruction pair");
        }
        return;
    }
    
    // Godunov-Kolgan-Rodionov requires at least P1
    if (solver == "godunov-kolgan-rodionov") {
        if (reconstruction.find("p0") != std::string::npos) {
            std::cerr << "\nError: Incompatible solver-reconstruction pair!\n";
            std::cerr << "  Solver: " << solver << "\n";
            std::cerr << "  Reconstruction: " << reconstruction << "\n\n";
            std::cerr << "The 'godunov-kolgan-rodionov' solver requires at least P1 reconstruction.\n";
            std::cerr << "It supports: P1, ENO (any order), WENO (any order)\n\n";
            std::cerr << "To see all compatible solver-reconstruction pairs, use:\n";
            std::cerr << "  ./cfd_numerical_solver --list-solvers\n\n";
            throw std::runtime_error("Incompatible solver-reconstruction pair");
        }
        return;
    }
    
    // Godunov-Kolgan supports all reconstruction schemes
    if (solver == "godunov-kolgan") {
        return;
    }
}
}  // namespace

auto SolverFactory::Create(Settings& settings) -> std::unique_ptr<Solver> {
    // Convert to lowercase for case-insensitive comparison
    std::string type_lower = settings.solver;
    std::transform(type_lower.begin(), type_lower.end(), type_lower.begin(),
                   [](unsigned char c) -> int { return std::tolower(c); });
    
    std::string recon_lower = settings.reconstruction;
    std::transform(recon_lower.begin(), recon_lower.end(), recon_lower.begin(),
                   [](unsigned char c) -> int { return std::tolower(c); });

    // Validate compatibility before creating solver
    ValidateSolverReconstructionCompatibility(type_lower, recon_lower);

    if (type_lower == "analytical") {
        return std::make_unique<AnalyticalSolver>(settings);
    }
    if (type_lower == "godunov") {
        return std::make_unique<GodunovSolver>(settings);
    }
    if (type_lower == "godunov-kolgan") {
        return std::make_unique<GodunovSolver>(settings);
    }
    if (type_lower == "godunov-kolgan-rodionov") {
        return std::make_unique<GodunovKolganRodionovSolver>(settings);
    }

    std::cerr << "\nError: Unknown solver type: " << settings.solver << "\n\n";
    std::cerr << "Supported solvers:\n";
    std::cerr << "  - godunov\n";
    std::cerr << "  - godunov-kolgan\n";
    std::cerr << "  - godunov-kolgan-rodionov\n";
    std::cerr << "  - analytical\n\n";
    std::cerr << "To see all solver options and their compatible reconstructions, use:\n";
    std::cerr << "  ./cfd_numerical_solver --list-solvers\n\n";
    
    throw std::runtime_error("Unknown solver type: " + settings.solver);
}
