#include "solver/SolverFactory.hpp"

#include <algorithm>
#include <cctype>
#include <iostream>
#include <stdexcept>

#include "solver/AnalyticalSolver.hpp"
#include "solver/FiniteVolumeSolver.hpp"
#include "solver/GodunovKolganRodionovSolver.hpp"
#include "solver/GodunovSolver.hpp"
#include "solver/MacCormackSolver.hpp"
#include "spatial/ForwardEulerSpatialOperator.hpp"
#include "spatial/GodunovKolganRodionovSpatialOperator.hpp"
#include "spatial/GodunovSpatialOperator.hpp"
#include "spatial/GodunovSpatialOperator2D.hpp"
#include "time/ForwardEulerTimeIntegrator.hpp"
#include "time/SSPRK3TimeIntegrator.hpp"

namespace {
void ValidateSolverReconstructionCompatibility(const std::string& solver,
                                               const std::string& reconstruction) {
    if (solver == "analytical") return;

    if (solver == "godunov") {
        if (reconstruction.find("p0") == std::string::npos) {
            std::cerr << "\nError: Incompatible solver-reconstruction pair!\n";
            std::cerr << "  Solver: " << solver << ", Reconstruction: " << reconstruction << "\n";
            throw std::runtime_error("Incompatible solver-reconstruction pair");
        }
        return;
    }

    if (solver == "godunov-kolgan-rodionov") {
        if (reconstruction.find("p0") != std::string::npos) {
            std::cerr << "\nError: godunov-kolgan-rodionov requires at least P1\n";
            throw std::runtime_error("Incompatible solver-reconstruction pair");
        }
        return;
    }
}
} // namespace

auto SolverFactory::Create(Settings& settings) -> std::unique_ptr<Solver> {
    std::string type_lower = settings.solver;
    std::string recon_lower = settings.reconstruction;

    ValidateSolverReconstructionCompatibility(type_lower, recon_lower);

    if (type_lower == "analytical") {
        return std::make_unique<AnalyticalSolver>(settings);
    }

    // For 2D, use GodunovSpatialOperator2D
    if (settings.dim >= 2) {
        if (type_lower == "godunov" || type_lower == "godunov-kolgan") {
            auto spatial_operator = std::make_shared<GodunovSpatialOperator2D>(settings);
            return std::make_unique<FiniteVolumeSolver>(settings, spatial_operator);
        }
        if (type_lower == "godunov-kolgan-rodionov") {
            // For now, use GodunovSpatialOperator2D for 2D (P0 reconstruction)
            // Higher-order 2D reconstruction can be added later
            auto spatial_operator = std::make_shared<GodunovSpatialOperator2D>(settings);
            return std::make_unique<FiniteVolumeSolver>(settings, spatial_operator);
        }
    }

    // 1D path (unchanged)
    if (type_lower == "godunov" || type_lower == "godunov-kolgan") {
        auto spatial_operator = std::make_shared<GodunovSpatialOperator>(settings);
        return std::make_unique<FiniteVolumeSolver>(settings, spatial_operator);
    }
    if (type_lower == "godunov-kolgan-rodionov") {
        auto spatial_operator =
            std::make_shared<GodunovKolganRodionovSpatialOperator>(settings);
        return std::make_unique<FiniteVolumeSolver>(settings, spatial_operator);
    }
    if (type_lower == "maccormack") {
        return std::make_unique<FiniteVolumeSolver>(settings, nullptr);
    }

    std::cerr << "\nError: Unknown solver type: " << settings.solver << "\n\n";
    throw std::runtime_error("Unknown solver type: " + settings.solver);
}
