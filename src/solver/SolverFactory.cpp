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

void Validate2DReconstructionSupport(const std::string& reconstruction, int dim) {
    if (dim < 2) return;

    // 2D currently only supports P0 reconstruction in the spatial operator.
    // Higher-order reconstruction objects are created but unused in 2D.
    std::string recon_lower = reconstruction;
    std::transform(recon_lower.begin(), recon_lower.end(), recon_lower.begin(),
                   [](unsigned char c) -> char { return static_cast<char>(std::tolower(c)); });

    if (recon_lower.find("p0") == std::string::npos) {
        std::cerr << "\nWarning: 2D spatial operator currently uses P0 reconstruction only.\n"
                  << "  Configured reconstruction '" << reconstruction
                  << "' will be ignored in 2D.\n"
                  << "  The solver will run with P0 (piecewise-constant) reconstruction.\n\n";
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
        Validate2DReconstructionSupport(recon_lower, settings.dim);

        if (settings.time_integrator != "euler") {
            std::cerr << "Warning: 2D currently uses Forward Euler only. "
                      << "Configured time_integrator '" << settings.time_integrator
                      << "' will be ignored.\n";
        }

        if (type_lower == "godunov" || type_lower == "godunov-kolgan") {
            auto spatial_operator = std::make_shared<GodunovSpatialOperator2D>(settings);
            return std::make_unique<FiniteVolumeSolver>(settings, spatial_operator);
        }
        if (type_lower == "godunov-kolgan-rodionov") {
            // GKR in 2D: use GodunovSpatialOperator2D (P0 fallback)
            // Higher-order 2D reconstruction can be added later
            std::cerr << "\nWarning: godunov-kolgan-rodionov in 2D uses P0 spatial operator.\n"
                      << "  Second-order GKR predictor-corrector is not yet implemented for 2D.\n"
                      << "  Running with first-order Godunov method.\n\n";
            auto spatial_operator = std::make_shared<GodunovSpatialOperator2D>(settings);
            return std::make_unique<FiniteVolumeSolver>(settings, spatial_operator);
        }
        if (type_lower == "maccormack") {
            std::cerr << "\nError: MacCormack solver is not supported in 2D.\n\n";
            throw std::runtime_error("MacCormack solver does not support 2D");
        }

        std::cerr << "\nError: Unknown solver type for 2D: " << settings.solver << "\n\n";
        throw std::runtime_error("Unknown solver type: " + settings.solver);
    }

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
