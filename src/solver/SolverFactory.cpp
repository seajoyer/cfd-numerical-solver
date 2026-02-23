#include "solver/SolverFactory.hpp"

#include <algorithm>
#include <cctype>
#include <iostream>
#include <stdexcept>
#include <string>

#include "solver/FiniteVolumeSolver.hpp"

#include "spatial/GodunovSpatialOperator.hpp"
#include "spatial/GodunovKolganRodionovSpatialOperator.hpp"
// #include "spatial/MacCormackSpatialOperator.hpp"

#include "time/ForwardEulerTimeIntegrator.hpp"
#include "time/SSPRK2TimeIntegrator.hpp"
#include "time/SSPRK3TimeIntegrator.hpp"
#include "time/MacCormackTimeIntegrator.hpp"


static std::string ToLowerCopy(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(),
                   [](unsigned char c) -> char {
                       return static_cast<char>(std::tolower(c));
                   });
    return s;
}

static void ValidateSolverReconstructionCompatibility(const std::string& solver_lower,
                                                      const std::string& reconstruction_lower) {
    if (solver_lower == "analytical") return;

    if (solver_lower == "godunov") {
        if (reconstruction_lower.find("p0") == std::string::npos) {
            std::cerr << "\nError: Incompatible solver-reconstruction pair!\n"
                << "  Solver: " << solver_lower
                << ", Reconstruction: " << reconstruction_lower << "\n";
            throw std::runtime_error("Incompatible solver-reconstruction pair");
        }
        return;
    }

    if (solver_lower == "godunov-kolgan-rodionov") {
        if (reconstruction_lower.find("p0") != std::string::npos) {
            std::cerr << "\nError: godunov-kolgan-rodionov requires at least P1\n";
            throw std::runtime_error("Incompatible solver-reconstruction pair");
        }
    }
}

static auto CreateTimeIntegrator(const Settings& settings,
                                 const std::shared_ptr<BoundaryManager>& boundary_manager) -> std::shared_ptr<
    TimeIntegrator> {
    const std::string ti = ToLowerCopy(settings.time_integrator);

    if (ti == "euler") return std::make_shared<ForwardEulerTimeIntegrator>();
    if (ti == "ssprk2") return std::make_shared<SSPRK2TimeIntegrator>();
    if (ti == "ssprk3") return std::make_shared<SSPRK3TimeIntegrator>();
    if (ti == "maccormack") return std::make_shared<MacCormackTimeIntegrator>(settings, boundary_manager);

    throw std::runtime_error("Unknown time integrator type: " + settings.time_integrator);
}

static auto CreateSpatialOperator(const Settings& settings,
                                  const std::shared_ptr<BoundaryManager>& boundary_manager)
    -> std::shared_ptr<SpatialOperator> {
    const std::string solver_lower = ToLowerCopy(settings.solver);

    if (solver_lower == "godunov" || solver_lower == "godunov-kolgan") {
        return std::make_shared<GodunovSpatialOperator>(settings, boundary_manager);
    }

    if (solver_lower == "godunov-kolgan-rodionov") {
        return std::make_shared<GodunovKolganRodionovSpatialOperator>(settings, boundary_manager);
    }

    throw std::runtime_error("Unknown solver type: " + settings.solver);
}

void SolverFactory::AddBoundary(const Axis axis,
                                std::shared_ptr<BoundaryCondition> left_bc,
                                std::shared_ptr<BoundaryCondition> right_bc) {
    boundary_manager_->Set(axis, std::move(left_bc), std::move(right_bc));
}

auto SolverFactory::Create(Settings& settings,
                           const std::shared_ptr<BoundaryManager>& boundary_manager) -> std::unique_ptr<Solver> {
    const std::string solver_lower = ToLowerCopy(settings.solver);
    const std::string recon_lower = ToLowerCopy(settings.reconstruction);

    ValidateSolverReconstructionCompatibility(solver_lower, recon_lower);

    auto spatial_operator = CreateSpatialOperator(settings, boundary_manager);
    auto time_integrator = CreateTimeIntegrator(settings, boundary_manager);

    return std::make_unique<FiniteVolumeSolver>(settings, spatial_operator, time_integrator);
}
