#include "solver/SolverFactory.hpp"

#include <cctype>
#include <memory>
#include <stdexcept>
#include <string>

#include "parallel/HaloExchange.hpp"
#include "solver/FiniteVolumeSolver.hpp"
#include "spatial/GodunovKolganRodionovSpatialOperator.hpp"
#include "spatial/GodunovSpatialOperator.hpp"
#include "spatial/SpatialOperator.hpp"
#include "time/ForwardEulerTimeIntegrator.hpp"
#include "time/SSPRK2TimeIntegrator.hpp"
#include "time/SSPRK3TimeIntegrator.hpp"
#include "time/TimeIntegrator.hpp"

std::unique_ptr<Solver> SolverFactory::Create(const Settings& settings,
                                              const std::shared_ptr<Mesh>& mesh,
                                              const std::shared_ptr<BoundaryManager>& boundary_manager,
                                              const MPIContext* mpi_context,
                                              const StateSynchronizer* halo_exchange) {
    if (!mesh) {
        throw std::runtime_error("SolverFactory::Create: mesh is null");
    }

    if (!boundary_manager) {
        throw std::runtime_error("SolverFactory::Create: boundary_manager is null");
    }

    std::string solver_name = settings.solver;
    for (char& c : solver_name) {
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }

    std::string time_integrator_name = settings.time_integrator;
    for (char& c : time_integrator_name) {
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }

    std::shared_ptr<SpatialOperator> spatial_operator;
    if (solver_name == "godunov" || solver_name == "godunov-kolgan") {
        spatial_operator = std::make_shared<GodunovSpatialOperator>(settings, boundary_manager, halo_exchange);
    }
    else if (solver_name == "godunov-kolgan-rodionov") {
        spatial_operator = std::make_shared<GodunovKolganRodionovSpatialOperator>(
            settings, boundary_manager, halo_exchange);
    }
    else {
        throw std::runtime_error(
            "SolverFactory::Create: unsupported solver '" + settings.solver + "'"
        );
    }

    std::shared_ptr<TimeIntegrator> time_integrator;
    if (time_integrator_name == "euler") {
        time_integrator = std::make_shared<ForwardEulerTimeIntegrator>();
    }
    else if (time_integrator_name == "ssprk2") {
        time_integrator = std::make_shared<SSPRK2TimeIntegrator>();
    }
    else if (time_integrator_name == "ssprk3") {
        time_integrator = std::make_shared<SSPRK3TimeIntegrator>();
    }
    else {
        throw std::runtime_error(
            "SolverFactory::Create: unsupported time integrator '" +
            settings.time_integrator + "'"
        );
    }

    return std::make_unique<FiniteVolumeSolver>(
        settings,
        mesh,
        std::move(spatial_operator),
        std::move(time_integrator),
        mpi_context,
        halo_exchange
    );
}
