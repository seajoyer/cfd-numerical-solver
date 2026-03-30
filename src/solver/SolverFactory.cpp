#include "solver/SolverFactory.hpp"

#include <algorithm>
#include <cctype>
#include <iostream>
#include <stdexcept>
#include <string>

#include "solver/FiniteVolumeSolver.hpp"
#include "spatial/FLICSpatialOperator.hpp"
#include "spatial/MaderSpatialOperator.hpp"
#include "spatial/GodunovKolganRodionovSpatialOperator.hpp"
#include "spatial/GodunovSpatialOperator.hpp"
#include "time/ForwardEulerTimeIntegrator.hpp"
#include "time/MacCormackTimeIntegrator.hpp"
#include "time/MaderTimeIntegrator.hpp"
#include "time/SSPRK2TimeIntegrator.hpp"
#include "time/SSPRK3TimeIntegrator.hpp"

#include "solver/SolverFactory.hpp"

#include <memory>
#include <stdexcept>
#include <utility>

#include "utils/StringUtils.hpp"

namespace {
    auto CreateTimeIntegrator(const Settings& settings,
                              const std::shared_ptr<BoundaryManager>& boundary_manager)
        -> std::shared_ptr<TimeIntegrator> {
        const std::string ti = utils::ToLower(settings.time_integrator);

        if (ti == "euler") {
            return std::make_shared<ForwardEulerTimeIntegrator>();
        }
        if (ti == "ssprk2") {
            return std::make_shared<SSPRK2TimeIntegrator>();
        }
        if (ti == "ssprk3") {
            return std::make_shared<SSPRK3TimeIntegrator>();
        }
        if (ti == "maccormack") {
            return std::make_shared<MacCormackTimeIntegrator>(settings, boundary_manager);
        }
        if (ti == "mader") {
            return std::make_shared<MaderTimeIntegrator>();
        }

        throw std::runtime_error("Unknown time integrator type: " + settings.time_integrator);
    }

    auto CreateSpatialOperator(const Settings& settings,
                               const std::shared_ptr<BoundaryManager>& boundary_manager)
        -> std::shared_ptr<SpatialOperator> {
        const std::string solver = utils::ToLower(settings.solver);

        if (solver == "godunov" || solver == "godunov-kolgan") {
            return std::make_shared<GodunovSpatialOperator>(settings, boundary_manager);
        }

        if (solver == "godunov-kolgan-rodionov") {
            return std::make_shared<GodunovKolganRodionovSpatialOperator>(settings, boundary_manager);
        }

        if (solver == "flic") {
            return std::make_shared<FLICSpatialOperator>(settings, boundary_manager);
        }

        if (solver == "mader") {
            return std::make_shared<MaderSpatialOperator>(settings, boundary_manager);
        }

        throw std::runtime_error("Unknown solver type: " + settings.solver);
    }
} // namespace

void SolverFactory::AddBoundary(const Axis axis,
                                std::shared_ptr<BoundaryCondition> left_bc,
                                std::shared_ptr<BoundaryCondition> right_bc) {
    boundary_manager_->Set(axis, std::move(left_bc), std::move(right_bc));
}

auto SolverFactory::Create(const Settings& settings,
                           Mesh mesh,
                           const std::shared_ptr<BoundaryManager>& boundary_manager,
                           const MPIContext* mpi_context) -> std::unique_ptr<Solver> {
    auto spatial_operator = CreateSpatialOperator(settings, boundary_manager);
    auto time_integrator = CreateTimeIntegrator(settings, boundary_manager);

    return std::make_unique<FiniteVolumeSolver>(settings,
                                                std::move(mesh),
                                                std::move(spatial_operator),
                                                std::move(time_integrator),
                                                mpi_context
    );
}
