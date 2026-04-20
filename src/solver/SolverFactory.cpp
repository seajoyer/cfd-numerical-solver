#include "solver/SolverFactory.hpp"

#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

#include "solver/FiniteVolumeSolver.hpp"
#include "pressure_velocity/PimpleSolver.hpp"
#include "pressure_velocity/PisoSolver.hpp"
#include "pressure_velocity/SimpleSolver.hpp"
#include "spatial/FLICSpatialOperator.hpp"
#include "spatial/GodunovKolganRodionovSpatialOperator.hpp"
#include "spatial/GodunovSpatialOperator.hpp"
#include "spatial/MaderSpatialOperator.hpp"
#include "time/ForwardEulerTimeIntegrator.hpp"
#include "time/MacCormackTimeIntegrator.hpp"
#include "time/MaderTimeIntegrator.hpp"
#include "time/SSPRK2TimeIntegrator.hpp"
#include "time/SSPRK3TimeIntegrator.hpp"
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

        throw std::runtime_error("Unknown finite-volume solver type: " + settings.solver);
    }

    [[nodiscard]] bool IsPressureVelocitySolver(const std::string& solver_name) {
        const std::string s = utils::ToLower(solver_name);
        return s == "simple" || s == "piso" || s == "pimple";
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
                           const MPIContext* mpi_context,
                           std::shared_ptr<EOS> eos) -> std::unique_ptr<Solver> {
    const std::string solver = utils::ToLower(settings.solver);

    if (IsPressureVelocitySolver(solver)) {
        if (solver == "simple") {
            return std::make_unique<SimpleSolver>(settings,
                                                  std::move(mesh),
                                                  boundary_manager,
                                                  mpi_context);
        }

        if (solver == "piso") {
            return std::make_unique<PisoSolver>(settings,
                                                std::move(mesh),
                                                boundary_manager,
                                                mpi_context);
        }

        if (solver == "pimple") {
            return std::make_unique<PimpleSolver>(settings,
                                                  std::move(mesh),
                                                  boundary_manager,
                                                  mpi_context);
        }

        throw std::runtime_error("Unknown pressure-velocity solver type: " + settings.solver);
    }

    auto spatial_operator = CreateSpatialOperator(settings, boundary_manager);

    if (auto mader_op = std::dynamic_pointer_cast<MaderSpatialOperator>(spatial_operator)) {
        mader_op->SetEos(eos);
    }

    auto time_integrator = CreateTimeIntegrator(settings, boundary_manager);

    return std::make_unique<FiniteVolumeSolver>(settings,
                                                std::move(mesh),
                                                std::move(spatial_operator),
                                                std::move(time_integrator),
                                                mpi_context,
                                                eos);
}
