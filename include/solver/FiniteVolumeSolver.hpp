// FiniteVolumeSolver.hpp
#ifndef FINITEVOLUMESOLVER_HPP
#define FINITEVOLUMESOLVER_HPP

#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

#include "solver/Solver.hpp"
#include "config/Settings.hpp"
#include "bc/BoundaryManager.hpp"
#include "spatial/SpatialOperator.hpp"
#include "time/TimeIntegrator.hpp"
#include "solver/TimeStepCalculator.hpp"
#include "data/DataLayer.hpp"
#include "data/Workspace.hpp"

class FiniteVolumeSolver final : public Solver {
public:
    FiniteVolumeSolver(const Settings& settings,
                       std::shared_ptr<SpatialOperator> spatial_operator,
                       std::shared_ptr<TimeIntegrator> time_integrator);

    auto Step(DataLayer& layer, double& t_cur) -> double override;
    void SetCfl(double cfl) override;

private:
    Settings settings_;

    std::shared_ptr<SpatialOperator> spatial_operator_;
    std::shared_ptr<TimeIntegrator> time_integrator_;

    Workspace workspace_;

    double rho_min_ = 1e-10;
    double p_min_ = 1e-10;

    void EnsureWorkspaceSized(const DataLayer& layer);
};

#endif  // FINITEVOLUMESOLVER_HPP
