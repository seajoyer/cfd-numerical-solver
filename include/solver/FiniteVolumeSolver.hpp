#ifndef FINITEVOLUMESOLVER_HPP
#define FINITEVOLUMESOLVER_HPP

#include <memory>

#include "Solver.hpp"
#include "bc/BoundaryManager.hpp"
#include "config/Settings.hpp"
#include "data/DataLayer.hpp"
#include "data/Variables.hpp"
#include "spatial/SpatialOperator.hpp"
#include "time/TimeIntegrator.hpp"
#include "limiter/GlobalLimiter.hpp"
#include "limiter/VacuumFixLimiter.hpp"
#include "filter/SolutionFilter.hpp"

/**
 * @class FiniteVolumeSolver
 * @brief Generic finite-volume solver supporting 1D and 2D via pluggable components.
 *
 * For 1D: delegates to TimeIntegrator + SpatialOperator as before.
 * For 2D: uses Step2D() with Forward Euler and the 2D spatial operator.
 */
class FiniteVolumeSolver : public Solver {
public:
    FiniteVolumeSolver(const Settings& settings,
                       std::shared_ptr<SpatialOperator> spatial_operator);

    auto Step(DataLayer& layer, double& t_cur) -> double override;
    void SetCfl(double cfl) override;
    void AddBoundary(int axis,
                     std::shared_ptr<BoundaryCondition> left_bc,
                     std::shared_ptr<BoundaryCondition> right_bc) override;

private:
    Settings settings_;
    std::shared_ptr<BoundaryManager> boundary_manager_;
    std::shared_ptr<SpatialOperator> spatial_operator_;
    std::shared_ptr<TimeIntegrator> time_integrator_;
    double rho_min_;
    double p_min_;

    std::unique_ptr<GlobalLimiter> global_limiter_;
    std::unique_ptr<SolutionFilter> diffusion_;
    std::unique_ptr<VacuumFixLimiter> vacuum_fix_limiter_;
    SolutionFilter solution_filter_;

    [[nodiscard]] auto ComputeDx(const DataLayer& layer) const -> double;
    [[nodiscard]] auto ComputeDy(const DataLayer& layer) const -> double;

    /**
     * @brief 2D time step using Forward Euler with the 2D spatial operator.
     */
    void Step2D(DataLayer& layer, double dt, double dx, double dy);

    void InitializeTimeIntegrator();
};

#endif  // FINITEVOLUMESOLVER_HPPP
