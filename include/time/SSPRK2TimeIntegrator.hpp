#ifndef SSPRK2TIMEINTEGRATOR_HPP
#define SSPRK2TIMEINTEGRATOR_HPP

#include "time/TimeIntegrator.hpp"
#include "data/DataLayer.hpp"
#include "data/Variables.hpp"
#include "config/Settings.hpp"
#include "bc/BoundaryManager.hpp"

/**
 * @class SSPRK2TimeIntegrator
 * @brief Strong-Stability-Preserving 2-stage Runge–Kutta (TVD RK2) integrator.
 *
 * This class implements the classical second-order SSPRK2 scheme:
 *
 *  Stage 1:
 *    U^(1) = U^n + dt * L(U^n)
 *
 *  Stage 2:
 *    U^{n+1} = 0.5 * U^n + 0.5 * ( U^(1) + dt * L(U^(1)) )
 *
 * where L(U) is the semi-discrete spatial operator provided by SpatialOperator.
 *
 * Design:
 *  - Conservative variables are built from DataLayer via EOS::PrimToCons().
 *  - The SpatialOperator computes rhs = dU/dt.
 *  - PositivityLimiter is applied only at the final stage (U^{n+1}).
 *  - Core cells are written back to DataLayer via StoreConservativeCell().
 */
class SSPRK2TimeIntegrator : public TimeIntegrator {
public:
    SSPRK2TimeIntegrator() = default;

    /// Constructor with small positivity thresholds.
    SSPRK2TimeIntegrator(std::shared_ptr<BoundaryManager> bm);

    /**
     * @brief Advances the solution by one SSPRK2 step.
     *
     * @param layer DataLayer with current primitive state (modified in-place).
     * @param dt    Time step size.
     * @param dx    Spatial step size.
     * @param settings Global simulation settings (gamma, etc.).
     * @param op    Spatial operator providing rhs = dU/dt.
     */
    void Advance(DataLayer& layer,
                 double dt,
                 double dx,
                 const Settings& settings,
                 const SpatialOperator& op) const override;

private:
    std::shared_ptr<BoundaryManager> boundary_manager_;
};

#endif  // SSPRK2TIMEINTEGRATOR_HPP