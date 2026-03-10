#ifndef SSPRK2TIMEINTEGRATOR_HPP
#define SSPRK2TIMEINTEGRATOR_HPP

#include "time/TimeIntegrator.hpp"

/**
 * @class SSPRK2TimeIntegrator
 * @brief Strong-Stability-Preserving 2-stage Runge–Kutta (TVD RK2).
 *
 * Scheme:
 *  Stage 1: U1 = U^n + dt * L(U^n)
 *  Stage 2: U^{n+1} = 0.5 * U^n + 0.5 * (U1 + dt * L(U1))
 *
 * Notes:
 *  - SpatialOperator handles halo + physical BC internally.
 *  - This integrator updates only fluid core cells.
 *  - PositivityLimiter is applied after the final stage.
 */
class SSPRK2TimeIntegrator final : public TimeIntegrator {
public:
    SSPRK2TimeIntegrator() = default;

    void Advance(DataLayer& layer,
                 const Mesh& mesh,
                 Workspace& workspace,
                 double dt,
                 double gamma,
                 const SpatialOperator& op) const override;
};

#endif  // SSPRK2TIMEINTEGRATOR_HPP