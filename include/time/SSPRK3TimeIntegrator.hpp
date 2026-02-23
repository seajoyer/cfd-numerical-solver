#ifndef SSPRK3TIMEINTEGRATOR_HPP
#define SSPRK3TIMEINTEGRATOR_HPP

#include "time/TimeIntegrator.hpp"

/**
 * @class SSPRK3TimeIntegrator
 * @brief Third-order strong-stability-preserving Runge–Kutta (SSPRK(3,3)).
 *
 * Scheme (Shu–Osher):
 *  U1 = U^n + dt * L(U^n)
 *  U2 = 3/4 U^n + 1/4 (U1 + dt * L(U1))
 *  U^{n+1} = 1/3 U^n + 2/3 (U2 + dt * L(U2))
 *
 * Notes:
 *  - SpatialOperator handles halo + physical BC internally.
 *  - Updates only core cells.
 *  - PositivityLimiter is applied after the final stage.
 */
class SSPRK3TimeIntegrator final : public TimeIntegrator {
public:
    SSPRK3TimeIntegrator() = default;

    void Advance(DataLayer& layer,
                 Workspace& workspace,
                 double dt,
                 double gamma,
                 const SpatialOperator& op) const override;
};

#endif  // SSPRK3TIMEINTEGRATOR_HPP