#ifndef SSPRK3TIMEINTEGRATOR_HPP
#define SSPRK3TIMEINTEGRATOR_HPP

#include "time/TimeIntegrator.hpp"

/**
 * @class SSPRK3TimeIntegrator
 * @brief Strong-Stability-Preserving 3-stage Runge-Kutta integrator.
 *
 * Scheme:
 *   Stage 1: U1       = U^n + dt * L(U^n)
 *   Stage 2: U2       = 3/4 U^n + 1/4 (U1 + dt * L(U1))
 *   Stage 3: U^{n+1}  = 1/3 U^n + 2/3 (U2 + dt * L(U2))
 *
 * Notes:
 * - Works on generic face-based meshes.
 * - SpatialOperator handles face fluxes and boundary conditions internally.
 * - PositivityLimiter is applied after the final stage.
 */
class SSPRK3TimeIntegrator final : public TimeIntegrator {
public:
    SSPRK3TimeIntegrator() = default;

    void Advance(DataLayer& layer,
                 const Mesh& mesh,
                 Workspace& workspace,
                 double dt,
                 double gamma,
                 const SpatialOperator& op,
                 const StateSynchronizer* halo_exchange) const override;
};

#endif  // SSPRK3TIMEINTEGRATOR_HPP
