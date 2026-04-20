#ifndef SSPRK2TIMEINTEGRATOR_HPP
#define SSPRK2TIMEINTEGRATOR_HPP

#include "time/TimeIntegrator.hpp"

/**
 * @class SSPRK2TimeIntegrator
 * @brief Strong-Stability-Preserving 2-stage Runge-Kutta integrator.
 *
 * Scheme:
 *   Stage 1: U1       = U^n + dt * L(U^n)
 *   Stage 2: U^{n+1}  = 0.5 * U^n + 0.5 * (U1 + dt * L(U1))
 *
 * Notes:
 * - Works on generic face-based meshes.
 * - SpatialOperator handles face fluxes and boundary conditions internally.
 * - PositivityLimiter is applied after the final stage.
 */
class SSPRK2TimeIntegrator final : public TimeIntegrator {
public:
    SSPRK2TimeIntegrator() = default;

    void Advance(DataLayer& layer,
                 const Mesh& mesh,
                 Workspace& workspace,
                 double dt,
                 double gamma,
                 const SpatialOperator& op,
                 const StateSynchronizer* halo_exchange) const override;
};

#endif  // SSPRK2TIMEINTEGRATOR_HPP
