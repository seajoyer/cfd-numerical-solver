#ifndef FORWARDEULERTIMEINTEGRATOR_HPP
#define FORWARDEULERTIMEINTEGRATOR_HPP

#include "time/TimeIntegrator.hpp"

/**
 * @class ForwardEulerTimeIntegrator
 * @brief First-order explicit Euler time integration for conservative FV systems.
 *
 * Implements:
 *   U^{n+1} = U^n + dt * RHS(U^n)
 */
class ForwardEulerTimeIntegrator final : public TimeIntegrator {
public:
    ForwardEulerTimeIntegrator() = default;

    void Advance(DataLayer& layer,
                 const Mesh& mesh,
                 Workspace& workspace,
                 double dt,
                 double gamma,
                 const SpatialOperator& op) const override;
};

#endif  // FORWARDEULERTIMEINTEGRATOR_HPP