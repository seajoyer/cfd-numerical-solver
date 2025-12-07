#ifndef FORWARDEULERTIMEINTEGRATOR_HPP
#define FORWARDEULERTIMEINTEGRATOR_HPP

#include "TimeIntegrator.hpp"

/**
 * @class ForwardEulerTimeIntegrator
 * @brief First-order explicit Euler time integration for FV systems.
 *
 * Implements:
 *
 *   U^{n+1} = U^n + dt * L(U^n),
 *
 * where L(U) is provided by a SpatialOperator. Positivity limiting
 * may be applied after the update.
 */
class ForwardEulerTimeIntegrator : public TimeIntegrator {
public:
    ForwardEulerTimeIntegrator();

    void Advance(DataLayer& layer,
                 double dt,
                 double dx,
                 const Settings& settings,
                 const SpatialOperator& op) const override;
};

#endif  // FORWARDEULERTIMEINTEGRATOR_HPP
