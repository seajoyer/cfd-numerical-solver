#ifndef MACCORMACKTIMEINTEGRATOR_HPP
#define MACCORMACKTIMEINTEGRATOR_HPP

#include <memory>

#include "time/TimeIntegrator.hpp"
#include "spatial/ForwardEulerSpatialOperator.hpp"
#include "spatial/BackwardEulerSpatialOperator.hpp"

class BoundaryManager;

/**
 * @class MacCormackTimeIntegrator
 * @brief Two-step predictor–corrector (MacCormack) for conservative FV systems.
 *
 * Predictor (forward difference):
 *   U* = U^n + dt * L_fwd(U^n)
 *
 * Corrector (backward difference):
 *   U^{n+1} = 0.5 * [ U^n + U* + dt * L_bwd(U*) ].
 *
 * Notes:
 *  - Forward/Backward spatial operators apply halo + physical BC internally.
 *  - Updates only core cells.
 *  - PositivityLimiter is applied after the final update.
 */
class MacCormackTimeIntegrator final : public TimeIntegrator {
public:
    MacCormackTimeIntegrator(const Settings& settings,
                             std::shared_ptr<BoundaryManager> boundary_manager);

    void Advance(DataLayer& layer,
                 Workspace& workspace,
                 double dt,
                 double gamma,
                 const SpatialOperator& op) const override;

private:
    // We ignore 'op' and use dedicated directional operators for MacCormack.
    // (Keeping the 'op' parameter for the common TimeIntegrator interface.)
    ForwardEulerSpatialOperator forward_op_;
    BackwardEulerSpatialOperator backward_op_;
};

#endif  // MACCORMACKTIMEINTEGRATOR_HPP