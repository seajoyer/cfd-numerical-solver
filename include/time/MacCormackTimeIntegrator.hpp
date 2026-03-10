#ifndef MACCORMACKTIMEINTEGRATOR_HPP
#define MACCORMACKTIMEINTEGRATOR_HPP

#include <memory>

#include "spatial/BackwardEulerSpatialOperator.hpp"
#include "spatial/ForwardEulerSpatialOperator.hpp"
#include "time/TimeIntegrator.hpp"

class BoundaryManager;

/**
 * @class MacCormackTimeIntegrator
 * @brief Two-step predictor-corrector (MacCormack) for conservative FV systems.
 *
 * Predictor (forward difference):
 *   U* = U^n + dt * L_fwd(U^n)
 *
 * Corrector (backward difference):
 *   U^{n+1} = 0.5 * [ U^n + U* + dt * L_bwd(U*) ].
 *
 * Notes:
 *  - Forward/Backward spatial operators apply halo + physical BC internally.
 *  - Updates only fluid core cells.
 *  - PositivityLimiter is applied after the final update.
 */
class MacCormackTimeIntegrator final : public TimeIntegrator {
public:
    MacCormackTimeIntegrator(const Settings& settings,
                             std::shared_ptr<BoundaryManager> boundary_manager);

    void Advance(DataLayer& layer,
                 const Mesh& mesh,
                 Workspace& workspace,
                 double dt,
                 double gamma,
                 const SpatialOperator& op) const override;

private:
    ForwardEulerSpatialOperator forward_op_;
    BackwardEulerSpatialOperator backward_op_;
};

#endif  // MACCORMACKTIMEINTEGRATOR_HPP