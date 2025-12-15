#ifndef MACCORMACKTIMEINTEGRATOR_HPP
#define MACCORMACKTIMEINTEGRATOR_HPP

#include "spatial/ForwardEulerSpatialOperator.hpp"
#include "spatial/BackwardEulerSpatialOperator.hpp"
#include "TimeIntegrator.hpp"
#include "bc/BoundaryManager.hpp"

/**
 * @class MacCormackTimeIntegrator
 * @brief Two-step predictor–corrector (MacCormack) time integration.
 *
 * Classical MacCormack scheme in conservative form:
 *
 *  Predictor (forward):
 *    U^* = U^n + dt * L_forward(U^n)
 *
 *  Corrector (backward):
 *    U^{n+1} = 0.5 * [ U^n + U^* + dt * L_backward(U^*) ].
 *
 * Here L_forward and L_backward can be constructed from the same
 * SpatialOperator with directional bias, or simply approximated
 * through multiple calls to a symmetric L(U) for a more generic
 * implementation. Positivity limiting can be applied at the end.
 */
class MacCormackTimeIntegrator : public TimeIntegrator {
public:
    MacCormackTimeIntegrator() = default;
    MacCormackTimeIntegrator(std::shared_ptr<BoundaryManager> bm, Settings &settings);

    void Advance(DataLayer& layer,
                 double dt,
                 double dx,
                 const Settings& settings,
                 const SpatialOperator& op) const override;

private:
    std::unique_ptr<ForwardEulerSpatialOperator> forward_operator_;
    std::unique_ptr<BackwardEulerSpatialOperator> backward_operator_;

    std::shared_ptr<BoundaryManager> boundary_manager_;
};

#endif  // MACCORMACKTIMEINTEGRATOR_HPP