#ifndef SSPRK3TIMEINTEGRATOR_HPP
#define SSPRK3TIMEINTEGRATOR_HPP

#include "TimeIntegrator.hpp"
#include "bc/BoundaryManager.hpp"

/**
 * @class SSPRK3TimeIntegrator
 * @brief Third-order strong-stability-preserving Runge–Kutta (SSPRK(3,3)).
 *
 * Implements the classical Shu–Osher SSPRK3 scheme:
 *
 *  U^{(1)} = U^n + dt * L(U^n)
 *  U^{(2)} = 3/4 U^n + 1/4 (U^{(1)} + dt * L(U^{(1)}))
 *  U^{n+1} = 1/3 U^n + 2/3 (U^{(2)} + dt * L(U^{(2)}))
 *
 * Positivity limiting may be applied after each stage or only at the end,
 * depending on the desired robustness.
 */
class SSPRK3TimeIntegrator : public TimeIntegrator {
public:
    SSPRK3TimeIntegrator() = default;

    SSPRK3TimeIntegrator(std::shared_ptr<BoundaryManager> bc);

    void Advance(DataLayer& layer,
                 double dt,
                 double dx,
                 const Settings& settings,
                 const SpatialOperator& op) const override;

private:
    std::shared_ptr<BoundaryManager> boundary_manager_;
};

#endif  // SSPRK3TIMEINTEGRATOR_HPP
