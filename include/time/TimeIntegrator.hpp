#ifndef TIMEINTEGRATOR_HPP
#define TIMEINTEGRATOR_HPP

#include "data/DataLayer.hpp"
#include "data/Workspace.hpp"

class SpatialOperator;

/**
 * @class TimeIntegrator
 * @brief Abstract explicit time integration scheme for semi-discrete FV systems.
 *
 * Integrators advance the conservative state U stored in DataLayer.
 * They do NOT apply boundary conditions or halos; those are handled inside SpatialOperator.
 *
 * Workspace is provided by the caller and reused to avoid allocations.
 */
class TimeIntegrator {
public:
    virtual ~TimeIntegrator() = default;

    /**
     * @brief Advances the solution by one time step dt.
     *
     * @param layer Conservative state owner (updated in-place).
     * @param workspace Scratch buffers (W, rhs).
     * @param dt Time step size.
     * @param gamma Ratio of specific heats.
     * @param op Spatial operator providing RHS evaluations.
     */
    virtual void Advance(DataLayer& layer,
                         Workspace& workspace,
                         double dt,
                         double gamma,
                         const SpatialOperator& op) const = 0;

    /** @brief Sets positivity thresholds used by integrator post-update limiter. */
    virtual void SetPositivityThresholds(double rho_min, double p_min) {
        rho_min_ = rho_min;
        p_min_ = p_min;
    }

protected:
    double rho_min_ = 1e-10;
    double p_min_ = 1e-10;
};

#endif  // TIMEINTEGRATOR_HPP
