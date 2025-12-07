#ifndef TIMEINTEGRATOR_HPP
#define TIMEINTEGRATOR_HPP

#include <xtensor.hpp>

#include "config/Settings.hpp"
#include "data/DataLayer.hpp"
#include "data/Variables.hpp"
#include "spatial/SpatialOperator.hpp"

/**
 * @class TimeIntegrator
 * @brief Abstract time integration scheme for semi-discrete FV systems.
 *
 * Given a semi-discrete system
 *
 *   dU/dt = L(U),
 *
 * this interface advances the conservative solution by a single time step
 * of size dt, using a provided SpatialOperator for repeated evaluations
 * of L(U).
 *
 * Implementations may internally:
 *  - extract U from DataLayer into xt::xarray<Conservative>,
 *  - perform one or more stages (Forward Euler, SSPRK2/3, MacCormack, ...),
 *  - apply local positivity limiters during or after stages,
 *  - write the final U^{n+1} back into DataLayer.
 */
class TimeIntegrator {
public:
    virtual ~TimeIntegrator() = default;

    /**
     * @brief Advances the solution by a single time step dt.
     *
     * The method operates in-place on the DataLayer:
     *  - reads the current state from DataLayer,
     *  - performs time integration using the spatial operator,
     *  - writes the updated state back into DataLayer.
     *
     * @param layer    DataLayer containing the current state (modified in-place).
     * @param dt       Time step to perform (assumed already CFL-limited).
     * @param dx       Spatial step size.
     * @param settings Global simulation settings (for gamma, etc.).
     */
    virtual void Advance(DataLayer& layer,
                         double dt,
                         double dx,
                         const Settings& settings,
                         const SpatialOperator& op) const = 0;

    /**
     * @brief Sets positivity thresholds used inside the integrator.
     *
     * Implementations may use these thresholds for PositivityLimiter
     * calls at intermediate or final stages.
     *
     * @param rho_min Minimal allowed density.
     * @param p_min   Minimal allowed pressure.
     */
    virtual void SetPositivityThresholds(double rho_min, double p_min);

    /**
     * @brief Writes a conservative state back into the DataLayer.
     *
     * Converts (rho, rhoU, E) into primitive and auxiliary fields:
     *  - rho, u, P
     *  - p, V, U, e, m
     *
     * using EOS::Pressure and the same semantics as GodunovSolver.
     */
    void StoreConservativeCell(const Conservative& uc,
                               int i,
                               double dx,
                               const Settings& settings,
                               DataLayer& layer) const;
protected:
    double rho_min_{1e-6};
    double p_min_{1e-6};
};

#endif  // TIMEINTEGRATOR_HPP