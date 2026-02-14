#ifndef SPATIALOPERATOR_HPP
#define SPATIALOPERATOR_HPP

#include <xtensor.hpp>

#include "config/Settings.hpp"
#include "data/DataLayer.hpp"
#include "data/Variables.hpp"

/**
 * @class SpatialOperator
 * @brief Abstract semi-discrete finite-volume operator for the Euler equations.
 *
 * This interface defines the spatial semi-discrete operator
 *
 *    dU_j/dt = L_j(U),
 *
 * in conservative form. Implementations are responsible for:
 *  - reconstructing interface or cell-center states,
 *  - invoking a Riemann solver or physical flux,
 *  - assembling conservative flux differences to obtain the RHS.
 *
 * Time integration is not handled here; it is delegated to TimeIntegrator.
 *
 * For 2D operators, the 2-argument overload ComputeRHS(layer, dx, dy, gamma, rhs)
 * should be overridden. The 1-argument dx overload serves as a dispatcher:
 * for 2D data it defaults to dy = dx (square cells).
 */
class SpatialOperator {
public:
    virtual ~SpatialOperator() = default;

    /**
     * @brief Computes the semi-discrete RHS dU/dt for all cells (1D or square-cell 2D).
     *
     * For 1D operators, this is the primary entry point.
     * For 2D operators, this may dispatch to ComputeRHS2D with dy = dx.
     *
     * @param layer Mesh and primitive state data.
     * @param dx    Uniform cell size (x-direction).
     * @param gamma Ratio of specific heats.
     * @param rhs   Output array for dU/dt.
     */
    virtual void ComputeRHS(const DataLayer& layer,
                            double dx,
                            double gamma,
                            xt::xarray<Conservative>& rhs) const = 0;

    /**
     * @brief Computes the semi-discrete RHS for 2D with separate dx, dy.
     *
     * Default implementation calls ComputeRHS(layer, dx, gamma, rhs),
     * which is correct for 1D operators. 2D operators should override this.
     *
     * @param layer Mesh and primitive state data.
     * @param dx    Cell size in x-direction.
     * @param dy    Cell size in y-direction.
     * @param gamma Ratio of specific heats.
     * @param rhs   Output array for dU/dt.
     */
    virtual void ComputeRHS2D(const DataLayer& layer,
                              double dx, double dy,
                              double gamma,
                              xt::xarray<Conservative>& rhs) const {
        (void)dy;  // 1D operators ignore dy
        ComputeRHS(layer, dx, gamma, rhs);
    }
};

#endif  // SPATIALOPERATOR_HPP
