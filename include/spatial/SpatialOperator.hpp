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
 * in conservative form, evaluated on a given 1D mesh. Implementations
 * are responsible for:
 *  - reconstructing interface or cell-center states,
 *  - invoking a Riemann solver or physical flux,
 *  - assembling conservative flux differences to obtain the RHS.
 *
 * Time integration is not handled here; it is delegated to TimeIntegrator.
 */
class SpatialOperator {
public:
    virtual ~SpatialOperator() = default;

    /**
     * @brief Computes the semi-discrete RHS dU/dt for all cells.
     *
     * The operator evaluates the conservative right-hand side L(U) based on:
     *  - current state stored in DataLayer (including ghost cells),
     *  - mesh geometry stored in DataLayer (e.g. xc),
     *  - uniform cell size dx,
     *  - thermodynamic parameter gamma (from Settings).
     *
     * Implementations must fill rhs with the same length as the 1D grid.
     * Ghost-cell entries in rhs may be left unused or set to zero.
     *
     * @param layer Mesh and primitive state data.
     * @param dx    Uniform cell size.
     * @param gamma Ratio of specific heats.
     * @param rhs   Output array for dU/dt (same 1D length as the grid).
     */
    virtual void ComputeRHS(const DataLayer& layer,
                            double dx,
                            double gamma,
                            xt::xarray<Conservative>& rhs) const = 0;
};

#endif  // SPATIALOPERATOR_HPP
