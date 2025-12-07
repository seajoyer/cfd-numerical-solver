#ifndef ARTIFICIALVISCOSITY_HPP
#define ARTIFICIALVISCOSITY_HPP

#include <xtensor.hpp>

#include "config/Settings.hpp"

struct DataLayer;

/**
 * @class ArtificialViscosity
 * @brief Abstract interface for artificial viscosity models.
 *
 * Artificial viscosity is modeled here as an additional scalar quantity
 * q_{j+1/2} defined on cell interfaces. Typical usage:
 *
 *  - Compute q on all interfaces via ComputeInterfaceQ().
 *  - For each interface j+1/2, add q_{j+1/2} to the pressure of left/right
 *    reconstructed states before calling the Riemann solver:
 *
 *      P_L_eff = P_L + q_{j+1/2},
 *      P_R_eff = P_R + q_{j+1/2}.
 *
 * This increases numerical dissipation in compressive regions (e.g. shocks)
 * while leaving expansion regions unaffected.
 */
class ArtificialViscosity {
public:
    virtual ~ArtificialViscosity() = default;

    /**
     * @brief Computes artificial viscosity q on all cell interfaces.
     *
     * The implementation must fill the array `q` with size equal to the
     * number of interfaces: n_interfaces = total_size - 1, where total_size
     * is the linear size of 1D arrays in DataLayer (including ghosts).
     *
     * Interfaces are indexed as:
     *  - interface i corresponds to the face between cells i and i+1.
     *
     * @param layer    Current solution state (primitive fields, geometry).
     * @param dx       Spatial step size.
     * @param q        Output array of length total_size - 1 with q_{i}.
     *                 On entry it may have any shape; implementation should
     *                 resize/assign it as needed.
     */
    virtual void ComputeInterfaceQ(const DataLayer& layer,
                                   double dx,
                                   xt::xarray<double>& q) const = 0;
};

#endif  // ARTIFICIALVISCOSITY_HPP