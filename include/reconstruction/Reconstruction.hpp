#ifndef RECONSTRUCTION_HPP
#define RECONSTRUCTION_HPP

#include "data/DataLayer.hpp"
#include "data/Variables.hpp"

/**
 * @class Reconstruction
 * @brief Interface for 1D spatial reconstruction schemes.
 *
 * Reconstruction schemes take the cell-centered primitive variables stored
 * in a DataLayer and build left/right primitive states at cell interfaces.
 *
 * Implementations in this project include:
 *  - P0Reconstruction: piecewise constant Godunov reconstruction,
 *  - P1Reconstruction: piecewise linear reconstruction with a minmod limiter.
 *
 * Design:
 *  - The interface is per-interface and stateless: given an interface index,
 *    it computes the corresponding left/right states using the current
 *    DataLayer.
 *  - This design avoids storing temporary vectors of reconstructed states
 *    and works naturally with SoA storage (DataLayer).
 */
class Reconstruction {
   public:
    virtual ~Reconstruction() = default;

    /**
     * @brief Computes left/right primitive states at a given interface.
     *
     * Interface i is understood as lying between cell i and cell i+1
     * in the 1D layout (including ghost cells). Implementations may
     * use neighboring cells (ghosts included) to construct higher-order
     * approximations.
     *
     * @param layer          DataLayer containing current primitive fields.
     * @param interfaceIndex Index of the interface (0 ≤ i < totalSize - 1).
     * @param leftState      Output: left primitive state at the interface.
     * @param rightState     Output: right primitive state at the interface.
     */
    virtual void ComputeInterfaceStates(const DataLayer& layer, int interface_index,
                                        Primitive& left_state,
                                        Primitive& right_state) const = 0;
};

#endif  // RECONSTRUCTION_HPP
