#ifndef RECONSTRUCTION_HPP
#define RECONSTRUCTION_HPP

#include "data/DataLayer.hpp"
#include "data/Variables.hpp"

/**
 * @class Reconstruction
 * @brief Interface for 1D spatial reconstruction schemes on all interfaces.
 *
 * Reconstruction schemes take the cell-centered primitive variables stored
 * in a DataLayer and build left/right primitive states at all cell interfaces
 * in a single call.
 *
 * The interface is batch-oriented and stateless: given the current DataLayer,
 * it fills two arrays of Primitive states corresponding to the left and right
 * states at every interface in the 1D layout (including ghosts).
 */
class Reconstruction {
public:
    virtual ~Reconstruction() = default;

    /**
     * @brief Reconstructs left/right primitive states at all interfaces.
     *
     * For a 1D grid with total_size cells (including ghost cells), there are
     * total_size - 1 interfaces. Interface i lies between cell i and cell i+1.
     *
     * Implementations must fill left_states(i) and right_states(i) with the
     * reconstructed primitive states at interface i for all
     * 0 ≤ i < total_size - 1. Implementations are allowed to resize the
     * provided arrays to the required length.
     *
     * @param layer        DataLayer containing current primitive fields.
     * @param left_states  Output array of left primitive states at interfaces.
     * @param right_states Output array of right primitive states at interfaces.
     */
    virtual void ReconstructStates(const DataLayer& layer,
                                   xt::xarray<Primitive>& left_states,
                                   xt::xarray<Primitive>& right_states) const = 0;
};

#endif  // RECONSTRUCTION_HPP