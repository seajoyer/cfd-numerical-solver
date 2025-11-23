#ifndef P0RECONSTRUCTION_HPP
#define P0RECONSTRUCTION_HPP

#include "reconstruction/Reconstruction.hpp"

/**
 * @class P0Reconstruction
 * @brief Piecewise-constant reconstruction (Godunov's original scheme).
 *
 * For each interface between cell i and cell i+1, this scheme assigns
 *  - left state  = cell-average in cell i,
 *  - right state = cell-average in cell i+1.
 *
 * This yields a first-order accurate Godunov method in space.
 */
class P0Reconstruction : public Reconstruction {
public:
    P0Reconstruction() = default;
    ~P0Reconstruction() override = default;

    /**
     * @brief Computes left/right primitive states at all interfaces.
     *
     * For a 1D layout with total_size cells (including ghosts), there are
     * total_size - 1 interfaces. Interface i takes the primitive state from
     * cell i on the left and from cell i+1 on the right.
     *
     * @param layer        DataLayer with primitive fields.
     * @param left_states  Output: left primitive states at all interfaces.
     * @param right_states Output: right primitive states at all interfaces.
     */
    void ReconstructStates(const DataLayer& layer,
                           xt::xarray<Primitive>& left_states,
                           xt::xarray<Primitive>& right_states) const override;
};

#endif  // P0RECONSTRUCTION_HPP