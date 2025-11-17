#ifndef P0RECONSTRUCTION_HPP
#define P0RECONSTRUCTION_HPP

#include "reconstruction/Reconstruction.hpp"

/**
 * @class P0Reconstruction
 * @brief Piecewise-constant reconstruction (Godunov's original scheme).
 *
 * For each interface between cell i and cell i+1, this scheme assigns
 *
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
     * @brief Computes left/right primitive states at interface i.
     *
     * Simply reads the primitive states from cell i (left) and cell i+1
     * (right) using DataLayer::GetPrimitive().
     *
     * @param layer          DataLayer with primitive fields.
     * @param interface_index Interface index (0 ≤ i < totalSize - 1).
     * @param left_state      Output: primitive state from cell i.
     * @param right_state     Output: primitive state from cell i+1.
     */
    void ComputeInterfaceStates(const DataLayer& layer,
                                int interface_index,
                                Primitive& left_state,
                                Primitive& right_state) const override;
};

#endif  // P0RECONSTRUCTION_HPP