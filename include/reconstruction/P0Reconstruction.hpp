#ifndef P0RECONSTRUCTION_HPP
#define P0RECONSTRUCTION_HPP

#include "reconstruction/Reconstruction.hpp"

/**
 * @class P0Reconstruction
 * @brief Piecewise-constant reconstruction (first order).
 *
 * For each face:
 *  - WL = W in left cell
 *  - WR = W in right cell
 */
class P0Reconstruction final : public Reconstruction {
public:
    P0Reconstruction() = default;
    ~P0Reconstruction() override = default;

    void ReconstructFace(const xt::xtensor<double, 4>& W,
                         Axis axis,
                         int i, int j, int k,
                         PrimitiveCell& WL,
                         PrimitiveCell& WR) const override;
};

#endif  // P0RECONSTRUCTION_HPP