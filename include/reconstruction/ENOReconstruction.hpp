#ifndef ENORECONSTRUCTION_HPP
#define ENORECONSTRUCTION_HPP

#include "reconstruction/Reconstruction.hpp"

/**
 * @class ENOReconstruction
 * @brief Essentially Non-Oscillatory reconstruction of arbitrary order r (runtime).
 *
 * Reconstructs primitive variables W=(rho,u,v,w,P) at a face using an ENO stencil
 * selected by the divided-difference smoothness criterion (uniform grid in index space).
 *
 * Face indexing convention:
 *  - (i,j,k) is the LEFT cell adjacent to the face.
 *  - RIGHT cell is +1 along axis.
 *
 * For each face:
 *  - WL is reconstructed from base cell (i,j,k) evaluated at x = base + 0.5
 *  - WR is reconstructed from base cell (i+1 along axis) evaluated at x = base - 0.5
 *
 * Contract:
 *  - No allocations in the hot path after first call (thread-local scratch buffers).
 */
class ENOReconstruction final : public Reconstruction {
public:
    /** @brief Constructs ENO of given order r (r >= 1). */
    explicit ENOReconstruction(int order);

    ~ENOReconstruction() override = default;

    [[nodiscard]] int GetOrder() const { return order_; }
    void SetOrder(int order);

    void ReconstructFace(const xt::xtensor<double, 4>& W,
                         Axis axis,
                         int i, int j, int k,
                         PrimitiveCell& WL,
                         PrimitiveCell& WR) const override;

private:
    int order_ = 1;
};

#endif  // ENORECONSTRUCTION_HPP