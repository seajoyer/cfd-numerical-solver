#ifndef ENORECONSTRUCTION_HPP
#define ENORECONSTRUCTION_HPP

#include "reconstruction/Reconstruction.hpp"
#include "data/DataLayer.hpp"
#include "data/Variables.hpp"

/**
 * @class ENOReconstruction
 * @brief Essentially Non-Oscillatory reconstruction of arbitrary order.
 *
 * This class implements an ENO reconstruction on a 1D uniform grid for
 * the primitive variables stored in a DataLayer. The formal order of
 * accuracy in space is controlled by the integer parameter @p order.
 *
 * For a given order r, interface values are reconstructed using stencils
 * of size r chosen adaptively according to the ENO smoothness criterion.
 * The same stencil is applied component-wise to all primitive variables.
 */
class ENOReconstruction : public Reconstruction {
public:
    /**
     * @brief Constructs an ENO reconstruction of the given order.
     *
     * @param order Formal ENO reconstruction order (r ≥ 1).
     */
    explicit ENOReconstruction(int order);

    /**
     * @brief Virtual destructor.
     */
    ~ENOReconstruction() override = default;

    /**
     * @brief Returns the current ENO order.
     *
     * @return Integer order r.
     */
    [[nodiscard]] auto GetOrder() const -> int { return order_; }

    /**
     * @brief Sets the ENO reconstruction order.
     *
     * @param order New integer order r.
     */
    void SetOrder(int order);

    /**
     * @brief Reconstructs left/right primitive states at all interfaces.
     *
     * For a 1D layer of total size @c M , this method fills
     *
     *  - @p left_states(i)  with the left state at interface i+1/2,
     *  - @p right_states(i) with the right state at interface i+1/2,
     *
     * for i = 0,1,...,M-2. The arrays @p left_states and @p right_states
     * are resized internally to match the number of interfaces.
     *
     * @param layer        DataLayer with primitive fields.
     * @param left_states  Output array of left interface states.
     * @param right_states Output array of right interface states.
     */
    void ReconstructStates(const DataLayer& layer,
                           xt::xarray<Primitive>& left_states,
                           xt::xarray<Primitive>& right_states) const override;

private:
    int order_;
};

#endif  // ENORECONSTRUCTION_HPP