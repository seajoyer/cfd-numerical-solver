#ifndef WENORECONSTRUCTION_HPP
#define WENORECONSTRUCTION_HPP

#include "reconstruction/Reconstruction.hpp"

/**
 * @class WENOReconstruction
 * @brief Weighted Essentially Non-Oscillatory (WENO) reconstruction.
 *
 * This class implements a family of WENO schemes for 1D uniform grids.
 * It reconstructs left and right primitive states at cell interfaces
 * from cell-centered primitive variables stored in a DataLayer.
 *
 * The nominal order of accuracy is controlled by an integer parameter
 * `order_`. In the initial implementation, the solver will support
 * orders 3 and 5 (WENO3, WENO5), using precomputed
 * linear weights and smoothness indicators. The implementation is
 * structured such that common logic for WENO stencils, smoothness
 * indicators and nonlinear weights can be factored into reusable
 * helpers, allowing future extension to arbitrary orders with
 * dynamically computed coefficients.
 *
 * The reconstruction operates on primitive variables (rho, u, P)
 * component-wise:
 *  - for each interface i+1/2, a left state w_{i+1/2}^- is built
 *    from data around cell i,
 *  - a right state w_{i+1/2}^+ is built from data around cell i+1,
 * using one-sided WENO polynomials.
 *
 * Ghost cells are assumed to be already filled by the caller; this
 * class does not perform any boundary handling itself.
 */
class WENOReconstruction : public Reconstruction {
public:
    /**
     * @brief Constructs a WENO reconstruction with a given nominal order.
     *
     * Currently supported nominal orders are 3 and 5.
     *
     * @param order Nominal order of accuracy (e.g. 3 or 5).
     */
    explicit WENOReconstruction(int order): order_(order) {}

    /**
     * @brief Virtual destructor.
     */
    ~WENOReconstruction() override = default;

    /**
     * @brief Reconstructs left/right primitive states at all interfaces.
     *
     * For a 1D grid with total_size cells (including ghost cells), there are
     * total_size - 1 interfaces. For each interface i in [0, total_size - 2],
     * this method computes:
     *
     *  - left_states(i)  ≈ w_{i+1/2}^-  (value approached from cell i),
     *  - right_states(i) ≈ w_{i+1/2}^+  (value approached from cell i+1),
     *
     * using a one-dimensional WENO reconstruction of the configured order.
     * The arrays @p left_states and @p right_states are resized internally
     * to match the number of interfaces.
     *
     * @param layer        DataLayer containing current primitive fields.
     * @param left_states  Output array of left primitive states at interfaces.
     * @param right_states Output array of right primitive states at interfaces.
     */
    void ReconstructStates(const DataLayer& layer,
                           xt::xarray<Primitive>& left_states,
                           xt::xarray<Primitive>& right_states) const override;

    /**
     * @brief Returns the current nominal WENO order.
     *
     * This is the order requested at construction (or via SetOrder())
     * after mapping to one of the internally supported orders (e.g. 3, 5).
     *
     * @return Effective nominal order used by the reconstruction.
     */
    [[nodiscard]] auto GetOrder() const -> int { return order_; }

    /**
     * @brief Sets the nominal WENO order.
     *
     * The requested order is internally mapped to one of the supported
     * implementations (e.g. WENO3, WENO5). The exact mapping and
     * clamping rules are defined in the implementation.
     *
     * @param order Nominal order of accuracy (e.g. 3 or 5).
     */
    void SetOrder(int order);

    /**
     * @brief Returns the small positive parameter ε used in nonlinear weights.
     *
     * The parameter ε appears in the definition of nonlinear WENO weights:
     *
     *  \f[
     *      \tilde{\omega}_k =
     *      \frac{d_k}{(\varepsilon + \beta_k)^p}, \quad
     *      \omega_k = \frac{\tilde{\omega}_k}{\sum_j \tilde{\omega}_j}.
     *  \f]
     *
     * It prevents division by zero and controls sensitivity to smoothness
     * indicators.
     *
     * @return Current ε value used in the weight formula.
     */
    [[nodiscard]] auto GetEpsilon() const -> double { return epsilon_; }

    /**
     * @brief Sets the small positive parameter ε used in nonlinear weights.
     *
     * @param epsilon New ε value (should be positive and typically small,
     *                e.g. 1e-6 or 1e-10).
     */
    void SetEpsilon(double epsilon);

    /**
     * @brief Returns the power p used in nonlinear WENO weights.
     *
     * The exponent p appears in the weight formula
     * (ε + β_k)^(-p). Typical choices are p = 2 or p = 1,
     * depending on the specific WENO variant.
     *
     * @return Current exponent p for nonlinear weights.
     */
    [[nodiscard]] auto GetNonlinearWeightPower() const -> int { return p_; }

    /**
     * @brief Sets the power p used in nonlinear WENO weights.
     *
     * @param p New exponent for nonlinear WENO weights. Usually p ≥ 1,
     *          with p = 2 in classical WENO-JS.
     */
    void SetNonlinearWeightPower(int p);

private:
    /**
     * @brief Nominal WENO order (mapped to one of the supported variants).
     *
     * This value determines which set of precomputed coefficients and
     * smoothness indicator formulas is used (e.g. WENO3, WENO5).
     */
    int order_{5};

    /**
     * @brief Small positive parameter ε in the nonlinear weight formula.
     */
    double epsilon_{1e-6};

    /**
     * @brief Exponent p used in the nonlinear weight formula.
     */
    int p_{2};
};

#endif  // WENORECONSTRUCTION_HPP
