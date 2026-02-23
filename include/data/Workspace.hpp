#ifndef WORKSPACE_HPP
#define WORKSPACE_HPP

#include <cstddef>
#include <xtensor.hpp>

class DataLayer;

/**
 * @class Workspace
 * @brief Reusable scratch buffers for one local subdomain.
 *
 * Workspace stores temporary arrays reused across time steps / stages.
 * Minimal buffers for the smoke-test pipeline:
 *  - W   : primitive variables W(var,i,j,k) = (rho,u,v,w,P)
 *  - rhs : RHS in conservative form rhs(var,i,j,k)
 *
 * Shapes are always (5, sx, sy, sz) matching DataLayer total sizes (including ghosts).
 * No other state (e.g., U_stage) is stored here.
 */
class Workspace final {
public:
    /** @brief Number of variables in W and rhs. */
    static constexpr std::size_t k_nvar = 5;

    Workspace() = default;

    /**
     * @brief Resize buffers to match DataLayer total sizes (including ghosts).
     * @details Allocates only if shape changed. If shape unchanged, does nothing.
     */
    void ResizeFrom(const DataLayer& layer);

    /** @brief Return primitive buffer W(var,i,j,k). */
    [[nodiscard]] xt::xtensor<double, 4>& W();
    [[nodiscard]] const xt::xtensor<double, 4>& W() const;

    /** @brief Return RHS buffer rhs(var,i,j,k). */
    [[nodiscard]] xt::xtensor<double, 4>& Rhs();
    [[nodiscard]] const xt::xtensor<double, 4>& Rhs() const;

    /** @brief Fill rhs with zeros (no allocations). */
    void ZeroRhs();

    /** @brief Fill W with zeros (no allocations). */
    void ZeroW();

    /** @brief Check whether buffers are allocated with correct rank (4). */
    [[nodiscard]] bool IsAllocated() const;

private:
    xt::xtensor<double, 4> W_;
    xt::xtensor<double, 4> rhs_;

    std::size_t sx_ = 0;
    std::size_t sy_ = 0;
    std::size_t sz_ = 0;

    void Allocate(std::size_t sx, std::size_t sy, std::size_t sz);
};

#endif  // WORKSPACE_HPP
