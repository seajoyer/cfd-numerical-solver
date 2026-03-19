#ifndef WORKSPACE_HPP
#define WORKSPACE_HPP

#include <cstddef>
#include <xtensor.hpp>

class Mesh;

/**
 * @class Workspace
 * @brief Reusable scratch buffers for one local subdomain.
 *
 * This version is prepared for a staggered Mader / large-particle scheme.
 *
 * Cell-centered fields:
 *   Wc(var,i,j,k) : auxiliary primitive-like cell data
 *                   (rho, u_cell, v_cell, w_cell, P)
 *   rhs(var,i,j,k): conservative RHS (kept for compatibility)
 *   T(i,j,k)      : temperature
 *   I(i,j,k)      : specific internal energy
 *   Reactant(i,j,k): mass fraction of unreacted component
 *   Q(i,j,k)      : cell-centered artificial viscosity
 *   D(type,i,j,k) : transport accumulators
 *
 * Face-centered velocity fields:
 *   Ux(i_face,j,k) : x-directed velocity on x-faces, shape (sx+1, sy, sz)
 *   Vy(i,j_face,k) : y-directed velocity on y-faces, shape (sx, sy+1, sz)
 *   Wz(i,j,k_face) : z-directed velocity on z-faces, shape (sx, sy, sz+1)
 *
 * Old face-centered velocity fields:
 *   UxOld, VyOld, WzOld
 */
class Workspace final {
public:
    static constexpr std::size_t k_nvar = 5;

    // Wc(var,i,j,k) indices
    static constexpr std::size_t k_rho = 0;
    static constexpr std::size_t k_u = 1;
    static constexpr std::size_t k_v = 2;
    static constexpr std::size_t k_w = 3;
    static constexpr std::size_t k_p = 4;

    // D(type,i,j,k) indices
    static constexpr std::size_t k_ndelta = 5;
    static constexpr std::size_t k_dm = 0;
    static constexpr std::size_t k_de = 1;
    static constexpr std::size_t k_dw = 2;
    static constexpr std::size_t k_dpu = 3;
    static constexpr std::size_t k_dpv = 4;

    Workspace() = default;

    /**
     * @brief Resize all buffers to match Mesh padded sizes.
     * @details Allocates only if shape changed.
     */
    void ResizeFrom(const Mesh& mesh);

    // -------------------- cell-centered arrays --------------------

    [[nodiscard]] xt::xtensor<double, 4>& W();
    [[nodiscard]] const xt::xtensor<double, 4>& W() const;

    [[nodiscard]] xt::xtensor<double, 4>& Rhs();
    [[nodiscard]] const xt::xtensor<double, 4>& Rhs() const;

    [[nodiscard]] xt::xtensor<double, 3>& Temperature();
    [[nodiscard]] const xt::xtensor<double, 3>& Temperature() const;

    [[nodiscard]] xt::xtensor<double, 3>& InternalEnergy();
    [[nodiscard]] const xt::xtensor<double, 3>& InternalEnergy() const;

    /**
     * @brief Cell-centered artificial viscosity q(i,j,k).
     */
    [[nodiscard]] xt::xtensor<double, 3>& Q();
    [[nodiscard]] const xt::xtensor<double, 3>& Q() const;

    /**
     * @brief Transport accumulators D(type,i,j,k).
     */
    [[nodiscard]] xt::xtensor<double, 4>& D();
    [[nodiscard]] const xt::xtensor<double, 4>& D() const;

    // -------------------- face-centered velocity arrays --------------------

    /**
     * @brief X-face velocity Ux(i_face,j,k), shape (sx+1, sy, sz).
     */
    [[nodiscard]] xt::xtensor<double, 3>& Ux();
    [[nodiscard]] const xt::xtensor<double, 3>& Ux() const;

    /**
     * @brief Y-face velocity Vy(i,j_face,k), shape (sx, sy+1, sz).
     */
    [[nodiscard]] xt::xtensor<double, 3>& Vy();
    [[nodiscard]] const xt::xtensor<double, 3>& Vy() const;

    /**
     * @brief Z-face velocity Wz(i,j,k_face), shape (sx, sy, sz+1).
     * @details Included for completeness; in 2D can remain zero.
     */
    [[nodiscard]] xt::xtensor<double, 3>& Wz();
    [[nodiscard]] const xt::xtensor<double, 3>& Wz() const;

    [[nodiscard]] xt::xtensor<double, 3>& UxOld();
    [[nodiscard]] const xt::xtensor<double, 3>& UxOld() const;

    [[nodiscard]] xt::xtensor<double, 3>& VyOld();
    [[nodiscard]] const xt::xtensor<double, 3>& VyOld() const;

    [[nodiscard]] xt::xtensor<double, 3>& WzOld();
    [[nodiscard]] const xt::xtensor<double, 3>& WzOld() const;

    // -------------------- zero helpers --------------------

    void ZeroWc();
    void ZeroRhs();
    void ZeroTemperature();
    void ZeroInternalEnergy();
    void ZeroQ();
    void ZeroD();

    void ZeroUx();
    void ZeroVy();
    void ZeroWz();

    void ZeroUxOld();
    void ZeroVyOld();
    void ZeroWzOld();

    void ZeroAll();

    [[nodiscard]] bool IsAllocated() const;

private:
    // cell-centered
    xt::xtensor<double, 4> W_;
    xt::xtensor<double, 4> rhs_;
    xt::xtensor<double, 3> temperature_;
    xt::xtensor<double, 3> internal_energy_;
    xt::xtensor<double, 3> q_;
    xt::xtensor<double, 4> D_;

    // face-centered velocities
    xt::xtensor<double, 3> ux_;
    xt::xtensor<double, 3> vy_;
    xt::xtensor<double, 3> wz_;

    xt::xtensor<double, 3> ux_old_;
    xt::xtensor<double, 3> vy_old_;
    xt::xtensor<double, 3> wz_old_;

    std::size_t sx_ = 0;
    std::size_t sy_ = 0;
    std::size_t sz_ = 0;

    void Allocate(std::size_t sx, std::size_t sy, std::size_t sz);
};

#endif  // WORKSPACE_HPP
