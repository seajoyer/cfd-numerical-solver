#ifndef WORKSPACE_HPP
#define WORKSPACE_HPP

#include <cstddef>

#include <xtensor.hpp>

class Mesh;

/**
 * @brief Reusable scratch buffers for one unstructured single-process mesh.
 *
 * Cell-centered fields:
 *   W(cell, var)          : auxiliary primitive-like cell data
 *                           (rho, u, v, w, P)
 *   rhs(cell, var)        : conservative RHS
 *   temperature(cell)     : temperature
 *   internal_energy(cell) : specific internal energy
 *   q(cell)               : cell-centered artificial viscosity
 *   D(cell, type)         : transport accumulators
 *
 * Face-centered fields:
 *   ux(face)              : x-component of face-centered velocity
 *   vy(face)              : y-component of face-centered velocity
 *   wz(face)              : z-component of face-centered velocity
 *
 * Old face-centered velocity fields:
 *   ux_old(face), vy_old(face), wz_old(face)
 */
class Workspace final {
public:
    static constexpr std::size_t k_nvar = 5;

    // W(cell,var) indices
    static constexpr std::size_t k_rho = 0;
    static constexpr std::size_t k_u = 1;
    static constexpr std::size_t k_v = 2;
    static constexpr std::size_t k_w = 3;
    static constexpr std::size_t k_p = 4;

    // D(cell,type) indices
    static constexpr std::size_t k_ndelta = 5;
    static constexpr std::size_t k_dm = 0;
    static constexpr std::size_t k_de = 1;
    static constexpr std::size_t k_dw = 2;
    static constexpr std::size_t k_dpu = 3;
    static constexpr std::size_t k_dpv = 4;

    Workspace() = default;

    /**
     * @brief Resize all buffers to match mesh cell and face counts.
     * @details Allocates only if shape changed.
     */
    void ResizeFrom(const Mesh& mesh);

    // -------------------- cell-centered arrays --------------------

    [[nodiscard]] xt::xtensor<double, 2>& W();
    [[nodiscard]] const xt::xtensor<double, 2>& W() const;

    [[nodiscard]] xt::xtensor<double, 2>& Rhs();
    [[nodiscard]] const xt::xtensor<double, 2>& Rhs() const;

    [[nodiscard]] xt::xtensor<double, 1>& Temperature();
    [[nodiscard]] const xt::xtensor<double, 1>& Temperature() const;

    [[nodiscard]] xt::xtensor<double, 1>& InternalEnergy();
    [[nodiscard]] const xt::xtensor<double, 1>& InternalEnergy() const;

    /**
     * @brief Cell-centered artificial viscosity q(cell).
     */
    [[nodiscard]] xt::xtensor<double, 1>& Q();
    [[nodiscard]] const xt::xtensor<double, 1>& Q() const;

    /**
     * @brief Transport accumulators D(cell,type).
     */
    [[nodiscard]] xt::xtensor<double, 2>& D();
    [[nodiscard]] const xt::xtensor<double, 2>& D() const;

    // -------------------- face-centered velocity arrays --------------------

    /**
     * @brief X-component of face-centered velocity. Shape (n_faces).
     */
    [[nodiscard]] xt::xtensor<double, 1>& Ux();
    [[nodiscard]] const xt::xtensor<double, 1>& Ux() const;

    /**
     * @brief Y-component of face-centered velocity. Shape (n_faces).
     */
    [[nodiscard]] xt::xtensor<double, 1>& Vy();
    [[nodiscard]] const xt::xtensor<double, 1>& Vy() const;

    /**
     * @brief Z-component of face-centered velocity. Shape (n_faces).
     */
    [[nodiscard]] xt::xtensor<double, 1>& Wz();
    [[nodiscard]] const xt::xtensor<double, 1>& Wz() const;

    [[nodiscard]] xt::xtensor<double, 1>& UxOld();
    [[nodiscard]] const xt::xtensor<double, 1>& UxOld() const;

    [[nodiscard]] xt::xtensor<double, 1>& VyOld();
    [[nodiscard]] const xt::xtensor<double, 1>& VyOld() const;

    [[nodiscard]] xt::xtensor<double, 1>& WzOld();
    [[nodiscard]] const xt::xtensor<double, 1>& WzOld() const;

    // -------------------- zero helpers --------------------

    void ZeroW();
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

    [[nodiscard]] std::size_t GetCellCount() const;
    [[nodiscard]] std::size_t GetFaceCount() const;

private:
    // cell-centered
    xt::xtensor<double, 2> W_;
    xt::xtensor<double, 2> rhs_;
    xt::xtensor<double, 1> temperature_;
    xt::xtensor<double, 1> internal_energy_;
    xt::xtensor<double, 1> q_;
    xt::xtensor<double, 2> D_;

    // face-centered velocities
    xt::xtensor<double, 1> ux_;
    xt::xtensor<double, 1> vy_;
    xt::xtensor<double, 1> wz_;

    xt::xtensor<double, 1> ux_old_;
    xt::xtensor<double, 1> vy_old_;
    xt::xtensor<double, 1> wz_old_;

    std::size_t n_cells_ = 0;
    std::size_t n_faces_ = 0;

    void Allocate(std::size_t n_cells, std::size_t n_faces);
};

#endif  // WORKSPACE_HPP
