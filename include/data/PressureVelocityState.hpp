#ifndef PRESSUREVELOCITYSTATE_HPP
#define PRESSUREVELOCITYSTATE_HPP

#include <cstddef>
#include <xtensor.hpp>

class Mesh;
// CHECK: STAGGERED_GRID
/**
 * @class PressureVelocityState
 * @brief Persistent solution fields for pressure-velocity coupling solvers
 *        on a structured staggered grid.
 *
 * Storage:
 *  - p(i,j,k)        : cell-centered pressure, shape (sx, sy, sz)
 *  - ux(i_face,j,k)  : x-face velocity,    shape (sx+1, sy, sz)
 *  - vy(i,j_face,k)  : y-face velocity,    shape (sx, sy+1, sz)
 *  - wz(i,j,k_face)  : z-face velocity,    shape (sx, sy, sz+1)
 *
 * Old fields are kept for transient algorithms (PISO / PIMPLE).
 */
class PressureVelocityState final {
public:
    PressureVelocityState() = default;

    /**
     * @brief Resize all fields to match padded mesh sizes.
     * @details Reallocates only if shapes changed.
     */
    void ResizeFrom(const Mesh& mesh);

    // -------------------- pressure --------------------

    [[nodiscard]] xt::xtensor<double, 3>& Pressure();
    [[nodiscard]] const xt::xtensor<double, 3>& Pressure() const;

    [[nodiscard]] xt::xtensor<double, 3>& PressureOld();
    [[nodiscard]] const xt::xtensor<double, 3>& PressureOld() const;

    // -------------------- face velocities --------------------

    [[nodiscard]] xt::xtensor<double, 3>& Ux();
    [[nodiscard]] const xt::xtensor<double, 3>& Ux() const;

    [[nodiscard]] xt::xtensor<double, 3>& Vy();
    [[nodiscard]] const xt::xtensor<double, 3>& Vy() const;

    [[nodiscard]] xt::xtensor<double, 3>& Wz();
    [[nodiscard]] const xt::xtensor<double, 3>& Wz() const;

    [[nodiscard]] xt::xtensor<double, 3>& UxOld();
    [[nodiscard]] const xt::xtensor<double, 3>& UxOld() const;

    [[nodiscard]] xt::xtensor<double, 3>& VyOld();
    [[nodiscard]] const xt::xtensor<double, 3>& VyOld() const;

    [[nodiscard]] xt::xtensor<double, 3>& WzOld();
    [[nodiscard]] const xt::xtensor<double, 3>& WzOld() const;

    // -------------------- helpers --------------------

    void CopyCurrentToOld();

    void ZeroPressure();
    void ZeroPressureOld();

    void ZeroUx();
    void ZeroVy();
    void ZeroWz();

    void ZeroUxOld();
    void ZeroVyOld();
    void ZeroWzOld();

    void ZeroAll();

    [[nodiscard]] bool IsAllocated() const;

    [[nodiscard]] std::size_t GetSx() const;
    [[nodiscard]] std::size_t GetSy() const;
    [[nodiscard]] std::size_t GetSz() const;

private:
    xt::xtensor<double, 3> pressure_;
    xt::xtensor<double, 3> pressure_old_;

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

#endif  // PRESSUREVELOCITYSTATE_HPP