#ifndef PRESSUREVELOCITYWORKSPACE_HPP
#define PRESSUREVELOCITYWORKSPACE_HPP

#include <cstddef>
#include <xtensor.hpp>

class Mesh;

/**
 * @class PressureVelocityWorkspace
 * @brief Reusable scratch buffers for pressure-velocity coupling solvers
 *        on a structured staggered grid.
 *
 * Cell-centered:
 *  - p_corr(i,j,k)       : pressure correction
 *  - pressure_rhs(i,j,k) : RHS of pressure correction equation
 *
 * Face-centered:
 *  - ux_star / vy_star / wz_star : predicted velocities
 *  - aPu / aPv / aPw             : diagonal momentum coefficients
 *
 * Optional helper fields for pressure equation assembly:
 *  - apE/apW/apN/apS/apT/apB/apP : pressure matrix coefficients
 */
class PressureVelocityWorkspace final {
public:
    PressureVelocityWorkspace() = default;

    /**
     * @brief Resize all buffers to match padded mesh sizes.
     * @details Reallocates only if shapes changed.
     */
    void ResizeFrom(const Mesh& mesh);

    // -------------------- cell-centered --------------------

    [[nodiscard]] xt::xtensor<double, 3>& PressureCorrection();
    [[nodiscard]] const xt::xtensor<double, 3>& PressureCorrection() const;

    [[nodiscard]] xt::xtensor<double, 3>& PressureRhs();
    [[nodiscard]] const xt::xtensor<double, 3>& PressureRhs() const;

    [[nodiscard]] xt::xtensor<double, 3>& ApE();
    [[nodiscard]] const xt::xtensor<double, 3>& ApE() const;

    [[nodiscard]] xt::xtensor<double, 3>& ApW();
    [[nodiscard]] const xt::xtensor<double, 3>& ApW() const;

    [[nodiscard]] xt::xtensor<double, 3>& ApN();
    [[nodiscard]] const xt::xtensor<double, 3>& ApN() const;

    [[nodiscard]] xt::xtensor<double, 3>& ApS();
    [[nodiscard]] const xt::xtensor<double, 3>& ApS() const;

    [[nodiscard]] xt::xtensor<double, 3>& ApT();
    [[nodiscard]] const xt::xtensor<double, 3>& ApT() const;

    [[nodiscard]] xt::xtensor<double, 3>& ApB();
    [[nodiscard]] const xt::xtensor<double, 3>& ApB() const;

    [[nodiscard]] xt::xtensor<double, 3>& ApP();
    [[nodiscard]] const xt::xtensor<double, 3>& ApP() const;

    // -------------------- face-centered --------------------

    [[nodiscard]] xt::xtensor<double, 3>& UxStar();
    [[nodiscard]] const xt::xtensor<double, 3>& UxStar() const;

    [[nodiscard]] xt::xtensor<double, 3>& VyStar();
    [[nodiscard]] const xt::xtensor<double, 3>& VyStar() const;

    [[nodiscard]] xt::xtensor<double, 3>& WzStar();
    [[nodiscard]] const xt::xtensor<double, 3>& WzStar() const;

    [[nodiscard]] xt::xtensor<double, 3>& APu();
    [[nodiscard]] const xt::xtensor<double, 3>& APu() const;

    [[nodiscard]] xt::xtensor<double, 3>& APv();
    [[nodiscard]] const xt::xtensor<double, 3>& APv() const;

    [[nodiscard]] xt::xtensor<double, 3>& APw();
    [[nodiscard]] const xt::xtensor<double, 3>& APw() const;

    // -------------------- helpers --------------------

    void ZeroPressureCorrection();
    void ZeroPressureRhs();

    void ZeroPressureCoefficients();

    void ZeroUxStar();
    void ZeroVyStar();
    void ZeroWzStar();

    void ZeroAPu();
    void ZeroAPv();
    void ZeroAPw();

    void ZeroAll();

    [[nodiscard]] bool IsAllocated() const;

    [[nodiscard]] std::size_t GetSx() const;
    [[nodiscard]] std::size_t GetSy() const;
    [[nodiscard]] std::size_t GetSz() const;

    [[nodiscard]] xt::xtensor<double, 3>& AuE();
    [[nodiscard]] const xt::xtensor<double, 3>& AuE() const;

    [[nodiscard]] xt::xtensor<double, 3>& AuW();
    [[nodiscard]] const xt::xtensor<double, 3>& AuW() const;

    [[nodiscard]] xt::xtensor<double, 3>& AuN();
    [[nodiscard]] const xt::xtensor<double, 3>& AuN() const;

    [[nodiscard]] xt::xtensor<double, 3>& AuS();
    [[nodiscard]] const xt::xtensor<double, 3>& AuS() const;

    [[nodiscard]] xt::xtensor<double, 3>& AvE();
    [[nodiscard]] const xt::xtensor<double, 3>& AvE() const;

    [[nodiscard]] xt::xtensor<double, 3>& AvW();
    [[nodiscard]] const xt::xtensor<double, 3>& AvW() const;

    [[nodiscard]] xt::xtensor<double, 3>& AvN();
    [[nodiscard]] const xt::xtensor<double, 3>& AvN() const;

    [[nodiscard]] xt::xtensor<double, 3>& AvS();
    [[nodiscard]] const xt::xtensor<double, 3>& AvS() const;

    void ZeroMomentumCoefficients();

private:
    xt::xtensor<double, 3> pressure_correction_;
    xt::xtensor<double, 3> pressure_rhs_;

    xt::xtensor<double, 3> ap_e_;
    xt::xtensor<double, 3> ap_w_;
    xt::xtensor<double, 3> ap_n_;
    xt::xtensor<double, 3> ap_s_;
    xt::xtensor<double, 3> ap_t_;
    xt::xtensor<double, 3> ap_b_;
    xt::xtensor<double, 3> ap_p_;

    xt::xtensor<double, 3> ux_star_;
    xt::xtensor<double, 3> vy_star_;
    xt::xtensor<double, 3> wz_star_;

    xt::xtensor<double, 3> a_pu_;
    xt::xtensor<double, 3> a_pv_;
    xt::xtensor<double, 3> a_pw_;

    xt::xtensor<double, 3> au_e_;
    xt::xtensor<double, 3> au_w_;
    xt::xtensor<double, 3> au_n_;
    xt::xtensor<double, 3> au_s_;

    xt::xtensor<double, 3> av_e_;
    xt::xtensor<double, 3> av_w_;
    xt::xtensor<double, 3> av_n_;
    xt::xtensor<double, 3> av_s_;

    std::size_t sx_ = 0;
    std::size_t sy_ = 0;
    std::size_t sz_ = 0;

    void Allocate(std::size_t sx, std::size_t sy, std::size_t sz);
};

#endif  // PRESSUREVELOCITYWORKSPACE_HPP
