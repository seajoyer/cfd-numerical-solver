#ifndef VNRARTIFICIALVISCOSITY_HPP
#define VNRARTIFICIALVISCOSITY_HPP

#include <cmath>

#include "viscosity/ArtificialViscosity.hpp"
#include "data/Variables.hpp"
#include "config/Settings.hpp"

/**
 * @class VNRArtificialViscosity
 * @brief von Neumann–Richtmyer artificial viscosity (axis-split, face-based).
 *
 * For each axis, on each face (indexed by its left cell):
 *   du = u_n(R) - u_n(L)
 *   if du < 0 (compression):
 *     q = rho_bar * ( C2 * du^2 + C1 * c_bar * |du| )
 *   else q = 0
 *
 * Viscosity flux correction along an axis:
 *   momentum_normal += q
 *   energy         += u_face * q
 *
 * and is added to RHS as a flux divergence.
 */
class VNRArtificialViscosity final : public ArtificialViscosity {
public:
    explicit VNRArtificialViscosity(const Settings& settings,
                                    double C1 = 1.5,
                                    double C2 = 6.0);

    void AddToRhs(const DataLayer& layer,
                  const xt::xtensor<double, 4>& W,
                  double gamma,
                  double dt,
                  xt::xtensor<double, 4>& rhs) const override;

    [[nodiscard]] int GetRequiredPadding() const override { return 1; }

private:
    double C1_;
    double C2_;
    Settings settings_;

    // Cached face buffers (mutable to allow resizing in const AddToRhs)
    mutable xt::xtensor<double, 3> qx_; // (sx-1, sy, sz)
    mutable xt::xtensor<double, 3> qy_; // (sx, sy-1, sz)
    mutable xt::xtensor<double, 3> qz_; // (sx, sy, sz-1)
    mutable int sx_ = 0;
    mutable int sy_ = 0;
    mutable int sz_ = 0;

    void ResizeFrom(const DataLayer& layer) const;

    [[nodiscard]] xt::xtensor<double, 3>& QFace(Axis axis) const;
    [[nodiscard]] const xt::xtensor<double, 1>& InvMetric(const DataLayer& layer, Axis axis) const;

    [[nodiscard]] static std::size_t MomVarIndex(Axis axis);
    [[nodiscard]] static std::size_t VelVarIndex(Axis axis);

    void ComputeQFaces(const DataLayer& layer,
                       const xt::xtensor<double, 4>& W,
                       double gamma,
                       Axis axis) const;

    void AddAxisContribution(const DataLayer& layer,
                             const xt::xtensor<double, 4>& W,
                             xt::xtensor<double, 4>& rhs,
                             Axis axis) const;

    static PrimitiveCell LoadPrimitive(const xt::xtensor<double, 4>& W, int i, int j, int k);
};

#endif  // VNRARTIFICIALVISCOSITY_HPP