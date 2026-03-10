#ifndef VNRARTIFICIALVISCOSITY_HPP
#define VNRARTIFICIALVISCOSITY_HPP

#include <cmath>

#include "config/Settings.hpp"
#include "data/Variables.hpp"
#include "viscosity/ArtificialViscosity.hpp"

/**
 * @class VNRArtificialViscosity
 * @brief von Neumann-Richtmyer artificial viscosity (axis-split, face-based).
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
class DataLayer;
class Mesh;

class VNRArtificialViscosity final : public ArtificialViscosity {
public:
    explicit VNRArtificialViscosity(const Settings& settings,
                                    double C1 = 1.5,
                                    double C2 = 6.0);

    void AddToRhs(const DataLayer& layer,
                  const Mesh& mesh,
                  const xt::xtensor<double, 4>& W,
                  double gamma,
                  double dt,
                  xt::xtensor<double, 4>& rhs) const override;

    [[nodiscard]] int GetRequiredPadding() const override { return 1; }

private:
    double C1_;
    double C2_;
    Settings settings_;

    mutable xt::xtensor<double, 3> qx_;
    mutable xt::xtensor<double, 3> qy_;
    mutable xt::xtensor<double, 3> qz_;
    mutable int sx_ = 0;
    mutable int sy_ = 0;
    mutable int sz_ = 0;

    void ResizeFrom(const Mesh& mesh) const;

    [[nodiscard]] xt::xtensor<double, 3>& QFace(Axis axis) const;
    [[nodiscard]] const xt::xtensor<double, 1>& InvMetric(const Mesh& mesh, Axis axis) const;

    [[nodiscard]] static std::size_t MomVarIndex(Axis axis);
    [[nodiscard]] static std::size_t VelVarIndex(Axis axis);

    void ComputeQFaces(const Mesh& mesh,
                       const xt::xtensor<double, 4>& W,
                       double gamma,
                       Axis axis) const;

    void AddAxisContribution(const Mesh& mesh,
                             const xt::xtensor<double, 4>& W,
                             xt::xtensor<double, 4>& rhs,
                             Axis axis) const;

    static PrimitiveCell LoadPrimitive(const xt::xtensor<double, 4>& W, int i, int j, int k);
};

#endif  // VNRARTIFICIALVISCOSITY_HPP