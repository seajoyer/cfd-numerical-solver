#ifndef HLLCRIEMANNSOLVER_HPP
#define HLLCRIEMANNSOLVER_HPP

#include "riemann/RiemannSolver.hpp"

/**
 * @class HLLCRiemannSolver
 * @brief HLLC approximate Riemann solver for generic face normals.
 *
 * Computes conservative numerical flux for Euler equations of an ideal gas
 * through one face with an arbitrary unit normal.
 */
class HLLCRiemannSolver final : public RiemannSolver {
public:
    HLLCRiemannSolver() = default;

    [[nodiscard]] ConservativeCell ComputeFlux(const PrimitiveCell& left,
                                               const PrimitiveCell& right,
                                               double gamma,
                                               const FaceNormal& normal) const override;

private:
    void SplitVelocity(const PrimitiveCell& state,
                       const FaceNormal& normal,
                       double& normal_velocity,
                       double& tangential_x,
                       double& tangential_y,
                       double& tangential_z) const;

    void BuildStarMomentum(double rho_star,
                           double star_normal_velocity,
                           double tangential_x,
                           double tangential_y,
                           double tangential_z,
                           const FaceNormal& normal,
                           double& rhoU_star,
                           double& rhoV_star,
                           double& rhoW_star) const;
};

#endif  // HLLCRIEMANNSOLVER_HPP
