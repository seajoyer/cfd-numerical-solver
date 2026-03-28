#ifndef ROERIEMANNSOLVER_HPP
#define ROERIEMANNSOLVER_HPP

#include "riemann/RiemannSolver.hpp"

/**
 * @class RoeRiemannSolver
 * @brief Roe approximate Riemann solver for Euler equations (ideal gas).
 *
 * Works for a generic face unit normal by rotating the interface problem into
 * a local orthonormal basis:
 * - e_n  : face normal
 * - e_t1 : first tangential direction
 * - e_t2 : second tangential direction
 *
 * Roe linearization is then applied in the local normal direction.
 */
class RoeRiemannSolver final : public RiemannSolver {
public:
    RoeRiemannSolver() = default;

    [[nodiscard]] ConservativeCell ComputeFlux(const PrimitiveCell& left,
                                               const PrimitiveCell& right,
                                               double gamma,
                                               const FaceNormal& normal) const override;

private:
    [[nodiscard]] double EntropyFix(double lambda,
                                    double lambda_left,
                                    double lambda_right) const;

    void BuildTangentialBasis(const FaceNormal& normal,
                              double& t1_x,
                              double& t1_y,
                              double& t1_z,
                              double& t2_x,
                              double& t2_y,
                              double& t2_z) const;

    void ProjectVelocityToLocalBasis(const PrimitiveCell& state,
                                     const FaceNormal& normal,
                                     double t1_x,
                                     double t1_y,
                                     double t1_z,
                                     double t2_x,
                                     double t2_y,
                                     double t2_z,
                                     double& u_n,
                                     double& u_t1,
                                     double& u_t2) const;

    void ComposeVectorFromLocalBasis(double v_n,
                                     double v_t1,
                                     double v_t2,
                                     const FaceNormal& normal,
                                     double t1_x,
                                     double t1_y,
                                     double t1_z,
                                     double t2_x,
                                     double t2_y,
                                     double t2_z,
                                     double& v_x,
                                     double& v_y,
                                     double& v_z) const;
};

#endif  // ROERIEMANNSOLVER_HPP
