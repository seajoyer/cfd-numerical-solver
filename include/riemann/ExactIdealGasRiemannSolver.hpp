#ifndef EXACTIDEALGASRIEMANNSOLVER_HPP
#define EXACTIDEALGASRIEMANNSOLVER_HPP

#include "riemann/RiemannSolver.hpp"

/**
 * @class ExactIdealGasRiemannSolver
 * @brief Exact Riemann solver for ideal-gas Euler equations on a generic face.
 *
 * Solves the exact 1D Riemann problem in the local normal direction of the face,
 * then reconstructs a full 3D primitive sample using tangential velocity
 * components taken from the appropriate side of the contact discontinuity.
 *
 * Intended mainly for validation and reference solutions.
 */
class ExactIdealGasRiemannSolver final : public RiemannSolver {
public:
    /**
     * @brief Construct exact solver with xi = 0 and Q = 2.
     */
    ExactIdealGasRiemannSolver();

    /**
     * @brief Construct exact solver with custom sampling coordinate and Q threshold.
     *
     * @param xi Similarity coordinate used for sampling.
     * @param Q_user Pressure-ratio threshold used in initial pressure guess selection.
     */
    ExactIdealGasRiemannSolver(double xi, double Q_user);

    /**
     * @brief Set similarity coordinate for exact-solution sampling.
     */
    void SetXi(double xi);

    /**
     * @brief Set pressure-ratio threshold used in initial pressure guess selection.
     */
    void SetQ(double Q);

    /**
     * @brief Sample exact solution at given similarity coordinate.
     *
     * @param left Owner-side primitive state.
     * @param right Neighbor-side or exterior primitive state.
     * @param gamma Ratio of specific heats.
     * @param xi Similarity coordinate.
     * @param normal Face unit normal.
     * @return Primitive sample of exact solution.
     */
    [[nodiscard]] PrimitiveCell Sample(const PrimitiveCell& left,
                                       const PrimitiveCell& right,
                                       double gamma,
                                       double xi,
                                       const FaceNormal& normal) const;

    /**
     * @brief Compute exact flux at the interface sample point xi = xi_.
     *
     * @param left Owner-side primitive state.
     * @param right Neighbor-side or exterior primitive state.
     * @param gamma Ratio of specific heats.
     * @param normal Face unit normal.
     * @return Conservative numerical flux through the face.
     */
    [[nodiscard]] ConservativeCell ComputeFlux(const PrimitiveCell& left,
                                               const PrimitiveCell& right,
                                               double gamma,
                                               const FaceNormal& normal) const override;

private:
    struct State1D final {
        double rho = 0.0;
        double un = 0.0;
        double p = 0.0;
        double a = 0.0;
    };

    struct Primitive1D final {
        double rho = 0.0;
        double un = 0.0;
        double p = 0.0;
    };

    double xi_ = 0.0;
    double Q_user_ = 2.0;

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

    void ComposeVelocityFromLocalBasis(double u_n,
                                       double u_t1,
                                       double u_t2,
                                       const FaceNormal& normal,
                                       double t1_x,
                                       double t1_y,
                                       double t1_z,
                                       double t2_x,
                                       double t2_y,
                                       double t2_z,
                                       double& u,
                                       double& v,
                                       double& w) const;

    [[nodiscard]] State1D MakeState1D(const PrimitiveCell& state,
                                      const FaceNormal& normal,
                                      double t1_x,
                                      double t1_y,
                                      double t1_z,
                                      double t2_x,
                                      double t2_y,
                                      double t2_z,
                                      double gamma) const;

    [[nodiscard]] double PhiRarefaction(double p,
                                        const State1D& state,
                                        double gamma) const;

    [[nodiscard]] double PhiRarefactionDerivative(double p,
                                                  const State1D& state,
                                                  double gamma) const;

    [[nodiscard]] double PhiShock(double p,
                                  const State1D& state,
                                  double gamma) const;

    [[nodiscard]] double PhiShockDerivative(double p,
                                            const State1D& state,
                                            double gamma) const;

    [[nodiscard]] double SidePhi(double p,
                                 const State1D& state,
                                 double gamma) const;

    [[nodiscard]] double SidePhiDerivative(double p,
                                           const State1D& state,
                                           double gamma) const;

    [[nodiscard]] double InitialGuess(const State1D& left,
                                      const State1D& right,
                                      double Q_user,
                                      double gamma) const;

    [[nodiscard]] double SolveStarPressure(const State1D& left,
                                           const State1D& right,
                                           double gamma,
                                           double Q_user) const;

    [[nodiscard]] Primitive1D SampleVacuum(double xi,
                                           const State1D& left,
                                           const State1D& right,
                                           double gamma) const;

    [[nodiscard]] Primitive1D SampleNonVacuum(double xi,
                                              double p_star,
                                              double u_star,
                                              const State1D& left,
                                              const State1D& right,
                                              double gamma) const;
};

#endif  // EXACTIDEALGASRIEMANNSOLVER_HPP
