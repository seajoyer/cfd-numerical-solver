#include "riemann/RoeRiemannSolver.hpp"

#include <algorithm>
#include <cmath>

#include "riemann/RiemannHelpers.hpp"
#include "data/Variables.hpp"  // NormalVelocity, PrimitiveToConservative, SoundSpeed, EulerFlux

namespace {

// Harten–Hyman entropy fix on acoustic eigenvalues (same spirit as your old code).
inline double EntropyFix(const double lambda, const double lambda_L, const double lambda_R) {
    const double delta = std::max(0.0, lambda_R - lambda_L);
    const double abs_lambda = std::abs(lambda);

    if (delta <= 0.0) return abs_lambda;
    if (abs_lambda >= delta) return abs_lambda;

    return 0.5 * (lambda * lambda / delta + delta);
}

}  // namespace

auto RoeRiemannSolver::ComputeFlux(const PrimitiveCell& left,
                                  const PrimitiveCell& right,
                                  const double gamma,
                                  const Axis axis) const -> FluxCell {
    // Work in a rotated 1D frame: (un, ut1, ut2) aligned with interface normal.
    double unL, ut1L, ut2L;
    double unR, ut1R, ut2R;
    riemann::SplitVelocity(left,  axis, unL, ut1L, ut2L);
    riemann::SplitVelocity(right, axis, unR, ut1R, ut2R);

    // Build "rotated" primitive states where:
    // u=un, v=ut1, w=ut2, and we compute flux along Axis::X (normal direction).
    PrimitiveCell WL{left.rho,  unL, ut1L, ut2L, left.P};
    PrimitiveCell WR{right.rho, unR, ut1R, ut2R, right.P};

    const FluxCell FL_rot = EulerFlux(WL, gamma, Axis::X);
    const FluxCell FR_rot = EulerFlux(WR, gamma, Axis::X);

    // Conservative states (rotated frame)
    double rhoL, rhoUnL, rhoUt1L, rhoUt2L, EL;
    double rhoR, rhoUnR, rhoUt1R, rhoUt2R, ER;
    PrimitiveToConservative(WL, gamma, rhoL, rhoUnL, rhoUt1L, rhoUt2L, EL);
    PrimitiveToConservative(WR, gamma, rhoR, rhoUnR, rhoUt1R, rhoUt2R, ER);

    // Floors to avoid nonphysical inputs exploding the solver.
    const double rho_floor = 1e-14;
    rhoL = std::max(rhoL, rho_floor);
    rhoR = std::max(rhoR, rho_floor);

    const double pL = WL.P;
    const double pR = WR.P;

    // Roe averages
    const double sqL = std::sqrt(rhoL);
    const double sqR = std::sqrt(rhoR);
    const double denom = sqL + sqR;

    const double un_tilde  = (sqL * unL  + sqR * unR ) / denom;
    const double ut1_tilde = (sqL * ut1L + sqR * ut1R) / denom;
    const double ut2_tilde = (sqL * ut2L + sqR * ut2R) / denom;

    const double HL = (EL + pL) / rhoL;
    const double HR = (ER + pR) / rhoR;
    const double H_tilde = (sqL * HL + sqR * HR) / denom;

    const double q2 = 0.5 * (un_tilde*un_tilde + ut1_tilde*ut1_tilde + ut2_tilde*ut2_tilde);
    double a2_tilde = (gamma - 1.0) * (H_tilde - q2);
    // floor to keep stable
    a2_tilde = std::max(a2_tilde, 1e-14);
    const double a_tilde = std::sqrt(a2_tilde);

    const double rho_tilde = sqL * sqR;

    // Jumps in primitive (rotated)
    const double drho = rhoR - rhoL;
    const double dun  = unR - unL;
    const double dut1 = ut1R - ut1L;
    const double dut2 = ut2R - ut2L;
    const double dp   = pR - pL;

    // Wave strengths (5-wave Roe decomposition: 2 acoustic + 1 entropy + 2 shear)
    const double alpha2 = drho - dp / a2_tilde;  // entropy/contact density part
    const double alpha1 = 0.5 / a2_tilde * (dp - rho_tilde * a_tilde * dun);
    const double alpha3 = 0.5 / a2_tilde * (dp + rho_tilde * a_tilde * dun);
    const double alpha4 = rho_tilde * dut1;      // shear-1
    const double alpha5 = rho_tilde * dut2;      // shear-2

    // Eigenvalues (with entropy fix on acoustic ones)
    const double lambda1 = un_tilde - a_tilde;
    const double lambda2 = un_tilde;
    const double lambda3 = un_tilde + a_tilde;
    const double lambda4 = un_tilde;
    const double lambda5 = un_tilde;

    const double aL = SoundSpeed(WL, gamma);
    const double aR = SoundSpeed(WR, gamma);

    const double lambda1_L = unL - aL;
    const double lambda1_R = unR - aR;
    const double lambda3_L = unL + aL;
    const double lambda3_R = unR + aR;

    const double a1 = EntropyFix(lambda1, lambda1_L, lambda1_R);
    const double a2 = std::abs(lambda2);
    const double a3 = EntropyFix(lambda3, lambda3_L, lambda3_R);
    const double a4 = std::abs(lambda4);
    const double a5 = std::abs(lambda5);

    // Right eigenvectors in rotated frame (rho, rho*un, rho*ut1, rho*ut2, E)
    const double r1_0 = 1.0;
    const double r1_1 = un_tilde - a_tilde;
    const double r1_2 = ut1_tilde;
    const double r1_3 = ut2_tilde;
    const double r1_4 = H_tilde - un_tilde * a_tilde;

    const double r2_0 = 1.0;
    const double r2_1 = un_tilde;
    const double r2_2 = ut1_tilde;
    const double r2_3 = ut2_tilde;
    const double r2_4 = q2;

    const double r3_0 = 1.0;
    const double r3_1 = un_tilde + a_tilde;
    const double r3_2 = ut1_tilde;
    const double r3_3 = ut2_tilde;
    const double r3_4 = H_tilde + un_tilde * a_tilde;

    // Shear eigenvectors (advected with u~)
    const double r4_0 = 0.0;
    const double r4_1 = 0.0;
    const double r4_2 = 1.0;
    const double r4_3 = 0.0;
    const double r4_4 = ut1_tilde;

    const double r5_0 = 0.0;
    const double r5_1 = 0.0;
    const double r5_2 = 0.0;
    const double r5_3 = 1.0;
    const double r5_4 = ut2_tilde;

    // Dissipation term dF = sum |lambda_k| alpha_k r_k
    const double dF0 = a1 * alpha1 * r1_0 +
                       a2 * alpha2 * r2_0 +
                       a3 * alpha3 * r3_0 +
                       a4 * alpha4 * r4_0 +
                       a5 * alpha5 * r5_0;

    const double dF1 = a1 * alpha1 * r1_1 +
                       a2 * alpha2 * r2_1 +
                       a3 * alpha3 * r3_1 +
                       a4 * alpha4 * r4_1 +
                       a5 * alpha5 * r5_1;

    const double dF2 = a1 * alpha1 * r1_2 +
                       a2 * alpha2 * r2_2 +
                       a3 * alpha3 * r3_2 +
                       a4 * alpha4 * r4_2 +
                       a5 * alpha5 * r5_2;

    const double dF3 = a1 * alpha1 * r1_3 +
                       a2 * alpha2 * r2_3 +
                       a3 * alpha3 * r3_3 +
                       a4 * alpha4 * r4_3 +
                       a5 * alpha5 * r5_3;

    const double dF4 = a1 * alpha1 * r1_4 +
                       a2 * alpha2 * r2_4 +
                       a3 * alpha3 * r3_4 +
                       a4 * alpha4 * r4_4 +
                       a5 * alpha5 * r5_4;

    // Roe flux in rotated frame
    FluxCell F_rot;
    F_rot.mass   = 0.5 * (FL_rot.mass   + FR_rot.mass)   - 0.5 * dF0;
    F_rot.mom_x  = 0.5 * (FL_rot.mom_x  + FR_rot.mom_x)  - 0.5 * dF1; // mom_n
    F_rot.mom_y  = 0.5 * (FL_rot.mom_y  + FR_rot.mom_y)  - 0.5 * dF2; // mom_t1
    F_rot.mom_z  = 0.5 * (FL_rot.mom_z  + FR_rot.mom_z)  - 0.5 * dF3; // mom_t2
    F_rot.energy = 0.5 * (FL_rot.energy + FR_rot.energy) - 0.5 * dF4;

    // Map rotated momentum flux components (mom_n, mom_t1, mom_t2) back to global (x,y,z)
    FluxCell F;
    F.mass = F_rot.mass;
    F.energy = F_rot.energy;

    riemann::ComposeVector(F_rot.mom_x, F_rot.mom_y, F_rot.mom_z, axis,
                           F.mom_x, F.mom_y, F.mom_z);

    return F;
}