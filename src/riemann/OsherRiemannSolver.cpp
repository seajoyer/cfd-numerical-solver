#include "riemann/OsherRiemannSolver.hpp"

#include <algorithm>
#include <cmath>

#include "data/Variables.hpp"
#include "riemann/RiemannHelpers.hpp"

bool OsherRiemannSolver::IsPhysical(const PrimitiveCell& w) {
    return std::isfinite(w.rho) && std::isfinite(w.P) &&
        (w.rho > 0.0) && (w.P > 0.0) &&
        std::isfinite(w.u) && std::isfinite(w.v) && std::isfinite(w.w);
}

PrimitiveCell OsherRiemannSolver::ToNormalFrame(const PrimitiveCell& w, const Axis axis) {
    double un = 0.0, ut1 = 0.0, ut2 = 0.0;
    riemann::SplitVelocity(w, axis, un, ut1, ut2);
    // store (un,ut1,ut2) into (u,v,w) so that we can treat as Axis::X everywhere
    return PrimitiveCell{w.rho, un, ut1, ut2, w.P};
}

ConservativeCell OsherRiemannSolver::Diff(const ConservativeCell& a, const ConservativeCell& b) {
    return a - b;
}

bool OsherRiemannSolver::Solve3x3(const double r1_0, const double r1_1, const double r1_2,
                                  const double r2_0, const double r2_1, const double r2_2,
                                  const double r3_0, const double r3_1, const double r3_2,
                                  const double b0, const double b1, const double b2,
                                  double& x1, double& x2, double& x3) {
    double a11 = r1_0, a12 = r2_0, a13 = r3_0;
    double a21 = r1_1, a22 = r2_1, a23 = r3_1;
    double a31 = r1_2, a32 = r2_2, a33 = r3_2;

    double bb1 = b0, bb2 = b1, bb3 = b2;

    const double eps = 1e-14;

    // Pivot for column 1
    int pivot_row = 0;
    double max_abs = std::fabs(a11);
    if (std::fabs(a21) > max_abs) {
        max_abs = std::fabs(a21);
        pivot_row = 1;
    }
    if (std::fabs(a31) > max_abs) {
        max_abs = std::fabs(a31);
        pivot_row = 2;
    }
    if (max_abs < eps) return false;

    if (pivot_row == 1) {
        std::swap(a11, a21);
        std::swap(a12, a22);
        std::swap(a13, a23);
        std::swap(bb1, bb2);
    }
    else if (pivot_row == 2) {
        std::swap(a11, a31);
        std::swap(a12, a32);
        std::swap(a13, a33);
        std::swap(bb1, bb3);
    }

    const double m21 = a21 / a11;
    const double m31 = a31 / a11;

    a21 -= m21 * a11;
    a22 -= m21 * a12;
    a23 -= m21 * a13;
    bb2 -= m21 * bb1;
    a31 -= m31 * a11;
    a32 -= m31 * a12;
    a33 -= m31 * a13;
    bb3 -= m31 * bb1;

    // Pivot for column 2
    pivot_row = 1;
    max_abs = std::fabs(a22);
    if (std::fabs(a32) > max_abs) {
        max_abs = std::fabs(a32);
        pivot_row = 2;
    }
    if (max_abs < eps) return false;

    if (pivot_row == 2) {
        std::swap(a22, a32);
        std::swap(a23, a33);
        std::swap(bb2, bb3);
    }

    const double m32 = a32 / a22;
    a32 -= m32 * a22;
    a33 -= m32 * a23;
    bb3 -= m32 * bb2;

    if (std::fabs(a33) < eps) return false;

    x3 = bb3 / a33;
    x2 = (bb2 - a23 * x3) / a22;
    x1 = (bb1 - a12 * x2 - a13 * x3) / a11;

    return std::isfinite(x1) && std::isfinite(x2) && std::isfinite(x3);
}

bool OsherRiemannSolver::ApplyAbsJacobian(const PrimitiveCell& w,
                                          const double gamma,
                                          const ConservativeCell& dU,
                                          ConservativeCell& out) {
    if (!IsPhysical(w)) return false;

    const double rho = w.rho;
    const double un = w.u; // normal velocity in rotated frame
    const double p = w.P;

    const double a2 = gamma * p / rho;
    if (!(a2 > 0.0) || !std::isfinite(a2)) return false;
    const double a = std::sqrt(a2);

    // Total enthalpy H = (E + p)/rho
    const double kinetic = 0.5 * (w.u * w.u + w.v * w.v + w.w * w.w);
    const double E = p / (gamma - 1.0) + rho * kinetic;
    const double H = (E + p) / rho;

    const double lambda1 = un - a;
    const double lambda2 = un;
    const double lambda3 = un + a;

    // Eigenvectors for (rho, rho*un, E) in normal frame
    const double r1_0 = 1.0;
    const double r1_1 = un - a;
    const double r1_2 = H - un * a;

    const double r2_0 = 1.0;
    const double r2_1 = un;
    const double r2_2 = 0.5 * un * un;

    const double r3_0 = 1.0;
    const double r3_1 = un + a;
    const double r3_2 = H + un * a;

    // Solve R * alpha = dU_1D where dU_1D = [d(rho), d(rho*un), dE]
    double alpha1 = 0.0, alpha2 = 0.0, alpha3 = 0.0;
    if (!Solve3x3(r1_0, r1_1, r1_2,
                  r2_0, r2_1, r2_2,
                  r3_0, r3_1, r3_2,
                  dU.rho, dU.rhoU, dU.E,
                  alpha1, alpha2, alpha3)) {
        return false;
    }

    const double s1 = std::fabs(lambda1) * alpha1;
    const double s2 = std::fabs(lambda2) * alpha2;
    const double s3 = std::fabs(lambda3) * alpha3;

    // 1D part in (rho, rho*un, E)
    out.rho = s1 * r1_0 + s2 * r2_0 + s3 * r3_0;
    out.rhoU = s1 * r1_1 + s2 * r2_1 + s3 * r3_1;
    out.E = s1 * r1_2 + s2 * r2_2 + s3 * r3_2;

    // Tangential momenta: passive advection with eigenvalue un
    const double sT = std::fabs(un);
    out.rhoV = sT * dU.rhoV; // corresponds to t1 momentum in rotated frame
    out.rhoW = sT * dU.rhoW; // corresponds to t2 momentum in rotated frame

    return std::isfinite(out.rho) && std::isfinite(out.rhoU) &&
        std::isfinite(out.rhoV) && std::isfinite(out.rhoW) &&
        std::isfinite(out.E);
}

FluxCell OsherRiemannSolver::RusanovFluxRot(const PrimitiveCell& WL,
                                            const PrimitiveCell& WR,
                                            const double gamma) {
    const FluxCell FL = EulerFlux(WL, gamma, Axis::X);
    const FluxCell FR = EulerFlux(WR, gamma, Axis::X);

    const double aL = SoundSpeed(WL, gamma);
    const double aR = SoundSpeed(WR, gamma);

    const double unL = WL.u;
    const double unR = WR.u;

    const double a_max = std::max(std::abs(unL) + aL, std::abs(unR) + aR);

    const ConservativeCell UL = ConservativeFromPrimitive(WL, gamma);
    const ConservativeCell UR = ConservativeFromPrimitive(WR, gamma);
    const ConservativeCell dU = Diff(UR, UL);

    FluxCell F;
    F.mass = 0.5 * (FL.mass + FR.mass) - 0.5 * a_max * dU.rho;
    F.mom_x = 0.5 * (FL.mom_x + FR.mom_x) - 0.5 * a_max * dU.rhoU;
    F.mom_y = 0.5 * (FL.mom_y + FR.mom_y) - 0.5 * a_max * dU.rhoV;
    F.mom_z = 0.5 * (FL.mom_z + FR.mom_z) - 0.5 * a_max * dU.rhoW;
    F.energy = 0.5 * (FL.energy + FR.energy) - 0.5 * a_max * dU.E;
    return F;
}

void OsherRiemannSolver::ComposeMomentumFluxToGlobal(const FluxCell& F_rot,
                                                     const Axis axis,
                                                     FluxCell& F) {
    F.mass = F_rot.mass;
    F.energy = F_rot.energy;
    riemann::ComposeVector(F_rot.mom_x, F_rot.mom_y, F_rot.mom_z,
                           axis, F.mom_x, F.mom_y, F.mom_z);
}

auto OsherRiemannSolver::ComputeFlux(const PrimitiveCell& left,
                                     const PrimitiveCell& right,
                                     const double gamma,
                                     const Axis axis) const -> FluxCell {
    // Rotate to normal frame; store (un,ut1,ut2) in (u,v,w) and treat as Axis::X.
    const PrimitiveCell WL = ToNormalFrame(left, axis);
    const PrimitiveCell WR = ToNormalFrame(right, axis);

    // Fallback for non-physical inputs
    if (!IsPhysical(WL) || !IsPhysical(WR)) {
        const FluxCell F_rot = RusanovFluxRot(WL, WR, gamma);
        FluxCell F;
        ComposeMomentumFluxToGlobal(F_rot, axis, F);
        return F;
    }

    const ConservativeCell UL = ConservativeFromPrimitive(WL, gamma);
    const ConservativeCell UR = ConservativeFromPrimitive(WR, gamma);
    const ConservativeCell dU = Diff(UR, UL);

    const double norm_dU =
        std::abs(dU.rho) + std::abs(dU.rhoU) + std::abs(dU.rhoV) +
        std::abs(dU.rhoW) + std::abs(dU.E);

    if (norm_dU < 1e-14) {
        const FluxCell F_rot = EulerFlux(WL, gamma, Axis::X);
        FluxCell F;
        ComposeMomentumFluxToGlobal(F_rot, axis, F);
        return F;
    }

    const FluxCell FL_rot = EulerFlux(WL, gamma, Axis::X);
    const FluxCell FR_rot = EulerFlux(WR, gamma, Axis::X);

    // 3-point Gauss–Legendre on [0,1]
    const double s[3] = {
        0.5 * (1.0 - std::sqrt(3.0 / 5.0)),
        0.5,
        0.5 * (1.0 + std::sqrt(3.0 / 5.0))
    };
    const double wq[3] = {5.0 / 18.0, 4.0 / 9.0, 5.0 / 18.0};

    ConservativeCell integral{}; // all zeros by default

    for (int m = 0; m < 3; ++m) {
        const double theta = s[m];

        PrimitiveCell Wt;
        Wt.rho = (1.0 - theta) * WL.rho + theta * WR.rho;
        Wt.u = (1.0 - theta) * WL.u + theta * WR.u;
        Wt.v = (1.0 - theta) * WL.v + theta * WR.v;
        Wt.w = (1.0 - theta) * WL.w + theta * WR.w;
        Wt.P = (1.0 - theta) * WL.P + theta * WR.P;

        ConservativeCell absA_dU;
        if (!ApplyAbsJacobian(Wt, gamma, dU, absA_dU)) {
            const FluxCell F_rot = RusanovFluxRot(WL, WR, gamma);
            FluxCell F;
            ComposeMomentumFluxToGlobal(F_rot, axis, F);
            return F;
        }

        integral += (wq[m] * absA_dU);
    }

    // Osher flux in normal frame
    FluxCell F_rot;
    F_rot.mass = 0.5 * (FL_rot.mass + FR_rot.mass) - 0.5 * integral.rho;
    F_rot.mom_x = 0.5 * (FL_rot.mom_x + FR_rot.mom_x) - 0.5 * integral.rhoU; // mom_n
    F_rot.mom_y = 0.5 * (FL_rot.mom_y + FR_rot.mom_y) - 0.5 * integral.rhoV; // mom_t1
    F_rot.mom_z = 0.5 * (FL_rot.mom_z + FR_rot.mom_z) - 0.5 * integral.rhoW; // mom_t2
    F_rot.energy = 0.5 * (FL_rot.energy + FR_rot.energy) - 0.5 * integral.E;

    // Map rotated momentum flux back to global (x,y,z)
    FluxCell F;
    ComposeMomentumFluxToGlobal(F_rot, axis, F);
    return F;
}
