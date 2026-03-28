#include "riemann/ExactIdealGasRiemannSolver.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

ExactIdealGasRiemannSolver::ExactIdealGasRiemannSolver() = default;

ExactIdealGasRiemannSolver::ExactIdealGasRiemannSolver(const double xi,
                                                       const double Q_user)
    : xi_(xi),
      Q_user_(Q_user) {}

void ExactIdealGasRiemannSolver::SetXi(const double xi) {
    xi_ = xi;
}

void ExactIdealGasRiemannSolver::SetQ(const double Q) {
    Q_user_ = Q;
}

void ExactIdealGasRiemannSolver::BuildTangentialBasis(const FaceNormal& normal,
                                                      double& t1_x,
                                                      double& t1_y,
                                                      double& t1_z,
                                                      double& t2_x,
                                                      double& t2_y,
                                                      double& t2_z) const {
    double ref_x = 0.0;
    double ref_y = 0.0;
    double ref_z = 0.0;

    if (std::abs(normal.x) < 0.9) {
        ref_x = 1.0;
    }
    else {
        ref_y = 1.0;
    }

    const double dot = ref_x * normal.x + ref_y * normal.y + ref_z * normal.z;

    t1_x = ref_x - dot * normal.x;
    t1_y = ref_y - dot * normal.y;
    t1_z = ref_z - dot * normal.z;

    const double t1_norm = std::sqrt(t1_x * t1_x + t1_y * t1_y + t1_z * t1_z);
    if (t1_norm <= 1e-14) {
        throw std::runtime_error(
            "ExactIdealGasRiemannSolver::BuildTangentialBasis: failed to build first tangential vector"
        );
    }

    t1_x /= t1_norm;
    t1_y /= t1_norm;
    t1_z /= t1_norm;

    t2_x = normal.y * t1_z - normal.z * t1_y;
    t2_y = normal.z * t1_x - normal.x * t1_z;
    t2_z = normal.x * t1_y - normal.y * t1_x;

    const double t2_norm = std::sqrt(t2_x * t2_x + t2_y * t2_y + t2_z * t2_z);
    if (t2_norm <= 1e-14) {
        throw std::runtime_error(
            "ExactIdealGasRiemannSolver::BuildTangentialBasis: failed to build second tangential vector"
        );
    }

    t2_x /= t2_norm;
    t2_y /= t2_norm;
    t2_z /= t2_norm;
}

void ExactIdealGasRiemannSolver::ProjectVelocityToLocalBasis(const PrimitiveCell& state,
                                                             const FaceNormal& normal,
                                                             const double t1_x,
                                                             const double t1_y,
                                                             const double t1_z,
                                                             const double t2_x,
                                                             const double t2_y,
                                                             const double t2_z,
                                                             double& u_n,
                                                             double& u_t1,
                                                             double& u_t2) const {
    u_n = state.u * normal.x + state.v * normal.y + state.w * normal.z;
    u_t1 = state.u * t1_x + state.v * t1_y + state.w * t1_z;
    u_t2 = state.u * t2_x + state.v * t2_y + state.w * t2_z;
}

void ExactIdealGasRiemannSolver::ComposeVelocityFromLocalBasis(const double u_n,
                                                               const double u_t1,
                                                               const double u_t2,
                                                               const FaceNormal& normal,
                                                               const double t1_x,
                                                               const double t1_y,
                                                               const double t1_z,
                                                               const double t2_x,
                                                               const double t2_y,
                                                               const double t2_z,
                                                               double& u,
                                                               double& v,
                                                               double& w) const {
    u = u_n * normal.x + u_t1 * t1_x + u_t2 * t2_x;
    v = u_n * normal.y + u_t1 * t1_y + u_t2 * t2_y;
    w = u_n * normal.z + u_t1 * t1_z + u_t2 * t2_z;
}

ExactIdealGasRiemannSolver::State1D ExactIdealGasRiemannSolver::MakeState1D(
    const PrimitiveCell& state,
    const FaceNormal& normal,
    const double t1_x,
    const double t1_y,
    const double t1_z,
    const double t2_x,
    const double t2_y,
    const double t2_z,
    const double gamma
) const {
    double u_n = 0.0;
    double u_t1 = 0.0;
    double u_t2 = 0.0;

    ProjectVelocityToLocalBasis(state,
                                normal,
                                t1_x, t1_y, t1_z,
                                t2_x, t2_y, t2_z,
                                u_n, u_t1, u_t2);

    State1D s;
    s.rho = state.rho;
    s.un = u_n;
    s.p = state.P;
    s.a = std::sqrt(std::max(gamma * s.p / s.rho, 0.0));

    return s;
}

double ExactIdealGasRiemannSolver::PhiRarefaction(const double p,
                                                  const State1D& state,
                                                  const double gamma) const {
    const double pressure_ratio = p / state.p;
    const double exponent = (gamma - 1.0) / (2.0 * gamma);
    return 2.0 * state.a / (gamma - 1.0) * (std::pow(pressure_ratio, exponent) - 1.0);
}

double ExactIdealGasRiemannSolver::PhiRarefactionDerivative(const double p,
                                                            const State1D& state,
                                                            const double gamma) const {
    const double pressure_ratio = p / state.p;
    const double exponent = (gamma - 1.0) / (2.0 * gamma);
    const double coefficient = (2.0 * state.a / (gamma - 1.0)) * (exponent / state.p);
    return coefficient * std::pow(pressure_ratio, exponent - 1.0);
}

double ExactIdealGasRiemannSolver::PhiShock(const double p,
                                            const State1D& state,
                                            const double gamma) const {
    const double A = 2.0 / ((gamma + 1.0) * state.rho);
    const double B = (gamma - 1.0) / (gamma + 1.0) * state.p;
    return (p - state.p) * std::sqrt(A / (p + B));
}

double ExactIdealGasRiemannSolver::PhiShockDerivative(const double p,
                                                      const State1D& state,
                                                      const double gamma) const {
    const double A = 2.0 / ((gamma + 1.0) * state.rho);
    const double B = (gamma - 1.0) / (gamma + 1.0) * state.p;
    const double sqrt_term = std::sqrt(A / (p + B));
    const double term = (p - state.p) / (2.0 * (p + B));
    return sqrt_term * (1.0 - term);
}

double ExactIdealGasRiemannSolver::SidePhi(const double p,
                                           const State1D& state,
                                           const double gamma) const {
    return (p <= state.p)
               ? PhiRarefaction(p, state, gamma)
               : PhiShock(p, state, gamma);
}

double ExactIdealGasRiemannSolver::SidePhiDerivative(const double p,
                                                     const State1D& state,
                                                     const double gamma) const {
    return (p <= state.p)
               ? PhiRarefactionDerivative(p, state, gamma)
               : PhiShockDerivative(p, state, gamma);
}

double ExactIdealGasRiemannSolver::InitialGuess(const State1D& left,
                                                const State1D& right,
                                                const double Q_user,
                                                const double gamma) const {
    const double p_pvrs =
        0.5 * (left.p + right.p) -
        0.125 * (right.un - left.un) * (left.rho + right.rho) * (left.a + right.a);

    const double p_min = std::min(left.p, right.p);
    const double p_max = std::max(left.p, right.p);
    const double Q = p_max / p_min;

    double p_star = std::max(p_pvrs, 1e-16);

    if (p_min < p_pvrs && p_pvrs < p_max && Q < Q_user) {
        return p_star;
    }

    if (p_pvrs <= p_min) {
        const double exponent = (gamma - 1.0) / (2.0 * gamma);
        const double numerator =
            left.a + right.a - 0.5 * (gamma - 1.0) * (right.un - left.un);
        const double denominator =
            left.a / std::pow(left.p, exponent) +
            right.a / std::pow(right.p, exponent);

        p_star = std::pow(numerator / denominator, 1.0 / exponent);
        return std::max(p_star, 1e-16);
    }

    const auto A = [gamma](const State1D& state) -> double {
        return 2.0 / ((gamma + 1.0) * state.rho);
    };

    const auto B = [gamma](const State1D& state) -> double {
        return (gamma - 1.0) / (gamma + 1.0) * state.p;
    };

    const auto g = [&](const State1D& state) -> double {
        return std::sqrt(A(state) / (p_pvrs + B(state)));
    };

    const double g_left = g(left);
    const double g_right = g(right);

    p_star =
        (g_left * left.p + g_right * right.p - (right.un - left.un)) /
        (g_left + g_right);

    return std::max(p_star, 1e-16);
}

double ExactIdealGasRiemannSolver::SolveStarPressure(const State1D& left,
                                                     const State1D& right,
                                                     const double gamma,
                                                     const double Q_user) const {
    const double delta_u = right.un - left.un;
    const double critical = 2.0 * (left.a + right.a) / (gamma - 1.0);

    if (delta_u >= critical) {
        return 0.0;
    }

    double p = InitialGuess(left, right, Q_user, gamma);
    const double p_min = 1e-16;

    for (int iter = 0; iter < 40; ++iter) {
        const double f_left = SidePhi(p, left, gamma);
        const double f_right = SidePhi(p, right, gamma);
        const double f = f_left + f_right + (right.un - left.un);

        if (std::fabs(f) < 1e-10) {
            break;
        }

        const double df_left = SidePhiDerivative(p, left, gamma);
        const double df_right = SidePhiDerivative(p, right, gamma);
        const double df = df_left + df_right;

        if (df == 0.0 || !std::isfinite(df)) {
            break;
        }

        double p_new = p - f / df;

        if (!std::isfinite(p_new) || p_new < p_min) {
            p_new = p_min;
        }

        if (std::fabs(p_new - p) <= 1e-8 * (p + p_min)) {
            p = p_new;
            break;
        }

        p = p_new;
    }

    if (!std::isfinite(p) || p <= 0.0) {
        return 0.0;
    }

    return p;
}

ExactIdealGasRiemannSolver::Primitive1D ExactIdealGasRiemannSolver::SampleVacuum(
    const double xi,
    const State1D& left,
    const State1D& right,
    const double gamma
) const {
    const double shl = left.un - left.a;
    const double svl = left.un + 2.0 * left.a / (gamma - 1.0);

    const double shr = right.un + right.a;
    const double svr = right.un - 2.0 * right.a / (gamma - 1.0);

    if (xi <= shl) {
        return Primitive1D{left.rho, left.un, left.p};
    }

    if (xi >= shr) {
        return Primitive1D{right.rho, right.un, right.p};
    }

    if (xi > shl && xi < svl) {
        const double u_n =
            2.0 / (gamma + 1.0) *
            (left.a + 0.5 * (gamma - 1.0) * left.un + xi);

        const double a =
            2.0 / (gamma + 1.0) *
            (left.a + 0.5 * (gamma - 1.0) * (left.un - xi));

        const double rho = left.rho * std::pow(a / left.a, 2.0 / (gamma - 1.0));
        const double p = left.p * std::pow(a / left.a, 2.0 * gamma / (gamma - 1.0));

        return Primitive1D{rho, u_n, p};
    }

    if (xi > svr && xi < shr) {
        const double u_n =
            2.0 / (gamma + 1.0) *
            (-right.a + 0.5 * (gamma - 1.0) * right.un + xi);

        const double a =
            2.0 / (gamma + 1.0) *
            (right.a - 0.5 * (gamma - 1.0) * right.un + 0.5 * (gamma - 1.0) * xi);

        const double rho = right.rho * std::pow(a / right.a, 2.0 / (gamma - 1.0));
        const double p = right.p * std::pow(a / right.a, 2.0 * gamma / (gamma - 1.0));

        return Primitive1D{rho, u_n, p};
    }

    return Primitive1D{0.0, 0.0, 0.0};
}

ExactIdealGasRiemannSolver::Primitive1D ExactIdealGasRiemannSolver::SampleNonVacuum(
    const double xi,
    const double p_star,
    const double u_star,
    const State1D& left,
    const State1D& right,
    const double gamma
) const {
    if (xi <= u_star) {
        if (p_star > left.p) {
            const double q = p_star / left.p;
            const double shock_speed =
                left.un -
                left.a * std::sqrt(
                    0.5 * ((gamma + 1.0) / gamma * q + (gamma - 1.0) / gamma)
                );

            if (xi <= shock_speed) {
                return Primitive1D{left.rho, left.un, left.p};
            }

            const double factor =
                (q + (gamma - 1.0) / (gamma + 1.0)) /
                ((gamma - 1.0) / (gamma + 1.0) * q + 1.0);

            const double rho_star_left = left.rho * factor;
            return Primitive1D{rho_star_left, u_star, p_star};
        }

        const double shl = left.un - left.a;
        const double rho_star_left = left.rho * std::pow(p_star / left.p, 1.0 / gamma);
        const double a_star_left = std::sqrt(std::max(gamma * p_star / rho_star_left, 0.0));
        const double stl = u_star - a_star_left;

        if (xi <= shl) {
            return Primitive1D{left.rho, left.un, left.p};
        }

        if (xi >= stl) {
            return Primitive1D{rho_star_left, u_star, p_star};
        }

        const double u_n =
            2.0 / (gamma + 1.0) *
            (left.a + 0.5 * (gamma - 1.0) * left.un + xi);

        const double a =
            2.0 / (gamma + 1.0) *
            (left.a + 0.5 * (gamma - 1.0) * (left.un - xi));

        const double rho = left.rho * std::pow(a / left.a, 2.0 / (gamma - 1.0));
        const double p = left.p * std::pow(a / left.a, 2.0 * gamma / (gamma - 1.0));

        return Primitive1D{rho, u_n, p};
    }

    if (p_star > right.p) {
        const double q = p_star / right.p;
        const double shock_speed =
            right.un +
            right.a * std::sqrt(
                0.5 * ((gamma + 1.0) / gamma * q + (gamma - 1.0) / gamma)
            );

        if (xi >= shock_speed) {
            return Primitive1D{right.rho, right.un, right.p};
        }

        const double factor =
            (q + (gamma - 1.0) / (gamma + 1.0)) /
            ((gamma - 1.0) / (gamma + 1.0) * q + 1.0);

        const double rho_star_right = right.rho * factor;
        return Primitive1D{rho_star_right, u_star, p_star};
    }

    const double shr = right.un + right.a;

    const double rho_star_right = right.rho * std::pow(p_star / right.p, 1.0 / gamma);
    const double a_star_right = std::sqrt(std::max(gamma * p_star / rho_star_right, 0.0));
    const double str = u_star + a_star_right;

    if (xi >= shr) {
        return Primitive1D{right.rho, right.un, right.p};
    }

    if (xi <= str) {
        return Primitive1D{rho_star_right, u_star, p_star};
    }

    const double u_n =
        2.0 / (gamma + 1.0) *
        (-right.a + 0.5 * (gamma - 1.0) * right.un + xi);

    const double a =
        2.0 / (gamma + 1.0) *
        (right.a - 0.5 * (gamma - 1.0) * right.un + 0.5 * (gamma - 1.0) * xi);

    const double rho = right.rho * std::pow(a / right.a, 2.0 / (gamma - 1.0));
    const double p = right.p * std::pow(a / right.a, 2.0 * gamma / (gamma - 1.0));

    return Primitive1D{rho, u_n, p};
}

PrimitiveCell ExactIdealGasRiemannSolver::Sample(const PrimitiveCell& left,
                                                 const PrimitiveCell& right,
                                                 const double gamma,
                                                 const double xi,
                                                 const FaceNormal& normal) const {
    if (!IsUnitNormal(normal)) {
        throw std::runtime_error(
            "ExactIdealGasRiemannSolver::Sample: face normal must be unit-length"
        );
    }

    double t1_x = 0.0;
    double t1_y = 0.0;
    double t1_z = 0.0;
    double t2_x = 0.0;
    double t2_y = 0.0;
    double t2_z = 0.0;
    BuildTangentialBasis(normal, t1_x, t1_y, t1_z, t2_x, t2_y, t2_z);

    const State1D left_state =
        MakeState1D(left, normal, t1_x, t1_y, t1_z, t2_x, t2_y, t2_z, gamma);

    const State1D right_state =
        MakeState1D(right, normal, t1_x, t1_y, t1_z, t2_x, t2_y, t2_z, gamma);

    const double p_star = SolveStarPressure(left_state, right_state, gamma, Q_user_);

    Primitive1D sample_1d{};
    double u_star = 0.0;

    if (p_star <= 0.0) {
        sample_1d = SampleVacuum(xi, left_state, right_state, gamma);
        u_star = 0.0;
    }
    else {
        const double f_left = SidePhi(p_star, left_state, gamma);
        const double f_right = SidePhi(p_star, right_state, gamma);
        u_star = 0.5 * (left_state.un + right_state.un + f_right - f_left);
        sample_1d = SampleNonVacuum(xi, p_star, u_star, left_state, right_state, gamma);
    }

    double un_left = 0.0;
    double ut1_left = 0.0;
    double ut2_left = 0.0;
    ProjectVelocityToLocalBasis(left,
                                normal,
                                t1_x, t1_y, t1_z,
                                t2_x, t2_y, t2_z,
                                un_left, ut1_left, ut2_left);

    double un_right = 0.0;
    double ut1_right = 0.0;
    double ut2_right = 0.0;
    ProjectVelocityToLocalBasis(right,
                                normal,
                                t1_x, t1_y, t1_z,
                                t2_x, t2_y, t2_z,
                                un_right, ut1_right, ut2_right);

    const bool use_left_tangential = (xi <= u_star);

    const double ut1 = use_left_tangential ? ut1_left : ut1_right;
    const double ut2 = use_left_tangential ? ut2_left : ut2_right;

    PrimitiveCell sample;
    sample.rho = sample_1d.rho;
    sample.P = sample_1d.p;

    ComposeVelocityFromLocalBasis(sample_1d.un,
                                  ut1,
                                  ut2,
                                  normal,
                                  t1_x, t1_y, t1_z,
                                  t2_x, t2_y, t2_z,
                                  sample.u, sample.v, sample.w);

    return sample;
}

ConservativeCell ExactIdealGasRiemannSolver::ComputeFlux(const PrimitiveCell& left,
                                                         const PrimitiveCell& right,
                                                         const double gamma,
                                                         const FaceNormal& normal) const {
    const PrimitiveCell sample = Sample(left, right, gamma, xi_, normal);
    return PhysicalFlux(sample, gamma, normal);
}
