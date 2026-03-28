#include "riemann/RoeRiemannSolver.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

double RoeRiemannSolver::EntropyFix(const double lambda,
                                    const double lambda_left,
                                    const double lambda_right) const {
    const double delta = std::max(0.0, lambda_right - lambda_left);
    const double abs_lambda = std::abs(lambda);

    if (delta <= 0.0) {
        return abs_lambda;
    }

    if (abs_lambda >= delta) {
        return abs_lambda;
    }

    return 0.5 * (lambda * lambda / delta + delta);
}

void RoeRiemannSolver::BuildTangentialBasis(const FaceNormal& normal,
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

    const double t1_norm =
        std::sqrt(t1_x * t1_x + t1_y * t1_y + t1_z * t1_z);

    if (t1_norm <= 1e-14) {
        throw std::runtime_error(
            "RoeRiemannSolver::BuildTangentialBasis: failed to build first tangential vector"
        );
    }

    t1_x /= t1_norm;
    t1_y /= t1_norm;
    t1_z /= t1_norm;

    t2_x = normal.y * t1_z - normal.z * t1_y;
    t2_y = normal.z * t1_x - normal.x * t1_z;
    t2_z = normal.x * t1_y - normal.y * t1_x;

    const double t2_norm =
        std::sqrt(t2_x * t2_x + t2_y * t2_y + t2_z * t2_z);

    if (t2_norm <= 1e-14) {
        throw std::runtime_error(
            "RoeRiemannSolver::BuildTangentialBasis: failed to build second tangential vector"
        );
    }

    t2_x /= t2_norm;
    t2_y /= t2_norm;
    t2_z /= t2_norm;
}

void RoeRiemannSolver::ProjectVelocityToLocalBasis(const PrimitiveCell& state,
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

void RoeRiemannSolver::ComposeVectorFromLocalBasis(const double v_n,
                                                   const double v_t1,
                                                   const double v_t2,
                                                   const FaceNormal& normal,
                                                   const double t1_x,
                                                   const double t1_y,
                                                   const double t1_z,
                                                   const double t2_x,
                                                   const double t2_y,
                                                   const double t2_z,
                                                   double& v_x,
                                                   double& v_y,
                                                   double& v_z) const {
    v_x = v_n * normal.x + v_t1 * t1_x + v_t2 * t2_x;
    v_y = v_n * normal.y + v_t1 * t1_y + v_t2 * t2_y;
    v_z = v_n * normal.z + v_t1 * t1_z + v_t2 * t2_z;
}

ConservativeCell RoeRiemannSolver::ComputeFlux(const PrimitiveCell& left,
                                               const PrimitiveCell& right,
                                               const double gamma,
                                               const FaceNormal& normal) const {
    if (!IsUnitNormal(normal)) {
        throw std::runtime_error(
            "RoeRiemannSolver::ComputeFlux: face normal must be unit-length"
        );
    }

    double t1_x = 0.0;
    double t1_y = 0.0;
    double t1_z = 0.0;
    double t2_x = 0.0;
    double t2_y = 0.0;
    double t2_z = 0.0;
    BuildTangentialBasis(normal, t1_x, t1_y, t1_z, t2_x, t2_y, t2_z);

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

    const PrimitiveCell left_local{
        left.rho,
        un_left,
        ut1_left,
        ut2_left,
        left.P
    };

    const PrimitiveCell right_local{
        right.rho,
        un_right,
        ut1_right,
        ut2_right,
        right.P
    };

    const FaceNormal x_normal{1.0, 0.0, 0.0};

    const ConservativeCell flux_left_local = PhysicalFlux(left_local, gamma, x_normal);
    const ConservativeCell flux_right_local = PhysicalFlux(right_local, gamma, x_normal);

    const ConservativeCell U_left_local = ConservativeFromPrimitive(left_local, gamma);
    const ConservativeCell U_right_local = ConservativeFromPrimitive(right_local, gamma);

    const double rho_floor = 1e-14;
    const double rho_left = std::max(U_left_local.rho, rho_floor);
    const double rho_right = std::max(U_right_local.rho, rho_floor);

    const double p_left = left_local.P;
    const double p_right = right_local.P;

    const double sqrt_rho_left = std::sqrt(rho_left);
    const double sqrt_rho_right = std::sqrt(rho_right);
    const double denom = sqrt_rho_left + sqrt_rho_right;

    const double un_tilde =
        (sqrt_rho_left * un_left + sqrt_rho_right * un_right) / denom;
    const double ut1_tilde =
        (sqrt_rho_left * ut1_left + sqrt_rho_right * ut1_right) / denom;
    const double ut2_tilde =
        (sqrt_rho_left * ut2_left + sqrt_rho_right * ut2_right) / denom;

    const double H_left = (U_left_local.E + p_left) / rho_left;
    const double H_right = (U_right_local.E + p_right) / rho_right;
    const double H_tilde =
        (sqrt_rho_left * H_left + sqrt_rho_right * H_right) / denom;

    const double q2 =
        0.5 * (un_tilde * un_tilde +
            ut1_tilde * ut1_tilde +
            ut2_tilde * ut2_tilde);

    double a2_tilde = (gamma - 1.0) * (H_tilde - q2);
    a2_tilde = std::max(a2_tilde, 1e-14);

    const double a_tilde = std::sqrt(a2_tilde);
    const double rho_tilde = sqrt_rho_left * sqrt_rho_right;

    const double drho = rho_right - rho_left;
    const double dun = un_right - un_left;
    const double dut1 = ut1_right - ut1_left;
    const double dut2 = ut2_right - ut2_left;
    const double dp = p_right - p_left;

    const double alpha2 = drho - dp / a2_tilde;
    const double alpha1 = 0.5 / a2_tilde * (dp - rho_tilde * a_tilde * dun);
    const double alpha3 = 0.5 / a2_tilde * (dp + rho_tilde * a_tilde * dun);
    const double alpha4 = rho_tilde * dut1;
    const double alpha5 = rho_tilde * dut2;

    const double lambda1 = un_tilde - a_tilde;
    const double lambda2 = un_tilde;
    const double lambda3 = un_tilde + a_tilde;
    const double lambda4 = un_tilde;
    const double lambda5 = un_tilde;

    const double a_left = SoundSpeed(left_local, gamma);
    const double a_right = SoundSpeed(right_local, gamma);

    const double lambda1_left = un_left - a_left;
    const double lambda1_right = un_right - a_right;
    const double lambda3_left = un_left + a_left;
    const double lambda3_right = un_right + a_right;

    const double abs_lambda1 = EntropyFix(lambda1, lambda1_left, lambda1_right);
    const double abs_lambda2 = std::abs(lambda2);
    const double abs_lambda3 = EntropyFix(lambda3, lambda3_left, lambda3_right);
    const double abs_lambda4 = std::abs(lambda4);
    const double abs_lambda5 = std::abs(lambda5);

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

    const double dF0 =
        abs_lambda1 * alpha1 * r1_0 +
        abs_lambda2 * alpha2 * r2_0 +
        abs_lambda3 * alpha3 * r3_0 +
        abs_lambda4 * alpha4 * r4_0 +
        abs_lambda5 * alpha5 * r5_0;

    const double dF1 =
        abs_lambda1 * alpha1 * r1_1 +
        abs_lambda2 * alpha2 * r2_1 +
        abs_lambda3 * alpha3 * r3_1 +
        abs_lambda4 * alpha4 * r4_1 +
        abs_lambda5 * alpha5 * r5_1;

    const double dF2 =
        abs_lambda1 * alpha1 * r1_2 +
        abs_lambda2 * alpha2 * r2_2 +
        abs_lambda3 * alpha3 * r3_2 +
        abs_lambda4 * alpha4 * r4_2 +
        abs_lambda5 * alpha5 * r5_2;

    const double dF3 =
        abs_lambda1 * alpha1 * r1_3 +
        abs_lambda2 * alpha2 * r2_3 +
        abs_lambda3 * alpha3 * r3_3 +
        abs_lambda4 * alpha4 * r4_3 +
        abs_lambda5 * alpha5 * r5_3;

    const double dF4 =
        abs_lambda1 * alpha1 * r1_4 +
        abs_lambda2 * alpha2 * r2_4 +
        abs_lambda3 * alpha3 * r3_4 +
        abs_lambda4 * alpha4 * r4_4 +
        abs_lambda5 * alpha5 * r5_4;

    const double flux_mass_local =
        0.5 * (flux_left_local.rho + flux_right_local.rho) - 0.5 * dF0;
    const double flux_mom_n_local =
        0.5 * (flux_left_local.rhoU + flux_right_local.rhoU) - 0.5 * dF1;
    const double flux_mom_t1_local =
        0.5 * (flux_left_local.rhoV + flux_right_local.rhoV) - 0.5 * dF2;
    const double flux_mom_t2_local =
        0.5 * (flux_left_local.rhoW + flux_right_local.rhoW) - 0.5 * dF3;
    const double flux_energy_local =
        0.5 * (flux_left_local.E + flux_right_local.E) - 0.5 * dF4;

    ConservativeCell flux;
    flux.rho = flux_mass_local;
    flux.E = flux_energy_local;

    ComposeVectorFromLocalBasis(flux_mom_n_local,
                                flux_mom_t1_local,
                                flux_mom_t2_local,
                                normal,
                                t1_x, t1_y, t1_z,
                                t2_x, t2_y, t2_z,
                                flux.rhoU, flux.rhoV, flux.rhoW);

    return flux;
}
