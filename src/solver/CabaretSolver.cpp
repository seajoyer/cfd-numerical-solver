#include "solver/CabaretSolver.hpp"

#include <algorithm>
#include <cmath>
#include <utility>

#include "solver/EOS.hpp"
#include "solver/PositivityLimiter.hpp"
#include "solver/TimeStepCalculator.hpp"

CabaretSolver::CabaretSolver(const Settings& settings)
    : settings_(settings), boundary_manager_(settings.dim) {
    cfl_ = settings_.cfl;
}

auto CabaretSolver::Step(DataLayer& layer, double& t_cur) -> double {
    boundary_manager_.ApplyAll(layer);

    const double dx = ComputeDx(layer);
    const int total_size = layer.GetTotalSize();
    const int core_start = layer.GetCoreStart(0);
    const int core_end = layer.GetCoreEndExclusive(0);
    const int n_core = core_end - core_start;

    if (n_core < 2 || total_size < 3) {
        return 0.0;
    }

    double dt = TimeStepCalculator::ComputeDt(layer, dx, cfl_, settings_.gamma);
    if (dt <= 0.0) {
        return 0.0;
    }

    if (settings_.t_end > 0.0 && t_cur + dt > settings_.t_end) {
        dt = settings_.t_end - t_cur;
        if (dt <= 0.0) {
            return 0.0;
        }
    }

    const double dt_over_dx = dt / dx;
    const int n_interfaces = total_size - 1;

    xt::xarray<Conservative> U =
        xt::xarray<Conservative>::from_shape({static_cast<std::size_t>(total_size)});
    for (int i = 0; i < total_size; ++i) {
        U(i) = EOS::PrimToCons(layer.GetPrimitive(i), settings_.gamma);
    }

    xt::xarray<Flux> F =
        xt::xarray<Flux>::from_shape({static_cast<std::size_t>(total_size)});
    for (int i = 0; i < total_size; ++i) {
        F(i) = EulerFlux(layer.GetPrimitive(i), settings_.gamma);
    }

    EnsurePsiStorage(U, F, dt_over_dx);

    xt::xarray<Conservative> psi_half_new =
        xt::xarray<Conservative>::from_shape({static_cast<std::size_t>(n_interfaces)});

    for (int k = 0; k < n_interfaces; ++k) {
        const Flux dF = Flux::Diff(F(k + 1), F(k));
        const Conservative dF_u(dF.mass, dF.momentum, dF.energy);
        psi_half_new(k) = psi_half_prev_(k) - dt_over_dx * dF_u;
    }

    xt::xarray<Conservative> du =
        xt::xarray<Conservative>::from_shape({static_cast<std::size_t>(total_size)});
    for (int i = 0; i < total_size; ++i) {
        du(i) = Conservative(0.0, 0.0, 0.0);
    }

    for (int k = 0; k < n_interfaces; ++k) {
        const Primitive wl = layer.GetPrimitive(k);
        const Primitive wr = layer.GetPrimitive(k + 1);

        const RoeMatrices M = BuildRoeMatrices(wl, wr, settings_.gamma);

        const CharVec cL = MatVec(M.L, U(k));
        const CharVec cR = MatVec(M.L, U(k + 1));
        const CharVec cH = MatVec(M.L, psi_half_new(k));

        CharVec dc_left = {0.0, 0.0, 0.0};
        CharVec dc_right = {0.0, 0.0, 0.0};

        const double smax = MaxSignalSpeed(M.u, M.a);
        const double sonic_eps = 1e-6 * std::max(1.0, smax);

        for (int q = 0; q < 3; ++q) {
            const double cmin = std::min(cL[q], cR[q]);
            const double cmax = std::max(cL[q], cR[q]);

            const double lambda = M.lambda[q];

            if (lambda > 0.0) {
                if (std::abs(lambda) <= sonic_eps) {
                    const double dc = -lambda * (cR[q] - cL[q]) * dt_over_dx;
                    dc_right[q] = dc;
                } else {
                    const double cand = 2.0 * cH[q] - cL[q];
                    const double cnew = CabaretLimit(cH[q], cL[q], cR[q], cmin, cmax);
                    dc_right[q] = cnew - cR[q];
                }
            } else if (lambda < 0.0) {
                if (std::abs(lambda) <= sonic_eps) {
                    const double dc = -lambda * (cR[q] - cL[q]) * dt_over_dx;
                    dc_left[q] = dc;
                } else {
                    const double cand = 2.0 * cH[q] - cR[q];
                    (void)cand;
                    const double cnew = CabaretLimit(cH[q], cR[q], cL[q], cmin, cmax);
                    dc_left[q] = cnew - cL[q];
                }
            }
        }

        const Conservative dU_left = MatVec(M.R, dc_left);
        const Conservative dU_right = MatVec(M.R, dc_right);

        du(k) += dU_left;
        du(k + 1) += dU_right;
    }

    for (int j = core_start; j < core_end; ++j) {
        Conservative Unew = U(j) + du(j);
        PositivityLimiter::Apply(Unew, settings_.gamma, rho_min_, p_min_);
        StoreConservativeCell(Unew, j, dx, layer);
    }

    psi_half_prev_ = std::move(psi_half_new);
    t_cur += dt;
    return dt;
}

void CabaretSolver::SetCfl(double cfl) {
    cfl_ = cfl;
}

void CabaretSolver::AddBoundary(int axis,
                                std::shared_ptr<BoundaryCondition> left_bc,
                                std::shared_ptr<BoundaryCondition> right_bc) {
    boundary_manager_.Set(axis, std::move(left_bc), std::move(right_bc));
}

auto CabaretSolver::ComputeDx(const DataLayer& layer) const -> double {
    if (settings_.N > 0 && settings_.L_x > 0.0) {
        return settings_.L_x / static_cast<double>(settings_.N);
    }

    const int core_start = layer.GetCoreStart(0);
    const int core_end = layer.GetCoreEndExclusive(0);
    if (core_end - core_start > 1) {
        return layer.xc(core_start + 1) - layer.xc(core_start);
    }

    return 1.0;
}

void CabaretSolver::StoreConservativeCell(const Conservative& uc,
                                          const int i,
                                          const double dx,
                                          DataLayer& layer) const {
    const double rho = uc.rho;
    const double rhoU = uc.rhoU;

    const double uvel = rho > 0.0 ? rhoU / rho : 0.0;
    const double P = (rho > 0.0) ? EOS::Pressure(uc, settings_.gamma) : 0.0;

    layer.rho(i) = rho;
    layer.u(i) = uvel;
    layer.P(i) = P;

    layer.p(i) = rhoU;
    layer.V(i) = rho > 0.0 ? 1.0 / rho : 0.0;

    const double kinetic = 0.5 * rho * uvel * uvel;
    const double Eint = uc.E - kinetic;
    const double eint = rho > 0.0 ? Eint / rho : 0.0;
    const double Etot = rho > 0.0 ? uc.E : 0.0;

    layer.U(i) = eint;
    layer.e(i) = Etot;
    layer.m(i) = rho * dx;
}

void CabaretSolver::EnsurePsiStorage(const xt::xarray<Conservative>& U,
                                     const xt::xarray<Flux>& F,
                                     const double dt_over_dx) {
    const int total_size = static_cast<int>(U.size());
    const int n_interfaces = total_size - 1;

    if (n_interfaces <= 0) {
        psi_initialized_ = true;
        return;
    }

    if (!psi_initialized_ || static_cast<int>(psi_half_prev_.size()) != n_interfaces) {
        psi_half_prev_ =
            xt::xarray<Conservative>::from_shape({static_cast<std::size_t>(n_interfaces)});

        for (int k = 0; k < n_interfaces; ++k) {
            const Flux dF = Flux::Diff(F(k + 1), F(k));
            const Conservative dF_u(dF.mass, dF.momentum, dF.energy);

            const Conservative avg = 0.5 * (U(k) + U(k + 1));

            psi_half_prev_(k) = avg + dt_over_dx * dF_u;
        }

        psi_initialized_ = true;
    }
}

auto CabaretSolver::BuildRoeMatrices(const Primitive& left,
                                     const Primitive& right,
                                     const double gamma) -> RoeMatrices {
    RoeMatrices M{};

    const double rho_l = std::max(left.rho, 1e-14);
    const double rho_r = std::max(right.rho, 1e-14);

    const double u_l = left.u;
    const double u_r = right.u;

    const double p_l = std::max(left.P, 1e-14);
    const double p_r = std::max(right.P, 1e-14);

    const Conservative Ul = EOS::PrimToCons({rho_l, u_l, p_l}, gamma);
    const Conservative Ur = EOS::PrimToCons({rho_r, u_r, p_r}, gamma);

    const double H_l = (Ul.E + p_l) / rho_l;
    const double H_r = (Ur.E + p_r) / rho_r;

    const double sr_l = std::sqrt(rho_l);
    const double sr_r = std::sqrt(rho_r);
    const double denom = sr_l + sr_r;

    const double u = (sr_l * u_l + sr_r * u_r) / denom;
    const double H = (sr_l * H_l + sr_r * H_r) / denom;

    const double a2 = (gamma - 1.0) * (H - 0.5 * u * u);
    const double a = std::sqrt(std::max(a2, 1e-14));

    M.u = u;
    M.a = a;

    M.lambda[0] = u - a;
    M.lambda[1] = u;
    M.lambda[2] = u + a;

    M.R[0][0] = 1.0;
    M.R[1][0] = u - a;
    M.R[2][0] = H - u * a;

    M.R[0][1] = 1.0;
    M.R[1][1] = u;
    M.R[2][1] = 0.5 * u * u;

    M.R[0][2] = 1.0;
    M.R[1][2] = u + a;
    M.R[2][2] = H + u * a;

    const double beta = (gamma - 1.0) / (a * a);

    M.L[0][0] = 0.5 * (beta * u * u + u / a);
    M.L[0][1] = -0.5 * (beta * u + 1.0 / a);
    M.L[0][2] = 0.5 * beta;

    M.L[1][0] = 1.0 - beta * u * u;
    M.L[1][1] = beta * u;
    M.L[1][2] = -beta;

    M.L[2][0] = 0.5 * (beta * u * u - u / a);
    M.L[2][1] = -0.5 * (beta * u - 1.0 / a);
    M.L[2][2] = 0.5 * beta;

    return M;
}

auto CabaretSolver::MatVec(const double A[3][3], const Conservative& u) -> CharVec {
    const double x0 = u.rho;
    const double x1 = u.rhoU;
    const double x2 = u.E;

    CharVec w{};
    w[0] = A[0][0] * x0 + A[0][1] * x1 + A[0][2] * x2;
    w[1] = A[1][0] * x0 + A[1][1] * x1 + A[1][2] * x2;
    w[2] = A[2][0] * x0 + A[2][1] * x1 + A[2][2] * x2;
    return w;
}

auto CabaretSolver::MatVec(const double A[3][3], const CharVec& w) -> Conservative {
    Conservative u;
    u.rho = A[0][0] * w[0] + A[0][1] * w[1] + A[0][2] * w[2];
    u.rhoU = A[1][0] * w[0] + A[1][1] * w[1] + A[1][2] * w[2];
    u.E = A[2][0] * w[0] + A[2][1] * w[1] + A[2][2] * w[2];
    return u;
}

auto CabaretSolver::CabaretLimit(const double w_half,
                                 const double w_upwind,
                                 const double w_downwind,
                                 const double w_min,
                                 const double w_max) -> double {
    const double cand = 2.0 * w_half - w_upwind;

    const double d = cand - w_half;
    if (std::abs(d) <= 1e-30) {
        return w_half;
    }

    double theta = 1.0;

    if (d > 0.0) {
        theta = (w_max - w_half) / d;
    } else {
        theta = (w_min - w_half) / d;
    }

    theta = std::min(1.0, std::max(0.0, theta));

    const double w_new = w_half + theta * d;

    (void)w_downwind;
    return w_new;
}

auto CabaretSolver::MaxSignalSpeed(const double u, const double a) -> double {
    return std::abs(u) + a;
}
