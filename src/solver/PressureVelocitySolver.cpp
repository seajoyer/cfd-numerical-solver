#include "solver/PressureVelocitySolver.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>

PressureVelocitySolver::PressureVelocitySolver(const Settings& settings,
                                               Mesh mesh,
                                               std::shared_ptr<BoundaryManager> boundary_manager,
                                               const MPIContext* mpi_context,
                                               const bool steady)
    : settings_(settings),
      mesh_(std::move(mesh)),
      boundary_manager_(std::move(boundary_manager)),
      mpi_context_(mpi_context),
      steady_(steady) {
    cfl_ = settings_.cfl;

    if (!boundary_manager_) {
        throw std::runtime_error("PressureVelocitySolver: boundary_manager is null");
    }

    internal_boundary_condition_ = std::make_unique<WallInternalBoundary>();
    EnsureStorageSized();
    state_.ZeroAll();
    workspace_.ZeroAll();
    ApplyHaloAndPhysicalBc();
    ApplyImmersedVelocityConstraints();
}

const Mesh& PressureVelocitySolver::GetMesh() const {
    return mesh_;
}

Mesh& PressureVelocitySolver::GetMesh() {
    return mesh_;
}

const PressureVelocityState& PressureVelocitySolver::GetState() const {
    return state_;
}

PressureVelocityState& PressureVelocitySolver::GetState() {
    return state_;
}

void PressureVelocitySolver::SetCfl(const double cfl) {
    cfl_ = cfl;
}

void PressureVelocitySolver::EnsureStorageSized() {
    state_.ResizeFrom(mesh_);
    workspace_.ResizeFrom(mesh_);
}

void PressureVelocitySolver::ApplyHaloAndPhysicalBc() {
    if (!boundary_manager_) {
        throw std::runtime_error("PressureVelocitySolver: boundary_manager is null");
    }

    boundary_manager_->UpdateHalo(state_, mesh_);
    boundary_manager_->ApplyPhysicalBc(state_, mesh_);
}

void PressureVelocitySolver::ApplyImmersedVelocityConstraints() {
    if (!settings_.immersed_enabled || !internal_boundary_condition_) {
        return;
    }

    internal_boundary_condition_->ApplyPressureVelocityVelocityConstraints(state_, workspace_, mesh_);
}

void PressureVelocitySolver::ApplyImmersedMomentumCorrections(const double dt, const double nu) {
    if (!settings_.immersed_enabled || !internal_boundary_condition_) {
        return;
    }

    internal_boundary_condition_->ApplyPressureVelocityMomentumCorrections(
        state_, workspace_, mesh_, steady_, dt, nu);
}

double PressureVelocitySolver::ComputeDt(const double t_cur) const {
    const auto& ux = state_.Ux();
    const auto& vy = state_.Vy();
    const auto& dx = mesh_.Dx();
    const auto& dy = mesh_.Dy();

    const int i0 = mesh_.GetCoreStartX();
    const int i1 = mesh_.GetCoreEndExclusiveX();
    const int j0 = mesh_.GetCoreStartY();
    const int j1 = mesh_.GetCoreEndExclusiveY();

    const double eps = 1e-14;
    const double u_char_fallback = 1.0;

    double dt_min = std::numeric_limits<double>::max();

    for (int j = j0; j < j1; ++j) {
        for (int i = i0; i < i1; ++i) {
            if (!mesh_.IsFluidCell(i, j, 0)) {
                continue;
            }

            const double dx_i = std::max(dx(i), eps);
            const double dy_j = std::max(dy(j), eps);

            const double u = 0.5 * (ux(i, j, 0) + ux(i + 1, j, 0));
            const double v = 0.5 * (vy(i, j, 0) + vy(i, j + 1, 0));

            double convective = std::abs(u) / dx_i + std::abs(v) / dy_j;
            if (convective < eps) {
                convective = u_char_fallback / std::min(dx_i, dy_j);
            }

            const double local_dt = cfl_ / convective;
            dt_min = std::min(dt_min, local_dt);
        }
    }

    if (mpi_context_) {
        dt_min = mpi_context_->GlobalMin(dt_min);
    }

    if (!std::isfinite(dt_min) || dt_min <= 0.0) {
        return 0.0;
    }

    if (settings_.t_end > 0.0 && t_cur + dt_min > settings_.t_end) {
        dt_min = settings_.t_end - t_cur;
    }

    return std::max(dt_min, 0.0);
}

void PressureVelocitySolver::BuildMomentumCoefficients(const double dt, const double nu) {
    ApplyHaloAndPhysicalBc();
    ApplyImmersedVelocityConstraints();

    workspace_.ZeroMomentumCoefficients();
    workspace_.ZeroAPu();
    workspace_.ZeroAPv();

    auto& au_e = workspace_.AuE();
    auto& au_w = workspace_.AuW();
    auto& au_n = workspace_.AuN();
    auto& au_s = workspace_.AuS();

    auto& av_e = workspace_.AvE();
    auto& av_w = workspace_.AvW();
    auto& av_n = workspace_.AvN();
    auto& av_s = workspace_.AvS();

    auto& a_pu = workspace_.APu();
    auto& a_pv = workspace_.APv();

    const auto& ux_old = state_.UxOld();
    const auto& vy_old = state_.VyOld();
    const auto& dx = mesh_.Dx();
    const auto& dy = mesh_.Dy();

    const bool x_periodic =
        boundary_manager_->IsPeriodic(Axis::X, Side::Left) &&
        boundary_manager_->IsPeriodic(Axis::X, Side::Right);

    const bool y_periodic =
        mesh_.GetDim() >= 2 &&
        boundary_manager_->IsPeriodic(Axis::Y, Side::Left) &&
        boundary_manager_->IsPeriodic(Axis::Y, Side::Right);

    const double rho = 1.0;
    const double inv_dt = (!steady_) ? 1.0 / std::max(dt, 1e-14) : 0.0;

    const int i0 = mesh_.GetCoreStartX();
    const int i1 = mesh_.GetCoreEndExclusiveX();
    const int j0 = mesh_.GetCoreStartY();
    const int j1 = mesh_.GetCoreEndExclusiveY();

    const int u_begin = x_periodic ? i0 : (i0 + 1);
    const int u_end = x_periodic ? i1 : i1;
    const int v_begin = y_periodic ? j0 : (j0 + 1);
    const int v_end = y_periodic ? j1 : j1;

    for (int j = j0; j < j1; ++j) {
        for (int i = u_begin; i <= u_end; ++i) {
            int left_cell = i - 1;
            int right_cell = i;

            if (x_periodic && (i == i0 || i == i1)) {
                left_cell = i1 - 1;
                right_cell = i0;
            }

            if (!(mesh_.IsFluidCell(left_cell, j, 0) && mesh_.IsFluidCell(right_cell, j, 0))) {
                continue;
            }

            const double dx_u = 0.5 * (dx(left_cell) + dx(right_cell));
            const double dy_u = dy(j);

            const double ue = 0.5 * (ux_old(i, j, 0) + ux_old(i + 1, j, 0));
            const double uw = 0.5 * (ux_old(i, j, 0) + ux_old(i - 1, j, 0));
            const double vn = 0.5 * (vy_old(left_cell, j + 1, 0) + vy_old(right_cell, j + 1, 0));
            const double vs = 0.5 * (vy_old(left_cell, j, 0) + vy_old(right_cell, j, 0));

            const double fe = rho * ue * dy_u;
            const double fw = rho * uw * dy_u;
            const double fn = rho * vn * dx_u;
            const double fs = rho * vs * dx_u;

            const double de = nu * dy_u / std::max(dx_u, 1e-14);
            const double dw = de;
            const double dn = nu * dx_u / std::max(dy_u, 1e-14);
            const double ds = dn;

            au_e(i, j, 0) = std::max(-fe, 0.0) + de;
            au_w(i, j, 0) = std::max(fw, 0.0) + dw;
            au_n(i, j, 0) = std::max(-fn, 0.0) + dn;
            au_s(i, j, 0) = std::max(fs, 0.0) + ds;

            const double transient = (!steady_) ? rho * dx_u * dy_u * inv_dt : 0.0;

            a_pu(i, j, 0) =
                au_e(i, j, 0) +
                au_w(i, j, 0) +
                au_n(i, j, 0) +
                au_s(i, j, 0) +
                transient;
        }
    }

    for (int j = v_begin; j <= v_end; ++j) {
        for (int i = i0; i < i1; ++i) {
            int bottom_cell = j - 1;
            int top_cell = j;

            if (y_periodic && (j == j0 || j == j1)) {
                bottom_cell = j1 - 1;
                top_cell = j0;
            }

            if (!(mesh_.IsFluidCell(i, bottom_cell, 0) && mesh_.IsFluidCell(i, top_cell, 0))) {
                continue;
            }

            const double dx_v = dx(i);
            const double dy_v = 0.5 * (dy(bottom_cell) + dy(top_cell));

            const double ue = 0.5 * (ux_old(i + 1, bottom_cell, 0) + ux_old(i + 1, top_cell, 0));
            const double uw = 0.5 * (ux_old(i, bottom_cell, 0) + ux_old(i, top_cell, 0));
            const double vn = 0.5 * (vy_old(i, j, 0) + vy_old(i, j + 1, 0));
            const double vs = 0.5 * (vy_old(i, j, 0) + vy_old(i, j - 1, 0));

            const double fe = rho * ue * dy_v;
            const double fw = rho * uw * dy_v;
            const double fn = rho * vn * dx_v;
            const double fs = rho * vs * dx_v;

            const double de = nu * dy_v / std::max(dx_v, 1e-14);
            const double dw = de;
            const double dn = nu * dx_v / std::max(dy_v, 1e-14);
            const double ds = dn;

            av_e(i, j, 0) = std::max(-fe, 0.0) + de;
            av_w(i, j, 0) = std::max(fw, 0.0) + dw;
            av_n(i, j, 0) = std::max(-fn, 0.0) + dn;
            av_s(i, j, 0) = std::max(fs, 0.0) + ds;

            const double transient = (!steady_) ? rho * dx_v * dy_v * inv_dt : 0.0;

            a_pv(i, j, 0) =
                av_e(i, j, 0) +
                av_w(i, j, 0) +
                av_n(i, j, 0) +
                av_s(i, j, 0) +
                transient;
        }
    }

    boundary_manager_->ApplyPressureVelocityBoundary(
        state_,
        workspace_,
        mesh_,
        PvAssemblyStage::MomentumCoefficients,
        steady_,
        dt,
        nu);

    ApplyImmersedMomentumCorrections(dt, nu);
}

void PressureVelocitySolver::SolveMomentumPredictor(const double dt, const double alpha_u) {
    auto& ux_star = workspace_.UxStar();
    auto& vy_star = workspace_.VyStar();

    ux_star = state_.UxOld();
    vy_star = state_.VyOld();

    state_.Ux() = ux_star;
    state_.Vy() = vy_star;
    ApplyHaloAndPhysicalBc();
    ApplyImmersedVelocityConstraints();
    ux_star = state_.Ux();
    vy_star = state_.Vy();

    const auto& p_old = state_.PressureOld();

    const auto& au_e = workspace_.AuE();
    const auto& au_w = workspace_.AuW();
    const auto& au_n = workspace_.AuN();
    const auto& au_s = workspace_.AuS();

    const auto& av_e = workspace_.AvE();
    const auto& av_w = workspace_.AvW();
    const auto& av_n = workspace_.AvN();
    const auto& av_s = workspace_.AvS();

    const auto& a_pu = workspace_.APu();
    const auto& a_pv = workspace_.APv();

    const auto& ux_old = state_.UxOld();
    const auto& vy_old = state_.VyOld();
    const auto& dx = mesh_.Dx();
    const auto& dy = mesh_.Dy();

    const bool x_periodic =
        boundary_manager_->IsPeriodic(Axis::X, Side::Left) &&
        boundary_manager_->IsPeriodic(Axis::X, Side::Right);

    const bool y_periodic =
        mesh_.GetDim() >= 2 &&
        boundary_manager_->IsPeriodic(Axis::Y, Side::Left) &&
        boundary_manager_->IsPeriodic(Axis::Y, Side::Right);

    const double rho = 1.0;
    const double inv_dt = (!steady_) ? 1.0 / std::max(dt, 1e-14) : 0.0;

    const int i0 = mesh_.GetCoreStartX();
    const int i1 = mesh_.GetCoreEndExclusiveX();
    const int j0 = mesh_.GetCoreStartY();
    const int j1 = mesh_.GetCoreEndExclusiveY();

    const int u_begin = x_periodic ? i0 : (i0 + 1);
    const int u_end = x_periodic ? i1 : (i1 - 1);
    const int v_begin = y_periodic ? j0 : (j0 + 1);
    const int v_end = y_periodic ? j1 : (j1 - 1);

    const int n_sweeps = steady_ ? 3 : 2;

    for (int sweep = 0; sweep < n_sweeps; ++sweep) {
        for (int j = j0; j < j1; ++j) {
            for (int i = u_begin; i <= u_end; ++i) {
                int left_cell = i - 1;
                int right_cell = i;

                if (x_periodic && (i == i0 || i == i1)) {
                    left_cell = i1 - 1;
                    right_cell = i0;
                }

                if (!(mesh_.IsFluidCell(left_cell, j, 0) && mesh_.IsFluidCell(right_cell, j, 0))) {
                    continue;
                }

                const double dx_u = 0.5 * (dx(left_cell) + dx(right_cell));
                const double dy_u = dy(j);

                double h =
                    au_e(i, j, 0) * ux_star(i + 1, j, 0) +
                    au_w(i, j, 0) * ux_star(i - 1, j, 0) +
                    au_n(i, j, 0) * ux_star(i, j + 1, 0) +
                    au_s(i, j, 0) * ux_star(i, j - 1, 0);

                if (!steady_) {
                    h += rho * dx_u * dy_u * inv_dt * ux_old(i, j, 0);
                }

                const double ap = std::max(a_pu(i, j, 0), 1e-14);
                const double raw_u =
                    (h - (p_old(right_cell, j, 0) - p_old(left_cell, j, 0)) * dy_u) / ap;

                ux_star(i, j, 0) =
                    (1.0 - alpha_u) * ux_old(i, j, 0) + alpha_u * raw_u;
            }
        }

        for (int j = v_begin; j <= v_end; ++j) {
            for (int i = i0; i < i1; ++i) {
                int bottom_cell = j - 1;
                int top_cell = j;

                if (y_periodic && (j == j0 || j == j1)) {
                    bottom_cell = j1 - 1;
                    top_cell = j0;
                }

                if (!(mesh_.IsFluidCell(i, bottom_cell, 0) && mesh_.IsFluidCell(i, top_cell, 0))) {
                    continue;
                }

                const double dx_v = dx(i);
                const double dy_v = 0.5 * (dy(bottom_cell) + dy(top_cell));

                double h =
                    av_e(i, j, 0) * vy_star(i + 1, j, 0) +
                    av_w(i, j, 0) * vy_star(i - 1, j, 0) +
                    av_n(i, j, 0) * vy_star(i, j + 1, 0) +
                    av_s(i, j, 0) * vy_star(i, j - 1, 0);

                if (!steady_) {
                    h += rho * dx_v * dy_v * inv_dt * vy_old(i, j, 0);
                }

                const double ap = std::max(a_pv(i, j, 0), 1e-14);
                const double raw_v =
                    (h - (p_old(i, top_cell, 0) - p_old(i, bottom_cell, 0)) * dx_v) / ap;

                vy_star(i, j, 0) =
                    (1.0 - alpha_u) * vy_old(i, j, 0) + alpha_u * raw_v;
            }
        }

        state_.Ux() = ux_star;
        state_.Vy() = vy_star;
        ApplyHaloAndPhysicalBc();
        ApplyImmersedVelocityConstraints();
        ux_star = state_.Ux();
        vy_star = state_.Vy();
    }
}

void PressureVelocitySolver::BuildPressureCorrectionEquation(const double dt, const double nu) {
    workspace_.ZeroPressureCorrection();
    workspace_.ZeroPressureRhs();
    workspace_.ZeroPressureCoefficients();

    auto& ap_e = workspace_.ApE();
    auto& ap_w = workspace_.ApW();
    auto& ap_n = workspace_.ApN();
    auto& ap_s = workspace_.ApS();
    auto& ap_p = workspace_.ApP();
    auto& rhs = workspace_.PressureRhs();

    const auto& ux_star = workspace_.UxStar();
    const auto& vy_star = workspace_.VyStar();
    const auto& a_pu = workspace_.APu();
    const auto& a_pv = workspace_.APv();
    const auto& dx = mesh_.Dx();
    const auto& dy = mesh_.Dy();

    const bool x_periodic =
        boundary_manager_->IsPeriodic(Axis::X, Side::Left) &&
        boundary_manager_->IsPeriodic(Axis::X, Side::Right);

    const bool y_periodic =
        mesh_.GetDim() >= 2 &&
        boundary_manager_->IsPeriodic(Axis::Y, Side::Left) &&
        boundary_manager_->IsPeriodic(Axis::Y, Side::Right);

    const int i0 = mesh_.GetCoreStartX();
    const int i1 = mesh_.GetCoreEndExclusiveX();
    const int j0 = mesh_.GetCoreStartY();
    const int j1 = mesh_.GetCoreEndExclusiveY();

    for (int j = j0; j < j1; ++j) {
        for (int i = i0; i < i1; ++i) {
            if (!mesh_.IsFluidCell(i, j, 0)) {
                continue;
            }

            const double dx_c = dx(i);
            const double dy_c = dy(j);

            double ae = 0.0;
            double aw = 0.0;
            double an = 0.0;
            double as = 0.0;

            if (i + 1 < i1) {
                if (mesh_.IsFluidCell(i + 1, j, 0)) {
                    ae = (dy_c * dy_c) / std::max(a_pu(i + 1, j, 0), 1e-14);
                }
            }
            else if (x_periodic && mesh_.IsFluidCell(i0, j, 0)) {
                ae = (dy_c * dy_c) / std::max(a_pu(i1, j, 0), 1e-14);
            }

            if (i - 1 >= i0) {
                if (mesh_.IsFluidCell(i - 1, j, 0)) {
                    aw = (dy_c * dy_c) / std::max(a_pu(i, j, 0), 1e-14);
                }
            }
            else if (x_periodic && mesh_.IsFluidCell(i1 - 1, j, 0)) {
                aw = (dy_c * dy_c) / std::max(a_pu(i0, j, 0), 1e-14);
            }

            if (j + 1 < j1) {
                if (mesh_.IsFluidCell(i, j + 1, 0)) {
                    an = (dx_c * dx_c) / std::max(a_pv(i, j + 1, 0), 1e-14);
                }
            }
            else if (y_periodic && mesh_.IsFluidCell(i, j0, 0)) {
                an = (dx_c * dx_c) / std::max(a_pv(i, j1, 0), 1e-14);
            }

            if (j - 1 >= j0) {
                if (mesh_.IsFluidCell(i, j - 1, 0)) {
                    as = (dx_c * dx_c) / std::max(a_pv(i, j, 0), 1e-14);
                }
            }
            else if (y_periodic && mesh_.IsFluidCell(i, j1 - 1, 0)) {
                as = (dx_c * dx_c) / std::max(a_pv(i, j0, 0), 1e-14);
            }

            ap_e(i, j, 0) = ae;
            ap_w(i, j, 0) = aw;
            ap_n(i, j, 0) = an;
            ap_s(i, j, 0) = as;
            ap_p(i, j, 0) = ae + aw + an + as;

            rhs(i, j, 0) =
                -((ux_star(i + 1, j, 0) - ux_star(i, j, 0)) * dy_c +
                    (vy_star(i, j + 1, 0) - vy_star(i, j, 0)) * dx_c);
        }
    }

    boundary_manager_->ApplyPressureVelocityBoundary(
        state_,
        workspace_,
        mesh_,
        PvAssemblyStage::PressureCorrectionEquation,
        steady_,
        dt,
        nu);

    double sum_rhs = 0.0;
    int n_fluid = 0;
    for (int j = j0; j < j1; ++j) {
        for (int i = i0; i < i1; ++i) {
            if (mesh_.IsFluidCell(i, j, 0)) {
                sum_rhs += rhs(i, j, 0);
                n_fluid++;
            }
        }
    }

    if (mpi_context_ && n_fluid > 0) {
        sum_rhs = mpi_context_->GlobalSum(sum_rhs);
        n_fluid = mpi_context_->GlobalSum(n_fluid);
    }

    if (n_fluid > 0) {
        const double mean_rhs = sum_rhs / static_cast<double>(n_fluid);
        for (int j = j0; j < j1; ++j) {
            for (int i = i0; i < i1; ++i) {
                if (mesh_.IsFluidCell(i, j, 0)) {
                    rhs(i, j, 0) -= mean_rhs;
                }
            }
        }
    }
}

int PressureVelocitySolver::PressureLinearIndex(const int i, const int j) const {
    const int nx = mesh_.GetCoreEndExclusiveX() - mesh_.GetCoreStartX();
    return (j - mesh_.GetCoreStartY()) * nx + (i - mesh_.GetCoreStartX());
}

void PressureVelocitySolver::SolvePressureCorrection() {
    const int i0 = mesh_.GetCoreStartX();
    const int i1 = mesh_.GetCoreEndExclusiveX();
    const int j0 = mesh_.GetCoreStartY();
    const int j1 = mesh_.GetCoreEndExclusiveY();

    const bool x_periodic =
        boundary_manager_->IsPeriodic(Axis::X, Side::Left) &&
        boundary_manager_->IsPeriodic(Axis::X, Side::Right);

    const bool y_periodic =
        mesh_.GetDim() >= 2 &&
        boundary_manager_->IsPeriodic(Axis::Y, Side::Left) &&
        boundary_manager_->IsPeriodic(Axis::Y, Side::Right);

    const int nx = i1 - i0;
    const int ny = j1 - j0;
    const int np = nx * ny;

    Eigen::SparseMatrix<double> a(np, np);
    Eigen::VectorXd b(np);
    b.setZero();

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(static_cast<std::size_t>(np) * 5);

    const auto& ap_e = workspace_.ApE();
    const auto& ap_w = workspace_.ApW();
    const auto& ap_n = workspace_.ApN();
    const auto& ap_s = workspace_.ApS();
    const auto& ap_p = workspace_.ApP();
    const auto& rhs = workspace_.PressureRhs();

    int gauge_i = -1;
    int gauge_j = -1;
    for (int j = j0; j < j1 && gauge_i < 0; ++j) {
        for (int i = i0; i < i1; ++i) {
            if (mesh_.IsFluidCell(i, j, 0)) {
                gauge_i = i;
                gauge_j = j;
                break;
            }
        }
    }

    if (gauge_i < 0) {
        throw std::runtime_error("PressureVelocitySolver: no fluid cells found for pressure gauge");
    }

    const int gauge_row = PressureLinearIndex(gauge_i, gauge_j);

    for (int j = j0; j < j1; ++j) {
        for (int i = i0; i < i1; ++i) {
            const int row = PressureLinearIndex(i, j);

            if (!mesh_.IsFluidCell(i, j, 0)) {
                triplets.emplace_back(row, row, 1.0);
                b(row) = 0.0;
                continue;
            }

            if (row == gauge_row) {
                triplets.emplace_back(row, row, 1.0);
                b(row) = 0.0;
                continue;
            }

            const double app = ap_p(i, j, 0);

            if (app <= 0.0) {
                triplets.emplace_back(row, row, 1.0);
                b(row) = 0.0;
                continue;
            }

            triplets.emplace_back(row, row, app);
            b(row) = rhs(i, j, 0);

            if (ap_e(i, j, 0) != 0.0) {
                int ie = i + 1;
                if (ie >= i1) {
                    ie = x_periodic ? i0 : -1;
                }
                if (ie >= i0 && ie < i1 && mesh_.IsFluidCell(ie, j, 0)) {
                    int target_row = PressureLinearIndex(ie, j);
                    if (target_row != gauge_row) {
                        // ВОТ ОНО! Сохраняем симметрию
                        triplets.emplace_back(row, target_row, -ap_e(i, j, 0));
                    }
                }
            }

            if (ap_w(i, j, 0) != 0.0) {
                int iw = i - 1;
                if (iw < i0) {
                    iw = x_periodic ? (i1 - 1) : -1;
                }
                if (iw >= i0 && iw < i1 && mesh_.IsFluidCell(iw, j, 0)) {
                    int target_row = PressureLinearIndex(iw, j);
                    if (target_row != gauge_row) {
                        triplets.emplace_back(row, target_row, -ap_w(i, j, 0));
                    }
                }
            }

            if (ap_n(i, j, 0) != 0.0) {
                int jn = j + 1;
                if (jn >= j1) {
                    jn = y_periodic ? j0 : -1;
                }
                if (jn >= j0 && jn < j1 && mesh_.IsFluidCell(i, jn, 0)) {
                    int target_row = PressureLinearIndex(i, jn);
                    if (target_row != gauge_row) {
                        triplets.emplace_back(row, target_row, -ap_n(i, j, 0));
                    }
                }
            }

            if (ap_s(i, j, 0) != 0.0) {
                int js = j - 1;
                if (js < j0) {
                    js = y_periodic ? (j1 - 1) : -1;
                }
                if (js >= j0 && js < j1 && mesh_.IsFluidCell(i, js, 0)) {
                    int target_row = PressureLinearIndex(i, js);
                    if (target_row != gauge_row) {
                        triplets.emplace_back(row, target_row, -ap_s(i, j, 0));
                    }
                }
            }
        }
    }

    a.setFromTriplets(triplets.begin(), triplets.end());

    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower | Eigen::Upper> solver;
    solver.setTolerance(1e-10);
    solver.setMaxIterations(std::max(500, np * 4));
    solver.compute(a);

    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("PressureVelocitySolver: pressure correction matrix factorization failed");
    }

    const Eigen::VectorXd x = solver.solve(b);

    if (solver.info() != Eigen::Success || !std::isfinite(solver.error())) {
        throw std::runtime_error("PressureVelocitySolver: pressure correction solve failed");
    }

    auto& p_corr = workspace_.PressureCorrection();
    p_corr.fill(0.0);

    for (int j = j0; j < j1; ++j) {
        for (int i = i0; i < i1; ++i) {
            if (!mesh_.IsFluidCell(i, j, 0)) {
                continue;
            }
            p_corr(i, j, 0) = x(PressureLinearIndex(i, j));
        }
    }
}

void PressureVelocitySolver::ApplyPressureCorrection(const double alpha_p) {
    auto& p = state_.Pressure();
    auto& ux = state_.Ux();
    auto& vy = state_.Vy();

    const auto& p_corr = workspace_.PressureCorrection();
    const auto& ux_star = workspace_.UxStar();
    const auto& vy_star = workspace_.VyStar();
    const auto& a_pu = workspace_.APu();
    const auto& a_pv = workspace_.APv();
    const auto& dx = mesh_.Dx();
    const auto& dy = mesh_.Dy();

    const bool x_periodic =
        boundary_manager_->IsPeriodic(Axis::X, Side::Left) &&
        boundary_manager_->IsPeriodic(Axis::X, Side::Right);

    const bool y_periodic =
        mesh_.GetDim() >= 2 &&
        boundary_manager_->IsPeriodic(Axis::Y, Side::Left) &&
        boundary_manager_->IsPeriodic(Axis::Y, Side::Right);

    const int i0 = mesh_.GetCoreStartX();
    const int i1 = mesh_.GetCoreEndExclusiveX();
    const int j0 = mesh_.GetCoreStartY();
    const int j1 = mesh_.GetCoreEndExclusiveY();

    for (int j = j0; j < j1; ++j) {
        for (int i = i0; i < i1; ++i) {
            if (!mesh_.IsFluidCell(i, j, 0)) {
                continue;
            }
            p(i, j, 0) += alpha_p * p_corr(i, j, 0);
        }
    }

    const int u_begin = x_periodic ? i0 : (i0 + 1);
    const int u_end = x_periodic ? i1 : (i1 - 1);

    for (int j = j0; j < j1; ++j) {
        for (int i = u_begin; i <= u_end; ++i) {
            int left_cell = i - 1;
            int right_cell = i;

            if (x_periodic && (i == i0 || i == i1)) {
                left_cell = i1 - 1;
                right_cell = i0;
            }

            if (!(mesh_.IsFluidCell(left_cell, j, 0) && mesh_.IsFluidCell(right_cell, j, 0))) {
                continue;
            }

            ux(i, j, 0) =
                ux_star(i, j, 0) -
                dy(j) / std::max(a_pu(i, j, 0), 1e-14) *
                (p_corr(right_cell, j, 0) - p_corr(left_cell, j, 0));
        }
    }

    const int v_begin = y_periodic ? j0 : (j0 + 1);
    const int v_end = y_periodic ? j1 : (j1 - 1);

    for (int j = v_begin; j <= v_end; ++j) {
        for (int i = i0; i < i1; ++i) {
            int bottom_cell = j - 1;
            int top_cell = j;

            if (y_periodic && (j == j0 || j == j1)) {
                bottom_cell = j1 - 1;
                top_cell = j0;
            }

            if (!(mesh_.IsFluidCell(i, bottom_cell, 0) && mesh_.IsFluidCell(i, top_cell, 0))) {
                continue;
            }

            vy(i, j, 0) =
                vy_star(i, j, 0) -
                dx(i) / std::max(a_pv(i, j, 0), 1e-14) *
                (p_corr(i, top_cell, 0) - p_corr(i, bottom_cell, 0));
        }
    }

    ApplyHaloAndPhysicalBc();
    ApplyImmersedVelocityConstraints();
}

void PressureVelocitySolver::CopyStarToState() {
    state_.Ux() = workspace_.UxStar();
    state_.Vy() = workspace_.VyStar();
    ApplyHaloAndPhysicalBc();
    ApplyImmersedVelocityConstraints();
}

double PressureVelocitySolver::ComputeMassResidual() const {
    const auto& ux = state_.Ux();
    const auto& vy = state_.Vy();
    const auto& dx = mesh_.Dx();
    const auto& dy = mesh_.Dy();

    const int i0 = mesh_.GetCoreStartX();
    const int i1 = mesh_.GetCoreEndExclusiveX();
    const int j0 = mesh_.GetCoreStartY();
    const int j1 = mesh_.GetCoreEndExclusiveY();

    double rmax = 0.0;

    for (int j = j0; j < j1; ++j) {
        for (int i = i0; i < i1; ++i) {
            if (!mesh_.IsFluidCell(i, j, 0)) {
                continue;
            }

            const double div =
                (ux(i + 1, j, 0) - ux(i, j, 0)) * dy(j) +
                (vy(i, j + 1, 0) - vy(i, j, 0)) * dx(i);

            rmax = std::max(rmax, std::abs(div));
        }
    }

    if (mpi_context_) {
        rmax = mpi_context_->GlobalMax(rmax);
    }

    return rmax;
}
