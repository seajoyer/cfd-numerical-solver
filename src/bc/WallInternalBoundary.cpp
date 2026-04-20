#include "bc/WallInternalBoundary.hpp"

#include <algorithm>
#include <cmath>

void WallInternalBoundary::BuildBoundaryState(const PrimitiveCell& fluid_state,
                                              const ImmersedFaceInfo& face_info,
                                              PrimitiveCell& boundary_state) const {
    if (!face_info.is_active) {
        boundary_state = fluid_state;
        return;
    }

    const double nx_in = face_info.normal_x;
    const double ny_in = face_info.normal_y;
    const double nz_in = face_info.normal_z;

    const double norm2 = nx_in * nx_in + ny_in * ny_in + nz_in * nz_in;
    if (norm2 <= k_eps_ * k_eps_) {
        boundary_state = fluid_state;
        return;
    }

    const double norm = std::sqrt(norm2);
    const double nx = nx_in / norm;
    const double ny = ny_in / norm;
    const double nz = nz_in / norm;

    const double un = fluid_state.u * nx + fluid_state.v * ny + fluid_state.w * nz;
    const double ut_x = fluid_state.u - un * nx;
    const double ut_y = fluid_state.v - un * ny;
    const double ut_z = fluid_state.w - un * nz;

    boundary_state.rho = fluid_state.rho;
    boundary_state.P = fluid_state.P;

    boundary_state.u = -ut_x - un * nx;
    boundary_state.v = -ut_y - un * ny;
    boundary_state.w = -ut_z - un * nz;
}

void WallInternalBoundary::ApplyPressureVelocityVelocityConstraints(PressureVelocityState& state,
                                                                    PressureVelocityWorkspace& workspace,
                                                                    const Mesh& mesh) const {
    auto& ux = state.Ux();
    auto& vy = state.Vy();
    auto& ux_star = workspace.UxStar();
    auto& vy_star = workspace.VyStar();

    const int i0 = mesh.GetCoreStartX();
    const int i1 = mesh.GetCoreEndExclusiveX();
    const int j0 = mesh.GetCoreStartY();
    const int j1 = mesh.GetCoreEndExclusiveY();

    for (int j = j0; j < j1; ++j) {
        for (int i = i0; i <= i1; ++i) {
            const bool has_left = (i - 1 >= 0 && i - 1 < mesh.GetSx());
            const bool has_right = (i >= 0 && i < mesh.GetSx());

            if (!has_left || !has_right) {
                continue;
            }

            const bool left_fluid = mesh.IsFluidCell(i - 1, j, 0);
            const bool right_fluid = mesh.IsFluidCell(i, j, 0);

            if (left_fluid != right_fluid) {
                ux(i, j, 0) = 0.0;
                ux_star(i, j, 0) = 0.0;
            }
        }
    }

    for (int j = j0; j <= j1; ++j) {
        for (int i = i0; i < i1; ++i) {
            const bool has_bottom = (j - 1 >= 0 && j - 1 < mesh.GetSy());
            const bool has_top = (j >= 0 && j < mesh.GetSy());

            if (!has_bottom || !has_top) {
                continue;
            }

            const bool bottom_fluid = mesh.IsFluidCell(i, j - 1, 0);
            const bool top_fluid = mesh.IsFluidCell(i, j, 0);

            if (bottom_fluid != top_fluid) {
                vy(i, j, 0) = 0.0;
                vy_star(i, j, 0) = 0.0;
            }
        }
    }
}

void WallInternalBoundary::ApplyPressureVelocityMomentumCorrections(PressureVelocityState& state,
                                                                    PressureVelocityWorkspace& workspace,
                                                                    const Mesh& mesh,
                                                                    const bool steady,
                                                                    const double dt,
                                                                    const double nu) const {
    (void)state;
    (void)steady;
    (void)dt;

    auto& au_e = workspace.AuE();
    auto& au_w = workspace.AuW();
    auto& au_n = workspace.AuN();
    auto& au_s = workspace.AuS();
    auto& av_e = workspace.AvE();
    auto& av_w = workspace.AvW();
    auto& av_n = workspace.AvN();
    auto& av_s = workspace.AvS();
    auto& a_pu = workspace.APu();
    auto& a_pv = workspace.APv();

    const auto& dx = mesh.Dx();
    const auto& dy = mesh.Dy();

    const int i0 = mesh.GetCoreStartX();
    const int i1 = mesh.GetCoreEndExclusiveX();
    const int j0 = mesh.GetCoreStartY();
    const int j1 = mesh.GetCoreEndExclusiveY();

    for (int j = j0; j < j1; ++j) {
        for (int i = i0 + 1; i < i1; ++i) {
            if (!(mesh.IsFluidCell(i - 1, j, 0) && mesh.IsFluidCell(i, j, 0))) {
                continue;
            }

            const double dx_u = 0.5 * (dx(i - 1) + dx(i));
            const double dy_u = dy(j);

            if (j + 1 < mesh.GetSy()) {
                const bool north_left_fluid = mesh.IsFluidCell(i - 1, j + 1, 0);
                const bool north_right_fluid = mesh.IsFluidCell(i, j + 1, 0);
                if (!(north_left_fluid && north_right_fluid)) {
                    const double dn = nu * dx_u / std::max(0.5 * dy_u, 1e-14);
                    a_pu(i, j, 0) += dn - au_n(i, j, 0);
                    au_n(i, j, 0) = dn;
                }
            }

            if (j - 1 >= 0) {
                const bool south_left_fluid = mesh.IsFluidCell(i - 1, j - 1, 0);
                const bool south_right_fluid = mesh.IsFluidCell(i, j - 1, 0);
                if (!(south_left_fluid && south_right_fluid)) {
                    const double ds = nu * dx_u / std::max(0.5 * dy_u, 1e-14);
                    a_pu(i, j, 0) += ds - au_s(i, j, 0);
                    au_s(i, j, 0) = ds;
                }
            }
        }
    }

    for (int j = j0 + 1; j < j1; ++j) {
        for (int i = i0; i < i1; ++i) {
            if (!(mesh.IsFluidCell(i, j - 1, 0) && mesh.IsFluidCell(i, j, 0))) {
                continue;
            }

            const double dx_v = dx(i);
            const double dy_v = 0.5 * (dy(j - 1) + dy(j));

            if (i + 1 < mesh.GetSx()) {
                const bool east_bottom_fluid = mesh.IsFluidCell(i + 1, j - 1, 0);
                const bool east_top_fluid = mesh.IsFluidCell(i + 1, j, 0);
                if (!(east_bottom_fluid && east_top_fluid)) {
                    const double de = nu * dy_v / std::max(0.5 * dx_v, 1e-14);
                    a_pv(i, j, 0) += de - av_e(i, j, 0);
                    av_e(i, j, 0) = de;
                }
            }

            if (i - 1 >= 0) {
                const bool west_bottom_fluid = mesh.IsFluidCell(i - 1, j - 1, 0);
                const bool west_top_fluid = mesh.IsFluidCell(i - 1, j, 0);
                if (!(west_bottom_fluid && west_top_fluid)) {
                    const double dw = nu * dy_v / std::max(0.5 * dx_v, 1e-14);
                    a_pv(i, j, 0) += dw - av_w(i, j, 0);
                    av_w(i, j, 0) = dw;
                }
            }
        }
    }
}
