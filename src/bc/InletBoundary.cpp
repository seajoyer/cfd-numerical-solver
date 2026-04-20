#include "bc/InletBoundary.hpp"

#include <algorithm>

#include "data/DataLayer.hpp"
#include "data/Mesh.hpp"
#include "data/PressureVelocityState.hpp"
#include "data/PressureVelocityWorkspace.hpp"
#include "data/Variables.hpp"

namespace {
    void SetPlaneToStateX(xt::xtensor<double, 4>& U, const int dst_i, const FarfieldConservative& s) {
        xt::view(U, DataLayer::k_rho, dst_i, xt::all(), xt::all()) = s.rho;
        xt::view(U, DataLayer::k_rhoU, dst_i, xt::all(), xt::all()) = s.rhoU;
        xt::view(U, DataLayer::k_rhoV, dst_i, xt::all(), xt::all()) = s.rhoV;
        xt::view(U, DataLayer::k_rhoW, dst_i, xt::all(), xt::all()) = s.rhoW;
        xt::view(U, DataLayer::k_E, dst_i, xt::all(), xt::all()) = s.E;
    }

    void SetPlaneToStateY(xt::xtensor<double, 4>& U, const int dst_j, const FarfieldConservative& s) {
        xt::view(U, DataLayer::k_rho, xt::all(), dst_j, xt::all()) = s.rho;
        xt::view(U, DataLayer::k_rhoU, xt::all(), dst_j, xt::all()) = s.rhoU;
        xt::view(U, DataLayer::k_rhoV, xt::all(), dst_j, xt::all()) = s.rhoV;
        xt::view(U, DataLayer::k_rhoW, xt::all(), dst_j, xt::all()) = s.rhoW;
        xt::view(U, DataLayer::k_E, xt::all(), dst_j, xt::all()) = s.E;
    }

    void SetPlaneToStateZ(xt::xtensor<double, 4>& U, const int dst_k, const FarfieldConservative& s) {
        xt::view(U, DataLayer::k_rho, xt::all(), xt::all(), dst_k) = s.rho;
        xt::view(U, DataLayer::k_rhoU, xt::all(), xt::all(), dst_k) = s.rhoU;
        xt::view(U, DataLayer::k_rhoV, xt::all(), xt::all(), dst_k) = s.rhoV;
        xt::view(U, DataLayer::k_rhoW, xt::all(), xt::all(), dst_k) = s.rhoW;
        xt::view(U, DataLayer::k_E, xt::all(), xt::all(), dst_k) = s.E;
    }
}

InletBoundary::InletBoundary(const FarfieldConservative& inflow_U)
    : inflow_U_(inflow_U) {}

InletBoundary::InletBoundary(const BoundaryStateSettings& primitive_state)
    : primitive_state_(primitive_state) {}

void InletBoundary::Apply(DataLayer& layer, const Mesh& mesh, const Axis axis, const Side side) const {
    const int ng = mesh.GetPadding();
    if (ng == 0) {
        return;
    }

    const int dim = mesh.GetDim();
    if (axis == Axis::Y && dim < 2) {
        return;
    }
    if (axis == Axis::Z && dim < 3) {
        return;
    }

    auto& U = layer.U();

    const int i0 = mesh.GetCoreStartX();
    const int i1 = mesh.GetCoreEndExclusiveX();
    const int j0 = mesh.GetCoreStartY();
    const int j1 = mesh.GetCoreEndExclusiveY();
    const int k0 = mesh.GetCoreStartZ();
    const int k1 = mesh.GetCoreEndExclusiveZ();

    for (int g = 0; g < ng; ++g) {
        if (axis == Axis::X) {
            const int dst_i = side == Side::Left ? (i0 - 1 - g) : (i1 + g);
            SetPlaneToStateX(U, dst_i, inflow_U_);
            continue;
        }

        if (axis == Axis::Y) {
            const int dst_j = side == Side::Left ? (j0 - 1 - g) : (j1 + g);
            SetPlaneToStateY(U, dst_j, inflow_U_);
            continue;
        }

        const int dst_k = side == Side::Left ? (k0 - 1 - g) : (k1 + g);
        SetPlaneToStateZ(U, dst_k, inflow_U_);
    }
}

void InletBoundary::Apply(PressureVelocityState& state,
                          const Mesh& mesh,
                          const Axis axis,
                          const Side side) const {
    const int ng = mesh.GetPadding();
    if (ng == 0) {
        return;
    }

    const int dim = mesh.GetDim();
    if (axis == Axis::Y && dim < 2) {
        return;
    }
    if (axis == Axis::Z && dim < 3) {
        return;
    }

    auto& p = state.Pressure();
    auto& ux = state.Ux();
    auto& vy = state.Vy();
    auto& wz = state.Wz();

    const int sx = mesh.GetSx();
    const int sy = mesh.GetSy();
    const int sz = mesh.GetSz();

    const int i0 = mesh.GetCoreStartX();
    const int i1 = mesh.GetCoreEndExclusiveX();
    const int j0 = mesh.GetCoreStartY();
    const int j1 = mesh.GetCoreEndExclusiveY();
    const int k0 = mesh.GetCoreStartZ();
    const int k1 = mesh.GetCoreEndExclusiveZ();

    const double u_in = primitive_state_.u;
    const double v_in = primitive_state_.v;
    const double w_in = primitive_state_.w;
    const double p_in = primitive_state_.p;

    if (axis == Axis::X && side == Side::Left) {
        for (int j = 0; j < sy; ++j) {
            for (int k = 0; k < sz; ++k) {
                const bool open_segment = mesh.IsFluidCell(i0, j, 0);

                for (int g = 0; g < ng; ++g) {
                    const int ig = i0 - 1 - g;
                    const int ii = i0 + g;

                    if (open_segment) {
                        p(ig, j, k) = 2.0 * p_in - p(ii, j, k);
                    }
                    else {
                        p(ig, j, k) = p(ii, j, k);
                    }
                }

                ux(i0, j, k) = open_segment ? u_in : 0.0;
                for (int g = 1; g <= ng; ++g) {
                    const int ig = i0 - g;
                    ux(ig, j, k) = open_segment ? u_in : 0.0;
                }
            }
        }

        if (dim >= 2) {
            for (int j = 0; j < sy + 1; ++j) {
                for (int k = 0; k < sz; ++k) {
                    const int jc = std::clamp(j, 0, sy - 1);
                    const bool open_segment = mesh.IsFluidCell(i0, jc, 0);

                    for (int g = 0; g < ng; ++g) {
                        const int ig = i0 - 1 - g;
                        const int ii = i0 + g;

                        if (open_segment) {
                            vy(ig, j, k) = 2.0 * v_in - vy(ii, j, k);
                        }
                        else {
                            vy(ig, j, k) = -vy(ii, j, k);
                        }
                    }
                }
            }
        }

        if (dim >= 3) {
            for (int j = 0; j < sy; ++j) {
                for (int k = 0; k < sz + 1; ++k) {
                    const bool open_segment = mesh.IsFluidCell(i0, j, 0);

                    for (int g = 0; g < ng; ++g) {
                        const int ig = i0 - 1 - g;
                        const int ii = i0 + g;

                        if (open_segment) {
                            wz(ig, j, k) = 2.0 * w_in - wz(ii, j, k);
                        }
                        else {
                            wz(ig, j, k) = -wz(ii, j, k);
                        }
                    }
                }
            }
        }

        return;
    }

    if (axis == Axis::X && side == Side::Right) {
        for (int j = 0; j < sy; ++j) {
            for (int k = 0; k < sz; ++k) {
                const bool open_segment = mesh.IsFluidCell(i1 - 1, j, 0);

                for (int g = 0; g < ng; ++g) {
                    const int ig = i1 + g;
                    const int ii = i1 - 1 - g;

                    if (open_segment) {
                        p(ig, j, k) = 2.0 * p_in - p(ii, j, k);
                    }
                    else {
                        p(ig, j, k) = p(ii, j, k);
                    }
                }

                ux(i1, j, k) = open_segment ? u_in : 0.0;
                for (int g = 1; g <= ng; ++g) {
                    const int ig = i1 + g;
                    ux(ig, j, k) = open_segment ? u_in : 0.0;
                }
            }
        }

        if (dim >= 2) {
            for (int j = 0; j < sy + 1; ++j) {
                for (int k = 0; k < sz; ++k) {
                    const int jc = std::clamp(j, 0, sy - 1);
                    const bool open_segment = mesh.IsFluidCell(i1 - 1, jc, 0);

                    for (int g = 0; g < ng; ++g) {
                        const int ig = i1 + g;
                        const int ii = i1 - 1 - g;

                        if (open_segment) {
                            vy(ig, j, k) = 2.0 * v_in - vy(ii, j, k);
                        }
                        else {
                            vy(ig, j, k) = -vy(ii, j, k);
                        }
                    }
                }
            }
        }

        if (dim >= 3) {
            for (int j = 0; j < sy; ++j) {
                for (int k = 0; k < sz + 1; ++k) {
                    const bool open_segment = mesh.IsFluidCell(i1 - 1, j, 0);

                    for (int g = 0; g < ng; ++g) {
                        const int ig = i1 + g;
                        const int ii = i1 - 1 - g;

                        if (open_segment) {
                            wz(ig, j, k) = 2.0 * w_in - wz(ii, j, k);
                        }
                        else {
                            wz(ig, j, k) = -wz(ii, j, k);
                        }
                    }
                }
            }
        }

        return;
    }

    if (axis == Axis::Y && side == Side::Left) {
        for (int i = 0; i < sx; ++i) {
            for (int k = 0; k < sz; ++k) {
                const bool open_segment = mesh.IsFluidCell(i, j0, 0);

                for (int g = 0; g < ng; ++g) {
                    const int jg = j0 - 1 - g;
                    const int ji = j0 + g;

                    if (open_segment) {
                        p(i, jg, k) = 2.0 * p_in - p(i, ji, k);
                    }
                    else {
                        p(i, jg, k) = p(i, ji, k);
                    }
                }

                vy(i, j0, k) = open_segment ? v_in : 0.0;
                for (int g = 1; g <= ng; ++g) {
                    const int jg = j0 - g;
                    vy(i, jg, k) = open_segment ? v_in : 0.0;
                }
            }
        }

        for (int i = 0; i < sx + 1; ++i) {
            for (int k = 0; k < sz; ++k) {
                const int ic = std::clamp(i, 0, sx - 1);
                const bool open_segment = mesh.IsFluidCell(ic, j0, 0);

                for (int g = 0; g < ng; ++g) {
                    const int jg = j0 - 1 - g;
                    const int ji = j0 + g;

                    if (open_segment) {
                        ux(i, jg, k) = 2.0 * u_in - ux(i, ji, k);
                    }
                    else {
                        ux(i, jg, k) = -ux(i, ji, k);
                    }
                }
            }
        }

        if (dim >= 3) {
            for (int i = 0; i < sx; ++i) {
                for (int k = 0; k < sz + 1; ++k) {
                    const bool open_segment = mesh.IsFluidCell(i, j0, 0);

                    for (int g = 0; g < ng; ++g) {
                        const int jg = j0 - 1 - g;
                        const int ji = j0 + g;

                        if (open_segment) {
                            wz(i, jg, k) = 2.0 * w_in - wz(i, ji, k);
                        }
                        else {
                            wz(i, jg, k) = -wz(i, ji, k);
                        }
                    }
                }
            }
        }

        return;
    }

    if (axis == Axis::Y && side == Side::Right) {
        for (int i = 0; i < sx; ++i) {
            for (int k = 0; k < sz; ++k) {
                const bool open_segment = mesh.IsFluidCell(i, j1 - 1, 0);

                for (int g = 0; g < ng; ++g) {
                    const int jg = j1 + g;
                    const int ji = j1 - 1 - g;

                    if (open_segment) {
                        p(i, jg, k) = 2.0 * p_in - p(i, ji, k);
                    }
                    else {
                        p(i, jg, k) = p(i, ji, k);
                    }
                }

                vy(i, j1, k) = open_segment ? v_in : 0.0;
                for (int g = 1; g <= ng; ++g) {
                    const int jg = j1 + g;
                    vy(i, jg, k) = open_segment ? v_in : 0.0;
                }
            }
        }

        for (int i = 0; i < sx + 1; ++i) {
            for (int k = 0; k < sz; ++k) {
                const int ic = std::clamp(i, 0, sx - 1);
                const bool open_segment = mesh.IsFluidCell(ic, j1 - 1, 0);

                for (int g = 0; g < ng; ++g) {
                    const int jg = j1 + g;
                    const int ji = j1 - 1 - g;

                    if (open_segment) {
                        ux(i, jg, k) = 2.0 * u_in - ux(i, ji, k);
                    }
                    else {
                        ux(i, jg, k) = -ux(i, ji, k);
                    }
                }
            }
        }

        if (dim >= 3) {
            for (int i = 0; i < sx; ++i) {
                for (int k = 0; k < sz + 1; ++k) {
                    const bool open_segment = mesh.IsFluidCell(i, j1 - 1, 0);

                    for (int g = 0; g < ng; ++g) {
                        const int jg = j1 + g;
                        const int ji = j1 - 1 - g;

                        if (open_segment) {
                            wz(i, jg, k) = 2.0 * w_in - wz(i, ji, k);
                        }
                        else {
                            wz(i, jg, k) = -wz(i, ji, k);
                        }
                    }
                }
            }
        }

        return;
    }
}

void InletBoundary::ApplyPressureVelocityBoundary(PressureVelocityState& state,
                                                  PressureVelocityWorkspace& workspace,
                                                  const Mesh& mesh,
                                                  const Axis axis,
                                                  const Side side,
                                                  const PvAssemblyStage stage,
                                                  const bool steady,
                                                  const double dt,
                                                  const double nu) const {
    (void)state;
    (void)steady;
    (void)dt;

    if (stage == PvAssemblyStage::PressureCorrectionEquation) {
        return;
    }

    const int dim = mesh.GetDim();
    const int i0 = mesh.GetCoreStartX();
    const int i1 = mesh.GetCoreEndExclusiveX();
    const int j0 = mesh.GetCoreStartY();
    const int j1 = mesh.GetCoreEndExclusiveY();

    const auto& dx = mesh.Dx();
    const auto& dy = mesh.Dy();

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

    if (axis == Axis::Y && dim >= 2) {
        if (side == Side::Left) {
            const int j = j0;
            for (int i = i0 + 1; i < i1; ++i) {
                if (!(mesh.IsFluidCell(i - 1, j, 0) && mesh.IsFluidCell(i, j, 0))) {
                    continue;
                }

                const double dx_u = 0.5 * (dx(i - 1) + dx(i));
                const double ds = nu * dx_u / std::max(0.5 * dy(j), 1e-14);

                a_pu(i, j, 0) += ds - au_s(i, j, 0);
                au_s(i, j, 0) = ds;
            }

            return;
        }

        {
            const int j = j1 - 1;
            for (int i = i0 + 1; i < i1; ++i) {
                if (!(mesh.IsFluidCell(i - 1, j, 0) && mesh.IsFluidCell(i, j, 0))) {
                    continue;
                }

                const double dx_u = 0.5 * (dx(i - 1) + dx(i));
                const double dn = nu * dx_u / std::max(0.5 * dy(j), 1e-14);

                a_pu(i, j, 0) += dn - au_n(i, j, 0);
                au_n(i, j, 0) = dn;
            }

            return;
        }
    }

    if (axis == Axis::X) {
        if (side == Side::Left) {
            const int i = i0;
            for (int j = j0 + 1; j < j1; ++j) {
                if (!(mesh.IsFluidCell(i, j - 1, 0) && mesh.IsFluidCell(i, j, 0))) {
                    continue;
                }

                const double dy_v = 0.5 * (dy(j - 1) + dy(j));
                const double dw = nu * dy_v / std::max(0.5 * dx(i), 1e-14);

                a_pv(i, j, 0) += dw - av_w(i, j, 0);
                av_w(i, j, 0) = dw;
            }

            return;
        }

        {
            const int i = i1 - 1;
            for (int j = j0 + 1; j < j1; ++j) {
                if (!(mesh.IsFluidCell(i, j - 1, 0) && mesh.IsFluidCell(i, j, 0))) {
                    continue;
                }

                const double dy_v = 0.5 * (dy(j - 1) + dy(j));
                const double de = nu * dy_v / std::max(0.5 * dx(i), 1e-14);

                a_pv(i, j, 0) += de - av_e(i, j, 0);
                av_e(i, j, 0) = de;
            }

            return;
        }
    }
}
