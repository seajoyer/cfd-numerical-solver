#include "bc/WallBoundary.hpp"

#include <algorithm>

#include "data/DataLayer.hpp"
#include "data/Mesh.hpp"
#include "data/PressureVelocityState.hpp"
#include "data/PressureVelocityWorkspace.hpp"
#include "data/Variables.hpp"

namespace {
    void CopyScalarsAndZeroMomentaX(xt::xtensor<double, 4>& U, const int dst_i, const int src_i) {
        xt::view(U, DataLayer::k_rho, dst_i, xt::all(), xt::all()) =
            xt::view(U, DataLayer::k_rho, src_i, xt::all(), xt::all());
        xt::view(U, DataLayer::k_E, dst_i, xt::all(), xt::all()) =
            xt::view(U, DataLayer::k_E, src_i, xt::all(), xt::all());

        xt::view(U, DataLayer::k_rhoU, dst_i, xt::all(), xt::all()) = 0.0;
        xt::view(U, DataLayer::k_rhoV, dst_i, xt::all(), xt::all()) = 0.0;
        xt::view(U, DataLayer::k_rhoW, dst_i, xt::all(), xt::all()) = 0.0;
    }

    void CopyScalarsAndZeroMomentaY(xt::xtensor<double, 4>& U, const int dst_j, const int src_j) {
        xt::view(U, DataLayer::k_rho, xt::all(), dst_j, xt::all()) =
            xt::view(U, DataLayer::k_rho, xt::all(), src_j, xt::all());
        xt::view(U, DataLayer::k_E, xt::all(), dst_j, xt::all()) =
            xt::view(U, DataLayer::k_E, xt::all(), src_j, xt::all());

        xt::view(U, DataLayer::k_rhoU, xt::all(), dst_j, xt::all()) = 0.0;
        xt::view(U, DataLayer::k_rhoV, xt::all(), dst_j, xt::all()) = 0.0;
        xt::view(U, DataLayer::k_rhoW, xt::all(), dst_j, xt::all()) = 0.0;
    }

    void CopyScalarsAndZeroMomentaZ(xt::xtensor<double, 4>& U, const int dst_k, const int src_k) {
        xt::view(U, DataLayer::k_rho, xt::all(), xt::all(), dst_k) =
            xt::view(U, DataLayer::k_rho, xt::all(), xt::all(), src_k);
        xt::view(U, DataLayer::k_E, xt::all(), xt::all(), dst_k) =
            xt::view(U, DataLayer::k_E, xt::all(), xt::all(), src_k);

        xt::view(U, DataLayer::k_rhoU, xt::all(), xt::all(), dst_k) = 0.0;
        xt::view(U, DataLayer::k_rhoV, xt::all(), xt::all(), dst_k) = 0.0;
        xt::view(U, DataLayer::k_rhoW, xt::all(), xt::all(), dst_k) = 0.0;
    }
}

void WallBoundary::Apply(DataLayer& layer, const Mesh& mesh, const Axis axis, const Side side) const {
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

    if (axis == Axis::X) {
        const int src_i = side == Side::Left ? i0 : (i1 - 1);
        for (int g = 0; g < ng; ++g) {
            const int dst_i = side == Side::Left ? (i0 - 1 - g) : (i1 + g);
            CopyScalarsAndZeroMomentaX(U, dst_i, src_i);
        }
        return;
    }

    if (axis == Axis::Y) {
        const int src_j = side == Side::Left ? j0 : (j1 - 1);
        for (int g = 0; g < ng; ++g) {
            const int dst_j = side == Side::Left ? (j0 - 1 - g) : (j1 + g);
            CopyScalarsAndZeroMomentaY(U, dst_j, src_j);
        }
        return;
    }

    const int src_k = side == Side::Left ? k0 : (k1 - 1);
    for (int g = 0; g < ng; ++g) {
        const int dst_k = side == Side::Left ? (k0 - 1 - g) : (k1 + g);
        CopyScalarsAndZeroMomentaZ(U, dst_k, src_k);
    }
}

void WallBoundary::Apply(PressureVelocityState& state,
                         const Mesh& mesh,
                         const Axis axis,
                         const Side side) const {
    const int ng = mesh.GetPadding();
    if (ng == 0) {
        return;
    }

    auto& p = state.Pressure();
    auto& ux = state.Ux();
    auto& vy = state.Vy();
    auto& wz = state.Wz();

    const int i0 = mesh.GetCoreStartX();
    const int i1 = mesh.GetCoreEndExclusiveX();
    const int j0 = mesh.GetCoreStartY();
    const int j1 = mesh.GetCoreEndExclusiveY();
    const int k0 = mesh.GetCoreStartZ();
    const int k1 = mesh.GetCoreEndExclusiveZ();

    if (axis == Axis::X) {
        if (side == Side::Left) {
            for (int g = 0; g < ng; ++g) {
                const int ig = i0 - 1 - g;
                const int ii = i0 + g;
                xt::view(p, ig, xt::all(), xt::all()) = xt::view(p, ii, xt::all(), xt::all());
            }

            xt::view(ux, i0, xt::all(), xt::all()) = 0.0;
            xt::view(ux, xt::range(i0 - ng, i0), xt::all(), xt::all()) = 0.0;

            if (mesh.GetDim() >= 2) {
                for (int g = 0; g < ng; ++g) {
                    const int ig = i0 - 1 - g;
                    const int ii = i0 + g;
                    xt::view(vy, ig, xt::all(), xt::all()) =
                        -xt::view(vy, ii, xt::all(), xt::all());
                }
            }

            if (mesh.GetDim() >= 3) {
                for (int g = 0; g < ng; ++g) {
                    const int ig = i0 - 1 - g;
                    const int ii = i0 + g;
                    xt::view(wz, ig, xt::all(), xt::all()) =
                        -xt::view(wz, ii, xt::all(), xt::all());
                }
            }

            return;
        }

        for (int g = 0; g < ng; ++g) {
            const int ig = i1 + g;
            const int ii = i1 - 1 - g;
            xt::view(p, ig, xt::all(), xt::all()) = xt::view(p, ii, xt::all(), xt::all());
        }

        xt::view(ux, i1, xt::all(), xt::all()) = 0.0;
        xt::view(ux, xt::range(i1 + 1, i1 + ng + 1), xt::all(), xt::all()) = 0.0;

        if (mesh.GetDim() >= 2) {
            for (int g = 0; g < ng; ++g) {
                const int ig = i1 + g;
                const int ii = i1 - 1 - g;
                xt::view(vy, ig, xt::all(), xt::all()) =
                    -xt::view(vy, ii, xt::all(), xt::all());
            }
        }

        if (mesh.GetDim() >= 3) {
            for (int g = 0; g < ng; ++g) {
                const int ig = i1 + g;
                const int ii = i1 - 1 - g;
                xt::view(wz, ig, xt::all(), xt::all()) =
                    -xt::view(wz, ii, xt::all(), xt::all());
            }
        }

        return;
    }

    if (axis == Axis::Y) {
        if (side == Side::Left) {
            for (int g = 0; g < ng; ++g) {
                const int jg = j0 - 1 - g;
                const int ji = j0 + g;
                xt::view(p, xt::all(), jg, xt::all()) = xt::view(p, xt::all(), ji, xt::all());
            }

            xt::view(vy, xt::all(), j0, xt::all()) = 0.0;
            xt::view(vy, xt::all(), xt::range(j0 - ng, j0), xt::all()) = 0.0;

            for (int g = 0; g < ng; ++g) {
                const int jg = j0 - 1 - g;
                const int ji = j0 + g;
                xt::view(ux, xt::all(), jg, xt::all()) =
                    -xt::view(ux, xt::all(), ji, xt::all());
            }

            if (mesh.GetDim() >= 3) {
                for (int g = 0; g < ng; ++g) {
                    const int jg = j0 - 1 - g;
                    const int ji = j0 + g;
                    xt::view(wz, xt::all(), jg, xt::all()) =
                        -xt::view(wz, xt::all(), ji, xt::all());
                }
            }

            return;
        }

        for (int g = 0; g < ng; ++g) {
            const int jg = j1 + g;
            const int ji = j1 - 1 - g;
            xt::view(p, xt::all(), jg, xt::all()) = xt::view(p, xt::all(), ji, xt::all());
        }

        xt::view(vy, xt::all(), j1, xt::all()) = 0.0;
        xt::view(vy, xt::all(), xt::range(j1 + 1, j1 + ng + 1), xt::all()) = 0.0;

        for (int g = 0; g < ng; ++g) {
            const int jg = j1 + g;
            const int ji = j1 - 1 - g;
            xt::view(ux, xt::all(), jg, xt::all()) =
                -xt::view(ux, xt::all(), ji, xt::all());
        }

        if (mesh.GetDim() >= 3) {
            for (int g = 0; g < ng; ++g) {
                const int jg = j1 + g;
                const int ji = j1 - 1 - g;
                xt::view(wz, xt::all(), jg, xt::all()) =
                    -xt::view(wz, xt::all(), ji, xt::all());
            }
        }

        return;
    }

    if (axis == Axis::Z) {
        if (side == Side::Left) {
            for (int g = 0; g < ng; ++g) {
                const int kg = k0 - 1 - g;
                const int ki = k0 + g;
                xt::view(p, xt::all(), xt::all(), kg) = xt::view(p, xt::all(), xt::all(), ki);
            }

            xt::view(wz, xt::all(), xt::all(), k0) = 0.0;
            xt::view(wz, xt::all(), xt::all(), xt::range(k0 - ng, k0)) = 0.0;

            xt::view(ux, xt::all(), xt::all(), xt::range(k0 - ng, k0)) =
                -xt::view(ux, xt::all(), xt::all(), xt::range(k0, k0 + ng));
            xt::view(vy, xt::all(), xt::all(), xt::range(k0 - ng, k0)) =
                -xt::view(vy, xt::all(), xt::all(), xt::range(k0, k0 + ng));

            return;
        }

        for (int g = 0; g < ng; ++g) {
            const int kg = k1 + g;
            const int ki = k1 - 1 - g;
            xt::view(p, xt::all(), xt::all(), kg) = xt::view(p, xt::all(), xt::all(), ki);
        }

        xt::view(wz, xt::all(), xt::all(), k1) = 0.0;
        xt::view(wz, xt::all(), xt::all(), xt::range(k1 + 1, k1 + ng + 1)) = 0.0;

        xt::view(ux, xt::all(), xt::all(), xt::range(k1, k1 + ng)) =
            -xt::view(ux, xt::all(), xt::all(), xt::range(k1 - ng, k1));
        xt::view(vy, xt::all(), xt::all(), xt::range(k1, k1 + ng)) =
            -xt::view(vy, xt::all(), xt::all(), xt::range(k1 - ng, k1));
    }
}

void WallBoundary::ApplyPressureVelocityBoundary(PressureVelocityState& state,
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
