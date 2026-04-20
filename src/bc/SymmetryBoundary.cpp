#include "bc/SymmetryBoundary.hpp"

#include "data/DataLayer.hpp"
#include "data/Mesh.hpp"
#include "data/PressureVelocityState.hpp"
#include "data/PressureVelocityWorkspace.hpp"
#include "data/Variables.hpp"

void SymmetryBoundary::Apply(DataLayer& layer, const Mesh& mesh, const Axis axis, const Side side) const {
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
            const int src_i = side == Side::Left ? (i0 + g) : (i1 - 1 - g);

            xt::view(U, xt::all(), dst_i, xt::all(), xt::all()) =
                xt::view(U, xt::all(), src_i, xt::all(), xt::all());

            xt::view(U, DataLayer::k_rhoU, dst_i, xt::all(), xt::all()) *= -1.0;
            continue;
        }

        if (axis == Axis::Y) {
            const int dst_j = side == Side::Left ? (j0 - 1 - g) : (j1 + g);
            const int src_j = side == Side::Left ? (j0 + g) : (j1 - 1 - g);

            xt::view(U, xt::all(), xt::all(), dst_j, xt::all()) =
                xt::view(U, xt::all(), xt::all(), src_j, xt::all());

            xt::view(U, DataLayer::k_rhoV, xt::all(), dst_j, xt::all()) *= -1.0;
            continue;
        }

        const int dst_k = side == Side::Left ? (k0 - 1 - g) : (k1 + g);
        const int src_k = side == Side::Left ? (k0 + g) : (k1 - 1 - g);

        xt::view(U, xt::all(), xt::all(), xt::all(), dst_k) =
            xt::view(U, xt::all(), xt::all(), xt::all(), src_k);

        xt::view(U, DataLayer::k_rhoW, xt::all(), xt::all(), dst_k) *= -1.0;
    }
}

void SymmetryBoundary::Apply(PressureVelocityState& state,
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

                xt::view(p, ig, xt::all(), xt::all()) =
                    xt::view(p, ii, xt::all(), xt::all());

                xt::view(vy, ig, xt::all(), xt::all()) =
                    xt::view(vy, ii, xt::all(), xt::all());

                if (dim >= 3) {
                    xt::view(wz, ig, xt::all(), xt::all()) =
                        xt::view(wz, ii, xt::all(), xt::all());
                }
            }

            xt::view(ux, i0, xt::all(), xt::all()) = 0.0;
            xt::view(ux, xt::range(i0 - ng, i0), xt::all(), xt::all()) = 0.0;
            return;
        }

        for (int g = 0; g < ng; ++g) {
            const int ig = i1 + g;
            const int ii = i1 - 1 - g;

            xt::view(p, ig, xt::all(), xt::all()) =
                xt::view(p, ii, xt::all(), xt::all());

            xt::view(vy, ig, xt::all(), xt::all()) =
                xt::view(vy, ii, xt::all(), xt::all());

            if (dim >= 3) {
                xt::view(wz, ig, xt::all(), xt::all()) =
                    xt::view(wz, ii, xt::all(), xt::all());
            }
        }

        xt::view(ux, i1, xt::all(), xt::all()) = 0.0;
        xt::view(ux, xt::range(i1 + 1, i1 + ng + 1), xt::all(), xt::all()) = 0.0;
        return;
    }

    if (axis == Axis::Y) {
        if (side == Side::Left) {
            for (int g = 0; g < ng; ++g) {
                const int jg = j0 - 1 - g;
                const int ji = j0 + g;

                xt::view(p, xt::all(), jg, xt::all()) =
                    xt::view(p, xt::all(), ji, xt::all());

                xt::view(ux, xt::all(), jg, xt::all()) =
                    xt::view(ux, xt::all(), ji, xt::all());

                if (dim >= 3) {
                    xt::view(wz, xt::all(), jg, xt::all()) =
                        xt::view(wz, xt::all(), ji, xt::all());
                }
            }

            xt::view(vy, xt::all(), j0, xt::all()) = 0.0;
            xt::view(vy, xt::all(), xt::range(j0 - ng, j0), xt::all()) = 0.0;
            return;
        }

        for (int g = 0; g < ng; ++g) {
            const int jg = j1 + g;
            const int ji = j1 - 1 - g;

            xt::view(p, xt::all(), jg, xt::all()) =
                xt::view(p, xt::all(), ji, xt::all());

            xt::view(ux, xt::all(), jg, xt::all()) =
                xt::view(ux, xt::all(), ji, xt::all());

            if (dim >= 3) {
                xt::view(wz, xt::all(), jg, xt::all()) =
                    xt::view(wz, xt::all(), ji, xt::all());
            }
        }

        xt::view(vy, xt::all(), j1, xt::all()) = 0.0;
        xt::view(vy, xt::all(), xt::range(j1 + 1, j1 + ng + 1), xt::all()) = 0.0;
        return;
    }

    if (axis == Axis::Z) {
        if (side == Side::Left) {
            for (int g = 0; g < ng; ++g) {
                const int kg = k0 - 1 - g;
                const int ki = k0 + g;

                xt::view(p, xt::all(), xt::all(), kg) =
                    xt::view(p, xt::all(), xt::all(), ki);

                xt::view(ux, xt::all(), xt::all(), kg) =
                    xt::view(ux, xt::all(), xt::all(), ki);

                xt::view(vy, xt::all(), xt::all(), kg) =
                    xt::view(vy, xt::all(), xt::all(), ki);
            }

            xt::view(wz, xt::all(), xt::all(), k0) = 0.0;
            xt::view(wz, xt::all(), xt::all(), xt::range(k0 - ng, k0)) = 0.0;
            return;
        }

        for (int g = 0; g < ng; ++g) {
            const int kg = k1 + g;
            const int ki = k1 - 1 - g;

            xt::view(p, xt::all(), xt::all(), kg) =
                xt::view(p, xt::all(), xt::all(), ki);

            xt::view(ux, xt::all(), xt::all(), kg) =
                xt::view(ux, xt::all(), xt::all(), ki);

            xt::view(vy, xt::all(), xt::all(), kg) =
                xt::view(vy, xt::all(), xt::all(), ki);
        }

        xt::view(wz, xt::all(), xt::all(), k1) = 0.0;
        xt::view(wz, xt::all(), xt::all(), xt::range(k1 + 1, k1 + ng + 1)) = 0.0;
    }
}

void SymmetryBoundary::ApplyPressureVelocityBoundary(PressureVelocityState& state,
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
    (void)nu;

    if (stage == PvAssemblyStage::PressureCorrectionEquation) {
        return;
    }

    const int dim = mesh.GetDim();
    const int i0 = mesh.GetCoreStartX();
    const int i1 = mesh.GetCoreEndExclusiveX();
    const int j0 = mesh.GetCoreStartY();
    const int j1 = mesh.GetCoreEndExclusiveY();

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

                a_pu(i, j, 0) -= au_s(i, j, 0);
                au_s(i, j, 0) = 0.0;
            }

            return;
        }

        {
            const int j = j1 - 1;
            for (int i = i0 + 1; i < i1; ++i) {
                if (!(mesh.IsFluidCell(i - 1, j, 0) && mesh.IsFluidCell(i, j, 0))) {
                    continue;
                }

                a_pu(i, j, 0) -= au_n(i, j, 0);
                au_n(i, j, 0) = 0.0;
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

                a_pv(i, j, 0) -= av_w(i, j, 0);
                av_w(i, j, 0) = 0.0;
            }

            return;
        }

        {
            const int i = i1 - 1;
            for (int j = j0 + 1; j < j1; ++j) {
                if (!(mesh.IsFluidCell(i, j - 1, 0) && mesh.IsFluidCell(i, j, 0))) {
                    continue;
                }

                a_pv(i, j, 0) -= av_e(i, j, 0);
                av_e(i, j, 0) = 0.0;
            }
        }
    }
}
