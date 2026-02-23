#include "spatial/BackwardEulerSpatialOperator.hpp"

#include <stdexcept>

#include "bc/BoundaryManager.hpp"
#include "data/Variables.hpp"
#include "viscosity/VNRArtificialViscosity.hpp"

BackwardEulerSpatialOperator::BackwardEulerSpatialOperator(
    const Settings& settings,
    std::shared_ptr<BoundaryManager> boundary_manager)
    : SpatialOperator(std::move(boundary_manager)) {
    viscosity_ = nullptr;
    if (settings.viscosity) {
        viscosity_ = std::make_shared<VNRArtificialViscosity>(settings);
    }
}

PrimitiveCell BackwardEulerSpatialOperator::LoadPrimitive(const xt::xtensor<double, 4>& W,
                                                          const int i, const int j, const int k) {
    PrimitiveCell w;
    w.rho = W(var::u_rho, i, j, k);
    w.u = W(var::u_u, i, j, k);
    w.v = W(var::u_v, i, j, k);
    w.w = W(var::u_w, i, j, k);
    w.P = W(var::u_P, i, j, k);
    return w;
}

[[nodiscard]] const xt::xtensor<double, 1>&
BackwardEulerSpatialOperator::InvMetric(const DataLayer& layer, const Axis axis) const {
    if (axis == Axis::X) return layer.InvDx();
    if (axis == Axis::Y) return layer.InvDy();
    return layer.InvDz();
}

void BackwardEulerSpatialOperator::AccumulateAxis(const DataLayer& layer,
                                                  const xt::xtensor<double, 4>& W,
                                                  xt::xtensor<double, 4>& rhs,
                                                  const double gamma,
                                                  const Axis axis) const {
    const AxisStride st = AxisStride::FromAxis(axis);
    const auto& inv_h = InvMetric(layer, axis);

    const int i0 = layer.GetCoreStartX();
    const int i1 = layer.GetCoreEndExclusiveX();
    const int j0 = layer.GetCoreStartY();
    const int j1 = layer.GetCoreEndExclusiveY();
    const int k0 = layer.GetCoreStartZ();
    const int k1 = layer.GetCoreEndExclusiveZ();

    for (int k = k0; k < k1; ++k) {
        for (int j = j0; j < j1; ++j) {
            for (int i = i0; i < i1; ++i) {
                const int im = i - st.di;
                const int jm = j - st.dj;
                const int km = k - st.dk;

                const PrimitiveCell w0 = LoadPrimitive(W, i, j, k);
                const PrimitiveCell wm = LoadPrimitive(W, im, jm, km);

                const FluxCell F0 = EulerFlux(w0, gamma, axis);
                const FluxCell Fm = EulerFlux(wm, gamma, axis);

                const int idx = (axis == Axis::X) ? i : (axis == Axis::Y) ? j : k;
                const double inv = inv_h(static_cast<std::size_t>(idx));

                rhs(DataLayer::k_rho, i, j, k) += -(F0.mass - Fm.mass) * inv;
                rhs(DataLayer::k_rhoU, i, j, k) += -(F0.mom_x - Fm.mom_x) * inv;
                rhs(DataLayer::k_rhoV, i, j, k) += -(F0.mom_y - Fm.mom_y) * inv;
                rhs(DataLayer::k_rhoW, i, j, k) += -(F0.mom_z - Fm.mom_z) * inv;
                rhs(DataLayer::k_E, i, j, k) += -(F0.energy - Fm.energy) * inv;
            }
        }
    }
}

void BackwardEulerSpatialOperator::ComputeRHS(DataLayer& layer, Workspace& workspace, const double gamma,
                                              const double dt) const {
    if (!boundary_manager_) {
        throw std::runtime_error("BackwardEulerSpatialOperator: boundary_manager_ is null");
    }

    const int ng = layer.GetPadding();
    if (ng < 1) {
        throw std::runtime_error("BackwardEulerSpatialOperator: requires at least 1 ghost cell (ng>=1)");
    }

    workspace.ResizeFrom(layer);

    boundary_manager_->UpdateHalo(layer);
    boundary_manager_->ApplyPhysicalBc(layer);

    ConvertUtoW(layer.U(), workspace.W(), gamma,
                0, layer.GetSx(),
                0, layer.GetSy(),
                0, layer.GetSz());

    workspace.ZeroRhs();

    auto& rhs = workspace.Rhs();
    const auto& W = workspace.W();

    AccumulateAxis(layer, W, rhs, gamma, Axis::X);
    if (layer.GetDim() >= 2) AccumulateAxis(layer, W, rhs, gamma, Axis::Y);
    if (layer.GetDim() >= 3) AccumulateAxis(layer, W, rhs, gamma, Axis::Z);

    if (viscosity_) {
        viscosity_->AddToRhs(layer, workspace.W(), gamma, dt, workspace.Rhs());
    }
}
