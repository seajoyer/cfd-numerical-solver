#include "spatial/ForwardEulerSpatialOperator.hpp"

#include <stdexcept>

#include "bc/BoundaryManager.hpp"
#include "data/Variables.hpp"
#include "viscosity/VNRArtificialViscosity.hpp"

ForwardEulerSpatialOperator::ForwardEulerSpatialOperator(
    const Settings& settings,
    std::shared_ptr<BoundaryManager> boundary_manager)
    : SpatialOperator(std::move(boundary_manager)) {
    viscosity_ = nullptr;
    if (settings.viscosity) {
        viscosity_ = std::make_shared<VNRArtificialViscosity>(settings);
    }
}

PrimitiveCell ForwardEulerSpatialOperator::LoadPrimitive(const xt::xtensor<double, 4>& W,
                                                         const int i,
                                                         const int j,
                                                         const int k) {
    PrimitiveCell w;
    w.rho = W(var::u_rho, i, j, k);
    w.u = W(var::u_u, i, j, k);
    w.v = W(var::u_v, i, j, k);
    w.w = W(var::u_w, i, j, k);
    w.P = W(var::u_P, i, j, k);
    return w;
}

const xt::xtensor<double, 1>& ForwardEulerSpatialOperator::InvMetric(const Mesh& mesh,
                                                                     const Axis axis) const {
    if (axis == Axis::X) {
        return mesh.InvDx();
    }
    if (axis == Axis::Y) {
        return mesh.InvDy();
    }
    return mesh.InvDz();
}

void ForwardEulerSpatialOperator::AccumulateAxis(const DataLayer& layer,
                                                 const Mesh& mesh,
                                                 const xt::xtensor<double, 4>& W,
                                                 xt::xtensor<double, 4>& rhs,
                                                 const double gamma,
                                                 const Axis axis) const {
    (void)layer;

    const AxisStride st = AxisStride::FromAxis(axis);
    const auto& inv_h = InvMetric(mesh, axis);

    const int i0 = mesh.GetCoreStartX();
    const int i1 = mesh.GetCoreEndExclusiveX();
    const int j0 = mesh.GetCoreStartY();
    const int j1 = mesh.GetCoreEndExclusiveY();
    const int k0 = mesh.GetCoreStartZ();
    const int k1 = mesh.GetCoreEndExclusiveZ();

    for (int k = k0; k < k1; ++k) {
        for (int j = j0; j < j1; ++j) {
            for (int i = i0; i < i1; ++i) {
                if (!mesh.IsFluidCell(i, j, k)) {
                    continue;
                }

                const int ip = i + st.di;
                const int jp = j + st.dj;
                const int kp = k + st.dk;

                if (!mesh.IsFluidCell(ip, jp, kp)) {
                    continue;
                }

                const PrimitiveCell w0 = LoadPrimitive(W, i, j, k);
                const PrimitiveCell w1 = LoadPrimitive(W, ip, jp, kp);

                const FluxCell F0 = EulerFlux(w0, gamma, axis);
                const FluxCell F1 = EulerFlux(w1, gamma, axis);

                const int idx = axis == Axis::X ? i : (axis == Axis::Y ? j : k);
                const double inv = inv_h(static_cast<std::size_t>(idx));

                rhs(DataLayer::k_rho, i, j, k) += -(F1.mass - F0.mass) * inv;
                rhs(DataLayer::k_rhoU, i, j, k) += -(F1.mom_x - F0.mom_x) * inv;
                rhs(DataLayer::k_rhoV, i, j, k) += -(F1.mom_y - F0.mom_y) * inv;
                rhs(DataLayer::k_rhoW, i, j, k) += -(F1.mom_z - F0.mom_z) * inv;
                rhs(DataLayer::k_E, i, j, k) += -(F1.energy - F0.energy) * inv;
            }
        }
    }
}

void ForwardEulerSpatialOperator::ComputeRHS(DataLayer& layer,
                                             const Mesh& mesh,
                                             Workspace& workspace,
                                             const double gamma,
                                             const double dt) const {
    if (!boundary_manager_) {
        throw std::runtime_error("ForwardEulerSpatialOperator: boundary_manager_ is null");
    }

    const int ng = mesh.GetPadding();
    if (ng < 1) {
        throw std::runtime_error("ForwardEulerSpatialOperator: requires at least 1 ghost cell (ng>=1)");
    }

    workspace.ResizeFrom(mesh);

    boundary_manager_->UpdateHalo(layer, mesh);
    boundary_manager_->ApplyPhysicalBc(layer, mesh);

    ConvertUtoW(layer.U(), workspace.W(), gamma,
                0, mesh.GetSx(),
                0, mesh.GetSy(),
                0, mesh.GetSz());

    workspace.ZeroRhs();

    auto& rhs = workspace.Rhs();
    const auto& W = workspace.W();

    AccumulateAxis(layer, mesh, W, rhs, gamma, Axis::X);
    if (mesh.GetDim() >= 2) {
        AccumulateAxis(layer, mesh, W, rhs, gamma, Axis::Y);
    }
    if (mesh.GetDim() >= 3) {
        AccumulateAxis(layer, mesh, W, rhs, gamma, Axis::Z);
    }

    if (viscosity_) {
        viscosity_->AddToRhs(layer, mesh, workspace.W(), gamma, dt, workspace.Rhs());
    }
}