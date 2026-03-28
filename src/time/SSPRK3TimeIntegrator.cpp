#include "time/SSPRK3TimeIntegrator.hpp"

#include "data/DataLayer.hpp"
#include "data/Workspace.hpp"
#include "geometry/Mesh.hpp"
#include "solver/PositivityLimiter.hpp"
#include "spatial/SpatialOperator.hpp"

void SSPRK3TimeIntegrator::Advance(DataLayer& layer,
                                   const Mesh& mesh,
                                   Workspace& workspace,
                                   const double dt,
                                   const double gamma,
                                   const SpatialOperator& op) const {
    if (dt <= 0.0) {
        return;
    }

    workspace.ResizeFrom(mesh);

    auto& U = layer.U();
    const xt::xtensor<double, 2> U0 = U;

    op.ComputeRHS(layer, mesh, workspace, gamma, dt);

    const auto& rhs_stage1 = workspace.Rhs();
    for (std::size_t cell_id = 0; cell_id < mesh.GetCellCount(); ++cell_id) {
        for (std::size_t var = 0; var < DataLayer::k_nvar; ++var) {
            U(cell_id, var) = U0(cell_id, var) + dt * rhs_stage1(cell_id, var);
        }
    }

    op.ComputeRHS(layer, mesh, workspace, gamma, dt);

    const auto& rhs_stage2 = workspace.Rhs();
    for (std::size_t cell_id = 0; cell_id < mesh.GetCellCount(); ++cell_id) {
        for (std::size_t var = 0; var < DataLayer::k_nvar; ++var) {
            U(cell_id, var) =
                3.0 / 4.0 * U0(cell_id, var) +
                1.0 / 4.0 * (U(cell_id, var) + dt * rhs_stage2(cell_id, var));
        }
    }

    op.ComputeRHS(layer, mesh, workspace, gamma, dt);

    const auto& rhs_stage3 = workspace.Rhs();
    for (std::size_t cell_id = 0; cell_id < mesh.GetCellCount(); ++cell_id) {
        for (std::size_t var = 0; var < DataLayer::k_nvar; ++var) {
            U(cell_id, var) =
                1.0 / 3.0 * U0(cell_id, var) +
                2.0 / 3.0 * (U(cell_id, var) + dt * rhs_stage3(cell_id, var));
        }
    }

    PositivityLimiter::Apply(layer, mesh, gamma, rho_min_, p_min_);
}
