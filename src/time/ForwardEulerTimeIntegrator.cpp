#include "time/ForwardEulerTimeIntegrator.hpp"

#include "data/DataLayer.hpp"
#include "data/Workspace.hpp"
#include "geometry/Mesh.hpp"
#include "solver/PositivityLimiter.hpp"
#include "spatial/SpatialOperator.hpp"

void ForwardEulerTimeIntegrator::Advance(DataLayer& layer,
                                         const Mesh& mesh,
                                         Workspace& workspace,
                                         const double dt,
                                         const double gamma,
                                         const SpatialOperator& op) const {
    if (dt <= 0.0) {
        return;
    }

    workspace.ResizeFrom(mesh);

    op.ComputeRHS(layer, mesh, workspace, gamma, dt);

    auto& U = layer.U();
    const auto& rhs = workspace.Rhs();

    for (std::size_t cell_id = 0; cell_id < mesh.GetCellCount(); ++cell_id) {
        for (std::size_t var = 0; var < DataLayer::k_nvar; ++var) {
            U(cell_id, var) += dt * rhs(cell_id, var);
        }
    }

    PositivityLimiter::Apply(layer, mesh, gamma, rho_min_, p_min_);
}
