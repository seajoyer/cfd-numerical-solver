#include "time/MaderTimeIntegrator.hpp"

#include <stdexcept>

#include "data/DataLayer.hpp"
#include "data/Mesh.hpp"
#include "data/Workspace.hpp"
#include "spatial/MaderSpatialOperator.hpp"
#include "spatial/SpatialOperator.hpp"

void MaderTimeIntegrator::Advance(DataLayer& layer,
                                  const Mesh& mesh,
                                  Workspace& workspace,
                                  const double dt,
                                  const double gamma,
                                  const SpatialOperator& op) const {
    const auto* mader_op = dynamic_cast<const MaderSpatialOperator*>(&op);
    if (!mader_op) {
        throw std::runtime_error(
            "MaderTimeIntegrator::Advance requires SpatialOperator of type MaderSpatialOperator");
    }

    if (dt <= 0.0) {
        return;
    }

    // Full split step of Mader 2DE:
    mader_op->Phase1_EosAndChemistry(layer, mesh, workspace, gamma, dt);
    mader_op->Phase2_PressureForces(layer, mesh, workspace, gamma, dt);
    mader_op->Phase3_ZipEnergy(layer, mesh, workspace, gamma, dt);
    mader_op->Phase4_Transport(layer, mesh, workspace, gamma, dt);
    mader_op->Phase5_Finalize(layer, mesh, workspace, gamma, dt);
}
