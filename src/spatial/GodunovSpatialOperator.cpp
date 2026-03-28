#include "spatial/GodunovSpatialOperator.hpp"

#include <cctype>
#include <memory>
#include <stdexcept>
#include <string>

#include "bc/BoundaryManager.hpp"
#include "data/DataLayer.hpp"
#include "data/Variables.hpp"
#include "data/Workspace.hpp"
#include "geometry/Cell.hpp"
#include "geometry/Face.hpp"
#include "geometry/Mesh.hpp"
#include "reconstruction/P0Reconstruction.hpp"
#include "reconstruction/P1Reconstruction.hpp"
#include "reconstruction/Reconstruction.hpp"
#include "riemann/RusanovRiemannSolver.hpp"
#include "riemann/HLLRiemannSolver.hpp"
#include "riemann/HLLCRiemannSolver.hpp"
#include "riemann/RoeRiemannSolver.hpp"
#include "riemann/ExactIdealGasRiemannSolver.hpp"
#include "riemann/RiemannSolver.hpp"

GodunovSpatialOperator::GodunovSpatialOperator(
    const Settings& settings,
    std::shared_ptr<BoundaryManager> boundary_manager
) : SpatialOperator(std::move(boundary_manager)) {
    InitializeReconstruction(settings);
    InitializeRiemannSolver(settings);

    if (!boundary_manager_) {
        throw std::runtime_error("GodunovSpatialOperator: boundary_manager is null");
    }
    if (!reconstruction_) {
        throw std::runtime_error("GodunovSpatialOperator: reconstruction is null");
    }
    if (!riemann_solver_) {
        throw std::runtime_error("GodunovSpatialOperator: riemann_solver is null");
    }
}

void GodunovSpatialOperator::InitializeReconstruction(const Settings& settings) {
    std::string name = settings.reconstruction;
    for (char& c : name) {
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }

    if (name == "p0") {
        reconstruction_ = std::make_shared<P0Reconstruction>();
        return;
    }

    if (name == "p1") {
        reconstruction_ = std::make_shared<P1Reconstruction>();
        return;
    }

    throw std::runtime_error(
        "GodunovSpatialOperator::InitializeReconstruction: unsupported reconstruction '" +
        settings.reconstruction + "'"
    );
}

void GodunovSpatialOperator::InitializeRiemannSolver(const Settings& settings) {
    std::string name = settings.riemann_solver;
    for (char& c : name) {
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }

    if (name == "rusanov") {
        riemann_solver_ = std::make_shared<RusanovRiemannSolver>();
        return;
    }

    if (name == "hll") {
        riemann_solver_ = std::make_shared<HLLCRiemannSolver>();
        return;
    }

    if (name == "hllc") {
        riemann_solver_ = std::make_shared<HLLCRiemannSolver>();
        return;
    }

    if (name == "roe") {
        riemann_solver_ = std::make_shared<RoeRiemannSolver>();
        return;
    }

    if (name == "exact") {
        riemann_solver_ = std::make_shared<ExactIdealGasRiemannSolver>(0.0, settings.Q_user);
        return;
    }

    throw std::runtime_error(
        "GodunovSpatialOperator::InitializeRiemannSolver: unsupported riemann solver '" +
        settings.riemann_solver + "'"
    );
}

void GodunovSpatialOperator::FillPrimitiveCache(const DataLayer& layer,
                                                const Mesh& mesh,
                                                Workspace& workspace,
                                                const double gamma) const {
    const auto& U = layer.U();
    auto& W = workspace.W();

    for (std::size_t cell_id = 0; cell_id < mesh.GetCellCount(); ++cell_id) {
        ConservativeCell U_cell;
        U_cell.rho = U(cell_id, DataLayer::k_rho);
        U_cell.rhoU = U(cell_id, DataLayer::k_rhoU);
        U_cell.rhoV = U(cell_id, DataLayer::k_rhoV);
        U_cell.rhoW = U(cell_id, DataLayer::k_rhoW);
        U_cell.E = U(cell_id, DataLayer::k_E);

        const PrimitiveCell w = PrimitiveFromConservativeCell(U_cell, gamma);

        W(cell_id, Workspace::k_rho) = w.rho;
        W(cell_id, Workspace::k_u) = w.u;
        W(cell_id, Workspace::k_v) = w.v;
        W(cell_id, Workspace::k_w) = w.w;
        W(cell_id, Workspace::k_p) = w.P;
    }
}

FaceNormal GodunovSpatialOperator::BuildFaceNormal(const Face& face) const {
    FaceNormal normal;
    normal.x = face.normal_x;
    normal.y = face.normal_y;
    normal.z = face.normal_z;
    return normal;
}

void GodunovSpatialOperator::AccumulateFluxToOwner(const Mesh& mesh,
                                                   const Face& face,
                                                   const ConservativeCell& flux,
                                                   Workspace& workspace) const {
    auto& rhs = workspace.Rhs();
    const Cell& owner = mesh.GetCell(face.owner_cell_id);

    const double scale = face.measure / owner.volume;

    rhs(face.owner_cell_id, DataLayer::k_rho) -= flux.rho * scale;
    rhs(face.owner_cell_id, DataLayer::k_rhoU) -= flux.rhoU * scale;
    rhs(face.owner_cell_id, DataLayer::k_rhoV) -= flux.rhoV * scale;
    rhs(face.owner_cell_id, DataLayer::k_rhoW) -= flux.rhoW * scale;
    rhs(face.owner_cell_id, DataLayer::k_E) -= flux.E * scale;
}

void GodunovSpatialOperator::AccumulateFluxToNeighbor(const Mesh& mesh,
                                                      const Face& face,
                                                      const ConservativeCell& flux,
                                                      Workspace& workspace) const {
    auto& rhs = workspace.Rhs();
    const Cell& neighbor = mesh.GetCell(face.neighbor_cell_id);

    const double scale = face.measure / neighbor.volume;

    rhs(face.neighbor_cell_id, DataLayer::k_rho) += flux.rho * scale;
    rhs(face.neighbor_cell_id, DataLayer::k_rhoU) += flux.rhoU * scale;
    rhs(face.neighbor_cell_id, DataLayer::k_rhoV) += flux.rhoV * scale;
    rhs(face.neighbor_cell_id, DataLayer::k_rhoW) += flux.rhoW * scale;
    rhs(face.neighbor_cell_id, DataLayer::k_E) += flux.E * scale;
}

void GodunovSpatialOperator::AccumulateInternalFace(const DataLayer& layer,
                                                    const Mesh& mesh,
                                                    const Face& face,
                                                    Workspace& workspace,
                                                    const double gamma) const {
    PrimitiveCell owner_state;
    PrimitiveCell neighbor_state;

    reconstruction_->ReconstructInteriorFace(mesh, workspace, face, owner_state, neighbor_state);

    const FaceNormal normal = BuildFaceNormal(face);
    const ConservativeCell flux =
        riemann_solver_->ComputeFlux(owner_state, neighbor_state, gamma, normal);

    AccumulateFluxToOwner(mesh, face, flux, workspace);
    AccumulateFluxToNeighbor(mesh, face, flux, workspace);
}

void GodunovSpatialOperator::AccumulateBoundaryFace(const DataLayer& layer,
                                                    const Mesh& mesh,
                                                    const Face& face,
                                                    Workspace& workspace,
                                                    const double gamma) const {
    PrimitiveCell interior_state;
    reconstruction_->ReconstructBoundaryFaceInterior(mesh, workspace, face, interior_state);

    const PrimitiveCell exterior_state =
        boundary_manager_->BuildExteriorState(layer, mesh, face, interior_state);

    const FaceNormal normal = BuildFaceNormal(face);
    const ConservativeCell flux =
        riemann_solver_->ComputeFlux(interior_state, exterior_state, gamma, normal);

    AccumulateFluxToOwner(mesh, face, flux, workspace);
}

void GodunovSpatialOperator::ComputeRHS(const DataLayer& layer,
                                        const Mesh& mesh,
                                        Workspace& workspace,
                                        const double gamma,
                                        const double dt) const {
    (void)dt;

    workspace.ResizeFrom(mesh);
    workspace.ZeroRhs();

    FillPrimitiveCache(layer, mesh, workspace, gamma);

    for (const Face& face : mesh.Faces()) {
        if (face.IsInternal()) {
            AccumulateInternalFace(layer, mesh, face, workspace, gamma);
            continue;
        }

        if (face.IsBoundary()) {
            AccumulateBoundaryFace(layer, mesh, face, workspace, gamma);
            continue;
        }

        throw std::runtime_error(
            "GodunovSpatialOperator::ComputeRHS: face has invalid topology classification"
        );
    }
}
