#include "spatial/GodunovKolganRodionovSpatialOperator.hpp"

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
#include "parallel/StateSynchronizer.hpp"
#include "reconstruction/P0Reconstruction.hpp"
#include "reconstruction/P1Reconstruction.hpp"
#include "reconstruction/Reconstruction.hpp"
#include "riemann/ExactIdealGasRiemannSolver.hpp"
#include "riemann/HLLCRiemannSolver.hpp"
#include "riemann/HLLRiemannSolver.hpp"
#include "riemann/RiemannSolver.hpp"
#include "riemann/RoeRiemannSolver.hpp"
#include "riemann/RusanovRiemannSolver.hpp"

GodunovKolganRodionovSpatialOperator::GodunovKolganRodionovSpatialOperator(
    const Settings& settings,
    std::shared_ptr<BoundaryManager> boundary_manager,
    const StateSynchronizer* synchronizer
) : SpatialOperator(std::move(boundary_manager), synchronizer) {
    InitializeReconstruction(settings);
    InitializeRiemannSolver(settings);

    if (!boundary_manager_) {
        throw std::runtime_error(
            "GodunovKolganRodionovSpatialOperator: boundary_manager is null"
        );
    }
    if (!reconstruction_) {
        throw std::runtime_error(
            "GodunovKolganRodionovSpatialOperator: reconstruction is null"
        );
    }
    if (!riemann_solver_) {
        throw std::runtime_error(
            "GodunovKolganRodionovSpatialOperator: riemann_solver is null"
        );
    }
}

void GodunovKolganRodionovSpatialOperator::InitializeReconstruction(const Settings& settings) {
    std::string name = settings.reconstruction;
    for (char& c : name) {
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }

    if (name == "p1") {
        reconstruction_ = std::make_shared<P1Reconstruction>();
        return;
    }

    if (name == "p0") {
        reconstruction_ = std::make_shared<P0Reconstruction>();
        return;
    }

    throw std::runtime_error(
        "GodunovKolganRodionovSpatialOperator::InitializeReconstruction: unsupported reconstruction '" +
        settings.reconstruction + "'"
    );
}

void GodunovKolganRodionovSpatialOperator::InitializeRiemannSolver(const Settings& settings) {
    std::string name = settings.riemann_solver;
    for (char& c : name) {
        c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    }

    if (name == "rusanov") {
        riemann_solver_ = std::make_shared<RusanovRiemannSolver>();
        return;
    }

    if (name == "hll") {
        riemann_solver_ = std::make_shared<HLLRiemannSolver>();
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
        "GodunovKolganRodionovSpatialOperator::InitializeRiemannSolver: unsupported riemann solver '" +
        settings.riemann_solver + "'"
    );
}

void GodunovKolganRodionovSpatialOperator::FillPrimitiveCacheFromConservative(const xt::xtensor<double, 2>& U,
                                                                              const Mesh& mesh,
                                                                              Workspace& workspace,
                                                                              const double gamma) const {
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

PrimitiveCell GodunovKolganRodionovSpatialOperator::LoadCellPrimitive(const Workspace& workspace,
                                                                      const std::size_t cell_id) const {
    const auto& W = workspace.W();

    PrimitiveCell state;
    state.rho = W(cell_id, Workspace::k_rho);
    state.u = W(cell_id, Workspace::k_u);
    state.v = W(cell_id, Workspace::k_v);
    state.w = W(cell_id, Workspace::k_w);
    state.P = W(cell_id, Workspace::k_p);

    return state;
}

FaceNormal GodunovKolganRodionovSpatialOperator::BuildFaceNormal(const Face& face) const {
    FaceNormal normal;
    normal.x = face.normal_x;
    normal.y = face.normal_y;
    normal.z = face.normal_z;
    return normal;
}

void GodunovKolganRodionovSpatialOperator::AccumulateFluxToOwner(const Mesh& mesh,
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

void GodunovKolganRodionovSpatialOperator::AccumulateFluxToNeighbor(const Mesh& mesh,
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

void GodunovKolganRodionovSpatialOperator::AccumulatePredictorInternalFace(
    const DataLayer& layer,
    const Mesh& mesh,
    const Face& face,
    Workspace& workspace,
    const double gamma
) const {
    const PrimitiveCell owner_state = LoadCellPrimitive(workspace, face.owner_cell_id);
    const PrimitiveCell neighbor_state = LoadCellPrimitive(workspace, face.neighbor_cell_id);

    const FaceNormal normal = BuildFaceNormal(face);
    const ConservativeCell flux =
        riemann_solver_->ComputeFlux(owner_state, neighbor_state, gamma, normal);

    AccumulateFluxToOwner(mesh, face, flux, workspace);
    AccumulateFluxToNeighbor(mesh, face, flux, workspace);
}

void GodunovKolganRodionovSpatialOperator::AccumulatePredictorBoundaryFace(
    const DataLayer& layer,
    const Mesh& mesh,
    const Face& face,
    Workspace& workspace,
    const double gamma
) const {
    const PrimitiveCell interior_state = LoadCellPrimitive(workspace, face.owner_cell_id);

    const PrimitiveCell exterior_state =
        boundary_manager_->BuildExteriorState(layer, mesh, face, interior_state);

    const FaceNormal normal = BuildFaceNormal(face);
    const ConservativeCell flux =
        riemann_solver_->ComputeFlux(interior_state, exterior_state, gamma, normal);

    AccumulateFluxToOwner(mesh, face, flux, workspace);
}

void GodunovKolganRodionovSpatialOperator::ComputePredictorRhs(const DataLayer& layer,
                                                               const Mesh& mesh,
                                                               Workspace& workspace,
                                                               const double gamma) const {
    workspace.ZeroRhs();

    for (const Face& face : mesh.Faces()) {
        if (face.IsInternal() || face.IsMPIBoundary()) {
            AccumulatePredictorInternalFace(layer, mesh, face, workspace, gamma);
            continue;
        }

        if (face.IsPhysicalBoundary()) {
            AccumulatePredictorBoundaryFace(layer, mesh, face, workspace, gamma);
            continue;
        }

        throw std::runtime_error(
            "GodunovKolganRodionovSpatialOperator::ComputePredictorRhs: invalid face topology"
        );
    }
}

void GodunovKolganRodionovSpatialOperator::AccumulateFinalInternalFace(
    const DataLayer& layer,
    const Mesh& mesh,
    const Face& face,
    Workspace& workspace,
    const double gamma
) const {
    PrimitiveCell owner_state;
    PrimitiveCell neighbor_state;

    reconstruction_->ReconstructInteriorFace(mesh, workspace, face, owner_state, neighbor_state);

    const FaceNormal normal = BuildFaceNormal(face);
    const ConservativeCell flux = riemann_solver_->ComputeFlux(owner_state, neighbor_state, gamma, normal);

    AccumulateFluxToOwner(mesh, face, flux, workspace);
    AccumulateFluxToNeighbor(mesh, face, flux, workspace);
}

void GodunovKolganRodionovSpatialOperator::AccumulateFinalBoundaryFace(
    const DataLayer& layer,
    const Mesh& mesh,
    const Face& face,
    Workspace& workspace,
    const double gamma
) const {
    PrimitiveCell interior_state;
    reconstruction_->ReconstructBoundaryFaceInterior(mesh, workspace, face, interior_state);

    const PrimitiveCell exterior_state =
        boundary_manager_->BuildExteriorState(layer, mesh, face, interior_state);

    const FaceNormal normal = BuildFaceNormal(face);
    const ConservativeCell flux =
        riemann_solver_->ComputeFlux(interior_state, exterior_state, gamma, normal);

    AccumulateFluxToOwner(mesh, face, flux, workspace);
}

void GodunovKolganRodionovSpatialOperator::ComputeFinalRhs(const DataLayer& layer,
                                                           const Mesh& mesh,
                                                           Workspace& workspace,
                                                           const double gamma) const {
    workspace.ZeroRhs();

    for (const Face& face : mesh.Faces()) {
        if (face.IsInternal() || face.IsMPIBoundary()) {
            AccumulateFinalInternalFace(layer, mesh, face, workspace, gamma);
            continue;
        }

        if (face.IsPhysicalBoundary()) {
            AccumulateFinalBoundaryFace(layer, mesh, face, workspace, gamma);
            continue;
        }

        throw std::runtime_error(
            "GodunovKolganRodionovSpatialOperator::ComputeFinalRhs: invalid face topology"
        );
    }
}

void GodunovKolganRodionovSpatialOperator::ComputeRHS(const DataLayer& layer,
                                                      const Mesh& mesh,
                                                      Workspace& workspace,
                                                      const double gamma,
                                                      const double dt) const {
    if (dt <= 0.0) {
        throw std::runtime_error(
            "GodunovKolganRodionovSpatialOperator::ComputeRHS: dt must be positive"
        );
    }

    workspace.ResizeFrom(mesh);

    // Step 1: primitive cache from U^n
    FillPrimitiveCacheFromConservative(layer.U(), mesh, workspace, gamma);

    // Step 2: predictor RHS from first-order face fluxes
    ComputePredictorRhs(layer, mesh, workspace, gamma);

    const auto& rhs_predictor = workspace.Rhs();

    DataLayer half_layer;
    half_layer.Resize(mesh.GetCellCount());

    half_layer.U() = layer.U();
    half_layer.ReactantMassFraction() = layer.ReactantMassFraction();

    for (std::size_t cell_id = 0; cell_id < mesh.GetOwnedCellCount(); ++cell_id) {
        for (std::size_t var = 0; var < DataLayer::k_nvar; ++var) {
            half_layer.U()(cell_id, var) += 0.5 * dt * rhs_predictor(cell_id, var);
        }
    }

    if (synchronizer_) {
        synchronizer_->Synchronize(half_layer);
    }

    // Step 3: primitive cache from synchronized U^{n+1/2}
    FillPrimitiveCacheFromConservative(half_layer.U(), mesh, workspace, gamma);

    // Step 4: final RHS with selected reconstruction on predicted state
    ComputeFinalRhs(layer, mesh, workspace, gamma);
}
