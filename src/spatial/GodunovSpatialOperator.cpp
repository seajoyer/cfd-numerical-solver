#include "spatial/GodunovSpatialOperator.hpp"

#include <algorithm>
#include <cctype>
#include <stdexcept>
#include <string>

#include "bc/BoundaryManager.hpp"
#include "data/Variables.hpp"

#include "reconstruction/P0Reconstruction.hpp"
#include "reconstruction/P1Reconstruction.hpp"
#include "reconstruction/ENOReconstruction.hpp"
#include "reconstruction/WENOReconstruction.hpp"

#include "riemann/AcousticRiemannSolver.hpp"
#include "riemann/RusanovRiemannSolver.hpp"
#include "riemann/HLLRiemannSolver.hpp"
#include "riemann/HLLCRiemannSolver.hpp"
#include "riemann/RoeRiemannSolver.hpp"
#include "riemann/OsherRiemannSolver.hpp"
#include "riemann/ExactIdealGasRiemannSolver.hpp"

#include "viscosity/VNRArtificialViscosity.hpp"

GodunovSpatialOperator::GodunovSpatialOperator(const Settings& settings,
                                               std::shared_ptr<BoundaryManager> boundary_manager)
    : SpatialOperator(std::move(boundary_manager)) {
    InitializeReconstruction(settings);
    InitializeRiemannSolver(settings);

    viscosity_ = nullptr;
    if (settings.viscosity) {
        viscosity_ = std::make_shared<VNRArtificialViscosity>(settings);
    }
}

void GodunovSpatialOperator::InitializeReconstruction(const Settings& settings) {
    std::string name = settings.reconstruction;

    if (name.find("p1") != std::string::npos) {
        reconstruction_ = std::make_shared<P1Reconstruction>();
    }
    else if (name.find("p0") != std::string::npos) {
        reconstruction_ = std::make_shared<P0Reconstruction>();
    }
    else if (name.rfind("eno", 0) == 0) {
        int order = 3;
        try {
            order = std::stoi(name.substr(3));
        }
        catch (...) {}
        reconstruction_ = std::make_shared<ENOReconstruction>(order);
    }
    else if (name.rfind("weno", 0) == 0) {
        int order = 5;
        try {
            order = std::stoi(name.substr(4));
        }
        catch (...) {}
        reconstruction_ = std::make_shared<WENOReconstruction>(order);
    }
    else {
        reconstruction_ = std::make_shared<P0Reconstruction>();
    }
}

void GodunovSpatialOperator::InitializeRiemannSolver(const Settings& settings) {
    std::string name = settings.riemann_solver;
    std::string lower(name.size(), '\0');

    std::transform(name.begin(), name.end(), lower.begin(),
                   [](unsigned char c) -> char {
                       return static_cast<char>(std::tolower(c));
                   });

    if (lower.find("rusanov") != std::string::npos) {
        riemann_solver_ = std::make_shared<RusanovRiemannSolver>();
    }
    else if (lower.find("hllc") != std::string::npos) {
        riemann_solver_ = std::make_shared<HLLCRiemannSolver>();
    }
    else if (lower.find("hll") != std::string::npos) {
        riemann_solver_ = std::make_shared<HLLRiemannSolver>();
    }
    else if (lower.find("roe") != std::string::npos) {
        riemann_solver_ = std::make_shared<RoeRiemannSolver>();
    }
    else if (lower.find("osher") != std::string::npos) {
        riemann_solver_ = std::make_shared<OsherRiemannSolver>();
    }
    else if (lower.find("exact") != std::string::npos) {
        riemann_solver_ = std::make_shared<ExactIdealGasRiemannSolver>(0.0, settings.Q_user);
    }
    else if (lower.find("acoustic") != std::string::npos) {
        riemann_solver_ = std::make_shared<AcousticRiemannSolver>();
    }
    else {
        riemann_solver_ = std::make_shared<HLLRiemannSolver>();
    }
}

void GodunovSpatialOperator::ComputeRHS(DataLayer& layer, Workspace& workspace, const double gamma,
                                        const double dt) const {
    if (!reconstruction_ || !riemann_solver_) {
        throw std::runtime_error("GodunovSpatialOperator: reconstruction_ or riemann_solver_ is null");
    }

    const int ng = layer.GetPadding();
    if (ng < 1) {
        throw std::runtime_error("GodunovSpatialOperator: requires at least 1 ghost cell (ng>=1)");
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

    AccumulateAxisFluxDivergence(layer, W, rhs, gamma, Axis::X);

    if (layer.GetDim() >= 2) {
        AccumulateAxisFluxDivergence(layer, W, rhs, gamma, Axis::Y);
    }
    if (layer.GetDim() >= 3) {
        AccumulateAxisFluxDivergence(layer, W, rhs, gamma, Axis::Z);
    }

    if (viscosity_) {
        viscosity_->AddToRhs(layer, workspace.W(), gamma, dt, workspace.Rhs());
    }
}

[[nodiscard]] const xt::xtensor<double, 1>&
GodunovSpatialOperator::InvMetric(const DataLayer& layer, const Axis axis) const {
    if (axis == Axis::X) return layer.InvDx();
    if (axis == Axis::Y) return layer.InvDy();
    return layer.InvDz();
}

void GodunovSpatialOperator::AccumulateAxisFluxDivergence(const DataLayer& layer,
                                                          const xt::xtensor<double, 4>& W,
                                                          xt::xtensor<double, 4>& rhs,
                                                          const double gamma,
                                                          const Axis axis) const {
    const AxisStride s = AxisStride::FromAxis(axis);
    const auto& inv_h = InvMetric(layer, axis);

    const int i0 = layer.GetCoreStartX();
    const int i1 = layer.GetCoreEndExclusiveX();
    const int j0 = layer.GetCoreStartY();
    const int j1 = layer.GetCoreEndExclusiveY();
    const int k0 = layer.GetCoreStartZ();
    const int k1 = layer.GetCoreEndExclusiveZ();

    PrimitiveCell WL;
    PrimitiveCell WR;
    FluxCell Fp;
    FluxCell Fm;

    for (int k = k0; k < k1; ++k) {
        for (int j = j0; j < j1; ++j) {
            for (int i = i0; i < i1; ++i) {
                // plus face: left cell = (i,j,k)
                reconstruction_->ReconstructFace(W, axis, i, j, k, WL, WR);
                Fp = riemann_solver_->ComputeFlux(WL, WR, gamma, axis);

                // minus face: left cell = (i,j,k) - stride
                reconstruction_->ReconstructFace(W, axis, i - s.di, j - s.dj, k - s.dk, WL, WR);
                Fm = riemann_solver_->ComputeFlux(WL, WR, gamma, axis);

                const int idx = (axis == Axis::X) ? i : (axis == Axis::Y) ? j : k;
                const double inv = inv_h(static_cast<std::size_t>(idx));

                rhs(DataLayer::k_rho, i, j, k) += -(Fp.mass - Fm.mass) * inv;
                rhs(DataLayer::k_rhoU, i, j, k) += -(Fp.mom_x - Fm.mom_x) * inv;
                rhs(DataLayer::k_rhoV, i, j, k) += -(Fp.mom_y - Fm.mom_y) * inv;
                rhs(DataLayer::k_rhoW, i, j, k) += -(Fp.mom_z - Fm.mom_z) * inv;
                rhs(DataLayer::k_E, i, j, k) += -(Fp.energy - Fm.energy) * inv;
            }
        }
    }
}
