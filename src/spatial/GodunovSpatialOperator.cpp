#include "spatial/GodunovSpatialOperator.hpp"

#include <algorithm>
#include <cctype>
#include <stdexcept>
#include <string>

#include "bc/BoundaryManager.hpp"
#include "bc/BoundaryCondition.hpp"
#include "bc/InternalBoundaryCondition.hpp"
#include "bc/SlipWallInternalBoundary.hpp"
#include "data/Variables.hpp"
#include "reconstruction/ENOReconstruction.hpp"
#include "reconstruction/P0Reconstruction.hpp"
#include "reconstruction/P1Reconstruction.hpp"
#include "reconstruction/WENOReconstruction.hpp"
#include "riemann/AcousticRiemannSolver.hpp"
#include "riemann/ExactIdealGasRiemannSolver.hpp"
#include "riemann/HLLCRiemannSolver.hpp"
#include "riemann/HLLRiemannSolver.hpp"
#include "riemann/OsherRiemannSolver.hpp"
#include "riemann/RoeRiemannSolver.hpp"
#include "riemann/RusanovRiemannSolver.hpp"
#include "viscosity/VNRArtificialViscosity.hpp"

namespace {

PrimitiveCell LoadPrimitiveCell(const xt::xtensor<double, 4>& W, const int i, const int j, const int k) {
    PrimitiveCell w;
    w.rho = W(var::u_rho, i, j, k);
    w.u = W(var::u_u, i, j, k);
    w.v = W(var::u_v, i, j, k);
    w.w = W(var::u_w, i, j, k);
    w.P = W(var::u_P, i, j, k);
    return w;
}

bool TryComputeFaceFlux(const Mesh& mesh,
                        const xt::xtensor<double, 4>& W,
                        const Reconstruction& reconstruction,
                        const RiemannSolver& riemann_solver,
                        const InternalBoundaryCondition& internal_bc,
                        const double gamma,
                        const Axis axis,
                        const int il,
                        const int jl,
                        const int kl,
                        FluxCell& flux) {
    const AxisStride s = AxisStride::FromAxis(axis);
    const int ir = il + s.di;
    const int jr = jl + s.dj;
    const int kr = kl + s.dk;

    const bool left_fluid = mesh.IsFluidCell(il, jl, kl);
    const bool right_fluid = mesh.IsFluidCell(ir, jr, kr);

    if (!left_fluid && !right_fluid) {
        return false;
    }

    PrimitiveCell WL;
    PrimitiveCell WR;

    if (left_fluid && right_fluid) {
        reconstruction.ReconstructFace(W, axis, il, jl, kl, WL, WR);
        flux = riemann_solver.ComputeFlux(WL, WR, gamma, axis);
        return true;
    }

    if (left_fluid && !right_fluid) {
        WL = LoadPrimitiveCell(W, il, jl, kl);

        const Side side = Side::Right;
        if (mesh.HasImmersedFace(axis, side, il, jl, kl)) {
            const ImmersedFaceInfo info = mesh.GetImmersedFaceInfo(axis, side, il, jl, kl);
            internal_bc.BuildBoundaryState(WL, info, WR);
        } else {
            WR = WL;
        }

        flux = riemann_solver.ComputeFlux(WL, WR, gamma, axis);
        return true;
    }

    WR = LoadPrimitiveCell(W, ir, jr, kr);

    const Side side = Side::Left;
    if (mesh.HasImmersedFace(axis, side, ir, jr, kr)) {
        const ImmersedFaceInfo info = mesh.GetImmersedFaceInfo(axis, side, ir, jr, kr);
        internal_bc.BuildBoundaryState(WR, info, WL);
    } else {
        WL = WR;
    }

    flux = riemann_solver.ComputeFlux(WL, WR, gamma, axis);
    return true;
}

}  // namespace

GodunovSpatialOperator::GodunovSpatialOperator(const Settings& settings,
                                               std::shared_ptr<BoundaryManager> boundary_manager)
    : SpatialOperator(std::move(boundary_manager)) {
    InitializeReconstruction(settings);
    InitializeRiemannSolver(settings);

    viscosity_ = nullptr;
    if (settings.viscosity) {
        viscosity_ = std::make_shared<VNRArtificialViscosity>(settings);
    }

    internal_boundary_condition_ = std::make_shared<SlipWallInternalBoundary>();
}

void GodunovSpatialOperator::InitializeReconstruction(const Settings& settings) {
    std::string name = settings.reconstruction;

    if (name.find("p1") != std::string::npos) {
        reconstruction_ = std::make_shared<P1Reconstruction>();
    } else if (name.find("p0") != std::string::npos) {
        reconstruction_ = std::make_shared<P0Reconstruction>();
    } else if (name.rfind("eno", 0) == 0) {
        int order = 3;
        try {
            order = std::stoi(name.substr(3));
        } catch (...) {
        }
        reconstruction_ = std::make_shared<ENOReconstruction>(order);
    } else if (name.rfind("weno", 0) == 0) {
        int order = 5;
        try {
            order = std::stoi(name.substr(4));
        } catch (...) {
        }
        reconstruction_ = std::make_shared<WENOReconstruction>(order);
    } else {
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
    } else if (lower.find("hllc") != std::string::npos) {
        riemann_solver_ = std::make_shared<HLLCRiemannSolver>();
    } else if (lower.find("hll") != std::string::npos) {
        riemann_solver_ = std::make_shared<HLLRiemannSolver>();
    } else if (lower.find("roe") != std::string::npos) {
        riemann_solver_ = std::make_shared<RoeRiemannSolver>();
    } else if (lower.find("osher") != std::string::npos) {
        riemann_solver_ = std::make_shared<OsherRiemannSolver>();
    } else if (lower.find("exact") != std::string::npos) {
        riemann_solver_ = std::make_shared<ExactIdealGasRiemannSolver>(0.0, settings.Q_user);
    } else if (lower.find("acoustic") != std::string::npos) {
        riemann_solver_ = std::make_shared<AcousticRiemannSolver>();
    } else {
        riemann_solver_ = std::make_shared<HLLRiemannSolver>();
    }
}

void GodunovSpatialOperator::ComputeRHS(DataLayer& layer,
                                        const Mesh& mesh,
                                        Workspace& workspace,
                                        const double gamma,
                                        const double dt) const {
    if (!reconstruction_ || !riemann_solver_ || !internal_boundary_condition_) {
        throw std::runtime_error("GodunovSpatialOperator: required component is null");
    }

    const int ng = mesh.GetPadding();
    if (ng < 1) {
        throw std::runtime_error("GodunovSpatialOperator: requires at least 1 ghost cell (ng>=1)");
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

    AccumulateAxisFluxDivergence(layer, mesh, W, rhs, gamma, Axis::X);

    if (mesh.GetDim() >= 2) {
        AccumulateAxisFluxDivergence(layer, mesh, W, rhs, gamma, Axis::Y);
    }
    if (mesh.GetDim() >= 3) {
        AccumulateAxisFluxDivergence(layer, mesh, W, rhs, gamma, Axis::Z);
    }

    if (viscosity_) {
        viscosity_->AddToRhs(layer, mesh, workspace.W(), gamma, dt, workspace.Rhs());
    }
}

const xt::xtensor<double, 1>& GodunovSpatialOperator::InvMetric(const Mesh& mesh, const Axis axis) const {
    if (axis == Axis::X) {
        return mesh.InvDx();
    }
    if (axis == Axis::Y) {
        return mesh.InvDy();
    }
    return mesh.InvDz();
}

void GodunovSpatialOperator::AccumulateAxisFluxDivergence(const DataLayer& layer,
                                                          const Mesh& mesh,
                                                          const xt::xtensor<double, 4>& W,
                                                          xt::xtensor<double, 4>& rhs,
                                                          const double gamma,
                                                          const Axis axis) const {
    const AxisStride s = AxisStride::FromAxis(axis);
    const auto& inv_h = InvMetric(mesh, axis);

    const int i0 = mesh.GetCoreStartX();
    const int i1 = mesh.GetCoreEndExclusiveX();
    const int j0 = mesh.GetCoreStartY();
    const int j1 = mesh.GetCoreEndExclusiveY();
    const int k0 = mesh.GetCoreStartZ();
    const int k1 = mesh.GetCoreEndExclusiveZ();

    FluxCell Fp{};
    FluxCell Fm{};

    for (int k = k0; k < k1; ++k) {
        for (int j = j0; j < j1; ++j) {
            for (int i = i0; i < i1; ++i) {
                if (!mesh.IsFluidCell(i, j, k)) {
                    continue;
                }

                const bool ok_plus = TryComputeFaceFlux(mesh,
                                                        W,
                                                        *reconstruction_,
                                                        *riemann_solver_,
                                                        *internal_boundary_condition_,
                                                        gamma,
                                                        axis,
                                                        i,
                                                        j,
                                                        k,
                                                        Fp);

                const bool ok_minus = TryComputeFaceFlux(mesh,
                                                         W,
                                                         *reconstruction_,
                                                         *riemann_solver_,
                                                         *internal_boundary_condition_,
                                                         gamma,
                                                         axis,
                                                         i - s.di,
                                                         j - s.dj,
                                                         k - s.dk,
                                                         Fm);

                if (!ok_plus && !ok_minus) {
                    continue;
                }
                if (!ok_plus) {
                    Fp = FluxCell{};
                }
                if (!ok_minus) {
                    Fm = FluxCell{};
                }

                const int idx = axis == Axis::X ? i : (axis == Axis::Y ? j : k);
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