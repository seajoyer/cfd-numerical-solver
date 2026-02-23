#include "spatial/GodunovKolganRodionovSpatialOperator.hpp"

#include <stdexcept>

#include "bc/BoundaryManager.hpp"
#include "data/DataLayer.hpp"
#include "data/Variables.hpp"
#include "data/Workspace.hpp"

#include "reconstruction/Reconstruction.hpp"
#include "riemann/RiemannSolver.hpp"

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

GodunovKolganRodionovSpatialOperator::GodunovKolganRodionovSpatialOperator(
    const Settings& settings,
    std::shared_ptr<BoundaryManager> boundary_manager)
    : SpatialOperator(std::move(boundary_manager)) {
    InitializeReconstruction(settings);
    InitializeRiemannSolver(settings);

    viscosity_ = nullptr;
    if (settings.viscosity) {
        viscosity_ = std::make_shared<VNRArtificialViscosity>(settings);
    }
}

void GodunovKolganRodionovSpatialOperator::InitializeReconstruction(const Settings& settings) {
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

void GodunovKolganRodionovSpatialOperator::InitializeRiemannSolver(const Settings& settings) {
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

void GodunovKolganRodionovSpatialOperator::ComputeRHS(DataLayer& layer,
                                                      Workspace& workspace,
                                                      const double gamma, const double dt) const {
    if (!reconstruction_ || !riemann_solver_) {
        throw std::runtime_error("GKR: reconstruction_ or riemann_solver_ is null");
    }

    const int ng = layer.GetPadding();
    if (ng < 1) {
        throw std::runtime_error("GKR: requires at least 1 ghost cell (ng>=1)");
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

    AccumulateAxis(layer, W, rhs, gamma, Axis::X, dt);
    if (layer.GetDim() >= 2) AccumulateAxis(layer, W, rhs, gamma, Axis::Y, dt);
    if (layer.GetDim() >= 3) AccumulateAxis(layer, W, rhs, gamma, Axis::Z, dt);

    if (viscosity_) {
        viscosity_->AddToRhs(layer, workspace.W(), gamma, dt, workspace.Rhs());
    }
}


[[nodiscard]] double GodunovKolganRodionovSpatialOperator::InvMetricAt(
    const DataLayer& layer, const Axis axis, const int i, const int j, const int k) const {
    if (axis == Axis::X) return layer.InvDx()(static_cast<std::size_t>(i));
    if (axis == Axis::Y) return layer.InvDy()(static_cast<std::size_t>(j));
    return layer.InvDz()(static_cast<std::size_t>(k));
}

void GodunovKolganRodionovSpatialOperator::MapCellIndices(
    const Axis axis, const int s, const int a, const int b, int& i, int& j, int& k) const {
    if (axis == Axis::X) {
        i = s;
        j = a;
        k = b;
        return;
    }
    if (axis == Axis::Y) {
        i = a;
        j = s;
        k = b;
        return;
    }
    // Axis::Z
    i = a;
    j = b;
    k = s;
}

void GodunovKolganRodionovSpatialOperator::ApplyPredictor(
    ConservativeCell& UL,
    ConservativeCell& UR,
    const FluxCell& FL,
    const FluxCell& FR,
    const double half_dt_over_d) {
    const double d_mass = FR.mass - FL.mass;
    const double d_mom_x = FR.mom_x - FL.mom_x;
    const double d_mom_y = FR.mom_y - FL.mom_y;
    const double d_mom_z = FR.mom_z - FL.mom_z;
    const double d_energy = FR.energy - FL.energy;

    UL.rho -= half_dt_over_d * d_mass;
    UL.rhoU -= half_dt_over_d * d_mom_x;
    UL.rhoV -= half_dt_over_d * d_mom_y;
    UL.rhoW -= half_dt_over_d * d_mom_z;
    UL.E -= half_dt_over_d * d_energy;

    UR.rho -= half_dt_over_d * d_mass;
    UR.rhoU -= half_dt_over_d * d_mom_x;
    UR.rhoV -= half_dt_over_d * d_mom_y;
    UR.rhoW -= half_dt_over_d * d_mom_z;
    UR.E -= half_dt_over_d * d_energy;
}

void GodunovKolganRodionovSpatialOperator::AccumulateCellRhs(
    xt::xtensor<double, 4>& rhs,
    const int i, const int j, const int k,
    const double invd,
    const FluxCell& Fp,
    const FluxCell& Fm) {
    rhs(DataLayer::k_rho, i, j, k) += -(Fp.mass - Fm.mass) * invd;
    rhs(DataLayer::k_rhoU, i, j, k) += -(Fp.mom_x - Fm.mom_x) * invd;
    rhs(DataLayer::k_rhoV, i, j, k) += -(Fp.mom_y - Fm.mom_y) * invd;
    rhs(DataLayer::k_rhoW, i, j, k) += -(Fp.mom_z - Fm.mom_z) * invd;
    rhs(DataLayer::k_E, i, j, k) += -(Fp.energy - Fm.energy) * invd;
}

void GodunovKolganRodionovSpatialOperator::ComputeStarForCell(const xt::xtensor<double, 4>& W,
                                                              const double gamma,
                                                              const Axis axis,
                                                              const AxisStride& st,
                                                              const DataLayer& layer,
                                                              const int ci, const int cj, const int ck,
                                                              PrimitiveCell& WL_face,
                                                              PrimitiveCell& WR_face,
                                                              PrimitiveCell& WLm_face,
                                                              PrimitiveCell& WRm_face,
                                                              ConservativeCell& UL_star,
                                                              ConservativeCell& UR_star,
                                                              double dt) const {
    // Left face of cell C: interface (C-1/2) has left cell index = C - stride
    reconstruction_->ReconstructFace(W, axis,
                                     ci - st.di, cj - st.dj, ck - st.dk,
                                     WLm_face, WRm_face);

    // Right face of cell C: interface (C+1/2) has left cell index = C
    reconstruction_->ReconstructFace(W, axis,
                                     ci, cj, ck,
                                     WL_face, WR_face);

    // States *inside cell C* at its left/right faces
    const PrimitiveCell W_L_cell = WRm_face; // right state at (C-1/2)
    const PrimitiveCell W_R_cell = WL_face; // left  state at (C+1/2)

    ConservativeCell U_L = ConservativeFromPrimitive(W_L_cell, gamma);
    ConservativeCell U_R = ConservativeFromPrimitive(W_R_cell, gamma);

    const FluxCell F_L = EulerFlux(W_L_cell, gamma, axis);
    const FluxCell F_R = EulerFlux(W_R_cell, gamma, axis);

    const double half_dt_over_d = 0.5 * dt * InvMetricAt(layer, axis, ci, cj, ck);
    ApplyPredictor(U_L, U_R, F_L, F_R, half_dt_over_d);

    UL_star = U_L;
    UR_star = U_R;
}

void GodunovKolganRodionovSpatialOperator::AccumulateAxis(DataLayer& layer,
                                                          const xt::xtensor<double, 4>& W,
                                                          xt::xtensor<double, 4>& rhs,
                                                          const double gamma,
                                                          const Axis axis,
                                                          double dt) const {
    const AxisStride st = AxisStride::FromAxis(axis);

    const int i0 = layer.GetCoreStartX();
    const int i1 = layer.GetCoreEndExclusiveX();
    const int j0 = layer.GetCoreStartY();
    const int j1 = layer.GetCoreEndExclusiveY();
    const int k0 = layer.GetCoreStartZ();
    const int k1 = layer.GetCoreEndExclusiveZ();

    // transverse ranges (a,b) and sweep range (s)
    int a0 = 0, a1 = 0, b0 = 0, b1 = 0;
    if (axis == Axis::X) {
        a0 = j0;
        a1 = j1;
        b0 = k0;
        b1 = k1;
    }
    if (axis == Axis::Y) {
        a0 = i0;
        a1 = i1;
        b0 = k0;
        b1 = k1;
    }
    if (axis == Axis::Z) {
        a0 = i0;
        a1 = i1;
        b0 = j0;
        b1 = j1;
    }

    int s_core0 = 0, s_core1 = 0;
    if (axis == Axis::X) {
        s_core0 = i0;
        s_core1 = i1;
    }
    if (axis == Axis::Y) {
        s_core0 = j0;
        s_core1 = j1;
    }
    if (axis == Axis::Z) {
        s_core0 = k0;
        s_core1 = k1;
    }

    // We need F_{s_core0-1/2} as "Fm" for the first core cell.
    const int s_left = s_core0 - 1;

    PrimitiveCell WL_face{}, WR_face{}, WLm_face{}, WRm_face{};
    ConservativeCell ULs{}, URs{}, ULn{}, URn{};

    for (int b = b0; b < b1; ++b) {
        for (int a = a0; a < a1; ++a) {
            int iL, jL, kL;
            int iR, jR, kR;

            // Left cell at s_left, right cell at s_left+1 (needed for initial Fm)
            MapCellIndices(axis, s_left, a, b, iL, jL, kL);
            MapCellIndices(axis, s_left + 1, a, b, iR, jR, kR);

            ComputeStarForCell(W, gamma, axis, st, layer,
                               iL, jL, kL,
                               WL_face, WR_face, WLm_face, WRm_face,
                               ULs, URs, dt);

            ComputeStarForCell(W, gamma, axis, st, layer,
                               iR, jR, kR,
                               WL_face, WR_face, WLm_face, WRm_face,
                               ULn, URn, dt);

            // Interface between cell(s_left) and cell(s_left+1):
            // left state = UR* of left cell, right state = UL* of right cell
            FluxCell Fm = riemann_solver_->ComputeFlux(
                PrimitiveFromConservativeCell(URs, gamma),
                PrimitiveFromConservativeCell(ULn, gamma),
                gamma, axis);

            for (int s = s_core0; s < s_core1; ++s) {
                // shift: previous "right cell" becomes current "left cell"
                ULs = ULn;
                URs = URn;

                // build next right cell at s+1
                MapCellIndices(axis, s + 1, a, b, iR, jR, kR);
                ComputeStarForCell(W, gamma, axis, st, layer,
                                   iR, jR, kR,
                                   WL_face, WR_face, WLm_face, WRm_face,
                                   ULn, URn, dt);

                const FluxCell Fp = riemann_solver_->ComputeFlux(
                    PrimitiveFromConservativeCell(URs, gamma),
                    PrimitiveFromConservativeCell(ULn, gamma),
                    gamma, axis);

                int ci, cj, ck;
                MapCellIndices(axis, s, a, b, ci, cj, ck);

                const double invd = InvMetricAt(layer, axis, ci, cj, ck);
                AccumulateCellRhs(rhs, ci, cj, ck, invd, Fp, Fm);

                Fm = Fp;
            }
        }
    }
}
