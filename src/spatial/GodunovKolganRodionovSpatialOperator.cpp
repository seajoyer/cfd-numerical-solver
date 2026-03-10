#include "spatial/GodunovKolganRodionovSpatialOperator.hpp"

#include <algorithm>
#include <cctype>
#include <stdexcept>
#include <string>

#include "bc/BoundaryManager.hpp"
#include "data/DataLayer.hpp"
#include "data/Mesh.hpp"
#include "data/Variables.hpp"
#include "data/Workspace.hpp"
#include "reconstruction/ENOReconstruction.hpp"
#include "reconstruction/P0Reconstruction.hpp"
#include "reconstruction/P1Reconstruction.hpp"
#include "reconstruction/Reconstruction.hpp"
#include "reconstruction/WENOReconstruction.hpp"
#include "riemann/AcousticRiemannSolver.hpp"
#include "riemann/ExactIdealGasRiemannSolver.hpp"
#include "riemann/HLLCRiemannSolver.hpp"
#include "riemann/HLLRiemannSolver.hpp"
#include "riemann/OsherRiemannSolver.hpp"
#include "riemann/RiemannSolver.hpp"
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

    bool IsFluidInterface(const Mesh& mesh,
                          const int iL, const int jL, const int kL,
                          const int iR, const int jR, const int kR) {
        return mesh.IsFluidCell(iL, jL, kL) && mesh.IsFluidCell(iR, jR, kR);
    }
} // namespace

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
                                                      const Mesh& mesh,
                                                      Workspace& workspace,
                                                      const double gamma,
                                                      const double dt) const {
    if (!reconstruction_ || !riemann_solver_) {
        throw std::runtime_error("GKR: reconstruction_ or riemann_solver_ is null");
    }

    const int ng = mesh.GetPadding();
    if (ng < 1) {
        throw std::runtime_error("GKR: requires at least 1 ghost cell (ng>=1)");
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

    AccumulateAxis(layer, mesh, W, rhs, gamma, Axis::X, dt);
    if (mesh.GetDim() >= 2) {
        AccumulateAxis(layer, mesh, W, rhs, gamma, Axis::Y, dt);
    }
    if (mesh.GetDim() >= 3) {
        AccumulateAxis(layer, mesh, W, rhs, gamma, Axis::Z, dt);
    }

    if (viscosity_) {
        viscosity_->AddToRhs(layer, mesh, workspace.W(), gamma, dt, workspace.Rhs());
    }
}

double GodunovKolganRodionovSpatialOperator::InvMetricAt(
    const Mesh& mesh,
    const Axis axis,
    const int i,
    const int j,
    const int k) const {
    if (axis == Axis::X) {
        return mesh.InvDx()(static_cast<std::size_t>(i));
    }
    if (axis == Axis::Y) {
        return mesh.InvDy()(static_cast<std::size_t>(j));
    }
    return mesh.InvDz()(static_cast<std::size_t>(k));
}

void GodunovKolganRodionovSpatialOperator::MapCellIndices(
    const Axis axis,
    const int s,
    const int a,
    const int b,
    int& i,
    int& j,
    int& k) const {
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
    const int i,
    const int j,
    const int k,
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
                                                              const Mesh& mesh,
                                                              const int ci,
                                                              const int cj,
                                                              const int ck,
                                                              PrimitiveCell& WL_face,
                                                              PrimitiveCell& WR_face,
                                                              PrimitiveCell& WLm_face,
                                                              PrimitiveCell& WRm_face,
                                                              ConservativeCell& UL_star,
                                                              ConservativeCell& UR_star,
                                                              const double dt) const {
    const int im = ci - st.di;
    const int jm = cj - st.dj;
    const int km = ck - st.dk;

    const int ip = ci + st.di;
    const int jp = cj + st.dj;
    const int kp = ck + st.dk;

    if (!mesh.IsFluidCell(ci, cj, ck)) {
        UL_star = ConservativeCell{};
        UR_star = ConservativeCell{};
        return;
    }

    if (!mesh.IsFluidCell(im, jm, km) || !mesh.IsFluidCell(ip, jp, kp)) {
        const PrimitiveCell wc = LoadPrimitiveCell(W, ci, cj, ck);
        const ConservativeCell Uc = ConservativeFromPrimitive(wc, gamma);
        UL_star = Uc;
        UR_star = Uc;
        return;
    }

    reconstruction_->ReconstructFace(W, axis,
                                     ci - st.di, cj - st.dj, ck - st.dk,
                                     WLm_face, WRm_face);

    reconstruction_->ReconstructFace(W, axis,
                                     ci, cj, ck,
                                     WL_face, WR_face);

    const PrimitiveCell W_L_cell = WRm_face;
    const PrimitiveCell W_R_cell = WL_face;

    ConservativeCell U_L = ConservativeFromPrimitive(W_L_cell, gamma);
    ConservativeCell U_R = ConservativeFromPrimitive(W_R_cell, gamma);

    const FluxCell F_L = EulerFlux(W_L_cell, gamma, axis);
    const FluxCell F_R = EulerFlux(W_R_cell, gamma, axis);

    const double half_dt_over_d = 0.5 * dt * InvMetricAt(mesh, axis, ci, cj, ck);
    ApplyPredictor(U_L, U_R, F_L, F_R, half_dt_over_d);

    UL_star = U_L;
    UR_star = U_R;
}

void GodunovKolganRodionovSpatialOperator::AccumulateAxis(DataLayer& layer,
                                                          const Mesh& mesh,
                                                          const xt::xtensor<double, 4>& W,
                                                          xt::xtensor<double, 4>& rhs,
                                                          const double gamma,
                                                          const Axis axis,
                                                          const double dt) const {
    (void)layer;

    const AxisStride st = AxisStride::FromAxis(axis);

    const int i0 = mesh.GetCoreStartX();
    const int i1 = mesh.GetCoreEndExclusiveX();
    const int j0 = mesh.GetCoreStartY();
    const int j1 = mesh.GetCoreEndExclusiveY();
    const int k0 = mesh.GetCoreStartZ();
    const int k1 = mesh.GetCoreEndExclusiveZ();

    int a0 = 0;
    int a1 = 0;
    int b0 = 0;
    int b1 = 0;

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

    int s_core0 = 0;
    int s_core1 = 0;

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

    PrimitiveCell WL_face{};
    PrimitiveCell WR_face{};
    PrimitiveCell WLm_face{};
    PrimitiveCell WRm_face{};
    ConservativeCell UL_left{};
    ConservativeCell UR_left{};
    ConservativeCell UL_right{};
    ConservativeCell UR_right{};

    for (int b = b0; b < b1; ++b) {
        for (int a = a0; a < a1; ++a) {
            FluxCell Fm{};

            {
                int iL = 0;
                int jL = 0;
                int kL = 0;
                int iR = 0;
                int jR = 0;
                int kR = 0;

                MapCellIndices(axis, s_core0 - 1, a, b, iL, jL, kL);
                MapCellIndices(axis, s_core0, a, b, iR, jR, kR);

                if (IsFluidInterface(mesh, iL, jL, kL, iR, jR, kR)) {
                    ComputeStarForCell(W, gamma, axis, st, mesh,
                                       iL, jL, kL,
                                       WL_face, WR_face, WLm_face, WRm_face,
                                       UL_left, UR_left, dt);

                    ComputeStarForCell(W, gamma, axis, st, mesh,
                                       iR, jR, kR,
                                       WL_face, WR_face, WLm_face, WRm_face,
                                       UL_right, UR_right, dt);

                    Fm = riemann_solver_->ComputeFlux(
                        PrimitiveFromConservativeCell(UR_left, gamma),
                        PrimitiveFromConservativeCell(UL_right, gamma),
                        gamma, axis);
                }
            }

            for (int s = s_core0; s < s_core1; ++s) {
                FluxCell Fp{};

                int ci = 0;
                int cj = 0;
                int ck = 0;
                MapCellIndices(axis, s, a, b, ci, cj, ck);

                int iL = 0;
                int jL = 0;
                int kL = 0;
                int iR = 0;
                int jR = 0;
                int kR = 0;

                MapCellIndices(axis, s, a, b, iL, jL, kL);
                MapCellIndices(axis, s + 1, a, b, iR, jR, kR);

                if (IsFluidInterface(mesh, iL, jL, kL, iR, jR, kR)) {
                    ComputeStarForCell(W, gamma, axis, st, mesh,
                                       iL, jL, kL,
                                       WL_face, WR_face, WLm_face, WRm_face,
                                       UL_left, UR_left, dt);

                    ComputeStarForCell(W, gamma, axis, st, mesh,
                                       iR, jR, kR,
                                       WL_face, WR_face, WLm_face, WRm_face,
                                       UL_right, UR_right, dt);

                    Fp = riemann_solver_->ComputeFlux(
                        PrimitiveFromConservativeCell(UR_left, gamma),
                        PrimitiveFromConservativeCell(UL_right, gamma),
                        gamma, axis);
                }

                if (mesh.IsFluidCell(ci, cj, ck)) {
                    const double invd = InvMetricAt(mesh, axis, ci, cj, ck);
                    AccumulateCellRhs(rhs, ci, cj, ck, invd, Fp, Fm);
                }

                Fm = Fp;
            }
        }
    }
}
