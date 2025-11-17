#include "solver/GodunovKolganRodionovSolver.hpp"

#include <cmath>

GodunovKolganRodionovSolver::GodunovKolganRodionovSolver(const Settings& settings)
    : settings_(settings),
      boundaryManager_(settings.dim),
      rhoMin_(1e-10),
      pMin_(1e-10) {
    cfl_ = settings_.cfl;
    InitializeReconstruction();
    InitializeRiemannSolver();
}

void GodunovKolganRodionovSolver::AddBoundary(
    int axis,
    std::shared_ptr<BoundaryCondition> left_bc,
    std::shared_ptr<BoundaryCondition> right_bc) {
    boundaryManager_.Set(axis, std::move(left_bc), std::move(right_bc));
}

void GodunovKolganRodionovSolver::InitializeReconstruction() {
    reconstruction_ = std::make_shared<P1Reconstruction>();
}

void GodunovKolganRodionovSolver::InitializeRiemannSolver() {
    std::string name = settings_.riemann_solver;
    std::string lower(name.size(), '\0');

    std::transform(name.begin(), name.end(), lower.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

    if (lower.find("hllc") != std::string::npos) {
        riemannSolver_ = std::make_shared<HLLCRiemannSolver>();
    } else if (lower.find("hll") != std::string::npos) {
        riemannSolver_ = std::make_shared<HLLRiemannSolver>();
    } else if (lower.find("exact") != std::string::npos) {
        riemannSolver_ = std::make_shared<ExactIdealGasRiemannSolver>(0., settings_.Q);
    } else {
        riemannSolver_ = std::make_shared<HLLRiemannSolver>();
    }
}

double GodunovKolganRodionovSolver::ComputeDx(const DataLayer& layer) const {
    if (settings_.N > 0 && settings_.L_x > 0.0) {
        return settings_.L_x / static_cast<double>(settings_.N);
    }

    const int coreStart = layer.GetCoreStart(0);
    const int coreEnd = layer.GetCoreEndExclusive(0);

    if (coreEnd - coreStart > 1) {
        return layer.xc(coreStart + 1) - layer.xc(coreStart);
    }

    return 1.0;
}

double GodunovKolganRodionovSolver::Step(DataLayer& layer, double& t_cur) {
    // 1. Apply boundary conditions (update ghost cells).
    boundaryManager_.ApplyAll(layer);

    const double dx = ComputeDx(layer);
    const int totalSize = layer.GetTotalSize();
    const int coreStart = layer.GetCoreStart(0);
    const int coreEnd = layer.GetCoreEndExclusive(0);
    const int nCore = coreEnd - coreStart;

    if (nCore < 2 || totalSize < 3) {
        return 0.0;
    }

    // 2. Compute dt from CFL condition using primitive fields in DataLayer.
    double dt = TimeStepCalculator::ComputeDt(layer, dx, cfl_, settings_.gamma);
    if (dt <= 0.0) {
        return 0.0;
    }

    // Clamp dt to not exceed t_end.
    if (t_cur + dt > settings_.t_end) {
        dt = settings_.t_end - t_cur;
        if (dt <= 0.0) {
            return 0.0;
        }
    }

    const double halfDtOverDx = 0.5 * dt / dx;
    const double dtOverDx = dt / dx;

    // 3. Allocate buffer for updated conservative states (core cells).
    xt::xarray<Conservative> updated =
        xt::xarray<Conservative>::from_shape({static_cast<std::size_t>(totalSize)});

    // 4. Main loop: update physical cells using MUSCL–Hancock + Riemann solver.
    for (int j = coreStart; j < coreEnd; ++j) {
        // Conservative state in cell j at time n.
        Primitive w_j = layer.GetPrimitive(j);
        Conservative Uj = ToConservative(w_j, settings_.gamma);

        // Predicted states for j-1, j, j+1.
        Conservative UL_star_left, UR_star_left;
        Conservative UL_star_center, UR_star_center;
        Conservative UL_star_right, UR_star_right;

        const int iL = j - 1;
        const int iC = j;
        const int iR = j + 1;

        ComputePredictedStatesAtCell(layer, iL, halfDtOverDx,
                                     UL_star_left, UR_star_left);
        ComputePredictedStatesAtCell(layer, iC, halfDtOverDx,
                                     UL_star_center, UR_star_center);
        ComputePredictedStatesAtCell(layer, iR, halfDtOverDx,
                                     UL_star_right, UR_star_right);

        // Interface j-1/2: right from cell (j-1), left from cell j.
        Primitive WL_minus = ToPrimitive(UR_star_left, settings_.gamma);
        Primitive WR_minus = ToPrimitive(UL_star_center, settings_.gamma);
        Flux Fminus = riemannSolver_->ComputeFlux(WL_minus, WR_minus, settings_.gamma);

        // Interface j+1/2: right from cell j, left from cell (j+1).
        Primitive WL_plus = ToPrimitive(UR_star_center, settings_.gamma);
        Primitive WR_plus = ToPrimitive(UL_star_right, settings_.gamma);
        Flux Fplus = riemannSolver_->ComputeFlux(WL_plus, WR_plus, settings_.gamma);

        // Finite-volume update for cell j.
        Uj -= dtOverDx * Flux::Diff(Fplus, Fminus);

        // Positivity limiter on updated state.
        PositivityLimiter::Apply(Uj, settings_.gamma, rhoMin_, pMin_);

        // Store into buffer; DataLayer is not modified yet.
        updated(j) = Uj;
    }

    // 5. Write updated conservative states back to DataLayer (core cells).
    for (int j = coreStart; j < coreEnd; ++j) {
        const Conservative& Uj_new = updated(j);
        StoreConservativeCell(Uj_new, j, dx, layer);
    }

    // 6. Advance time.
    t_cur += dt;
    return dt;
}

void GodunovKolganRodionovSolver::StoreConservativeCell(const Conservative& uc,
                                                        const int i,
                                                        const double dx,
                                                        DataLayer& layer) const {
    const double rho = uc.rho;
    const double rhoU = uc.rhoU;

    const double uvel = (rho > 0.0) ? (rhoU / rho) : 0.0;
    const double P = EOS::Pressure(uc, settings_.gamma);

    layer.rho(i) = rho;
    layer.u(i) = uvel;
    layer.P(i) = P;

    layer.p(i) = rhoU;
    layer.V(i) = (rho > 0.0) ? (1.0 / rho) : 0.0;

    const double kinetic = 0.5 * rho * uvel * uvel;
    const double Eint = uc.E - kinetic;
    const double eint = (rho > 0.0) ? (Eint / rho) : 0.0;
    const double Etot = (rho > 0.0) ? uc.E : 0.0;

    layer.U(i) = eint;
    layer.e(i) = Etot;
    layer.m(i) = rho * dx;
}

void GodunovKolganRodionovSolver::ComputePredictedStatesAtCell(
    const DataLayer& layer,
    const int i,
    const double halfDtOverDx,
    Conservative& U_L_star_out,
    Conservative& U_R_star_out) const {
    const int totalSize = layer.GetTotalSize();
    const int nInterfaces = totalSize - 1;

    auto clampInterface = [nInterfaces](int idx) {
        if (idx < 0) {
            return 0;
        }
        if (idx >= nInterfaces) {
            return nInterfaces - 1;
        }
        return idx;
    };

    const int leftInterface = clampInterface(i - 1); // i-1/2
    const int rightInterface = clampInterface(i); // i+1/2

    Primitive WL_left, WR_left;
    Primitive WL_right, WR_right;

    // States at interface i-1/2: (WL_left, WR_left)
    reconstruction_->ComputeInterfaceStates(layer, leftInterface,
                                            WL_left, WR_left);
    // States at interface i+1/2: (WL_right, WR_right)
    reconstruction_->ComputeInterfaceStates(layer, rightInterface,
                                            WL_right, WR_right);

    // Cell-internal states at time n:
    //  - left state in cell i:  right state from interface i-1/2
    //  - right state in cell i: left state  from interface i+1/2
    Primitive W_L_cell = WR_left;
    Primitive W_R_cell = WL_right;

    // Convert to conservative form.
    Conservative U_L = ToConservative(W_L_cell, settings_.gamma);
    Conservative U_R = ToConservative(W_R_cell, settings_.gamma);

    // Euler fluxes from these states.
    Flux F_L = EulerFlux(W_L_cell, settings_.gamma);
    Flux F_R = EulerFlux(W_R_cell, settings_.gamma);

    // MUSCL–Hancock predictor to half time.
    U_L_star_out = U_L;
    U_R_star_out = U_R;

    U_L_star_out -= halfDtOverDx * Flux::Diff(F_R, F_L);
    U_R_star_out -= halfDtOverDx * Flux::Diff(F_R, F_L);
}