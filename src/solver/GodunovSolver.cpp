#include "solver/GodunovSolver.hpp"


GodunovSolver::GodunovSolver(const Settings& settings)
    : settings_(settings),
      boundaryManager_(settings.dim),
      rhoMin_(1e-10),
      pMin_(1e-10) {
    cfl_ = settings_.cfl;
    InitializeReconstruction();
    InitializeRiemannSolver();
}

double GodunovSolver::Step(DataLayer& layer, double& t_cur) {
    boundaryManager_.ApplyAll(layer);

    const double dx = ComputeDx(layer);
    const int totalSize = layer.GetTotalSize();
    const int coreStart = layer.GetCoreStart(0);
    const int coreEnd = layer.GetCoreEndExclusive(0);
    const int nCore = coreEnd - coreStart;

    if (nCore < 2 || totalSize < 3) {
        return 0.0;
    }

    double dt = TimeStepCalculator::ComputeDt(layer, dx, cfl_, settings_.gamma);
    if (dt <= 0.0) {
        return 0.0;
    }

    if (t_cur + dt > settings_.t_end) {
        dt = settings_.t_end - t_cur;
        if (dt <= 0.0) {
            return 0.0;
        }
    }

    const double dtOverDx = dt / dx;

    const int nInterfaces = totalSize - 1;
    xt::xarray<Flux> fluxes =
        xt::xarray<Flux>::from_shape({static_cast<std::size_t>(nInterfaces)});

    for (int i = 0; i < nInterfaces; ++i) {
        Primitive WL, WR;
        reconstruction_->ComputeInterfaceStates(layer, i, WL, WR);
        fluxes(i) = riemannSolver_->ComputeFlux(WL, WR, settings_.gamma);
    }

    for (int j = coreStart; j < coreEnd; ++j) {
        const Flux& Fminus = fluxes(j - 1); // F_{j-1/2}
        const Flux& Fplus = fluxes(j); // F_{j+1/2}

        Primitive w = layer.GetPrimitive(j);
        Conservative U = ToConservative(w, settings_.gamma);

        U -= dtOverDx * Flux::Diff(Fplus, Fminus);

        PositivityLimiter::Apply(U, settings_.gamma, rhoMin_, pMin_);
        StoreConservativeCell(U, j, dx, layer);
    }

    t_cur += dt;
    return dt;
}


void GodunovSolver::SetCfl(double cfl) {
    cfl_ = cfl;
}

void GodunovSolver::AddBoundary(
    int axis,
    std::shared_ptr<BoundaryCondition> left_bc,
    std::shared_ptr<BoundaryCondition> right_bc) {
    boundaryManager_.Set(axis, std::move(left_bc), std::move(right_bc));
}

void GodunovSolver::InitializeReconstruction() {
    std::string name = settings_.reconstruction;
    std::string lower(name.size(), '\0');
    std::transform(name.begin(), name.end(), lower.begin(),
                   [](unsigned char c) {
                       return static_cast<char>(std::tolower(c));
                   });

    if (lower.find("p1") != std::string::npos) {
        reconstruction_ = std::make_shared<P1Reconstruction>();
    } else if (lower.find("p0") != std::string::npos) {
        reconstruction_ = std::make_shared<P0Reconstruction>();
    } else {
        reconstruction_ = std::make_shared<P0Reconstruction>();
    }
}

void GodunovSolver::InitializeRiemannSolver() {
    std::string name = settings_.riemann_solver;
    std::string lower(name.size(), '\0');
    std::transform(name.begin(), name.end(), lower.begin(),
                   [](unsigned char c) {
                       return static_cast<char>(std::tolower(c));
                   });

    if (lower.find("hllc") != std::string::npos) {
        riemannSolver_ = std::make_shared<HLLCRiemannSolver>();
    } else if (lower.find("hll") != std::string::npos) {
        riemannSolver_ = std::make_shared<HLLRiemannSolver>();
    } else if (lower.find("exact") != std::string::npos) {
        riemannSolver_ = std::make_shared<ExactIdealGasRiemannSolver>(0, settings_.Q);
    } else if (lower.find("acoustic") != std::string::npos) {
        riemannSolver_ = std::make_shared<AcousticRiemannSolver>();
    } else {
        riemannSolver_ = std::make_shared<HLLRiemannSolver>();
    }
}

double GodunovSolver::ComputeDx(const DataLayer& layer) const {
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

void GodunovSolver::StoreConservativeCell(const Conservative& uc,
                                          const int i,
                                          const double dx,
                                          DataLayer& layer) const {
    const double rho = uc.rho;
    const double rhoU = uc.rhoU;

    const double uvel = rho > 0.0 ? rhoU / rho : 0.0;
    const double P = EOS::Pressure(uc, settings_.gamma);

    layer.rho(i) = rho;
    layer.u(i) = uvel;
    layer.P(i) = P;

    layer.p(i) = rhoU;
    layer.V(i) = rho > 0.0 ? 1.0 / rho : 0.0;

    const double kinetic = 0.5 * rho * uvel * uvel;
    const double Eint = uc.E - kinetic;
    const double eint = rho > 0.0 ? Eint / rho : 0.0;
    const double Etot = rho > 0.0 ? uc.E : 0.0;

    layer.U(i) = eint;
    layer.e(i) = Etot;
    layer.m(i) = rho * dx;
}