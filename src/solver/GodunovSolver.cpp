#include "solver/GodunovSolver.hpp"

#include "reconstruction/ENOReconstruction.hpp"
#include "reconstruction/WENOReconstruction.hpp"
#include "reconstruction/P0Reconstruction.hpp"
#include "reconstruction/P1Reconstruction.hpp"
#include "riemann/AcousticRiemannSolver.hpp"
#include "riemann/ExactIdealGasRiemannSolver.hpp"
#include "riemann/HLLCRiemannSolver.hpp"
#include "riemann/HLLRiemannSolver.hpp"
#include "solver/EOS.hpp"
#include "solver/PositivityLimiter.hpp"
#include "solver/TimeStepCalculator.hpp"

GodunovSolver::GodunovSolver(const Settings& settings)
    : settings_(settings),
      boundary_manager_(settings.dim),
      rho_min_(1e-10),
      p_min_(1e-10) {
    cfl_ = settings_.cfl;
    InitializeReconstruction();
    InitializeRiemannSolver();
}

auto GodunovSolver::Step(DataLayer& layer, double& t_cur) -> double {
    boundary_manager_.ApplyAll(layer);

    const double dx = ComputeDx(layer);
    const int total_size = layer.GetTotalSize();
    const int core_start = layer.GetCoreStart(0);
    const int core_end = layer.GetCoreEndExclusive(0);
    const int n_core = core_end - core_start;

    if (n_core < 2 || total_size < 3) {
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

    const double dt_over_dx = dt / dx;

    const int n_interfaces = total_size - 1;

    xt::xarray<Primitive> left_states;
    xt::xarray<Primitive> right_states;
    reconstruction_->ReconstructStates(layer, left_states, right_states);

    xt::xarray<Flux> fluxes =
        xt::xarray<Flux>::from_shape({static_cast<std::size_t>(n_interfaces)});

    for (int i = 0; i < n_interfaces; ++i) {
        fluxes(i) = riemann_solver_->ComputeFlux(left_states(i), right_states(i),
                                                 settings_.gamma);
    }

    for (int j = core_start; j < core_end; ++j) {
        const Flux& Fminus = fluxes(j - 1);  // F_{j-1/2}
        const Flux& Fplus = fluxes(j);       // F_{j+1/2}

        Primitive w = layer.GetPrimitive(j);
        Conservative U = EOS::PrimToCons(w, settings_.gamma);

        U -= dt_over_dx * Flux::Diff(Fplus, Fminus);

        PositivityLimiter::Apply(U, settings_.gamma, rho_min_, p_min_);
        StoreConservativeCell(U, j, dx, layer);
    }

    t_cur += dt;
    return dt;
}

void GodunovSolver::SetCfl(double cfl) { cfl_ = cfl; }

void GodunovSolver::AddBoundary(int axis, std::shared_ptr<BoundaryCondition> left_bc,
                                std::shared_ptr<BoundaryCondition> right_bc) {
    boundary_manager_.Set(axis, std::move(left_bc), std::move(right_bc));
}

void GodunovSolver::InitializeReconstruction() {
    std::string name = settings_.reconstruction;
    if (name.find("p1") != std::string::npos) {
        reconstruction_ = std::make_shared<P1Reconstruction>();
    } else if (name.find("p0") != std::string::npos) {
        reconstruction_ = std::make_shared<P0Reconstruction>();
    } else if (name.starts_with("eno")) {
        int order = 3;
        try {
            order = std::stoi(name.substr(3, std::string::npos));
        } catch (...) {
            std::cout << "Order of ENO don't found. Set order to 3" << std::endl;
        }
        reconstruction_ = std::make_shared<ENOReconstruction>(order);
    } else if (name.starts_with("weno")) {
        int order = 5;
        try {
            order = std::stoi(name.substr(4, std::string::npos));
        } catch (...) {
            std::cout << "Order of WENO don't found. Set order to 5" << std::endl;
        }
        if (order != 3 and order != 5) {
            std::cout << "WENO supports only orders 3 or 5 for now. Set order to 5" <<
                std::endl;
        }
        reconstruction_ = std::make_shared<WENOReconstruction>(order);
    } else {
        reconstruction_ = std::make_shared<P0Reconstruction>();
    }
}

void GodunovSolver::InitializeRiemannSolver() {
    std::string name = settings_.riemann_solver;
    std::string lower(name.size(), '\0');
    std::transform(name.begin(), name.end(), lower.begin(), [](unsigned char c) -> char {
        return static_cast<char>(std::tolower(c));
    });

    if (lower.find("hllc") != std::string::npos) {
        riemann_solver_ = std::make_shared<HLLCRiemannSolver>();
    } else if (lower.find("hll") != std::string::npos) {
        riemann_solver_ = std::make_shared<HLLRiemannSolver>();
    } else if (lower.find("exact") != std::string::npos) {
        riemann_solver_ =
            std::make_shared<ExactIdealGasRiemannSolver>(0, settings_.Q_user);
    } else if (lower.find("acoustic") != std::string::npos) {
        riemann_solver_ = std::make_shared<AcousticRiemannSolver>();
    } else {
        riemann_solver_ = std::make_shared<HLLRiemannSolver>();
    }
}

auto GodunovSolver::ComputeDx(const DataLayer& layer) const -> double {
    if (settings_.N > 0 && settings_.L_x > 0.0) {
        return settings_.L_x / static_cast<double>(settings_.N);
    }

    const int core_start = layer.GetCoreStart(0);
    const int core_end = layer.GetCoreEndExclusive(0);
    if (core_end - core_start > 1) {
        return layer.xc(core_start + 1) - layer.xc(core_start);
    }

    return 1.0;
}

void GodunovSolver::StoreConservativeCell(const Conservative& uc, const int i,
                                          const double dx, DataLayer& layer) const {
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