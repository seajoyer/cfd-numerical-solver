#include "solver/GodunovSolver.hpp"

#include <algorithm>
#include <cctype>
#include <string>
#include <vector>

#include "solver/EOS.hpp"
#include "solver/Variables.hpp"
#include "reconstruction/Reconstruction.hpp"
#include "reconstruction/P0Reconstruction.hpp"
#include "riemann/RiemannSolver.hpp"
#include "riemann/HLLRiemannSolver.hpp"
#include "riemann/HLLCRiemannSolver.hpp"
#include "riemann/ExactIdealGasRiemannSolver.hpp"
#include "solver/TimeStepCalculator.hpp"
#include "solver/PositivityLimiter.hpp"

GodunovSolver::GodunovSolver(const Settings &settings)
    : settings_(settings),
      boundaryManager_(settings.dim),
      rhoMin_(1e-10),
      pMin_(1e-10) {
    cfl_ = settings_.cfl;
    InitializeReconstruction();
    InitializeRiemannSolver();
}

double GodunovSolver::Step(DataLayer &layer, double &t_cur) {
    boundaryManager_.ApplyAll(layer);

    const double dx = ComputeDx(layer);
    const int totalSize = layer.GetTotalSize();
    const int coreStart = layer.GetCoreStart(0);
    const int coreEnd = layer.GetCoreEndExclusive(0);

    if (coreEnd - coreStart < 2 || totalSize < 3) {
        return 0.0;
    }

    std::vector<Primitive> prim(static_cast<std::size_t>(totalSize));
    for (int i = 0; i < totalSize; ++i) {
        Primitive w;
        w.rho = layer.rho(i);
        w.u = layer.u(i);
        w.P = layer.P(i);
        prim[static_cast<std::size_t>(i)] = w;
    }

    std::vector<Primitive> primCore;
    primCore.reserve(static_cast<std::size_t>(coreEnd - coreStart));
    for (int i = coreStart; i < coreEnd; ++i) {
        primCore.push_back(prim[static_cast<std::size_t>(i)]);
    }

    double dt = TimeStepCalculator::ComputeDt(primCore, dx, cfl_, settings_.gamma);
    if (dt <= 0.0) {
        return 0.0;
    }

    if (t_cur + dt > settings_.t_end) {
        dt = settings_.t_end - t_cur;
        if (dt <= 0.0) {
            return 0.0;
        }
    }

    const int nInterfaces = totalSize - 1;
    std::vector<Primitive> leftStates(static_cast<std::size_t>(nInterfaces));
    std::vector<Primitive> rightStates(static_cast<std::size_t>(nInterfaces));
    reconstruction_->Reconstruct(prim, leftStates, rightStates);

    std::vector<Flux> fluxes(static_cast<std::size_t>(nInterfaces));
    for (int i = 0; i < nInterfaces; ++i) {
        fluxes[static_cast<std::size_t>(i)] =
            riemannSolver_->ComputeFlux(leftStates[static_cast<std::size_t>(i)],
                                        rightStates[static_cast<std::size_t>(i)],
                                        settings_.gamma);
    }

    std::vector<Conservative> cons(static_cast<std::size_t>(totalSize));
    for (int i = 0; i < totalSize; ++i) {
        const Primitive &w = prim[static_cast<std::size_t>(i)];

        Conservative u;
        u.rho = w.rho;
        u.rhoU = w.rho * w.u;
        u.E = w.P / (settings_.gamma - 1.0)
              + 0.5 * w.rho * w.u * w.u;

        cons[static_cast<std::size_t>(i)] = u;
    }

    std::vector<Conservative> updated = cons;
    const double dtOverDx = dt / dx;

    for (int j = coreStart; j < coreEnd; ++j) {
        const Flux &Fminus = fluxes[static_cast<std::size_t>(j - 1)];
        const Flux &Fplus = fluxes[static_cast<std::size_t>(j)];

        Conservative u = cons[static_cast<std::size_t>(j)];

        u.rho -= dtOverDx * (Fplus.mass - Fminus.mass);
        u.rhoU -= dtOverDx * (Fplus.momentum - Fminus.momentum);
        u.E -= dtOverDx * (Fplus.energy - Fminus.energy);

        PositivityLimiter::Apply(u, settings_.gamma, rhoMin_, pMin_);
        updated[static_cast<std::size_t>(j)] = u;
    }

    StoreConservativeArray(updated, layer);

    t_cur += dt;
    return dt;
}

void GodunovSolver::SetCfl(double cfl) {
    cfl_ = cfl;
}

void GodunovSolver::AddBoundary(int axis,
                                std::shared_ptr<BoundaryCondition> left_bc,
                                std::shared_ptr<BoundaryCondition> right_bc) {
    boundaryManager_.Set(axis, std::move(left_bc), std::move(right_bc));
}

void GodunovSolver::InitializeReconstruction() {
    std::string name = settings_.reconstruction;
    std::string lower(name.size(), '\0');
    std::transform(name.begin(), name.end(), lower.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    if (lower.find("p0") != std::string::npos) {
        reconstruction_ = std::make_shared<P0Reconstruction>();
    } else {
        reconstruction_ = std::make_shared<P0Reconstruction>();
    }

}

void GodunovSolver::InitializeRiemannSolver() {
    std::string name = settings_.riemann_solver;
    std::string lower(name.size(), '\0');
    std::transform(name.begin(), name.end(), lower.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });

    if (lower.find("hllc") != std::string::npos) {
        riemannSolver_ = std::make_shared<HLLCRiemannSolver>();
    } else if (lower.find("hll") != std::string::npos) {
        riemannSolver_ = std::make_shared<HLLRiemannSolver>();
    } else if (lower.find("exact") != std::string::npos) {
        riemannSolver_ = std::make_shared<ExactIdealGasRiemannSolver>();
    } else {
        riemannSolver_ = std::make_shared<HLLRiemannSolver>();
    }
}

double GodunovSolver::ComputeDx(const DataLayer &layer) const {
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

void GodunovSolver::BuildPrimitiveArray(const DataLayer &layer,
                                        std::vector<Primitive> &primitives) const {
    const int totalSize = layer.GetTotalSize();
    primitives.resize(static_cast<std::size_t>(totalSize));

    for (int i = 0; i < totalSize; ++i) {
        Primitive w;
        w.rho = layer.rho(i);
        w.u = layer.u(i);
        w.P = layer.P(i);
        primitives[static_cast<std::size_t>(i)] = w;
    }
}

void GodunovSolver::StoreConservativeArray(const std::vector<Conservative> &updatedConservative,
                                           DataLayer &layer) const {
    const int coreStart = layer.GetCoreStart(0);
    const int coreEnd = layer.GetCoreEndExclusive(0);
    const int totalSize = layer.GetTotalSize();

    if (static_cast<int>(updatedConservative.size()) < totalSize) {
        return;
    }

    const double dx = ComputeDx(layer);

    for (int i = coreStart; i < coreEnd; ++i) {
        const Conservative &uc = updatedConservative[static_cast<std::size_t>(i)];

        const double rho = uc.rho;
        const double uvel = uc.rhoU / rho;
        const double P = EOS::Pressure(uc, settings_.gamma);

        layer.rho(i) = rho;
        layer.u(i) = uvel;
        layer.P(i) = P;

        layer.p(i) = uc.rhoU;

        layer.V(i) = rho > 0.0 ? 1.0 / rho : 0.0;

        const double kinetic = 0.5 * rho * uvel * uvel;
        const double Eint = uc.E - kinetic;
        const double eint = rho > 0.0 ? Eint / rho : 0.0;
        const double etot = rho > 0.0 ? uc.E : 0.0;

        layer.U(i) = eint;
        layer.e(i) = etot;
        layer.m(i) = rho * dx;
    }
}
