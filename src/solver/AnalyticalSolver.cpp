#include "solver/AnalyticalSolver.hpp"

#include <algorithm>
#include <cmath>

AnalyticalSolver::AnalyticalSolver(const Settings& settings) :
    settings_(settings),
    gamma_(settings.gamma),
    x0_(settings.x0),
    has_initial_states_(false),
    dt_external_(-1.0) {
    exact_.SetQ(settings_.Q);
}

void AnalyticalSolver::SetDt(double dt) {
    dt_external_ = dt;
}

void AnalyticalSolver::InitializeStatesFromLayer(const DataLayer& layer) {
    const int coreStart = layer.GetCoreStart(0);
    const int coreEnd = layer.GetCoreEndExclusive(0);

    // Find interface index near x0_
    int idx = coreStart;
    while (idx < coreEnd && layer.xc(idx) < x0_) {
        ++idx;
    }

    const int iL = std::max(coreStart, idx - 1);
    const int iR = std::min(coreEnd - 1, idx);

    Primitive wl;
    wl.rho = layer.rho(iL);
    wl.u = layer.u(iL);
    wl.P = layer.P(iL);

    Primitive wr;
    wr.rho = layer.rho(iR);
    wr.u = layer.u(iR);
    wr.P = layer.P(iR);

    left_ = wl;
    right_ = wr;
    has_initial_states_ = true;
}

double AnalyticalSolver::Step(DataLayer& layer, double& t_cur) {
    const double t_end = settings_.t_end;
    if (t_cur >= t_end) {
        return 0.0;
    }

    // On first call, read left/right states from the initialized DataLayer
    if (!has_initial_states_) {
        InitializeStatesFromLayer(layer);
    }

    double dt = dt_external_;
    if (dt <= 0.0) {
        dt = t_end - t_cur;
    }
    if (t_cur + dt > t_end) {
        dt = t_end - t_cur;
    }
    if (dt <= 0.0) {
        return 0.0;
    }

    const double t = t_cur + dt;

    const int coreStart = layer.GetCoreStart(0);
    const int coreEnd = layer.GetCoreEndExclusive(0);

    double dx = 0.0;
    if (coreEnd - coreStart > 1) {
        dx = layer.xc(coreStart + 1) - layer.xc(coreStart);
    }

    for (int i = coreStart; i < coreEnd; ++i) {
        const double x = layer.xc(i);
        const double xi = t > 0.0 ? (x - x0_) / t : 0.0;

        Primitive w = exact_.Sample(left_, right_, gamma_, xi);

        const double rho = w.rho;
        const double u = w.u;
        const double P = w.P;

        layer.rho(i) = rho;
        layer.u(i) = u;
        layer.P(i) = P;

        // Momentum density
        layer.p(i) = rho * u;

        // Specific volume
        layer.V(i) = rho > 0.0 ? 1.0 / rho : 0.0;

        // Energies: match GodunovSolver semantics
        // Total energy density:
        const double kinetic = 0.5 * rho * u * u;
        const double Eint = P / (gamma_ - 1.0);
        const double Etot = Eint + kinetic;

        // specific internal energy
        const double eint = rho > 0.0 ? Eint / rho : 0.0;

        // In your GodunovSolver:
        // U = specific internal energy, e = total energy density.
        layer.U(i) = eint;
        layer.e(i) = Etot;

        // Cell mass
        layer.m(i) = rho * dx;
    }

    return dt;
}

void AnalyticalSolver::SetCfl(double cfl) {
    (void)cfl;
}

void AnalyticalSolver::AddBoundary(int axis,
                                   std::shared_ptr<BoundaryCondition> left_bc,
                                   std::shared_ptr<BoundaryCondition> right_bc) {
    (void)axis;
    (void)left_bc;
    (void)right_bc;
}