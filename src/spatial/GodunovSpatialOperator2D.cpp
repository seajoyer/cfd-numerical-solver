#include "spatial/GodunovSpatialOperator2D.hpp"

#include <algorithm>
#include <cctype>
#include <string>
#include <vector>

#include "data/DataLayer.hpp"
#include "reconstruction/ENOReconstruction.hpp"
#include "reconstruction/P0Reconstruction.hpp"
#include "reconstruction/P1Reconstruction.hpp"
#include "reconstruction/WENOReconstruction.hpp"
#include "riemann/AcousticRiemannSolver.hpp"
#include "riemann/ExactIdealGasRiemannSolver.hpp"
#include "riemann/RusanovRiemannSolver.hpp"
#include "riemann/HLLCRiemannSolver.hpp"
#include "riemann/HLLRiemannSolver.hpp"
#include "riemann/RoeRiemannSolver.hpp"
#include "riemann/OsherRiemannSolver.hpp"
#include "viscosity/VNRArtificialViscosity.hpp"

GodunovSpatialOperator2D::GodunovSpatialOperator2D(const Settings& settings)
    : dim_(settings.dim) {
    InitializeReconstruction(settings);
    InitializeRiemannSolver(settings);
    viscosity_ = nullptr;
    if (settings.viscosity) {
        viscosity_ = std::make_unique<VNRArtificialViscosity>(settings);
    }
}

void GodunovSpatialOperator2D::InitializeReconstruction(const Settings& settings) {
    std::string name = settings.reconstruction;

    if (name.find("p1") != std::string::npos) {
        reconstruction_ = std::make_shared<P1Reconstruction>();
    } else if (name.find("p0") != std::string::npos) {
        reconstruction_ = std::make_shared<P0Reconstruction>();
    } else if (name.starts_with("eno")) {
        int order = 3;
        try { order = std::stoi(name.substr(3)); } catch (...) {}
        reconstruction_ = std::make_shared<ENOReconstruction>(order);
    } else if (name.starts_with("weno")) {
        int order = 5;
        try { order = std::stoi(name.substr(4)); } catch (...) {}
        reconstruction_ = std::make_shared<WENOReconstruction>(order);
    } else {
        reconstruction_ = std::make_shared<P0Reconstruction>();
    }
}

void GodunovSpatialOperator2D::InitializeRiemannSolver(const Settings& settings) {
    std::string name = settings.riemann_solver;
    std::string lower(name.size(), '\0');
    std::transform(name.begin(), name.end(), lower.begin(),
                   [](unsigned char c) -> char { return static_cast<char>(std::tolower(c)); });

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

// ============================================================================
// ComputeRHS dispatcher
// ============================================================================

void GodunovSpatialOperator2D::ComputeRHS(const DataLayer& layer,
                                           double dx,
                                           double gamma,
                                           xt::xarray<Conservative>& rhs) const {
    if (dim_ <= 1 || layer.GetDim() <= 1) {
        ComputeRHS1D(layer, dx, gamma, rhs);
    } else {
        // For the 1-argument dx overload in 2D, assume dy = dx
        ComputeRHS2D(layer, dx, dx, gamma, rhs);
    }
}

void GodunovSpatialOperator2D::ComputeRHS2D(const DataLayer& layer,
                                             double dx, double dy,
                                             double gamma,
                                             xt::xarray<Conservative>& rhs) const {
    const int tx = layer.GetTotalSize(0);
    const int ty = layer.GetTotalSize(1);

    // Initialize 2D RHS to zero
    rhs = xt::xarray<Conservative>::from_shape(
        {static_cast<std::size_t>(tx), static_cast<std::size_t>(ty)});
    for (int i = 0; i < tx; ++i) {
        for (int j = 0; j < ty; ++j) {
            rhs(i, j) = Conservative{};
        }
    }

    // Accumulate x-direction flux contributions
    ComputeXFluxes(layer, dx, gamma, rhs);

    // Accumulate y-direction flux contributions
    ComputeYFluxes(layer, dy, gamma, rhs);
}

// ============================================================================
// X-direction flux sweep
// ============================================================================

void GodunovSpatialOperator2D::ComputeXFluxes(const DataLayer& layer,
                                               double dx, double gamma,
                                               xt::xarray<Conservative>& rhs) const {
    const int tx = layer.GetTotalSize(0);
    const int ty = layer.GetTotalSize(1);
    const int core_start_x = layer.GetCoreStart(0);
    const int core_end_x = layer.GetCoreEndExclusive(0);
    const int core_start_y = layer.GetCoreStart(1);
    const int core_end_y = layer.GetCoreEndExclusive(1);
    const int n_x_interfaces = tx - 1;

    // For each row j (including ghost rows so we have complete flux coverage)
    for (int j = 0; j < ty; ++j) {
        // Build a temporary 1D DataLayer for this row
        // Instead, we directly build left/right states from 2D data
        const int n_iface = n_x_interfaces;

        // Reconstruct interface states along x for this row
        // For P0: left_state(i) = cell(i,j), right_state(i) = cell(i+1,j)
        std::vector<Primitive> left_states(n_iface);
        std::vector<Primitive> right_states(n_iface);

        for (int i = 0; i < n_iface; ++i) {
            // For P0 reconstruction: piecewise constant
            // The existing Reconstruction interface works on a 1D DataLayer,
            // so we manually do P0 here for 2D. For higher-order, we'd need
            // to extract a 1D slice.
            Primitive wL = layer.GetPrimitive2D(i, j);
            Primitive wR = layer.GetPrimitive2D(i + 1, j);

            // For the 1D Riemann solver, we pass the x-normal velocity as 'u'
            // and store the transverse velocity separately.
            left_states[i] = wL;
            right_states[i] = wR;
        }

        // Compute fluxes at each x-interface
        std::vector<Flux> fluxes(n_iface);
        for (int i = 0; i < n_iface; ++i) {
            // Create 1D Primitive with x-velocity as normal velocity
            // The 1D Riemann solver uses .u as the normal velocity
            Primitive wL_1d(left_states[i].rho, left_states[i].u, left_states[i].P);
            Primitive wR_1d(right_states[i].rho, right_states[i].u, right_states[i].P);

            // Solve 1D Riemann problem (returns mass, momentum_x, energy fluxes)
            Flux f1d = riemann_solver_->ComputeFlux(wL_1d, wR_1d, gamma);

            // Build the full 2D flux in x-direction:
            // F = (rho*u, rho*u^2+P, rho*u*v, u*(E+P))
            // The 1D solver gives us mass, x-momentum, and energy correctly.
            // We need to add the y-momentum flux = mass_flux * v_interface
            // For the Godunov method, v_interface is upwinded by the mass flux:
            double v_at_interface;
            if (f1d.mass >= 0.0) {
                v_at_interface = left_states[i].v;
            } else {
                v_at_interface = right_states[i].v;
            }

            fluxes[i].mass = f1d.mass;
            fluxes[i].momentum_x = f1d.momentum_x;  // rho*u^2+P from Riemann solver
            fluxes[i].momentum_y = f1d.mass * v_at_interface;  // rho*u*v (passive advection)
            fluxes[i].energy = f1d.energy;

            // Note: The energy flux from the 1D solver already accounts for
            // the total energy E = e_int + 0.5*rho*u^2, but in 2D the total
            // energy is E = e_int + 0.5*rho*(u^2+v^2). The 1D Riemann solver
            // computed the flux using E_1d = P/(gamma-1) + 0.5*rho*u^2.
            // We need to correct for the missing v^2 kinetic energy.
            // Energy flux correction: add 0.5*rho*v^2 * u to the energy flux
            // Actually, the proper way: the 1D solver's energy flux is
            // u*(E_1d + P) where E_1d misses 0.5*rho*v^2.
            // Correction: energy_flux += mass_flux * 0.5 * v^2
            double v_energy = v_at_interface;
            fluxes[i].energy += f1d.mass * 0.5 * v_energy * v_energy;
        }

        // Accumulate flux differences into RHS for core cells in this row
        if (j >= core_start_y && j < core_end_y) {
            for (int i = core_start_x; i < core_end_x; ++i) {
                const Flux& Fminus = fluxes[i - 1];  // F_{i-1/2}
                const Flux& Fplus = fluxes[i];        // F_{i+1/2}

                rhs(i, j).rho  -= (Fplus.mass - Fminus.mass) / dx;
                rhs(i, j).rhoU -= (Fplus.momentum_x - Fminus.momentum_x) / dx;
                rhs(i, j).rhoV -= (Fplus.momentum_y - Fminus.momentum_y) / dx;
                rhs(i, j).E    -= (Fplus.energy - Fminus.energy) / dx;
            }
        }
    }
}

// ============================================================================
// Y-direction flux sweep
// ============================================================================

void GodunovSpatialOperator2D::ComputeYFluxes(const DataLayer& layer,
                                               double dy, double gamma,
                                               xt::xarray<Conservative>& rhs) const {
    const int tx = layer.GetTotalSize(0);
    const int ty = layer.GetTotalSize(1);
    const int core_start_x = layer.GetCoreStart(0);
    const int core_end_x = layer.GetCoreEndExclusive(0);
    const int core_start_y = layer.GetCoreStart(1);
    const int core_end_y = layer.GetCoreEndExclusive(1);
    const int n_y_interfaces = ty - 1;

    // For each column i
    for (int i = 0; i < tx; ++i) {
        const int n_iface = n_y_interfaces;

        // Reconstruct interface states along y for this column
        std::vector<Primitive> bottom_states(n_iface);
        std::vector<Primitive> top_states(n_iface);

        for (int j = 0; j < n_iface; ++j) {
            bottom_states[j] = layer.GetPrimitive2D(i, j);
            top_states[j] = layer.GetPrimitive2D(i, j + 1);
        }

        // Compute fluxes at each y-interface using rotated velocities
        std::vector<Flux> fluxes(n_iface);
        for (int j = 0; j < n_iface; ++j) {
            // Rotate: y-velocity becomes the normal velocity for the 1D solver
            Primitive wB_1d(bottom_states[j].rho, bottom_states[j].v, bottom_states[j].P);
            Primitive wT_1d(top_states[j].rho, top_states[j].v, top_states[j].P);

            // Solve 1D Riemann problem with v as normal velocity
            Flux g1d = riemann_solver_->ComputeFlux(wB_1d, wT_1d, gamma);

            // Build the full 2D flux in y-direction:
            // G = (rho*v, rho*u*v, rho*v^2+P, v*(E+P))
            // The 1D solver gives: mass=rho*v, momentum=rho*v^2+P, energy=v*(E_1d+P)
            // We need: x-momentum flux = mass_flux * u (passive advection of u by v)
            double u_at_interface;
            if (g1d.mass >= 0.0) {
                u_at_interface = bottom_states[j].u;
            } else {
                u_at_interface = top_states[j].u;
            }

            fluxes[j].mass = g1d.mass;                           // rho*v
            fluxes[j].momentum_x = g1d.mass * u_at_interface;    // rho*v*u (passive)
            fluxes[j].momentum_y = g1d.momentum_x;               // rho*v^2+P (from solver)
            fluxes[j].energy = g1d.energy;

            // Energy correction: add missing u^2 kinetic energy contribution
            double u_energy = u_at_interface;
            fluxes[j].energy += g1d.mass * 0.5 * u_energy * u_energy;
        }

        // Accumulate flux differences into RHS for core cells in this column
        if (i >= core_start_x && i < core_end_x) {
            for (int j = core_start_y; j < core_end_y; ++j) {
                const Flux& Gminus = fluxes[j - 1];  // G_{j-1/2}
                const Flux& Gplus = fluxes[j];        // G_{j+1/2}

                rhs(i, j).rho  -= (Gplus.mass - Gminus.mass) / dy;
                rhs(i, j).rhoU -= (Gplus.momentum_x - Gminus.momentum_x) / dy;
                rhs(i, j).rhoV -= (Gplus.momentum_y - Gminus.momentum_y) / dy;
                rhs(i, j).E    -= (Gplus.energy - Gminus.energy) / dy;
            }
        }
    }
}

// ============================================================================
// 1D fallback (identical to original GodunovSpatialOperator)
// ============================================================================

void GodunovSpatialOperator2D::ComputeRHS1D(const DataLayer& layer,
                                             double dx, double gamma,
                                             xt::xarray<Conservative>& rhs) const {
    const int total_size = layer.GetTotalSize();
    if (total_size < 3 || dx <= 0.0) {
        rhs = xt::xarray<Conservative>::from_shape(
            {static_cast<std::size_t>(std::max(total_size, 0))});
        for (int j = 0; j < total_size; ++j) rhs(j) = Conservative{};
        return;
    }

    const int core_start = layer.GetCoreStart(0);
    const int core_end = layer.GetCoreEndExclusive(0);
    const int n_core = core_end - core_start;
    if (n_core < 2) {
        rhs = xt::xarray<Conservative>::from_shape({static_cast<std::size_t>(total_size)});
        for (int j = 0; j < total_size; ++j) rhs(j) = Conservative{};
        return;
    }

    const int n_interfaces = total_size - 1;

    rhs = xt::xarray<Conservative>::from_shape({static_cast<std::size_t>(total_size)});

    xt::xarray<Primitive> left_states;
    xt::xarray<Primitive> right_states;
    reconstruction_->ReconstructStates(layer, left_states, right_states);

    xt::xarray<Flux> fluxes =
        xt::xarray<Flux>::from_shape({static_cast<std::size_t>(n_interfaces)});

    for (int i = 0; i < n_interfaces; ++i) {
        Primitive WL = left_states(i);
        Primitive WR = right_states(i);
        fluxes(i) = riemann_solver_->ComputeFlux(WL, WR, gamma);
    }

    for (int j = 0; j < total_size; ++j) rhs(j) = Conservative{};

    for (int j = core_start; j < core_end; ++j) {
        const Flux& Fminus = fluxes(j - 1);
        const Flux& Fplus = fluxes(j);
        const Flux dF = Flux::Diff(Fplus, Fminus);
        rhs(j) = dF * (-1.0 / dx);
    }
}
