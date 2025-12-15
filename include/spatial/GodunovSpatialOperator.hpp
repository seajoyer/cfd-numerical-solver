#ifndef GODUNOVSPATIALOPERATOR_HPP
#define GODUNOVSPATIALOPERATOR_HPP

#include <memory>

#include "SpatialOperator.hpp"
#include "viscosity/ArtificialViscosity.hpp"
#include "config/Settings.hpp"
#include "reconstruction/Reconstruction.hpp"
#include "riemann/RiemannSolver.hpp"

/**
 * @class GodunovSpatialOperator
 * @brief First-order Godunov finite-volume spatial operator.
 *
 * This operator:
 *  - reconstructs P0 (piecewise-constant) interface states,
 *  - evaluates numerical fluxes with a chosen Riemann solver,
 *  - assembles flux differences to compute dU/dt.
 *
 * It corresponds to the spatial part of the existing GodunovSolver,
 * but without any time integration or positivity limiting.
 */
class GodunovSpatialOperator : public SpatialOperator {
public:
    /**
     * @brief Constructs a Godunov spatial operator from Settings.
     *
     * Uses:
     *  - settings.reconstruction to choose P0/P1/ENO/WENO (P0 by default),
     *  - settings.riemann_solver to choose HLL/HLLC/Exact/Acoustic.
     *
     * @param settings Global simulation settings.
     */
    explicit GodunovSpatialOperator(const Settings& settings);

    /**
     * @brief Computes dU/dt using a first-order Godunov flux difference.
     *
     * Uses the same reconstruction and Riemann solver selection
     * as the current GodunovSolver, but operates purely on the
     * provided conservative array U and fills rhs with
     *
     *     rhs_j = -(F_{j+1/2} - F_{j-1/2}) / dx.
     */
    void ComputeRHS(const DataLayer& layer,
                    double dx,
                    double gamma,
                    xt::xarray<Conservative>& rhs) const override;

private:
    /** @brief Selected reconstruction strategy (P0/P1/ENO/WENO). */
    std::shared_ptr<Reconstruction> reconstruction_;

    /** @brief Selected Riemann solver (HLL, HLLC, Exact, Acoustic). */
    std::shared_ptr<RiemannSolver> riemann_solver_;

    /** @brief Von Neumann Richtmyer viscosity. */
    std::shared_ptr<ArtificialViscosity> viscosity_;

    /** @brief Initializes reconstruction_ based on settings.reconstruction. */
    void InitializeReconstruction(const Settings& settings);

    /** @brief Initializes riemann_solver_ based on settings.riemann_solver. */
    void InitializeRiemannSolver(const Settings& settings);
};

#endif  // GODUNOVSPATIALOPERATOR_HPP
