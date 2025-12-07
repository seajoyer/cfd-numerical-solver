#ifndef GODUNOVKOLGANRADIONOVSPATIALOPERATOR
#define GODUNOVKOLGANRADIONOVSPATIALOPERATOR

#include <memory>

#include "SpatialOperator.hpp"
#include "config/Settings.hpp"
#include "reconstruction/Reconstruction.hpp"
#include "riemann/RiemannSolver.hpp"

/**
 * @class GodunovKolganRodionovSpatialOperator
 * @brief Second-order Godunov-Kolgan-Rodionov spatial operator.
 *
 * This operator:
 *  - performs piecewise-linear (or higher-order ENO/WENO) reconstruction
 *    of primitive variables,
 *  - advances reconstructed edge values by half a time step (predictor),
 *  - builds half-time interface states and calls a Riemann solver,
 *  - assembles the resulting flux differences into dU/dt.
 *
 * Time integration (including the choice of dt) is handled by the
 * enclosing TimeIntegrator / FiniteVolumeSolver.
 */
class GodunovKolganRodionovSpatialOperator : public SpatialOperator {
public:
    /**
     * @brief Constructs a Godunov-Kolgan-Rodionov spatial operator from Settings.
     *
     * Uses:
     *  - settings.reconstruction to choose P1/ENO/WENO (forcing ≥ P1),
     *  - settings.riemann_solver to choose HLL/HLLC/Exact.
     *
     * @param settings Global simulation settings.
     */
    explicit GodunovKolganRodionovSpatialOperator(const Settings& settings);

    /**
     * @brief Computes dU/dt using a Godunov-Kolgan-Rodionov predictor and Riemann fluxes.
     *
     * The method internally uses a half-time predictor based on dt and dx
     * (dt is provided by the TimeIntegrator via its calling pattern).
     *
     * @param layer Mesh and auxiliary data.
     * @param dx    Cell size.
     * @param gamma Ratio of specific heats.
     * @param rhs   Output array for dU/dt.
     */
    void ComputeRHS(const DataLayer& layer,
                    double dx,
                    double gamma,
                    xt::xarray<Conservative>& rhs) const override;

    /**
     * @brief Sets the local time step used for the half-step predictor.
     *
     * Many time integrators will call this before evaluating dU/dt
     * if a predictor requires knowledge of the current dt.
     *
     * @param dt Time step size used inside the predictor.
     */
    void SetLocalTimeStep(double dt);

private:
    /** @brief Reconstruction strategy (typically P1, ENO, or WENO). */
    std::shared_ptr<Reconstruction> reconstruction_;

    /** @brief Riemann solver. */
    std::shared_ptr<RiemannSolver> riemann_solver_;

    /** @brief Local dt used for the half-time predictor. */
    double dt_local_;

    void InitializeReconstruction(const Settings& settings);
    void InitializeRiemannSolver(const Settings& settings);
};

#endif  // GODUNOVKOLGANRADIONOVSPATIALOPERATOR
