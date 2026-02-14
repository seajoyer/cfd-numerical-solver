#ifndef GODUNOVSPATIALOPERATOR2D_HPP
#define GODUNOVSPATIALOPERATOR2D_HPP

#include <memory>

#include "SpatialOperator.hpp"
#include "viscosity/ArtificialViscosity.hpp"
#include "config/Settings.hpp"
#include "reconstruction/Reconstruction.hpp"
#include "riemann/RiemannSolver.hpp"

/**
 * @class GodunovSpatialOperator2D
 * @brief 2D Godunov spatial operator using unsplit finite volume method.
 *
 * This operator computes the full 2D semi-discrete RHS:
 *
 *   dU/dt = -dF/dx - dG/dy
 *
 * by performing:
 *   1. X-sweep: for each row j, reconstruct states along x, solve Riemann
 *      problems at x-interfaces, compute x-flux differences.
 *   2. Y-sweep: for each column i, reconstruct states along y (with rotated
 *      velocities), solve Riemann problems at y-interfaces, compute y-flux
 *      differences.
 *   3. Sum both contributions for the total RHS.
 *
 * Existing 1D Riemann solvers are reused via velocity rotation:
 * - For x-sweep: Primitive = (rho, u, P), v is passively advected.
 * - For y-sweep: Primitive = (rho, v, P), u is passively advected.
 *
 * @note Currently only supports P0 (piecewise constant) reconstruction in 2D.
 *       The reconstruction_ member is only used for 1D fallback.
 * @note Falls back to 1D behavior when dim=1.
 */
class GodunovSpatialOperator2D : public SpatialOperator {
public:
    explicit GodunovSpatialOperator2D(const Settings& settings);

    /**
     * @brief Computes RHS. For dim=1, delegates to 1D logic.
     *        For dim=2, calls ComputeRHS2D with dy = dx (square cells).
     */
    void ComputeRHS(const DataLayer& layer,
                    double dx,
                    double gamma,
                    xt::xarray<Conservative>& rhs) const override;

    /**
     * @brief Computes 2D RHS with separate dx, dy (the primary 2D entry point).
     */
    void ComputeRHS2D(const DataLayer& layer,
                      double dx, double dy,
                      double gamma,
                      xt::xarray<Conservative>& rhs) const override;

private:
    std::shared_ptr<Reconstruction> reconstruction_;
    std::shared_ptr<RiemannSolver> riemann_solver_;
    std::shared_ptr<ArtificialViscosity> viscosity_;
    int dim_;

    void InitializeReconstruction(const Settings& settings);
    void InitializeRiemannSolver(const Settings& settings);

    /**
     * @brief Computes x-direction flux contributions for all cells.
     */
    void ComputeXFluxes(const DataLayer& layer,
                        double dx, double gamma,
                        xt::xarray<Conservative>& rhs) const;

    /**
     * @brief Computes y-direction flux contributions for all cells.
     */
    void ComputeYFluxes(const DataLayer& layer,
                        double dy, double gamma,
                        xt::xarray<Conservative>& rhs) const;

    /**
     * @brief 1D RHS computation (backward compatible fallback).
     */
    void ComputeRHS1D(const DataLayer& layer,
                      double dx, double gamma,
                      xt::xarray<Conservative>& rhs) const;
};

#endif  // GODUNOVSPATIALOPERATOR2D_HPP
