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
 * @note This operator works with 2D DataLayer arrays.
 * @note Falls back to 1D behavior when dim=1.
 */
class GodunovSpatialOperator2D : public SpatialOperator {
public:
    explicit GodunovSpatialOperator2D(const Settings& settings);

    /**
     * @brief Computes 2D RHS. For dim=1, delegates to 1D logic.
     *
     * For dim=2:
     *  - rhs is a 2D xarray<Conservative> of shape (total_x, total_y)
     *  - Accumulates x-direction flux differences / dx
     *  - Accumulates y-direction flux differences / dy
     */
    void ComputeRHS(const DataLayer& layer,
                    double dx,
                    double gamma,
                    xt::xarray<Conservative>& rhs) const override;

    /**
     * @brief 2D version with separate dx, dy.
     */
    void ComputeRHS2D(const DataLayer& layer,
                      double dx, double dy,
                      double gamma,
                      xt::xarray<Conservative>& rhs) const;

private:
    std::shared_ptr<Reconstruction> reconstruction_;
    std::shared_ptr<RiemannSolver> riemann_solver_;
    std::shared_ptr<ArtificialViscosity> viscosity_;
    int dim_;

    void InitializeReconstruction(const Settings& settings);
    void InitializeRiemannSolver(const Settings& settings);

    /**
     * @brief Computes x-direction flux contributions for all cells.
     *
     * For each row j (from core_start_y to core_end_y), extracts a 1D slice
     * of primitive states, runs reconstruction and Riemann solver, and
     * accumulates -dF/dx into the rhs array.
     */
    void ComputeXFluxes(const DataLayer& layer,
                        double dx, double gamma,
                        xt::xarray<Conservative>& rhs) const;

    /**
     * @brief Computes y-direction flux contributions for all cells.
     *
     * For each column i (from core_start_x to core_end_x), extracts a 1D slice
     * with rotated velocities (v -> normal, u -> transverse), runs reconstruction
     * and Riemann solver, and accumulates -dG/dy into the rhs array.
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
