#ifndef BACKWARDEULERSPATIALOPERATOR_HPP
#define BACKWARDEULERSPATIALOPERATOR_HPP

#include "spatial/SpatialOperator.hpp"

/**
 * @class BackwardEulerSpatialOperator
 * @brief Backward-Euler-style spatial operator for MacCormack corrector.
 *
 * Implements the backward-difference form
 *
 *   L_j(U) = - ( F_j - F_{j-1} ) / dx,
 *
 * where F_j is the physical Euler flux evaluated at cell centers.
 * Ghost cells are assumed to be already filled by boundary conditions.
 */
class BackwardEulerSpatialOperator : public SpatialOperator {
public:
    BackwardEulerSpatialOperator() = default;
    ~BackwardEulerSpatialOperator() override = default;

    /**
     * @brief Computes backward-difference RHS using cell-center Euler fluxes.
     *
     * On core cells j = [core_start, core_end):
     *
     *   rhs_j.rho  = - (F_j.mass     - F_{j-1}.mass)     / dx
     *   rhs_j.rhoU = - (F_j.momentum - F_{j-1}.momentum) / dx
     *   rhs_j.E    = - (F_j.energy   - F_{j-1}.energy)   / dx
     *
     * Non-core entries in rhs are set to zero.
     */
    void ComputeRHS(const DataLayer& layer,
                    double dx,
                    double gamma,
                    xt::xarray<Conservative>& rhs) const override;
};

#endif  // BACKWARDEULERSPATIALOPERATOR_HPP