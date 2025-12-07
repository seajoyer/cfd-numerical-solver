#ifndef FORWARDEULERSPATIALOPERATOR_HPP
#define FORWARDEULERSPATIALOPERATOR_HPP

#include "spatial/SpatialOperator.hpp"

/**
 * @class ForwardEulerSpatialOperator
 * @brief Forward-Euler-style spatial operator for MacCormack predictor.
 *
 * Implements the forward-difference form
 *
 *   L_j(U) = - ( F_{j+1} - F_j ) / dx,
 *
 * where F_j is the physical Euler flux evaluated at cell centers.
 * Ghost cells are assumed to be already filled by boundary conditions.
 */
class ForwardEulerSpatialOperator : public SpatialOperator {
public:
    ForwardEulerSpatialOperator() = default;
    ~ForwardEulerSpatialOperator() override = default;

    /**
     * @brief Computes forward-difference RHS using cell-center Euler fluxes.
     *
     * On core cells j = [core_start, core_end):
     *
     *   rhs_j.rho  = - (F_{j+1}.mass     - F_j.mass)     / dx
     *   rhs_j.rhoU = - (F_{j+1}.momentum - F_j.momentum) / dx
     *   rhs_j.E    = - (F_{j+1}.energy   - F_j.energy)   / dx
     *
     * Non-core entries in rhs are set to zero.
     */
    void ComputeRHS(const DataLayer& layer,
                    double dx,
                    double gamma,
                    xt::xarray<Conservative>& rhs) const override;
};

#endif  // FORWARDEULERSPATIALOPERATOR_HPP
