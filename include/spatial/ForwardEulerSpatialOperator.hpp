#ifndef FORWARDEULERSPATIALOPERATOR_HPP
#define FORWARDEULERSPATIALOPERATOR_HPP

#include <memory>

#include "spatial/SpatialOperator.hpp"
#include "viscosity/ArtificialViscosity.hpp"
#include "config/Settings.hpp"

/**
 * @class ForwardEulerSpatialOperator
 * @brief Forward-difference RHS operator for MacCormack predictor with optional artificial viscosity.
 *
 * Computes:
 *   L_j(U) = - (F_{j+1} - F_j) / dx,
 * where F_j is the physical Euler flux evaluated at cell centers.
 *
 * If viscosity is enabled, an effective pressure is used:
 *   P_eff(j) = P(j) + q_cell(j),
 * where q_cell is obtained from interface viscosity q_{j+1/2} by averaging.
 */
class ForwardEulerSpatialOperator : public SpatialOperator {
public:
    ForwardEulerSpatialOperator() = default;

    /**
     * @brief Constructs operator with artificial viscosity model.
     *
     * @param settings Global settings used to construct viscosity model.
     */
    explicit ForwardEulerSpatialOperator(const Settings& settings);

    ~ForwardEulerSpatialOperator() override = default;

    void ComputeRHS(const DataLayer& layer,
                    double dx,
                    double gamma,
                    xt::xarray<Conservative>& rhs) const override;

private:
    std::shared_ptr<ArtificialViscosity> viscosity_;
};

#endif  // FORWARDEULERSPATIALOPERATOR_HPP
