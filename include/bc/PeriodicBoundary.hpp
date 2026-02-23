#ifndef PERIODICBOUNDARY_HPP
#define PERIODICBOUNDARY_HPP

#include "bc/BoundaryCondition.hpp"
#include "data/Variables.hpp"

/**
 * @class PeriodicBoundary
 * @brief Periodic boundary condition for conservative state U.
 *
 * Periodicity is enforced by copying conservative state U from the opposite
 * side of the core domain into ghost layers along the selected axis.
 *
 * Works only with conservative state U(var,i,j,k).
 */
class PeriodicBoundary final : public BoundaryCondition {
public:
    void Apply(DataLayer& layer, Axis axis, Side side) const override;
};

#endif  // PERIODICBOUNDARY_HPP