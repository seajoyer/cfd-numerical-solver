#ifndef NONREFLECTIVEBOUNDARY_HPP
#define NONREFLECTIVEBOUNDARY_HPP

#include "bc/BoundaryCondition.hpp"
#include "data/Variables.hpp"

/**
 * @class NonReflectiveBoundary
 * @brief Placeholder non-reflective boundary condition for conservative Euler state U.
 *
 * Current implementation enforces a zero-gradient (outlet-like) condition by copying
 * the nearest interior core layer into ghost layers along the selected axis.
 *
 * This is a simplified fallback; characteristic-based radiation conditions can be
 * added later without changing the interface.
 */
class NonReflectiveBoundary final : public BoundaryCondition {
public:
    void Apply(DataLayer& layer, Axis axis, Side side) const override;
};

#endif  // NONREFLECTIVEBOUNDARY_HPP