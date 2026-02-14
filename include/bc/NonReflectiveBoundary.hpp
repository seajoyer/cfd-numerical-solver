#ifndef NONREFLECTIVEBOUNDARY_HPP
#define NONREFLECTIVEBOUNDARY_HPP

#include "bc/BoundaryCondition.hpp"

/**
 * @class NonReflectiveBoundary
 * @brief Zero-gradient (extrapolation) boundary condition for 1D and 2D.
 *
 * Copies the nearest core cell state into all ghost cells, allowing
 * waves to exit the domain without reflection.
 *
 * For 2D, fills ghost columns (axis=0) or ghost rows (axis=1) by
 * extrapolating from the nearest core cell along the boundary-normal direction.
 */
class NonReflectiveBoundary : public BoundaryCondition {
public:
    void Apply(DataLayer& layer, int axis, Side side) const override;

private:
    void Apply2D(DataLayer& layer, int axis, Side side) const;
};

#endif  // NONREFLECTIVEBOUNDARY_HPP
