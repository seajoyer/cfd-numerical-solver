#ifndef REFLECTIVEBOUNDARY_HPP
#define REFLECTIVEBOUNDARY_HPP

#include "bc/BoundaryCondition.hpp"

/**
 * @class ReflectiveBoundary
 * @brief Reflecting wall boundary condition for 1D and 2D.
 *
 * Mirrors the state across the boundary and reverses the normal velocity component.
 * In 2D, only the velocity component normal to the boundary is reflected.
 */
class ReflectiveBoundary : public BoundaryCondition {
public:
    void Apply(DataLayer& layer, int axis, Side side) const override;
private:
    void Apply2D(DataLayer& layer, int axis, Side side) const;
};

#endif  // REFLECTIVEBOUNDARY_HPP
