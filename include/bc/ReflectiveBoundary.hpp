#ifndef REFLECTIVEBOUNDARY_HPP
#define REFLECTIVEBOUNDARY_HPP

#include "bc/BoundaryCondition.hpp"
#include "data/DataLayer.hpp"
#include "data/Variables.hpp"

/**
 * @class ReflectiveBoundary
 * @brief Reflective (slip-wall) boundary condition for conservative Euler state U.
 *
 * Ghost cells are filled by mirroring interior cells across the boundary.
 * The normal momentum component is inverted, tangential components are copied.
 *
 * For axis:
 *  - X: rhoU changes sign
 *  - Y: rhoV changes sign
 *  - Z: rhoW changes sign
 */
class ReflectiveBoundary final : public BoundaryCondition {
public:
    void Apply(DataLayer& layer, Axis axis, Side side) const override;
};

#endif  // REFLECTIVEBOUNDARY_HPP