#ifndef REFLECTIVEBOUNDARY_HPP
#define REFLECTIVEBOUNDARY_HPP

#include "bc/BoundaryCondition.hpp"

/**
 * @class ReflectiveBoundary
 * @brief Slip reflective boundary.
 *
 * Normal velocity component changes sign.
 * Tangential velocity components are preserved.
 */
class ReflectiveBoundary final : public BoundaryCondition {
public:
    [[nodiscard]] PrimitiveCell BuildExteriorState(const DataLayer& layer,
                                                   const Mesh& mesh,
                                                   const Face& face,
                                                   const PrimitiveCell& interior_state) const override;
};

#endif  // REFLECTIVEBOUNDARY_HPP
