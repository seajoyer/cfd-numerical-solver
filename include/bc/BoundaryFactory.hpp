#ifndef BOUNDARYFACTORY_HPP
#define BOUNDARYFACTORY_HPP

#include <memory>
#include <string>

#include "config/Settings.hpp"
#include "data/Variables.hpp"

class BoundaryCondition;

/**
 * @class BoundaryFactory
 * @brief Factory for creating boundary conditions by type.
 */
class BoundaryFactory final {
public:
    /**
     * @brief Create boundary condition from one boundary-condition settings block.
     */
    [[nodiscard]] static std::shared_ptr<BoundaryCondition> Create(
        const BoundaryConditionSettings& boundary_settings,
        const Settings& settings
    );

private:
    [[nodiscard]] static PrimitiveCell PrimitiveFromBoundaryState(
        const BoundaryStateSettings& state
    );
};

#endif  // BOUNDARYFACTORY_HPP
