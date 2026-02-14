#ifndef BOUNDARYFACTORY_HPP
#define BOUNDARYFACTORY_HPP

#include <memory>
#include <string>

#include "bc/BoundaryCondition.hpp"

/**
 * @class BoundaryFactory
 * @brief Factory for creating boundary condition instances.
 */
class BoundaryFactory {
public:
    /**
     * @brief Creates a 1D boundary condition.
     */
    static auto Create(const std::string& type, double rho_inf, double u_inf, double p_inf)
        -> std::shared_ptr<BoundaryCondition>;

    /**
     * @brief Creates a 2D boundary condition with both u and v far-field values.
     */
    static auto Create2D(const std::string& type, double rho_inf,
                         double u_inf, double v_inf, double p_inf)
        -> std::shared_ptr<BoundaryCondition>;
};

#endif  // BOUNDARYFACTORY_HPP
