#ifndef BOUNDARYFACTORY_HPP
#define BOUNDARYFACTORY_HPP

#include <memory>
#include <string>

#include "bc/BoundaryCondition.hpp"

/**
 * @class BoundaryFactory
 * @brief Factory class for creating boundary condition instances.
 *
 * This class provides a centralized way to create different boundary condition types
 * without coupling the Simulation class to specific boundary implementations.
 */
class BoundaryFactory {
   public:
    /**
     * @brief Creates a boundary condition instance based on the type string.
     *
     * Supported boundary types:
     * - "free_stream": Free stream boundary condition
     * - "inlet": Inlet boundary condition
     * - "outlet": Outlet boundary condition
     * - "reflective": Reflective boundary condition
     * - "non_reflective": Non-reflective boundary condition
     * - "periodic": Periodic boundary condition
     * - "symmetry": Symmetry boundary condition
     * - "wall": Wall boundary condition
     *
     * For "free_stream" and "inlet", the provided parameters (rho_inf, u_inf, p_inf,
     * gamma) are used to initialize the boundary condition. For other types, these
     * parameters are ignored, assuming those boundary conditions use default or
     * no-argument constructors.
     *
     * @param boundary_type String identifier for the boundary condition.
     * @param rho_inf External/inlet density (default 0.0).
     * @param u_inf External/inlet velocity (default 0.0).
     * @param p_inf External/inlet pressure (default 0.0).
     * @return Shared pointer to the created boundary condition.
     * @throws std::runtime_error if boundary type is not recognized.
     */
    static auto Create(const std::string& boundary_type, double rho_inf = 0.0,
                       double u_inf = 0.0, double p_inf = 0.0)
        -> std::shared_ptr<BoundaryCondition>;
};

#endif  // BOUNDARYFACTORY_HPP
