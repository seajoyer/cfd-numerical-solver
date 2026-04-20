#ifndef BOUNDARYFACTORY_HPP
#define BOUNDARYFACTORY_HPP

#include <memory>
#include <string>

#include "config/Settings.hpp"

class BoundaryCondition;
class MPIContext;

/**
 * @struct FarfieldConservative
 * @brief Conservative far-field state (rho, rhoU, rhoV, rhoW, E).
 */
struct FarfieldConservative final {
    double rho = 0.0;
    double rhoU = 0.0;
    double rhoV = 0.0;
    double rhoW = 0.0;
    double E = 0.0;
};

/**
 * @class BoundaryFactory
 * @brief Factory for creating boundary condition instances by type.
 */
class BoundaryFactory {
public:
    static auto Create(const std::string& boundary_type) -> std::shared_ptr<BoundaryCondition>;

    static auto Create(const std::string& boundary_type,
                       const FarfieldConservative& farfield_U,
                       const Settings& settings,
                       int mpi_size) -> std::shared_ptr<BoundaryCondition>;

    static auto Create(const std::string& boundary_type,
                       const BoundaryStateSettings& primitive_state,
                       const Settings& settings,
                       int mpi_size) -> std::shared_ptr<BoundaryCondition>;
};

#endif  // BOUNDARYFACTORY_HPP
