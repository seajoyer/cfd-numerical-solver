#ifndef BOUNDARYFACTORY_HPP
#define BOUNDARYFACTORY_HPP

#include <memory>
#include <string>

class BoundaryCondition;

/**
 * @struct FarfieldConservative
 * @brief Conservative far-field state (rho, rhoU, rhoV, rhoW, E).
 *
 * Used by boundary conditions that require a prescribed state (e.g., free-stream, inlet).
 */
struct FarfieldConservative final {
    double rho  = 0.0;
    double rhoU = 0.0;
    double rhoV = 0.0;
    double rhoW = 0.0;
    double E    = 0.0;
};

/**
 * @class BoundaryFactory
 * @brief Factory for creating boundary condition instances by string type.
 *
 * Types are expected to come from config parser as strings, e.g.:
 * - "outlet"
 * - "free_stream" (future)
 * - "inlet"       (future)
 * etc.
 *
 * For now only "outlet" is guaranteed to be supported in the smoke-test pipeline.
 */
class BoundaryFactory {
public:
    /**
     * @brief Creates a boundary condition instance by type.
     *
     * @param boundary_type Boundary condition type string.
     * @return Shared pointer to boundary condition instance.
     * @throws std::runtime_error if boundary_type is unknown.
     */
    static auto Create(const std::string& boundary_type) -> std::shared_ptr<BoundaryCondition>;

    /**
     * @brief Creates a boundary condition instance by type with far-field conservative state.
     *
     * @param boundary_type Boundary condition type string.
     * @param farfield_U Far-field conservative state.
     * @return Shared pointer to boundary condition instance.
     * @throws std::runtime_error if boundary_type is unknown or type does not use farfield.
     */
    static auto Create(const std::string& boundary_type, const FarfieldConservative& farfield_U)
        -> std::shared_ptr<BoundaryCondition>;
};

#endif  // BOUNDARYFACTORY_HPP