#ifndef BOUNDARYMANAGER_HPP
#define BOUNDARYMANAGER_HPP

#include <memory>
#include <unordered_map>

#include "data/Variables.hpp"

class BoundaryCondition;
class DataLayer;
class Face;
class Mesh;

/**
 * @class BoundaryManager
 * @brief Stores boundary conditions indexed by boundary tag.
 */
class BoundaryManager final {
public:
    BoundaryManager() = default;

    /**
     * @brief Register boundary condition for one boundary tag.
     */
    void Register(int boundary_tag, std::shared_ptr<BoundaryCondition> boundary_condition);

    /**
     * @brief Check whether a boundary tag has a registered boundary condition.
     */
    [[nodiscard]] bool Has(int boundary_tag) const;

    /**
     * @brief Get boundary condition by boundary tag.
     * @throws std::runtime_error if boundary tag is not registered.
     */
    [[nodiscard]] const BoundaryCondition& Get(int boundary_tag) const;

    /**
     * @brief Build external primitive state for one boundary face.
     * @throws std::runtime_error if face is not boundary or no BC is registered for its tag.
     */
    [[nodiscard]] PrimitiveCell BuildExteriorState(const DataLayer& layer,
                                                   const Mesh& mesh,
                                                   const Face& face,
                                                   const PrimitiveCell& interior_state) const;

private:
    std::unordered_map<int, std::shared_ptr<BoundaryCondition>> by_tag_;
};

#endif  // BOUNDARYMANAGER_HPP
