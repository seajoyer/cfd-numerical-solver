#ifndef BOUNDARYMANAGER_HPP
#define BOUNDARYMANAGER_HPP

#include <memory>
#include <vector>

struct DataLayer;
class BoundaryCondition;


struct AxisBc {
    std::shared_ptr<BoundaryCondition> min;
    std::shared_ptr<BoundaryCondition> max;
};

/**
 * @class BoundaryManager
 * @brief Manages and applies boundary conditions for each spatial axis.
 *
 * This class holds pairs of boundary condition objects for each dimension
 * (minimum and maximum sides) and provides a simple interface to apply
 * them all at once to a DataLayer.
 *
 * It allows flexible configuration of boundary types for each axis
 * without hardcoding logic into the solver.
 *
 * @note BoundaryManager enables extension to 2D/3D cases by maintaining
 *       independent boundary pairs for each coordinate direction.
 */
class BoundaryManager {
public:
    /**
     * @brief Constructs BoundaryManager for a given dimensionality.
     *
     * @param dim Number of spatial dimensions (1, 2, or 3).
     */
    explicit BoundaryManager(int dim);

    /**
     * @brief Assigns boundary conditions for a specific axis.
     *
     * @param axis Index of spatial axis (0 for X, 1 for Y, 2 for Z).
     * @param min Pointer to boundary condition for the lower side (Side::Min).
     * @param max Pointer to boundary condition for the upper side (Side::Max).
     *
     * @note Passing nullptr disables the corresponding boundary.
     */
    void Set(int axis,
             std::shared_ptr<BoundaryCondition> min,
             std::shared_ptr<BoundaryCondition> max);

    /**
     * @brief Applies all configured boundary conditions to the given DataLayer.
     *
     * Iterates over all axes and applies each pair (min/max) in order.
     *
     * @param layer Reference to the DataLayer whose ghost zones should be updated.
     */
    void ApplyAll(DataLayer &layer) const;

private:
    std::vector<AxisBc> axes;
};

#endif // BOUNDARYMANAGER_HPP
