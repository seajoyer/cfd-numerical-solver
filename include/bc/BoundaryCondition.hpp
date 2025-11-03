#ifndef BOUNDARYCONDITION_HPP
#define BOUNDARYCONDITION_HPP

#include <cstddef>

struct DataLayer;

enum class Side { Min = 0, Max = 1 };

/**
 * @class BoundaryCondition
 * @brief Abstract base class for all boundary condition types.
 *
 * Provides a unified interface for applying boundary conditions to
 * a DataLayer object. Each derived class implements its own logic
 * for updating ghost cells depending on the physical type of the boundary.
 *
 * Examples of derived classes include:
 *  - PeriodicBoundary
 *  - ReflectiveBoundary
 *  - InletBoundary / OutletBoundary
 *  - SymmetryBoundary
 *  - FreeStreamBoundary
 *
 * @note Boundary conditions are typically created and managed
 *       through BoundaryManager rather than instantiated directly.
 */
class BoundaryCondition {
public:
    /**
     * @brief Virtual destructor for safe polymorphic deletion.
     */
    virtual ~BoundaryCondition() = default;

    /**
     * @brief Applies the boundary condition to a given DataLayer along a specific axis.
     *
     * This method updates ghost cells corresponding to the specified side (Min or Max)
     * of the computational domain.
     *
     * @param layer Reference to the current data layer being modified.
     * @param axis Spatial axis along which the condition is applied (0 for 1D).
     * @param side Which boundary to apply: Side::Min (left) or Side::Max (right).
     *
     * @note Derived classes define the concrete logic of this method.
     */
    virtual void Apply(DataLayer &layer, int axis, Side side) const = 0;
};

#endif // BOUNDARYCONDITION_HPP
