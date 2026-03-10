#ifndef BOUNDARYCONDITION_HPP
#define BOUNDARYCONDITION_HPP

#include <cstdint>

class DataLayer;
class Mesh;

/**
 * @enum Side
 * @brief Boundary side along a selected axis.
 *
 * The boundary is defined by a pair (axis, side):
 * - Side::Left  : lower/min side of the axis
 * - Side::Right : upper/max side of the axis
 */
enum class Side : std::uint8_t { Left = 0, Right = 1 };

/**
 * @brief Spatial axis.
 *
 * Axis is defined in data/Variables.hpp (project-wide enum).
 * Forward-declared here to avoid heavy includes.
 */
enum class Axis : std::uint8_t;

/**
 * @class BoundaryCondition
 * @brief Abstract base class for physical boundary conditions on a structured grid.
 *
 * Boundary conditions fill ghost cells of the conservative state U(var,i,j,k)
 * for a given (axis, side).
 *
 * Contract:
 *  - Must modify only ghost cells.
 *  - Must not modify core cells.
 *  - Intended for global external boundaries (not for MPI internal interfaces).
 */
class BoundaryCondition {
public:
    /** @brief Virtual destructor for safe polymorphic deletion. */
    virtual ~BoundaryCondition() = default;

    /**
     * @brief Apply boundary condition along a specified axis and side.
     *
     * @param layer Data layer to modify (ghost cells of U will be written).
     * @param mesh Structured mesh with ranges and metadata.
     * @param axis Spatial axis (X,Y,Z).
     * @param side Boundary side (Left or Right).
     */
    virtual void Apply(DataLayer& layer, const Mesh& mesh, Axis axis, Side side) const = 0;
};

#endif  // BOUNDARYCONDITION_HPP
