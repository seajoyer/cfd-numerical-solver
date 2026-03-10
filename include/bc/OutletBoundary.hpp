#ifndef OUTLETBOUNDARY_HPP
#define OUTLETBOUNDARY_HPP

#include "bc/BoundaryCondition.hpp"

/**
 * @class OutletBoundary
 * @brief Outlet (zero-gradient / Neumann) boundary condition for conservative state U.
 *
 * The outlet condition copies the nearest interior (core) layer into ghost layers
 * along the selected axis, enforcing zero normal gradient at the boundary.
 *
 * Works only with conservative state U(var,i,j,k).
 */
class OutletBoundary final : public BoundaryCondition {
public:
    /**
     * @brief Apply outlet BC along the specified axis and side.
     * @param layer Data layer to modify (ghost cells of U will be written).
     * @param mesh Structured mesh with ranges and metadata.
     * @param axis Axis (X/Y/Z).
     * @param side Side (Left/Right).
     */
    void Apply(DataLayer& layer, const Mesh& mesh, Axis axis, Side side) const override;
};

#endif  // OUTLETBOUNDARY_HPP