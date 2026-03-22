#ifndef INLETBOUNDARY_HPP
#define INLETBOUNDARY_HPP

#include "bc/BoundaryCondition.hpp"
#include "bc/BoundaryFactory.hpp"

/**
 * @class InletBoundary
 * @brief Conditional inlet boundary with prescribed conservative inflow state.
 *
 * If the local normal velocity at the nearest interior core cell is directed
 * into the domain, ghost cells are set to the prescribed inflow state.
 * Otherwise, behaves as outlet (zero-gradient): copies the nearest interior layer
 * into ghost cells.
 *
 * Works only with conservative state U(var,i,j,k).
 */
class InletBoundary final : public BoundaryCondition {
public:
    /**
     * @brief Constructs inlet boundary with prescribed inflow conservative state.
     * @param inflow_U Conservative state imposed during inflow.
     */
    explicit InletBoundary(const FarfieldConservative& inflow_U);

    /**
     * @brief Apply inlet BC along the specified axis and side.
     * @param layer Data layer to modify (ghost cells of U will be written).
     * @param mesh Structured mesh with ranges and metadata.
     * @param axis Axis (X/Y/Z).
     * @param side Side (Left/Right).
     */
    void Apply(DataLayer& layer, const Mesh& mesh, Axis axis, Side side) const override;

private:
    FarfieldConservative inflow_U_;
};

#endif  // INLETBOUNDARY_HPP