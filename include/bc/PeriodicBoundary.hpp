#ifndef PERIODICBOUNDARY_HPP
#define PERIODICBOUNDARY_HPP

#include "bc/BoundaryCondition.hpp"

/**
 * @class PeriodicBoundary
 * @brief Periodic boundary condition for conservative state U.
 *
 * Ghost cells are filled from the opposite side of the local core domain.
 *
 * This implementation is intended for:
 * - single-process runs
 * - MPI runs with a single rank
 *
 * For multi-rank MPI, periodicity should be handled by domain decomposition
 * and halo exchange, not by physical boundary application.
 */
class PeriodicBoundary final : public BoundaryCondition {
public:
    void Apply(DataLayer& layer, const Mesh& mesh, Axis axis, Side side) const override;
};

#endif  // PERIODICBOUNDARY_HPP