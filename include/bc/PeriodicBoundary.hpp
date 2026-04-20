#ifndef PERIODICBOUNDARY_HPP
#define PERIODICBOUNDARY_HPP

#include "bc/BoundaryCondition.hpp"

/**
 * @class PeriodicBoundary
 * @brief Periodic boundary condition for conservative and pressure-velocity storage.
 *
 * Ghost cells / ghost faces are filled from the opposite side of the local core domain.
 *
 * Intended for:
 *  - single-process runs
 *  - fallback local periodic wrapping when no MPI periodic neighbor exchange is used
 *
 * For multi-rank MPI periodic runs, periodicity should be handled by domain
 * decomposition neighbor ranks and HaloExchange.
 */
class PeriodicBoundary final : public BoundaryCondition {
public:
    void Apply(DataLayer& layer, const Mesh& mesh, Axis axis, Side side) const override;

    void Apply(PressureVelocityState& state,
               const Mesh& mesh,
               Axis axis,
               Side side) const override;

    void ApplyPressureVelocityBoundary(PressureVelocityState& state,
                                       PressureVelocityWorkspace& workspace,
                                       const Mesh& mesh,
                                       Axis axis,
                                       Side side,
                                       PvAssemblyStage stage,
                                       bool steady,
                                       double dt,
                                       double nu) const override;

    [[nodiscard]] bool IsPeriodic() const override;
};

#endif  // PERIODICBOUNDARY_HPP
