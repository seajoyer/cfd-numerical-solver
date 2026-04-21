#ifndef BOUNDARYCONDITION_HPP
#define BOUNDARYCONDITION_HPP

#include <cstdint>
#include <stdexcept>

class DataLayer;
class Mesh;
class PressureVelocityState;
class PressureVelocityWorkspace;

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
 * @enum PvAssemblyStage
 * @brief Assembly stage for pressure-velocity boundary contributions.
 */
enum class PvAssemblyStage : std::uint8_t {
    MomentumCoefficients = 0,
    PressureCorrectionEquation = 1
};
// CHECK: TG_BC
/**
 * @class BoundaryCondition
 * @brief Abstract base class for physical boundary conditions on a structured grid.
 */
class BoundaryCondition {
public:
    virtual ~BoundaryCondition() = default;

    /**
     * @brief Apply BC to conservative state storage.
     */
    virtual void Apply(DataLayer& layer, const Mesh& mesh, Axis axis, Side side) const = 0;

    /**
     * @brief Apply BC to pressure-velocity staggered state.
     */
    virtual void Apply(PressureVelocityState& state,
                       const Mesh& mesh,
                       Axis axis,
                       Side side) const {
        (void)state;
        (void)mesh;
        (void)axis;
        (void)side;
        throw std::runtime_error(
            "BoundaryCondition: pressure-velocity Apply() is not implemented for this BC");
    }

    /**
     * @brief Apply BC-specific contributions to pressure-velocity equation assembly.
     *
     * Default implementation does nothing.
     */
    virtual void ApplyPressureVelocityBoundary(PressureVelocityState& state,
                                              PressureVelocityWorkspace& workspace,
                                              const Mesh& mesh,
                                              Axis axis,
                                              Side side,
                                              PvAssemblyStage stage,
                                              bool steady,
                                              double dt,
                                              double nu) const {
        (void)state;
        (void)workspace;
        (void)mesh;
        (void)axis;
        (void)side;
        (void)stage;
        (void)steady;
        (void)dt;
        (void)nu;
    }

    /**
     * @brief Whether this BC represents periodic topology.
     *
     * @details
     * Used by pressure-velocity solvers to distinguish
     * "no local neighbor because this is a real boundary"
     * from "no local neighbor because this is a periodic seam".
     */
    [[nodiscard]] virtual bool IsPeriodic() const {
        return false;
    }
};

#endif  // BOUNDARYCONDITION_HPP