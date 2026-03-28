#ifndef TIMESTEPCALCULATOR_HPP
#define TIMESTEPCALCULATOR_HPP

class DataLayer;
class Mesh;

/**
 * @class TimeStepCalculator
 * @brief CFL-based timestep selection for generic face-based finite-volume meshes.
 *
 * For each cell, timestep is estimated from face contributions:
 *
 *   dt_cell = cfl * V / sum_faces( (|u_n| + c) * S )
 *
 * where:
 * - V is cell volume (or area in 2D)
 * - S is face measure
 * - u_n is velocity projected onto the face normal
 * - c is sound speed
 *
 * Returns the minimum stable timestep over all mesh cells.
 */
class TimeStepCalculator final {
public:
    /**
     * @brief Compute stable explicit timestep.
     *
     * @param layer Conservative state storage U(cell,var).
     * @param mesh Mesh with cells, faces, geometry, and connectivity.
     * @param gamma Ratio of specific heats.
     * @param cfl CFL number.
     * @return Stable timestep, or 0.0 if it cannot be computed.
     */
    static double ComputeDt(const DataLayer& layer,
                            const Mesh& mesh,
                            double gamma,
                            double cfl);
};

#endif  // TIMESTEPCALCULATOR_HPP
