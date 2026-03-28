#ifndef POSITIVITYLIMITER_HPP
#define POSITIVITYLIMITER_HPP

class DataLayer;
class Mesh;

/**
 * @class PositivityLimiter
 * @brief Ensures physically admissible conservative Euler state on generic meshes.
 *
 * Applies in-place corrections on conservative cell-centered state U(cell,var):
 * - density is clamped to rho_min while preserving velocity when possible
 * - pressure is clamped to p_min by adjusting total energy
 */
class PositivityLimiter final {
public:
    /**
     * @brief Apply positivity corrections in-place on all mesh cells.
     *
     * @param layer DataLayer containing conservative state U(cell,var).
     * @param mesh Mesh defining the number of real cells.
     * @param gamma Ratio of specific heats.
     * @param rho_min Density floor.
     * @param p_min Pressure floor.
     */
    static void Apply(DataLayer& layer,
                      const Mesh& mesh,
                      double gamma,
                      double rho_min,
                      double p_min);
};

#endif  // POSITIVITYLIMITER_HPP