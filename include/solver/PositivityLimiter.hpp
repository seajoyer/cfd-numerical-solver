#ifndef POSITIVITYLIMITER_HPP
#define POSITIVITYLIMITER_HPP

class DataLayer;
class Mesh;

/**
 * @class PositivityLimiter
 * @brief Ensures physically admissible conservative Euler state U (rho and P floors).
 *
 * Applies in-place corrections on conservative U(var,i,j,k):
 *  - rho is clamped to rho_min (preserving velocities when possible)
 *  - pressure P(U) is clamped to p_min by adjusting total energy E
 *
 * Corrections are applied only on fluid core cells.
 */
class PositivityLimiter final {
public:
    /**
     * @brief Apply positivity corrections in-place on fluid core cells.
     * @param layer DataLayer containing U.
     * @param mesh Structured mesh with core ranges and cell types.
     * @param gamma Ratio of specific heats.
     * @param rho_min Minimum density.
     * @param p_min Minimum pressure.
     */
    static void Apply(DataLayer& layer, const Mesh& mesh, double gamma, double rho_min, double p_min);
};

#endif  // POSITIVITYLIMITER_HPP