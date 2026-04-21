#ifndef TIMESTEPCALCULATOR_HPP
#define TIMESTEPCALCULATOR_HPP
#include <memory>

class DataLayer;
class Mesh;
struct Settings;
class EOS;

/**
 * @class TimeStepCalculator
 * @brief CFL-based timestep selection from conservative state U on structured mesh.
 *
 * Uses conservative U(var,i,j,k) and the active Equation of State (EOS)
 * to compute the correct local sound speed for physical timestep restriction.
 */
class TimeStepCalculator final {
public:
    // CHECK: CFL_2D
    // CHECK: FLIC_CFL
    /**
     * @brief Computes stable explicit timestep dt.
     * @param layer DataLayer with conservative state U.
     * @param mesh Structured mesh with geometry, metrics, and cell types.
     * @param settings Global simulation settings (contains CFL and chemistry flags).
     * @param eos Initialized Equation of State object to evaluate true sound speed.
     * @return dt > 0 on success; 0.0 if dt cannot be computed.
     */
    static auto ComputeDt(const DataLayer& layer,
                          const Mesh& mesh,
                          const Settings& settings,
                          std::shared_ptr<EOS> eos) -> double;
};

#endif  // TIMESTEPCALCULATOR_HPP
