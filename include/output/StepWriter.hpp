#ifndef STEPWRITER_HPP
#define STEPWRITER_HPP

#include <cstddef>
#include <string>

#include "config/Settings.hpp"

class DataLayer;
class Mesh;

/**
 * @class StepWriter
 * @brief Abstract interface for writing simulation data to disk.
 *
 * This class defines a minimal I/O interface that allows saving
 * the state of the computational domain at each time step.
 *
 * Concrete implementations define specific file formats and data serialization
 * strategies.
 *
 * Writers that accumulate data over multiple steps should override Finalize().
 */
class StepWriter {
public:
    virtual ~StepWriter() = default;

    /**
     * @brief Write simulation data to disk.
     *
     * @param layer The numerical solution data layer.
     * @param mesh Structured mesh with geometry and ranges.
     * @param settings Solver settings for output file name construction.
     * @param step Current simulation step number.
     * @param time Current simulation time.
     */
    virtual void Write(const DataLayer& layer,
                       const Mesh& mesh,
                       const Settings& settings,
                       std::size_t step,
                       double time) const = 0;

    /**
     * @brief Write simulation data with optional analytical comparison.
     *
     * Default implementation ignores analytical data and delegates to Write().
     *
     * @param layer The numerical solution data layer.
     * @param analytical_layer Optional analytical solution data.
     * @param mesh Structured mesh for numerical solution.
     * @param analytical_mesh Optional analytical mesh.
     * @param settings Solver settings for output file name construction.
     * @param step Current simulation step number.
     * @param time Current simulation time.
     */
    virtual void Write(const DataLayer& layer,
                       const DataLayer* analytical_layer,
                       const Mesh& mesh,
                       const Mesh* analytical_mesh,
                       const Settings& settings,
                       std::size_t step,
                       double time) const {
        (void)analytical_layer;
        (void)analytical_mesh;
        Write(layer, mesh, settings, step, time);
    }

    /**
     * @brief Finalize output and write any accumulated data.
     *
     * @param settings Solver settings for filename construction.
     * @return Path to the generated output file.
     */
    virtual auto Finalize(const Settings& settings) -> std::string {
        (void)settings;
        return "";
    }

    /**
     * @brief Whether this writer requires Finalize() to be called.
     * @return true if Finalize() should be called at end of simulation.
     */
    [[nodiscard]] virtual auto RequiresFinalization() const -> bool {
        return false;
    }
};

#endif  // STEPWRITER_HPP
