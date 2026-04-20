#ifndef STEPWRITER_HPP
#define STEPWRITER_HPP

#include <cstddef>
#include <stdexcept>
#include <string>

#include "config/Settings.hpp"

class DataLayer;
class Mesh;
class EOS;
class PressureVelocityState;

/**
 * @class StepWriter
 * @brief Abstract interface for writing simulation data to disk.
 *
 * Supports both:
 *  - conservative DataLayer output
 *  - pressure-velocity state output
 *
 * Writers that accumulate data over multiple steps should override Finalize().
 */
class StepWriter {
public:
    virtual ~StepWriter() = default;

    /**
     * @brief Write simulation data to disk from conservative storage.
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
     * @brief Write simulation data from pressure-velocity staggered storage.
     *
     * Default implementation throws. Concrete writers that support
     * SIMPLE / PISO / PIMPLE output should override this method.
     */
    virtual void Write(const PressureVelocityState& state,
                       const Mesh& mesh,
                       const Settings& settings,
                       std::size_t step,
                       double time) const {
        (void)state;
        (void)mesh;
        (void)settings;
        (void)step;
        (void)time;
        throw std::runtime_error(
            "StepWriter: PressureVelocityState writing is not implemented");
    }

    /**
     * @brief Finalize output and write any accumulated data.
     *
     * @param settings Solver settings for filename construction.
     */
    virtual void Finalize(const Settings& settings) {
        (void)settings;
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
