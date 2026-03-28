#ifndef STEPWRITER_HPP
#define STEPWRITER_HPP

#include <cstddef>

class DataLayer;
class Settings;
class Mesh;

/**
 * @class StepWriter
 * @brief Abstract interface for writing simulation data to disk.
 *
 * This interface defines minimal output functionality for one simulation state.
 * Concrete implementations define specific file formats and serialization logic.
 *
 * Writers that accumulate data over multiple steps may override Finalize().
 */
class StepWriter {
public:
    virtual ~StepWriter() = default;

    /**
     * @brief Write simulation data to disk.
     *
     * @param layer Numerical solution data.
     * @param mesh Generic mesh with geometry and connectivity.
     * @param settings Solver settings used for metadata and file naming.
     * @param step Current simulation step number.
     * @param time Current simulation time.
     */
    virtual void Write(const DataLayer& layer,
                       const Mesh& mesh,
                       const Settings& settings,
                       std::size_t step,
                       double time) const = 0;

    /**
     * @brief Finalize output and write any accumulated data.
     *
     * Default implementation does nothing.
     */
    virtual void Finalize(const Settings& settings) {
        (void)settings;
    }

    /**
     * @brief Whether this writer requires Finalize() to be called.
     */
    [[nodiscard]] virtual bool RequiresFinalization() const {
        return false;
    }
};

#endif  // STEPWRITER_HPP
