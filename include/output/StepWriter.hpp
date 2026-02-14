#ifndef STEPWRITER_HPP
#define STEPWRITER_HPP

#include <cstddef>
#include <string>
#include "config/Settings.hpp"

struct DataLayer;

/**
 * @class StepWriter
 * @brief Abstract interface for writing simulation data to disk.
 *
 * This class defines a minimal I/O interface that allows saving
 * the state of the computational domain at each time step.
 *
 * Concrete implementations (e.g., VTKWriter, PNGWriter, GIFWriter) define 
 * specific file formats and data serialization strategies.
 *
 * Writers that accumulate data over multiple steps (e.g., GIFWriter) should
 * override the Finalize() method to perform final output operations.
 *
 * @note StepWriter is used by Simulation to record results periodically
 *       without depending on the underlying file format.
 */
class StepWriter {
   public:
    virtual ~StepWriter() = default;

    /**
     * @brief Write 1D simulation data to disk
     * 
     * @param layer The numerical solution data layer
     * @param settings Solver settings for output file name construction
     * @param step Current simulation step number
     * @param time Current simulation time
     */
    virtual void Write(const DataLayer& layer, const Settings& settings, std::size_t step,
                       double time) const = 0;

    /**
     * @brief Write simulation data with optional analytical comparison
     * 
     * This overload allows writers to include analytical solution data
     * for comparison plots (e.g., PNG writer, GIF writer).
     * 
     * @param layer The numerical solution data layer
     * @param analytical_layer Optional analytical solution data (may be nullptr)
     * @param settings Solver settings for output file name construction
     * @param step Current simulation step number
     * @param time Current simulation time
     */
    virtual void Write(const DataLayer& layer, const DataLayer* analytical_layer,
                       const Settings& settings, std::size_t step, double time) const {
        // Default implementation ignores analytical data
        Write(layer, settings, step, time);
    }

    /**
     * @brief Write 2D simulation data to disk
     * 
     * @param layer The numerical solution data layer
     * @param settings Solver settings for output file name construction
     * @param step Current simulation step number
     * @param time Current simulation time
     */
    virtual void Write2D(const DataLayer& layer, const Settings& settings, std::size_t step,
                       double time) const = 0;


    /**
     * @brief Finalize output and write any accumulated data
     * 
     * This method is called at the end of simulation to allow writers
     * that accumulate data (like GIFWriter) to perform final output.
     * 
     * @param settings Solver settings for output file name construction
     * @return Path to the finalized output file (empty if not applicable)
     */
    virtual auto Finalize(const Settings& settings) -> std::string {
        // Default implementation does nothing
        (void)settings;
        return "";
    }
    
    /**
     * @brief Check if this writer requires finalization
     * 
     * @return true if Finalize() should be called at end of simulation
     */
    [[nodiscard]] virtual auto RequiresFinalization() const -> bool {
        return false;
    }
};

#endif  // STEPWRITER_HPP
