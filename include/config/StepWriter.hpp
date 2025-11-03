#ifndef STEPWRITER_HPP
#define STEPWRITER_HPP

#include <cstddef>

struct DataLayer;

/**
 * @class StepWriter
 * @brief Abstract interface for writing simulation data to disk.
 *
 * This class defines a minimal I/O interface that allows saving
 * the state of the computational domain at each time step.
 *
 * Concrete implementations (e.g., VtkWriter) define specific file
 * formats and data serialization strategies.
 *
 * @note StepWriter is used by Simulation to record results periodically
 *       without depending on the underlying file format.
 */
class StepWriter {
public:
    virtual ~StepWriter() = default;

    virtual void Write(const DataLayer &layer, std::size_t step, double time) const = 0;
};

#endif // STEPWRITER_HPP