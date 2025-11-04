#ifndef WRITERFACTORY_H_
#define WRITERFACTORY_H_

#include <memory>
#include <string>

#include "output/StepWriter.hpp"

/**
 * @class WriterFactory
 * @brief Factory class for creating StepWriter instances based on output format.
 *
 * This class provides a static factory method to instantiate concrete StepWriter
 * implementations (e.g., VTKWriter) without exposing the underlying creation logic.
 * It promotes extensibility by allowing new output formats to be added centrally
 * without modifying client code (e.g., Simulation).
 *
 * @note Currently supports "vtk" format; additional formats can be added in the
 *       implementation file.
 */
class WriterFactory {
   public:
    /**
     * @brief Creates a StepWriter instance for the specified output format.
     *
     * @param output_format String identifier for the desired output format (e.g., "vtk").
     * @param output_dir Directory path where output files will be written.
     * @return Unique pointer to the created StepWriter.
     * @throws std::runtime_error If the output format is unrecognized.
     */
    static auto Create(const std::string& output_format, const std::string& output_dir)
        -> std::unique_ptr<StepWriter>;
};

#endif  // WRITERFACTORY_H_
