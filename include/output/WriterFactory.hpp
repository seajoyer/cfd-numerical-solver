#ifndef WRITERFACTORY_H_
#define WRITERFACTORY_H_

#include <memory>
#include <string>
#include <vector>

#include "output/StepWriter.hpp"

/**
 * @class WriterFactory
 * @brief Factory class for creating StepWriter instances based on output format.
 *
 * This class provides static factory methods to instantiate concrete StepWriter
 * implementations (e.g., VTKWriter, PNGWriter) without exposing the underlying 
 * creation logic. It promotes extensibility by allowing new output formats to be 
 * added centrally without modifying client code (e.g., Simulation).
 *
 * Supported formats:
 * - "vtk": VTK structured grid format
 * - "png" or "png<width>x<height>": PNG image with 4 subplot panels
 *   - Default: png (1200x900)
 *   - Custom: png1920x1080, png800x600, png3840x2160, etc.
 *   - Fonts and line widths scale automatically with resolution
 */
class WriterFactory {
   public:
    /**
     * @brief Creates a StepWriter instance for the specified output format.
     *
     * @param output_format String identifier for the desired output format (e.g., "vtk", "png").
     * @param output_dir Directory path where output files will be written.
     * @param is_analytical If true, creates writer for analytical solution output.
     * @return Unique pointer to the created StepWriter.
     * @throws std::runtime_error If the output format is unrecognized.
     */
    static auto Create(const std::string& output_format, 
                       const std::string& output_dir,
                       bool is_analytical = false)
        -> std::unique_ptr<StepWriter>;

    /**
     * @brief Creates multiple StepWriter instances for the specified output formats.
     *
     * @param output_formats Vector of format identifiers (e.g., {"vtk", "png"}).
     * @param output_dir Base directory path for output files.
     * @param is_analytical If true, creates writers for analytical solution output.
     * @return Vector of unique pointers to created StepWriters.
     * @throws std::runtime_error If any output format is unrecognized.
     */
    static auto CreateMultiple(const std::vector<std::string>& output_formats,
                               const std::string& output_dir,
                               bool is_analytical = false)
        -> std::vector<std::unique_ptr<StepWriter>>;

    /**
     * @brief Check if a format is supported
     * @param format Format string to check
     * @return true if the format is supported
     */
    static auto IsFormatSupported(const std::string& format) -> bool;

    /**
     * @brief Get list of all supported formats
     * @return Vector of supported format strings
     */
    static auto GetSupportedFormats() -> std::vector<std::string>;
};

#endif  // WRITERFACTORY_H_
