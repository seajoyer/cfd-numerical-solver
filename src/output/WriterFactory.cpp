#include "output/WriterFactory.hpp"

#include <stdexcept>

#include "output/VTKWriter.hpp"

/**
 * @namespace output
 * @brief Namespace containing classes and factories for simulation output handling.
 *
 * This namespace encapsulates all I/O-related components, ensuring separation of
 * concerns between simulation logic and file writing.
 */

// Static member function implementation
auto WriterFactory::Create(const std::string& output_format,
                           const std::string& output_dir) -> std::unique_ptr<StepWriter> {
    if (output_format == "vtk") {
        /**
         * @class VTKWriter
         * @brief Concrete StepWriter implementation for VTK file format.
         *
         * Writes simulation data (e.g., density, velocity, pressure) to VTK files,
         * which are compatible with visualization tools like ParaView or VisIt.
         * Supports 1D data layers by generating legacy VTK unstructured grid files.
         *
         * @note Assumes DataLayer provides cell-centered data (xc, xb) and
         * primitive/conservative variables.
         */
        return std::make_unique<VTKWriter>(output_dir);
    }

    // Future formats can be added here:
    //
    // if (output_format == "plt") {
    //     /**
    //      * @class PLTWriter
    //      * @brief Concrete StepWriter for Tecplot PLT format.
    //      *
    //      * ... (documentation for future class)
    //      */
    //     return std::make_unique<PLTWriter>(output_dir);
    // }
    //
    // etc ...

    throw std::runtime_error("Unknown output format: " + output_format);
}
