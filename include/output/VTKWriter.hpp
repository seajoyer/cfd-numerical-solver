#ifndef VTKWRITER_HPP
#define VTKWRITER_HPP

#include "StepWriter.hpp"
#include <string>
#include <memory>

/**
 * @class VTKWriter
 * @brief Concrete implementation of StepWriter for VTK format output.
 *
 * This class writes simulation data to VTK structured grid files.
 * It supports 1D, 2D, and 3D data output with automatic conversion
 * from DataLayer's internal representation to VTK format.
 *
 * Files are named as: output_dir/N_<gridsize>__step_<stepnum>.vtk
 */
class VTKWriter : public StepWriter {
public:
    /**
     * @brief Constructs a VTKWriter with specified output directory.
     * @param output_dir Directory where VTK files will be written
     */
    explicit VTKWriter(const std::string& output_dir);
    ~VTKWriter() override;

    /**
     * @brief Writes the current state of DataLayer to a VTK file.
     * @param layer The data layer containing simulation state
     * @param step Current simulation step number
     * @param time Current simulation time
     */
    void Write(const DataLayer& layer, std::size_t step, double time) const override;

private:
    std::string output_dir_;
    
    // PIMPL to hide VTK implementation details
    class Impl;
    std::unique_ptr<Impl> pimpl_;

    /**
     * @brief Generates filename based on grid size and step number.
     * @param N Grid size
     * @param step Step number
     * @return Full path to output file
     */
    [[nodiscard]] auto GenerateFilename(int N, std::size_t step) const -> std::string;

    /**
     * @brief Writes 1D data to VTK format.
     */
    void Write1D(const DataLayer& layer, std::size_t step, double time) const;

    /**
     * @brief Writes 2D data to VTK format (future implementation).
     */
    void Write2D(const DataLayer& layer, std::size_t step, double time) const;

    /**
     * @brief Writes 3D data to VTK format (future implementation).
     */
    void Write3D(const DataLayer& layer, std::size_t step, double time) const;
};

#endif // VTKWRITER_HPP
