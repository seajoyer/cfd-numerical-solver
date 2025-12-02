#ifndef VTKWRITER_HPP
#define VTKWRITER_HPP

#include <memory>
#include <string>

#include "StepWriter.hpp"

/**
 * @class VTKWriter
 * @brief Concrete implementation of StepWriter for VTK format output.
 *
 * This class writes simulation data to VTK structured grid files.
 * It supports 1D, 2D, and 3D data output with automatic conversion
 * from DataLayer's internal representation to VTK format.
 *
 * Output structure:
 * ```
 * output_dir/
 *   ├── analytical/              (if analytical solution enabled)
 *   │   └── step_NNNN.vtk
 *   └── solver__R_recon__N_size__CFL_value/
 *       └── solver__...__step_NNNN.vtk
 * ```
 */
class VTKWriter : public StepWriter {
   public:
    /**
     * @brief Constructs a VTKWriter with specified output directory.
     * @param output_dir Directory where VTK files will be written
     * @param is_analytical If true, writes to analytical subdirectory with simpler naming
     */
    explicit VTKWriter(const std::string& output_dir, bool is_analytical = false);
    ~VTKWriter() override;

    /**
     * @brief Writes the current state of DataLayer to a VTK file.
     * @param layer The data layer containing simulation state
     * @param settings Solver settings for output file name construction
     * @param step Current simulation step number
     * @param time Current simulation time
     */
    void Write(const DataLayer& layer, const Settings& settings, std::size_t step,
               double time) const override;

   private:
    std::string output_dir_;
    bool is_analytical_;

    // PIMPL to hide VTK implementation details
    class Impl;
    std::unique_ptr<Impl> pimpl_;

    /**
     * @brief Generates filename based on grid size and step number.
     * @param N Grid size
     * @param step Step number
     * @param settings Solver settings for output file name construction
     * @return Full path to output file
     */
    [[nodiscard]] auto GenerateFilename(int N, std::size_t step,
                                        const Settings& settings) const -> std::string;

    /**
     * @brief Writes 1D data to VTK format.
     */
    void Write1D(const DataLayer& layer, const Settings& settings, std::size_t step,
                 double time) const;

    /**
     * @brief Writes 2D data to VTK format (future implementation).
     */
    void Write2D(const DataLayer& layer, const Settings& settings, std::size_t step,
                 double time) const;

    /**
     * @brief Writes 3D data to VTK format (future implementation).
     */
    void Write3D(const DataLayer& layer, const Settings& settings, std::size_t step,
                 double time) const;
};

#endif  // VTKWRITER_HPP
