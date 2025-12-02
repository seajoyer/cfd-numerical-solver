#ifndef PNGWRITER_HPP
#define PNGWRITER_HPP

#include <memory>
#include <string>

#include "StepWriter.hpp"

/**
 * @class PNGWriter
 * @brief Concrete implementation of StepWriter for PNG image output.
 *
 * This class generates PNG images containing 4 subplots showing:
 * - Density vs X (top-left)
 * - Velocity vs X (top-right)  
 * - Pressure vs X (bottom-left)
 * - Specific Internal Energy vs X (bottom-right)
 *
 * Numerical solutions are plotted as red lines, analytical solutions
 * (when provided) are plotted as black lines.
 *
 * Uses VTK's charting capabilities (vtkChartXY, vtkContextView) for
 * rendering 2D plots to off-screen images.
 *
 * Files are named as: output_dir/step_<stepnum>.png
 */
class PNGWriter : public StepWriter {
   public:
    /**
     * @brief Constructs a PNGWriter with specified output directory.
     * @param output_dir Directory where PNG files will be written
     * @param width Image width in pixels (default: 1200)
     * @param height Image height in pixels (default: 900)
     */
    explicit PNGWriter(const std::string& output_dir, 
                       int width = 1200, 
                       int height = 900);
    
    ~PNGWriter() override;

    /**
     * @brief Writes the current state of DataLayer to a PNG file.
     * 
     * Creates a 2x2 grid of plots showing density, velocity, pressure,
     * and specific internal energy as functions of X coordinate.
     * Only numerical solution is plotted (red line).
     * 
     * @param layer The data layer containing simulation state
     * @param settings Solver settings for output file name construction
     * @param step Current simulation step number
     * @param time Current simulation time
     */
    void Write(const DataLayer& layer, const Settings& settings, std::size_t step,
               double time) const override;

    /**
     * @brief Writes simulation data with analytical comparison.
     * 
     * Creates a 2x2 grid of plots showing density, velocity, pressure,
     * and specific internal energy. Numerical solution is plotted as
     * red line, analytical solution as black line.
     * 
     * @param layer The numerical solution data layer
     * @param analytical_layer The analytical solution data (may be nullptr)
     * @param settings Solver settings for output file name construction
     * @param step Current simulation step number
     * @param time Current simulation time
     */
    void Write(const DataLayer& layer, const DataLayer* analytical_layer,
               const Settings& settings, std::size_t step, double time) const override;

   private:
    std::string output_dir_;
    int width_;
    int height_;

    // PIMPL to hide VTK implementation details
    class Impl;
    std::unique_ptr<Impl> pimpl_;

    /**
     * @brief Generates filename based on step number.
     * @param step Step number
     * @return Full path to output file
     */
    [[nodiscard]] auto GenerateFilename(std::size_t step) const -> std::string;
};

#endif  // PNGWRITER_HPP
