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
 * Resolution can be customized and all visual elements (fonts, line widths,
 * margins) automatically scale based on the specified dimensions.
 *
 * Currently supports 1D data only. For 2D simulations, Write() logs a
 * warning and skips (use VTK format for 2D output).
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
     * @brief Writes simulation data to a PNG file.
     * 
     * For 1D data: creates a 2x2 grid of plots.
     * For 2D data: logs a warning and skips (not yet supported).
     */
    void Write(const DataLayer& layer, const Settings& settings, std::size_t step,
               double time) const override;

    /**
     * @brief Writes simulation data with analytical comparison.
     * 
     * For 1D data: creates plots with both numerical and analytical solutions.
     * For 2D data: logs a warning and skips.
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

    [[nodiscard]] auto GenerateFilename(std::size_t step) const -> std::string;
};

#endif  // PNGWRITER_HPP
