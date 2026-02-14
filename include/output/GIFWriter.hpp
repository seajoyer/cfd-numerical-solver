#ifndef GIFWRITER_HPP
#define GIFWRITER_HPP

#include <memory>
#include <string>
#include <vector>

#include "StepWriter.hpp"

/**
 * @class GIFWriter
 * @brief Concrete implementation of StepWriter for animated GIF output.
 *
 * This class generates an animated GIF showing the evolution of the simulation
 * over time. Each frame contains 4 subplots showing:
 * - Density vs X (top-left)
 * - Velocity vs X (top-right)
 * - Pressure vs X (bottom-left)
 * - Specific Internal Energy vs X (bottom-right)
 *
 * Numerical solutions are plotted as red lines, analytical solutions
 * (when provided) are plotted as black lines.
 *
 * The GIF is written to the case directory with the naming convention:
 *   <solver>__R_<reconstruction>__N_<cells>__CFL_<cfl>.gif
 *
 * Resolution can be customized via the format string (e.g., "gif1920x1080").
 * Visual elements scale automatically based on resolution, using 1200x900
 * as the baseline.
 *
 * ## Implementation Details
 * 
 * The writer buffers frames during the simulation and generates the final
 * GIF when Finalize() is called. Each frame is rendered using VTK's charting
 * capabilities (same as PNGWriter) and stored in memory as raw pixel data.
 *
 * Uses the gif.h single-header library for GIF encoding.
 *
 * ## Memory Considerations
 *
 * Each frame uses (width * height * 4) bytes of memory for RGBA data.
 * For a 1200x900 GIF with 1000 frames, this would require ~4.3 GB.
 * For long simulations, consider:
 * - Reducing output frequency (output_every_steps)
 * - Using smaller resolution
 * - Using VTK format for full data and PNG for selected frames
 *
 * @see PNGWriter for similar single-frame output
 */
class GIFWriter : public StepWriter {
   public:
    /**
     * @brief Constructs a GIFWriter with specified output directory and resolution.
     * 
     * @param output_dir Base directory for output (GIF will be written here directly)
     * @param width Frame width in pixels (default: 1200)
     * @param height Frame height in pixels (default: 900)
     * @param delay_centiseconds Delay between frames in centiseconds (default: 10 = 0.1s)
     */
    explicit GIFWriter(const std::string& output_dir,
                       int width = 1200,
                       int height = 900,
                       int delay_centiseconds = 10);

    ~GIFWriter() override;

    // Prevent copying (writer holds frame buffer)
    GIFWriter(const GIFWriter&) = delete;
    auto operator=(const GIFWriter&) -> GIFWriter& = delete;

    // Allow moving
    GIFWriter(GIFWriter&&) noexcept;
    auto operator=(GIFWriter&&) noexcept -> GIFWriter&;

    /**
     * @brief Captures the current state as a frame for the GIF.
     * 
     * Renders the 2x2 grid of plots and stores the frame in memory.
     * Only numerical solution is included.
     * 
     * @param layer The data layer containing simulation state
     * @param settings Solver settings (used for filename generation in Finalize)
     * @param step Current simulation step number
     * @param time Current simulation time
     */
    void Write(const DataLayer& layer, const Settings& settings, std::size_t step,
               double time) const override;

    /**
     * @brief Captures the current state with analytical comparison as a frame.
     * 
     * Renders the 2x2 grid of plots with both numerical (red) and
     * analytical (black) solutions and stores the frame in memory.
     * 
     * @param layer The numerical solution data layer
     * @param analytical_layer The analytical solution data (may be nullptr)
     * @param settings Solver settings
     * @param step Current simulation step number
     * @param time Current simulation time
     */
    void Write(const DataLayer& layer, const DataLayer* analytical_layer,
               const Settings& settings, std::size_t step, double time) const override;

    /**
     * @brief Write 2D simulation data (not implemented for GIF writer).
     * 
     * GIFWriter currently only supports 1D data visualization.
     * This method is required by the StepWriter interface.
     * 
     * @param layer The data layer containing simulation state
     * @param settings Solver settings
     * @param step Current simulation step number
     * @param time Current simulation time
     */
    void Write2D(const DataLayer& layer, const Settings& settings, std::size_t step,
                 double time) const override;

    /**
     * @brief Writes all accumulated frames to the final GIF file.
     * 
     * This must be called at the end of simulation to generate the GIF.
     * The filename is constructed from the solver settings:
     *   <output_dir>/<solver>__R_<reconstruction>__N_<N>__CFL_<cfl>.gif
     * 
     * @param settings Solver settings for filename construction
     * @return Path to the generated GIF file
     */
    auto Finalize(const Settings& settings) -> std::string override;

    /**
     * @brief Indicates that this writer requires finalization.
     * @return Always returns true
     */
    [[nodiscard]] auto RequiresFinalization() const -> bool override {
        return true;
    }

    /**
     * @brief Returns the number of frames currently buffered.
     * @return Frame count
     */
    [[nodiscard]] auto GetFrameCount() const -> std::size_t;

    /**
     * @brief Sets the delay between frames.
     * @param centiseconds Delay in centiseconds (100ths of a second)
     */
    void SetFrameDelay(int centiseconds);

    /**
     * @brief Gets the current frame delay.
     * @return Delay in centiseconds
     */
    [[nodiscard]] auto GetFrameDelay() const -> int;

   private:
    std::string output_dir_;
    int width_;
    int height_;
    int delay_centiseconds_;

    // PIMPL to hide implementation details and frame buffer
    class Impl;
    std::unique_ptr<Impl> pimpl_;

    /**
     * @brief Generates the output filename based on settings.
     * @param settings Solver settings
     * @return Full path to the GIF file
     */
    [[nodiscard]] auto GenerateFilename(const Settings& settings) const -> std::string;

    /**
     * @brief Renders a single frame and adds it to the buffer.
     * 
     * @param layer Numerical solution data
     * @param analytical_layer Analytical solution data (may be nullptr)
     * @param settings Solver settings
     * @param step Step number for title
     * @param time Current time for title
     */
    void RenderFrame(const DataLayer& layer, const DataLayer* analytical_layer,
                     const Settings& settings, std::size_t step, double time) const;
};

#endif  // GIFWRITER_HPP
