#ifndef GIFWRITER_HPP
#define GIFWRITER_HPP

#include <memory>
#include <string>

#include "StepWriter.hpp"

/**
 * @class GIFWriter
 * @brief Concrete implementation of StepWriter for animated GIF output.
 *
 * This class generates animated GIFs by capturing frames during
 * the simulation and assembling them when Finalize() is called.
 *
 * Currently supports 1D data only. For 2D simulations, Write() logs a
 * warning and skips (use VTK format for 2D output).
 *
 * Requires Finalize() to be called at the end of simulation.
 */
class GIFWriter : public StepWriter {
   public:
    explicit GIFWriter(const std::string& output_dir,
                       int width = 1200,
                       int height = 900,
                       int delay_centiseconds = 10);

    ~GIFWriter() override;

    GIFWriter(GIFWriter&&) noexcept;
    auto operator=(GIFWriter&&) noexcept -> GIFWriter&;

    /**
     * @brief Captures current state as a GIF frame.
     * For 2D data: logs a warning and skips.
     */
    void Write(const DataLayer& layer, const Settings& settings, std::size_t step,
               double time) const override;

    /**
     * @brief Captures current state with analytical comparison as a frame.
     * For 2D data: logs a warning and skips.
     */
    void Write(const DataLayer& layer, const DataLayer* analytical_layer,
               const Settings& settings, std::size_t step, double time) const override;

    /**
     * @brief Writes all accumulated frames to the final GIF file.
     */
    auto Finalize(const Settings& settings) -> std::string override;

    [[nodiscard]] auto RequiresFinalization() const -> bool override {
        return true;
    }

    [[nodiscard]] auto GetFrameCount() const -> std::size_t;
    void SetFrameDelay(int centiseconds);
    [[nodiscard]] auto GetFrameDelay() const -> int;

   private:
    std::string output_dir_;
    int width_;
    int height_;
    int delay_centiseconds_;

    class Impl;
    std::unique_ptr<Impl> pimpl_;

    [[nodiscard]] auto GenerateFilename(const Settings& settings) const -> std::string;
    void RenderFrame(const DataLayer& layer, const DataLayer* analytical_layer,
                     const Settings& settings, std::size_t step, double time) const;
};

#endif  // GIFWRITER_HPP
