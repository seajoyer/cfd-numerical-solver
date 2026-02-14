#ifndef VTKWRITER_HPP
#define VTKWRITER_HPP

#include <memory>
#include <string>
#include "output/StepWriter.hpp"
#include "data/DataLayer.hpp"
#include "config/Settings.hpp"

/**
 * @class VTKWriter
 * @brief VTK structured grid output for 1D, 2D, and 3D simulations.
 *
 * Write() checks layer.GetDim() and dispatches internally:
 * - 1D: writes VTK structured grid (expanded to 2D slab for visualization)
 *       with density, velocity, pressure, momentum, energy fields
 * - 2D: writes VTK structured grid with density, velocity components,
 *       pressure, energy, and velocity vectors
 * - 3D: placeholder for future implementation
 */
class VTKWriter : public StepWriter {
public:
    VTKWriter(const std::string& output_dir, bool is_analytical);
    ~VTKWriter();

    void Write(const DataLayer& layer, const Settings& settings,
               std::size_t step, double time) const override;

    void Write(const DataLayer& layer, const DataLayer* analytical_layer,
               const Settings& settings, std::size_t step, double time) const override;

    [[nodiscard]] auto RequiresFinalization() const -> bool override;
    auto Finalize(const Settings& settings) -> std::string override;

private:
    std::string output_dir_;
    bool is_analytical_;

    // PIMPL idiom to hide VTK library dependencies from header
    class Impl;
    std::unique_ptr<Impl> pimpl_;

    /** @brief Generate output filename based on settings and step */
    auto GenerateFilename(int N, std::size_t step, const Settings& settings) const
        -> std::string;

    /** @brief Internal 1D VTK writing using VTK library objects */
    void Write1D(const DataLayer& layer, const Settings& settings,
                 std::size_t step, double time) const;

    /** @brief Internal 2D VTK writing using VTK library objects */
    void Write2D(const DataLayer& layer, const Settings& settings,
                 std::size_t step, double time) const;

    /** @brief Internal 3D VTK writing (placeholder) */
    void Write3D(const DataLayer& layer, const Settings& settings,
                 std::size_t step, double time) const;
};

#endif  // VTKWRITER_HPP
