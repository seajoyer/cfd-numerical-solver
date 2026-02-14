#ifndef VTKWRITER_HPP
#define VTKWRITER_HPP

#include <string>
#include "output/StepWriter.hpp"
#include "data/DataLayer.hpp"
#include "config/Settings.hpp"

/**
 * @class VTKWriter
 * @brief VTK structured grid output for 1D and 2D simulations.
 */
class VTKWriter : public StepWriter {
public:
    VTKWriter(const std::string& output_dir, bool is_analytical);

    void Write(const DataLayer& layer, const Settings& settings,
               std::size_t step, double time) const override;

    void Write(const DataLayer& layer, const DataLayer* analytical_layer,
               const Settings& settings, std::size_t step, double time) const override;

    void Write2D(const DataLayer& layer, const Settings& settings,
                 std::size_t step, double time) const override;

    [[nodiscard]] auto RequiresFinalization() const -> bool override;
    auto Finalize(const Settings& settings) -> std::string override;

private:
    std::string output_dir_;
    bool is_analytical_;
};

#endif  // VTKWRITER_HPP
