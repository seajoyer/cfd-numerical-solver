#ifndef VTKWRITER_HPP
#define VTKWRITER_HPP

#include <string>

#include "output/StepWriter.hpp"

class DataLayer;
class Settings;
class Mesh;

/**
 * @class VTKWriter
 * @brief VTK unstructured-grid output for cell-centered simulations on generic meshes.
 */
class VTKWriter final : public StepWriter {
public:
    explicit VTKWriter(std::string output_dir);

    void Write(const DataLayer& layer,
               const Mesh& mesh,
               const Settings& settings,
               std::size_t step,
               double time) const override;

    [[nodiscard]] bool RequiresFinalization() const override;
    void Finalize(const Settings& settings) override;

private:
    std::string output_dir_;

    [[nodiscard]] std::string GenerateRankDirectory() const;
    [[nodiscard]] std::string GenerateFilename(std::size_t step,
                                               const Settings& settings) const;

    void EnsureDirectoriesExist() const;

    void WriteUnstructuredGrid(const DataLayer& layer,
                               const Mesh& mesh,
                               const Settings& settings,
                               std::size_t step,
                               double time) const;
};

#endif  // VTKWRITER_HPP
