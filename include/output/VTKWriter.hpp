#ifndef VTKWRITER_HPP
#define VTKWRITER_HPP

#include <memory>
#include <string>

#include "config/Settings.hpp"
#include "output/StepWriter.hpp"

class DataLayer;
class Mesh;

/**
 * @class VTKWriter
 * @brief VTK structured grid output for structured Cartesian simulations.
 */
class VTKWriter : public StepWriter {
public:
    VTKWriter(const std::string& output_dir, bool is_analytical, int rank = 0, int size = 1);
    ~VTKWriter();

    void Write(const DataLayer& layer,
               const Mesh& mesh,
               const Settings& settings,
               std::size_t step,
               double time) const override;

    void Write(const DataLayer& layer,
               const DataLayer* analytical_layer,
               const Mesh& mesh,
               const Mesh* analytical_mesh,
               const Settings& settings,
               std::size_t step,
               double time) const override;

    [[nodiscard]] auto RequiresFinalization() const -> bool override;
    auto Finalize(const Settings& settings) -> std::string override;

private:
    std::string output_dir_;
    std::string rank_output_dir_;
    bool is_analytical_;
    int rank_ = 0;
    int size_ = 1;

    class Impl;
    std::unique_ptr<Impl> pimpl_;

    auto GenerateFilename(int N, std::size_t step, const Settings& settings) const
        -> std::string;

    void Write3D(const DataLayer& layer,
                 const Mesh& mesh,
                 const Settings& settings,
                 std::size_t step,
                 double time) const;
};

#endif  // VTKWRITER_HPP
