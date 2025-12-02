#include "output/VTKWriter.hpp"

#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkStructuredGrid.h>
#include <vtkStructuredGridWriter.h>

#include <cstddef>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "data/DataLayer.hpp"
#include "utils/StringUtils.hpp"

// PIMPL implementation
class VTKWriter::Impl {
   public:
    vtkSmartPointer<vtkStructuredGrid> structured_grid;
    vtkSmartPointer<vtkPoints> points;

    Impl() {
        structured_grid = vtkSmartPointer<vtkStructuredGrid>::New();
        points = vtkSmartPointer<vtkPoints>::New();
    }
};

VTKWriter::VTKWriter(const std::string& output_dir, bool is_analytical)
    : output_dir_(output_dir), is_analytical_(is_analytical), pimpl_(std::make_unique<Impl>()) {
    // Create output directory if it doesn't exist
    std::filesystem::create_directories(output_dir_);
}

VTKWriter::~VTKWriter() = default;

auto VTKWriter::GenerateFilename(int N, std::size_t step, const Settings& settings) const
    -> std::string {
    std::ostringstream oss;

    if (is_analytical_) {
        // Analytical solution uses simple naming
        oss << output_dir_ << "/step_" << std::setw(4) << std::setfill('0') << step
            << ".vtk";
    } else {
        // Numerical solution includes solver parameters in filename
        oss << output_dir_ << "/" << settings.solver << "__R_" << settings.reconstruction
            << "__N_" << N << "__CFL_" << utils::DoubleWithoutDot(settings.cfl) <<
            "__step_" << std::setw(4) << std::setfill('0') << step << ".vtk";
    }

    return oss.str();
}

void VTKWriter::Write(const DataLayer& layer, const Settings& settings, std::size_t step,
                      double time) const {
    const int dim = layer.GetDim();

    switch (dim) {
        case 1:
            Write1D(layer, settings, step, time);
            break;
        case 2:
            Write2D(layer, settings, step, time);
            break;
        case 3:
            Write3D(layer, settings, step, time);
            break;
        default:
            throw std::runtime_error("Unsupported dimension: " + std::to_string(dim));
    }
}

void VTKWriter::Write1D(const DataLayer& layer, const Settings& settings,
                        std::size_t step, double time) const {
    const int N = layer.GetN();
    const int start = layer.GetCoreStart();
    const int end = layer.GetCoreEndExclusive();

    // Fix: Use actual core size for nx
    const int nx = end - start;
    if (nx <= 0) {
        throw std::runtime_error("Invalid core range: start=" + std::to_string(start) +
                                 ", end=" + std::to_string(end));
    }

    // Generate filename
    std::string filename = GenerateFilename(N, step, settings);

    // Create fresh VTK objects for this write
    vtkSmartPointer<vtkStructuredGrid> grid = vtkSmartPointer<vtkStructuredGrid>::New();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    // Expansion along Y for better visualization: represent 1D as thin slab in XY plane
    const int ny = 20;  // Number of points along Y (hardcoded for visualization)
    const int nz = 1;
    grid->SetDimensions(nx, ny, nz);

    // Compute domain length in X for scaling Y extent
    double min_x = layer.xc(start);
    double max_x = layer.xc(start + nx - 1);
    double Lx = max_x - min_x;
    double Ly = Lx;  // Match Y extent to X for consistent scaling
    double dy = (ny > 1) ? Ly / (ny - 1) : 0.0;

    // Total number of points
    vtkIdType num_points = static_cast<vtkIdType>(nx) * ny * nz;
    points->SetNumberOfPoints(num_points);

    // Set points: structured grid ordering (i fastest, then j, then k)
    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            double y = 0.0 + static_cast<double>(j) * dy;
            for (int i = 0; i < nx; ++i) {
                double x = layer.xc(start + i);
                double z = 0.0;
                auto pid = static_cast<vtkIdType>(i + j * nx) +
                           static_cast<vtkIdType>(k * nx * ny);
                points->SetPoint(pid, x, y, z);
            }
        }
    }
    grid->SetPoints(points);

    // Helper lambda to add scalar field (repeated along Y)
    auto add_scalar_field = [&](const xt::xarray<double>& data,
                                const char* name) -> void {
        vtkSmartPointer<vtkDoubleArray> array = vtkSmartPointer<vtkDoubleArray>::New();
        array->SetName(name);
        array->SetNumberOfComponents(1);
        array->SetNumberOfTuples(num_points);
        for (int k = 0; k < nz; ++k) {
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    const int idx = start + i;
                    double val = data(idx);
                    auto pid = static_cast<vtkIdType>(i + j * nx) +
                               static_cast<vtkIdType>(k * nx * ny);
                    array->SetValue(pid, val);
                }
            }
        }
        grid->GetPointData()->AddArray(array);
    };

    // Add all scalar fields
    add_scalar_field(layer.rho, "density");
    add_scalar_field(layer.u, "velocity");
    add_scalar_field(layer.P, "pressure");
    add_scalar_field(layer.p, "momentum");
    add_scalar_field(layer.e, "specific_internal_energy");
    add_scalar_field(layer.U, "conserved_energy");
    add_scalar_field(layer.V, "volume");
    add_scalar_field(layer.m, "mass");

    // Add time as field data
    vtkSmartPointer<vtkDoubleArray> time_array = vtkSmartPointer<vtkDoubleArray>::New();
    time_array->SetName("TimeValue");
    time_array->SetNumberOfComponents(1);
    time_array->InsertNextValue(time);
    grid->GetFieldData()->AddArray(time_array);

    // Write to file
    vtkSmartPointer<vtkStructuredGridWriter> writer =
        vtkSmartPointer<vtkStructuredGridWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData(grid);
    writer->SetFileTypeToBinary();
    writer->Write();
}

void VTKWriter::Write2D(const DataLayer& layer, const Settings& settings,
                        std::size_t step, double time) const {
    // Placeholder for 2D implementation
    throw std::runtime_error("2D VTK output not yet implemented");
}

void VTKWriter::Write3D(const DataLayer& layer, const Settings& settings,
                        std::size_t step, double time) const {
    // Placeholder for 3D implementation
    throw std::runtime_error("3D VTK output not yet implemented");
}
