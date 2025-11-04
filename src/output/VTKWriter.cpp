#include "output/VTKWriter.hpp"
#include "data/DataLayer.hpp"

#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkStructuredGrid.h>
#include <vtkStructuredGridWriter.h>

#include <cstddef>
#include <filesystem>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <stdexcept>

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

VTKWriter::VTKWriter(const std::string& output_dir)
    : output_dir_(output_dir), pimpl_(std::make_unique<Impl>()) {
    // Create output directory if it doesn't exist
    std::filesystem::create_directories(output_dir_);
}

VTKWriter::~VTKWriter() = default;

auto VTKWriter::GenerateFilename(int N, std::size_t step) const -> std::string {
    std::ostringstream oss;
    oss << output_dir_ << "/N_" << N << "__step_" 
        << std::setw(4) << std::setfill('0') << step << ".vtk";
    return oss.str();
}

void VTKWriter::Write(const DataLayer& layer, std::size_t step, double time) const {
    const int dim = layer.GetDim();
    
    switch (dim) {
        case 1:
            Write1D(layer, step, time);
            break;
        case 2:
            Write2D(layer, step, time);
            break;
        case 3:
            Write3D(layer, step, time);
            break;
        default:
            throw std::runtime_error("Unsupported dimension: " + std::to_string(dim));
    }
}

void VTKWriter::Write1D(const DataLayer& layer, std::size_t step, double time) const {
    const int N = layer.GetN();
    const int start = layer.GetCoreStart();
    const int end = layer.GetCoreEndExclusive();
    
    // Generate filename
    std::string filename = GenerateFilename(N, step);

    // Create fresh VTK objects for this write
    vtkSmartPointer<vtkStructuredGrid> grid = vtkSmartPointer<vtkStructuredGrid>::New();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    // Set dimensions: 1D data represented as line (nx points, 1x1 in y,z)
    const int nx = N;
    const int ny = 1;
    const int nz = 1;
    grid->SetDimensions(nx, ny, nz);

    // Create points from cell centers
    points->SetNumberOfPoints(static_cast<vtkIdType>(nx * ny) * nz);
    for (int i = 0; i < nx; ++i) {
        const int idx = start + i;
        const double x = layer.xc(idx);
        points->SetPoint(i, x, 0.0, 0.0);
    }
    grid->SetPoints(points);

    // Helper lambda to add scalar field
    auto add_scalar_field = [&](const xt::xarray<double>& data, const char* name) -> void {
        vtkSmartPointer<vtkDoubleArray> array = vtkSmartPointer<vtkDoubleArray>::New();
        array->SetName(name);
        array->SetNumberOfComponents(1);
        array->SetNumberOfValues(nx);
        for (int i = 0; i < nx; ++i) {
            const int idx = start + i;
            array->SetValue(i, data(idx));
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
    // add_scalar_field(layer.xb, "cell_boundary_coords");
    // add_scalar_field(layer.xc, "cell_center_coords");

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

    std::cout << "Wrote VTK file: " << filename 
              << " (step=" << step << ", time=" << time << ")" << '\n';
}

void VTKWriter::Write2D(const DataLayer& layer, std::size_t step, double time) const {
    // Placeholder for 2D implementation
    // In the future, this will handle 2D grids with proper indexing
    throw std::runtime_error("2D VTK output not yet implemented");
}

void VTKWriter::Write3D(const DataLayer& layer, std::size_t step, double time) const {
    // Placeholder for 3D implementation
    // In the future, this will handle 3D grids with proper indexing
    throw std::runtime_error("3D VTK output not yet implemented");
}
