#include "visualization/VTKWriter.hpp"

#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkStructuredGrid.h>
#include <vtkStructuredGridWriter.h>

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

VTKWriter::VTKWriter()
    : nx_(0), ny_(0), nz_(0), pimpl_(std::make_unique<Impl>()) {}

VTKWriter::~VTKWriter() = default;

void VTKWriter::SetDimensions(int nx, int ny, int nz) {
    if (nx <= 0 || ny <= 0 || nz <= 0) {
        throw std::invalid_argument("Dimensions must be positive");
    }
    nx_ = nx;
    ny_ = ny;
    nz_ = nz;
}

void VTKWriter::SetCoordinates(const std::vector<double>& x_coords,
                               const std::vector<double>& y_coords,
                               const std::vector<double>& z_coords) {
    if (x_coords.size() != static_cast<size_t>(nx_) ||
        y_coords.size() != static_cast<size_t>(ny_) ||
        z_coords.size() != static_cast<size_t>(nz_)) {
        throw std::invalid_argument("Coordinate sizes don't match dimensions");
    }

    x_coords_ = x_coords;
    y_coords_ = y_coords;
    z_coords_ = z_coords;
}

void VTKWriter::AddScalarField(
    const std::string& field_name,
    const std::vector<std::vector<std::vector<double>>>& data) {
    if (!ValidateDimensions()) {
        throw std::runtime_error("Dimensions not set properly");
    }

    // Validate data dimensions
    if (data.size() != static_cast<size_t>(nx_) ||
        data[0].size() != static_cast<size_t>(ny_) ||
        data[0][0].size() != static_cast<size_t>(nz_)) {
        throw std::invalid_argument(
            "Scalar field dimensions don't match mesh dimensions");
    }

    // Create VTK array
    vtkSmartPointer<vtkDoubleArray> scalar_array =
        vtkSmartPointer<vtkDoubleArray>::New();
    scalar_array->SetName(field_name.c_str());
    scalar_array->SetNumberOfComponents(1);
    scalar_array->SetNumberOfValues(nx_ * ny_ * nz_);

    // Fill the array (VTK uses Fortran ordering: x fastest, then y, then z)
    int index = 0;
    for (int k = 0; k < nz_; ++k) {
        for (int j = 0; j < ny_; ++j) {
            for (int i = 0; i < nx_; ++i) {
                scalar_array->SetValue(index++, data[i][j][k]);
            }
        }
    }

    pimpl_->structured_grid->GetPointData()->AddArray(scalar_array);
}

void VTKWriter::AddVectorField(
    const std::string& field_name,
    const std::vector<std::vector<std::vector<double>>>& data_x,
    const std::vector<std::vector<std::vector<double>>>& data_y,
    const std::vector<std::vector<std::vector<double>>>& data_z) {
    if (!ValidateDimensions()) {
        throw std::runtime_error("Dimensions not set properly");
    }

    // Validate data dimensions
    if (data_x.size() != static_cast<size_t>(nx_) ||
        data_x[0].size() != static_cast<size_t>(ny_) ||
        data_x[0][0].size() != static_cast<size_t>(nz_) ||
        data_y.size() != static_cast<size_t>(nx_) ||
        data_y[0].size() != static_cast<size_t>(ny_) ||
        data_y[0][0].size() != static_cast<size_t>(nz_) ||
        data_z.size() != static_cast<size_t>(nx_) ||
        data_z[0].size() != static_cast<size_t>(ny_) ||
        data_z[0][0].size() != static_cast<size_t>(nz_)) {
        throw std::invalid_argument(
            "Vector field dimensions don't match mesh dimensions");
    }

    // Create VTK array
    vtkSmartPointer<vtkDoubleArray> vector_array =
        vtkSmartPointer<vtkDoubleArray>::New();
    vector_array->SetName(field_name.c_str());
    vector_array->SetNumberOfComponents(3);
    vector_array->SetNumberOfTuples(nx_ * ny_ * nz_);

    // Fill the array
    int index = 0;
    for (int k = 0; k < nz_; ++k) {
        for (int j = 0; j < ny_; ++j) {
            for (int i = 0; i < nx_; ++i) {
                double vector[3] = {data_x[i][j][k], data_y[i][j][k],
                                    data_z[i][j][k]};
                vector_array->SetTuple(index++, vector);
            }
        }
    }

    pimpl_->structured_grid->GetPointData()->AddArray(vector_array);
}

auto VTKWriter::WriteToFile(const std::string& filename) -> bool {
    if (!ValidateDimensions()) {
        std::cerr << "Error: Dimensions not set properly" << '\n';
        return false;
    }

    try {
        // Create the structured grid if not already created
        if (pimpl_->points->GetNumberOfPoints() == 0) {
            CreateStructuredGrid();
        }

        // Create writer
        vtkSmartPointer<vtkStructuredGridWriter> writer =
            vtkSmartPointer<vtkStructuredGridWriter>::New();
        writer->SetFileName(filename.c_str());
        writer->SetInputData(pimpl_->structured_grid);

        // Write binary format for better performance
        writer->SetFileTypeToBinary();

        // Write the file
        writer->Write();

        std::cout << "Successfully wrote VTK file: " << filename << '\n';
        return true;

    } catch (const std::exception& e) {
        std::cerr << "Error writing VTK file: " << e.what() << '\n';
        return false;
    }
}

void VTKWriter::CreateStructuredGrid() {
    // Set grid dimensions
    pimpl_->structured_grid->SetDimensions(nx_, ny_, nz_);

    // Create points
    pimpl_->points->SetNumberOfPoints(nx_ * ny_ * nz_);

    int point_id = 0;
    for (int k = 0; k < nz_; ++k) {
        for (int j = 0; j < ny_; ++j) {
            for (int i = 0; i < nx_; ++i) {
                double x =
                    (x_coords_.empty()) ? static_cast<double>(i) : x_coords_[i];
                double y =
                    (y_coords_.empty()) ? static_cast<double>(j) : y_coords_[j];
                double z =
                    (z_coords_.empty()) ? static_cast<double>(k) : z_coords_[k];
                pimpl_->points->SetPoint(point_id++, x, y, z);
            }
        }
    }

    pimpl_->structured_grid->SetPoints(pimpl_->points);
}

[[nodiscard]] auto VTKWriter::ValidateDimensions() const -> bool {
    return (nx_ > 0 && ny_ > 0 && nz_ > 0);
}

void VTKWriter::Clear() {
    nx_ = ny_ = nz_ = 0;
    x_coords_.clear();
    y_coords_.clear();
    z_coords_.clear();
    pimpl_->structured_grid = vtkSmartPointer<vtkStructuredGrid>::New();
    pimpl_->points = vtkSmartPointer<vtkPoints>::New();
}
