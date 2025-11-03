#include "VtkWriter.hpp"
#include <vtkStructuredGrid.h>
#include <vtkPoints.h>
#include <vtkDoubleArray.h>
#include <vtkStructuredGridWriter.h>
#include <vtkPointData.h>

// PIMPL implementation
class VtkWriter::Impl {
public:
    vtkSmartPointer<vtkStructuredGrid> structuredGrid;
    vtkSmartPointer<vtkPoints> points;
    
    Impl() {
        structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();
        points = vtkSmartPointer<vtkPoints>::New();
    }
};

VtkWriter::VtkWriter() : nx_(0), ny_(0), nz_(0), pimpl_(std::make_unique<Impl>()) {}

VtkWriter::~VtkWriter() = default;

void VtkWriter::SetDimensions(int nx, int ny, int nz) {
    if (nx <= 0 || ny <= 0 || nz <= 0) {
        throw std::invalid_argument("Dimensions must be positive");
    }
    nx_ = nx;
    ny_ = ny;
    nz_ = nz;
}

void VtkWriter::SetCoordinates(const std::vector<double>& x_coords, 
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

void VtkWriter::AddScalarField(const std::string& field_name, 
                               const std::vector<std::vector<std::vector<double>>>& data) {
    if (!ValidateDimensions()) {
        throw std::runtime_error("Dimensions not set properly");
    }
    
    // Validate data dimensions
    if (data.size() != static_cast<size_t>(nx_) || 
        data[0].size() != static_cast<size_t>(ny_) || 
        data[0][0].size() != static_cast<size_t>(nz_)) {
        throw std::invalid_argument("Scalar field dimensions don't match mesh dimensions");
    }
    
    // Create VTK array
    vtkSmartPointer<vtkDoubleArray> scalarArray = vtkSmartPointer<vtkDoubleArray>::New();
    scalarArray->SetName(field_name.c_str());
    scalarArray->SetNumberOfComponents(1);
    scalarArray->SetNumberOfValues(nx_ * ny_ * nz_);
    
    // Fill the array (VTK uses Fortran ordering: x fastest, then y, then z)
    int index = 0;
    for (int k = 0; k < nz_; ++k) {
        for (int j = 0; j < ny_; ++j) {
            for (int i = 0; i < nx_; ++i) {
                scalarArray->SetValue(index++, data[i][j][k]);
            }
        }
    }
    
    pimpl_->structuredGrid->GetPointData()->AddArray(scalarArray);
}

void VtkWriter::AddVectorField(const std::string& field_name,
                               const std::vector<std::vector<std::vector<double>>>& data_x,
                               const std::vector<std::vector<std::vector<double>>>& data_y,
                               const std::vector<std::vector<std::vector<double>>>& data_z) {
    if (!ValidateDimensions()) {
        throw std::runtime_error("Dimensions not set properly");
    }
    
    // Validate data dimensions
    if (data_x.size() != static_cast<size_t>(nx_) || data_x[0].size() != static_cast<size_t>(ny_) || data_x[0][0].size() != static_cast<size_t>(nz_) ||
        data_y.size() != static_cast<size_t>(nx_) || data_y[0].size() != static_cast<size_t>(ny_) || data_y[0][0].size() != static_cast<size_t>(nz_) ||
        data_z.size() != static_cast<size_t>(nx_) || data_z[0].size() != static_cast<size_t>(ny_) || data_z[0][0].size() != static_cast<size_t>(nz_)) {
        throw std::invalid_argument("Vector field dimensions don't match mesh dimensions");
    }
    
    // Create VTK array
    vtkSmartPointer<vtkDoubleArray> vectorArray = vtkSmartPointer<vtkDoubleArray>::New();
    vectorArray->SetName(field_name.c_str());
    vectorArray->SetNumberOfComponents(3);
    vectorArray->SetNumberOfTuples(nx_ * ny_ * nz_);
    
    // Fill the array
    int index = 0;
    for (int k = 0; k < nz_; ++k) {
        for (int j = 0; j < ny_; ++j) {
            for (int i = 0; i < nx_; ++i) {
                double vector[3] = {data_x[i][j][k], data_y[i][j][k], data_z[i][j][k]};
                vectorArray->SetTuple(index++, vector);
            }
        }
    }
    
    pimpl_->structuredGrid->GetPointData()->AddArray(vectorArray);
}

bool VtkWriter::WriteToFile(const std::string& filename) {
    if (!ValidateDimensions()) {
        std::cerr << "Error: Dimensions not set properly" << std::endl;
        return false;
    }
    
    try {
        // Create the structured grid if not already created
        if (pimpl_->points->GetNumberOfPoints() == 0) {
            CreateStructuredGrid();
        }
        
        // Create writer
        vtkSmartPointer<vtkStructuredGridWriter> writer = vtkSmartPointer<vtkStructuredGridWriter>::New();
        writer->SetFileName(filename.c_str());
        writer->SetInputData(pimpl_->structuredGrid);
        
        // Write binary format for better performance
        writer->SetFileTypeToBinary();
        
        // Write the file
        writer->Write();
        
        std::cout << "Successfully wrote VTK file: " << filename << std::endl;
        return true;
        
    } catch (const std::exception& e) {
        std::cerr << "Error writing VTK file: " << e.what() << std::endl;
        return false;
    }
}

void VtkWriter::CreateStructuredGrid() {
    // Set grid dimensions
    pimpl_->structuredGrid->SetDimensions(nx_, ny_, nz_);
    
    // Create points
    pimpl_->points->SetNumberOfPoints(nx_ * ny_ * nz_);
    
    int pointId = 0;
    for (int k = 0; k < nz_; ++k) {
        for (int j = 0; j < ny_; ++j) {
            for (int i = 0; i < nx_; ++i) {
                double x = (x_coords_.empty()) ? static_cast<double>(i) : x_coords_[i];
                double y = (y_coords_.empty()) ? static_cast<double>(j) : y_coords_[j];
                double z = (z_coords_.empty()) ? static_cast<double>(k) : z_coords_[k];
                pimpl_->points->SetPoint(pointId++, x, y, z);
            }
        }
    }
    
    pimpl_->structuredGrid->SetPoints(pimpl_->points);
}

bool VtkWriter::ValidateDimensions() const {
    return (nx_ > 0 && ny_ > 0 && nz_ > 0);
}

void VtkWriter::Clear() {
    nx_ = ny_ = nz_ = 0;
    x_coords_.clear();
    y_coords_.clear();
    z_coords_.clear();
    pimpl_->structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();
    pimpl_->points = vtkSmartPointer<vtkPoints>::New();
}
