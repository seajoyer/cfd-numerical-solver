#ifndef VTK_WRITER_H_
#define VTK_WRITER_H_

#include <vector>
#include <string>
#include <memory>

// Forward declarations to avoid including VTK headers in the header
class vtkStructuredGrid;
class vtkDoubleArray;
class vtkPoints;

class VTKWriter {
public:
    // Constructor
    VTKWriter();
    
    // Destructor
    ~VTKWriter();
    
    /**
     * @brief Set the mesh dimensions (nx, ny, nz)
     * @param nx Number of points in x-direction
     * @param ny Number of points in y-direction
     * @param nz Number of points in z-direction
     */
    void SetDimensions(int nx, int ny, int nz);
    
    /**
     * @brief Set the physical coordinates of the mesh points
     * @param x_coords 1D vector of x-coordinates (size should be nx)
     * @param y_coords 1D vector of y-coordinates (size should be ny)
     * @param z_coords 1D vector of z-coordinates (size should be nz)
     */
    void SetCoordinates(const std::vector<double>& x_coords, 
                       const std::vector<double>& y_coords, 
                       const std::vector<double>& z_coords);
    
    /**
     * @brief Add a scalar field to the mesh
     * @param field_name Name of the field (e.g., "Temperature", "Pressure")
     * @param data 3D matrix of scalar values (nx × ny × nz)
     */
    void AddScalarField(const std::string& field_name, 
                       const std::vector<std::vector<std::vector<double>>>& data);
    
    /**
     * @brief Add a vector field to the mesh
     * @param field_name Name of the field (e.g., "Velocity", "Force")
     * @param data_x 3D matrix of x-components (nx × ny × nz)
     * @param data_y 3D matrix of y-components (nx × ny × nz)
     * @param data_z 3D matrix of z-components (nx × ny × nz)
     */
    void AddVectorField(const std::string& field_name,
                       const std::vector<std::vector<std::vector<double>>>& data_x,
                       const std::vector<std::vector<std::vector<double>>>& data_y,
                       const std::vector<std::vector<std::vector<double>>>& data_z);
    
    /**
     * @brief Write the mesh and fields to a VTK file
     * @param filename Output filename (should end with .vtk)
     * @return true if successful, false otherwise
     */
    bool WriteToFile(const std::string& filename);
    
    /**
     * @brief Clear all data from the class
     */
    void Clear();
    
private:
    int nx_, ny_, nz_;
    std::vector<double> x_coords_, y_coords_, z_coords_;
    
    // PIMPL idiom to hide VTK implementation details
    class Impl;
    std::unique_ptr<Impl> pimpl_;
    
    // Helper methods
    bool ValidateDimensions() const;
    void CreateStructuredGrid();
};

#endif // VTK_WRITER_H_
