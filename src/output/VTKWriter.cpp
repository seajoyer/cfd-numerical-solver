#include "output/VTKWriter.hpp"

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkFieldData.h>
#include <vtkIdList.h>
#include <vtkIntArray.h>
#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <cstddef>
#include <filesystem>
#include <iomanip>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

#include "config/Settings.hpp"
#include "data/DataLayer.hpp"
#include "data/Variables.hpp"
#include "geometry/Cell.hpp"
#include "geometry/Mesh.hpp"
#include "geometry/Node.hpp"
#include "utils/StringUtils.hpp"

namespace {
    [[nodiscard]] int DetectVtkCellType(const Mesh& mesh, const Cell& cell) {
        const std::size_t node_count = cell.node_ids.size();
        const int dim = mesh.GetDim();

        if (dim == 2) {
            if (node_count == 3) {
                return VTK_TRIANGLE;
            }
            if (node_count == 4) {
                return VTK_QUAD;
            }
            if (node_count > 4) {
                return VTK_POLYGON;
            }

            throw std::runtime_error("VTKWriter: unsupported 2D cell node count");
        }

        if (dim == 3) {
            if (node_count == 4) {
                return VTK_TETRA;
            }
            if (node_count == 5) {
                return VTK_PYRAMID;
            }
            if (node_count == 6) {
                return VTK_WEDGE;
            }
            if (node_count == 8) {
                return VTK_HEXAHEDRON;
            }

            throw std::runtime_error("VTKWriter: unsupported 3D cell node count");
        }

        throw std::runtime_error("VTKWriter: only 2D and 3D meshes are supported");
    }

    [[nodiscard]] PrimitiveCell ConservativeRowToPrimitive(const xt::xtensor<double, 2>& U,
                                                           const std::size_t cell_id,
                                                           const double gamma) {
        ConservativeCell conservative;
        conservative.rho = U(cell_id, DataLayer::k_rho);
        conservative.rhoU = U(cell_id, DataLayer::k_rhoU);
        conservative.rhoV = U(cell_id, DataLayer::k_rhoV);
        conservative.rhoW = U(cell_id, DataLayer::k_rhoW);
        conservative.E = U(cell_id, DataLayer::k_E);

        return PrimitiveFromConservativeCell(conservative, gamma);
    }
} // namespace

VTKWriter::VTKWriter(std::string output_dir)
    : output_dir_(std::move(output_dir)) {}

void VTKWriter::Write(const DataLayer& layer,
                      const Mesh& mesh,
                      const Settings& settings,
                      const std::size_t step,
                      const double time) const {
    WriteUnstructuredGrid(layer, mesh, settings, step, time);
}

bool VTKWriter::RequiresFinalization() const {
    return false;
}

void VTKWriter::Finalize(const Settings& settings) {
    (void)settings;
}

std::string VTKWriter::GenerateFilename(const std::size_t step,
                                        const Settings& settings) const {
    std::ostringstream oss;
    oss << output_dir_
        << "/"
        << settings.solver
        << "__R_" << settings.reconstruction
        << "__step_" << std::setw(4) << std::setfill('0') << step
        << ".vtu";
    return oss.str();
}

void VTKWriter::WriteUnstructuredGrid(const DataLayer& layer,
                                      const Mesh& mesh,
                                      const Settings& settings,
                                      const std::size_t step,
                                      const double time) const {
    if (mesh.GetCellCount() == 0) {
        throw std::runtime_error("VTKWriter: mesh has zero cells");
    }
    if (mesh.GetNodeCount() == 0) {
        throw std::runtime_error("VTKWriter: mesh has zero nodes");
    }
    if (!layer.IsAllocated()) {
        throw std::runtime_error("VTKWriter: DataLayer is not allocated");
    }
    if (layer.GetCellCount() != mesh.GetCellCount()) {
        throw std::runtime_error("VTKWriter: DataLayer cell count does not match mesh cell count");
    }

    std::filesystem::create_directories(output_dir_);
    const std::string filename = GenerateFilename(step, settings);

    vtkSmartPointer<vtkUnstructuredGrid> grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    points->SetNumberOfPoints(static_cast<vtkIdType>(mesh.GetNodeCount()));

    for (std::size_t node_id = 0; node_id < mesh.GetNodeCount(); ++node_id) {
        const Node& node = mesh.GetNode(node_id);
        points->SetPoint(static_cast<vtkIdType>(node_id), node.x, node.y, node.z);
    }

    grid->SetPoints(points);

    for (std::size_t cell_id = 0; cell_id < mesh.GetCellCount(); ++cell_id) {
        const Cell& cell = mesh.GetCell(cell_id);
        const int vtk_cell_type = DetectVtkCellType(mesh, cell);

        vtkSmartPointer<vtkIdList> ids = vtkSmartPointer<vtkIdList>::New();
        ids->SetNumberOfIds(static_cast<vtkIdType>(cell.node_ids.size()));

        for (std::size_t local_node = 0; local_node < cell.node_ids.size(); ++local_node) {
            ids->SetId(static_cast<vtkIdType>(local_node),
                       static_cast<vtkIdType>(cell.node_ids[local_node]));
        }

        grid->InsertNextCell(vtk_cell_type, ids);
    }

    const vtkIdType n_cells = static_cast<vtkIdType>(mesh.GetCellCount());
    const auto& U = layer.U();
    const double gamma = settings.gamma;

    vtkSmartPointer<vtkDoubleArray> arr_rho = vtkSmartPointer<vtkDoubleArray>::New();
    arr_rho->SetName("density");
    arr_rho->SetNumberOfComponents(1);
    arr_rho->SetNumberOfTuples(n_cells);

    vtkSmartPointer<vtkDoubleArray> arr_vel = vtkSmartPointer<vtkDoubleArray>::New();
    arr_vel->SetName("velocity");
    arr_vel->SetNumberOfComponents(3);
    arr_vel->SetNumberOfTuples(n_cells);

    vtkSmartPointer<vtkDoubleArray> arr_p = vtkSmartPointer<vtkDoubleArray>::New();
    arr_p->SetName("pressure");
    arr_p->SetNumberOfComponents(1);
    arr_p->SetNumberOfTuples(n_cells);

    vtkSmartPointer<vtkDoubleArray> arr_e = vtkSmartPointer<vtkDoubleArray>::New();
    arr_e->SetName("conserved_energy");
    arr_e->SetNumberOfComponents(1);
    arr_e->SetNumberOfTuples(n_cells);

    vtkSmartPointer<vtkDoubleArray> arr_eint = vtkSmartPointer<vtkDoubleArray>::New();
    arr_eint->SetName("internal_energy");
    arr_eint->SetNumberOfComponents(1);
    arr_eint->SetNumberOfTuples(n_cells);

    vtkSmartPointer<vtkDoubleArray> arr_lambda = vtkSmartPointer<vtkDoubleArray>::New();
    arr_lambda->SetName("reactant_mass_fraction");
    arr_lambda->SetNumberOfComponents(1);
    arr_lambda->SetNumberOfTuples(n_cells);

    vtkSmartPointer<vtkIntArray> arr_cell_id = vtkSmartPointer<vtkIntArray>::New();
    arr_cell_id->SetName("cell_id");
    arr_cell_id->SetNumberOfComponents(1);
    arr_cell_id->SetNumberOfTuples(n_cells);

    for (std::size_t cell_id = 0; cell_id < mesh.GetCellCount(); ++cell_id) {
        const PrimitiveCell primitive = ConservativeRowToPrimitive(U, cell_id, gamma);

        const double rho = primitive.rho;
        const double u = primitive.u;
        const double v = primitive.v;
        const double w = primitive.w;
        const double P = primitive.P;
        const double E = U(cell_id, DataLayer::k_E);
        const double kinetic = 0.5 * rho * (u * u + v * v + w * w);
        const double eint = rho > 0.0 ? (E - kinetic) / rho : 0.0;

        arr_rho->SetValue(static_cast<vtkIdType>(cell_id), rho);

        double velocity[3] = {u, v, w};
        arr_vel->SetTuple(static_cast<vtkIdType>(cell_id), velocity);

        arr_p->SetValue(static_cast<vtkIdType>(cell_id), P);
        arr_e->SetValue(static_cast<vtkIdType>(cell_id), E);
        arr_eint->SetValue(static_cast<vtkIdType>(cell_id), eint);
        arr_lambda->SetValue(static_cast<vtkIdType>(cell_id),
                             layer.ReactantMassFraction()(cell_id));
        arr_cell_id->SetValue(static_cast<vtkIdType>(cell_id),
                              static_cast<int>(cell_id));
    }

    grid->GetCellData()->AddArray(arr_rho);
    grid->GetCellData()->AddArray(arr_vel);
    grid->GetCellData()->AddArray(arr_p);
    grid->GetCellData()->AddArray(arr_e);
    grid->GetCellData()->AddArray(arr_eint);
    grid->GetCellData()->AddArray(arr_lambda);
    grid->GetCellData()->AddArray(arr_cell_id);

    vtkSmartPointer<vtkDoubleArray> time_array = vtkSmartPointer<vtkDoubleArray>::New();
    time_array->SetName("TimeValue");
    time_array->SetNumberOfComponents(1);
    time_array->InsertNextValue(time);
    grid->GetFieldData()->AddArray(time_array);

    vtkSmartPointer<vtkIntArray> dim_array = vtkSmartPointer<vtkIntArray>::New();
    dim_array->SetName("MeshDimension");
    dim_array->SetNumberOfComponents(1);
    dim_array->InsertNextValue(mesh.GetDim());
    grid->GetFieldData()->AddArray(dim_array);

    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
        vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData(grid);
    writer->SetDataModeToBinary();

    if (!writer->Write()) {
        throw std::runtime_error("VTKWriter: failed to write VTU file");
    }
}
