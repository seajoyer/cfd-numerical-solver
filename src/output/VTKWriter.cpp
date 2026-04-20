#include "output/VTKWriter.hpp"

#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkFieldData.h>
#include <vtkIdList.h>
#include <vtkIntArray.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnsignedCharArray.h>

#include <cstddef>
#include <filesystem>
#include <iomanip>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_set>

#include "config/Settings.hpp"
#include "data/DataLayer.hpp"
#include "data/Variables.hpp"
#include "geometry/Cell.hpp"
#include "geometry/Mesh.hpp"
#include "geometry/Node.hpp"
#include "parallel/MPIContext.hpp"
#include "output/VTKRecomposer.hpp"

namespace {
    vtkSmartPointer<vtkDoubleArray> CreateDoubleArray(const char* name, vtkIdType n_tuples, int n_comp = 1) {
        auto arr = vtkSmartPointer<vtkDoubleArray>::New();
        arr->SetName(name);
        arr->SetNumberOfComponents(n_comp);
        arr->SetNumberOfTuples(n_tuples);
        return arr;
    }

    vtkSmartPointer<vtkIntArray> CreateIntArray(const char* name, vtkIdType n_tuples, int n_comp = 1) {
        auto arr = vtkSmartPointer<vtkIntArray>::New();
        arr->SetName(name);
        arr->SetNumberOfComponents(n_comp);
        arr->SetNumberOfTuples(n_tuples);
        return arr;
    }

    template <typename T>
    void AddFieldDataArray(vtkUnstructuredGrid* grid, const char* name, const T& value) {
        if constexpr (std::is_same_v<T, double> || std::is_same_v<T, float>) {
            auto arr = vtkSmartPointer<vtkDoubleArray>::New();
            arr->SetName(name);
            arr->SetNumberOfComponents(1);
            arr->InsertNextValue(static_cast<double>(value));
            grid->GetFieldData()->AddArray(arr);
        }
        else {
            auto arr = vtkSmartPointer<vtkIntArray>::New();
            arr->SetName(name);
            arr->SetNumberOfComponents(1);
            arr->InsertNextValue(static_cast<int>(value));
            grid->GetFieldData()->AddArray(arr);
        }
    }

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
    return true;
}

void VTKWriter::Finalize(const Settings& settings) {
    MPIContext mpi;

    mpi.Barrier();

    if (settings.mpi_enabled && mpi.Size() > 1 && mpi.IsRoot()) {
        VTKRecomposer::RecomposeCaseDirectory(output_dir_, mpi.Size());
    }

    mpi.Barrier();
}

std::string VTKWriter::GenerateRankDirectory() const {
    MPIContext mpi;

    std::ostringstream oss;
    oss << output_dir_
        << "/rank_"
        << std::setw(3) << std::setfill('0') << mpi.Rank();

    return oss.str();
}

std::string VTKWriter::GenerateFilename(const std::size_t step,
                                        const Settings& settings) const {
    std::ostringstream oss;
    oss << GenerateRankDirectory()
        << "/"
        << settings.solver
        << "__R_" << settings.reconstruction
        << "__step_" << std::setw(6) << std::setfill('0') << step
        << ".vtu";
    return oss.str();
}

void VTKWriter::EnsureDirectoriesExist() const {
    MPIContext mpi;

    if (mpi.IsRoot()) {
        std::filesystem::create_directories(output_dir_);
    }

    mpi.Barrier();

    std::filesystem::create_directories(GenerateRankDirectory());

    mpi.Barrier();
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

    EnsureDirectoriesExist();
    const std::string filename = GenerateFilename(step, settings);

    const std::size_t owned_cell_count = mesh.GetOwnedCellCount();

    std::unordered_set<std::size_t> used_nodes_set;
    for (std::size_t cell_id = 0; cell_id < owned_cell_count; ++cell_id) {
        const Cell& cell = mesh.GetCell(cell_id);
        for (std::size_t node_id : cell.node_ids) {
            used_nodes_set.insert(node_id);
        }
    }

    std::vector used_nodes(used_nodes_set.begin(), used_nodes_set.end());

    std::unordered_map<std::size_t, vtkIdType> old_to_new_node_id;
    for (std::size_t i = 0; i < used_nodes.size(); ++i) {
        old_to_new_node_id[used_nodes[i]] = static_cast<vtkIdType>(i);
    }

    vtkSmartPointer<vtkUnstructuredGrid> grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    const vtkIdType n_nodes = static_cast<vtkIdType>(used_nodes.size());
    points->SetNumberOfPoints(n_nodes);

    for (std::size_t i = 0; i < used_nodes.size(); ++i) {
        const Node& node = mesh.GetNode(used_nodes[i]);
        points->SetPoint(static_cast<vtkIdType>(i), node.x, node.y, node.z);
    }
    grid->SetPoints(points);

    for (std::size_t cell_id = 0; cell_id < owned_cell_count; ++cell_id) {
        const Cell& cell = mesh.GetCell(cell_id);
        const int vtk_cell_type = DetectVtkCellType(mesh, cell);

        vtkSmartPointer<vtkIdList> ids = vtkSmartPointer<vtkIdList>::New();
        ids->SetNumberOfIds(static_cast<vtkIdType>(cell.node_ids.size()));

        for (std::size_t local_node = 0; local_node < cell.node_ids.size(); ++local_node) {
            std::size_t old_node_id = cell.node_ids[local_node];
            vtkIdType new_id = old_to_new_node_id.at(old_node_id);
            ids->SetId(static_cast<vtkIdType>(local_node), new_id);
        }

        grid->InsertNextCell(vtk_cell_type, ids);
    }

    const vtkIdType n_cells = static_cast<vtkIdType>(owned_cell_count);
    const auto& U = layer.U();
    const double gamma = settings.gamma;

    vtkSmartPointer<vtkDoubleArray> arr_rho = CreateDoubleArray("density", n_cells);
    vtkSmartPointer<vtkDoubleArray> arr_vel = CreateDoubleArray("velocity", n_cells, 3);
    vtkSmartPointer<vtkDoubleArray> arr_p = CreateDoubleArray("pressure", n_cells);
    vtkSmartPointer<vtkDoubleArray> arr_e = CreateDoubleArray("conserved_energy", n_cells);
    vtkSmartPointer<vtkDoubleArray> arr_eint = CreateDoubleArray("internal_energy", n_cells);
    vtkSmartPointer<vtkDoubleArray> arr_lambda = CreateDoubleArray("reactant_mass_fraction", n_cells);

    vtkSmartPointer<vtkIntArray> arr_cell_local_id = CreateIntArray("cell_local_id", n_cells);
    vtkSmartPointer<vtkIntArray> arr_cell_original_id = CreateIntArray("cell_original_id", n_cells);

    for (std::size_t cell_id = 0; cell_id < owned_cell_count; ++cell_id) {
        const PrimitiveCell primitive = ConservativeRowToPrimitive(U, cell_id, gamma);

        const double rho = primitive.rho;
        const double u = primitive.u;
        const double v = primitive.v;
        const double w = primitive.w;
        const double P = primitive.P;
        const double E = U(cell_id, DataLayer::k_E);
        const double kinetic = 0.5 * rho * (u * u + v * v + w * w);
        const double eint = rho > 0.0 ? (E - kinetic) / rho : 0.0;

        arr_rho->SetValue(cell_id, rho);
        double velocity[3] = {u, v, w};
        arr_vel->SetTuple(cell_id, velocity);
        arr_p->SetValue(cell_id, P);
        arr_e->SetValue(cell_id, E);
        arr_eint->SetValue(cell_id, eint);
        arr_lambda->SetValue(cell_id, layer.ReactantMassFraction()(cell_id));

        const Cell& cell = mesh.GetCell(cell_id);
        arr_cell_local_id->SetValue(cell_id, static_cast<int>(cell.local_id));
        arr_cell_original_id->SetValue(cell_id, static_cast<int>(cell.id));
    }

    grid->GetCellData()->AddArray(arr_rho);
    grid->GetCellData()->AddArray(arr_vel);
    grid->GetCellData()->AddArray(arr_p);
    grid->GetCellData()->AddArray(arr_e);
    grid->GetCellData()->AddArray(arr_eint);
    grid->GetCellData()->AddArray(arr_lambda);
    grid->GetCellData()->AddArray(arr_cell_local_id);
    grid->GetCellData()->AddArray(arr_cell_original_id);

    AddFieldDataArray<double>(grid, "TimeValue", time);
    AddFieldDataArray<int>(grid, "MeshDimension", mesh.GetDim());
    AddFieldDataArray<int>(grid, "OwnedCellCount", static_cast<int>(owned_cell_count));
    AddFieldDataArray<int>(grid, "GhostCellCount", 0);
    MPIContext mpi;
    AddFieldDataArray<int>(grid, "Rank", mpi.Rank());

    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
        vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData(grid);
    writer->SetDataModeToBinary();

    if (!writer->Write()) {
        throw std::runtime_error("VTKWriter: failed to write VTU file");
    }
}
