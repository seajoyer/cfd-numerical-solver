#include "output/VTKRecomposer.hpp"

#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkFieldData.h>
#include <vtkIdList.h>
#include <vtkIntArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkXMLUnstructuredGridWriter.h>

#include <algorithm>
#include <filesystem>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace fs = std::filesystem;

namespace {
    constexpr double MERGE_TOLERANCE = 1e-8;

    struct PointKey {
        double x, y, z;

        bool operator==(const PointKey& other) const {
            return std::abs(x - other.x) < MERGE_TOLERANCE &&
                std::abs(y - other.y) < MERGE_TOLERANCE &&
                std::abs(z - other.z) < MERGE_TOLERANCE;
        }
    };

    struct PointKeyHash {
        std::size_t operator()(const PointKey& p) const {
            auto hash_double = [](double val) {
                long long rounded = static_cast<long long>(val / MERGE_TOLERANCE);
                return std::hash<long long>{}(rounded);
            };
            return hash_double(p.x) ^ (hash_double(p.y) << 1) ^ (hash_double(p.z) << 2);
        }
    };

    void AppendDataArray(vtkDataArray* src, vtkDataArray* dst, vtkIdType offset) {
        if (!src || !dst) return;
        vtkIdType n_tuples = src->GetNumberOfTuples();
        int n_comps = src->GetNumberOfComponents();
        for (vtkIdType i = 0; i < n_tuples; ++i) {
            double* tuple = src->GetTuple(i);
            dst->SetTuple(offset + i, tuple);
        }
    }

    vtkSmartPointer<vtkDataArray> CreateCompatibleArray(vtkDataArray* src, vtkIdType total_tuples) {
        vtkSmartPointer<vtkDataArray> arr;
        arr.TakeReference(src->NewInstance());
        arr->SetName(src->GetName());
        arr->SetNumberOfComponents(src->GetNumberOfComponents());
        arr->SetNumberOfTuples(total_tuples);
        return arr;
    }

    bool AreFieldDataArraysConsistent(const std::vector<vtkSmartPointer<vtkUnstructuredGrid>>& grids) {
        if (grids.empty()) return true;
        vtkFieldData* ref_fd = grids[0]->GetFieldData();
        for (size_t i = 1; i < grids.size(); ++i) {
            vtkFieldData* fd = grids[i]->GetFieldData();
            if (fd->GetNumberOfArrays() != ref_fd->GetNumberOfArrays()) return false;
            for (int j = 0; j < ref_fd->GetNumberOfArrays(); ++j) {
                vtkDataArray* arr_ref = ref_fd->GetArray(j);
                vtkDataArray* arr = fd->GetArray(arr_ref->GetName());
                if (!arr) return false;
                if (arr->GetNumberOfTuples() != arr_ref->GetNumberOfTuples() ||
                    arr->GetNumberOfComponents() != arr_ref->GetNumberOfComponents())
                    return false;
            }
        }
        return true;
    }
} // namespace

std::vector<std::string> VTKRecomposer::FindRankDirectories(const std::string& case_dir,
                                                            const int expected_rank_count) {
    if (expected_rank_count <= 0) {
        throw std::runtime_error("VTKRecomposer: expected_rank_count must be positive");
    }
    if (!fs::exists(case_dir) || !fs::is_directory(case_dir)) {
        throw std::runtime_error("VTKRecomposer: case directory does not exist or is not a directory: " + case_dir);
    }

    std::vector<std::string> rank_dirs;
    rank_dirs.reserve(expected_rank_count);
    for (int rank = 0; rank < expected_rank_count; ++rank) {
        std::ostringstream oss;
        oss << "rank_" << std::setw(3) << std::setfill('0') << rank;
        fs::path rank_dir = fs::path(case_dir) / oss.str();
        if (!fs::exists(rank_dir) || !fs::is_directory(rank_dir)) {
            throw std::runtime_error("VTKRecomposer: missing rank directory: " + rank_dir.string());
        }
        rank_dirs.push_back(rank_dir.string());
    }
    return rank_dirs;
}

void VTKRecomposer::RecomposeCaseDirectory(const std::string& case_dir,
                                           const int expected_rank_count) {
    const std::vector<std::string> rank_dirs = FindRankDirectories(case_dir, expected_rank_count);

    std::map<std::string, std::vector<fs::path>> files_by_step;
    for (const std::string& rank_dir_str : rank_dirs) {
        fs::path rank_dir(rank_dir_str);
        for (const auto& entry : fs::directory_iterator(rank_dir)) {
            if (!entry.is_regular_file()) continue;
            if (entry.path().extension() != ".vtu") continue;
            std::string filename = entry.path().filename().string();
            files_by_step[filename].push_back(entry.path());
        }
    }

    if (files_by_step.empty()) {
        throw std::runtime_error("VTKRecomposer: no .vtu files found in rank directories");
    }

    const fs::path recomposed_dir = fs::path(case_dir) / "recomposed";
    fs::create_directories(recomposed_dir);

    for (const auto& [step_filename, paths] : files_by_step) {
        if (static_cast<int>(paths.size()) != expected_rank_count) {
            throw std::runtime_error("VTKRecomposer: step '" + step_filename +
                "' does not have data from all expected ranks");
        }

        std::vector<vtkSmartPointer<vtkUnstructuredGrid>> grids;
        grids.reserve(paths.size());
        for (const auto& p : paths) {
            auto reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
            reader->SetFileName(p.string().c_str());
            reader->Update();
            grids.push_back(reader->GetOutput());
        }

        vtkSmartPointer<vtkPoints> merged_points = vtkSmartPointer<vtkPoints>::New();
        merged_points->SetDataType(VTK_DOUBLE);
        std::unordered_map<PointKey, vtkIdType, PointKeyHash> point_map;
        std::vector<std::vector<vtkIdType>> global_cell_nodes;
        global_cell_nodes.resize(grids.size());

        for (size_t g = 0; g < grids.size(); ++g) {
            vtkUnstructuredGrid* grid = grids[g];
            vtkPoints* points = grid->GetPoints();
            if (!points) continue;
            vtkIdType n_points = points->GetNumberOfPoints();
            vtkIdType n_cells = grid->GetNumberOfCells();

            global_cell_nodes[g].resize(n_cells * 8);

            for (vtkIdType cell_id = 0; cell_id < n_cells; ++cell_id) {
                vtkSmartPointer<vtkIdList> pt_ids = vtkSmartPointer<vtkIdList>::New();
                grid->GetCellPoints(cell_id, pt_ids);
                vtkIdType n_pts = pt_ids->GetNumberOfIds();
                for (vtkIdType j = 0; j < n_pts; ++j) {
                    global_cell_nodes[g][cell_id * 8 + j] = pt_ids->GetId(j);
                }
            }

            std::vector<vtkIdType> local_to_global(n_points, -1);
            for (vtkIdType pt_id = 0; pt_id < n_points; ++pt_id) {
                double xyz[3];
                points->GetPoint(pt_id, xyz);
                PointKey key{xyz[0], xyz[1], xyz[2]};
                auto it = point_map.find(key);
                if (it == point_map.end()) {
                    vtkIdType new_id = merged_points->GetNumberOfPoints();
                    merged_points->InsertNextPoint(xyz);
                    point_map[key] = new_id;
                    local_to_global[pt_id] = new_id;
                }
                else {
                    local_to_global[pt_id] = it->second;
                }
            }

            for (vtkIdType cell_id = 0; cell_id < n_cells; ++cell_id) {
                vtkSmartPointer<vtkIdList> pt_ids = vtkSmartPointer<vtkIdList>::New();
                grid->GetCellPoints(cell_id, pt_ids);
                vtkIdType n_pts = pt_ids->GetNumberOfIds();
                for (vtkIdType j = 0; j < n_pts; ++j) {
                    vtkIdType local_id = global_cell_nodes[g][cell_id * 8 + j];
                    global_cell_nodes[g][cell_id * 8 + j] = local_to_global[local_id];
                }
            }
        }

        vtkSmartPointer<vtkUnstructuredGrid> merged_grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        merged_grid->SetPoints(merged_points);

        vtkIdType total_cells = 0;
        for (auto& grid : grids) {
            total_cells += grid->GetNumberOfCells();
        }

        vtkCellData* ref_cd = grids[0]->GetCellData();
        std::vector<vtkSmartPointer<vtkDataArray>> merged_cell_arrays;
        for (int i = 0; i < ref_cd->GetNumberOfArrays(); ++i) {
            vtkDataArray* arr = ref_cd->GetArray(i);
            auto new_arr = CreateCompatibleArray(arr, total_cells);
            merged_cell_arrays.push_back(new_arr);
        }

        vtkIdType cell_offset = 0;
        for (size_t g = 0; g < grids.size(); ++g) {
            vtkUnstructuredGrid* grid = grids[g];
            vtkIdType n_cells = grid->GetNumberOfCells();

            for (vtkIdType cell_id = 0; cell_id < n_cells; ++cell_id) {
                vtkSmartPointer<vtkIdList> pt_ids = vtkSmartPointer<vtkIdList>::New();
                grid->GetCellPoints(cell_id, pt_ids);
                vtkIdType n_pts = pt_ids->GetNumberOfIds();
                vtkSmartPointer<vtkIdList> new_ids = vtkSmartPointer<vtkIdList>::New();
                new_ids->SetNumberOfIds(n_pts);
                for (vtkIdType j = 0; j < n_pts; ++j) {
                    new_ids->SetId(j, global_cell_nodes[g][cell_id * 8 + j]);
                }
                merged_grid->InsertNextCell(grid->GetCellType(cell_id), new_ids);
            }

            vtkCellData* cd = grid->GetCellData();
            for (size_t arr_idx = 0; arr_idx < merged_cell_arrays.size(); ++arr_idx) {
                vtkDataArray* src_arr = cd->GetArray(merged_cell_arrays[arr_idx]->GetName());
                if (src_arr) {
                    AppendDataArray(src_arr, merged_cell_arrays[arr_idx], cell_offset);
                }
            }
            cell_offset += n_cells;
        }

        for (auto& arr : merged_cell_arrays) {
            merged_grid->GetCellData()->AddArray(arr);
        }

        if (AreFieldDataArraysConsistent(grids)) {
            vtkFieldData* fd0 = grids[0]->GetFieldData();
            for (int i = 0; i < fd0->GetNumberOfArrays(); ++i) {
                vtkDataArray* arr = fd0->GetArray(i);
                vtkSmartPointer<vtkDataArray> copy;
                copy.TakeReference(arr->NewInstance());
                copy->DeepCopy(arr);
                merged_grid->GetFieldData()->AddArray(copy);
            }
        }

        fs::path input_name(step_filename);
        std::string stem = input_name.stem().string();
        fs::path out_path = recomposed_dir / (stem + ".vtu");

        auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        writer->SetFileName(out_path.string().c_str());
        writer->SetInputData(merged_grid);
        writer->SetDataModeToBinary();
        if (!writer->Write()) {
            throw std::runtime_error("VTKRecomposer: failed to write merged VTU file: " + out_path.string());
        }
    }
}
