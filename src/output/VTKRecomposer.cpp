#include "output/VTKRecomposer.hpp"

#include <filesystem>
#include <iomanip>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <vtkDataArray.h>
#include <vtkDataSetAttributes.h>
#include <vtkDoubleArray.h>
#include <vtkFieldData.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include <vtkStructuredGridReader.h>
#include <vtkStructuredGridWriter.h>

#include "utils/StringUtils.hpp"

namespace fs = std::filesystem;

namespace {
    struct PieceInfo {
        vtkSmartPointer<vtkStructuredGrid> grid;
        int offset_x = 0;
        int offset_y = 0;
        int offset_z = 0;
        int local_nx = 0;
        int local_ny = 0;
        int local_nz = 0;
    };

    auto ReadScalarFieldData(vtkStructuredGrid* grid, const char* name) -> int {
        if (!grid) {
            throw std::runtime_error("VTKRecomposer: null grid in ReadScalarFieldData");
        }

        vtkFieldData* field_data = grid->GetFieldData();
        if (!field_data) {
            throw std::runtime_error("VTKRecomposer: missing field data");
        }

        vtkDataArray* arr = field_data->GetArray(name);
        if (!arr || arr->GetNumberOfTuples() < 1) {
            throw std::runtime_error(std::string("VTKRecomposer: missing field-data array: ") + name);
        }

        return static_cast<int>(arr->GetTuple1(0));
    }

    auto ReadGrid(const std::string& filename) -> vtkSmartPointer<vtkStructuredGrid> {
        vtkNew<vtkStructuredGridReader> reader;
        reader->SetFileName(filename.c_str());
        reader->Update();

        vtkStructuredGrid* raw = reader->GetOutput();
        if (!raw) {
            throw std::runtime_error("VTKRecomposer: failed to read file: " + filename);
        }

        vtkSmartPointer<vtkStructuredGrid> grid = vtkSmartPointer<vtkStructuredGrid>::New();
        grid->DeepCopy(raw);
        return grid;
    }

    auto MakeStepFilename(const int step, const Settings& settings) -> std::string {
        std::ostringstream oss;
        oss << settings.solver
            << "__R_" << settings.reconstruction
            << "__N_" << settings.GetNx() << "x" << settings.GetNy() << "x" << settings.GetNz()
            << "__CFL_" << utils::DoubleWithoutDot(settings.cfl)
            << "__step_" << std::setw(4) << std::setfill('0') << step << ".vtk";
        return oss.str();
    }

    auto IsGhostArrayName(const std::string& name) -> bool {
        const char* ghost_name = vtkDataSetAttributes::GhostArrayName();
        if (ghost_name != nullptr && name == ghost_name) {
            return true;
        }

        // Дополнительная защита, если в какой-то сборке/формате имя будет нестандартно обрабатываться
        if (name == "vtkGhostType") {
            return true;
        }

        return false;
    }

    void CopyPointDataArray(vtkDataArray* src_array,
                            vtkDataArray* dst_array,
                            const int local_nx,
                            const int local_ny,
                            const int local_nz,
                            const int offset_x,
                            const int offset_y,
                            const int offset_z,
                            const int global_nx,
                            const int global_ny,
                            const int global_nz) {
        if (!src_array || !dst_array) {
            throw std::runtime_error("VTKRecomposer: null array in CopyPointDataArray");
        }

        const int num_components = src_array->GetNumberOfComponents();
        if (dst_array->GetNumberOfComponents() != num_components) {
            throw std::runtime_error("VTKRecomposer: component count mismatch while copying point data");
        }

        std::vector<double> tuple(static_cast<std::size_t>(num_components), 0.0);

        for (int k = 0; k < local_nz; ++k) {
            for (int j = 0; j < local_ny; ++j) {
                for (int i = 0; i < local_nx; ++i) {
                    const vtkIdType local_id =
                        static_cast<vtkIdType>(i) +
                        static_cast<vtkIdType>(j) * local_nx +
                        static_cast<vtkIdType>(k) * local_nx * local_ny;

                    const int gi = offset_x + i;
                    const int gj = offset_y + j;
                    const int gk = offset_z + k;

                    if (gi < 0 || gi >= global_nx ||
                        gj < 0 || gj >= global_ny ||
                        gk < 0 || gk >= global_nz) {
                        throw std::runtime_error("VTKRecomposer: global index out of range while copying point data");
                    }

                    const vtkIdType global_id =
                        static_cast<vtkIdType>(gi) +
                        static_cast<vtkIdType>(gj) * global_nx +
                        static_cast<vtkIdType>(gk) * global_nx * global_ny;

                    src_array->GetTuple(local_id, tuple.data());
                    dst_array->SetTuple(global_id, tuple.data());
                }
            }
        }
    }

    void CopyBlanking(vtkStructuredGrid* local_grid,
                      vtkStructuredGrid* global_grid,
                      const int local_nx,
                      const int local_ny,
                      const int local_nz,
                      const int offset_x,
                      const int offset_y,
                      const int offset_z,
                      const int global_nx,
                      const int global_ny,
                      const int global_nz) {
        if (!local_grid || !global_grid) {
            throw std::runtime_error("VTKRecomposer: null grid in CopyBlanking");
        }

        for (int k = 0; k < local_nz; ++k) {
            for (int j = 0; j < local_ny; ++j) {
                for (int i = 0; i < local_nx; ++i) {
                    const vtkIdType local_id =
                        static_cast<vtkIdType>(i) +
                        static_cast<vtkIdType>(j) * local_nx +
                        static_cast<vtkIdType>(k) * local_nx * local_ny;

                    const int gi = offset_x + i;
                    const int gj = offset_y + j;
                    const int gk = offset_z + k;

                    if (gi < 0 || gi >= global_nx ||
                        gj < 0 || gj >= global_ny ||
                        gk < 0 || gk >= global_nz) {
                        throw std::runtime_error("VTKRecomposer: global index out of range while copying blanking");
                    }

                    const vtkIdType global_id =
                        static_cast<vtkIdType>(gi) +
                        static_cast<vtkIdType>(gj) * global_nx +
                        static_cast<vtkIdType>(gk) * global_nx * global_ny;

                    if (!local_grid->IsPointVisible(local_id)) {
                        global_grid->BlankPoint(global_id);
                    }
                }
            }
        }
    }

    void CopyPoints(vtkPoints* local_points,
                    vtkPoints* global_points,
                    const int local_nx,
                    const int local_ny,
                    const int local_nz,
                    const int offset_x,
                    const int offset_y,
                    const int offset_z,
                    const int global_nx,
                    const int global_ny,
                    const int global_nz) {
        if (!local_points || !global_points) {
            throw std::runtime_error("VTKRecomposer: null points in CopyPoints");
        }

        for (int k = 0; k < local_nz; ++k) {
            for (int j = 0; j < local_ny; ++j) {
                for (int i = 0; i < local_nx; ++i) {
                    const vtkIdType local_id =
                        static_cast<vtkIdType>(i) +
                        static_cast<vtkIdType>(j) * local_nx +
                        static_cast<vtkIdType>(k) * local_nx * local_ny;

                    const int gi = offset_x + i;
                    const int gj = offset_y + j;
                    const int gk = offset_z + k;

                    if (gi < 0 || gi >= global_nx ||
                        gj < 0 || gj >= global_ny ||
                        gk < 0 || gk >= global_nz) {
                        throw std::runtime_error("VTKRecomposer: global index out of range while copying points");
                    }

                    const vtkIdType global_id =
                        static_cast<vtkIdType>(gi) +
                        static_cast<vtkIdType>(gj) * global_nx +
                        static_cast<vtkIdType>(gk) * global_nx * global_ny;

                    double xyz[3] = {0.0, 0.0, 0.0};
                    local_points->GetPoint(local_id, xyz);
                    global_points->SetPoint(global_id, xyz);
                }
            }
        }
    }

    auto BuildRankDirectories(const fs::path& base_dir, const int size) -> std::vector<fs::path> {
        std::vector<fs::path> rank_dirs;
        rank_dirs.reserve(static_cast<std::size_t>(size));

        for (int r = 0; r < size; ++r) {
            std::ostringstream oss;
            oss << "rank_" << std::setw(4) << std::setfill('0') << r;
            const fs::path rank_dir = base_dir / oss.str();

            if (!fs::exists(rank_dir)) {
                throw std::runtime_error("VTKRecomposer: missing rank directory: " + rank_dir.string());
            }
            if (!fs::is_directory(rank_dir)) {
                throw std::runtime_error("VTKRecomposer: path is not a directory: " + rank_dir.string());
            }

            rank_dirs.push_back(rank_dir);
        }

        return rank_dirs;
    }

    auto CollectAllSteps(const std::vector<fs::path>& rank_dirs) -> std::vector<int> {
        std::set<int> all_steps;

        for (const auto& rank_dir : rank_dirs) {
            for (const auto& entry : fs::directory_iterator(rank_dir)) {
                if (!entry.is_regular_file()) {
                    continue;
                }
                if (entry.path().extension() != ".vtk") {
                    continue;
                }

                const int step = VTKRecomposer::ExtractStepNumber(entry.path().filename().string());
                if (step >= 0) {
                    all_steps.insert(step);
                }
            }
        }

        return std::vector<int>(all_steps.begin(), all_steps.end());
    }

    auto ReadPiecesForStep(const std::vector<fs::path>& rank_dirs,
                           const int step,
                           const Settings& settings,
                           int& global_nx,
                           int& global_ny,
                           int& global_nz) -> std::vector<PieceInfo> {
        std::vector<PieceInfo> pieces;
        pieces.reserve(rank_dirs.size());

        const std::string step_filename = MakeStepFilename(step, settings);

        global_nx = -1;
        global_ny = -1;
        global_nz = -1;

        for (const auto& rank_dir : rank_dirs) {
            const fs::path file = rank_dir / step_filename;
            if (!fs::exists(file)) {
                throw std::runtime_error("VTKRecomposer: missing file: " + file.string());
            }

            vtkSmartPointer<vtkStructuredGrid> grid = ReadGrid(file.string());

            int dims[3] = {0, 0, 0};
            grid->GetDimensions(dims);

            PieceInfo piece;
            piece.grid = grid;
            piece.offset_x = ReadScalarFieldData(grid, "OffsetX");
            piece.offset_y = ReadScalarFieldData(grid, "OffsetY");
            piece.offset_z = ReadScalarFieldData(grid, "OffsetZ");
            piece.local_nx = dims[0];
            piece.local_ny = dims[1];
            piece.local_nz = dims[2];

            const int gx = ReadScalarFieldData(grid, "GlobalNx");
            const int gy = ReadScalarFieldData(grid, "GlobalNy");
            const int gz = ReadScalarFieldData(grid, "GlobalNz");

            if (piece.local_nx <= 0 || piece.local_ny <= 0 || piece.local_nz <= 0) {
                throw std::runtime_error("VTKRecomposer: invalid local grid dimensions in file: " + file.string());
            }

            if (global_nx < 0) {
                global_nx = gx;
                global_ny = gy;
                global_nz = gz;
            }
            else if (global_nx != gx || global_ny != gy || global_nz != gz) {
                throw std::runtime_error("VTKRecomposer: inconsistent global sizes across rank files");
            }

            pieces.push_back(std::move(piece));
        }

        if (global_nx <= 0 || global_ny <= 0 || global_nz <= 0) {
            throw std::runtime_error("VTKRecomposer: invalid global dimensions for recomposed grid");
        }

        return pieces;
    }

    auto CreateGlobalArrays(vtkPointData* ref_pd,
                            const vtkIdType global_num_points)
        -> std::unordered_map<std::string, vtkSmartPointer<vtkDataArray>> {
        if (!ref_pd) {
            throw std::runtime_error("VTKRecomposer: missing reference point data");
        }

        std::unordered_map<std::string, vtkSmartPointer<vtkDataArray>> global_arrays;

        for (int a = 0; a < ref_pd->GetNumberOfArrays(); ++a) {
            vtkDataArray* ref_array = ref_pd->GetArray(a);
            if (!ref_array) {
                continue;
            }

            const std::string name = ref_array->GetName() ? ref_array->GetName() : "";
            if (name.empty()) {
                continue;
            }

            // Критично: ghost/blanking массивы не переносим как обычные данные.
            if (IsGhostArrayName(name)) {
                continue;
            }

            vtkSmartPointer<vtkDataArray> arr;
            arr.TakeReference(ref_array->NewInstance());
            arr->SetName(name.c_str());
            arr->SetNumberOfComponents(ref_array->GetNumberOfComponents());
            arr->SetNumberOfTuples(global_num_points);

            // Явная инициализация
            arr->FillComponent(0, 0.0);
            for (int c = 1; c < ref_array->GetNumberOfComponents(); ++c) {
                arr->FillComponent(c, 0.0);
            }

            global_arrays.emplace(name, arr);
        }

        return global_arrays;
    }

    void CopyAllPointData(const std::vector<PieceInfo>& pieces,
                          std::unordered_map<std::string, vtkSmartPointer<vtkDataArray>>& global_arrays,
                          const int global_nx,
                          const int global_ny,
                          const int global_nz) {
        for (const auto& piece : pieces) {
            vtkStructuredGrid* local_grid = piece.grid;
            vtkPointData* local_pd = local_grid->GetPointData();
            if (!local_pd) {
                throw std::runtime_error("VTKRecomposer: missing point data in one of the pieces");
            }

            for (int a = 0; a < local_pd->GetNumberOfArrays(); ++a) {
                vtkDataArray* src_array = local_pd->GetArray(a);
                if (!src_array) {
                    continue;
                }

                const std::string name = src_array->GetName() ? src_array->GetName() : "";
                if (name.empty()) {
                    continue;
                }

                if (IsGhostArrayName(name)) {
                    continue;
                }

                auto it = global_arrays.find(name);
                if (it == global_arrays.end()) {
                    throw std::runtime_error("VTKRecomposer: destination array not found for: " + name);
                }

                vtkDataArray* dst_array = it->second;
                CopyPointDataArray(src_array,
                                   dst_array,
                                   piece.local_nx, piece.local_ny, piece.local_nz,
                                   piece.offset_x, piece.offset_y, piece.offset_z,
                                   global_nx, global_ny, global_nz);
            }
        }
    }

    void CopyFieldDataExceptTopologyMetadata(vtkStructuredGrid* src_grid,
                                             vtkStructuredGrid* dst_grid) {
        if (!src_grid || !dst_grid) {
            throw std::runtime_error("VTKRecomposer: null grid in CopyFieldDataExceptTopologyMetadata");
        }

        vtkFieldData* local_field = src_grid->GetFieldData();
        if (!local_field) {
            return;
        }

        for (int a = 0; a < local_field->GetNumberOfArrays(); ++a) {
            vtkDataArray* arr = local_field->GetArray(a);
            if (!arr) {
                continue;
            }

            const std::string name = arr->GetName() ? arr->GetName() : "";
            if (name == "GlobalNx" || name == "GlobalNy" || name == "GlobalNz" ||
                name == "OffsetX" || name == "OffsetY" || name == "OffsetZ") {
                continue;
            }

            vtkSmartPointer<vtkDataArray> copy;
            copy.TakeReference(arr->NewInstance());
            copy->DeepCopy(arr);
            dst_grid->GetFieldData()->AddArray(copy);
        }
    }
} // namespace

VTKRecomposer::VTKRecomposer(std::string output_dir, const int rank, const int size)
    : output_dir_(std::move(output_dir)), rank_(rank), size_(size) {}

auto VTKRecomposer::ExtractStepNumber(const std::string& filename) -> int {
    const std::size_t pos = filename.find("step_");
    if (pos == std::string::npos) {
        return -1;
    }

    const std::size_t start = pos + 5;
    const std::size_t end = filename.find(".vtk", start);
    if (end == std::string::npos) {
        return -1;
    }

    try {
        return std::stoi(filename.substr(start, end - start));
    }
    catch (...) {
        return -1;
    }
}

void VTKRecomposer::RecomposeAssigned(const Settings& settings) const {
    if (size_ <= 1) {
        return;
    }

    const fs::path base_dir(output_dir_);
    if (!fs::exists(base_dir) || !fs::is_directory(base_dir)) {
        throw std::runtime_error("VTKRecomposer: invalid base directory: " + base_dir.string());
    }

    const std::vector<fs::path> rank_dirs = BuildRankDirectories(base_dir, size_);
    const std::vector<int> steps = CollectAllSteps(rank_dirs);

    if (steps.empty()) {
        return;
    }

    const fs::path recomposed_dir = base_dir / "recomposed";
    fs::create_directories(recomposed_dir);

    for (std::size_t idx = 0; idx < steps.size(); ++idx) {
        if (static_cast<int>(idx % static_cast<std::size_t>(size_)) != rank_) {
            continue;
        }

        const int step = steps[idx];

        int global_nx = -1;
        int global_ny = -1;
        int global_nz = -1;

        const std::vector<PieceInfo> pieces =
            ReadPiecesForStep(rank_dirs, step, settings, global_nx, global_ny, global_nz);

        vtkNew<vtkStructuredGrid> global_grid;
        global_grid->SetDimensions(global_nx, global_ny, global_nz);

        const vtkIdType global_num_points =
            static_cast<vtkIdType>(global_nx) *
            static_cast<vtkIdType>(global_ny) *
            static_cast<vtkIdType>(global_nz);

        vtkNew<vtkPoints> global_points;
        global_points->SetNumberOfPoints(global_num_points);

        vtkPointData* ref_pd = pieces.front().grid->GetPointData();
        auto global_arrays = CreateGlobalArrays(ref_pd, global_num_points);

        for (const auto& piece : pieces) {
            vtkStructuredGrid* local_grid = piece.grid;
            vtkPoints* local_points = local_grid->GetPoints();
            if (!local_points) {
                throw std::runtime_error("VTKRecomposer: missing points in one of the pieces");
            }

            CopyPoints(local_points,
                       global_points,
                       piece.local_nx, piece.local_ny, piece.local_nz,
                       piece.offset_x, piece.offset_y, piece.offset_z,
                       global_nx, global_ny, global_nz);
        }

        CopyAllPointData(pieces, global_arrays, global_nx, global_ny, global_nz);

        global_grid->SetPoints(global_points);

        for (const auto& piece : pieces) {
            CopyBlanking(piece.grid,
                         global_grid,
                         piece.local_nx, piece.local_ny, piece.local_nz,
                         piece.offset_x, piece.offset_y, piece.offset_z,
                         global_nx, global_ny, global_nz);
        }

        for (const auto& [name, arr] : global_arrays) {
            (void)name;
            global_grid->GetPointData()->AddArray(arr);
        }

        CopyFieldDataExceptTopologyMetadata(pieces.front().grid, global_grid);

        vtkNew<vtkStructuredGridWriter> writer;
        const std::string out_name = (recomposed_dir / MakeStepFilename(step, settings)).string();

        writer->SetFileName(out_name.c_str());
        writer->SetInputData(global_grid);
        writer->SetFileTypeToBinary();

        if (writer->Write() == 0) {
            throw std::runtime_error("VTKRecomposer: failed to write recomposed file: " + out_name);
        }
    }
}
