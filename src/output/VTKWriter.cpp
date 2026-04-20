#include "output/VTKWriter.hpp"

#include <vtkDoubleArray.h>
#include <vtkFieldData.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkStructuredGrid.h>
#include <vtkStructuredGridWriter.h>

#include <cstddef>
#include <filesystem>
#include <iomanip>
#include <memory>
#include <sstream>
#include <format>
#include <stdexcept>
#include <algorithm>

#include "data/DataLayer.hpp"
#include "data/Mesh.hpp"
#include "data/PressureVelocityState.hpp"
#include "utils/StringUtils.hpp"
#include "output/VTKRecomposer.hpp"
#include "solver/EOS.hpp"


class VTKWriter::Impl {
public:
    vtkSmartPointer<vtkStructuredGrid> structured_grid;
    vtkSmartPointer<vtkPoints> points;

    Impl() {
        structured_grid = vtkSmartPointer<vtkStructuredGrid>::New();
        points = vtkSmartPointer<vtkPoints>::New();
    }
};

static inline void PressureVelocityCellCenteredVelocity(
    const PressureVelocityState& state,
    const int i,
    const int j,
    const int k,
    double& u,
    double& v,
    double& w
) {
    const auto& ux = state.Ux();
    const auto& vy = state.Vy();
    const auto& wz = state.Wz();

    u = 0.5 * (ux(i, j, k) + ux(i + 1, j, k));
    v = 0.5 * (vy(i, j, k) + vy(i, j + 1, k));
    w = 0.5 * (wz(i, j, k) + wz(i, j, k + 1));
}

VTKWriter::VTKWriter(const std::string& output_dir, std::shared_ptr<EOS> eos, const bool is_analytical, const int rank,
                     const int size) : output_dir_(output_dir),
                                       is_analytical_(is_analytical),
                                       rank_(rank),
                                       size_(size),
                                       eos_(std::move(eos)),
                                       pimpl_(std::make_unique<Impl>()) {
    if (size_ > 1) {
        rank_output_dir_ = output_dir_ + "/" + std::format("rank_{:04d}", rank_);
    }
    else {
        rank_output_dir_ = output_dir_;
    }
}

VTKWriter::~VTKWriter() = default;

auto VTKWriter::RequiresFinalization() const -> bool {
    return size_ > 1;
}

void VTKWriter::Finalize(const Settings& settings) {
    if (size_ <= 1) {
        return;
    }

    VTKRecomposer recomposer(output_dir_, rank_, size_);
    recomposer.RecomposeAssigned(settings);
}

auto VTKWriter::GenerateFilename(const int N, const std::size_t step, const Settings& settings) const
    -> std::string {
    std::ostringstream oss;

    if (is_analytical_) {
        oss << rank_output_dir_ << "/step_" << std::setw(4) << std::setfill('0') << step << ".vtk";
    }
    else {
        oss << rank_output_dir_ << "/" << settings.solver << "__R_" << settings.reconstruction
            << "__N_" << settings.GetNx() << "x" << settings.GetNy() << "x" << settings.GetNz()
            << "__CFL_" << utils::DoubleWithoutDot(settings.cfl)
            << "__step_" << std::setw(4) << std::setfill('0') << step << ".vtk";
    }

    return oss.str();
}

void VTKWriter::Write(const DataLayer& layer,
                      const Mesh& mesh,
                      const Settings& settings,
                      const std::size_t step,
                      const double time) const {
    Write3D(layer, mesh, settings, step, time);
}

void VTKWriter::Write(const PressureVelocityState& state,
                      const Mesh& mesh,
                      const Settings& settings,
                      const std::size_t step,
                      const double time) const {
    Write3D(state, mesh, settings, step, time);
}

void VTKWriter::Write(const DataLayer& layer,
                      const DataLayer* analytical_layer,
                      const Mesh& mesh,
                      const Mesh* analytical_mesh,
                      const Settings& settings,
                      const std::size_t step,
                      const double time) const {
    (void)analytical_layer;
    (void)analytical_mesh;
    Write(layer, mesh, settings, step, time);
}

static inline void ConservativeToPrimitive(
    const xt::xtensor<double, 4>& U,
    const int i, const int j, const int k,
    const std::shared_ptr<EOS>& eos,
    const double lambda,
    double& rho, double& u, double& v, double& w, double& P, double& T) {
    rho = U(DataLayer::k_rho, i, j, k);
    if (rho <= 0.0) {
        rho = 0.0;
        u = 0.0;
        v = 0.0;
        w = 0.0;
        P = 0.0;
        T = 0.0;
        return;
    }

    const double inv_rho = 1.0 / rho;
    u = U(DataLayer::k_rhoU, i, j, k) * inv_rho;
    v = U(DataLayer::k_rhoV, i, j, k) * inv_rho;
    w = U(DataLayer::k_rhoW, i, j, k) * inv_rho;

    const double E = U(DataLayer::k_E, i, j, k);
    const double kinetic = 0.5 * rho * (u * u + v * v + w * w);

    const double I_cell = std::max((E - kinetic) * inv_rho, 0.0);

    EosCellInput eos_in{rho, I_cell, lambda};
    EosCellOutput eos_out = eos->Evaluate(eos_in);

    P = eos_out.P;
    T = eos_out.T;
}

void VTKWriter::Write3D(const DataLayer& layer,
                        const Mesh& mesh,
                        const Settings& settings,
                        const std::size_t step,
                        const double time) const {
    const int cs_x = mesh.GetCoreStartX();
    const int ce_x = mesh.GetCoreEndExclusiveX();
    const int cs_y = mesh.GetCoreStartY();
    const int ce_y = mesh.GetCoreEndExclusiveY();
    const int cs_z = mesh.GetCoreStartZ();
    const int ce_z = mesh.GetCoreEndExclusiveZ();

    const int nx = ce_x - cs_x;
    const int ny = ce_y - cs_y;
    const int nz = ce_z - cs_z;

    if (nx <= 0 || ny <= 0 || nz <= 0) {
        throw std::runtime_error("Invalid core range");
    }

    const std::string filename =
        is_analytical_
            ? (rank_output_dir_ + "/step_" + (static_cast<std::ostringstream&&>(
                std::ostringstream() << std::setw(4) << std::setfill('0') << step)).str() + ".vtk")
            : GenerateFilename(settings.GetNx(), step, settings);

    vtkSmartPointer<vtkStructuredGrid> grid = vtkSmartPointer<vtkStructuredGrid>::New();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    grid->SetDimensions(nx, ny, nz);

    const vtkIdType num_points = static_cast<vtkIdType>(nx) * ny * nz;
    points->SetNumberOfPoints(num_points);

    const auto& xc = mesh.Xc();
    const auto& yc = mesh.Yc();
    const auto& zc = mesh.Zc();

    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                const double x = xc(static_cast<std::size_t>(cs_x + i));
                const double y = yc(static_cast<std::size_t>(cs_y + j));
                const double z = zc(static_cast<std::size_t>(cs_z + k));
                const vtkIdType pid = i + j * nx + k * nx * ny;
                points->SetPoint(pid, x, y, z);
            }
        }
    }
    grid->SetPoints(points);

    const auto& U = layer.U();

    auto arr_rho = vtkSmartPointer<vtkDoubleArray>::New();
    arr_rho->SetName("density");
    arr_rho->SetNumberOfComponents(1);
    arr_rho->SetNumberOfTuples(num_points);

    auto arr_vel = vtkSmartPointer<vtkDoubleArray>::New();
    arr_vel->SetName("velocity");
    arr_vel->SetNumberOfComponents(3);
    arr_vel->SetNumberOfTuples(num_points);

    auto arr_p = vtkSmartPointer<vtkDoubleArray>::New();
    arr_p->SetName("pressure");
    arr_p->SetNumberOfComponents(1);
    arr_p->SetNumberOfTuples(num_points);

    // Добавим вывод температуры, раз она у нас теперь честно считается
    auto arr_t = vtkSmartPointer<vtkDoubleArray>::New();
    arr_t->SetName("temperature");
    arr_t->SetNumberOfComponents(1);
    arr_t->SetNumberOfTuples(num_points);

    auto arr_e = vtkSmartPointer<vtkDoubleArray>::New();
    arr_e->SetName("conserved_energy");
    arr_e->SetNumberOfComponents(1);
    arr_e->SetNumberOfTuples(num_points);

    auto arr_eint = vtkSmartPointer<vtkDoubleArray>::New();
    arr_eint->SetName("internal_energy");
    arr_eint->SetNumberOfComponents(1);
    arr_eint->SetNumberOfTuples(num_points);

    const bool write_lambda = utils::ToLower(settings.solver) == "mader";

    vtkSmartPointer<vtkDoubleArray> arr_lambda;
    if (write_lambda) {
        arr_lambda = vtkSmartPointer<vtkDoubleArray>::New();
        arr_lambda->SetName("reactant_mass_fraction");
        arr_lambda->SetNumberOfComponents(1);
        arr_lambda->SetNumberOfTuples(num_points);
    }

    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                const int ii = cs_x + i;
                const int jj = cs_y + j;
                const int kk = cs_z + k;
                const vtkIdType pid = i + j * nx + k * nx * ny;

                if (mesh.IsSolidCell(ii, jj, kk)) {
                    grid->BlankPoint(pid);
                    const double nan = std::numeric_limits<double>::quiet_NaN();

                    arr_rho->SetValue(pid, nan);

                    double vec[3] = {nan, nan, nan};
                    arr_vel->SetTuple(pid, vec);

                    arr_p->SetValue(pid, nan);
                    arr_t->SetValue(pid, nan);
                    arr_e->SetValue(pid, nan);
                    arr_eint->SetValue(pid, nan);

                    if (write_lambda) {
                        arr_lambda->SetValue(pid, nan);
                    }
                    continue;
                }

                double lambda = 1.0;
                if (write_lambda) {
                    lambda = layer.ReactantMassFraction()(ii, jj, kk);
                }

                double rho = 0.0;
                double u = 0.0;
                double v = 0.0;
                double w = 0.0;
                double P = 0.0;
                double T = 0.0;
                ConservativeToPrimitive(U, ii, jj, kk, eos_, lambda, rho, u, v, w, P, T);

                const double E = U(DataLayer::k_E, ii, jj, kk);
                const double kinetic = 0.5 * rho * (u * u + v * v + w * w);
                const double eint = rho > 0.0 ? (E - kinetic) / rho : 0.0;

                arr_rho->SetValue(pid, rho);

                double vec[3] = {u, v, w};
                arr_vel->SetTuple(pid, vec);

                arr_p->SetValue(pid, P);
                arr_t->SetValue(pid, T);
                arr_e->SetValue(pid, E);
                arr_eint->SetValue(pid, eint);

                if (write_lambda) {
                    arr_lambda->SetValue(pid, lambda);
                }
            }
        }
    }

    grid->GetPointData()->AddArray(arr_rho);
    grid->GetPointData()->AddArray(arr_p);
    grid->GetPointData()->AddArray(arr_t);
    grid->GetPointData()->AddArray(arr_eint);
    grid->GetPointData()->AddArray(arr_e);
    grid->GetPointData()->AddArray(arr_vel);

    if (write_lambda) {
        grid->GetPointData()->AddArray(arr_lambda);
    }

    vtkSmartPointer<vtkDoubleArray> time_array = vtkSmartPointer<vtkDoubleArray>::New();
    time_array->SetName("TimeValue");
    time_array->SetNumberOfComponents(1);
    time_array->InsertNextValue(time);
    grid->GetFieldData()->AddArray(time_array);

    vtkSmartPointer<vtkStructuredGridWriter> writer = vtkSmartPointer<vtkStructuredGridWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData(grid);
    writer->SetFileTypeToBinary();
    vtkSmartPointer<vtkDoubleArray> arr_global_nx = vtkSmartPointer<vtkDoubleArray>::New();
    arr_global_nx->SetName("GlobalNx");
    arr_global_nx->InsertNextValue(mesh.GetGlobalNx());
    grid->GetFieldData()->AddArray(arr_global_nx);

    vtkSmartPointer<vtkDoubleArray> arr_global_ny = vtkSmartPointer<vtkDoubleArray>::New();
    arr_global_ny->SetName("GlobalNy");
    arr_global_ny->InsertNextValue(mesh.GetGlobalNy());
    grid->GetFieldData()->AddArray(arr_global_ny);

    vtkSmartPointer<vtkDoubleArray> arr_global_nz = vtkSmartPointer<vtkDoubleArray>::New();
    arr_global_nz->SetName("GlobalNz");
    arr_global_nz->InsertNextValue(mesh.GetGlobalNz());
    grid->GetFieldData()->AddArray(arr_global_nz);

    vtkSmartPointer<vtkDoubleArray> arr_offset_x = vtkSmartPointer<vtkDoubleArray>::New();
    arr_offset_x->SetName("OffsetX");
    arr_offset_x->InsertNextValue(mesh.GetOffsetX());
    grid->GetFieldData()->AddArray(arr_offset_x);

    vtkSmartPointer<vtkDoubleArray> arr_offset_y = vtkSmartPointer<vtkDoubleArray>::New();
    arr_offset_y->SetName("OffsetY");
    arr_offset_y->InsertNextValue(mesh.GetOffsetY());
    grid->GetFieldData()->AddArray(arr_offset_y);

    vtkSmartPointer<vtkDoubleArray> arr_offset_z = vtkSmartPointer<vtkDoubleArray>::New();
    arr_offset_z->SetName("OffsetZ");
    arr_offset_z->InsertNextValue(mesh.GetOffsetZ());
    grid->GetFieldData()->AddArray(arr_offset_z);

    writer->Write();
}


void VTKWriter::Write3D(const PressureVelocityState& state,
                        const Mesh& mesh,
                        const Settings& settings,
                        const std::size_t step,
                        const double time) const {
    const int cs_x = mesh.GetCoreStartX();
    const int ce_x = mesh.GetCoreEndExclusiveX();
    const int cs_y = mesh.GetCoreStartY();
    const int ce_y = mesh.GetCoreEndExclusiveY();
    const int cs_z = mesh.GetCoreStartZ();
    const int ce_z = mesh.GetCoreEndExclusiveZ();

    const int nx = ce_x - cs_x;
    const int ny = ce_y - cs_y;
    const int nz = ce_z - cs_z;

    if (nx <= 0 || ny <= 0 || nz <= 0) {
        throw std::runtime_error("Invalid core range");
    }

    const std::string filename =
        is_analytical_
            ? (rank_output_dir_ + "/step_" + (static_cast<std::ostringstream&&>(
                std::ostringstream() << std::setw(4) << std::setfill('0') << step)).str() + ".vtk")
            : GenerateFilename(settings.GetNx(), step, settings);

    vtkSmartPointer<vtkStructuredGrid> grid = vtkSmartPointer<vtkStructuredGrid>::New();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    grid->SetDimensions(nx, ny, nz);

    const vtkIdType num_points = static_cast<vtkIdType>(nx) * ny * nz;
    points->SetNumberOfPoints(num_points);

    const auto& xc = mesh.Xc();
    const auto& yc = mesh.Yc();
    const auto& zc = mesh.Zc();

    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                const double x = xc(static_cast<std::size_t>(cs_x + i));
                const double y = yc(static_cast<std::size_t>(cs_y + j));
                const double z = zc(static_cast<std::size_t>(cs_z + k));
                const vtkIdType pid = i + j * nx + k * nx * ny;
                points->SetPoint(pid, x, y, z);
            }
        }
    }
    grid->SetPoints(points);

    const auto& p_field = state.Pressure();

    auto arr_p = vtkSmartPointer<vtkDoubleArray>::New();
    arr_p->SetName("pressure");
    arr_p->SetNumberOfComponents(1);
    arr_p->SetNumberOfTuples(num_points);

    auto arr_vel = vtkSmartPointer<vtkDoubleArray>::New();
    arr_vel->SetName("velocity");
    arr_vel->SetNumberOfComponents(3);
    arr_vel->SetNumberOfTuples(num_points);

    auto arr_u = vtkSmartPointer<vtkDoubleArray>::New();
    arr_u->SetName("u");
    arr_u->SetNumberOfComponents(1);
    arr_u->SetNumberOfTuples(num_points);

    auto arr_v = vtkSmartPointer<vtkDoubleArray>::New();
    arr_v->SetName("v");
    arr_v->SetNumberOfComponents(1);
    arr_v->SetNumberOfTuples(num_points);

    auto arr_w = vtkSmartPointer<vtkDoubleArray>::New();
    arr_w->SetName("w");
    arr_w->SetNumberOfComponents(1);
    arr_w->SetNumberOfTuples(num_points);

    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                const int ii = cs_x + i;
                const int jj = cs_y + j;
                const int kk = cs_z + k;
                const vtkIdType pid = i + j * nx + k * nx * ny;

                if (mesh.IsSolidCell(ii, jj, kk)) {
                    grid->BlankPoint(pid);
                    const double nan = std::numeric_limits<double>::quiet_NaN();

                    arr_p->SetValue(pid, nan);

                    double vec[3] = {nan, nan, nan};
                    arr_vel->SetTuple(pid, vec);

                    arr_u->SetValue(pid, nan);
                    arr_v->SetValue(pid, nan);
                    arr_w->SetValue(pid, nan);
                    continue;
                }

                const double p = p_field(ii, jj, kk);

                double u = 0.0;
                double v = 0.0;
                double w = 0.0;
                PressureVelocityCellCenteredVelocity(state, ii, jj, kk, u, v, w);

                arr_p->SetValue(pid, p);

                double vec[3] = {u, v, w};
                arr_vel->SetTuple(pid, vec);

                arr_u->SetValue(pid, u);
                arr_v->SetValue(pid, v);
                arr_w->SetValue(pid, w);
            }
        }
    }

    grid->GetPointData()->AddArray(arr_p);
    grid->GetPointData()->AddArray(arr_vel);
    grid->GetPointData()->AddArray(arr_u);
    grid->GetPointData()->AddArray(arr_v);
    grid->GetPointData()->AddArray(arr_w);

    vtkSmartPointer<vtkDoubleArray> time_array = vtkSmartPointer<vtkDoubleArray>::New();
    time_array->SetName("TimeValue");
    time_array->SetNumberOfComponents(1);
    time_array->InsertNextValue(time);
    grid->GetFieldData()->AddArray(time_array);

    vtkSmartPointer<vtkDoubleArray> arr_global_nx = vtkSmartPointer<vtkDoubleArray>::New();
    arr_global_nx->SetName("GlobalNx");
    arr_global_nx->InsertNextValue(mesh.GetGlobalNx());
    grid->GetFieldData()->AddArray(arr_global_nx);

    vtkSmartPointer<vtkDoubleArray> arr_global_ny = vtkSmartPointer<vtkDoubleArray>::New();
    arr_global_ny->SetName("GlobalNy");
    arr_global_ny->InsertNextValue(mesh.GetGlobalNy());
    grid->GetFieldData()->AddArray(arr_global_ny);

    vtkSmartPointer<vtkDoubleArray> arr_global_nz = vtkSmartPointer<vtkDoubleArray>::New();
    arr_global_nz->SetName("GlobalNz");
    arr_global_nz->InsertNextValue(mesh.GetGlobalNz());
    grid->GetFieldData()->AddArray(arr_global_nz);

    vtkSmartPointer<vtkDoubleArray> arr_offset_x = vtkSmartPointer<vtkDoubleArray>::New();
    arr_offset_x->SetName("OffsetX");
    arr_offset_x->InsertNextValue(mesh.GetOffsetX());
    grid->GetFieldData()->AddArray(arr_offset_x);

    vtkSmartPointer<vtkDoubleArray> arr_offset_y = vtkSmartPointer<vtkDoubleArray>::New();
    arr_offset_y->SetName("OffsetY");
    arr_offset_y->InsertNextValue(mesh.GetOffsetY());
    grid->GetFieldData()->AddArray(arr_offset_y);

    vtkSmartPointer<vtkDoubleArray> arr_offset_z = vtkSmartPointer<vtkDoubleArray>::New();
    arr_offset_z->SetName("OffsetZ");
    arr_offset_z->InsertNextValue(mesh.GetOffsetZ());
    grid->GetFieldData()->AddArray(arr_offset_z);

    vtkSmartPointer<vtkStructuredGridWriter> writer = vtkSmartPointer<vtkStructuredGridWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData(grid);
    writer->SetFileTypeToBinary();
    writer->Write();
}
