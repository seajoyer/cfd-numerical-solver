#include "output/VTKWriter.hpp"

#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkFieldData.h>
#include <vtkPoints.h>
#include <vtkStructuredGrid.h>
#include <vtkStructuredGridWriter.h>

#include <cmath>
#include <cstddef>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include "data/DataLayer.hpp"

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
    std::filesystem::create_directories(output_dir_);
}

VTKWriter::~VTKWriter() = default;

auto VTKWriter::RequiresFinalization() const -> bool {
    return false;
}

auto VTKWriter::Finalize(const Settings&) -> std::string {
    return "";
}

auto VTKWriter::GenerateFilename(int N, std::size_t step, const Settings& settings) const
    -> std::string {
    std::ostringstream oss;

    if (is_analytical_) {
        oss << output_dir_ << "/step_" << std::setw(4) << std::setfill('0') << step << ".vtk";
    } else {
        oss << output_dir_ << "/" << settings.solver << "__R_" << settings.reconstruction
            << "__N_" << N << "__CFL_" << std::fixed << std::setprecision(1) << settings.cfl
            << "__step_" << std::setw(4) << std::setfill('0') << step << ".vtk";
    }

    return oss.str();
}

void VTKWriter::Write(const DataLayer& layer, const Settings& settings, std::size_t step,
                      double time) const {
    Write3D(layer, settings, step, time);
}

void VTKWriter::Write(const DataLayer& layer, const DataLayer* analytical_layer,
                      const Settings& settings, std::size_t step, double time) const {
    Write(layer, settings, step, time);
}

static inline void ConservativeToPrimitive(
    const xt::xtensor<double, 4>& U, int i, int j, int k, double gamma,
    double& rho, double& u, double& v, double& w, double& P
) {
    rho = U(DataLayer::k_rho, i, j, k);
    if (rho <= 0.0) {
        rho = 0.0;
        u = v = w = 0.0;
        P = 0.0;
        return;
    }

    const double inv_rho = 1.0 / rho;
    u = U(DataLayer::k_rhoU, i, j, k) * inv_rho;
    v = U(DataLayer::k_rhoV, i, j, k) * inv_rho;
    w = U(DataLayer::k_rhoW, i, j, k) * inv_rho;

    const double E = U(DataLayer::k_E, i, j, k);
    const double kinetic = 0.5 * rho * (u*u + v*v + w*w);
    const double eint = E - kinetic;
    P = (gamma - 1.0) * eint;
}

void VTKWriter::Write3D(const DataLayer& layer, const Settings& settings,
                        std::size_t step, double time) const {
    const int cs_x = layer.GetCoreStartX();
    const int ce_x = layer.GetCoreEndExclusiveX();
    const int cs_y = layer.GetCoreStartY();
    const int ce_y = layer.GetCoreEndExclusiveY();
    const int cs_z = layer.GetCoreStartZ();
    const int ce_z = layer.GetCoreEndExclusiveZ();

    const int nx = ce_x - cs_x;
    const int ny = ce_y - cs_y;
    const int nz = ce_z - cs_z;

    if (nx <= 0 || ny <= 0 || nz <= 0) {
        throw std::runtime_error("Invalid core range");
    }

    std::ostringstream oss;
    if (is_analytical_) {
        oss << output_dir_ << "/step_" << std::setw(4) << std::setfill('0') << step << ".vtk";
    } else {
        oss << output_dir_ << "/"
            << settings.solver << "__R_" << settings.reconstruction
            << "__N_" << settings.GetNx() << "x" << settings.GetNy() << "x" << settings.GetNz()
            << "__CFL_" << std::fixed << std::setprecision(1) << settings.cfl
            << "__step_" << std::setw(4) << std::setfill('0') << step << ".vtk";
    }
    const std::string filename = oss.str();

    vtkSmartPointer<vtkStructuredGrid> grid = vtkSmartPointer<vtkStructuredGrid>::New();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    grid->SetDimensions(nx, ny, nz);

    const vtkIdType num_points = static_cast<vtkIdType>(nx) * ny * nz;
    points->SetNumberOfPoints(num_points);

    const auto& xc = layer.Xc();
    const auto& yc = layer.Yc();
    const auto& zc = layer.Zc();

    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                const double x = xc(static_cast<std::size_t>(cs_x + i));
                const double y = yc(static_cast<std::size_t>(cs_y + j));
                const double z = zc(static_cast<std::size_t>(cs_z + k));
                const vtkIdType pid = static_cast<vtkIdType>(i + j * nx + k * nx * ny);
                points->SetPoint(pid, x, y, z);
            }
        }
    }
    grid->SetPoints(points);

    const double gamma = settings.gamma;
    const auto& U = layer.U();

    auto add_scalar_field = [&](auto&& accessor, const char* name) {
        auto array = vtkSmartPointer<vtkDoubleArray>::New();
        array->SetName(name);
        array->SetNumberOfComponents(1);
        array->SetNumberOfTuples(num_points);

        for (int k = 0; k < nz; ++k) {
            for (int j = 0; j < ny; ++j) {
                for (int i = 0; i < nx; ++i) {
                    const int ii = cs_x + i;
                    const int jj = cs_y + j;
                    const int kk = cs_z + k;

                    const double val = accessor(ii, jj, kk);
                    const vtkIdType pid = static_cast<vtkIdType>(i + j * nx + k * nx * ny);
                    array->SetValue(pid, val);
                }
            }
        }
        grid->GetPointData()->AddArray(array);
    };

    add_scalar_field([&](int i, int j, int k) {
        return U(DataLayer::k_rho, i, j, k);
    }, "density");

    add_scalar_field([&](int i, int j, int k) {
        double rho, u, v, w, P;
        ConservativeToPrimitive(U, i, j, k, gamma, rho, u, v, w, P);
        return u;
    }, "velocity_x");

    add_scalar_field([&](int i, int j, int k) {
        double rho, u, v, w, P;
        ConservativeToPrimitive(U, i, j, k, gamma, rho, u, v, w, P);
        return v;
    }, "velocity_y");

    add_scalar_field([&](int i, int j, int k) {
        double rho, u, v, w, P;
        ConservativeToPrimitive(U, i, j, k, gamma, rho, u, v, w, P);
        return w;
    }, "velocity_z");

    add_scalar_field([&](int i, int j, int k) {
        double rho, u, v, w, P;
        ConservativeToPrimitive(U, i, j, k, gamma, rho, u, v, w, P);
        return std::sqrt(u*u + v*v + w*w);
    }, "velocity_magnitude");

    add_scalar_field([&](int i, int j, int k) {
        double rho, u, v, w, P;
        ConservativeToPrimitive(U, i, j, k, gamma, rho, u, v, w, P);
        return P;
    }, "pressure");

    add_scalar_field([&](int i, int j, int k) {
        return U(DataLayer::k_E, i, j, k);
    }, "energy");

    vtkSmartPointer<vtkDoubleArray> velocity_vectors = vtkSmartPointer<vtkDoubleArray>::New();
    velocity_vectors->SetName("velocity");
    velocity_vectors->SetNumberOfComponents(3);
    velocity_vectors->SetNumberOfTuples(num_points);

    for (int k = 0; k < nz; ++k) {
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                const int ii = cs_x + i;
                const int jj = cs_y + j;
                const int kk = cs_z + k;

                double rho, u, v, w, P;
                ConservativeToPrimitive(U, ii, jj, kk, gamma, rho, u, v, w, P);

                double vec[3] = {u, v, w};
                const vtkIdType pid = static_cast<vtkIdType>(i + j * nx + k * nx * ny);
                velocity_vectors->SetTuple(pid, vec);
            }
        }
    }
    grid->GetPointData()->AddArray(velocity_vectors);

    vtkSmartPointer<vtkDoubleArray> time_array = vtkSmartPointer<vtkDoubleArray>::New();
    time_array->SetName("TimeValue");
    time_array->SetNumberOfComponents(1);
    time_array->InsertNextValue(time);
    grid->GetFieldData()->AddArray(time_array);

    vtkSmartPointer<vtkStructuredGridWriter> writer = vtkSmartPointer<vtkStructuredGridWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData(grid);
    writer->SetFileTypeToBinary();
    writer->Write();
}