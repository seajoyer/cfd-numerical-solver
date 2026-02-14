#include "output/VTKWriter.hpp"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

VTKWriter::VTKWriter(const std::string& output_dir, bool is_analytical)
    : output_dir_(output_dir), is_analytical_(is_analytical) {
    std::filesystem::create_directories(output_dir);
}

auto VTKWriter::RequiresFinalization() const -> bool { return false; }
auto VTKWriter::Finalize(const Settings&) -> std::string { return ""; }

void VTKWriter::Write(const DataLayer& layer, const Settings& settings,
                      std::size_t step, double time) const {
    if (settings.dim >= 2) {
        Write2D(layer, settings, step, time);
        return;
    }

    // --- Original 1D VTK writing ---
    std::ostringstream oss;
    if (is_analytical_) {
        oss << output_dir_ << "/step_" << std::setw(4) << std::setfill('0') << step << ".vtk";
    } else {
        oss << output_dir_ << "/"
            << settings.solver << "__R_" << settings.reconstruction
            << "__N_" << settings.N
            << "__CFL_" << std::fixed << std::setprecision(1) << settings.cfl
            << "__step_" << std::setw(4) << std::setfill('0') << step << ".vtk";
    }
    std::string filepath = oss.str();

    std::ofstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "VTKWriter: Failed to open " << filepath << '\n';
        return;
    }

    const int core_start = layer.GetCoreStart(0);
    const int core_end = layer.GetCoreEndExclusive(0);
    const int n = core_end - core_start;

    file << "# vtk DataFile Version 3.0\n";
    file << "CFD Simulation step=" << step << " time=" << time << "\n";
    file << "ASCII\n";
    file << "DATASET STRUCTURED_POINTS\n";
    file << "DIMENSIONS " << n << " 1 1\n";
    file << "ORIGIN " << layer.xc(core_start) << " 0 0\n";
    double dx = (n > 1) ? (layer.xc(core_start + 1) - layer.xc(core_start)) : 1.0;
    file << "SPACING " << dx << " 1 1\n";
    file << "POINT_DATA " << n << "\n";

    // Density
    file << "SCALARS density double 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int i = core_start; i < core_end; ++i) file << layer.rho(i) << "\n";

    // Velocity
    file << "SCALARS velocity double 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int i = core_start; i < core_end; ++i) file << layer.u(i) << "\n";

    // Pressure
    file << "SCALARS pressure double 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int i = core_start; i < core_end; ++i) file << layer.P(i) << "\n";

    // Energy
    file << "SCALARS energy double 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int i = core_start; i < core_end; ++i) file << layer.e(i) << "\n";

    // Internal Energy
    file << "SCALARS internal_energy double 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int i = core_start; i < core_end; ++i) file << layer.U(i) << "\n";

    file.close();
}

void VTKWriter::Write(const DataLayer& layer, const DataLayer* analytical_layer,
                      const Settings& settings, std::size_t step, double time) const {
    Write(layer, settings, step, time);
}

void VTKWriter::Write2D(const DataLayer& layer, const Settings& settings,
                        std::size_t step, double time) const {
    std::ostringstream oss;
    oss << output_dir_ << "/"
        << settings.solver << "__step_" << std::setw(4) << std::setfill('0') << step << ".vtk";
    std::string filepath = oss.str();

    std::ofstream file(filepath);
    if (!file.is_open()) {
        std::cerr << "VTKWriter: Failed to open " << filepath << '\n';
        return;
    }

    const int cs_x = layer.GetCoreStart(0);
    const int ce_x = layer.GetCoreEndExclusive(0);
    const int cs_y = layer.GetCoreStart(1);
    const int ce_y = layer.GetCoreEndExclusive(1);
    const int nx = ce_x - cs_x;
    const int ny = ce_y - cs_y;
    const int total_points = nx * ny;

    const double dx = settings.L_x / static_cast<double>(nx);
    const double dy = settings.L_y / static_cast<double>(ny);

    file << "# vtk DataFile Version 3.0\n";
    file << "CFD 2D Simulation step=" << step << " time=" << time << "\n";
    file << "ASCII\n";
    file << "DATASET STRUCTURED_POINTS\n";
    file << "DIMENSIONS " << nx << " " << ny << " 1\n";
    file << "ORIGIN " << layer.xc(cs_x) << " " << layer.yc(cs_y) << " 0\n";
    file << "SPACING " << dx << " " << dy << " 1\n";
    file << "POINT_DATA " << total_points << "\n";

    // Density
    file << "SCALARS density double 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int j = cs_y; j < ce_y; ++j)
        for (int i = cs_x; i < ce_x; ++i)
            file << layer.rho(i, j) << "\n";

    // Velocity-X
    file << "SCALARS velocity_x double 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int j = cs_y; j < ce_y; ++j)
        for (int i = cs_x; i < ce_x; ++i)
            file << layer.u(i, j) << "\n";

    // Velocity-Y
    file << "SCALARS velocity_y double 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int j = cs_y; j < ce_y; ++j)
        for (int i = cs_x; i < ce_x; ++i)
            file << layer.v(i, j) << "\n";

    // Velocity magnitude
    file << "SCALARS velocity_magnitude double 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int j = cs_y; j < ce_y; ++j)
        for (int i = cs_x; i < ce_x; ++i) {
            double vmag = std::sqrt(layer.u(i,j)*layer.u(i,j) + layer.v(i,j)*layer.v(i,j));
            file << vmag << "\n";
        }

    // Pressure
    file << "SCALARS pressure double 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int j = cs_y; j < ce_y; ++j)
        for (int i = cs_x; i < ce_x; ++i)
            file << layer.P(i, j) << "\n";

    // Energy
    file << "SCALARS energy double 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int j = cs_y; j < ce_y; ++j)
        for (int i = cs_x; i < ce_x; ++i)
            file << layer.e(i, j) << "\n";

    // Velocity vector
    file << "VECTORS velocity double\n";
    for (int j = cs_y; j < ce_y; ++j)
        for (int i = cs_x; i < ce_x; ++i)
            file << layer.u(i, j) << " " << layer.v(i, j) << " 0\n";

    file.close();
}
