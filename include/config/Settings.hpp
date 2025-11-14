#ifndef SETTINGS_HPP
#define SETTINGS_HPP

#include <cstddef>
#include <string>

struct Settings {
    std::string solver = "godunov";
    std::string riemann_solver = "exact";
    std::string reconstruction = "p0";
    std::string left_boundary  = "free_stream";
    std::string right_boundary = "free_stream";
    double Q_user = 2.0;

    int N = 200;
    double cfl = 0.5;
    int padding = 2;
    double c = 340;
    double gamma = 1.4;
    int dim = 1;
    double L_x = 10;
    double L_y = 10;
    double L_z = 10;

    double t_end = 0;
    std::size_t step_end = 10000;

    std::size_t log_every_steps = 1;
    double      log_every_time  = 0.0;

    std::size_t output_every_steps = 1;
    double      output_every_time  = 0.0;

    std::string output_format = "vtk";
    std::string output_dir = "data/output";

    int sod_test_num = 1;
};

#endif // SETTINGS_HPP
