#ifndef SETTINGS_HPP
#define SETTINGS_HPP

#include <string>
#include <vector>

struct Settings {
    std::string solver = "godunov";

    int N = 0;
    double CFL = 0.5;
    double t_end = 0.2;
    int padding = 2;
    double c = 340;
    double gamma = 1.4;
    int dim = 1;
    double L_x = 10;
    double L_y = 10;
    double L_z = 10;

    std::size_t output_every_steps = 1;
    std::vector<std::string> output_types = {"vtk"};
    std::string output_dir = "data/output";
};

#endif // SETTINGS_HPP
