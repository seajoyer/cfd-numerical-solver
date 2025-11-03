#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <string>

struct Config {
    int N = 0;
    int padding = 2;
    int dim = 1;
    double CFL = 0.5;
    double t_end = 0.2;
    std::string output_dir = "data/output";
    std::size_t output_every = 1;
};

#endif // CONFIG_HPP
