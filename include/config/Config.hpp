#ifndef CONFIG_HPP
#define CONFIG_HPP

#include <string>

struct Config {
    int n = 0;
    int padding = 2;
    int dim = 1;
    double cfl = 0.5;
    double tEnd = 0.2;
    std::string outputDir = "data/output";
    std::size_t outputEvery = 1;
};

#endif // CONFIG_HPP