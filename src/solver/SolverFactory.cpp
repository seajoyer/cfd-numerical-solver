#include "solver/SolverFactory.hpp"
#include "solver/GodunovSolver.hpp"
#include <algorithm>
#include <cctype>
#include <stdexcept>

auto SolverFactory::Create(Settings &settings) -> std::unique_ptr<Solver> {
    // Convert to lowercase for case-insensitive comparison
    std::string type_lower = settings.solver;
    std::transform(type_lower.begin(), type_lower.end(), type_lower.begin(),
                   [](unsigned char c) -> int { return std::tolower(c); });

    if (type_lower == "godunov") {
        return std::make_unique<GodunovSolver>(settings);
    }

    // Future solvers can be added here:
    // 
    // if (type_lower == "muscl") {
    //     return std::make_unique<MUSCLSolver>(dim);
    // }
    // 
    // if (type_lower == "weno") {
    //     return std::make_unique<WENOSolver>(dim);
    // }
    //
    // etc ...

    throw std::runtime_error("Unknown solver type: " + settings.solver);
}
