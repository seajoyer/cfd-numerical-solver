#include "solver/SolverFactory.hpp"
#include "solver/GodunovKolganRodionovSolver.hpp"
#include "solver/GodunovSolver.hpp"
#include "solver/AnalyticalSolver.hpp"
#include <algorithm>
#include <cctype>
#include <stdexcept>

auto SolverFactory::Create(Settings &settings) -> std::unique_ptr<Solver> {
    // Convert to lowercase for case-insensitive comparison
    std::string type_lower = settings.solver;
    std::transform(type_lower.begin(), type_lower.end(), type_lower.begin(),
                   [](unsigned char c) -> int { return std::tolower(c); });

    if (type_lower == "analytical") {
        return std::make_unique<AnalyticalSolver>(settings);
    }
    if (type_lower == "godunov") {
        settings.reconstruction = "p0";
        return std::make_unique<GodunovSolver>(settings);
    }
    if (type_lower == "godunov-kolgan") {
        return std::make_unique<GodunovSolver>(settings);
    }
    if (type_lower == "godunov-kolgan-rodionov") {
        return std::make_unique<GodunovKolganRodionovSolver>(settings);
    }

    throw std::runtime_error("Unknown solver type: " + settings.solver);
}
