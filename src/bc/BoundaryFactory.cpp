#include "bc/BoundaryFactory.hpp"

#include <algorithm>
#include <cctype>
#include <iostream>
#include <stdexcept>

#include "bc/FreeStreamBoundary.hpp"
#include "bc/InletBoundary.hpp"
#include "bc/NonReflectiveBoundary.hpp"
#include "bc/OutletBoundary.hpp"
#include "bc/PeriodicBoundary.hpp"
#include "bc/ReflectiveBoundary.hpp"
#include "bc/SymmetryBoundary.hpp"
#include "bc/WallBoundary.hpp"

auto BoundaryFactory::Create(const std::string& type, double rho_inf,
                              double u_inf, double p_inf)
    -> std::shared_ptr<BoundaryCondition> {

    std::string lower = type;
    std::transform(lower.begin(), lower.end(), lower.begin(),
                   [](unsigned char c) { return std::tolower(c); });

    if (lower == "free_stream" || lower == "freestream") {
        return std::make_shared<FreeStreamBoundary>(rho_inf, u_inf, p_inf);
    }
    if (lower == "reflective") {
        return std::make_shared<ReflectiveBoundary>();
    }
    if (lower == "outlet") {
        return std::make_shared<OutletBoundary>();
    }
    if (lower == "inlet") {
        return std::make_shared<InletBoundary>(rho_inf, u_inf, p_inf);
    }
    if (lower == "periodic") {
        return std::make_shared<PeriodicBoundary>();
    }
    if (lower == "symmetry") {
        return std::make_shared<SymmetryBoundary>();
    }
    if (lower == "wall") {
        return std::make_shared<WallBoundary>();
    }
    if (lower == "non_reflective" || lower == "nonreflective") {
        return std::make_shared<NonReflectiveBoundary>();
    }

    std::cerr << "Warning: Unknown boundary type '" << type
              << "', defaulting to free_stream\n";
    return std::make_shared<FreeStreamBoundary>(rho_inf, u_inf, p_inf);
}

auto BoundaryFactory::Create2D(const std::string& type, double rho_inf,
                                double u_inf, double v_inf, double p_inf)
    -> std::shared_ptr<BoundaryCondition> {

    std::string lower = type;
    std::transform(lower.begin(), lower.end(), lower.begin(),
                   [](unsigned char c) { return std::tolower(c); });

    if (lower == "free_stream" || lower == "freestream") {
        return std::make_shared<FreeStreamBoundary>(rho_inf, u_inf, v_inf, p_inf);
    }
    if (lower == "reflective") {
        return std::make_shared<ReflectiveBoundary>();
    }
    if (lower == "outlet") {
        return std::make_shared<OutletBoundary>();
    }
    if (lower == "inlet") {
        return std::make_shared<InletBoundary>(rho_inf, u_inf, v_inf, p_inf);
    }
    if (lower == "periodic") {
        return std::make_shared<PeriodicBoundary>();
    }
    if (lower == "symmetry") {
        return std::make_shared<SymmetryBoundary>();
    }
    if (lower == "wall") {
        return std::make_shared<WallBoundary>();
    }
    if (lower == "non_reflective" || lower == "nonreflective") {
        return std::make_shared<NonReflectiveBoundary>();
    }

    std::cerr << "Warning: Unknown boundary type '" << type
              << "', defaulting to free_stream\n";
    return std::make_shared<FreeStreamBoundary>(rho_inf, u_inf, v_inf, p_inf);
}
