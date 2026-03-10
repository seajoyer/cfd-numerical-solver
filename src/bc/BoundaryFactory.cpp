#include "bc/BoundaryFactory.hpp"

#include "bc/BoundaryCondition.hpp"
#include "bc/OutletBoundary.hpp"
// #include "bc/ReflectiveBoundary.hpp"
// #include "bc/InletBoundary.hpp"
// #include "bc/NonReflectiveBoundary.hpp"
// #include "bc/SymmetryBoundary.hpp"
// #include "bc/PeriodicBoundary.hpp"
// #include "bc/WallBoundary.hpp"
#include "bc/FreeStreamBoundary.hpp"

#include <stdexcept>


auto BoundaryFactory::Create(const std::string& boundary_type) -> std::shared_ptr<BoundaryCondition> {
    if (boundary_type == "outlet") return std::make_shared<OutletBoundary>();
    // if (boundary_type == "reflective") return std::make_shared<ReflectiveBoundary>();
    // if (boundary_type == "non_reflective") return std::make_shared<NonReflectiveBoundary>();
    // if (boundary_type == "periodic") return std::make_shared<PeriodicBoundary>();
    // if (boundary_type == "symmetry") return std::make_shared<SymmetryBoundary>();
    // if (boundary_type == "wall") return std::make_shared<WallBoundary>();

    throw std::runtime_error("Unknown boundary condition type: " + boundary_type);
}

auto BoundaryFactory::Create(const std::string& boundary_type, const FarfieldConservative& farfield_U)
    -> std::shared_ptr<BoundaryCondition> {
    if (boundary_type == "free_stream") return std::make_shared<FreeStreamBoundary>(farfield_U);
    // if (boundary_type == "inlet") return std::make_shared<InletBoundary>(farfield_U);
    return Create(boundary_type);
}
