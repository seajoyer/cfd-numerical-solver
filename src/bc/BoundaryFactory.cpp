#include "bc/BoundaryFactory.hpp"

#include <memory>
#include <stdexcept>

#include "bc/BoundaryCondition.hpp"
#include "bc/FreeStreamBoundary.hpp"
#include "bc/InletBoundary.hpp"
#include "bc/NonReflectiveBoundary.hpp"
#include "bc/OutletBoundary.hpp"
#include "bc/PeriodicBoundary.hpp"
#include "bc/ReflectiveBoundary.hpp"
#include "bc/SymmetryBoundary.hpp"
#include "bc/WallBoundary.hpp"
#include "parallel/MPIContext.hpp"
#include "utils/StringUtils.hpp"

auto BoundaryFactory::Create(const std::string& boundary_type) -> std::shared_ptr<BoundaryCondition> {
    const std::string type = utils::ToLower(boundary_type);

    if (type == "outlet") return std::make_shared<OutletBoundary>();
    if (type == "reflective") return std::make_shared<ReflectiveBoundary>();
    if (type == "symmetry") return std::make_shared<SymmetryBoundary>();
    if (type == "wall") return std::make_shared<WallBoundary>();
    if (type == "periodic") return std::make_shared<PeriodicBoundary>();

    throw std::runtime_error("Unknown boundary condition type: " + boundary_type);
}

auto BoundaryFactory::Create(const std::string& boundary_type,
                             const FarfieldConservative& farfield_U,
                             const Settings& settings,
                             const int mpi_size) -> std::shared_ptr<BoundaryCondition> {
    const std::string type = utils::ToLower(boundary_type);
    const bool multi_rank_mpi = settings.mpi_enabled && mpi_size > 1;

    if (type == "free_stream") {
        return std::make_shared<FreeStreamBoundary>(farfield_U);
    }

    if (type == "inlet") {
        return std::make_shared<InletBoundary>(farfield_U);
    }

    if (type == "non_reflective") {
        return std::make_shared<NonReflectiveBoundary>(farfield_U, settings.gamma);
    }

    if (type == "periodic") {
        if (multi_rank_mpi) {
            return nullptr;
        }
        return std::make_shared<PeriodicBoundary>();
    }

    return Create(type);
}
