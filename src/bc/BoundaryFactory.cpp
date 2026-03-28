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
#include "utils/StringUtils.hpp"

PrimitiveCell BoundaryFactory::PrimitiveFromBoundaryState(
    const BoundaryStateSettings& state
) {
    PrimitiveCell primitive;
    primitive.rho = state.rho;
    primitive.u = state.u;
    primitive.v = state.v;
    primitive.w = state.w;
    primitive.P = state.p;
    return primitive;
}

std::shared_ptr<BoundaryCondition> BoundaryFactory::Create(
    const BoundaryConditionSettings& boundary_settings,
    const Settings& settings
) {
    (void)settings;

    const std::string type = utils::ToLower(boundary_settings.type);

    if (type == "outlet") {
        return std::make_shared<OutletBoundary>();
    }

    if (type == "reflective") {
        return std::make_shared<ReflectiveBoundary>();
    }

    if (type == "symmetry") {
        return std::make_shared<SymmetryBoundary>();
    }

    if (type == "wall") {
        return std::make_shared<WallBoundary>();
    }

    if (type == "periodic") {
        return std::make_shared<PeriodicBoundary>();
    }

    if (type == "free_stream") {
        if (!boundary_settings.state.has_value()) {
            throw std::runtime_error("BoundaryFactory: free_stream requires boundary state");
        }
        return std::make_shared<FreeStreamBoundary>(
                                                    PrimitiveFromBoundaryState(*boundary_settings.state)
                                                   );
    }

    if (type == "inlet") {
        if (!boundary_settings.state.has_value()) {
            throw std::runtime_error("BoundaryFactory: inlet requires boundary state");
        }
        return std::make_shared<InletBoundary>(
                                               PrimitiveFromBoundaryState(*boundary_settings.state)
                                              );
    }

    if (type == "non_reflective") {
        if (!boundary_settings.state.has_value()) {
            throw std::runtime_error("BoundaryFactory: non_reflective requires boundary state");
        }
        return std::make_shared<NonReflectiveBoundary>(
                                                       PrimitiveFromBoundaryState(*boundary_settings.state),
                                                       settings.gamma
                                                      );
    }

    throw std::runtime_error("BoundaryFactory: unknown boundary condition type: " + type);
}
