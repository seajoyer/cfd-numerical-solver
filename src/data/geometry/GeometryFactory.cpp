#ifndef GEOMETRYFACTORY_HPP
#define GEOMETRYFACTORY_HPP

#include <memory>
#include <stdexcept>
#include <string>

#include "config/Settings.hpp"
#include "data/geometry/CircleGeometry.hpp"
#include "data/geometry/GeometryPrimitive.hpp"
#include "data/geometry/RectangleGeometry.hpp"

inline auto CreateGeometryPrimitive(const ImmersedObjectSettings& object)
    -> std::shared_ptr<GeometryPrimitive> {
    if (object.type == "circle") {
        return std::make_shared<CircleGeometry>(object.cx, object.cy, object.radius);
    }

    if (object.type == "rectangle") {
        return std::make_shared<RectangleGeometry>(object.cx, object.cy, object.size_x, object.size_y);
    }

    throw std::runtime_error("Unknown immersed object type: " + object.type);
}

#endif  // GEOMETRYFACTORY_HPP