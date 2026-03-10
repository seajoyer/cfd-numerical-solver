#include "data/geometry/CircleGeometry.hpp"

#include <cmath>
#include <stdexcept>

CircleGeometry::CircleGeometry(const double cx, const double cy, const double radius)
    : cx_(cx), cy_(cy), radius_(radius) {
    if (!(radius_ > 0.0)) {
        throw std::invalid_argument("CircleGeometry: radius must be > 0");
    }
}

bool CircleGeometry::Contains(const GeometryPoint& p) const {
    const double dx = p.x - cx_;
    const double dy = p.y - cy_;
    return dx * dx + dy * dy < radius_ * radius_;
}

double CircleGeometry::SignedDistance(const GeometryPoint& p) const {
    const double dx = p.x - cx_;
    const double dy = p.y - cy_;
    return std::sqrt(dx * dx + dy * dy) - radius_;
}

GeometryVector CircleGeometry::OutwardNormal(const GeometryPoint& p) const {
    const double dx = p.x - cx_;
    const double dy = p.y - cy_;
    const double norm = std::sqrt(dx * dx + dy * dy);

    if (norm > 0.0) {
        return GeometryVector{dx / norm, dy / norm, 0.0};
    }

    return GeometryVector{1.0, 0.0, 0.0};
}

int CircleGeometry::GetDimension() const {
    return 2;
}

double CircleGeometry::GetCenterX() const {
    return cx_;
}

double CircleGeometry::GetCenterY() const {
    return cy_;
}

double CircleGeometry::GetRadius() const {
    return radius_;
}
