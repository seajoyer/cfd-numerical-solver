#include "data/geometry/RectangleGeometry.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

RectangleGeometry::RectangleGeometry(const double cx,
                                     const double cy,
                                     const double size_x,
                                     const double size_y)
    : cx_(cx), cy_(cy), hx_(0.5 * size_x), hy_(0.5 * size_y) {
    if (!(size_x > 0.0)) {
        throw std::invalid_argument("RectangleGeometry: size_x must be > 0");
    }
    if (!(size_y > 0.0)) {
        throw std::invalid_argument("RectangleGeometry: size_y must be > 0");
    }
}

bool RectangleGeometry::Contains(const GeometryPoint& p) const {
    const double qx = std::abs(p.x - cx_);
    const double qy = std::abs(p.y - cy_);
    return qx < hx_ && qy < hy_;
}

double RectangleGeometry::SignedDistance(const GeometryPoint& p) const {
    const double qx = std::abs(p.x - cx_) - hx_;
    const double qy = std::abs(p.y - cy_) - hy_;

    const double ox = std::max(qx, 0.0);
    const double oy = std::max(qy, 0.0);
    const double outside = std::sqrt(ox * ox + oy * oy);
    const double inside = std::min(std::max(qx, qy), 0.0);

    return outside + inside;
}

GeometryVector RectangleGeometry::OutwardNormal(const GeometryPoint& p) const {
    const double dx = p.x - cx_;
    const double dy = p.y - cy_;

    const double ax = std::abs(dx);
    const double ay = std::abs(dy);

    const double qx = ax - hx_;
    const double qy = ay - hy_;

    if (qx > 0.0 || qy > 0.0) {
        double nx = qx > 0.0 ? (dx >= 0.0 ? qx : -qx) : 0.0;
        double ny = qy > 0.0 ? (dy >= 0.0 ? qy : -qy) : 0.0;
        const double norm = std::sqrt(nx * nx + ny * ny);

        if (norm > 0.0) {
            return GeometryVector{nx / norm, ny / norm, 0.0};
        }
    }

    if ((hx_ - ax) <= (hy_ - ay)) {
        return GeometryVector{dx >= 0.0 ? 1.0 : -1.0, 0.0, 0.0};
    }

    return GeometryVector{0.0, dy >= 0.0 ? 1.0 : -1.0, 0.0};
}

int RectangleGeometry::GetDimension() const {
    return 2;
}

double RectangleGeometry::GetCenterX() const {
    return cx_;
}

double RectangleGeometry::GetCenterY() const {
    return cy_;
}

double RectangleGeometry::GetSizeX() const {
    return 2.0 * hx_;
}

double RectangleGeometry::GetSizeY() const {
    return 2.0 * hy_;
}