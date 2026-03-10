#ifndef CIRCLEGEOMETRY_HPP
#define CIRCLEGEOMETRY_HPP

#include "data/geometry/GeometryPrimitive.hpp"

/**
 * @class CircleGeometry
 * @brief 2D circle primitive embedded in the XY plane.
 *
 * The circle is defined by center (cx, cy) and radius r.
 * Coordinate z is ignored in all geometric queries.
 *
 * Sign convention for signed distance:
 *  - negative: point is inside circle
 *  - zero:     point is on boundary
 *  - positive: point is outside circle
 *
 * Normal convention:
 *  - outward unit normal from circle interior to exterior
 */
class CircleGeometry final : public GeometryPrimitive {
public:
    /**
     * @brief Construct circle from center and radius.
     * @param cx Circle center x-coordinate.
     * @param cy Circle center y-coordinate.
     * @param radius Circle radius, must be > 0.
     */
    CircleGeometry(double cx, double cy, double radius);

    [[nodiscard]] bool Contains(const GeometryPoint& p) const override;
    [[nodiscard]] double SignedDistance(const GeometryPoint& p) const override;
    [[nodiscard]] GeometryVector OutwardNormal(const GeometryPoint& p) const override;
    [[nodiscard]] int GetDimension() const override;

    [[nodiscard]] double GetCenterX() const;
    [[nodiscard]] double GetCenterY() const;
    [[nodiscard]] double GetRadius() const;

private:
    double cx_ = 0.0;
    double cy_ = 0.0;
    double radius_ = 1.0;
};

#endif  // CIRCLEGEOMETRY_HPP