#ifndef RECTANGLEGEOMETRY_HPP
#define RECTANGLEGEOMETRY_HPP

#include "data/geometry/GeometryPrimitive.hpp"

/**
 * @class RectangleGeometry
 * @brief Axis-aligned 2D rectangle primitive embedded in the XY plane.
 *
 * The rectangle is defined by center (cx, cy) and full sizes (size_x, size_y).
 * Coordinate z is ignored in all geometric queries.
 *
 * Sign convention for signed distance:
 *  - negative: point is inside rectangle
 *  - zero:     point is on boundary
 *  - positive: point is outside rectangle
 *
 * Normal convention:
 *  - outward unit normal from rectangle interior to exterior
 */
class RectangleGeometry final : public GeometryPrimitive {
public:
    /**
     * @brief Construct rectangle from center and full side lengths.
     * @param cx Rectangle center x-coordinate.
     * @param cy Rectangle center y-coordinate.
     * @param size_x Full rectangle size along x, must be > 0.
     * @param size_y Full rectangle size along y, must be > 0.
     */
    RectangleGeometry(double cx, double cy, double size_x, double size_y);

    [[nodiscard]] bool Contains(const GeometryPoint& p) const override;
    [[nodiscard]] double SignedDistance(const GeometryPoint& p) const override;
    [[nodiscard]] GeometryVector OutwardNormal(const GeometryPoint& p) const override;
    [[nodiscard]] int GetDimension() const override;

    [[nodiscard]] double GetCenterX() const;
    [[nodiscard]] double GetCenterY() const;
    [[nodiscard]] double GetSizeX() const;
    [[nodiscard]] double GetSizeY() const;

private:
    double cx_ = 0.0;
    double cy_ = 0.0;
    double hx_ = 0.5;
    double hy_ = 0.5;
};

#endif  // RECTANGLEGEOMETRY_HPP