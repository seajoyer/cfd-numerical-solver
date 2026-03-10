#ifndef GEOMETRYPRIMITIVE_HPP
#define GEOMETRYPRIMITIVE_HPP

#include <cstdint>

/**
 * @struct GeometryPoint
 * @brief Cartesian point in 3D space.
 *
 * For 1D and 2D use cases, unused coordinates remain valid and may be zero.
 */
struct GeometryPoint final {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

/**
 * @struct GeometryVector
 * @brief Cartesian vector in 3D space.
 */
struct GeometryVector final {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

/**
 * @class GeometryPrimitive
 * @brief Abstract base class for embedded geometric primitives on Cartesian mesh.
 *
 * A primitive is defined in global physical coordinates and is responsible for:
 *  - point-in-solid query
 *  - signed-distance evaluation
 *  - outward-normal evaluation on or near the boundary
 *
 * Sign convention for signed distance:
 *  - negative: point is inside solid
 *  - zero:     point lies on boundary
 *  - positive: point is outside solid
 *
 * Normal convention:
 *  - outward unit normal from solid to fluid/outside region
 *
 * Implementations must be allocation-free in hot-path queries.
 */
class GeometryPrimitive {
public:
    virtual ~GeometryPrimitive() = default;

    /**
     * @brief Check whether a point lies inside the solid primitive.
     * @param p Query point in global physical coordinates.
     * @return true if point is inside solid.
     */
    [[nodiscard]] virtual bool Contains(const GeometryPoint& p) const = 0;

    /**
     * @brief Signed distance from point to primitive boundary.
     * @param p Query point in global physical coordinates.
     * @return Signed distance: negative inside, positive outside.
     */
    [[nodiscard]] virtual double SignedDistance(const GeometryPoint& p) const = 0;

    /**
     * @brief Outward unit normal at or near the primitive boundary.
     * @param p Query point in global physical coordinates.
     * @return Outward normal vector. Implementations should return a best-effort
     *         normalized vector even slightly away from the boundary.
     */
    [[nodiscard]] virtual GeometryVector OutwardNormal(const GeometryPoint& p) const = 0;

    /**
     * @brief Spatial dimension supported by this primitive.
     * @return 2 for planar primitive, 3 for volumetric primitive.
     */
    [[nodiscard]] virtual int GetDimension() const = 0;
};

#endif  // GEOMETRYPRIMITIVE_HPP