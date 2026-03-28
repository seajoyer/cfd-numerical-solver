#ifndef STRUCTUREDMESHBUILDER_HPP
#define STRUCTUREDMESHBUILDER_HPP

#include "geometry/Mesh.hpp"

/**
 * @brief Builder of uniform Cartesian meshes in generic unstructured format.
 *
 * Resulting mesh is stored as:
 *  - nodes
 *  - faces
 *  - cells
 *
 * Boundary tags:
 *  - 0 : xmin
 *  - 1 : xmax
 *  - 2 : ymin
 *  - 3 : ymax
 *  - 4 : zmin
 *  - 5 : zmax
 */
class StructuredMeshBuilder final {
public:
    static constexpr int k_xmin_tag = 0;
    static constexpr int k_xmax_tag = 1;
    static constexpr int k_ymin_tag = 2;
    static constexpr int k_ymax_tag = 3;
    static constexpr int k_zmin_tag = 4;
    static constexpr int k_zmax_tag = 5;

    /**
     * @brief Build uniform Cartesian 2D mesh in generic mesh format.
     * @param nx Number of cells in x.
     * @param ny Number of cells in y.
     * @param x_min Domain minimum x.
     * @param x_max Domain maximum x.
     * @param y_min Domain minimum y.
     * @param y_max Domain maximum y.
     */
    [[nodiscard]] static Mesh BuildUniformCartesian2D(int nx, int ny,
                                                      double x_min, double x_max,
                                                      double y_min, double y_max);

    /**
     * @brief Build uniform Cartesian 3D mesh in generic mesh format.
     * @param nx Number of cells in x.
     * @param ny Number of cells in y.
     * @param nz Number of cells in z.
     * @param x_min Domain minimum x.
     * @param x_max Domain maximum x.
     * @param y_min Domain minimum y.
     * @param y_max Domain maximum y.
     * @param z_min Domain minimum z.
     * @param z_max Domain maximum z.
     */
    [[nodiscard]] static Mesh BuildUniformCartesian3D(int nx, int ny, int nz,
                                                      double x_min, double x_max,
                                                      double y_min, double y_max,
                                                      double z_min, double z_max);

private:
    static void Validate2DInput(int nx, int ny,
                                double x_min, double x_max,
                                double y_min, double y_max);

    static void Validate3DInput(int nx, int ny, int nz,
                                double x_min, double x_max,
                                double y_min, double y_max,
                                double z_min, double z_max);
};

#endif  // STRUCTUREDMESHBUILDER_HPP
