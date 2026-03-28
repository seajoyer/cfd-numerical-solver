#ifndef GMSHMESHBUILDER_HPP
#define GMSHMESHBUILDER_HPP

#include <cstddef>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

#include "geometry/Mesh.hpp"

/**
 * @brief Builder that converts a Gmsh mesh into internal generic mesh format.
 *
 * Supported cell element types:
 *  - 2D: triangle (type 2), quadrangle (type 3)
 *  - 3D: tetrahedron (type 4), hexahedron (type 5),
 *        prism (type 6), pyramid (type 7)
 *
 * Supported boundary element types:
 *  - 2D boundaries: line (type 1)
 *  - 3D boundaries: triangle (type 2), quadrangle (type 3)
 *
 * Boundary tags are taken from Gmsh physical groups of dimension dim - 1.
 */
class GmshMeshBuilder final {
public:
    /**
     * @brief Build mesh from an already prepared current Gmsh model.
     * @param dim Spatial dimension of target mesh (2 or 3).
     */
    [[nodiscard]] static Mesh BuildFromCurrentModel(int dim);

    /**
     * @brief Open mesh file in Gmsh and build internal mesh from it.
     * @param file_path Path to .msh file or another Gmsh-readable mesh file.
     * @param dim Spatial dimension of target mesh (2 or 3).
     */
    [[nodiscard]] static Mesh BuildFromFile(const std::string& file_path, int dim);

    /**
     * @brief Open .geo file in Gmsh, generate mesh, and build internal mesh.
     * @param file_path Path to .geo file.
     * @param dim Spatial dimension of target mesh (2 or 3).
     */
    [[nodiscard]] static Mesh BuildFromGeoFile(const std::string& file_path, int dim);

private:
    /**
     * @brief Small 3D vector utility.
     */
    struct Vec3 final {
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
    };

    /**
     * @brief Sorted face-node key used for face deduplication.
     */
    struct FaceKey final {
        std::vector<std::size_t> node_ids_sorted;

        [[nodiscard]] bool operator<(const FaceKey& other) const;
    };

    /**
     * @brief Boundary metadata attached to one boundary face key.
     */
    struct BoundaryFaceRecord final {
        int boundary_tag = -1;
    };

    /**
     * @brief RAII helper for local Gmsh session ownership.
     */
    class GmshSessionGuard final {
    public:
        explicit GmshSessionGuard(bool finalize_on_destroy);
        ~GmshSessionGuard();

        GmshSessionGuard(const GmshSessionGuard&) = delete;
        auto operator=(const GmshSessionGuard&) -> GmshSessionGuard& = delete;

    private:
        bool finalize_on_destroy_ = false;
    };

    static void ValidateInputDimension(int dim);
    static void EnsureGmshInitialized();

    [[nodiscard]] static auto IsSupportedCellElementType(int dim, int element_type) -> bool;
    [[nodiscard]] static auto IsSupportedBoundaryElementType(int dim, int element_type) -> bool;
    [[nodiscard]] static auto GetExpectedNodeCountForElementType(int element_type) -> int;

    [[nodiscard]] static std::vector<std::vector<int>> GetLocalFacesForElementType(int element_type);
    [[nodiscard]] static std::vector<std::size_t> ConvertElementConnectivity(
        const std::vector<std::size_t>& element_nodes,
        const std::vector<int>& local_face_pattern
    );

    [[nodiscard]] static FaceKey MakeFaceKey(const std::vector<std::size_t>& node_ids);

    [[nodiscard]] static Vec3 Add(const Vec3& a, const Vec3& b);
    [[nodiscard]] static Vec3 Subtract(const Vec3& a, const Vec3& b);
    [[nodiscard]] static Vec3 Multiply(double scalar, const Vec3& v);
    [[nodiscard]] static Vec3 Divide(const Vec3& v, double scalar);
    [[nodiscard]] static double Dot(const Vec3& a, const Vec3& b);
    [[nodiscard]] static Vec3 Cross(const Vec3& a, const Vec3& b);
    [[nodiscard]] static double Norm(const Vec3& v);
    [[nodiscard]] static Vec3 Normalize(const Vec3& v);
    [[nodiscard]] static Vec3 ToVec3(const Node& node);

    [[nodiscard]] static auto ComputePolygonCenter(
        const std::vector<std::size_t>& node_ids,
        const std::vector<Node>& nodes
    ) -> Vec3;

    [[nodiscard]] static double ComputeEdgeLength(
        const std::vector<std::size_t>& node_ids,
        const std::vector<Node>& nodes
    );

    [[nodiscard]] static double ComputePolygonArea3D(
        const std::vector<std::size_t>& node_ids,
        const std::vector<Node>& nodes
    );

    [[nodiscard]] static double ComputeCellArea2D(
        const std::vector<std::size_t>& node_ids,
        const std::vector<Node>& nodes
    );

    [[nodiscard]] static double ComputeTetraVolume(
        const Vec3& a,
        const Vec3& b,
        const Vec3& c,
        const Vec3& d
    );

    [[nodiscard]] static double ComputeCellVolume3D(
        int element_type,
        const std::vector<std::size_t>& node_ids,
        const std::vector<Node>& nodes
    );

    [[nodiscard]] static std::map<FaceKey, BoundaryFaceRecord> BuildBoundaryFaceMap(
        int dim,
        const std::unordered_map<std::size_t, std::size_t>& gmsh_to_internal_node
    );

    static void BuildNodes(
        Mesh& mesh,
        std::unordered_map<std::size_t, std::size_t>& gmsh_to_internal_node
    );

    static void BuildCellsAndFaces(
        Mesh& mesh,
        int dim,
        const std::unordered_map<std::size_t, std::size_t>& gmsh_to_internal_node,
        const std::map<FaceKey, BoundaryFaceRecord>& boundary_faces
    );

    static void FinalizeFaceNormals(Mesh& mesh);
};

#endif  // GMSHMESHBUILDER_HPP
