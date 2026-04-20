#ifndef P1RECONSTRUCTION_HPP
#define P1RECONSTRUCTION_HPP

#include <cstddef>
#include <cstdint>
#include <vector>

#include "reconstruction/Reconstruction.hpp"

class Mesh;
class Workspace;
class Face;
struct Cell;

/**
 * @enum LimiterType
 * @brief Limiter type for generic-mesh piecewise-linear reconstruction.
 */
enum class LimiterType : std::uint8_t {
    kNone = 0,
    kBarthJespersen = 1
};

/**
 * @class P1Reconstruction
 * @brief Piecewise-linear reconstruction on generic face-based meshes.
 *
 * Reconstruction model:
 * - compute one cell gradient from neighboring cell-centered primitive values
 * - optionally limit reconstructed variation
 * - extrapolate from cell center to face center
 *
 * For internal or MPI face:
 * - reconstruct owner-side state at face center
 * - reconstruct neighbor-side state at face center
 *
 * For physical boundary face:
 * - reconstruct owner-side interior state at face center
 *
 * Notes:
 * - Supports local meshes with owned and ghost cells
 * - Uses local face connectivity
 * - Ghost cells may participate in gradients and stencils
 */
class P1Reconstruction final : public Reconstruction {
public:
    P1Reconstruction() = default;
    ~P1Reconstruction() override = default;

    void ReconstructInteriorFace(const Mesh& mesh,
                                 const Workspace& workspace,
                                 const Face& face,
                                 PrimitiveCell& owner_state,
                                 PrimitiveCell& neighbor_state) const override;

    void ReconstructBoundaryFaceInterior(const Mesh& mesh,
                                         const Workspace& workspace,
                                         const Face& face,
                                         PrimitiveCell& interior_state) const override;

    void SetLimiter(LimiterType type);

private:
    struct PrimitiveGradient final {
        PrimitiveCell dx;
        PrimitiveCell dy;
        PrimitiveCell dz;
    };

    LimiterType limiter_type_ = LimiterType::kBarthJespersen;

    [[nodiscard]] PrimitiveCell LoadCellPrimitive(const Workspace& workspace,
                                                  std::size_t cell_id) const;

    void CollectNeighborCellIds(const Mesh& mesh,
                                const Cell& cell,
                                std::vector<std::size_t>& neighbor_cell_ids) const;

    [[nodiscard]] PrimitiveGradient ComputeGradient(const Mesh& mesh,
                                                    const Workspace& workspace,
                                                    const Cell& cell) const;

    [[nodiscard]] double ComputeLimiter(const Mesh& mesh,
                                        const Workspace& workspace,
                                        const Cell& cell,
                                        const PrimitiveGradient& gradient,
                                        double face_center_x,
                                        double face_center_y,
                                        double face_center_z) const;

    [[nodiscard]] PrimitiveCell ExtrapolateToPoint(const Cell& cell,
                                                   const PrimitiveCell& cell_state,
                                                   const PrimitiveGradient& gradient,
                                                   double limiter,
                                                   double x,
                                                   double y,
                                                   double z) const;

    double SolveLeastSquaresComponent(double a11, double a12, double a13,
                                      double a22, double a23, double a33,
                                      double b1, double b2, double b3,
                                      int dim, double& gx, double& gy,
                                      double& gz) const;

    [[nodiscard]] double ComputeBarthJespersenPhi(double w_cell, double w_min, double w_max,
                                                  double w_face_candidate) const;
};

#endif  // P1RECONSTRUCTION_HPP
