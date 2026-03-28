#ifndef P0RECONSTRUCTION_HPP
#define P0RECONSTRUCTION_HPP

#include "reconstruction/Reconstruction.hpp"

/**
 * @class P0Reconstruction
 * @brief Piecewise-constant face reconstruction on generic meshes.
 *
 * For internal face:
 * - owner_state    = primitive state in owner cell
 * - neighbor_state = primitive state in neighbor cell
 *
 * For boundary face:
 * - interior_state = primitive state in owner cell
 */
class P0Reconstruction final : public Reconstruction {
public:
    P0Reconstruction() = default;
    ~P0Reconstruction() override = default;

    void ReconstructInteriorFace(const Mesh& mesh,
                                 const Workspace& workspace,
                                 const Face& face,
                                 PrimitiveCell& owner_state,
                                 PrimitiveCell& neighbor_state) const override;

    void ReconstructBoundaryFaceInterior(const Mesh& mesh,
                                         const Workspace& workspace,
                                         const Face& face,
                                         PrimitiveCell& interior_state) const override;

private:
    [[nodiscard]] PrimitiveCell LoadCellPrimitive(const Workspace& workspace,
                                                  std::size_t cell_id) const;
};

#endif  // P0RECONSTRUCTION_HPP
