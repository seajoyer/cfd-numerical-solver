#include "reconstruction/P1Reconstruction.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "data/Workspace.hpp"
#include "geometry/Cell.hpp"
#include "geometry/Face.hpp"
#include "geometry/Mesh.hpp"

void P1Reconstruction::SetLimiter(const LimiterType type) {
    limiter_type_ = type;
}

PrimitiveCell P1Reconstruction::LoadCellPrimitive(const Workspace& workspace,
                                                  const std::size_t cell_id) const {
    const auto& W = workspace.W();

    PrimitiveCell state;
    state.rho = W(cell_id, Workspace::k_rho);
    state.u = W(cell_id, Workspace::k_u);
    state.v = W(cell_id, Workspace::k_v);
    state.w = W(cell_id, Workspace::k_w);
    state.P = W(cell_id, Workspace::k_p);

    return state;
}

void P1Reconstruction::CollectNeighborCellIds(const Mesh& mesh,
                                              const Cell& cell,
                                              std::vector<std::size_t>& neighbor_cell_ids) const {
    neighbor_cell_ids.clear();
    neighbor_cell_ids.reserve(cell.face_ids.size());

    for (const std::size_t face_id : cell.face_ids) {
        const Face& face = mesh.GetFace(face_id);

        if (face.IsPhysicalBoundary()) {
            continue;
        }

        const std::size_t neighbor_id =
            (face.owner_cell_id == cell.local_id) ? face.neighbor_cell_id : face.owner_cell_id;

        if (neighbor_id == Face::k_invalid_cell_id || neighbor_id == cell.local_id) {
            continue;
        }

        neighbor_cell_ids.push_back(neighbor_id);
    }
}

double P1Reconstruction::SolveLeastSquaresComponent(const double a11,
                                                    const double a12,
                                                    const double a13,
                                                    const double a22,
                                                    const double a23,
                                                    const double a33,
                                                    const double b1,
                                                    const double b2,
                                                    const double b3,
                                                    const int dim,
                                                    double& gx,
                                                    double& gy,
                                                    double& gz) const {
    gx = 0.0;
    gy = 0.0;
    gz = 0.0;

    constexpr double eps = 1e-14;

    if (dim == 1) {
        if (std::abs(a11) <= eps) {
            return 0.0;
        }
        gx = b1 / a11;
        return 1.0;
    }

    if (dim == 2) {
        const double det = a11 * a22 - a12 * a12;
        if (std::abs(det) <= eps) {
            return 0.0;
        }

        gx = (b1 * a22 - b2 * a12) / det;
        gy = (-b1 * a12 + b2 * a11) / det;
        return 1.0;
    }

    const double det =
        a11 * (a22 * a33 - a23 * a23)
        - a12 * (a12 * a33 - a13 * a23)
        + a13 * (a12 * a23 - a13 * a22);

    if (std::abs(det) <= eps) {
        return 0.0;
    }

    gx =
    (b1 * (a22 * a33 - a23 * a23)
        - a12 * (b2 * a33 - a23 * b3)
        + a13 * (b2 * a23 - a22 * b3)) / det;

    gy =
    (a11 * (b2 * a33 - a23 * b3)
        - b1 * (a12 * a33 - a13 * a23)
        + a13 * (a12 * b3 - b2 * a13)) / det;

    gz =
    (a11 * (a22 * b3 - b2 * a23)
        - a12 * (a12 * b3 - b2 * a13)
        + b1 * (a12 * a23 - a22 * a13)) / det;

    return 1.0;
}

P1Reconstruction::PrimitiveGradient P1Reconstruction::ComputeGradient(const Mesh& mesh,
                                                                      const Workspace& workspace,
                                                                      const Cell& cell) const {
    PrimitiveGradient gradient{};

    std::vector<std::size_t> neighbor_ids;
    CollectNeighborCellIds(mesh, cell, neighbor_ids);

    if (neighbor_ids.empty()) {
        return gradient;
    }

    const PrimitiveCell wc = LoadCellPrimitive(workspace, cell.local_id);

    double a11 = 0.0;
    double a12 = 0.0;
    double a13 = 0.0;
    double a22 = 0.0;
    double a23 = 0.0;
    double a33 = 0.0;

    double b1_rho = 0.0, b2_rho = 0.0, b3_rho = 0.0;
    double b1_u = 0.0, b2_u = 0.0, b3_u = 0.0;
    double b1_v = 0.0, b2_v = 0.0, b3_v = 0.0;
    double b1_w = 0.0, b2_w = 0.0, b3_w = 0.0;
    double b1_P = 0.0, b2_P = 0.0, b3_P = 0.0;

    for (const std::size_t neighbor_id : neighbor_ids) {
        const Cell& neighbor = mesh.GetCell(neighbor_id);
        const PrimitiveCell wn = LoadCellPrimitive(workspace, neighbor.local_id);

        const double dx = neighbor.center_x - cell.center_x;
        const double dy = neighbor.center_y - cell.center_y;
        const double dz = neighbor.center_z - cell.center_z;

        a11 += dx * dx;
        a12 += dx * dy;
        a13 += dx * dz;
        a22 += dy * dy;
        a23 += dy * dz;
        a33 += dz * dz;

        const double drho = wn.rho - wc.rho;
        const double du = wn.u - wc.u;
        const double dv = wn.v - wc.v;
        const double dw = wn.w - wc.w;
        const double dP = wn.P - wc.P;

        b1_rho += dx * drho;
        b2_rho += dy * drho;
        b3_rho += dz * drho;

        b1_u += dx * du;
        b2_u += dy * du;
        b3_u += dz * du;

        b1_v += dx * dv;
        b2_v += dy * dv;
        b3_v += dz * dv;

        b1_w += dx * dw;
        b2_w += dy * dw;
        b3_w += dz * dw;

        b1_P += dx * dP;
        b2_P += dy * dP;
        b3_P += dz * dP;
    }

    const int dim = mesh.GetDim();

    SolveLeastSquaresComponent(a11, a12, a13, a22, a23, a33,
                               b1_rho, b2_rho, b3_rho, dim,
                               gradient.dx.rho, gradient.dy.rho, gradient.dz.rho);

    SolveLeastSquaresComponent(a11, a12, a13, a22, a23, a33,
                               b1_u, b2_u, b3_u, dim,
                               gradient.dx.u, gradient.dy.u, gradient.dz.u);

    SolveLeastSquaresComponent(a11, a12, a13, a22, a23, a33,
                               b1_v, b2_v, b3_v, dim,
                               gradient.dx.v, gradient.dy.v, gradient.dz.v);

    SolveLeastSquaresComponent(a11, a12, a13, a22, a23, a33,
                               b1_w, b2_w, b3_w, dim,
                               gradient.dx.w, gradient.dy.w, gradient.dz.w);

    SolveLeastSquaresComponent(a11, a12, a13, a22, a23, a33,
                               b1_P, b2_P, b3_P, dim,
                               gradient.dx.P, gradient.dy.P, gradient.dz.P);

    return gradient;
}

double P1Reconstruction::ComputeBarthJespersenPhi(const double w_cell,
                                                  const double w_min,
                                                  const double w_max,
                                                  const double w_face_candidate) const {
    constexpr double eps = 1e-14;

    const double delta = w_face_candidate - w_cell;

    if (delta > eps) {
        return std::min(1.0, (w_max - w_cell) / delta);
    }

    if (delta < -eps) {
        return std::min(1.0, (w_min - w_cell) / delta);
    }

    return 1.0;
}

double P1Reconstruction::ComputeLimiter(const Mesh& mesh,
                                        const Workspace& workspace,
                                        const Cell& cell,
                                        const PrimitiveGradient& gradient,
                                        const double face_center_x,
                                        const double face_center_y,
                                        const double face_center_z) const {
    if (limiter_type_ == LimiterType::kNone) {
        return 1.0;
    }

    const PrimitiveCell wc = LoadCellPrimitive(workspace, cell.local_id);

    PrimitiveCell w_min = wc;
    PrimitiveCell w_max = wc;

    std::vector<std::size_t> neighbor_ids;
    CollectNeighborCellIds(mesh, cell, neighbor_ids);

    for (const std::size_t neighbor_id : neighbor_ids) {
        const PrimitiveCell wn = LoadCellPrimitive(workspace, neighbor_id);

        w_min.rho = std::min(w_min.rho, wn.rho);
        w_min.u = std::min(w_min.u, wn.u);
        w_min.v = std::min(w_min.v, wn.v);
        w_min.w = std::min(w_min.w, wn.w);
        w_min.P = std::min(w_min.P, wn.P);

        w_max.rho = std::max(w_max.rho, wn.rho);
        w_max.u = std::max(w_max.u, wn.u);
        w_max.v = std::max(w_max.v, wn.v);
        w_max.w = std::max(w_max.w, wn.w);
        w_max.P = std::max(w_max.P, wn.P);
    }

    const double dx = face_center_x - cell.center_x;
    const double dy = face_center_y - cell.center_y;
    const double dz = face_center_z - cell.center_z;

    const PrimitiveCell w_face_candidate{
        wc.rho + gradient.dx.rho * dx + gradient.dy.rho * dy + gradient.dz.rho * dz,
        wc.u + gradient.dx.u * dx + gradient.dy.u * dy + gradient.dz.u * dz,
        wc.v + gradient.dx.v * dx + gradient.dy.v * dy + gradient.dz.v * dz,
        wc.w + gradient.dx.w * dx + gradient.dy.w * dy + gradient.dz.w * dz,
        wc.P + gradient.dx.P * dx + gradient.dy.P * dy + gradient.dz.P * dz
    };

    double phi = 1.0;

    phi = std::min(phi, ComputeBarthJespersenPhi(wc.rho, w_min.rho, w_max.rho, w_face_candidate.rho));
    phi = std::min(phi, ComputeBarthJespersenPhi(wc.u, w_min.u, w_max.u, w_face_candidate.u));
    phi = std::min(phi, ComputeBarthJespersenPhi(wc.v, w_min.v, w_max.v, w_face_candidate.v));
    phi = std::min(phi, ComputeBarthJespersenPhi(wc.w, w_min.w, w_max.w, w_face_candidate.w));
    phi = std::min(phi, ComputeBarthJespersenPhi(wc.P, w_min.P, w_max.P, w_face_candidate.P));

    return std::clamp(phi, 0.0, 1.0);
}

PrimitiveCell P1Reconstruction::ExtrapolateToPoint(const Cell& cell,
                                                   const PrimitiveCell& cell_state,
                                                   const PrimitiveGradient& gradient,
                                                   const double limiter,
                                                   const double x,
                                                   const double y,
                                                   const double z) const {
    const double dx = x - cell.center_x;
    const double dy = y - cell.center_y;
    const double dz = z - cell.center_z;

    PrimitiveCell state;
    state.rho = cell_state.rho + limiter * (gradient.dx.rho * dx + gradient.dy.rho * dy + gradient.dz.rho * dz);
    state.u = cell_state.u + limiter * (gradient.dx.u * dx + gradient.dy.u * dy + gradient.dz.u * dz);
    state.v = cell_state.v + limiter * (gradient.dx.v * dx + gradient.dy.v * dy + gradient.dz.v * dz);
    state.w = cell_state.w + limiter * (gradient.dx.w * dx + gradient.dy.w * dy + gradient.dz.w * dz);
    state.P = cell_state.P + limiter * (gradient.dx.P * dx + gradient.dy.P * dy + gradient.dz.P * dz);

    return state;
}

void P1Reconstruction::ReconstructInteriorFace(const Mesh& mesh,
                                               const Workspace& workspace,
                                               const Face& face,
                                               PrimitiveCell& owner_state,
                                               PrimitiveCell& neighbor_state) const {
    if (!(face.IsInternal() || face.IsMPIBoundary())) {
        throw std::runtime_error(
            "P1Reconstruction::ReconstructInteriorFace: face is not internal or MPI boundary"
        );
    }

    const Cell& owner = mesh.GetCell(face.owner_cell_id);
    const Cell& neighbor = mesh.GetCell(face.neighbor_cell_id);

    const PrimitiveCell owner_cell_state = LoadCellPrimitive(workspace, owner.local_id);
    const PrimitiveCell neighbor_cell_state = LoadCellPrimitive(workspace, neighbor.local_id);

    const PrimitiveGradient owner_gradient = ComputeGradient(mesh, workspace, owner);
    const PrimitiveGradient neighbor_gradient = ComputeGradient(mesh, workspace, neighbor);

    const double owner_limiter =
        ComputeLimiter(mesh, workspace, owner, owner_gradient,
                       face.center_x, face.center_y, face.center_z);

    const double neighbor_limiter =
        ComputeLimiter(mesh, workspace, neighbor, neighbor_gradient,
                       face.center_x, face.center_y, face.center_z);

    owner_state = ExtrapolateToPoint(owner,
                                     owner_cell_state,
                                     owner_gradient,
                                     owner_limiter,
                                     face.center_x,
                                     face.center_y,
                                     face.center_z);

    neighbor_state = ExtrapolateToPoint(neighbor,
                                        neighbor_cell_state,
                                        neighbor_gradient,
                                        neighbor_limiter,
                                        face.center_x,
                                        face.center_y,
                                        face.center_z);
}

void P1Reconstruction::ReconstructBoundaryFaceInterior(const Mesh& mesh,
                                                       const Workspace& workspace,
                                                       const Face& face,
                                                       PrimitiveCell& interior_state) const {
    if (!face.IsPhysicalBoundary()) {
        throw std::runtime_error(
            "P1Reconstruction::ReconstructBoundaryFaceInterior: face is not physical boundary"
        );
    }

    const Cell& owner = mesh.GetCell(face.owner_cell_id);
    const PrimitiveCell owner_cell_state = LoadCellPrimitive(workspace, owner.local_id);
    const PrimitiveGradient owner_gradient = ComputeGradient(mesh, workspace, owner);

    const double owner_limiter =
        ComputeLimiter(mesh, workspace, owner, owner_gradient,
                       face.center_x, face.center_y, face.center_z);

    interior_state = ExtrapolateToPoint(owner,
                                        owner_cell_state,
                                        owner_gradient,
                                        owner_limiter,
                                        face.center_x,
                                        face.center_y,
                                        face.center_z);
}
