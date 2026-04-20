#include "solver/TimeStepCalculator.hpp"

#include <cmath>
#include <limits>
#include <vector>

#include "data/DataLayer.hpp"
#include "data/Variables.hpp"
#include "geometry/Cell.hpp"
#include "geometry/Face.hpp"
#include "geometry/Mesh.hpp"

double TimeStepCalculator::ComputeDt(const DataLayer& layer,
                                     const Mesh& mesh,
                                     const double gamma,
                                     const double cfl) {
    if (cfl <= 0.0) {
        return 0.0;
    }

    if (mesh.GetOwnedCellCount() == 0 || mesh.GetFaceCount() == 0) {
        return 0.0;
    }

    const auto& U = layer.U();

    std::vector<PrimitiveCell> primitive_by_cell(mesh.GetCellCount());

    for (std::size_t cell_id = 0; cell_id < mesh.GetCellCount(); ++cell_id) {
        ConservativeCell U_cell;
        U_cell.rho = U(cell_id, DataLayer::k_rho);
        U_cell.rhoU = U(cell_id, DataLayer::k_rhoU);
        U_cell.rhoV = U(cell_id, DataLayer::k_rhoV);
        U_cell.rhoW = U(cell_id, DataLayer::k_rhoW);
        U_cell.E = U(cell_id, DataLayer::k_E);

        primitive_by_cell[cell_id] = PrimitiveFromConservativeCell(U_cell, gamma);
    }

    double dt_min = std::numeric_limits<double>::infinity();
    bool has_valid_dt = false;

    for (std::size_t cell_id = 0; cell_id < mesh.GetOwnedCellCount(); ++cell_id) {
        const Cell& cell = mesh.GetCell(cell_id);

        double spectral_sum = 0.0;

        for (const std::size_t face_id : cell.face_ids) {
            const Face& face = mesh.GetFace(face_id);

            FaceNormal normal;
            normal.x = face.normal_x;
            normal.y = face.normal_y;
            normal.z = face.normal_z;

            if (face.owner_cell_id != cell.local_id) {
                normal.x = -normal.x;
                normal.y = -normal.y;
                normal.z = -normal.z;
            }

            const PrimitiveCell& w = primitive_by_cell[cell.local_id];
            const double un = NormalVelocity(w, normal);
            const double c = SoundSpeed(w, gamma);

            spectral_sum += (std::abs(un) + c) * face.measure;
        }

        if (spectral_sum <= 0.0) {
            continue;
        }

        const double dt_cell = cfl * cell.volume / spectral_sum;

        if (std::isfinite(dt_cell) && dt_cell > 0.0) {
            dt_min = std::min(dt_min, dt_cell);
            has_valid_dt = true;
        }
    }

    if (!has_valid_dt || !std::isfinite(dt_min) || dt_min <= 0.0) {
        return 0.0;
    }

    return dt_min;
}
