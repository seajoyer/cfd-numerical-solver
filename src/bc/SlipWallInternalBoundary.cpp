#include "bc/SlipWallInternalBoundary.hpp"

#include <cmath>

void SlipWallInternalBoundary::BuildBoundaryState(const PrimitiveCell& fluid_state,
                                                  const ImmersedFaceInfo& face_info,
                                                  PrimitiveCell& boundary_state) const {
    const double nx_in = face_info.normal_x;
    const double ny_in = face_info.normal_y;
    const double nz_in = face_info.normal_z;

    const double norm = std::sqrt(nx_in * nx_in + ny_in * ny_in + nz_in * nz_in);
    if (!(norm > 0.0)) {
        boundary_state = fluid_state;
        return;
    }

    const double nx = nx_in / norm;
    const double ny = ny_in / norm;
    const double nz = nz_in / norm;

    const double un = fluid_state.u * nx + fluid_state.v * ny + fluid_state.w * nz;

    boundary_state.rho = fluid_state.rho;
    boundary_state.P = fluid_state.P;

    boundary_state.u = fluid_state.u - 2.0 * un * nx;
    boundary_state.v = fluid_state.v - 2.0 * un * ny;
    boundary_state.w = fluid_state.w - 2.0 * un * nz;
}