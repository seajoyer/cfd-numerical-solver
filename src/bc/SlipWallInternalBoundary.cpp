#include "bc/SlipWallInternalBoundary.hpp"

#include <cmath>

void SlipWallInternalBoundary::BuildBoundaryState(const PrimitiveCell& fluid_state,
                                                  const ImmersedFaceInfo& face_info,
                                                  PrimitiveCell& boundary_state) const {
    if (!face_info.is_active) {
        boundary_state = fluid_state;
        return;
    }

    const double nx_in = face_info.normal_x;
    const double ny_in = face_info.normal_y;
    const double nz_in = face_info.normal_z;

    const double norm2 = nx_in * nx_in + ny_in * ny_in + nz_in * nz_in;
    if (norm2 <= k_eps * k_eps) {
        boundary_state = fluid_state;
        return;
    }

    const double norm = std::sqrt(norm2);
    const double nx = nx_in / norm;
    const double ny = ny_in / norm;
    const double nz = nz_in / norm;

    const double un = fluid_state.u * nx + fluid_state.v * ny + fluid_state.w * nz;

    const double ut_x = fluid_state.u - un * nx;
    const double ut_y = fluid_state.v - un * ny;
    const double ut_z = fluid_state.w - un * nz;

    boundary_state.rho = fluid_state.rho;
    boundary_state.P = fluid_state.P;

    boundary_state.u = ut_x - un * nx;
    boundary_state.v = ut_y - un * ny;
    boundary_state.w = ut_z - un * nz;
}
