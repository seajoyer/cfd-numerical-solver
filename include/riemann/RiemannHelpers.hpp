#ifndef RIEMANNHELPERS_HPP
#define RIEMANNHELPERS_HPP

#include <cstdint>

#include "data/Variables.hpp"

/**
 * @file RiemannHelpers.hpp
 * @brief Helper utilities for axis-aligned Riemann solvers (no allocations).
 *
 * This module contains small utilities that depend on the interface normal axis,
 * such as splitting velocity into normal/tangential components and composing
 * momentum components back into global (x,y,z) ordering.
 *
 * Add further Riemann-solver-specific helpers here later (Roe averages, entropy fixes, etc.).
 */

namespace riemann {
    /**
     * @brief Split velocity into normal and two tangential components for a face normal along axis.
     *
     * Mapping:
     *  - Axis::X: (un, ut1, ut2) = (u, v, w)
     *  - Axis::Y: (un, ut1, ut2) = (v, u, w)
     *  - Axis::Z: (un, ut1, ut2) = (w, u, v)
     */
    inline void SplitVelocity(const PrimitiveCell& w, const Axis axis, double& un, double& ut1, double& ut2) {
        if (axis == Axis::X) {
            un = w.u;
            ut1 = w.v;
            ut2 = w.w;
            return;
        }
        if (axis == Axis::Y) {
            un = w.v;
            ut1 = w.u;
            ut2 = w.w;
            return;
        }
        // Axis::Z
        un = w.w;
        ut1 = w.u;
        ut2 = w.v;
    }

    /**
     * @brief Compose conservative momentum components (rhoU,rhoV,rhoW) from (rho, un, ut1, ut2).
     *
     * Inverse mapping to SplitVelocity.
     */
    inline void ComposeMomentum(const double rho, const double un, const double ut1, const double ut2, const Axis axis,
                                double& rhoU, double& rhoV, double& rhoW) {
        if (axis == Axis::X) {
            rhoU = rho * un;
            rhoV = rho * ut1;
            rhoW = rho * ut2;
            return;
        }
        if (axis == Axis::Y) {
            rhoU = rho * ut1;
            rhoV = rho * un;
            rhoW = rho * ut2;
            return;
        }
        // Axis::Z
        rhoU = rho * ut1;
        rhoV = rho * ut2;
        rhoW = rho * un;
    }

    /**
     * @brief Compose (u,v,w) from (un, ut1, ut2) for a face normal along axis.
     *
     * Inverse mapping to SplitVelocity.
     */
    inline void ComposeVelocity(const double un, const double ut1, const double ut2, const Axis axis,
                                double& u, double& v, double& w) {
        if (axis == Axis::X) {
            u = un;
            v = ut1;
            w = ut2;
            return;
        }
        if (axis == Axis::Y) {
            u = ut1;
            v = un;
            w = ut2;
            return;
        }
        // Axis::Z
        u = ut1;
        v = ut2;
        w = un;
    }

    inline void ComposeVector(const double vn, const double vt1, const double vt2, const Axis axis,
                              double& vx, double& vy, double& vz) {
        if (axis == Axis::X) {
            vx = vn;
            vy = vt1;
            vz = vt2;
            return;
        }
        if (axis == Axis::Y) {
            vx = vt1;
            vy = vn;
            vz = vt2;
            return;
        }
        // Axis::Z
        vx = vt1;
        vy = vt2;
        vz = vn;
    }
} // namespace riemann

#endif  // RIEMANNHELPERS_HPP
