#ifndef TIMESTEPCALCULATOR_HPP
#define TIMESTEPCALCULATOR_HPP

#include "data/DataLayer.hpp"

/**
 * @class TimeStepCalculator
 * @brief CFL-based timestep selection from conservative state U (3D-by-default).
 *
 * Uses only conservative U(var,i,j,k) with var=(rho, rhoU, rhoV, rhoW, E).
 * Mesh metrics are taken from DataLayer::Dx/Dy/Dz (ready for nonuniform grids).
 *
 * For each core cell and each active axis:
 *   s_axis = |v_axis| + c,  c = sqrt(gamma * P / rho),  P = (gamma-1)*(E - 0.5*rho*|v|^2)
 * Local dt candidates:
 *   dt_x = cfl * dx(i) / s_x
 *   dt_y = cfl * dy(j) / s_y
 *   dt_z = cfl * dz(k) / s_z
 * Returns global minimum over core cells and active axes.
 */
class TimeStepCalculator final {
public:
    /**
     * @brief Computes stable explicit timestep dt.
     * @param layer DataLayer with conservative state U and grid metrics.
     * @param gamma Ratio of specific heats.
     * @param cfl CFL number (0 < cfl <= 1).
     * @return dt > 0 on success; 0.0 if dt cannot be computed.
     */
    static auto ComputeDt(const DataLayer& layer, double gamma, double cfl) -> double;
};

#endif  // TIMESTEPCALCULATOR_HPP
