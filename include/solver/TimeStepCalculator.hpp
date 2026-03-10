#ifndef TIMESTEPCALCULATOR_HPP
#define TIMESTEPCALCULATOR_HPP

class DataLayer;
class Mesh;

/**
 * @class TimeStepCalculator
 * @brief CFL-based timestep selection from conservative state U on structured mesh.
 *
 * Uses only conservative U(var,i,j,k) with var=(rho, rhoU, rhoV, rhoW, E).
 * Mesh metrics are taken from Mesh::Dx/Dy/Dz.
 *
 * For each fluid core cell and each active axis:
 *   s_axis = |v_axis| + c,  c = sqrt(gamma * P / rho),  P = (gamma-1)*(E - 0.5*rho*|v|^2)
 * Local dt candidates:
 *   dt_x = cfl * dx(i) / s_x
 *   dt_y = cfl * dy(j) / s_y
 *   dt_z = cfl * dz(k) / s_z
 * Returns global minimum over fluid core cells and active axes.
 */
class TimeStepCalculator final {
public:
    /**
     * @brief Computes stable explicit timestep dt.
     * @param layer DataLayer with conservative state U.
     * @param mesh Structured mesh with geometry, metrics, and cell types.
     * @param gamma Ratio of specific heats.
     * @param cfl CFL number (0 < cfl <= 1).
     * @return dt > 0 on success; 0.0 if dt cannot be computed.
     */
    static auto ComputeDt(const DataLayer& layer, const Mesh& mesh, double gamma, double cfl) -> double;
};

#endif  // TIMESTEPCALCULATOR_HPP