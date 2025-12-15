#ifndef VNRARTIFICIALVISCOSITY_HPP
#define VNRARTIFICIALVISCOSITY_HPP

#include "viscosity/ArtificialViscosity.hpp"

/**
 * @class VNRArtificialViscosity
 * @brief von Neumann–Richtmyer artificial viscosity model.
 *
 * Implements a classic compressive artificial viscosity of the form:
 *
 *   q_{j+1/2} = 0                      if Δu >= 0  (no compression),
 *              C2 * rho_bar * (Δu)^2 +
 *              C1 * rho_bar * c_bar * |Δu|   otherwise,
 *
 * where:
 *   Δu      = u_{j+1} - u_j,
 *   rho_bar = 0.5 * (rho_j + rho_{j+1}),
 *   c_bar   = 0.5 * (c_j + c_{j+1}),  c = sound speed,
 *
 * and C1, C2 are user-tunable coefficients.
 *
 * The resulting q_{j+1/2} is meant to be added to the pressure on both
 * sides of the interface prior to the Riemann solve.
 */
class VNRArtificialViscosity : public ArtificialViscosity {
public:
    /**
     * @brief Constructs the VNR viscosity model with given coefficients.
     *
     * @param C1 Linear term coefficient (typical 0.1–0.5).
     * @param C2 Quadratic term coefficient (typical 1.0–2.0).
     */
    explicit VNRArtificialViscosity(const Settings& settings, double C1 = 1.5,
                                    double C2 = 6.0)
        : settings_(settings), C1_(C1), C2_(C2) {
    }

    void ComputeInterfaceQ(const DataLayer& layer,
                           double dx,
                           xt::xarray<double>& q) const override;

private:
    double C1_; ///< Linear coefficient in q
    double C2_; ///< Quadratic coefficient in q
    Settings settings_;
};

#endif  // VNRARTIFICIALVISCOSITY_HPP