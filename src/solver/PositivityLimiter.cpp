#include "solver/PositivityLimiter.hpp"

#include "solver/EOS.hpp"

void PositivityLimiter::Apply(Conservative& u, double gamma, double rho_min,
                              double p_min) {
    if (u.rho < rho_min) {
        double u_vel = 0.0;
        if (u.rho > 0.0) {
            u_vel = u.rhoU / u.rho;
        }
        u.rho = rho_min;
        u.rhoU = u.rho * u_vel;
    }

    double p = EOS::Pressure(u, gamma);

    if (p < p_min) {
        const double rho = u.rho;
        double u_vel = 0.0;
        if (rho > 0.0) {
            u_vel = u.rhoU / rho;
        }

        const double kinetic = 0.5 * rho * u_vel * u_vel;
        const double internal = p_min / (gamma - 1.0);
        u.E = internal + kinetic;
    }
}
