#include "solver/PositivityLimiter.hpp"

void PositivityLimiter::Apply(Conservative &u,
                              double gamma,
                              double rhoMin,
                              double pMin) {
    if (u.rho < rhoMin) {
        double uVel = 0.0;
        if (u.rho > 0.0) {
            uVel = u.rhoU / u.rho;
        }
        u.rho = rhoMin;
        u.rhoU = u.rho * uVel;
    }

    double p = EOS::Pressure(u, gamma);

    if (p < pMin) {
        const double rho = u.rho;
        double uVel = 0.0;
        if (rho > 0.0) {
            uVel = u.rhoU / rho;
        }

        const double kinetic = 0.5 * rho * uVel * uVel;
        const double internal = pMin / (gamma - 1.0);
        u.E = internal + kinetic;
    }
}
