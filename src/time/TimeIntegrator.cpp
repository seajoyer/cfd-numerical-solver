#include "time/TimeIntegrator.hpp"
#include "solver/EOS.hpp"


void TimeIntegrator::SetPositivityThresholds(double rho_min, double p_min) {
    rho_min_ = rho_min;
    p_min_ = p_min;
}


void TimeIntegrator::StoreConservativeCell(const Conservative& uc,
                                           const int i,
                                           const double dx,
                                           const Settings& settings,
                                           DataLayer& layer) const {
    const double rho = uc.rho;
    const double rhoU = uc.rhoU;
    const double uvel = rho > 0.0 ? rhoU / rho : 0.0;
    const double P = EOS::Pressure(uc, settings.gamma);

    layer.rho(i) = rho;
    layer.u(i) = uvel;
    layer.P(i) = P;

    layer.p(i) = rhoU;
    layer.V(i) = rho > 0.0 ? 1.0 / rho : 0.0;

    const double kinetic = 0.5 * rho * uvel * uvel;
    const double Eint = uc.E - kinetic;
    const double eint = rho > 0.0 ? Eint / rho : 0.0;
    const double Etot = rho > 0.0 ? uc.E : 0.0;

    layer.U(i) = eint;
    layer.e(i) = Etot;
    layer.m(i) = rho * dx;
}