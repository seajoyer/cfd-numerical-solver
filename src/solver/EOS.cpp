#include "solver/EOS.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace {
    constexpr double k_eps = 1e-14;
}

void EOS::SetType(const EosType type) {
    type_ = type;
}

EosType EOS::GetType() const {
    return type_;
}

void EOS::SetIdealGasParameters(const IdealGasEosParameters& params) {
    if (params.gamma <= 1.0) {
        throw std::invalid_argument("EOS: ideal-gas gamma must be > 1");
    }
    if (params.R_specific <= 0.0) {
        throw std::invalid_argument("EOS: ideal-gas R_specific must be > 0");
    }
    if (params.rho_floor <= 0.0) {
        throw std::invalid_argument("EOS: ideal-gas rho_floor must be > 0");
    }
    if (params.p_floor < 0.0) {
        throw std::invalid_argument("EOS: ideal-gas p_floor must be >= 0");
    }
    if (params.T_floor < 0.0) {
        throw std::invalid_argument("EOS: ideal-gas T_floor must be >= 0");
    }

    ideal_gas_params_ = params;
}

const IdealGasEosParameters& EOS::GetIdealGasParameters() const {
    return ideal_gas_params_;
}

void EOS::SetHugoniotGruneisenParameters(const HugoniotGruneisenEosParameters& params) {
    if (params.rho0 <= 0.0) {
        throw std::invalid_argument("EOS: Hugoniot-Gruneisen rho0 must be > 0");
    }
    if (params.C <= 0.0) {
        throw std::invalid_argument("EOS: Hugoniot-Gruneisen C must be > 0");
    }
    if (params.S <= 0.0) {
        throw std::invalid_argument("EOS: Hugoniot-Gruneisen S must be > 0");
    }
    if (params.gamma_s < 0.0) {
        throw std::invalid_argument("EOS: Hugoniot-Gruneisen gamma_s must be >= 0");
    }
    if (params.c_v <= 0.0) {
        throw std::invalid_argument("EOS: Hugoniot-Gruneisen c_v must be > 0");
    }
    if (params.rho_floor <= 0.0) {
        throw std::invalid_argument("EOS: Hugoniot-Gruneisen rho_floor must be > 0");
    }
    if (params.p_floor < 0.0) {
        throw std::invalid_argument("EOS: Hugoniot-Gruneisen p_floor must be >= 0");
    }
    if (params.T_floor < 0.0) {
        throw std::invalid_argument("EOS: Hugoniot-Gruneisen T_floor must be >= 0");
    }

    hugoniot_gruneisen_params_ = params;
}

const HugoniotGruneisenEosParameters& EOS::GetHugoniotGruneisenParameters() const {
    return hugoniot_gruneisen_params_;
}

EosCellOutput EOS::Evaluate(const EosCellInput& in) const {
    switch (type_) {
    case EosType::IdealGas:
        return EvaluateIdealGas(in);

    case EosType::HugoniotGruneisen:
        return EvaluateHugoniotGruneisen(in);

    case EosType::Jwl:
        throw std::logic_error("EOS::Evaluate: JWL is not implemented yet");
    }

    throw std::logic_error("EOS::Evaluate: unknown EOS type");
}

EosCellOutput EOS::EvaluateIdealGas(const EosCellInput& in) const {
    const double rho = std::max(in.rho, ideal_gas_params_.rho_floor);
    const double I = in.I;

    const double P_raw = (ideal_gas_params_.gamma - 1.0) * rho * I;
    const double P = std::max(P_raw, ideal_gas_params_.p_floor);

    double T = P / (rho * ideal_gas_params_.R_specific);
    T = std::max(T, ideal_gas_params_.T_floor);

    const double c2 = ideal_gas_params_.gamma * P / rho;
    const double c = std::sqrt(std::max(c2, 0.0));

    EosCellOutput out;
    out.P = P;
    out.T = T;
    out.c = c;
    return out;
}

EosCellOutput EOS::EvaluateHugoniotGruneisen(const EosCellInput& in) const {
    const double rho = std::max(in.rho, hugoniot_gruneisen_params_.rho_floor);
    const double I = in.I;

    const double V = 1.0 / rho;
    const double V0 = 1.0 / hugoniot_gruneisen_params_.rho0;
    const double dV = V0 - V;

    const double denom =
        V0 - hugoniot_gruneisen_params_.S * dV;

    const double denom_safe =
        (std::abs(denom) > k_eps)
            ? denom
            : (denom >= 0.0 ? k_eps : -k_eps);

    double P_H = 0.0;
    double I_H = 0.0;

    if (dV > 0.0) {
        P_H = hugoniot_gruneisen_params_.C * hugoniot_gruneisen_params_.C * dV
            / (denom_safe * denom_safe);
        I_H = 0.5 * P_H * dV;
    }

    const double P_raw =
        P_H + (hugoniot_gruneisen_params_.gamma_s / V) * (I - I_H);
    const double P = std::max(P_raw, hugoniot_gruneisen_params_.p_floor);

    double T =
        hugoniot_gruneisen_params_.T_ref
        + (I - hugoniot_gruneisen_params_.I_ref) / hugoniot_gruneisen_params_.c_v;
    T = std::max(T, hugoniot_gruneisen_params_.T_floor);

    const double drho = std::max(1e-6 * rho, 1e-8);
    const double rho_p = rho + drho;
    const double V_p = 1.0 / rho_p;
    const double dV_p = V0 - V_p;
    const double denom_p = V0 - hugoniot_gruneisen_params_.S * dV_p;
    const double denom_p_safe =
        (std::abs(denom_p) > k_eps)
            ? denom_p
            : (denom_p >= 0.0 ? k_eps : -k_eps);

    double P_H_p = 0.0;
    double I_H_p = 0.0;
    if (dV_p > 0.0) {
        P_H_p = hugoniot_gruneisen_params_.C * hugoniot_gruneisen_params_.C * dV_p
            / (denom_p_safe * denom_p_safe);
        I_H_p = 0.5 * P_H_p * dV_p;
    }
    const double P_p =
        P_H_p + (hugoniot_gruneisen_params_.gamma_s / V_p) * (I - I_H_p);

    const double dPdrho = (P_p - P_raw) / drho;
    const double c = std::sqrt(std::max(dPdrho, 0.0));

    EosCellOutput out;
    out.P = P;
    out.T = T;
    out.c = c;
    return out;
}


double EOS::ComputeInternalEnergy(const double rho, const double P, const double lambda) const {
    (void)lambda;
    switch (type_) {
    case EosType::IdealGas:
        return ComputeInternalEnergyIdealGas(rho, P);

    case EosType::HugoniotGruneisen:
        return ComputeInternalEnergyHugoniotGruneisen(rho, P);

    case EosType::Jwl:
        throw std::logic_error("EOS::ComputeInternalEnergy: JWL is not implemented yet");
    }

    throw std::logic_error("EOS::ComputeInternalEnergy: unknown EOS type");
}

double EOS::ComputeInternalEnergyIdealGas(const double rho, const double P) const {
    const double safe_rho = std::max(rho, ideal_gas_params_.rho_floor);
    const double I = P / ((ideal_gas_params_.gamma - 1.0) * safe_rho);
    return std::max(I, 0.0);
}

double EOS::ComputeInternalEnergyHugoniotGruneisen(const double rho, const double P) const {
    const double safe_rho = std::max(rho, hugoniot_gruneisen_params_.rho_floor);
    const double V = 1.0 / safe_rho;
    const double V0 = 1.0 / hugoniot_gruneisen_params_.rho0;
    const double dV = V0 - V;

    const double denom = V0 - hugoniot_gruneisen_params_.S * dV;
    const double denom_safe = (std::abs(denom) > 1e-14) ? denom : (denom >= 0.0 ? 1e-14 : -1e-14);

    double P_H = 0.0;
    double I_H = 0.0;

    if (dV > 0.0) {
        P_H = hugoniot_gruneisen_params_.C * hugoniot_gruneisen_params_.C * dV / (denom_safe * denom_safe);
        I_H = 0.5 * P_H * dV;
    }

    const double gamma_s = std::max(hugoniot_gruneisen_params_.gamma_s, 1e-6); // Защита от деления на 0
    const double I = I_H + (P - P_H) * V / gamma_s;

    return std::max(I, 0.0);
}