#ifndef EOS_HPP
#define EOS_HPP

#include <stdexcept>

/**
 * @enum EosType
 * @brief Supported equation-of-state models.
 */
enum class EosType {
    IdealGas = 0,
    HugoniotGruneisen = 1,
    Jwl = 2
};

/**
 * @struct EosCellInput
 * @brief Input thermodynamic state for one cell.
 *
 * rho    : density
 * I      : specific internal energy
 * lambda : mass fraction of unreacted component
 */
struct EosCellInput final {
    double rho = 0.0;
    double I = 0.0;
    double lambda = 0.0;
};

/**
 * @struct EosCellOutput
 * @brief Output thermodynamic state for one cell.
 */
struct EosCellOutput final {
    double P = 0.0;
    double T = 0.0;
    double c = 0.0;
};

/**
 * @struct IdealGasEosParameters
 * @brief Parameters for ideal-gas EOS.
 *
 * P = (gamma - 1) rho I
 * T = P / (rho R_specific)
 */
struct IdealGasEosParameters final {
    double gamma = 1.4;
    double R_specific = 1.0;
    double rho_floor = 1e-14;
    double p_floor = 1e-14;
    double T_floor = 0.0;
};

/**
 * @struct HugoniotGruneisenEosParameters
 * @brief Parameters for Hugoniot + Gruneisen EOS.
 *
 * Given:
 *   U_s = C + S U_p
 *   P_H = C^2 (V0 - V) / [V0 - S (V0 - V)]^2
 *   I_H = 0.5 P_H (V0 - V)
 *   P   = P_H + (gamma_s / V) (I - I_H)
 *
 * where:
 *   V = 1 / rho
 *   V0 = 1 / rho0
 *
 * Temperature is not defined by these relations alone, so here we use
 * a simple caloric closure:
 *   T = T_ref + (I - I_ref) / c_v
 */
struct HugoniotGruneisenEosParameters final {
    double rho0 = 1.0; // reference density
    double C = 1.0; // shock Hugoniot parameter
    double S = 1.0; // shock Hugoniot parameter
    double gamma_s = 1.0; // Gruneisen coefficient

    double c_v = 1.0; // caloric closure for temperature
    double T_ref = 300.0;
    double I_ref = 0.0;

    double rho_floor = 1e-14;
    double p_floor = 0.0;
    double T_floor = 0.0;
};

/**
 * @class EOS
 * @brief Cell-wise equation of state evaluator.
 */
class EOS final {
public:
    EOS() = default;

    void SetType(EosType type);
    [[nodiscard]] EosType GetType() const;

    void SetIdealGasParameters(const IdealGasEosParameters& params);
    [[nodiscard]] const IdealGasEosParameters& GetIdealGasParameters() const;

    void SetHugoniotGruneisenParameters(const HugoniotGruneisenEosParameters& params);
    [[nodiscard]] const HugoniotGruneisenEosParameters& GetHugoniotGruneisenParameters() const;

    /**
     * @brief Evaluate EOS for one cell.
     */
    [[nodiscard]] EosCellOutput Evaluate(const EosCellInput& in) const;

private:
    EosType type_ = EosType::IdealGas;

    IdealGasEosParameters ideal_gas_params_;
    HugoniotGruneisenEosParameters hugoniot_gruneisen_params_;

    [[nodiscard]] EosCellOutput EvaluateIdealGas(const EosCellInput& in) const;
    [[nodiscard]] EosCellOutput EvaluateHugoniotGruneisen(const EosCellInput& in) const;
};

#endif  // EOS_HPP
