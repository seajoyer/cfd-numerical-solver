#ifndef MADERSPATIALOPERATOR_HPP
#define MADERSPATIALOPERATOR_HPP

#include <cstddef>
#include <memory>

#include "config/Settings.hpp"
#include "solver/EOS.hpp"
#include "spatial/SpatialOperator.hpp"

class DataLayer;
class Mesh;
class Workspace;
class BoundaryManager;

/**
 * @brief Chemistry controls for Mader phase I.
 */
struct MaderChemistryParameters final {
    bool enabled = true;

    double z_freq = 5.0e15; // Arrhenius pre-exponential factor
    double activation_energy = 2.0e4; // E*
    double gas_constant = 1.0; // R_g

    double heat_release = 4.84e6; // Q_chem

    double min_temperature = 1500.0; // MINWT
    double min_reactant = 1.0e-6; // GASW
    std::size_t delay_cycles = 0; // VCNT

    double reactant_floor = 0.0;
    double reactant_initial = 1.0;
};

/**
 * @brief Artificial viscosity controls.
 *
 * Current target:
 *  - slab geometry only
 *  - large-particle / Mader-style staggered update
 */
struct MaderViscosityParameters final {
    bool enabled = true;
    double coefficient = 0.5; // analogous to VISC / K_B depending on chosen form
};

/**
 * @brief Transport method selector for phase IV.
 */
enum class MaderTransportMethod : std::size_t {
    DonorAcceptor = 0,
    Alternative = 1
};

enum class MaderCompositionTransportMode : std::size_t {
    Standard = 0,
    Shargatov = 1
};

struct MaderTransportParameters final {
    MaderTransportMethod method = MaderTransportMethod::DonorAcceptor;

    bool slab_geometry = true;
    MaderCompositionTransportMode composition_mode =
        MaderCompositionTransportMode::Standard;
};

/**
 * @brief Transport controls for phase IV.
 */
// struct MaderTransportParameters final {
//     MaderTransportMethod method = MaderTransportMethod::DonorAcceptor;
//
//     bool use_shargatov_correction = false; // reserved for later
//     bool slab_geometry = true; // current implementation target
// };

/**
 * @class MaderSpatialOperator
 * @brief Five-phase Mader / large-particle operator.
 *
 * Important:
 *  - This operator is phase-split.
 *  - It is not naturally represented as a single FV RHS.
 *  - ComputeRHS() remains only for base-interface compatibility.
 *
 * Variable placement intended by this class:
 *  - cell-centered: rho, I, P, T, lambda, transport accumulators
 *  - face-centered: Ux, Vy (, Wz later if needed)
 */
class MaderSpatialOperator final : public SpatialOperator {
public:
    explicit MaderSpatialOperator(const Settings& settings, std::shared_ptr<BoundaryManager> boundary_manager);

    /**
     * @brief Compatibility stub.
     *
     * Mader scheme is implemented by explicit phases, not by a single RHS operator.
     */
    void ComputeRHS(DataLayer& layer,
                    const Mesh& mesh,
                    Workspace& workspace,
                    double gamma,
                    double dt) const override;
    // CHECK: MADER_ARRHENIUS
    // ------------------------------------------------------------
    // Phase I: EOS + chemistry
    // ------------------------------------------------------------
    void Phase1_EosAndChemistry(DataLayer& layer,
                                const Mesh& mesh,
                                Workspace& workspace,
                                double gamma,
                                double dt) const;
    // CHECK: MADER_VISC
    // CHECK: MADER_VELOCITY
    // ------------------------------------------------------------
    // Phase II: pressure forces + artificial viscosity
    // updates face-centered velocities
    // ------------------------------------------------------------
    void Phase2_PressureForces(DataLayer& layer,
                               const Mesh& mesh,
                               Workspace& workspace,
                               double gamma,
                               double dt) const;
    // CHECK: MADER_ZIP
    // ------------------------------------------------------------
    // Phase III: ZIP energy + intermediate density
    // ------------------------------------------------------------
    void Phase3_ZipEnergy(DataLayer& layer,
                          const Mesh& mesh,
                          Workspace& workspace,
                          double gamma,
                          double dt) const;
    // CHECK: MADER_DONOR
    // ------------------------------------------------------------
    // Phase IV: transport accumulation
    // ------------------------------------------------------------
    void Phase4_Transport(DataLayer& layer,
                          const Mesh& mesh,
                          Workspace& workspace,
                          double gamma,
                          double dt) const;
    // CHECK: MADER_REPARTITION
    // ------------------------------------------------------------
    // Phase V: repartition / apply accumulated increments
    // ------------------------------------------------------------
    void Phase5_Finalize(DataLayer& layer,
                         const Mesh& mesh,
                         Workspace& workspace,
                         double gamma,
                         double dt) const;

    // -------------------- configuration --------------------

    void SetEos(std::shared_ptr<EOS> eos);
    [[nodiscard]] std::shared_ptr<EOS> GetEos() const;

    void SetChemistryParameters(const MaderChemistryParameters& params);
    [[nodiscard]] const MaderChemistryParameters& GetChemistryParameters() const;

    void SetViscosityParameters(const MaderViscosityParameters& params);
    [[nodiscard]] const MaderViscosityParameters& GetViscosityParameters() const;

    void SetTransportParameters(const MaderTransportParameters& params);
    [[nodiscard]] const MaderTransportParameters& GetTransportParameters() const;

private:
    // -------------------- generic helpers --------------------

    void ApplyBoundaryConditions(DataLayer& layer, const Mesh& mesh) const;

    // -------------------- phase I helpers --------------------

    /**
     * @brief Fill cell-centered thermodynamic fields on local padded mesh.
     *
     * Expected outputs in Workspace:
     *  - W()(rho, u_cell, v_cell, w_cell, P)
     *  - InternalEnergy()
     *  - Temperature()
     *  - ReactantMassFraction()
     */
    void ComputeCellCenteredThermodynamics(DataLayer& layer,
                                           const Mesh& mesh,
                                           Workspace& workspace,
                                           double gamma,
                                           double dt) const;

    /**
     * @brief Initialize face-centered velocities from conservative state if needed.
     *
     * For staggered Mader scheme this usually means reconstructing/interface-assigning
     * Ux, Vy from current state at the beginning of the step or first cycle.
     */
    void InitializeFaceVelocitiesIfNeeded(DataLayer& layer,
                                          const Mesh& mesh,
                                          Workspace& workspace) const;

    // -------------------- phase II helpers --------------------

    /**
     * @brief Save old face-centered velocities before pressure update.
     */
    void SaveOldFaceVelocities(Workspace& workspace) const;

    /**
     * @brief Compute artificial viscosity in the chosen large-particle form.
     *
     * Current storage target:
     *  - Workspace::Q() is cell-centered.
     */
    void ComputeArtificialViscosity(const Mesh& mesh,
                                    Workspace& workspace) const;

    /**
     * @brief Update x-face velocities Ux by pressure/viscosity forces.
     */
    void UpdateUxFaces(const Mesh& mesh,
                       Workspace& workspace,
                       double dt) const;

    /**
     * @brief Update y-face velocities Vy by pressure/viscosity forces.
     */
    void UpdateVyFaces(const Mesh& mesh,
                       Workspace& workspace,
                       double dt) const;

    // -------------------- phase III helpers --------------------

    /**
     * @brief Compute intermediate density on cells from staggered velocity divergence.
     */
    void ComputeIntermediateDensity(const Mesh& mesh,
                                    Workspace& workspace,
                                    double dt) const;

    /**
     * @brief Update internal energy using ZIP energy equation.
     */
    void UpdateInternalEnergyZip(const Mesh& mesh,
                                 Workspace& workspace,
                                 double dt) const;

    /**
     * @brief Rebuild total conservative energy from updated I and staggered/cell fields.
     */
    void RebuildConservativeEnergy(DataLayer& layer,
                                   const Mesh& mesh,
                                   Workspace& workspace) const;

    // -------------------- phase IV helpers --------------------

    /**
     * @brief Reset transport accumulators DM/DE/DW/DPU/DPV.
     */
    void ZeroTransportAccumulators(Workspace& workspace) const;

    /**
     * @brief Main donor-acceptor transport dispatcher.
     */
    void TransportDonorAcceptor(const Mesh& mesh,
                                Workspace& workspace,
                                DataLayer& layer,
                                double dt) const;

    /**
     * @brief Donor-acceptor transport through x/R interfaces.
     */
    void TransportDonorAcceptorR(const Mesh& mesh,
                                 Workspace& workspace,
                                 DataLayer& layer,
                                 double dt) const;

    /**
     * @brief Donor-acceptor transport through y/Z interfaces.
     */
    void TransportDonorAcceptorZ(const Mesh& mesh,
                                 Workspace& workspace,
                                 DataLayer& layer,
                                 double dt) const;

    // -------------------- phase V helpers --------------------

    /**
     * @brief Apply accumulated transport increments to cell-centered state.
     */
    void ApplyTransportAccumulators(DataLayer& layer,
                                    const Mesh& mesh,
                                    Workspace& workspace,
                                    double gamma) const;

    /**
     * @brief Reconstruct cell-centered primitive fields from updated conservative state.
     */
    void ReconstructCellFieldsFromConservative(DataLayer& layer,
                                               const Mesh& mesh,
                                               Workspace& workspace,
                                               double gamma) const;

    Settings settings_;
    std::shared_ptr<EOS> eos_;
    MaderChemistryParameters chemistry_params_;
    MaderViscosityParameters viscosity_params_;
    MaderTransportParameters transport_params_;

    mutable std::size_t cycle_counter_ = 0;

    [[nodiscard]] bool IsPureA(double w) const;
    [[nodiscard]] bool IsPureB(double w) const;
    [[nodiscard]] bool IsMixed(double w) const;

    [[nodiscard]] double ComputeTransferredReactantFractionStandard(
        double donor_w) const;
    // CHECK: SHARGATOV
    [[nodiscard]] double ComputeTransferredReactantFractionShargatovR(
        const Mesh& mesh,
        const DataLayer& layer,
        int donor_i, int donor_j, int donor_k,
        int accept_i, int accept_j, int accept_k) const;

    [[nodiscard]] double ComputeTransferredReactantFractionShargatovZ(
        const Mesh& mesh,
        const DataLayer& layer,
        int donor_i, int donor_j, int donor_k,
        int accept_i, int accept_j, int accept_k) const;
};

#endif  // MADERSPATIALOPERATOR_HPP
