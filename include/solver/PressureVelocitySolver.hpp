#ifndef PRESSUREVELOCITYSOLVER_HPP
#define PRESSUREVELOCITYSOLVER_HPP

#include <memory>

#include "bc/BoundaryManager.hpp"
#include "bc/InternalBoundaryCondition.hpp"
#include "config/Settings.hpp"
#include "data/Mesh.hpp"
#include "data/PressureVelocityState.hpp"
#include "data/PressureVelocityWorkspace.hpp"
#include "bc/InternalBoundaryCondition.hpp"
#include "bc/WallInternalBoundary.hpp"
#include "solver/Solver.hpp"

class MPIContext;

class PressureVelocitySolver : public Solver {
public:
    PressureVelocitySolver(const Settings& settings,
                           Mesh mesh,
                           std::shared_ptr<BoundaryManager> boundary_manager,
                           const MPIContext* mpi_context,
                           bool steady = false);

    ~PressureVelocitySolver() override = default;

    [[nodiscard]] const Mesh& GetMesh() const;
    [[nodiscard]] Mesh& GetMesh();

    [[nodiscard]] const PressureVelocityState& GetState() const;
    [[nodiscard]] PressureVelocityState& GetState();

    void SetCfl(double cfl) override;

    void ApplyImmersedVelocityConstraints();

protected:
    Settings settings_;
    Mesh mesh_;
    std::shared_ptr<BoundaryManager> boundary_manager_;
    const MPIContext* mpi_context_ = nullptr;
    bool steady_;

    PressureVelocityState state_;
    PressureVelocityWorkspace workspace_;

    std::unique_ptr<InternalBoundaryCondition> internal_boundary_condition_;

    void EnsureStorageSized();
    void ApplyHaloAndPhysicalBc();
    void ApplyImmersedMomentumCorrections(double dt, double nu);

    [[nodiscard]] double ComputeDt(double t_cur) const;

    void BuildMomentumCoefficients(double dt, double nu);
    void SolveMomentumPredictor(double dt, double alpha_u);
    void BuildPressureCorrectionEquation(double dt, double nu);
    void SolvePressureCorrection();
    void ApplyPressureCorrection(double alpha_p);
    void CopyStarToState();

    [[nodiscard]] double ComputeMassResidual() const;
    [[nodiscard]] int PressureLinearIndex(int i, int j) const;
};

#endif
