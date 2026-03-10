#ifndef FINITEVOLUMESOLVER_HPP
#define FINITEVOLUMESOLVER_HPP

#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

#include "bc/BoundaryManager.hpp"
#include "config/Settings.hpp"
#include "data/DataLayer.hpp"
#include "data/Mesh.hpp"
#include "data/Workspace.hpp"
#include "filter/SolutionFilter.hpp"
#include "solver/Solver.hpp"
#include "solver/TimeStepCalculator.hpp"
#include "spatial/SpatialOperator.hpp"
#include "time/TimeIntegrator.hpp"

class FiniteVolumeSolver final : public Solver {
public:
    FiniteVolumeSolver(const Settings& settings,
                       Mesh mesh,
                       std::shared_ptr<SpatialOperator> spatial_operator,
                       std::shared_ptr<TimeIntegrator> time_integrator,
                       const MPIContext* mpi_context);

    auto Step(DataLayer& layer, double& t_cur) -> double override;
    void SetCfl(double cfl) override;

    [[nodiscard]] const Mesh& GetMesh() const;
    [[nodiscard]] Mesh& GetMesh();

private:
    Settings settings_;
    Mesh mesh_;

    std::shared_ptr<SpatialOperator> spatial_operator_;
    std::shared_ptr<TimeIntegrator> time_integrator_;
    std::unique_ptr<SolutionFilter> diffusion_;

    const MPIContext* mpi_context_ = nullptr;

    Workspace workspace_;

    double rho_min_ = 1e-10;
    double p_min_ = 1e-10;

    void EnsureWorkspaceSized();
};

#endif  // FINITEVOLUMESOLVER_HPP
