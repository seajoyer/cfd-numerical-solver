#ifndef FINITEVOLUMESOLVER_HPP
#define FINITEVOLUMESOLVER_HPP

#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

#include "config/Settings.hpp"
#include "data/DataLayer.hpp"
#include "data/Workspace.hpp"
#include "geometry/Mesh.hpp"
#include "solver/Solver.hpp"
#include "solver/TimeStepCalculator.hpp"
#include "spatial/SpatialOperator.hpp"
#include "time/TimeIntegrator.hpp"

class MPIContext;
class StateSynchronizer;

class FiniteVolumeSolver final : public Solver {
public:
    FiniteVolumeSolver(const Settings& settings,
                       std::shared_ptr<Mesh> mesh,
                       std::shared_ptr<SpatialOperator> spatial_operator,
                       std::shared_ptr<TimeIntegrator> time_integrator,
                       const MPIContext* mpi_context,
                       const StateSynchronizer* halo_exchange);

    auto Step(DataLayer& layer, double& t_cur) -> double override;
    void SetCfl(double cfl) override;

    [[nodiscard]] const Mesh& GetMesh() const;
    [[nodiscard]] Mesh& GetMesh();

private:
    Settings settings_;
    std::shared_ptr<Mesh> mesh_;

    std::shared_ptr<SpatialOperator> spatial_operator_;
    std::shared_ptr<TimeIntegrator> time_integrator_;
    // std::unique_ptr<SolutionFilter> diffusion_;

    const MPIContext* mpi_context_ = nullptr;
    const StateSynchronizer* halo_exchange_ = nullptr;

    Workspace workspace_;

    double rho_min_ = 1e-10;
    double p_min_ = 1e-10;

    void EnsureWorkspaceSized();
};

#endif  // FINITEVOLUMESOLVER_HPP
