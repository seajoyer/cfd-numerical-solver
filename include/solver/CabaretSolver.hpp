#ifndef CABARETSOLVER_HPP
#define CABARETSOLVER_HPP

#include <array>
#include <memory>

#include "Solver.hpp"
#include "bc/BoundaryCondition.hpp"
#include "bc/BoundaryManager.hpp"
#include "config/Settings.hpp"
#include "data/DataLayer.hpp"
#include "data/Variables.hpp"

class CabaretSolver : public Solver {
public:
    explicit CabaretSolver(const Settings& settings);

    auto Step(DataLayer& layer, double& t_cur) -> double override;

    void SetCfl(double cfl) override;

    void AddBoundary(int axis,
                     std::shared_ptr<BoundaryCondition> left_bc,
                     std::shared_ptr<BoundaryCondition> right_bc) override;

private:
    using CharVec = std::array<double, 3>;

    struct RoeMatrices {
        double L[3][3];
        double R[3][3];
        double lambda[3];
        double a;
        double u;
    };

    Settings settings_;
    BoundaryManager boundary_manager_;

    xt::xarray<Conservative> psi_half_prev_;
    bool psi_initialized_ = false;

    double rho_min_ = 1e-10;
    double p_min_ = 1e-10;

    [[nodiscard]] auto ComputeDx(const DataLayer& layer) const -> double;

    void StoreConservativeCell(const Conservative& uc, int i, double dx, DataLayer& layer) const;

    void EnsurePsiStorage(const xt::xarray<Conservative>& U,
                          const xt::xarray<Flux>& F,
                          double dt_over_dx);

    [[nodiscard]] static auto BuildRoeMatrices(const Primitive& left,
                                              const Primitive& right,
                                              double gamma) -> RoeMatrices;

    [[nodiscard]] static auto MatVec(const double A[3][3], const Conservative& u) -> CharVec;
    [[nodiscard]] static auto MatVec(const double A[3][3], const CharVec& w) -> Conservative;

    [[nodiscard]] static auto CabaretLimit(double w_half,
                                          double w_upwind,
                                          double w_downwind,
                                          double w_min,
                                          double w_max) -> double;

    [[nodiscard]] static auto MaxSignalSpeed(double u, double a) -> double;
};

#endif  // CABARETSOLVER_HPP
