#ifndef GODUNOVKOLGANRODIONOVSPATIALOPERATOR_HPP
#define GODUNOVKOLGANRODIONOVSPATIALOPERATOR_HPP

#include <memory>

#include "config/Settings.hpp"
#include "spatial/SpatialOperator.hpp"
#include "data/Variables.hpp"

class Reconstruction;
class RiemannSolver;
class ArtificialViscosity;

/**
 * @class GodunovKolganRodionovSpatialOperator
 * @brief Second-order Godunov–Kolgan–Rodionov (MUSCL–Hancock) spatial operator (axis-aligned).
 *
 * Pipeline (per ComputeRHS):
 *  1) UpdateHalo(U) -> ApplyPhysicalBc(U)
 *  2) ConvertUtoW(U -> workspace.W) on full padded domain
 *  3) rhs = 0
 *  4) For each axis:
 *      - reconstruct face states (WL/WR)
 *      - predictor: evolve each cell’s left/right edge conservative states by dt/2 using local flux difference
 *      - solve Riemann on predicted interface states
 *      - accumulate flux divergence into rhs
 *
 * dt is provided externally via SetLocalTimeStep(dt) (required).
 */
class GodunovKolganRodionovSpatialOperator final : public SpatialOperator {
public:
    GodunovKolganRodionovSpatialOperator(const Settings& settings,
                                         std::shared_ptr<BoundaryManager> boundary_manager);

    void ComputeRHS(DataLayer& layer, Workspace& workspace, double gamma, double dt) const override;

private:
    std::shared_ptr<Reconstruction> reconstruction_;
    std::shared_ptr<RiemannSolver> riemann_solver_;
    std::shared_ptr<ArtificialViscosity> viscosity_; // optional; may be null

    void InitializeReconstruction(const Settings& settings);
    void InitializeRiemannSolver(const Settings& settings);

    void AccumulateAxis(DataLayer& layer,
                        const xt::xtensor<double, 4>& W,
                        xt::xtensor<double, 4>& rhs,
                        double gamma,
                        Axis axis,
                        double dt) const;

    [[nodiscard]] double InvMetricAt(const DataLayer& layer, Axis axis, int i, int j, int k) const;

    static void ApplyPredictor(ConservativeCell& UL,
                               ConservativeCell& UR,
                               const FluxCell& FL,
                               const FluxCell& FR,
                               double half_dt_over_d);

    void MapCellIndices(Axis axis, int s, int a, int b, int& i, int& j, int& k) const;

    void ComputeStarForCell(const xt::xtensor<double, 4>& W,
                            double gamma,
                            Axis axis,
                            const AxisStride& st,
                            const DataLayer& layer,
                            int ci, int cj, int ck,
                            PrimitiveCell& WL_face,
                            PrimitiveCell& WR_face,
                            PrimitiveCell& WLm_face,
                            PrimitiveCell& WRm_face,
                            ConservativeCell& UL_star,
                            ConservativeCell& UR_star,
                            double dt) const;

    static void AccumulateCellRhs(xt::xtensor<double, 4>& rhs,
                                  int i, int j, int k,
                                  double invd,
                                  const FluxCell& Fp,
                                  const FluxCell& Fm);
};

#endif  // GODUNOVKOLGANRODIONOVSPATIALOPERATOR_HPP
