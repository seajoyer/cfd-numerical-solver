#ifndef WENORECONSTRUCTION_HPP
#define WENORECONSTRUCTION_HPP

#include "reconstruction/Reconstruction.hpp"

/**
 * @class WENOReconstruction
 * @brief WENO-JS reconstruction (WENO3 or WENO5) for axis-aligned faces on primitives W.
 *
 * Supports nominal orders:
 *  - 3 (WENO3-JS)
 *  - 5 (WENO5-JS)
 *
 * Input:
 *  - W(var,i,j,k) = (rho,u,v,w,P), shape (5,sx,sy,sz)
 *
 * For each face with LEFT cell index (i,j,k):
 *  - WL = left-biased reconstruction at the face i+1/2 (from cell i)
 *  - WR = right-biased reconstruction at the same face (from cell i+1)
 *
 * Contract:
 *  - No allocations in hot path.
 *  - Ghost cells provide required stencil support.
 */
class WENOReconstruction final : public Reconstruction {
public:
    explicit WENOReconstruction(int order);

    ~WENOReconstruction() override = default;

    [[nodiscard]] int GetOrder() const {
        return order_;
    }

    void SetOrder(int order);

    [[nodiscard]] double GetEpsilon() const {
        return epsilon_;
    }

    void SetEpsilon(double epsilon);

    [[nodiscard]] int GetNonlinearWeightPower() const {
        return p_;
    }

    void SetNonlinearWeightPower(int p);

    void ReconstructFace(const xt::xtensor<double, 4>& W,
                         Axis axis,
                         int i, int j, int k,
                         PrimitiveCell& WL,
                         PrimitiveCell& WR) const override;

private:
    int order_{5}; // 3 or 5
    double epsilon_{1e-6};
    int p_{2};

    static double PowInt(double x, int p);
    static void ComputeNonlinearWeights(const double* beta, const double* d, int r, double eps, int p, double* omega);

    static double LoadAt(const xt::xtensor<double, 4>& W, std::size_t v,
                         int i, int j, int k, const AxisStride& st, int offset);

    static void StorePrimitiveComponent(PrimitiveCell& w, std::size_t v, double val);

    static double Weno3Left(double fim1, double fi, double fip1, double eps, int p);
    static double Weno3Right(double fi, double fip1, double fip2, double eps, int p);

    static double Weno5Left(double fim2, double fim1, double fi, double fip1, double fip2, double eps, int p);
    static double Weno5Right(double fim1, double fi, double fip1, double fip2, double fip3, double eps, int p);
};

#endif  // WENORECONSTRUCTION_HPP
