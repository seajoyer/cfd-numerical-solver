#ifndef HALOEXCHANGE_HPP
#define HALOEXCHANGE_HPP

#include <vector>

#include "data/DataLayer.hpp"
#include "data/Variables.hpp"
#include "bc/BoundaryCondition.hpp"
#include "data/Mesh.hpp"
#include "data/PressureVelocityState.hpp"
#include "data/Variables.hpp"
#include "parallel/MPIContext.hpp"
// CHECK: HALO_EXCHANGE
/**
 * @class HaloExchange
 * @brief MPI halo exchange for structured Cartesian subdomains.
 *
 * Supported exchanged data:
 *  - conservative DataLayer U(var,i,j,k)
 *  - pressure-velocity state:
 *      * pressure(i,j,k)
 *      * ux(i_face,j,k)
 *      * vy(i,j_face,k)
 *      * wz(i,j,k_face)
 */
class HaloExchange final {
public:
    HaloExchange(MPI_Comm comm, int size, bool exchange_reactant_mass_fraction);

    void Exchange(DataLayer& layer, const Mesh& mesh) const;
    void Exchange(PressureVelocityState& state, const Mesh& mesh) const;

private:
    MPI_Comm comm_ = MPI_COMM_NULL;
    int size_ = 1;

    bool exchange_reactant_mass_fraction_ = false;

    void ExchangeX(DataLayer& layer, const Mesh& mesh) const;
    void ExchangeY(DataLayer& layer, const Mesh& mesh) const;
    void ExchangeZ(DataLayer& layer, const Mesh& mesh) const;

    [[nodiscard]] std::vector<double> PackX(const DataLayer& layer, const Mesh& mesh, int i_begin) const;
    [[nodiscard]] std::vector<double> PackY(const DataLayer& layer, const Mesh& mesh, int j_begin) const;
    [[nodiscard]] std::vector<double> PackZ(const DataLayer& layer, const Mesh& mesh, int k_begin) const;

    void UnpackX(DataLayer& layer, const Mesh& mesh, int i_begin, const std::vector<double>& buffer) const;
    void UnpackY(DataLayer& layer, const Mesh& mesh, int j_begin, const std::vector<double>& buffer) const;
    void UnpackZ(DataLayer& layer, const Mesh& mesh, int k_begin, const std::vector<double>& buffer) const;

    void ExchangeX(PressureVelocityState& state, const Mesh& mesh) const;
    void ExchangeY(PressureVelocityState& state, const Mesh& mesh) const;
    void ExchangeZ(PressureVelocityState& state, const Mesh& mesh) const;

    [[nodiscard]] std::vector<double> PackX(const PressureVelocityState& state, const Mesh& mesh, int i_begin) const;
    [[nodiscard]] std::vector<double> PackY(const PressureVelocityState& state, const Mesh& mesh, int j_begin) const;
    [[nodiscard]] std::vector<double> PackZ(const PressureVelocityState& state, const Mesh& mesh, int k_begin) const;

    void UnpackX(PressureVelocityState& state, const Mesh& mesh, int i_begin, const std::vector<double>& buffer) const;
    void UnpackY(PressureVelocityState& state, const Mesh& mesh, int j_begin, const std::vector<double>& buffer) const;
    void UnpackZ(PressureVelocityState& state, const Mesh& mesh, int k_begin, const std::vector<double>& buffer) const;
};

#endif  // HALOEXCHANGE_HPP
