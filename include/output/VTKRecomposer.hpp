#ifndef VTKRECOMPOSER_HPP
#define VTKRECOMPOSER_HPP

#include <string>

#include "config/Settings.hpp"

/**
 * @class VTKRecomposer
 * @brief Parallel recomposition of per-rank VTK structured-grid files.
 *
 * Each MPI rank processes only its assigned subset of step files:
 *   step_index % mpi_size == mpi_rank
 *
 * Input layout:
 *   output_dir/
 *     rank_0000/
 *     rank_0001/
 *     ...
 *
 * Output layout:
 *   output_dir/recomposed/
 *     <solver>__R_<reconstruction>__N_<Nx>x<Ny>x<Nz>__CFL_<...>__step_XXXX.vtk
 */
class VTKRecomposer final {
public:
    VTKRecomposer(std::string output_dir, int rank, int size);

    /**
     * @brief Recompose only the subset of step files assigned to this rank.
     */
    void RecomposeAssigned(const Settings& settings) const;

    static auto ExtractStepNumber(const std::string& filename) -> int;
private:
    std::string output_dir_;
    int rank_ = 0;
    int size_ = 1;

};

#endif  // VTKRECOMPOSER_HPP