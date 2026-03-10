#ifndef MPICONTEXT_HPP
#define MPICONTEXT_HPP

#include <stdexcept>

#include <mpi.h>

/**
 * @class MPIContext
 * @brief Small RAII wrapper around MPI communicator and basic collective operations.
 *
 * This class is responsible for:
 *  - optional MPI initialization/finalization ownership
 *  - storing communicator handle
 *  - exposing rank/size helpers
 *  - providing basic collectives needed by the solver
 *
 * By default it uses MPI_COMM_WORLD.
 */
class MPIContext final {
public:
    /**
     * @brief Construct MPI context on a communicator.
     * @param comm MPI communicator, default is MPI_COMM_WORLD.
     * @param owns_lifetime If true, Initialize/Finalize are managed by this object.
     */
    explicit MPIContext(MPI_Comm comm = MPI_COMM_WORLD, bool owns_lifetime = false);

    ~MPIContext();

    MPIContext(const MPIContext&) = delete;
    MPIContext& operator=(const MPIContext&) = delete;

    MPIContext(MPIContext&&) = delete;
    MPIContext& operator=(MPIContext&&) = delete;

    /**
     * @brief Initialize MPI if not initialized yet.
     * @param argc Command-line argc.
     * @param argv Command-line argv.
     */
    static void Initialize(int& argc, char**& argv);

    /**
     * @brief Finalize MPI if not finalized yet.
     */
    static void Finalize();

    /**
     * @brief Check whether MPI is initialized.
     */
    [[nodiscard]] static bool IsInitialized();

    /**
     * @brief Check whether MPI is finalized.
     */
    [[nodiscard]] static bool IsFinalized();

    /**
     * @brief Get communicator.
     */
    [[nodiscard]] MPI_Comm Comm() const;

    /**
     * @brief Get rank in communicator.
     */
    [[nodiscard]] int Rank() const;

    /**
     * @brief Get communicator size.
     */
    [[nodiscard]] int Size() const;

    /**
     * @brief Check whether this rank is root.
     * @param root Root rank id, default 0.
     */
    [[nodiscard]] bool IsRoot(int root = 0) const;

    /**
     * @brief Global barrier on communicator.
     */
    void Barrier() const;

    /**
     * @brief Global minimum for double.
     */
    [[nodiscard]] double GlobalMin(double value) const;

    /**
     * @brief Global maximum for double.
     */
    [[nodiscard]] double GlobalMax(double value) const;

    /**
     * @brief Global sum for double.
     */
    [[nodiscard]] double GlobalSum(double value) const;

    /**
     * @brief Global minimum for int.
     */
    [[nodiscard]] int GlobalMin(int value) const;

    /**
     * @brief Global maximum for int.
     */
    [[nodiscard]] int GlobalMax(int value) const;

    /**
     * @brief Global sum for int.
     */
    [[nodiscard]] int GlobalSum(int value) const;

private:
    MPI_Comm comm_ = MPI_COMM_WORLD;
    int rank_ = 0;
    int size_ = 1;
    bool owns_lifetime_ = false;

    void RefreshRankAndSize();
};

#endif  // MPICONTEXT_HPP