#include "parallel/MPIContext.hpp"

MPIContext::MPIContext(const MPI_Comm comm, const bool owns_lifetime)
    : comm_(comm), owns_lifetime_(owns_lifetime) {
    if (!IsInitialized()) {
        if (owns_lifetime_) {
            int argc = 0;
            char** argv = nullptr;
            Initialize(argc, argv);
        }
        else {
            throw std::runtime_error("MPIContext: MPI is not initialized");
        }
    }

    if (IsFinalized()) {
        throw std::runtime_error("MPIContext: MPI is already finalized");
    }

    RefreshRankAndSize();
}

MPIContext::~MPIContext() {
    if (owns_lifetime_ && IsInitialized() && !IsFinalized()) {
        Finalize();
    }
}

void MPIContext::Initialize(int& argc, char**& argv) {
    int initialized = 0;
    MPI_Initialized(&initialized);
    if (!initialized) {
        MPI_Init(&argc, &argv);
    }
}

void MPIContext::Finalize() {
    int initialized = 0;
    MPI_Initialized(&initialized);
    if (!initialized) {
        return;
    }

    int finalized = 0;
    MPI_Finalized(&finalized);
    if (!finalized) {
        MPI_Finalize();
    }
}

bool MPIContext::IsInitialized() {
    int initialized = 0;
    MPI_Initialized(&initialized);
    return initialized != 0;
}

bool MPIContext::IsFinalized() {
    int finalized = 0;
    MPI_Finalized(&finalized);
    return finalized != 0;
}

MPI_Comm MPIContext::Comm() const {
    return comm_;
}

int MPIContext::Rank() const {
    return rank_;
}

int MPIContext::Size() const {
    return size_;
}

bool MPIContext::IsRoot(const int root) const {
    return rank_ == root;
}

void MPIContext::Barrier() const {
    MPI_Barrier(comm_);
}

double MPIContext::GlobalMin(const double value) const {
    double result = value;
    MPI_Allreduce(&value, &result, 1, MPI_DOUBLE, MPI_MIN, comm_);
    return result;
}

double MPIContext::GlobalMax(const double value) const {
    double result = value;
    MPI_Allreduce(&value, &result, 1, MPI_DOUBLE, MPI_MAX, comm_);
    return result;
}

double MPIContext::GlobalSum(const double value) const {
    double result = value;
    MPI_Allreduce(&value, &result, 1, MPI_DOUBLE, MPI_SUM, comm_);
    return result;
}

int MPIContext::GlobalMin(const int value) const {
    int result = value;
    MPI_Allreduce(&value, &result, 1, MPI_INT, MPI_MIN, comm_);
    return result;
}

int MPIContext::GlobalMax(const int value) const {
    int result = value;
    MPI_Allreduce(&value, &result, 1, MPI_INT, MPI_MAX, comm_);
    return result;
}

int MPIContext::GlobalSum(const int value) const {
    int result = value;
    MPI_Allreduce(&value, &result, 1, MPI_INT, MPI_SUM, comm_);
    return result;
}

void MPIContext::RefreshRankAndSize() {
    MPI_Comm_rank(comm_, &rank_);
    MPI_Comm_size(comm_, &size_);
}
