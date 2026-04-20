#ifndef STATESYNCHRONIZER_HPP
#define STATESYNCHRONIZER_HPP

class DataLayer;

/**
 * @class StateSynchronizer
 * @brief Abstract interface for synchronizing ghost-cell state before flux evaluations.
 */
class StateSynchronizer {
public:
    virtual ~StateSynchronizer() = default;

    /**
     * @brief Synchronize ghost-cell state in-place.
     */
    virtual void Synchronize(DataLayer& layer) const = 0;
};

#endif  // STATESYNCHRONIZER_HPP
