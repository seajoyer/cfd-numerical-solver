#ifndef INITIALCONDITIONINITIALIZER_HPP
#define INITIALCONDITIONINITIALIZER_HPP

#include "config/InitialConditions.hpp"
#include "config/Settings.hpp"

class DataLayer;
class Mesh;

/**
 * @class InitialConditionInitializer
 * @brief Applies parsed initial-condition description to cell-centered solution storage.
 */
class InitialConditionInitializer final {
public:
    InitialConditionInitializer(const Settings& settings,
                                const InitialConditions& initial_conditions);

    /**
     * @brief Fill conservative solution fields on the given mesh.
     */
    void Apply(const Mesh& mesh, DataLayer& layer) const;

private:
    const Settings& settings_;
    const InitialConditions& initial_conditions_;

    void ApplyStructuredRegions(const Mesh& mesh, DataLayer& layer) const;
    void ApplyConstant(const Mesh& mesh, DataLayer& layer) const;

    [[nodiscard]] static std::size_t RegionIndex(double coord,
                                                 const std::vector<double>& interfaces);

    static void ValidateStructuredRegionShape(const StructuredRegionInitialCondition& ic);
};

#endif  // INITIALCONDITIONINITIALIZER_HPP
