#ifndef INITIALCONDITIONS_HPP
#define INITIALCONDITIONS_HPP

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

#include "config/Settings.hpp"

/**
 * @brief Supported initial-condition representations.
 */
enum class InitialConditionType {
    StructuredRegions,
    Constant,
    RegionMarkers
};

/**
 * @brief Structured region values stored in normalized 3D form.
 *
 * Layout:
 * - 1D: [nx][1][1]
 * - 2D: [nx][ny][1]
 * - 3D: [nx][ny][nz]
 */
struct Field3DValues final {
    std::vector<std::vector<std::vector<double>>> values;

    [[nodiscard]] std::size_t Nx() const {
        return values.size();
    }

    [[nodiscard]] std::size_t Ny() const {
        return values.empty() ? 0 : values[0].size();
    }

    [[nodiscard]] std::size_t Nz() const {
        return values.empty() || values[0].empty() ? 0 : values[0][0].size();
    }

    [[nodiscard]] double At(std::size_t ix, std::size_t iy, std::size_t iz) const {
        return values.at(ix).at(iy).at(iz);
    }

    [[nodiscard]] bool Empty() const {
        return values.empty();
    }
};

/**
 * @brief Region-based structured initial-condition description.
 *
 * Region counts:
 * - x regions = interfaces_x.size() + 1
 * - y regions = interfaces_y.size() + 1
 * - z regions = interfaces_z.size() + 1
 */
struct StructuredRegionInitialCondition final {
    std::vector<double> interfaces_x;
    std::vector<double> interfaces_y;
    std::vector<double> interfaces_z;

    Field3DValues rho;
    Field3DValues u;
    Field3DValues v;
    Field3DValues w;
    Field3DValues p;

    std::optional<Field3DValues> reactant_mass_fraction;

    [[nodiscard]] std::size_t RegionCountX() const {
        return interfaces_x.size() + 1;
    }

    [[nodiscard]] std::size_t RegionCountY() const {
        return interfaces_y.size() + 1;
    }

    [[nodiscard]] std::size_t RegionCountZ() const {
        return interfaces_z.size() + 1;
    }
};

/**
 * @brief Constant primitive initial state over the whole mesh.
 */
struct ConstantInitialCondition final {
    double rho = 0.0;
    double u = 0.0;
    double v = 0.0;
    double w = 0.0;
    double p = 0.0;

    std::optional<double> reactant_mass_fraction;
};

/**
 * @brief Initial-condition container for one case.
 *
 * At the current stage only one representation is expected to be populated
 * according to the selected type.
 */
struct InitialConditions final {
    InitialConditionType type = InitialConditionType::StructuredRegions;

    std::optional<StructuredRegionInitialCondition> structured_regions;
    std::optional<ConstantInitialCondition> constant;

    /**
     * @brief Case-local runtime overrides.
     */
    CaseSettings overrides;
};

#endif  // INITIALCONDITIONS_HPP
