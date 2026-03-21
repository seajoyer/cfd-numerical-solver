#ifndef INITIALCONDITIONS_HPP
#define INITIALCONDITIONS_HPP

#include <cstddef>
#include <string>
#include <vector>

#include "config/Settings.hpp"

/**
 * @file InitialConditions.hpp
 * @brief Structured initial conditions in region-based form
 */

struct Field3DValues {
    std::vector<std::vector<std::vector<double>>> values;

    [[nodiscard]] auto Nx() const -> std::size_t {
        return values.size();
    }

    [[nodiscard]] auto Ny() const -> std::size_t {
        return values.empty() ? 0 : values[0].size();
    }

    [[nodiscard]] auto Nz() const -> std::size_t {
        return values.empty() || values[0].empty() ? 0 : values[0][0].size();
    }

    [[nodiscard]] auto At(std::size_t ix, std::size_t iy, std::size_t iz) const -> double {
        return values.at(ix).at(iy).at(iz);
    }

    [[nodiscard]] auto Empty() const -> bool {
        return values.empty();
    }
};

struct InitialConditions {
    /**
     * Current supported values:
     * - structured_regions
     * - region_markers (reserved for future unstructured support)
     */
    std::string ic_type = "structured_regions";

    /**
     * Region interfaces along axes.
     * Number of regions along each axis is interfaces.size() + 1.
     */
    std::vector<double> interfaces_x;
    std::vector<double> interfaces_y;
    std::vector<double> interfaces_z;

    /**
     * Structured region values stored in normalized 3D form:
     * - 1D: [nx][1][1]
     * - 2D: [nx][ny][1]
     * - 3D: [nx][ny][nz]
     */
    Field3DValues rho;
    Field3DValues u;
    Field3DValues v;
    Field3DValues w;
    Field3DValues p;

    /**
     * Case-local runtime overrides.
     */
    CaseSettings overrides;

    [[nodiscard]] auto RegionCountX() const -> std::size_t {
        return interfaces_x.size() + 1;
    }

    [[nodiscard]] auto RegionCountY() const -> std::size_t {
        return interfaces_y.size() + 1;
    }

    [[nodiscard]] auto RegionCountZ() const -> std::size_t {
        return interfaces_z.size() + 1;
    }
};

#endif  // INITIALCONDITIONS_HPP