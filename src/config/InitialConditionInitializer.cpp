#include "config/InitialConditionInitializer.hpp"

#include <algorithm>
#include <stdexcept>

#include "data/DataLayer.hpp"
#include "data/Variables.hpp"
#include "geometry/Cell.hpp"
#include "geometry/Mesh.hpp"

InitialConditionInitializer::InitialConditionInitializer(
    const Settings& settings,
    const InitialConditions& initial_conditions
) : settings_(settings),
    initial_conditions_(initial_conditions) {}

void InitialConditionInitializer::Apply(const Mesh& mesh, DataLayer& layer) const {
    if (!layer.IsAllocated()) {
        throw std::runtime_error("InitialConditionInitializer: DataLayer is not allocated");
    }

    if (layer.GetCellCount() != mesh.GetCellCount()) {
        throw std::runtime_error("InitialConditionInitializer: DataLayer cell count does not match mesh");
    }

    switch (initial_conditions_.type) {
    case InitialConditionType::StructuredRegions:
        ApplyStructuredRegions(mesh, layer);
        return;

    case InitialConditionType::Constant:
        ApplyConstant(mesh, layer);
        return;

    case InitialConditionType::RegionMarkers:
        throw std::runtime_error(
                                 "InitialConditionInitializer: region_markers initial condition is not implemented yet"
                                );
    }

    throw std::runtime_error("InitialConditionInitializer: unsupported initial condition type");
}

void InitialConditionInitializer::ApplyStructuredRegions(const Mesh& mesh, DataLayer& layer) const {
    if (!initial_conditions_.structured_regions.has_value()) {
        throw std::runtime_error(
                                 "InitialConditionInitializer: structured_regions data is missing"
                                );
    }

    const StructuredRegionInitialCondition& ic = *initial_conditions_.structured_regions;
    ValidateStructuredRegionShape(ic);

    auto& U = layer.U();
    auto& lambda = layer.ReactantMassFraction();

    for (std::size_t cell_id = 0; cell_id < mesh.GetCellCount(); ++cell_id) {
        const Cell& cell = mesh.GetCell(cell_id);

        const std::size_t ix = RegionIndex(cell.center_x, ic.interfaces_x);
        const std::size_t iy = (settings_.mesh.dim >= 2)
                                   ? RegionIndex(cell.center_y, ic.interfaces_y)
                                   : 0;
        const std::size_t iz = (settings_.mesh.dim >= 3)
                                   ? RegionIndex(cell.center_z, ic.interfaces_z)
                                   : 0;

        PrimitiveCell primitive;
        primitive.rho = ic.rho.At(ix, iy, iz);
        primitive.u = ic.u.At(ix, iy, iz);
        primitive.v = (settings_.mesh.dim >= 2) ? ic.v.At(ix, iy, iz) : 0.0;
        primitive.w = (settings_.mesh.dim >= 3) ? ic.w.At(ix, iy, iz) : 0.0;
        primitive.P = ic.p.At(ix, iy, iz);

        const ConservativeCell conservative =
            ConservativeFromPrimitive(primitive, settings_.gamma);

        U(cell_id, DataLayer::k_rho) = conservative.rho;
        U(cell_id, DataLayer::k_rhoU) = conservative.rhoU;
        U(cell_id, DataLayer::k_rhoV) = conservative.rhoV;
        U(cell_id, DataLayer::k_rhoW) = conservative.rhoW;
        U(cell_id, DataLayer::k_E) = conservative.E;

        if (ic.reactant_mass_fraction.has_value()) {
            lambda(cell_id) = ic.reactant_mass_fraction->At(ix, iy, iz);
        } else {
            lambda(cell_id) = 0.0;
        }
    }
}

void InitialConditionInitializer::ApplyConstant(const Mesh& mesh, DataLayer& layer) const {
    (void)mesh;

    if (!initial_conditions_.constant.has_value()) {
        throw std::runtime_error("InitialConditionInitializer: constant initial condition is missing");
    }

    const ConstantInitialCondition& ic = *initial_conditions_.constant;

    PrimitiveCell primitive;
    primitive.rho = ic.rho;
    primitive.u = ic.u;
    primitive.v = (settings_.mesh.dim >= 2) ? ic.v : 0.0;
    primitive.w = (settings_.mesh.dim >= 3) ? ic.w : 0.0;
    primitive.P = ic.p;

    const ConservativeCell conservative =
        ConservativeFromPrimitive(primitive, settings_.gamma);

    auto& U = layer.U();
    auto& lambda = layer.ReactantMassFraction();

    for (std::size_t cell_id = 0; cell_id < layer.GetCellCount(); ++cell_id) {
        U(cell_id, DataLayer::k_rho) = conservative.rho;
        U(cell_id, DataLayer::k_rhoU) = conservative.rhoU;
        U(cell_id, DataLayer::k_rhoV) = conservative.rhoV;
        U(cell_id, DataLayer::k_rhoW) = conservative.rhoW;
        U(cell_id, DataLayer::k_E) = conservative.E;

        lambda(cell_id) = ic.reactant_mass_fraction.value_or(0.0);
    }
}

std::size_t InitialConditionInitializer::RegionIndex(const double coord,
                                                     const std::vector<double>& interfaces) {
    return static_cast<std::size_t>(
        std::upper_bound(interfaces.begin(), interfaces.end(), coord) - interfaces.begin()
    );
}

void InitialConditionInitializer::ValidateStructuredRegionShape(
    const StructuredRegionInitialCondition& ic
) {
    const std::size_t nx = ic.RegionCountX();
    const std::size_t ny = ic.RegionCountY();
    const std::size_t nz = ic.RegionCountZ();

    auto validate = [&](const Field3DValues& field, const char* name) {
        if (field.Nx() != nx) {
            throw std::runtime_error(
                                     std::string("InitialConditionInitializer: ") + name +
                                     " x-shape does not match interface count"
                                    );
        }
        if (field.Ny() != ny) {
            throw std::runtime_error(
                                     std::string("InitialConditionInitializer: ") + name +
                                     " y-shape does not match interface count"
                                    );
        }
        if (field.Nz() != nz) {
            throw std::runtime_error(
                                     std::string("InitialConditionInitializer: ") + name +
                                     " z-shape does not match interface count"
                                    );
        }
    };

    validate(ic.rho, "rho");
    validate(ic.u, "u");
    validate(ic.v, "v");
    validate(ic.w, "w");
    validate(ic.p, "p");

    if (ic.reactant_mass_fraction.has_value()) {
        validate(*ic.reactant_mass_fraction, "reactant_mass_fraction");
    }
}
