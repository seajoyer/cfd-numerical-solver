#include "parallel/DomainDecomposition.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "geometry/Cell.hpp"
#include "geometry/Face.hpp"
#include "geometry/Mesh.hpp"
#include "geometry/Node.hpp"
#include "parallel/MPIContext.hpp"

namespace {
    struct HaloAccumulator final {
        int remote_rank = -1;

        std::vector<std::size_t> send_local_ids;
        std::vector<std::size_t> recv_local_ids;

        std::unordered_set<std::size_t> send_seen;
        std::unordered_set<std::size_t> recv_seen;
    };

    [[nodiscard]] std::vector<std::size_t> MakeAllCellIds(const Mesh& mesh) {
        std::vector<std::size_t> ids(mesh.GetCellCount());
        for (std::size_t i = 0; i < ids.size(); ++i) {
            ids[i] = i;
        }
        return ids;
    }

    [[nodiscard]] double CellCoordByAxis(const Cell& cell, const int axis) {
        if (axis == 0) {
            return cell.center_x;
        }
        if (axis == 1) {
            return cell.center_y;
        }
        return cell.center_z;
    }

    [[nodiscard]] int ChooseSplitAxis(const Mesh& mesh,
                                      const std::vector<std::size_t>& cell_ids) {
        double min_x = 0.0, max_x = 0.0;
        double min_y = 0.0, max_y = 0.0;
        double min_z = 0.0, max_z = 0.0;

        bool first = true;
        for (const std::size_t cell_id : cell_ids) {
            const Cell& cell = mesh.GetCell(cell_id);

            if (first) {
                min_x = max_x = cell.center_x;
                min_y = max_y = cell.center_y;
                min_z = max_z = cell.center_z;
                first = false;
                continue;
            }

            min_x = std::min(min_x, cell.center_x);
            max_x = std::max(max_x, cell.center_x);

            min_y = std::min(min_y, cell.center_y);
            max_y = std::max(max_y, cell.center_y);

            min_z = std::min(min_z, cell.center_z);
            max_z = std::max(max_z, cell.center_z);
        }

        const double span_x = max_x - min_x;
        const double span_y = (mesh.GetDim() >= 2) ? (max_y - min_y) : -1.0;
        const double span_z = (mesh.GetDim() >= 3) ? (max_z - min_z) : -1.0;

        if (mesh.GetDim() == 1) {
            return 0;
        }

        if (mesh.GetDim() == 2) {
            return (span_x >= span_y) ? 0 : 1;
        }

        if (span_x >= span_y && span_x >= span_z) {
            return 0;
        }
        if (span_y >= span_z) {
            return 1;
        }
        return 2;
    }

    void BuildRCBRecursive(const Mesh& mesh,
                           const std::vector<std::size_t>& cell_ids,
                           const int proc_begin,
                           const int proc_end,
                           std::vector<int>& part) {
        const int proc_count = proc_end - proc_begin;
        if (proc_count <= 0) {
            throw std::runtime_error("DomainDecomposition: invalid processor interval");
        }

        if (proc_count == 1) {
            for (const std::size_t cell_id : cell_ids) {
                part[cell_id] = proc_begin;
            }
            return;
        }

        if (cell_ids.empty()) {
            return;
        }

        const int axis = ChooseSplitAxis(mesh, cell_ids);

        std::vector<std::size_t> sorted = cell_ids;
        std::stable_sort(sorted.begin(), sorted.end(),
                         [&](const std::size_t lhs, const std::size_t rhs) {
                             const Cell& a = mesh.GetCell(lhs);
                             const Cell& b = mesh.GetCell(rhs);
                             const double ca = CellCoordByAxis(a, axis);
                             const double cb = CellCoordByAxis(b, axis);
                             if (ca != cb) {
                                 return ca < cb;
                             }
                             return lhs < rhs;
                         });

        const int left_proc_count = proc_count / 2;
        const int right_proc_count = proc_count - left_proc_count;

        const std::size_t n = sorted.size();
        const std::size_t cut =
            static_cast<std::size_t>(std::llround(
                static_cast<double>(n) * static_cast<double>(left_proc_count) /
                static_cast<double>(proc_count)
            ));

        const std::vector<std::size_t> left(
            sorted.begin(),
            sorted.begin() + static_cast<std::ptrdiff_t>(cut)
        );
        const std::vector<std::size_t> right(
            sorted.begin() + static_cast<std::ptrdiff_t>(cut),
            sorted.end()
        );

        BuildRCBRecursive(mesh, left, proc_begin, proc_begin + left_proc_count, part);
        BuildRCBRecursive(mesh, right, proc_begin + left_proc_count, proc_end, part);
    }

    [[nodiscard]] std::vector<int> BuildRCBPartition(const Mesh& global_mesh,
                                                     const int nproc) {
        if (nproc <= 0) {
            throw std::runtime_error("DomainDecomposition: MPI size must be positive");
        }

        std::vector<int> part(global_mesh.GetCellCount(), -1);
        const std::vector<std::size_t> all_ids = MakeAllCellIds(global_mesh);

        BuildRCBRecursive(global_mesh, all_ids, 0, nproc, part);

        for (std::size_t i = 0; i < part.size(); ++i) {
            if (part[i] < 0 || part[i] >= nproc) {
                throw std::runtime_error("DomainDecomposition: invalid partition result");
            }
        }

        return part;
    }

    [[nodiscard]] std::vector<std::size_t> CollectOwnedGlobalIds(const Mesh& global_mesh,
                                                                 const std::vector<int>& part,
                                                                 const int rank) {
        std::vector<std::size_t> owned;
        owned.reserve(global_mesh.GetCellCount());

        for (std::size_t global_cell_id = 0; global_cell_id < global_mesh.GetCellCount(); ++global_cell_id) {
            if (part[global_cell_id] == rank) {
                owned.push_back(global_cell_id);
            }
        }

        return owned;
    }

    [[nodiscard]] std::set<std::size_t> CollectGhostGlobalIds(const Mesh& global_mesh,
                                                              const std::vector<int>& part,
                                                              const int rank) {
        std::set<std::size_t> ghosts;

        for (const Face& face : global_mesh.Faces()) {
            if (!face.IsInternal()) {
                continue;
            }

            const int owner_rank = part[face.owner_cell_id];
            const int neighbor_rank = part[face.neighbor_cell_id];

            if (owner_rank == rank && neighbor_rank != rank) {
                ghosts.insert(face.neighbor_cell_id);
            }
            else if (neighbor_rank == rank && owner_rank != rank) {
                ghosts.insert(face.owner_cell_id);
            }
        }

        return ghosts;
    }

    void CopyLocalCells(const Mesh& global_mesh,
                        const std::vector<std::size_t>& owned_global_ids,
                        const std::vector<std::size_t>& ghost_global_ids,
                        DomainDecomposition::Result& result) {
        auto& local_cells = result.local_mesh.Cells();

        result.local_to_global_cell.clear();
        result.global_to_local_cell.clear();

        local_cells.reserve(owned_global_ids.size() + ghost_global_ids.size());
        result.local_to_global_cell.reserve(owned_global_ids.size() + ghost_global_ids.size());

        auto copy_one = [&](const std::size_t global_cell_id) {
            Cell cell = global_mesh.GetCell(global_cell_id);
            cell.local_id = local_cells.size();
            cell.face_ids.clear();

            result.global_to_local_cell[global_cell_id] = cell.local_id;
            result.local_to_global_cell.push_back(global_cell_id);

            local_cells.push_back(std::move(cell));
        };

        for (const std::size_t global_cell_id : owned_global_ids) {
            copy_one(global_cell_id);
        }
        for (const std::size_t global_cell_id : ghost_global_ids) {
            copy_one(global_cell_id);
        }

        result.n_owned_cells = owned_global_ids.size();
        result.n_ghost_cells = ghost_global_ids.size();

        result.local_mesh.SetOwnedCellCount(result.n_owned_cells);
        result.local_mesh.SetGhostCellCount(result.n_ghost_cells);
    }

    void CopyLocalNodesAndRemapConnectivity(DomainDecomposition::Result& result,
                                            const Mesh& global_mesh) {
        std::set<std::size_t> used_global_node_ids;

        for (const Cell& local_cell : result.local_mesh.Cells()) {
            const std::size_t global_cell_id = result.local_to_global_cell[local_cell.local_id];
            const Cell& global_cell = global_mesh.GetCell(global_cell_id);

            for (const std::size_t global_node_id : global_cell.node_ids) {
                used_global_node_ids.insert(global_node_id);
            }
        }

        std::unordered_map<std::size_t, std::size_t> global_to_local_node;
        auto& local_nodes = result.local_mesh.Nodes();
        local_nodes.clear();
        local_nodes.reserve(used_global_node_ids.size());

        for (const std::size_t global_node_id : used_global_node_ids) {
            Node node = global_mesh.GetNode(global_node_id);
            node.id = local_nodes.size();
            global_to_local_node[global_node_id] = node.id;
            local_nodes.push_back(std::move(node));
        }

        for (Cell& local_cell : result.local_mesh.Cells()) {
            const std::size_t global_cell_id = result.local_to_global_cell[local_cell.local_id];
            const Cell& global_cell = global_mesh.GetCell(global_cell_id);

            local_cell.node_ids.clear();
            local_cell.node_ids.reserve(global_cell.node_ids.size());

            for (const std::size_t global_node_id : global_cell.node_ids) {
                local_cell.node_ids.push_back(global_to_local_node.at(global_node_id));
            }
        }
    }

    void AppendHaloId(std::vector<std::size_t>& ids,
                      std::unordered_set<std::size_t>& seen,
                      const std::size_t value) {
        if (seen.insert(value).second) {
            ids.push_back(value);
        }
    }

    void BuildLocalFacesAndHalos(const Mesh& global_mesh,
                                 const std::vector<int>& part,
                                 const int rank,
                                 DomainDecomposition::Result& result) {
        auto& local_faces = result.local_mesh.Faces();
        auto& local_cells = result.local_mesh.Cells();

        std::map<int, HaloAccumulator> halos_by_rank;

        std::unordered_map<std::size_t, std::size_t> global_to_local_node;

        for (const Cell& local_cell : result.local_mesh.Cells()) {
            const std::size_t global_cell_id = result.local_to_global_cell[local_cell.local_id];
            const Cell& global_cell = global_mesh.GetCell(global_cell_id);

            for (std::size_t k = 0; k < global_cell.node_ids.size(); ++k) {
                global_to_local_node[global_cell.node_ids[k]] = local_cell.node_ids[k];
            }
        }

        auto remap_face_nodes_to_local = [&](Face& face) {
            for (std::size_t& node_id : face.node_ids) {
                node_id = global_to_local_node.at(node_id);
            }
        };

        auto add_local_face = [&](Face face) {
            face.id = local_faces.size();
            local_faces.push_back(face);

            local_cells[face.owner_cell_id].face_ids.push_back(face.id);
            if (face.IsInternal() || face.IsMPIBoundary()) {
                local_cells[face.neighbor_cell_id].face_ids.push_back(face.id);
            }
        };

        for (const Face& global_face : global_mesh.Faces()) {
            if (global_face.IsPhysicalBoundary()) {
                const int owner_rank = part[global_face.owner_cell_id];
                if (owner_rank != rank) {
                    continue;
                }

                Face local_face = global_face;
                local_face.owner_cell_id = result.global_to_local_cell.at(global_face.owner_cell_id);
                local_face.neighbor_cell_id = Face::k_invalid_cell_id;
                local_face.kind = FaceKind::PhysicalBoundary;
                local_face.remote_rank = -1;
                local_face.remote_cell_id = Face::k_invalid_cell_id;

                remap_face_nodes_to_local(local_face);
                add_local_face(std::move(local_face));
                continue;
            }

            if (!global_face.IsInternal()) {
                continue;
            }

            const std::size_t global_owner = global_face.owner_cell_id;
            const std::size_t global_neighbor = global_face.neighbor_cell_id;

            const int owner_rank = part[global_owner];
            const int neighbor_rank = part[global_neighbor];

            if (owner_rank == rank && neighbor_rank == rank) {
                Face local_face = global_face;
                local_face.owner_cell_id = result.global_to_local_cell.at(global_owner);
                local_face.neighbor_cell_id = result.global_to_local_cell.at(global_neighbor);
                local_face.kind = FaceKind::Interior;
                local_face.remote_rank = -1;
                local_face.remote_cell_id = Face::k_invalid_cell_id;

                remap_face_nodes_to_local(local_face);
                add_local_face(std::move(local_face));
                continue;
            }

            if (owner_rank == rank && neighbor_rank != rank) {
                Face local_face = global_face;
                local_face.owner_cell_id = result.global_to_local_cell.at(global_owner);
                local_face.neighbor_cell_id = result.global_to_local_cell.at(global_neighbor);
                local_face.kind = FaceKind::MPIBoundary;
                local_face.boundary_tag = -1;
                local_face.remote_rank = neighbor_rank;
                local_face.remote_cell_id = global_neighbor;

                remap_face_nodes_to_local(local_face);
                add_local_face(local_face);

                HaloAccumulator& halo = halos_by_rank[neighbor_rank];
                halo.remote_rank = neighbor_rank;
                AppendHaloId(halo.send_local_ids, halo.send_seen, local_face.owner_cell_id);
                AppendHaloId(halo.recv_local_ids, halo.recv_seen, local_face.neighbor_cell_id);
                continue;
            }

            if (neighbor_rank == rank && owner_rank != rank) {
                Face local_face = global_face;
                local_face.owner_cell_id = result.global_to_local_cell.at(global_neighbor);
                local_face.neighbor_cell_id = result.global_to_local_cell.at(global_owner);
                local_face.kind = FaceKind::MPIBoundary;
                local_face.boundary_tag = -1;
                local_face.remote_rank = owner_rank;
                local_face.remote_cell_id = global_owner;

                local_face.normal_x *= -1.0;
                local_face.normal_y *= -1.0;
                local_face.normal_z *= -1.0;

                remap_face_nodes_to_local(local_face);
                add_local_face(local_face);

                HaloAccumulator& halo = halos_by_rank[owner_rank];
                halo.remote_rank = owner_rank;
                AppendHaloId(halo.send_local_ids, halo.send_seen, local_face.owner_cell_id);
                AppendHaloId(halo.recv_local_ids, halo.recv_seen, local_face.neighbor_cell_id);
                continue;
            }
        }

        result.halos.clear();
        result.halos.reserve(halos_by_rank.size());

        for (auto& [remote_rank, halo_acc] : halos_by_rank) {
            (void)remote_rank;
            DomainDecomposition::NeighborHalo halo;
            halo.remote_rank = halo_acc.remote_rank;
            halo.send_local_ids = std::move(halo_acc.send_local_ids);
            halo.recv_local_ids = std::move(halo_acc.recv_local_ids);
            result.halos.push_back(std::move(halo));
        }

        std::sort(result.halos.begin(), result.halos.end(),
                  [](const DomainDecomposition::NeighborHalo& a,
                     const DomainDecomposition::NeighborHalo& b) {
                      return a.remote_rank < b.remote_rank;
                  });
    }
} // namespace


std::vector<int> DomainDecomposition::BuildRCBPartition(const Mesh& global_mesh,
                                                        const int nproc) {
    return ::BuildRCBPartition(global_mesh, nproc);
}

DomainDecomposition::Result DomainDecomposition::BuildLocalResultForRank(
    const Mesh& global_mesh,
    const std::vector<int>& global_part,
    const int rank
) {
    Result result;
    result.local_mesh = Mesh(global_mesh.GetDim());
    result.global_part = global_part;
    result.owned_global_ids = CollectOwnedGlobalIds(global_mesh, result.global_part, rank);

    {
        const std::set<std::size_t> ghost_set =
            CollectGhostGlobalIds(global_mesh, result.global_part, rank);
        result.ghost_global_ids.assign(ghost_set.begin(), ghost_set.end());
    }

    CopyLocalCells(global_mesh,
                   result.owned_global_ids,
                   result.ghost_global_ids,
                   result);

    CopyLocalNodesAndRemapConnectivity(result, global_mesh);
    BuildLocalFacesAndHalos(global_mesh, result.global_part, rank, result);

    result.local_mesh.Validate();
    return result;
}

DomainDecomposition::Result DomainDecomposition::BuildRCB(const Mesh& global_mesh,
                                                          const MPIContext& mpi) {
    const std::vector<int> global_part = BuildRCBPartition(global_mesh, mpi.Size());
    return BuildLocalResultForRank(global_mesh, global_part, mpi.Rank());
}
