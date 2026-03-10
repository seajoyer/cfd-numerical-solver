#include "data/Mesh.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

#include "bc/BoundaryCondition.hpp"
#include "data/Variables.hpp"
#include "data/geometry/GeometryPrimitive.hpp"

Mesh::Mesh(const int nx, const int ny, const int nz, const int padding, const int dim)
    : nx_(nx), ny_(ny), nz_(nz), ng_(padding), dim_(dim) {
    Validate();

    if (dim_ < 2) {
        ny_ = 1;
    }
    if (dim_ < 3) {
        nz_ = 1;
    }

    RecomputeSizes();

    global_nx_ = nx_;
    global_ny_ = ny_;
    global_nz_ = nz_;

    Allocate();
}

int Mesh::GetDim() const {
    return dim_;
}

int Mesh::GetPadding() const {
    return ng_;
}

int Mesh::GetNx() const {
    return nx_;
}

int Mesh::GetNy() const {
    return ny_;
}

int Mesh::GetNz() const {
    return nz_;
}

int Mesh::GetSx() const {
    return sx_;
}

int Mesh::GetSy() const {
    return sy_;
}

int Mesh::GetSz() const {
    return sz_;
}

int Mesh::GetCoreStartX() const {
    return ng_;
}

int Mesh::GetCoreStartY() const {
    return dim_ >= 2 ? ng_ : 0;
}

int Mesh::GetCoreStartZ() const {
    return dim_ >= 3 ? ng_ : 0;
}

int Mesh::GetCoreEndExclusiveX() const {
    return ng_ + nx_;
}

int Mesh::GetCoreEndExclusiveY() const {
    return dim_ >= 2 ? ng_ + ny_ : 1;
}

int Mesh::GetCoreEndExclusiveZ() const {
    return dim_ >= 3 ? ng_ + nz_ : 1;
}

int Mesh::GetCoreNx() const {
    return nx_;
}

int Mesh::GetCoreNy() const {
    return ny_;
}

int Mesh::GetCoreNz() const {
    return nz_;
}

void Mesh::SetGlobalDecomposition(const int global_nx,
                                  const int global_ny,
                                  const int global_nz,
                                  const int offset_x,
                                  const int offset_y,
                                  const int offset_z) {
    if (global_nx <= 0) {
        throw std::invalid_argument("global_nx must be > 0");
    }
    if (dim_ >= 2 && global_ny <= 0) {
        throw std::invalid_argument("global_ny must be > 0 for dim>=2");
    }
    if (dim_ >= 3 && global_nz <= 0) {
        throw std::invalid_argument("global_nz must be > 0 for dim>=3");
    }
    if (offset_x < 0 || offset_y < 0 || offset_z < 0) {
        throw std::invalid_argument("global offsets must be >= 0");
    }

    global_nx_ = global_nx;
    global_ny_ = dim_ >= 2 ? global_ny : 1;
    global_nz_ = dim_ >= 3 ? global_nz : 1;

    offset_x_ = offset_x;
    offset_y_ = dim_ >= 2 ? offset_y : 0;
    offset_z_ = dim_ >= 3 ? offset_z : 0;
}

int Mesh::GetGlobalNx() const {
    return global_nx_;
}

int Mesh::GetGlobalNy() const {
    return global_ny_;
}

int Mesh::GetGlobalNz() const {
    return global_nz_;
}

int Mesh::GetOffsetX() const {
    return offset_x_;
}

int Mesh::GetOffsetY() const {
    return offset_y_;
}

int Mesh::GetOffsetZ() const {
    return offset_z_;
}

void Mesh::SetNeighborRank(const Axis axis, const Side side, const int rank) {
    ValidateAxisActive(axis);

    const auto a = static_cast<std::size_t>(axis);
    const auto s = static_cast<std::size_t>(side);
    neighbor_rank_[a][s] = rank;
}

int Mesh::GetNeighborRank(const Axis axis, const Side side) const {
    ValidateAxisActive(axis);

    const auto a = static_cast<std::size_t>(axis);
    const auto s = static_cast<std::size_t>(side);
    return neighbor_rank_[a][s];
}

bool Mesh::HasNeighborRank(const Axis axis, const Side side) const {
    return GetNeighborRank(axis, side) >= 0;
}

xt::xtensor<double, 1>& Mesh::Xb() {
    return xb_;
}

const xt::xtensor<double, 1>& Mesh::Xb() const {
    return xb_;
}

xt::xtensor<double, 1>& Mesh::Xc() {
    return xc_;
}

const xt::xtensor<double, 1>& Mesh::Xc() const {
    return xc_;
}

xt::xtensor<double, 1>& Mesh::Yb() {
    return yb_;
}

const xt::xtensor<double, 1>& Mesh::Yb() const {
    return yb_;
}

xt::xtensor<double, 1>& Mesh::Yc() {
    return yc_;
}

const xt::xtensor<double, 1>& Mesh::Yc() const {
    return yc_;
}

xt::xtensor<double, 1>& Mesh::Zb() {
    return zb_;
}

const xt::xtensor<double, 1>& Mesh::Zb() const {
    return zb_;
}

xt::xtensor<double, 1>& Mesh::Zc() {
    return zc_;
}

const xt::xtensor<double, 1>& Mesh::Zc() const {
    return zc_;
}

const xt::xtensor<double, 1>& Mesh::Dx() const {
    return dx_;
}

const xt::xtensor<double, 1>& Mesh::Dy() const {
    return dy_;
}

const xt::xtensor<double, 1>& Mesh::Dz() const {
    return dz_;
}

const xt::xtensor<double, 1>& Mesh::InvDx() const {
    return inv_dx_;
}

const xt::xtensor<double, 1>& Mesh::InvDy() const {
    return inv_dy_;
}

const xt::xtensor<double, 1>& Mesh::InvDz() const {
    return inv_dz_;
}

void Mesh::UpdateMetricsFromCoordinates() {
    for (int i = 0; i < sx_; ++i) {
        const double xl = xb_(static_cast<std::size_t>(i));
        const double xr = xb_(static_cast<std::size_t>(i + 1));
        const double d = xr - xl;

        xc_(static_cast<std::size_t>(i)) = 0.5 * (xl + xr);
        dx_(static_cast<std::size_t>(i)) = d;
        inv_dx_(static_cast<std::size_t>(i)) = d != 0.0 ? 1.0 / d : 0.0;
    }

    for (int j = 0; j < sy_; ++j) {
        const double yl = yb_(static_cast<std::size_t>(j));
        const double yr = yb_(static_cast<std::size_t>(j + 1));
        const double d = yr - yl;

        yc_(static_cast<std::size_t>(j)) = 0.5 * (yl + yr);
        dy_(static_cast<std::size_t>(j)) = d;
        inv_dy_(static_cast<std::size_t>(j)) = d != 0.0 ? 1.0 / d : 0.0;
    }

    for (int k = 0; k < sz_; ++k) {
        const double zl = zb_(static_cast<std::size_t>(k));
        const double zr = zb_(static_cast<std::size_t>(k + 1));
        const double d = zr - zl;

        zc_(static_cast<std::size_t>(k)) = 0.5 * (zl + zr);
        dz_(static_cast<std::size_t>(k)) = d;
        inv_dz_(static_cast<std::size_t>(k)) = d != 0.0 ? 1.0 / d : 0.0;
    }
}

bool Mesh::IsGlobalBoundary(const Axis axis, const Side side) const {
    const auto a = static_cast<std::size_t>(axis);
    const auto s = static_cast<std::size_t>(side);

    if (a == 1 && dim_ < 2) {
        return false;
    }
    if (a == 2 && dim_ < 3) {
        return false;
    }

    return is_global_boundary_[a][s];
}

void Mesh::SetGlobalBoundary(const Axis axis, const Side side, const bool is_global) {
    ValidateAxisActive(axis);

    const auto a = static_cast<std::size_t>(axis);
    const auto s = static_cast<std::size_t>(side);
    is_global_boundary_[a][s] = is_global;
}

void Mesh::SetAllGlobalBoundaries(const bool is_global) {
    for (int a = 0; a < 3; ++a) {
        is_global_boundary_[a][0] = is_global;
        is_global_boundary_[a][1] = is_global;
    }
}

xt::xtensor<std::uint8_t, 3>& Mesh::CellTypes() {
    return cell_type_;
}

const xt::xtensor<std::uint8_t, 3>& Mesh::CellTypes() const {
    return cell_type_;
}

void Mesh::SetAllCellsFluid() {
    cell_type_.fill(static_cast<std::uint8_t>(CellType::Fluid));
}

void Mesh::SetAllCellsSolid() {
    cell_type_.fill(static_cast<std::uint8_t>(CellType::Solid));
}

void Mesh::SetCellType(const int i, const int j, const int k, const CellType type) {
    ValidateCellIndex(i, j, k);
    cell_type_(static_cast<std::size_t>(i),
               static_cast<std::size_t>(j),
               static_cast<std::size_t>(k)) = static_cast<std::uint8_t>(type);
}

CellType Mesh::GetCellType(const int i, const int j, const int k) const {
    ValidateCellIndex(i, j, k);
    return static_cast<CellType>(cell_type_(static_cast<std::size_t>(i),
                                            static_cast<std::size_t>(j),
                                            static_cast<std::size_t>(k)));
}

bool Mesh::IsFluidCell(const int i, const int j, const int k) const {
    return GetCellType(i, j, k) == CellType::Fluid;
}

bool Mesh::IsSolidCell(const int i, const int j, const int k) const {
    return GetCellType(i, j, k) == CellType::Solid;
}

void Mesh::AddPrimitive(std::shared_ptr<GeometryPrimitive> primitive) {
    if (!primitive) {
        throw std::invalid_argument("Mesh::AddPrimitive: primitive is null");
    }
    if (primitive->GetDimension() > dim_) {
        throw std::invalid_argument("Mesh::AddPrimitive: primitive dimension exceeds mesh dimension");
    }
    primitives_.push_back(std::move(primitive));
}

void Mesh::ClearPrimitives() {
    primitives_.clear();
}

const std::vector<std::shared_ptr<GeometryPrimitive>>& Mesh::Primitives() const {
    return primitives_;
}

void Mesh::BuildCellTypesFromPrimitives() {
    SetAllCellsFluid();

    if (primitives_.empty()) {
        return;
    }

    for (int k = 0; k < sz_; ++k) {
        const double z = zc_(static_cast<std::size_t>(k));
        for (int j = 0; j < sy_; ++j) {
            const double y = yc_(static_cast<std::size_t>(j));
            for (int i = 0; i < sx_; ++i) {
                const double x = xc_(static_cast<std::size_t>(i));

                GeometryPoint p{x, y, z};
                for (const auto& primitive : primitives_) {
                    if (primitive->Contains(p)) {
                        SetCellType(i, j, k, CellType::Solid);
                        break;
                    }
                }
            }
        }
    }
}

void Mesh::BuildImmersedFaces() {
    ClearAllImmersedFaces();

    if (primitives_.empty()) {
        return;
    }

    for (int k = 0; k < sz_; ++k) {
        const double zc = zc_(static_cast<std::size_t>(k));
        for (int j = 0; j < sy_; ++j) {
            const double yc = yc_(static_cast<std::size_t>(j));
            for (int i = 0; i < sx_; ++i) {
                if (!IsFluidCell(i, j, k)) {
                    continue;
                }

                const double xc = xc_(static_cast<std::size_t>(i));
                const GeometryPoint cell_center{xc, yc, zc};

                for (int axis_index = 0; axis_index < 3; ++axis_index) {
                    if (axis_index == 1 && dim_ < 2) {
                        continue;
                    }
                    if (axis_index == 2 && dim_ < 3) {
                        continue;
                    }

                    const Axis axis = axis_index == 0 ? Axis::X : (axis_index == 1 ? Axis::Y : Axis::Z);
                    const AxisStride stride = AxisStride::FromAxis(axis);

                    for (int side_index = 0; side_index < 2; ++side_index) {
                        const Side side = side_index == 0 ? Side::Left : Side::Right;

                        const int sign = side == Side::Left ? -1 : 1;
                        const int in = i + sign * stride.di;
                        const int jn = j + sign * stride.dj;
                        const int kn = k + sign * stride.dk;

                        if (in < 0 || in >= sx_ || jn < 0 || jn >= sy_ || kn < 0 || kn >= sz_) {
                            continue;
                        }

                        if (!IsSolidCell(in, jn, kn)) {
                            continue;
                        }

                        const double xn = xc_(static_cast<std::size_t>(in));
                        const double yn = yc_(static_cast<std::size_t>(jn));
                        const double zn = zc_(static_cast<std::size_t>(kn));

                        const GeometryPoint face_mid{
                            0.5 * (xc + xn),
                            0.5 * (yc + yn),
                            0.5 * (zc + zn)
                        };

                        const GeometryPrimitive* primitive = FindClosestPrimitive(
                            face_mid.x, face_mid.y, face_mid.z);

                        if (!primitive) {
                            continue;
                        }

                        const GeometryVector normal = primitive->OutwardNormal(face_mid);
                        const double distance = std::max(primitive->SignedDistance(cell_center), 0.0);

                        SetImmersedFace(axis, side, i, j, k,
                                        normal.x, normal.y, normal.z,
                                        distance);
                    }
                }
            }
        }
    }
}

void Mesh::ClearAllImmersedFaces() {
    immersed_face_active_.fill(0);
    immersed_face_normal_.fill(0.0);
    immersed_face_distance_.fill(0.0);
}

void Mesh::ClearImmersedFace(const Axis axis, const Side side, const int i, const int j, const int k) {
    ValidateAxisActive(axis);
    ValidateCellIndex(i, j, k);

    const auto a = static_cast<std::size_t>(AxisIndex(axis));
    const auto s = static_cast<std::size_t>(SideIndex(side));
    const auto I = static_cast<std::size_t>(i);
    const auto J = static_cast<std::size_t>(j);
    const auto K = static_cast<std::size_t>(k);

    immersed_face_active_(a, s, I, J, K) = 0;
    immersed_face_normal_(a, s, 0, I, J, K) = 0.0;
    immersed_face_normal_(a, s, 1, I, J, K) = 0.0;
    immersed_face_normal_(a, s, 2, I, J, K) = 0.0;
    immersed_face_distance_(a, s, I, J, K) = 0.0;
}

void Mesh::SetImmersedFace(const Axis axis,
                           const Side side,
                           const int i,
                           const int j,
                           const int k,
                           const double normal_x,
                           const double normal_y,
                           const double normal_z,
                           const double distance) {
    ValidateAxisActive(axis);
    ValidateCellIndex(i, j, k);

    const auto a = static_cast<std::size_t>(AxisIndex(axis));
    const auto s = static_cast<std::size_t>(SideIndex(side));
    const auto I = static_cast<std::size_t>(i);
    const auto J = static_cast<std::size_t>(j);
    const auto K = static_cast<std::size_t>(k);

    immersed_face_active_(a, s, I, J, K) = 1;
    immersed_face_normal_(a, s, 0, I, J, K) = normal_x;
    immersed_face_normal_(a, s, 1, I, J, K) = normal_y;
    immersed_face_normal_(a, s, 2, I, J, K) = normal_z;
    immersed_face_distance_(a, s, I, J, K) = distance;
}

bool Mesh::HasImmersedFace(const Axis axis, const Side side, const int i, const int j, const int k) const {
    ValidateAxisActive(axis);
    ValidateCellIndex(i, j, k);

    const auto a = static_cast<std::size_t>(AxisIndex(axis));
    const auto s = static_cast<std::size_t>(SideIndex(side));
    return immersed_face_active_(a, s,
                                 static_cast<std::size_t>(i),
                                 static_cast<std::size_t>(j),
                                 static_cast<std::size_t>(k)) != 0;
}

ImmersedFaceInfo Mesh::GetImmersedFaceInfo(const Axis axis,
                                           const Side side,
                                           const int i,
                                           const int j,
                                           const int k) const {
    ValidateAxisActive(axis);
    ValidateCellIndex(i, j, k);

    const auto a = static_cast<std::size_t>(AxisIndex(axis));
    const auto s = static_cast<std::size_t>(SideIndex(side));
    const auto I = static_cast<std::size_t>(i);
    const auto J = static_cast<std::size_t>(j);
    const auto K = static_cast<std::size_t>(k);

    ImmersedFaceInfo info;
    info.is_active = immersed_face_active_(a, s, I, J, K) != 0;
    info.normal_x = immersed_face_normal_(a, s, 0, I, J, K);
    info.normal_y = immersed_face_normal_(a, s, 1, I, J, K);
    info.normal_z = immersed_face_normal_(a, s, 2, I, J, K);
    info.distance = immersed_face_distance_(a, s, I, J, K);
    return info;
}

bool Mesh::HasAnyImmersedFace(const int i, const int j, const int k) const {
    ValidateCellIndex(i, j, k);

    const auto I = static_cast<std::size_t>(i);
    const auto J = static_cast<std::size_t>(j);
    const auto K = static_cast<std::size_t>(k);

    for (int a = 0; a < 3; ++a) {
        if (a == 1 && dim_ < 2) {
            continue;
        }
        if (a == 2 && dim_ < 3) {
            continue;
        }
        for (int s = 0; s < 2; ++s) {
            if (immersed_face_active_(static_cast<std::size_t>(a),
                                      static_cast<std::size_t>(s),
                                      I, J, K) != 0) {
                return true;
            }
        }
    }

    return false;
}

void Mesh::Validate() const {
    if (nx_ <= 0) {
        throw std::invalid_argument("nx must be > 0");
    }
    if (dim_ < 1 || dim_ > 3) {
        throw std::invalid_argument("dim must be 1..3");
    }
    if (ng_ < 0) {
        throw std::invalid_argument("padding must be >= 0");
    }
    if (dim_ >= 2 && ny_ <= 0) {
        throw std::invalid_argument("ny must be > 0 for dim>=2");
    }
    if (dim_ >= 3 && nz_ <= 0) {
        throw std::invalid_argument("nz must be > 0 for dim>=3");
    }
}

void Mesh::RecomputeSizes() {
    sx_ = nx_ + 2 * ng_;
    sy_ = dim_ >= 2 ? ny_ + 2 * ng_ : 1;
    sz_ = dim_ >= 3 ? nz_ + 2 * ng_ : 1;
}

void Mesh::Allocate() {
    xb_ = xt::zeros<double>({static_cast<std::size_t>(sx_ + 1)});
    xc_ = xt::zeros<double>({static_cast<std::size_t>(sx_)});
    dx_ = xt::zeros<double>({static_cast<std::size_t>(sx_)});
    inv_dx_ = xt::zeros<double>({static_cast<std::size_t>(sx_)});

    yb_ = xt::zeros<double>({static_cast<std::size_t>(sy_ + 1)});
    yc_ = xt::zeros<double>({static_cast<std::size_t>(sy_)});
    dy_ = xt::zeros<double>({static_cast<std::size_t>(sy_)});
    inv_dy_ = xt::zeros<double>({static_cast<std::size_t>(sy_)});

    zb_ = xt::zeros<double>({static_cast<std::size_t>(sz_ + 1)});
    zc_ = xt::zeros<double>({static_cast<std::size_t>(sz_)});
    dz_ = xt::zeros<double>({static_cast<std::size_t>(sz_)});
    inv_dz_ = xt::zeros<double>({static_cast<std::size_t>(sz_)});

    cell_type_ = xt::zeros<std::uint8_t>({
        static_cast<std::size_t>(sx_),
        static_cast<std::size_t>(sy_),
        static_cast<std::size_t>(sz_)
    });

    immersed_face_active_ = xt::zeros<std::uint8_t>({
        std::size_t{3},
        std::size_t{2},
        static_cast<std::size_t>(sx_),
        static_cast<std::size_t>(sy_),
        static_cast<std::size_t>(sz_)
    });

    immersed_face_normal_ = xt::zeros<double>({
        std::size_t{3},
        std::size_t{2},
        std::size_t{3},
        static_cast<std::size_t>(sx_),
        static_cast<std::size_t>(sy_),
        static_cast<std::size_t>(sz_)
    });

    immersed_face_distance_ = xt::zeros<double>({
        std::size_t{3},
        std::size_t{2},
        static_cast<std::size_t>(sx_),
        static_cast<std::size_t>(sy_),
        static_cast<std::size_t>(sz_)
    });

    SetAllCellsFluid();
    ClearAllImmersedFaces();
}

void Mesh::ValidateCellIndex(const int i, const int j, const int k) const {
    if (i < 0 || i >= sx_) {
        throw std::out_of_range("Mesh: i is out of range");
    }
    if (j < 0 || j >= sy_) {
        throw std::out_of_range("Mesh: j is out of range");
    }
    if (k < 0 || k >= sz_) {
        throw std::out_of_range("Mesh: k is out of range");
    }
}

void Mesh::ValidateAxisActive(const Axis axis) const {
    const auto a = static_cast<std::size_t>(axis);

    if (a == 1 && dim_ < 2) {
        throw std::invalid_argument("Mesh: Axis Y is not active for dim<2");
    }
    if (a == 2 && dim_ < 3) {
        throw std::invalid_argument("Mesh: Axis Z is not active for dim<3");
    }
}

int Mesh::AxisIndex(const Axis axis) const {
    return static_cast<int>(static_cast<std::uint8_t>(axis));
}

int Mesh::SideIndex(const Side side) const {
    return static_cast<int>(static_cast<std::uint8_t>(side));
}

const GeometryPrimitive* Mesh::FindClosestPrimitive(const double x, const double y, const double z) const {
    if (primitives_.empty()) {
        return nullptr;
    }

    const GeometryPoint p{x, y, z};

    const GeometryPrimitive* best = nullptr;
    double best_abs_distance = std::numeric_limits<double>::infinity();

    for (const auto& primitive : primitives_) {
        const double d = std::abs(primitive->SignedDistance(p));
        if (d < best_abs_distance) {
            best_abs_distance = d;
            best = primitive.get();
        }
    }

    return best;
}
