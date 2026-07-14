#include "tabi_calderon_aux_surface.h"

#include "tabi_analytic_ses.h"

#include "../titan_executor.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace TitanPBE::Tabi::calderon {
namespace {

// Intersections computed from adjacent tetrahedra can differ by a few ulps.
// Welding at 1e-5 A is far below any supported grid spacing (>= 0.25 A), but
// prevents duplicate seam vertices and removes solver-hostile micro-slivers.
constexpr double kQuantize = 1.0e5;
constexpr std::size_t kMaxGridVertices = 4000000;

struct Vec3 {
  double x{0.0};
  double y{0.0};
  double z{0.0};
};

struct CellKey {
  int x{0};
  int y{0};
  int z{0};

  bool operator==(const CellKey &) const = default;
};

struct CellKeyHash {
  std::size_t operator()(const CellKey &key) const noexcept {
    std::size_t hash = static_cast<std::size_t>(key.x) * 73856093U;
    hash ^= static_cast<std::size_t>(key.y) * 19349663U;
    hash ^= static_cast<std::size_t>(key.z) * 83492791U;
    return hash;
  }
};

struct VertexKey {
  std::int64_t x{0};
  std::int64_t y{0};
  std::int64_t z{0};

  bool operator==(const VertexKey &) const = default;
};

struct VertexKeyHash {
  std::size_t operator()(const VertexKey &key) const noexcept {
    std::size_t hash = static_cast<std::size_t>(key.x);
    hash ^= static_cast<std::size_t>(key.y) + 0x9e3779b9U + (hash << 6U) +
            (hash >> 2U);
    hash ^= static_cast<std::size_t>(key.z) + 0x9e3779b9U + (hash << 6U) +
            (hash >> 2U);
    return hash;
  }
};

struct EdgeKey {
  std::size_t lo{0};
  std::size_t hi{0};

  bool operator==(const EdgeKey &) const = default;
};

struct EdgeKeyHash {
  std::size_t operator()(const EdgeKey &key) const noexcept {
    std::size_t hash = key.lo;
    hash ^= key.hi + 0x9e3779b9U + (hash << 6U) + (hash >> 2U);
    return hash;
  }
};

struct FieldSample {
  double value{1.0};
  std::size_t atom_index{0};
};

struct AtomGrid {
  Vec3 origin;
  double cell_size{1.0};
  double max_radius{0.0};
  std::unordered_map<CellKey, std::vector<std::size_t>, CellKeyHash> cells;
};

[[nodiscard]] Vec3 operator+(const Vec3 &left, const Vec3 &right) noexcept {
  return {left.x + right.x, left.y + right.y, left.z + right.z};
}

[[nodiscard]] Vec3 operator-(const Vec3 &left, const Vec3 &right) noexcept {
  return {left.x - right.x, left.y - right.y, left.z - right.z};
}

[[nodiscard]] Vec3 operator*(double scale, const Vec3 &value) noexcept {
  return {scale * value.x, scale * value.y, scale * value.z};
}

[[nodiscard]] double dot(const Vec3 &left, const Vec3 &right) noexcept {
  return left.x * right.x + left.y * right.y + left.z * right.z;
}

[[nodiscard]] Vec3 cross(const Vec3 &left, const Vec3 &right) noexcept {
  return {left.y * right.z - left.z * right.y,
          left.z * right.x - left.x * right.z,
          left.x * right.y - left.y * right.x};
}

[[nodiscard]] double norm(const Vec3 &value) noexcept {
  return std::sqrt(dot(value, value));
}

[[nodiscard]] Vec3 normalized(const Vec3 &value) noexcept {
  const double length = norm(value);
  return length > 1.0e-300 ? (1.0 / length) * value : Vec3{1.0, 0.0, 0.0};
}

[[nodiscard]] CellKey cell_for(const AtomGrid &grid,
                               const Vec3 &point) noexcept {
  return {
      static_cast<int>(std::floor((point.x - grid.origin.x) / grid.cell_size)),
      static_cast<int>(std::floor((point.y - grid.origin.y) / grid.cell_size)),
      static_cast<int>(std::floor((point.z - grid.origin.z) / grid.cell_size))};
}

[[nodiscard]] AtomGrid build_atom_grid(std::span<const Atom> atoms,
                                       const Vec3 &origin) {
  AtomGrid grid;
  grid.origin = origin;
  for (const Atom &atom : atoms) {
    grid.max_radius = std::max(grid.max_radius, std::max(0.0, atom.radius));
  }
  grid.cell_size = std::max(1.0e-6, grid.max_radius);
  grid.cells.reserve(atoms.size() * 2);
  for (std::size_t index = 0; index < atoms.size(); ++index) {
    const Atom &atom = atoms[index];
    if (!(atom.radius > 0.0) || !std::isfinite(atom.radius)) {
      continue;
    }
    grid.cells[cell_for(grid, {atom.x, atom.y, atom.z})].push_back(index);
  }
  return grid;
}

[[nodiscard]] FieldSample sample_union_sdf(const AtomGrid &grid,
                                           std::span<const Atom> atoms,
                                           const Vec3 &point) noexcept {
  const CellKey center = cell_for(grid, point);
  FieldSample best;
  best.value = std::numeric_limits<double>::infinity();
  // A sphere that contains a point must have its centre in this 3x3x3 cell
  // neighborhood because the cell width equals the largest radius.
  for (int dz = -1; dz <= 1; ++dz) {
    for (int dy = -1; dy <= 1; ++dy) {
      for (int dx = -1; dx <= 1; ++dx) {
        const auto found =
            grid.cells.find({center.x + dx, center.y + dy, center.z + dz});
        if (found == grid.cells.end()) {
          continue;
        }
        for (const std::size_t index : found->second) {
          const Atom &atom = atoms[index];
          const double rx = point.x - atom.x;
          const double ry = point.y - atom.y;
          const double rz = point.z - atom.z;
          const double value =
              std::sqrt(rx * rx + ry * ry + rz * rz) - atom.radius;
          if (value < best.value) {
            best = {value, index};
          }
        }
      }
    }
  }
  if (!std::isfinite(best.value)) {
    best.value = grid.max_radius;
    best.atom_index = 0;
  }
  return best;
}

[[nodiscard]] bool point_inside_union(const AtomGrid &grid,
                                      std::span<const Atom> atoms,
                                      const Vec3 &point) noexcept {
  const CellKey center = cell_for(grid, point);
  for (int dz = -1; dz <= 1; ++dz) {
    for (int dy = -1; dy <= 1; ++dy) {
      for (int dx = -1; dx <= 1; ++dx) {
        const auto found =
            grid.cells.find({center.x + dx, center.y + dy, center.z + dz});
        if (found == grid.cells.end())
          continue;
        for (const std::size_t index : found->second) {
          const Atom &atom = atoms[index];
          const double rx = point.x - atom.x;
          const double ry = point.y - atom.y;
          const double rz = point.z - atom.z;
          if (rx * rx + ry * ry + rz * rz < atom.radius * atom.radius)
            return true;
        }
      }
    }
  }
  return false;
}

[[nodiscard]] const std::array<Vec3, 27> &topology_sample_directions() {
  // Construct once rather than normalizing the same 26 directions at every
  // grid vertex. The zero entry represents the unshifted probe centre.
  static const std::array<Vec3, 27> directions = [] {
    std::array<Vec3, 27> result{};
    std::size_t index = 0;
    for (int dz = -1; dz <= 1; ++dz) {
      for (int dy = -1; dy <= 1; ++dy) {
        for (int dx = -1; dx <= 1; ++dx) {
          if (dx == 0 && dy == 0 && dz == 0) {
            result[index++] = {};
          } else {
            result[index++] =
                normalized({static_cast<double>(dx), static_cast<double>(dy),
                            static_cast<double>(dz)});
          }
        }
      }
    }
    return result;
  }();
  return directions;
}

[[nodiscard]] VertexKey vertex_key(const Vec3 &point) noexcept {
  return {static_cast<std::int64_t>(std::llround(point.x * kQuantize)),
          static_cast<std::int64_t>(std::llround(point.y * kQuantize)),
          static_cast<std::int64_t>(std::llround(point.z * kQuantize))};
}

std::size_t
add_vertex(SurfaceMesh &mesh,
           std::unordered_map<VertexKey, std::size_t, VertexKeyHash> &lookup,
           const AtomGrid &atom_grid, std::span<const Atom> atoms,
           const Vec3 &point) {
  const VertexKey key = vertex_key(point);
  if (const auto found = lookup.find(key); found != lookup.end()) {
    return found->second;
  }
  const FieldSample sample = sample_union_sdf(atom_grid, atoms, point);
  const Atom &atom = atoms[sample.atom_index];
  const Vec3 normal =
      normalized({point.x - atom.x, point.y - atom.y, point.z - atom.z});
  const std::size_t index = mesh.vertices.size();
  mesh.vertices.push_back({point.x, point.y, point.z, normal.x, normal.y,
                           normal.z, sample.atom_index});
  lookup.emplace(key, index);
  return index;
}

void add_triangle(
    SurfaceMesh &mesh,
    std::unordered_map<VertexKey, std::size_t, VertexKeyHash> &lookup,
    const AtomGrid &atom_grid, std::span<const Atom> atoms, Vec3 a, Vec3 b,
    Vec3 c) {
  const Vec3 centroid = (1.0 / 3.0) * (a + b + c);
  const FieldSample field = sample_union_sdf(atom_grid, atoms, centroid);
  const Atom &normal_atom = atoms[field.atom_index];
  const Vec3 outward =
      normalized({centroid.x - normal_atom.x, centroid.y - normal_atom.y,
                  centroid.z - normal_atom.z});
  Vec3 normal_vector = cross(b - a, c - a);
  if (dot(normal_vector, outward) < 0.0) {
    std::swap(b, c);
    normal_vector = cross(b - a, c - a);
  }
  const double twice_area = norm(normal_vector);
  if (!(twice_area > 1.0e-16) || !std::isfinite(twice_area)) {
    return;
  }
  const Vec3 normal = (1.0 / twice_area) * normal_vector;
  const std::size_t ia = add_vertex(mesh, lookup, atom_grid, atoms, a);
  const std::size_t ib = add_vertex(mesh, lookup, atom_grid, atoms, b);
  const std::size_t ic = add_vertex(mesh, lookup, atom_grid, atoms, c);
  if (ia == ib || ib == ic || ic == ia) {
    return;
  }
  const double area = 0.5 * twice_area;
  Triangle triangle;
  triangle.v0 = ia;
  triangle.v1 = ib;
  triangle.v2 = ic;
  triangle.atom_index = field.atom_index;
  triangle.cx = centroid.x;
  triangle.cy = centroid.y;
  triangle.cz = centroid.z;
  triangle.nx = normal.x;
  triangle.ny = normal.y;
  triangle.nz = normal.z;
  triangle.area = area;
  mesh.triangles.push_back(triangle);

  Panel panel;
  panel.x = centroid.x;
  panel.y = centroid.y;
  panel.z = centroid.z;
  panel.nx = normal.x;
  panel.ny = normal.y;
  panel.nz = normal.z;
  panel.area = area;
  panel.atom_index = field.atom_index;
  panel.v0x = a.x;
  panel.v0y = a.y;
  panel.v0z = a.z;
  panel.v1x = b.x;
  panel.v1y = b.y;
  panel.v1z = b.z;
  panel.v2x = c.x;
  panel.v2y = c.y;
  panel.v2z = c.z;
  panel.v0_index = ia;
  panel.v1_index = ib;
  panel.v2_index = ic;
  panel.has_triangle_geometry = true;
  mesh.panels.push_back(panel);
}

[[nodiscard]] Vec3 interpolate(const Vec3 &a, const Vec3 &b, double fa,
                               double fb) noexcept {
  const double denominator = fa - fb;
  const double t = std::abs(denominator) > 1.0e-300
                       ? std::clamp(fa / denominator, 0.0, 1.0)
                       : 0.5;
  return a + t * (b - a);
}

void polygonize_tetrahedron(
    SurfaceMesh &mesh,
    std::unordered_map<VertexKey, std::size_t, VertexKeyHash> &vertex_lookup,
    const AtomGrid &atom_grid, std::span<const Atom> atoms,
    const std::array<Vec3, 4> &points, const std::array<double, 4> &values) {
  std::array<int, 4> inside{};
  std::array<int, 4> outside{};
  int inside_count = 0;
  int outside_count = 0;
  for (int i = 0; i < 4; ++i) {
    if (values[static_cast<std::size_t>(i)] < 0.0) {
      inside[static_cast<std::size_t>(inside_count++)] = i;
    } else {
      outside[static_cast<std::size_t>(outside_count++)] = i;
    }
  }
  if (inside_count == 0 || inside_count == 4) {
    return;
  }
  const auto edge = [&](int a, int b) {
    return interpolate(points[static_cast<std::size_t>(a)],
                       points[static_cast<std::size_t>(b)],
                       values[static_cast<std::size_t>(a)],
                       values[static_cast<std::size_t>(b)]);
  };
  if (inside_count == 1) {
    const int i = inside[0];
    add_triangle(mesh, vertex_lookup, atom_grid, atoms, edge(i, outside[0]),
                 edge(i, outside[1]), edge(i, outside[2]));
  } else if (inside_count == 3) {
    const int o = outside[0];
    add_triangle(mesh, vertex_lookup, atom_grid, atoms, edge(o, inside[0]),
                 edge(o, inside[2]), edge(o, inside[1]));
  } else {
    const int i0 = inside[0];
    const int i1 = inside[1];
    const int o0 = outside[0];
    const int o1 = outside[1];
    const Vec3 p_i0_o0 = edge(i0, o0);
    const Vec3 p_i0_o1 = edge(i0, o1);
    const Vec3 p_i1_o0 = edge(i1, o0);
    const Vec3 p_i1_o1 = edge(i1, o1);
    add_triangle(mesh, vertex_lookup, atom_grid, atoms, p_i0_o0, p_i0_o1,
                 p_i1_o1);
    add_triangle(mesh, vertex_lookup, atom_grid, atoms, p_i0_o0, p_i1_o1,
                 p_i1_o0);
  }
}

struct DirectedTriangleEdge {
  std::size_t triangle{0};
  int direction{1};
};

[[nodiscard]] EdgeKey edge_key(std::size_t a, std::size_t b) noexcept {
  return {std::min(a, b), std::max(a, b)};
}

[[nodiscard]] int edge_direction(std::size_t a, std::size_t b) noexcept {
  return a < b ? 1 : -1;
}

[[nodiscard]] Vec3 mesh_vertex(const SurfaceMesh &mesh,
                               std::size_t index) noexcept {
  const auto &vertex = mesh.vertices[index];
  return {vertex.x, vertex.y, vertex.z};
}

void update_triangle_geometry(SurfaceMesh &mesh, std::size_t triangle_index) {
  Triangle &triangle = mesh.triangles[triangle_index];
  Panel &panel = mesh.panels[triangle_index];
  const Vec3 a = mesh_vertex(mesh, triangle.v0);
  const Vec3 b = mesh_vertex(mesh, triangle.v1);
  const Vec3 c = mesh_vertex(mesh, triangle.v2);
  const Vec3 normal_vector = cross(b - a, c - a);
  const double twice_area = norm(normal_vector);
  if (!(twice_area > 1.0e-16)) {
    return;
  }
  const Vec3 normal = (1.0 / twice_area) * normal_vector;
  const Vec3 centroid = (1.0 / 3.0) * (a + b + c);
  triangle.cx = centroid.x;
  triangle.cy = centroid.y;
  triangle.cz = centroid.z;
  triangle.nx = normal.x;
  triangle.ny = normal.y;
  triangle.nz = normal.z;
  triangle.area = 0.5 * twice_area;
  panel.x = centroid.x;
  panel.y = centroid.y;
  panel.z = centroid.z;
  panel.nx = normal.x;
  panel.ny = normal.y;
  panel.nz = normal.z;
  panel.area = triangle.area;
  panel.v0x = a.x;
  panel.v0y = a.y;
  panel.v0z = a.z;
  panel.v1x = b.x;
  panel.v1y = b.y;
  panel.v1z = b.z;
  panel.v2x = c.x;
  panel.v2y = c.y;
  panel.v2z = c.z;
  panel.v0_index = triangle.v0;
  panel.v1_index = triangle.v1;
  panel.v2_index = triangle.v2;
}

void orient_mesh_consistently(SurfaceMesh &mesh);

[[nodiscard]] std::size_t collapse_short_manifold_edges(SurfaceMesh &mesh,
                                                        double threshold) {
  if (mesh.vertices.empty() || mesh.triangles.empty() || !(threshold > 0.0))
    return 0;
  struct EdgeRecord {
    std::vector<std::size_t> triangles;
    double length{0.0};
  };
  std::unordered_map<EdgeKey, EdgeRecord, EdgeKeyHash> edges;
  std::vector<std::unordered_set<std::size_t>> neighbors(mesh.vertices.size());
  edges.reserve(mesh.triangles.size() * 3U);
  for (std::size_t triangle_index = 0; triangle_index < mesh.triangles.size();
       ++triangle_index) {
    const Triangle &triangle = mesh.triangles[triangle_index];
    const std::array<std::pair<std::size_t, std::size_t>, 3> local_edges = {
        std::pair{triangle.v0, triangle.v1},
        std::pair{triangle.v1, triangle.v2},
        std::pair{triangle.v2, triangle.v0}};
    for (const auto &[a, b] : local_edges) {
      EdgeRecord &record = edges[edge_key(a, b)];
      record.triangles.push_back(triangle_index);
      neighbors[a].insert(b);
      neighbors[b].insert(a);
    }
  }
  struct Candidate {
    EdgeKey edge;
    double length{0.0};
  };
  std::vector<Candidate> candidates;
  for (auto &[key, record] : edges) {
    const Vec3 a = mesh_vertex(mesh, key.lo);
    const Vec3 b = mesh_vertex(mesh, key.hi);
    record.length = norm(b - a);
    if (record.triangles.size() == 2U && record.length < threshold)
      candidates.push_back({key, record.length});
  }
  std::sort(candidates.begin(), candidates.end(),
            [](const Candidate &a, const Candidate &b) {
              return a.length < b.length;
            });
  std::vector<std::size_t> remap(mesh.vertices.size());
  for (std::size_t index = 0; index < remap.size(); ++index)
    remap[index] = index;
  std::vector<bool> touched(mesh.vertices.size(), false);
  std::size_t collapsed = 0;
  for (const Candidate &candidate : candidates) {
    const std::size_t keep = candidate.edge.lo;
    const std::size_t remove = candidate.edge.hi;
    if (touched[keep] || touched[remove])
      continue;
    std::size_t common_neighbors = 0;
    for (const std::size_t neighbor : neighbors[keep])
      common_neighbors += neighbors[remove].contains(neighbor) ? 1U : 0U;
    // Link condition for collapsing an interior edge of a closed 2-manifold.
    if (common_neighbors != 2U)
      continue;
    Vertex &kept = mesh.vertices[keep];
    const Vertex &removed = mesh.vertices[remove];
    kept.x = 0.5 * (kept.x + removed.x);
    kept.y = 0.5 * (kept.y + removed.y);
    kept.z = 0.5 * (kept.z + removed.z);
    const Vec3 blended_normal = normalized(
        {kept.nx + removed.nx, kept.ny + removed.ny, kept.nz + removed.nz});
    kept.nx = blended_normal.x;
    kept.ny = blended_normal.y;
    kept.nz = blended_normal.z;
    remap[remove] = keep;
    touched[keep] = true;
    touched[remove] = true;
    for (const std::size_t neighbor : neighbors[keep])
      touched[neighbor] = true;
    for (const std::size_t neighbor : neighbors[remove])
      touched[neighbor] = true;
    ++collapsed;
  }
  if (collapsed == 0)
    return 0;

  std::vector<Triangle> triangles;
  std::vector<Panel> panels;
  triangles.reserve(mesh.triangles.size());
  panels.reserve(mesh.panels.size());
  for (std::size_t index = 0; index < mesh.triangles.size(); ++index) {
    Triangle triangle = mesh.triangles[index];
    triangle.v0 = remap[triangle.v0];
    triangle.v1 = remap[triangle.v1];
    triangle.v2 = remap[triangle.v2];
    if (triangle.v0 == triangle.v1 || triangle.v1 == triangle.v2 ||
        triangle.v2 == triangle.v0)
      continue;
    triangles.push_back(triangle);
    panels.push_back(mesh.panels[index]);
  }
  std::vector<bool> referenced(mesh.vertices.size(), false);
  for (const Triangle &triangle : triangles) {
    referenced[triangle.v0] = true;
    referenced[triangle.v1] = true;
    referenced[triangle.v2] = true;
  }
  std::vector<std::size_t> compact(mesh.vertices.size(), 0U);
  std::vector<Vertex> vertices;
  vertices.reserve(mesh.vertices.size() - collapsed);
  for (std::size_t index = 0; index < mesh.vertices.size(); ++index) {
    if (!referenced[index])
      continue;
    compact[index] = vertices.size();
    vertices.push_back(mesh.vertices[index]);
  }
  for (Triangle &triangle : triangles) {
    triangle.v0 = compact[triangle.v0];
    triangle.v1 = compact[triangle.v1];
    triangle.v2 = compact[triangle.v2];
  }
  mesh.vertices = std::move(vertices);
  mesh.triangles = std::move(triangles);
  mesh.panels = std::move(panels);
  for (std::size_t index = 0; index < mesh.triangles.size(); ++index)
    update_triangle_geometry(mesh, index);
  return collapsed;
}

[[nodiscard]] double normalized_triangle_quality(const Vec3 &a, const Vec3 &b,
                                                 const Vec3 &c) noexcept {
  const double ab2 = dot(b - a, b - a);
  const double bc2 = dot(c - b, c - b);
  const double ca2 = dot(a - c, a - c);
  const double denominator = ab2 + bc2 + ca2;
  if (!(denominator > 1.0e-30))
    return 0.0;
  const double twice_area = norm(cross(b - a, c - a));
  return 2.0 * std::sqrt(3.0) * twice_area / denominator;
}

[[nodiscard]] std::size_t flip_edges_for_quality(SurfaceMesh &mesh) {
  struct Owner {
    std::size_t triangle{0};
    std::size_t from{0};
    std::size_t to{0};
    std::size_t opposite{0};
  };
  std::unordered_map<EdgeKey, std::vector<Owner>, EdgeKeyHash> owners;
  owners.reserve(mesh.triangles.size() * 3U);
  for (std::size_t index = 0; index < mesh.triangles.size(); ++index) {
    const Triangle &triangle = mesh.triangles[index];
    owners[edge_key(triangle.v0, triangle.v1)].push_back(
        {index, triangle.v0, triangle.v1, triangle.v2});
    owners[edge_key(triangle.v1, triangle.v2)].push_back(
        {index, triangle.v1, triangle.v2, triangle.v0});
    owners[edge_key(triangle.v2, triangle.v0)].push_back(
        {index, triangle.v2, triangle.v0, triangle.v1});
  }
  struct Candidate {
    EdgeKey edge;
    Owner first;
    Owner second;
    double old_quality{0.0};
  };
  std::vector<Candidate> candidates;
  for (const auto &[edge, edge_owners] : owners) {
    if (edge_owners.size() != 2U)
      continue;
    const Owner &first = edge_owners[0];
    const Owner &second = edge_owners[1];
    if (first.from != second.to || first.to != second.from ||
        first.opposite == second.opposite)
      continue;
    const Vec3 a = mesh_vertex(mesh, first.from);
    const Vec3 b = mesh_vertex(mesh, first.to);
    const Vec3 c = mesh_vertex(mesh, first.opposite);
    const Vec3 d = mesh_vertex(mesh, second.opposite);
    const double old_quality = std::min(normalized_triangle_quality(a, b, c),
                                        normalized_triangle_quality(b, a, d));
    candidates.push_back({edge, first, second, old_quality});
  }
  std::sort(candidates.begin(), candidates.end(),
            [](const Candidate &a, const Candidate &b) {
              return a.old_quality < b.old_quality;
            });
  std::vector<bool> touched(mesh.triangles.size(), false);
  std::size_t flipped = 0;
  for (const Candidate &candidate : candidates) {
    const Owner &first = candidate.first;
    const Owner &second = candidate.second;
    if (touched[first.triangle] || touched[second.triangle])
      continue;
    const EdgeKey replacement = edge_key(first.opposite, second.opposite);
    if (owners.contains(replacement))
      continue;
    const Vec3 a = mesh_vertex(mesh, first.from);
    const Vec3 b = mesh_vertex(mesh, first.to);
    const Vec3 c = mesh_vertex(mesh, first.opposite);
    const Vec3 d = mesh_vertex(mesh, second.opposite);
    const double new_first_quality = normalized_triangle_quality(c, d, b);
    const double new_second_quality = normalized_triangle_quality(d, c, a);
    const double new_quality = std::min(new_first_quality, new_second_quality);
    if (!(new_quality > candidate.old_quality * 1.05 + 1.0e-8))
      continue;
    const Vec3 old_normal = cross(b - a, c - a) + cross(a - b, d - b);
    const Vec3 new_normal = cross(d - c, b - c) + cross(c - d, a - d);
    if (dot(old_normal, new_normal) <= 0.0)
      continue;
    Triangle &triangle_first = mesh.triangles[first.triangle];
    Triangle &triangle_second = mesh.triangles[second.triangle];
    triangle_first.v0 = first.opposite;
    triangle_first.v1 = second.opposite;
    triangle_first.v2 = first.to;
    triangle_second.v0 = second.opposite;
    triangle_second.v1 = first.opposite;
    triangle_second.v2 = first.from;
    update_triangle_geometry(mesh, first.triangle);
    update_triangle_geometry(mesh, second.triangle);
    touched[first.triangle] = true;
    touched[second.triangle] = true;
    ++flipped;
  }
  return flipped;
}

void repair_micro_slivers(SurfaceMesh &mesh, double grid_spacing) {
  orient_mesh_consistently(mesh);
  for (int round = 0; round < 8; ++round) {
    SurfaceMesh candidate = mesh;
    const std::size_t collapsed =
        collapse_short_manifold_edges(candidate, 1.0e-2 * grid_spacing);
    if (collapsed == 0)
      break;
    orient_mesh_consistently(candidate);
    const TraceSpaces audit = build_closed_surface_trace_spaces(candidate);
    if (!audit.quality.valid)
      break;
    mesh = std::move(candidate);
  }
  for (int round = 0; round < 8; ++round) {
    SurfaceMesh candidate = mesh;
    if (flip_edges_for_quality(candidate) == 0)
      return;
    orient_mesh_consistently(candidate);
    const TraceSpaces audit = build_closed_surface_trace_spaces(candidate);
    if (!audit.quality.valid)
      return;
    mesh = std::move(candidate);
  }
}

[[nodiscard]] double mesh_max_edge_ratio(const SurfaceMesh &mesh) noexcept {
  double result = 0.0;
  for (const Triangle &triangle : mesh.triangles) {
    const Vec3 a = mesh_vertex(mesh, triangle.v0);
    const Vec3 b = mesh_vertex(mesh, triangle.v1);
    const Vec3 c = mesh_vertex(mesh, triangle.v2);
    const double ab = norm(b - a);
    const double bc = norm(c - b);
    const double ca = norm(a - c);
    const double shortest = std::min({ab, bc, ca});
    if (shortest > 1.0e-14)
      result = std::max(result, std::max({ab, bc, ca}) / shortest);
  }
  return result;
}

void smooth_analytic_slivers(SurfaceMesh &mesh,
                             const analytic_ses::PatchComplex &patches,
                             double grid_spacing) {
  for (int round = 0; round < 4; ++round) {
    std::vector<std::unordered_set<std::size_t>> neighbors(
        mesh.vertices.size());
    std::vector<Vec3> normal_sums(mesh.vertices.size());
    std::vector<bool> poor(mesh.vertices.size(), false);
    for (const Triangle &triangle : mesh.triangles) {
      const Vec3 a = mesh_vertex(mesh, triangle.v0);
      const Vec3 b = mesh_vertex(mesh, triangle.v1);
      const Vec3 c = mesh_vertex(mesh, triangle.v2);
      const Vec3 face_normal = cross(b - a, c - a);
      normal_sums[triangle.v0] = normal_sums[triangle.v0] + face_normal;
      normal_sums[triangle.v1] = normal_sums[triangle.v1] + face_normal;
      normal_sums[triangle.v2] = normal_sums[triangle.v2] + face_normal;
      neighbors[triangle.v0].insert(triangle.v1);
      neighbors[triangle.v0].insert(triangle.v2);
      neighbors[triangle.v1].insert(triangle.v0);
      neighbors[triangle.v1].insert(triangle.v2);
      neighbors[triangle.v2].insert(triangle.v0);
      neighbors[triangle.v2].insert(triangle.v1);
      const double quality = normalized_triangle_quality(a, b, c);
      const double ab = norm(b - a);
      const double bc = norm(c - b);
      const double ca = norm(a - c);
      const double shortest = std::min({ab, bc, ca});
      const double ratio = shortest > 1.0e-14
                               ? std::max({ab, bc, ca}) / shortest
                               : std::numeric_limits<double>::infinity();
      if (quality < 0.02 || ratio > 50.0) {
        poor[triangle.v0] = true;
        poor[triangle.v1] = true;
        poor[triangle.v2] = true;
      }
    }
    SurfaceMesh candidate = mesh;
    std::size_t moved = 0;
    for (std::size_t index = 0; index < mesh.vertices.size(); ++index) {
      if (!poor[index] || neighbors[index].size() < 3U)
        continue;
      const Vec3 point = mesh_vertex(mesh, index);
      Vec3 average;
      for (const std::size_t neighbor : neighbors[index])
        average = average + mesh_vertex(mesh, neighbor);
      average = (1.0 / static_cast<double>(neighbors[index].size())) * average;
      const Vec3 normal = normalized(normal_sums[index]);
      Vec3 displacement = average - point;
      displacement = displacement - dot(displacement, normal) * normal;
      const double displacement_norm = norm(displacement);
      const double cap = 0.15 * grid_spacing;
      if (displacement_norm > cap)
        displacement = (cap / displacement_norm) * displacement;
      Vec3 proposed = point + 0.25 * displacement;
      const auto hits =
          analytic_ses::intersect_ray_x(patches, proposed.y, proposed.z);
      if (hits.empty())
        continue;
      const auto closest = std::min_element(
          hits.begin(), hits.end(), [&](const auto &left, const auto &right) {
            return std::abs(left.x - proposed.x) <
                   std::abs(right.x - proposed.x);
          });
      if (std::abs(closest->x - point.x) > 0.25 * grid_spacing)
        continue;
      proposed.x = closest->x;
      Vertex &vertex = candidate.vertices[index];
      vertex.x = proposed.x;
      vertex.y = proposed.y;
      vertex.z = proposed.z;
      vertex.nx = closest->normal.x;
      vertex.ny = closest->normal.y;
      vertex.nz = closest->normal.z;
      ++moved;
    }
    if (moved == 0)
      return;
    for (std::size_t index = 0; index < candidate.triangles.size(); ++index)
      update_triangle_geometry(candidate, index);
    orient_mesh_consistently(candidate);
    const TraceSpaces audit = build_closed_surface_trace_spaces(candidate);
    if (!audit.quality.valid)
      return;
    const double old_ratio = mesh_max_edge_ratio(mesh);
    const double new_ratio = mesh_max_edge_ratio(candidate);
    if (!(new_ratio < old_ratio * 0.995))
      return;
    mesh = std::move(candidate);
  }
}

[[nodiscard]] bool
point_inside_triangle_component(const SurfaceMesh &mesh, const Vec3 &point,
                                std::span<const std::size_t> triangle_indices) {
  const Vec3 direction =
      normalized({1.0, 0.3713906763541037, 0.2171827465516321});
  std::vector<double> intersections;
  intersections.reserve(16);
  for (const std::size_t triangle_index : triangle_indices) {
    const Triangle &triangle = mesh.triangles[triangle_index];
    const Vec3 a = mesh_vertex(mesh, triangle.v0);
    const Vec3 b = mesh_vertex(mesh, triangle.v1);
    const Vec3 c = mesh_vertex(mesh, triangle.v2);
    const Vec3 edge1 = b - a;
    const Vec3 edge2 = c - a;
    const Vec3 pvec = cross(direction, edge2);
    const double determinant = dot(edge1, pvec);
    if (std::abs(determinant) <= 1.0e-13)
      continue;
    const double inverse = 1.0 / determinant;
    const Vec3 tvec = point - a;
    const double u = dot(tvec, pvec) * inverse;
    if (u < -1.0e-10 || u > 1.0 + 1.0e-10)
      continue;
    const Vec3 qvec = cross(tvec, edge1);
    const double v = dot(direction, qvec) * inverse;
    if (v < -1.0e-10 || u + v > 1.0 + 1.0e-10)
      continue;
    const double distance = dot(edge2, qvec) * inverse;
    if (distance > 1.0e-9)
      intersections.push_back(distance);
  }
  std::sort(intersections.begin(), intersections.end());
  intersections.erase(
      std::unique(intersections.begin(), intersections.end(),
                  [](double a, double b) {
                    return std::abs(a - b) <=
                           1.0e-8 * std::max({1.0, std::abs(a), std::abs(b)});
                  }),
      intersections.end());
  return (intersections.size() & 1U) != 0U;
}

[[nodiscard]] std::size_t
remove_nonphysical_components(SurfaceMesh &mesh, std::span<const Atom> atoms) {
  if (mesh.triangles.empty())
    return 0;
  std::unordered_map<EdgeKey, std::vector<std::size_t>, EdgeKeyHash> owners;
  owners.reserve(mesh.triangles.size() * 3U);
  for (std::size_t index = 0; index < mesh.triangles.size(); ++index) {
    const Triangle &triangle = mesh.triangles[index];
    owners[edge_key(triangle.v0, triangle.v1)].push_back(index);
    owners[edge_key(triangle.v1, triangle.v2)].push_back(index);
    owners[edge_key(triangle.v2, triangle.v0)].push_back(index);
  }
  std::vector<std::vector<std::size_t>> adjacency(mesh.triangles.size());
  for (const auto &[edge, edge_owners] : owners) {
    (void)edge;
    if (edge_owners.size() != 2U)
      continue;
    adjacency[edge_owners[0]].push_back(edge_owners[1]);
    adjacency[edge_owners[1]].push_back(edge_owners[0]);
  }
  std::vector<int> component_of(mesh.triangles.size(), -1);
  std::vector<std::vector<std::size_t>> components;
  for (std::size_t root = 0; root < mesh.triangles.size(); ++root) {
    if (component_of[root] >= 0)
      continue;
    const int component = static_cast<int>(components.size());
    components.emplace_back();
    std::vector<std::size_t> stack = {root};
    component_of[root] = component;
    while (!stack.empty()) {
      const std::size_t current = stack.back();
      stack.pop_back();
      components.back().push_back(current);
      for (const std::size_t neighbor : adjacency[current]) {
        if (component_of[neighbor] < 0) {
          component_of[neighbor] = component;
          stack.push_back(neighbor);
        }
      }
    }
  }
  if (components.size() <= 1U)
    return 0;

  struct Bounds {
    Vec3 lower{std::numeric_limits<double>::infinity(),
               std::numeric_limits<double>::infinity(),
               std::numeric_limits<double>::infinity()};
    Vec3 upper{-std::numeric_limits<double>::infinity(),
               -std::numeric_limits<double>::infinity(),
               -std::numeric_limits<double>::infinity()};
    Vec3 sample;
  };
  std::vector<Bounds> bounds(components.size());
  for (std::size_t component = 0; component < components.size(); ++component) {
    for (const std::size_t triangle_index : components[component]) {
      const Triangle &triangle = mesh.triangles[triangle_index];
      for (const std::size_t vertex_index :
           {triangle.v0, triangle.v1, triangle.v2}) {
        const Vec3 point = mesh_vertex(mesh, vertex_index);
        bounds[component].lower.x =
            std::min(bounds[component].lower.x, point.x);
        bounds[component].lower.y =
            std::min(bounds[component].lower.y, point.y);
        bounds[component].lower.z =
            std::min(bounds[component].lower.z, point.z);
        bounds[component].upper.x =
            std::max(bounds[component].upper.x, point.x);
        bounds[component].upper.y =
            std::max(bounds[component].upper.y, point.y);
        bounds[component].upper.z =
            std::max(bounds[component].upper.z, point.z);
      }
    }
    const Triangle &triangle = mesh.triangles[components[component].front()];
    bounds[component].sample = (1.0 / 3.0) * (mesh_vertex(mesh, triangle.v0) +
                                              mesh_vertex(mesh, triangle.v1) +
                                              mesh_vertex(mesh, triangle.v2));
  }
  std::vector<bool> keep(components.size(), false);
  for (std::size_t component = 0; component < components.size(); ++component) {
    const Bounds &box = bounds[component];
    for (const Atom &atom : atoms) {
      const Vec3 center{atom.x, atom.y, atom.z};
      if (center.x <= box.lower.x || center.x >= box.upper.x ||
          center.y <= box.lower.y || center.y >= box.upper.y ||
          center.z <= box.lower.z || center.z >= box.upper.z)
        continue;
      if (point_inside_triangle_component(mesh, center,
                                          components[component])) {
        keep[component] = true;
        break;
      }
    }
  }
  for (std::size_t inner = 0; inner < components.size(); ++inner) {
    if (!keep[inner])
      continue;
    const Vec3 point = bounds[inner].sample;
    for (std::size_t outer = 0; outer < components.size(); ++outer) {
      if (inner == outer || !keep[outer])
        continue;
      const Bounds &box = bounds[outer];
      if (point.x <= box.lower.x || point.x >= box.upper.x ||
          point.y <= box.lower.y || point.y >= box.upper.y ||
          point.z <= box.lower.z || point.z >= box.upper.z)
        continue;
      if (point_inside_triangle_component(mesh, point, components[outer])) {
        keep[inner] = false;
        break;
      }
    }
  }
  const std::size_t removed_components =
      static_cast<std::size_t>(std::count(keep.begin(), keep.end(), false));
  if (removed_components == 0)
    return 0;
  std::vector<Triangle> triangles;
  std::vector<Panel> panels;
  for (std::size_t index = 0; index < mesh.triangles.size(); ++index) {
    if (keep[static_cast<std::size_t>(component_of[index])]) {
      triangles.push_back(mesh.triangles[index]);
      panels.push_back(mesh.panels[index]);
    }
  }
  std::vector<bool> referenced(mesh.vertices.size(), false);
  for (const Triangle &triangle : triangles) {
    referenced[triangle.v0] = true;
    referenced[triangle.v1] = true;
    referenced[triangle.v2] = true;
  }
  std::vector<std::size_t> compact(mesh.vertices.size(), 0U);
  std::vector<Vertex> vertices;
  for (std::size_t index = 0; index < mesh.vertices.size(); ++index) {
    if (!referenced[index])
      continue;
    compact[index] = vertices.size();
    vertices.push_back(mesh.vertices[index]);
  }
  for (Triangle &triangle : triangles) {
    triangle.v0 = compact[triangle.v0];
    triangle.v1 = compact[triangle.v1];
    triangle.v2 = compact[triangle.v2];
  }
  mesh.vertices = std::move(vertices);
  mesh.triangles = std::move(triangles);
  mesh.panels = std::move(panels);
  for (std::size_t index = 0; index < mesh.triangles.size(); ++index)
    update_triangle_geometry(mesh, index);
  return removed_components;
}

// Marching-tetrahedra cases are locally oriented, but a shared diagonal can
// still give neighboring tetrahedra the same edge direction. Repair the mesh
// deterministically before the closed-surface trace-space check.
void orient_mesh_consistently(SurfaceMesh &mesh) {
  if (mesh.triangles.empty() || mesh.panels.size() != mesh.triangles.size()) {
    return;
  }
  std::unordered_map<EdgeKey, std::vector<DirectedTriangleEdge>, EdgeKeyHash>
      edge_owners;
  edge_owners.reserve(mesh.triangles.size() * 3);
  for (std::size_t index = 0; index < mesh.triangles.size(); ++index) {
    const Triangle &triangle = mesh.triangles[index];
    const std::array<std::pair<std::size_t, std::size_t>, 3> edges = {
        std::pair{triangle.v0, triangle.v1},
        std::pair{triangle.v1, triangle.v2},
        std::pair{triangle.v2, triangle.v0}};
    for (const auto &[a, b] : edges) {
      edge_owners[edge_key(a, b)].push_back({index, edge_direction(a, b)});
    }
  }

  struct Neighbor {
    std::size_t index{0};
    int direction_product{1};
  };
  std::vector<std::vector<Neighbor>> adjacency(mesh.triangles.size());
  for (const auto &[key, owners] : edge_owners) {
    (void)key;
    if (owners.size() != 2) {
      continue;
    }
    adjacency[owners[0].triangle].push_back(
        {owners[1].triangle, owners[0].direction * owners[1].direction});
    adjacency[owners[1].triangle].push_back(
        {owners[0].triangle, owners[0].direction * owners[1].direction});
  }

  std::vector<int> orientation(mesh.triangles.size(), 0);
  std::vector<std::vector<std::size_t>> components;
  for (std::size_t root = 0; root < mesh.triangles.size(); ++root) {
    if (orientation[root] != 0) {
      continue;
    }
    orientation[root] = 1;
    components.emplace_back();
    std::vector<std::size_t> stack = {root};
    while (!stack.empty()) {
      const std::size_t current = stack.back();
      stack.pop_back();
      components.back().push_back(current);
      for (const Neighbor neighbor : adjacency[current]) {
        const int required = -orientation[current] * neighbor.direction_product;
        if (orientation[neighbor.index] == 0) {
          orientation[neighbor.index] = required;
          stack.push_back(neighbor.index);
        }
      }
    }
  }

  for (const auto &component : components) {
    double signed_volume_six = 0.0;
    for (const std::size_t index : component) {
      const Triangle &triangle = mesh.triangles[index];
      const Vec3 a = mesh_vertex(mesh, triangle.v0);
      const Vec3 b = mesh_vertex(mesh, triangle.v1);
      const Vec3 c = mesh_vertex(mesh, triangle.v2);
      signed_volume_six +=
          static_cast<double>(orientation[index]) * dot(a, cross(b, c));
    }
    if (signed_volume_six < 0.0) {
      for (const std::size_t index : component) {
        orientation[index] = -orientation[index];
      }
    }
  }

  for (std::size_t index = 0; index < mesh.triangles.size(); ++index) {
    if (orientation[index] < 0) {
      Triangle &triangle = mesh.triangles[index];
      std::swap(triangle.v1, triangle.v2);
    }
    update_triangle_geometry(mesh, index);
  }
}

void polygonize_grid(SurfaceMesh &mesh, const AtomGrid &ownership_atom_grid,
                     std::span<const Atom> atoms, const Vec3 &lower,
                     double grid_spacing, std::size_t nx, std::size_t ny,
                     std::size_t nz, std::span<const double> values) {
  const auto grid_index = [nx, ny](std::size_t ix, std::size_t iy,
                                   std::size_t iz) {
    return ix + nx * (iy + ny * iz);
  };
  static constexpr std::array<std::array<int, 4>, 6> tetrahedra = {
      std::array<int, 4>{0, 5, 1, 6}, std::array<int, 4>{0, 1, 2, 6},
      std::array<int, 4>{0, 2, 3, 6}, std::array<int, 4>{0, 3, 7, 6},
      std::array<int, 4>{0, 7, 4, 6}, std::array<int, 4>{0, 4, 5, 6}};
  std::unordered_map<VertexKey, std::size_t, VertexKeyHash> vertex_lookup;
  vertex_lookup.reserve((nx - 1) * (ny - 1));
  for (std::size_t iz = 0; iz + 1 < nz; ++iz) {
    for (std::size_t iy = 0; iy + 1 < ny; ++iy) {
      for (std::size_t ix = 0; ix + 1 < nx; ++ix) {
        const std::array<Vec3, 8> cube_points = {
            Vec3{lower.x + grid_spacing * static_cast<double>(ix),
                 lower.y + grid_spacing * static_cast<double>(iy),
                 lower.z + grid_spacing * static_cast<double>(iz)},
            Vec3{lower.x + grid_spacing * static_cast<double>(ix + 1),
                 lower.y + grid_spacing * static_cast<double>(iy),
                 lower.z + grid_spacing * static_cast<double>(iz)},
            Vec3{lower.x + grid_spacing * static_cast<double>(ix + 1),
                 lower.y + grid_spacing * static_cast<double>(iy + 1),
                 lower.z + grid_spacing * static_cast<double>(iz)},
            Vec3{lower.x + grid_spacing * static_cast<double>(ix),
                 lower.y + grid_spacing * static_cast<double>(iy + 1),
                 lower.z + grid_spacing * static_cast<double>(iz)},
            Vec3{lower.x + grid_spacing * static_cast<double>(ix),
                 lower.y + grid_spacing * static_cast<double>(iy),
                 lower.z + grid_spacing * static_cast<double>(iz + 1)},
            Vec3{lower.x + grid_spacing * static_cast<double>(ix + 1),
                 lower.y + grid_spacing * static_cast<double>(iy),
                 lower.z + grid_spacing * static_cast<double>(iz + 1)},
            Vec3{lower.x + grid_spacing * static_cast<double>(ix + 1),
                 lower.y + grid_spacing * static_cast<double>(iy + 1),
                 lower.z + grid_spacing * static_cast<double>(iz + 1)},
            Vec3{lower.x + grid_spacing * static_cast<double>(ix),
                 lower.y + grid_spacing * static_cast<double>(iy + 1),
                 lower.z + grid_spacing * static_cast<double>(iz + 1)}};
        const std::array<double, 8> cube_values = {
            values[grid_index(ix, iy, iz)],
            values[grid_index(ix + 1, iy, iz)],
            values[grid_index(ix + 1, iy + 1, iz)],
            values[grid_index(ix, iy + 1, iz)],
            values[grid_index(ix, iy, iz + 1)],
            values[grid_index(ix + 1, iy, iz + 1)],
            values[grid_index(ix + 1, iy + 1, iz + 1)],
            values[grid_index(ix, iy + 1, iz + 1)]};
        for (const auto &tetrahedron : tetrahedra) {
          std::array<Vec3, 4> points{};
          std::array<double, 4> tetra_values{};
          for (std::size_t k = 0; k < tetrahedron.size(); ++k) {
            const std::size_t vertex = static_cast<std::size_t>(tetrahedron[k]);
            points[k] = cube_points[vertex];
            tetra_values[k] = cube_values[vertex];
          }
          polygonize_tetrahedron(mesh, vertex_lookup, ownership_atom_grid,
                                 atoms, points, tetra_values);
        }
      }
    }
  }
}

} // namespace

AuxiliarySurfaceResult
build_implicit_surface(std::span<const Atom> atoms, double grid_spacing,
                       double probe_radius, TitanExecutor *executor,
                       std::size_t max_grid_vertices, bool validate_topology) {
  const auto begin = std::chrono::steady_clock::now();
  AuxiliarySurfaceResult result;
  if (atoms.empty() || !(grid_spacing > 0.0) || !std::isfinite(grid_spacing) ||
      !(probe_radius >= 0.0) || !std::isfinite(probe_radius)) {
    result.status = "invalid_input";
    return result;
  }

  Vec3 lower{std::numeric_limits<double>::infinity(),
             std::numeric_limits<double>::infinity(),
             std::numeric_limits<double>::infinity()};
  Vec3 upper{-std::numeric_limits<double>::infinity(),
             -std::numeric_limits<double>::infinity(),
             -std::numeric_limits<double>::infinity()};
  std::size_t valid_atoms = 0;
  for (const Atom &atom : atoms) {
    if (!(atom.radius > 0.0) || !std::isfinite(atom.radius)) {
      continue;
    }
    ++valid_atoms;
    const double expanded_radius = atom.radius + probe_radius;
    lower.x = std::min(lower.x, atom.x - expanded_radius);
    lower.y = std::min(lower.y, atom.y - expanded_radius);
    lower.z = std::min(lower.z, atom.z - expanded_radius);
    upper.x = std::max(upper.x, atom.x + expanded_radius);
    upper.y = std::max(upper.y, atom.y + expanded_radius);
    upper.z = std::max(upper.z, atom.z + expanded_radius);
  }
  if (valid_atoms == 0) {
    result.status = "no_valid_atoms";
    return result;
  }
  const double margin = 2.0 * grid_spacing;
  lower = lower - Vec3{margin, margin, margin};
  upper = upper + Vec3{margin, margin, margin};
  const auto axis_count = [&](double lo, double hi) {
    return static_cast<std::size_t>(std::ceil((hi - lo) / grid_spacing)) + 1;
  };
  const std::size_t nx = axis_count(lower.x, upper.x);
  const std::size_t ny = axis_count(lower.y, upper.y);
  const std::size_t nz = axis_count(lower.z, upper.z);
  const std::size_t grid_vertex_cap =
      std::max(kMaxGridVertices, max_grid_vertices);
  if (nx < 2 || ny < 2 || nz < 2 || nx > grid_vertex_cap / ny ||
      nx * ny > grid_vertex_cap / nz) {
    result.status = "grid_too_large";
    return result;
  }
  result.grid_vertices = nx * ny * nz;
  std::vector<Atom> field_atoms(atoms.begin(), atoms.end());
  if (probe_radius > 0.0) {
    for (Atom &atom : field_atoms) {
      atom.radius += probe_radius;
    }
  }
  const AtomGrid field_atom_grid = build_atom_grid(field_atoms, lower);
  const AtomGrid ownership_atom_grid = build_atom_grid(atoms, lower);
  std::vector<double> values(result.grid_vertices, 1.0);
  const auto grid_index = [nx, ny](std::size_t ix, std::size_t iy,
                                   std::size_t iz) {
    return ix + nx * (iy + ny * iz);
  };
  const auto evaluate_slices = [&](std::size_t begin_z, std::size_t end_z) {
    for (std::size_t iz = begin_z; iz < end_z; ++iz) {
      for (std::size_t iy = 0; iy < ny; ++iy) {
        for (std::size_t ix = 0; ix < nx; ++ix) {
          const Vec3 point{lower.x + grid_spacing * static_cast<double>(ix),
                           lower.y + grid_spacing * static_cast<double>(iy),
                           lower.z + grid_spacing * static_cast<double>(iz)};
          if (probe_radius == 0.0) {
            values[grid_index(ix, iy, iz)] =
                sample_union_sdf(field_atom_grid, field_atoms, point).value;
            continue;
          }

          // The solvent-excluded molecular volume is the morphological closing
          // (VDW union dilated by the probe, then eroded by the same probe).
          // Evaluate the erosion with a symmetric 26-direction spherical
          // stencil at every grid vertex. The zero isosurface is subsequently
          // welded and topology-audited, so no open patch reaches the solver.
          double eroded_value =
              sample_union_sdf(field_atom_grid, field_atoms, point).value;
          for (int dz = -1; dz <= 1; ++dz) {
            for (int dy = -1; dy <= 1; ++dy) {
              for (int dx = -1; dx <= 1; ++dx) {
                if (dx == 0 && dy == 0 && dz == 0) {
                  continue;
                }
                const Vec3 direction = normalized({static_cast<double>(dx),
                                                   static_cast<double>(dy),
                                                   static_cast<double>(dz)});
                const Vec3 probe_point = point + probe_radius * direction;
                eroded_value = std::max(
                    eroded_value,
                    sample_union_sdf(field_atom_grid, field_atoms, probe_point)
                        .value);
              }
            }
          }
          values[grid_index(ix, iy, iz)] = eroded_value;
        }
      }
    }
  };
  if (executor != nullptr && executor->capacity() > 1 && nz > 1) {
    executor->parallel_for(
        nz, std::min<std::size_t>(executor->capacity(), nz), 1,
        [&](std::size_t, std::size_t begin_z, std::size_t end_z) {
          evaluate_slices(begin_z, end_z);
        });
  } else {
    evaluate_slices(0, nz);
  }

  if (probe_radius > 0.0) {
    // Keep only solvent connected to the grid boundary. Positive pockets that
    // cannot be reached by a probe from the exterior are part of the molecular
    // volume, not additional dielectric interfaces.
    std::vector<std::uint8_t> exterior(result.grid_vertices, 0U);
    std::vector<std::size_t> queue;
    queue.reserve(result.grid_vertices / 8U);
    const auto seed = [&](std::size_t ix, std::size_t iy, std::size_t iz) {
      const std::size_t index = grid_index(ix, iy, iz);
      if (values[index] >= 0.0 && exterior[index] == 0U) {
        exterior[index] = 1U;
        queue.push_back(index);
      }
    };
    for (std::size_t iz = 0; iz < nz; ++iz) {
      for (std::size_t iy = 0; iy < ny; ++iy) {
        seed(0, iy, iz);
        seed(nx - 1, iy, iz);
      }
    }
    for (std::size_t iz = 0; iz < nz; ++iz) {
      for (std::size_t ix = 0; ix < nx; ++ix) {
        seed(ix, 0, iz);
        seed(ix, ny - 1, iz);
      }
    }
    for (std::size_t iy = 0; iy < ny; ++iy) {
      for (std::size_t ix = 0; ix < nx; ++ix) {
        seed(ix, iy, 0);
        seed(ix, iy, nz - 1);
      }
    }
    for (std::size_t head = 0; head < queue.size(); ++head) {
      const std::size_t index = queue[head];
      const std::size_t iz = index / (nx * ny);
      const std::size_t remainder = index - iz * nx * ny;
      const std::size_t iy = remainder / nx;
      const std::size_t ix = remainder - iy * nx;
      const auto visit = [&](std::size_t neighbor) {
        if (values[neighbor] >= 0.0 && exterior[neighbor] == 0U) {
          exterior[neighbor] = 1U;
          queue.push_back(neighbor);
        }
      };
      if (ix > 0)
        visit(index - 1);
      if (ix + 1 < nx)
        visit(index + 1);
      if (iy > 0)
        visit(index - nx);
      if (iy + 1 < ny)
        visit(index + nx);
      if (iz > 0)
        visit(index - nx * ny);
      if (iz + 1 < nz)
        visit(index + nx * ny);
    }
    for (std::size_t index = 0; index < values.size(); ++index) {
      if (values[index] >= 0.0 && exterior[index] == 0U) {
        values[index] = -grid_spacing;
      }
    }
  }

  polygonize_grid(result.mesh, ownership_atom_grid, atoms, lower, grid_spacing,
                  nx, ny, nz, values);
  repair_micro_slivers(result.mesh, grid_spacing);

  if (result.mesh.triangles.empty()) {
    result.status = "empty_implicit_surface";
  } else if (!validate_topology) {
    result.trace_spaces.quality.vertex_count = result.mesh.vertices.size();
    result.trace_spaces.quality.triangle_count = result.mesh.triangles.size();
    result.trace_spaces.quality.failure_reason = "surface_only_not_audited";
    result.status = "pass";
  } else {
    orient_mesh_consistently(result.mesh);
    result.trace_spaces = build_closed_surface_trace_spaces(result.mesh);
    result.status = result.trace_spaces.quality.failure_reason;
  }
  const auto end = std::chrono::steady_clock::now();
  result.build_ms =
      std::chrono::duration<double, std::milli>(end - begin).count();
  return result;
}

AuxiliarySurfaceResult build_implicit_vdw_auxiliary_surface(
    std::span<const Atom> atoms, double grid_spacing, TitanExecutor *executor) {
  return build_implicit_surface(atoms, grid_spacing, 0.0, executor,
                                kMaxGridVertices, true);
}

AuxiliarySurfaceResult build_implicit_ses_surface(std::span<const Atom> atoms,
                                                  double grid_spacing,
                                                  double probe_radius,
                                                  TitanExecutor *executor,
                                                  std::size_t max_grid_vertices,
                                                  bool validate_topology) {
  return build_implicit_surface(atoms, grid_spacing, probe_radius, executor,
                                max_grid_vertices, validate_topology);
}

AuxiliarySurfaceResult build_analytic_ses_surface(std::span<const Atom> atoms,
                                                  double grid_spacing,
                                                  double probe_radius,
                                                  TitanExecutor *executor,
                                                  std::size_t max_grid_vertices,
                                                  bool validate_topology) {
  const auto begin = std::chrono::steady_clock::now();
  AuxiliarySurfaceResult result;
  if (atoms.empty() || !(grid_spacing > 0.0) || !std::isfinite(grid_spacing) ||
      !(probe_radius > 0.0) || !std::isfinite(probe_radius)) {
    result.status = "invalid_input";
    return result;
  }
  const analytic_ses::PatchComplex patches =
      analytic_ses::build_patch_complex(atoms, probe_radius);
  if (patches.status != "pass") {
    result.status = "analytic_patch_" + patches.status;
    return result;
  }

  Vec3 lower{std::numeric_limits<double>::infinity(),
             std::numeric_limits<double>::infinity(),
             std::numeric_limits<double>::infinity()};
  Vec3 upper{-std::numeric_limits<double>::infinity(),
             -std::numeric_limits<double>::infinity(),
             -std::numeric_limits<double>::infinity()};
  for (const Atom &atom : patches.atoms) {
    const double radius = atom.radius + probe_radius;
    lower.x = std::min(lower.x, atom.x - radius);
    lower.y = std::min(lower.y, atom.y - radius);
    lower.z = std::min(lower.z, atom.z - radius);
    upper.x = std::max(upper.x, atom.x + radius);
    upper.y = std::max(upper.y, atom.y + radius);
    upper.z = std::max(upper.z, atom.z + radius);
  }
  const double margin = 2.0 * grid_spacing;
  lower = lower - Vec3{margin, margin, margin};
  upper = upper + Vec3{margin, margin, margin};
  const auto axis_count = [&](double lo, double hi) {
    return static_cast<std::size_t>(std::ceil((hi - lo) / grid_spacing)) + 1U;
  };
  const std::size_t nx = axis_count(lower.x, upper.x);
  const std::size_t ny = axis_count(lower.y, upper.y);
  const std::size_t nz = axis_count(lower.z, upper.z);
  const std::size_t cap = std::max<std::size_t>(1U, max_grid_vertices);
  if (nx < 2 || ny < 2 || nz < 2 || nx > cap / ny || nx * ny > cap / nz) {
    result.status = "grid_too_large";
    return result;
  }
  result.grid_vertices = nx * ny * nz;
  const auto grid_index = [nx, ny](std::size_t ix, std::size_t iy,
                                   std::size_t iz) {
    return ix + nx * (iy + ny * iz);
  };
  // Use the flood-filled configuration-space closing only as an exterior
  // topology mask. Intersection locations and final vertex projection still
  // come from the analytic contact/toroidal/reentrant patches below.
  std::vector<Atom> expanded_atoms = patches.atoms;
  for (Atom &atom : expanded_atoms)
    atom.radius += probe_radius;
  const AtomGrid expanded_grid = build_atom_grid(expanded_atoms, lower);
  // This mask only carries topology/sign. Exact surface distances below come
  // from analytic ray crossings, so avoid 27 square roots and minimum-distance
  // reductions per grid vertex here.
  std::vector<std::uint8_t> topology_inside(result.grid_vertices, 0U);
  const auto &sample_directions = topology_sample_directions();
  const auto evaluate_topology_slices = [&](std::size_t begin_z,
                                            std::size_t end_z) {
    for (std::size_t iz = begin_z; iz < end_z; ++iz) {
      for (std::size_t iy = 0; iy < ny; ++iy) {
        for (std::size_t ix = 0; ix < nx; ++ix) {
          const Vec3 point{lower.x + grid_spacing * static_cast<double>(ix),
                           lower.y + grid_spacing * static_cast<double>(iy),
                           lower.z + grid_spacing * static_cast<double>(iz)};
          bool inside = true;
          for (const Vec3 &direction : sample_directions) {
            if (!point_inside_union(expanded_grid, expanded_atoms,
                                    point + probe_radius * direction)) {
              inside = false;
              break;
            }
          }
          topology_inside[grid_index(ix, iy, iz)] = inside ? 1U : 0U;
        }
      }
    }
  };
  if (executor != nullptr && executor->capacity() > 1 && nz > 1) {
    executor->parallel_for(
        nz, std::min<std::size_t>(executor->capacity(), nz), 1,
        [&](std::size_t, std::size_t begin_z, std::size_t end_z) {
          evaluate_topology_slices(begin_z, end_z);
        });
  } else {
    evaluate_topology_slices(0, nz);
  }
  std::vector<std::uint8_t> exterior(result.grid_vertices, 0U);
  std::vector<std::uint8_t> cavity_fill(result.grid_vertices, 0U);
  std::vector<std::size_t> queue;
  queue.reserve(result.grid_vertices / 8U);
  const auto seed = [&](std::size_t ix, std::size_t iy, std::size_t iz) {
    const std::size_t index = grid_index(ix, iy, iz);
    if (topology_inside[index] == 0U && exterior[index] == 0U) {
      exterior[index] = 1U;
      queue.push_back(index);
    }
  };
  for (std::size_t iz = 0; iz < nz; ++iz) {
    for (std::size_t iy = 0; iy < ny; ++iy) {
      seed(0, iy, iz);
      seed(nx - 1, iy, iz);
    }
  }
  for (std::size_t iz = 0; iz < nz; ++iz) {
    for (std::size_t ix = 0; ix < nx; ++ix) {
      seed(ix, 0, iz);
      seed(ix, ny - 1, iz);
    }
  }
  for (std::size_t iy = 0; iy < ny; ++iy) {
    for (std::size_t ix = 0; ix < nx; ++ix) {
      seed(ix, iy, 0);
      seed(ix, iy, nz - 1);
    }
  }
  for (std::size_t head = 0; head < queue.size(); ++head) {
    const std::size_t index = queue[head];
    const std::size_t iz = index / (nx * ny);
    const std::size_t remainder = index - iz * nx * ny;
    const std::size_t iy = remainder / nx;
    const std::size_t ix = remainder - iy * nx;
    const auto visit = [&](std::size_t neighbor) {
      if (topology_inside[neighbor] == 0U && exterior[neighbor] == 0U) {
        exterior[neighbor] = 1U;
        queue.push_back(neighbor);
      }
    };
    if (ix > 0)
      visit(index - 1);
    if (ix + 1 < nx)
      visit(index + 1);
    if (iy > 0)
      visit(index - nx);
    if (iy + 1 < ny)
      visit(index + nx);
    if (iz > 0)
      visit(index - nx * ny);
    if (iz + 1 < nz)
      visit(index + nx * ny);
  }
  for (std::size_t index = 0; index < topology_inside.size(); ++index) {
    if (topology_inside[index] == 0U && exterior[index] == 0U) {
      topology_inside[index] = 1U;
      cavity_fill[index] = 1U;
    }
  }

  // Marching tetrahedra only uses field magnitudes on edges whose endpoints
  // have opposite signs. Reconstruct the original eroded-SDF magnitude just
  // for that one-vertex boundary band; all bulk vertices retain a cheap
  // grid-spacing fallback. This preserves the old interpolation geometry
  // without paying for 27 SDF evaluations throughout the volume.
  std::vector<double> topology_magnitude(result.grid_vertices, grid_spacing);
  const auto evaluate_boundary_magnitudes = [&](std::size_t begin_z,
                                                std::size_t end_z) {
    for (std::size_t iz = begin_z; iz < end_z; ++iz) {
      for (std::size_t iy = 0; iy < ny; ++iy) {
        for (std::size_t ix = 0; ix < nx; ++ix) {
          const std::size_t index = grid_index(ix, iy, iz);
          if (cavity_fill[index] != 0U)
            continue;
          const bool inside = topology_inside[index] != 0U;
          bool boundary = false;
          for (int dz = -1; dz <= 1 && !boundary; ++dz) {
            for (int dy = -1; dy <= 1 && !boundary; ++dy) {
              for (int dx = -1; dx <= 1; ++dx) {
                if (dx == 0 && dy == 0 && dz == 0)
                  continue;
                const auto neighbor_x = static_cast<std::ptrdiff_t>(ix) + dx;
                const auto neighbor_y = static_cast<std::ptrdiff_t>(iy) + dy;
                const auto neighbor_z = static_cast<std::ptrdiff_t>(iz) + dz;
                if (neighbor_x < 0 || neighbor_y < 0 || neighbor_z < 0 ||
                    neighbor_x >= static_cast<std::ptrdiff_t>(nx) ||
                    neighbor_y >= static_cast<std::ptrdiff_t>(ny) ||
                    neighbor_z >= static_cast<std::ptrdiff_t>(nz)) {
                  continue;
                }
                const std::size_t neighbor =
                    grid_index(static_cast<std::size_t>(neighbor_x),
                               static_cast<std::size_t>(neighbor_y),
                               static_cast<std::size_t>(neighbor_z));
                if ((topology_inside[neighbor] != 0U) != inside) {
                  boundary = true;
                  break;
                }
              }
            }
          }
          if (!boundary)
            continue;
          const Vec3 point{lower.x + grid_spacing * static_cast<double>(ix),
                           lower.y + grid_spacing * static_cast<double>(iy),
                           lower.z + grid_spacing * static_cast<double>(iz)};
          double eroded = -std::numeric_limits<double>::infinity();
          for (const Vec3 &direction : sample_directions) {
            eroded = std::max(eroded,
                              sample_union_sdf(expanded_grid, expanded_atoms,
                                               point + probe_radius * direction)
                                  .value);
          }
          topology_magnitude[index] = std::abs(eroded);
        }
      }
    }
  };
  if (executor != nullptr && executor->capacity() > 1 && nz > 1) {
    executor->parallel_for(
        nz, std::min<std::size_t>(executor->capacity(), nz), 1,
        [&](std::size_t, std::size_t begin_z, std::size_t end_z) {
          evaluate_boundary_magnitudes(begin_z, end_z);
        });
  } else {
    evaluate_boundary_magnitudes(0, nz);
  }

  std::vector<double> values(result.grid_vertices, grid_spacing);
  std::atomic_size_t odd_ray_rows{0};
  const auto evaluate_rows = [&](std::size_t begin_row, std::size_t end_row) {
    for (std::size_t row = begin_row; row < end_row; ++row) {
      const std::size_t iz = row / ny;
      const std::size_t iy = row - iz * ny;
      const double y = lower.y + grid_spacing * static_cast<double>(iy);
      const double z = lower.z + grid_spacing * static_cast<double>(iz);
      const std::vector<analytic_ses::RayHit> hits =
          analytic_ses::intersect_ray_x(patches, y, z);
      std::vector<double> crossings;
      crossings.reserve(hits.size());
      for (const auto &hit : hits) {
        if (std::abs(hit.normal.x) > 1.0e-8)
          crossings.push_back(hit.x);
      }
      if ((crossings.size() & 1U) != 0U)
        odd_ray_rows.fetch_add(1U, std::memory_order_relaxed);
      for (std::size_t ix = 0; ix < nx; ++ix) {
        const double x = lower.x + grid_spacing * static_cast<double>(ix);
        const auto next =
            std::lower_bound(crossings.begin(), crossings.end(), x);
        double distance = topology_magnitude[grid_index(ix, iy, iz)];
        if (next != crossings.end())
          distance = std::min(distance, std::abs(*next - x));
        if (next != crossings.begin())
          distance = std::min(distance, std::abs(*std::prev(next) - x));
        distance = std::max(distance, 1.0e-9 * grid_spacing);
        values[grid_index(ix, iy, iz)] =
            topology_inside[grid_index(ix, iy, iz)] != 0U ? -distance
                                                          : distance;
      }
    }
  };
  const std::size_t row_count = ny * nz;
  if (executor != nullptr && executor->capacity() > 1 && row_count > 1) {
    const std::size_t worker_count =
        std::min<std::size_t>(executor->capacity(), row_count);
    // Central rays intersect many more patches than peripheral rays. Static
    // contiguous z/y ranges therefore leave workers idle near the end of this
    // phase. Keep at least 64 tiles per worker, with bounds that amortize the
    // atomic fetch without making tiles too coarse for small molecules.
    const std::size_t target_tiles = 64U * worker_count;
    const std::size_t tile_rows =
        std::clamp((row_count + target_tiles - 1U) / target_tiles,
                   std::size_t{2U}, std::size_t{8U});
    std::atomic_size_t next_row{0U};
    executor->parallel_do(worker_count, [&](std::size_t) {
      for (;;) {
        const std::size_t row_begin =
            next_row.fetch_add(tile_rows, std::memory_order_relaxed);
        if (row_begin >= row_count)
          break;
        evaluate_rows(row_begin, std::min(row_count, row_begin + tile_rows));
      }
    });
  } else {
    evaluate_rows(0, row_count);
  }
  result.odd_ray_rows = odd_ray_rows.load(std::memory_order_relaxed);

  const AtomGrid ownership_grid = build_atom_grid(patches.atoms, lower);
  polygonize_grid(result.mesh, ownership_grid, patches.atoms, lower,
                  grid_spacing, nx, ny, nz, values);
  repair_micro_slivers(result.mesh, grid_spacing);
  if (remove_nonphysical_components(result.mesh, patches.atoms) > 0U)
    orient_mesh_consistently(result.mesh);
  smooth_analytic_slivers(result.mesh, patches, grid_spacing);
  repair_micro_slivers(result.mesh, grid_spacing);
  if (result.mesh.triangles.empty()) {
    result.status = "empty_analytic_surface";
  } else {
    orient_mesh_consistently(result.mesh);
    result.trace_spaces = build_closed_surface_trace_spaces(result.mesh);
    result.status = result.trace_spaces.quality.failure_reason;
    if (!validate_topology && result.status == "pass")
      result.trace_spaces.quality.failure_reason = "pass";
  }
  result.build_ms = std::chrono::duration<double, std::milli>(
                        std::chrono::steady_clock::now() - begin)
                        .count();
  return result;
}

} // namespace TitanPBE::Tabi::calderon
