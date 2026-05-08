#include "tabi_types.h"

#include "../titan_surface.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <unordered_map>
#include <vector>

namespace TitanPBE::Tabi {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kTwoPi = 2.0 * kPi;
constexpr double kVertexQuantize = 1.0e8;

struct UnitVertex {
  double x{0.0};
  double y{0.0};
  double z{0.0};
};

struct UnitTriangle {
  std::size_t a{0};
  std::size_t b{0};
  std::size_t c{0};
};

struct UnitSphereMesh {
  std::vector<UnitVertex> vertices;
  std::vector<UnitTriangle> triangles;
};

struct AtomCellGrid {
  double min_x{0.0};
  double min_y{0.0};
  double min_z{0.0};
  double cell{1.0};
  std::size_t nx{1};
  std::size_t ny{1};
  std::size_t nz{1};
  std::size_t cell_count{1};
};

struct QuantizedKey {
  std::int64_t x{0};
  std::int64_t y{0};
  std::int64_t z{0};

  bool operator==(const QuantizedKey &other) const noexcept {
    return x == other.x && y == other.y && z == other.z;
  }
};

struct QuantizedKeyHash {
  std::size_t operator()(const QuantizedKey &key) const noexcept {
    std::uint64_t h = 1469598103934665603ull;
    auto mix = [&](std::int64_t value) {
      h ^= static_cast<std::uint64_t>(value) + 0x9e3779b97f4a7c15ull +
           (h << 6u) + (h >> 2u);
      h *= 1099511628211ull;
    };
    mix(key.x);
    mix(key.y);
    mix(key.z);
    return static_cast<std::size_t>(h);
  }
};

struct EdgeKey {
  std::size_t a{0};
  std::size_t b{0};

  bool operator==(const EdgeKey &other) const noexcept {
    return a == other.a && b == other.b;
  }
};

struct EdgeKeyHash {
  std::size_t operator()(const EdgeKey &edge) const noexcept {
    std::uint64_t h = 1469598103934665603ull;
    h ^= static_cast<std::uint64_t>(edge.a) + 0x9e3779b97f4a7c15ull +
         (h << 6u) + (h >> 2u);
    h *= 1099511628211ull;
    h ^= static_cast<std::uint64_t>(edge.b) + 0x9e3779b97f4a7c15ull +
         (h << 6u) + (h >> 2u);
    h *= 1099511628211ull;
    return static_cast<std::size_t>(h);
  }
};

bool finite_atom(const Atom &atom) {
  return std::isfinite(atom.x) && std::isfinite(atom.y) &&
         std::isfinite(atom.z) && std::isfinite(atom.radius);
}

double elapsed_ms(std::chrono::steady_clock::time_point begin,
                  std::chrono::steady_clock::time_point end) noexcept {
  return std::chrono::duration<double, std::milli>(end - begin).count();
}

double norm3(double x, double y, double z) noexcept {
  return std::sqrt(std::max(0.0, x * x + y * y + z * z));
}

double spherical_triangle_area_unit(const UnitVertex &a, const UnitVertex &b,
                                    const UnitVertex &c) noexcept {
  const double bx_cx = b.y * c.z - b.z * c.y;
  const double by_cx = b.z * c.x - b.x * c.z;
  const double bz_cx = b.x * c.y - b.y * c.x;
  const double det = a.x * bx_cx + a.y * by_cx + a.z * bz_cx;
  const double ab = a.x * b.x + a.y * b.y + a.z * b.z;
  const double bc = b.x * c.x + b.y * c.y + b.z * c.z;
  const double ca = c.x * a.x + c.y * a.y + c.z * a.z;
  const double denom = 1.0 + ab + bc + ca;
  return std::abs(2.0 * std::atan2(det, denom));
}

QuantizedKey quantized_key(double x, double y, double z) noexcept {
  return {static_cast<std::int64_t>(std::llround(x * kVertexQuantize)),
          static_cast<std::int64_t>(std::llround(y * kVertexQuantize)),
          static_cast<std::int64_t>(std::llround(z * kVertexQuantize))};
}

int mesh_rings_for_point_count(int point_count) noexcept {
  const int safe_points = std::max(24, point_count);
  const int rings =
      static_cast<int>(std::lround(std::sqrt(static_cast<double>(safe_points) *
                                             0.75))) +
      1;
  return std::clamp(rings, 6, 18);
}

UnitSphereMesh make_unit_sphere_mesh(int point_count) {
  const int rings = mesh_rings_for_point_count(point_count);
  const int segments = std::max(12, rings * 2);
  UnitSphereMesh mesh;
  mesh.vertices.reserve(2 + static_cast<std::size_t>(rings - 1) *
                                static_cast<std::size_t>(segments));
  mesh.vertices.push_back({0.0, 0.0, 1.0});
  for (int ring = 1; ring < rings; ++ring) {
    const double theta = kPi * static_cast<double>(ring) /
                         static_cast<double>(rings);
    const double z = std::cos(theta);
    const double r = std::sin(theta);
    for (int seg = 0; seg < segments; ++seg) {
      const double phi = kTwoPi * static_cast<double>(seg) /
                         static_cast<double>(segments);
      mesh.vertices.push_back({r * std::cos(phi), r * std::sin(phi), z});
    }
  }
  const std::size_t south = mesh.vertices.size();
  mesh.vertices.push_back({0.0, 0.0, -1.0});

  auto ring_index = [&](int ring, int seg) -> std::size_t {
    const int wrapped = (seg % segments + segments) % segments;
    return 1 + static_cast<std::size_t>(ring - 1) *
                   static_cast<std::size_t>(segments) +
           static_cast<std::size_t>(wrapped);
  };

  for (int seg = 0; seg < segments; ++seg)
    mesh.triangles.push_back({0, ring_index(1, seg),
                              ring_index(1, seg + 1)});

  for (int ring = 1; ring < rings - 1; ++ring) {
    for (int seg = 0; seg < segments; ++seg) {
      const std::size_t a = ring_index(ring, seg);
      const std::size_t b = ring_index(ring, seg + 1);
      const std::size_t c = ring_index(ring + 1, seg);
      const std::size_t d = ring_index(ring + 1, seg + 1);
      mesh.triangles.push_back({a, c, b});
      mesh.triangles.push_back({b, c, d});
    }
  }

  for (int seg = 0; seg < segments; ++seg)
    mesh.triangles.push_back({ring_index(rings - 1, seg), south,
                              ring_index(rings - 1, seg + 1)});
  return mesh;
}

const UnitSphereMesh &cached_unit_sphere_mesh(int point_count) {
  static const UnitSphereMesh default_mesh =
      make_unit_sphere_mesh(TitanSurface::kDefaultPointCount);
  if (point_count == TitanSurface::kDefaultPointCount)
    return default_mesh;
  static thread_local UnitSphereMesh local_mesh;
  static thread_local int local_points = 0;
  if (local_points != point_count) {
    local_mesh = make_unit_sphere_mesh(point_count);
    local_points = point_count;
  }
  return local_mesh;
}

AtomCellGrid make_atom_cell_grid(const std::vector<Atom> &atoms,
                                 const std::vector<double> &surface_radii) {
  AtomCellGrid grid;
  if (atoms.empty()) {
    return grid;
  }

  double max_radius = 0.0;
  double max_x = -std::numeric_limits<double>::infinity();
  double max_y = -std::numeric_limits<double>::infinity();
  double max_z = -std::numeric_limits<double>::infinity();
  grid.min_x = std::numeric_limits<double>::infinity();
  grid.min_y = std::numeric_limits<double>::infinity();
  grid.min_z = std::numeric_limits<double>::infinity();
  for (std::size_t i = 0; i < atoms.size(); ++i) {
    max_radius = std::max(max_radius, surface_radii[i]);
    grid.min_x = std::min(grid.min_x, atoms[i].x);
    grid.min_y = std::min(grid.min_y, atoms[i].y);
    grid.min_z = std::min(grid.min_z, atoms[i].z);
    max_x = std::max(max_x, atoms[i].x);
    max_y = std::max(max_y, atoms[i].y);
    max_z = std::max(max_z, atoms[i].z);
  }

  grid.cell = std::max(1.0e-6, 2.0 * max_radius);
  grid.nx = std::max<std::size_t>(
      1, static_cast<std::size_t>(
             std::floor((max_x - grid.min_x) / grid.cell)) +
             1);
  grid.ny = std::max<std::size_t>(
      1, static_cast<std::size_t>(
             std::floor((max_y - grid.min_y) / grid.cell)) +
             1);
  grid.nz = std::max<std::size_t>(
      1, static_cast<std::size_t>(
             std::floor((max_z - grid.min_z) / grid.cell)) +
             1);
  grid.cell_count = grid.nx * grid.ny * grid.nz;
  return grid;
}

std::size_t atom_cell_index(const AtomCellGrid &grid, double x, double y,
                            double z) noexcept {
  const std::size_t ix = std::min<std::size_t>(
      grid.nx - 1, static_cast<std::size_t>(
                       std::max(0.0, std::floor((x - grid.min_x) / grid.cell))));
  const std::size_t iy = std::min<std::size_t>(
      grid.ny - 1, static_cast<std::size_t>(
                       std::max(0.0, std::floor((y - grid.min_y) / grid.cell))));
  const std::size_t iz = std::min<std::size_t>(
      grid.nz - 1, static_cast<std::size_t>(
                       std::max(0.0, std::floor((z - grid.min_z) / grid.cell))));
  return (iz * grid.ny + iy) * grid.nx + ix;
}

void atom_cell_coords(const AtomCellGrid &grid, std::size_t cell,
                      std::size_t &ix, std::size_t &iy,
                      std::size_t &iz) noexcept {
  ix = cell % grid.nx;
  const std::size_t yz = cell / grid.nx;
  iy = yz % grid.ny;
  iz = yz / grid.ny;
}

std::vector<std::vector<std::size_t>>
build_atom_neighbor_lists(const std::vector<Atom> &atoms,
                          const std::vector<double> &surface_radii,
                          SurfaceMesh &mesh) {
  const auto begin = std::chrono::steady_clock::now();
  std::vector<std::vector<std::size_t>> neighbors(atoms.size());
  if (atoms.empty()) {
    mesh.neighbor_build_ms = elapsed_ms(begin, std::chrono::steady_clock::now());
    return neighbors;
  }

  const AtomCellGrid grid = make_atom_cell_grid(atoms, surface_radii);
  mesh.cell_count = grid.cell_count;
  std::vector<std::size_t> cell_offsets(grid.cell_count + 1, 0);
  std::vector<std::size_t> atom_cell(atoms.size(), 0);
  for (std::size_t i = 0; i < atoms.size(); ++i) {
    atom_cell[i] = atom_cell_index(grid, atoms[i].x, atoms[i].y, atoms[i].z);
    ++cell_offsets[atom_cell[i] + 1];
  }
  for (std::size_t c = 0; c < grid.cell_count; ++c)
    cell_offsets[c + 1] += cell_offsets[c];

  std::vector<std::size_t> fill = cell_offsets;
  std::vector<std::size_t> cell_atoms(atoms.size(), 0);
  for (std::size_t i = 0; i < atoms.size(); ++i)
    cell_atoms[fill[atom_cell[i]]++] = i;

  std::size_t neighbor_sum = 0;
  for (std::size_t i = 0; i < atoms.size(); ++i) {
    std::size_t ix = 0;
    std::size_t iy = 0;
    std::size_t iz = 0;
    atom_cell_coords(grid, atom_cell[i], ix, iy, iz);
    const std::size_t x0 = ix == 0 ? 0 : ix - 1;
    const std::size_t y0 = iy == 0 ? 0 : iy - 1;
    const std::size_t z0 = iz == 0 ? 0 : iz - 1;
    const std::size_t x1 = std::min(grid.nx - 1, ix + 1);
    const std::size_t y1 = std::min(grid.ny - 1, iy + 1);
    const std::size_t z1 = std::min(grid.nz - 1, iz + 1);
    auto &list = neighbors[i];
    for (std::size_t zz = z0; zz <= z1; ++zz) {
      for (std::size_t yy = y0; yy <= y1; ++yy) {
        for (std::size_t xx = x0; xx <= x1; ++xx) {
          const std::size_t cell = (zz * grid.ny + yy) * grid.nx + xx;
          for (std::size_t pos = cell_offsets[cell];
               pos < cell_offsets[cell + 1]; ++pos) {
            const std::size_t j = cell_atoms[pos];
            if (i == j)
              continue;
            const double cutoff = surface_radii[i] + surface_radii[j];
            const double dx = atoms[i].x - atoms[j].x;
            const double dy = atoms[i].y - atoms[j].y;
            const double dz = atoms[i].z - atoms[j].z;
            if (dx * dx + dy * dy + dz * dz < cutoff * cutoff)
              list.push_back(j);
          }
        }
      }
    }
    std::sort(list.begin(), list.end());
    list.erase(std::unique(list.begin(), list.end()), list.end());
    neighbor_sum += list.size();
    mesh.max_neighbor_count = std::max(mesh.max_neighbor_count, list.size());
  }

  mesh.mean_neighbor_count =
      atoms.empty() ? 0.0
                    : static_cast<double>(neighbor_sum) /
                          static_cast<double>(atoms.size());
  mesh.neighbor_build_ms = elapsed_ms(begin, std::chrono::steady_clock::now());
  return neighbors;
}

bool point_buried(double x, double y, double z,
                  const std::vector<Atom> &atoms,
                  const std::vector<double> &surface_radii,
                  const std::vector<std::size_t> &neighbors) noexcept {
  for (const std::size_t j : neighbors) {
    const Atom &other = atoms[j];
    const double dx = x - other.x;
    const double dy = y - other.y;
    const double dz = z - other.z;
    const double r = surface_radii[j];
    if (dx * dx + dy * dy + dz * dz < r * r)
      return true;
  }
  return false;
}

double triangle_solid_angle(const Vertex &a, const Vertex &b, const Vertex &c,
                            double px, double py, double pz) noexcept {
  const double ax = a.x - px;
  const double ay = a.y - py;
  const double az = a.z - pz;
  const double bx = b.x - px;
  const double by = b.y - py;
  const double bz = b.z - pz;
  const double cx = c.x - px;
  const double cy = c.y - py;
  const double cz = c.z - pz;
  const double la = norm3(ax, ay, az);
  const double lb = norm3(bx, by, bz);
  const double lc = norm3(cx, cy, cz);
  if (la <= 1.0e-14 || lb <= 1.0e-14 || lc <= 1.0e-14)
    return 0.0;
  const double cross_x = by * cz - bz * cy;
  const double cross_y = bz * cx - bx * cz;
  const double cross_z = bx * cy - by * cx;
  const double numerator = ax * cross_x + ay * cross_y + az * cross_z;
  const double ab = ax * bx + ay * by + az * bz;
  const double bc = bx * cx + by * cy + bz * cz;
  const double ca = cx * ax + cy * ay + cz * az;
  const double denominator = la * lb * lc + ab * lc + bc * la + ca * lb;
  return 2.0 * std::atan2(numerator, denominator);
}

std::size_t add_vertex(
    SurfaceMesh &mesh,
    std::unordered_map<QuantizedKey, std::size_t, QuantizedKeyHash> &lookup,
    const Atom &atom, std::size_t atom_index, const UnitVertex &unit,
    double radius) {
  const double x = atom.x + unit.x * radius;
  const double y = atom.y + unit.y * radius;
  const double z = atom.z + unit.z * radius;
  const QuantizedKey key = quantized_key(x, y, z);
  auto found = lookup.find(key);
  if (found != lookup.end())
    return found->second;
  const std::size_t index = mesh.vertices.size();
  mesh.vertices.push_back({x, y, z, unit.x, unit.y, unit.z, atom_index});
  lookup.emplace(key, index);
  return index;
}

void add_edge(std::unordered_map<EdgeKey, std::size_t, EdgeKeyHash> &edges,
              std::size_t a, std::size_t b) {
  if (a > b)
    std::swap(a, b);
  ++edges[{a, b}];
}

void finalize_solid_angle_diagonal(SurfaceMesh &mesh,
                                   const std::string &diagonal_mode) {
  const auto begin = std::chrono::steady_clock::now();
  if (mesh.panels.empty()) {
    mesh.diag_min = 1.0;
    mesh.diag_max = 1.0;
    mesh.solid_angle_ms = elapsed_ms(begin, std::chrono::steady_clock::now());
    return;
  }

  std::vector<std::vector<std::size_t>> triangles_by_vertex(
      mesh.vertices.size());
  std::unordered_map<EdgeKey, std::size_t, EdgeKeyHash> edge_counts;
  edge_counts.reserve(mesh.triangles.size() * 3);
  for (std::size_t t = 0; t < mesh.triangles.size(); ++t) {
    const Triangle &tri = mesh.triangles[t];
    triangles_by_vertex[tri.v0].push_back(t);
    triangles_by_vertex[tri.v1].push_back(t);
    triangles_by_vertex[tri.v2].push_back(t);
    add_edge(edge_counts, tri.v0, tri.v1);
    add_edge(edge_counts, tri.v1, tri.v2);
    add_edge(edge_counts, tri.v2, tri.v0);
  }

  std::size_t valence_sum = 0;
  std::size_t active_vertices = 0;
  for (const auto &incident : triangles_by_vertex) {
    if (!incident.empty()) {
      valence_sum += incident.size();
      ++active_vertices;
    }
  }
  mesh.avg_valence = active_vertices > 0
                         ? static_cast<double>(valence_sum) /
                               static_cast<double>(active_vertices)
                         : 0.0;
  mesh.nonmanifold_edges = 0;
  for (const auto &[edge, count] : edge_counts) {
    (void)edge;
    if (count > 2)
      ++mesh.nonmanifold_edges;
  }

  mesh.diag_min = std::numeric_limits<double>::infinity();
  mesh.diag_max = -std::numeric_limits<double>::infinity();
  std::vector<std::vector<std::size_t>> patches(mesh.triangles.size());
  if (diagonal_mode == "solid-angle") {
    std::vector<std::size_t> patch;
    for (std::size_t t = 0; t < mesh.triangles.size(); ++t) {
      const Triangle &tri = mesh.triangles[t];
      patch.clear();
      patch.push_back(t);
      for (const std::size_t vertex : {tri.v0, tri.v1, tri.v2}) {
        const auto &incident = triangles_by_vertex[vertex];
        patch.insert(patch.end(), incident.begin(), incident.end());
      }
      std::sort(patch.begin(), patch.end());
      patch.erase(std::unique(patch.begin(), patch.end()), patch.end());
      patches[t] = patch;
    }
  }
  for (std::size_t t = 0; t < mesh.triangles.size(); ++t) {
    Panel &panel = mesh.panels[t];
    double diagonal = 1.0;
    double solid_angle = kTwoPi;
    if (diagonal_mode == "off") {
      diagonal = 0.0;
      solid_angle = 0.0;
    } else if (diagonal_mode == "solid-angle") {
      const double offset = std::max(1.0e-7, std::sqrt(panel.area) * 1.0e-5);
      const double px = panel.x + panel.nx * offset;
      const double py = panel.y + panel.ny * offset;
      const double pz = panel.z + panel.nz * offset;
      double omega = 0.0;
      for (const std::size_t patch_tri_index : patches[t]) {
        const Triangle &patch_tri = mesh.triangles[patch_tri_index];
        omega += triangle_solid_angle(mesh.vertices[patch_tri.v0],
                                      mesh.vertices[patch_tri.v1],
                                      mesh.vertices[patch_tri.v2], px, py, pz);
      }
      solid_angle = std::abs(omega);
      diagonal = solid_angle / kTwoPi;
      if (!std::isfinite(diagonal) || diagonal < 0.10 || diagonal > 2.50) {
        ++mesh.solid_angle_fallbacks;
        diagonal = 1.0;
        solid_angle = kTwoPi;
      } else {
        diagonal = std::clamp(diagonal, 0.25, 1.75);
      }
    }
    panel.diagonal = diagonal;
    panel.solid_angle = solid_angle;
    mesh.diag_min = std::min(mesh.diag_min, diagonal);
    mesh.diag_max = std::max(mesh.diag_max, diagonal);
  }
  if (!std::isfinite(mesh.diag_min))
    mesh.diag_min = 1.0;
  if (!std::isfinite(mesh.diag_max))
    mesh.diag_max = 1.0;
  mesh.solid_angle_ms = elapsed_ms(begin, std::chrono::steady_clock::now());
}

SurfaceMesh build_triangle_mesh(const std::vector<Atom> &atoms,
                                const std::vector<std::size_t> &original_indices,
                                const Options &options) {
  const double probe =
      options.surface_mode == "sas" ? std::max(0.0, options.probe_radius) : 0.0;
  const int point_count = std::max(24, TitanSurface::kDefaultPointCount);
  const UnitSphereMesh &unit_mesh = cached_unit_sphere_mesh(point_count);

  std::vector<double> surface_radii(atoms.size(), 0.0);
  for (std::size_t i = 0; i < atoms.size(); ++i)
    surface_radii[i] = std::max(0.0, atoms[i].radius) + probe;

  SurfaceMesh mesh;
  const std::vector<std::vector<std::size_t>> neighbors =
      build_atom_neighbor_lists(atoms, surface_radii, mesh);
  const auto mesh_begin = std::chrono::steady_clock::now();
  mesh.vertices.reserve(atoms.size() * unit_mesh.vertices.size());
  mesh.triangles.reserve(atoms.size() * unit_mesh.triangles.size());
  mesh.panels.reserve(atoms.size() * unit_mesh.triangles.size());
  std::unordered_map<QuantizedKey, std::size_t, QuantizedKeyHash> vertex_lookup;
  vertex_lookup.reserve(atoms.size() * unit_mesh.vertices.size());

  for (std::size_t atom_index = 0; atom_index < atoms.size(); ++atom_index) {
    const Atom &atom = atoms[atom_index];
    const double radius = surface_radii[atom_index];
    if (radius <= 0.0)
      continue;
    for (const UnitTriangle &unit_tri : unit_mesh.triangles) {
      const UnitVertex &ua = unit_mesh.vertices[unit_tri.a];
      const UnitVertex &ub = unit_mesh.vertices[unit_tri.b];
      const UnitVertex &uc = unit_mesh.vertices[unit_tri.c];
      const double ax = atom.x + ua.x * radius;
      const double ay = atom.y + ua.y * radius;
      const double az = atom.z + ua.z * radius;
      const double bx = atom.x + ub.x * radius;
      const double by = atom.y + ub.y * radius;
      const double bz = atom.z + ub.z * radius;
      const double cx = atom.x + uc.x * radius;
      const double cy = atom.y + uc.y * radius;
      const double cz = atom.z + uc.z * radius;
      double panel_ux = (ua.x + ub.x + uc.x) / 3.0;
      double panel_uy = (ua.y + ub.y + uc.y) / 3.0;
      double panel_uz = (ua.z + ub.z + uc.z) / 3.0;
      const double panel_len = norm3(panel_ux, panel_uy, panel_uz);
      if (panel_len <= 1.0e-14)
        continue;
      panel_ux /= panel_len;
      panel_uy /= panel_len;
      panel_uz /= panel_len;
      const double px = atom.x + panel_ux * radius;
      const double py = atom.y + panel_uy * radius;
      const double pz = atom.z + panel_uz * radius;

      const bool centroid_buried =
          point_buried(px, py, pz, atoms, surface_radii, neighbors[atom_index]);
      const int buried_vertices =
          static_cast<int>(point_buried(ax, ay, az, atoms, surface_radii,
                                        neighbors[atom_index])) +
          static_cast<int>(point_buried(bx, by, bz, atoms, surface_radii,
                                        neighbors[atom_index])) +
          static_cast<int>(point_buried(cx, cy, cz, atoms, surface_radii,
                                        neighbors[atom_index]));
      if (centroid_buried || buried_vertices == 3) {
        ++mesh.buried_candidates;
        continue;
      }
      ++mesh.exposed_candidates;

      std::size_t ia =
          add_vertex(mesh, vertex_lookup, atom, original_indices[atom_index], ua,
                     radius);
      std::size_t ib =
          add_vertex(mesh, vertex_lookup, atom, original_indices[atom_index], ub,
                     radius);
      std::size_t ic =
          add_vertex(mesh, vertex_lookup, atom, original_indices[atom_index], uc,
                     radius);

      const Vertex &va = mesh.vertices[ia];
      const Vertex &vb = mesh.vertices[ib];
      const Vertex &vc = mesh.vertices[ic];
      double ux = vb.x - va.x;
      double uy = vb.y - va.y;
      double uz = vb.z - va.z;
      double vx = vc.x - va.x;
      double vy = vc.y - va.y;
      double vz = vc.z - va.z;
      double nx = uy * vz - uz * vy;
      double ny = uz * vx - ux * vz;
      double nz = ux * vy - uy * vx;
      double normal_len = norm3(nx, ny, nz);
      if (normal_len <= 1.0e-14)
        continue;
      double area = spherical_triangle_area_unit(ua, ub, uc) * radius * radius;
      if (!std::isfinite(area) || area <= 1.0e-14)
        area = 0.5 * normal_len;
      nx /= normal_len;
      ny /= normal_len;
      nz /= normal_len;
      const double radial_x = px - atom.x;
      const double radial_y = py - atom.y;
      const double radial_z = pz - atom.z;
      if (nx * radial_x + ny * radial_y + nz * radial_z < 0.0) {
        std::swap(ib, ic);
        nx = -nx;
        ny = -ny;
        nz = -nz;
      }
      nx = panel_ux;
      ny = panel_uy;
      nz = panel_uz;

      Triangle tri;
      tri.v0 = ia;
      tri.v1 = ib;
      tri.v2 = ic;
      tri.atom_index = original_indices[atom_index];
      tri.cx = px;
      tri.cy = py;
      tri.cz = pz;
      tri.nx = nx;
      tri.ny = ny;
      tri.nz = nz;
      tri.area = area;
      mesh.triangles.push_back(tri);
      mesh.panels.push_back(
          Panel{px, py, pz, nx, ny, nz, area, 1.0, kTwoPi, tri.atom_index});
    }
  }

  mesh.mesh_build_ms = elapsed_ms(mesh_begin, std::chrono::steady_clock::now());
  finalize_solid_angle_diagonal(mesh, options.diagonal_mode);
  return mesh;
}

std::vector<Panel> build_point_panel_surface(const std::vector<Atom> &atoms,
                                             const Options &options,
                                             Result *result,
                                             double elapsed_surface_ms,
                                             TitanExecutor *executor) {
  std::vector<TitanSurface::Atom> surface_atoms;
  std::vector<std::size_t> original_indices;
  surface_atoms.reserve(atoms.size());
  original_indices.reserve(atoms.size());
  for (std::size_t i = 0; i < atoms.size(); ++i) {
    const Atom &atom = atoms[i];
    if (!finite_atom(atom))
      continue;
    surface_atoms.push_back({atom.x, atom.y, atom.z, atom.radius});
    original_indices.push_back(i);
  }

  TitanSurface::Options surface_options;
  surface_options.probe_radius = options.probe_radius;
  surface_options.min_cell_size = 1.0;
  surface_options.point_count = TitanSurface::kDefaultPointCount;
  surface_options.surface_mode =
      options.surface_mode == "sas" ? std::string("sas") : std::string("vdw");
  surface_options.emit_panels = true;
  surface_options.emit_neighbor_workspace = false;

  TitanSurface::Workspace workspace =
      TitanSurface::build_workspace(surface_atoms, surface_options, {},
                                    executor);

  std::vector<Panel> panels;
  panels.reserve(workspace.panels.size());
  for (const TitanSurface::Panel &surface_panel : workspace.panels) {
    Panel panel;
    panel.x = surface_panel.x;
    panel.y = surface_panel.y;
    panel.z = surface_panel.z;
    panel.nx = surface_panel.nx;
    panel.ny = surface_panel.ny;
    panel.nz = surface_panel.nz;
    panel.area = surface_panel.area;
    panel.diagonal = options.diagonal_mode == "off" ? 0.0 : 1.0;
    panel.solid_angle = panel.diagonal > 0.0 ? kTwoPi : 0.0;
    panel.atom_index = surface_panel.atom_index < original_indices.size()
                           ? original_indices[surface_panel.atom_index]
                           : surface_panel.atom_index;
    panels.push_back(panel);
  }

  if (result != nullptr) {
    result->surface_geometry = "point-panel";
    result->panel_count = panels.size();
    result->surface_area_a2 = workspace.total_area;
    result->surface_build_ms = elapsed_surface_ms + workspace.diagnostics.build_ms;
    result->surface_neighbor_ms =
        workspace.diagnostics.neighbor_build_ms;
    result->surface_mesh_ms = workspace.diagnostics.exposure_ms;
    result->solid_angle_ms = 0.0;
    result->exposed_points = workspace.diagnostics.exposed_points;
    result->buried_points = workspace.diagnostics.buried_points;
    result->boundary_nodes = panels.size();
    result->diag_min = options.diagonal_mode == "off" ? 0.0 : 1.0;
    result->diag_max = options.diagonal_mode == "off" ? 0.0 : 1.0;
  }
  return panels;
}

} // namespace

std::vector<Panel> build_surface(const std::vector<Atom> &atoms,
                                 const Options &options, Result *result,
                                 TitanExecutor *executor) {
  const auto surface_begin = std::chrono::steady_clock::now();
  if (options.surface_geometry == "point-panel") {
    return build_point_panel_surface(atoms, options, result, 0.0, executor);
  }

  std::vector<Atom> finite_atoms;
  std::vector<std::size_t> original_indices;
  finite_atoms.reserve(atoms.size());
  original_indices.reserve(atoms.size());
  for (std::size_t i = 0; i < atoms.size(); ++i) {
    if (!finite_atom(atoms[i]))
      continue;
    finite_atoms.push_back(atoms[i]);
    original_indices.push_back(i);
  }

  SurfaceMesh mesh = build_triangle_mesh(finite_atoms, original_indices, options);
  if (result != nullptr) {
    result->surface_geometry = "triangle";
    result->panel_count = mesh.panels.size();
    result->surface_area_a2 = std::accumulate(
        mesh.panels.begin(), mesh.panels.end(), 0.0,
        [](double sum, const Panel &panel) { return sum + panel.area; });
    result->surface_build_ms =
        elapsed_ms(surface_begin, std::chrono::steady_clock::now());
    result->surface_neighbor_ms = mesh.neighbor_build_ms;
    result->surface_mesh_ms = mesh.mesh_build_ms;
    result->solid_angle_ms = mesh.solid_angle_ms;
    result->exposed_points = mesh.exposed_candidates;
    result->buried_points = mesh.buried_candidates;
    result->mesh_vertices = mesh.vertices.size();
    result->mesh_triangles = mesh.triangles.size();
    result->boundary_nodes = mesh.panels.size();
    result->avg_valence = mesh.avg_valence;
    result->nonmanifold_edges = mesh.nonmanifold_edges;
    result->solid_angle_fallbacks = mesh.solid_angle_fallbacks;
    result->diag_min = mesh.diag_min;
    result->diag_max = mesh.diag_max;
  }
  return std::move(mesh.panels);
}

} // namespace TitanPBE::Tabi
