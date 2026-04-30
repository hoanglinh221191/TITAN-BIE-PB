#include "titan_sasa.h"

#include "titan_executor.h"
#include "titan_simd_dispatch.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>

namespace TitanSASA {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

struct SpherePoint {
  double x{0.0};
  double y{0.0};
  double z{0.0};
};

const std::vector<SpherePoint> &sphere_points(int point_count) {
  static const std::vector<SpherePoint> default_points = [] {
    std::vector<SpherePoint> out(kDefaultPointCount);
    const double inc = kPi * (3.0 - std::sqrt(5.0));
    const double off = 2.0 / static_cast<double>(kDefaultPointCount);
    for (int k = 0; k < kDefaultPointCount; ++k) {
      const double y = static_cast<double>(k) * off - 1.0 + 0.5 * off;
      const double r = std::sqrt(std::max(0.0, 1.0 - y * y));
      const double phi = static_cast<double>(k) * inc;
      out[static_cast<std::size_t>(k)] = {std::cos(phi) * r, y,
                                          std::sin(phi) * r};
    }
    return out;
  }();
  if (point_count == kDefaultPointCount)
    return default_points;

  static thread_local std::vector<SpherePoint> custom_points;
  static thread_local int custom_count = 0;
  if (custom_count != point_count) {
    custom_points.assign(static_cast<std::size_t>(point_count), {});
    const double inc = kPi * (3.0 - std::sqrt(5.0));
    const double off = 2.0 / static_cast<double>(point_count);
    for (int k = 0; k < point_count; ++k) {
      const double y = static_cast<double>(k) * off - 1.0 + 0.5 * off;
      const double r = std::sqrt(std::max(0.0, 1.0 - y * y));
      const double phi = static_cast<double>(k) * inc;
      custom_points[static_cast<std::size_t>(k)] = {std::cos(phi) * r, y,
                                                    std::sin(phi) * r};
    }
    custom_count = point_count;
  }
  return custom_points;
}

std::size_t effective_workers(const TitanExecutor *executor,
                              std::size_t work_items) noexcept {
  if (executor == nullptr || work_items == 0)
    return 1;
  return std::max<std::size_t>(
      1, std::min<std::size_t>(executor->capacity(), work_items));
}

} // namespace

std::vector<double>
compute_per_atom_sasa(std::span<const Atom> atoms, const Options &options,
                      std::span<const std::size_t> targets,
                      TitanExecutor *executor) {
  std::vector<double> area_by_atom(atoms.size(), 0.0);
  if (atoms.empty())
    return area_by_atom;

  const int point_count =
      options.point_count > 0 ? options.point_count : kDefaultPointCount;
  const auto &sphere = sphere_points(point_count);
  const double probe = std::max(0.0, options.probe_radius);

  double min_x = std::numeric_limits<double>::max();
  double min_y = std::numeric_limits<double>::max();
  double min_z = std::numeric_limits<double>::max();
  double max_x = std::numeric_limits<double>::lowest();
  double max_y = std::numeric_limits<double>::lowest();
  double max_z = std::numeric_limits<double>::lowest();
  double max_surface_radius = 0.0;
  for (const Atom &atom : atoms) {
    const double surface_radius = std::max(0.0, atom.radius) + probe;
    max_surface_radius = std::max(max_surface_radius, surface_radius);
    min_x = std::min(min_x, atom.x);
    min_y = std::min(min_y, atom.y);
    min_z = std::min(min_z, atom.z);
    max_x = std::max(max_x, atom.x);
    max_y = std::max(max_y, atom.y);
    max_z = std::max(max_z, atom.z);
  }

  const double cell =
      std::max(options.min_cell_size, 2.0 * max_surface_radius);
  const std::size_t bx = std::max<std::size_t>(
      1, static_cast<std::size_t>(std::floor((max_x - min_x) / cell)) + 1);
  const std::size_t by = std::max<std::size_t>(
      1, static_cast<std::size_t>(std::floor((max_y - min_y) / cell)) + 1);
  const std::size_t bz = std::max<std::size_t>(
      1, static_cast<std::size_t>(std::floor((max_z - min_z) / cell)) + 1);

  std::vector<std::vector<std::size_t>> bins(bx * by * bz);
  auto bin_index = [&](double x, double y, double z) {
    const std::size_t ix = std::min<std::size_t>(
        bx - 1, static_cast<std::size_t>(
                    std::max(0.0, std::floor((x - min_x) / cell))));
    const std::size_t iy = std::min<std::size_t>(
        by - 1, static_cast<std::size_t>(
                    std::max(0.0, std::floor((y - min_y) / cell))));
    const std::size_t iz = std::min<std::size_t>(
        bz - 1, static_cast<std::size_t>(
                    std::max(0.0, std::floor((z - min_z) / cell))));
    return (iz * by + iy) * bx + ix;
  };

  for (std::size_t i = 0; i < atoms.size(); ++i)
    bins[bin_index(atoms[i].x, atoms[i].y, atoms[i].z)].push_back(i);

  std::vector<std::size_t> all_targets;
  if (targets.empty()) {
    all_targets.resize(atoms.size());
    std::iota(all_targets.begin(), all_targets.end(), std::size_t{0});
    targets = std::span<const std::size_t>(all_targets.data(),
                                           all_targets.size());
  }

  auto compute_range = [&](std::size_t begin, std::size_t end) {
    for (std::size_t target_pos = begin; target_pos < end; ++target_pos) {
      const std::size_t i = targets[target_pos];
      if (i >= atoms.size())
        continue;
      const Atom &atom = atoms[i];
      const double ri = std::max(0.0, atom.radius) + probe;
      if (ri <= 0.0)
        continue;

      const double point_area =
          4.0 * kPi * ri * ri / static_cast<double>(point_count);
      const std::size_t base = bin_index(atom.x, atom.y, atom.z);
      const std::size_t bix = base % bx;
      const std::size_t biy = (base / bx) % by;
      const std::size_t biz = base / (bx * by);
      const std::size_t x0 = (bix == 0) ? 0 : bix - 1;
      const std::size_t y0 = (biy == 0) ? 0 : biy - 1;
      const std::size_t z0 = (biz == 0) ? 0 : biz - 1;
      const std::size_t x1 = std::min(bx - 1, bix + 1);
      const std::size_t y1 = std::min(by - 1, biy + 1);
      const std::size_t z1 = std::min(bz - 1, biz + 1);

      std::vector<std::size_t> neighbors;
      neighbors.reserve(64);
      for (std::size_t kk = z0; kk <= z1; ++kk) {
        for (std::size_t jj = y0; jj <= y1; ++jj) {
          for (std::size_t ii = x0; ii <= x1; ++ii) {
            const auto &cell_atoms = bins[(kk * by + jj) * bx + ii];
            for (std::size_t j : cell_atoms) {
              if (j == i)
                continue;
              const Atom &other = atoms[j];
              const double rj = std::max(0.0, other.radius) + probe;
              const double dx = atom.x - other.x;
              const double dy = atom.y - other.y;
              const double dz = atom.z - other.z;
              if (dx * dx + dy * dy + dz * dz < (ri + rj) * (ri + rj))
                neighbors.push_back(j);
            }
          }
        }
      }

      int exposed = 0;
      const std::size_t full_groups = sphere.size() / 4;
      std::vector<unsigned> buried_mask(full_groups, 0u);
      std::vector<std::array<double, 4>> px(full_groups);
      std::vector<std::array<double, 4>> py(full_groups);
      std::vector<std::array<double, 4>> pz(full_groups);
      for (std::size_t g = 0; g < full_groups; ++g) {
        for (int lane = 0; lane < 4; ++lane) {
          const SpherePoint &pt = sphere[g * 4 + static_cast<std::size_t>(lane)];
          px[g][static_cast<std::size_t>(lane)] = atom.x + pt.x * ri;
          py[g][static_cast<std::size_t>(lane)] = atom.y + pt.y * ri;
          pz[g][static_cast<std::size_t>(lane)] = atom.z + pt.z * ri;
        }
      }

      for (std::size_t j : neighbors) {
        const Atom &other = atoms[j];
        const double rj = std::max(0.0, other.radius) + probe;
        const double cutoff_sq = rj * rj;
        for (std::size_t g = 0; g < full_groups; ++g) {
          buried_mask[g] |= titan_simd_dispatch::point_vs_four_cutoff_mask(
              other.x, other.y, other.z, px[g].data(), py[g].data(),
              pz[g].data(), cutoff_sq, -1.0, false);
        }
      }

      for (std::size_t g = 0; g < full_groups; ++g) {
        const unsigned mask = buried_mask[g];
        if ((mask & 1u) == 0u)
          ++exposed;
        if ((mask & 2u) == 0u)
          ++exposed;
        if ((mask & 4u) == 0u)
          ++exposed;
        if ((mask & 8u) == 0u)
          ++exposed;
      }

      const std::size_t remainder_begin = full_groups * 4;
      for (std::size_t p = remainder_begin; p < sphere.size(); ++p) {
        const SpherePoint &pt = sphere[p];
        const double sx = atom.x + pt.x * ri;
        const double sy = atom.y + pt.y * ri;
        const double sz = atom.z + pt.z * ri;
        bool buried = false;
        for (std::size_t j : neighbors) {
          const Atom &other = atoms[j];
          const double rj = std::max(0.0, other.radius) + probe;
          const double dx = sx - other.x;
          const double dy = sy - other.y;
          const double dz = sz - other.z;
          if (dx * dx + dy * dy + dz * dz < rj * rj) {
            buried = true;
            break;
          }
        }
        if (!buried)
          ++exposed;
      }

      area_by_atom[i] = static_cast<double>(exposed) * point_area;
    }
  };

  if (executor != nullptr && targets.size() > 1) {
    const std::size_t workers = effective_workers(executor, targets.size());
    executor->parallel_for(targets.size(), workers, 8,
                           [&](std::size_t worker, std::size_t begin,
                               std::size_t end) {
                             (void)worker;
                             compute_range(begin, end);
                           });
  } else {
    compute_range(0, targets.size());
  }

  return area_by_atom;
}

double sum_sasa(std::span<const double> sasa) noexcept {
  return std::accumulate(sasa.begin(), sasa.end(), 0.0);
}

double nonpolar_energy_from_delta_sasa_kj(double delta_sasa_a2) noexcept {
  return kNonpolarGammaKJPerMolA2 * delta_sasa_a2;
}

} // namespace TitanSASA
