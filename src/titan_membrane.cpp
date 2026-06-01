#include "titan_membrane.h"

#include "BoxMatrix.h"
#include "info.h"
#include "titan_executor.h"
#include "titan_ui.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <functional>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unordered_set>

namespace {

constexpr double kEpsilon = 1.0e-12;

struct Point2D {
  double x{0.0};
  double y{0.0};
};

std::string join_output_path(const std::string &base_dir,
                             const std::string &filename) {
  if (base_dir.empty() || base_dir == ".")
    return filename;
  return (std::filesystem::path(base_dir) / filename).string();
}

std::string trim_copy(std::string value) {
  const auto first = value.find_first_not_of(" \t\r\n");
  if (first == std::string::npos)
    return {};
  const auto last = value.find_last_not_of(" \t\r\n");
  return value.substr(first, last - first + 1);
}

std::string uppercase_copy(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(),
                 [](unsigned char c) {
                   return static_cast<char>(std::toupper(c));
                 });
  return value;
}

std::vector<std::string> split(const std::string &text, char delimiter) {
  std::vector<std::string> values;
  std::stringstream stream(text);
  std::string value;
  while (std::getline(stream, value, delimiter)) {
    value = trim_copy(value);
    if (!value.empty())
      values.push_back(value);
  }
  return values;
}

double dot(const std::array<double, 3> &a, const std::array<double, 3> &b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

std::array<double, 3> cross(const std::array<double, 3> &a,
                            const std::array<double, 3> &b) {
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
          a[0] * b[1] - a[1] * b[0]};
}

double norm(const std::array<double, 3> &a) { return std::sqrt(dot(a, a)); }

std::array<double, 3> normalized(std::array<double, 3> a) {
  const double length = norm(a);
  if (!(length > kEpsilon))
    throw std::runtime_error("Membrane: degenerate normal vector.");
  for (double &value : a)
    value /= length;
  return a;
}

void apply_mic(const MOL &frame, double &dx, double &dy, double &dz) {
  if (BoxMatrix::valid_box(frame.Box)) {
    BoxMatrix matrix(frame.Box);
    matrix.apply_mic(dx, dy, dz);
  }
}

double polygon_area(const std::vector<Point2D> &polygon) {
  if (polygon.size() < 3)
    return 0.0;
  double twice_area = 0.0;
  for (std::size_t i = 0; i < polygon.size(); ++i) {
    const auto &a = polygon[i];
    const auto &b = polygon[(i + 1) % polygon.size()];
    twice_area += a.x * b.y - a.y * b.x;
  }
  return 0.5 * std::abs(twice_area);
}

std::vector<Point2D>
clip_bisector(const std::vector<Point2D> &polygon,
              double dx, double dy) {
  if (polygon.empty())
    return {};
  const double rhs = 0.5 * (dx * dx + dy * dy);
  auto inside = [&](const Point2D &point) {
    return point.x * dx + point.y * dy <= rhs + 1.0e-10;
  };
  auto intersection = [&](const Point2D &a, const Point2D &b) {
    const double da = a.x * dx + a.y * dy - rhs;
    const double db = b.x * dx + b.y * dy - rhs;
    const double denom = da - db;
    const double t = std::abs(denom) <= kEpsilon ? 0.0 : da / denom;
    return Point2D{a.x + t * (b.x - a.x), a.y + t * (b.y - a.y)};
  };

  std::vector<Point2D> result;
  result.reserve(polygon.size() + 1);
  for (std::size_t i = 0; i < polygon.size(); ++i) {
    const auto &current = polygon[i];
    const auto &next = polygon[(i + 1) % polygon.size()];
    const bool current_inside = inside(current);
    const bool next_inside = inside(next);
    if (current_inside)
      result.push_back(current);
    if (current_inside != next_inside)
      result.push_back(intersection(current, next));
  }
  return result;
}

std::array<double, 3>
smallest_covariance_axis(const AlignedDoubleVector &x,
                         const AlignedDoubleVector &y,
                         const AlignedDoubleVector &z) {
  if (x.size() < 4 || y.size() != x.size() || z.size() != x.size())
    throw std::runtime_error(
        "Membrane: at least four lipid headgroups are required.");
  const double inv_n = 1.0 / static_cast<double>(x.size());
  const double cx = std::accumulate(x.begin(), x.end(), 0.0) * inv_n;
  const double cy = std::accumulate(y.begin(), y.end(), 0.0) * inv_n;
  const double cz = std::accumulate(z.begin(), z.end(), 0.0) * inv_n;
  double a[3][3]{};
  for (std::size_t i = 0; i < x.size(); ++i) {
    const double dx = x[i] - cx;
    const double dy = y[i] - cy;
    const double dz = z[i] - cz;
    a[0][0] += dx * dx;
    a[0][1] += dx * dy;
    a[0][2] += dx * dz;
    a[1][1] += dy * dy;
    a[1][2] += dy * dz;
    a[2][2] += dz * dz;
  }
  a[1][0] = a[0][1];
  a[2][0] = a[0][2];
  a[2][1] = a[1][2];
  double vectors[3][3]{{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
  for (int iter = 0; iter < 24; ++iter) {
    int p = 0;
    int q = 1;
    double largest = std::abs(a[p][q]);
    for (int i = 0; i < 3; ++i) {
      for (int j = i + 1; j < 3; ++j) {
        if (std::abs(a[i][j]) > largest) {
          largest = std::abs(a[i][j]);
          p = i;
          q = j;
        }
      }
    }
    if (largest <= 1.0e-12)
      break;
    const double angle = 0.5 * std::atan2(2.0 * a[p][q], a[q][q] - a[p][p]);
    const double c = std::cos(angle);
    const double s = std::sin(angle);
    for (int k = 0; k < 3; ++k) {
      const double apk = a[p][k];
      const double aqk = a[q][k];
      a[p][k] = c * apk - s * aqk;
      a[q][k] = s * apk + c * aqk;
    }
    for (int k = 0; k < 3; ++k) {
      const double akp = a[k][p];
      const double akq = a[k][q];
      a[k][p] = c * akp - s * akq;
      a[k][q] = s * akp + c * akq;
      const double vkp = vectors[k][p];
      const double vkq = vectors[k][q];
      vectors[k][p] = c * vkp - s * vkq;
      vectors[k][q] = s * vkp + c * vkq;
    }
  }
  int smallest = 0;
  if (a[1][1] < a[smallest][smallest])
    smallest = 1;
  if (a[2][2] < a[smallest][smallest])
    smallest = 2;
  return normalized(
      {vectors[0][smallest], vectors[1][smallest], vectors[2][smallest]});
}

std::pair<double, double> tangent_lengths(const MOL &frame,
                                         const std::array<double, 3> &normal,
                                         const std::array<double, 3> &u,
                                         const std::array<double, 3> &v) {
  BoxMatrix matrix(frame.Box);
  const std::array<std::array<double, 3>, 3> edges{{
      {matrix.m[0][0], matrix.m[0][1], matrix.m[0][2]},
      {matrix.m[1][0], matrix.m[1][1], matrix.m[1][2]},
      {matrix.m[2][0], matrix.m[2][1], matrix.m[2][2]},
  }};
  double surface_area = 0.0;
  for (int i = 0; i < 3; ++i)
    for (int j = i + 1; j < 3; ++j)
      surface_area = std::max(surface_area, std::abs(dot(cross(edges[i], edges[j]),
                                                    normal)));
  double lu = 0.0;
  double lv = 0.0;
  for (const auto &edge : edges) {
    lu = std::max(lu, std::abs(dot(edge, u)));
    lv = std::max(lv, std::abs(dot(edge, v)));
  }
  if (!(surface_area > kEpsilon) || !(lu > kEpsilon) || !(lv > kEpsilon))
    throw std::runtime_error("Membrane: invalid projected periodic surface.");
  const double scale = std::sqrt(surface_area / (lu * lv));
  return {lu * scale, lv * scale};
}

} // namespace

void TitanMembraneAnalyzer::configure(const MembraneArgs &args) {
  args_ = args;
  requested_ = args.flag_membrane;
  enabled_ = requested_ && initialized_;
}

void TitanMembraneAnalyzer::init(const Multi_MOL &lipids, std::fstream &log) {
  if (!requested_)
    return;
  build_topology_cache(lipids, log);
  if (lipid_plans_.size() < 4)
    throw std::runtime_error(
        "Membrane: fewer than four supported lipids were selected.");
  const MOL *reference_frame =
      lipids.m_mol.empty() ? nullptr : &lipids.m_mol.front();
  std::tie(grid_x_, grid_y_) = resolved_grid(reference_frame);
  if (args_.leaflet == "user")
    load_user_leaflets();
  flip_states_.assign(lipid_plans_.size(), {});
  thickness_cells_.assign(static_cast<std::size_t>(grid_x_ * grid_y_), {});
  initialized_ = true;
  enabled_ = true;
  reset_trajectory_state();
  if (log.is_open()) {
    TitanUI::log_kv(log, "Membrane Lipids",
                    std::to_string(lipid_plans_.size()));
    TitanUI::log_kv(log, "Membrane Grid",
                    std::to_string(grid_x_) + "," + std::to_string(grid_y_));
    TitanUI::log_kv(log, "Membrane Normal", args_.normal);
  }
}

void TitanMembraneAnalyzer::reset_trajectory_state() {
  frame_summaries_.clear();
  apl_rows_.clear();
  leaflet_rows_.clear();
  order_accumulators_.clear();
  flip_events_.clear();
  std::fill(flip_states_.begin(), flip_states_.end(), FlipState{});
  std::fill(thickness_cells_.begin(), thickness_cells_.end(), ThicknessCell{});
  skipped_frames_ = 0;
  have_previous_normal_ = false;
}

void TitanMembraneAnalyzer::build_topology_cache(const Multi_MOL &lipids,
                                                 std::fstream &log) {
  lipid_plans_.clear();
  unknown_lipids_ = 0;
  const std::unordered_set<std::string> supported{
      "DOPA", "DOPC", "DOPE", "DOPG", "DOPS", "DPPC", "DPPE", "DPPG",
      "DPPS", "POPC", "POPE", "POPG", "POPS", "DLPC", "DMPC", "CHOL",
      "CHL1"};
  std::unordered_set<std::string> selected;
  if (args_.lipids.rfind("resname:", 0) == 0) {
    for (std::string value : split(args_.lipids.substr(8), ','))
      selected.insert(uppercase_copy(value));
  }
  const std::vector<std::string> explicit_heads =
      args_.head.rfind("atomname:", 0) == 0 ? split(args_.head.substr(9), ',')
                                            : std::vector<std::string>{};
  const std::vector<std::string> explicit_tail_groups =
      args_.tail.rfind("atomname:", 0) == 0 ? split(args_.tail.substr(9), ';')
                                            : std::vector<std::string>{};

  struct ResidueAtoms {
    std::string chain;
    std::string resname;
    std::size_t resnr{0};
    std::vector<std::size_t> atoms;
  };
  std::vector<ResidueAtoms> residues;
  std::unordered_map<std::string, std::size_t> residue_slots;
  for (std::size_t atom_idx = 0; atom_idx < lipids.info.size(); ++atom_idx) {
    const mol_info &atom = lipids.info[atom_idx];
    const std::string resname = uppercase_copy(atom.get_original_res_name());
    const std::string key =
        atom.get_chain() + ":" + std::to_string(atom.get_resnr()) + ":" + resname;
    auto [it, inserted] = residue_slots.emplace(key, residues.size());
    if (inserted)
      residues.push_back({atom.get_chain(), resname, atom.get_resnr(), {}});
    residues[it->second].atoms.push_back(atom_idx);
  }

  for (const ResidueAtoms &residue : residues) {
    const bool accepted = selected.empty() ? supported.contains(residue.resname)
                                           : selected.contains(residue.resname);
    if (!accepted) {
      ++unknown_lipids_;
      continue;
    }
    LipidPlan plan;
    plan.chain = residue.chain;
    plan.resname = residue.resname;
    plan.resnr = residue.resnr;
    plan.id = plan.chain + ":" + std::to_string(plan.resnr) + ":" + plan.resname;
    plan.head_atom = std::numeric_limits<std::size_t>::max();

    auto atom_index_named = [&](std::string_view name) {
      for (std::size_t idx : residue.atoms)
        if (uppercase_copy(lipids.info[idx].get_atom()) ==
            uppercase_copy(std::string(name)))
          return idx;
      return std::numeric_limits<std::size_t>::max();
    };
    if (!explicit_heads.empty()) {
      for (const std::string &name : explicit_heads) {
        plan.head_atom = atom_index_named(name);
        if (plan.head_atom != std::numeric_limits<std::size_t>::max())
          break;
      }
    } else {
      plan.head_atom = atom_index_named("PO4");
      if (plan.head_atom == std::numeric_limits<std::size_t>::max())
        plan.head_atom = atom_index_named("ROH");
    }
    if (plan.head_atom == std::numeric_limits<std::size_t>::max()) {
      ++unknown_lipids_;
      continue;
    }

    if (!explicit_tail_groups.empty()) {
      for (const std::string &group : explicit_tail_groups) {
        std::vector<std::size_t> tail;
        for (const std::string &name : split(group, ',')) {
          const std::size_t idx = atom_index_named(name);
          if (idx != std::numeric_limits<std::size_t>::max())
            tail.push_back(idx);
        }
        if (tail.size() >= 2)
          plan.tails.push_back(std::move(tail));
      }
    } else {
      std::map<char, std::vector<std::pair<int, std::size_t>>> tail_atoms;
      for (std::size_t idx : residue.atoms) {
        const std::string name = uppercase_copy(lipids.info[idx].get_atom());
        if (name.size() < 3)
          continue;
        const char suffix = name.back();
        if (suffix != 'A' && suffix != 'B')
          continue;
        int position = -1;
        for (char c : name)
          if (std::isdigit(static_cast<unsigned char>(c)) != 0) {
            position = c - '0';
            break;
          }
        if (position >= 0)
          tail_atoms[suffix].push_back({position, idx});
      }
      for (auto &[suffix, atoms] : tail_atoms) {
        (void)suffix;
        std::sort(atoms.begin(), atoms.end());
        std::vector<std::size_t> tail;
        for (const auto &[position, idx] : atoms) {
          (void)position;
          tail.push_back(idx);
        }
        if (tail.size() >= 2)
          plan.tails.push_back(std::move(tail));
      }
    }
    lipid_plans_.push_back(std::move(plan));
  }
  if (unknown_lipids_ > 0) {
    const std::string warning =
        "Membrane: skipped " + std::to_string(unknown_lipids_) +
        " residue(s) without a supported lipid template/head bead.";
    TitanUI::print_warning(warning);
    if (log.is_open())
      TitanUI::log_kv(log, "Membrane Warning", warning);
  }
}

std::array<double, 3>
TitanMembraneAnalyzer::resolve_normal(const MOL &) const {
  std::array<double, 3> normal{};
  if (args_.normal == "x")
    normal = {1.0, 0.0, 0.0};
  else if (args_.normal == "y")
    normal = {0.0, 1.0, 0.0};
  else if (args_.normal == "z")
    normal = {0.0, 0.0, 1.0};
  else
    normal = smallest_covariance_axis(head_x_, head_y_, head_z_);
  if (have_previous_normal_ && dot(normal, previous_normal_) < 0.0)
    for (double &value : normal)
      value = -value;
  previous_normal_ = normal;
  have_previous_normal_ = true;
  return normal;
}

void TitanMembraneAnalyzer::build_tangent_basis(
    const std::array<double, 3> &normal, std::array<double, 3> &u,
    std::array<double, 3> &v) const noexcept {
  const std::array<double, 3> seed =
      std::abs(normal[2]) < 0.9 ? std::array<double, 3>{0.0, 0.0, 1.0}
                                : std::array<double, 3>{1.0, 0.0, 0.0};
  u = cross(seed, normal);
  const double u_norm = norm(u);
  for (double &value : u)
    value /= u_norm;
  v = cross(normal, u);
}

void TitanMembraneAnalyzer::classify_leaflets(
    const MOL &frame, const std::array<double, 3> &normal,
    std::vector<int> &leaflets, AlignedDoubleVector &signed_distance) const {
  std::array<double, 3> center{0.0, 0.0, 0.0};
  for (std::size_t i = 0; i < head_x_.size(); ++i) {
    center[0] += head_x_[i];
    center[1] += head_y_[i];
    center[2] += head_z_[i];
  }
  for (double &value : center)
    value /= static_cast<double>(head_x_.size());
  leaflets.resize(lipid_plans_.size());
  signed_distance.resize(lipid_plans_.size());

  const bool use_local_midplane =
      args_.leaflet == "auto" && args_.normal == "auto";
  std::array<double, 3> u{}, v{};
  double lu = 0.0;
  double lv = 0.0;
  std::vector<std::vector<std::size_t>> cells;
  if (use_local_midplane) {
    build_tangent_basis(normal, u, v);
    std::tie(lu, lv) = tangent_lengths(frame, normal, u, v);
    cells.resize(static_cast<std::size_t>(grid_x_ * grid_y_));
    auto periodic_cell = [](double value, double length, int bins) {
      value = std::fmod(value, length);
      if (value < 0.0)
        value += length;
      return std::min(bins - 1,
                      static_cast<int>(value / length * bins));
    };
    for (std::size_t i = 0; i < lipid_plans_.size(); ++i) {
      const std::array<double, 3> head{head_x_[i], head_y_[i], head_z_[i]};
      const int ix = periodic_cell(dot(head, u), lu, grid_x_);
      const int iy = periodic_cell(dot(head, v), lv, grid_y_);
      cells[static_cast<std::size_t>(iy * grid_x_ + ix)].push_back(i);
    }
  }

  for (std::size_t i = 0; i < lipid_plans_.size(); ++i) {
    const auto user = user_leaflets_.find(lipid_plans_[i].id);
    if (args_.leaflet == "user" && user != user_leaflets_.end()) {
      leaflets[i] = user->second;
      signed_distance[i] = 0.0;
      continue;
    }
    std::array<double, 3> local_center = center;
    std::array<double, 3> local_normal = normal;
    if (use_local_midplane) {
      auto periodic_cell = [](double value, double length, int bins) {
        value = std::fmod(value, length);
        if (value < 0.0)
          value += length;
        return std::min(bins - 1,
                        static_cast<int>(value / length * bins));
      };
      const std::array<double, 3> anchor{head_x_[i], head_y_[i], head_z_[i]};
      const int ix = periodic_cell(dot(anchor, u), lu, grid_x_);
      const int iy = periodic_cell(dot(anchor, v), lv, grid_y_);
      AlignedDoubleVector local_x, local_y, local_z;
      std::vector<bool> included(lipid_plans_.size(), false);
      for (int dy = -1; dy <= 1; ++dy) {
        const int cy = (iy + dy + grid_y_) % grid_y_;
        for (int dx = -1; dx <= 1; ++dx) {
          const int cx = (ix + dx + grid_x_) % grid_x_;
          for (const std::size_t j :
               cells[static_cast<std::size_t>(cy * grid_x_ + cx)]) {
            if (included[j])
              continue;
            included[j] = true;
            double rx = head_x_[j] - anchor[0];
            double ry = head_y_[j] - anchor[1];
            double rz = head_z_[j] - anchor[2];
            apply_mic(frame, rx, ry, rz);
            local_x.push_back(anchor[0] + rx);
            local_y.push_back(anchor[1] + ry);
            local_z.push_back(anchor[2] + rz);
          }
        }
      }
      if (local_x.size() >= 4) {
        local_center = {
            std::accumulate(local_x.begin(), local_x.end(), 0.0) /
                static_cast<double>(local_x.size()),
            std::accumulate(local_y.begin(), local_y.end(), 0.0) /
                static_cast<double>(local_y.size()),
            std::accumulate(local_z.begin(), local_z.end(), 0.0) /
                static_cast<double>(local_z.size())};
        local_normal = smallest_covariance_axis(local_x, local_y, local_z);
        if (dot(local_normal, normal) < 0.0)
          for (double &value : local_normal)
            value = -value;
      }
    }
    std::array<double, 3> local_delta{head_x_[i] - local_center[0],
                                      head_y_[i] - local_center[1],
                                      head_z_[i] - local_center[2]};
    apply_mic(frame, local_delta[0], local_delta[1], local_delta[2]);
    signed_distance[i] = dot(local_delta, local_normal);
    leaflets[i] = signed_distance[i] >= 0.0 ? 1 : -1;
  }
}

void TitanMembraneAnalyzer::compute_apl(
    const MOL &frame, const std::array<double, 3> &normal,
    const std::vector<int> &leaflets, TitanExecutor &executor,
    std::vector<double> &apl, double &surface_area) const {
  std::array<double, 3> u{}, v{};
  build_tangent_basis(normal, u, v);
  const auto [lu, lv] = tangent_lengths(frame, normal, u, v);
  surface_area = lu * lv;
  std::vector<Point2D> projected(lipid_plans_.size());
  std::size_t projected_count = 0;
#if TITAN_COMPILED_WITH_AVX
  const __m256d ux = _mm256_set1_pd(u[0]);
  const __m256d uy = _mm256_set1_pd(u[1]);
  const __m256d uz = _mm256_set1_pd(u[2]);
  const __m256d vx = _mm256_set1_pd(v[0]);
  const __m256d vy = _mm256_set1_pd(v[1]);
  const __m256d vz = _mm256_set1_pd(v[2]);
  for (; projected_count + 4 <= lipid_plans_.size(); projected_count += 4) {
    const __m256d x = _mm256_load_pd(head_x_.data() + projected_count);
    const __m256d y = _mm256_load_pd(head_y_.data() + projected_count);
    const __m256d z = _mm256_load_pd(head_z_.data() + projected_count);
    const __m256d pu = simd_compat::fmadd_pd(
        ux, x, simd_compat::fmadd_pd(uy, y, _mm256_mul_pd(uz, z)));
    const __m256d pv = simd_compat::fmadd_pd(
        vx, x, simd_compat::fmadd_pd(vy, y, _mm256_mul_pd(vz, z)));
    alignas(32) double pu_lane[4];
    alignas(32) double pv_lane[4];
    _mm256_store_pd(pu_lane, pu);
    _mm256_store_pd(pv_lane, pv);
    for (std::size_t lane = 0; lane < 4; ++lane)
      projected[projected_count + lane] = {pu_lane[lane], pv_lane[lane]};
  }
#endif
  for (std::size_t i = projected_count; i < lipid_plans_.size(); ++i)
    projected[i] = {dot({head_x_[i], head_y_[i], head_z_[i]}, u),
                    dot({head_x_[i], head_y_[i], head_z_[i]}, v)};
  apl.assign(lipid_plans_.size(), 0.0);
  executor.parallel_for(lipid_plans_.size(), executor.capacity(), 8,
                        [&](std::size_t, std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      std::vector<Point2D> polygon{{-0.5 * lu, -0.5 * lv},
                                   {0.5 * lu, -0.5 * lv},
                                   {0.5 * lu, 0.5 * lv},
                                   {-0.5 * lu, 0.5 * lv}};
      for (std::size_t j = 0; j < projected.size() && !polygon.empty(); ++j) {
        if (leaflets[j] != leaflets[i])
          continue;
        for (int sx = -1; sx <= 1 && !polygon.empty(); ++sx) {
          for (int sy = -1; sy <= 1 && !polygon.empty(); ++sy) {
            const double dx =
                projected[j].x - projected[i].x + static_cast<double>(sx) * lu;
            const double dy =
                projected[j].y - projected[i].y + static_cast<double>(sy) * lv;
            if (dx * dx + dy * dy <= kEpsilon)
              continue;
            polygon = clip_bisector(polygon, dx, dy);
          }
        }
      }
      apl[i] = polygon_area(polygon);
    }
  });
}

void TitanMembraneAnalyzer::update_order_parameters(
    const MOL &frame, const std::array<double, 3> &normal) {
  for (const LipidPlan &lipid : lipid_plans_) {
    for (std::size_t tail_idx = 0; tail_idx < lipid.tails.size(); ++tail_idx) {
      const auto &tail = lipid.tails[tail_idx];
      for (std::size_t bond = 1; bond < tail.size(); ++bond) {
        const std::size_t a = tail[bond - 1];
        const std::size_t b = tail[bond];
        double dx = frame.pos.x[b] - frame.pos.x[a];
        double dy = frame.pos.y[b] - frame.pos.y[a];
        double dz = frame.pos.z[b] - frame.pos.z[a];
        apply_mic(frame, dx, dy, dz);
        const double value =
            order_parameter(dx, dy, dz, normal[0], normal[1], normal[2]);
        if (!std::isfinite(value))
          continue;
        const std::string key = lipid.resname + ",tail" +
                                std::to_string(tail_idx + 1) + ",bond" +
                                std::to_string(bond);
        auto &acc = order_accumulators_[key];
        acc.sum += value;
        acc.sum_sq += value * value;
        ++acc.count;
      }
    }
  }
}

void TitanMembraneAnalyzer::compute_interdigitation(
    const MOL &frame, const std::array<double, 3> &normal,
    const std::vector<int> &leaflets, double &overlap, double &penetration) const {
  std::vector<double> upper;
  std::vector<double> lower;
  for (std::size_t lipid_idx = 0; lipid_idx < lipid_plans_.size(); ++lipid_idx) {
    for (const auto &tail : lipid_plans_[lipid_idx].tails) {
      for (std::size_t atom : tail) {
        const double projected =
            frame.pos.x[atom] * normal[0] + frame.pos.y[atom] * normal[1] +
            frame.pos.z[atom] * normal[2];
        (leaflets[lipid_idx] > 0 ? upper : lower).push_back(projected);
      }
    }
  }
  if (upper.empty() || lower.empty()) {
    overlap = 0.0;
    penetration = 0.0;
    return;
  }
  const auto [upper_min, upper_max] =
      std::minmax_element(upper.begin(), upper.end());
  const auto [lower_min, lower_max] =
      std::minmax_element(lower.begin(), lower.end());
  penetration = std::max(0.0, std::min(*upper_max, *lower_max) -
                                  std::max(*upper_min, *lower_min));
  const double total_span =
      std::max(*upper_max, *lower_max) - std::min(*upper_min, *lower_min);
  overlap = total_span > kEpsilon ? std::clamp(penetration / total_span, 0.0, 1.0)
                                  : 0.0;
}

void TitanMembraneAnalyzer::update_thickness_map(
    const MOL &frame, const std::array<double, 3> &normal,
    const std::vector<int> &leaflets, double &mean_thickness) {
  std::array<double, 3> u{}, v{};
  build_tangent_basis(normal, u, v);
  const auto [lu, lv] = tangent_lengths(frame, normal, u, v);
  const std::size_t cell_count = static_cast<std::size_t>(grid_x_ * grid_y_);
  std::vector<double> upper_sum(cell_count, 0.0), lower_sum(cell_count, 0.0);
  std::vector<std::size_t> upper_count(cell_count, 0), lower_count(cell_count, 0);
  for (std::size_t i = 0; i < lipid_plans_.size(); ++i) {
    const std::array<double, 3> point{head_x_[i], head_y_[i], head_z_[i]};
    double pu = dot(point, u);
    double pv = dot(point, v);
    pu -= std::floor(pu / lu) * lu;
    pv -= std::floor(pv / lv) * lv;
    const int ix = std::min(grid_x_ - 1, static_cast<int>(pu / lu * grid_x_));
    const int iy = std::min(grid_y_ - 1, static_cast<int>(pv / lv * grid_y_));
    const std::size_t cell = static_cast<std::size_t>(iy * grid_x_ + ix);
    const double height = dot(point, normal);
    if (leaflets[i] > 0) {
      upper_sum[cell] += height;
      ++upper_count[cell];
    } else {
      lower_sum[cell] += height;
      ++lower_count[cell];
    }
  }
  double thickness_sum = 0.0;
  std::size_t thickness_count = 0;
  for (std::size_t cell = 0; cell < cell_count; ++cell) {
    if (upper_count[cell] == 0 || lower_count[cell] == 0)
      continue;
    const double thickness = upper_sum[cell] / upper_count[cell] -
                             lower_sum[cell] / lower_count[cell];
    thickness_cells_[cell].sum += thickness;
    thickness_cells_[cell].sum_sq += thickness * thickness;
    ++thickness_cells_[cell].count;
    thickness_sum += thickness;
    ++thickness_count;
  }
  if (thickness_count == 0) {
    double upper = 0.0;
    double lower = 0.0;
    std::size_t nu = 0;
    std::size_t nl = 0;
    for (std::size_t i = 0; i < lipid_plans_.size(); ++i) {
      const double height = head_x_[i] * normal[0] + head_y_[i] * normal[1] +
                            head_z_[i] * normal[2];
      if (leaflets[i] > 0) {
        upper += height;
        ++nu;
      } else {
        lower += height;
        ++nl;
      }
    }
    if (nu == 0 || nl == 0)
      throw std::runtime_error("Membrane: could not resolve both leaflets.");
    mean_thickness = upper / nu - lower / nl;
    return;
  }
  mean_thickness = thickness_sum / static_cast<double>(thickness_count);
}

void TitanMembraneAnalyzer::update_flipflop(const std::vector<int> &leaflets,
                                            std::size_t frame,
                                            double time_ps) {
  for (std::size_t i = 0; i < leaflets.size(); ++i) {
    FlipState &state = flip_states_[i];
    const int current = leaflets[i];
    if (state.stable == 0) {
      state.stable = current;
      continue;
    }
    if (current == state.stable) {
      state.pending = 0;
      state.pending_count = 0;
      continue;
    }
    if (state.pending != current) {
      state.pending = current;
      state.pending_count = 1;
    } else {
      ++state.pending_count;
    }
    if (state.pending_count >= args_.flipflop_persistence) {
      flip_events_.push_back(
          {frame, time_ps, lipid_plans_[i].id, leaflet_name(state.stable),
           leaflet_name(current)});
      state.stable = current;
      state.pending = 0;
      state.pending_count = 0;
    }
  }
}

void TitanMembraneAnalyzer::process_frame(const Multi_MOL &lipids,
                                          std::size_t snap_idx,
                                          TitanExecutor &executor) {
  if (!enabled_ || snap_idx >= lipids.m_mol.size())
    return;
  const MOL &frame = lipids.m_mol[snap_idx];
  if (!frame.coordinates_consistent() || !BoxMatrix::valid_box(frame.Box)) {
    ++skipped_frames_;
    return;
  }
  head_x_.resize(lipid_plans_.size());
  head_y_.resize(lipid_plans_.size());
  head_z_.resize(lipid_plans_.size());
  for (std::size_t i = 0; i < lipid_plans_.size(); ++i) {
    const std::size_t head = lipid_plans_[i].head_atom;
    if (head >= frame.size())
      throw std::runtime_error("Membrane: head-bead index exceeds frame size.");
    head_x_[i] = frame.pos.x[head];
    head_y_[i] = frame.pos.y[head];
    head_z_[i] = frame.pos.z[head];
  }
  const std::array<double, 3> normal = resolve_normal(frame);
  std::vector<int> leaflets;
  AlignedDoubleVector signed_distance;
  classify_leaflets(frame, normal, leaflets, signed_distance);
  std::size_t upper_count = 0;
  std::size_t lower_count = 0;
  for (int leaflet : leaflets)
    leaflet > 0 ? ++upper_count : ++lower_count;
  if (upper_count == 0 || lower_count == 0)
    throw std::runtime_error("Membrane: could not resolve both leaflets.");

  std::vector<double> apl;
  double surface_area = 0.0;
  compute_apl(frame, normal, leaflets, executor, apl, surface_area);
  const double mean_apl =
      std::accumulate(apl.begin(), apl.end(), 0.0) / apl.size();
  double thickness = 0.0;
  update_thickness_map(frame, normal, leaflets, thickness);
  double interdigitation = 0.0;
  double penetration = 0.0;
  compute_interdigitation(frame, normal, leaflets, interdigitation, penetration);
  update_order_parameters(frame, normal);
  update_flipflop(leaflets, snap_idx, frame.get_time());

  frame_summaries_.push_back({snap_idx, frame.get_time(), mean_apl, thickness,
                              interdigitation, penetration, upper_count,
                              lower_count});
  for (std::size_t i = 0; i < lipid_plans_.size(); ++i) {
    const std::string leaflet = leaflet_name(leaflets[i]);
    apl_rows_.push_back({snap_idx, frame.get_time(), lipid_plans_[i].id,
                         leaflet, apl[i]});
    leaflet_rows_.push_back({snap_idx, frame.get_time(), lipid_plans_[i].id,
                             leaflet, signed_distance[i]});
  }
  (void)surface_area;
}

double TitanMembraneAnalyzer::order_parameter(double dx, double dy, double dz,
                                              double nx, double ny,
                                              double nz) noexcept {
  const double bond_norm = std::sqrt(dx * dx + dy * dy + dz * dz);
  const double normal_norm = std::sqrt(nx * nx + ny * ny + nz * nz);
  if (!(bond_norm > kEpsilon) || !(normal_norm > kEpsilon))
    return std::numeric_limits<double>::quiet_NaN();
  const double cosine = (dx * nx + dy * ny + dz * nz) / (bond_norm * normal_norm);
  return 0.5 * (3.0 * cosine * cosine - 1.0);
}

std::pair<int, int>
TitanMembraneAnalyzer::resolved_grid(const MOL *frame) const {
  if (args_.grid != "auto") {
    const auto values = split(args_.grid, ',');
    return {std::stoi(values.at(0)), std::stoi(values.at(1))};
  }
  const int cells = std::max(
      1, static_cast<int>(std::ceil(lipid_plans_.size() / 16.0)));
  double aspect = 1.0;
  if (frame != nullptr && BoxMatrix::valid_box(frame->Box)) {
    BoxMatrix matrix(frame->Box);
    std::array<double, 3> lengths{
        norm({matrix.m[0][0], matrix.m[0][1], matrix.m[0][2]}),
        norm({matrix.m[1][0], matrix.m[1][1], matrix.m[1][2]}),
        norm({matrix.m[2][0], matrix.m[2][1], matrix.m[2][2]})};
    std::sort(lengths.begin(), lengths.end(), std::greater<>());
    if (lengths[1] > kEpsilon)
      aspect = lengths[0] / lengths[1];
  }
  const int nx = std::clamp(
      static_cast<int>(std::ceil(std::sqrt(cells * aspect))), 4, 64);
  const int ny = std::clamp(
      static_cast<int>(std::ceil(static_cast<double>(cells) / nx)), 4, 64);
  return {nx, ny};
}

std::string TitanMembraneAnalyzer::leaflet_name(int leaflet) const {
  return leaflet > 0 ? "upper" : "lower";
}

void TitanMembraneAnalyzer::load_user_leaflets() {
  user_leaflets_.clear();
  std::ifstream input(args_.leaflet_file);
  if (!input.is_open())
    throw std::runtime_error("Membrane: cannot open leaflet CSV: " +
                             args_.leaflet_file);
  std::string line;
  while (std::getline(input, line)) {
    const auto values = split(line, ',');
    if (values.size() < 2 || values[0] == "lipid_id")
      continue;
    const std::string leaflet = uppercase_copy(values[1]);
    if (leaflet != "UPPER" && leaflet != "LOWER")
      throw std::runtime_error("Membrane: invalid user leaflet value.");
    user_leaflets_[values[0]] = leaflet == "UPPER" ? 1 : -1;
  }
  for (const auto &lipid : lipid_plans_)
    if (!user_leaflets_.contains(lipid.id))
      throw std::runtime_error("Membrane: user leaflet CSV misses lipid " +
                               lipid.id);
}

void TitanMembraneAnalyzer::export_leaflets(const std::string &filename) const {
  std::ofstream out(filename);
  out << "frame,time_ps,lipid_id,leaflet,signed_distance_A\n";
  for (const auto &row : leaflet_rows_)
    out << row.frame << ',' << row.time_ps << ',' << row.lipid_id << ','
        << row.leaflet << ',' << row.signed_distance_a << '\n';
}

void TitanMembraneAnalyzer::export_apl_timeseries(
    const std::string &filename) const {
  std::ofstream out(filename);
  out << "frame,time_ps,mean_apl_A2,upper_lipids,lower_lipids\n";
  for (const auto &row : frame_summaries_)
    out << row.frame << ',' << row.time_ps << ',' << row.mean_apl_a2 << ','
        << row.upper_count << ',' << row.lower_count << '\n';
}

void TitanMembraneAnalyzer::export_apl_per_lipid(
    const std::string &filename) const {
  std::ofstream out(filename);
  out << "frame,time_ps,lipid_id,leaflet,apl_A2\n";
  for (const auto &row : apl_rows_)
    out << row.frame << ',' << row.time_ps << ',' << row.lipid_id << ','
        << row.leaflet << ',' << row.apl_a2 << '\n';
}

void TitanMembraneAnalyzer::export_thickness_timeseries(
    const std::string &filename) const {
  std::ofstream out(filename);
  out << "frame,time_ps,mean_thickness_A\n";
  for (const auto &row : frame_summaries_)
    out << row.frame << ',' << row.time_ps << ',' << row.mean_thickness_a
        << '\n';
}

void TitanMembraneAnalyzer::export_thickness_map(
    const std::string &filename) const {
  std::ofstream out(filename);
  out << "grid_x,grid_y,mean_thickness_A,samples\n";
  for (int iy = 0; iy < grid_y_; ++iy) {
    for (int ix = 0; ix < grid_x_; ++ix) {
      const auto &cell = thickness_cells_[static_cast<std::size_t>(
          iy * grid_x_ + ix)];
      out << ix << ',' << iy << ',';
      if (cell.count == 0)
        out << "nan,0\n";
      else
        out << cell.sum / cell.count << ',' << cell.count << '\n';
    }
  }
}

void TitanMembraneAnalyzer::export_order_parameters(
    const std::string &filename) const {
  std::ofstream out(filename);
  out << "resname,tail,bond,mean_S,std_S,samples\n";
  for (const auto &[key, acc] : order_accumulators_) {
    const double mean = acc.count == 0 ? 0.0 : acc.sum / acc.count;
    const double variance =
        acc.count == 0 ? 0.0 : std::max(0.0, acc.sum_sq / acc.count - mean * mean);
    out << key << ',' << mean << ',' << std::sqrt(variance) << ',' << acc.count
        << '\n';
  }
}

void TitanMembraneAnalyzer::export_interdigitation(
    const std::string &filename) const {
  std::ofstream out(filename);
  out << "frame time_ps overlap_coefficient penetration_depth_A\n";
  for (const auto &row : frame_summaries_)
    out << row.frame << ' ' << row.time_ps << ' ' << row.interdigitation << ' '
        << row.penetration_depth_a << '\n';
}

void TitanMembraneAnalyzer::export_flipflop_events(
    const std::string &filename) const {
  std::ofstream out(filename);
  out << "frame,time_ps,lipid_id,from_leaflet,to_leaflet\n";
  for (const auto &row : flip_events_)
    out << row.frame << ',' << row.time_ps << ',' << row.lipid_id << ','
        << row.from_leaflet << ',' << row.to_leaflet << '\n';
}

void TitanMembraneAnalyzer::export_summary(const std::string &filename) const {
  std::ofstream out(filename);
  out << std::fixed << std::setprecision(8);
  out << "module: MARTINI membrane analysis\n";
  out << "reference_note: Leaflets use head-bead geometry; auto normal uses PCA "
         "with sign continuity between analyzed frames, and auto leaflets use "
         "periodic local-neighborhood PCA midplanes.\n";
  out << "apl_note: APL is a periodic local Voronoi estimate per leaflet.\n";
  out << "order_parameter_note: MARTINI tail order uses 0.5*(3*cos(theta)^2-1) "
         "for consecutive tail-bead vectors.\n";
  out << "frames_analyzed: " << frame_summaries_.size() << '\n';
  out << "frames_skipped: " << skipped_frames_ << '\n';
  out << "lipids: " << lipid_plans_.size() << '\n';
  out << "unknown_lipids_skipped: " << unknown_lipids_ << '\n';
  out << "normal_policy: " << args_.normal << '\n';
  out << "leaflet_policy: " << args_.leaflet << '\n';
  out << "grid: " << grid_x_ << ',' << grid_y_ << '\n';
  out << "flipflop_persistence_frames: " << args_.flipflop_persistence << '\n';
  if (!frame_summaries_.empty()) {
    double apl = 0.0;
    double thickness = 0.0;
    double interdigitation = 0.0;
    for (const auto &row : frame_summaries_) {
      apl += row.mean_apl_a2;
      thickness += row.mean_thickness_a;
      interdigitation += row.interdigitation;
    }
    const double count = static_cast<double>(frame_summaries_.size());
    out << "mean_apl_A2: " << apl / count << '\n';
    out << "mean_thickness_A: " << thickness / count << '\n';
    out << "mean_interdigitation: " << interdigitation / count << '\n';
  }
  out << "flipflop_events: " << flip_events_.size() << '\n';
}

void TitanMembraneAnalyzer::accumulate_trajectory() {
  if (frame_summaries_.empty())
    return;
  auto append = [](ThicknessCell &acc, double value) {
    acc.sum += value;
    acc.sum_sq += value * value;
    ++acc.count;
  };
  double apl = 0.0;
  double thickness = 0.0;
  double interdigitation = 0.0;
  for (const auto &row : frame_summaries_) {
    apl += row.mean_apl_a2;
    thickness += row.mean_thickness_a;
    interdigitation += row.interdigitation;
  }
  const double frame_count = static_cast<double>(frame_summaries_.size());
  append(aggregate_apl_, apl / frame_count);
  append(aggregate_thickness_, thickness / frame_count);
  append(aggregate_interdigitation_, interdigitation / frame_count);

  if (aggregate_thickness_cells_.size() < thickness_cells_.size())
    aggregate_thickness_cells_.resize(thickness_cells_.size());
  for (std::size_t i = 0; i < thickness_cells_.size(); ++i) {
    if (thickness_cells_[i].count == 0)
      continue;
    append(aggregate_thickness_cells_[i],
           thickness_cells_[i].sum / thickness_cells_[i].count);
  }
  for (const auto &[key, acc] : order_accumulators_) {
    if (acc.count == 0)
      continue;
    append(aggregate_order_parameters_[key], acc.sum / acc.count);
  }
  ++aggregate_trajectories_;
}

void TitanMembraneAnalyzer::export_average(const std::string &base_dir,
                                           const info &info_output) const {
  auto mean_sem = [](const ThicknessCell &acc) {
    if (acc.count == 0)
      return std::pair{std::numeric_limits<double>::quiet_NaN(),
                       std::numeric_limits<double>::quiet_NaN()};
    const double mean = acc.sum / acc.count;
    const double variance =
        std::max(0.0, acc.sum_sq / acc.count - mean * mean);
    const double sem = acc.count > 1 ? std::sqrt(variance / acc.count) : 0.0;
    return std::pair{mean, sem};
  };

  {
    std::ofstream out(join_output_path(
        base_dir, info_output.get_membrane_summary_average_filename()));
    out << std::fixed << std::setprecision(8);
    out << "module: MARTINI membrane replicate aggregate\n";
    out << "reference_note: Aggregate values are replicate means with SEM; "
           "trajectory timeseries remain separate.\n";
    out << "replicates: " << aggregate_trajectories_ << '\n';
    const auto [apl, apl_sem] = mean_sem(aggregate_apl_);
    const auto [thickness, thickness_sem] = mean_sem(aggregate_thickness_);
    const auto [interdigitation, interdigitation_sem] =
        mean_sem(aggregate_interdigitation_);
    out << "mean_apl_A2: " << apl << '\n';
    out << "sem_apl_A2: " << apl_sem << '\n';
    out << "mean_thickness_A: " << thickness << '\n';
    out << "sem_thickness_A: " << thickness_sem << '\n';
    out << "mean_interdigitation: " << interdigitation << '\n';
    out << "sem_interdigitation: " << interdigitation_sem << '\n';
  }
  {
    std::ofstream out(join_output_path(
        base_dir, info_output.get_membrane_thickness_map_average_filename()));
    out << "grid_x,grid_y,mean_thickness_A,sem_thickness_A,replicates\n";
    for (int iy = 0; iy < grid_y_; ++iy) {
      for (int ix = 0; ix < grid_x_; ++ix) {
        const auto &acc = aggregate_thickness_cells_[static_cast<std::size_t>(
            iy * grid_x_ + ix)];
        const auto [mean, sem] = mean_sem(acc);
        out << ix << ',' << iy << ',' << mean << ',' << sem << ',' << acc.count
            << '\n';
      }
    }
  }
  {
    std::ofstream out(join_output_path(
        base_dir, info_output.get_membrane_order_parameters_average_filename()));
    out << "resname,tail,bond,mean_S,sem_S,replicates\n";
    for (const auto &[key, acc] : aggregate_order_parameters_) {
      const auto [mean, sem] = mean_sem(acc);
      out << key << ',' << mean << ',' << sem << ',' << acc.count << '\n';
    }
  }
}

void TitanMembraneAnalyzer::finalize_trajectory(
    const std::string &base_dir, const info &info_output, const arguments &args,
    const std::vector<std::string> &input_line, int input_name,
    std::fstream &log) {
  if (!enabled_)
    return;
  export_leaflets(join_output_path(
      base_dir, info_output.get_membrane_leaflets_filename(args, input_line,
                                                           input_name)));
  export_apl_timeseries(join_output_path(
      base_dir, info_output.get_membrane_apl_timeseries_filename(
                    args, input_line, input_name)));
  export_apl_per_lipid(join_output_path(
      base_dir, info_output.get_membrane_apl_per_lipid_filename(
                    args, input_line, input_name)));
  export_thickness_timeseries(join_output_path(
      base_dir, info_output.get_membrane_thickness_timeseries_filename(
                    args, input_line, input_name)));
  export_thickness_map(join_output_path(
      base_dir, info_output.get_membrane_thickness_map_filename(
                    args, input_line, input_name)));
  export_order_parameters(join_output_path(
      base_dir, info_output.get_membrane_order_parameters_filename(
                    args, input_line, input_name)));
  export_interdigitation(join_output_path(
      base_dir, info_output.get_membrane_interdigitation_filename(
                    args, input_line, input_name)));
  export_flipflop_events(join_output_path(
      base_dir, info_output.get_membrane_flipflop_events_filename(
                    args, input_line, input_name)));
  export_summary(join_output_path(
      base_dir,
      info_output.get_membrane_summary_filename(args, input_line, input_name)));
  accumulate_trajectory();
  const bool multi_input = input_line.size() > 1 || args.flag_listfile ||
                           args.flag_pattern_input || args.flag_xdr;
  if (multi_input)
    export_average(base_dir, info_output);
  if (log.is_open())
    TitanUI::log_kv(log, "Membrane Frames",
                    std::to_string(frame_summaries_.size()));
}
