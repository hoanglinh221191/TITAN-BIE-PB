#include "titan_pbe.h"

#include "titan_energy_decomposition.h"
#include "titan_executor.h"
#include "titan_interaction_entropy.h"
#include "titan_sasa.h"
#include "titan_ui.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cctype>
#include <cstdint>
#include <filesystem>
#include <iomanip>
#include <limits>
#include <map>
#include <numeric>
#include <span>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <utility>

namespace TitanPBE {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kCoulombKcalPerMol = 332.063713299;
constexpr double kKcalToKJ = 4.184;
constexpr double kKJToKcal = 1.0 / kKcalToKJ;
constexpr double kKBKJPerMolK = 0.00831446261815324;
constexpr double kKBKcalPerMolK = kKBKJPerMolK * kKJToKcal;
constexpr double kDelphiEpk = 167100.9162872952;
constexpr double kDelphiReactionFieldEnergyFactor = 0.9549296586;

std::size_t idx3(std::size_t i, std::size_t j, std::size_t k,
                 std::size_t nx, std::size_t ny) noexcept {
  return (k * ny + j) * nx + i;
}

std::size_t effective_workers(const TitanExecutor &executor,
                              std::size_t work_items) noexcept {
  if (work_items == 0)
    return 1;
  return std::max<std::size_t>(
      1, std::min<std::size_t>(executor.capacity(), work_items));
}

double harmonic_eps(double a, double b) noexcept {
  if (a <= 0.0 || b <= 0.0)
    return 0.5 * (a + b);
  if (std::abs(a - b) <= 1.0e-14)
    return a;
  return 2.0 * a * b / (a + b);
}

double debye_kappa_inv_a(double salt_m, double temp_k) noexcept {
  if (salt_m <= 0.0)
    return 0.0;
  const double temp_scale = std::sqrt(298.15 / std::max(1.0, temp_k));
  return 0.3291 * std::sqrt(salt_m) * temp_scale;
}

std::string format_fixed(double value, int digits = 6) {
  std::ostringstream oss;
  oss << std::fixed << std::setprecision(digits) << value;
  return oss.str();
}

std::string format_suffix(std::size_t value, int width) {
  std::ostringstream oss;
  oss << std::setw(width) << std::setfill('0') << value;
  return oss.str();
}

double elapsed_ms(std::chrono::steady_clock::time_point begin,
                  std::chrono::steady_clock::time_point end) {
  return std::chrono::duration<double, std::milli>(end - begin).count();
}

std::string csv_escape(const std::string &value) {
  if (value.find_first_of(",\"\n\r") == std::string::npos)
    return value;
  std::string out = "\"";
  for (char ch : value) {
    if (ch == '"')
      out += "\"\"";
    else
      out += ch;
  }
  out += '"';
  return out;
}

template <typename Shape>
bool is_boundary(std::size_t i, std::size_t j, std::size_t k,
                 const Shape &shape) noexcept {
  return i == 0 || j == 0 || k == 0 || i + 1 == shape.nx ||
         j + 1 == shape.ny || k + 1 == shape.nz;
}

double clamp_fraction(double value) noexcept {
  if (value < 0.0)
    return 0.0;
  if (value > 1.0)
    return 1.0;
  return value;
}

bool parse_double_token(const std::string &token, double &value) {
  try {
    std::size_t used = 0;
    value = std::stod(token, &used);
    return used == token.size() && std::isfinite(value);
  } catch (...) {
    return false;
  }
}

std::string trim_copy(const std::string &value) {
  std::size_t first = 0;
  while (first < value.size() &&
         std::isspace(static_cast<unsigned char>(value[first])))
    ++first;
  std::size_t last = value.size();
  while (last > first &&
         std::isspace(static_cast<unsigned char>(value[last - 1])))
    --last;
  std::string out = value.substr(first, last - first);
  if (out.size() >= 2 && out.front() == '"' && out.back() == '"')
    out = out.substr(1, out.size() - 2);
  return out;
}

std::string lowercase_copy(std::string value) {
  for (char &ch : value)
    ch = static_cast<char>(std::tolower(static_cast<unsigned char>(ch)));
  return value;
}

std::vector<std::string> split_csv_simple(const std::string &line) {
  std::vector<std::string> fields;
  std::string field;
  bool quoted = false;
  for (std::size_t i = 0; i < line.size(); ++i) {
    const char ch = line[i];
    if (ch == '"') {
      if (quoted && i + 1 < line.size() && line[i + 1] == '"') {
        field.push_back('"');
        ++i;
      } else {
        quoted = !quoted;
      }
    } else if (ch == ',' && !quoted) {
      fields.push_back(trim_copy(field));
      field.clear();
    } else {
      field.push_back(ch);
    }
  }
  fields.push_back(trim_copy(field));
  return fields;
}

std::unordered_map<std::string, std::size_t>
csv_header_map(const std::vector<std::string> &header) {
  std::unordered_map<std::string, std::size_t> map;
  for (std::size_t i = 0; i < header.size(); ++i)
    map[lowercase_copy(header[i])] = i;
  return map;
}

bool csv_double(const std::vector<std::string> &fields,
                const std::unordered_map<std::string, std::size_t> &header,
                const std::vector<std::string> &names, double &value) {
  for (const std::string &name : names) {
    const auto it = header.find(lowercase_copy(name));
    if (it != header.end() && it->second < fields.size() &&
        parse_double_token(fields[it->second], value))
      return true;
  }
  return false;
}

bool csv_index(const std::vector<std::string> &fields,
               const std::unordered_map<std::string, std::size_t> &header,
               const std::string &name, std::size_t &value) {
  double parsed = 0.0;
  if (!csv_double(fields, header, {name}, parsed) || parsed < 0.0)
    return false;
  value = static_cast<std::size_t>(std::llround(parsed));
  return true;
}

bool parse_boolish(const std::string &value, bool &out) {
  const std::string v = lowercase_copy(trim_copy(value));
  if (v == "1" || v == "true" || v == "yes" || v == "solute" ||
      v == "inside") {
    out = true;
    return true;
  }
  if (v == "0" || v == "false" || v == "no" || v == "solvent" ||
      v == "outside") {
    out = false;
    return true;
  }
  double numeric = 0.0;
  if (parse_double_token(v, numeric)) {
    out = std::abs(numeric) > 0.5;
    return true;
  }
  return false;
}

bool surface_uses_probe_exclusion(const std::string &surface) noexcept {
  return surface == "sas" || surface == "ses" || surface == "delphi-ses";
}

bool surface_uses_ses_flood_fill(const std::string &surface) noexcept {
  return surface == "ses" || surface == "delphi-ses";
}

bool surface_uses_face_refinement(const std::string &surface) noexcept {
  return surface == "delphi-ses" || surface == "delphi-vdw" ||
         surface == "gaussian" || surface == "gaussian-cutoff";
}

bool surface_uses_gaussian(const std::string &surface) noexcept {
  return surface == "gaussian" || surface == "gaussian-cutoff";
}

std::string uppercase_letters(std::string value) {
  std::string out;
  out.reserve(value.size());
  for (char ch : value) {
    if (std::isalpha(static_cast<unsigned char>(ch)))
      out.push_back(static_cast<char>(
          std::toupper(static_cast<unsigned char>(ch))));
  }
  return out;
}

struct AtomSpatialBins {
  double min_x{0.0};
  double min_y{0.0};
  double min_z{0.0};
  double cell{1.0};
  std::size_t bx{1};
  std::size_t by{1};
  std::size_t bz{1};
  double max_radius{0.0};
  std::vector<std::vector<std::size_t>> bins;
  bool empty{true};
};

template <typename AtomRange>
AtomSpatialBins build_atom_bins(const AtomRange &atoms, double inflate,
                                double min_cell) {
  AtomSpatialBins out;
  if (atoms.empty())
    return out;

  out.empty = false;
  out.min_x = std::numeric_limits<double>::max();
  out.min_y = std::numeric_limits<double>::max();
  out.min_z = std::numeric_limits<double>::max();
  double max_x = std::numeric_limits<double>::lowest();
  double max_y = std::numeric_limits<double>::lowest();
  double max_z = std::numeric_limits<double>::lowest();
  for (const auto &atom : atoms) {
    out.max_radius = std::max(out.max_radius, atom.radius + inflate);
    out.min_x = std::min(out.min_x, atom.x);
    out.min_y = std::min(out.min_y, atom.y);
    out.min_z = std::min(out.min_z, atom.z);
    max_x = std::max(max_x, atom.x);
    max_y = std::max(max_y, atom.y);
    max_z = std::max(max_z, atom.z);
  }
  out.cell = std::max({out.max_radius, min_cell, 1.0e-6});
  out.bx = std::max<std::size_t>(
      1, static_cast<std::size_t>(std::floor((max_x - out.min_x) / out.cell)) +
             1);
  out.by = std::max<std::size_t>(
      1, static_cast<std::size_t>(std::floor((max_y - out.min_y) / out.cell)) +
             1);
  out.bz = std::max<std::size_t>(
      1, static_cast<std::size_t>(std::floor((max_z - out.min_z) / out.cell)) +
             1);
  out.bins.assign(out.bx * out.by * out.bz, {});

  auto bin_index = [&](double x, double y, double z) {
    const std::size_t ix = std::min<std::size_t>(
        out.bx - 1, static_cast<std::size_t>(
                        std::max(0.0, std::floor((x - out.min_x) / out.cell))));
    const std::size_t iy = std::min<std::size_t>(
        out.by - 1, static_cast<std::size_t>(
                        std::max(0.0, std::floor((y - out.min_y) / out.cell))));
    const std::size_t iz = std::min<std::size_t>(
        out.bz - 1, static_cast<std::size_t>(
                        std::max(0.0, std::floor((z - out.min_z) / out.cell))));
    return (iz * out.by + iy) * out.bx + ix;
  };
  for (std::size_t atom_idx = 0; atom_idx < atoms.size(); ++atom_idx)
    out.bins[bin_index(atoms[atom_idx].x, atoms[atom_idx].y,
                       atoms[atom_idx].z)]
        .push_back(atom_idx);
  return out;
}

template <typename AtomRange>
double max_atom_radius(const AtomRange &atoms) {
  double radius = 0.0;
  for (const auto &atom : atoms)
    radius = std::max(radius, atom.radius);
  return radius;
}

template <typename AtomRange>
bool point_inside_inflated_atom(const AtomSpatialBins &bins,
                                const AtomRange &atoms, double x, double y,
                                double z, double inflate) {
  if (bins.empty)
    return false;
  const std::size_t bix = std::min<std::size_t>(
      bins.bx - 1, static_cast<std::size_t>(
                       std::max(0.0, std::floor((x - bins.min_x) / bins.cell))));
  const std::size_t biy = std::min<std::size_t>(
      bins.by - 1, static_cast<std::size_t>(
                       std::max(0.0, std::floor((y - bins.min_y) / bins.cell))));
  const std::size_t biz = std::min<std::size_t>(
      bins.bz - 1, static_cast<std::size_t>(
                       std::max(0.0, std::floor((z - bins.min_z) / bins.cell))));
  const std::size_t x0 = (bix == 0) ? 0 : bix - 1;
  const std::size_t y0 = (biy == 0) ? 0 : biy - 1;
  const std::size_t z0 = (biz == 0) ? 0 : biz - 1;
  const std::size_t x1 = std::min(bins.bx - 1, bix + 1);
  const std::size_t y1 = std::min(bins.by - 1, biy + 1);
  const std::size_t z1 = std::min(bins.bz - 1, biz + 1);
  for (std::size_t kk = z0; kk <= z1; ++kk) {
    for (std::size_t jj = y0; jj <= y1; ++jj) {
      for (std::size_t ii = x0; ii <= x1; ++ii) {
        const auto &cell_atoms = bins.bins[(kk * bins.by + jj) * bins.bx + ii];
        for (std::size_t atom_idx : cell_atoms) {
          const auto &atom = atoms[atom_idx];
          const double radius = atom.radius + inflate;
          const double dx = x - atom.x;
          const double dy = y - atom.y;
          const double dz = z - atom.z;
          if (dx * dx + dy * dy + dz * dz <= radius * radius)
            return true;
        }
      }
    }
  }
  return false;
}

template <typename AtomRange>
bool point_inside_vdw_half_grid(const AtomSpatialBins &bins,
                                const AtomRange &atoms, double x, double y,
                                double z) {
  if (bins.empty)
    return false;
  const std::size_t bix = std::min<std::size_t>(
      bins.bx - 1, static_cast<std::size_t>(
                       std::max(0.0, std::floor((x - bins.min_x) / bins.cell))));
  const std::size_t biy = std::min<std::size_t>(
      bins.by - 1, static_cast<std::size_t>(
                       std::max(0.0, std::floor((y - bins.min_y) / bins.cell))));
  const std::size_t biz = std::min<std::size_t>(
      bins.bz - 1, static_cast<std::size_t>(
                       std::max(0.0, std::floor((z - bins.min_z) / bins.cell))));
  const std::size_t x0 = (bix == 0) ? 0 : bix - 1;
  const std::size_t y0 = (biy == 0) ? 0 : biy - 1;
  const std::size_t z0 = (biz == 0) ? 0 : biz - 1;
  const std::size_t x1 = std::min(bins.bx - 1, bix + 1);
  const std::size_t y1 = std::min(bins.by - 1, biy + 1);
  const std::size_t z1 = std::min(bins.bz - 1, biz + 1);
  for (std::size_t kk = z0; kk <= z1; ++kk) {
    for (std::size_t jj = y0; jj <= y1; ++jj) {
      for (std::size_t ii = x0; ii <= x1; ++ii) {
        const auto &cell_atoms = bins.bins[(kk * bins.by + jj) * bins.bx + ii];
        for (std::size_t atom_idx : cell_atoms) {
          const auto &atom = atoms[atom_idx];
          const double radius = atom.radius;
          const double dx = x - atom.x;
          const double dy = y - atom.y;
          const double dz = z - atom.z;
          if (dx * dx + dy * dy + dz * dz < radius * radius)
            return true;
        }
      }
    }
  }
  return false;
}

double dielectric_from_inside_fraction(double inside_fraction, double eps_in,
                                       double eps_out) noexcept {
  const double f = std::clamp(inside_fraction, 0.0, 1.0);
  const double denom = f / std::max(1.0e-12, eps_in) +
                       (1.0 - f) / std::max(1.0e-12, eps_out);
  return denom > 0.0 ? 1.0 / denom : 0.5 * (eps_in + eps_out);
}

template <typename AtomRange>
double segment_inside_fraction(const AtomSpatialBins &bins,
                               const AtomRange &atoms, double x0, double y0,
                               double z0, double x1, double y1, double z1,
                               double inflate) {
  if (bins.empty)
    return 0.0;

  const double sx = x1 - x0;
  const double sy = y1 - y0;
  const double sz = z1 - z0;
  const double a = sx * sx + sy * sy + sz * sz;
  if (a <= 0.0)
    return 0.0;

  const double pad = bins.max_radius + 1.0e-12;
  auto clamp_cell = [](double value, std::size_t n) {
    return static_cast<std::size_t>(std::clamp<double>(
        std::floor(value), 0.0, static_cast<double>(n - 1)));
  };
  const std::size_t ix0 = clamp_cell(
      (std::min(x0, x1) - pad - bins.min_x) / bins.cell, bins.bx);
  const std::size_t ix1 = clamp_cell(
      (std::max(x0, x1) + pad - bins.min_x) / bins.cell, bins.bx);
  const std::size_t iy0 = clamp_cell(
      (std::min(y0, y1) - pad - bins.min_y) / bins.cell, bins.by);
  const std::size_t iy1 = clamp_cell(
      (std::max(y0, y1) + pad - bins.min_y) / bins.cell, bins.by);
  const std::size_t iz0 = clamp_cell(
      (std::min(z0, z1) - pad - bins.min_z) / bins.cell, bins.bz);
  const std::size_t iz1 = clamp_cell(
      (std::max(z0, z1) + pad - bins.min_z) / bins.cell, bins.bz);

  std::vector<std::pair<double, double>> intervals;
  for (std::size_t iz = iz0; iz <= iz1; ++iz) {
    for (std::size_t iy = iy0; iy <= iy1; ++iy) {
      for (std::size_t ix = ix0; ix <= ix1; ++ix) {
        const auto &cell_atoms = bins.bins[(iz * bins.by + iy) * bins.bx + ix];
        for (std::size_t atom_idx : cell_atoms) {
          const auto &atom = atoms[atom_idx];
          const double radius = atom.radius + inflate;
          const double px = x0 - atom.x;
          const double py = y0 - atom.y;
          const double pz = z0 - atom.z;
          const double b = 2.0 * (px * sx + py * sy + pz * sz);
          const double c = px * px + py * py + pz * pz - radius * radius;
          const double disc = b * b - 4.0 * a * c;
          if (disc < 0.0)
            continue;
          const double root = std::sqrt(std::max(0.0, disc));
          double t0 = (-b - root) / (2.0 * a);
          double t1 = (-b + root) / (2.0 * a);
          if (t1 < 0.0 || t0 > 1.0)
            continue;
          t0 = std::clamp(t0, 0.0, 1.0);
          t1 = std::clamp(t1, 0.0, 1.0);
          if (t1 > t0)
            intervals.emplace_back(t0, t1);
        }
      }
    }
  }
  if (intervals.empty())
    return 0.0;

  std::sort(intervals.begin(), intervals.end());
  double covered = 0.0;
  double begin = intervals.front().first;
  double end = intervals.front().second;
  for (std::size_t i = 1; i < intervals.size(); ++i) {
    if (intervals[i].first <= end) {
      end = std::max(end, intervals[i].second);
    } else {
      covered += end - begin;
      begin = intervals[i].first;
      end = intervals[i].second;
    }
  }
  covered += end - begin;
  return std::clamp(covered, 0.0, 1.0);
}

double fallback_radius_from_mass(double mass) noexcept {
  if (!std::isfinite(mass) || mass <= 0.0)
    return 0.0;
  if (mass < 2.5)
    return 1.20; // H
  if (mass < 13.5)
    return 1.70; // C
  if (mass < 15.5)
    return 1.55; // N
  if (mass < 18.5)
    return 1.52; // O/F
  if (mass < 32.5)
    return 1.80; // P/S
  return 1.80;
}

} // namespace

void Analyzer::configure(const arguments &args) {
  args_ = args.pbe;
  if (args_.pbsa) {
    args_.flag_pbe = true;
    args_.dgsolv = true;
  }
}

void Analyzer::begin_trajectory(std::size_t trajectory_id,
                                const std::string &source_label) {
  if (!enabled())
    return;
  current_trajectory_id_ = trajectory_id;
  current_source_label_ = source_label;
}

void Analyzer::end_trajectory(std::fstream &log) {
  if (!enabled())
    return;
  if (log.is_open()) {
    TitanUI::log_kv(log, "PBE Trajectory Complete",
                    std::to_string(current_trajectory_id_) + " " +
                        current_source_label_);
  }
}

void Analyzer::collect_atoms(const MOL &mol, const Multi_MOL_parameters &params,
                             std::vector<Atom> &atoms,
                             const char *label) const {
  const std::size_t n = mol.size();
  if (params.size() < n) {
    throw std::runtime_error(std::string("PBE: ") + label +
                             " parameter table is shorter than coordinates.");
  }
  const auto xs = mol.get_pos_x();
  const auto ys = mol.get_pos_y();
  const auto zs = mol.get_pos_z();
  for (std::size_t i = 0; i < n; ++i) {
    const double q = params.charge[i];
    double r = params.radius[i];
    if ((!std::isfinite(r) || r <= 0.0) && i < params.mass.size())
      r = fallback_radius_from_mass(params.mass[i]);
    if (!std::isfinite(q) || !std::isfinite(r) || r <= 0.0) {
      std::ostringstream msg;
      msg << "PBE: missing/non-finite charge or VDW radius in " << label
          << " atom " << (i + 1)
          << ". Check lib_ff assignment before running -pbe.";
      throw std::runtime_error(msg.str());
    }
    if (!std::isfinite(xs[i]) || !std::isfinite(ys[i]) ||
        !std::isfinite(zs[i])) {
      std::ostringstream msg;
      msg << "PBE: non-finite coordinate in " << label << " atom " << (i + 1)
          << ".";
      throw std::runtime_error(msg.str());
    }
    atoms.push_back({xs[i], ys[i], zs[i], q, r, {}, {}, {}, 0, {}});
  }
}

void Analyzer::collect_test_pqr_atoms(std::vector<Atom> &atoms) const {
  std::ifstream in(args_.test_pqr);
  if (!in.is_open()) {
    throw std::runtime_error("PBE: failed to open test PQR fixture: " +
                             args_.test_pqr);
  }

  std::string line;
  std::size_t line_no = 0;
  while (std::getline(in, line)) {
    ++line_no;
    const std::size_t first = line.find_first_not_of(" \t\r\n");
    if (first == std::string::npos)
      continue;
    const bool atom_record = line.compare(first, 4, "ATOM") == 0 ||
                             line.compare(first, 6, "HETATM") == 0;
    if (!atom_record)
      continue;

    std::istringstream iss(line);
    std::vector<double> numeric_tokens;
    std::string token;
    while (iss >> token) {
      double value = 0.0;
      if (parse_double_token(token, value))
        numeric_tokens.push_back(value);
    }
    if (numeric_tokens.size() < 5) {
      std::ostringstream msg;
      msg << "PBE: malformed PQR fixture line " << line_no << " in "
          << args_.test_pqr
          << ". Expected x y z charge radius numeric fields.";
      throw std::runtime_error(msg.str());
    }

    const std::size_t n = numeric_tokens.size();
    const double x = numeric_tokens[n - 5];
    const double y = numeric_tokens[n - 4];
    const double z = numeric_tokens[n - 3];
    const double q = numeric_tokens[n - 2];
    const double r = numeric_tokens[n - 1];
    if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z) ||
        !std::isfinite(q) || !std::isfinite(r) || r <= 0.0) {
      std::ostringstream msg;
      msg << "PBE: invalid coordinate, charge, or radius on PQR fixture line "
          << line_no << " in " << args_.test_pqr << ".";
      throw std::runtime_error(msg.str());
    }
    std::string element;
    std::string atom_name;
    std::string resname;
    std::string chain;
    int resnr = 0;
    {
      std::istringstream token_stream(line);
      std::vector<std::string> tokens;
      while (token_stream >> token)
        tokens.push_back(token);
      if (tokens.size() >= 3)
        atom_name = tokens[2];
      if (tokens.size() >= 4)
        resname = tokens[3];
      if (tokens.size() >= 5) {
        double maybe_resnr = 0.0;
        if (parse_double_token(tokens[4], maybe_resnr)) {
          resnr = static_cast<int>(std::lround(maybe_resnr));
        } else {
          chain = tokens[4];
          if (tokens.size() >= 6 && parse_double_token(tokens[5], maybe_resnr))
            resnr = static_cast<int>(std::lround(maybe_resnr));
        }
      }
      if (!tokens.empty()) {
        double maybe_element_number = 0.0;
        if (!parse_double_token(tokens.back(), maybe_element_number))
          element = tokens.back();
      }
    }
    atoms.push_back({x, y, z, q, r, element, atom_name, resname, resnr, chain});
  }

  if (atoms.empty()) {
    throw std::runtime_error("PBE: no ATOM/HETATM records found in test PQR "
                             "fixture: " +
                             args_.test_pqr);
  }
}

std::string Analyzer::inferred_element_symbol(const Atom &atom) const {
  const std::string explicit_element = uppercase_letters(atom.element);
  if (!explicit_element.empty()) {
    if (explicit_element.size() == 1 || explicit_element == "CL" ||
        explicit_element == "BR" || explicit_element == "NA" ||
        explicit_element == "MG" || explicit_element == "ZN" ||
        explicit_element == "FE" || explicit_element == "CA")
      return explicit_element;
    return std::string(1, explicit_element.front());
  }

  const std::string name = uppercase_letters(atom.atom_name);
  if (name.empty())
    return "UNK";
  const std::string res = uppercase_letters(atom.resname);
  if (name == "CL" || name == "BR")
    return name;
  if (name == "NA" && (res == "NA" || res == "SOD"))
    return "NA";
  if (name == "MG" && res == "MG")
    return "MG";
  if (name == "ZN" && res == "ZN")
    return "ZN";
  if (name == "FE" && res == "FE")
    return "FE";
  if (name == "CA" && (res == "CA" || res == "CAL" || res == "CA2"))
    return "CA";
  return std::string(1, name.front());
}

double Analyzer::effective_surface_radius(const Atom &atom) const {
  if (effective_radius_convention() != "delphi-bondi")
    return atom.radius;
  if (args_.surface == "delphi-vdw" && !args_.test_pqr.empty() &&
      atom.radius > 0.0) {
    return atom.radius;
  }
  const std::string element = inferred_element_symbol(atom);
  if (element == "H")
    return std::max(atom.radius, 1.30);
  if (element == "C")
    return std::max(atom.radius, 1.80);
  if (element == "N")
    return std::max(atom.radius, 1.60);
  if (element == "O")
    return std::max(atom.radius, 1.50);
  if (element == "S" || element == "P")
    return std::max(atom.radius, 1.80);
  if (element == "F")
    return std::max(atom.radius, 1.47);
  if (element == "CL")
    return std::max(atom.radius, 1.75);
  if (element == "BR")
    return std::max(atom.radius, 1.85);
  if (element == "I")
    return std::max(atom.radius, 1.98);
  return atom.radius;
}

std::string Analyzer::radius_media_source(const Atom &atom) const {
  if (effective_radius_convention() != "delphi-bondi")
    return "native";
  if (args_.surface == "delphi-vdw" && !args_.test_pqr.empty() &&
      atom.radius > 0.0)
    return "pqr_siz_radius";
  if (effective_surface_radius(atom) > atom.radius + 1.0e-12)
    return "bondi_floor";
  return "input_radius";
}

std::string Analyzer::effective_radius_convention() const {
  if (args_.radius_convention != "auto")
    return args_.radius_convention;
  if (args_.surface.rfind("delphi-", 0) == 0 || uses_gaussian_model() ||
      args_.biomodel == "rpbe" || args_.charge_assignment == "delphi" ||
      args_.charge_assignment == "gaussian")
    return "delphi-bondi";
  return "native";
}

bool Analyzer::uses_gaussian_model() const {
  return args_.dielectric_model == "gaussian" ||
         surface_uses_gaussian(args_.surface) ||
         args_.charge_assignment == "gaussian" || args_.biomodel == "rpbe";
}

double Analyzer::gaussian_density_at(const std::vector<Atom> &atoms, double x,
                                     double y, double z) const {
  const double exponent = std::max(1.0e-12, args_.gaussian_exp);
  const double sigma = std::max(1.0e-12, args_.gaussian_sigma);
  double density = 0.0;
  for (const Atom &atom : atoms) {
    const double r0 = std::max(1.0e-12, effective_surface_radius(atom));
    const double cutoff = 3.0 * sigma * r0;
    const double dx = x - atom.x;
    const double dy = y - atom.y;
    const double dz = z - atom.z;
    const double r2 = dx * dx + dy * dy + dz * dz;
    if (r2 > cutoff * cutoff)
      continue;
    const double scaled = std::sqrt(r2) / (sigma * r0);
    density += std::exp(-std::pow(scaled, exponent));
    if (density >= 1.0)
      return 1.0;
  }
  return std::clamp(density, 0.0, 1.0);
}

double Analyzer::gaussian_dielectric(double density, double eps_in,
                                     double eps_out) const {
  const double d = std::clamp(density, 0.0, 1.0);
  if (args_.surface == "gaussian-cutoff") {
    const double cutoff =
        args_.surface_cutoff == "auto" ? args_.density_cutoff
                                       : args_.surface_cutoff_value;
    if (d >= cutoff)
      return eps_in;
    if (d > 0.0 && args_.gap_diel > eps_in && args_.gap_diel < eps_out)
      return args_.gap_diel;
  }
  return eps_out - (eps_out - eps_in) * d;
}

std::vector<Analyzer::Atom>
Analyzer::adjusted_surface_atoms(const std::vector<Atom> &atoms) const {
  std::vector<Atom> adjusted = atoms;
  if (effective_radius_convention() != "delphi-bondi")
    return adjusted;
  for (Atom &atom : adjusted)
    atom.radius = effective_surface_radius(atom);
  return adjusted;
}

Analyzer::GridShape Analyzer::build_shape(const std::vector<Atom> &atoms) const {
  return build_shape(atoms, args_.scale, args_.margin);
}

Analyzer::GridShape Analyzer::build_shape(const std::vector<Atom> &atoms,
                                          double scale, double margin) const {
  if (atoms.empty())
    throw std::runtime_error("PBE: no atoms were available for the frame.");

  double min_x = std::numeric_limits<double>::max();
  double min_y = std::numeric_limits<double>::max();
  double min_z = std::numeric_limits<double>::max();
  double max_x = std::numeric_limits<double>::lowest();
  double max_y = std::numeric_limits<double>::lowest();
  double max_z = std::numeric_limits<double>::lowest();
  double coord_min_x = std::numeric_limits<double>::max();
  double coord_min_y = std::numeric_limits<double>::max();
  double coord_min_z = std::numeric_limits<double>::max();
  double coord_max_x = std::numeric_limits<double>::lowest();
  double coord_max_y = std::numeric_limits<double>::lowest();
  double coord_max_z = std::numeric_limits<double>::lowest();
  const double surface_inflate =
      surface_uses_probe_exclusion(args_.surface) ? args_.probe_radius : 0.0;
  for (const Atom &atom : atoms) {
    coord_min_x = std::min(coord_min_x, atom.x);
    coord_min_y = std::min(coord_min_y, atom.y);
    coord_min_z = std::min(coord_min_z, atom.z);
    coord_max_x = std::max(coord_max_x, atom.x);
    coord_max_y = std::max(coord_max_y, atom.y);
    coord_max_z = std::max(coord_max_z, atom.z);
    const double r = effective_surface_radius(atom) + surface_inflate;
    min_x = std::min(min_x, atom.x - r);
    min_y = std::min(min_y, atom.y - r);
    min_z = std::min(min_z, atom.z - r);
    max_x = std::max(max_x, atom.x + r);
    max_y = std::max(max_y, atom.y + r);
    max_z = std::max(max_z, atom.z + r);
  }

  min_x -= margin;
  min_y -= margin;
  min_z -= margin;
  max_x += margin;
  max_y += margin;
  max_z += margin;

  if (args_.gridbox == "cubic") {
    const double cx = 0.5 * (min_x + max_x);
    const double cy = 0.5 * (min_y + max_y);
    const double cz = 0.5 * (min_z + max_z);
    const double span =
        std::max({max_x - min_x, max_y - min_y, max_z - min_z});
    min_x = cx - 0.5 * span;
    max_x = cx + 0.5 * span;
    min_y = cy - 0.5 * span;
    max_y = cy + 0.5 * span;
    min_z = cz - 0.5 * span;
    max_z = cz + 0.5 * span;
  }

  GridShape shape;
  shape.h = 1.0 / scale;
  if (!args_.test_pqr.empty() && args_.test_grid_size > 0) {
    const std::size_t n =
        static_cast<std::size_t>(args_.test_grid_size);
    double cx = 0.5 * (coord_min_x + coord_max_x);
    double cy = 0.5 * (coord_min_y + coord_max_y);
    double cz = 0.5 * (coord_min_z + coord_max_z);
    if (args_.test_center_user) {
      cx = args_.test_center_x;
      cy = args_.test_center_y;
      cz = args_.test_center_z;
    } else if (args_.focus_center_user) {
      cx = args_.focus_center_x;
      cy = args_.focus_center_y;
      cz = args_.focus_center_z;
    }
    shape.nx = n;
    shape.ny = n;
    shape.nz = n;
    shape.ox = cx - 0.5 * static_cast<double>(n - 1) * shape.h;
    shape.oy = cy - 0.5 * static_cast<double>(n - 1) * shape.h;
    shape.oz = cz - 0.5 * static_cast<double>(n - 1) * shape.h;
    if (shape.nx > std::numeric_limits<std::size_t>::max() / shape.ny ||
        shape.nx * shape.ny >
            std::numeric_limits<std::size_t>::max() / shape.nz) {
      throw std::runtime_error("PBE: grid dimensions overflow size_t.");
    }
    shape.n = shape.nx * shape.ny * shape.nz;
    shape.memory_bytes =
        shape.n * (7 * sizeof(double) + 2 * sizeof(unsigned char));
    if (shape.n > args_.max_grid_points) {
      std::ostringstream msg;
      msg << "PBE: requested grid " << shape.nx << "x" << shape.ny << "x"
          << shape.nz << " = " << shape.n
          << " points exceeds -pbe-max-grid " << args_.max_grid_points << ".";
      throw std::runtime_error(msg.str());
    }
    return shape;
  }
  shape.ox = min_x;
  shape.oy = min_y;
  shape.oz = min_z;
  auto count_axis = [&](double span) {
    std::size_t n =
        static_cast<std::size_t>(std::ceil(std::max(0.0, span) / shape.h)) + 1;
    n = std::max<std::size_t>(3, n);
    if ((n & 1ULL) == 0)
      ++n;
    return n;
  };
  shape.nx = count_axis(max_x - min_x);
  shape.ny = count_axis(max_y - min_y);
  shape.nz = count_axis(max_z - min_z);
  if (args_.focus_center_user && args_.focus == "delphi") {
    shape.ox =
        args_.focus_center_x - 0.5 * static_cast<double>(shape.nx - 1) * shape.h;
    shape.oy =
        args_.focus_center_y - 0.5 * static_cast<double>(shape.ny - 1) * shape.h;
    shape.oz =
        args_.focus_center_z - 0.5 * static_cast<double>(shape.nz - 1) * shape.h;
  }

  if (shape.nx > std::numeric_limits<std::size_t>::max() / shape.ny ||
      shape.nx * shape.ny >
          std::numeric_limits<std::size_t>::max() / shape.nz) {
    throw std::runtime_error("PBE: grid dimensions overflow size_t.");
  }
  shape.n = shape.nx * shape.ny * shape.nz;
  shape.memory_bytes =
      shape.n * (7 * sizeof(double) + 2 * sizeof(unsigned char));
  if (shape.n > args_.max_grid_points) {
    std::ostringstream msg;
    msg << "PBE: requested grid " << shape.nx << "x" << shape.ny << "x"
        << shape.nz << " = " << shape.n << " points exceeds -pbe-max-grid "
        << args_.max_grid_points << ".";
    throw std::runtime_error(msg.str());
  }
  return shape;
}

void Analyzer::ensure_workspace(const GridShape &shape) {
  phi_.assign(shape.n, 0.0);
  rho_.assign(shape.n, 0.0);
  eps_x_.assign(shape.n, args_.ext_diel);
  eps_y_.assign(shape.n, args_.ext_diel);
  eps_z_.assign(shape.n, args_.ext_diel);
  denom_.assign(shape.n, 1.0);
  screen_.assign(shape.n, 0.0);
  inside_.assign(shape.n, 0);
  ion_accessible_.assign(shape.n, 1);
  last_gaussian_density_sum_ = 0.0;
  last_gaussian_cutoff_points_ = 0;
  last_gaussian_eps_min_ = 0.0;
  last_gaussian_eps_max_ = 0.0;
  last_gaussian_eps_mean_ = 0.0;
  last_vdwms_audit_ = VdwmsAudit{};
}

void Analyzer::classify_surface(const GridShape &shape,
                                const std::vector<Atom> &atoms,
                                TitanExecutor &executor) {
  const std::vector<Atom> surface_atoms = adjusted_surface_atoms(atoms);
  if (uses_gaussian_model()) {
    double max_base_radius = shape.h;
    for (const Atom &atom : surface_atoms)
      max_base_radius = std::max(max_base_radius, atom.radius);
    const double sigma = std::max(1.0e-12, args_.gaussian_sigma);
    const double influence = 3.0 * sigma * max_base_radius;
    const AtomSpatialBins gaussian_bins =
        build_atom_bins(surface_atoms, std::max(0.0, influence), shape.h);
    auto density_at = [&](double x, double y, double z) {
      const double exponent = std::max(1.0e-12, args_.gaussian_exp);
      const int radius_bins = static_cast<int>(
          std::ceil(std::max(shape.h, influence) / gaussian_bins.cell)) + 1;
      const int bx0 = static_cast<int>(std::floor((x - gaussian_bins.min_x) /
                                                  gaussian_bins.cell));
      const int by0 = static_cast<int>(std::floor((y - gaussian_bins.min_y) /
                                                  gaussian_bins.cell));
      const int bz0 = static_cast<int>(std::floor((z - gaussian_bins.min_z) /
                                                  gaussian_bins.cell));
      double density = 0.0;
      for (int dz = -radius_bins; dz <= radius_bins; ++dz) {
        const int bz = bz0 + dz;
        if (bz < 0 || bz >= static_cast<int>(gaussian_bins.bz))
          continue;
        for (int dy = -radius_bins; dy <= radius_bins; ++dy) {
          const int by = by0 + dy;
          if (by < 0 || by >= static_cast<int>(gaussian_bins.by))
            continue;
          for (int dx = -radius_bins; dx <= radius_bins; ++dx) {
            const int bx = bx0 + dx;
            if (bx < 0 || bx >= static_cast<int>(gaussian_bins.bx))
              continue;
            const auto &bucket = gaussian_bins.bins
                [(static_cast<std::size_t>(bz) * gaussian_bins.by +
                  static_cast<std::size_t>(by)) *
                     gaussian_bins.bx +
                 static_cast<std::size_t>(bx)];
            for (std::size_t atom_idx : bucket) {
              const Atom &atom = surface_atoms[atom_idx];
              const double r0 = std::max(1.0e-12, atom.radius);
              const double cutoff = 3.0 * sigma * r0;
              const double ddx = x - atom.x;
              const double ddy = y - atom.y;
              const double ddz = z - atom.z;
              const double r2 = ddx * ddx + ddy * ddy + ddz * ddz;
              if (r2 > cutoff * cutoff)
                continue;
              const double scaled = std::sqrt(r2) / (sigma * r0);
              density += std::exp(-std::pow(scaled, exponent));
              if (density >= 1.0)
                return 1.0;
            }
          }
        }
      }
      return std::clamp(density, 0.0, 1.0);
    };
    const double cutoff =
        args_.surface_cutoff == "auto" ? args_.density_cutoff
                                       : args_.surface_cutoff_value;
    const std::size_t workers = effective_workers(executor, shape.nz);
    std::vector<double> local_density_sum(workers, 0.0);
    std::vector<std::size_t> local_cutoff_points(workers, 0);
    executor.parallel_for(shape.nz, workers, 1,
                          [&](std::size_t worker, std::size_t z_begin,
                              std::size_t z_end) {
      double density_sum = 0.0;
      std::size_t cutoff_points = 0;
      for (std::size_t k = z_begin; k < z_end; ++k) {
        const double z = shape.oz + static_cast<double>(k) * shape.h;
        for (std::size_t j = 0; j < shape.ny; ++j) {
          const double y = shape.oy + static_cast<double>(j) * shape.h;
          for (std::size_t i = 0; i < shape.nx; ++i) {
            const double x = shape.ox + static_cast<double>(i) * shape.h;
            const double density = density_at(x, y, z);
            density_sum += density;
            const bool inside = density >= cutoff;
            if (inside)
              ++cutoff_points;
            inside_[idx3(i, j, k, shape.nx, shape.ny)] = inside ? 1U : 0U;
          }
        }
      }
      local_density_sum[worker] += density_sum;
      local_cutoff_points[worker] += cutoff_points;
    });
    last_gaussian_density_sum_ =
        std::accumulate(local_density_sum.begin(), local_density_sum.end(),
                        0.0);
    last_gaussian_cutoff_points_ = std::accumulate(
        local_cutoff_points.begin(), local_cutoff_points.end(),
        static_cast<std::size_t>(0));
    return;
  }
  double max_radius = 0.0;
  double min_x = std::numeric_limits<double>::max();
  double min_y = std::numeric_limits<double>::max();
  double min_z = std::numeric_limits<double>::max();
  double max_x = std::numeric_limits<double>::lowest();
  double max_y = std::numeric_limits<double>::lowest();
  double max_z = std::numeric_limits<double>::lowest();
  for (const Atom &atom : surface_atoms) {
    const double surface_radius =
        atom.radius +
        (surface_uses_probe_exclusion(args_.surface) ? args_.probe_radius
                                                     : 0.0);
    max_radius = std::max(max_radius, surface_radius);
    min_x = std::min(min_x, atom.x);
    min_y = std::min(min_y, atom.y);
    min_z = std::min(min_z, atom.z);
    max_x = std::max(max_x, atom.x);
    max_y = std::max(max_y, atom.y);
    max_z = std::max(max_z, atom.z);
  }
  const double cell = std::max(max_radius, shape.h);
  const std::size_t bx = std::max<std::size_t>(
      1, static_cast<std::size_t>(std::floor((max_x - min_x) / cell)) + 1);
  const std::size_t by = std::max<std::size_t>(
      1, static_cast<std::size_t>(std::floor((max_y - min_y) / cell)) + 1);
  const std::size_t bz = std::max<std::size_t>(
      1, static_cast<std::size_t>(std::floor((max_z - min_z) / cell)) + 1);

  std::vector<std::vector<std::size_t>> bins(bx * by * bz);
  auto bin_index = [&](double x, double y, double z) {
    const std::size_t ix = std::min<std::size_t>(
        bx - 1, static_cast<std::size_t>(std::max(0.0, std::floor((x - min_x) / cell))));
    const std::size_t iy = std::min<std::size_t>(
        by - 1, static_cast<std::size_t>(std::max(0.0, std::floor((y - min_y) / cell))));
    const std::size_t iz = std::min<std::size_t>(
        bz - 1, static_cast<std::size_t>(std::max(0.0, std::floor((z - min_z) / cell))));
    return (iz * by + iy) * bx + ix;
  };
  for (std::size_t atom_idx = 0; atom_idx < surface_atoms.size();
       ++atom_idx) {
    const Atom &atom = surface_atoms[atom_idx];
    bins[bin_index(atom.x, atom.y, atom.z)].push_back(atom_idx);
  }

  const std::size_t workers = effective_workers(executor, shape.nz);
  executor.parallel_for(shape.nz, workers, 1,
                        [&](std::size_t, std::size_t z_begin,
                            std::size_t z_end) {
    for (std::size_t k = z_begin; k < z_end; ++k) {
      const double z = shape.oz + static_cast<double>(k) * shape.h;
      for (std::size_t j = 0; j < shape.ny; ++j) {
        const double y = shape.oy + static_cast<double>(j) * shape.h;
        for (std::size_t i = 0; i < shape.nx; ++i) {
          const double x = shape.ox + static_cast<double>(i) * shape.h;
          const std::size_t bix = std::min<std::size_t>(
              bx - 1, static_cast<std::size_t>(
                          std::max(0.0, std::floor((x - min_x) / cell))));
          const std::size_t biy = std::min<std::size_t>(
              by - 1, static_cast<std::size_t>(
                          std::max(0.0, std::floor((y - min_y) / cell))));
          const std::size_t biz = std::min<std::size_t>(
              bz - 1, static_cast<std::size_t>(
                          std::max(0.0, std::floor((z - min_z) / cell))));

          bool hit = false;
          const std::size_t x0 = (bix == 0) ? 0 : bix - 1;
          const std::size_t y0 = (biy == 0) ? 0 : biy - 1;
          const std::size_t z0 = (biz == 0) ? 0 : biz - 1;
          const std::size_t x1 = std::min(bx - 1, bix + 1);
          const std::size_t y1 = std::min(by - 1, biy + 1);
          const std::size_t z1 = std::min(bz - 1, biz + 1);
          for (std::size_t kk = z0; kk <= z1 && !hit; ++kk) {
            for (std::size_t jj = y0; jj <= y1 && !hit; ++jj) {
              for (std::size_t ii = x0; ii <= x1 && !hit; ++ii) {
                const auto &cell_atoms = bins[(kk * by + jj) * bx + ii];
                for (std::size_t atom_idx : cell_atoms) {
                  const Atom &atom = surface_atoms[atom_idx];
                  const double dx = x - atom.x;
                  const double dy = y - atom.y;
                  const double dz = z - atom.z;
                  const double surface_radius =
                      atom.radius +
                      (surface_uses_probe_exclusion(args_.surface)
                           ? args_.probe_radius
                           : 0.0);
                  if (dx * dx + dy * dy + dz * dz <=
                      surface_radius * surface_radius) {
                    hit = true;
                    break;
                  }
                }
              }
            }
          }
          inside_[idx3(i, j, k, shape.nx, shape.ny)] = hit ? 1U : 0U;
        }
      }
    }
  });

  if (!surface_uses_ses_flood_fill(args_.surface))
    return;

  if (shape.n > static_cast<std::size_t>(
                    std::numeric_limits<std::uint32_t>::max())) {
    throw std::runtime_error(
        "PBE: SES surface flood fill currently requires <= 2^32 grid points.");
  }

  std::vector<unsigned char> solvent_accessible(shape.n, 0);
  std::vector<std::uint32_t> queue;
  queue.reserve(std::min<std::size_t>(shape.n, 1048576));
  auto enqueue = [&](std::size_t id) {
    if (inside_[id] || solvent_accessible[id])
      return;
    solvent_accessible[id] = 1U;
    queue.push_back(static_cast<std::uint32_t>(id));
  };

  for (std::size_t k = 0; k < shape.nz; ++k) {
    for (std::size_t j = 0; j < shape.ny; ++j) {
      enqueue(idx3(0, j, k, shape.nx, shape.ny));
      enqueue(idx3(shape.nx - 1, j, k, shape.nx, shape.ny));
    }
    for (std::size_t i = 0; i < shape.nx; ++i) {
      enqueue(idx3(i, 0, k, shape.nx, shape.ny));
      enqueue(idx3(i, shape.ny - 1, k, shape.nx, shape.ny));
    }
  }
  for (std::size_t j = 0; j < shape.ny; ++j) {
    for (std::size_t i = 0; i < shape.nx; ++i) {
      enqueue(idx3(i, j, 0, shape.nx, shape.ny));
      enqueue(idx3(i, j, shape.nz - 1, shape.nx, shape.ny));
    }
  }

  for (std::size_t head = 0; head < queue.size(); ++head) {
    const std::size_t id = static_cast<std::size_t>(queue[head]);
    const std::size_t i = id % shape.nx;
    const std::size_t j = (id / shape.nx) % shape.ny;
    const std::size_t k = id / (shape.nx * shape.ny);
    if (i > 0)
      enqueue(id - 1);
    if (i + 1 < shape.nx)
      enqueue(id + 1);
    if (j > 0)
      enqueue(id - shape.nx);
    if (j + 1 < shape.ny)
      enqueue(id + shape.nx);
    if (k > 0)
      enqueue(id - shape.nx * shape.ny);
    if (k + 1 < shape.nz)
      enqueue(id + shape.nx * shape.ny);
  }

  for (std::size_t id = 0; id < shape.n; ++id)
    inside_[id] = solvent_accessible[id] ? 0U : 1U;
}

void Analyzer::classify_ion_accessibility(const GridShape &shape,
                                          const std::vector<Atom> &atoms,
                                          TitanExecutor &executor) {
  ion_accessible_.assign(shape.n, 0);
  if (args_.ion_access == "dielectric") {
    for (std::size_t id = 0; id < shape.n; ++id)
      ion_accessible_[id] = inside_[id] ? 0U : 1U;
    return;
  }

  if (shape.n > static_cast<std::size_t>(
                    std::numeric_limits<std::uint32_t>::max())) {
    throw std::runtime_error(
        "PBE: Stern ion flood fill currently requires <= 2^32 grid points.");
  }

  const AtomSpatialBins bins =
      build_atom_bins(atoms, args_.ion_radius, shape.h);
  std::vector<unsigned char> blocked(shape.n, 0);
  const std::size_t workers = effective_workers(executor, shape.nz);
  executor.parallel_for(
      shape.nz, workers, 1,
      [&](std::size_t, std::size_t z_begin, std::size_t z_end) {
        for (std::size_t k = z_begin; k < z_end; ++k) {
          const double z = shape.oz + static_cast<double>(k) * shape.h;
          for (std::size_t j = 0; j < shape.ny; ++j) {
            const double y = shape.oy + static_cast<double>(j) * shape.h;
            for (std::size_t i = 0; i < shape.nx; ++i) {
              const double x = shape.ox + static_cast<double>(i) * shape.h;
              const std::size_t id = idx3(i, j, k, shape.nx, shape.ny);
              blocked[id] =
                  point_inside_inflated_atom(bins, atoms, x, y, z,
                                             args_.ion_radius)
                      ? 1U
                      : 0U;
            }
          }
        }
      });

  std::vector<std::uint32_t> queue;
  queue.reserve(std::min<std::size_t>(shape.n, 1048576));
  auto enqueue = [&](std::size_t id) {
    if (blocked[id] || ion_accessible_[id])
      return;
    ion_accessible_[id] = 1U;
    queue.push_back(static_cast<std::uint32_t>(id));
  };

  for (std::size_t k = 0; k < shape.nz; ++k) {
    for (std::size_t j = 0; j < shape.ny; ++j) {
      enqueue(idx3(0, j, k, shape.nx, shape.ny));
      enqueue(idx3(shape.nx - 1, j, k, shape.nx, shape.ny));
    }
    for (std::size_t i = 0; i < shape.nx; ++i) {
      enqueue(idx3(i, 0, k, shape.nx, shape.ny));
      enqueue(idx3(i, shape.ny - 1, k, shape.nx, shape.ny));
    }
  }
  for (std::size_t j = 0; j < shape.ny; ++j) {
    for (std::size_t i = 0; i < shape.nx; ++i) {
      enqueue(idx3(i, j, 0, shape.nx, shape.ny));
      enqueue(idx3(i, j, shape.nz - 1, shape.nx, shape.ny));
    }
  }

  for (std::size_t head = 0; head < queue.size(); ++head) {
    const std::size_t id = static_cast<std::size_t>(queue[head]);
    const std::size_t i = id % shape.nx;
    const std::size_t j = (id / shape.nx) % shape.ny;
    const std::size_t k = id / (shape.nx * shape.ny);
    if (i > 0)
      enqueue(id - 1);
    if (i + 1 < shape.nx)
      enqueue(id + 1);
    if (j > 0)
      enqueue(id - shape.nx);
    if (j + 1 < shape.ny)
      enqueue(id + shape.nx);
    if (k > 0)
      enqueue(id - shape.nx * shape.ny);
    if (k + 1 < shape.nz)
      enqueue(id + shape.nx * shape.ny);
  }
}

Analyzer::ChargeDiagnostics
Analyzer::assign_charges(const GridShape &shape,
                         const std::vector<Atom> &atoms) {
  std::fill(rho_.begin(), rho_.end(), 0.0);
  ChargeDiagnostics diag;
  double original_mx = 0.0;
  double original_my = 0.0;
  double original_mz = 0.0;
  double assigned_mx = 0.0;
  double assigned_my = 0.0;
  double assigned_mz = 0.0;
  auto add_charge = [&](std::size_t i, std::size_t j, std::size_t k,
                        double q) {
    if (q == 0.0)
      return;
    rho_[idx3(i, j, k, shape.nx, shape.ny)] += q;
    diag.assigned_charge += q;
    assigned_mx += q * (shape.ox + static_cast<double>(i) * shape.h);
    assigned_my += q * (shape.oy + static_cast<double>(j) * shape.h);
    assigned_mz += q * (shape.oz + static_cast<double>(k) * shape.h);
  };
  for (const Atom &atom : atoms) {
    original_mx += atom.q * atom.x;
    original_my += atom.q * atom.y;
    original_mz += atom.q * atom.z;
    if (args_.charge_assignment == "gaussian" || args_.biomodel == "rpbe") {
      const double r0 = std::max(1.0e-12, effective_surface_radius(atom));
      const double sigma = std::max(1.0e-12, args_.gaussian_sigma);
      const double influence = std::max(shape.h, 3.0 * sigma * r0);
      const int imin = std::max(
          0, static_cast<int>(std::floor((atom.x - influence - shape.ox) /
                                         shape.h)));
      const int imax = std::min(
          static_cast<int>(shape.nx) - 1,
          static_cast<int>(std::ceil((atom.x + influence - shape.ox) /
                                     shape.h)));
      const int jmin = std::max(
          0, static_cast<int>(std::floor((atom.y - influence - shape.oy) /
                                         shape.h)));
      const int jmax = std::min(
          static_cast<int>(shape.ny) - 1,
          static_cast<int>(std::ceil((atom.y + influence - shape.oy) /
                                     shape.h)));
      const int kmin = std::max(
          0, static_cast<int>(std::floor((atom.z - influence - shape.oz) /
                                         shape.h)));
      const int kmax = std::min(
          static_cast<int>(shape.nz) - 1,
          static_cast<int>(std::ceil((atom.z + influence - shape.oz) /
                                     shape.h)));
      double norm = 0.0;
      for (int kk = kmin; kk <= kmax; ++kk) {
        const double z = shape.oz + static_cast<double>(kk) * shape.h;
        for (int jj = jmin; jj <= jmax; ++jj) {
          const double y = shape.oy + static_cast<double>(jj) * shape.h;
          for (int ii = imin; ii <= imax; ++ii) {
            const double x = shape.ox + static_cast<double>(ii) * shape.h;
            const double dx = x - atom.x;
            const double dy = y - atom.y;
            const double dz = z - atom.z;
            const double scaled =
                std::sqrt(dx * dx + dy * dy + dz * dz) / (sigma * r0);
            norm += std::exp(-std::pow(scaled, args_.gaussian_exp));
          }
        }
      }
      if (norm > 0.0) {
        for (int kk = kmin; kk <= kmax; ++kk) {
          const double z = shape.oz + static_cast<double>(kk) * shape.h;
          for (int jj = jmin; jj <= jmax; ++jj) {
            const double y = shape.oy + static_cast<double>(jj) * shape.h;
            for (int ii = imin; ii <= imax; ++ii) {
              const double x = shape.ox + static_cast<double>(ii) * shape.h;
              const double dx = x - atom.x;
              const double dy = y - atom.y;
              const double dz = z - atom.z;
              const double scaled =
                  std::sqrt(dx * dx + dy * dy + dz * dz) / (sigma * r0);
              const double weight =
                  std::exp(-std::pow(scaled, args_.gaussian_exp)) / norm;
              add_charge(static_cast<std::size_t>(ii),
                         static_cast<std::size_t>(jj),
                         static_cast<std::size_t>(kk), atom.q * weight);
            }
          }
        }
        continue;
      }
    }
    double gx = (atom.x - shape.ox) / shape.h;
    double gy = (atom.y - shape.oy) / shape.h;
    double gz = (atom.z - shape.oz) / shape.h;

    if (args_.charge_assignment == "trilinear") {
      gx = std::clamp(gx, 0.0, static_cast<double>(shape.nx - 1));
      gy = std::clamp(gy, 0.0, static_cast<double>(shape.ny - 1));
      gz = std::clamp(gz, 0.0, static_cast<double>(shape.nz - 1));
    }
    const double fgx = std::floor(gx);
    const double fgy = std::floor(gy);
    const double fgz = std::floor(gz);
    if (args_.charge_assignment == "delphi" &&
        (fgx < 0.0 || fgy < 0.0 || fgz < 0.0 ||
         fgx > static_cast<double>(shape.nx - 2) ||
         fgy > static_cast<double>(shape.ny - 2) ||
         fgz > static_cast<double>(shape.nz - 2))) {
      continue;
    }
    std::size_t i0 = std::min<std::size_t>(
        shape.nx - 2, static_cast<std::size_t>(std::max(0.0, fgx)));
    std::size_t j0 = std::min<std::size_t>(
        shape.ny - 2, static_cast<std::size_t>(std::max(0.0, fgy)));
    std::size_t k0 = std::min<std::size_t>(
        shape.nz - 2, static_cast<std::size_t>(std::max(0.0, fgz)));
    const double fx = clamp_fraction(gx - static_cast<double>(i0));
    const double fy = clamp_fraction(gy - static_cast<double>(j0));
    const double fz = clamp_fraction(gz - static_cast<double>(k0));
    for (std::size_t dz = 0; dz <= 1; ++dz) {
      const double wz = dz ? fz : (1.0 - fz);
      for (std::size_t dy = 0; dy <= 1; ++dy) {
        const double wy = dy ? fy : (1.0 - fy);
        for (std::size_t dx = 0; dx <= 1; ++dx) {
          const double wx = dx ? fx : (1.0 - fx);
          add_charge(i0 + dx, j0 + dy, k0 + dz, atom.q * wx * wy * wz);
        }
      }
    }
  }
  const double dmx = assigned_mx - original_mx;
  const double dmy = assigned_my - original_my;
  const double dmz = assigned_mz - original_mz;
  diag.moment_error = std::sqrt(dmx * dmx + dmy * dmy + dmz * dmz);
  return diag;
}

bool Analyzer::oracle_replay_applies(const GridShape &shape) const {
  if (args_.test_pqr.empty())
    return false;
  if (args_.test_grid_size <= 0)
    return true;
  const std::size_t n = static_cast<std::size_t>(args_.test_grid_size);
  return shape.nx == n && shape.ny == n && shape.nz == n;
}

void Analyzer::apply_oracle_surface_map(const GridShape &shape) {
  if (args_.test_surface_map != "oracle" || !oracle_replay_applies(shape))
    return;
  std::ifstream in("pbe_oracle_surface_map.csv");
  if (!in.is_open()) {
    throw std::runtime_error(
        "PBE: -pbe-test-surface-map=oracle requires "
        "pbe_oracle_surface_map.csv in the run directory.");
  }
  std::string line;
  if (!std::getline(in, line))
    throw std::runtime_error("PBE: empty pbe_oracle_surface_map.csv.");
  const auto header = csv_header_map(split_csv_simple(line));
  const bool has_inside = header.find("inside") != header.end();
  if (!has_inside) {
    throw std::runtime_error(
        "PBE: pbe_oracle_surface_map.csv must contain i,j,k,inside columns.");
  }
  std::fill(inside_.begin(), inside_.end(), 0U);
  std::size_t rows = 0;
  while (std::getline(in, line)) {
    if (trim_copy(line).empty())
      continue;
    const auto fields = split_csv_simple(line);
    std::size_t i = 0, j = 0, k = 0;
    if (!csv_index(fields, header, "i", i) ||
        !csv_index(fields, header, "j", j) ||
        !csv_index(fields, header, "k", k)) {
      throw std::runtime_error(
          "PBE: pbe_oracle_surface_map.csv row is missing i,j,k.");
    }
    if (i >= shape.nx || j >= shape.ny || k >= shape.nz)
      throw std::runtime_error(
          "PBE: pbe_oracle_surface_map.csv index is outside the grid.");
    bool inside = false;
    const auto inside_it = header.find("inside");
    if (inside_it == header.end() || inside_it->second >= fields.size() ||
        !parse_boolish(fields[inside_it->second], inside)) {
      throw std::runtime_error(
          "PBE: pbe_oracle_surface_map.csv inside value is invalid.");
    }
    inside_[idx3(i, j, k, shape.nx, shape.ny)] = inside ? 1U : 0U;
    ++rows;
  }
  if (rows == 0)
    throw std::runtime_error("PBE: pbe_oracle_surface_map.csv has no rows.");
}

Analyzer::ChargeDiagnostics
Analyzer::assign_oracle_charge_grid(const GridShape &shape,
                                    const std::vector<Atom> &atoms) {
  if (args_.test_charge_grid != "oracle" || !oracle_replay_applies(shape))
    return assign_charges(shape, atoms);
  std::ifstream in("pbe_oracle_charge_grid.csv");
  if (!in.is_open()) {
    throw std::runtime_error(
        "PBE: -pbe-test-charge-grid=oracle requires "
        "pbe_oracle_charge_grid.csv in the run directory.");
  }
  std::string line;
  if (!std::getline(in, line))
    throw std::runtime_error("PBE: empty pbe_oracle_charge_grid.csv.");
  const auto header = csv_header_map(split_csv_simple(line));
  std::fill(rho_.begin(), rho_.end(), 0.0);

  ChargeDiagnostics diag;
  double original_mx = 0.0;
  double original_my = 0.0;
  double original_mz = 0.0;
  for (const Atom &atom : atoms) {
    original_mx += atom.q * atom.x;
    original_my += atom.q * atom.y;
    original_mz += atom.q * atom.z;
  }

  double assigned_mx = 0.0;
  double assigned_my = 0.0;
  double assigned_mz = 0.0;
  std::size_t rows = 0;
  while (std::getline(in, line)) {
    if (trim_copy(line).empty())
      continue;
    const auto fields = split_csv_simple(line);
    std::size_t i = 0, j = 0, k = 0;
    double q = 0.0;
    if (!csv_index(fields, header, "i", i) ||
        !csv_index(fields, header, "j", j) ||
        !csv_index(fields, header, "k", k)) {
      throw std::runtime_error(
          "PBE: pbe_oracle_charge_grid.csv row is missing i,j,k.");
    }
    if (!csv_double(fields, header, {"rho", "q", "charge"}, q)) {
      throw std::runtime_error(
          "PBE: pbe_oracle_charge_grid.csv row is missing rho/q/charge.");
    }
    if (i >= shape.nx || j >= shape.ny || k >= shape.nz)
      throw std::runtime_error(
          "PBE: pbe_oracle_charge_grid.csv index is outside the grid.");
    const std::size_t id = idx3(i, j, k, shape.nx, shape.ny);
    rho_[id] += q;
    diag.assigned_charge += q;
    assigned_mx += q * (shape.ox + static_cast<double>(i) * shape.h);
    assigned_my += q * (shape.oy + static_cast<double>(j) * shape.h);
    assigned_mz += q * (shape.oz + static_cast<double>(k) * shape.h);
    ++rows;
  }
  if (rows == 0)
    throw std::runtime_error("PBE: pbe_oracle_charge_grid.csv has no rows.");
  const double dmx = assigned_mx - original_mx;
  const double dmy = assigned_my - original_my;
  const double dmz = assigned_mz - original_mz;
  diag.moment_error = std::sqrt(dmx * dmx + dmy * dmy + dmz * dmz);
  return diag;
}

void Analyzer::apply_oracle_midpoint_epsilon_map(const GridShape &shape) {
  if (args_.test_midpoint_epsilon != "oracle" ||
      !oracle_replay_applies(shape))
    return;
  std::ifstream in("pbe_oracle_midpoint_epsilon.csv");
  if (!in.is_open()) {
    throw std::runtime_error(
        "PBE: -pbe-test-midpoint-epsilon=oracle requires "
        "pbe_oracle_midpoint_epsilon.csv in the run directory.");
  }
  std::string line;
  if (!std::getline(in, line))
    throw std::runtime_error("PBE: empty pbe_oracle_midpoint_epsilon.csv.");
  const auto header = csv_header_map(split_csv_simple(line));
  std::vector<unsigned char> seen(shape.n, 0U);
  std::size_t rows = 0;
  while (std::getline(in, line)) {
    if (trim_copy(line).empty())
      continue;
    const auto fields = split_csv_simple(line);
    std::size_t i = 0, j = 0, k = 0;
    double ex = 0.0, ey = 0.0, ez = 0.0;
    if (!csv_index(fields, header, "i", i) ||
        !csv_index(fields, header, "j", j) ||
        !csv_index(fields, header, "k", k)) {
      throw std::runtime_error(
          "PBE: pbe_oracle_midpoint_epsilon.csv row is missing i,j,k.");
    }
    if (!csv_double(fields, header, {"eps_x", "epsilon_x", "ex"}, ex) ||
        !csv_double(fields, header, {"eps_y", "epsilon_y", "ey"}, ey) ||
        !csv_double(fields, header, {"eps_z", "epsilon_z", "ez"}, ez)) {
      throw std::runtime_error(
          "PBE: pbe_oracle_midpoint_epsilon.csv row is missing eps_x/eps_y/eps_z.");
    }
    if (i >= shape.nx || j >= shape.ny || k >= shape.nz)
      throw std::runtime_error(
          "PBE: pbe_oracle_midpoint_epsilon.csv index is outside the grid.");
    if (!std::isfinite(ex) || !std::isfinite(ey) || !std::isfinite(ez))
      throw std::runtime_error(
          "PBE: pbe_oracle_midpoint_epsilon.csv contains non-finite epsilon.");
    const std::size_t id = idx3(i, j, k, shape.nx, shape.ny);
    eps_x_[id] = ex;
    eps_y_[id] = ey;
    eps_z_[id] = ez;
    seen[id] = 1U;
    ++rows;
  }
  if (rows == 0)
    throw std::runtime_error(
        "PBE: pbe_oracle_midpoint_epsilon.csv has no rows.");
  if (std::any_of(seen.begin(), seen.end(),
                  [](unsigned char value) { return value == 0U; })) {
    throw std::runtime_error(
        "PBE: pbe_oracle_midpoint_epsilon.csv must provide one row for every grid point.");
  }
}

Analyzer::DelphiVdwMediaMap Analyzer::build_delphi_vdw_media_map(
    const GridShape &shape, const std::vector<Atom> &surface_atoms,
    TitanExecutor &executor) const {
  DelphiVdwMediaMap out;
  out.face_media.assign(shape.n * 3, 0);
  out.face_owner.assign(shape.n * 3, 0);
  out.audit.enabled = "raw";
  if (shape.n == 0 || surface_atoms.empty())
    return out;

  const double grid_spacing_extension = 0.5 * shape.h * std::sqrt(3.0);
  const double delphi_voxel_side =
      max_atom_radius(surface_atoms) + args_.ion_radius + grid_spacing_extension;
  const AtomSpatialBins raw_bins =
      build_atom_bins(surface_atoms, 0.0,
                      std::max(shape.h, delphi_voxel_side));

  auto face_index = [](std::size_t id, int axis) {
    return id * 3 + static_cast<std::size_t>(axis);
  };
  auto find_vdw_face_owner = [&](double x, double y,
                                 double z) -> std::int32_t {
    if (raw_bins.empty)
      return 0;
    const std::size_t bix = std::min<std::size_t>(
        raw_bins.bx - 1,
        static_cast<std::size_t>(std::max(
            0.0, std::floor((x - raw_bins.min_x) / raw_bins.cell))));
    const std::size_t biy = std::min<std::size_t>(
        raw_bins.by - 1,
        static_cast<std::size_t>(std::max(
            0.0, std::floor((y - raw_bins.min_y) / raw_bins.cell))));
    const std::size_t biz = std::min<std::size_t>(
        raw_bins.bz - 1,
        static_cast<std::size_t>(std::max(
            0.0, std::floor((z - raw_bins.min_z) / raw_bins.cell))));
    const std::size_t x0 = (bix == 0) ? 0 : bix - 1;
    const std::size_t y0 = (biy == 0) ? 0 : biy - 1;
    const std::size_t z0 = (biz == 0) ? 0 : biz - 1;
    const std::size_t x1 = std::min(raw_bins.bx - 1, bix + 1);
    const std::size_t y1 = std::min(raw_bins.by - 1, biy + 1);
    const std::size_t z1 = std::min(raw_bins.bz - 1, biz + 1);
    for (std::size_t kk = z0; kk <= z1; ++kk) {
      for (std::size_t jj = y0; jj <= y1; ++jj) {
        for (std::size_t ii = x0; ii <= x1; ++ii) {
          const auto &bucket =
              raw_bins.bins[(kk * raw_bins.by + jj) * raw_bins.bx + ii];
          for (std::size_t atom_idx : bucket) {
            const Atom &atom = surface_atoms[atom_idx];
            const double dx = x - atom.x;
            const double dy = y - atom.y;
            const double dz = z - atom.z;
            if (dx * dx + dy * dy + dz * dz < atom.radius * atom.radius)
              return static_cast<std::int32_t>(atom_idx + 2);
          }
        }
      }
    }
    return 0;
  };
  const std::size_t workers = effective_workers(executor, shape.nz);
  executor.parallel_for(
      shape.nz, workers, 1,
      [&](std::size_t, std::size_t z_begin, std::size_t z_end) {
        for (std::size_t k = z_begin; k < z_end; ++k) {
          const double z = shape.oz + static_cast<double>(k) * shape.h;
          for (std::size_t j = 0; j < shape.ny; ++j) {
            const double y = shape.oy + static_cast<double>(j) * shape.h;
            for (std::size_t i = 0; i < shape.nx; ++i) {
              const std::size_t id = idx3(i, j, k, shape.nx, shape.ny);
              const double x = shape.ox + static_cast<double>(i) * shape.h;
              if (i + 1 < shape.nx) {
                const std::int32_t owner =
                    find_vdw_face_owner(x + 0.5 * shape.h, y, z);
                if (owner != 0) {
                  out.face_media[face_index(id, 0)] = 1;
                  out.face_owner[face_index(id, 0)] = owner;
                }
              }
              if (j + 1 < shape.ny) {
                const std::int32_t owner =
                    find_vdw_face_owner(x, y + 0.5 * shape.h, z);
                if (owner != 0) {
                  out.face_media[face_index(id, 1)] = 1;
                  out.face_owner[face_index(id, 1)] = owner;
                }
              }
              if (k + 1 < shape.nz) {
                const std::int32_t owner =
                    find_vdw_face_owner(x, y, z + 0.5 * shape.h);
                if (owner != 0) {
                  out.face_media[face_index(id, 2)] = 1;
                  out.face_owner[face_index(id, 2)] = owner;
                }
              }
            }
          }
        }
      });

  const std::size_t plane = shape.nx * shape.ny;
  auto positive_face = [&](std::size_t i, std::size_t j, std::size_t k,
                           int axis) -> std::int32_t {
    if (i >= shape.nx || j >= shape.ny || k >= shape.nz)
      return 0;
    const std::size_t id = idx3(i, j, k, shape.nx, shape.ny);
    return out.face_media[face_index(id, axis)];
  };
  auto face_media6 = [&](std::size_t i, std::size_t j, std::size_t k,
                         std::array<std::int32_t, 6> &faces) {
    const std::size_t id = idx3(i, j, k, shape.nx, shape.ny);
    faces[0] = out.face_media[face_index(id, 0)];
    faces[1] = out.face_media[face_index(id, 1)];
    faces[2] = out.face_media[face_index(id, 2)];
    faces[3] = i > 0 ? out.face_media[face_index(id - 1, 0)] : 0;
    faces[4] = j > 0 ? out.face_media[face_index(id - shape.nx, 1)] : 0;
    faces[5] = k > 0 ? out.face_media[face_index(id - plane, 2)] : 0;
  };
  auto collect_boundary_voxels = [&]() {
    std::vector<std::size_t> voxels;
    if (shape.nx < 3 || shape.ny < 3 || shape.nz < 3)
      return voxels;
    std::array<std::int32_t, 6> faces{};
    for (std::size_t k = 1; k + 1 < shape.nz; ++k) {
      for (std::size_t j = 1; j + 1 < shape.ny; ++j) {
        for (std::size_t i = 1; i + 1 < shape.nx; ++i) {
          face_media6(i, j, k, faces);
          bool has_solute = false;
          bool has_solvent = false;
          for (std::int32_t value : faces) {
            has_solute = has_solute || value != 0;
            has_solvent = has_solvent || value == 0;
          }
          if (has_solute && has_solvent)
            voxels.push_back(idx3(i, j, k, shape.nx, shape.ny));
        }
      }
    }
    return voxels;
  };

  std::vector<std::size_t> boundary_voxels = collect_boundary_voxels();
  out.audit.boundary_initial = boundary_voxels.size();
  if (args_.test_vdwms_remap == "off" || args_.probe_radius <= 0.0) {
    out.audit.enabled = args_.test_vdwms_remap == "off" ? "off" : "zero_probe";
    out.audit.boundary_final = out.audit.boundary_initial;
    return out;
  }

  out.audit.enabled = "auto";
  const double shell = std::min(args_.probe_radius, 0.5 * shape.h);
  const double shell_sq_eps = 1.0e-24;
  const AtomSpatialBins shell_bins =
      build_atom_bins(surface_atoms, shell, std::max(shape.h, delphi_voxel_side));
  const AtomSpatialBins probe_bins =
      build_atom_bins(surface_atoms, args_.probe_radius,
                      std::max(shape.h, delphi_voxel_side));
  auto candidate_shell_atom = [&](double x, double y,
                                  double z) -> std::ptrdiff_t {
    if (shell_bins.empty)
      return -1;
    const int radius_bins = static_cast<int>(
        std::ceil(std::max(shell_bins.max_radius, shape.h) /
                  std::max(1.0e-12, shell_bins.cell))) +
                            1;
    const int bx0 = static_cast<int>(
        std::floor((x - shell_bins.min_x) / shell_bins.cell));
    const int by0 = static_cast<int>(
        std::floor((y - shell_bins.min_y) / shell_bins.cell));
    const int bz0 = static_cast<int>(
        std::floor((z - shell_bins.min_z) / shell_bins.cell));
    std::ptrdiff_t best_atom = -1;
    double best_surface_distance = std::numeric_limits<double>::infinity();
    for (int dz = -radius_bins; dz <= radius_bins; ++dz) {
      const int bz = bz0 + dz;
      if (bz < 0 || bz >= static_cast<int>(shell_bins.bz))
        continue;
      for (int dy = -radius_bins; dy <= radius_bins; ++dy) {
        const int by = by0 + dy;
        if (by < 0 || by >= static_cast<int>(shell_bins.by))
          continue;
        for (int dx = -radius_bins; dx <= radius_bins; ++dx) {
          const int bx = bx0 + dx;
          if (bx < 0 || bx >= static_cast<int>(shell_bins.bx))
            continue;
          const auto &bucket =
              shell_bins.bins[(static_cast<std::size_t>(bz) * shell_bins.by +
                               static_cast<std::size_t>(by)) *
                                  shell_bins.bx +
                              static_cast<std::size_t>(bx)];
          for (std::size_t atom_idx : bucket) {
            const Atom &atom = surface_atoms[atom_idx];
            const double radius = atom.radius + shell;
            const double ddx = x - atom.x;
            const double ddy = y - atom.y;
            const double ddz = z - atom.z;
            const double r2 = ddx * ddx + ddy * ddy + ddz * ddz;
            if (r2 < radius * radius - shell_sq_eps) {
              const double surface_distance = std::sqrt(r2) - atom.radius;
              if (surface_distance < best_surface_distance) {
                best_surface_distance = surface_distance;
                best_atom = static_cast<std::ptrdiff_t>(atom_idx);
              }
            }
          }
        }
      }
    }
    return best_atom;
  };
  const double core_shell =
      std::min(args_.probe_radius, 0.25 * std::sqrt(3.0) * shape.h);
  const AtomSpatialBins core_bins =
      build_atom_bins(surface_atoms, core_shell,
                      std::max(shape.h, delphi_voxel_side));
  auto in_core_shell = [&](double x, double y, double z) {
    if (core_bins.empty)
      return false;
    const int radius_bins = static_cast<int>(
        std::ceil(std::max(core_bins.max_radius, shape.h) /
                  std::max(1.0e-12, core_bins.cell))) +
                            1;
    const int bx0 = static_cast<int>(
        std::floor((x - core_bins.min_x) / core_bins.cell));
    const int by0 = static_cast<int>(
        std::floor((y - core_bins.min_y) / core_bins.cell));
    const int bz0 = static_cast<int>(
        std::floor((z - core_bins.min_z) / core_bins.cell));
    for (int dz = -radius_bins; dz <= radius_bins; ++dz) {
      const int bz = bz0 + dz;
      if (bz < 0 || bz >= static_cast<int>(core_bins.bz))
        continue;
      for (int dy = -radius_bins; dy <= radius_bins; ++dy) {
        const int by = by0 + dy;
        if (by < 0 || by >= static_cast<int>(core_bins.by))
          continue;
        for (int dx = -radius_bins; dx <= radius_bins; ++dx) {
          const int bx = bx0 + dx;
          if (bx < 0 || bx >= static_cast<int>(core_bins.bx))
            continue;
          const auto &bucket =
              core_bins.bins[(static_cast<std::size_t>(bz) * core_bins.by +
                              static_cast<std::size_t>(by)) *
                                 core_bins.bx +
                             static_cast<std::size_t>(bx)];
          for (std::size_t atom_idx : bucket) {
            const Atom &atom = surface_atoms[atom_idx];
            const double radius = atom.radius + core_shell;
            const double ddx = x - atom.x;
            const double ddy = y - atom.y;
            const double ddz = z - atom.z;
            const double r2 = ddx * ddx + ddy * ddy + ddz * ddz;
            if (r2 < radius * radius - shell_sq_eps)
              return true;
          }
        }
      }
    }
    return false;
  };
  auto probe_surface_occluded = [&](double x, double y, double z,
                                    std::size_t closest_atom_idx) {
    if (probe_bins.empty || closest_atom_idx >= surface_atoms.size())
      return false;
    const Atom &closest = surface_atoms[closest_atom_idx];
    const double vx = x - closest.x;
    const double vy = y - closest.y;
    const double vz = z - closest.z;
    const double distance = std::sqrt(vx * vx + vy * vy + vz * vz);
    if (distance <= 1.0e-12)
      return true;
    const double surface_radius = closest.radius + args_.probe_radius;
    const double sx = closest.x + vx * (surface_radius / distance);
    const double sy = closest.y + vy * (surface_radius / distance);
    const double sz = closest.z + vz * (surface_radius / distance);
    const int radius_bins = static_cast<int>(
        std::ceil(std::max(probe_bins.max_radius, shape.h) /
                  std::max(1.0e-12, probe_bins.cell))) +
                            1;
    const int bx0 = static_cast<int>(
        std::floor((sx - probe_bins.min_x) / probe_bins.cell));
    const int by0 = static_cast<int>(
        std::floor((sy - probe_bins.min_y) / probe_bins.cell));
    const int bz0 = static_cast<int>(
        std::floor((sz - probe_bins.min_z) / probe_bins.cell));
    constexpr double kSasShrink = 0.99999;
    for (int dz = -radius_bins; dz <= radius_bins; ++dz) {
      const int bz = bz0 + dz;
      if (bz < 0 || bz >= static_cast<int>(probe_bins.bz))
        continue;
      for (int dy = -radius_bins; dy <= radius_bins; ++dy) {
        const int by = by0 + dy;
        if (by < 0 || by >= static_cast<int>(probe_bins.by))
          continue;
        for (int dx = -radius_bins; dx <= radius_bins; ++dx) {
          const int bx = bx0 + dx;
          if (bx < 0 || bx >= static_cast<int>(probe_bins.bx))
            continue;
          const auto &bucket =
              probe_bins.bins[(static_cast<std::size_t>(bz) * probe_bins.by +
                               static_cast<std::size_t>(by)) *
                                  probe_bins.bx +
                              static_cast<std::size_t>(bx)];
          for (std::size_t atom_idx : bucket) {
            if (atom_idx == closest_atom_idx)
              continue;
            const Atom &atom = surface_atoms[atom_idx];
            const double radius = atom.radius + args_.probe_radius;
            const double ddx = sx - atom.x;
            const double ddy = sy - atom.y;
            const double ddz = sz - atom.z;
            const double r2 = ddx * ddx + ddy * ddy + ddz * ddz;
            if (r2 < radius * radius * kSasShrink)
              return true;
          }
        }
      }
    }
    return false;
  };
  struct ExposedPointBins {
    bool empty{true};
    double min_x{0.0};
    double min_y{0.0};
    double min_z{0.0};
    double cell{1.0};
    std::size_t bx{0};
    std::size_t by{0};
    std::size_t bz{0};
    std::vector<std::array<double, 3>> points;
    std::vector<std::vector<std::size_t>> bins;
  };
  auto build_exposed_point_bins = [&]() {
    ExposedPointBins bins;
    if (surface_atoms.empty() || args_.probe_radius <= 0.0)
      return bins;
    constexpr std::size_t kSurfaceSamples = 520;
    constexpr double kGoldenAngle = 2.3999632297286533222;
    bins.points.reserve(surface_atoms.size() * 64);
    for (std::size_t atom_idx = 0; atom_idx < surface_atoms.size(); ++atom_idx) {
      const Atom &atom = surface_atoms[atom_idx];
      const double radius = atom.radius + args_.probe_radius;
      if (radius <= 0.0)
        continue;
      for (std::size_t sample = 0; sample < kSurfaceSamples; ++sample) {
        const double t =
            (static_cast<double>(sample) + 0.5) /
            static_cast<double>(kSurfaceSamples);
        const double z_unit = 1.0 - 2.0 * t;
        const double radial =
            std::sqrt(std::max(0.0, 1.0 - z_unit * z_unit));
        const double theta = kGoldenAngle * static_cast<double>(sample);
        const double sx = atom.x + radius * radial * std::cos(theta);
        const double sy = atom.y + radius * radial * std::sin(theta);
        const double sz = atom.z + radius * z_unit;
        if (!probe_surface_occluded(sx, sy, sz, atom_idx))
          bins.points.push_back({sx, sy, sz});
      }
    }
    if (bins.points.empty())
      return bins;
    double max_x = bins.points.front()[0];
    double max_y = bins.points.front()[1];
    double max_z = bins.points.front()[2];
    bins.min_x = max_x;
    bins.min_y = max_y;
    bins.min_z = max_z;
    for (const auto &point : bins.points) {
      bins.min_x = std::min(bins.min_x, point[0]);
      bins.min_y = std::min(bins.min_y, point[1]);
      bins.min_z = std::min(bins.min_z, point[2]);
      max_x = std::max(max_x, point[0]);
      max_y = std::max(max_y, point[1]);
      max_z = std::max(max_z, point[2]);
    }
    bins.cell = std::max(shape.h, args_.probe_radius);
    bins.bx = std::max<std::size_t>(
        1, static_cast<std::size_t>(
               std::floor((max_x - bins.min_x) / bins.cell)) +
               1);
    bins.by = std::max<std::size_t>(
        1, static_cast<std::size_t>(
               std::floor((max_y - bins.min_y) / bins.cell)) +
               1);
    bins.bz = std::max<std::size_t>(
        1, static_cast<std::size_t>(
               std::floor((max_z - bins.min_z) / bins.cell)) +
               1);
    bins.bins.assign(bins.bx * bins.by * bins.bz, {});
    for (std::size_t point_idx = 0; point_idx < bins.points.size();
         ++point_idx) {
      const auto &point = bins.points[point_idx];
      const std::size_t ix = std::min<std::size_t>(
          bins.bx - 1,
          static_cast<std::size_t>(
              std::max(0.0, std::floor((point[0] - bins.min_x) / bins.cell))));
      const std::size_t iy = std::min<std::size_t>(
          bins.by - 1,
          static_cast<std::size_t>(
              std::max(0.0, std::floor((point[1] - bins.min_y) / bins.cell))));
      const std::size_t iz = std::min<std::size_t>(
          bins.bz - 1,
          static_cast<std::size_t>(
              std::max(0.0, std::floor((point[2] - bins.min_z) / bins.cell))));
      bins.bins[(iz * bins.by + iy) * bins.bx + ix].push_back(point_idx);
    }
    bins.empty = false;
    return bins;
  };
  const ExposedPointBins exposed_bins = build_exposed_point_bins();
  auto near_exposed_point = [&](double x, double y, double z) {
    if (exposed_bins.empty)
      return false;
    const double probe_sq = args_.probe_radius * args_.probe_radius;
    const int radius_bins = static_cast<int>(
        std::ceil(args_.probe_radius / std::max(1.0e-12, exposed_bins.cell))) +
                            1;
    const int bx0 = static_cast<int>(
        std::floor((x - exposed_bins.min_x) / exposed_bins.cell));
    const int by0 = static_cast<int>(
        std::floor((y - exposed_bins.min_y) / exposed_bins.cell));
    const int bz0 = static_cast<int>(
        std::floor((z - exposed_bins.min_z) / exposed_bins.cell));
    for (int dz = -radius_bins; dz <= radius_bins; ++dz) {
      const int bz = bz0 + dz;
      if (bz < 0 || bz >= static_cast<int>(exposed_bins.bz))
        continue;
      for (int dy = -radius_bins; dy <= radius_bins; ++dy) {
        const int by = by0 + dy;
        if (by < 0 || by >= static_cast<int>(exposed_bins.by))
          continue;
        for (int dx = -radius_bins; dx <= radius_bins; ++dx) {
          const int bx = bx0 + dx;
          if (bx < 0 || bx >= static_cast<int>(exposed_bins.bx))
            continue;
          const auto &bucket =
              exposed_bins.bins[(static_cast<std::size_t>(bz) *
                                     exposed_bins.by +
                                 static_cast<std::size_t>(by)) *
                                    exposed_bins.bx +
                                static_cast<std::size_t>(bx)];
          for (std::size_t point_idx : bucket) {
            const auto &point = exposed_bins.points[point_idx];
            const double ddx = x - point[0];
            const double ddy = y - point[1];
            const double ddz = z - point[2];
            const double r2 = ddx * ddx + ddy * ddy + ddz * ddz;
            if (r2 < probe_sq)
              return true;
          }
        }
      }
    }
    return false;
  };

  std::size_t added = 0;
  std::size_t exposed = 0;
  std::vector<unsigned char> processed_faces(shape.n * 3, 0U);
  auto process_face = [&](std::size_t face_id, int axis, double x, double y,
                          double z) {
    if (face_id >= shape.n || axis < 0 || axis > 2)
      return;
    const std::size_t fid = face_index(face_id, axis);
    if (processed_faces[fid] != 0U)
      return;
    processed_faces[fid] = 1U;
    if (out.face_media[fid] != 0)
      return;
    const std::ptrdiff_t atom_idx = candidate_shell_atom(x, y, z);
    if (atom_idx >= 0 && (in_core_shell(x, y, z) || !near_exposed_point(x, y, z))) {
      out.face_media[fid] = 1;
      out.face_owner[fid] = static_cast<std::int32_t>(atom_idx + 2);
      ++added;
    } else {
      ++exposed;
    }
  };
  for (std::size_t center_id : boundary_voxels) {
    const std::size_t k = center_id / plane;
    const std::size_t rem = center_id - k * plane;
    const std::size_t j = rem / shape.nx;
    const std::size_t i = rem - j * shape.nx;
    const double x = shape.ox + static_cast<double>(i) * shape.h;
    const double y = shape.oy + static_cast<double>(j) * shape.h;
    const double z = shape.oz + static_cast<double>(k) * shape.h;
    if (i + 1 < shape.nx)
      process_face(center_id, 0, x + 0.5 * shape.h, y, z);
    if (j + 1 < shape.ny)
      process_face(center_id, 1, x, y + 0.5 * shape.h, z);
    if (k + 1 < shape.nz)
      process_face(center_id, 2, x, y, z + 0.5 * shape.h);
    if (i > 0)
      process_face(center_id - 1, 0, x - 0.5 * shape.h, y, z);
    if (j > 0)
      process_face(center_id - shape.nx, 1, x, y - 0.5 * shape.h, z);
    if (k > 0)
      process_face(center_id - plane, 2, x, y, z - 0.5 * shape.h);
  }
  out.audit.added = added;
  out.audit.removed = 0;
  out.audit.exposed_points = exposed_bins.points.empty()
                                  ? exposed
                                  : exposed_bins.points.size();
  boundary_voxels = collect_boundary_voxels();
  out.audit.boundary_final = boundary_voxels.size();
  const double face_count =
      static_cast<double>(std::max<std::size_t>(1, shape.n * 3));
  out.audit.midpoint_rms =
      (args_.ext_diel - args_.int_diel) *
      std::sqrt(static_cast<double>(added) / face_count);
  (void)positive_face;
  return out;
}

void Analyzer::apply_delphi_vdw_media_map(
    const GridShape &shape, const DelphiVdwMediaMap &media_map, double eps_in,
    double eps_out) {
  if (media_map.face_media.size() != shape.n * 3)
    throw std::runtime_error("PBE: internal DelPhi VDW media map has wrong size.");
  for (std::size_t id = 0; id < shape.n; ++id) {
    eps_x_[id] = media_map.face_media[id * 3] != 0 ? eps_in : eps_out;
    eps_y_[id] = media_map.face_media[id * 3 + 1] != 0 ? eps_in : eps_out;
    eps_z_[id] = media_map.face_media[id * 3 + 2] != 0 ? eps_in : eps_out;
  }
  last_vdwms_audit_ = media_map.audit;
}

Analyzer::BoundaryDiagnostics
Analyzer::apply_oracle_boundary_map(const GridShape &shape) {
  BoundaryDiagnostics diag;
  if (args_.test_boundary_map != "oracle" || !oracle_replay_applies(shape))
    return diag;
  std::ifstream in("pbe_oracle_boundary_map.csv");
  if (!in.is_open()) {
    throw std::runtime_error(
        "PBE: -pbe-test-boundary-map=oracle requires "
        "pbe_oracle_boundary_map.csv in the run directory.");
  }
  std::string line;
  if (!std::getline(in, line))
    throw std::runtime_error("PBE: empty pbe_oracle_boundary_map.csv.");
  const auto header = csv_header_map(split_csv_simple(line));
  std::vector<unsigned char> seen(shape.n, 0U);
  std::size_t rows = 0;
  while (std::getline(in, line)) {
    if (trim_copy(line).empty())
      continue;
    const auto fields = split_csv_simple(line);
    std::size_t i = 0, j = 0, k = 0;
    double value = 0.0;
    if (!csv_index(fields, header, "i", i) ||
        !csv_index(fields, header, "j", j) ||
        !csv_index(fields, header, "k", k)) {
      throw std::runtime_error(
          "PBE: pbe_oracle_boundary_map.csv row is missing i,j,k.");
    }
    if (!csv_double(fields, header, {"phi", "potential", "boundary_phi"},
                    value)) {
      throw std::runtime_error(
          "PBE: pbe_oracle_boundary_map.csv row is missing phi/potential.");
    }
    if (i >= shape.nx || j >= shape.ny || k >= shape.nz)
      throw std::runtime_error(
          "PBE: pbe_oracle_boundary_map.csv index is outside the grid.");
    if (!is_boundary(i, j, k, shape))
      throw std::runtime_error(
          "PBE: pbe_oracle_boundary_map.csv contains non-boundary indices.");
    const std::size_t id = idx3(i, j, k, shape.nx, shape.ny);
    phi_[id] = value;
    seen[id] = 1U;
    ++rows;
  }
  std::size_t expected = 0;
  double sumsq = 0.0;
  double max_abs = 0.0;
  for (std::size_t k = 0; k < shape.nz; ++k) {
    for (std::size_t j = 0; j < shape.ny; ++j) {
      for (std::size_t i = 0; i < shape.nx; ++i) {
        if (!is_boundary(i, j, k, shape))
          continue;
        ++expected;
        const std::size_t id = idx3(i, j, k, shape.nx, shape.ny);
        if (!seen[id]) {
          throw std::runtime_error(
              "PBE: pbe_oracle_boundary_map.csv must provide every "
              "boundary grid point for replay.");
        }
        sumsq += phi_[id] * phi_[id];
        max_abs = std::max(max_abs, std::abs(phi_[id]));
      }
    }
  }
  if (rows == 0 || expected == 0)
    throw std::runtime_error("PBE: pbe_oracle_boundary_map.csv has no rows.");
  diag.rms = std::sqrt(sumsq / static_cast<double>(expected));
  diag.max_abs = max_abs;
  return diag;
}

void Analyzer::build_coefficients(const GridShape &shape,
                                  const std::vector<Atom> &atoms,
                                  double eps_in, double eps_out, double kappa,
                                  TitanExecutor &executor) {
  const double kappa2h2 = kappa * kappa * shape.h * shape.h;
  const std::vector<Atom> surface_atoms = adjusted_surface_atoms(atoms);
  auto center_eps = [&](std::size_t id) {
    return inside_[id] ? eps_in : eps_out;
  };
  const bool refine_faces = surface_uses_face_refinement(args_.surface);
  const bool delphi_vdw_faces = args_.surface == "delphi-vdw";
  const double surface_inflate =
      surface_uses_probe_exclusion(args_.surface) ? args_.probe_radius : 0.0;
  const AtomSpatialBins face_bins =
      refine_faces ? build_atom_bins(surface_atoms, surface_inflate, shape.h)
                   : AtomSpatialBins{};
  auto segment_face_eps = [&](double x0, double y0, double z0, double x1,
                              double y1, double z1) {
    const double inside_fraction = segment_inside_fraction(
        face_bins, surface_atoms, x0, y0, z0, x1, y1, z1, surface_inflate);
    return dielectric_from_inside_fraction(inside_fraction, eps_in, eps_out);
  };
  auto sampled_face_eps = [&](double x, double y, double z, int axis) {
    int inside_samples = 0;
    const double q = 0.25 * shape.h;
    for (int a = -1; a <= 1; a += 2) {
      for (int b = -1; b <= 1; b += 2) {
        double sx = x;
        double sy = y;
        double sz = z;
        if (axis == 0) {
          sy += static_cast<double>(a) * q;
          sz += static_cast<double>(b) * q;
        } else if (axis == 1) {
          sx += static_cast<double>(a) * q;
          sz += static_cast<double>(b) * q;
        } else {
          sx += static_cast<double>(a) * q;
          sy += static_cast<double>(b) * q;
        }
        if (point_inside_inflated_atom(face_bins, surface_atoms, sx, sy, sz,
                                       surface_inflate))
          ++inside_samples;
      }
    }
    const double frac_inside = static_cast<double>(inside_samples) / 4.0;
    return frac_inside * eps_in + (1.0 - frac_inside) * eps_out;
  };

  if (uses_gaussian_model()) {
    double max_base_radius = shape.h;
    for (const Atom &atom : surface_atoms)
      max_base_radius = std::max(max_base_radius, atom.radius);
    const double sigma = std::max(1.0e-12, args_.gaussian_sigma);
    const double influence = 3.0 * sigma * max_base_radius;
    const AtomSpatialBins gaussian_bins =
        build_atom_bins(surface_atoms, std::max(0.0, influence), shape.h);
    auto density_at = [&](double x, double y, double z) {
      const double exponent = std::max(1.0e-12, args_.gaussian_exp);
      const int radius_bins = static_cast<int>(
          std::ceil(std::max(shape.h, influence) / gaussian_bins.cell)) + 1;
      const int bx0 = static_cast<int>(std::floor((x - gaussian_bins.min_x) /
                                                  gaussian_bins.cell));
      const int by0 = static_cast<int>(std::floor((y - gaussian_bins.min_y) /
                                                  gaussian_bins.cell));
      const int bz0 = static_cast<int>(std::floor((z - gaussian_bins.min_z) /
                                                  gaussian_bins.cell));
      double density = 0.0;
      for (int dz = -radius_bins; dz <= radius_bins; ++dz) {
        const int bz = bz0 + dz;
        if (bz < 0 || bz >= static_cast<int>(gaussian_bins.bz))
          continue;
        for (int dy = -radius_bins; dy <= radius_bins; ++dy) {
          const int by = by0 + dy;
          if (by < 0 || by >= static_cast<int>(gaussian_bins.by))
            continue;
          for (int dx = -radius_bins; dx <= radius_bins; ++dx) {
            const int bx = bx0 + dx;
            if (bx < 0 || bx >= static_cast<int>(gaussian_bins.bx))
              continue;
            const auto &bucket = gaussian_bins.bins
                [(static_cast<std::size_t>(bz) * gaussian_bins.by +
                  static_cast<std::size_t>(by)) *
                     gaussian_bins.bx +
                 static_cast<std::size_t>(bx)];
            for (std::size_t atom_idx : bucket) {
              const Atom &atom = surface_atoms[atom_idx];
              const double r0 = std::max(1.0e-12, atom.radius);
              const double cutoff = 3.0 * sigma * r0;
              const double ddx = x - atom.x;
              const double ddy = y - atom.y;
              const double ddz = z - atom.z;
              const double r2 = ddx * ddx + ddy * ddy + ddz * ddz;
              if (r2 > cutoff * cutoff)
                continue;
              const double scaled = std::sqrt(r2) / (sigma * r0);
              density += std::exp(-std::pow(scaled, exponent));
              if (density >= 1.0)
                return 1.0;
            }
          }
        }
      }
      return std::clamp(density, 0.0, 1.0);
    };
    for (std::size_t k = 0; k < shape.nz; ++k) {
      for (std::size_t j = 0; j < shape.ny; ++j) {
        for (std::size_t i = 0; i < shape.nx; ++i) {
          const std::size_t id = idx3(i, j, k, shape.nx, shape.ny);
          const double x = shape.ox + static_cast<double>(i) * shape.h;
          const double y = shape.oy + static_cast<double>(j) * shape.h;
          const double z = shape.oz + static_cast<double>(k) * shape.h;
          const double e0 =
              gaussian_dielectric(density_at(x, y, z), eps_in, eps_out);
          eps_x_[id] = e0;
          eps_y_[id] = e0;
          eps_z_[id] = e0;
          if (i + 1 < shape.nx) {
            eps_x_[id] = gaussian_dielectric(
                density_at(x + 0.5 * shape.h, y, z),
                eps_in, eps_out);
          }
          if (j + 1 < shape.ny) {
            eps_y_[id] = gaussian_dielectric(
                density_at(x, y + 0.5 * shape.h, z),
                eps_in, eps_out);
          }
          if (k + 1 < shape.nz) {
            eps_z_[id] = gaussian_dielectric(
                density_at(x, y, z + 0.5 * shape.h),
                eps_in, eps_out);
          }
        }
      }
    }
    double eps_min = std::numeric_limits<double>::max();
    double eps_max = std::numeric_limits<double>::lowest();
    double eps_sum = 0.0;
    std::size_t eps_count = 0;
    for (std::size_t id = 0; id < shape.n; ++id) {
      const double values[3] = {eps_x_[id], eps_y_[id], eps_z_[id]};
      for (double eps_value : values) {
        if (!std::isfinite(eps_value))
          continue;
        eps_min = std::min(eps_min, eps_value);
        eps_max = std::max(eps_max, eps_value);
        eps_sum += eps_value;
        ++eps_count;
      }
    }
    last_gaussian_eps_min_ = eps_count > 0 ? eps_min : 0.0;
    last_gaussian_eps_max_ = eps_count > 0 ? eps_max : 0.0;
    last_gaussian_eps_mean_ =
        eps_count > 0 ? eps_sum / static_cast<double>(eps_count) : 0.0;
    for (std::size_t k = 1; k + 1 < shape.nz; ++k) {
      for (std::size_t j = 1; j + 1 < shape.ny; ++j) {
        for (std::size_t i = 1; i + 1 < shape.nx; ++i) {
          const std::size_t id = idx3(i, j, k, shape.nx, shape.ny);
          const double screen = ion_accessible_[id] ? kappa2h2 : 0.0;
          screen_[id] = screen;
          denom_[id] = eps_x_[id] + eps_x_[id - 1] + eps_y_[id] +
                       eps_y_[id - shape.nx] + eps_z_[id] +
                       eps_z_[id - shape.nx * shape.ny] + screen;
        }
      }
    }
    return;
  }

  if (delphi_vdw_faces) {
    const DelphiVdwMediaMap media_map =
        build_delphi_vdw_media_map(shape, surface_atoms, executor);
    apply_delphi_vdw_media_map(shape, media_map, eps_in, eps_out);

    apply_oracle_midpoint_epsilon_map(shape);
    for (std::size_t k = 1; k + 1 < shape.nz; ++k) {
      for (std::size_t j = 1; j + 1 < shape.ny; ++j) {
        for (std::size_t i = 1; i + 1 < shape.nx; ++i) {
          const std::size_t id = idx3(i, j, k, shape.nx, shape.ny);
          const double screen = ion_accessible_[id] ? kappa2h2 : 0.0;
          screen_[id] = screen;
          denom_[id] = eps_x_[id] + eps_x_[id - 1] + eps_y_[id] +
                       eps_y_[id - shape.nx] + eps_z_[id] +
                       eps_z_[id - shape.nx * shape.ny] + screen;
        }
      }
    }
    return;
  }

  for (std::size_t k = 0; k < shape.nz; ++k) {
    for (std::size_t j = 0; j < shape.ny; ++j) {
      for (std::size_t i = 0; i < shape.nx; ++i) {
        const std::size_t id = idx3(i, j, k, shape.nx, shape.ny);
        const double e0 = center_eps(id);
        eps_x_[id] = e0;
        eps_y_[id] = e0;
        eps_z_[id] = e0;
        if (i + 1 < shape.nx) {
          const double e1 = center_eps(id + 1);
          const double x_face =
              shape.ox + (static_cast<double>(i) + 0.5) * shape.h;
          const double y_face = shape.oy + static_cast<double>(j) * shape.h;
          const double z_face = shape.oz + static_cast<double>(k) * shape.h;
          eps_x_[id] =
              (delphi_vdw_faces && inside_[id] != inside_[id + 1])
                  ? segment_face_eps(shape.ox + static_cast<double>(i) *
                                                    shape.h,
                                     y_face, z_face,
                                     shape.ox + static_cast<double>(i + 1) *
                                                    shape.h,
                                     y_face, z_face)
              : (refine_faces && inside_[id] != inside_[id + 1])
                  ? sampled_face_eps(x_face, y_face, z_face, 0)
                  : harmonic_eps(e0, e1);
        }
        if (j + 1 < shape.ny) {
          const double e1 = center_eps(id + shape.nx);
          const double x_face = shape.ox + static_cast<double>(i) * shape.h;
          const double y_face =
              shape.oy + (static_cast<double>(j) + 0.5) * shape.h;
          const double z_face = shape.oz + static_cast<double>(k) * shape.h;
          eps_y_[id] =
              (delphi_vdw_faces && inside_[id] != inside_[id + shape.nx])
                  ? segment_face_eps(x_face,
                                     shape.oy + static_cast<double>(j) *
                                                    shape.h,
                                     z_face, x_face,
                                     shape.oy + static_cast<double>(j + 1) *
                                                    shape.h,
                                     z_face)
              : (refine_faces && inside_[id] != inside_[id + shape.nx])
                  ? sampled_face_eps(x_face, y_face, z_face, 1)
                  : harmonic_eps(e0, e1);
        }
        if (k + 1 < shape.nz) {
          const std::size_t dz = shape.nx * shape.ny;
          const double e1 = center_eps(id + dz);
          const double x_face = shape.ox + static_cast<double>(i) * shape.h;
          const double y_face = shape.oy + static_cast<double>(j) * shape.h;
          const double z_face =
              shape.oz + (static_cast<double>(k) + 0.5) * shape.h;
          eps_z_[id] =
              (delphi_vdw_faces && inside_[id] != inside_[id + dz])
                  ? segment_face_eps(x_face, y_face,
                                     shape.oz + static_cast<double>(k) *
                                                    shape.h,
                                     x_face, y_face,
                                     shape.oz + static_cast<double>(k + 1) *
                                                    shape.h)
              : (refine_faces && inside_[id] != inside_[id + dz])
                  ? sampled_face_eps(x_face, y_face, z_face, 2)
                  : harmonic_eps(e0, e1);
        }
      }
    }
  }

  apply_oracle_midpoint_epsilon_map(shape);
  for (std::size_t k = 1; k + 1 < shape.nz; ++k) {
    for (std::size_t j = 1; j + 1 < shape.ny; ++j) {
      for (std::size_t i = 1; i + 1 < shape.nx; ++i) {
        const std::size_t id = idx3(i, j, k, shape.nx, shape.ny);
        const double screen = ion_accessible_[id] ? kappa2h2 : 0.0;
        screen_[id] = screen;
        denom_[id] = eps_x_[id] + eps_x_[id - 1] + eps_y_[id] +
                     eps_y_[id - shape.nx] + eps_z_[id] +
                     eps_z_[id - shape.nx * shape.ny] + screen;
      }
    }
  }
}

Analyzer::BoundaryDiagnostics
Analyzer::set_boundary(const GridShape &shape, const std::vector<Atom> &atoms,
                       double eps_out, double kappa) {
  std::fill(phi_.begin(), phi_.end(), 0.0);
  BoundaryDiagnostics diag;
  if (args_.boundary == "zero")
    return diag;

  double total_charge = 0.0;
  double cx = 0.0;
  double cy = 0.0;
  double cz = 0.0;
  for (const Atom &atom : atoms) {
    total_charge += atom.q;
    cx += atom.q * atom.x;
    cy += atom.q * atom.y;
    cz += atom.q * atom.z;
  }
  const double center_x = shape.ox + 0.5 * static_cast<double>(shape.nx - 1) * shape.h;
  const double center_y = shape.oy + 0.5 * static_cast<double>(shape.ny - 1) * shape.h;
  const double center_z = shape.oz + 0.5 * static_cast<double>(shape.nz - 1) * shape.h;
  if (std::abs(total_charge) > 1.0e-12) {
    cx /= total_charge;
    cy /= total_charge;
    cz /= total_charge;
  } else {
    cx = center_x;
    cy = center_y;
    cz = center_z;
  }
  double px = 0.0;
  double py = 0.0;
  double pz = 0.0;
  for (const Atom &atom : atoms) {
    px += atom.q * (atom.x - cx);
    py += atom.q * (atom.y - cy);
    pz += atom.q * (atom.z - cz);
  }
  double sumsq = 0.0;
  std::size_t count = 0;
  for (std::size_t k = 0; k < shape.nz; ++k) {
    const double z = shape.oz + static_cast<double>(k) * shape.h;
    for (std::size_t j = 0; j < shape.ny; ++j) {
      const double y = shape.oy + static_cast<double>(j) * shape.h;
      for (std::size_t i = 0; i < shape.nx; ++i) {
        if (!is_boundary(i, j, k, shape))
          continue;
        const double x = shape.ox + static_cast<double>(i) * shape.h;
        double potential = 0.0;
        if (args_.boundary == "dipolar") {
          const double dx = x - cx;
          const double dy = y - cy;
          const double dz = z - cz;
          const double r2 = dx * dx + dy * dy + dz * dz;
          const double r = std::sqrt(r2);
          if (r > 1.0e-10) {
            const double screened = std::exp(-kappa * r);
            const double monopole = total_charge / r;
            const double dipole = (px * dx + py * dy + pz * dz) /
                                  std::max(1.0e-30, r2 * r);
            potential = screened * (monopole + dipole) /
                        std::max(1.0e-12, eps_out);
          }
        } else {
          for (const Atom &atom : atoms) {
            const double dx = x - atom.x;
            const double dy = y - atom.y;
            const double dz = z - atom.z;
            const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
            if (r > 1.0e-10) {
              potential += atom.q * std::exp(-kappa * r) /
                           (std::max(1.0e-12, eps_out) * r);
            }
          }
        }
        phi_[idx3(i, j, k, shape.nx, shape.ny)] = potential;
        sumsq += potential * potential;
        diag.max_abs = std::max(diag.max_abs, std::abs(potential));
        ++count;
      }
    }
  }
  diag.rms = count > 0 ? std::sqrt(sumsq / static_cast<double>(count)) : 0.0;
  return diag;
}

Analyzer::SolveStats Analyzer::solve_sor(const GridShape &shape,
                                          TitanExecutor &executor) {
  SolveStats stats;
  const double omega =
      args_.omega == "auto"
          ? std::clamp(2.0 / (1.0 + std::sin(kPi / static_cast<double>(
                                             std::max({shape.nx, shape.ny,
                                                       shape.nz})))),
                       1.0, 1.95)
          : args_.omega_value;
  const std::size_t interior_z = shape.nz > 2 ? shape.nz - 2 : 0;
  const std::size_t workers = effective_workers(executor, interior_z);
  const double rhs = rhs_scale(shape);

  std::vector<double> local_sumsq(workers, 0.0);
  std::vector<double> local_max(workers, 0.0);
  std::vector<std::size_t> local_count(workers, 0);

  for (int iter = 1; iter <= args_.max_iter; ++iter) {
    std::fill(local_sumsq.begin(), local_sumsq.end(), 0.0);
    std::fill(local_max.begin(), local_max.end(), 0.0);
    std::fill(local_count.begin(), local_count.end(), 0);

    for (int color = 0; color < 2; ++color) {
      executor.parallel_for(
          interior_z, workers, 1,
          [&](std::size_t worker, std::size_t begin, std::size_t end) {
            double sumsq = 0.0;
            double max_delta = 0.0;
            std::size_t count = 0;
            for (std::size_t kk = begin; kk < end; ++kk) {
              const std::size_t k = kk + 1;
              for (std::size_t j = 1; j + 1 < shape.ny; ++j) {
                for (std::size_t i = 1; i + 1 < shape.nx; ++i) {
                  if (static_cast<int>((i + j + k) & 1ULL) != color)
                    continue;
                  const std::size_t id = idx3(i, j, k, shape.nx, shape.ny);
                  const double neighbor_sum =
                      eps_x_[id] * phi_[id + 1] + eps_x_[id - 1] * phi_[id - 1] +
                      eps_y_[id] * phi_[id + shape.nx] +
                      eps_y_[id - shape.nx] * phi_[id - shape.nx] +
                      eps_z_[id] * phi_[id + shape.nx * shape.ny] +
                      eps_z_[id - shape.nx * shape.ny] *
                          phi_[id - shape.nx * shape.ny];
                  const double target =
                      (neighbor_sum + rhs * rho_[id]) / denom_[id];
                  const double updated = phi_[id] + omega * (target - phi_[id]);
                  const double delta = updated - phi_[id];
                  phi_[id] = updated;
                  const double abs_delta = std::abs(delta);
                  sumsq += delta * delta;
                  max_delta = std::max(max_delta, abs_delta);
                  ++count;
                }
              }
            }
            local_sumsq[worker] += sumsq;
            local_max[worker] = std::max(local_max[worker], max_delta);
            local_count[worker] += count;
          });
    }

    double sumsq = 0.0;
    double max_delta = 0.0;
    std::size_t count = 0;
    for (std::size_t worker = 0; worker < workers; ++worker) {
      sumsq += local_sumsq[worker];
      max_delta = std::max(max_delta, local_max[worker]);
      count += local_count[worker];
    }

    stats.iterations = iter;
    stats.rms = count > 0 ? std::sqrt(sumsq / static_cast<double>(count)) : 0.0;
    stats.max_dphi = max_delta;
    if (stats.rms <= args_.rms_tol || stats.max_dphi <= args_.dphi_tol) {
      stats.converged = true;
      break;
    }
  }
  return stats;
}

double Analyzer::compute_linear_residual(const GridShape &shape) const {
  if (shape.n == 0 || phi_.size() != shape.n || rho_.size() != shape.n)
    return 0.0;
  const double rhs = rhs_scale(shape);
  double sumsq = 0.0;
  std::size_t count = 0;
  for (std::size_t k = 1; k + 1 < shape.nz; ++k) {
    for (std::size_t j = 1; j + 1 < shape.ny; ++j) {
      for (std::size_t i = 1; i + 1 < shape.nx; ++i) {
        const std::size_t id = idx3(i, j, k, shape.nx, shape.ny);
        const double neighbor_sum =
            eps_x_[id] * phi_[id + 1] + eps_x_[id - 1] * phi_[id - 1] +
            eps_y_[id] * phi_[id + shape.nx] +
            eps_y_[id - shape.nx] * phi_[id - shape.nx] +
            eps_z_[id] * phi_[id + shape.nx * shape.ny] +
            eps_z_[id - shape.nx * shape.ny] *
                phi_[id - shape.nx * shape.ny];
        const double residual =
            rhs * rho_[id] - (denom_[id] * phi_[id] - neighbor_sum);
        sumsq += residual * residual;
        ++count;
      }
    }
  }
  return count > 0 ? std::sqrt(sumsq / static_cast<double>(count)) : 0.0;
}

double Analyzer::compute_nonlinear_residual(
    const GridShape &shape, const AlignedDoubleVector &source,
    double screen_scale) const {
  return compute_nonlinear_residual_from_phi(shape, phi_, source,
                                             screen_scale);
}

double Analyzer::compute_nonlinear_residual_from_phi(
    const GridShape &shape, const AlignedDoubleVector &field,
    const AlignedDoubleVector &source, double screen_scale) const {
  if (shape.n == 0 || field.size() != shape.n || rho_.size() != shape.n ||
      source.size() != shape.n)
    return 0.0;
  const double rhs = nonlinear_rhs_scale(shape, nonlinear_alpha());
  double sumsq = 0.0;
  std::size_t count = 0;
  for (std::size_t k = 1; k + 1 < shape.nz; ++k) {
    for (std::size_t j = 1; j + 1 < shape.ny; ++j) {
      for (std::size_t i = 1; i + 1 < shape.nx; ++i) {
        const std::size_t id = idx3(i, j, k, shape.nx, shape.ny);
        const double neighbor_sum =
            eps_x_[id] * field[id + 1] + eps_x_[id - 1] * field[id - 1] +
            eps_y_[id] * field[id + shape.nx] +
            eps_y_[id - shape.nx] * field[id - shape.nx] +
            eps_z_[id] * field[id + shape.nx * shape.ny] +
            eps_z_[id - shape.nx * shape.ny] *
                field[id - shape.nx * shape.ny];
        const double eps_sum = std::max(1.0e-12, denom_[id] - screen_[id]);
        const double residual =
            rhs * rho_[id] -
            (eps_sum * field[id] - neighbor_sum +
             (screen_[id] * screen_scale) * source[id]);
        sumsq += residual * residual;
        ++count;
      }
    }
  }
  return count > 0 ? std::sqrt(sumsq / static_cast<double>(count)) : 0.0;
}

double Analyzer::compute_nonlinear_residual_vector(
    const GridShape &shape, const AlignedDoubleVector &field,
    const AlignedDoubleVector &source, double screen_scale,
    AlignedDoubleVector &residual, TitanExecutor &executor) const {
  residual.assign(shape.n, 0.0);
  if (shape.n == 0 || field.size() != shape.n || rho_.size() != shape.n ||
      source.size() != shape.n)
    return 0.0;

  const std::size_t interior_z = shape.nz > 2 ? shape.nz - 2 : 0;
  const std::size_t workers = effective_workers(executor, interior_z);
  std::vector<double> local_sumsq(workers, 0.0);
  std::vector<std::size_t> local_count(workers, 0);
  const double rhs = nonlinear_rhs_scale(shape, nonlinear_alpha());
  const std::size_t xy = shape.nx * shape.ny;

  executor.parallel_for(
      interior_z, workers, 1,
      [&](std::size_t worker, std::size_t begin, std::size_t end) {
        double sumsq = 0.0;
        std::size_t count = 0;
        for (std::size_t kk = begin; kk < end; ++kk) {
          const std::size_t k = kk + 1;
          for (std::size_t j = 1; j + 1 < shape.ny; ++j) {
            for (std::size_t i = 1; i + 1 < shape.nx; ++i) {
              const std::size_t id = idx3(i, j, k, shape.nx, shape.ny);
              const double neighbor_sum =
                  eps_x_[id] * field[id + 1] +
                  eps_x_[id - 1] * field[id - 1] +
                  eps_y_[id] * field[id + shape.nx] +
                  eps_y_[id - shape.nx] * field[id - shape.nx] +
                  eps_z_[id] * field[id + xy] +
                  eps_z_[id - xy] * field[id - xy];
              const double eps_sum =
                  std::max(1.0e-12, denom_[id] - screen_[id]);
              const double r =
                  rhs * rho_[id] -
                  (eps_sum * field[id] - neighbor_sum +
                   (screen_[id] * screen_scale) * source[id]);
              residual[id] = r;
              sumsq += r * r;
              ++count;
            }
          }
        }
        local_sumsq[worker] = sumsq;
        local_count[worker] = count;
      });

  double sumsq = 0.0;
  std::size_t count = 0;
  for (std::size_t worker = 0; worker < workers; ++worker) {
    sumsq += local_sumsq[worker];
    count += local_count[worker];
  }
  return count > 0 ? std::sqrt(sumsq / static_cast<double>(count)) : 0.0;
}

Analyzer::NewtonCorrectionStats Analyzer::solve_newton_correction_sor(
    const GridShape &shape, const AlignedDoubleVector &residual,
    double residual_rms, double alpha, double clip, double screen_scale,
    AlignedDoubleVector &delta, int max_sweeps,
    TitanExecutor &executor) const {
  NewtonCorrectionStats stats;
  delta.assign(shape.n, 0.0);
  if (shape.n == 0 || residual.size() != shape.n)
    return stats;

  const int requested_sweeps =
      std::clamp(args_.nl_inner_iter, 1, std::max(1, args_.nl_inner_iter));
  const int sweeps =
      max_sweeps > 0 ? std::min(requested_sweeps, max_sweeps)
                     : requested_sweeps;
  const double omega =
      args_.omega == "auto"
          ? std::clamp(2.0 / (1.0 + std::sin(kPi / static_cast<double>(
                                             std::max({shape.nx, shape.ny,
                                                       shape.nz})))),
                       1.0, 1.95)
          : args_.omega_value;
  const double max_correction = 20.0 * std::max(1.0e-12, shape.h);
  const std::size_t interior_z = shape.nz > 2 ? shape.nz - 2 : 0;
  const std::size_t workers = effective_workers(executor, interior_z);
  const std::size_t xy = shape.nx * shape.ny;
  std::vector<double> local_sumsq(workers, 0.0);
  std::vector<double> local_max(workers, 0.0);
  std::vector<std::size_t> local_count(workers, 0);
  std::vector<double> local_resid_sumsq(workers, 0.0);
  std::vector<double> local_resid_max(workers, 0.0);
  std::vector<std::size_t> local_resid_count(workers, 0);
  std::vector<double> local_jac_min(workers,
                                    std::numeric_limits<double>::max());
  std::vector<double> local_jac_max(workers, 0.0);
  std::vector<std::size_t> local_screen_points(workers, 0);
  std::vector<std::size_t> local_active_points(workers, 0);

  auto compute_correction_residual = [&]() {
    std::fill(local_resid_sumsq.begin(), local_resid_sumsq.end(), 0.0);
    std::fill(local_resid_max.begin(), local_resid_max.end(), 0.0);
    std::fill(local_resid_count.begin(), local_resid_count.end(), 0);
    std::fill(local_jac_min.begin(), local_jac_min.end(),
              std::numeric_limits<double>::max());
    std::fill(local_jac_max.begin(), local_jac_max.end(), 0.0);
    std::fill(local_screen_points.begin(), local_screen_points.end(), 0);
    std::fill(local_active_points.begin(), local_active_points.end(), 0);

    executor.parallel_for(
        interior_z, workers, 1,
        [&](std::size_t worker, std::size_t begin, std::size_t end) {
          double sumsq = 0.0;
          double max_abs = 0.0;
          double jac_min = std::numeric_limits<double>::max();
          double jac_max = 0.0;
          std::size_t count = 0;
          std::size_t screen_points = 0;
          std::size_t active_points = 0;
          for (std::size_t kk = begin; kk < end; ++kk) {
            const std::size_t k = kk + 1;
            for (std::size_t j = 1; j + 1 < shape.ny; ++j) {
              for (std::size_t i = 1; i + 1 < shape.nx; ++i) {
                const std::size_t id = idx3(i, j, k, shape.nx, shape.ny);
                const double eps_sum =
                    std::max(1.0e-12, denom_[id] - screen_[id]);
                const double diag =
                    std::max(1.0e-12,
                             eps_sum + (screen_[id] * screen_scale) *
                                           nonlinear_source_derivative(
                                               phi_[id], alpha, clip));
                const double neighbor_delta_sum =
                    eps_x_[id] * delta[id + 1] +
                    eps_x_[id - 1] * delta[id - 1] +
                    eps_y_[id] * delta[id + shape.nx] +
                    eps_y_[id - shape.nx] * delta[id - shape.nx] +
                    eps_z_[id] * delta[id + xy] +
                    eps_z_[id - xy] * delta[id - xy];
                const double j_delta = diag * delta[id] - neighbor_delta_sum;
                const double r = residual[id] - j_delta;
                const double abs_r = std::abs(r);
                sumsq += r * r;
                max_abs = std::max(max_abs, abs_r);
                jac_min = std::min(jac_min, diag);
                jac_max = std::max(jac_max, diag);
                if (screen_[id] > 0.0)
                  ++screen_points;
                if (screen_[id] * screen_scale > 0.0)
                  ++active_points;
                ++count;
              }
            }
          }
          local_resid_sumsq[worker] = sumsq;
          local_resid_max[worker] = max_abs;
          local_resid_count[worker] = count;
          local_jac_min[worker] = jac_min;
          local_jac_max[worker] = jac_max;
          local_screen_points[worker] = screen_points;
          local_active_points[worker] = active_points;
        });

    double sumsq = 0.0;
    double max_abs = 0.0;
    double jac_min = std::numeric_limits<double>::max();
    double jac_max = 0.0;
    std::size_t count = 0;
    std::size_t screen_points = 0;
    std::size_t active_points = 0;
    for (std::size_t worker = 0; worker < workers; ++worker) {
      sumsq += local_resid_sumsq[worker];
      max_abs = std::max(max_abs, local_resid_max[worker]);
      count += local_resid_count[worker];
      jac_min = std::min(jac_min, local_jac_min[worker]);
      jac_max = std::max(jac_max, local_jac_max[worker]);
      screen_points += local_screen_points[worker];
      active_points += local_active_points[worker];
    }
    stats.residual_rms =
        count > 0 ? std::sqrt(sumsq / static_cast<double>(count)) : 0.0;
    stats.residual_max = max_abs;
    stats.jac_diag_min = count > 0 ? jac_min : 0.0;
    stats.jac_diag_max = count > 0 ? jac_max : 0.0;
    stats.screen_points = screen_points;
    stats.active_nonlinear_points = active_points;
  };

  compute_correction_residual();
  const double initial_correction_residual =
      std::isfinite(stats.residual_rms) ? stats.residual_rms : residual_rms;
  double best_residual = stats.residual_rms;
  int stale_iterations = 0;

  for (int iter = 1; iter <= sweeps; ++iter) {
    std::fill(local_sumsq.begin(), local_sumsq.end(), 0.0);
    std::fill(local_max.begin(), local_max.end(), 0.0);
    std::fill(local_count.begin(), local_count.end(), 0);

    for (int color = 0; color < 2; ++color) {
      executor.parallel_for(
          interior_z, workers, 1,
          [&](std::size_t worker, std::size_t begin, std::size_t end) {
            double sumsq = 0.0;
            double max_update = 0.0;
            std::size_t count = 0;
            for (std::size_t kk = begin; kk < end; ++kk) {
              const std::size_t k = kk + 1;
              for (std::size_t j = 1; j + 1 < shape.ny; ++j) {
                for (std::size_t i = 1; i + 1 < shape.nx; ++i) {
                  if (static_cast<int>((i + j + k) & 1ULL) != color)
                    continue;
                  const std::size_t id = idx3(i, j, k, shape.nx, shape.ny);
                  const double eps_sum =
                      std::max(1.0e-12, denom_[id] - screen_[id]);
                  const double diag =
                      std::max(1.0e-12,
                               eps_sum + (screen_[id] * screen_scale) *
                                             nonlinear_source_derivative(
                                                 phi_[id], alpha, clip));
                  const double neighbor_delta_sum =
                      eps_x_[id] * delta[id + 1] +
                      eps_x_[id - 1] * delta[id - 1] +
                      eps_y_[id] * delta[id + shape.nx] +
                      eps_y_[id - shape.nx] * delta[id - shape.nx] +
                      eps_z_[id] * delta[id + xy] +
                      eps_z_[id - xy] * delta[id - xy];
                  const double target = std::clamp(
                      (residual[id] + neighbor_delta_sum) / diag,
                      -max_correction, max_correction);
                  const double updated =
                      delta[id] + omega * (target - delta[id]);
                  const double d = updated - delta[id];
                  delta[id] = updated;
                  const double abs_d = std::abs(d);
                  sumsq += d * d;
                  max_update = std::max(max_update, abs_d);
                  ++count;
                }
              }
            }
            local_sumsq[worker] += sumsq;
            local_max[worker] = std::max(local_max[worker], max_update);
            local_count[worker] += count;
          });
    }

    double sumsq = 0.0;
    double max_update = 0.0;
    std::size_t count = 0;
    for (std::size_t worker = 0; worker < workers; ++worker) {
      sumsq += local_sumsq[worker];
      max_update = std::max(max_update, local_max[worker]);
      count += local_count[worker];
    }
    stats.sweeps.iterations = iter;
    stats.sweeps.rms =
        count > 0 ? std::sqrt(sumsq / static_cast<double>(count)) : 0.0;
    stats.sweeps.max_dphi = max_update;

    compute_correction_residual();
    const double correction_target =
        std::max(args_.nl_residual, 0.05 * initial_correction_residual);
    if (stats.residual_rms <= correction_target) {
      stats.sweeps.converged = true;
      stats.residual_decreased = true;
      break;
    }
    if (std::isfinite(stats.residual_rms) &&
        stats.residual_rms < best_residual * 0.999) {
      best_residual = stats.residual_rms;
      stale_iterations = 0;
    } else if (++stale_iterations >= std::min(4, args_.nl_stagnation)) {
      break;
    }
  }
  stats.residual_decreased =
      std::isfinite(stats.residual_rms) &&
      stats.residual_rms < initial_correction_residual;
  stats.inner_levels = 1;
  stats.inner_cycles = stats.sweeps.iterations;
  stats.inner_residual_before = initial_correction_residual;
  stats.inner_residual_after = stats.residual_rms;
  stats.inner_effective = "sor";
  return stats;
}

void Analyzer::prepare_newton_correction_hierarchy(
    const GridShape &shape, const AlignedDoubleVector &residual, double alpha,
    double clip, double screen_scale,
    std::vector<NewtonCorrectionLevel> &levels) const {
  const int level_count = effective_vcycle_levels(shape);
  levels.clear();
  levels.resize(static_cast<std::size_t>(level_count));

  auto allocate_level = [](NewtonCorrectionLevel &level,
                           const GridShape &level_shape,
                           double max_correction) {
    level.shape = level_shape;
    level.delta.assign(level_shape.n, 0.0);
    level.rhs.assign(level_shape.n, 0.0);
    level.eps_x.assign(level_shape.n, 0.0);
    level.eps_y.assign(level_shape.n, 0.0);
    level.eps_z.assign(level_shape.n, 0.0);
    level.diag.assign(level_shape.n, 1.0);
    level.screen.assign(level_shape.n, 0.0);
    level.base_phi.assign(level_shape.n, 0.0);
    level.max_correction = max_correction;
  };

  auto fill_diag = [&](NewtonCorrectionLevel &level) {
    const GridShape &s = level.shape;
    const std::size_t xy = s.nx * s.ny;
    std::fill(level.diag.begin(), level.diag.end(), 1.0);
    if (s.nx < 3 || s.ny < 3 || s.nz < 3)
      return;
    for (std::size_t k = 1; k + 1 < s.nz; ++k) {
      for (std::size_t j = 1; j + 1 < s.ny; ++j) {
        for (std::size_t i = 1; i + 1 < s.nx; ++i) {
          const std::size_t id = idx3(i, j, k, s.nx, s.ny);
          const double eps_sum =
              level.eps_x[id] + level.eps_x[id - 1] +
              level.eps_y[id] + level.eps_y[id - s.nx] +
              level.eps_z[id] + level.eps_z[id - xy];
          level.diag[id] =
              std::max(1.0e-12, eps_sum +
                                      (level.screen[id] * screen_scale) *
                                          nonlinear_source_derivative(
                                              level.base_phi[id], alpha,
                                              clip));
        }
      }
    }
  };

  const double correction_limit = 20.0 * std::max(1.0e-12, shape.h);
  allocate_level(levels.front(), shape, correction_limit);
  levels.front().rhs = residual;
  levels.front().eps_x = eps_x_;
  levels.front().eps_y = eps_y_;
  levels.front().eps_z = eps_z_;
  levels.front().screen = screen_;
  levels.front().base_phi = phi_;
  for (std::size_t k = 0; k < shape.nz; ++k) {
    for (std::size_t j = 0; j < shape.ny; ++j) {
      for (std::size_t i = 0; i < shape.nx; ++i) {
        if (i == 0 || j == 0 || k == 0 || i + 1 == shape.nx ||
            j + 1 == shape.ny || k + 1 == shape.nz) {
          levels.front().rhs[idx3(i, j, k, shape.nx, shape.ny)] = 0.0;
        }
      }
    }
  }
  fill_diag(levels.front());

  for (std::size_t level_index = 1; level_index < levels.size();
       ++level_index) {
    const NewtonCorrectionLevel &fine = levels[level_index - 1];
    NewtonCorrectionLevel &coarse = levels[level_index];
    allocate_level(coarse, coarsen_shape(fine.shape), correction_limit);
    for (std::size_t k = 0; k < coarse.shape.nz; ++k) {
      const double z = coarse.shape.oz + static_cast<double>(k) * coarse.shape.h;
      for (std::size_t j = 0; j < coarse.shape.ny; ++j) {
        const double y =
            coarse.shape.oy + static_cast<double>(j) * coarse.shape.h;
        for (std::size_t i = 0; i < coarse.shape.nx; ++i) {
          const double x =
              coarse.shape.ox + static_cast<double>(i) * coarse.shape.h;
          const std::size_t id =
              idx3(i, j, k, coarse.shape.nx, coarse.shape.ny);
          double value = 0.0;
          if (interpolate_field_value(fine.shape, fine.eps_x, x, y, z, value))
            coarse.eps_x[id] = value;
          if (interpolate_field_value(fine.shape, fine.eps_y, x, y, z, value))
            coarse.eps_y[id] = value;
          if (interpolate_field_value(fine.shape, fine.eps_z, x, y, z, value))
            coarse.eps_z[id] = value;
          if (interpolate_field_value(fine.shape, fine.screen, x, y, z, value))
            coarse.screen[id] = 4.0 * value;
          if (interpolate_field_value(fine.shape, fine.base_phi, x, y, z,
                                      value))
            coarse.base_phi[id] = value;
        }
      }
    }
    fill_diag(coarse);
  }
}

Analyzer::SolveStats Analyzer::smooth_newton_correction_level(
    NewtonCorrectionLevel &level, int iterations,
    bool stop_on_convergence) const {
  SolveStats stats;
  if (iterations <= 0)
    return stats;
  const GridShape &shape = level.shape;
  const double omega =
      args_.omega == "auto"
          ? std::clamp(2.0 / (1.0 + std::sin(kPi / static_cast<double>(
                                             std::max({shape.nx, shape.ny,
                                                       shape.nz})))),
                       1.0, 1.95)
          : args_.omega_value;
  const double max_correction =
      std::max(1.0e-12, level.max_correction);
  for (int iter = 1; iter <= iterations; ++iter) {
    double sumsq = 0.0;
    double max_delta = 0.0;
    std::size_t count = 0;
    for (int color = 0; color < 2; ++color) {
      for (std::size_t k = 1; k + 1 < shape.nz; ++k) {
        for (std::size_t j = 1; j + 1 < shape.ny; ++j) {
          for (std::size_t i = 1; i + 1 < shape.nx; ++i) {
            if (static_cast<int>((i + j + k) & 1ULL) != color)
              continue;
            const std::size_t id = idx3(i, j, k, shape.nx, shape.ny);
            const double neighbor_sum =
                level.eps_x[id] * level.delta[id + 1] +
                level.eps_x[id - 1] * level.delta[id - 1] +
                level.eps_y[id] * level.delta[id + shape.nx] +
                level.eps_y[id - shape.nx] * level.delta[id - shape.nx] +
                level.eps_z[id] * level.delta[id + shape.nx * shape.ny] +
                level.eps_z[id - shape.nx * shape.ny] *
                    level.delta[id - shape.nx * shape.ny];
            const double target = std::clamp(
                (level.rhs[id] + neighbor_sum) /
                    std::max(1.0e-12, level.diag[id]),
                -max_correction, max_correction);
            const double updated =
                std::clamp(level.delta[id] + omega * (target - level.delta[id]),
                           -max_correction, max_correction);
            const double d = updated - level.delta[id];
            level.delta[id] = updated;
            sumsq += d * d;
            max_delta = std::max(max_delta, std::abs(d));
            ++count;
          }
        }
      }
    }
    stats.iterations = iter;
    stats.rms = count > 0 ? std::sqrt(sumsq / static_cast<double>(count)) : 0.0;
    stats.max_dphi = max_delta;
    if (stop_on_convergence &&
        (stats.rms <= args_.nl_residual || stats.max_dphi <= args_.nl_dphi)) {
      stats.converged = true;
      break;
    }
  }
  return stats;
}

double Analyzer::compute_newton_correction_level_residual(
    const NewtonCorrectionLevel &level, AlignedDoubleVector &residual) const {
  const GridShape &shape = level.shape;
  residual.assign(shape.n, 0.0);
  double sumsq = 0.0;
  std::size_t count = 0;
  for (std::size_t k = 1; k + 1 < shape.nz; ++k) {
    for (std::size_t j = 1; j + 1 < shape.ny; ++j) {
      for (std::size_t i = 1; i + 1 < shape.nx; ++i) {
        const std::size_t id = idx3(i, j, k, shape.nx, shape.ny);
        const double neighbor_sum =
            level.eps_x[id] * level.delta[id + 1] +
            level.eps_x[id - 1] * level.delta[id - 1] +
            level.eps_y[id] * level.delta[id + shape.nx] +
            level.eps_y[id - shape.nx] * level.delta[id - shape.nx] +
            level.eps_z[id] * level.delta[id + shape.nx * shape.ny] +
            level.eps_z[id - shape.nx * shape.ny] *
                level.delta[id - shape.nx * shape.ny];
        const double j_delta = level.diag[id] * level.delta[id] - neighbor_sum;
        const double r = level.rhs[id] - j_delta;
        residual[id] = r;
        sumsq += r * r;
        ++count;
      }
    }
  }
  return count > 0 ? std::sqrt(sumsq / static_cast<double>(count)) : 0.0;
}

void Analyzer::restrict_newton_correction_residual(
    const NewtonCorrectionLevel &fine,
    const AlignedDoubleVector &fine_residual,
    NewtonCorrectionLevel &coarse) const {
  std::fill(coarse.rhs.begin(), coarse.rhs.end(), 0.0);
  std::fill(coarse.delta.begin(), coarse.delta.end(), 0.0);
  const int weights[3] = {1, 2, 1};
  for (std::size_t kc = 1; kc + 1 < coarse.shape.nz; ++kc) {
    for (std::size_t jc = 1; jc + 1 < coarse.shape.ny; ++jc) {
      for (std::size_t ic = 1; ic + 1 < coarse.shape.nx; ++ic) {
        const std::size_t fi = std::min(fine.shape.nx - 2, 2 * ic);
        const std::size_t fj = std::min(fine.shape.ny - 2, 2 * jc);
        const std::size_t fk = std::min(fine.shape.nz - 2, 2 * kc);
        double sum = 0.0;
        int wsum = 0;
        for (int dz = -1; dz <= 1; ++dz) {
          const std::size_t iz = static_cast<std::size_t>(
              std::clamp<int>(static_cast<int>(fk) + dz, 1,
                              static_cast<int>(fine.shape.nz - 2)));
          for (int dy = -1; dy <= 1; ++dy) {
            const std::size_t iy = static_cast<std::size_t>(
                std::clamp<int>(static_cast<int>(fj) + dy, 1,
                                static_cast<int>(fine.shape.ny - 2)));
            for (int dx = -1; dx <= 1; ++dx) {
              const std::size_t ix = static_cast<std::size_t>(
                  std::clamp<int>(static_cast<int>(fi) + dx, 1,
                                  static_cast<int>(fine.shape.nx - 2)));
              const int w = weights[dx + 1] * weights[dy + 1] *
                            weights[dz + 1];
              sum += static_cast<double>(w) *
                     fine_residual[idx3(ix, iy, iz, fine.shape.nx,
                                        fine.shape.ny)];
              wsum += w;
            }
          }
        }
        coarse.rhs[idx3(ic, jc, kc, coarse.shape.nx, coarse.shape.ny)] =
            wsum > 0 ? sum / static_cast<double>(wsum) : 0.0;
      }
    }
  }
}

void Analyzer::prolongate_add_newton_correction(
    const NewtonCorrectionLevel &coarse, NewtonCorrectionLevel &fine) const {
  const double max_correction =
      std::max(1.0e-12, fine.max_correction);
  const double correction_damping = 0.25;
  for (std::size_t k = 1; k + 1 < fine.shape.nz; ++k) {
    const double z = fine.shape.oz + static_cast<double>(k) * fine.shape.h;
    for (std::size_t j = 1; j + 1 < fine.shape.ny; ++j) {
      const double y = fine.shape.oy + static_cast<double>(j) * fine.shape.h;
      for (std::size_t i = 1; i + 1 < fine.shape.nx; ++i) {
        const double x = fine.shape.ox + static_cast<double>(i) * fine.shape.h;
        double correction = 0.0;
        if (interpolate_field_value(coarse.shape, coarse.delta, x, y, z,
                                    correction)) {
          const std::size_t id = idx3(i, j, k, fine.shape.nx, fine.shape.ny);
          fine.delta[id] = std::clamp(fine.delta[id] +
                                          correction_damping * correction,
                                      -max_correction, max_correction);
        }
      }
    }
  }
}

void Analyzer::run_newton_correction_vcycle(
    std::vector<NewtonCorrectionLevel> &levels,
    std::size_t level_index) const {
  NewtonCorrectionLevel &level = levels[level_index];
  if (level_index + 1 >= levels.size()) {
    smooth_newton_correction_level(level, args_.mg_coarse_iter, true);
    return;
  }

  smooth_newton_correction_level(level, args_.mg_pre, false);
  AlignedDoubleVector residual;
  compute_newton_correction_level_residual(level, residual);
  NewtonCorrectionLevel &coarse = levels[level_index + 1];
  restrict_newton_correction_residual(level, residual, coarse);
  run_newton_correction_vcycle(levels, level_index + 1);
  prolongate_add_newton_correction(coarse, level);
  smooth_newton_correction_level(level, args_.mg_post, false);
}

Analyzer::NewtonCorrectionStats Analyzer::solve_newton_correction_vcycle(
    const GridShape &shape, const AlignedDoubleVector &residual,
    double residual_rms, double alpha, double clip, double screen_scale,
    AlignedDoubleVector &delta, int max_cycles) const {
  NewtonCorrectionStats stats;
  delta.assign(shape.n, 0.0);
  if (shape.n == 0 || residual.size() != shape.n)
    return stats;

  std::vector<NewtonCorrectionLevel> levels;
  prepare_newton_correction_hierarchy(shape, residual, alpha, clip,
                                      screen_scale, levels);
  stats.inner_effective = "vcycle";
  stats.inner_levels = static_cast<int>(levels.size());

  AlignedDoubleVector correction_residual;
  stats.residual_rms =
      compute_newton_correction_level_residual(levels.front(),
                                               correction_residual);
  stats.residual_max = 0.0;
  stats.jac_diag_min = std::numeric_limits<double>::max();
  stats.jac_diag_max = 0.0;
  stats.screen_points = 0;
  stats.active_nonlinear_points = 0;
  for (std::size_t k = 1; k + 1 < shape.nz; ++k) {
    for (std::size_t j = 1; j + 1 < shape.ny; ++j) {
      for (std::size_t i = 1; i + 1 < shape.nx; ++i) {
        const std::size_t id = idx3(i, j, k, shape.nx, shape.ny);
        stats.residual_max =
            std::max(stats.residual_max, std::abs(correction_residual[id]));
        stats.jac_diag_min = std::min(stats.jac_diag_min, levels.front().diag[id]);
        stats.jac_diag_max = std::max(stats.jac_diag_max, levels.front().diag[id]);
        if (levels.front().screen[id] > 0.0)
          ++stats.screen_points;
        if (levels.front().screen[id] * screen_scale > 0.0)
          ++stats.active_nonlinear_points;
      }
    }
  }
  if (stats.jac_diag_min == std::numeric_limits<double>::max())
    stats.jac_diag_min = 0.0;

  const double initial_correction_residual =
      std::isfinite(stats.residual_rms) ? stats.residual_rms : residual_rms;
  stats.inner_residual_before = initial_correction_residual;
  double best_residual = stats.residual_rms;
  int stale_cycles = 0;
  const int cycles =
      max_cycles > 0 ? max_cycles : std::max(1, args_.nl_inner_iter);

  for (int cycle = 1; cycle <= cycles; ++cycle) {
    const AlignedDoubleVector before = levels.front().delta;
    run_newton_correction_vcycle(levels, 0);
    double sumsq = 0.0;
    double max_update = 0.0;
    std::size_t count = 0;
    for (std::size_t k = 1; k + 1 < shape.nz; ++k) {
      for (std::size_t j = 1; j + 1 < shape.ny; ++j) {
        for (std::size_t i = 1; i + 1 < shape.nx; ++i) {
          const std::size_t id = idx3(i, j, k, shape.nx, shape.ny);
          const double d = levels.front().delta[id] - before[id];
          sumsq += d * d;
          max_update = std::max(max_update, std::abs(d));
          ++count;
        }
      }
    }
    stats.sweeps.iterations = cycle;
    stats.sweeps.rms =
        count > 0 ? std::sqrt(sumsq / static_cast<double>(count)) : 0.0;
    stats.sweeps.max_dphi = max_update;
    stats.residual_rms =
        compute_newton_correction_level_residual(levels.front(),
                                                 correction_residual);
    stats.residual_max = 0.0;
    for (double r : correction_residual)
      stats.residual_max = std::max(stats.residual_max, std::abs(r));
    const double correction_target =
        std::max(args_.nl_residual, 0.05 * initial_correction_residual);
    if (stats.residual_rms <= correction_target) {
      stats.sweeps.converged = true;
      break;
    }
    if (std::isfinite(stats.residual_rms) &&
        stats.residual_rms < best_residual * 0.999) {
      best_residual = stats.residual_rms;
      stale_cycles = 0;
    } else if (++stale_cycles >= std::min(4, args_.nl_stagnation)) {
      break;
    }
  }

  delta = levels.front().delta;
  stats.inner_cycles = stats.sweeps.iterations;
  stats.inner_residual_after = stats.residual_rms;
  stats.residual_decreased =
      std::isfinite(stats.residual_rms) &&
      stats.residual_rms < initial_correction_residual;
  return stats;
}

Analyzer::NewtonCorrectionStats Analyzer::solve_newton_correction(
    const GridShape &shape, const AlignedDoubleVector &residual,
    double residual_rms, double alpha, double clip, double screen_scale,
    AlignedDoubleVector &delta, int max_iterations,
    const std::string &inner_solver, TitanExecutor &executor) const {
  if (inner_solver == "vcycle") {
    NewtonCorrectionStats vcycle =
        solve_newton_correction_vcycle(shape, residual, residual_rms, alpha,
                                       clip, screen_scale, delta,
                                       max_iterations);
    if (vcycle.residual_decreased)
      return vcycle;
    NewtonCorrectionStats fallback =
        solve_newton_correction_sor(shape, residual, residual_rms, alpha, clip,
                                    screen_scale, delta, max_iterations,
                                    executor);
    fallback.inner_effective = "vcycle_sor_fallback";
    fallback.inner_levels = vcycle.inner_levels;
    fallback.inner_cycles = vcycle.inner_cycles;
    fallback.inner_residual_before = vcycle.inner_residual_before;
    fallback.inner_residual_after = fallback.residual_rms;
    return fallback;
  }
  return solve_newton_correction_sor(shape, residual, residual_rms, alpha,
                                     clip, screen_scale, delta,
                                     max_iterations, executor);
}

std::size_t Analyzer::update_nonlinear_source(
    const GridShape &shape, AlignedDoubleVector &source, double alpha,
    double clip) const {
  if (source.size() != shape.n)
    source.assign(shape.n, 0.0);
  std::size_t clipped = 0;
  for (std::size_t id = 0; id < shape.n; ++id) {
    if (screen_[id] <= 0.0) {
      source[id] = 0.0;
      continue;
    }
    const double raw_arg = nonlinear_argument_value(phi_[id], alpha);
    if (raw_arg > clip) {
      ++clipped;
    } else if (raw_arg < -clip) {
      ++clipped;
    }
    source[id] = nonlinear_source_value(phi_[id], alpha, clip);
  }
  return clipped;
}

double Analyzer::grid_reaction_energy_kcal(
    const GridShape &shape, const AlignedDoubleVector &solvent_phi,
    const AlignedDoubleVector &vacuum_phi) const {
  if (shape.n == 0 || solvent_phi.size() != shape.n ||
      vacuum_phi.size() != shape.n || rho_.size() != shape.n)
    return 0.0;
  double sum = 0.0;
  for (std::size_t id = 0; id < shape.n; ++id)
    sum += rho_[id] * (solvent_phi[id] - vacuum_phi[id]);
  return 0.5 * kCoulombKcalPerMol * sum;
}

double Analyzer::grid_phase_energy_kcal(
    const GridShape &shape, const AlignedDoubleVector &field) const {
  if (shape.n == 0 || field.size() != shape.n || rho_.size() != shape.n)
    return 0.0;
  double sum = 0.0;
  for (std::size_t id = 0; id < shape.n; ++id)
    sum += rho_[id] * field[id];
  return 0.5 * kCoulombKcalPerMol * sum;
}

double Analyzer::atom_phase_energy_kcal(
    const std::vector<Atom> &atoms, const std::vector<double> &atom_phi) const {
  if (atoms.size() != atom_phi.size())
    return 0.0;
  double sum = 0.0;
  for (std::size_t i = 0; i < atoms.size(); ++i)
    sum += atoms[i].q * atom_phi[i];
  return 0.5 * kCoulombKcalPerMol * sum;
}

double Analyzer::analytic_coulomb_energy_kcal(
    const std::vector<Atom> &atoms, double eps) const {
  const double dielectric = std::max(1.0e-12, eps);
  double energy = 0.0;
  for (std::size_t i = 0; i < atoms.size(); ++i) {
    for (std::size_t j = i + 1; j < atoms.size(); ++j) {
      const double dx = atoms[i].x - atoms[j].x;
      const double dy = atoms[i].y - atoms[j].y;
      const double dz = atoms[i].z - atoms[j].z;
      const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
      if (r > 1.0e-12)
        energy += kCoulombKcalPerMol * atoms[i].q * atoms[j].q /
                  (dielectric * r);
    }
  }
  return energy;
}

double Analyzer::effective_nl_relax() const noexcept {
  if (args_.nl_relax == "auto")
    return 0.5;
  return std::clamp(args_.nl_relax_value, 0.05, 1.0);
}

double Analyzer::nonlinear_alpha() const noexcept {
  return kCoulombKcalPerMol /
         (kKBKcalPerMolK * std::max(1.0, args_.temp));
}

double Analyzer::rhs_scale(const GridShape &shape) const noexcept {
  return 4.0 * kPi / std::max(1.0e-12, shape.h);
}

bool Analyzer::nonlinear_dimensionless_mode() const noexcept {
  return args_.test_nl_source == "delphi-dimensionless";
}

double Analyzer::nonlinear_rhs_scale(const GridShape &shape,
                                     double alpha) const noexcept {
  const double base = rhs_scale(shape);
  return nonlinear_dimensionless_mode() ? alpha * base : base;
}

double Analyzer::nonlinear_argument_value(double field_value,
                                          double alpha) const noexcept {
  return nonlinear_dimensionless_mode() ? field_value : alpha * field_value;
}

double Analyzer::nonlinear_source_value(double phi, double alpha,
                                        double clip) const noexcept {
  if (alpha <= 0.0)
    return 0.0;
  const double arg =
      std::clamp(nonlinear_argument_value(phi, alpha), -clip, clip);
  if (args_.test_nl_source == "delphi-dimensionless")
    return std::sinh(arg);
  return std::sinh(arg) / alpha;
}

double Analyzer::nonlinear_source_derivative(double phi, double alpha,
                                             double clip) const noexcept {
  if (alpha <= 0.0)
    return 0.0;
  const double arg =
      std::clamp(nonlinear_argument_value(phi, alpha), -clip, clip);
  return std::cosh(arg);
}

void Analyzer::scale_nonlinear_field_to_solve_units(AlignedDoubleVector &field,
                                                    double alpha) const {
  if (!nonlinear_dimensionless_mode() || alpha <= 0.0)
    return;
  for (double &value : field)
    value *= alpha;
}

void Analyzer::scale_nonlinear_field_to_physical_units(
    AlignedDoubleVector &field, double alpha) const {
  if (!nonlinear_dimensionless_mode() || alpha <= 0.0)
    return;
  const double inv_alpha = 1.0 / alpha;
  for (double &value : field)
    value *= inv_alpha;
}

Analyzer::ParityAudit Analyzer::compute_parity_audit(
    const GridShape &shape, const AlignedDoubleVector *source,
    const AlignedDoubleVector *residual) const {
  ParityAudit audit;
  if (shape.n == 0 || denom_.size() != shape.n || screen_.size() != shape.n)
    return audit;

  double eps_min = std::numeric_limits<double>::max();
  double eps_max = 0.0;
  double eps_sum_total = 0.0;
  double screen_sumsq = 0.0;
  double source_sumsq = 0.0;
  std::size_t count = 0;
  for (std::size_t id = 0; id < shape.n; ++id) {
    const double eps_sum = std::max(1.0e-12, denom_[id] - screen_[id]);
    eps_min = std::min(eps_min, eps_sum);
    eps_max = std::max(eps_max, eps_sum);
    eps_sum_total += eps_sum;
    if (screen_[id] > 0.0) {
      ++audit.screen_points;
      audit.screen_sum += screen_[id];
      screen_sumsq += screen_[id] * screen_[id];
      audit.screen_max = std::max(audit.screen_max, std::abs(screen_[id]));
    }
    if (source != nullptr && source->size() == shape.n) {
      const double value = (*source)[id];
      audit.source_sum += value;
      source_sumsq += value * value;
      audit.source_max = std::max(audit.source_max, std::abs(value));
    }
    if (residual != nullptr && residual->size() == shape.n)
      audit.residual_max_abs =
          std::max(audit.residual_max_abs, std::abs((*residual)[id]));
    ++count;
  }
  audit.eps_sum_min = count > 0 ? eps_min : 0.0;
  audit.eps_sum_max = count > 0 ? eps_max : 0.0;
  audit.eps_sum_mean =
      count > 0 ? eps_sum_total / static_cast<double>(count) : 0.0;
  audit.screen_rms =
      audit.screen_points > 0
          ? std::sqrt(screen_sumsq / static_cast<double>(audit.screen_points))
          : 0.0;
  audit.source_rms =
      count > 0 ? std::sqrt(source_sumsq / static_cast<double>(count)) : 0.0;
  if (eps_x_.size() == shape.n && eps_y_.size() == shape.n &&
      eps_z_.size() == shape.n) {
    const std::size_t xy = shape.nx * shape.ny;
    double midpoint_sum_total = 0.0;
    double midpoint_interior_total = 0.0;
    std::size_t midpoint_count = 0;
    std::size_t midpoint_interior_count = 0;
    double midpoint_min = std::numeric_limits<double>::max();
    double midpoint_max = std::numeric_limits<double>::lowest();
    for (std::size_t k = 0; k < shape.nz; ++k) {
      for (std::size_t j = 0; j < shape.ny; ++j) {
        for (std::size_t i = 0; i < shape.nx; ++i) {
          const std::size_t id = idx3(i, j, k, shape.nx, shape.ny);
          double midpoint_sum = eps_x_[id] + eps_y_[id] + eps_z_[id];
          if (i > 0)
            midpoint_sum += eps_x_[id - 1];
          if (j > 0)
            midpoint_sum += eps_y_[id - shape.nx];
          if (k > 0)
            midpoint_sum += eps_z_[id - xy];
          if (!std::isfinite(midpoint_sum))
            continue;
          midpoint_min = std::min(midpoint_min, midpoint_sum);
          midpoint_max = std::max(midpoint_max, midpoint_sum);
          midpoint_sum_total += midpoint_sum;
          ++midpoint_count;
          if (i > 0 && i + 1 < shape.nx && j > 0 && j + 1 < shape.ny &&
              k > 0 && k + 1 < shape.nz) {
            midpoint_interior_total += midpoint_sum;
            ++midpoint_interior_count;
          }
        }
      }
    }
    audit.midpoint_eps_sum_min =
        midpoint_count > 0 ? midpoint_min : 0.0;
    audit.midpoint_eps_sum_max =
        midpoint_count > 0 ? midpoint_max : 0.0;
    audit.midpoint_eps_sum_mean =
        midpoint_count > 0
            ? midpoint_sum_total / static_cast<double>(midpoint_count)
            : 0.0;
    audit.midpoint_eps_sum_interior_mean =
        midpoint_interior_count > 0
            ? midpoint_interior_total /
                  static_cast<double>(midpoint_interior_count)
            : 0.0;
  }
  return audit;
}

Analyzer::ForensicAudit
Analyzer::compute_forensic_audit(const GridShape &shape) const {
  ForensicAudit audit;
  if (shape.n == 0 || inside_.size() != shape.n || rho_.size() != shape.n)
    return audit;

  const double voxel_volume = shape.h * shape.h * shape.h;
  std::size_t inside_count = 0;
  std::size_t face_count = 0;
  double eps_face_sum = 0.0;
  double eps_face_min = std::numeric_limits<double>::max();
  double eps_face_max = 0.0;

  auto add_face_eps = [&](double value) {
    if (!std::isfinite(value))
      return;
    eps_face_min = std::min(eps_face_min, value);
    eps_face_max = std::max(eps_face_max, value);
    eps_face_sum += value;
    ++face_count;
  };
  double eps_axis_sum[3] = {0.0, 0.0, 0.0};
  std::size_t eps_axis_count[3] = {0, 0, 0};
  std::size_t internal_axis_count[3] = {0, 0, 0};
  const double internal_threshold = 0.5 * (args_.int_diel + args_.ext_diel);
  auto add_axis_face_eps = [&](int axis, double value) {
    if (!std::isfinite(value) || axis < 0 || axis > 2)
      return;
    eps_axis_sum[axis] += value;
    ++eps_axis_count[axis];
    if (value <= internal_threshold)
      ++internal_axis_count[axis];
  };

  const std::size_t xy = shape.nx * shape.ny;
  for (std::size_t k = 0; k < shape.nz; ++k) {
    for (std::size_t j = 0; j < shape.ny; ++j) {
      for (std::size_t i = 0; i < shape.nx; ++i) {
        const std::size_t id = idx3(i, j, k, shape.nx, shape.ny);
        if (inside_[id])
          ++inside_count;
        if (i + 1 < shape.nx) {
          if (inside_[id] != inside_[id + 1])
            ++audit.face_cross_x;
          if (eps_x_.size() == shape.n) {
            add_face_eps(eps_x_[id]);
            add_axis_face_eps(0, eps_x_[id]);
          }
        }
        if (j + 1 < shape.ny) {
          if (inside_[id] != inside_[id + shape.nx])
            ++audit.face_cross_y;
          if (eps_y_.size() == shape.n) {
            add_face_eps(eps_y_[id]);
            add_axis_face_eps(1, eps_y_[id]);
          }
        }
        if (k + 1 < shape.nz) {
          if (inside_[id] != inside_[id + xy])
            ++audit.face_cross_z;
          if (eps_z_.size() == shape.n) {
            add_face_eps(eps_z_[id]);
            add_axis_face_eps(2, eps_z_[id]);
          }
        }

        const double rho = rho_[id];
        const double x = shape.ox + static_cast<double>(i) * shape.h;
        const double y = shape.oy + static_cast<double>(j) * shape.h;
        const double z = shape.oz + static_cast<double>(k) * shape.h;
        audit.rho_sum += rho;
        audit.rho_l1 += std::abs(rho);
        audit.rho_l2 += rho * rho;
        audit.rho_mx += rho * x;
        audit.rho_my += rho * y;
        audit.rho_mz += rho * z;
      }
    }
  }
  audit.vdw_volume_estimate_A3 = static_cast<double>(inside_count) * voxel_volume;
  audit.eps_face_min = face_count > 0 ? eps_face_min : 0.0;
  audit.eps_face_max = face_count > 0 ? eps_face_max : 0.0;
  audit.eps_face_mean =
      face_count > 0 ? eps_face_sum / static_cast<double>(face_count) : 0.0;
  audit.eps_face_count_x = eps_axis_count[0];
  audit.eps_face_count_y = eps_axis_count[1];
  audit.eps_face_count_z = eps_axis_count[2];
  audit.internal_face_count_x = internal_axis_count[0];
  audit.internal_face_count_y = internal_axis_count[1];
  audit.internal_face_count_z = internal_axis_count[2];
  audit.eps_face_mean_x =
      eps_axis_count[0] > 0
          ? eps_axis_sum[0] / static_cast<double>(eps_axis_count[0])
          : 0.0;
  audit.eps_face_mean_y =
      eps_axis_count[1] > 0
          ? eps_axis_sum[1] / static_cast<double>(eps_axis_count[1])
          : 0.0;
  audit.eps_face_mean_z =
      eps_axis_count[2] > 0
          ? eps_axis_sum[2] / static_cast<double>(eps_axis_count[2])
          : 0.0;
  audit.rho_l2 = std::sqrt(audit.rho_l2);

  double denom_sum = 0.0;
  std::size_t denom_count = 0;
  double denom_min = std::numeric_limits<double>::max();
  double denom_max = 0.0;
  if (denom_.size() == shape.n) {
    for (double value : denom_) {
      if (!std::isfinite(value))
        continue;
      denom_min = std::min(denom_min, value);
      denom_max = std::max(denom_max, value);
      denom_sum += value;
      ++denom_count;
    }
  }
  audit.denom_min = denom_count > 0 ? denom_min : 0.0;
  audit.denom_max = denom_count > 0 ? denom_max : 0.0;
  audit.denom_mean =
      denom_count > 0 ? denom_sum / static_cast<double>(denom_count) : 0.0;

  if (screen_.size() == shape.n) {
    for (double value : screen_) {
      if (!std::isfinite(value))
        continue;
      audit.screen_sum += value;
      audit.screen_max = std::max(audit.screen_max, std::abs(value));
    }
  }
  return audit;
}

Analyzer::IscmAudit Analyzer::compute_iscm_audit(
    const GridShape &shape, const std::vector<Atom> &atoms,
    const AlignedDoubleVector &solvent_phi,
    const AlignedDoubleVector &solvent_screen,
    const std::vector<unsigned char> &solvent_inside,
    const AlignedDoubleVector &solvent_eps_x,
    const AlignedDoubleVector &solvent_eps_y,
    const AlignedDoubleVector &solvent_eps_z,
    double grid_vacuum_phase_kcal_mol) const {
  IscmAudit audit;
  if (args_.test_pqr.empty()) {
    audit.status = "requires_test_pqr";
    return audit;
  }
  if (args_.dielectric_model != "two-dielectric" || uses_gaussian_model()) {
    audit.status = "unsupported_dielectric_model";
    return audit;
  }
  if (args_.model != "linear") {
    audit.status = "linear_only";
    return audit;
  }
  if (shape.n == 0 || shape.h <= 0.0 || atoms.empty() ||
      solvent_phi.size() != shape.n || solvent_inside.size() != shape.n ||
      rho_.size() != shape.n) {
    audit.status = "missing_field";
    return audit;
  }

  const double kt_kj = kKBKJPerMolK * std::max(1.0, args_.temp);
  const double kt_kcal = kKBKcalPerMolK * std::max(1.0, args_.temp);
  const double epkt = kDelphiEpk / std::max(1.0, args_.temp);
  const double potential_scale = nonlinear_alpha();
  const double grid_scale = 1.0 / std::max(1.0e-12, shape.h);
  const double induced_scale =
      kDelphiReactionFieldEnergyFactor /
      (2.0 * std::max(1.0e-12, grid_scale) * std::max(1.0e-12, epkt));
  const std::size_t plane = shape.nx * shape.ny;
  const bool delphi_vdw_boundary =
      args_.surface == "delphi-vdw" && solvent_eps_x.size() == shape.n &&
      solvent_eps_y.size() == shape.n && solvent_eps_z.size() == shape.n;
  auto face_is_solute = [&](double value) {
    return std::abs(value - args_.int_diel) <=
           std::abs(value - args_.ext_diel);
  };

  std::vector<IscmSample> deterministic_samples;
  std::vector<IscmSample> top_samples;
  auto add_sample = [&](IscmSample sample) {
    if (deterministic_samples.size() < 16) {
      sample.kind = "boundary_sample";
      sample.rank = deterministic_samples.size() + 1;
      deterministic_samples.push_back(sample);
    }
    sample.kind = "top_energy";
    top_samples.push_back(sample);
    std::sort(top_samples.begin(), top_samples.end(),
              [](const IscmSample &a, const IscmSample &b) {
                return std::abs(a.energy_contribution_kt) >
                       std::abs(b.energy_contribution_kt);
              });
    if (top_samples.size() > 32)
      top_samples.resize(32);
  };

  for (std::size_t k = 1; k + 1 < shape.nz; ++k) {
    for (std::size_t j = 1; j + 1 < shape.ny; ++j) {
      for (std::size_t i = 1; i + 1 < shape.nx; ++i) {
        const std::size_t id = idx3(i, j, k, shape.nx, shape.ny);
        const unsigned char center_inside = solvent_inside[id];
        bool dielectric_boundary = false;
        if (delphi_vdw_boundary) {
          const double faces[6] = {solvent_eps_x[id],
                                   solvent_eps_x[id - 1],
                                   solvent_eps_y[id],
                                   solvent_eps_y[id - shape.nx],
                                   solvent_eps_z[id],
                                   solvent_eps_z[id - plane]};
          int solute_faces = 0;
          int solvent_faces = 0;
          for (double face_eps : faces) {
            if (!std::isfinite(face_eps))
              continue;
            if (face_is_solute(face_eps))
              ++solute_faces;
            else
              ++solvent_faces;
          }
          dielectric_boundary = solute_faces > 0 && solvent_faces > 0;
          if (dielectric_boundary) {
            ++audit.delphi_face_crossings;
            audit.delphi_eps_face_matches +=
                static_cast<std::size_t>(solute_faces);
          }
        } else {
          dielectric_boundary =
              solvent_inside[id - 1] != center_inside ||
              solvent_inside[id + 1] != center_inside ||
              solvent_inside[id - shape.nx] != center_inside ||
              solvent_inside[id + shape.nx] != center_inside ||
              solvent_inside[id - plane] != center_inside ||
              solvent_inside[id + plane] != center_inside;
        }
        if (!dielectric_boundary)
          continue;

        const double phi_center = potential_scale * solvent_phi[id];
        const double phi_neigh =
            potential_scale *
            (solvent_phi[id - 1] + solvent_phi[id + 1] +
             solvent_phi[id - shape.nx] + solvent_phi[id + shape.nx] +
             solvent_phi[id - plane] + solvent_phi[id + plane]) /
            6.0;
        double charge_adjustment = 0.0;
        if (std::abs(rho_[id]) > 1.0e-14) {
          ++audit.charged_boundary_points;
          const double solvent_flag = center_inside ? 0.0 : 1.0;
          const double salt_diag =
              solvent_screen.size() == shape.n ? solvent_screen[id] * solvent_flag
                                                : 0.0;
          const double denom =
              6.0 * args_.int_diel / std::max(1.0e-12, epkt) + salt_diag;
          if (std::abs(denom) > 1.0e-30)
            charge_adjustment = rho_[id] * 4.0 * kPi * grid_scale / denom;
        }

        const double phi_effective = phi_center - phi_neigh - charge_adjustment;
        const double induced_charge = phi_effective * induced_scale;
        const double x = shape.ox + static_cast<double>(i) * shape.h;
        const double y = shape.oy + static_cast<double>(j) * shape.h;
        const double z = shape.oz + static_cast<double>(k) * shape.h;
        double atom_potential_sum = 0.0;
        for (const Atom &atom : atoms) {
          const double dx = x - atom.x;
          const double dy = y - atom.y;
          const double dz = z - atom.z;
          const double r2 = dx * dx + dy * dy + dz * dz;
          if (r2 > 1.0e-24)
            atom_potential_sum += atom.q / std::sqrt(r2);
        }
        const double contribution_kt =
            0.5 * epkt * induced_charge * atom_potential_sum;
        audit.water_phase_kt += contribution_kt;
        audit.induced_charge_sum += induced_charge;
        audit.induced_charge_l1 += std::abs(induced_charge);
        audit.induced_charge_l2 += induced_charge * induced_charge;
        audit.induced_charge_max =
            std::max(audit.induced_charge_max, std::abs(induced_charge));
        ++audit.boundary_points;

        IscmSample sample;
        sample.i = i;
        sample.j = j;
        sample.k = k;
        sample.x = x;
        sample.y = y;
        sample.z = z;
        sample.phi_center = phi_center;
        sample.phi_neighbor_average = phi_neigh;
        sample.charge_adjustment = charge_adjustment;
        sample.induced_charge = induced_charge;
        sample.atom_potential_sum = atom_potential_sum;
        sample.energy_contribution_kt = contribution_kt;
        add_sample(sample);
      }
    }
  }

  audit.induced_charge_l2 = std::sqrt(audit.induced_charge_l2);
  audit.water_phase_kj_mol = audit.water_phase_kt * kt_kj;
  audit.water_phase_kcal_mol = audit.water_phase_kt * kt_kcal;
  audit.vacuum_phase_kcal_mol = grid_vacuum_phase_kcal_mol;
  audit.vacuum_phase_kj_mol = grid_vacuum_phase_kcal_mol * kKcalToKJ;
  audit.vacuum_phase_kt =
      kt_kcal > 0.0 ? grid_vacuum_phase_kcal_mol / kt_kcal : 0.0;
  // DelPhi's two-dielectric VDW ISCM path finalizes corrected RF from the
  // water induced-surface-charge phase alone. Keep the vacuum grid phase in
  // the audit so convention differences remain inspectable.
  audit.corrected_kt = audit.water_phase_kt;
  audit.corrected_kj_mol = audit.corrected_kt * kt_kj;
  audit.corrected_kcal_mol = audit.corrected_kt * kt_kcal;
  audit.status = audit.boundary_points > 0 ? "ok" : "no_boundary_points";

  audit.samples = deterministic_samples;
  for (std::size_t idx = 0; idx < top_samples.size(); ++idx) {
    top_samples[idx].rank = idx + 1;
    audit.samples.push_back(top_samples[idx]);
  }
  return audit;
}

void Analyzer::copy_iscm_audit_to_solve(const IscmAudit &audit,
                                        ReactionSolve &solve) const {
  solve.iscm_boundary_points = audit.boundary_points;
  solve.iscm_charged_boundary_points = audit.charged_boundary_points;
  solve.iscm_induced_charge_sum = audit.induced_charge_sum;
  solve.iscm_induced_charge_l1 = audit.induced_charge_l1;
  solve.iscm_induced_charge_l2 = audit.induced_charge_l2;
  solve.iscm_induced_charge_max = audit.induced_charge_max;
  solve.iscm_water_phase_kcal_mol = audit.water_phase_kcal_mol;
  solve.iscm_water_phase_kj_mol = audit.water_phase_kj_mol;
  solve.iscm_water_phase_kt = audit.water_phase_kt;
  solve.iscm_vacuum_phase_kcal_mol = audit.vacuum_phase_kcal_mol;
  solve.iscm_vacuum_phase_kj_mol = audit.vacuum_phase_kj_mol;
  solve.iscm_vacuum_phase_kt = audit.vacuum_phase_kt;
  solve.iscm_corrected_kcal_mol = audit.corrected_kcal_mol;
  solve.iscm_corrected_kj_mol = audit.corrected_kj_mol;
  solve.iscm_corrected_kt = audit.corrected_kt;
  solve.iscm_status = audit.status;
  solve.delphi_face_crossings = audit.delphi_face_crossings;
  solve.delphi_eps_face_matches = audit.delphi_eps_face_matches;
  solve.delphi_boundary_grid_points = audit.boundary_points;
  solve.delphi_charged_boundary_points = audit.charged_boundary_points;
  solve.delphi_iscm_status = audit.status;
}

void Analyzer::update_newton_derivative_audit(
    const GridShape &shape, const AlignedDoubleVector &source,
    const AlignedDoubleVector &delta, double alpha, double clip,
    double screen_scale, ParityAudit &audit) const {
  if (!args_.test_newton_derivative || shape.n == 0 ||
      phi_.size() != shape.n || source.size() != shape.n ||
      delta.size() != shape.n || rho_.size() != shape.n)
    return;

  double max_delta = 0.0;
  for (double value : delta)
    max_delta = std::max(max_delta, std::abs(value));
  if (max_delta <= 1.0e-14)
    return;

  const std::size_t xy = shape.nx * shape.ny;
  const double rhs = nonlinear_rhs_scale(shape, alpha);
  const double eps_step = 1.0e-6 / std::max(1.0, max_delta);
  const std::size_t interior =
      shape.nx > 2 && shape.ny > 2 && shape.nz > 2
          ? (shape.nx - 2) * (shape.ny - 2) * (shape.nz - 2)
          : 0;
  const std::size_t stride = std::max<std::size_t>(1, interior / 4096);
  std::size_t ordinal = 0;
  for (std::size_t k = 1; k + 1 < shape.nz; ++k) {
    for (std::size_t j = 1; j + 1 < shape.ny; ++j) {
      for (std::size_t i = 1; i + 1 < shape.nx; ++i, ++ordinal) {
        if (ordinal % stride != 0)
          continue;
        const std::size_t id = idx3(i, j, k, shape.nx, shape.ny);
        const double eps_sum = std::max(1.0e-12, denom_[id] - screen_[id]);
        const double diag =
            std::max(1.0e-12, eps_sum + (screen_[id] * screen_scale) *
                                           nonlinear_source_derivative(
                                               phi_[id], alpha, clip));
        const double neighbor_delta_sum =
            eps_x_[id] * delta[id + 1] + eps_x_[id - 1] * delta[id - 1] +
            eps_y_[id] * delta[id + shape.nx] +
            eps_y_[id - shape.nx] * delta[id - shape.nx] +
            eps_z_[id] * delta[id + xy] + eps_z_[id - xy] * delta[id - xy];
        const double analytic = diag * delta[id] - neighbor_delta_sum;

        const double neighbor_sum =
            eps_x_[id] * phi_[id + 1] + eps_x_[id - 1] * phi_[id - 1] +
            eps_y_[id] * phi_[id + shape.nx] +
            eps_y_[id - shape.nx] * phi_[id - shape.nx] +
            eps_z_[id] * phi_[id + xy] + eps_z_[id - xy] * phi_[id - xy];
        const double f0 =
            eps_sum * phi_[id] - neighbor_sum +
            (screen_[id] * screen_scale) * source[id] - rhs * rho_[id];

        const auto trial_phi = [&](std::size_t probe) {
          return phi_[probe] + eps_step * delta[probe];
        };
        const double trial_neighbor_sum =
            eps_x_[id] * trial_phi(id + 1) +
            eps_x_[id - 1] * trial_phi(id - 1) +
            eps_y_[id] * trial_phi(id + shape.nx) +
            eps_y_[id - shape.nx] * trial_phi(id - shape.nx) +
            eps_z_[id] * trial_phi(id + xy) +
            eps_z_[id - xy] * trial_phi(id - xy);
        const double trial_center = trial_phi(id);
        const double trial_source =
            nonlinear_source_value(trial_center, alpha, clip);
        const double f1 =
            eps_sum * trial_center - trial_neighbor_sum +
            (screen_[id] * screen_scale) * trial_source - rhs * rho_[id];
        const double finite_difference = (f1 - f0) / eps_step;
        const double abs_error = std::abs(analytic - finite_difference);
        const double rel_error =
            abs_error / std::max(1.0, std::abs(finite_difference));
        audit.derivative_max_abs =
            std::max(audit.derivative_max_abs, abs_error);
        audit.derivative_max_rel =
            std::max(audit.derivative_max_rel, rel_error);
        ++audit.derivative_samples;
      }
    }
  }
}

Analyzer::SolveStats Analyzer::solve_nonlinear_sor_block(
    const GridShape &shape, const AlignedDoubleVector &source, int sweeps,
    double relax, double screen_scale, TitanExecutor &executor,
    std::size_t &clipped_points) {
  SolveStats stats;
  const double omega =
      args_.omega == "auto"
          ? std::clamp(2.0 / (1.0 + std::sin(kPi / static_cast<double>(
                                             std::max({shape.nx, shape.ny,
                                                       shape.nz})))),
                       1.0, 1.95)
          : args_.omega_value;
  const double effective_omega = std::clamp(omega * relax, 0.05, 1.95);
  const std::size_t interior_z = shape.nz > 2 ? shape.nz - 2 : 0;
  const std::size_t workers = effective_workers(executor, interior_z);
  const double rhs = nonlinear_rhs_scale(shape, nonlinear_alpha());

  std::vector<double> local_sumsq(workers, 0.0);
  std::vector<double> local_max(workers, 0.0);
  std::vector<std::size_t> local_count(workers, 0);
  clipped_points = 0;

  for (int iter = 1; iter <= sweeps; ++iter) {
    std::fill(local_sumsq.begin(), local_sumsq.end(), 0.0);
    std::fill(local_max.begin(), local_max.end(), 0.0);
    std::fill(local_count.begin(), local_count.end(), 0);

    for (int color = 0; color < 2; ++color) {
      executor.parallel_for(
          interior_z, workers, 1,
          [&](std::size_t worker, std::size_t begin, std::size_t end) {
            double sumsq = 0.0;
            double max_delta = 0.0;
            std::size_t count = 0;
            for (std::size_t kk = begin; kk < end; ++kk) {
              const std::size_t k = kk + 1;
              for (std::size_t j = 1; j + 1 < shape.ny; ++j) {
                for (std::size_t i = 1; i + 1 < shape.nx; ++i) {
                  if (static_cast<int>((i + j + k) & 1ULL) != color)
                    continue;
                  const std::size_t id = idx3(i, j, k, shape.nx, shape.ny);
                  const double neighbor_sum =
                      eps_x_[id] * phi_[id + 1] +
                      eps_x_[id - 1] * phi_[id - 1] +
                      eps_y_[id] * phi_[id + shape.nx] +
                      eps_y_[id - shape.nx] * phi_[id - shape.nx] +
                      eps_z_[id] * phi_[id + shape.nx * shape.ny] +
                      eps_z_[id - shape.nx * shape.ny] *
                          phi_[id - shape.nx * shape.ny];
                  const double eps_sum =
                      std::max(1.0e-12, denom_[id] - screen_[id]);
                  const double target =
                      (neighbor_sum + rhs * rho_[id] -
                       (screen_[id] * screen_scale) * source[id]) /
                      eps_sum;
                  const double updated =
                      phi_[id] + effective_omega * (target - phi_[id]);
                  const double delta = updated - phi_[id];
                  phi_[id] = updated;
                  const double abs_delta = std::abs(delta);
                  sumsq += delta * delta;
                  max_delta = std::max(max_delta, abs_delta);
                  ++count;
                }
              }
            }
            local_sumsq[worker] += sumsq;
            local_max[worker] = std::max(local_max[worker], max_delta);
            local_count[worker] += count;
          });
    }

    double sumsq = 0.0;
    double max_delta = 0.0;
    std::size_t count = 0;
    for (std::size_t worker = 0; worker < workers; ++worker) {
      sumsq += local_sumsq[worker];
      max_delta = std::max(max_delta, local_max[worker]);
      count += local_count[worker];
    }
    stats.iterations = iter;
    stats.rms = count > 0 ? std::sqrt(sumsq / static_cast<double>(count)) : 0.0;
    stats.max_dphi = max_delta;
  }
  stats.converged = stats.max_dphi <= args_.nl_dphi;
  return stats;
}

int Analyzer::effective_focus_levels() const noexcept {
  if (args_.focus != "auto" && args_.focus != "delphi")
    return 1;
  return std::clamp(args_.focus_levels <= 1 ? 2 : args_.focus_levels, 2, 4);
}

std::string Analyzer::effective_solver_mode() const {
  if (args_.multigrid == "vcycle")
    return "vcycle";
  if (args_.multigrid == "focus" || args_.multigrid == "auto" ||
      args_.focus == "auto")
    return "nested_focus";
  if (args_.focus == "delphi")
    return "delphi_focus";
  return "sor";
}

bool Analyzer::interpolate_field_value(const GridShape &shape,
                                       const AlignedDoubleVector &field,
                                       double x, double y, double z,
                                       double &value) const {
  if (shape.nx < 2 || shape.ny < 2 || shape.nz < 2 ||
      field.size() != shape.n)
    return false;
  const double gx = (x - shape.ox) / shape.h;
  const double gy = (y - shape.oy) / shape.h;
  const double gz = (z - shape.oz) / shape.h;
  if (gx < 0.0 || gy < 0.0 || gz < 0.0 ||
      gx > static_cast<double>(shape.nx - 1) ||
      gy > static_cast<double>(shape.ny - 1) ||
      gz > static_cast<double>(shape.nz - 1))
    return false;

  const std::size_t i0 = std::min<std::size_t>(
      shape.nx - 2, static_cast<std::size_t>(std::floor(gx)));
  const std::size_t j0 = std::min<std::size_t>(
      shape.ny - 2, static_cast<std::size_t>(std::floor(gy)));
  const std::size_t k0 = std::min<std::size_t>(
      shape.nz - 2, static_cast<std::size_t>(std::floor(gz)));
  const double fx = clamp_fraction(gx - static_cast<double>(i0));
  const double fy = clamp_fraction(gy - static_cast<double>(j0));
  const double fz = clamp_fraction(gz - static_cast<double>(k0));

  value = 0.0;
  for (std::size_t dz = 0; dz <= 1; ++dz) {
    const double wz = dz ? fz : (1.0 - fz);
    for (std::size_t dy = 0; dy <= 1; ++dy) {
      const double wy = dy ? fy : (1.0 - fy);
      for (std::size_t dx = 0; dx <= 1; ++dx) {
        const double wx = dx ? fx : (1.0 - fx);
        value += wx * wy * wz *
                 field[idx3(i0 + dx, j0 + dy, k0 + dz, shape.nx, shape.ny)];
      }
    }
  }
  return true;
}

void Analyzer::seed_from_focus_source(const GridShape &shape,
                                      const GridShape &source_shape,
                                      const AlignedDoubleVector &source_phi,
                                      TitanExecutor &executor) {
  seed_from_focus_source(shape, source_shape, source_phi, executor, false);
}

void Analyzer::seed_from_focus_source(const GridShape &shape,
                                      const GridShape &source_shape,
                                      const AlignedDoubleVector &source_phi,
                                      TitanExecutor &executor,
                                      bool boundary_only) {
  if (source_phi.empty())
    return;
  const std::size_t workers = effective_workers(executor, shape.nz);
  executor.parallel_for(shape.nz, workers, 1,
                        [&](std::size_t, std::size_t z_begin,
                            std::size_t z_end) {
    for (std::size_t k = z_begin; k < z_end; ++k) {
      const double z = shape.oz + static_cast<double>(k) * shape.h;
      for (std::size_t j = 0; j < shape.ny; ++j) {
        const double y = shape.oy + static_cast<double>(j) * shape.h;
        for (std::size_t i = 0; i < shape.nx; ++i) {
          if (boundary_only && !is_boundary(i, j, k, shape))
            continue;
          const double x = shape.ox + static_cast<double>(i) * shape.h;
          double seeded = 0.0;
          if (interpolate_field_value(source_shape, source_phi, x, y, z,
                                      seeded)) {
            phi_[idx3(i, j, k, shape.nx, shape.ny)] = seeded;
          }
        }
      }
    }
  });
}

Analyzer::FieldSolve Analyzer::solve_focused_field(
    const std::vector<Atom> &atoms, double eps_in, double eps_out, double kappa,
    TitanExecutor &executor, const FrameResult *map_result) {
  FieldSolve field;
  const int levels = effective_focus_levels();
  field.levels = levels;
  field.solver_mode =
      levels > 1 ? (args_.focus == "delphi" ? "delphi_focus"
                                             : "nested_focus")
                 : "sor";
  field.mg_levels = levels;

  GridShape previous_shape;
  AlignedDoubleVector previous_phi;
  bool has_previous = false;
  for (int level = 0; level < levels; ++level) {
    const int remaining = levels - level - 1;
    double scale = args_.scale;
    if (remaining > 0) {
      const double divisor = std::pow(2.0, static_cast<double>(remaining));
      scale = std::max(0.25, args_.scale / divisor);
      scale = std::min(scale, args_.scale);
    }
    double margin = args_.margin;
    if (remaining > 0 && levels > 1) {
      margin += args_.focus_margin *
                (static_cast<double>(remaining) /
                 static_cast<double>(std::max(1, levels - 1)));
    }
    if (level + 1 == levels) {
      scale = args_.scale;
      margin = args_.margin;
    }
    if (level == 0) {
      field.initial_scale = scale;
      field.initial_margin = margin;
    }

    const GridShape shape = build_shape(atoms, scale, margin);
    const std::size_t peak_bytes =
        shape.memory_bytes + previous_phi.size() * sizeof(double);
    field.peak_memory_bytes = std::max(field.peak_memory_bytes, peak_bytes);
    ensure_workspace(shape);
    classify_surface(shape, atoms, executor);
    apply_oracle_surface_map(shape);
    classify_ion_accessibility(shape, atoms, executor);
    const ChargeDiagnostics charge_diag =
        args_.test_charge_grid == "oracle"
            ? assign_oracle_charge_grid(shape, atoms)
            : assign_charges(shape, atoms);
    build_coefficients(shape, atoms, eps_in, eps_out, kappa, executor);
    if (level + 1 == levels && map_result != nullptr)
      maybe_write_maps(shape, *map_result, false);
    BoundaryDiagnostics boundary_diag =
        set_boundary(shape, atoms, eps_out, kappa);
    if (args_.test_boundary_map == "oracle" && oracle_replay_applies(shape))
      boundary_diag = apply_oracle_boundary_map(shape);
    if (has_previous)
      seed_from_focus_source(shape, previous_shape, previous_phi, executor,
                             args_.focus == "delphi");

    SolveStats stats = solve_sor(shape, executor);
    field.total_iterations += stats.iterations;
    if (level + 1 == levels) {
      if (map_result != nullptr)
        maybe_write_maps(shape, *map_result, true);
      field.shape = shape;
      field.stats = stats;
      if (has_previous) {
        field.focus_parent_grid = std::to_string(previous_shape.nx) + "x" +
                                  std::to_string(previous_shape.ny) + "x" +
                                  std::to_string(previous_shape.nz);
      }
      field.focus_child_grid = std::to_string(shape.nx) + "x" +
                               std::to_string(shape.ny) + "x" +
                               std::to_string(shape.nz);
      field.rhs_scale = rhs_scale(shape);
      field.charge_assigned = charge_diag.assigned_charge;
      field.charge_moment_error = charge_diag.moment_error;
      field.gaussian_density_sum = last_gaussian_density_sum_;
      field.gaussian_cutoff_points = last_gaussian_cutoff_points_;
      field.gaussian_eps_min = last_gaussian_eps_min_;
      field.gaussian_eps_max = last_gaussian_eps_max_;
      field.gaussian_eps_mean = last_gaussian_eps_mean_;
      field.vdwms_enabled = last_vdwms_audit_.enabled;
      field.vdwms_boundary_initial = last_vdwms_audit_.boundary_initial;
      field.vdwms_boundary_final = last_vdwms_audit_.boundary_final;
      field.vdwms_added = last_vdwms_audit_.added;
      field.vdwms_removed = last_vdwms_audit_.removed;
      field.vdwms_exposed_points = last_vdwms_audit_.exposed_points;
      field.vdwms_midpoint_rms = last_vdwms_audit_.midpoint_rms;
      field.gaussian_charge_norm_error =
          std::abs(charge_diag.assigned_charge -
                   std::accumulate(atoms.begin(), atoms.end(), 0.0,
                                   [](double sum, const Atom &atom) {
                                     return sum + atom.q;
                                   }));
      field.boundary_rms = boundary_diag.rms;
      field.boundary_max = boundary_diag.max_abs;
      field.linear_residual = compute_linear_residual(shape);
      field.parity = compute_parity_audit(shape, nullptr, nullptr);
      field.mg_pre = 0;
      field.mg_post = 0;
      field.mg_coarse_iter = 0;
      field.solute_points = static_cast<std::size_t>(
          std::count(inside_.begin(), inside_.end(), 1U));
      field.solvent_points = shape.n - field.solute_points;
      field.ion_accessible_points = static_cast<std::size_t>(
          std::count(ion_accessible_.begin(), ion_accessible_.end(), 1U));
      field.stern_excluded_points =
          args_.ion_access == "stern" ? shape.n - field.ion_accessible_points
                                      : 0;
      field.atom_phi = interpolate_atoms(shape, atoms);
      const bool keep_field_grids =
          parity_map_dump_enabled() ||
          args_.test_energy_convention == "delphi-iscm" ||
          args_.test_energy_convention == "delphi-corrected" ||
          barnase_forensic_enabled();
      if (args_.charge_assignment == "delphi" ||
          args_.charge_assignment == "gaussian" || keep_field_grids)
        field.phi_grid = phi_;
      if (keep_field_grids) {
        field.rho_grid = rho_;
        field.screen_grid = screen_;
        field.inside_grid = inside_;
        field.eps_x_grid = eps_x_;
        field.eps_y_grid = eps_y_;
        field.eps_z_grid = eps_z_;
        field.eps_sum_grid.assign(shape.n, 0.0);
        for (std::size_t id = 0; id < shape.n; ++id) {
          field.eps_sum_grid[id] =
              denom_.size() == shape.n && screen_.size() == shape.n
                  ? std::max(1.0e-12, denom_[id] - screen_[id])
                  : 0.0;
        }
      }
    } else {
      previous_shape = shape;
      previous_phi = phi_;
      has_previous = true;
    }
  }

  return field;
}

Analyzer::FieldSolve Analyzer::solve_linear_field(
    const std::vector<Atom> &atoms, double eps_in, double eps_out, double kappa,
    TitanExecutor &executor, const FrameResult *map_result) {
  if (effective_solver_mode() == "vcycle")
    return solve_vcycle_field(atoms, eps_in, eps_out, kappa, executor,
                              map_result);
  return solve_focused_field(atoms, eps_in, eps_out, kappa, executor,
                             map_result);
}

Analyzer::FieldSolve Analyzer::solve_newton_nonlinear_field(
    const std::vector<Atom> &atoms, double eps_in, double eps_out, double kappa,
    TitanExecutor &executor, const FrameResult *map_result) {
  FieldSolve field =
      solve_linear_field(atoms, eps_in, eps_out, kappa, executor, nullptr);
  field.linear_seed_solver = field.solver_mode;
  field.pbe_model_effective = "nonlinear";
  field.nl_solver = "newton";
  field.nl_inner_solver =
      args_.nl_inner == "auto"
          ? (field.shape.n > 250000 ? std::string("vcycle") : std::string("sor"))
          : args_.nl_inner;
  field.nl_mg_mode = args_.nl_multigrid;
  field.nl_relax = args_.nl_relax == "auto" ? 1.0 : effective_nl_relax();

  if (kappa <= 0.0 || args_.salt <= 0.0) {
    field.pbe_model_effective = "linear";
    field.nl_status = "linearized_no_salt";
    if (map_result != nullptr) {
      maybe_write_maps(field.shape, *map_result, false);
      maybe_write_maps(field.shape, *map_result, true);
    }
    return field;
  }

  if (map_result != nullptr)
    maybe_write_maps(field.shape, *map_result, false);

  const double alpha = nonlinear_alpha();
  const double clip = std::max(1.0e-12, args_.nl_clip);
  const bool dimensionless_mode = nonlinear_dimensionless_mode();
  scale_nonlinear_field_to_solve_units(phi_, alpha);
  AlignedDoubleVector source(field.shape.n, 0.0);
  AlignedDoubleVector residual_vec(field.shape.n, 0.0);
  AlignedDoubleVector delta(field.shape.n, 0.0);
  SolveStats last_stats;
  int total_iterations = 0;
  int blocks = 0;
  int backtracks = 0;
  int continuation_steps = 0;
  int inner_iterations = 0;
  int rejected_steps = 0;
  std::size_t max_clipped = 0;
  double max_step = 0.0;
  double last_step_norm = 0.0;
  double last_residual_before = 0.0;
  double last_residual_after = 0.0;
  double last_correction_residual = 0.0;
  double last_correction_max_residual = 0.0;
  double last_accepted_lambda = 0.0;
  double last_jac_diag_min = 0.0;
  double last_jac_diag_max = 0.0;
  std::size_t last_screen_points = 0;
  std::size_t last_active_nonlinear_points = 0;
  int last_inner_levels = 1;
  int last_inner_cycles = 0;
  double last_inner_residual_before = 0.0;
  double last_inner_residual_after = 0.0;
  std::string last_inner_effective = field.nl_inner_solver;
  double residual = std::numeric_limits<double>::infinity();
  double final_screen_scale = 1.0;
  double trust_psi = 0.5;
  bool finite = true;
  bool converged = false;
  std::string stop_reason = "work_budget";
  std::string line_search_status =
      args_.nl_multigrid == "fas" ? "fas_sor_correction" : "not_started";

  std::vector<double> screen_scales;
  if (args_.nl_continuation == "auto" && args_.salt > 0.5)
    screen_scales = {0.1, 0.25, 0.5, 0.75, 1.0};
  else if (args_.nl_continuation == "auto")
    screen_scales = {0.25, 0.5, 0.75, 1.0};
  else
    screen_scales = {1.0};

  for (double screen_scale : screen_scales) {
    final_screen_scale = screen_scale;
    ++continuation_steps;
    std::size_t clipped = update_nonlinear_source(field.shape, source, alpha,
                                                  clip);
    max_clipped = std::max(max_clipped, clipped);
    residual = compute_nonlinear_residual_vector(
        field.shape, phi_, source, screen_scale, residual_vec, executor);
    last_residual_after = residual;
    const double stage_initial_residual =
        std::isfinite(residual) ? residual : std::numeric_limits<double>::max();
    const double stage_residual_target =
        screen_scale < 1.0
            ? std::max({args_.nl_residual, 1.0e-4,
                        0.20 * stage_initial_residual})
            : args_.nl_residual;
    bool step_done = false;
    int stagnant_steps = 0;
    double best_stage_residual = residual;
    while (total_iterations < args_.nl_max_iter && !step_done) {
      const AlignedDoubleVector before = phi_;
      last_residual_before = residual;
      if (screen_scale >= 1.0 && residual <= args_.nl_residual) {
        converged = true;
        last_stats.converged = true;
        stop_reason = "residual_tol";
        step_done = true;
        break;
      }

      const NewtonCorrectionStats correction_stats =
          solve_newton_correction(field.shape, residual_vec, residual, alpha,
                                  clip, screen_scale, delta,
                                  args_.nl_inner_iter, field.nl_inner_solver,
                                  executor);
      last_stats = correction_stats.sweeps;
      const int used_inner = std::max(1, correction_stats.sweeps.iterations);
      inner_iterations += used_inner;
      ++total_iterations;
      ++blocks;
      last_correction_residual = correction_stats.residual_rms;
      last_correction_max_residual = correction_stats.residual_max;
      last_jac_diag_min = correction_stats.jac_diag_min;
      last_jac_diag_max = correction_stats.jac_diag_max;
      last_screen_points = correction_stats.screen_points;
      last_active_nonlinear_points =
          correction_stats.active_nonlinear_points;
      last_inner_levels = correction_stats.inner_levels;
      last_inner_cycles = correction_stats.inner_cycles;
      last_inner_residual_before = correction_stats.inner_residual_before;
      last_inner_residual_after = correction_stats.inner_residual_after;
      last_inner_effective = correction_stats.inner_effective;

      bool accepted = false;
      double accepted_residual = std::numeric_limits<double>::infinity();
      double accepted_step_norm = 0.0;
      double accepted_max_step = 0.0;
      std::size_t accepted_clipped = clipped;
      double trial_lambda = field.nl_relax;
      while (trial_lambda >= 1.0e-4) {
        phi_ = before;
        double sumsq = 0.0;
        double local_max_step = 0.0;
        std::size_t count = 0;
        const double max_field_step =
            args_.nl_line_search == "auto"
                ? trust_psi /
                      (dimensionless_mode ? 1.0 : std::max(1.0e-12, alpha))
                : std::numeric_limits<double>::infinity();
        for (std::size_t k = 1; k + 1 < field.shape.nz; ++k) {
          for (std::size_t j = 1; j + 1 < field.shape.ny; ++j) {
            for (std::size_t i = 1; i + 1 < field.shape.nx; ++i) {
              const std::size_t id =
                  idx3(i, j, k, field.shape.nx, field.shape.ny);
              const double raw_step = trial_lambda * delta[id];
              const double step =
                  std::clamp(raw_step, -max_field_step, max_field_step);
              phi_[id] = before[id] + step;
              sumsq += step * step;
              local_max_step = std::max(local_max_step, std::abs(step));
              ++count;
            }
          }
        }
        finite = std::all_of(phi_.begin(), phi_.end(),
                             [](double v) { return std::isfinite(v); });
        if (finite) {
          const std::size_t trial_clipped =
              update_nonlinear_source(field.shape, source, alpha, clip);
          const double trial_residual = compute_nonlinear_residual_vector(
              field.shape, phi_, source, screen_scale, residual_vec, executor);
          const bool residual_reduced = std::isfinite(trial_residual) &&
                                        trial_residual <=
                                            residual * (1.0 - 1.0e-4 *
                                                                 trial_lambda);
          if (args_.nl_line_search != "auto" || residual_reduced ||
              trial_residual <= args_.nl_residual) {
            accepted = true;
            accepted_residual = trial_residual;
            accepted_step_norm =
                count > 0 ? std::sqrt(sumsq / static_cast<double>(count)) : 0.0;
            accepted_max_step = local_max_step;
            accepted_clipped = trial_clipped;
            line_search_status =
                trial_lambda < field.nl_relax ? "accepted_backtracked"
                                             : "accepted_full";
            last_accepted_lambda = trial_lambda;
            break;
          }
        }
        if (args_.nl_line_search != "auto")
          break;
        trial_lambda *= 0.5;
        ++backtracks;
      }

      if (!finite) {
        trust_psi = std::max(0.05, trust_psi * 0.5);
        line_search_status = "nonfinite";
        stop_reason = "nonfinite";
        break;
      }
      if (!accepted) {
        trust_psi = std::max(0.05, trust_psi * 0.5);
        phi_ = before;
        update_nonlinear_source(field.shape, source, alpha, clip);
        std::size_t fallback_clipped = 0;
        SolveStats fallback_stats = solve_nonlinear_sor_block(
            field.shape, source, std::min(args_.nl_block, 4), 0.25,
            screen_scale, executor, fallback_clipped);
        const std::size_t fallback_source_clipped =
            update_nonlinear_source(field.shape, source, alpha, clip);
        const double fallback_residual = compute_nonlinear_residual_vector(
            field.shape, phi_, source, screen_scale, residual_vec, executor);
        if (std::isfinite(fallback_residual) && fallback_residual < residual) {
          accepted = true;
          residual = fallback_residual;
          last_residual_after = residual;
          last_step_norm = fallback_stats.rms;
          max_step = std::max(max_step, fallback_stats.max_dphi);
          max_clipped =
              std::max(max_clipped,
                       std::max(fallback_clipped, fallback_source_clipped));
          line_search_status = "fallback_picard";
          last_stats = fallback_stats;
          last_stats.rms = residual;
          ++total_iterations;
          ++blocks;
          if (residual < best_stage_residual * 0.999) {
            best_stage_residual = residual;
            stagnant_steps = 0;
          } else if (++stagnant_steps >= args_.nl_stagnation) {
            stop_reason = "stagnation";
            step_done = true;
          }
          if (screen_scale >= 1.0 && residual <= args_.nl_residual) {
            converged = true;
            stop_reason = "residual_tol";
            step_done = true;
          } else if (screen_scale < 1.0 &&
                     (residual <= stage_residual_target ||
                      fallback_stats.max_dphi <= args_.nl_dphi)) {
            stop_reason = "stage_target";
            step_done = true;
          }
          continue;
        }
        phi_ = before;
        update_nonlinear_source(field.shape, source, alpha, clip);
        compute_nonlinear_residual_vector(field.shape, phi_, source,
                                          screen_scale, residual_vec,
                                          executor);
        if (screen_scale < 1.0 &&
            (residual <= stage_residual_target ||
             residual < stage_initial_residual)) {
          line_search_status = "stage_target_line_search_floor";
          stop_reason = "stage_target";
          step_done = true;
        } else {
          line_search_status = "line_search_floor";
          stop_reason = "line_search_floor";
        }
        ++rejected_steps;
        ++backtracks;
        break;
      }

      if (args_.nl_line_search == "auto" && last_accepted_lambda >= field.nl_relax &&
          accepted_residual < 0.95 * last_residual_before) {
        field.nl_relax = std::min(1.0, field.nl_relax * 1.25);
        trust_psi = std::min(2.0, trust_psi * 1.25);
      }
      residual = accepted_residual;
      last_residual_after = residual;
      last_step_norm = accepted_step_norm;
      max_step = std::max(max_step, accepted_max_step);
      max_clipped = std::max(max_clipped, accepted_clipped);
      last_stats.rms = residual;
      last_stats.max_dphi = accepted_max_step;
      last_stats.converged = false;
      if (accepted_residual < best_stage_residual * 0.999) {
        best_stage_residual = accepted_residual;
        stagnant_steps = 0;
      } else if (++stagnant_steps >= args_.nl_stagnation) {
        stop_reason = "stagnation";
        step_done = true;
      }
      if (screen_scale >= 1.0 && residual <= args_.nl_residual) {
        last_stats.converged = true;
        converged = true;
        stop_reason = "residual_tol";
        step_done = true;
      } else if (screen_scale < 1.0 &&
                 (residual <= stage_residual_target ||
                  accepted_max_step <= args_.nl_dphi)) {
        stop_reason = "stage_target";
        step_done = true;
      }
    }
    if (!finite || total_iterations >= args_.nl_max_iter ||
        stop_reason == "line_search_floor" ||
        (screen_scale >= 1.0 && (converged || stop_reason == "stagnation")))
      break;
  }
  if (!converged && total_iterations >= args_.nl_max_iter &&
      stop_reason != "line_search_floor" && stop_reason != "nonfinite" &&
      stop_reason != "stagnation")
    stop_reason = "work_budget";

  field.nl_iterations = total_iterations;
  field.nl_blocks = blocks;
  field.nl_max_dphi = max_step;
  field.nl_residual = std::isfinite(residual) ? residual : 0.0;
  field.nl_clipped_points = max_clipped;
  field.nl_continuation_steps = continuation_steps;
  field.nl_backtracks = backtracks;
  field.newton_inner_iterations = inner_iterations;
  field.newton_step_norm = last_step_norm;
  field.newton_residual_before =
      std::isfinite(last_residual_before) ? last_residual_before : 0.0;
  field.newton_residual_after =
      std::isfinite(last_residual_after) ? last_residual_after : 0.0;
  field.newton_line_search_status = line_search_status;
  field.newton_correction_residual =
      std::isfinite(last_correction_residual) ? last_correction_residual : 0.0;
  field.newton_correction_max_residual =
      std::isfinite(last_correction_max_residual)
          ? last_correction_max_residual
          : 0.0;
  field.newton_accepted_lambda = last_accepted_lambda;
  field.newton_rejected_steps = rejected_steps;
  field.newton_jac_diag_min =
      std::isfinite(last_jac_diag_min) ? last_jac_diag_min : 0.0;
  field.newton_jac_diag_max =
      std::isfinite(last_jac_diag_max) ? last_jac_diag_max : 0.0;
  field.newton_screen_points = last_screen_points;
  field.newton_active_nonlinear_points = last_active_nonlinear_points;
  field.newton_stop_reason = stop_reason;
  field.newton_inner_levels = last_inner_levels;
  field.newton_inner_cycles = last_inner_cycles;
  field.newton_inner_residual_before =
      std::isfinite(last_inner_residual_before) ? last_inner_residual_before
                                                : 0.0;
  field.newton_inner_residual_after =
      std::isfinite(last_inner_residual_after) ? last_inner_residual_after
                                               : 0.0;
  field.newton_inner_effective = last_inner_effective;
  field.nl_status =
      finite ? (converged ? "ok"
                          : (stop_reason == "stagnation"
                                 ? "nonlinear_not_converged_stagnated"
                                 : (stop_reason == "line_search_floor"
                                        ? "nonlinear_not_converged_line_search_floor"
                                        : "nonlinear_not_converged")))
             : "nonlinear_nonfinite";
  field.stats = last_stats;
  field.stats.iterations = total_iterations;
  field.stats.converged = finite && converged;
  field.total_iterations += total_iterations;
  field.parity = compute_parity_audit(field.shape, &source, &residual_vec);
  update_newton_derivative_audit(field.shape, source, delta, alpha, clip,
                                 final_screen_scale, field.parity);
  if (map_result != nullptr &&
      (args_.parity_dump == "points" || args_.parity_dump == "all")) {
    maybe_write_parity_points(*map_result, field.shape, source, residual_vec,
                              delta, alpha, clip, final_screen_scale);
  }
  scale_nonlinear_field_to_physical_units(phi_, alpha);
  field.atom_phi = interpolate_atoms(field.shape, atoms);
  if (args_.charge_assignment == "delphi" ||
      args_.charge_assignment == "gaussian")
    field.phi_grid = phi_;

  if (map_result != nullptr)
    maybe_write_maps(field.shape, *map_result, true);
  return field;
}

Analyzer::FieldSolve Analyzer::solve_nonlinear_field(
    const std::vector<Atom> &atoms, double eps_in, double eps_out, double kappa,
    TitanExecutor &executor, const FrameResult *map_result) {
  if (args_.nl_solver == "newton" || args_.nl_multigrid == "fas")
    return solve_newton_nonlinear_field(atoms, eps_in, eps_out, kappa,
                                        executor, map_result);

  FieldSolve field =
      solve_linear_field(atoms, eps_in, eps_out, kappa, executor, nullptr);
  field.linear_seed_solver = field.solver_mode;
  field.pbe_model_effective = "nonlinear";
  field.nl_solver = "picard";
  field.nl_inner_solver = "sor";
  field.nl_mg_mode = "off";
  field.nl_relax = effective_nl_relax();

  if (kappa <= 0.0 || args_.salt <= 0.0) {
    field.pbe_model_effective = "linear";
    field.nl_status = "linearized_no_salt";
    if (map_result != nullptr) {
      maybe_write_maps(field.shape, *map_result, false);
      maybe_write_maps(field.shape, *map_result, true);
    }
    return field;
  }

  if (map_result != nullptr)
    maybe_write_maps(field.shape, *map_result, false);

  const double alpha = nonlinear_alpha();
  const double clip = std::max(1.0e-12, args_.nl_clip);
  scale_nonlinear_field_to_solve_units(phi_, alpha);
  AlignedDoubleVector source(field.shape.n, 0.0);
  SolveStats last_stats;
  int total_iterations = 0;
  int blocks = 0;
  double max_block_dphi = 0.0;
  double residual = std::numeric_limits<double>::infinity();
  bool finite = true;
  std::size_t max_clipped = 0;
  int continuation_steps = 0;
  int backtracks = 0;
  std::vector<double> screen_scales;
  if (args_.nl_continuation == "auto" && args_.salt > 0.5)
    screen_scales = {0.25, 0.5, 1.0};
  else if (args_.nl_continuation == "auto")
    screen_scales = {0.5, 1.0};
  else
    screen_scales = {1.0};

  for (double screen_scale : screen_scales) {
    ++continuation_steps;
    double previous_residual = std::numeric_limits<double>::infinity();
    bool step_done = false;
    while (total_iterations < args_.nl_max_iter && !step_done) {
      const std::size_t clipped_this =
          update_nonlinear_source(field.shape, source, alpha, clip);
      max_clipped = std::max(max_clipped, clipped_this);

      const int sweeps =
          std::min(args_.nl_block, args_.nl_max_iter - total_iterations);
      const AlignedDoubleVector before =
          args_.nl_line_search == "auto" ? phi_ : AlignedDoubleVector{};
      std::size_t unused_clipped = 0;
      last_stats = solve_nonlinear_sor_block(
          field.shape, source, sweeps, field.nl_relax, screen_scale, executor,
          unused_clipped);
      total_iterations += sweeps;
      ++blocks;
      max_block_dphi = std::max(max_block_dphi, last_stats.max_dphi);

      finite = std::all_of(phi_.begin(), phi_.end(),
                           [](double v) { return std::isfinite(v); });
      if (!finite)
        break;

      update_nonlinear_source(field.shape, source, alpha, clip);
      residual = compute_nonlinear_residual(field.shape, source, screen_scale);
      if (args_.nl_line_search == "auto" && std::isfinite(previous_residual) &&
          residual > previous_residual * 1.05 && field.nl_relax > 0.051 &&
          before.size() == phi_.size()) {
        phi_ = before;
        field.nl_relax = std::max(0.05, 0.5 * field.nl_relax);
        ++backtracks;
        residual = previous_residual;
        continue;
      }
      previous_residual = residual;

      if (last_stats.max_dphi <= args_.nl_dphi ||
          residual <= args_.nl_residual) {
        if (screen_scale >= 1.0) {
          last_stats.converged = true;
        }
        step_done = true;
      }
    }
    if (!finite || total_iterations >= args_.nl_max_iter)
      break;
  }

  field.nl_iterations = total_iterations;
  field.nl_blocks = blocks;
  field.nl_max_dphi = max_block_dphi;
  field.nl_residual = std::isfinite(residual) ? residual : 0.0;
  field.nl_clipped_points = max_clipped;
  field.nl_continuation_steps = continuation_steps;
  field.nl_backtracks = backtracks;
  field.nl_status =
      finite ? (last_stats.converged ? "ok" : "nonlinear_not_converged")
             : "nonlinear_nonfinite";
  field.stats = last_stats;
  field.stats.iterations = total_iterations;
  field.stats.converged = finite && last_stats.converged;
  field.total_iterations += total_iterations;
  scale_nonlinear_field_to_physical_units(phi_, alpha);
  field.atom_phi = interpolate_atoms(field.shape, atoms);
  if (args_.charge_assignment == "delphi" ||
      args_.charge_assignment == "gaussian")
    field.phi_grid = phi_;

  if (map_result != nullptr)
    maybe_write_maps(field.shape, *map_result, true);
  return field;
}

Analyzer::FieldSolve Analyzer::solve_physical_field(
    const std::vector<Atom> &atoms, double eps_in, double eps_out, double kappa,
    TitanExecutor &executor, const FrameResult *map_result) {
  if (args_.model == "nonlinear")
    return solve_nonlinear_field(atoms, eps_in, eps_out, kappa, executor,
                                 map_result);
  return solve_linear_field(atoms, eps_in, eps_out, kappa, executor,
                            map_result);
}

Analyzer::GridShape Analyzer::coarsen_shape(const GridShape &fine) const {
  GridShape coarse;
  coarse.ox = fine.ox;
  coarse.oy = fine.oy;
  coarse.oz = fine.oz;
  coarse.h = fine.h * 2.0;
  coarse.nx = std::max<std::size_t>(3, ((fine.nx - 1) / 2) + 1);
  coarse.ny = std::max<std::size_t>(3, ((fine.ny - 1) / 2) + 1);
  coarse.nz = std::max<std::size_t>(3, ((fine.nz - 1) / 2) + 1);
  if (coarse.nx > std::numeric_limits<std::size_t>::max() / coarse.ny ||
      coarse.nx * coarse.ny >
          std::numeric_limits<std::size_t>::max() / coarse.nz) {
    throw std::runtime_error("PBE: multigrid dimensions overflow size_t.");
  }
  coarse.n = coarse.nx * coarse.ny * coarse.nz;
  coarse.memory_bytes =
      coarse.n * (7 * sizeof(double) + 2 * sizeof(unsigned char));
  if (coarse.n > args_.max_grid_points) {
    throw std::runtime_error(
        "PBE: multigrid coarse grid exceeds -pbe-max-grid.");
  }
  return coarse;
}

int Analyzer::effective_vcycle_levels(const GridShape &fine_shape) const {
  int requested = 4;
  if (args_.mg_levels != "auto") {
    requested = std::clamp(std::stoi(args_.mg_levels), 1, 4);
  }
  int levels = 1;
  GridShape current = fine_shape;
  while (levels < requested) {
    const std::size_t next_nx = ((current.nx - 1) / 2) + 1;
    const std::size_t next_ny = ((current.ny - 1) / 2) + 1;
    const std::size_t next_nz = ((current.nz - 1) / 2) + 1;
    if (next_nx < 3 || next_ny < 3 || next_nz < 3)
      break;
    ++levels;
    current = coarsen_shape(current);
    if (current.nx == 3 && current.ny == 3 && current.nz == 3)
      break;
  }
  return levels;
}

void Analyzer::prepare_multigrid_level(
    MultigridLevel &level, const GridShape &shape,
    const std::vector<Atom> &atoms, double eps_in, double eps_out, double kappa,
    bool physical_rhs, bool physical_boundary, TitanExecutor &executor) {
  ensure_workspace(shape);
  classify_surface(shape, atoms, executor);
  apply_oracle_surface_map(shape);
  classify_ion_accessibility(shape, atoms, executor);
  const ChargeDiagnostics charge_diag =
      args_.test_charge_grid == "oracle"
          ? assign_oracle_charge_grid(shape, atoms)
          : assign_charges(shape, atoms);
  build_coefficients(shape, atoms, eps_in, eps_out, kappa, executor);
  BoundaryDiagnostics boundary_diag;
  if (physical_boundary) {
    boundary_diag = set_boundary(shape, atoms, eps_out, kappa);
    if (args_.test_boundary_map == "oracle" && oracle_replay_applies(shape))
      boundary_diag = apply_oracle_boundary_map(shape);
  } else {
    std::fill(phi_.begin(), phi_.end(), 0.0);
  }

  level.shape = shape;
  level.phi = phi_;
  level.rho = rho_;
  level.eps_x = eps_x_;
  level.eps_y = eps_y_;
  level.eps_z = eps_z_;
  level.denom = denom_;
  level.screen = screen_;
  level.inside = inside_;
  level.ion_accessible = ion_accessible_;
  level.charge_diag = charge_diag;
  level.boundary_diag = boundary_diag;
  level.rhs.assign(shape.n, 0.0);
  if (physical_rhs) {
    const double rhs = rhs_scale(shape);
    for (std::size_t id = 0; id < shape.n; ++id)
      level.rhs[id] = rhs * rho_[id];
  }
}

void Analyzer::copy_level_to_workspace(const MultigridLevel &level) {
  phi_ = level.phi;
  rho_ = level.rho;
  eps_x_ = level.eps_x;
  eps_y_ = level.eps_y;
  eps_z_ = level.eps_z;
  denom_ = level.denom;
  screen_ = level.screen;
  inside_ = level.inside;
  ion_accessible_ = level.ion_accessible;
}

Analyzer::SolveStats
Analyzer::smooth_multigrid_level(MultigridLevel &level, int iterations,
                                 bool stop_on_convergence) const {
  SolveStats stats;
  if (iterations <= 0)
    return stats;
  const GridShape &shape = level.shape;
  const double omega =
      args_.omega == "auto"
          ? std::clamp(2.0 / (1.0 + std::sin(kPi / static_cast<double>(
                                             std::max({shape.nx, shape.ny,
                                                       shape.nz})))),
                       1.0, 1.95)
          : args_.omega_value;
  for (int iter = 1; iter <= iterations; ++iter) {
    double sumsq = 0.0;
    double max_delta = 0.0;
    std::size_t count = 0;
    for (int color = 0; color < 2; ++color) {
      for (std::size_t k = 1; k + 1 < shape.nz; ++k) {
        for (std::size_t j = 1; j + 1 < shape.ny; ++j) {
          for (std::size_t i = 1; i + 1 < shape.nx; ++i) {
            if (static_cast<int>((i + j + k) & 1ULL) != color)
              continue;
            const std::size_t id = idx3(i, j, k, shape.nx, shape.ny);
            const double neighbor_sum =
                level.eps_x[id] * level.phi[id + 1] +
                level.eps_x[id - 1] * level.phi[id - 1] +
                level.eps_y[id] * level.phi[id + shape.nx] +
                level.eps_y[id - shape.nx] * level.phi[id - shape.nx] +
                level.eps_z[id] * level.phi[id + shape.nx * shape.ny] +
                level.eps_z[id - shape.nx * shape.ny] *
                    level.phi[id - shape.nx * shape.ny];
            const double target =
                (neighbor_sum + level.rhs[id]) / level.denom[id];
            const double updated =
                level.phi[id] + omega * (target - level.phi[id]);
            const double delta = updated - level.phi[id];
            level.phi[id] = updated;
            const double abs_delta = std::abs(delta);
            sumsq += delta * delta;
            max_delta = std::max(max_delta, abs_delta);
            ++count;
          }
        }
      }
    }
    stats.iterations = iter;
    stats.rms = count > 0 ? std::sqrt(sumsq / static_cast<double>(count)) : 0.0;
    stats.max_dphi = max_delta;
    if (stop_on_convergence &&
        (stats.rms <= args_.rms_tol || stats.max_dphi <= args_.dphi_tol)) {
      stats.converged = true;
      break;
    }
  }
  return stats;
}

double Analyzer::compute_multigrid_residual(
    const MultigridLevel &level, AlignedDoubleVector &residual) const {
  const GridShape &shape = level.shape;
  residual.assign(shape.n, 0.0);
  double sumsq = 0.0;
  std::size_t count = 0;
  for (std::size_t k = 1; k + 1 < shape.nz; ++k) {
    for (std::size_t j = 1; j + 1 < shape.ny; ++j) {
      for (std::size_t i = 1; i + 1 < shape.nx; ++i) {
        const std::size_t id = idx3(i, j, k, shape.nx, shape.ny);
        const double neighbor_sum =
            level.eps_x[id] * level.phi[id + 1] +
            level.eps_x[id - 1] * level.phi[id - 1] +
            level.eps_y[id] * level.phi[id + shape.nx] +
            level.eps_y[id - shape.nx] * level.phi[id - shape.nx] +
            level.eps_z[id] * level.phi[id + shape.nx * shape.ny] +
            level.eps_z[id - shape.nx * shape.ny] *
                level.phi[id - shape.nx * shape.ny];
        const double a_phi = level.denom[id] * level.phi[id] - neighbor_sum;
        const double r = level.rhs[id] - a_phi;
        residual[id] = r;
        sumsq += r * r;
        ++count;
      }
    }
  }
  return count > 0 ? std::sqrt(sumsq / static_cast<double>(count)) : 0.0;
}

void Analyzer::restrict_residual_to_coarse(
    const MultigridLevel &fine, const AlignedDoubleVector &fine_residual,
    MultigridLevel &coarse) const {
  std::fill(coarse.rhs.begin(), coarse.rhs.end(), 0.0);
  std::fill(coarse.phi.begin(), coarse.phi.end(), 0.0);
  const int weights[3] = {1, 2, 1};
  for (std::size_t kc = 1; kc + 1 < coarse.shape.nz; ++kc) {
    for (std::size_t jc = 1; jc + 1 < coarse.shape.ny; ++jc) {
      for (std::size_t ic = 1; ic + 1 < coarse.shape.nx; ++ic) {
        const std::size_t fi = std::min(fine.shape.nx - 2, 2 * ic);
        const std::size_t fj = std::min(fine.shape.ny - 2, 2 * jc);
        const std::size_t fk = std::min(fine.shape.nz - 2, 2 * kc);
        double sum = 0.0;
        int wsum = 0;
        for (int dz = -1; dz <= 1; ++dz) {
          const std::size_t iz = static_cast<std::size_t>(
              std::clamp<int>(static_cast<int>(fk) + dz, 1,
                              static_cast<int>(fine.shape.nz - 2)));
          for (int dy = -1; dy <= 1; ++dy) {
            const std::size_t iy = static_cast<std::size_t>(
                std::clamp<int>(static_cast<int>(fj) + dy, 1,
                                static_cast<int>(fine.shape.ny - 2)));
            for (int dx = -1; dx <= 1; ++dx) {
              const std::size_t ix = static_cast<std::size_t>(
                  std::clamp<int>(static_cast<int>(fi) + dx, 1,
                                  static_cast<int>(fine.shape.nx - 2)));
              const int w = weights[dx + 1] * weights[dy + 1] *
                            weights[dz + 1];
              sum += static_cast<double>(w) *
                     fine_residual[idx3(ix, iy, iz, fine.shape.nx,
                                        fine.shape.ny)];
              wsum += w;
            }
          }
        }
        coarse.rhs[idx3(ic, jc, kc, coarse.shape.nx, coarse.shape.ny)] =
            wsum > 0 ? sum / static_cast<double>(wsum) : 0.0;
      }
    }
  }
}

void Analyzer::prolongate_add_correction(const MultigridLevel &coarse,
                                         MultigridLevel &fine) const {
  for (std::size_t k = 1; k + 1 < fine.shape.nz; ++k) {
    const double z = fine.shape.oz + static_cast<double>(k) * fine.shape.h;
    for (std::size_t j = 1; j + 1 < fine.shape.ny; ++j) {
      const double y = fine.shape.oy + static_cast<double>(j) * fine.shape.h;
      for (std::size_t i = 1; i + 1 < fine.shape.nx; ++i) {
        const double x = fine.shape.ox + static_cast<double>(i) * fine.shape.h;
        double correction = 0.0;
        if (interpolate_field_value(coarse.shape, coarse.phi, x, y, z,
                                    correction)) {
          fine.phi[idx3(i, j, k, fine.shape.nx, fine.shape.ny)] += correction;
        }
      }
    }
  }
}

void Analyzer::run_multigrid_vcycle(std::vector<MultigridLevel> &levels,
                                    std::size_t level_index) {
  MultigridLevel &level = levels[level_index];
  if (level_index + 1 >= levels.size()) {
    smooth_multigrid_level(level, args_.mg_coarse_iter, true);
    return;
  }

  smooth_multigrid_level(level, args_.mg_pre, false);
  AlignedDoubleVector residual;
  compute_multigrid_residual(level, residual);
  MultigridLevel &coarse = levels[level_index + 1];
  restrict_residual_to_coarse(level, residual, coarse);
  run_multigrid_vcycle(levels, level_index + 1);
  prolongate_add_correction(coarse, level);
  smooth_multigrid_level(level, args_.mg_post, false);
}

std::size_t Analyzer::multigrid_memory_bytes(
    const std::vector<MultigridLevel> &levels) const {
  std::size_t bytes = 0;
  for (const MultigridLevel &level : levels) {
    bytes += level.phi.size() * sizeof(double);
    bytes += level.rho.size() * sizeof(double);
    bytes += level.rhs.size() * sizeof(double);
    bytes += level.eps_x.size() * sizeof(double);
    bytes += level.eps_y.size() * sizeof(double);
    bytes += level.eps_z.size() * sizeof(double);
    bytes += level.denom.size() * sizeof(double);
    bytes += level.screen.size() * sizeof(double);
    bytes += level.inside.size() * sizeof(unsigned char);
    bytes += level.ion_accessible.size() * sizeof(unsigned char);
  }
  return bytes;
}

Analyzer::FieldSolve Analyzer::solve_vcycle_field(
    const std::vector<Atom> &atoms, double eps_in, double eps_out, double kappa,
    TitanExecutor &executor, const FrameResult *map_result) {
  FieldSolve field;
  field.solver_mode = "vcycle";
  field.mg_pre = args_.mg_pre;
  field.mg_post = args_.mg_post;
  field.mg_coarse_iter = args_.mg_coarse_iter;

  const GridShape fine_shape = build_shape(atoms);
  const int levels_count = effective_vcycle_levels(fine_shape);
  field.levels = levels_count;
  field.mg_levels = levels_count;
  field.initial_scale = args_.scale;
  field.initial_margin = args_.margin;

  std::vector<MultigridLevel> levels(static_cast<std::size_t>(levels_count));
  GridShape shape = fine_shape;
  for (int level = 0; level < levels_count; ++level) {
    prepare_multigrid_level(levels[static_cast<std::size_t>(level)], shape,
                            atoms, eps_in, eps_out, kappa, level == 0,
                            level == 0, executor);
    if (level + 1 < levels_count)
      shape = coarsen_shape(shape);
  }
  field.peak_memory_bytes = multigrid_memory_bytes(levels);

  copy_level_to_workspace(levels.front());
  if (map_result != nullptr)
    maybe_write_maps(fine_shape, *map_result, false);

  const int cycle_work =
      std::max(1, args_.mg_pre + args_.mg_post + args_.mg_coarse_iter);
  const int max_cycles = std::max(1, args_.max_iter / cycle_work);
  SolveStats cycle_stats;
  for (int cycle = 1; cycle <= max_cycles; ++cycle) {
    const AlignedDoubleVector before = levels.front().phi;
    run_multigrid_vcycle(levels, 0);
    double sumsq = 0.0;
    double max_delta = 0.0;
    std::size_t count = 0;
    for (std::size_t k = 1; k + 1 < fine_shape.nz; ++k) {
      for (std::size_t j = 1; j + 1 < fine_shape.ny; ++j) {
        for (std::size_t i = 1; i + 1 < fine_shape.nx; ++i) {
          const std::size_t id = idx3(i, j, k, fine_shape.nx, fine_shape.ny);
          const double delta = levels.front().phi[id] - before[id];
          sumsq += delta * delta;
          max_delta = std::max(max_delta, std::abs(delta));
          ++count;
        }
      }
    }
    cycle_stats.iterations = cycle;
    cycle_stats.rms =
        count > 0 ? std::sqrt(sumsq / static_cast<double>(count)) : 0.0;
    cycle_stats.max_dphi = max_delta;
    if (cycle_stats.rms <= args_.rms_tol ||
        cycle_stats.max_dphi <= args_.dphi_tol) {
      cycle_stats.converged = true;
      break;
    }
  }

  field.mg_cycles = cycle_stats.iterations;
  copy_level_to_workspace(levels.front());
  SolveStats final_stats = solve_sor(fine_shape, executor);
  field.total_iterations =
      field.mg_cycles * cycle_work + final_stats.iterations;
  field.stats = final_stats;
  field.rhs_scale = rhs_scale(fine_shape);
  field.charge_assigned = levels.front().charge_diag.assigned_charge;
  field.charge_moment_error = levels.front().charge_diag.moment_error;
  field.gaussian_density_sum = last_gaussian_density_sum_;
  field.gaussian_cutoff_points = last_gaussian_cutoff_points_;
  field.gaussian_eps_min = last_gaussian_eps_min_;
  field.gaussian_eps_max = last_gaussian_eps_max_;
  field.gaussian_eps_mean = last_gaussian_eps_mean_;
  field.vdwms_enabled = last_vdwms_audit_.enabled;
  field.vdwms_boundary_initial = last_vdwms_audit_.boundary_initial;
  field.vdwms_boundary_final = last_vdwms_audit_.boundary_final;
  field.vdwms_added = last_vdwms_audit_.added;
  field.vdwms_removed = last_vdwms_audit_.removed;
  field.vdwms_exposed_points = last_vdwms_audit_.exposed_points;
  field.vdwms_midpoint_rms = last_vdwms_audit_.midpoint_rms;
  field.gaussian_charge_norm_error =
      std::abs(levels.front().charge_diag.assigned_charge -
               std::accumulate(atoms.begin(), atoms.end(), 0.0,
                               [](double sum, const Atom &atom) {
                                 return sum + atom.q;
                               }));
  field.boundary_rms = levels.front().boundary_diag.rms;
  field.boundary_max = levels.front().boundary_diag.max_abs;
  field.linear_residual = compute_linear_residual(fine_shape);
  field.parity = compute_parity_audit(fine_shape, nullptr, nullptr);
  if (map_result != nullptr)
    maybe_write_maps(fine_shape, *map_result, true);
  field.shape = fine_shape;
  field.solute_points = static_cast<std::size_t>(
      std::count(inside_.begin(), inside_.end(), 1U));
  field.solvent_points = fine_shape.n - field.solute_points;
  field.ion_accessible_points = static_cast<std::size_t>(
      std::count(ion_accessible_.begin(), ion_accessible_.end(), 1U));
  field.stern_excluded_points =
      args_.ion_access == "stern" ? fine_shape.n - field.ion_accessible_points
                                  : 0;
  field.atom_phi = interpolate_atoms(fine_shape, atoms);
  const bool keep_field_grids =
      parity_map_dump_enabled() ||
      args_.test_energy_convention == "delphi-iscm" ||
      args_.test_energy_convention == "delphi-corrected" ||
      barnase_forensic_enabled();
  if (args_.charge_assignment == "delphi" ||
      args_.charge_assignment == "gaussian" || keep_field_grids)
    field.phi_grid = phi_;
  if (keep_field_grids) {
    field.rho_grid = rho_;
    field.screen_grid = screen_;
    field.inside_grid = inside_;
    field.eps_x_grid = eps_x_;
    field.eps_y_grid = eps_y_;
    field.eps_z_grid = eps_z_;
    field.eps_sum_grid.assign(fine_shape.n, 0.0);
    for (std::size_t id = 0; id < fine_shape.n; ++id) {
      field.eps_sum_grid[id] =
          denom_.size() == fine_shape.n && screen_.size() == fine_shape.n
              ? std::max(1.0e-12, denom_[id] - screen_[id])
              : 0.0;
    }
  }
  return field;
}

std::vector<double>
Analyzer::interpolate_atoms(const GridShape &shape,
                            const std::vector<Atom> &atoms) const {
  std::vector<double> potentials;
  potentials.reserve(atoms.size());
  for (const Atom &atom : atoms) {
    double gx = (atom.x - shape.ox) / shape.h;
    double gy = (atom.y - shape.oy) / shape.h;
    double gz = (atom.z - shape.oz) / shape.h;
    gx = std::clamp(gx, 0.0, static_cast<double>(shape.nx - 1));
    gy = std::clamp(gy, 0.0, static_cast<double>(shape.ny - 1));
    gz = std::clamp(gz, 0.0, static_cast<double>(shape.nz - 1));
    std::size_t i0 = std::min<std::size_t>(
        shape.nx - 2, static_cast<std::size_t>(std::floor(gx)));
    std::size_t j0 = std::min<std::size_t>(
        shape.ny - 2, static_cast<std::size_t>(std::floor(gy)));
    std::size_t k0 = std::min<std::size_t>(
        shape.nz - 2, static_cast<std::size_t>(std::floor(gz)));
    const double fx = clamp_fraction(gx - static_cast<double>(i0));
    const double fy = clamp_fraction(gy - static_cast<double>(j0));
    const double fz = clamp_fraction(gz - static_cast<double>(k0));
    double value = 0.0;
    for (std::size_t dz = 0; dz <= 1; ++dz) {
      const double wz = dz ? fz : (1.0 - fz);
      for (std::size_t dy = 0; dy <= 1; ++dy) {
        const double wy = dy ? fy : (1.0 - fy);
        for (std::size_t dx = 0; dx <= 1; ++dx) {
          const double wx = dx ? fx : (1.0 - fx);
          value += wx * wy * wz *
                   phi_[idx3(i0 + dx, j0 + dy, k0 + dz, shape.nx, shape.ny)];
        }
      }
    }
    potentials.push_back(value);
  }
  return potentials;
}

std::vector<double>
Analyzer::compute_sasa_per_atom(const std::vector<Atom> &atoms,
                                TitanExecutor &executor) const {
  std::vector<TitanSASA::Atom> sasa_atoms;
  sasa_atoms.reserve(atoms.size());
  for (const Atom &atom : atoms)
    sasa_atoms.push_back({atom.x, atom.y, atom.z, atom.radius});

  TitanSASA::Options options;
  options.probe_radius = args_.probe_radius;
  options.min_cell_size = 1.0;
  options.point_count = TitanSASA::kDefaultPointCount;
  return TitanSASA::compute_per_atom_sasa(
      std::span<const TitanSASA::Atom>(sasa_atoms.data(), sasa_atoms.size()),
      options, {}, &executor);
}

double Analyzer::compute_sasa_total(const std::vector<Atom> &atoms,
                                    TitanExecutor &executor) const {
  const std::vector<double> area_by_atom =
      compute_sasa_per_atom(atoms, executor);
  return TitanSASA::sum_sasa(
      std::span<const double>(area_by_atom.data(), area_by_atom.size()));
}

Analyzer::PbsaComponents Analyzer::compute_pbsa_components(
    const MOL &mol1, const Multi_MOL_parameters &params1, const MOL &mol2,
    const Multi_MOL_parameters &params2, const std::vector<Atom> &mol1_atoms,
    const std::vector<Atom> &mol2_atoms,
    const std::vector<Atom> &complex_atoms, double deltaG_solv_kj_mol,
    TitanExecutor &executor,
    const EnergyDecompositionProvider *energy_provider) const {
  (void)params1;
  (void)params2;
  if (energy_provider == nullptr || !energy_provider->ready()) {
    throw std::runtime_error(
        "PBSA: -pbsa requires initialized two-molecule force-field "
        "parameters for elec/vdw interaction energy.");
  }

  PbsaComponents comp;
  const std::size_t n1 = mol1.size();
  const std::size_t n2 = mol2.size();
  const std::size_t workers = effective_workers(executor, n1);
  std::vector<double> local_elec(workers, 0.0);
  std::vector<double> local_vdw(workers, 0.0);
  executor.parallel_for(
      n1, workers, 1,
      [&](std::size_t worker, std::size_t begin, std::size_t end) {
        double elec = 0.0;
        double vdw = 0.0;
        for (std::size_t i = begin; i < end; ++i) {
          for (std::size_t j = 0; j < n2; ++j) {
            const EnergyPairDecomposition pair =
                energy_provider->pair_energy(mol1, mol2, i, j);
            elec += pair.elec;
            vdw += pair.vdw;
          }
        }
        local_elec[worker] = elec;
        local_vdw[worker] = vdw;
      });
  for (std::size_t worker = 0; worker < workers; ++worker) {
    comp.elec_kj_mol += local_elec[worker];
    comp.vdw_kj_mol += local_vdw[worker];
  }
  comp.mm_kj_mol = comp.elec_kj_mol + comp.vdw_kj_mol;

  comp.sasa_mol1_a2 = compute_sasa_total(mol1_atoms, executor);
  comp.sasa_mol2_a2 = compute_sasa_total(mol2_atoms, executor);
  comp.sasa_complex_a2 = compute_sasa_total(complex_atoms, executor);
  comp.delta_sasa_a2 =
      comp.sasa_complex_a2 - comp.sasa_mol1_a2 - comp.sasa_mol2_a2;
  comp.nonpolar_kj_mol =
      TitanSASA::nonpolar_energy_from_delta_sasa_kj(comp.delta_sasa_a2);
  comp.delta_solv_kj_mol = deltaG_solv_kj_mol;
  comp.no_entropy_kj_mol =
      comp.mm_kj_mol + comp.delta_solv_kj_mol + comp.nonpolar_kj_mol;
  return comp;
}

std::vector<Analyzer::ResidueFrameRow> Analyzer::compute_pbsa_residue_rows(
    const MOL &mol1, const MOL &mol2, const std::vector<Atom> &mol1_atoms,
    const std::vector<Atom> &mol2_atoms,
    const std::vector<Atom> &complex_atoms,
    const ReactionSolve &complex_solve, const ReactionSolve &mol1_solve,
    const ReactionSolve &mol2_solve, const FrameResult &result,
    TitanExecutor &executor,
    const EnergyDecompositionProvider *energy_provider) const {
  std::vector<ResidueFrameRow> rows;
  if (args_.pbsa_decomp != "residue" || energy_provider == nullptr ||
      !energy_provider->ready())
    return rows;

  const std::size_t n1 = mol1_atoms.size();
  const std::size_t n2 = mol2_atoms.size();
  if (complex_solve.atom_reaction_kj_mol.size() < n1 + n2 ||
      mol1_solve.atom_reaction_kj_mol.size() < n1 ||
      mol2_solve.atom_reaction_kj_mol.size() < n2)
    return rows;

  auto make_key = [](const char *molecule, const mol_info *info,
                     std::size_t atom_index) {
    ResidueKey key;
    key.molecule = molecule;
    if (info != nullptr) {
      key.chainnum = info->get_chainnum();
      key.chain = info->get_chain();
      key.resnr = info->get_resnr();
      key.respos = info->get_respos();
      key.resname = info->get_res_name();
    } else {
      key.resnr = atom_index + 1;
      key.respos = atom_index + 1;
      key.resname = "UNK";
    }
    return key;
  };
  auto key_id = [](const ResidueKey &key) {
    return key.molecule + "|" + std::to_string(key.chainnum) + "|" +
           key.chain + "|" + std::to_string(key.resnr) + "|" +
           std::to_string(key.respos) + "|" + key.resname;
  };
  auto add = [&](std::map<std::string, ResidueContribution> &accum,
                 const ResidueKey &key, double elec, double vdw, double polar,
                 double nonpolar) {
    ResidueContribution &entry = accum[key_id(key)];
    entry.key = key;
    entry.elec_kj_mol += elec;
    entry.vdw_kj_mol += vdw;
    entry.mm_kj_mol += elec + vdw;
    entry.polar_kj_mol += polar;
    entry.nonpolar_kj_mol += nonpolar;
  };

  std::map<std::string, ResidueContribution> accum;
  for (std::size_t i = 0; i < n1; ++i) {
    const ResidueKey rec_key =
        make_key("receptor", energy_provider->rec_info(i), i);
    for (std::size_t j = 0; j < n2; ++j) {
      const ResidueKey lig_key =
          make_key("ligand", energy_provider->lig_info(j), j);
      const EnergyPairDecomposition pair =
          energy_provider->pair_energy(mol1, mol2, i, j);
      add(accum, rec_key, 0.5 * pair.elec, 0.5 * pair.vdw, 0.0, 0.0);
      add(accum, lig_key, 0.5 * pair.elec, 0.5 * pair.vdw, 0.0, 0.0);
    }
  }

  const std::vector<double> sasa_complex =
      compute_sasa_per_atom(complex_atoms, executor);
  const std::vector<double> sasa_mol1 =
      compute_sasa_per_atom(mol1_atoms, executor);
  const std::vector<double> sasa_mol2 =
      compute_sasa_per_atom(mol2_atoms, executor);

  for (std::size_t i = 0; i < n1; ++i) {
    const ResidueKey rec_key =
        make_key("receptor", energy_provider->rec_info(i), i);
    const double polar =
        complex_solve.atom_reaction_kj_mol[i] -
        mol1_solve.atom_reaction_kj_mol[i];
    const double delta_sasa =
        (i < sasa_complex.size() ? sasa_complex[i] : 0.0) -
        (i < sasa_mol1.size() ? sasa_mol1[i] : 0.0);
    add(accum, rec_key, 0.0, 0.0, polar,
        TitanSASA::nonpolar_energy_from_delta_sasa_kj(delta_sasa));
  }
  for (std::size_t j = 0; j < n2; ++j) {
    const ResidueKey lig_key =
        make_key("ligand", energy_provider->lig_info(j), j);
    const std::size_t complex_index = n1 + j;
    const double polar =
        complex_solve.atom_reaction_kj_mol[complex_index] -
        mol2_solve.atom_reaction_kj_mol[j];
    const double delta_sasa =
        (complex_index < sasa_complex.size() ? sasa_complex[complex_index]
                                             : 0.0) -
        (j < sasa_mol2.size() ? sasa_mol2[j] : 0.0);
    add(accum, lig_key, 0.0, 0.0, polar,
        TitanSASA::nonpolar_energy_from_delta_sasa_kj(delta_sasa));
  }

  rows.reserve(accum.size());
  for (auto &[id, contribution] : accum) {
    (void)id;
    contribution.total_kj_mol = contribution.mm_kj_mol +
                                contribution.polar_kj_mol +
                                contribution.nonpolar_kj_mol;
    ResidueFrameRow row;
    static_cast<ResidueContribution &>(row) = contribution;
    row.trajectory_id = result.trajectory_id;
    row.frame_index = result.frame_index;
    row.global_index = result.global_index;
    row.time_ps = result.time_ps;
    rows.push_back(row);
  }
  return rows;
}

void Analyzer::accumulate_pbsa_residue_rows(
    const FrameResult &result, const std::vector<ResidueFrameRow> &rows) {
  if (rows.empty())
    return;

  auto key_id = [](const ResidueKey &key) {
    return key.molecule + "|" + std::to_string(key.chainnum) + "|" +
           key.chain + "|" + std::to_string(key.resnr) + "|" +
           std::to_string(key.respos) + "|" + key.resname;
  };
  auto add_stats = [&](ResidueStats &stats, const ResidueFrameRow &row) {
    stats.key = row.key;
    ++stats.frames;
    add_scalar(stats.elec, row.elec_kj_mol);
    add_scalar(stats.vdw, row.vdw_kj_mol);
    add_scalar(stats.mm, row.mm_kj_mol);
    add_scalar(stats.polar, row.polar_kj_mol);
    add_scalar(stats.nonpolar, row.nonpolar_kj_mol);
    add_scalar(stats.total, row.total_kj_mol);
  };

  double elec = 0.0;
  double vdw = 0.0;
  double polar = 0.0;
  double nonpolar = 0.0;
  double total = 0.0;
  for (const ResidueFrameRow &row : rows) {
    const std::string id = key_id(row.key);
    add_stats(pbsa_residue_stats_[id], row);
    add_stats(pbsa_residue_traj_stats_[std::to_string(row.trajectory_id) +
                                       "|" + id],
              row);
    elec += row.elec_kj_mol;
    vdw += row.vdw_kj_mol;
    polar += row.polar_kj_mol;
    nonpolar += row.nonpolar_kj_mol;
    total += row.total_kj_mol;
  }

  pbsa_residue_closure_elec_max_abs_kj_ =
      std::max(pbsa_residue_closure_elec_max_abs_kj_,
               std::abs(elec - result.pbsa_elec_kj_mol));
  pbsa_residue_closure_vdw_max_abs_kj_ =
      std::max(pbsa_residue_closure_vdw_max_abs_kj_,
               std::abs(vdw - result.pbsa_vdw_kj_mol));
  pbsa_residue_closure_polar_max_abs_kj_ =
      std::max(pbsa_residue_closure_polar_max_abs_kj_,
               std::abs(polar - result.pbsa_delta_solv_kj_mol));
  pbsa_residue_closure_nonpolar_max_abs_kj_ =
      std::max(pbsa_residue_closure_nonpolar_max_abs_kj_,
               std::abs(nonpolar - result.pbsa_nonpolar_kj_mol));
  pbsa_residue_closure_max_abs_kj_ =
      std::max(pbsa_residue_closure_max_abs_kj_,
               std::abs(total - result.pbsa_no_entropy_kj_mol));

  if (args_.pbsa_residue_frames == "all") {
    pbsa_residue_frame_rows_.insert(pbsa_residue_frame_rows_.end(),
                                    rows.begin(), rows.end());
  }
}

void Analyzer::write_energy_header_if_needed() {
  if (!energy_out_.is_open()) {
    energy_out_.open("pbe_energy.csv", std::ios::out | std::ios::trunc);
    if (!energy_out_.is_open())
      throw std::runtime_error("PBE: failed to open pbe_energy.csv.");
  }
  if (!energy_header_written_) {
    energy_out_
        << "trajectory_id,frame_index,global_index,time_ps,atom_count,"
           "total_charge,grid_nx,grid_ny,grid_nz,grid_points,memory_bytes,"
           "surface,probe_A,solute_points,solvent_points,"
           "salt_M,iterations_solvent,iterations_vacuum,rms_solvent,"
           "rms_vacuum,max_dphi_solvent,max_dphi_vacuum,status,"
           "E_rxn_kT,E_rxn_kJmol,E_rxn_kcalmol,deltaG_solv_kT,"
           "deltaG_solv_kJmol,deltaG_solv_kcalmol,G_solv_complex_kJmol,"
           "G_solv_mol1_kJmol,G_solv_mol2_kJmol,G_solv_complex_kcalmol,"
           "G_solv_mol1_kcalmol,G_solv_mol2_kcalmol,dgsolv_status,"
           "focus_enabled,focus_levels,focus_initial_scale,"
           "focus_initial_margin,iterations_solvent_final,"
           "iterations_vacuum_final,solver_mode,mg_levels,mg_cycles,mg_pre,"
           "mg_post,mg_coarse_iter,surface_effective,ion_access,"
           "ion_radius_A,ion_accessible_points,stern_excluded_points,"
           "face_refine,pbe_model,pbe_model_effective,linear_seed_solver,"
           "nl_iterations,nl_blocks,nl_relax,nl_max_dphi,"
           "nl_clipped_points,nl_status,boundary_mode,charge_assignment,"
           "rhs_scale,charge_assigned,charge_moment_error,boundary_rms,"
           "linear_residual,nl_residual,nl_continuation_steps,nl_backtracks,"
           "grid_origin_x_A,grid_origin_y_A,grid_origin_z_A,grid_center_x_A,"
           "grid_center_y_A,grid_center_z_A,grid_spacing_A,boundary_max,"
           "reaction_energy_source,solvent_phase_rf_kT,"
           "solvent_phase_rf_kJmol,solvent_phase_rf_kcalmol,"
           "vacuum_phase_rf_kT,vacuum_phase_rf_kJmol,"
           "vacuum_phase_rf_kcalmol,atom_corrected_rxn_kT,"
           "atom_corrected_rxn_kJmol,atom_corrected_rxn_kcalmol,"
           "grid_corrected_rxn_kT,grid_corrected_rxn_kJmol,"
           "grid_corrected_rxn_kcalmol,coulomb_energy_kT,"
           "coulomb_energy_kJmol,coulomb_energy_kcalmol,pbe_biomodel,"
           "dielectric_model,nl_solver,nl_inner_solver,nl_mg_mode,"
           "gaussian_exp,gaussian_sigma,density_cutoff,surface_cutoff,"
           "gap_diel,radius_convention,focus_parent_grid,focus_child_grid,"
           "probe_status,newton_inner_iterations,newton_step_norm,"
           "newton_residual_before,newton_residual_after,"
           "newton_line_search_status,newton_correction_residual,"
           "newton_correction_max_residual,newton_accepted_lambda,"
           "newton_rejected_steps,newton_jac_diag_min,newton_jac_diag_max,"
           "newton_screen_points,newton_active_nonlinear_points,"
           "newton_stop_reason,gaussian_density_sum,"
           "gaussian_cutoff_points,gaussian_eps_min,gaussian_eps_max,"
           "gaussian_eps_mean,gaussian_charge_norm_error,"
           "newton_inner_levels,newton_inner_cycles,"
           "newton_inner_residual_before,newton_inner_residual_after,"
           "newton_inner_effective,iscm_boundary_points,"
           "iscm_charged_boundary_points,iscm_induced_charge_sum,"
           "iscm_induced_charge_l1,iscm_induced_charge_l2,"
           "iscm_induced_charge_max,iscm_water_rf_kT,"
           "iscm_water_rf_kJmol,iscm_water_rf_kcalmol,iscm_vacuum_rf_kT,"
           "iscm_vacuum_rf_kJmol,iscm_vacuum_rf_kcalmol,"
           "iscm_corrected_rxn_kT,iscm_corrected_rxn_kJmol,"
           "iscm_corrected_rxn_kcalmol,iscm_status,delphi_face_crossings,"
           "delphi_eps_face_matches,delphi_boundary_grid_points,"
           "delphi_charged_boundary_points,delphi_iscm_status,"
           "vdwms_enabled,vdwms_boundary_initial,vdwms_boundary_final,"
           "vdwms_added,vdwms_removed,vdwms_exposed_points,"
           "vdwms_midpoint_rms\n";
    energy_header_written_ = true;
  }
}

void Analyzer::append_energy_row(const FrameResult &result) {
  write_energy_header_if_needed();
  energy_out_ << result.trajectory_id << ',' << result.frame_index << ','
              << result.global_index << ',' << std::setprecision(10)
              << result.time_ps << ',' << result.atom_count << ','
              << result.total_charge << ',' << result.nx << ',' << result.ny
              << ',' << result.nz << ',' << result.grid_points << ','
              << result.memory_bytes << ',' << args_.surface << ','
              << args_.probe_radius << ',' << result.solute_points << ','
              << result.solvent_points << ',' << args_.salt << ','
              << result.solvent_iterations << ',' << result.vacuum_iterations
              << ',' << result.solvent_rms << ',' << result.vacuum_rms << ','
              << result.solvent_max_dphi << ',' << result.vacuum_max_dphi
              << ',' << result.status << ',' << result.energy_kt << ','
              << result.energy_kj_mol << ',' << result.energy_kcal_mol << ',';
  if (result.dgsolv_enabled) {
    energy_out_ << result.deltaG_solv_kt << ',' << result.deltaG_solv_kj_mol
                << ',' << result.deltaG_solv_kcal_mol << ','
                << result.complex_energy_kj_mol << ','
                << result.mol1_energy_kj_mol << ','
                << result.mol2_energy_kj_mol << ','
                << result.complex_energy_kcal_mol << ','
                << result.mol1_energy_kcal_mol << ','
                << result.mol2_energy_kcal_mol << ','
                << result.dgsolv_status << ',';
  } else {
    energy_out_ << ",,,,,,,,,off,";
  }
  energy_out_ << (result.focus_enabled ? "true" : "false") << ','
              << result.focus_levels << ',' << result.focus_initial_scale << ','
              << result.focus_initial_margin << ','
              << result.solvent_iterations_final << ','
              << result.vacuum_iterations_final << ',' << result.solver_mode
              << ',' << result.mg_levels << ',' << result.mg_cycles << ','
              << result.mg_pre << ',' << result.mg_post << ','
              << result.mg_coarse_iter << ',' << result.surface_effective
              << ',' << result.ion_access << ',' << result.ion_radius_A << ','
              << result.ion_accessible_points << ','
              << result.stern_excluded_points << ','
              << (result.face_refine ? "true" : "false") << ','
              << result.pbe_model << ',' << result.pbe_model_effective << ','
              << result.linear_seed_solver << ',' << result.nl_iterations
              << ',' << result.nl_blocks << ',' << result.nl_relax << ','
              << result.nl_max_dphi << ',' << result.nl_clipped_points << ','
              << result.nl_status << ',' << result.boundary_mode << ','
              << result.charge_assignment << ',' << result.rhs_scale << ','
              << result.charge_assigned << ',' << result.charge_moment_error
              << ',' << result.boundary_rms << ',' << result.linear_residual
              << ',' << result.nl_residual << ','
              << result.nl_continuation_steps << ',' << result.nl_backtracks
              << ',' << result.grid_origin_x_A << ',' << result.grid_origin_y_A
              << ',' << result.grid_origin_z_A << ',' << result.grid_center_x_A
              << ',' << result.grid_center_y_A << ',' << result.grid_center_z_A
              << ',' << result.grid_spacing_A << ',' << result.boundary_max
              << ',' << result.reaction_energy_source << ','
              << result.solvent_phase_kt << ',' << result.solvent_phase_kj_mol
              << ',' << result.solvent_phase_kcal_mol << ','
              << result.vacuum_phase_kt << ',' << result.vacuum_phase_kj_mol
              << ',' << result.vacuum_phase_kcal_mol << ','
              << result.atom_corrected_kt << ','
              << result.atom_corrected_kj_mol << ','
              << result.atom_corrected_kcal_mol << ','
              << result.grid_corrected_kt << ','
              << result.grid_corrected_kj_mol << ','
              << result.grid_corrected_kcal_mol << ','
              << result.coulomb_energy_kt << ','
              << result.coulomb_energy_kj_mol << ','
              << result.coulomb_energy_kcal_mol << ','
              << result.pbe_biomodel << ',' << result.dielectric_model << ','
              << result.nl_solver << ',' << result.nl_inner_solver << ','
              << result.nl_mg_mode << ',' << result.gaussian_exp << ','
              << result.gaussian_sigma << ',' << result.density_cutoff << ','
              << result.surface_cutoff << ',' << result.gap_diel << ','
              << result.radius_convention << ','
              << csv_escape(result.focus_parent_grid) << ','
              << csv_escape(result.focus_child_grid) << ','
              << csv_escape(result.probe_status) << ','
              << result.newton_inner_iterations << ','
              << result.newton_step_norm << ','
              << result.newton_residual_before << ','
              << result.newton_residual_after << ','
              << csv_escape(result.newton_line_search_status) << ','
              << result.newton_correction_residual << ','
              << result.newton_correction_max_residual << ','
              << result.newton_accepted_lambda << ','
              << result.newton_rejected_steps << ','
              << result.newton_jac_diag_min << ','
              << result.newton_jac_diag_max << ','
              << result.newton_screen_points << ','
              << result.newton_active_nonlinear_points << ','
              << csv_escape(result.newton_stop_reason) << ','
              << result.gaussian_density_sum << ','
              << result.gaussian_cutoff_points << ','
              << result.gaussian_eps_min << ','
              << result.gaussian_eps_max << ','
              << result.gaussian_eps_mean << ','
              << result.gaussian_charge_norm_error << ','
              << result.newton_inner_levels << ','
              << result.newton_inner_cycles << ','
              << result.newton_inner_residual_before << ','
              << result.newton_inner_residual_after << ','
              << csv_escape(result.newton_inner_effective) << ','
              << result.iscm_boundary_points << ','
              << result.iscm_charged_boundary_points << ','
              << result.iscm_induced_charge_sum << ','
              << result.iscm_induced_charge_l1 << ','
              << result.iscm_induced_charge_l2 << ','
              << result.iscm_induced_charge_max << ','
              << result.iscm_water_phase_kt << ','
              << result.iscm_water_phase_kj_mol << ','
              << result.iscm_water_phase_kcal_mol << ','
              << result.iscm_vacuum_phase_kt << ','
              << result.iscm_vacuum_phase_kj_mol << ','
              << result.iscm_vacuum_phase_kcal_mol << ','
              << result.iscm_corrected_kt << ','
              << result.iscm_corrected_kj_mol << ','
              << result.iscm_corrected_kcal_mol << ','
              << csv_escape(result.iscm_status) << ','
              << result.delphi_face_crossings << ','
              << result.delphi_eps_face_matches << ','
              << result.delphi_boundary_grid_points << ','
              << result.delphi_charged_boundary_points << ','
              << csv_escape(result.delphi_iscm_status) << ','
              << csv_escape(result.vdwms_enabled) << ','
              << result.vdwms_boundary_initial << ','
              << result.vdwms_boundary_final << ',' << result.vdwms_added
              << ',' << result.vdwms_removed << ','
              << result.vdwms_exposed_points << ','
              << result.vdwms_midpoint_rms
              << '\n';
}

void Analyzer::maybe_write_parity_stage(const FrameResult &result,
                                        const ReactionSolve &solve) const {
  if (args_.test_pqr.empty() || args_.parity_dump == "off")
    return;
  std::ofstream out;
  if (!parity_stage_header_written_) {
    out.open("pbe_parity_stage.csv", std::ios::out | std::ios::trunc);
    if (!out.is_open())
      throw std::runtime_error("PBE: failed to open pbe_parity_stage.csv.");
    out << "trajectory_id,frame_index,global_index,source_label,pbe_model,"
           "pbe_model_effective,solver_mode,nl_solver,nl_inner_solver,"
           "nl_inner_effective,nl_status,boundary_mode,charge_assignment,"
           "surface_effective,ion_access,test_nl_source,"
           "test_energy_convention,grid_nx,grid_ny,grid_nz,grid_points,"
           "grid_origin_x_A,grid_origin_y_A,grid_origin_z_A,"
           "grid_center_x_A,grid_center_y_A,grid_center_z_A,grid_spacing_A,"
           "total_charge,charge_assigned,charge_moment_error,solute_points,"
           "solvent_points,ion_accessible_points,stern_excluded_points,"
           "eps_sum_min,eps_sum_max,eps_sum_mean,screen_points,screen_sum,"
           "screen_rms,screen_max,boundary_rms,boundary_max,rhs_scale,"
           "linear_residual,nl_residual,nl_continuation_steps,nl_backtracks,"
           "nl_clipped_points,source_sum,source_rms,source_max,"
           "nonlinear_residual_max_abs,newton_residual_before,"
           "newton_residual_after,newton_correction_residual,"
           "newton_correction_max_residual,newton_inner_residual_before,"
           "newton_inner_residual_after,newton_inner_iterations,"
           "newton_inner_levels,newton_inner_cycles,newton_jac_diag_min,"
           "newton_jac_diag_max,newton_screen_points,"
           "newton_active_nonlinear_points,newton_accepted_lambda,"
           "newton_rejected_steps,newton_line_search_status,"
           "newton_stop_reason,newton_derivative_samples,"
           "newton_derivative_max_abs,newton_derivative_max_rel,"
           "reaction_energy_source,solvent_phase_rf_kT,"
           "solvent_phase_rf_kJmol,vacuum_phase_rf_kT,"
           "vacuum_phase_rf_kJmol,atom_corrected_rxn_kT,"
           "atom_corrected_rxn_kJmol,grid_corrected_rxn_kT,"
           "grid_corrected_rxn_kJmol,coulomb_energy_kT,"
           "coulomb_energy_kJmol,final_E_rxn_kT,final_E_rxn_kJmol,"
           "iscm_boundary_points,iscm_charged_boundary_points,"
           "iscm_induced_charge_sum,iscm_induced_charge_l1,"
           "iscm_induced_charge_l2,iscm_induced_charge_max,"
           "iscm_water_rf_kT,iscm_water_rf_kJmol,iscm_vacuum_rf_kT,"
           "iscm_vacuum_rf_kJmol,iscm_corrected_rxn_kT,"
           "iscm_corrected_rxn_kJmol,iscm_status,parity_map_dump,"
           "test_charge_grid,test_surface_map,test_boundary_map,"
           "status,test_midpoint_epsilon,delphi_face_crossings,"
           "delphi_eps_face_matches,delphi_boundary_grid_points,"
            "delphi_charged_boundary_points,delphi_iscm_status,"
            "test_vdwms_remap,vdwms_enabled,vdwms_boundary_initial,"
            "vdwms_boundary_final,vdwms_added,vdwms_removed,"
            "vdwms_exposed_points,vdwms_midpoint_rms,"
            "midpoint_eps_sum_min,midpoint_eps_sum_max,"
            "midpoint_eps_sum_mean,midpoint_eps_sum_interior_mean\n";
    parity_stage_header_written_ = true;
  } else {
    out.open("pbe_parity_stage.csv", std::ios::out | std::ios::app);
    if (!out.is_open())
      throw std::runtime_error("PBE: failed to append pbe_parity_stage.csv.");
  }

  out << std::setprecision(12);
  const ParityAudit &audit = solve.parity;
  out << result.trajectory_id << ',' << result.frame_index << ','
      << result.global_index << ',' << csv_escape(current_source_label_) << ','
      << result.pbe_model << ',' << result.pbe_model_effective << ','
      << result.solver_mode << ',' << result.nl_solver << ','
      << result.nl_inner_solver << ',' << result.newton_inner_effective << ','
      << result.nl_status << ',' << result.boundary_mode << ','
      << result.charge_assignment << ',' << result.surface_effective << ','
      << result.ion_access << ',' << args_.test_nl_source << ','
      << args_.test_energy_convention << ',' << result.nx << ','
      << result.ny << ',' << result.nz << ',' << result.grid_points << ','
      << result.grid_origin_x_A << ',' << result.grid_origin_y_A << ','
      << result.grid_origin_z_A << ',' << result.grid_center_x_A << ','
      << result.grid_center_y_A << ',' << result.grid_center_z_A << ','
      << result.grid_spacing_A << ',' << result.total_charge << ','
      << result.charge_assigned << ',' << result.charge_moment_error << ','
      << result.solute_points << ',' << result.solvent_points << ','
      << result.ion_accessible_points << ',' << result.stern_excluded_points
      << ',' << audit.eps_sum_min << ',' << audit.eps_sum_max << ','
      << audit.eps_sum_mean << ',' << audit.screen_points << ','
      << audit.screen_sum << ',' << audit.screen_rms << ','
      << audit.screen_max << ',' << result.boundary_rms << ','
      << result.boundary_max << ',' << result.rhs_scale << ','
      << result.linear_residual << ',' << result.nl_residual << ','
      << result.nl_continuation_steps << ',' << result.nl_backtracks << ','
      << result.nl_clipped_points << ',' << audit.source_sum << ','
      << audit.source_rms << ',' << audit.source_max << ','
      << audit.residual_max_abs << ',' << result.newton_residual_before << ','
      << result.newton_residual_after << ','
      << result.newton_correction_residual << ','
      << result.newton_correction_max_residual << ','
      << result.newton_inner_residual_before << ','
      << result.newton_inner_residual_after << ','
      << result.newton_inner_iterations << ',' << result.newton_inner_levels
      << ',' << result.newton_inner_cycles << ','
      << result.newton_jac_diag_min << ',' << result.newton_jac_diag_max << ','
      << result.newton_screen_points << ','
      << result.newton_active_nonlinear_points << ','
      << result.newton_accepted_lambda << ',' << result.newton_rejected_steps
      << ',' << csv_escape(result.newton_line_search_status) << ','
      << csv_escape(result.newton_stop_reason) << ','
      << audit.derivative_samples << ',' << audit.derivative_max_abs << ','
      << audit.derivative_max_rel << ','
      << csv_escape(result.reaction_energy_source) << ','
      << result.solvent_phase_kt << ',' << result.solvent_phase_kj_mol << ','
      << result.vacuum_phase_kt << ',' << result.vacuum_phase_kj_mol << ','
      << result.atom_corrected_kt << ',' << result.atom_corrected_kj_mol
      << ',' << result.grid_corrected_kt << ','
      << result.grid_corrected_kj_mol << ',' << result.coulomb_energy_kt
      << ',' << result.coulomb_energy_kj_mol << ',' << result.energy_kt
      << ',' << result.energy_kj_mol << ',' << result.iscm_boundary_points
      << ',' << result.iscm_charged_boundary_points << ','
      << result.iscm_induced_charge_sum << ',' << result.iscm_induced_charge_l1
      << ',' << result.iscm_induced_charge_l2 << ','
      << result.iscm_induced_charge_max << ',' << result.iscm_water_phase_kt
      << ',' << result.iscm_water_phase_kj_mol << ','
      << result.iscm_vacuum_phase_kt << ','
      << result.iscm_vacuum_phase_kj_mol << ','
      << result.iscm_corrected_kt << ',' << result.iscm_corrected_kj_mol
      << ',' << csv_escape(result.iscm_status) << ','
      << csv_escape(args_.parity_map_dump) << ','
      << csv_escape(args_.test_charge_grid) << ','
      << csv_escape(args_.test_surface_map) << ','
      << csv_escape(args_.test_boundary_map) << ',' << result.status << ','
      << csv_escape(args_.test_midpoint_epsilon) << ','
      << result.delphi_face_crossings << ','
      << result.delphi_eps_face_matches << ','
      << result.delphi_boundary_grid_points << ','
      << result.delphi_charged_boundary_points << ','
      << csv_escape(result.delphi_iscm_status) << ','
      << csv_escape(args_.test_vdwms_remap) << ','
      << csv_escape(result.vdwms_enabled) << ','
      << result.vdwms_boundary_initial << ',' << result.vdwms_boundary_final
      << ',' << result.vdwms_added << ',' << result.vdwms_removed << ','
      << result.vdwms_exposed_points << ',' << result.vdwms_midpoint_rms
      << ',' << audit.midpoint_eps_sum_min << ','
      << audit.midpoint_eps_sum_max << ','
      << audit.midpoint_eps_sum_mean << ','
      << audit.midpoint_eps_sum_interior_mean
      << '\n';
}

bool Analyzer::barnase_forensic_enabled() const {
  if (args_.test_pqr.empty() || args_.parity_dump == "off")
    return false;
  std::string probe = args_.test_pqr;
  try {
    probe += " ";
    probe += std::filesystem::path(args_.test_pqr).filename().string();
  } catch (...) {
  }
  std::transform(probe.begin(), probe.end(), probe.begin(), [](unsigned char c) {
    return static_cast<char>(std::tolower(c));
  });
  return probe.find("barnase") != std::string::npos;
}

void Analyzer::maybe_write_iscm_surface_charge_samples(
    const FrameResult &result, const IscmAudit &audit) const {
  if (args_.test_pqr.empty() || args_.parity_dump != "all" ||
      audit.samples.empty())
    return;

  std::ofstream out;
  if (!iscm_surface_samples_header_written_) {
    out.open("pbe_iscm_surface_charge_samples.csv",
             std::ios::out | std::ios::trunc);
    if (!out.is_open())
      throw std::runtime_error(
          "PBE: failed to open pbe_iscm_surface_charge_samples.csv.");
    out << "trajectory_id,frame_index,global_index,kind,rank,i,j,k,x_A,y_A,"
           "z_A,phi_center,phi_neighbor_average,charge_adjustment,"
           "induced_charge,atom_potential_sum,energy_contribution_kT,"
           "energy_contribution_kJmol,iscm_status\n";
    iscm_surface_samples_header_written_ = true;
  } else {
    out.open("pbe_iscm_surface_charge_samples.csv",
             std::ios::out | std::ios::app);
    if (!out.is_open())
      throw std::runtime_error(
          "PBE: failed to append pbe_iscm_surface_charge_samples.csv.");
  }

  const double kt_kj = kKBKJPerMolK * std::max(1.0, args_.temp);
  out << std::setprecision(12);
  for (const IscmSample &sample : audit.samples) {
    out << result.trajectory_id << ',' << result.frame_index << ','
        << result.global_index << ',' << csv_escape(sample.kind) << ','
        << sample.rank << ',' << sample.i << ',' << sample.j << ','
        << sample.k << ',' << sample.x << ',' << sample.y << ',' << sample.z
        << ',' << sample.phi_center << ',' << sample.phi_neighbor_average
        << ',' << sample.charge_adjustment << ',' << sample.induced_charge
        << ',' << sample.atom_potential_sum << ','
        << sample.energy_contribution_kt << ','
        << sample.energy_contribution_kt * kt_kj << ','
        << csv_escape(audit.status) << '\n';
  }
}

void Analyzer::maybe_write_barnase_forensic_audits(
    const FrameResult &result, const std::vector<Atom> &atoms,
    const FieldSolve &solvent, const FieldSolve &vacuum,
    const ReactionSolve &solve, const ForensicAudit &solvent_audit) const {
  if (!barnase_forensic_enabled())
    return;

  struct GroupStats {
    std::size_t count{0};
    double charge_sum{0.0};
    double radius_sum{0.0};
    double radius_min{std::numeric_limits<double>::max()};
    double radius_max{0.0};
    double media_radius_sum{0.0};
    double media_radius_min{std::numeric_limits<double>::max()};
    double media_radius_max{0.0};
    double x_min{std::numeric_limits<double>::max()};
    double x_max{std::numeric_limits<double>::lowest()};
    double y_min{std::numeric_limits<double>::max()};
    double y_max{std::numeric_limits<double>::lowest()};
    double z_min{std::numeric_limits<double>::max()};
    double z_max{std::numeric_limits<double>::lowest()};
    std::size_t radius_floor_count{0};
    std::size_t radius_changed_count{0};
  };

  auto update_group = [&](GroupStats &stats, const Atom &atom) {
    ++stats.count;
    const double media_radius = effective_surface_radius(atom);
    stats.charge_sum += atom.q;
    stats.radius_sum += atom.radius;
    stats.radius_min = std::min(stats.radius_min, atom.radius);
    stats.radius_max = std::max(stats.radius_max, atom.radius);
    stats.media_radius_sum += media_radius;
    stats.media_radius_min = std::min(stats.media_radius_min, media_radius);
    stats.media_radius_max = std::max(stats.media_radius_max, media_radius);
    stats.x_min = std::min(stats.x_min, atom.x);
    stats.x_max = std::max(stats.x_max, atom.x);
    stats.y_min = std::min(stats.y_min, atom.y);
    stats.y_max = std::max(stats.y_max, atom.y);
    stats.z_min = std::min(stats.z_min, atom.z);
    stats.z_max = std::max(stats.z_max, atom.z);
    if (media_radius > atom.radius + 1.0e-12)
      ++stats.radius_floor_count;
    if (std::abs(media_radius - atom.radius) > 1.0e-12)
      ++stats.radius_changed_count;
  };
  auto finish_minmax = [](double value, std::size_t count) {
    return count > 0 ? value : 0.0;
  };
  auto atom_element = [&](const Atom &atom) {
    return inferred_element_symbol(atom);
  };
  auto residue_key = [](const Atom &atom) {
    std::ostringstream oss;
    oss << (atom.chain.empty() ? "_" : atom.chain) << ':'
        << (atom.resname.empty() ? "UNK" : atom.resname) << ':'
        << atom.resnr;
    return oss.str();
  };

  GroupStats all;
  std::map<std::string, GroupStats> by_element;
  std::map<std::string, GroupStats> by_residue;
  std::size_t h_floor = 0;
  std::size_t c_floor = 0;
  std::size_t other_floor = 0;
  for (const Atom &atom : atoms) {
    update_group(all, atom);
    update_group(by_element[atom_element(atom)], atom);
    update_group(by_residue[residue_key(atom)], atom);
    if (effective_surface_radius(atom) > atom.radius + 1.0e-12) {
      const std::string element = atom_element(atom);
      if (element == "H")
        ++h_floor;
      else if (element == "C")
        ++c_floor;
      else
        ++other_floor;
    }
  }

  const GridShape &shape = solvent.shape;
  auto write_input_header = [&](std::ofstream &out) {
    out << "trajectory_id,frame_index,global_index,section,key,atom_count,"
           "total_charge,radius_sum_A,radius_min_A,radius_max_A,x_min_A,"
           "x_max_A,y_min_A,y_max_A,z_min_A,z_max_A,grid_nx,grid_ny,"
           "grid_nz,grid_origin_x_A,grid_origin_y_A,grid_origin_z_A,"
           "grid_center_x_A,grid_center_y_A,grid_center_z_A,grid_spacing_A,"
           "pqr_source,radius_convention,radius_floor_count,"
           "media_radius_sum_A,media_radius_min_A,media_radius_max_A,"
           "radius_changed_count\n";
  };
  auto write_group_row = [&](std::ofstream &out, const std::string &section,
                             const std::string &key, const GroupStats &stats) {
    out << result.trajectory_id << ',' << result.frame_index << ','
        << result.global_index << ',' << csv_escape(section) << ','
        << csv_escape(key) << ',' << stats.count << ',' << stats.charge_sum
        << ',' << stats.radius_sum << ','
        << finish_minmax(stats.radius_min, stats.count) << ','
        << finish_minmax(stats.radius_max, stats.count) << ','
        << finish_minmax(stats.x_min, stats.count) << ','
        << finish_minmax(stats.x_max, stats.count) << ','
        << finish_minmax(stats.y_min, stats.count) << ','
        << finish_minmax(stats.y_max, stats.count) << ','
        << finish_minmax(stats.z_min, stats.count) << ','
        << finish_minmax(stats.z_max, stats.count) << ',' << shape.nx << ','
        << shape.ny << ',' << shape.nz << ',' << shape.ox << ','
        << shape.oy << ',' << shape.oz << ','
        << shape.ox + 0.5 * static_cast<double>(shape.nx - 1) * shape.h
        << ',' << shape.oy + 0.5 * static_cast<double>(shape.ny - 1) * shape.h
        << ',' << shape.oz + 0.5 * static_cast<double>(shape.nz - 1) * shape.h
        << ',' << shape.h << ',' << csv_escape(args_.test_pqr) << ','
        << csv_escape(effective_radius_convention()) << ','
        << stats.radius_floor_count << ',' << stats.media_radius_sum << ','
        << finish_minmax(stats.media_radius_min, stats.count) << ','
        << finish_minmax(stats.media_radius_max, stats.count) << ','
        << stats.radius_changed_count << '\n';
  };

  std::ofstream inputs;
  if (!forensic_inputs_header_written_) {
    inputs.open("pbe_barnase_forensic_inputs.csv",
                std::ios::out | std::ios::trunc);
    if (!inputs.is_open())
      throw std::runtime_error(
          "PBE: failed to open pbe_barnase_forensic_inputs.csv.");
    write_input_header(inputs);
    forensic_inputs_header_written_ = true;
  } else {
    inputs.open("pbe_barnase_forensic_inputs.csv",
                std::ios::out | std::ios::app);
    if (!inputs.is_open())
      throw std::runtime_error(
          "PBE: failed to append pbe_barnase_forensic_inputs.csv.");
  }
  inputs << std::setprecision(12);
  write_group_row(inputs, "summary", "all", all);
  for (const auto &[key, stats] : by_element)
    write_group_row(inputs, "element", key, stats);
  for (const auto &[key, stats] : by_residue)
    write_group_row(inputs, "residue", key, stats);

  std::ofstream surface;
  if (!forensic_surface_header_written_) {
    surface.open("pbe_barnase_surface_audit.csv",
                 std::ios::out | std::ios::trunc);
    if (!surface.is_open())
      throw std::runtime_error(
          "PBE: failed to open pbe_barnase_surface_audit.csv.");
    surface
        << "trajectory_id,frame_index,global_index,surface,"
           "charge_assignment,radius_convention,solute_points,solvent_points,"
           "ion_accessible_points,stern_excluded_points,face_cross_x,"
           "face_cross_y,face_cross_z,eps_face_min,eps_face_max,"
           "eps_face_mean,eps_sum_min,eps_sum_max,eps_sum_mean,"
           "vdw_volume_estimate_A3,h_radius_floor_count,c_radius_floor_count,"
           "other_radius_floor_count,eps_face_count_x,eps_face_count_y,"
           "eps_face_count_z,internal_face_count_x,internal_face_count_y,"
            "internal_face_count_z,eps_face_mean_x,eps_face_mean_y,"
            "eps_face_mean_z,vdwms_enabled,vdwms_boundary_initial,"
            "vdwms_boundary_final,vdwms_added,vdwms_removed,"
            "vdwms_exposed_points,vdwms_midpoint_rms,"
            "midpoint_eps_sum_min,midpoint_eps_sum_max,"
            "midpoint_eps_sum_mean,midpoint_eps_sum_interior_mean\n";
    forensic_surface_header_written_ = true;
  } else {
    surface.open("pbe_barnase_surface_audit.csv",
                 std::ios::out | std::ios::app);
    if (!surface.is_open())
      throw std::runtime_error(
          "PBE: failed to append pbe_barnase_surface_audit.csv.");
  }
  surface << std::setprecision(12);
  surface << result.trajectory_id << ',' << result.frame_index << ','
          << result.global_index << ',' << csv_escape(args_.surface) << ','
          << csv_escape(args_.charge_assignment) << ','
          << csv_escape(effective_radius_convention()) << ','
          << solve.solute_points << ',' << solve.solvent_points << ','
          << solve.ion_accessible_points << ',' << solve.stern_excluded_points
          << ',' << solvent_audit.face_cross_x << ','
          << solvent_audit.face_cross_y << ',' << solvent_audit.face_cross_z
          << ',' << solvent_audit.eps_face_min << ','
          << solvent_audit.eps_face_max << ','
          << solvent_audit.eps_face_mean << ',' << solve.parity.eps_sum_min
          << ',' << solve.parity.eps_sum_max << ','
          << solve.parity.eps_sum_mean << ','
          << solvent_audit.vdw_volume_estimate_A3 << ',' << h_floor << ','
          << c_floor << ',' << other_floor << ','
          << solvent_audit.eps_face_count_x << ','
          << solvent_audit.eps_face_count_y << ','
          << solvent_audit.eps_face_count_z << ','
          << solvent_audit.internal_face_count_x << ','
          << solvent_audit.internal_face_count_y << ','
          << solvent_audit.internal_face_count_z << ','
          << solvent_audit.eps_face_mean_x << ','
          << solvent_audit.eps_face_mean_y << ','
          << solvent_audit.eps_face_mean_z << ','
          << csv_escape(solve.vdwms_enabled) << ','
          << solve.vdwms_boundary_initial << ',' << solve.vdwms_boundary_final
          << ',' << solve.vdwms_added << ',' << solve.vdwms_removed << ','
          << solve.vdwms_exposed_points << ',' << solve.vdwms_midpoint_rms
          << ',' << solve.parity.midpoint_eps_sum_min << ','
          << solve.parity.midpoint_eps_sum_max << ','
          << solve.parity.midpoint_eps_sum_mean << ','
          << solve.parity.midpoint_eps_sum_interior_mean
          << '\n';

  std::ofstream radius_media;
  if (!radius_media_audit_header_written_) {
    radius_media.open("pbe_barnase_radius_media_audit.csv",
                      std::ios::out | std::ios::trunc);
    if (!radius_media.is_open())
      throw std::runtime_error(
          "PBE: failed to open pbe_barnase_radius_media_audit.csv.");
    radius_media
        << "trajectory_id,frame_index,global_index,atom_index,atom_name,"
           "raw_element,inferred_element,resname,resnr,chain,q,raw_radius_A,"
           "media_radius_A,radius_source,radius_changed,radius_convention,"
           "pqr_source\n";
    radius_media_audit_header_written_ = true;
  } else {
    radius_media.open("pbe_barnase_radius_media_audit.csv",
                      std::ios::out | std::ios::app);
    if (!radius_media.is_open())
      throw std::runtime_error(
          "PBE: failed to append pbe_barnase_radius_media_audit.csv.");
  }
  radius_media << std::setprecision(12);
  for (std::size_t i = 0; i < atoms.size(); ++i) {
    const Atom &atom = atoms[i];
    const double media_radius = effective_surface_radius(atom);
    radius_media << result.trajectory_id << ',' << result.frame_index << ','
                 << result.global_index << ',' << (i + 1) << ','
                 << csv_escape(atom.atom_name) << ','
                 << csv_escape(atom.element.empty() ? "NA" : atom.element)
                 << ',' << csv_escape(inferred_element_symbol(atom)) << ','
                 << csv_escape(atom.resname) << ',' << atom.resnr << ','
                 << csv_escape(atom.chain) << ',' << atom.q << ','
                 << atom.radius << ',' << media_radius << ','
                 << csv_escape(radius_media_source(atom)) << ','
                 << (std::abs(media_radius - atom.radius) > 1.0e-12 ? 1 : 0)
                 << ',' << csv_escape(effective_radius_convention()) << ','
                 << csv_escape(args_.test_pqr) << '\n';
  }

  std::ofstream op;
  if (!forensic_operator_header_written_) {
    op.open("pbe_barnase_operator_audit.csv",
            std::ios::out | std::ios::trunc);
    if (!op.is_open())
      throw std::runtime_error(
          "PBE: failed to open pbe_barnase_operator_audit.csv.");
    op << "trajectory_id,frame_index,global_index,rhs_scale,rho_sum,rho_l1,"
          "rho_l2,rho_mx,rho_my,rho_mz,boundary_rms,boundary_max,"
          "linear_residual,denom_min,denom_max,denom_mean,screen_sum,"
          "screen_max,solvent_phase_rf_kJmol,vacuum_phase_rf_kJmol,"
          "atom_corrected_rxn_kJmol,grid_corrected_rxn_kJmol,"
          "coulomb_energy_kJmol,final_E_rxn_kJmol,reaction_energy_source,"
          "iscm_boundary_points,iscm_charged_boundary_points,"
          "iscm_induced_charge_sum,iscm_induced_charge_l1,"
          "iscm_induced_charge_l2,iscm_induced_charge_max,"
          "iscm_water_rf_kJmol,iscm_vacuum_rf_kJmol,"
          "iscm_corrected_rxn_kJmol,iscm_status,"
          "vdwms_enabled,vdwms_boundary_initial,vdwms_boundary_final,"
          "vdwms_added,vdwms_removed,vdwms_exposed_points,"
          "vdwms_midpoint_rms,"
          "status\n";
    forensic_operator_header_written_ = true;
  } else {
    op.open("pbe_barnase_operator_audit.csv",
            std::ios::out | std::ios::app);
    if (!op.is_open())
      throw std::runtime_error(
          "PBE: failed to append pbe_barnase_operator_audit.csv.");
  }
  op << std::setprecision(12);
  op << result.trajectory_id << ',' << result.frame_index << ','
     << result.global_index << ',' << solve.rhs_scale << ','
     << solvent_audit.rho_sum << ',' << solvent_audit.rho_l1 << ','
     << solvent_audit.rho_l2 << ',' << solvent_audit.rho_mx << ','
     << solvent_audit.rho_my << ',' << solvent_audit.rho_mz << ','
     << solve.boundary_rms << ',' << solve.boundary_max << ','
     << solve.linear_residual << ',' << solvent_audit.denom_min << ','
     << solvent_audit.denom_max << ',' << solvent_audit.denom_mean << ','
     << solvent_audit.screen_sum << ',' << solvent_audit.screen_max << ','
     << solve.solvent_phase_kj_mol << ',' << solve.vacuum_phase_kj_mol
     << ',' << solve.atom_corrected_kj_mol << ','
     << solve.grid_corrected_kj_mol << ',' << solve.coulomb_energy_kj_mol
     << ',' << solve.energy_kj_mol << ','
     << csv_escape(solve.reaction_energy_source) << ','
     << solve.iscm_boundary_points << ',' << solve.iscm_charged_boundary_points
     << ',' << solve.iscm_induced_charge_sum << ','
     << solve.iscm_induced_charge_l1 << ',' << solve.iscm_induced_charge_l2
     << ',' << solve.iscm_induced_charge_max << ','
     << solve.iscm_water_phase_kj_mol << ','
     << solve.iscm_vacuum_phase_kj_mol << ','
     << solve.iscm_corrected_kj_mol << ',' << csv_escape(solve.iscm_status)
     << ',' << csv_escape(solve.vdwms_enabled) << ','
     << solve.vdwms_boundary_initial << ',' << solve.vdwms_boundary_final
     << ',' << solve.vdwms_added << ',' << solve.vdwms_removed << ','
     << solve.vdwms_exposed_points << ',' << solve.vdwms_midpoint_rms << ','
     << (solve.converged ? "ok" : "not_converged") << '\n';

  maybe_write_atom_phi_audit(result, atoms, solvent, vacuum, solve);
  maybe_write_linear_forensic_points(result, atoms, solvent, vacuum, solve);
}

void Analyzer::maybe_write_atom_phi_audit(
    const FrameResult &result, const std::vector<Atom> &atoms,
    const FieldSolve &solvent, const FieldSolve &vacuum,
    const ReactionSolve &solve) const {
  if (!barnase_forensic_enabled() || solvent.atom_phi.size() != atoms.size() ||
      vacuum.atom_phi.size() != atoms.size())
    return;
  const GridShape &shape = solvent.shape;
  auto grid_index = [&](double value, double origin, std::size_t n) {
    if (n == 0 || shape.h <= 0.0 || !std::isfinite(value))
      return std::size_t{0};
    const double scaled = std::floor((value - origin) / shape.h);
    const double clamped =
        std::clamp(scaled, 0.0, static_cast<double>(n > 0 ? n - 1 : 0));
    return static_cast<std::size_t>(clamped);
  };

  std::ofstream out;
  if (!atom_phi_audit_header_written_) {
    out.open("pbe_atom_phi_audit.csv", std::ios::out | std::ios::trunc);
    if (!out.is_open())
      throw std::runtime_error("PBE: failed to open pbe_atom_phi_audit.csv.");
    out << "trajectory_id,frame_index,global_index,atom_index,atom_name,"
           "element,resname,resnr,chain,x_A,y_A,z_A,q,radius_A,grid_i,"
           "grid_j,grid_k,solvent_phi,vacuum_phi,reaction_phi,"
           "contribution_kJmol\n";
    atom_phi_audit_header_written_ = true;
  } else {
    out.open("pbe_atom_phi_audit.csv", std::ios::out | std::ios::app);
    if (!out.is_open())
      throw std::runtime_error("PBE: failed to append pbe_atom_phi_audit.csv.");
  }
  out << std::setprecision(12);

  for (std::size_t i = 0; i < atoms.size(); ++i) {
    const Atom &atom = atoms[i];
    const double solvent_phi = solvent.atom_phi[i];
    const double vacuum_phi = vacuum.atom_phi[i];
    const double contribution =
        i < solve.atom_reaction_kj_mol.size() ? solve.atom_reaction_kj_mol[i]
                                              : 0.0;
    out << result.trajectory_id << ',' << result.frame_index << ','
        << result.global_index << ',' << (i + 1) << ','
        << csv_escape(atom.atom_name) << ','
        << csv_escape(inferred_element_symbol(atom))
        << ',' << csv_escape(atom.resname) << ',' << atom.resnr << ','
        << csv_escape(atom.chain) << ',' << atom.x << ',' << atom.y << ','
        << atom.z << ',' << atom.q << ',' << atom.radius << ','
        << grid_index(atom.x, shape.ox, shape.nx) << ','
        << grid_index(atom.y, shape.oy, shape.ny) << ','
        << grid_index(atom.z, shape.oz, shape.nz) << ',' << solvent_phi
        << ',' << vacuum_phi << ',' << (solvent_phi - vacuum_phi) << ','
        << contribution << '\n';
  }
}

void Analyzer::maybe_write_linear_forensic_points(
    const FrameResult &result, const std::vector<Atom> &atoms,
    const FieldSolve &solvent, const FieldSolve &vacuum,
    const ReactionSolve &solve) const {
  if (!barnase_forensic_enabled() ||
      (args_.parity_dump != "points" && args_.parity_dump != "all") ||
      solvent.shape.n == 0 || solvent.phi_grid.size() != solvent.shape.n ||
      vacuum.phi_grid.size() != solvent.shape.n)
    return;
  const GridShape &shape = solvent.shape;
  auto point_id = [&](double x, double y, double z) {
    auto to_index = [&](double value, double origin, std::size_t n) {
      if (n <= 2 || shape.h <= 0.0 || !std::isfinite(value))
        return std::size_t{1};
      const double scaled = std::floor((value - origin) / shape.h);
      return static_cast<std::size_t>(
          std::clamp(scaled, 1.0, static_cast<double>(n - 2)));
    };
    return idx3(to_index(x, shape.ox, shape.nx), to_index(y, shape.oy, shape.ny),
                to_index(z, shape.oz, shape.nz), shape.nx, shape.ny);
  };
  std::vector<std::pair<std::string, std::size_t>> rows;
  std::vector<std::pair<double, std::size_t>> ranked_atoms;
  ranked_atoms.reserve(atoms.size());
  for (std::size_t i = 0; i < atoms.size(); ++i) {
    const double contribution =
        i < solve.atom_reaction_kj_mol.size()
            ? std::abs(solve.atom_reaction_kj_mol[i])
            : std::abs(atoms[i].q);
    ranked_atoms.push_back({contribution, i});
  }
  std::sort(ranked_atoms.begin(), ranked_atoms.end(),
            [](const auto &a, const auto &b) { return a.first > b.first; });
  for (std::size_t rank = 0; rank < std::min<std::size_t>(16, ranked_atoms.size());
       ++rank) {
    const Atom &atom = atoms[ranked_atoms[rank].second];
    rows.push_back({"top_atom_energy_" + std::to_string(rank + 1),
                    point_id(atom.x, atom.y, atom.z)});
  }
  std::vector<std::pair<double, std::size_t>> top_rho;
  top_rho.reserve(16);
  auto consider_rho = [&](std::size_t id) {
    const double value =
        rho_.size() == shape.n ? std::abs(rho_[id]) : 0.0;
    if (value <= 0.0)
      return;
    if (top_rho.size() < 16) {
      top_rho.push_back({value, id});
      return;
    }
    auto min_it = std::min_element(top_rho.begin(), top_rho.end());
    if (min_it != top_rho.end() && value > min_it->first)
      *min_it = {value, id};
  };
  if (rho_.size() == shape.n) {
    for (std::size_t id = 0; id < shape.n; ++id)
      consider_rho(id);
  }
  std::sort(top_rho.begin(), top_rho.end(),
            [](const auto &a, const auto &b) { return a.first > b.first; });
  for (std::size_t rank = 0; rank < top_rho.size(); ++rank)
    rows.push_back({"top_rho_" + std::to_string(rank + 1),
                    top_rho[rank].second});

  auto add_surface = [&](const std::string &label, std::size_t id) {
    if (rows.size() >= 48)
      return;
    rows.push_back({label, id});
  };
  if (inside_.size() == shape.n) {
    const std::size_t xy = shape.nx * shape.ny;
    std::size_t added = 0;
    for (std::size_t k = 1; k + 1 < shape.nz && added < 16; ++k) {
      for (std::size_t j = 1; j + 1 < shape.ny && added < 16; ++j) {
        for (std::size_t i = 1; i + 1 < shape.nx && added < 16; ++i) {
          const std::size_t id = idx3(i, j, k, shape.nx, shape.ny);
          const bool crossing =
              inside_[id] != inside_[id + 1] ||
              inside_[id] != inside_[id + shape.nx] ||
              inside_[id] != inside_[id + xy];
          if (crossing) {
            add_surface("surface_cross_" + std::to_string(++added), id);
          }
        }
      }
    }
  }
  if (rows.empty())
    return;

  std::ofstream out;
  if (!parity_points_header_written_) {
    out.open("pbe_parity_points.csv", std::ios::out | std::ios::trunc);
    if (!out.is_open())
      throw std::runtime_error("PBE: failed to open pbe_parity_points.csv.");
    out << "trajectory_id,frame_index,global_index,point_kind,rank,i,j,k,"
           "x_A,y_A,z_A,rho,eps_sum,screen,phi,source,rhs,F,Jdiag,Jdelta,"
           "residual\n";
    parity_points_header_written_ = true;
  } else {
    out.open("pbe_parity_points.csv", std::ios::out | std::ios::app);
    if (!out.is_open())
      throw std::runtime_error("PBE: failed to append pbe_parity_points.csv.");
  }
  out << std::setprecision(12);

  const double rhs = rhs_scale(shape);
  const std::size_t xy = shape.nx * shape.ny;
  std::size_t rank = 0;
  for (const auto &[kind, id] : rows) {
    if (id >= shape.n)
      continue;
    const std::size_t k = id / xy;
    const std::size_t rem = id - k * xy;
    const std::size_t j = rem / shape.nx;
    const std::size_t i = rem - j * shape.nx;
    const double eps_sum =
        denom_.size() == shape.n && screen_.size() == shape.n
            ? std::max(1.0e-12, denom_[id] - screen_[id])
            : 0.0;
    const double screen = screen_.size() == shape.n ? screen_[id] : 0.0;
    const double rho = rho_.size() == shape.n ? rho_[id] : 0.0;
    const double solvent_phi = solvent.phi_grid[id];
    const double vacuum_phi = vacuum.phi_grid[id];
    out << result.trajectory_id << ',' << result.frame_index << ','
        << result.global_index << ',' << csv_escape(kind) << ',' << ++rank
        << ',' << i << ',' << j << ',' << k << ','
        << shape.ox + static_cast<double>(i) * shape.h << ','
        << shape.oy + static_cast<double>(j) * shape.h << ','
        << shape.oz + static_cast<double>(k) * shape.h << ',' << rho << ','
        << eps_sum << ',' << screen << ',' << solvent_phi << ','
        << vacuum_phi << ',' << rhs * rho << ','
        << (solvent_phi - vacuum_phi) << ',' << eps_sum << ",0,0\n";
  }
}

bool Analyzer::parity_map_dump_enabled() const noexcept {
  return !args_.test_pqr.empty() && args_.parity_map_dump != "off";
}

std::vector<Analyzer::MapSamplePoint>
Analyzer::collect_map_sample_points(const GridShape &shape,
                                    const std::vector<Atom> &atoms,
                                    const ReactionSolve &solve) const {
  std::vector<MapSamplePoint> rows;
  rows.reserve(128);
  std::vector<unsigned char> used(shape.n, 0U);
  auto to_index = [&](double value, double origin, std::size_t n) {
    if (n <= 2 || shape.h <= 0.0 || !std::isfinite(value))
      return std::size_t{1};
    const double scaled = std::floor((value - origin) / shape.h);
    return static_cast<std::size_t>(
        std::clamp(scaled, 1.0, static_cast<double>(n - 2)));
  };
  auto add = [&](std::string kind, std::size_t id, double hint) {
    if (id >= shape.n || used[id])
      return;
    used[id] = 1U;
    rows.push_back({std::move(kind), rows.size() + 1, id, hint});
  };
  auto atom_id = [&](const Atom &atom) {
    return idx3(to_index(atom.x, shape.ox, shape.nx),
                to_index(atom.y, shape.oy, shape.ny),
                to_index(atom.z, shape.oz, shape.nz), shape.nx, shape.ny);
  };

  add("center", idx3(shape.nx / 2, shape.ny / 2, shape.nz / 2, shape.nx,
                     shape.ny),
      0.0);
  add("quartile_low",
      idx3(std::max<std::size_t>(1, shape.nx / 4),
           std::max<std::size_t>(1, shape.ny / 4),
           std::max<std::size_t>(1, shape.nz / 4), shape.nx, shape.ny),
      0.0);
  add("quartile_high",
      idx3(std::min<std::size_t>(shape.nx - 2, (shape.nx * 3) / 4),
           std::min<std::size_t>(shape.ny - 2, (shape.ny * 3) / 4),
           std::min<std::size_t>(shape.nz - 2, (shape.nz * 3) / 4), shape.nx,
           shape.ny),
      0.0);

  std::vector<std::pair<double, std::size_t>> ranked_atoms;
  ranked_atoms.reserve(atoms.size());
  for (std::size_t i = 0; i < atoms.size(); ++i) {
    const double hint =
        i < solve.atom_reaction_kj_mol.size()
            ? std::abs(solve.atom_reaction_kj_mol[i])
            : std::abs(atoms[i].q);
    ranked_atoms.push_back({hint, i});
  }
  std::sort(ranked_atoms.begin(), ranked_atoms.end(),
            [](const auto &a, const auto &b) { return a.first > b.first; });
  for (std::size_t rank = 0;
       rank < std::min<std::size_t>(32, ranked_atoms.size()); ++rank) {
    const std::size_t atom_index = ranked_atoms[rank].second;
    add("charge_site_" + std::to_string(rank + 1), atom_id(atoms[atom_index]),
        ranked_atoms[rank].first);
  }

  std::vector<std::pair<double, std::size_t>> top_rho;
  auto consider_top = [&](std::vector<std::pair<double, std::size_t>> &top,
                          double value, std::size_t id,
                          std::size_t limit) {
    if (value <= 0.0)
      return;
    if (top.size() < limit) {
      top.push_back({value, id});
      return;
    }
    auto min_it = std::min_element(top.begin(), top.end());
    if (min_it != top.end() && value > min_it->first)
      *min_it = {value, id};
  };
  if (rho_.size() == shape.n) {
    for (std::size_t id = 0; id < shape.n; ++id)
      consider_top(top_rho, std::abs(rho_[id]), id, 32);
  }
  std::sort(top_rho.begin(), top_rho.end(),
            [](const auto &a, const auto &b) { return a.first > b.first; });
  for (std::size_t rank = 0; rank < top_rho.size(); ++rank)
    add("top_rho_" + std::to_string(rank + 1), top_rho[rank].second,
        top_rho[rank].first);

  if (inside_.size() == shape.n) {
    const std::size_t xy = shape.nx * shape.ny;
    std::size_t added = 0;
    for (std::size_t k = 1; k + 1 < shape.nz && added < 32; ++k) {
      for (std::size_t j = 1; j + 1 < shape.ny && added < 32; ++j) {
        for (std::size_t i = 1; i + 1 < shape.nx && added < 32; ++i) {
          const std::size_t id = idx3(i, j, k, shape.nx, shape.ny);
          const bool crossing =
              inside_[id] != inside_[id + 1] ||
              inside_[id] != inside_[id + shape.nx] ||
              inside_[id] != inside_[id + xy] ||
              inside_[id] != inside_[id - 1] ||
              inside_[id] != inside_[id - shape.nx] ||
              inside_[id] != inside_[id - xy];
          if (crossing)
            add("surface_cross_" + std::to_string(++added), id, 0.0);
        }
      }
    }
  }
  for (std::size_t i = 0; i < rows.size(); ++i)
    rows[i].rank = i + 1;
  return rows;
}

void Analyzer::maybe_write_parity_map_dumps(
    const FrameResult &result, const std::vector<Atom> &atoms,
    const FieldSolve &solvent, const FieldSolve &vacuum,
    const ReactionSolve &solve) const {
  if (!parity_map_dump_enabled() || solvent.shape.n == 0 ||
      solvent.phi_grid.size() != solvent.shape.n ||
      vacuum.phi_grid.size() != solvent.shape.n)
    return;
  const GridShape &shape = solvent.shape;
  const bool full = args_.parity_map_dump == "full";
  const auto samples = collect_map_sample_points(shape, atoms, solve);
  auto coords = [&](std::size_t id, std::size_t &i, std::size_t &j,
                    std::size_t &k) {
    const std::size_t xy = shape.nx * shape.ny;
    k = id / xy;
    const std::size_t rem = id - k * xy;
    j = rem / shape.nx;
    i = rem - j * shape.nx;
  };
  auto eps_sum_at = [&](std::size_t id) {
    if (solvent.eps_sum_grid.size() == shape.n)
      return solvent.eps_sum_grid[id];
    return denom_.size() == shape.n && screen_.size() == shape.n
               ? std::max(1.0e-12, denom_[id] - screen_[id])
               : 0.0;
  };
  auto screen_at = [&](std::size_t id) {
    if (solvent.screen_grid.size() == shape.n)
      return solvent.screen_grid[id];
    return screen_.size() == shape.n ? screen_[id] : 0.0;
  };
  auto rho_at = [&](std::size_t id) {
    if (solvent.rho_grid.size() == shape.n)
      return solvent.rho_grid[id];
    return rho_.size() == shape.n ? rho_[id] : 0.0;
  };
  auto inside_at = [&](std::size_t id) {
    if (solvent.inside_grid.size() == shape.n)
      return solvent.inside_grid[id] ? 1 : 0;
    return inside_.size() == shape.n && inside_[id] ? 1 : 0;
  };
  auto eps_x_at = [&](std::size_t id) {
    if (solvent.eps_x_grid.size() == shape.n)
      return solvent.eps_x_grid[id];
    return eps_x_.size() == shape.n ? eps_x_[id] : 0.0;
  };
  auto eps_y_at = [&](std::size_t id) {
    if (solvent.eps_y_grid.size() == shape.n)
      return solvent.eps_y_grid[id];
    return eps_y_.size() == shape.n ? eps_y_[id] : 0.0;
  };
  auto eps_z_at = [&](std::size_t id) {
    if (solvent.eps_z_grid.size() == shape.n)
      return solvent.eps_z_grid[id];
    return eps_z_.size() == shape.n ? eps_z_[id] : 0.0;
  };

  {
    std::ofstream out;
    if (!field_samples_header_written_) {
      out.open("pbe_field_samples.csv", std::ios::out | std::ios::trunc);
      if (!out.is_open())
        throw std::runtime_error("PBE: failed to open pbe_field_samples.csv.");
      out << "trajectory_id,frame_index,global_index,point_kind,rank,i,j,k,"
             "x_A,y_A,z_A,rho,inside,eps_sum,eps_x,eps_y,eps_z,screen,solvent_phi,"
             "vacuum_phi,reaction_phi,atom_energy_hint_kJmol,is_boundary\n";
      field_samples_header_written_ = true;
    } else {
      out.open("pbe_field_samples.csv", std::ios::out | std::ios::app);
      if (!out.is_open())
        throw std::runtime_error("PBE: failed to append pbe_field_samples.csv.");
    }
    out << std::setprecision(12);
    auto write_field = [&](const std::string &kind, std::size_t rank,
                           std::size_t id, double hint) {
      if (id >= shape.n)
        return;
      std::size_t i = 0, j = 0, k = 0;
      coords(id, i, j, k);
      out << result.trajectory_id << ',' << result.frame_index << ','
          << result.global_index << ',' << csv_escape(kind) << ',' << rank
          << ',' << i << ',' << j << ',' << k << ','
          << shape.ox + static_cast<double>(i) * shape.h << ','
          << shape.oy + static_cast<double>(j) * shape.h << ','
          << shape.oz + static_cast<double>(k) * shape.h << ',' << rho_at(id)
          << ',' << inside_at(id) << ',' << eps_sum_at(id) << ','
          << eps_x_at(id) << ',' << eps_y_at(id) << ',' << eps_z_at(id) << ','
          << screen_at(id) << ',' << solvent.phi_grid[id] << ','
          << vacuum.phi_grid[id] << ','
          << solvent.phi_grid[id] - vacuum.phi_grid[id] << ',' << hint << ','
          << (is_boundary(i, j, k, shape) ? 1 : 0) << '\n';
    };
    if (full) {
      for (std::size_t id = 0; id < shape.n; ++id)
        write_field("full_grid", id + 1, id, 0.0);
    } else {
      for (const MapSamplePoint &sample : samples)
        write_field(sample.kind, sample.rank, sample.id, sample.hint);
    }
  }

  {
    std::ofstream out;
    if (!boundary_samples_header_written_) {
      out.open("pbe_boundary_samples.csv", std::ios::out | std::ios::trunc);
      if (!out.is_open())
        throw std::runtime_error(
            "PBE: failed to open pbe_boundary_samples.csv.");
      out << "trajectory_id,frame_index,global_index,boundary_kind,rank,i,j,k,"
             "x_A,y_A,z_A,solvent_phi,vacuum_phi,rho,inside,eps_sum,"
             "eps_x,eps_y,eps_z,screen\n";
      boundary_samples_header_written_ = true;
    } else {
      out.open("pbe_boundary_samples.csv", std::ios::out | std::ios::app);
      if (!out.is_open())
        throw std::runtime_error(
            "PBE: failed to append pbe_boundary_samples.csv.");
    }
    out << std::setprecision(12);
    auto write_boundary = [&](const std::string &kind, std::size_t rank,
                              std::size_t id) {
      std::size_t i = 0, j = 0, k = 0;
      coords(id, i, j, k);
      if (!is_boundary(i, j, k, shape))
        return;
      out << result.trajectory_id << ',' << result.frame_index << ','
          << result.global_index << ',' << csv_escape(kind) << ',' << rank
          << ',' << i << ',' << j << ',' << k << ','
          << shape.ox + static_cast<double>(i) * shape.h << ','
          << shape.oy + static_cast<double>(j) * shape.h << ','
          << shape.oz + static_cast<double>(k) * shape.h << ','
          << solvent.phi_grid[id] << ',' << vacuum.phi_grid[id] << ','
          << rho_at(id) << ',' << inside_at(id) << ',' << eps_sum_at(id)
          << ',' << eps_x_at(id) << ',' << eps_y_at(id) << ','
          << eps_z_at(id) << ',' << screen_at(id) << '\n';
    };
    if (full) {
      std::size_t rank = 0;
      for (std::size_t k = 0; k < shape.nz; ++k) {
        for (std::size_t j = 0; j < shape.ny; ++j) {
          for (std::size_t i = 0; i < shape.nx; ++i) {
            if (is_boundary(i, j, k, shape))
              write_boundary("full_boundary", ++rank,
                             idx3(i, j, k, shape.nx, shape.ny));
          }
        }
      }
    } else {
      std::vector<std::size_t> ids = {
          idx3(0, 0, 0, shape.nx, shape.ny),
          idx3(shape.nx - 1, 0, 0, shape.nx, shape.ny),
          idx3(0, shape.ny - 1, 0, shape.nx, shape.ny),
          idx3(0, 0, shape.nz - 1, shape.nx, shape.ny),
          idx3(shape.nx / 2, shape.ny / 2, 0, shape.nx, shape.ny),
          idx3(shape.nx / 2, shape.ny / 2, shape.nz - 1, shape.nx,
               shape.ny),
          idx3(shape.nx / 2, 0, shape.nz / 2, shape.nx, shape.ny),
          idx3(shape.nx / 2, shape.ny - 1, shape.nz / 2, shape.nx,
               shape.ny),
          idx3(0, shape.ny / 2, shape.nz / 2, shape.nx, shape.ny),
          idx3(shape.nx - 1, shape.ny / 2, shape.nz / 2, shape.nx,
               shape.ny)};
      std::sort(ids.begin(), ids.end());
      ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
      for (std::size_t rank = 0; rank < ids.size(); ++rank)
        write_boundary("sample_boundary", rank + 1, ids[rank]);
    }
  }

  {
    std::ofstream out;
    if (!charge_grid_audit_header_written_) {
      out.open("pbe_charge_grid_audit.csv", std::ios::out | std::ios::trunc);
      if (!out.is_open())
        throw std::runtime_error(
            "PBE: failed to open pbe_charge_grid_audit.csv.");
      out << "trajectory_id,frame_index,global_index,charge_kind,rank,i,j,k,"
             "x_A,y_A,z_A,rho,abs_rho,inside,eps_sum,eps_x,eps_y,eps_z,screen\n";
      charge_grid_audit_header_written_ = true;
    } else {
      out.open("pbe_charge_grid_audit.csv", std::ios::out | std::ios::app);
      if (!out.is_open())
        throw std::runtime_error(
            "PBE: failed to append pbe_charge_grid_audit.csv.");
    }
    out << std::setprecision(12);
    auto write_charge = [&](const std::string &kind, std::size_t rank,
                            std::size_t id) {
      const double rho = rho_at(id);
      if (rho == 0.0 && full)
        return;
      std::size_t i = 0, j = 0, k = 0;
      coords(id, i, j, k);
      out << result.trajectory_id << ',' << result.frame_index << ','
          << result.global_index << ',' << csv_escape(kind) << ',' << rank
          << ',' << i << ',' << j << ',' << k << ','
          << shape.ox + static_cast<double>(i) * shape.h << ','
          << shape.oy + static_cast<double>(j) * shape.h << ','
          << shape.oz + static_cast<double>(k) * shape.h << ',' << rho << ','
          << std::abs(rho) << ',' << inside_at(id) << ',' << eps_sum_at(id)
          << ',' << eps_x_at(id) << ',' << eps_y_at(id) << ','
          << eps_z_at(id) << ',' << screen_at(id) << '\n';
    };
    if (full) {
      std::size_t rank = 0;
      for (std::size_t id = 0; id < shape.n; ++id) {
        if (rho_at(id) != 0.0)
          write_charge("full_nonzero_rho", ++rank, id);
      }
    } else {
      std::size_t rank = 0;
      for (const MapSamplePoint &sample : samples) {
        if (std::abs(rho_at(sample.id)) > 0.0)
          write_charge("sample_rho", ++rank, sample.id);
      }
    }
  }
}

void Analyzer::maybe_write_parity_points(
    const FrameResult &result, const GridShape &shape,
    const AlignedDoubleVector &source, const AlignedDoubleVector &residual,
    const AlignedDoubleVector &delta, double alpha, double clip,
    double screen_scale) const {
  if (args_.test_pqr.empty() ||
      (args_.parity_dump != "points" && args_.parity_dump != "all") ||
      shape.n == 0 || source.size() != shape.n || residual.size() != shape.n ||
      delta.size() != shape.n || phi_.size() != shape.n)
    return;

  std::vector<std::pair<double, std::size_t>> top;
  top.reserve(16);
  auto consider_top = [&](std::size_t id) {
    const double value = std::abs(residual[id]);
    if (top.size() < 16) {
      top.push_back({value, id});
      return;
    }
    auto min_it = std::min_element(top.begin(), top.end());
    if (min_it != top.end() && value > min_it->first)
      *min_it = {value, id};
  };
  for (std::size_t k = 1; k + 1 < shape.nz; ++k) {
    for (std::size_t j = 1; j + 1 < shape.ny; ++j) {
      for (std::size_t i = 1; i + 1 < shape.nx; ++i)
        consider_top(idx3(i, j, k, shape.nx, shape.ny));
    }
  }
  std::sort(top.begin(), top.end(),
            [](const auto &a, const auto &b) { return a.first > b.first; });

  std::vector<std::pair<std::string, std::size_t>> rows;
  for (std::size_t rank = 0; rank < top.size(); ++rank)
    rows.push_back({"top_residual_" + std::to_string(rank + 1),
                    top[rank].second});
  auto add_fixed = [&](const std::string &label, std::size_t i, std::size_t j,
                       std::size_t k) {
    i = std::clamp<std::size_t>(i, 1, shape.nx > 1 ? shape.nx - 2 : 1);
    j = std::clamp<std::size_t>(j, 1, shape.ny > 1 ? shape.ny - 2 : 1);
    k = std::clamp<std::size_t>(k, 1, shape.nz > 1 ? shape.nz - 2 : 1);
    const std::size_t id = idx3(i, j, k, shape.nx, shape.ny);
    const auto duplicate = std::find_if(
        rows.begin(), rows.end(),
        [&](const auto &row) { return row.first == label || row.second == id; });
    if (duplicate == rows.end())
      rows.push_back({label, id});
  };
  add_fixed("center", shape.nx / 2, shape.ny / 2, shape.nz / 2);
  add_fixed("quartile_low", shape.nx / 4, shape.ny / 4, shape.nz / 4);
  add_fixed("quartile_high", (shape.nx * 3) / 4, (shape.ny * 3) / 4,
            (shape.nz * 3) / 4);

  std::ofstream out;
  if (!parity_points_header_written_) {
    out.open("pbe_parity_points.csv", std::ios::out | std::ios::trunc);
    if (!out.is_open())
      throw std::runtime_error("PBE: failed to open pbe_parity_points.csv.");
    out << "trajectory_id,frame_index,global_index,point_kind,rank,i,j,k,"
           "x_A,y_A,z_A,rho,eps_sum,screen,phi,source,rhs,F,Jdiag,Jdelta,"
           "residual\n";
    parity_points_header_written_ = true;
  } else {
    out.open("pbe_parity_points.csv", std::ios::out | std::ios::app);
    if (!out.is_open())
      throw std::runtime_error("PBE: failed to append pbe_parity_points.csv.");
  }

  const double rhs = nonlinear_rhs_scale(shape, alpha);
  const std::size_t xy = shape.nx * shape.ny;
  std::size_t rank = 0;
  for (const auto &[kind, id] : rows) {
    const std::size_t k = id / xy;
    const std::size_t rem = id - k * xy;
    const std::size_t j = rem / shape.nx;
    const std::size_t i = rem - j * shape.nx;
    const double eps_sum = std::max(1.0e-12, denom_[id] - screen_[id]);
    const double neighbor_sum =
        eps_x_[id] * phi_[id + 1] + eps_x_[id - 1] * phi_[id - 1] +
        eps_y_[id] * phi_[id + shape.nx] +
        eps_y_[id - shape.nx] * phi_[id - shape.nx] +
        eps_z_[id] * phi_[id + xy] + eps_z_[id - xy] * phi_[id - xy];
    const double f_value =
        eps_sum * phi_[id] - neighbor_sum +
        (screen_[id] * screen_scale) * source[id] - rhs * rho_[id];
    const double jdiag =
        std::max(1.0e-12,
                 eps_sum + (screen_[id] * screen_scale) *
                               nonlinear_source_derivative(phi_[id], alpha,
                                                           clip));
    const double neighbor_delta_sum =
        eps_x_[id] * delta[id + 1] + eps_x_[id - 1] * delta[id - 1] +
        eps_y_[id] * delta[id + shape.nx] +
        eps_y_[id - shape.nx] * delta[id - shape.nx] +
        eps_z_[id] * delta[id + xy] + eps_z_[id - xy] * delta[id - xy];
    const double jdelta = jdiag * delta[id] - neighbor_delta_sum;
    out << result.trajectory_id << ',' << result.frame_index << ','
        << result.global_index << ',' << csv_escape(kind) << ',' << ++rank
        << ',' << i << ',' << j << ',' << k << ','
        << shape.ox + static_cast<double>(i) * shape.h << ','
        << shape.oy + static_cast<double>(j) * shape.h << ','
        << shape.oz + static_cast<double>(k) * shape.h << ',' << rho_[id]
        << ',' << eps_sum << ',' << screen_[id] << ',' << phi_[id] << ','
        << source[id] << ',' << rhs * rho_[id] << ',' << f_value << ','
        << jdiag << ',' << jdelta << ',' << residual[id] << '\n';
  }
}

bool Analyzer::should_write_map() const {
  return args_.map != "off" && args_.map_frames != "none";
}

std::string Analyzer::map_suffix(const FrameResult &result) const {
  if (args_.map_frames == "last")
    return "_last";
  return "_traj" + format_suffix(result.trajectory_id, 4) + "_frame" +
         format_suffix(result.frame_index, 6);
}

void Analyzer::maybe_write_maps(const GridShape &shape,
                                const FrameResult &result,
                                bool write_phi_now) {
  if (!should_write_map())
    return;
  if (args_.map_frames == "first" && wrote_first_map_)
    return;

  const std::string suffix = map_suffix(result);
  if (write_phi_now && (args_.map == "phi" || args_.map == "all")) {
    write_cube("pbe_phi" + suffix + ".cube", shape, "phi", false);
  }
  if (!write_phi_now && (args_.map == "epsilon" || args_.map == "all")) {
    write_cube("pbe_epsilon" + suffix + ".cube", shape, "epsilon", true);
  }
  if (!write_phi_now && (args_.map == "ion" || args_.map == "all")) {
    write_cube("pbe_ion" + suffix + ".cube", shape, "ion_accessible", false,
               true);
  }
  const bool completed_first =
      (write_phi_now && (args_.map == "phi" || args_.map == "all")) ||
      (!write_phi_now && (args_.map == "epsilon" || args_.map == "ion"));
  if (args_.map_frames == "first" && completed_first)
    wrote_first_map_ = true;
}

std::string Analyzer::maybe_write_probe_sites(const FrameResult &result,
                                              const std::vector<Atom> &atoms,
                                              const ReactionSolve &solve) {
  if (args_.probe_residue.empty() && args_.probe_site.empty() &&
      args_.focus != "delphi")
    return "off";
  auto uppercase_alnum = [](const std::string &value) {
    std::string out;
    out.reserve(value.size());
    for (char ch : value) {
      if (std::isalnum(static_cast<unsigned char>(ch)))
        out.push_back(static_cast<char>(
            std::toupper(static_cast<unsigned char>(ch))));
    }
    return out;
  };
  auto digits_only = [](const std::string &value) {
    std::string out;
    for (char ch : value) {
      if (std::isdigit(static_cast<unsigned char>(ch)))
        out.push_back(ch);
    }
    return out;
  };
  auto site_tokens = [&]() {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream input(args_.probe_site);
    while (std::getline(input, token, ',')) {
      const std::string cleaned = uppercase_alnum(token);
      if (!cleaned.empty())
        tokens.push_back(cleaned);
    }
    return tokens;
  }();
  const std::string requested_residue = uppercase_alnum(args_.probe_residue);
  const std::string requested_residue_name = uppercase_letters(args_.probe_residue);
  const std::string requested_residue_number = digits_only(args_.probe_residue);
  auto residue_matches = [&](const Atom &atom) {
    if (requested_residue.empty())
      return true;
    const std::string atom_resname = uppercase_letters(atom.resname);
    const std::string atom_resid = std::to_string(atom.resnr);
    const std::string atom_key = atom_resname + atom_resid;
    if (requested_residue == atom_key || requested_residue == atom_resname)
      return true;
    if (!requested_residue_name.empty() &&
        requested_residue_name == atom_resname) {
      if (requested_residue_number.empty() ||
          requested_residue_number == atom_resid)
        return true;
      // DelPhi probe fixture labels can carry source-PDB numbering after
      // conversion to compact PQR residue indices; keep a resname fallback but
      // retain the exact key checks above for native numbering.
      return true;
    }
    return false;
  };
  auto site_matches = [&](const Atom &atom) {
    if (site_tokens.empty())
      return true;
    const std::string atom_name = uppercase_alnum(atom.atom_name);
    const std::string element = uppercase_letters(atom.element);
    for (const std::string &token : site_tokens) {
      if (token == "A" || token == "F")
        return true;
      if (token == "P") {
        if (element == "N" || element == "O" || element == "S" ||
            element == "P")
          return true;
      } else if (token == "Q") {
        if (std::abs(atom.q) >= 0.25)
          return true;
      } else if (atom_name == token || atom_name.find(token) != std::string::npos) {
        return true;
      }
    }
    return false;
  };
  double energy_kj_mol = 0.0;
  std::size_t matched = 0;
  for (std::size_t i = 0; i < atoms.size(); ++i) {
    const Atom &atom = atoms[i];
    if (!residue_matches(atom) || !site_matches(atom))
      continue;
    if (i < solve.atom_reaction_kj_mol.size())
      energy_kj_mol += solve.atom_reaction_kj_mol[i];
    ++matched;
  }
  const std::string status = matched > 0 ? "ok" : "no_matching_atoms";
  const std::string filename = "pbe_probe_sites.csv";
  const bool need_header =
      !std::filesystem::exists(filename) ||
      std::filesystem::file_size(filename) == 0;
  std::ofstream out(filename, std::ios::out | std::ios::app);
  if (!out.is_open())
    return "write_failed";
  if (need_header) {
    out << "trajectory_id,frame_index,global_index,time_ps,probe_residue,"
           "probe_site,matched_atoms,E_probe_kJmol,E_probe_kcalmol,"
           "E_probe_kT,status,solver_mode,focus_parent_grid,focus_child_grid\n";
  }
  const double kt = kKBKJPerMolK * std::max(1.0, args_.temp);
  out << result.trajectory_id << ',' << result.frame_index << ','
      << result.global_index << ',' << std::setprecision(10) << result.time_ps
      << ',' << csv_escape(args_.probe_residue) << ','
      << csv_escape(args_.probe_site) << ',' << matched << ','
      << energy_kj_mol << ',' << energy_kj_mol * kKJToKcal << ','
      << (kt > 0.0 ? energy_kj_mol / kt : 0.0) << ',' << status << ','
      << result.solver_mode << ',' << csv_escape(result.focus_parent_grid)
      << ',' << csv_escape(result.focus_child_grid) << '\n';
  return status;
}

void Analyzer::write_cube(const std::string &filename, const GridShape &shape,
                          const std::string &field_name,
                          bool epsilon_field, bool ion_field) const {
  std::ofstream out(filename, std::ios::out | std::ios::trunc);
  if (!out.is_open()) {
    TitanUI::print_warning("PBE: failed to write " + filename);
    return;
  }
  out << "TITAN PBE " << field_name << "\n";
  out << "units Angstrom; values are TITAN internal PBE units\n";
  out << std::setw(5) << 0 << " " << std::setw(12) << std::fixed
      << std::setprecision(6) << shape.ox << " " << std::setw(12) << shape.oy
      << " " << std::setw(12) << shape.oz << "\n";
  out << std::setw(5) << shape.nx << " " << std::setw(12) << shape.h << " "
      << std::setw(12) << 0.0 << " " << std::setw(12) << 0.0 << "\n";
  out << std::setw(5) << shape.ny << " " << std::setw(12) << 0.0 << " "
      << std::setw(12) << shape.h << " " << std::setw(12) << 0.0 << "\n";
  out << std::setw(5) << shape.nz << " " << std::setw(12) << 0.0 << " "
      << std::setw(12) << 0.0 << " " << std::setw(12) << shape.h << "\n";

  out << std::scientific << std::setprecision(6);
  int per_line = 0;
  for (std::size_t k = 0; k < shape.nz; ++k) {
    for (std::size_t j = 0; j < shape.ny; ++j) {
      for (std::size_t i = 0; i < shape.nx; ++i) {
        const std::size_t id = idx3(i, j, k, shape.nx, shape.ny);
        const double value =
            ion_field
                ? (ion_accessible_[id] ? 1.0 : 0.0)
                : (epsilon_field ? (inside_[id] ? args_.int_diel
                                                : args_.ext_diel)
                                 : phi_[id]);
        out << std::setw(13) << value;
        if (++per_line == 6) {
          out << "\n";
          per_line = 0;
        }
      }
    }
  }
  if (per_line != 0)
    out << "\n";
}

Analyzer::ReactionSolve Analyzer::solve_reaction_field(
    const std::vector<Atom> &atoms, TitanExecutor &executor,
    const FrameResult *map_result) {
  ReactionSolve solve;
  const double solvent_kappa = debye_kappa_inv_a(args_.salt, args_.temp);
  const FieldSolve solvent = solve_physical_field(
      atoms, args_.int_diel, args_.ext_diel, solvent_kappa, executor,
      map_result);
  const ForensicAudit solvent_forensic =
      barnase_forensic_enabled() ? compute_forensic_audit(solvent.shape)
                                 : ForensicAudit{};
  const bool wants_iscm =
      args_.test_energy_convention == "delphi-iscm" ||
      args_.test_energy_convention == "delphi-corrected" ||
      barnase_forensic_enabled();
  const AlignedDoubleVector solvent_screen_for_iscm =
      wants_iscm ? screen_ : AlignedDoubleVector{};
  const FieldSolve vacuum =
      solve_linear_field(atoms, args_.int_diel, args_.int_diel, 0.0, executor,
                         nullptr);

  solve.shape = solvent.shape;
  solve.solvent = solvent.stats;
  solve.vacuum = vacuum.stats;
  solve.solvent_total_iterations = solvent.total_iterations;
  solve.vacuum_total_iterations = vacuum.total_iterations;
  solve.focus_levels = solvent.levels;
  solve.focus_initial_scale = solvent.initial_scale;
  solve.focus_initial_margin = solvent.initial_margin;
  solve.solver_mode = solvent.solver_mode;
  solve.mg_levels = solvent.mg_levels;
  solve.mg_cycles = solvent.mg_cycles;
  solve.mg_pre = solvent.mg_pre;
  solve.mg_post = solvent.mg_post;
  solve.mg_coarse_iter = solvent.mg_coarse_iter;
  solve.solute_points = solvent.solute_points;
  solve.solvent_points = solvent.solvent_points;
  solve.ion_accessible_points = solvent.ion_accessible_points;
  solve.stern_excluded_points = solvent.stern_excluded_points;
  solve.peak_memory_bytes =
      std::max(solvent.peak_memory_bytes, vacuum.peak_memory_bytes);
  solve.pbe_model_effective = solvent.pbe_model_effective;
  solve.linear_seed_solver = solvent.linear_seed_solver;
  solve.nl_solver = solvent.nl_solver;
  solve.nl_inner_solver = solvent.nl_inner_solver;
  solve.nl_mg_mode = solvent.nl_mg_mode;
  solve.focus_parent_grid = solvent.focus_parent_grid;
  solve.focus_child_grid = solvent.focus_child_grid;
  solve.probe_status = solvent.probe_status;
  solve.newton_inner_iterations = solvent.newton_inner_iterations;
  solve.newton_step_norm = solvent.newton_step_norm;
  solve.newton_residual_before = solvent.newton_residual_before;
  solve.newton_residual_after = solvent.newton_residual_after;
  solve.newton_line_search_status = solvent.newton_line_search_status;
  solve.newton_correction_residual = solvent.newton_correction_residual;
  solve.newton_correction_max_residual =
      solvent.newton_correction_max_residual;
  solve.newton_accepted_lambda = solvent.newton_accepted_lambda;
  solve.newton_rejected_steps = solvent.newton_rejected_steps;
  solve.newton_jac_diag_min = solvent.newton_jac_diag_min;
  solve.newton_jac_diag_max = solvent.newton_jac_diag_max;
  solve.newton_screen_points = solvent.newton_screen_points;
  solve.newton_active_nonlinear_points =
      solvent.newton_active_nonlinear_points;
  solve.newton_stop_reason = solvent.newton_stop_reason;
  solve.newton_inner_levels = solvent.newton_inner_levels;
  solve.newton_inner_cycles = solvent.newton_inner_cycles;
  solve.newton_inner_residual_before = solvent.newton_inner_residual_before;
  solve.newton_inner_residual_after = solvent.newton_inner_residual_after;
  solve.newton_inner_effective = solvent.newton_inner_effective;
  solve.gaussian_density_sum = solvent.gaussian_density_sum;
  solve.gaussian_cutoff_points = solvent.gaussian_cutoff_points;
  solve.gaussian_eps_min = solvent.gaussian_eps_min;
  solve.gaussian_eps_max = solvent.gaussian_eps_max;
  solve.gaussian_eps_mean = solvent.gaussian_eps_mean;
  solve.gaussian_charge_norm_error = solvent.gaussian_charge_norm_error;
  solve.rhs_scale = solvent.rhs_scale;
  solve.charge_assigned = solvent.charge_assigned;
  solve.charge_moment_error = solvent.charge_moment_error;
  solve.boundary_rms = solvent.boundary_rms;
  solve.boundary_max = solvent.boundary_max;
  solve.linear_residual = solvent.linear_residual;
  solve.nl_iterations = solvent.nl_iterations;
  solve.nl_blocks = solvent.nl_blocks;
  solve.nl_relax = solvent.nl_relax;
  solve.nl_max_dphi = solvent.nl_max_dphi;
  solve.nl_residual = solvent.nl_residual;
  solve.nl_clipped_points = solvent.nl_clipped_points;
  solve.nl_continuation_steps = solvent.nl_continuation_steps;
  solve.nl_backtracks = solvent.nl_backtracks;
  solve.nl_status = solvent.nl_status;
  solve.parity = solvent.parity;
  solve.vdwms_enabled = solvent.vdwms_enabled;
  solve.vdwms_boundary_initial = solvent.vdwms_boundary_initial;
  solve.vdwms_boundary_final = solvent.vdwms_boundary_final;
  solve.vdwms_added = solvent.vdwms_added;
  solve.vdwms_removed = solvent.vdwms_removed;
  solve.vdwms_exposed_points = solvent.vdwms_exposed_points;
  solve.vdwms_midpoint_rms = solvent.vdwms_midpoint_rms;

  solve.atom_reaction_kj_mol.assign(atoms.size(), 0.0);
  const double kt_kcal =
      kKBKcalPerMolK * std::max(1.0, args_.temp);
  solve.solvent_phase_kcal_mol =
      atom_phase_energy_kcal(atoms, solvent.atom_phi);
  solve.vacuum_phase_kcal_mol =
      atom_phase_energy_kcal(atoms, vacuum.atom_phi);
  solve.atom_corrected_kcal_mol =
      solve.solvent_phase_kcal_mol - solve.vacuum_phase_kcal_mol;
  solve.coulomb_energy_kcal_mol =
      analytic_coulomb_energy_kcal(atoms, args_.int_diel);
  const double atom_solvent_phase_kcal = solve.solvent_phase_kcal_mol;
  const double atom_vacuum_phase_kcal = solve.vacuum_phase_kcal_mol;
  const bool has_grid_energy =
      solvent.phi_grid.size() == solvent.shape.n &&
      vacuum.phi_grid.size() == solvent.shape.n;
  IscmAudit iscm_audit;
  if (has_grid_energy) {
    const double grid_solvent_phase_kcal =
        grid_phase_energy_kcal(solvent.shape, solvent.phi_grid);
    const double grid_vacuum_phase_kcal =
        grid_phase_energy_kcal(solvent.shape, vacuum.phi_grid);
    solve.grid_corrected_kcal_mol =
        grid_solvent_phase_kcal - grid_vacuum_phase_kcal;
    if (wants_iscm) {
      iscm_audit = compute_iscm_audit(solvent.shape, atoms, solvent.phi_grid,
                                      solvent_screen_for_iscm,
                                      solvent.inside_grid, solvent.eps_x_grid,
                                      solvent.eps_y_grid, solvent.eps_z_grid,
                                      grid_vacuum_phase_kcal);
      copy_iscm_audit_to_solve(iscm_audit, solve);
    }
    const bool legacy_grid =
        args_.test_energy_convention == "auto" &&
        (args_.charge_assignment == "delphi" ||
         args_.charge_assignment == "gaussian");
    const bool forced_grid = args_.test_energy_convention == "grid";
    const bool forced_iscm = iscm_audit.status == "ok" &&
                             (args_.test_energy_convention ==
                                  "delphi-iscm" ||
                              args_.test_energy_convention ==
                                  "delphi-corrected");
    if (forced_iscm) {
      solve.solvent_phase_kcal_mol = solve.iscm_water_phase_kcal_mol;
      solve.vacuum_phase_kcal_mol = solve.iscm_vacuum_phase_kcal_mol;
      solve.energy_kcal_mol = solve.iscm_corrected_kcal_mol;
      solve.reaction_energy_source =
          args_.test_energy_convention == "delphi-corrected"
              ? "delphi-corrected-iscm"
              : "delphi-iscm";
    } else if (legacy_grid || forced_grid ||
               args_.test_energy_convention == "delphi-corrected" ||
               args_.test_energy_convention == "delphi-iscm") {
      solve.solvent_phase_kcal_mol = grid_solvent_phase_kcal;
      solve.vacuum_phase_kcal_mol = grid_vacuum_phase_kcal;
      solve.energy_kcal_mol = solve.grid_corrected_kcal_mol;
      if (args_.test_energy_convention == "delphi-corrected")
        solve.reaction_energy_source = "delphi-corrected-grid-fallback";
      else if (args_.test_energy_convention == "delphi-iscm")
        solve.reaction_energy_source = "delphi-iscm-grid-fallback";
      else
        solve.reaction_energy_source = "grid";
    } else {
      solve.solvent_phase_kcal_mol = atom_solvent_phase_kcal;
      solve.vacuum_phase_kcal_mol = atom_vacuum_phase_kcal;
      solve.energy_kcal_mol = solve.atom_corrected_kcal_mol;
      solve.reaction_energy_source = "atom";
    }
  } else {
    if (args_.test_energy_convention == "delphi-iscm" ||
        args_.test_energy_convention == "delphi-corrected")
      solve.iscm_status = "missing_grid_field";
    solve.grid_corrected_kcal_mol = solve.atom_corrected_kcal_mol;
    solve.energy_kcal_mol = solve.atom_corrected_kcal_mol;
    solve.reaction_energy_source = "atom";
  }
  if (args_.test_energy_convention == "atom") {
    solve.solvent_phase_kcal_mol = atom_solvent_phase_kcal;
    solve.vacuum_phase_kcal_mol = atom_vacuum_phase_kcal;
    solve.energy_kcal_mol = solve.atom_corrected_kcal_mol;
    solve.reaction_energy_source = "atom";
  }
  for (std::size_t i = 0; i < atoms.size(); ++i) {
    const double atom_kcal_mol =
        0.5 * kCoulombKcalPerMol * atoms[i].q *
        (solvent.atom_phi[i] - vacuum.atom_phi[i]);
    solve.atom_reaction_kj_mol[i] = atom_kcal_mol * kKcalToKJ;
  }
  solve.energy_kj_mol = solve.energy_kcal_mol * kKcalToKJ;
  solve.energy_kt =
      solve.energy_kj_mol / (kKBKJPerMolK * std::max(1.0, args_.temp));
  solve.solvent_phase_kj_mol = solve.solvent_phase_kcal_mol * kKcalToKJ;
  solve.solvent_phase_kt = kt_kcal > 0.0 ? solve.solvent_phase_kcal_mol / kt_kcal : 0.0;
  solve.vacuum_phase_kj_mol = solve.vacuum_phase_kcal_mol * kKcalToKJ;
  solve.vacuum_phase_kt = kt_kcal > 0.0 ? solve.vacuum_phase_kcal_mol / kt_kcal : 0.0;
  solve.atom_corrected_kj_mol = solve.atom_corrected_kcal_mol * kKcalToKJ;
  solve.atom_corrected_kt =
      kt_kcal > 0.0 ? solve.atom_corrected_kcal_mol / kt_kcal : 0.0;
  solve.grid_corrected_kj_mol = solve.grid_corrected_kcal_mol * kKcalToKJ;
  solve.grid_corrected_kt =
      kt_kcal > 0.0 ? solve.grid_corrected_kcal_mol / kt_kcal : 0.0;
  solve.coulomb_energy_kj_mol = solve.coulomb_energy_kcal_mol * kKcalToKJ;
  solve.coulomb_energy_kt =
      kt_kcal > 0.0 ? solve.coulomb_energy_kcal_mol / kt_kcal : 0.0;
  solve.converged = solve.solvent.converged && solve.vacuum.converged;
  if (map_result != nullptr) {
    maybe_write_iscm_surface_charge_samples(*map_result, iscm_audit);
    maybe_write_barnase_forensic_audits(*map_result, atoms, solvent, vacuum,
                                        solve, solvent_forensic);
    maybe_write_parity_map_dumps(*map_result, atoms, solvent, vacuum, solve);
  }
  return solve;
}

FrameResult Analyzer::consume_frame(const MOL &mol1,
                                    const Multi_MOL_parameters &params1,
                                    const MOL *mol2,
                                    const Multi_MOL_parameters *params2,
                                    std::size_t trajectory_id,
                                    std::size_t frame_index,
                                    std::size_t global_index,
                                    TitanExecutor &executor,
                                    std::fstream &log,
                                    const EnergyDecompositionProvider
                                        *energy_provider) {
  if (!enabled())
    return {};

  const auto frame_time_begin = std::chrono::steady_clock::now();
  atoms_.clear();
  std::vector<Atom> mol1_atoms;
  std::vector<Atom> mol2_atoms;
  if (!args_.test_pqr.empty()) {
    collect_test_pqr_atoms(atoms_);
  } else {
    collect_atoms(mol1, params1, mol1_atoms, "mol1");
    atoms_ = mol1_atoms;
    if (mol2 != nullptr && params2 != nullptr && mol2->size() > 0) {
      collect_atoms(*mol2, *params2, mol2_atoms, "mol2");
      atoms_.insert(atoms_.end(), mol2_atoms.begin(), mol2_atoms.end());
    }
  }
  if (args_.dgsolv && mol2_atoms.empty()) {
    throw std::runtime_error(
        "PBE: -pbe-dgsolv requires a two-molecule input with mol2 atoms.");
  }

  FrameResult result;
  result.trajectory_id = trajectory_id;
  result.frame_index = frame_index + 1;
  result.global_index = global_index + 1;
  result.time_ps = mol1.time_of_snap;
  result.atom_count = atoms_.size();
  result.surface_effective = args_.surface;
  result.ion_access = args_.ion_access;
  result.ion_radius_A = args_.ion_radius;
  result.face_refine = surface_uses_face_refinement(args_.surface);
  result.pbe_model = args_.model;
  result.pbe_biomodel = args_.biomodel;
  result.dielectric_model =
      uses_gaussian_model() ? std::string("gaussian") : args_.dielectric_model;
  result.nl_solver = args_.nl_solver;
  result.nl_inner_solver = args_.nl_inner;
  result.nl_mg_mode = args_.nl_multigrid;
  result.gaussian_exp = args_.gaussian_exp;
  result.gaussian_sigma = args_.gaussian_sigma;
  result.density_cutoff = args_.density_cutoff;
  result.surface_cutoff = args_.surface_cutoff == "auto"
                              ? args_.density_cutoff
                              : args_.surface_cutoff_value;
  result.gap_diel = args_.gap_diel;
  result.radius_convention = effective_radius_convention();
  result.boundary_mode = args_.boundary;
  result.charge_assignment = args_.charge_assignment;
  for (const Atom &atom : atoms_)
    result.total_charge += atom.q;

  const auto complex_time_begin = std::chrono::steady_clock::now();
  const ReactionSolve complex_solve =
      solve_reaction_field(atoms_, executor, &result);
  const auto complex_time_end = std::chrono::steady_clock::now();
  result.nx = complex_solve.shape.nx;
  result.ny = complex_solve.shape.ny;
  result.nz = complex_solve.shape.nz;
  result.grid_points = complex_solve.shape.n;
  result.memory_bytes = std::max(complex_solve.shape.memory_bytes,
                                 complex_solve.peak_memory_bytes);
  result.grid_origin_x_A = complex_solve.shape.ox;
  result.grid_origin_y_A = complex_solve.shape.oy;
  result.grid_origin_z_A = complex_solve.shape.oz;
  result.grid_spacing_A = complex_solve.shape.h;
  result.grid_center_x_A =
      complex_solve.shape.ox +
      0.5 * static_cast<double>(complex_solve.shape.nx - 1) *
          complex_solve.shape.h;
  result.grid_center_y_A =
      complex_solve.shape.oy +
      0.5 * static_cast<double>(complex_solve.shape.ny - 1) *
          complex_solve.shape.h;
  result.grid_center_z_A =
      complex_solve.shape.oz +
      0.5 * static_cast<double>(complex_solve.shape.nz - 1) *
          complex_solve.shape.h;
  result.solute_points = complex_solve.solute_points;
  result.solvent_points = complex_solve.solvent_points;
  result.ion_accessible_points = complex_solve.ion_accessible_points;
  result.stern_excluded_points = complex_solve.stern_excluded_points;
  result.pbe_model_effective = complex_solve.pbe_model_effective;
  result.linear_seed_solver = complex_solve.linear_seed_solver;
  result.nl_solver = complex_solve.nl_solver;
  result.nl_inner_solver = complex_solve.nl_inner_solver;
  result.nl_mg_mode = complex_solve.nl_mg_mode;
  result.focus_parent_grid = complex_solve.focus_parent_grid;
  result.focus_child_grid = complex_solve.focus_child_grid;
  result.probe_status = complex_solve.probe_status;
  result.newton_inner_iterations = complex_solve.newton_inner_iterations;
  result.newton_step_norm = complex_solve.newton_step_norm;
  result.newton_residual_before = complex_solve.newton_residual_before;
  result.newton_residual_after = complex_solve.newton_residual_after;
  result.newton_line_search_status = complex_solve.newton_line_search_status;
  result.newton_correction_residual =
      complex_solve.newton_correction_residual;
  result.newton_correction_max_residual =
      complex_solve.newton_correction_max_residual;
  result.newton_accepted_lambda = complex_solve.newton_accepted_lambda;
  result.newton_rejected_steps = complex_solve.newton_rejected_steps;
  result.newton_jac_diag_min = complex_solve.newton_jac_diag_min;
  result.newton_jac_diag_max = complex_solve.newton_jac_diag_max;
  result.newton_screen_points = complex_solve.newton_screen_points;
  result.newton_active_nonlinear_points =
      complex_solve.newton_active_nonlinear_points;
  result.newton_stop_reason = complex_solve.newton_stop_reason;
  result.newton_inner_levels = complex_solve.newton_inner_levels;
  result.newton_inner_cycles = complex_solve.newton_inner_cycles;
  result.newton_inner_residual_before =
      complex_solve.newton_inner_residual_before;
  result.newton_inner_residual_after =
      complex_solve.newton_inner_residual_after;
  result.newton_inner_effective = complex_solve.newton_inner_effective;
  result.gaussian_density_sum = complex_solve.gaussian_density_sum;
  result.gaussian_cutoff_points = complex_solve.gaussian_cutoff_points;
  result.gaussian_eps_min = complex_solve.gaussian_eps_min;
  result.gaussian_eps_max = complex_solve.gaussian_eps_max;
  result.gaussian_eps_mean = complex_solve.gaussian_eps_mean;
  result.gaussian_charge_norm_error = complex_solve.gaussian_charge_norm_error;
  result.rhs_scale = complex_solve.rhs_scale;
  result.charge_assigned = complex_solve.charge_assigned;
  result.charge_moment_error = complex_solve.charge_moment_error;
  result.boundary_rms = complex_solve.boundary_rms;
  result.boundary_max = complex_solve.boundary_max;
  result.linear_residual = complex_solve.linear_residual;
  result.nl_iterations = complex_solve.nl_iterations;
  result.nl_blocks = complex_solve.nl_blocks;
  result.nl_relax = complex_solve.nl_relax;
  result.nl_max_dphi = complex_solve.nl_max_dphi;
  result.nl_residual = complex_solve.nl_residual;
  result.nl_clipped_points = complex_solve.nl_clipped_points;
  result.nl_continuation_steps = complex_solve.nl_continuation_steps;
  result.nl_backtracks = complex_solve.nl_backtracks;
  result.nl_status = complex_solve.nl_status;
  result.vdwms_enabled = complex_solve.vdwms_enabled;
  result.vdwms_boundary_initial = complex_solve.vdwms_boundary_initial;
  result.vdwms_boundary_final = complex_solve.vdwms_boundary_final;
  result.vdwms_added = complex_solve.vdwms_added;
  result.vdwms_removed = complex_solve.vdwms_removed;
  result.vdwms_exposed_points = complex_solve.vdwms_exposed_points;
  result.vdwms_midpoint_rms = complex_solve.vdwms_midpoint_rms;
  result.focus_enabled =
      (complex_solve.solver_mode == "nested_focus" ||
       complex_solve.solver_mode == "delphi_focus") &&
      complex_solve.focus_levels > 1;
  result.focus_levels = complex_solve.focus_levels;
  result.focus_initial_scale = complex_solve.focus_initial_scale;
  result.focus_initial_margin = complex_solve.focus_initial_margin;
  result.solver_mode = complex_solve.solver_mode;
  result.mg_levels = complex_solve.mg_levels;
  result.mg_cycles = complex_solve.mg_cycles;
  result.mg_pre = complex_solve.mg_pre;
  result.mg_post = complex_solve.mg_post;
  result.mg_coarse_iter = complex_solve.mg_coarse_iter;
  result.solvent_iterations = complex_solve.solvent_total_iterations;
  result.vacuum_iterations = complex_solve.vacuum_total_iterations;
  result.solvent_iterations_final = complex_solve.solvent.iterations;
  result.vacuum_iterations_final = complex_solve.vacuum.iterations;
  result.solvent_rms = complex_solve.solvent.rms;
  result.vacuum_rms = complex_solve.vacuum.rms;
  result.solvent_max_dphi = complex_solve.solvent.max_dphi;
  result.vacuum_max_dphi = complex_solve.vacuum.max_dphi;
  result.solvent_converged = complex_solve.solvent.converged;
  result.vacuum_converged = complex_solve.vacuum.converged;
  result.energy_kcal_mol = complex_solve.energy_kcal_mol;
  result.energy_kj_mol = complex_solve.energy_kj_mol;
  result.energy_kt = complex_solve.energy_kt;
  result.solvent_phase_kcal_mol = complex_solve.solvent_phase_kcal_mol;
  result.solvent_phase_kj_mol = complex_solve.solvent_phase_kj_mol;
  result.solvent_phase_kt = complex_solve.solvent_phase_kt;
  result.vacuum_phase_kcal_mol = complex_solve.vacuum_phase_kcal_mol;
  result.vacuum_phase_kj_mol = complex_solve.vacuum_phase_kj_mol;
  result.vacuum_phase_kt = complex_solve.vacuum_phase_kt;
  result.atom_corrected_kcal_mol = complex_solve.atom_corrected_kcal_mol;
  result.atom_corrected_kj_mol = complex_solve.atom_corrected_kj_mol;
  result.atom_corrected_kt = complex_solve.atom_corrected_kt;
  result.grid_corrected_kcal_mol = complex_solve.grid_corrected_kcal_mol;
  result.grid_corrected_kj_mol = complex_solve.grid_corrected_kj_mol;
  result.grid_corrected_kt = complex_solve.grid_corrected_kt;
  result.iscm_water_phase_kcal_mol =
      complex_solve.iscm_water_phase_kcal_mol;
  result.iscm_water_phase_kj_mol = complex_solve.iscm_water_phase_kj_mol;
  result.iscm_water_phase_kt = complex_solve.iscm_water_phase_kt;
  result.iscm_vacuum_phase_kcal_mol =
      complex_solve.iscm_vacuum_phase_kcal_mol;
  result.iscm_vacuum_phase_kj_mol = complex_solve.iscm_vacuum_phase_kj_mol;
  result.iscm_vacuum_phase_kt = complex_solve.iscm_vacuum_phase_kt;
  result.iscm_corrected_kcal_mol = complex_solve.iscm_corrected_kcal_mol;
  result.iscm_corrected_kj_mol = complex_solve.iscm_corrected_kj_mol;
  result.iscm_corrected_kt = complex_solve.iscm_corrected_kt;
  result.iscm_boundary_points = complex_solve.iscm_boundary_points;
  result.iscm_charged_boundary_points =
      complex_solve.iscm_charged_boundary_points;
  result.iscm_induced_charge_sum = complex_solve.iscm_induced_charge_sum;
  result.iscm_induced_charge_l1 = complex_solve.iscm_induced_charge_l1;
  result.iscm_induced_charge_l2 = complex_solve.iscm_induced_charge_l2;
  result.iscm_induced_charge_max = complex_solve.iscm_induced_charge_max;
  result.iscm_status = complex_solve.iscm_status;
  result.delphi_face_crossings = complex_solve.delphi_face_crossings;
  result.delphi_eps_face_matches = complex_solve.delphi_eps_face_matches;
  result.delphi_boundary_grid_points =
      complex_solve.delphi_boundary_grid_points;
  result.delphi_charged_boundary_points =
      complex_solve.delphi_charged_boundary_points;
  result.delphi_iscm_status = complex_solve.delphi_iscm_status;
  result.vdwms_enabled = complex_solve.vdwms_enabled;
  result.vdwms_boundary_initial = complex_solve.vdwms_boundary_initial;
  result.vdwms_boundary_final = complex_solve.vdwms_boundary_final;
  result.vdwms_added = complex_solve.vdwms_added;
  result.vdwms_removed = complex_solve.vdwms_removed;
  result.vdwms_exposed_points = complex_solve.vdwms_exposed_points;
  result.vdwms_midpoint_rms = complex_solve.vdwms_midpoint_rms;
  result.coulomb_energy_kcal_mol = complex_solve.coulomb_energy_kcal_mol;
  result.coulomb_energy_kj_mol = complex_solve.coulomb_energy_kj_mol;
  result.coulomb_energy_kt = complex_solve.coulomb_energy_kt;
  result.reaction_energy_source = complex_solve.reaction_energy_source;
  result.complex_energy_kcal_mol = complex_solve.energy_kcal_mol;
  result.complex_energy_kj_mol = complex_solve.energy_kj_mol;

  ReactionSolve mol1_solve;
  ReactionSolve mol2_solve;
  bool have_endpoint_solve = false;
  if (args_.dgsolv) {
    const auto mol1_time_begin = std::chrono::steady_clock::now();
    mol1_solve = solve_reaction_field(mol1_atoms, executor, nullptr);
    const auto mol1_time_end = std::chrono::steady_clock::now();
    const auto mol2_time_begin = std::chrono::steady_clock::now();
    mol2_solve = solve_reaction_field(mol2_atoms, executor, nullptr);
    const auto mol2_time_end = std::chrono::steady_clock::now();
    result.pbsa_time_receptor_ms = elapsed_ms(mol1_time_begin, mol1_time_end);
    result.pbsa_time_ligand_ms = elapsed_ms(mol2_time_begin, mol2_time_end);
    have_endpoint_solve = true;
    result.dgsolv_enabled = true;
    result.mol1_energy_kcal_mol = mol1_solve.energy_kcal_mol;
    result.mol1_energy_kj_mol = mol1_solve.energy_kj_mol;
    result.mol2_energy_kcal_mol = mol2_solve.energy_kcal_mol;
    result.mol2_energy_kj_mol = mol2_solve.energy_kj_mol;
    result.deltaG_solv_kcal_mol =
        complex_solve.energy_kcal_mol - mol1_solve.energy_kcal_mol -
        mol2_solve.energy_kcal_mol;
    result.deltaG_solv_kj_mol =
        complex_solve.energy_kj_mol - mol1_solve.energy_kj_mol -
        mol2_solve.energy_kj_mol;
    result.deltaG_solv_kt =
        result.deltaG_solv_kj_mol /
        (kKBKJPerMolK * std::max(1.0, args_.temp));
    result.dgsolv_status = (complex_solve.converged && mol1_solve.converged &&
                            mol2_solve.converged)
                               ? "ok"
                               : "not_converged";
  }

  if (args_.pbsa) {
    result.pbsa_time_complex_ms =
        elapsed_ms(complex_time_begin, complex_time_end);
    const auto components_time_begin = std::chrono::steady_clock::now();
    const PbsaComponents pbsa = compute_pbsa_components(
        mol1, params1, *mol2, *params2, mol1_atoms, mol2_atoms, atoms_,
        result.deltaG_solv_kj_mol, executor, energy_provider);
    const auto components_time_end = std::chrono::steady_clock::now();
    result.pbsa_time_components_ms =
        elapsed_ms(components_time_begin, components_time_end);
    result.pbsa_enabled = true;
    result.pbsa_status =
        (have_endpoint_solve && complex_solve.converged &&
         mol1_solve.converged && mol2_solve.converged &&
         pbsa.status == "ok")
            ? "ok"
            : "not_converged";
    result.pbsa_elec_kj_mol = pbsa.elec_kj_mol;
    result.pbsa_vdw_kj_mol = pbsa.vdw_kj_mol;
    result.pbsa_mm_kj_mol = pbsa.mm_kj_mol;
    result.pbsa_sasa_complex_a2 = pbsa.sasa_complex_a2;
    result.pbsa_sasa_mol1_a2 = pbsa.sasa_mol1_a2;
    result.pbsa_sasa_mol2_a2 = pbsa.sasa_mol2_a2;
    result.pbsa_delta_sasa_a2 = pbsa.delta_sasa_a2;
    result.pbsa_nonpolar_kj_mol = pbsa.nonpolar_kj_mol;
    result.pbsa_delta_solv_kj_mol = pbsa.delta_solv_kj_mol;
    result.pbsa_delta_solv_kcal_mol = pbsa.delta_solv_kj_mol * kKJToKcal;
    result.pbsa_no_entropy_kj_mol = pbsa.no_entropy_kj_mol;
    result.pbsa_no_entropy_kcal_mol = pbsa.no_entropy_kj_mol * kKJToKcal;
    result.pbsa_receptor_atom_count = mol1_atoms.size();
    result.pbsa_ligand_atom_count = mol2_atoms.size();
    result.pbsa_complex_status =
        complex_solve.converged ? "ok" : "not_converged";
    result.pbsa_receptor_status =
        have_endpoint_solve
            ? (mol1_solve.converged ? "ok" : "not_converged")
            : "missing_endpoint";
    result.pbsa_ligand_status =
        have_endpoint_solve
            ? (mol2_solve.converged ? "ok" : "not_converged")
            : "missing_endpoint";
    if (have_endpoint_solve) {
      result.pbsa_receptor_solvent_iterations =
          mol1_solve.solvent_total_iterations;
      result.pbsa_receptor_vacuum_iterations =
          mol1_solve.vacuum_total_iterations;
      result.pbsa_ligand_solvent_iterations =
          mol2_solve.solvent_total_iterations;
      result.pbsa_ligand_vacuum_iterations =
          mol2_solve.vacuum_total_iterations;
      result.pbsa_receptor_mg_cycles = mol1_solve.mg_cycles;
      result.pbsa_ligand_mg_cycles = mol2_solve.mg_cycles;
      result.pbsa_receptor_solvent_rms = mol1_solve.solvent.rms;
      result.pbsa_receptor_vacuum_rms = mol1_solve.vacuum.rms;
      result.pbsa_ligand_solvent_rms = mol2_solve.solvent.rms;
      result.pbsa_ligand_vacuum_rms = mol2_solve.vacuum.rms;
      result.pbsa_receptor_solvent_max_dphi = mol1_solve.solvent.max_dphi;
      result.pbsa_receptor_vacuum_max_dphi = mol1_solve.vacuum.max_dphi;
      result.pbsa_ligand_solvent_max_dphi = mol2_solve.solvent.max_dphi;
      result.pbsa_ligand_vacuum_max_dphi = mol2_solve.vacuum.max_dphi;
    }
    if (args_.pbsa_decomp == "residue") {
      const auto residue_time_begin = std::chrono::steady_clock::now();
      const std::vector<ResidueFrameRow> residue_rows =
          compute_pbsa_residue_rows(mol1, *mol2, mol1_atoms, mol2_atoms,
                                    atoms_, complex_solve, mol1_solve,
                                    mol2_solve, result, executor,
                                    energy_provider);
      const auto residue_time_end = std::chrono::steady_clock::now();
      result.pbsa_time_residue_ms =
          elapsed_ms(residue_time_begin, residue_time_end);
      accumulate_pbsa_residue_rows(result, residue_rows);
    }
  }

  result.status =
      (result.solvent_converged && result.vacuum_converged) ? "ok"
                                                            : "not_converged";
  if (result.dgsolv_enabled && result.dgsolv_status != "ok")
    result.status = "not_converged";
  if (result.pbsa_enabled && result.pbsa_status != "ok")
    result.status = "not_converged";
  result.pbsa_time_frame_ms =
      elapsed_ms(frame_time_begin, std::chrono::steady_clock::now());
  result.probe_status = maybe_write_probe_sites(result, atoms_, complex_solve);

  maybe_write_parity_stage(result, complex_solve);
  append_energy_row(result);

  ++frames_total_;
  if (result.status == "ok")
    ++frames_converged_;
  else
    ++frames_failed_;
  max_grid_points_seen_ = std::max(max_grid_points_seen_, result.grid_points);
  max_memory_bytes_seen_ = std::max(max_memory_bytes_seen_, result.memory_bytes);
  max_solute_points_seen_ =
      std::max(max_solute_points_seen_, result.solute_points);
  max_solvent_points_seen_ =
      std::max(max_solvent_points_seen_, result.solvent_points);
  max_ion_accessible_points_seen_ =
      std::max(max_ion_accessible_points_seen_, result.ion_accessible_points);
  max_stern_excluded_points_seen_ =
      std::max(max_stern_excluded_points_seen_, result.stern_excluded_points);
  max_focus_levels_seen_ = std::max(max_focus_levels_seen_, result.focus_levels);
  max_mg_levels_seen_ = std::max(max_mg_levels_seen_, result.mg_levels);
  max_mg_cycles_seen_ = std::max(max_mg_cycles_seen_, result.mg_cycles);
  if (result.pbe_model == "nonlinear") {
    ++nonlinear_frames_;
    if (result.nl_status == "ok" || result.nl_status == "linearized_no_salt")
      ++nonlinear_converged_;
    else
      ++nonlinear_failed_;
    max_nl_blocks_seen_ = std::max(max_nl_blocks_seen_, result.nl_blocks);
    max_nl_iterations_seen_ =
        std::max(max_nl_iterations_seen_, result.nl_iterations);
    max_nl_clipped_points_seen_ =
        std::max(max_nl_clipped_points_seen_, result.nl_clipped_points);
    max_nl_dphi_seen_ = std::max(max_nl_dphi_seen_, result.nl_max_dphi);
    max_nl_residual_seen_ =
        std::max(max_nl_residual_seen_, result.nl_residual);
    max_newton_inner_iterations_seen_ = std::max(
        max_newton_inner_iterations_seen_, result.newton_inner_iterations);
    max_newton_step_norm_seen_ =
        std::max(max_newton_step_norm_seen_, result.newton_step_norm);
    max_newton_residual_before_seen_ = std::max(
        max_newton_residual_before_seen_, result.newton_residual_before);
    max_newton_residual_after_seen_ = std::max(
        max_newton_residual_after_seen_, result.newton_residual_after);
    max_newton_correction_residual_seen_ = std::max(
        max_newton_correction_residual_seen_,
        result.newton_correction_residual);
    max_newton_correction_max_residual_seen_ = std::max(
        max_newton_correction_max_residual_seen_,
        result.newton_correction_max_residual);
    max_newton_rejected_steps_seen_ =
        std::max(max_newton_rejected_steps_seen_, result.newton_rejected_steps);
    max_newton_inner_levels_seen_ =
        std::max(max_newton_inner_levels_seen_, result.newton_inner_levels);
    max_newton_inner_cycles_seen_ =
        std::max(max_newton_inner_cycles_seen_, result.newton_inner_cycles);
    max_newton_inner_residual_before_seen_ = std::max(
        max_newton_inner_residual_before_seen_,
        result.newton_inner_residual_before);
    max_newton_inner_residual_after_seen_ = std::max(
        max_newton_inner_residual_after_seen_,
        result.newton_inner_residual_after);
  }
  max_linear_residual_seen_ =
      std::max(max_linear_residual_seen_, result.linear_residual);
  max_boundary_rms_seen_ = std::max(max_boundary_rms_seen_, result.boundary_rms);
  max_boundary_max_seen_ = std::max(max_boundary_max_seen_, result.boundary_max);

  const double frame_n = static_cast<double>(frames_total_);
  const double delta = result.energy_kj_mol - energy_mean_;
  energy_mean_ += delta / frame_n;
  energy_m2_ += delta * (result.energy_kj_mol - energy_mean_);
  if (frames_total_ == 1) {
    energy_min_ = result.energy_kj_mol;
    energy_max_ = result.energy_kj_mol;
  } else {
    energy_min_ = std::min(energy_min_, result.energy_kj_mol);
    energy_max_ = std::max(energy_max_, result.energy_kj_mol);
  }
  iter_mean_ +=
      (static_cast<double>(result.solvent_iterations) - iter_mean_) / frame_n;
  mg_cycle_mean_ +=
      (static_cast<double>(result.mg_cycles) - mg_cycle_mean_) / frame_n;

  if (result.dgsolv_enabled) {
    ++dgsolv_frames_;
    const double dg_n = static_cast<double>(dgsolv_frames_);
    const double dg_delta = result.deltaG_solv_kj_mol - dgsolv_mean_;
    dgsolv_mean_ += dg_delta / dg_n;
    dgsolv_m2_ += dg_delta * (result.deltaG_solv_kj_mol - dgsolv_mean_);
    if (dgsolv_frames_ == 1) {
      dgsolv_min_ = result.deltaG_solv_kj_mol;
      dgsolv_max_ = result.deltaG_solv_kj_mol;
    } else {
      dgsolv_min_ = std::min(dgsolv_min_, result.deltaG_solv_kj_mol);
      dgsolv_max_ = std::max(dgsolv_max_, result.deltaG_solv_kj_mol);
    }
  }
  accumulate_pbsa(result);

  if (log.is_open() && result.status != "ok") {
    TitanUI::log_kv(log, "PBE Warning",
                    "frame " + std::to_string(result.global_index) +
                        " did not meet convergence thresholds");
  }
  return result;
}

void Analyzer::add_scalar(ScalarStats &stats, double value) {
  if (!std::isfinite(value))
    return;
  ++stats.count;
  const double n = static_cast<double>(stats.count);
  const double delta = value - stats.mean;
  stats.mean += delta / n;
  stats.m2 += delta * (value - stats.mean);
  if (stats.count == 1) {
    stats.min = value;
    stats.max = value;
  } else {
    stats.min = std::min(stats.min, value);
    stats.max = std::max(stats.max, value);
  }
}

double Analyzer::scalar_stddev(const ScalarStats &stats) {
  if (stats.count <= 1)
    return 0.0;
  return std::sqrt(stats.m2 / static_cast<double>(stats.count - 1));
}

void Analyzer::accumulate_pbsa(const FrameResult &result) {
  if (!result.pbsa_enabled)
    return;
  ++pbsa_frames_;
  if (result.pbsa_status != "ok")
    ++pbsa_failed_;
  add_scalar(pbsa_elec_stats_, result.pbsa_elec_kj_mol);
  add_scalar(pbsa_vdw_stats_, result.pbsa_vdw_kj_mol);
  add_scalar(pbsa_mm_stats_, result.pbsa_mm_kj_mol);
  add_scalar(pbsa_delta_sasa_stats_, result.pbsa_delta_sasa_a2);
  add_scalar(pbsa_nonpolar_stats_, result.pbsa_nonpolar_kj_mol);
  add_scalar(pbsa_polar_stats_, result.pbsa_delta_solv_kj_mol);
  add_scalar(pbsa_no_entropy_stats_, result.pbsa_no_entropy_kj_mol);
  pbsa_rows_.push_back(result);
  if (args_.pbsa_ie && result.pbsa_status == "ok")
    pbsa_ie_mm_energies_kj_.push_back(result.pbsa_mm_kj_mol);
}

Analyzer::PbsaEntropyResult
Analyzer::compute_pbsa_entropy_for_values(
    const std::vector<double> &values) const {
  PbsaEntropyResult out;
  out.samples = values.size();
  if (!args_.pbsa_ie) {
    out.status = "disabled";
    return out;
  }

  std::size_t block_size = 0;
  if (args_.pbsa_ie_block == "auto") {
    block_size = std::max<std::size_t>(
        1, static_cast<std::size_t>(std::floor(std::sqrt(values.size()))));
  } else {
    block_size =
        std::max<std::size_t>(1, static_cast<std::size_t>(
                                     std::stoi(args_.pbsa_ie_block)));
  }

  TitanEntropy::InteractionEntropyOptions options;
  options.rt_energy = kKBKJPerMolK * std::max(1.0, args_.temp);
  options.filter_sigma = args_.pbsa_ie_filter_sigma;
  options.min_samples =
      static_cast<std::size_t>(std::max(1, args_.pbsa_ie_min_samples));
  options.block_size = block_size;
  const TitanEntropy::InteractionEntropyResult entropy =
      TitanEntropy::compute_interaction_entropy(
          std::span<const double>(values.data(), values.size()), options);
  out.samples = entropy.samples;
  out.filtered_samples = entropy.filtered_samples;
  out.ess = entropy.ess;
  out.minus_tds_kj_mol = entropy.minus_tds_energy;
  out.block_size = entropy.block_size;
  out.block_count = entropy.block_count;
  out.block_mean_kj_mol = entropy.block_mean_energy;
  out.block_std_kj_mol = entropy.block_std_energy;
  out.block_drift_kj_mol = entropy.block_drift_energy;
  out.status = entropy.status;
  return out;
}

Analyzer::PbsaEntropyResult Analyzer::compute_pbsa_entropy() const {
  return compute_pbsa_entropy_for_values(pbsa_ie_mm_energies_kj_);
}

void Analyzer::write_pbsa_energy_csv(const PbsaEntropyResult &entropy) const {
  if (!args_.pbsa)
    return;

  std::ofstream out("mm_pbsa_energy.csv", std::ios::out | std::ios::trunc);
  if (!out.is_open()) {
    TitanUI::print_warning("PBSA: failed to write mm_pbsa_energy.csv.");
    return;
  }

  out << "trajectory_id,frame_index,global_index,time_ps,backend,"
         "atom_count_complex,atom_count_receptor,atom_count_ligand,"
         "grid_nx,grid_ny,grid_nz,solver_mode,mg_levels,"
         "mg_cycles_complex,mg_cycles_receptor,mg_cycles_ligand,"
         "status,complex_status,receptor_status,ligand_status,"
         "complex_solvent_iter,complex_vacuum_iter,"
         "receptor_solvent_iter,receptor_vacuum_iter,"
         "ligand_solvent_iter,ligand_vacuum_iter,"
         "complex_solvent_rms,complex_vacuum_rms,"
         "receptor_solvent_rms,receptor_vacuum_rms,"
         "ligand_solvent_rms,ligand_vacuum_rms,"
         "complex_solvent_max_dphi,complex_vacuum_max_dphi,"
         "receptor_solvent_max_dphi,receptor_vacuum_max_dphi,"
         "ligand_solvent_max_dphi,ligand_vacuum_max_dphi,"
         "E_elec_kJmol,E_elec_kcalmol,E_vdw_kJmol,E_vdw_kcalmol,"
         "E_MM_kJmol,E_MM_kcalmol,"
         "G_solv_complex_kJmol,G_solv_complex_kcalmol,"
         "G_solv_receptor_kJmol,G_solv_receptor_kcalmol,"
         "G_solv_ligand_kJmol,G_solv_ligand_kcalmol,"
         "DeltaDeltaG_solv_kJmol,DeltaDeltaG_solv_kcalmol,"
         "DeltaG_nonpolar_kJmol,DeltaG_nonpolar_kcalmol,"
         "deltaSASA_A2,minus_TDeltaS_IE_kJmol,"
         "minus_TDeltaS_IE_kcalmol,DeltaG_bind_kJmol,"
         "DeltaG_bind_kcalmol,ie_status\n";

  out << std::setprecision(10);
  const double minus_tds =
      (args_.pbsa_ie && entropy.status == "ok") ? entropy.minus_tds_kj_mol
                                                 : 0.0;
  for (const FrameResult &row : pbsa_rows_) {
    const double dg_bind = row.pbsa_no_entropy_kj_mol + minus_tds;
    out << row.trajectory_id << ',' << row.frame_index << ','
        << row.global_index << ',' << row.time_ps << ",pbe,"
        << row.atom_count << ',' << row.pbsa_receptor_atom_count << ','
        << row.pbsa_ligand_atom_count << ',' << row.nx << ',' << row.ny
        << ',' << row.nz << ',' << row.solver_mode << ',' << row.mg_levels
        << ',' << row.mg_cycles << ',' << row.pbsa_receptor_mg_cycles << ','
        << row.pbsa_ligand_mg_cycles << ',' << row.pbsa_status << ','
        << row.pbsa_complex_status << ',' << row.pbsa_receptor_status << ','
        << row.pbsa_ligand_status << ',' << row.solvent_iterations << ','
        << row.vacuum_iterations << ','
        << row.pbsa_receptor_solvent_iterations << ','
        << row.pbsa_receptor_vacuum_iterations << ','
        << row.pbsa_ligand_solvent_iterations << ','
        << row.pbsa_ligand_vacuum_iterations << ',' << row.solvent_rms << ','
        << row.vacuum_rms << ',' << row.pbsa_receptor_solvent_rms << ','
        << row.pbsa_receptor_vacuum_rms << ','
        << row.pbsa_ligand_solvent_rms << ','
        << row.pbsa_ligand_vacuum_rms << ','
        << row.solvent_max_dphi << ',' << row.vacuum_max_dphi << ','
        << row.pbsa_receptor_solvent_max_dphi << ','
        << row.pbsa_receptor_vacuum_max_dphi << ','
        << row.pbsa_ligand_solvent_max_dphi << ','
        << row.pbsa_ligand_vacuum_max_dphi << ',' << row.pbsa_elec_kj_mol
        << ',' << row.pbsa_elec_kj_mol * kKJToKcal << ','
        << row.pbsa_vdw_kj_mol << ',' << row.pbsa_vdw_kj_mol * kKJToKcal
        << ',' << row.pbsa_mm_kj_mol << ','
        << row.pbsa_mm_kj_mol * kKJToKcal << ','
        << row.complex_energy_kj_mol << ','
        << row.complex_energy_kcal_mol << ',' << row.mol1_energy_kj_mol
        << ',' << row.mol1_energy_kcal_mol << ','
        << row.mol2_energy_kj_mol << ',' << row.mol2_energy_kcal_mol << ','
        << row.pbsa_delta_solv_kj_mol << ','
        << row.pbsa_delta_solv_kcal_mol << ','
        << row.pbsa_nonpolar_kj_mol << ','
        << row.pbsa_nonpolar_kj_mol * kKJToKcal << ','
        << row.pbsa_delta_sasa_a2 << ',' << minus_tds << ','
        << minus_tds * kKJToKcal << ',' << dg_bind << ','
        << dg_bind * kKJToKcal << ',' << entropy.status << '\n';
  }
}

void Analyzer::write_pbsa_decomposition(
    const PbsaEntropyResult &entropy) const {
  if (!args_.pbsa || args_.pbsa_decomp == "off")
    return;

  std::ofstream out("mm_pbsa_decomposition.dat",
                    std::ios::out | std::ios::trunc);
  if (!out.is_open()) {
    TitanUI::print_warning(
        "PBSA: failed to write mm_pbsa_decomposition.dat.");
    return;
  }

  out << "# TITAN MM/PBSA endpoint decomposition\n";
  out << "requested_decomp " << args_.pbsa_decomp << "\n";
  out << "effective_decomp "
      << (args_.pbsa_decomp == "residue" ? "residue" : "molecule") << "\n";
  out << "residue_summary "
      << (args_.pbsa_decomp == "residue" ? "mm_pbsa_residue_summary.csv"
                                          : "off")
      << "\n";
  out << "residue_by_trajectory "
      << (args_.pbsa_decomp == "residue"
              ? "mm_pbsa_residue_by_trajectory.csv"
              : "off")
      << "\n";
  out << "residue_frame "
      << (args_.pbsa_decomp == "residue" &&
                  args_.pbsa_residue_frames == "all"
              ? "mm_pbsa_residue_frame.csv"
              : "off")
      << "\n";
  out << "residue_closure_total_max_abs_kJmol "
      << format_fixed(pbsa_residue_closure_max_abs_kj_, 10) << "\n";
  out << "residue_closure_elec_max_abs_kJmol "
      << format_fixed(pbsa_residue_closure_elec_max_abs_kj_, 10) << "\n";
  out << "residue_closure_vdw_max_abs_kJmol "
      << format_fixed(pbsa_residue_closure_vdw_max_abs_kj_, 10) << "\n";
  out << "residue_closure_polar_max_abs_kJmol "
      << format_fixed(pbsa_residue_closure_polar_max_abs_kj_, 10) << "\n";
  out << "residue_closure_nonpolar_max_abs_kJmol "
      << format_fixed(pbsa_residue_closure_nonpolar_max_abs_kj_, 10) << "\n";
  out << "interaction_entropy_status " << entropy.status << "\n";
  out << "component,mean_kJmol,mean_kcalmol,std_kJmol,std_kcalmol\n";
  auto write_stat_row = [&](const char *name, const ScalarStats &stats) {
    const double sd = scalar_stddev(stats);
    out << name << ',' << format_fixed(stats.mean) << ','
        << format_fixed(stats.mean * kKJToKcal) << ',' << format_fixed(sd)
        << ',' << format_fixed(sd * kKJToKcal) << "\n";
  };
  write_stat_row("DeltaE_elec", pbsa_elec_stats_);
  write_stat_row("DeltaE_vdw", pbsa_vdw_stats_);
  write_stat_row("DeltaE_MM", pbsa_mm_stats_);
  write_stat_row("DeltaDeltaG_solv_PB", pbsa_polar_stats_);
  write_stat_row("DeltaG_nonpolar_SASA", pbsa_nonpolar_stats_);
  write_stat_row("DeltaG_no_entropy", pbsa_no_entropy_stats_);
}

void Analyzer::write_pbsa_residue_outputs() const {
  if (!args_.pbsa || args_.pbsa_decomp != "residue")
    return;

  auto write_key = [](std::ofstream &out, const ResidueKey &key) {
    out << csv_escape(key.molecule) << ',' << key.chainnum << ','
        << csv_escape(key.chain) << ',' << key.resnr << ',' << key.respos
        << ',' << csv_escape(key.resname);
  };
  auto write_stats_values = [&](std::ofstream &out, const ScalarStats &stats) {
    out << ',' << format_fixed(stats.mean, 10) << ','
        << format_fixed(scalar_stddev(stats), 10) << ','
        << format_fixed(stats.min, 10) << ',' << format_fixed(stats.max, 10);
  };
  auto write_stats_header = [](std::ofstream &out) {
    const std::array<const char *, 6> names = {
        "E_elec", "E_vdw", "E_MM", "DeltaDeltaG_solv_PB",
        "DeltaG_nonpolar_SASA", "DeltaG_residue"};
    for (const char *name : names) {
      out << ',' << name << "_kJmol_mean," << name << "_kJmol_std," << name
          << "_kJmol_min," << name << "_kJmol_max";
    }
    out << ",DeltaG_residue_kcalmol_mean\n";
  };
  auto write_stats_row = [&](std::ofstream &out, const ResidueStats &stats) {
    write_key(out, stats.key);
    out << ',' << stats.frames;
    write_stats_values(out, stats.elec);
    write_stats_values(out, stats.vdw);
    write_stats_values(out, stats.mm);
    write_stats_values(out, stats.polar);
    write_stats_values(out, stats.nonpolar);
    write_stats_values(out, stats.total);
    out << ',' << format_fixed(stats.total.mean * kKJToKcal, 10) << "\n";
  };

  std::ofstream summary("mm_pbsa_residue_summary.csv",
                        std::ios::out | std::ios::trunc);
  if (summary.is_open()) {
    summary << "molecule,chainnum,chain,resnr,respos,resname,frames";
    write_stats_header(summary);
    for (const auto &[id, stats] : pbsa_residue_stats_) {
      (void)id;
      write_stats_row(summary, stats);
    }
  } else {
    TitanUI::print_warning(
        "PBSA: failed to write mm_pbsa_residue_summary.csv.");
  }

  std::ofstream by_traj("mm_pbsa_residue_by_trajectory.csv",
                        std::ios::out | std::ios::trunc);
  if (by_traj.is_open()) {
    by_traj << "trajectory_id,molecule,chainnum,chain,resnr,respos,resname,"
               "frames";
    write_stats_header(by_traj);
    for (const auto &[id, stats] : pbsa_residue_traj_stats_) {
      const std::size_t cut = id.find('|');
      by_traj << (cut == std::string::npos ? id : id.substr(0, cut)) << ',';
      write_stats_row(by_traj, stats);
    }
  } else {
    TitanUI::print_warning(
        "PBSA: failed to write mm_pbsa_residue_by_trajectory.csv.");
  }

  if (args_.pbsa_residue_frames == "all") {
    std::ofstream frame("mm_pbsa_residue_frame.csv",
                        std::ios::out | std::ios::trunc);
    if (frame.is_open()) {
      frame << "trajectory_id,frame_index,global_index,time_ps,molecule,"
               "chainnum,chain,resnr,respos,resname,E_elec_kJmol,"
               "E_vdw_kJmol,E_MM_kJmol,DeltaDeltaG_solv_PB_kJmol,"
               "DeltaG_nonpolar_SASA_kJmol,DeltaG_residue_kJmol,"
               "DeltaG_residue_kcalmol\n";
      frame << std::setprecision(10);
      for (const ResidueFrameRow &row : pbsa_residue_frame_rows_) {
        frame << row.trajectory_id << ',' << row.frame_index << ','
              << row.global_index << ',' << row.time_ps << ',';
        write_key(frame, row.key);
        frame << ',' << row.elec_kj_mol << ',' << row.vdw_kj_mol << ','
              << row.mm_kj_mol << ',' << row.polar_kj_mol << ','
              << row.nonpolar_kj_mol << ',' << row.total_kj_mol << ','
              << row.total_kj_mol * kKJToKcal << "\n";
      }
    } else {
      TitanUI::print_warning(
          "PBSA: failed to write mm_pbsa_residue_frame.csv.");
    }
  }
}

void Analyzer::write_pbsa_by_trajectory() const {
  if (!args_.pbsa)
    return;
  struct TrajStats {
    std::size_t frames{0};
    std::size_t failed{0};
    ScalarStats elec;
    ScalarStats vdw;
    ScalarStats mm;
    ScalarStats polar;
    ScalarStats nonpolar;
    ScalarStats no_entropy;
    std::vector<double> mm_values;
  };
  std::map<std::size_t, TrajStats> stats_by_traj;
  for (const FrameResult &row : pbsa_rows_) {
    TrajStats &stats = stats_by_traj[row.trajectory_id];
    ++stats.frames;
    if (row.pbsa_status != "ok")
      ++stats.failed;
    add_scalar(stats.elec, row.pbsa_elec_kj_mol);
    add_scalar(stats.vdw, row.pbsa_vdw_kj_mol);
    add_scalar(stats.mm, row.pbsa_mm_kj_mol);
    add_scalar(stats.polar, row.pbsa_delta_solv_kj_mol);
    add_scalar(stats.nonpolar, row.pbsa_nonpolar_kj_mol);
    add_scalar(stats.no_entropy, row.pbsa_no_entropy_kj_mol);
    if (args_.pbsa_ie && row.pbsa_status == "ok")
      stats.mm_values.push_back(row.pbsa_mm_kj_mol);
  }

  std::ofstream out("mm_pbsa_by_trajectory.csv",
                    std::ios::out | std::ios::trunc);
  if (!out.is_open()) {
    TitanUI::print_warning("PBSA: failed to write mm_pbsa_by_trajectory.csv.");
    return;
  }
  out << "trajectory_id,frames,failed,E_elec_kJmol_mean,E_elec_kJmol_std,"
         "E_vdw_kJmol_mean,E_vdw_kJmol_std,E_MM_kJmol_mean,"
         "E_MM_kJmol_std,DeltaDeltaG_solv_PB_kJmol_mean,"
         "DeltaDeltaG_solv_PB_kJmol_std,DeltaG_nonpolar_SASA_kJmol_mean,"
         "DeltaG_nonpolar_SASA_kJmol_std,DeltaG_bind_noIE_kJmol,"
         "DeltaG_bind_noIE_kcalmol,IE_status,IE_samples,IE_ESS,"
         "minus_TDeltaS_IE_kJmol,minus_TDeltaS_IE_kcalmol,"
         "DeltaG_bind_IE_kJmol,DeltaG_bind_IE_kcalmol\n";
  for (const auto &[trajectory_id, stats] : stats_by_traj) {
    const PbsaEntropyResult entropy =
        compute_pbsa_entropy_for_values(stats.mm_values);
    const double minus_tds =
        entropy.status == "ok" ? entropy.minus_tds_kj_mol : 0.0;
    const double dg_ie = stats.no_entropy.mean + minus_tds;
    out << trajectory_id << ',' << stats.frames << ',' << stats.failed << ','
        << format_fixed(stats.elec.mean, 10) << ','
        << format_fixed(scalar_stddev(stats.elec), 10) << ','
        << format_fixed(stats.vdw.mean, 10) << ','
        << format_fixed(scalar_stddev(stats.vdw), 10) << ','
        << format_fixed(stats.mm.mean, 10) << ','
        << format_fixed(scalar_stddev(stats.mm), 10) << ','
        << format_fixed(stats.polar.mean, 10) << ','
        << format_fixed(scalar_stddev(stats.polar), 10) << ','
        << format_fixed(stats.nonpolar.mean, 10) << ','
        << format_fixed(scalar_stddev(stats.nonpolar), 10) << ','
        << format_fixed(stats.no_entropy.mean, 10) << ','
        << format_fixed(stats.no_entropy.mean * kKJToKcal, 10) << ','
        << entropy.status << ',' << entropy.samples << ','
        << format_fixed(entropy.ess, 4) << ',' << format_fixed(minus_tds, 10)
        << ',' << format_fixed(minus_tds * kKJToKcal, 10) << ','
        << format_fixed(dg_ie, 10) << ','
        << format_fixed(dg_ie * kKJToKcal, 10) << "\n";
  }
}

void Analyzer::write_pbsa_ie_diagnostics(
    const PbsaEntropyResult &global_entropy) const {
  if (!args_.pbsa)
    return;
  std::ofstream out("mm_pbsa_ie_diagnostics.csv",
                    std::ios::out | std::ios::trunc);
  if (!out.is_open()) {
    TitanUI::print_warning("PBSA: failed to write mm_pbsa_ie_diagnostics.csv.");
    return;
  }
  out << "scope,trajectory_id,block_index,samples,filtered_samples,ess,"
         "status,block_size,block_count,minus_TDeltaS_IE_kJmol,"
         "minus_TDeltaS_IE_kcalmol,block_mean_kJmol,block_std_kJmol,"
         "block_drift_kJmol\n";
  auto write_entropy = [&](const char *scope, std::size_t traj,
                           int block_index,
                           const PbsaEntropyResult &entropy) {
    out << scope << ',' << traj << ',' << block_index << ','
        << entropy.samples << ',' << entropy.filtered_samples << ','
        << format_fixed(entropy.ess, 4) << ',' << entropy.status << ','
        << entropy.block_size << ',' << entropy.block_count << ','
        << format_fixed(entropy.minus_tds_kj_mol, 10) << ','
        << format_fixed(entropy.minus_tds_kj_mol * kKJToKcal, 10) << ','
        << format_fixed(entropy.block_mean_kj_mol, 10) << ','
        << format_fixed(entropy.block_std_kj_mol, 10) << ','
        << format_fixed(entropy.block_drift_kj_mol, 10) << "\n";
  };
  write_entropy("global", 0, -1, global_entropy);

  std::map<std::size_t, std::vector<double>> values_by_traj;
  for (const FrameResult &row : pbsa_rows_) {
    if (args_.pbsa_ie && row.pbsa_status == "ok")
      values_by_traj[row.trajectory_id].push_back(row.pbsa_mm_kj_mol);
  }
  for (const auto &[trajectory_id, values] : values_by_traj)
    write_entropy("trajectory", trajectory_id, -1,
                  compute_pbsa_entropy_for_values(values));

  const std::size_t block_size =
      global_entropy.block_size > 0 ? global_entropy.block_size : 1;
  for (std::size_t begin = 0, block_index = 0;
       begin < pbsa_ie_mm_energies_kj_.size();
       begin += block_size, ++block_index) {
    const std::size_t end =
        std::min(pbsa_ie_mm_energies_kj_.size(), begin + block_size);
    std::vector<double> block(
        pbsa_ie_mm_energies_kj_.begin() + static_cast<std::ptrdiff_t>(begin),
        pbsa_ie_mm_energies_kj_.begin() + static_cast<std::ptrdiff_t>(end));
    write_entropy("block", 0, static_cast<int>(block_index),
                  compute_pbsa_entropy_for_values(block));
  }
}

void Analyzer::write_pbsa_performance() const {
  if (!args_.pbsa || args_.pbsa_perf == "off")
    return;
  std::ofstream out("mm_pbsa_performance.dat", std::ios::out | std::ios::trunc);
  if (!out.is_open()) {
    TitanUI::print_warning("PBSA: failed to write mm_pbsa_performance.dat.");
    return;
  }
  auto add_time = [](ScalarStats &stats, double value) {
    Analyzer::add_scalar(stats, value);
  };
  ScalarStats complex_t;
  ScalarStats receptor_t;
  ScalarStats ligand_t;
  ScalarStats components_t;
  ScalarStats residue_t;
  ScalarStats frame_t;
  double slowest_ms = 0.0;
  std::size_t slowest_global = 0;
  for (const FrameResult &row : pbsa_rows_) {
    add_time(complex_t, row.pbsa_time_complex_ms);
    add_time(receptor_t, row.pbsa_time_receptor_ms);
    add_time(ligand_t, row.pbsa_time_ligand_ms);
    add_time(components_t, row.pbsa_time_components_ms);
    add_time(residue_t, row.pbsa_time_residue_ms);
    add_time(frame_t, row.pbsa_time_frame_ms);
    if (row.pbsa_time_frame_ms > slowest_ms) {
      slowest_ms = row.pbsa_time_frame_ms;
      slowest_global = row.global_index;
    }
  }
  auto write_time = [&](const char *name, const ScalarStats &stats) {
    out << name << "_ms_mean " << format_fixed(stats.mean, 6) << "\n";
    out << name << "_ms_max " << format_fixed(stats.max, 6) << "\n";
  };
  out << "# TITAN MM/PBSA performance summary\n";
  out << "frames " << pbsa_rows_.size() << "\n";
  out << "endpoint_solves " << (pbsa_rows_.size() * 6) << "\n";
  out << "max_grid_points " << max_grid_points_seen_ << "\n";
  out << "max_memory_bytes " << max_memory_bytes_seen_ << "\n";
  out << "residue_decomp " << args_.pbsa_decomp << "\n";
  write_time("complex_solve", complex_t);
  write_time("receptor_solve", receptor_t);
  write_time("ligand_solve", ligand_t);
  write_time("components", components_t);
  write_time("residue_decomp", residue_t);
  write_time("frame_total", frame_t);
  out << "slowest_global_frame " << slowest_global << "\n";
  out << "slowest_frame_ms " << format_fixed(slowest_ms, 6) << "\n";

  if (args_.pbsa_perf == "detail") {
    std::ofstream detail("mm_pbsa_performance_frames.csv",
                         std::ios::out | std::ios::trunc);
    if (detail.is_open()) {
      detail << "trajectory_id,frame_index,global_index,time_ps,"
                "complex_solve_ms,receptor_solve_ms,ligand_solve_ms,"
                "components_ms,residue_decomp_ms,frame_total_ms,"
                "grid_points,memory_bytes,status\n";
      detail << std::setprecision(10);
      for (const FrameResult &row : pbsa_rows_) {
        detail << row.trajectory_id << ',' << row.frame_index << ','
               << row.global_index << ',' << row.time_ps << ','
               << row.pbsa_time_complex_ms << ','
               << row.pbsa_time_receptor_ms << ','
               << row.pbsa_time_ligand_ms << ','
               << row.pbsa_time_components_ms << ','
               << row.pbsa_time_residue_ms << ','
               << row.pbsa_time_frame_ms << ',' << row.grid_points << ','
               << row.memory_bytes << ',' << row.pbsa_status << "\n";
      }
    } else {
      TitanUI::print_warning(
          "PBSA: failed to write mm_pbsa_performance_frames.csv.");
    }
  }
}

void Analyzer::write_pbsa_summary(std::fstream &log) {
  if (!args_.pbsa)
    return;

  const PbsaEntropyResult entropy = compute_pbsa_entropy();
  write_pbsa_energy_csv(entropy);
  write_pbsa_decomposition(entropy);
  write_pbsa_residue_outputs();
  write_pbsa_by_trajectory();
  write_pbsa_ie_diagnostics(entropy);
  write_pbsa_performance();

  const double dg_no_entropy_kj = pbsa_no_entropy_stats_.mean;
  const double minus_tds_ie_kj =
      (args_.pbsa_ie && entropy.status == "ok") ? entropy.minus_tds_kj_mol
                                                 : 0.0;
  const double dg_binding_ie_kj = dg_no_entropy_kj + minus_tds_ie_kj;

  std::ofstream out("mm_pbsa_summary.dat", std::ios::out | std::ios::trunc);
  if (!out.is_open()) {
    TitanUI::print_warning("PBSA: failed to write mm_pbsa_summary.dat.");
    return;
  }

  auto write_energy_stats = [&](const char *name, const ScalarStats &stats) {
    out << name << "_kJmol_mean " << format_fixed(stats.mean) << "\n";
    out << name << "_kJmol_std " << format_fixed(scalar_stddev(stats))
        << "\n";
    out << name << "_kJmol_min " << format_fixed(stats.min) << "\n";
    out << name << "_kJmol_max " << format_fixed(stats.max) << "\n";
    out << name << "_kcalmol_mean " << format_fixed(stats.mean * kKJToKcal)
        << "\n";
    out << name << "_kcalmol_std "
        << format_fixed(scalar_stddev(stats) * kKJToKcal) << "\n";
    out << name << "_kcalmol_min " << format_fixed(stats.min * kKJToKcal)
        << "\n";
    out << name << "_kcalmol_max " << format_fixed(stats.max * kKJToKcal)
        << "\n";
  };

  out << "# TITAN MM/PBSA summary\n";
  out << "model MM/PBSA_PBE_linear\n";
  out << "pbsa_enabled true\n";
  out << "pbsa_components_requested " << args_.pbsa_components << "\n";
  out << "pbsa_decomp_requested " << args_.pbsa_decomp << "\n";
  out << "frames_total " << pbsa_frames_ << "\n";
  out << "frames_failed " << pbsa_failed_ << "\n";
  out << "temperature_K " << format_fixed(args_.temp) << "\n";
  out << "pbe_salt_M " << format_fixed(args_.salt) << "\n";
  out << "pbe_surface " << args_.surface << "\n";
  out << "sasa_probe_A " << format_fixed(args_.probe_radius) << "\n";
  out << "sasa_points_per_atom " << TitanSASA::kDefaultPointCount << "\n";
  out << "nonpolar_gamma_kJmol_A2 "
      << format_fixed(TitanSASA::kNonpolarGammaKJPerMolA2, 8) << "\n";
  out << "nonpolar_gamma_kcalmol_A2 "
      << format_fixed(TitanSASA::kNonpolarGammaKcalPerMolA2, 8) << "\n";
  out << "endpoint_solves_per_frame 6\n";
  out << "energy_timeseries mm_pbsa_energy.csv\n";
  if (args_.pbsa_decomp != "off")
    out << "decomposition mm_pbsa_decomposition.dat\n";
  out << "interaction_entropy_enabled " << (args_.pbsa_ie ? "true" : "false")
      << "\n";
  out << "interaction_entropy_status " << entropy.status << "\n";
  out << "interaction_entropy_min_samples " << args_.pbsa_ie_min_samples
      << "\n";
  out << "interaction_entropy_filter_sigma "
      << format_fixed(args_.pbsa_ie_filter_sigma) << "\n";
  out << "interaction_entropy_block_config " << args_.pbsa_ie_block << "\n";
  out << "interaction_entropy_source MM_elec_plus_vdw\n";
  out << "interaction_entropy_samples " << entropy.samples
      << "\n";
  out << "interaction_entropy_filtered_samples " << entropy.filtered_samples
      << "\n";
  out << "interaction_entropy_ESS " << format_fixed(entropy.ess, 4) << "\n";
  out << "interaction_entropy_block_size " << entropy.block_size << "\n";
  out << "interaction_entropy_block_count " << entropy.block_count << "\n";
  out << "interaction_entropy_block_mean_kJmol "
      << format_fixed(entropy.block_mean_kj_mol) << "\n";
  out << "interaction_entropy_block_std_kJmol "
      << format_fixed(entropy.block_std_kj_mol) << "\n";
  out << "interaction_entropy_block_drift_kJmol "
      << format_fixed(entropy.block_drift_kj_mol) << "\n";
  out << "interaction_entropy_diagnostics mm_pbsa_ie_diagnostics.csv\n";
  out << "per_trajectory_summary mm_pbsa_by_trajectory.csv\n";
  if (args_.pbsa_perf != "off")
    out << "performance_summary mm_pbsa_performance.dat\n";
  if (args_.pbsa_decomp == "residue") {
    out << "residue_summary mm_pbsa_residue_summary.csv\n";
    out << "residue_by_trajectory mm_pbsa_residue_by_trajectory.csv\n";
    out << "residue_closure_total_max_abs_kJmol "
        << format_fixed(pbsa_residue_closure_max_abs_kj_, 10) << "\n";
  }
  write_energy_stats("E_elec", pbsa_elec_stats_);
  write_energy_stats("E_vdw", pbsa_vdw_stats_);
  write_energy_stats("E_MM", pbsa_mm_stats_);
  write_energy_stats("G_polar_PB", pbsa_polar_stats_);
  write_energy_stats("DeltaDeltaG_solv_PB", pbsa_polar_stats_);
  write_energy_stats("G_nonpolar_SASA", pbsa_nonpolar_stats_);
  out << "deltaSASA_A2_mean " << format_fixed(pbsa_delta_sasa_stats_.mean)
      << "\n";
  out << "deltaSASA_A2_std "
      << format_fixed(scalar_stddev(pbsa_delta_sasa_stats_)) << "\n";
  out << "deltaSASA_A2_min " << format_fixed(pbsa_delta_sasa_stats_.min)
      << "\n";
  out << "deltaSASA_A2_max " << format_fixed(pbsa_delta_sasa_stats_.max)
      << "\n";
  write_energy_stats("DeltaG_no_entropy", pbsa_no_entropy_stats_);
  out << "minus_TdS_IE_kJmol " << format_fixed(minus_tds_ie_kj) << "\n";
  out << "minus_TdS_IE_kcalmol "
      << format_fixed(minus_tds_ie_kj * kKJToKcal) << "\n";
  out << "DeltaG_binding_IE_kJmol " << format_fixed(dg_binding_ie_kj)
      << "\n";
  out << "DeltaG_binding_IE_kcalmol "
      << format_fixed(dg_binding_ie_kj * kKJToKcal) << "\n";
  std::map<std::size_t, std::vector<const FrameResult *>> rows_by_traj;
  for (const FrameResult &row : pbsa_rows_)
    rows_by_traj[row.trajectory_id].push_back(&row);
  ScalarStats traj_equal_noie;
  ScalarStats traj_equal_ie;
  for (const auto &[trajectory_id, rows] : rows_by_traj) {
    (void)trajectory_id;
    ScalarStats noie;
    std::vector<double> mm_values;
    for (const FrameResult *row : rows) {
      add_scalar(noie, row->pbsa_no_entropy_kj_mol);
      if (args_.pbsa_ie && row->pbsa_status == "ok")
        mm_values.push_back(row->pbsa_mm_kj_mol);
    }
    const PbsaEntropyResult traj_entropy =
        compute_pbsa_entropy_for_values(mm_values);
    const double traj_minus_tds =
        traj_entropy.status == "ok" ? traj_entropy.minus_tds_kj_mol : 0.0;
    add_scalar(traj_equal_noie, noie.mean);
    add_scalar(traj_equal_ie, noie.mean + traj_minus_tds);
  }
  out << "trajectory_equal_DeltaG_no_entropy_kJmol_mean "
      << format_fixed(traj_equal_noie.mean) << "\n";
  out << "trajectory_equal_DeltaG_no_entropy_kJmol_std "
      << format_fixed(scalar_stddev(traj_equal_noie)) << "\n";
  out << "trajectory_equal_DeltaG_binding_IE_kJmol_mean "
      << format_fixed(traj_equal_ie.mean) << "\n";
  out << "trajectory_equal_DeltaG_binding_IE_kJmol_std "
      << format_fixed(scalar_stddev(traj_equal_ie)) << "\n";

  if (log.is_open()) {
    TitanUI::log_section(log, "MM/PBSA OUTPUT AUDIT");
    TitanUI::log_kv(log, "Frames", std::to_string(pbsa_frames_));
    TitanUI::log_kv(log, "Failed", std::to_string(pbsa_failed_));
    TitanUI::log_kv(log, "Endpoint solves/frame", "6");
    TitanUI::log_kv(log, "IE Status", entropy.status);
    TitanUI::log_kv(log, "Interaction Entropy",
                    format_fixed(minus_tds_ie_kj) + " kJ/mol / " +
                        format_fixed(minus_tds_ie_kj * kKJToKcal) +
                        " kcal/mol");
    TitanUI::log_kv(log, "DeltaG Binding IE",
                    format_fixed(dg_binding_ie_kj) + " kJ/mol / " +
                        format_fixed(dg_binding_ie_kj * kKJToKcal) +
                        " kcal/mol");
    TitanUI::log_kv(log, "Summary", "mm_pbsa_summary.dat");
    TitanUI::log_kv(log, "Energy CSV", "mm_pbsa_energy.csv");
    TitanUI::log_kv(log, "By Trajectory", "mm_pbsa_by_trajectory.csv");
    TitanUI::log_kv(log, "IE Diagnostics", "mm_pbsa_ie_diagnostics.csv");
    if (args_.pbsa_perf != "off")
      TitanUI::log_kv(log, "Performance", "mm_pbsa_performance.dat");
    if (args_.pbsa_decomp != "off")
      TitanUI::log_kv(log, "Decomposition", "mm_pbsa_decomposition.dat");
    if (args_.pbsa_decomp == "residue")
      TitanUI::log_kv(log, "Residue Summary",
                      "mm_pbsa_residue_summary.csv");
  }
}

void Analyzer::finalize(std::fstream &log) {
  if (!enabled())
    return;
  if (energy_out_.is_open())
    energy_out_.flush();

  std::ofstream summary("pbe_summary.dat", std::ios::out | std::ios::trunc);
  if (summary.is_open()) {
    const double variance =
        frames_total_ > 1 ? energy_m2_ / static_cast<double>(frames_total_ - 1)
                          : 0.0;
    const double energy_kcal_mean = energy_mean_ / kKcalToKJ;
    const double energy_kcal_std = std::sqrt(variance) / kKcalToKJ;
    summary << "# TITAN PBE v1 summary\n";
    summary << "model " << (args_.model == "nonlinear" ? "nonlinear_pb"
                                                        : "linear_pb")
            << "\n";
    summary << "pbe_model " << args_.model << "\n";
    summary << "pbe_model_effective "
            << (args_.model == "nonlinear" ? "nonlinear_or_linearized"
                                            : "linear")
            << "\n";
    summary << "pbe_biomodel " << args_.biomodel << "\n";
    summary << "dielectric_model " << args_.dielectric_model << "\n";
    summary << "surface " << args_.surface << "\n";
    summary << "probe_radius_A " << format_fixed(args_.probe_radius) << "\n";
    summary << "surface_effective " << args_.surface << "\n";
    summary << "face_refine "
            << (surface_uses_face_refinement(args_.surface) ? "true" : "false")
            << "\n";
    summary << "boundary_mode " << args_.boundary << "\n";
    summary << "charge_assignment " << args_.charge_assignment << "\n";
    summary << "radius_convention " << effective_radius_convention() << "\n";
    summary << "gaussian_exp " << format_fixed(args_.gaussian_exp, 6) << "\n";
    summary << "gaussian_sigma " << format_fixed(args_.gaussian_sigma, 6)
            << "\n";
    summary << "density_cutoff " << format_fixed(args_.density_cutoff, 6)
            << "\n";
    summary << "surface_cutoff "
            << (args_.surface_cutoff == "auto"
                    ? std::string("auto")
                    : format_fixed(args_.surface_cutoff_value, 6))
            << "\n";
    summary << "gap_diel " << format_fixed(args_.gap_diel, 6) << "\n";
    summary << "ion_access " << args_.ion_access << "\n";
    summary << "ion_radius_A " << format_fixed(args_.ion_radius) << "\n";
    summary << "precision cpu_double\n";
    summary << "scale_points_per_A " << format_fixed(args_.scale) << "\n";
    summary << "margin_A " << format_fixed(args_.margin) << "\n";
    summary << "gridbox " << args_.gridbox << "\n";
    summary << "int_diel " << format_fixed(args_.int_diel) << "\n";
    summary << "ext_diel " << format_fixed(args_.ext_diel) << "\n";
    summary << "salt_M " << format_fixed(args_.salt) << "\n";
    summary << "temperature_K " << format_fixed(args_.temp) << "\n";
    summary << "kernel_requested " << args_.kernel << "\n";
    summary << "kernel_effective scalar_threaded\n";
    summary << "focus_mode " << args_.focus << "\n";
    summary << "multigrid_requested " << args_.multigrid << "\n";
    summary << "multigrid_effective " << effective_solver_mode() << "\n";
    summary << "focus_levels_effective " << effective_focus_levels() << "\n";
    summary << "focus_margin_extra_A " << format_fixed(args_.focus_margin)
            << "\n";
    summary << "mg_levels_config " << args_.mg_levels << "\n";
    summary << "mg_pre " << args_.mg_pre << "\n";
    summary << "mg_post " << args_.mg_post << "\n";
    summary << "mg_coarse_iter " << args_.mg_coarse_iter << "\n";
    summary << "nl_max_iter " << args_.nl_max_iter << "\n";
    summary << "nl_inner_iter " << args_.nl_inner_iter << "\n";
    summary << "nl_stagnation " << args_.nl_stagnation << "\n";
    summary << "nl_solver " << args_.nl_solver << "\n";
    summary << "nl_inner " << args_.nl_inner << "\n";
    summary << "nl_multigrid " << args_.nl_multigrid << "\n";
    summary << "nl_block " << args_.nl_block << "\n";
    summary << "nl_dphi_tol " << format_fixed(args_.nl_dphi, 10) << "\n";
    summary << "nl_relax " << format_fixed(effective_nl_relax(), 6) << "\n";
    summary << "nl_clip " << format_fixed(args_.nl_clip, 6) << "\n";
    summary << "nl_continuation " << args_.nl_continuation << "\n";
    summary << "nl_line_search " << args_.nl_line_search << "\n";
    summary << "nl_residual_tol " << format_fixed(args_.nl_residual, 10)
            << "\n";
    if (!args_.test_pqr.empty()) {
      summary << "test_pqr_adapter " << args_.test_pqr << "\n";
      summary << "test_grid_size " << args_.test_grid_size << "\n";
      summary << "test_center "
              << (args_.test_center_user
                      ? format_fixed(args_.test_center_x, 6) + "," +
                            format_fixed(args_.test_center_y, 6) + "," +
                            format_fixed(args_.test_center_z, 6)
                      : std::string("auto"))
              << "\n";
    }
    summary << "frames_total " << frames_total_ << "\n";
    summary << "frames_converged " << frames_converged_ << "\n";
    summary << "frames_failed " << frames_failed_ << "\n";
    summary << "energy_kJmol_mean " << format_fixed(energy_mean_) << "\n";
    summary << "energy_kJmol_std " << format_fixed(std::sqrt(variance)) << "\n";
    summary << "energy_kJmol_min " << format_fixed(energy_min_) << "\n";
    summary << "energy_kJmol_max " << format_fixed(energy_max_) << "\n";
    summary << "energy_kcalmol_mean " << format_fixed(energy_kcal_mean)
            << "\n";
    summary << "energy_kcalmol_std " << format_fixed(energy_kcal_std) << "\n";
    summary << "energy_kcalmol_min " << format_fixed(energy_min_ / kKcalToKJ)
            << "\n";
    summary << "energy_kcalmol_max " << format_fixed(energy_max_ / kKcalToKJ)
            << "\n";
    summary << "deltaG_solv_enabled " << (args_.dgsolv ? "true" : "false")
            << "\n";
    summary << "deltaG_solv_frames " << dgsolv_frames_ << "\n";
    if (dgsolv_frames_ > 0) {
      const double dg_variance =
          dgsolv_frames_ > 1
              ? dgsolv_m2_ / static_cast<double>(dgsolv_frames_ - 1)
              : 0.0;
      summary << "deltaG_solv_kJmol_mean " << format_fixed(dgsolv_mean_)
              << "\n";
      summary << "deltaG_solv_kJmol_std "
              << format_fixed(std::sqrt(dg_variance)) << "\n";
      summary << "deltaG_solv_kJmol_min " << format_fixed(dgsolv_min_)
              << "\n";
      summary << "deltaG_solv_kJmol_max " << format_fixed(dgsolv_max_)
              << "\n";
      summary << "deltaG_solv_kcalmol_mean "
              << format_fixed(dgsolv_mean_ / kKcalToKJ) << "\n";
      summary << "deltaG_solv_kcalmol_std "
              << format_fixed(std::sqrt(dg_variance) / kKcalToKJ) << "\n";
      summary << "deltaG_solv_kcalmol_min "
              << format_fixed(dgsolv_min_ / kKcalToKJ) << "\n";
      summary << "deltaG_solv_kcalmol_max "
              << format_fixed(dgsolv_max_ / kKcalToKJ) << "\n";
    }
    summary << "iterations_solvent_mean " << format_fixed(iter_mean_, 3)
            << "\n";
    summary << "max_grid_points " << max_grid_points_seen_ << "\n";
    summary << "max_memory_bytes " << max_memory_bytes_seen_ << "\n";
    summary << "max_solute_points " << max_solute_points_seen_ << "\n";
    summary << "max_solvent_points " << max_solvent_points_seen_ << "\n";
    summary << "max_ion_accessible_points " << max_ion_accessible_points_seen_
            << "\n";
    summary << "max_stern_excluded_points " << max_stern_excluded_points_seen_
            << "\n";
    summary << "max_focus_levels " << max_focus_levels_seen_ << "\n";
    summary << "max_mg_levels " << max_mg_levels_seen_ << "\n";
    summary << "mg_cycles_mean " << format_fixed(mg_cycle_mean_, 3) << "\n";
    summary << "max_mg_cycles " << max_mg_cycles_seen_ << "\n";
    summary << "max_boundary_rms " << format_fixed(max_boundary_rms_seen_, 10)
            << "\n";
    summary << "max_boundary_max " << format_fixed(max_boundary_max_seen_, 10)
            << "\n";
    summary << "max_linear_residual "
            << format_fixed(max_linear_residual_seen_, 10) << "\n";
    summary << "nonlinear_frames " << nonlinear_frames_ << "\n";
    summary << "nonlinear_converged " << nonlinear_converged_ << "\n";
    summary << "nonlinear_failed " << nonlinear_failed_ << "\n";
    summary << "max_nl_blocks " << max_nl_blocks_seen_ << "\n";
    summary << "max_nl_iterations " << max_nl_iterations_seen_ << "\n";
    summary << "max_nl_dphi " << format_fixed(max_nl_dphi_seen_, 10) << "\n";
    summary << "max_nl_residual " << format_fixed(max_nl_residual_seen_, 10)
            << "\n";
    summary << "max_nl_clipped_points " << max_nl_clipped_points_seen_
            << "\n";
    summary << "max_newton_inner_iterations "
            << max_newton_inner_iterations_seen_ << "\n";
    summary << "max_newton_step_norm "
            << format_fixed(max_newton_step_norm_seen_, 10) << "\n";
    summary << "max_newton_residual_before "
            << format_fixed(max_newton_residual_before_seen_, 10) << "\n";
    summary << "max_newton_residual_after "
            << format_fixed(max_newton_residual_after_seen_, 10) << "\n";
    summary << "max_newton_correction_residual "
            << format_fixed(max_newton_correction_residual_seen_, 10) << "\n";
    summary << "max_newton_correction_max_residual "
            << format_fixed(max_newton_correction_max_residual_seen_, 10)
            << "\n";
    summary << "max_newton_rejected_steps "
            << max_newton_rejected_steps_seen_ << "\n";
    summary << "max_newton_inner_levels "
            << max_newton_inner_levels_seen_ << "\n";
    summary << "max_newton_inner_cycles "
            << max_newton_inner_cycles_seen_ << "\n";
    summary << "max_newton_inner_residual_before "
            << format_fixed(max_newton_inner_residual_before_seen_, 10)
            << "\n";
    summary << "max_newton_inner_residual_after "
            << format_fixed(max_newton_inner_residual_after_seen_, 10)
            << "\n";
    summary << "energy_csv pbe_energy.csv\n";
    if (!args_.test_pqr.empty() && args_.parity_dump != "off") {
      summary << "parity_dump " << args_.parity_dump << "\n";
      summary << "test_energy_convention " << args_.test_energy_convention
              << "\n";
      summary << "test_midpoint_epsilon " << args_.test_midpoint_epsilon
              << "\n";
      summary << "parity_stage_csv pbe_parity_stage.csv\n";
      if (args_.parity_dump == "points" || args_.parity_dump == "all")
        summary << "parity_points_csv pbe_parity_points.csv\n";
      if (args_.parity_dump == "all" &&
          (args_.test_energy_convention == "delphi-iscm" ||
           args_.test_energy_convention == "delphi-corrected"))
        summary << "iscm_surface_charge_samples_csv "
                   "pbe_iscm_surface_charge_samples.csv\n";
    }
    if (!args_.test_pqr.empty() && args_.parity_map_dump != "off") {
      summary << "parity_map_dump " << args_.parity_map_dump << "\n";
      summary << "test_charge_grid " << args_.test_charge_grid << "\n";
      summary << "test_surface_map " << args_.test_surface_map << "\n";
      summary << "test_boundary_map " << args_.test_boundary_map << "\n";
      summary << "field_samples_csv pbe_field_samples.csv\n";
      summary << "boundary_samples_csv pbe_boundary_samples.csv\n";
      summary << "charge_grid_audit_csv pbe_charge_grid_audit.csv\n";
    }
  } else {
    TitanUI::print_warning("PBE: failed to write pbe_summary.dat.");
  }

  write_pbsa_summary(log);

  if (log.is_open()) {
    TitanUI::log_section(log, "PBE OUTPUT AUDIT");
    TitanUI::log_kv(log, "Frames", std::to_string(frames_total_));
    TitanUI::log_kv(log, "Converged", std::to_string(frames_converged_));
    TitanUI::log_kv(log, "Failed/Not Converged", std::to_string(frames_failed_));
    TitanUI::log_kv(log, "PBE Model", args_.model);
    TitanUI::log_kv(log, "Bio Model", args_.biomodel);
    TitanUI::log_kv(log, "Dielectric Model", args_.dielectric_model);
    TitanUI::log_kv(log, "Radius Convention", effective_radius_convention());
    TitanUI::log_kv(log, "Boundary", args_.boundary);
    TitanUI::log_kv(log, "Charge Assignment", args_.charge_assignment);
    TitanUI::log_kv(log, "Nonlinear Solver", args_.nl_solver);
    TitanUI::log_kv(log, "Nonlinear Inner", args_.nl_inner);
    TitanUI::log_kv(log, "Nonlinear MG", args_.nl_multigrid);
    TitanUI::log_kv(log, "Newton Inner Iter Cap",
                    std::to_string(args_.nl_inner_iter));
    TitanUI::log_kv(log, "Newton Stagnation Window",
                    std::to_string(args_.nl_stagnation));
    TitanUI::log_kv(log, "Nonlinear Frames",
                    std::to_string(nonlinear_frames_));
    TitanUI::log_kv(log, "Nonlinear Converged",
                    std::to_string(nonlinear_converged_));
    TitanUI::log_kv(log, "Nonlinear Failed",
                    std::to_string(nonlinear_failed_));
    TitanUI::log_kv(log, "Max Nonlinear dPhi",
                    format_fixed(max_nl_dphi_seen_, 10));
    TitanUI::log_kv(log, "Max Nonlinear Residual",
                    format_fixed(max_nl_residual_seen_, 10));
    TitanUI::log_kv(log, "Max Newton Inner Iterations",
                    std::to_string(max_newton_inner_iterations_seen_));
    TitanUI::log_kv(log, "Max Newton Residual Before/After",
                    format_fixed(max_newton_residual_before_seen_, 10) +
                        " / " +
                        format_fixed(max_newton_residual_after_seen_, 10));
    TitanUI::log_kv(log, "Max Newton Correction Residual",
                    format_fixed(max_newton_correction_residual_seen_, 10));
    TitanUI::log_kv(log, "Max Newton Correction Max Residual",
                    format_fixed(max_newton_correction_max_residual_seen_, 10));
    TitanUI::log_kv(log, "Max Newton Rejected Steps",
                    std::to_string(max_newton_rejected_steps_seen_));
    TitanUI::log_kv(log, "Max Newton Inner Levels/Cycles",
                    std::to_string(max_newton_inner_levels_seen_) + " / " +
                        std::to_string(max_newton_inner_cycles_seen_));
    TitanUI::log_kv(log, "Max Newton Inner Residual Before/After",
                    format_fixed(max_newton_inner_residual_before_seen_, 10) +
                        " / " +
                        format_fixed(max_newton_inner_residual_after_seen_, 10));
    TitanUI::log_kv(log, "Max Linear Residual",
                    format_fixed(max_linear_residual_seen_, 10));
    TitanUI::log_kv(log, "Max Boundary Potential",
                    format_fixed(max_boundary_max_seen_, 10));
    TitanUI::log_kv(log, "Max Nonlinear Clipped Points",
                    std::to_string(max_nl_clipped_points_seen_));
    TitanUI::log_kv(log, "Surface", args_.surface);
    if (uses_gaussian_model()) {
      TitanUI::log_kv(log, "Gaussian Exp", format_fixed(args_.gaussian_exp));
      TitanUI::log_kv(log, "Gaussian Sigma",
                      format_fixed(args_.gaussian_sigma));
      TitanUI::log_kv(log, "Density Cutoff",
                      format_fixed(args_.density_cutoff));
      TitanUI::log_kv(log, "Gap Dielectric", format_fixed(args_.gap_diel));
    }
    TitanUI::log_kv(log, "Probe Radius",
                    format_fixed(args_.probe_radius) + " A");
    TitanUI::log_kv(log, "Face Refinement",
                    surface_uses_face_refinement(args_.surface) ? "enabled"
                                                                : "off");
    TitanUI::log_kv(log, "Ion Access", args_.ion_access);
    TitanUI::log_kv(log, "Ion Radius", format_fixed(args_.ion_radius) + " A");
    TitanUI::log_kv(log, "Focus Mode", args_.focus);
    TitanUI::log_kv(log, "Multigrid Requested", args_.multigrid);
    TitanUI::log_kv(log, "Multigrid Effective", effective_solver_mode());
    TitanUI::log_kv(log, "Focus Levels",
                    std::to_string(effective_focus_levels()));
    TitanUI::log_kv(log, "MG Cycles Mean", format_fixed(mg_cycle_mean_, 3));
    TitanUI::log_kv(log, "DeltaG Solv", args_.dgsolv ? "enabled" : "off");
    if (dgsolv_frames_ > 0) {
      TitanUI::log_kv(log, "DeltaG Solv Frames",
                      std::to_string(dgsolv_frames_));
      TitanUI::log_kv(log, "DeltaG Solv Mean",
                      format_fixed(dgsolv_mean_) + " kJ/mol / " +
                          format_fixed(dgsolv_mean_ / kKcalToKJ) +
                          " kcal/mol");
    }
    if (!args_.test_pqr.empty()) {
      TitanUI::log_kv(log, "Test PQR Adapter", args_.test_pqr);
      TitanUI::log_kv(log, "Parity Dump", args_.parity_dump);
      TitanUI::log_kv(log, "Test Energy Convention",
                      args_.test_energy_convention);
      if (args_.parity_dump != "off") {
        TitanUI::log_kv(log, "Parity Stage CSV", "pbe_parity_stage.csv");
        if (args_.parity_dump == "points" || args_.parity_dump == "all")
          TitanUI::log_kv(log, "Parity Points CSV", "pbe_parity_points.csv");
        if (args_.parity_dump == "all" &&
            (args_.test_energy_convention == "delphi-iscm" ||
             args_.test_energy_convention == "delphi-corrected"))
          TitanUI::log_kv(log, "ISCM Surface Samples CSV",
                          "pbe_iscm_surface_charge_samples.csv");
      }
      TitanUI::log_kv(log, "Parity Map Dump", args_.parity_map_dump);
      TitanUI::log_kv(log, "Charge Grid Replay", args_.test_charge_grid);
      TitanUI::log_kv(log, "Surface Map Replay", args_.test_surface_map);
      TitanUI::log_kv(log, "Midpoint Epsilon Replay",
                      args_.test_midpoint_epsilon);
      TitanUI::log_kv(log, "Boundary Map Replay", args_.test_boundary_map);
      if (args_.parity_map_dump != "off") {
        TitanUI::log_kv(log, "Field Samples CSV", "pbe_field_samples.csv");
        TitanUI::log_kv(log, "Boundary Samples CSV",
                        "pbe_boundary_samples.csv");
        TitanUI::log_kv(log, "Charge Grid Audit CSV",
                        "pbe_charge_grid_audit.csv");
      }
    }
    TitanUI::log_kv(log, "Energy CSV", "pbe_energy.csv");
    TitanUI::log_kv(log, "Summary", "pbe_summary.dat");
    TitanUI::log_kv(log, "Max Grid Points",
                    std::to_string(max_grid_points_seen_));
    TitanUI::log_kv(log, "Max Memory Bytes",
                    std::to_string(max_memory_bytes_seen_));
    TitanUI::log_kv(log, "Max Solute Points",
                    std::to_string(max_solute_points_seen_));
    TitanUI::log_kv(log, "Max Solvent Points",
                    std::to_string(max_solvent_points_seen_));
    TitanUI::log_kv(log, "Max Ion Accessible Points",
                    std::to_string(max_ion_accessible_points_seen_));
    TitanUI::log_kv(log, "Max Stern Excluded Points",
                    std::to_string(max_stern_excluded_points_seen_));
  }
}

} // namespace TitanPBE
