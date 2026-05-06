#include "analyzer_impl.h"

#include "../titan_executor.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace TitanPBE {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kDelphiEpk = 167100.9162872952;

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

double delphi_debye_length_a(double salt_m, double temp_k,
                             double eps_out) noexcept {
  if (salt_m <= 0.0)
    return 1.0e6;
  return 0.01990076478 *
         std::sqrt(std::max(1.0, temp_k) * std::max(1.0e-12, eps_out)) /
         std::sqrt(salt_m);
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

std::string read_oracle_fixture_filter() {
  std::ifstream in("pbe_oracle_fixture.txt");
  if (!in.is_open())
    return {};
  std::string line;
  if (!std::getline(in, line))
    return {};
  return trim_copy(line);
}

bool oracle_fixture_matches(
    const std::vector<std::string> &fields,
    const std::unordered_map<std::string, std::size_t> &header,
    const std::string &fixture_filter) {
  if (fixture_filter.empty())
    return true;
  const auto it = header.find("fixture");
  if (it == header.end() || it->second >= fields.size())
    return true;
  return trim_copy(fields[it->second]) == fixture_filter;
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

} // namespace

double AnalyzerImpl::delphi_epkt() const noexcept {
  return kDelphiEpk / std::max(1.0, args_.temp);
}

bool AnalyzerImpl::delphi_nwt_parity_active() const noexcept {
  return !args_.test_pqr.empty() && args_.model == "nonlinear" &&
         args_.nl_solver == "newton" && args_.surface == "delphi-vdw" &&
         args_.charge_assignment == "delphi";
}

bool AnalyzerImpl::delphi_field_units_active() const noexcept {
  return args_.test_field_units == "delphi-dimensionless" ||
         (args_.test_operator_convention == "delphi" &&
          !args_.test_pqr.empty()) ||
         args_.focus == "delphi" ||
         args_.biomodel == "rpbe" ||
         delphi_nwt_parity_active();
}

double AnalyzerImpl::field_to_native_factor() const noexcept {
  return delphi_field_units_active() ? 1.0 / std::max(1.0e-30, delphi_epkt())
                                     : 1.0;
}

double AnalyzerImpl::effective_rhs_scale(const GridShape &shape) const noexcept {
  const double scale = rhs_scale(shape);
  return delphi_field_units_active() ? scale * delphi_epkt() : scale;
}

double AnalyzerImpl::effective_screen_term(const GridShape &shape, double eps_out,
                                       double kappa,
                                       std::size_t id) const noexcept {
  if (id >= ion_accessible_.size() || !ion_accessible_[id] ||
      args_.salt <= 0.0)
    return 0.0;
  if ((args_.test_operator_convention == "delphi" && !args_.test_pqr.empty()) ||
      args_.focus == "delphi" || delphi_nwt_parity_active()) {
    const double debye = delphi_debye_length_a(args_.salt, args_.temp, eps_out);
    if (debye <= 0.0 || !std::isfinite(debye))
      return 0.0;
    return eps_out * shape.h * shape.h / (debye * debye);
  }
  return kappa * kappa * shape.h * shape.h;
}

double AnalyzerImpl::delphi_relaxation_omega(const GridShape &shape) const {
  if (shape.n == 0 || shape.nx < 3 || shape.ny < 3 || shape.nz < 3 ||
      eps_x_.size() != shape.n || eps_y_.size() != shape.n ||
      eps_z_.size() != shape.n || denom_.size() != shape.n) {
    return 1.0;
  }
  const std::size_t nx = shape.nx;
  const std::size_t ny = shape.ny;
  const std::size_t nz = shape.nz;
  const std::size_t y_stride = nz;
  const std::size_t x_stride = ny * y_stride;
  const std::size_t half_n = (shape.n + 1) / 2;
  std::vector<double> even(half_n, 0.0);
  std::vector<double> odd(half_n, 0.0);
  std::vector<double> current(shape.n, 0.0);
  std::vector<double> sin_x(nx, 0.0), sin_y(ny, 0.0), sin_z(nz, 0.0);
  const double norm_x = std::sqrt(std::max(1.0, static_cast<double>(nx - 1)));
  const double norm_y = std::sqrt(std::max(1.0, static_cast<double>(ny - 1)));
  const double norm_z = std::sqrt(std::max(1.0, static_cast<double>(nz - 1)));
  const double sqrt2 = std::sqrt(2.0);
  for (std::size_t i = 0; i < nx; ++i)
    sin_x[i] = sqrt2 * std::sin(kPi * static_cast<double>(i) /
                                static_cast<double>(nx - 1)) /
               norm_x;
  for (std::size_t j = 0; j < ny; ++j)
    sin_y[j] = sqrt2 * std::sin(kPi * static_cast<double>(j) /
                                static_cast<double>(ny - 1)) /
               norm_y;
  for (std::size_t k = 0; k < nz; ++k)
    sin_z[k] = sqrt2 * std::sin(kPi * static_cast<double>(k) /
                                static_cast<double>(nz - 1)) /
               norm_z;
  for (std::size_t i = 0; i < nx; ++i) {
    for (std::size_t j = 0; j < ny; ++j) {
      for (std::size_t k = 0; k < nz; ++k) {
        const std::size_t flat = (i * ny + j) * nz + k;
        const double value = sin_x[i] * sin_y[j] * sin_z[k];
        current[flat] = value;
        if ((flat & 1ULL) == 0)
          even[flat / 2] = value;
        else
          odd[flat / 2] = value;
      }
    }
  }

  auto titan_id_from_flat = [&](std::size_t flat) {
    const std::size_t i = flat / x_stride;
    const std::size_t rem = flat - i * x_stride;
    const std::size_t j = rem / y_stride;
    const std::size_t k = rem - j * y_stride;
    return idx3(i, j, k, nx, ny);
  };
  auto half_value = [&](const std::vector<double> &values, std::size_t flat) {
    return values[flat / 2];
  };
  auto iterate_once = [&](int parity) {
    const std::vector<double> &neighbors = parity == 0 ? odd : even;
    std::vector<double> &target = parity == 0 ? even : odd;
    for (std::size_t flat = static_cast<std::size_t>(parity); flat < shape.n;
         flat += 2) {
      const std::size_t i = flat / x_stride;
      const std::size_t rem = flat - i * x_stride;
      const std::size_t j = rem / y_stride;
      const std::size_t k = rem - j * y_stride;
      if (is_boundary(i, j, k, shape))
        continue;
      const std::size_t id = titan_id_from_flat(flat);
      const std::size_t tid_xy = nx * ny;
      const double numerator =
          half_value(neighbors, flat - 1) * eps_z_[id - tid_xy] +
          half_value(neighbors, flat + 1) * eps_z_[id] +
          half_value(neighbors, flat - y_stride) * eps_y_[id - nx] +
          half_value(neighbors, flat + y_stride) * eps_y_[id] +
          half_value(neighbors, flat - x_stride) * eps_x_[id - 1] +
          half_value(neighbors, flat + x_stride) * eps_x_[id];
      target[flat / 2] = numerator / std::max(1.0e-30, denom_[id]);
    }
  };

  constexpr int kDelphiRelaxBlock = 1;
  for (int iter = 1; iter < kDelphiRelaxBlock; ++iter) {
    iterate_once(0);
    iterate_once(1);
  }
  for (std::size_t flat = 1; flat < shape.n; flat += 2)
    current[flat] = odd[flat / 2];
  iterate_once(0);
  iterate_once(1);

  double temp = 0.0;
  for (std::size_t flat = 1; flat < shape.n; flat += 2)
    temp += current[flat] * odd[flat / 2];
  double spectral_radius = 2.0 * temp;
  if (!std::isfinite(spectral_radius))
    return 1.0;
  spectral_radius = std::clamp(spectral_radius, 0.0, 1.0);
  return std::clamp(2.0 / (1.0 + std::sqrt(1.0 - spectral_radius)), 1.0,
                    1.95);
}

double AnalyzerImpl::debye_length_from_kappa(double kappa) const noexcept {
  if (kappa > 0.0)
    return 1.0 / kappa;
  return 1.0e6;
}

AnalyzerImpl::BoundaryFormulaSample AnalyzerImpl::boundary_formula_sample(
    const GridShape &shape, const std::vector<Atom> &atoms, std::size_t i,
    std::size_t j, std::size_t k, double eps_out, double kappa) const {
  BoundaryFormulaSample sample;
  const bool delphi_boundary = delphi_boundary_convention_active();
  sample.epkt = delphi_epkt();
  sample.debye_length =
      delphi_boundary
          ? delphi_debye_length_a(args_.salt, args_.temp, eps_out)
          : debye_length_from_kappa(kappa);
  sample.epsilon_temp = std::max(1.0e-30, eps_out / sample.epkt);
  const double x = shape.ox + static_cast<double>(i) * shape.h;
  const double y = shape.oy + static_cast<double>(j) * shape.h;
  const double z = shape.oz + static_cast<double>(k) * shape.h;
  double nearest_grid_dist2 = std::numeric_limits<double>::infinity();
  double nearest_real_dist2 = std::numeric_limits<double>::infinity();
  for (const Atom &atom : atoms) {
    if (std::abs(atom.q) <= 1.0e-13)
      continue;
    const double dx = x - atom.x;
    const double dy = y - atom.y;
    const double dz = z - atom.z;
    const double r2 = dx * dx + dy * dy + dz * dz;
    const double r = std::sqrt(r2);
    if (r > 1.0e-12) {
      sample.native_real_phi +=
          atom.q * std::exp(-kappa * r) / (std::max(1.0e-12, eps_out) * r);
      if (r2 < nearest_real_dist2)
        nearest_real_dist2 = r2;
    }

    const double gx = (atom.x - shape.ox) / std::max(1.0e-30, shape.h);
    const double gy = (atom.y - shape.oy) / std::max(1.0e-30, shape.h);
    const double gz = (atom.z - shape.oz) / std::max(1.0e-30, shape.h);
    const double gdx = static_cast<double>(i) - gx;
    const double gdy = static_cast<double>(j) - gy;
    const double gdz = static_cast<double>(k) - gz;
    const double grid_dist2 = gdx * gdx + gdy * gdy + gdz * gdz;
    const double distance_a = std::sqrt(grid_dist2) * shape.h;
    if (distance_a > 1.0e-12) {
      double q = atom.q;
      if (kappa > 0.0)
        q *= std::exp(-distance_a / sample.debye_length);
      sample.delphi_grid_phi_raw += q / (distance_a * sample.epsilon_temp);
      if (grid_dist2 < nearest_grid_dist2) {
        nearest_grid_dist2 = grid_dist2;
        sample.atom_grid_x = gx;
        sample.atom_grid_y = gy;
        sample.atom_grid_z = gz;
      }
    }
  }
  sample.nearest_atom_grid_distance =
      std::isfinite(nearest_grid_dist2) ? std::sqrt(nearest_grid_dist2) : 0.0;
  sample.nearest_atom_distance_a =
      std::isfinite(nearest_real_dist2) ? std::sqrt(nearest_real_dist2) : 0.0;
  const double native_from_delphi =
      sample.delphi_grid_phi_raw / std::max(1.0e-30, sample.epkt);
  if (delphi_boundary) {
    sample.selected_phi = delphi_field_units_active()
                              ? sample.delphi_grid_phi_raw
                              : native_from_delphi;
  } else {
    sample.selected_phi = delphi_field_units_active()
                              ? sample.native_real_phi * sample.epkt
                              : sample.native_real_phi;
  }
  return sample;
}

AnalyzerImpl::BoundaryDiagnostics
AnalyzerImpl::apply_oracle_boundary_map(const GridShape &shape, double eps_out,
                                    double kappa) {
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
  const bool want_vacuum =
      kappa <= 0.0 &&
      std::abs(eps_out - args_.int_diel) <=
          1.0e-8 * std::max(1.0, std::abs(args_.int_diel));
  auto row_matches_phase = [&](const std::vector<std::string> &fields) {
    const auto it = header.find("vacuum");
    if (it == header.end() || it->second >= fields.size())
      return true;
    const std::string value = trim_copy(fields[it->second]);
    if (value.empty() || value == "NA")
      return true;
    const bool row_vacuum =
        value == "1" || value == "true" || value == "True" ||
        value == "TRUE" || value == "vacuum";
    return row_vacuum == want_vacuum;
  };
  while (std::getline(in, line)) {
    if (trim_copy(line).empty())
      continue;
    const auto fields = split_csv_simple(line);
    if (!row_matches_phase(fields))
      continue;
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
    phi_[id] = delphi_field_units_active()
                   ? value
                   : value / std::max(1.0e-30, delphi_epkt());
    seen[id] = 1U;
    ++rows;
  }
  std::size_t expected = 0;
  double sumsq = 0.0;
  double max_abs = 0.0;
  if (rows == 0 && want_vacuum) {
    for (std::size_t kk = 0; kk < shape.nz; ++kk) {
      for (std::size_t jj = 0; jj < shape.ny; ++jj) {
        for (std::size_t ii = 0; ii < shape.nx; ++ii) {
          if (!is_boundary(ii, jj, kk, shape))
            continue;
          ++expected;
          const std::size_t id = idx3(ii, jj, kk, shape.nx, shape.ny);
          sumsq += phi_[id] * phi_[id];
          max_abs = std::max(max_abs, std::abs(phi_[id]));
        }
      }
    }
    diag.rms = expected > 0 ? std::sqrt(sumsq / static_cast<double>(expected))
                            : 0.0;
    diag.max_abs = max_abs;
    return diag;
  }
  expected = 0;
  sumsq = 0.0;
  max_abs = 0.0;
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

void AnalyzerImpl::maybe_write_boundary_diff(const GridShape &shape,
                                         const std::vector<Atom> &atoms,
                                         double eps_out, double kappa) const {
  if (args_.test_pqr.empty() || !oracle_replay_applies(shape))
    return;
  std::ifstream in("pbe_oracle_boundary_map.csv");
  if (!in.is_open())
    return;
  std::ofstream out;
  if (!boundary_diff_header_written_) {
    out.open("pbe_boundary_diff.csv", std::ios::out | std::ios::trunc);
    if (!out.is_open())
      throw std::runtime_error("PBE: failed to open pbe_boundary_diff.csv.");
    out << "i,j,k,flat_id,vacuum,oracle_phi_raw,oracle_phi_native,"
           "native_real_phi,delphi_grid_phi_raw,delphi_grid_phi_native,"
           "selected_phi,selected_phi_raw,selected_phi_native,"
           "delta_raw,delta_native,atom_grid_x,atom_grid_y,atom_grid_z,"
           "nearest_atom_grid_distance,nearest_atom_distance_A,epsilon_temp,"
           "debye_length,epkt,boundary_convention,field_units\n";
    boundary_diff_header_written_ = true;
  } else {
    out.open("pbe_boundary_diff.csv", std::ios::out | std::ios::app);
    if (!out.is_open())
      throw std::runtime_error("PBE: failed to append pbe_boundary_diff.csv.");
  }
  out << std::setprecision(12);
  std::string line;
  if (!std::getline(in, line))
    return;
  const auto header = csv_header_map(split_csv_simple(line));
  const bool want_vacuum =
      kappa <= 0.0 &&
      std::abs(eps_out - args_.int_diel) <=
          1.0e-8 * std::max(1.0, std::abs(args_.int_diel));
  auto row_matches_phase = [&](const std::vector<std::string> &fields) {
    const auto it = header.find("vacuum");
    if (it == header.end() || it->second >= fields.size())
      return true;
    const std::string value = trim_copy(fields[it->second]);
    if (value.empty() || value == "NA")
      return true;
    const bool row_vacuum =
        value == "1" || value == "true" || value == "True" ||
        value == "TRUE" || value == "vacuum";
    return row_vacuum == want_vacuum;
  };
  while (std::getline(in, line)) {
    if (trim_copy(line).empty())
      continue;
    const auto fields = split_csv_simple(line);
    if (!row_matches_phase(fields))
      continue;
    std::size_t i = 0, j = 0, k = 0;
    double oracle_raw = 0.0;
    if (!csv_index(fields, header, "i", i) ||
        !csv_index(fields, header, "j", j) ||
        !csv_index(fields, header, "k", k) ||
        !csv_double(fields, header, {"phi", "potential", "boundary_phi"},
                    oracle_raw))
      continue;
    if (i >= shape.nx || j >= shape.ny || k >= shape.nz ||
        !is_boundary(i, j, k, shape))
      continue;
    std::string vacuum_value = "NA";
    if (const auto it = header.find("vacuum"); it != header.end() &&
        it->second < fields.size())
      vacuum_value = fields[it->second];
    const BoundaryFormulaSample sample =
        boundary_formula_sample(shape, atoms, i, j, k, eps_out, kappa);
    const double epkt = std::max(1.0e-30, sample.epkt);
    const double oracle_native = oracle_raw / epkt;
    const double delphi_native = sample.delphi_grid_phi_raw / epkt;
    double selected_field = sample.selected_phi;
    double selected_raw =
        delphi_field_units_active() ? selected_field : selected_field * epkt;
    double selected_native = selected_field * field_to_native_factor();
    if (args_.test_boundary_map == "oracle" ||
        args_.test_boundary_convention == "oracle") {
      selected_raw = oracle_raw;
      selected_native = oracle_native;
      selected_field =
          delphi_field_units_active() ? oracle_raw : oracle_native;
    }
    out << i << ',' << j << ',' << k << ','
        << idx3(i, j, k, shape.nx, shape.ny) << ','
        << csv_escape(vacuum_value) << ',' << oracle_raw << ','
        << oracle_native << ',' << sample.native_real_phi << ','
        << sample.delphi_grid_phi_raw << ',' << delphi_native << ','
        << selected_field << ',' << selected_raw << ','
        << selected_native << ',' << (selected_raw - oracle_raw) << ','
        << (selected_native - oracle_native) << ',' << sample.atom_grid_x << ','
        << sample.atom_grid_y << ',' << sample.atom_grid_z << ','
        << sample.nearest_atom_grid_distance << ','
        << sample.nearest_atom_distance_a << ',' << sample.epsilon_temp << ','
        << sample.debye_length << ',' << sample.epkt << ','
        << csv_escape(args_.test_boundary_convention) << ','
        << csv_escape(args_.test_field_units) << '\n';
  }
}

AnalyzerImpl::BoundaryDiagnostics
AnalyzerImpl::set_boundary(const GridShape &shape, const std::vector<Atom> &atoms,
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
    for (std::size_t j = 0; j < shape.ny; ++j) {
      for (std::size_t i = 0; i < shape.nx; ++i) {
        if (!is_boundary(i, j, k, shape))
          continue;
        double potential = 0.0;
        if (args_.biomodel == "rpbe") {
          const bool vacuum_phase = eps_out <= args_.int_diel + 1.0e-9;
          const double x = shape.ox + static_cast<double>(i) * shape.h;
          const double y = shape.oy + static_cast<double>(j) * shape.h;
          const double z = shape.oz + static_cast<double>(k) * shape.h;
          const double epkt = delphi_epkt();
          const double indi_scaled =
              std::max(1.0e-30, args_.int_diel / std::max(1.0e-30, epkt));
          const double media_scaled =
              std::max(1.0e-30,
                       (vacuum_phase ? args_.int_diel : eps_out) /
                           std::max(1.0e-30, epkt));
          double coulomb_map = 0.0;
          double media_map = 0.0;
          for (const Atom &atom : atoms) {
            if (std::abs(atom.q) <= 0.0)
              continue;
            const double dx = x - atom.x;
            const double dy = y - atom.y;
            const double dz = z - atom.z;
            const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
            if (r <= 1.0e-10)
              continue;
            coulomb_map += atom.q / (r * indi_scaled);
            double q = atom.q;
            if (!vacuum_phase)
              q *= std::exp(-kappa * r);
            media_map += q / (r * media_scaled);
          }
          potential = media_map - coulomb_map;
        } else if (delphi_boundary_convention_active()) {
          potential =
              boundary_formula_sample(shape, atoms, i, j, k, eps_out, kappa)
                  .selected_phi;
        } else if (args_.boundary == "dipolar") {
          const double x = shape.ox + static_cast<double>(i) * shape.h;
          const double y = shape.oy + static_cast<double>(j) * shape.h;
          const double z = shape.oz + static_cast<double>(k) * shape.h;
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
          if (delphi_field_units_active())
            potential *= delphi_epkt();
        } else {
          const double x = shape.ox + static_cast<double>(i) * shape.h;
          const double y = shape.oy + static_cast<double>(j) * shape.h;
          const double z = shape.oz + static_cast<double>(k) * shape.h;
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
          if (delphi_field_units_active())
            potential *= delphi_epkt();
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

AnalyzerImpl::SolveStats AnalyzerImpl::solve_sor(const GridShape &shape,
                                          TitanExecutor &executor) {
  SolveStats stats;
  double omega =
      args_.omega == "auto"
          ? std::clamp(2.0 / (1.0 + std::sin(kPi / static_cast<double>(
                                             std::max({shape.nx, shape.ny,
                                                       shape.nz})))),
                       1.0, 1.95)
          : args_.omega_value;
  if (args_.omega == "auto" &&
      ((args_.test_operator_convention == "delphi" && !args_.test_pqr.empty()) ||
       args_.focus == "delphi")) {
    omega = delphi_relaxation_omega(shape);
  }
  const std::size_t interior_z = shape.nz > 2 ? shape.nz - 2 : 0;
  const std::size_t workers = effective_workers(executor, interior_z);
  const double rhs = effective_rhs_scale(shape);

  std::vector<double> local_sumsq(workers, 0.0);
  std::vector<double> local_max(workers, 0.0);
  std::vector<std::size_t> local_count(workers, 0);
  const bool write_trace =
      !args_.test_pqr.empty() && args_.parity_dump == "all";
  const bool delphi_block_stop =
      (!args_.test_pqr.empty() && args_.test_operator_convention == "delphi") ||
      args_.focus == "delphi";
  std::vector<std::size_t> trace_ids;
  std::unordered_map<std::string, double> oracle_trace_charge_density;
  auto trace_key = [](int color, std::size_t i, std::size_t j,
                      std::size_t k) {
    return std::to_string(color) + "|" + std::to_string(i) + "|" +
           std::to_string(j) + "|" + std::to_string(k);
  };
  if (write_trace && shape.nx > 2 && shape.ny > 2 && shape.nz > 2) {
    auto add_trace_id = [&](std::size_t i, std::size_t j, std::size_t k) {
      if (i == 0 || j == 0 || k == 0 || i + 1 >= shape.nx ||
          j + 1 >= shape.ny || k + 1 >= shape.nz)
        return;
      trace_ids.push_back(idx3(i, j, k, shape.nx, shape.ny));
    };
    {
      std::ifstream oracle_trace("pbe_oracle_sor_trace.csv");
      if (oracle_trace.is_open()) {
        std::string line;
        if (std::getline(oracle_trace, line)) {
          const auto header = csv_header_map(split_csv_simple(line));
          const std::string fixture_filter = read_oracle_fixture_filter();
          while (std::getline(oracle_trace, line) && trace_ids.size() < 64) {
            if (trim_copy(line).empty())
              continue;
            const auto fields = split_csv_simple(line);
            if (!oracle_fixture_matches(fields, header, fixture_filter))
              continue;
            std::size_t i = 0, j = 0, k = 0;
            if (csv_index(fields, header, "i", i) &&
                csv_index(fields, header, "j", j) &&
                csv_index(fields, header, "k", k)) {
              add_trace_id(i, j, k);
              double even_odd = 0.0;
              double charge_density = 0.0;
              if (csv_double(fields, header, {"even_odd", "color"},
                             even_odd) &&
                  csv_double(fields, header,
                             {"charge_density", "charge_map_scaled", "rhs"},
                             charge_density)) {
                oracle_trace_charge_density[trace_key(
                    static_cast<int>(std::llround(even_odd)), i, j, k)] =
                    charge_density;
              }
            }
          }
        }
      }
    }
    if (trace_ids.empty()) {
      trace_ids = {
          idx3(shape.nx / 2, shape.ny / 2, shape.nz / 2, shape.nx, shape.ny),
          idx3(1, 1, 1, shape.nx, shape.ny),
          idx3(shape.nx - 2, shape.ny - 2, shape.nz - 2, shape.nx, shape.ny),
          idx3(shape.nx / 2, shape.ny / 2, 1, shape.nx, shape.ny),
          idx3(shape.nx / 2, shape.ny / 2, shape.nz - 2, shape.nx,
               shape.ny)};
    }
    std::sort(trace_ids.begin(), trace_ids.end());
    trace_ids.erase(std::unique(trace_ids.begin(), trace_ids.end()),
                    trace_ids.end());
  }
  std::ofstream trace_out;
  if (write_trace) {
    if (!sor_trace_header_written_) {
      trace_out.open("pbe_sor_trace.csv", std::ios::out | std::ios::trunc);
      if (!trace_out.is_open())
        throw std::runtime_error("PBE: failed to open pbe_sor_trace.csv.");
      trace_out
          << "iteration,color,i,j,k,omega,old_phi,neighbor_sum,"
             "charge_density,denom,new_phi,dphi,rms,max_dphi,field_units,"
             "operator_convention,rhs_delta_vs_oracle,charge_grid_source\n";
      sor_trace_header_written_ = true;
    } else {
      trace_out.open("pbe_sor_trace.csv", std::ios::out | std::ios::app);
      if (!trace_out.is_open())
        throw std::runtime_error("PBE: failed to append pbe_sor_trace.csv.");
    }
    trace_out << std::setprecision(12);
  }
  for (int iter = 1; iter <= args_.max_iter; ++iter) {
    std::fill(local_sumsq.begin(), local_sumsq.end(), 0.0);
    std::fill(local_max.begin(), local_max.end(), 0.0);
    std::fill(local_count.begin(), local_count.end(), 0);

    for (int color = 0; color < 2; ++color) {
      std::vector<std::vector<std::string>> trace_lines(workers);
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
                  const double old_phi = phi_[id];
                  const double updated = phi_[id] + omega * (target - phi_[id]);
                  const double delta = updated - phi_[id];
                  phi_[id] = updated;
                  if (write_trace && iter <= 4 &&
                      std::binary_search(trace_ids.begin(), trace_ids.end(),
                                         id)) {
                    const auto oracle_charge_it =
                        oracle_trace_charge_density.find(trace_key(color, i, j, k));
                    const bool has_oracle_charge =
                        oracle_charge_it != oracle_trace_charge_density.end();
                    const double charge_density = rhs * rho_[id];
                    std::ostringstream row;
                    row << std::setprecision(12) << iter << ',' << color << ','
                        << i << ',' << j << ',' << k << ',' << omega << ','
                        << old_phi << ',' << neighbor_sum << ','
                        << charge_density << ',' << denom_[id] << ','
                        << updated << ',' << delta << ",0,0,"
                        << csv_escape(args_.test_field_units) << ','
                        << csv_escape(args_.test_operator_convention) << ',';
                    if (has_oracle_charge)
                      row << (charge_density - oracle_charge_it->second);
                    else
                      row << "NA";
                    row << ',' << csv_escape(args_.test_charge_grid) << '\n';
                    trace_lines[worker].push_back(row.str());
                  }
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
      if (write_trace) {
        for (const auto &worker_lines : trace_lines) {
          for (const std::string &row : worker_lines)
            trace_out << row;
        }
      }
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
    if (write_trace) {
      trace_out << iter << ",summary,0,0,0," << omega
                << ",0,0,0,0,0,0," << stats.rms << ',' << stats.max_dphi
                << ',' << csv_escape(args_.test_field_units) << ','
                << csv_escape(args_.test_operator_convention)
                << ",NA," << csv_escape(args_.test_charge_grid) << '\n';
    }
    const bool can_stop_this_iter =
        !delphi_block_stop || (iter % 50 == 0) || iter == args_.max_iter;
    if (can_stop_this_iter &&
        (stats.rms <= args_.rms_tol || stats.max_dphi <= args_.dphi_tol)) {
      stats.converged = true;
      break;
    }
  }
  return stats;
}

double AnalyzerImpl::compute_linear_residual(const GridShape &shape) const {
  if (shape.n == 0 || phi_.size() != shape.n || rho_.size() != shape.n)
    return 0.0;
  const double rhs = effective_rhs_scale(shape);
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

int AnalyzerImpl::effective_focus_levels() const noexcept {
  if (args_.focus != "auto" && args_.focus != "delphi")
    return 1;
  return std::clamp(args_.focus_levels <= 1 ? 2 : args_.focus_levels, 2, 4);
}

std::string AnalyzerImpl::effective_solver_mode() const {
  if (args_.multigrid == "vcycle")
    return "vcycle";
  if (args_.multigrid == "focus" || args_.multigrid == "auto" ||
      args_.focus == "auto")
    return "nested_focus";
  if (args_.focus == "delphi")
    return "delphi_focus";
  return "sor";
}

bool AnalyzerImpl::delphi_boundary_convention_active() const noexcept {
  return args_.test_boundary_convention == "delphi-grid" ||
         args_.test_boundary_convention == "oracle" ||
         delphi_nwt_parity_active() || args_.focus == "delphi";
}

std::string AnalyzerImpl::effective_boundary_convention() const {
  if (args_.test_boundary_convention == "oracle")
    return "oracle";
  if (delphi_boundary_convention_active())
    return "delphi-grid";
  return args_.test_boundary_convention;
}

AnalyzerImpl::FieldSolve AnalyzerImpl::solve_focused_field(
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

    auto build_delphi_focus_shape = [&]() {
      GridShape shape;
      shape.h = 1.0 / scale;
      std::size_t n = static_cast<std::size_t>(
          std::max(3, args_.test_grid_size));
      for (int r = 0; r < remaining; ++r)
        n = n * 2 + 7;
      if ((n & 1ULL) == 0)
        ++n;

      double cx = 0.0, cy = 0.0, cz = 0.0;
      if (level + 1 == levels) {
        cx = args_.test_center_user ? args_.test_center_x
                                    : args_.focus_center_x;
        cy = args_.test_center_user ? args_.test_center_y
                                    : args_.focus_center_y;
        cz = args_.test_center_user ? args_.test_center_z
                                    : args_.focus_center_z;
      } else {
        double min_x = std::numeric_limits<double>::infinity();
        double min_y = std::numeric_limits<double>::infinity();
        double min_z = std::numeric_limits<double>::infinity();
        double max_x = -std::numeric_limits<double>::infinity();
        double max_y = -std::numeric_limits<double>::infinity();
        double max_z = -std::numeric_limits<double>::infinity();
        for (const Atom &atom : atoms) {
          const double r = std::max(0.0, atom.radius);
          min_x = std::min(min_x, atom.x - r);
          min_y = std::min(min_y, atom.y - r);
          min_z = std::min(min_z, atom.z - r);
          max_x = std::max(max_x, atom.x + r);
          max_y = std::max(max_y, atom.y + r);
          max_z = std::max(max_z, atom.z + r);
        }
        cx = 0.5 * (min_x + max_x);
        cy = 0.5 * (min_y + max_y);
        cz = 0.5 * (min_z + max_z);
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
            << " points exceeds -pbe-max-grid " << args_.max_grid_points
            << ".";
        throw std::runtime_error(msg.str());
      }
      return shape;
    };

    const bool delphi_test_focus_shape =
        args_.focus == "delphi" && !args_.test_pqr.empty() &&
        args_.test_grid_size > 0 &&
        (args_.focus_center_user || args_.test_center_user);
    const GridShape shape =
        delphi_test_focus_shape ? build_delphi_focus_shape()
                                : build_shape(atoms, scale, margin);
    std::vector<Atom> focus_atoms;
    const std::vector<Atom> *level_atoms = &atoms;
    if (args_.focus == "delphi" && has_previous) {
      const double half_len =
          0.5 * static_cast<double>(std::max({shape.nx, shape.ny, shape.nz}) -
                                    1) *
          shape.h;
      const double cx = shape.ox + 0.5 * static_cast<double>(shape.nx - 1) *
                                       shape.h;
      const double cy = shape.oy + 0.5 * static_cast<double>(shape.ny - 1) *
                                       shape.h;
      const double cz = shape.oz + 0.5 * static_cast<double>(shape.nz - 1) *
                                       shape.h;
      constexpr double kDelphiFocusAtomMarginA = 3.5;
      const double lo_x = cx - half_len - kDelphiFocusAtomMarginA;
      const double lo_y = cy - half_len - kDelphiFocusAtomMarginA;
      const double lo_z = cz - half_len - kDelphiFocusAtomMarginA;
      const double hi_x = cx + half_len + kDelphiFocusAtomMarginA;
      const double hi_y = cy + half_len + kDelphiFocusAtomMarginA;
      const double hi_z = cz + half_len + kDelphiFocusAtomMarginA;
      focus_atoms.reserve(atoms.size());
      for (const Atom &atom : atoms) {
        if (atom.x >= lo_x && atom.x <= hi_x && atom.y >= lo_y &&
            atom.y <= hi_y && atom.z >= lo_z && atom.z <= hi_z) {
          focus_atoms.push_back(atom);
        }
      }
      if (!focus_atoms.empty())
        level_atoms = &focus_atoms;
    }
    const std::size_t peak_bytes =
        shape.memory_bytes + previous_phi.size() * sizeof(double);
    field.peak_memory_bytes = std::max(field.peak_memory_bytes, peak_bytes);
    ensure_workspace(shape);
    classify_surface(shape, *level_atoms, executor);
    apply_oracle_surface_map(shape);
    classify_ion_accessibility(shape, *level_atoms, executor);
    auto assign_focus_clipped_charges = [&]() {
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
      for (const Atom &atom : *level_atoms) {
        original_mx += atom.q * atom.x;
        original_my += atom.q * atom.y;
        original_mz += atom.q * atom.z;
        const double gx = (atom.x - shape.ox) / shape.h;
        const double gy = (atom.y - shape.oy) / shape.h;
        const double gz = (atom.z - shape.oz) / shape.h;
        const double focus_hi = static_cast<double>(
            std::max({shape.nx, shape.ny, shape.nz}) - 1);
        if (!(gx > 0.0 && gy > 0.0 && gz > 0.0 && gx < focus_hi &&
              gy < focus_hi && gz < focus_hi)) {
          continue;
        }
        const std::size_t i0 = std::min<std::size_t>(
            shape.nx - 2, static_cast<std::size_t>(std::floor(gx)));
        const std::size_t j0 = std::min<std::size_t>(
            shape.ny - 2, static_cast<std::size_t>(std::floor(gy)));
        const std::size_t k0 = std::min<std::size_t>(
            shape.nz - 2, static_cast<std::size_t>(std::floor(gz)));
        const double fx = clamp_fraction(gx - static_cast<double>(i0));
        const double fy = clamp_fraction(gy - static_cast<double>(j0));
        const double fz = clamp_fraction(gz - static_cast<double>(k0));
        for (std::size_t dz = 0; dz <= 1; ++dz) {
          const double wz = dz ? fz : (1.0 - fz);
          for (std::size_t dy = 0; dy <= 1; ++dy) {
            const double wy = dy ? fy : (1.0 - fy);
            for (std::size_t dx = 0; dx <= 1; ++dx) {
              const double wx = dx ? fx : (1.0 - fx);
              add_charge(i0 + dx, j0 + dy, k0 + dz,
                         atom.q * wx * wy * wz);
            }
          }
        }
      }
      const double dmx = assigned_mx - original_mx;
      const double dmy = assigned_my - original_my;
      const double dmz = assigned_mz - original_mz;
      diag.moment_error = std::sqrt(dmx * dmx + dmy * dmy + dmz * dmz);
      return diag;
    };
    const bool delphi_focus_charge_clip =
        args_.focus == "delphi" && has_previous &&
        args_.charge_assignment == "delphi";
    const ChargeDiagnostics charge_diag =
        args_.test_charge_grid == "oracle"
            ? assign_oracle_charge_grid(shape, *level_atoms)
            : (delphi_focus_charge_clip ? assign_focus_clipped_charges()
                                        : assign_charges(shape, *level_atoms));
    maybe_write_charge_grid_diff(shape, *level_atoms);
    AlignedDoubleVector energy_rho;
    build_coefficients(shape, *level_atoms, eps_in, eps_out, kappa, executor);
    if (args_.biomodel == "rpbe") {
      energy_rho = rho_;
      AlignedDoubleVector rpb_source;
      AlignedDoubleVector rpb_coulomb;
      (void)compute_rpb_source_audit(shape, *level_atoms, eps_x_, eps_y_, eps_z_,
                                     &rpb_source, &rpb_coulomb, eps_out);
      if (rpb_source.size() == shape.n) {
        const double source_to_native_rho =
            (shape.h * shape.h) /
            std::max(1.0e-30, effective_rhs_scale(shape));
        const bool solvent_phase = kappa > 0.0 && eps_out > eps_in + 1.0e-9;
        double kappa_square = 0.0;
        if (solvent_phase) {
          const double debye =
              delphi_debye_length_a(args_.salt, args_.temp, eps_out);
          if (debye > 0.0 && std::isfinite(debye))
            kappa_square = eps_out / (debye * debye);
        }
        for (std::size_t id = 0; id < shape.n; ++id)
          {
            double source_value = rpb_source[id];
            if (solvent_phase && rpb_coulomb.size() == shape.n &&
                last_gaussian_surface_.size() == shape.n) {
              source_value -=
                  (1.0 - std::clamp(last_gaussian_surface_[id], 0.0, 1.0)) *
                  kappa_square * rpb_coulomb[id];
            }
            rho_[id] = source_value * source_to_native_rho;
          }
      }
    }
    if (level + 1 == levels && map_result != nullptr)
      maybe_write_maps(shape, *map_result, false);
    BoundaryDiagnostics boundary_diag =
        set_boundary(shape, *level_atoms, eps_out, kappa);
    if (args_.test_boundary_map == "oracle" && oracle_replay_applies(shape))
      boundary_diag = apply_oracle_boundary_map(shape, eps_out, kappa);
    if (level + 1 == levels)
      maybe_write_boundary_diff(shape, *level_atoms, eps_out, kappa);
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
      if (args_.focus == "delphi" && map_result != nullptr && has_previous) {
        const std::string filename = "pbe_focus_boundary_samples.csv";
        const bool need_header =
            !std::filesystem::exists(filename) ||
            std::filesystem::file_size(filename) == 0;
        std::ofstream focus_out(filename, std::ios::out | std::ios::app);
        if (focus_out.is_open()) {
          if (need_header) {
            focus_out
                << "trajectory_id,frame_index,global_index,time_ps,i,j,k,"
                   "x_A,y_A,z_A,parent_phi,child_phi,parent_grid,child_grid,"
                   "interpolation_status\n";
          }
          std::size_t emitted = 0;
          constexpr std::size_t kMaxFocusBoundarySamples = 512;
          for (std::size_t k = 0;
               k < shape.nz && emitted < kMaxFocusBoundarySamples; ++k) {
            const double z = shape.oz + static_cast<double>(k) * shape.h;
            for (std::size_t j = 0;
                 j < shape.ny && emitted < kMaxFocusBoundarySamples; ++j) {
              const double y = shape.oy + static_cast<double>(j) * shape.h;
              for (std::size_t i = 0;
                   i < shape.nx && emitted < kMaxFocusBoundarySamples; ++i) {
                if (!is_boundary(i, j, k, shape))
                  continue;
                const double x = shape.ox + static_cast<double>(i) * shape.h;
                double parent_phi = 0.0;
                const bool ok = interpolate_field_value(
                    previous_shape, previous_phi, x, y, z, parent_phi);
                const std::size_t id = idx3(i, j, k, shape.nx, shape.ny);
                focus_out << map_result->trajectory_id << ','
                          << map_result->frame_index << ','
                          << map_result->global_index << ','
                          << std::setprecision(10) << map_result->time_ps
                          << ',' << i << ',' << j << ',' << k << ',' << x
                          << ',' << y << ',' << z << ',' << parent_phi << ','
                          << phi_[id] << ','
                          << csv_escape(field.focus_parent_grid) << ','
                          << csv_escape(field.focus_child_grid) << ','
                          << (ok ? "ok" : "outside_parent") << '\n';
                ++emitted;
              }
            }
          }
        }
      }
      field.boundary_convention = effective_boundary_convention();
      field.field_units = args_.test_field_units;
      field.rhs_scale = effective_rhs_scale(shape);
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
                   std::accumulate(level_atoms->begin(), level_atoms->end(), 0.0,
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
          !args_.probe_file.empty() ||
          args_.biomodel == "rpbe" ||
          barnase_forensic_enabled();
      if (args_.charge_assignment == "delphi" ||
          args_.charge_assignment == "gaussian" || keep_field_grids)
        field.phi_grid = phi_;
    if (keep_field_grids) {
      field.rho_grid = rho_;
      field.rho_energy_grid =
          (!energy_rho.empty() && energy_rho.size() == shape.n) ? energy_rho
                                                                : rho_;
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

AnalyzerImpl::FieldSolve AnalyzerImpl::solve_linear_field(
    const std::vector<Atom> &atoms, double eps_in, double eps_out, double kappa,
    TitanExecutor &executor, const FrameResult *map_result) {
  if (effective_solver_mode() == "vcycle")
    return solve_vcycle_field(atoms, eps_in, eps_out, kappa, executor,
                              map_result);
  return solve_focused_field(atoms, eps_in, eps_out, kappa, executor,
                             map_result);
}

} // namespace TitanPBE
