#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#include <psapi.h>
#else
#include <cstdio>
#endif

#include "analyzer_impl.h"

#include "titan_energy_decomposition.h"
#include "titan_executor.h"
#include "titan_interaction_entropy.h"
#include "titan_sasa.h"
#include "titan_simd_dispatch.h"
#include "titan_ui.h"
#include "tabi_types.h"

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
#include <thread>
#include <tuple>
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

struct ProcessMemorySample {
  std::size_t working_set_bytes{0};
  std::size_t peak_working_set_bytes{0};
};

ProcessMemorySample sample_process_memory() noexcept {
  ProcessMemorySample sample;
#if defined(_WIN32)
  PROCESS_MEMORY_COUNTERS counters{};
  counters.cb = sizeof(counters);
  if (K32GetProcessMemoryInfo(GetCurrentProcess(), &counters,
                              sizeof(counters))) {
    sample.working_set_bytes =
        static_cast<std::size_t>(counters.WorkingSetSize);
    sample.peak_working_set_bytes =
        static_cast<std::size_t>(counters.PeakWorkingSetSize);
  }
#else
  FILE *file = std::fopen("/proc/self/status", "r");
  if (file != nullptr) {
    char line[256];
    char label[64];
    long value_kb = 0;
    while (std::fgets(line, sizeof(line), file) != nullptr) {
      if (std::sscanf(line, "%63s %ld", label, &value_kb) != 2)
        continue;
      if (std::string(label) == "VmRSS:")
        sample.working_set_bytes =
            static_cast<std::size_t>(std::max<long>(0, value_kb)) * 1024u;
      else if (std::string(label) == "VmHWM:")
        sample.peak_working_set_bytes =
            static_cast<std::size_t>(std::max<long>(0, value_kb)) * 1024u;
    }
    std::fclose(file);
  }
#endif
  return sample;
}

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

double elapsed_ms(std::chrono::steady_clock::time_point begin,
                  std::chrono::steady_clock::time_point end) {
  return std::chrono::duration<double, std::milli>(end - begin).count();
}

struct NwtHomoOffsets {
  std::ptrdiff_t xp{0};
  std::ptrdiff_t xm{0};
  std::ptrdiff_t yp{0};
  std::ptrdiff_t ym{0};
  std::ptrdiff_t zp{0};
  std::ptrdiff_t zm{0};
};

NwtHomoOffsets nwt_homo_offsets(std::size_t nx_value,
                                std::size_t ny_value,
                                int color) noexcept {
  const std::ptrdiff_t nx = static_cast<std::ptrdiff_t>(nx_value);
  const std::ptrdiff_t xy =
      static_cast<std::ptrdiff_t>(nx_value * ny_value);
  if (color == 0) {
    return {0,
            -1,
            nx / 2,
            -((nx + 1) / 2),
            xy / 2,
            -((xy + 1) / 2)};
  }
  return {1, 0, (nx + 1) / 2, -(nx / 2), (xy + 1) / 2, -(xy / 2)};
}

#if TITAN_SIMD_X86
TITAN_TARGET_AVX2
void update_nwt_homo_linear_span_avx2(
    double *write_phi, const double *read_phi, const double *rhs,
    const double *eps_sum, const double *eps, std::uint32_t begin,
    std::uint32_t count, const NwtHomoOffsets offsets, double omega,
    bool accumulate_rms, double &sumsq, double &max_delta,
    std::size_t &sample_count) noexcept {
  const __m256d omega_v = _mm256_set1_pd(omega);
  const __m256d min_diag = _mm256_set1_pd(1.0e-12);
  const __m256d sign_mask = _mm256_set1_pd(-0.0);
  __m256d sumsq_v = _mm256_setzero_pd();
  __m256d max_v = _mm256_setzero_pd();

  std::size_t h = begin;
  const std::size_t end = static_cast<std::size_t>(begin) + count;
  for (; h + 4 <= end; h += 4) {
    const auto base = static_cast<std::ptrdiff_t>(h);
    const __m256d old_phi = _mm256_loadu_pd(write_phi + h);
    __m256d phi_sum =
        _mm256_add_pd(_mm256_loadu_pd(read_phi + base + offsets.xp),
                      _mm256_loadu_pd(read_phi + base + offsets.xm));
    phi_sum =
        _mm256_add_pd(phi_sum,
                      _mm256_loadu_pd(read_phi + base + offsets.yp));
    phi_sum =
        _mm256_add_pd(phi_sum,
                      _mm256_loadu_pd(read_phi + base + offsets.ym));
    phi_sum =
        _mm256_add_pd(phi_sum,
                      _mm256_loadu_pd(read_phi + base + offsets.zp));
    phi_sum =
        _mm256_add_pd(phi_sum,
                      _mm256_loadu_pd(read_phi + base + offsets.zm));
    const __m256d numerator =
        _mm256_add_pd(_mm256_mul_pd(_mm256_loadu_pd(eps + h), phi_sum),
                      _mm256_loadu_pd(rhs + h));
    const __m256d diagonal =
        _mm256_max_pd(_mm256_loadu_pd(eps_sum + h), min_diag);
    const __m256d target = _mm256_div_pd(numerator, diagonal);
    const __m256d delta =
        _mm256_mul_pd(omega_v, _mm256_sub_pd(target, old_phi));
    const __m256d updated = _mm256_add_pd(old_phi, delta);
    _mm256_storeu_pd(write_phi + h, updated);
    if (accumulate_rms) {
      sumsq_v = _mm256_add_pd(sumsq_v, _mm256_mul_pd(delta, delta));
      const __m256d abs_delta = _mm256_andnot_pd(sign_mask, delta);
      max_v = _mm256_max_pd(max_v, abs_delta);
    }
  }

  if (accumulate_rms) {
    alignas(32) double lanes[4];
    _mm256_store_pd(lanes, sumsq_v);
    sumsq += lanes[0] + lanes[1] + lanes[2] + lanes[3];
    _mm256_store_pd(lanes, max_v);
    max_delta = std::max(
        max_delta, std::max(std::max(lanes[0], lanes[1]),
                            std::max(lanes[2], lanes[3])));
    sample_count += h - static_cast<std::size_t>(begin);
  }

  for (; h < end; ++h) {
    const auto base = static_cast<std::ptrdiff_t>(h);
    const double old_phi = write_phi[h];
    const double neighbor_sum =
        eps[h] * (read_phi[base + offsets.xp] +
                  read_phi[base + offsets.xm] +
                  read_phi[base + offsets.yp] +
                  read_phi[base + offsets.ym] +
                  read_phi[base + offsets.zp] +
                  read_phi[base + offsets.zm]);
    const double target =
        (neighbor_sum + rhs[h]) / std::max(1.0e-12, eps_sum[h]);
    const double updated = old_phi + omega * (target - old_phi);
    const double delta = updated - old_phi;
    write_phi[h] = updated;
    if (accumulate_rms) {
      sumsq += delta * delta;
      max_delta = std::max(max_delta, std::abs(delta));
      ++sample_count;
    }
  }
}
#endif

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

class TabiPrepareProgressWriter {
public:
  TabiPrepareProgressWriter(const PbeArgs &args, std::size_t atom_count,
                            std::size_t threads)
      : enabled_(args.tabi_pre_solver_progress == "on"),
        start_(std::chrono::steady_clock::now()), atom_count_(atom_count),
        surface_points_(args.tabi_surface_points), threads_(threads),
        pqr_path_(args.test_pqr), tabi_sum_(args.tabi_sum),
        formulation_(args.tabi_formulation),
        tree_builder_(args.tabi_tree_builder),
        node_bounds_(args.tabi_node_bounds),
        mac_safety_(args.tabi_mac_safety),
        expansion_center_(args.tabi_expansion_center),
        fragmentation_safety_(args.tabi_fragmentation_safety),
        convergence_estimator_(args.tabi_convergence_estimator),
        convergence_estimator_model_(args.tabi_convergence_estimator_model),
        convergence_channel_calibration_(
            args.tabi_convergence_channel_calibration),
        convergence_split_policy_(args.tabi_convergence_split_policy),
        channel_missed_semantics_(args.tabi_channel_missed_semantics),
        progress_mode_(args.tabi_progress) {
    if (!enabled_)
      return;
    out_.open("pbe_tabi_prepare_progress.csv",
              std::ios::out | std::ios::trunc);
    if (!out_.is_open()) {
      enabled_ = false;
      return;
    }
    out_ << "elapsed_seconds,phase,atom_count,pqr_path,pqr_status,"
            "surface_points_requested,tabi_sum,tabi_formulation,"
            "tabi_tree_builder,tabi_node_bounds,tabi_mac_safety,"
            "tabi_expansion_center,tabi_fragmentation_safety,"
            "tabi_convergence_estimator,tabi_convergence_estimator_model,"
            "tabi_convergence_channel_calibration,"
            "tabi_convergence_split_policy,"
            "tabi_channel_missed_semantics,tabi_progress,threads,"
            "rss_bytes,status\n";
    out_.flush();
  }

  void write(const std::string &phase, const std::string &status) {
    if (!enabled_)
      return;
    const auto now = std::chrono::steady_clock::now();
    const ProcessMemorySample memory = sample_process_memory();
    const std::string pqr_status = pqr_path_.empty() ? "none" : "configured";
    out_ << std::fixed << std::setprecision(6)
         << std::chrono::duration<double>(now - start_).count() << ','
         << csv_escape(phase) << ',' << atom_count_ << ','
         << csv_escape(pqr_path_) << ',' << csv_escape(pqr_status) << ','
         << surface_points_ << ',' << csv_escape(tabi_sum_) << ','
         << csv_escape(formulation_) << ',' << csv_escape(tree_builder_)
         << ',' << csv_escape(node_bounds_) << ',' << csv_escape(mac_safety_)
         << ',' << csv_escape(expansion_center_) << ','
         << csv_escape(fragmentation_safety_) << ','
         << csv_escape(convergence_estimator_) << ','
         << csv_escape(convergence_estimator_model_) << ','
         << csv_escape(convergence_channel_calibration_) << ','
         << csv_escape(convergence_split_policy_) << ','
         << csv_escape(channel_missed_semantics_) << ','
         << csv_escape(progress_mode_) << ',' << threads_ << ','
         << memory.working_set_bytes << ',' << csv_escape(status) << '\n';
    out_.flush();
  }

private:
  bool enabled_{false};
  std::chrono::steady_clock::time_point start_{};
  std::size_t atom_count_{0};
  int surface_points_{96};
  std::size_t threads_{1};
  std::string pqr_path_;
  std::string tabi_sum_;
  std::string formulation_;
  std::string tree_builder_;
  std::string node_bounds_;
  std::string mac_safety_;
  std::string expansion_center_;
  std::string fragmentation_safety_;
  std::string convergence_estimator_;
  std::string convergence_estimator_model_;
  std::string convergence_channel_calibration_;
  std::string convergence_split_policy_;
  std::string channel_missed_semantics_;
  std::string progress_mode_;
  std::ofstream out_;
};

enum DelphiFaceAction : std::uint16_t {
  kDelphiFaceUnchanged = 0,
  kDelphiFaceRawVdw = 1,
  kDelphiFaceExposedSolvent = 2,
  kDelphiFaceContactSolvent = 3,
  kDelphiFaceCavityMedia = 4,
  kDelphiFaceSignRemap = 5,
  kDelphiFaceDelphiContactSign = 6,
  kDelphiFaceDelphiCavityRemap = 7,
  kDelphiFaceSignPreserve = 8,
};

const char *delphi_face_action_name(std::uint16_t action) noexcept {
  switch (action) {
  case kDelphiFaceRawVdw:
    return "raw_vdw";
  case kDelphiFaceExposedSolvent:
    return "exposed_solvent";
  case kDelphiFaceContactSolvent:
    return "contact_solvent";
  case kDelphiFaceCavityMedia:
    return "cavity_media";
  case kDelphiFaceSignRemap:
    return "sign_remap";
  case kDelphiFaceDelphiContactSign:
    return "delphi_contact_sign";
  case kDelphiFaceDelphiCavityRemap:
    return "delphi_cavity_remap";
  case kDelphiFaceSignPreserve:
    return "sign_preserve";
  default:
    return "unchanged";
  }
}

template <typename Shape>
bool is_boundary(std::size_t i, std::size_t j, std::size_t k,
                 const Shape &shape) noexcept {
  return i == 0 || j == 0 || k == 0 || i + 1 == shape.nx ||
         j + 1 == shape.ny || k + 1 == shape.nz;
}

std::pair<double, double> delphi_sinh_cosh_tail(double phi,
                                                double clip) noexcept {
  const double effective_clip = std::max(1.0e-12, clip);
  const double sinh_clip = std::sinh(effective_clip);
  if (phi > effective_clip) {
    const double sinh_phi = sinh_clip + sinh_clip * (phi - effective_clip);
    return {sinh_phi, std::sqrt(1.0 + sinh_phi * sinh_phi)};
  }
  if (phi < -effective_clip) {
    const double sinh_phi = -sinh_clip + sinh_clip * (phi + effective_clip);
    return {sinh_phi, std::sqrt(1.0 + sinh_phi * sinh_phi)};
  }
  return {std::sinh(phi), std::cosh(phi)};
}

std::pair<double, double>
delphi_sinh_cosh_tail_precomputed(double phi, double effective_clip,
                                  double sinh_clip) noexcept {
  if (phi > effective_clip) {
    const double sinh_phi = sinh_clip + sinh_clip * (phi - effective_clip);
    return {sinh_phi, std::sqrt(1.0 + sinh_phi * sinh_phi)};
  }
  if (phi < -effective_clip) {
    const double sinh_phi = -sinh_clip + sinh_clip * (phi + effective_clip);
    return {sinh_phi, std::sqrt(1.0 + sinh_phi * sinh_phi)};
  }
  return {std::sinh(phi), std::cosh(phi)};
}

std::pair<double, double>
delphi_sinh_cosh_tail_single_sinh(double phi, double effective_clip,
                                  double sinh_clip) noexcept {
  if (phi > effective_clip) {
    const double sinh_phi = sinh_clip + sinh_clip * (phi - effective_clip);
    return {sinh_phi, std::sqrt(1.0 + sinh_phi * sinh_phi)};
  }
  if (phi < -effective_clip) {
    const double sinh_phi = -sinh_clip + sinh_clip * (phi + effective_clip);
    return {sinh_phi, std::sqrt(1.0 + sinh_phi * sinh_phi)};
  }
  const double abs_phi = std::abs(phi);
  if (abs_phi <= 1.5) {
    const double x2 = phi * phi;
    const double sinh_poly =
        phi *
        (1.0 +
         x2 * (1.0 / 6.0 +
               x2 * (1.0 / 120.0 +
                     x2 * (1.0 / 5040.0 +
                           x2 * (1.0 / 362880.0 +
                                 x2 * (1.0 / 39916800.0 +
                                       x2 * (1.0 / 6227020800.0)))))));
    const double cosh_poly =
        1.0 +
        x2 * (0.5 +
              x2 * (1.0 / 24.0 +
                    x2 * (1.0 / 720.0 +
                          x2 * (1.0 / 40320.0 +
                                x2 * (1.0 / 3628800.0 +
                                      x2 * (1.0 / 479001600.0))))));
    return {sinh_poly, cosh_poly};
  }
  const double sinh_phi = std::sinh(phi);
  return {sinh_phi, std::sqrt(1.0 + sinh_phi * sinh_phi)};
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

bool surface_uses_face_refinement(const std::string &surface) noexcept {
  return surface == "delphi-ses" || surface == "delphi-vdw" ||
         surface == "gaussian" || surface == "gaussian-cutoff";
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

void AnalyzerImpl::configure(const arguments &args) {
  args_ = args.pbe;
  args_.pbsa_mol1_type = args.mol_1_type;
  args_.pbsa_mol2_type = args.mol_2_type;
  configured_threads_ = std::max<std::size_t>(1, args.nt);
  reset_grid_anchor();
  last_stability_state_ = StabilityFrameState{};
  last_stability_medium_state_ = StabilityMediumState{};
  last_stability_ensemble_ = StabilityEnsembleState{};
  stability_hysteresis_cells_kept_frame_ = 0;
  stability_hysteresis_cells_released_frame_ = 0;
  if (args_.pbsa) {
    args_.flag_pbe = true;
    args_.dgsolv = true;
  }
}

void AnalyzerImpl::begin_trajectory(std::size_t trajectory_id,
                                const std::string &source_label) {
  if (!enabled())
    return;
  current_trajectory_id_ = trajectory_id;
  current_source_label_ = source_label;
  reset_grid_anchor();
  last_stability_state_ = StabilityFrameState{};
  last_stability_medium_state_ = StabilityMediumState{};
  last_stability_ensemble_ = StabilityEnsembleState{};
  stability_hysteresis_cells_kept_frame_ = 0;
  stability_hysteresis_cells_released_frame_ = 0;
}

void AnalyzerImpl::end_trajectory(std::fstream &log) {
  if (!enabled())
    return;
  if (log.is_open()) {
    TitanUI::log_kv(log, "PBE Trajectory Complete",
                    std::to_string(current_trajectory_id_) + " " +
                        current_source_label_);
  }
}

void AnalyzerImpl::collect_atoms(const MOL &mol, const Multi_MOL_parameters &params,
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

void AnalyzerImpl::collect_test_pqr_atoms(std::vector<Atom> &atoms) const {
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

AnalyzerImpl::DelphiVdwMediaMap AnalyzerImpl::build_delphi_vdw_media_map(
    const GridShape &shape, const std::vector<Atom> &surface_atoms,
    TitanExecutor &executor) const {
  DelphiVdwMediaMap out;
  out.face_encoded.assign(shape.n * 3, 0);
  out.face_media.assign(shape.n * 3, 0);
  out.face_owner.assign(shape.n * 3, 0);
  out.face_debug.assign(shape.n * 3, {});
  out.audit.enabled = "raw";
  if (shape.n == 0 || surface_atoms.empty())
    return out;

  const std::int32_t delphi_epsilon_dim =
      static_cast<std::int32_t>(surface_atoms.size() + 2);
  auto face_index = [](std::size_t id, int axis) {
    return id * 3 + static_cast<std::size_t>(axis);
  };
  auto decode_media_from_value = [&](std::int32_t value) -> std::int32_t {
    return std::abs(value) / delphi_epsilon_dim;
  };
  auto decode_entity_from_value = [&](std::int32_t value) -> std::int32_t {
    if (value == 0)
      return 0;
    return value % delphi_epsilon_dim;
  };
  auto encode_delphi_face = [&](std::int32_t entity,
                                std::int32_t media) -> std::int32_t {
    const std::int32_t sign = entity < 0 ? -1 : 1;
    return sign * (std::abs(entity) + std::abs(media) * delphi_epsilon_dim);
  };
  auto sync_face_from_encoded = [&](std::size_t fid) {
    if (fid >= out.face_encoded.size())
      return;
    const std::int32_t value = out.face_encoded[fid];
    out.face_media[fid] = decode_media_from_value(value);
    out.face_owner[fid] = decode_entity_from_value(value);
    if (fid < out.face_debug.size())
      out.face_debug[fid].final_encoded = value;
  };
  auto set_face_encoded = [&](std::size_t fid, std::int32_t value) {
    if (fid >= out.face_encoded.size())
      return;
    out.face_encoded[fid] = value;
    sync_face_from_encoded(fid);
  };

  constexpr std::array<std::array<int, 3>, 27> kDelphiNeighborOffsets{{
      {{0, 0, 0}},   {{-1, 0, 0}},  {{1, 0, 0}},   {{0, -1, 0}},
      {{0, 1, 0}},   {{0, 0, -1}},  {{0, 0, 1}},   {{1, 0, 1}},
      {{-1, 0, 1}},  {{0, 1, 1}},   {{0, -1, 1}},  {{-1, -1, 0}},
      {{1, -1, 0}},  {{1, 1, 0}},   {{-1, 1, 0}},  {{-1, 0, -1}},
      {{1, 0, -1}},  {{0, 1, -1}},  {{0, -1, -1}}, {{-1, -1, -1}},
      {{1, -1, -1}}, {{1, 1, -1}},  {{-1, 1, -1}}, {{-1, 1, 1}},
      {{1, 1, 1}},   {{1, -1, 1}},  {{-1, -1, 1}},
  }};
  const double grid_spacing_extension = 0.5 * shape.h * std::sqrt(3.0);
  const double delphi_voxel_side =
      max_atom_radius(surface_atoms) + args_.ion_radius + grid_spacing_extension;
  double atom_min_x = surface_atoms.front().x;
  double atom_min_y = surface_atoms.front().y;
  double atom_min_z = surface_atoms.front().z;
  double atom_max_x = surface_atoms.front().x;
  double atom_max_y = surface_atoms.front().y;
  double atom_max_z = surface_atoms.front().z;
  for (const Atom &atom : surface_atoms) {
    atom_min_x = std::min(atom_min_x, atom.x);
    atom_min_y = std::min(atom_min_y, atom.y);
    atom_min_z = std::min(atom_min_z, atom.z);
    atom_max_x = std::max(atom_max_x, atom.x);
    atom_max_y = std::max(atom_max_y, atom.y);
    atom_max_z = std::max(atom_max_z, atom.z);
  }
  const double delphi_voxel_padding = delphi_voxel_side + 0.1;
  const double raw_voxel_origin_x = atom_min_x - delphi_voxel_padding;
  const double raw_voxel_origin_y = atom_min_y - delphi_voxel_padding;
  const double raw_voxel_origin_z = atom_min_z - delphi_voxel_padding;
  const auto raw_voxel_count = [&](double lo, double hi) -> std::size_t {
    return static_cast<std::size_t>(
        std::max<double>(
            0.0, std::floor((hi + delphi_voxel_padding -
                             (lo - delphi_voxel_padding)) /
                            delphi_voxel_side)));
  };
  const std::size_t raw_vx = raw_voxel_count(atom_min_x, atom_max_x);
  const std::size_t raw_vy = raw_voxel_count(atom_min_y, atom_max_y);
  const std::size_t raw_vz = raw_voxel_count(atom_min_z, atom_max_z);
  const std::size_t raw_sx = raw_vx + 1;
  const std::size_t raw_sy = raw_vy + 1;
  const std::size_t raw_sz = raw_vz + 1;
  const double raw_voxel_scale = 1.0 / std::max(1.0e-12, delphi_voxel_side);
  auto raw_voxel_flat = [&](std::size_t vx, std::size_t vy,
                            std::size_t vz) {
    return (vz * raw_sy + vy) * raw_sx + vx;
  };
  auto raw_atom_central = [&](const Atom &atom, std::size_t &vx,
                              std::size_t &vy, std::size_t &vz) {
    auto coord = [&](double value, double origin,
                     std::size_t high) -> std::size_t {
      const double scaled = (value - origin) * raw_voxel_scale;
      if (scaled <= 0.0)
        return 0;
      return std::min<std::size_t>(high, static_cast<std::size_t>(scaled));
    };
    vx = coord(atom.x, raw_voxel_origin_x, raw_vx);
    vy = coord(atom.y, raw_voxel_origin_y, raw_vy);
    vz = coord(atom.z, raw_voxel_origin_z, raw_vz);
  };
  std::vector<std::vector<std::size_t>> raw_voxel_atoms(raw_sx * raw_sy * raw_sz);
  for (const auto &offset : kDelphiNeighborOffsets) {
    for (std::size_t atom_idx = 0; atom_idx < surface_atoms.size(); ++atom_idx) {
      const Atom &atom = surface_atoms[atom_idx];
      if (atom.radius <= 0.0)
        continue;
      std::size_t cx = 0, cy = 0, cz = 0;
      raw_atom_central(atom, cx, cy, cz);
      const int tx = static_cast<int>(cx) + offset[0];
      const int ty = static_cast<int>(cy) + offset[1];
      const int tz = static_cast<int>(cz) + offset[2];
      if (tx < 0 || ty < 0 || tz < 0 ||
          tx > static_cast<int>(raw_vx) || ty > static_cast<int>(raw_vy) ||
          tz > static_cast<int>(raw_vz))
        continue;
      raw_voxel_atoms[raw_voxel_flat(static_cast<std::size_t>(tx),
                                     static_cast<std::size_t>(ty),
                                     static_cast<std::size_t>(tz))]
          .push_back(atom_idx);
    }
  }
  auto raw_grid_voxel = [&](double x, double y, double z, std::size_t &vx,
                            std::size_t &vy, std::size_t &vz) {
    auto coord = [&](double value, double origin,
                     std::size_t high) -> std::size_t {
      const double scaled = (value - origin) * raw_voxel_scale;
      if (scaled <= 0.0)
        return 0;
      return std::min<std::size_t>(high, static_cast<std::size_t>(scaled));
    };
    vx = coord(x, raw_voxel_origin_x, raw_vx);
    vy = coord(y, raw_voxel_origin_y, raw_vy);
    vz = coord(z, raw_voxel_origin_z, raw_vz);
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
              const bool valid_x = i + 1 < shape.nx;
              const bool valid_y = j + 1 < shape.ny;
              const bool valid_z = k + 1 < shape.nz;
              bool done_x = !valid_x;
              bool done_y = !valid_y;
              bool done_z = !valid_z;
              std::size_t gvx = 0, gvy = 0, gvz = 0;
              raw_grid_voxel(x, y, z, gvx, gvy, gvz);
              for (const auto &offset : kDelphiNeighborOffsets) {
                const int tx = static_cast<int>(gvx) + offset[0];
                const int ty = static_cast<int>(gvy) + offset[1];
                const int tz = static_cast<int>(gvz) + offset[2];
                if (tx < 0 || ty < 0 || tz < 0 ||
                    tx > static_cast<int>(raw_vx) ||
                    ty > static_cast<int>(raw_vy) ||
                    tz > static_cast<int>(raw_vz))
                  continue;
                const auto &bucket = raw_voxel_atoms[raw_voxel_flat(
                    static_cast<std::size_t>(tx),
                    static_cast<std::size_t>(ty),
                    static_cast<std::size_t>(tz))];
                for (std::size_t atom_idx : bucket) {
                  const Atom &atom = surface_atoms[atom_idx];
                  const double rdx = x - atom.x;
                  const double rdy = y - atom.y;
                  const double rdz = z - atom.z;
                  const double radius2 = atom.radius * atom.radius;
                  const std::int32_t owner =
                      static_cast<std::int32_t>(atom_idx + 2);
                  if (!done_x) {
                    const double dx = rdx + 0.5 * shape.h;
                    if (dx * dx + rdy * rdy + rdz * rdz < radius2) {
                      out.face_encoded[face_index(id, 0)] =
                          encode_delphi_face(owner, std::int32_t{1});
                      done_x = true;
                    }
                  }
                  if (!done_y) {
                    const double dy = rdy + 0.5 * shape.h;
                    if (rdx * rdx + dy * dy + rdz * rdz < radius2) {
                      out.face_encoded[face_index(id, 1)] =
                          encode_delphi_face(owner, std::int32_t{1});
                      done_y = true;
                    }
                  }
                  if (!done_z) {
                    const double dz = rdz + 0.5 * shape.h;
                    if (rdx * rdx + rdy * rdy + dz * dz < radius2) {
                      out.face_encoded[face_index(id, 2)] =
                          encode_delphi_face(owner, std::int32_t{1});
                      done_z = true;
                    }
                  }
                  if (done_x && done_y && done_z)
                    break;
                }
                if (done_x && done_y && done_z)
                  break;
              }
            }
          }
        }
      });

  for (std::size_t fid = 0; fid < out.face_media.size(); ++fid) {
    sync_face_from_encoded(fid);
    DelphiFaceDebug &debug = out.face_debug[fid];
    debug.raw_media = out.face_media[fid];
    debug.raw_owner = out.face_owner[fid];
    debug.raw_encoded = out.face_encoded[fid];
    debug.action = debug.raw_media != 0 ? kDelphiFaceRawVdw
                                        : kDelphiFaceUnchanged;
  }

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
  const AtomSpatialBins probe_bins =
      build_atom_bins(surface_atoms, args_.probe_radius,
                      std::max(shape.h, delphi_voxel_side));
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
  std::vector<unsigned char> vdwms_accessible_atoms(surface_atoms.size(), 1U);
  struct ClosestSurfaceAtom {
    std::ptrdiff_t index{-1};
    double surface_dist{std::numeric_limits<double>::infinity()};
  };
  auto closest_surface_atom = [&](double x, double y,
                                  double z) -> ClosestSurfaceAtom {
    ClosestSurfaceAtom best;
    if (probe_bins.empty)
      return best;
    const int radius_bins = static_cast<int>(
        std::ceil(std::max(probe_bins.max_radius, shape.h) /
                  std::max(1.0e-12, probe_bins.cell))) +
                            1;
    const int bx0 = static_cast<int>(
        std::floor((x - probe_bins.min_x) / probe_bins.cell));
    const int by0 = static_cast<int>(
        std::floor((y - probe_bins.min_y) / probe_bins.cell));
    const int bz0 = static_cast<int>(
        std::floor((z - probe_bins.min_z) / probe_bins.cell));
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
              probe_bins.bins[(static_cast<std::size_t>(bz) *
                                   probe_bins.by +
                               static_cast<std::size_t>(by)) *
                                  probe_bins.bx +
                              static_cast<std::size_t>(bx)];
          for (std::size_t atom_idx : bucket) {
            if (atom_idx < vdwms_accessible_atoms.size() &&
                vdwms_accessible_atoms[atom_idx] == 0U)
              continue;
            const Atom &atom = surface_atoms[atom_idx];
            const double ddx = x - atom.x;
            const double ddy = y - atom.y;
            const double ddz = z - atom.z;
            const double r2 = ddx * ddx + ddy * ddy + ddz * ddz;
            const double surface_distance = std::sqrt(r2) - atom.radius;
            if (surface_distance < best.surface_dist) {
              best.surface_dist = surface_distance;
              best.index = static_cast<std::ptrdiff_t>(atom_idx);
            }
          }
        }
      }
    }
    return best;
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
    constexpr int kInitialVertexCount = 12;
    constexpr int kRefinementLevel = 5;
    constexpr int kMaxVertices = 520;
    constexpr int kMaxEdges = 1040;
    constexpr double kSasShrink = 0.99999;
    auto probe_point_occupier = [&](double sx, double sy,
                                    double sz) -> std::int32_t {
      for (std::size_t atom_idx = 0; atom_idx < surface_atoms.size();
           ++atom_idx) {
        const Atom &atom = surface_atoms[atom_idx];
        const double radius = atom.radius + args_.probe_radius;
        const double ddx = sx - atom.x;
        const double ddy = sy - atom.y;
        const double ddz = sz - atom.z;
        const double r2 = ddx * ddx + ddy * ddy + ddz * ddz;
        if (r2 < radius * radius * kSasShrink)
          return static_cast<std::int32_t>(atom_idx + 1);
      }
      return 0;
    };
    struct SasEdge {
      int first{0};
      int second{0};
    };
    std::vector<std::array<double, 3>> vertex_array(kMaxVertices + 1);
    std::vector<SasEdge> edge_array(kMaxEdges + 1);
    std::vector<int> edges(kMaxEdges + 1, 0);
    const double angle_increment = 2.0 * kPi /
                                   static_cast<double>(kInitialVertexCount);
    for (int vertex_index = 1; vertex_index <= kInitialVertexCount;
         ++vertex_index) {
      const double angle =
          static_cast<double>(vertex_index - 1) * angle_increment;
      vertex_array[static_cast<std::size_t>(vertex_index)] = {
          std::cos(angle), std::sin(angle), 0.0};
      edge_array[static_cast<std::size_t>(vertex_index)] = {
          vertex_index,
          vertex_index != kInitialVertexCount ? vertex_index + 1 : 1};
    }
    int current_vertex_count = kInitialVertexCount;
    int current_edge_count = kInitialVertexCount;
    int start_edge_index = 1;
    int end_edge_index = 0;
    for (int level = 1; level <= kRefinementLevel; ++level) {
      start_edge_index = end_edge_index + 1;
      end_edge_index = current_edge_count;
      for (int edge_index = start_edge_index; edge_index <= end_edge_index;
           ++edge_index) {
        const int v1 = edge_array[static_cast<std::size_t>(edge_index)].first;
        const int v2 = edge_array[static_cast<std::size_t>(edge_index)].second;
        std::array<double, 3> midpoint{
            vertex_array[static_cast<std::size_t>(v1)][0] +
                vertex_array[static_cast<std::size_t>(v2)][0],
            vertex_array[static_cast<std::size_t>(v1)][1] +
                vertex_array[static_cast<std::size_t>(v2)][1],
            vertex_array[static_cast<std::size_t>(v1)][2] +
                vertex_array[static_cast<std::size_t>(v2)][2]};
        const double magnitude =
            std::sqrt(midpoint[0] * midpoint[0] +
                      midpoint[1] * midpoint[1] +
                      midpoint[2] * midpoint[2]);
        if (magnitude <= 1.0e-30)
          continue;
        ++current_vertex_count;
        vertex_array[static_cast<std::size_t>(current_vertex_count)] = {
            midpoint[0] / magnitude, midpoint[1] / magnitude,
            midpoint[2] / magnitude};
        ++current_edge_count;
        edges[static_cast<std::size_t>(edge_index)] = current_edge_count;
        edge_array[static_cast<std::size_t>(current_edge_count)] = {
            v1, current_vertex_count};
        ++current_edge_count;
        edge_array[static_cast<std::size_t>(current_edge_count)] = {
            current_vertex_count, v2};
      }
    }
    current_edge_count = end_edge_index;
    for (int edge_index = start_edge_index; edge_index <= current_edge_count;
         ++edge_index)
      edges[static_cast<std::size_t>(edge_index)] = -1;

    bins.points.reserve(surface_atoms.size() * 64);
    std::vector<unsigned char> has_pair(surface_atoms.size(), 0U);
    std::vector<unsigned char> accessible(surface_atoms.size(), 0U);
    for (std::size_t a = 0; a < surface_atoms.size(); ++a) {
      const Atom &atom_a = surface_atoms[a];
      const double ra = atom_a.radius + args_.probe_radius;
      if (atom_a.radius <= 0.0 || ra <= 0.0)
        continue;
      for (std::size_t b = a + 1; b < surface_atoms.size(); ++b) {
        const Atom &atom_b = surface_atoms[b];
        const double rb = atom_b.radius + args_.probe_radius;
        if (atom_b.radius <= 0.0 || rb <= 0.0)
          continue;
        const double dx = atom_b.x - atom_a.x;
        const double dy = atom_b.y - atom_a.y;
        const double dz = atom_b.z - atom_a.z;
        const double d2 = dx * dx + dy * dy + dz * dz;
        const double sum_r = ra + rb;
        const double diff_r = std::abs(ra - rb);
        if (sum_r * sum_r - d2 <= 0.01 || d2 <= diff_r * diff_r ||
            d2 <= 1.0e-16)
          continue;
        has_pair[a] = 1U;
        has_pair[b] = 1U;
        const double d = std::sqrt(d2);
        const double prefactor = 1.0 + (ra * ra - rb * rb) / d2;
        const double cx = atom_a.x + 0.5 * prefactor * dx;
        const double cy = atom_a.y + 0.5 * prefactor * dy;
        const double cz = atom_a.z + 0.5 * prefactor * dz;
        const double term1 = std::max(0.0, sum_r * sum_r - d2);
        const double term2 = std::max(0.0, d2 - diff_r * diff_r);
        const double circle_r = 0.5 * std::sqrt(term1) * std::sqrt(term2) / d;
        if (circle_r <= 0.0)
          continue;
        std::array<double, 3> circle_x{1.0, 0.0, 0.0};
        std::array<double, 3> circle_y{0.0, 1.0, 0.0};
        const double xy_projection = std::sqrt(dx * dx + dy * dy);
        if (xy_projection > 1.0e-8) {
          const double tangential_x = -dy / xy_projection;
          const double tangential_y = dx / xy_projection;
          const double cosine_theta = dz / d;
          const double sine_theta =
              std::sqrt(std::max(0.0, 1.0 - cosine_theta * cosine_theta));
          const double one_minus_cosine = 1.0 - cosine_theta;
          const double temp = one_minus_cosine * tangential_x;
          const double sine_tangential_x = sine_theta * tangential_x;
          const double sine_tangential_y = sine_theta * tangential_y;
          circle_x = {temp * tangential_x + cosine_theta,
                      temp * tangential_y, -sine_tangential_y};
          circle_y = {
              temp * tangential_y,
              one_minus_cosine * tangential_y * tangential_y + cosine_theta,
              sine_tangential_x};
        }
        auto circle_point = [&](const std::array<double, 3> &unit_point) {
          return std::array<double, 3>{
              cx + circle_r *
                       (unit_point[0] * circle_x[0] +
                        unit_point[1] * circle_y[0] + unit_point[2] * 0.0),
              cy + circle_r *
                       (unit_point[0] * circle_x[1] +
                        unit_point[1] * circle_y[1] + unit_point[2] * 0.0),
              cz + circle_r *
                       (unit_point[0] * circle_x[2] +
                        unit_point[1] * circle_y[2] + unit_point[2] * 1.0)};
        };
        std::vector<std::int32_t> vertex_occupation(kMaxVertices + 1, 0);
        std::vector<int> edge_state(kMaxEdges + 1, 0);
        int num_visible_vertices = 0;
        for (int vertex_index = 1; vertex_index <= kInitialVertexCount;
             ++vertex_index) {
          const auto point = circle_point(
              vertex_array[static_cast<std::size_t>(vertex_index)]);
          const std::int32_t occupation =
              probe_point_occupier(point[0], point[1], point[2]);
          vertex_occupation[static_cast<std::size_t>(vertex_index)] =
              occupation;
          if (occupation == 0) {
            ++num_visible_vertices;
            bins.points.push_back(point);
            vertex_occupation[static_cast<std::size_t>(vertex_index)] = 0;
          }
        }

        int edge_state_count = 0;
        for (int edge_index = kInitialVertexCount; edge_index >= 1;
             --edge_index) {
          const int v1 = edge_array[static_cast<std::size_t>(edge_index)].first;
          const int v2 = edge_array[static_cast<std::size_t>(edge_index)].second;
          const std::int32_t occ1 =
              vertex_occupation[static_cast<std::size_t>(v1)];
          const std::int32_t occ2 =
              vertex_occupation[static_cast<std::size_t>(v2)];
          if (occ1 > 0 && occ1 == occ2)
            continue;
          ++edge_state_count;
          edge_state[static_cast<std::size_t>(edge_state_count)] = edge_index;
        }

        while (edge_state_count > 0) {
          const int edge_index =
              edge_state[static_cast<std::size_t>(edge_state_count)];
          --edge_state_count;
          const int endpoint1 =
              edge_array[static_cast<std::size_t>(edge_index)].first;
          const int endpoint2 =
              edge_array[static_cast<std::size_t>(edge_index)].second;
          const std::int32_t vertex_index1 =
              vertex_occupation[static_cast<std::size_t>(endpoint1)];
          const std::int32_t vertex_index2 =
              vertex_occupation[static_cast<std::size_t>(endpoint2)];
          const int child_edge =
              edges[static_cast<std::size_t>(edge_index)];
          if (vertex_index1 > static_cast<std::int32_t>(surface_atoms.size()) ||
              vertex_index2 > static_cast<std::int32_t>(surface_atoms.size()))
            continue;
          const int refined_vertex_index = edge_index + kInitialVertexCount;
          if (refined_vertex_index < 1 ||
              refined_vertex_index > current_vertex_count)
            continue;
          const auto point = circle_point(
              vertex_array[static_cast<std::size_t>(refined_vertex_index)]);
          const std::int32_t occ_mid =
              probe_point_occupier(point[0], point[1], point[2]);
          bool handled = false;
          if (vertex_index1 != 0 &&
              occ_mid == vertex_index1) {
            vertex_occupation[static_cast<std::size_t>(refined_vertex_index)] =
                vertex_index1;
            if (child_edge > 0) {
              ++edge_state_count;
              edge_state[static_cast<std::size_t>(edge_state_count)] =
                  child_edge + 1;
            }
            handled = true;
          }
          if (!handled && vertex_index2 != 0 &&
              occ_mid == vertex_index2) {
            vertex_occupation[static_cast<std::size_t>(refined_vertex_index)] =
                vertex_index2;
            if (child_edge > 0) {
              ++edge_state_count;
              edge_state[static_cast<std::size_t>(edge_state_count)] =
                  child_edge;
            }
            handled = true;
          }
          if (!handled && occ_mid != 0) {
            vertex_occupation[static_cast<std::size_t>(refined_vertex_index)] =
                occ_mid;
            if (child_edge > 0) {
              ++edge_state_count;
              edge_state[static_cast<std::size_t>(edge_state_count)] =
                  child_edge + 1;
              ++edge_state_count;
              edge_state[static_cast<std::size_t>(edge_state_count)] =
                  child_edge;
            }
            handled = true;
          }
          if (handled)
            continue;

          ++num_visible_vertices;
          vertex_occupation[static_cast<std::size_t>(refined_vertex_index)] =
              0;
          bins.points.push_back(point);
          if (child_edge > 0) {
            const int child_a = child_edge;
            const int child_b = child_edge + 1;
            if ((edges[static_cast<std::size_t>(child_b)] > 0) ||
                vertex_index2 > 0) {
              ++edge_state_count;
              edge_state[static_cast<std::size_t>(edge_state_count)] =
                  child_b;
            }
            if ((edges[static_cast<std::size_t>(child_a)] > 0) ||
                vertex_index1 > 0) {
              ++edge_state_count;
              edge_state[static_cast<std::size_t>(edge_state_count)] =
                  child_a;
            }
          }
        }
        if (num_visible_vertices > 0) {
          accessible[a] = 1U;
          accessible[b] = 1U;
        }
      }
    }
    for (std::size_t atom_idx = 0; atom_idx < surface_atoms.size(); ++atom_idx) {
      if (has_pair[atom_idx] == 0U)
        accessible[atom_idx] = 1U;
    }
    vdwms_accessible_atoms = std::move(accessible);
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
  auto nearest_exposed_distance = [&](double x, double y, double z) {
    if (exposed_bins.empty)
      return std::numeric_limits<double>::infinity();
    const int radius_bins = static_cast<int>(
        std::ceil(args_.probe_radius / std::max(1.0e-12, exposed_bins.cell))) +
                            1;
    const int bx0 = static_cast<int>(
        std::floor((x - exposed_bins.min_x) / exposed_bins.cell));
    const int by0 = static_cast<int>(
        std::floor((y - exposed_bins.min_y) / exposed_bins.cell));
    const int bz0 = static_cast<int>(
        std::floor((z - exposed_bins.min_z) / exposed_bins.cell));
    double best = std::numeric_limits<double>::infinity();
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
            if (r2 < best)
              best = r2;
          }
        }
      }
    }
    return std::sqrt(best);
  };
  auto nearest_exposed_point = [&](double x, double y, double z,
                                   std::array<double, 3> &best_point) {
    if (exposed_bins.empty)
      return std::numeric_limits<double>::infinity();
    const int radius_bins = static_cast<int>(
        std::ceil(args_.probe_radius / std::max(1.0e-12, exposed_bins.cell))) +
                            1;
    const int bx0 = static_cast<int>(
        std::floor((x - exposed_bins.min_x) / exposed_bins.cell));
    const int by0 = static_cast<int>(
        std::floor((y - exposed_bins.min_y) / exposed_bins.cell));
    const int bz0 = static_cast<int>(
        std::floor((z - exposed_bins.min_z) / exposed_bins.cell));
    double best = std::numeric_limits<double>::infinity();
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
            if (r2 < best) {
              best = r2;
              best_point = point;
            }
          }
        }
      }
    }
    return std::sqrt(best);
  };

  struct DelphiFaceRef {
    std::size_t id{0};
    int axis{0};
    bool valid{false};
  };
  auto face_ref_for_neighbor = [&](std::size_t i, std::size_t j,
                                   std::size_t k,
                                   int neighbor) -> DelphiFaceRef {
    const std::size_t id = idx3(i, j, k, shape.nx, shape.ny);
    switch (neighbor) {
    case 1:
      return {id, 0, i + 1 < shape.nx};
    case 2:
      return {id, 1, j + 1 < shape.ny};
    case 3:
      return {id, 2, k + 1 < shape.nz};
    case 4:
      return {i > 0 ? id - 1 : id, 0, i > 0};
    case 5:
      return {j > 0 ? id - shape.nx : id, 1, j > 0};
    case 6:
      return {k > 0 ? id - plane : id, 2, k > 0};
    default:
      return {};
    }
  };
  auto face_midpoint = [&](std::size_t i, std::size_t j, std::size_t k,
                           int neighbor) {
    std::array<double, 3> point{
        shape.ox + static_cast<double>(i) * shape.h,
        shape.oy + static_cast<double>(j) * shape.h,
        shape.oz + static_cast<double>(k) * shape.h};
    if (neighbor == 1)
      point[0] += 0.5 * shape.h;
    else if (neighbor == 2)
      point[1] += 0.5 * shape.h;
    else if (neighbor == 3)
      point[2] += 0.5 * shape.h;
    else if (neighbor == 4)
      point[0] -= 0.5 * shape.h;
    else if (neighbor == 5)
      point[1] -= 0.5 * shape.h;
    else if (neighbor == 6)
      point[2] -= 0.5 * shape.h;
    return point;
  };
  auto face_flat = [&](const DelphiFaceRef &ref) {
    return face_index(ref.id, ref.axis);
  };
  auto decode_index = [&](std::size_t id, std::size_t &i, std::size_t &j,
                          std::size_t &k) {
    k = id / plane;
    const std::size_t rem = id - k * plane;
    j = rem / shape.nx;
    i = rem - j * shape.nx;
  };
  struct BoundaryState {
    bool boundary{false};
    bool external{false};
  };
  auto classify_boundary = [&](std::size_t i, std::size_t j,
                               std::size_t k) -> BoundaryState {
    BoundaryState state;
    if (i == 0 || j == 0 || k == 0 || i + 1 >= shape.nx ||
        j + 1 >= shape.ny || k + 1 >= shape.nz)
      return state;
    std::array<std::int32_t, 6> faces{};
    face_media6(i, j, k, faces);
    state.external = std::any_of(faces.begin(), faces.end(),
                                 [](std::int32_t v) { return v == 0; });
    state.boundary = faces[0] != faces[5];
    for (std::size_t n = 1; n < faces.size(); ++n)
      state.boundary = state.boundary || faces[n] != faces[n - 1];
    return state;
  };

  std::vector<std::uint32_t> boundary_rank(shape.n, 0U);
  std::vector<unsigned char> boundary_external(shape.n, 0U);
  std::vector<std::size_t> boundary_queue;
  boundary_queue.reserve(boundary_voxels.size() * 2);
  std::size_t discovered_boundary = 0;
  std::size_t external_boundary = 0;
  if (shape.nx >= 3 && shape.ny >= 3 && shape.nz >= 3) {
    for (std::size_t k = 1; k + 1 < shape.nz; ++k) {
      for (std::size_t j = 1; j + 1 < shape.ny; ++j) {
        for (std::size_t i = 1; i + 1 < shape.nx; ++i) {
          const BoundaryState state = classify_boundary(i, j, k);
          if (!state.boundary)
            continue;
          const std::size_t id = idx3(i, j, k, shape.nx, shape.ny);
          boundary_rank[id] = static_cast<std::uint32_t>(++discovered_boundary);
          boundary_external[id] = state.external ? 1U : 0U;
          if (state.external) {
            ++external_boundary;
            boundary_queue.push_back(id);
          }
        }
      }
    }
  }
  out.audit.boundary_initial = discovered_boundary;

  std::size_t added = 0;
  std::size_t removed = 0;
  std::size_t exposed = 0;
  std::size_t contact = 0;
  auto update_boundary_status = [&](std::size_t id, bool enqueue_new) {
    std::size_t i = 0, j = 0, k = 0;
    decode_index(id, i, j, k);
    const BoundaryState state = classify_boundary(i, j, k);
    if (!state.boundary) {
      if (boundary_rank[id] != 0U) {
        if (boundary_external[id] != 0U && external_boundary > 0)
          --external_boundary;
        boundary_rank[id] = 0U;
        boundary_external[id] = 0U;
        ++removed;
      }
      return;
    }
    if (boundary_rank[id] == 0U) {
      boundary_rank[id] = static_cast<std::uint32_t>(++discovered_boundary);
      boundary_external[id] = state.external ? 1U : 0U;
      if (state.external)
        ++external_boundary;
      if (enqueue_new)
        boundary_queue.push_back(id);
      return;
    }
    if (boundary_external[id] != (state.external ? 1U : 0U)) {
      if (state.external)
        ++external_boundary;
      else if (external_boundary > 0)
        --external_boundary;
      boundary_external[id] = state.external ? 1U : 0U;
    }
  };
  auto update_adjacent_boundary_status = [&](std::size_t i, std::size_t j,
                                             std::size_t k, int neighbor) {
    const std::size_t center_id = idx3(i, j, k, shape.nx, shape.ny);
    update_boundary_status(center_id, false);
    std::size_t ni = i, nj = j, nk = k;
    if (neighbor == 1)
      ++ni;
    else if (neighbor == 2)
      ++nj;
    else if (neighbor == 3)
      ++nk;
    else if (neighbor == 4)
      --ni;
    else if (neighbor == 5)
      --nj;
    else if (neighbor == 6)
      --nk;
    if (ni > 0 && nj > 0 && nk > 0 && ni + 1 < shape.nx &&
        nj + 1 < shape.ny && nk + 1 < shape.nz)
      update_boundary_status(idx3(ni, nj, nk, shape.nx, shape.ny), true);
  };
  auto set_contact_solvent = [&](const DelphiFaceRef &ref,
                                 std::int32_t entity,
                                 std::uint16_t action,
                                 std::uint32_t iteration,
                                 const ClosestSurfaceAtom &nearest_atom,
                                 double nearest_exposed_dist) {
    (void)entity;
    const std::size_t fid = face_flat(ref);
    const std::int32_t original = out.face_encoded[fid];
    std::int32_t original_entity = decode_entity_from_value(original);
    const std::int32_t original_media = decode_media_from_value(original);
    if (original_entity >= 0) {
      if (original_entity == 0)
        original_entity = 1;
      set_face_encoded(fid, encode_delphi_face(-original_entity,
                                               original_media));
    }
    if (fid < out.face_debug.size()) {
      DelphiFaceDebug &debug = out.face_debug[fid];
      debug.action = action;
      debug.iteration = iteration;
      debug.nearest_exposed_dist = nearest_exposed_dist;
      if (nearest_atom.index >= 0) {
        debug.nearest_atom =
            static_cast<std::int32_t>(nearest_atom.index + 1);
        debug.nearest_atom_surface_dist = nearest_atom.surface_dist;
      }
    }
  };
  auto set_cavity_media = [&](const DelphiFaceRef &ref,
                              std::uint32_t iteration,
                              const ClosestSurfaceAtom &nearest_atom,
                              double nearest_exposed_dist) {
    const std::size_t fid = face_flat(ref);
    if (decode_media_from_value(out.face_encoded[fid]) == 0)
      ++added;
    const std::int32_t media =
        nearest_atom.index >= 0 ? std::int32_t{1} : std::int32_t{1};
    set_face_encoded(fid, encode_delphi_face(std::int32_t{1}, media));
    if (fid < out.face_debug.size()) {
      DelphiFaceDebug &debug = out.face_debug[fid];
      debug.action = kDelphiFaceDelphiCavityRemap;
      debug.iteration = iteration;
      debug.nearest_exposed_dist = nearest_exposed_dist;
      if (nearest_atom.index >= 0) {
        debug.nearest_atom =
            static_cast<std::int32_t>(nearest_atom.index + 1);
        debug.nearest_atom_surface_dist = nearest_atom.surface_dist;
      }
    }
  };

  const std::size_t max_boundary_iterations =
      std::max<std::size_t>(32, shape.nx + shape.ny + shape.nz);
  std::size_t batch_begin = 0;
  std::size_t iterations = 0;
  while (batch_begin < boundary_queue.size() &&
         iterations < max_boundary_iterations) {
    const std::size_t batch_end = boundary_queue.size();
    ++iterations;
    for (std::size_t q = batch_begin; q < batch_end; ++q) {
      const std::size_t center_id = boundary_queue[q];
      if (center_id >= shape.n || boundary_rank[center_id] == 0U)
        continue;
      std::size_t i = 0, j = 0, k = 0;
      decode_index(center_id, i, j, k);
      std::array<std::int32_t, 7> neighbor_entities{};
      std::array<std::int32_t, 7> neighbor_media{};
      bool rm_boundary_pt = false;
      for (int neighbor = 1; neighbor <= 6; ++neighbor) {
        const DelphiFaceRef ref = face_ref_for_neighbor(i, j, k, neighbor);
        if (!ref.valid)
          continue;
        const std::size_t fid = face_flat(ref);
        neighbor_entities[static_cast<std::size_t>(neighbor)] =
            decode_entity_from_value(out.face_encoded[fid]);
        neighbor_media[static_cast<std::size_t>(neighbor)] =
            decode_media_from_value(out.face_encoded[fid]);
        const std::int32_t entity =
            neighbor_entities[static_cast<std::size_t>(neighbor)];
        if (entity > 1 &&
            entity <= static_cast<std::int32_t>(surface_atoms.size() + 1))
          rm_boundary_pt = true;
      }
      for (int neighbor = 1; neighbor <= 6; ++neighbor) {
        const DelphiFaceRef ref = face_ref_for_neighbor(i, j, k, neighbor);
        if (!ref.valid)
          continue;
        const std::int32_t entity =
            neighbor_entities[static_cast<std::size_t>(neighbor)];
        const std::int32_t media =
            neighbor_media[static_cast<std::size_t>(neighbor)];
        const bool should_process =
            entity == 0 ||
            (rm_boundary_pt &&
             entity > static_cast<std::int32_t>(surface_atoms.size() + 1)) ||
            (media == 0 && entity > 0);
        if (!should_process)
          continue;
        const auto midpoint = face_midpoint(i, j, k, neighbor);
        const double exposed_dist =
            nearest_exposed_distance(midpoint[0], midpoint[1], midpoint[2]);
        if (exposed_dist < args_.probe_radius) {
          set_contact_solvent(ref, 1, kDelphiFaceExposedSolvent,
                              static_cast<std::uint32_t>(iterations),
                              ClosestSurfaceAtom{}, exposed_dist);
          ++exposed;
          continue;
        }
        const ClosestSurfaceAtom contact_atom =
            closest_surface_atom(midpoint[0], midpoint[1], midpoint[2]);
        if (contact_atom.index >= 0 &&
            !probe_surface_occluded(midpoint[0], midpoint[1], midpoint[2],
                                    static_cast<std::size_t>(
                                        contact_atom.index))) {
          set_contact_solvent(ref,
                              static_cast<std::int32_t>(contact_atom.index + 2),
                              kDelphiFaceContactSolvent,
                              static_cast<std::uint32_t>(iterations),
                              contact_atom, exposed_dist);
          ++contact;
        } else {
          set_cavity_media(ref, static_cast<std::uint32_t>(iterations),
                           contact_atom, exposed_dist);
          update_adjacent_boundary_status(i, j, k, neighbor);
        }
      }
      update_boundary_status(center_id, false);
    }
    batch_begin = batch_end;
  }

  for (std::size_t fid = 0; fid < out.face_encoded.size(); ++fid) {
    const std::int32_t encoded = out.face_encoded[fid];
    if (encoded < 0) {
      const std::int32_t signed_entity = encoded % delphi_epsilon_dim;
      std::int32_t final_value = -encoded;
      if (signed_entity == -1)
        final_value -= 1;
      set_face_encoded(fid, final_value);
      if (fid < out.face_debug.size() &&
          out.face_debug[fid].action == kDelphiFaceUnchanged)
        out.face_debug[fid].action = kDelphiFaceSignPreserve;
    } else {
      sync_face_from_encoded(fid);
    }
  }

  out.audit.added = added;
  out.audit.removed = removed;
  out.audit.exposed_points = exposed_bins.points.empty()
                                  ? exposed + contact
                                  : exposed_bins.points.size();
  std::size_t final_boundary = 0;
  for (std::uint32_t rank : boundary_rank) {
    if (rank != 0U)
      ++final_boundary;
  }
  out.audit.boundary_final = final_boundary;
  const double face_count =
      static_cast<double>(std::max<std::size_t>(1, shape.n * 3));
  out.audit.midpoint_rms =
      (args_.ext_diel - args_.int_diel) *
      std::sqrt(static_cast<double>(added + removed) / face_count);
  out.iscm_surface_valid.assign(shape.n, 0U);
  out.iscm_surface_positions.assign(shape.n, {0.0, 0.0, 0.0});
  for (std::size_t id = 0; id < shape.n; ++id) {
    if (boundary_rank[id] == 0U)
      continue;
    std::size_t i = 0, j = 0, k = 0;
    decode_index(id, i, j, k);
    const std::array<double, 3> grid_point{
        shape.ox + static_cast<double>(i) * shape.h,
        shape.oy + static_cast<double>(j) * shape.h,
        shape.oz + static_cast<double>(k) * shape.h};
    std::array<double, 3> surface_point = grid_point;
    bool scaled = false;
    if (id < boundary_external.size() && boundary_external[id] == 0U &&
        surface_atoms.size() >= 2) {
      std::size_t closest1 = std::numeric_limits<std::size_t>::max();
      std::size_t closest2 = std::numeric_limits<std::size_t>::max();
      double min1 = std::numeric_limits<double>::max();
      double min2 = std::numeric_limits<double>::max();
      for (std::size_t atom_idx = 0; atom_idx < surface_atoms.size();
           ++atom_idx) {
        const Atom &atom = surface_atoms[atom_idx];
        const double dx = grid_point[0] - atom.x;
        const double dy = grid_point[1] - atom.y;
        const double dz = grid_point[2] - atom.z;
        const double distance_sq =
            dx * dx + dy * dy + dz * dz - atom.radius * atom.radius;
        if (distance_sq < min1) {
          closest2 = closest1;
          min2 = min1;
          closest1 = atom_idx;
          min1 = distance_sq;
        } else if (distance_sq <= min2) {
          closest2 = atom_idx;
          min2 = distance_sq;
        }
      }
      if (closest1 != std::numeric_limits<std::size_t>::max() &&
          closest2 != std::numeric_limits<std::size_t>::max() &&
          closest1 != closest2) {
        const Atom &atom1 = surface_atoms[closest1];
        const Atom &atom2 = surface_atoms[closest2];
        const double dx12 = atom2.x - atom1.x;
        const double dy12 = atom2.y - atom1.y;
        const double dz12 = atom2.z - atom1.z;
        const double dist12_sq =
            dx12 * dx12 + dy12 * dy12 + dz12 * dz12;
        if (dist12_sq > 1.0e-30) {
          const double factor = 0.5 * (min2 - min1) / dist12_sq;
          surface_point = {grid_point[0] + factor * dx12,
                           grid_point[1] + factor * dy12,
                           grid_point[2] + factor * dz12};
          scaled = true;
        }
      }
    }
    if (scaled) {
      out.iscm_surface_valid[id] = 1U;
      out.iscm_surface_positions[id] = surface_point;
      continue;
    }
    const ClosestSurfaceAtom nearest_atom =
        closest_surface_atom(grid_point[0], grid_point[1], grid_point[2]);
    if (nearest_atom.index >= 0 &&
        static_cast<std::size_t>(nearest_atom.index) < surface_atoms.size()) {
      const Atom &atom = surface_atoms[static_cast<std::size_t>(
          nearest_atom.index)];
      const double dx = grid_point[0] - atom.x;
      const double dy = grid_point[1] - atom.y;
      const double dz = grid_point[2] - atom.z;
      const double distance = std::sqrt(dx * dx + dy * dy + dz * dz);
      bool solvent_exposed = true;
      if (args_.probe_radius > 0.0) {
        solvent_exposed =
            !probe_surface_occluded(grid_point[0], grid_point[1],
                                    grid_point[2],
                                    static_cast<std::size_t>(
                                        nearest_atom.index));
      }
      if (distance > 1.0e-30 && solvent_exposed) {
        surface_point = {atom.x + dx * (atom.radius / distance),
                         atom.y + dy * (atom.radius / distance),
                         atom.z + dz * (atom.radius / distance)};
        scaled = true;
      }
    }
    if (!scaled && args_.probe_radius > 0.0 && !exposed_bins.empty) {
      std::array<double, 3> exposed_point{0.0, 0.0, 0.0};
      const double exposed_distance =
          nearest_exposed_point(grid_point[0], grid_point[1], grid_point[2],
                                exposed_point);
      if (std::isfinite(exposed_distance) && exposed_distance > 1.0e-30) {
        const double factor = args_.probe_radius / exposed_distance;
        surface_point = {
            exposed_point[0] + (grid_point[0] - exposed_point[0]) * factor,
            exposed_point[1] + (grid_point[1] - exposed_point[1]) * factor,
            exposed_point[2] + (grid_point[2] - exposed_point[2]) * factor};
      }
    }
    out.iscm_surface_valid[id] = 1U;
    out.iscm_surface_positions[id] = surface_point;
  }
  (void)positive_face;
  return out;
}

void AnalyzerImpl::apply_delphi_vdw_media_map(
    const GridShape &shape, const DelphiVdwMediaMap &media_map, double eps_in,
    double eps_out) {
  if (media_map.face_media.size() != shape.n * 3 ||
      media_map.face_encoded.size() != shape.n * 3)
    throw std::runtime_error("PBE: internal DelPhi VDW media map has wrong size.");
  for (std::size_t id = 0; id < shape.n; ++id) {
    eps_x_[id] = media_map.face_media[id * 3] != 0 ? eps_in : eps_out;
    eps_y_[id] = media_map.face_media[id * 3 + 1] != 0 ? eps_in : eps_out;
    eps_z_[id] = media_map.face_media[id * 3 + 2] != 0 ? eps_in : eps_out;
  }
  maybe_write_midpoint_face_diff(shape, media_map, eps_in, eps_out);
  last_vdwms_audit_ = media_map.audit;
  last_iscm_surface_valid_ = media_map.iscm_surface_valid;
  last_iscm_surface_positions_ = media_map.iscm_surface_positions;
}

void AnalyzerImpl::maybe_write_midpoint_face_diff(
    const GridShape &shape, const DelphiVdwMediaMap &media_map, double eps_in,
    double eps_out) const {
  if (args_.test_midpoint_diff == "off")
    return;
  if (args_.test_pqr.empty())
    throw std::runtime_error(
        "PBE: -pbe-test-midpoint-diff requires -pbe-test-pqr.");
  if (!oracle_replay_applies(shape))
    return;
  if (media_map.face_media.size() != shape.n * 3)
    throw std::runtime_error(
        "PBE: internal DelPhi VDW media map has wrong size for midpoint diff.");
  if (std::abs(eps_in - args_.int_diel) > 1.0e-10 ||
      std::abs(eps_out - args_.ext_diel) > 1.0e-10) {
    return;
  }

  std::ifstream in("pbe_oracle_midpoint_epsilon.csv");
  if (!in.is_open()) {
    throw std::runtime_error(
        "PBE: -pbe-test-midpoint-diff requires "
        "pbe_oracle_midpoint_epsilon.csv in the run directory.");
  }
  std::string line;
  if (!std::getline(in, line))
    throw std::runtime_error("PBE: empty pbe_oracle_midpoint_epsilon.csv.");
  const auto header = csv_header_map(split_csv_simple(line));
  const auto fixture_it = header.find("fixture");
  auto fixture_matches = [&](const std::vector<std::string> &fields) {
    if (fixture_it == header.end() || fixture_it->second >= fields.size())
      return true;
    const std::string fixture = fields[fixture_it->second];
    if (fixture.empty() || current_source_label_.empty())
      return true;
    if (fixture == current_source_label_)
      return true;
    if (current_source_label_.rfind(fixture + "__", 0) == 0 ||
        fixture.rfind(current_source_label_ + "__", 0) == 0)
      return true;
    // TITAN trajectory labels are often input filenames, while pyDelPhi oracle
    // rows use fixture names. The oracle file is copied per run, so keep the
    // rows comparable instead of dropping the entire file on label mismatch.
    return true;
  };

  auto face_index = [](std::size_t id, int axis) {
    return id * 3 + static_cast<std::size_t>(axis);
  };
  auto axis_name = [](int axis) -> const char * {
    if (axis == 0)
      return "x";
    if (axis == 1)
      return "y";
    return "z";
  };

  if (args_.test_midpoint_diff == "full") {
    std::ofstream state_out;
    if (!native_midpoint_state_header_written_) {
      state_out.open("pbe_native_midpoint_state.csv",
                     std::ios::out | std::ios::trunc);
      if (!state_out.is_open())
        throw std::runtime_error(
            "PBE: failed to create pbe_native_midpoint_state.csv.");
      state_out << "fixture,i,j,k,raw_idx_x,raw_idx_y,raw_idx_z,idx_x,idx_y,"
                   "idx_z,raw_media_x,raw_media_y,raw_media_z,media_x,"
                   "media_y,media_z,eps_x,eps_y,eps_z\n";
      native_midpoint_state_header_written_ = true;
    } else {
      state_out.open("pbe_native_midpoint_state.csv",
                     std::ios::out | std::ios::app);
      if (!state_out.is_open())
        throw std::runtime_error(
            "PBE: failed to append pbe_native_midpoint_state.csv.");
    }
    state_out << std::setprecision(17);
    const std::string fixture =
        current_source_label_.empty() ? "unknown" : current_source_label_;
    for (std::size_t i = 0; i < shape.nx; ++i) {
      for (std::size_t j = 0; j < shape.ny; ++j) {
        for (std::size_t k = 0; k < shape.nz; ++k) {
          const std::size_t id = idx3(i, j, k, shape.nx, shape.ny);
          std::int32_t raw_idx[3]{0, 0, 0};
          std::int32_t final_idx[3]{0, 0, 0};
          std::int32_t raw_media[3]{0, 0, 0};
          std::int32_t final_media[3]{0, 0, 0};
          double eps[3]{eps_out, eps_out, eps_out};
          for (int axis = 0; axis < 3; ++axis) {
            const std::size_t fid = face_index(id, axis);
            if (fid < media_map.face_debug.size()) {
              raw_idx[axis] = media_map.face_debug[fid].raw_encoded;
              raw_media[axis] = media_map.face_debug[fid].raw_media;
            }
            if (fid < media_map.face_encoded.size())
              final_idx[axis] = media_map.face_encoded[fid];
            if (fid < media_map.face_media.size())
              final_media[axis] = media_map.face_media[fid];
            eps[axis] = final_media[axis] != 0 ? eps_in : eps_out;
          }
          state_out << csv_escape(fixture) << ',' << i << ',' << j << ','
                    << k << ',' << raw_idx[0] << ',' << raw_idx[1] << ','
                    << raw_idx[2] << ',' << final_idx[0] << ','
                    << final_idx[1] << ',' << final_idx[2] << ','
                    << raw_media[0] << ',' << raw_media[1] << ','
                    << raw_media[2] << ',' << final_media[0] << ','
                    << final_media[1] << ',' << final_media[2] << ','
                    << eps[0] << ',' << eps[1] << ',' << eps[2] << '\n';
        }
      }
    }

    std::map<std::tuple<std::uint16_t, std::int32_t, std::int32_t>,
             std::size_t>
        histogram;
    for (std::size_t fid = 0; fid < media_map.face_debug.size(); ++fid) {
      const auto &debug = media_map.face_debug[fid];
      const std::int32_t final_media =
          fid < media_map.face_media.size() ? media_map.face_media[fid] : 0;
      ++histogram[std::make_tuple(debug.action, debug.raw_media,
                                  final_media)];
    }
    std::ofstream hist_out;
    if (!vdwms_action_histogram_header_written_) {
      hist_out.open("pbe_vdwms_action_histogram.csv",
                    std::ios::out | std::ios::trunc);
      if (!hist_out.is_open())
        throw std::runtime_error(
            "PBE: failed to create pbe_vdwms_action_histogram.csv.");
      hist_out << "fixture,action,raw_media,final_media,count\n";
      vdwms_action_histogram_header_written_ = true;
    } else {
      hist_out.open("pbe_vdwms_action_histogram.csv",
                    std::ios::out | std::ios::app);
      if (!hist_out.is_open())
        throw std::runtime_error(
            "PBE: failed to append pbe_vdwms_action_histogram.csv.");
    }
    for (const auto &[key, count] : histogram) {
      const auto [action, raw_media, final_media] = key;
      hist_out << csv_escape(fixture) << ','
               << delphi_face_action_name(action) << ',' << raw_media << ','
               << final_media << ',' << count << '\n';
    }
  }

  std::ofstream out;
  if (!midpoint_face_diff_header_written_) {
    out.open("pbe_midpoint_face_diff.csv", std::ios::out | std::ios::trunc);
    if (!out.is_open())
      throw std::runtime_error(
          "PBE: failed to create pbe_midpoint_face_diff.csv.");
    out << "i,j,k,axis,native_eps,oracle_eps,native_media,oracle_media,"
           "native_owner,native_encoded,oracle_encoded,face_x_A,face_y_A,"
           "face_z_A,queue_iteration,"
           "remap_action,nearest_atom,nearest_atom_surface_dist,"
           "nearest_exposed_dist,raw_media,raw_owner,raw_encoded,final_encoded\n";
    midpoint_face_diff_header_written_ = true;
  } else {
    out.open("pbe_midpoint_face_diff.csv", std::ios::out | std::ios::app);
    if (!out.is_open())
      throw std::runtime_error(
          "PBE: failed to append pbe_midpoint_face_diff.csv.");
  }
  out << std::setprecision(17);

  auto face_coord = [&](std::size_t i, std::size_t j, std::size_t k,
                        int axis) {
    std::array<double, 3> point{
        shape.ox + static_cast<double>(i) * shape.h,
        shape.oy + static_cast<double>(j) * shape.h,
        shape.oz + static_cast<double>(k) * shape.h};
    if (axis == 0)
      point[0] += 0.5 * shape.h;
    else if (axis == 1)
      point[1] += 0.5 * shape.h;
    else
      point[2] += 0.5 * shape.h;
    return point;
  };
  auto oracle_media = [&](double eps_value) {
    if (std::abs(eps_value - eps_in) <= 1.0e-9)
      return 1;
    if (std::abs(eps_value - eps_out) <= 1.0e-9)
      return 0;
    return -1;
  };
  auto write_optional_double = [&](double value) {
    if (std::isfinite(value) && std::abs(value) < 1.0e299)
      out << value;
    else
      out << "NA";
  };

  const bool full = args_.test_midpoint_diff == "full";
  const std::size_t sample_limit = full ? std::numeric_limits<std::size_t>::max()
                                        : std::size_t{50};
  std::size_t compared_rows = 0;
  std::size_t written_rows = 0;
  std::size_t skipped_fixture_rows = 0;
  while (std::getline(in, line)) {
    if (trim_copy(line).empty())
      continue;
    const auto fields = split_csv_simple(line);
    if (!fixture_matches(fields)) {
      ++skipped_fixture_rows;
      continue;
    }
    std::size_t i = 0, j = 0, k = 0;
    double oracle_eps[3]{0.0, 0.0, 0.0};
    double oracle_encoded[3]{std::numeric_limits<double>::quiet_NaN(),
                             std::numeric_limits<double>::quiet_NaN(),
                             std::numeric_limits<double>::quiet_NaN()};
    if (!csv_index(fields, header, "i", i) ||
        !csv_index(fields, header, "j", j) ||
        !csv_index(fields, header, "k", k)) {
      throw std::runtime_error(
          "PBE: pbe_oracle_midpoint_epsilon.csv row is missing i,j,k.");
    }
    if (!csv_double(fields, header, {"eps_x", "epsilon_x", "ex"},
                    oracle_eps[0]) ||
        !csv_double(fields, header, {"eps_y", "epsilon_y", "ey"},
                    oracle_eps[1]) ||
        !csv_double(fields, header, {"eps_z", "epsilon_z", "ez"},
                    oracle_eps[2])) {
      throw std::runtime_error(
          "PBE: pbe_oracle_midpoint_epsilon.csv row is missing eps_x/eps_y/eps_z.");
    }
    if (i >= shape.nx || j >= shape.ny || k >= shape.nz)
      throw std::runtime_error(
          "PBE: pbe_oracle_midpoint_epsilon.csv index is outside the grid.");
    (void)csv_double(fields, header, {"idx_x", "index_x", "encoded_x"},
                     oracle_encoded[0]);
    (void)csv_double(fields, header, {"idx_y", "index_y", "encoded_y"},
                     oracle_encoded[1]);
    (void)csv_double(fields, header, {"idx_z", "index_z", "encoded_z"},
                     oracle_encoded[2]);
    const std::size_t id = idx3(i, j, k, shape.nx, shape.ny);
    ++compared_rows;
    for (int axis = 0; axis < 3; ++axis) {
      const std::size_t fid = face_index(id, axis);
      const int native_media = media_map.face_media[fid] != 0 ? 1 : 0;
      const double native_eps = native_media != 0 ? eps_in : eps_out;
      if (std::abs(native_eps - oracle_eps[axis]) <= 1.0e-12)
        continue;
      if (written_rows >= sample_limit)
        continue;
      const auto point = face_coord(i, j, k, axis);
      const DelphiFaceDebug debug =
          fid < media_map.face_debug.size() ? media_map.face_debug[fid]
                                            : DelphiFaceDebug{};
      out << i << ',' << j << ',' << k << ',' << axis_name(axis) << ','
          << native_eps << ',' << oracle_eps[axis] << ',' << native_media
          << ',' << oracle_media(oracle_eps[axis]) << ','
          << media_map.face_owner[fid] << ','
          << (fid < media_map.face_encoded.size()
                  ? media_map.face_encoded[fid]
                  : std::int32_t{0})
          << ',';
      write_optional_double(oracle_encoded[axis]);
      out << ',' << point[0] << ',' << point[1]
          << ',' << point[2] << ',' << debug.iteration << ','
          << delphi_face_action_name(debug.action) << ','
          << debug.nearest_atom << ',';
      write_optional_double(debug.nearest_atom_surface_dist);
      out << ',';
      write_optional_double(debug.nearest_exposed_dist);
      out << ',' << debug.raw_media << ',' << debug.raw_owner << ','
          << debug.raw_encoded << ',' << debug.final_encoded << '\n';
      ++written_rows;
    }
  }
  if (compared_rows == 0) {
    const std::string skipped =
        skipped_fixture_rows > 0 ? " fixture filter removed all rows." : "";
    throw std::runtime_error(
        "PBE: pbe_oracle_midpoint_epsilon.csv has no comparable rows for "
        "midpoint diff." +
        skipped);
  }
}

double AnalyzerImpl::compute_nonlinear_residual(
    const GridShape &shape, const AlignedDoubleVector &source,
    double screen_scale) const {
  return compute_nonlinear_residual_from_phi(shape, phi_, source,
                                             screen_scale);
}

double AnalyzerImpl::compute_nonlinear_residual_from_phi(
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

double AnalyzerImpl::compute_nonlinear_residual_vector(
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

AnalyzerImpl::NewtonCorrectionStats AnalyzerImpl::solve_newton_correction_sor(
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
                const std::size_t first_i =
                    1 + (static_cast<int>((1 + j + k) & 1ULL) == color ? 0U
                                                                        : 1U);
                for (std::size_t i = first_i; i + 1 < shape.nx; i += 2) {
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

void AnalyzerImpl::prepare_newton_correction_hierarchy(
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

AnalyzerImpl::SolveStats AnalyzerImpl::smooth_newton_correction_level(
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
          const std::size_t first_i =
              1 + (static_cast<int>((1 + j + k) & 1ULL) == color ? 0U : 1U);
          for (std::size_t i = first_i; i + 1 < shape.nx; i += 2) {
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

double AnalyzerImpl::compute_newton_correction_level_residual(
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

void AnalyzerImpl::restrict_newton_correction_residual(
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

void AnalyzerImpl::prolongate_add_newton_correction(
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

void AnalyzerImpl::run_newton_correction_vcycle(
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

AnalyzerImpl::NewtonCorrectionStats AnalyzerImpl::solve_newton_correction_vcycle(
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

AnalyzerImpl::NewtonCorrectionStats AnalyzerImpl::solve_newton_correction(
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

std::size_t AnalyzerImpl::update_nonlinear_source(
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

double AnalyzerImpl::grid_reaction_energy_kcal(
    const GridShape &shape, const AlignedDoubleVector &solvent_phi,
    const AlignedDoubleVector &vacuum_phi) const {
  if (shape.n == 0 || solvent_phi.size() != shape.n ||
      vacuum_phi.size() != shape.n || rho_.size() != shape.n)
    return 0.0;
  double sum = 0.0;
  for (std::size_t id = 0; id < shape.n; ++id)
    sum += rho_[id] * (solvent_phi[id] - vacuum_phi[id]);
  return 0.5 * kCoulombKcalPerMol * field_to_native_factor() * sum;
}

double AnalyzerImpl::grid_phase_energy_kcal(
    const GridShape &shape, const AlignedDoubleVector &field) const {
  return grid_phase_energy_kcal(shape, field, rho_);
}

double AnalyzerImpl::grid_phase_energy_kcal(
    const GridShape &shape, const AlignedDoubleVector &field,
    const AlignedDoubleVector &rho) const {
  if (shape.n == 0 || field.size() != shape.n || rho.size() != shape.n)
    return 0.0;
  double sum = 0.0;
  for (std::size_t id = 0; id < shape.n; ++id)
    sum += rho[id] * field[id];
  return 0.5 * kCoulombKcalPerMol * field_to_native_factor() * sum;
}

double AnalyzerImpl::atom_phase_energy_kcal(
    const std::vector<Atom> &atoms, const std::vector<double> &atom_phi) const {
  if (atoms.size() != atom_phi.size())
    return 0.0;
  double sum = 0.0;
  for (std::size_t i = 0; i < atoms.size(); ++i)
    sum += atoms[i].q * atom_phi[i];
  return 0.5 * kCoulombKcalPerMol * field_to_native_factor() * sum;
}

double AnalyzerImpl::analytic_coulomb_energy_kcal(
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

double AnalyzerImpl::effective_nl_relax() const noexcept {
  if (args_.nl_relax == "auto")
    return 0.5;
  return std::clamp(args_.nl_relax_value, 0.05, 1.0);
}

double AnalyzerImpl::nonlinear_alpha() const noexcept {
  return kCoulombKcalPerMol /
         (kKBKcalPerMolK * std::max(1.0, args_.temp));
}

double AnalyzerImpl::effective_nonlinear_clip() const noexcept {
  if (delphi_nwt_parity_active() && delphi_field_units_active())
    return 50.0;
  return std::max(1.0e-12, args_.nl_clip);
}

double AnalyzerImpl::rhs_scale(const GridShape &shape) const noexcept {
  return 4.0 * kPi / std::max(1.0e-12, shape.h);
}

bool AnalyzerImpl::nonlinear_dimensionless_mode() const noexcept {
  return args_.test_nl_source == "delphi-dimensionless";
}

double AnalyzerImpl::nonlinear_rhs_scale(const GridShape &shape,
                                     double alpha) const noexcept {
  const double base = effective_rhs_scale(shape);
  if (delphi_field_units_active())
    return base;
  return nonlinear_dimensionless_mode() ? alpha * base : base;
}

double AnalyzerImpl::nonlinear_argument_value(double field_value,
                                          double alpha) const noexcept {
  if (delphi_field_units_active())
    return field_value;
  return nonlinear_dimensionless_mode() ? field_value : alpha * field_value;
}

double AnalyzerImpl::nonlinear_source_value(double phi, double alpha,
                                        double clip) const noexcept {
  if (alpha <= 0.0)
    return 0.0;
  const double arg =
      nonlinear_argument_value(phi, alpha);
  const auto [sinh_arg, unused_cosh] = delphi_sinh_cosh_tail(arg, clip);
  (void)unused_cosh;
  if (delphi_field_units_active() ||
      args_.test_nl_source == "delphi-dimensionless")
    return sinh_arg;
  return sinh_arg / alpha;
}

double AnalyzerImpl::nonlinear_source_derivative(double phi, double alpha,
                                             double clip) const noexcept {
  if (alpha <= 0.0)
    return 0.0;
  const double arg =
      nonlinear_argument_value(phi, alpha);
  const auto [unused_sinh, cosh_arg] = delphi_sinh_cosh_tail(arg, clip);
  (void)unused_sinh;
  return cosh_arg;
}

void AnalyzerImpl::scale_nonlinear_field_to_solve_units(AlignedDoubleVector &field,
                                                    double alpha) const {
  if (!nonlinear_dimensionless_mode() || delphi_field_units_active() ||
      alpha <= 0.0)
    return;
  for (double &value : field)
    value *= alpha;
}

void AnalyzerImpl::scale_nonlinear_field_to_physical_units(
    AlignedDoubleVector &field, double alpha) const {
  if (!nonlinear_dimensionless_mode() || delphi_field_units_active() ||
      alpha <= 0.0)
    return;
  const double inv_alpha = 1.0 / alpha;
  for (double &value : field)
    value *= inv_alpha;
}

AnalyzerImpl::ParityAudit AnalyzerImpl::compute_parity_audit(
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

AnalyzerImpl::ForensicAudit
AnalyzerImpl::compute_forensic_audit(const GridShape &shape) const {
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

AnalyzerImpl::IscmAudit AnalyzerImpl::compute_iscm_audit(
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
  if (args_.model != "linear" && !delphi_nwt_parity_active()) {
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
  const double field_native_scale = field_to_native_factor();
  const double grid_scale = 1.0 / std::max(1.0e-12, shape.h);
  const double induced_scale =
      kDelphiReactionFieldEnergyFactor /
      (2.0 * std::max(1.0e-12, grid_scale) * std::max(1.0e-12, epkt));
  const std::size_t plane = shape.nx * shape.ny;
  const bool delphi_vdw_boundary =
      args_.surface == "delphi-vdw" && solvent_eps_x.size() == shape.n &&
      solvent_eps_y.size() == shape.n && solvent_eps_z.size() == shape.n;
  const bool native_iscm_positions =
      args_.surface == "delphi-vdw" &&
      last_iscm_surface_valid_.size() == shape.n &&
      last_iscm_surface_positions_.size() == shape.n;
  std::unordered_map<std::size_t, std::array<double, 3>>
      oracle_iscm_positions;
  if (std::filesystem::exists("pbe_oracle_iscm_full.csv")) {
    std::ifstream oracle("pbe_oracle_iscm_full.csv");
    std::string line;
    if (oracle.is_open() && std::getline(oracle, line)) {
      const auto header = csv_header_map(split_csv_simple(line));
      while (std::getline(oracle, line)) {
        if (trim_copy(line).empty())
          continue;
        const auto fields = split_csv_simple(line);
        std::size_t oi = 0, oj = 0, ok = 0;
        double sx = 0.0, sy = 0.0, sz = 0.0;
        if (!csv_index(fields, header, "i", oi) ||
            !csv_index(fields, header, "j", oj) ||
            !csv_index(fields, header, "k", ok) ||
            !csv_double(fields, header, {"surface_x", "surface_x_A"}, sx) ||
            !csv_double(fields, header, {"surface_y", "surface_y_A"}, sy) ||
            !csv_double(fields, header, {"surface_z", "surface_z_A"}, sz))
          continue;
        if (oi >= shape.nx || oj >= shape.ny || ok >= shape.nz)
          continue;
        oracle_iscm_positions[idx3(oi, oj, ok, shape.nx, shape.ny)] = {
            sx, sy, sz};
      }
    }
  }
  audit.position_source =
      native_iscm_positions ? "native_delphi_vdw_scaled" : "grid_point";
  audit.native_independent = true;
  double position_delta_sq = 0.0;
  std::size_t position_delta_count = 0;
  struct IscmPositionDiff {
    std::size_t id{0};
    std::size_t i{0};
    std::size_t j{0};
    std::size_t k{0};
    double native_x{0.0};
    double native_y{0.0};
    double native_z{0.0};
    double oracle_x{0.0};
    double oracle_y{0.0};
    double oracle_z{0.0};
    double delta{0.0};
  };
  std::vector<IscmPositionDiff> position_diffs;
  if (args_.parity_dump == "all" && !oracle_iscm_positions.empty())
    position_diffs.reserve(oracle_iscm_positions.size());
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

        const double phi_center =
            potential_scale * solvent_phi[id] * field_native_scale;
        const double phi_neigh =
            potential_scale *
            (solvent_phi[id - 1] + solvent_phi[id + 1] +
             solvent_phi[id - shape.nx] + solvent_phi[id + shape.nx] +
             solvent_phi[id - plane] + solvent_phi[id + plane]) /
            6.0 * field_native_scale;
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
        double x = shape.ox + static_cast<double>(i) * shape.h;
        double y = shape.oy + static_cast<double>(j) * shape.h;
        double z = shape.oz + static_cast<double>(k) * shape.h;
        if (native_iscm_positions && last_iscm_surface_valid_[id] != 0U) {
          const auto &native_pos = last_iscm_surface_positions_[id];
          x = native_pos[0];
          y = native_pos[1];
          z = native_pos[2];
        }
        const auto oracle_pos_it = oracle_iscm_positions.find(id);
        if (oracle_pos_it != oracle_iscm_positions.end()) {
          const double dx = x - oracle_pos_it->second[0];
          const double dy = y - oracle_pos_it->second[1];
          const double dz = z - oracle_pos_it->second[2];
          const double delta_sq = dx * dx + dy * dy + dz * dz;
          position_delta_sq += delta_sq;
          ++position_delta_count;
          if (args_.parity_dump == "all") {
            position_diffs.push_back({id,
                                      i,
                                      j,
                                      k,
                                      x,
                                      y,
                                      z,
                                      oracle_pos_it->second[0],
                                      oracle_pos_it->second[1],
                                      oracle_pos_it->second[2],
                                      std::sqrt(delta_sq)});
          }
        }
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
  audit.position_oracle_matches = position_delta_count;
  audit.position_rms_vs_oracle =
      position_delta_count > 0
          ? std::sqrt(position_delta_sq /
                      static_cast<double>(position_delta_count))
          : 0.0;
  if (!position_diffs.empty()) {
    std::sort(position_diffs.begin(), position_diffs.end(),
              [](const IscmPositionDiff &a, const IscmPositionDiff &b) {
                return a.delta > b.delta;
              });
    std::ofstream diff_out("pbe_iscm_position_diff.csv",
                           std::ios::out | std::ios::trunc);
    if (diff_out.is_open()) {
      diff_out << "rank,id,i,j,k,native_x,native_y,native_z,oracle_x,"
                  "oracle_y,oracle_z,delta_A\n";
      diff_out << std::setprecision(12);
      std::size_t rank = 0;
      for (const auto &diff : position_diffs) {
        if (++rank > 512 && diff.delta <= 1.0e-10)
          break;
        diff_out << rank << ',' << diff.id << ',' << diff.i << ','
                 << diff.j << ',' << diff.k << ',' << diff.native_x << ','
                 << diff.native_y << ',' << diff.native_z << ','
                 << diff.oracle_x << ',' << diff.oracle_y << ','
                 << diff.oracle_z << ',' << diff.delta << '\n';
      }
    }
  }
  audit.water_phase_kj_mol = audit.water_phase_kt * kt_kj;
  audit.water_phase_kcal_mol = audit.water_phase_kt * kt_kcal;
  audit.vacuum_phase_kcal_mol = grid_vacuum_phase_kcal_mol;
  audit.vacuum_phase_kj_mol = grid_vacuum_phase_kcal_mol * kKcalToKJ;
  audit.vacuum_phase_kt =
      kt_kcal > 0.0 ? grid_vacuum_phase_kcal_mol / kt_kcal : 0.0;
  audit.water_minus_vacuum_kt = audit.water_phase_kt - audit.vacuum_phase_kt;
  audit.water_minus_vacuum_kj_mol =
      audit.water_phase_kj_mol - audit.vacuum_phase_kj_mol;
  audit.water_minus_vacuum_kcal_mol =
      audit.water_phase_kcal_mol - audit.vacuum_phase_kcal_mol;
  // DelPhi's two-dielectric VDW ISCM path finalizes corrected RF from the
  // water induced-surface-charge phase alone. Keep the vacuum grid phase in
  // the audit so convention differences remain inspectable.
  audit.corrected_kt = audit.water_phase_kt;
  audit.corrected_kj_mol = audit.corrected_kt * kt_kj;
  audit.corrected_kcal_mol = audit.corrected_kt * kt_kcal;
  audit.finalization_mode = "water_only";
  audit.status = audit.boundary_points > 0 ? "ok" : "no_boundary_points";

  audit.samples = deterministic_samples;
  for (std::size_t idx = 0; idx < top_samples.size(); ++idx) {
    top_samples[idx].rank = idx + 1;
    audit.samples.push_back(top_samples[idx]);
  }
  return audit;
}

void AnalyzerImpl::copy_iscm_audit_to_solve(const IscmAudit &audit,
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
  solve.iscm_water_minus_vacuum_kcal_mol = audit.water_minus_vacuum_kcal_mol;
  solve.iscm_water_minus_vacuum_kj_mol = audit.water_minus_vacuum_kj_mol;
  solve.iscm_water_minus_vacuum_kt = audit.water_minus_vacuum_kt;
  solve.iscm_finalization_mode = audit.finalization_mode;
  solve.iscm_position_source = audit.position_source;
  solve.iscm_position_rms_vs_oracle = audit.position_rms_vs_oracle;
  solve.iscm_position_oracle_matches = audit.position_oracle_matches;
  solve.iscm_native_independent = audit.native_independent;
  solve.iscm_status = audit.status;
  solve.delphi_face_crossings = audit.delphi_face_crossings;
  solve.delphi_eps_face_matches = audit.delphi_eps_face_matches;
  solve.delphi_boundary_grid_points = audit.boundary_points;
  solve.delphi_charged_boundary_points = audit.charged_boundary_points;
  solve.delphi_iscm_status = audit.status;
}

AnalyzerImpl::NonlinearEnergyAudit AnalyzerImpl::compute_nonlinear_energy_audit(
    const GridShape &shape, const AlignedDoubleVector &solvent_phi,
    const AlignedDoubleVector &solvent_screen) const {
  NonlinearEnergyAudit audit;
  if (args_.model != "nonlinear") {
    audit.status = "linear_model";
    return audit;
  }
  if (args_.salt <= 0.0) {
    audit.status = "linearized_no_salt";
    return audit;
  }
  if (shape.n == 0 || solvent_phi.size() != shape.n) {
    audit.status = "missing_field";
    return audit;
  }

  const double ion_concentration = std::max(0.0, args_.salt);
  const double coeffs[5] = {
      -2.0 * ion_concentration,
      0.0,
      -ion_concentration / 3.0,
      0.0,
      -ion_concentration / 60.0,
  };
  const double d_coeffs[5] = {
      -coeffs[0] / 2.0,
      -coeffs[1] / 3.0,
      -coeffs[2] / 4.0,
      -coeffs[3] / 5.0,
      -coeffs[4] / 6.0,
  };
  const double alpha = nonlinear_alpha();
  const double cell_volume = shape.h * shape.h * shape.h;
  constexpr double kDelphiEnergyUnitConversion = 0.0006023;
  double rho_phi_accum = 0.0;
  double osmotic_accum = 0.0;

  for (std::size_t k = 0; k < shape.nz; ++k) {
    const double cut_x_dummy = 1.0;
    (void)cut_x_dummy;
    const double cut_z = (k == 0 || k + 1 == shape.nz) ? 0.5 : 1.0;
    for (std::size_t j = 0; j < shape.ny; ++j) {
      const double cut_y = cut_z * ((j == 0 || j + 1 == shape.ny) ? 0.5 : 1.0);
      for (std::size_t i = 0; i < shape.nx; ++i) {
        const std::size_t id = idx3(i, j, k, shape.nx, shape.ny);
        const bool ion_accessible =
            solvent_screen.size() == shape.n
                ? solvent_screen[id] > 0.0
                : (id < ion_accessible_.size() && ion_accessible_[id] != 0U);
        if (!ion_accessible)
          continue;
        const double cut_edge =
            cut_y * ((i == 0 || i + 1 == shape.nx) ? 0.5 : 1.0);
        const double phi =
            delphi_field_units_active() ? solvent_phi[id] : alpha * solvent_phi[id];
        const double charge_density =
            ((((coeffs[4] * phi + coeffs[3]) * phi + coeffs[2]) * phi +
              coeffs[1]) *
                 phi +
             coeffs[0]);
        const double osmotic_density =
            ((((d_coeffs[4] * phi + d_coeffs[3]) * phi + d_coeffs[2]) * phi +
              d_coeffs[1]) *
                 phi +
             d_coeffs[0]);
        const double ionic_charge =
            cut_edge * charge_density * phi * cell_volume;
        rho_phi_accum -= ionic_charge * phi;
        osmotic_accum += cut_edge * osmotic_density * phi * phi;
        ++audit.sampled_points;
      }
    }
  }

  double electrostatic_stress =
      rho_phi_accum * 0.5 * kDelphiEnergyUnitConversion;
  double osmotic_pressure =
      osmotic_accum * (-kDelphiEnergyUnitConversion * cell_volume);
  // pyDelPhi negates both nonlinear terms before adding them to the final
  // corrected RF energy.
  electrostatic_stress *= -1.0;
  osmotic_pressure *= -1.0;

  audit.electrostatic_stress_kt = electrostatic_stress;
  audit.osmotic_pressure_kt = osmotic_pressure;
  audit.total_correction_kt = electrostatic_stress + osmotic_pressure;
  const double kt_kj = kKBKJPerMolK * std::max(1.0, args_.temp);
  const double kt_kcal = kKBKcalPerMolK * std::max(1.0, args_.temp);
  audit.electrostatic_stress_kj_mol = audit.electrostatic_stress_kt * kt_kj;
  audit.electrostatic_stress_kcal_mol =
      audit.electrostatic_stress_kt * kt_kcal;
  audit.osmotic_pressure_kj_mol = audit.osmotic_pressure_kt * kt_kj;
  audit.osmotic_pressure_kcal_mol = audit.osmotic_pressure_kt * kt_kcal;
  audit.total_correction_kj_mol = audit.total_correction_kt * kt_kj;
  audit.total_correction_kcal_mol = audit.total_correction_kt * kt_kcal;
  audit.status = audit.sampled_points > 0 ? "ok" : "no_ion_accessible_points";
  return audit;
}

void AnalyzerImpl::copy_nonlinear_energy_audit_to_solve(
    const NonlinearEnergyAudit &audit, ReactionSolve &solve) const {
  solve.nonlinear_stress_kt = audit.electrostatic_stress_kt;
  solve.nonlinear_stress_kj_mol = audit.electrostatic_stress_kj_mol;
  solve.nonlinear_stress_kcal_mol = audit.electrostatic_stress_kcal_mol;
  solve.nonlinear_osmotic_kt = audit.osmotic_pressure_kt;
  solve.nonlinear_osmotic_kj_mol = audit.osmotic_pressure_kj_mol;
  solve.nonlinear_osmotic_kcal_mol = audit.osmotic_pressure_kcal_mol;
  solve.nonlinear_energy_correction_kt = audit.total_correction_kt;
  solve.nonlinear_energy_correction_kj_mol = audit.total_correction_kj_mol;
  solve.nonlinear_energy_correction_kcal_mol =
      audit.total_correction_kcal_mol;
  solve.nonlinear_energy_status = audit.status;
}

void AnalyzerImpl::update_newton_derivative_audit(
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

AnalyzerImpl::SolveStats AnalyzerImpl::solve_nonlinear_sor_block(
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
                const std::size_t first_i =
                    1 + (static_cast<int>((1 + j + k) & 1ULL) == color ? 0U
                                                                        : 1U);
                for (std::size_t i = first_i; i + 1 < shape.nx; i += 2) {
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

bool AnalyzerImpl::nwt_half_shape_matches(const GridShape &shape) const
    noexcept {
  const GridShape &cached = nwt_half_workspace_.shape;
  return nwt_half_workspace_.active && cached.n == shape.n &&
         cached.nx == shape.nx && cached.ny == shape.ny &&
         cached.nz == shape.nz;
}

void AnalyzerImpl::invalidate_nwt_half_workspace() noexcept {
  nwt_half_workspace_.active = false;
  nwt_half_workspace_.coeff_ready = false;
  nwt_half_workspace_.packed_ready = false;
  nwt_half_workspace_.row_span_ready = false;
  nwt_half_workspace_.salt_tile_ready = false;
  nwt_half_workspace_.rhs_scale = 0.0;
}

void AnalyzerImpl::prepare_nwt_half_workspace(const GridShape &shape,
                                              double rhs_scale) {
  const std::size_t even_count = (shape.n + 1) / 2;
  const std::size_t odd_count = shape.n / 2;
  const bool shape_match = nwt_half_workspace_.shape.n == shape.n &&
                           nwt_half_workspace_.shape.nx == shape.nx &&
                           nwt_half_workspace_.shape.ny == shape.ny &&
                           nwt_half_workspace_.shape.nz == shape.nz;
  if (!shape_match ||
      nwt_half_workspace_.even_phi.size() != even_count ||
      nwt_half_workspace_.odd_phi.size() != odd_count) {
    nwt_half_workspace_.shape = shape;
    nwt_half_workspace_.even_phi.assign(even_count, 0.0);
    nwt_half_workspace_.odd_phi.assign(odd_count, 0.0);
    nwt_half_workspace_.rhs_even.assign(even_count, 0.0);
    nwt_half_workspace_.rhs_odd.assign(odd_count, 0.0);
    nwt_half_workspace_.eps_sum_even.assign(even_count, 1.0);
    nwt_half_workspace_.eps_sum_odd.assign(odd_count, 1.0);
    nwt_half_workspace_.screen_even.assign(even_count, 0.0);
    nwt_half_workspace_.screen_odd.assign(odd_count, 0.0);
    nwt_half_workspace_.homo_eps_even.assign(even_count, 0.0);
    nwt_half_workspace_.homo_eps_odd.assign(odd_count, 0.0);
    nwt_half_workspace_.coeff_ready = false;
    nwt_half_workspace_.packed_ready = false;
    nwt_half_workspace_.row_span_ready = false;
    nwt_half_workspace_.salt_tile_ready = false;
  }
  if (!nwt_half_workspace_.active) {
    for (std::size_t id = 0; id < shape.n; ++id) {
      if ((id & 1U) == 0U)
        nwt_half_workspace_.even_phi[id / 2] = phi_[id];
      else
        nwt_half_workspace_.odd_phi[id / 2] = phi_[id];
    }
    nwt_half_workspace_.active = true;
  }
  if (!nwt_half_workspace_.coeff_ready ||
      std::abs(nwt_half_workspace_.rhs_scale - rhs_scale) > 1.0e-30) {
    const auto pack_begin = std::chrono::steady_clock::now();
    for (std::size_t id = 0; id < shape.n; ++id) {
      const std::size_t h = id / 2;
      const double eps_sum = std::max(1.0e-12, denom_[id] - screen_[id]);
      if ((id & 1U) == 0U) {
        nwt_half_workspace_.rhs_even[h] = rhs_scale * rho_[id];
        nwt_half_workspace_.eps_sum_even[h] = eps_sum;
        nwt_half_workspace_.screen_even[h] = screen_[id];
      } else {
        nwt_half_workspace_.rhs_odd[h] = rhs_scale * rho_[id];
        nwt_half_workspace_.eps_sum_odd[h] = eps_sum;
        nwt_half_workspace_.screen_odd[h] = screen_[id];
      }
    }

    nwt_half_workspace_.even_homo_cells.clear();
    nwt_half_workspace_.even_hetero_cells.clear();
    nwt_half_workspace_.odd_homo_cells.clear();
    nwt_half_workspace_.odd_hetero_cells.clear();
    nwt_half_workspace_.even_homo_linear_spans.clear();
    nwt_half_workspace_.even_homo_salt_spans.clear();
    nwt_half_workspace_.odd_homo_linear_spans.clear();
    nwt_half_workspace_.odd_homo_salt_spans.clear();
    nwt_half_workspace_.even_homo_linear_tiles.clear();
    nwt_half_workspace_.even_homo_salt_tiles.clear();
    nwt_half_workspace_.odd_homo_linear_tiles.clear();
    nwt_half_workspace_.odd_homo_salt_tiles.clear();
    nwt_half_workspace_.homo_cells = 0;
    nwt_half_workspace_.homo_linear_cells = 0;
    nwt_half_workspace_.homo_salt_cells = 0;
    if (shape.n <= static_cast<std::size_t>(
                       std::numeric_limits<std::uint32_t>::max())) {
      const std::size_t xy = shape.nx * shape.ny;
      const std::size_t reserve_each =
          shape.n > 0 ? std::max<std::size_t>(1, shape.n / 4) : 0;
      nwt_half_workspace_.even_homo_cells.reserve(reserve_each);
      nwt_half_workspace_.even_hetero_cells.reserve(reserve_each);
      nwt_half_workspace_.odd_homo_cells.reserve(reserve_each);
      nwt_half_workspace_.odd_hetero_cells.reserve(reserve_each);
      const auto append_span =
          [](std::vector<NwtHalfWorkspace::HomoSpan> &spans,
             std::uint32_t h) {
            if (!spans.empty()) {
              auto &last = spans.back();
              if (last.begin + last.count == h) {
                ++last.count;
                return;
              }
            }
            spans.push_back({h, 1});
          };
      auto add_cell = [&](int color, std::size_t id) {
        const double ex_p = eps_x_[id];
        const double ex_m = eps_x_[id - 1];
        const double ey_p = eps_y_[id];
        const double ey_m = eps_y_[id - shape.nx];
        const double ez_p = eps_z_[id];
        const double ez_m = eps_z_[id - xy];
        NwtHalfWorkspace::PackedCell cell{};
        cell.h = static_cast<std::uint32_t>(id / 2);
        cell.xp = static_cast<std::uint32_t>((id + 1) / 2);
        cell.xm = static_cast<std::uint32_t>((id - 1) / 2);
        cell.yp = static_cast<std::uint32_t>((id + shape.nx) / 2);
        cell.ym = static_cast<std::uint32_t>((id - shape.nx) / 2);
        cell.zp = static_cast<std::uint32_t>((id + xy) / 2);
        cell.zm = static_cast<std::uint32_t>((id - xy) / 2);
        cell.rhs = rhs_scale * rho_[id];
        cell.eps_sum = std::max(1.0e-12, denom_[id] - screen_[id]);
        cell.screen = screen_[id];
        cell.ex_p = ex_p;
        cell.ex_m = ex_m;
        cell.ey_p = ey_p;
        cell.ey_m = ey_m;
        cell.ez_p = ez_p;
        cell.ez_m = ez_m;
        const bool homogeneous = ex_p == ex_m && ex_p == ey_p &&
                                 ex_p == ey_m && ex_p == ez_p &&
                                 ex_p == ez_m;
        if (color == 0) {
          (homogeneous ? nwt_half_workspace_.even_homo_cells
                       : nwt_half_workspace_.even_hetero_cells)
              .push_back(cell);
          if (homogeneous) {
            nwt_half_workspace_.homo_eps_even[cell.h] = ex_p;
            if (cell.screen > 0.0) {
              append_span(nwt_half_workspace_.even_homo_salt_spans, cell.h);
              ++nwt_half_workspace_.homo_salt_cells;
            } else {
              append_span(nwt_half_workspace_.even_homo_linear_spans, cell.h);
              ++nwt_half_workspace_.homo_linear_cells;
            }
            ++nwt_half_workspace_.homo_cells;
          }
        } else {
          (homogeneous ? nwt_half_workspace_.odd_homo_cells
                       : nwt_half_workspace_.odd_hetero_cells)
              .push_back(cell);
          if (homogeneous) {
            nwt_half_workspace_.homo_eps_odd[cell.h] = ex_p;
            if (cell.screen > 0.0) {
              append_span(nwt_half_workspace_.odd_homo_salt_spans, cell.h);
              ++nwt_half_workspace_.homo_salt_cells;
            } else {
              append_span(nwt_half_workspace_.odd_homo_linear_spans, cell.h);
              ++nwt_half_workspace_.homo_linear_cells;
            }
            ++nwt_half_workspace_.homo_cells;
          }
        }
      };
      for (std::size_t k = 1; k + 1 < shape.nz; ++k) {
        for (std::size_t j = 1; j + 1 < shape.ny; ++j) {
          for (int color = 0; color < 2; ++color) {
            const std::size_t first_i =
                1 + (static_cast<int>((1 + j + k) & 1ULL) == color ? 0U
                                                                    : 1U);
            for (std::size_t i = first_i; i + 1 < shape.nx; i += 2) {
              add_cell(color, idx3(i, j, k, shape.nx, shape.ny));
            }
          }
        }
      }
      const auto build_tiles =
          [](const std::vector<NwtHalfWorkspace::HomoSpan> &spans,
             std::vector<NwtHalfWorkspace::HomoTile> &tiles) {
            constexpr std::uint32_t kTargetTileCells = 2048;
            tiles.clear();
            for (const auto &span : spans) {
              std::uint32_t begin = span.begin;
              std::uint32_t remaining = span.count;
              while (remaining > 0) {
                const std::uint32_t take =
                    std::min<std::uint32_t>(remaining, kTargetTileCells);
                if (!tiles.empty()) {
                  auto &last = tiles.back();
                  const bool contiguous = last.begin + last.count == begin;
                  const bool has_room =
                      last.count < kTargetTileCells &&
                      take <= kTargetTileCells - last.count;
                  if (contiguous && has_room) {
                    last.count += take;
                  } else {
                    tiles.push_back({begin, take});
                  }
                } else {
                  tiles.push_back({begin, take});
                }
                begin += take;
                remaining -= take;
              }
            }
          };
      build_tiles(nwt_half_workspace_.even_homo_linear_spans,
                  nwt_half_workspace_.even_homo_linear_tiles);
      build_tiles(nwt_half_workspace_.even_homo_salt_spans,
                  nwt_half_workspace_.even_homo_salt_tiles);
      build_tiles(nwt_half_workspace_.odd_homo_linear_spans,
                  nwt_half_workspace_.odd_homo_linear_tiles);
      build_tiles(nwt_half_workspace_.odd_homo_salt_spans,
                  nwt_half_workspace_.odd_homo_salt_tiles);
      nwt_half_workspace_.packed_ready = true;
      nwt_half_workspace_.row_span_ready = true;
      nwt_half_workspace_.salt_tile_ready = true;
    } else {
      nwt_half_workspace_.packed_ready = false;
      nwt_half_workspace_.row_span_ready = false;
      nwt_half_workspace_.salt_tile_ready = false;
    }
    nwt_half_workspace_.rhs_scale = rhs_scale;
    nwt_half_workspace_.coeff_ready = true;
    add_timing_ms(TimingBucket::NwtCoeffPack,
                  std::chrono::duration<double, std::milli>(
                      std::chrono::steady_clock::now() - pack_begin)
                      .count());
  }
}

void AnalyzerImpl::flush_nwt_half_workspace(const GridShape &shape) {
  if (!nwt_half_shape_matches(shape))
    return;
  for (std::size_t id = 0; id < shape.n; ++id) {
    phi_[id] = (id & 1U) == 0U
                   ? nwt_half_workspace_.even_phi[id / 2]
                   : nwt_half_workspace_.odd_phi[id / 2];
  }
  nwt_half_workspace_.active = false;
}

AnalyzerImpl::SolveStats AnalyzerImpl::solve_newton_nwt_sor_block(
    const GridShape &shape, int sweeps, double screen_scale,
    TitanExecutor &executor, std::size_t &clipped_points) {
  const auto update_time_begin = std::chrono::steady_clock::now();
  SolveStats stats;
  if (shape.n == 0 || sweeps <= 0)
    return stats;

  double omega =
      args_.omega == "auto"
          ? std::clamp(2.0 / (1.0 + std::sin(kPi / static_cast<double>(
                                             std::max({shape.nx, shape.ny,
                                                       shape.nz})))), 1.0,
                       1.95)
          : args_.omega_value;
  if (args_.omega == "auto" && args_.test_operator_convention == "delphi" &&
      !args_.test_pqr.empty()) {
    omega = delphi_relaxation_omega(shape);
  }
  if (args_.omega == "auto" && delphi_nwt_parity_active()) {
    // pyDelPhi NWT defaults to nwt_adaptive_omega=1.0; keep this parity path
    // separate from the linear DelPhi SOR relaxation convention.
    omega = 1.0;
  }

  const std::size_t interior_z = shape.nz > 2 ? shape.nz - 2 : 0;
  const std::size_t workers = effective_workers(executor, interior_z);
  const double rhs = nonlinear_rhs_scale(shape, nonlinear_alpha());
  const double alpha = nonlinear_alpha();
  const double clip = effective_nonlinear_clip();
  const std::size_t xy = shape.nx * shape.ny;
  std::vector<double> local_sumsq(workers, 0.0);
  std::vector<double> local_max(workers, 0.0);
  std::vector<std::size_t> local_count(workers, 0);
  std::vector<std::size_t> local_clipped(workers, 0);
  clipped_points = 0;
  const bool delphi_block_stop = delphi_nwt_parity_active();
  const std::size_t delphi_rms_denominator =
      delphi_block_stop ? (shape.n + 1) / 2 : 0;
  const bool write_trace = !args_.test_pqr.empty() &&
                           args_.parity_dump == "all" &&
                           delphi_nwt_parity_active();
  const bool use_auto_half_kernel =
      args_.nwt_kernel == "auto" && args_.nwt_progress == "off" &&
      delphi_nwt_parity_active();
  const bool use_forced_half_kernel = args_.nwt_kernel == "delphi-half";
  if ((use_auto_half_kernel || use_forced_half_kernel) &&
      args_.nwt_progress == "off" && delphi_nwt_parity_active() && !write_trace &&
      (shape.nx & 1U) != 0U && (shape.ny & 1U) != 0U &&
      (shape.nz & 1U) != 0U) {
    return solve_newton_nwt_delphi_half_block(shape, sweeps, screen_scale,
                                              executor, clipped_points);
  }

  std::vector<std::size_t> trace_ids;
  if (write_trace && shape.nx > 2 && shape.ny > 2 && shape.nz > 2) {
    trace_ids = {
        idx3(shape.nx / 2, shape.ny / 2, shape.nz / 2, shape.nx, shape.ny),
        idx3(shape.nx / 2 - 1, shape.ny / 2, shape.nz / 2, shape.nx,
             shape.ny),
        idx3(shape.nx / 2, shape.ny / 2 - 1, shape.nz / 2, shape.nx,
             shape.ny),
        idx3(shape.nx / 2, shape.ny / 2, shape.nz / 2 - 1, shape.nx,
             shape.ny),
        idx3(1, 1, 1, shape.nx, shape.ny),
        idx3(shape.nx - 2, shape.ny - 2, shape.nz - 2, shape.nx, shape.ny),
        idx3(1, shape.ny / 2, shape.nz / 2, shape.nx, shape.ny),
        idx3(shape.nx - 2, shape.ny / 2, shape.nz / 2, shape.nx,
             shape.ny),
        idx3(shape.nx / 2, shape.ny / 2, 1, shape.nx, shape.ny),
        idx3(shape.nx / 2, shape.ny / 2, shape.nz - 2, shape.nx,
             shape.ny)};
    for (std::size_t id = 0; id < rho_.size(); ++id) {
      if (std::abs(rho_[id]) > 0.0) {
        trace_ids.push_back(id);
      }
    }
    std::sort(trace_ids.begin(), trace_ids.end());
    trace_ids.erase(std::unique(trace_ids.begin(), trace_ids.end()),
                    trace_ids.end());
  }
  std::ofstream trace_out;
  if (write_trace) {
    if (!nwt_trace_header_written_) {
      trace_out.open("pbe_nwt_trace.csv", std::ios::out | std::ios::trunc);
      if (!trace_out.is_open())
        throw std::runtime_error("PBE: failed to open pbe_nwt_trace.csv.");
      trace_out
          << "iteration,color,i,j,k,omega,old_phi,neighbor_sum,"
             "charge_density,eps_sum,screen,sinh_phi,cosh_phi,"
             "linearized_rhs,denom,target_phi,new_phi,dphi,rms,max_dphi,"
             "field_units,operator_convention,charge_grid_source\n";
      nwt_trace_header_written_ = true;
    } else {
      trace_out.open("pbe_nwt_trace.csv", std::ios::out | std::ios::app);
      if (!trace_out.is_open())
        throw std::runtime_error("PBE: failed to append pbe_nwt_trace.csv.");
    }
    trace_out << std::setprecision(12);
  }

  for (int iter = 1; iter <= sweeps; ++iter) {
    std::fill(local_sumsq.begin(), local_sumsq.end(), 0.0);
    std::fill(local_max.begin(), local_max.end(), 0.0);
    std::fill(local_count.begin(), local_count.end(), 0);
    std::fill(local_clipped.begin(), local_clipped.end(), 0);

    for (int color = 0; color < 2; ++color) {
      std::vector<std::vector<std::string>> trace_lines;
      if (write_trace && iter <= 4)
        trace_lines.resize(workers);
      executor.parallel_for(
          interior_z, workers, 1,
          [&](std::size_t worker, std::size_t begin, std::size_t end) {
            double sumsq = 0.0;
            double max_delta = 0.0;
            std::size_t count = 0;
            std::size_t clipped = 0;
            for (std::size_t kk = begin; kk < end; ++kk) {
              const std::size_t k = kk + 1;
              for (std::size_t j = 1; j + 1 < shape.ny; ++j) {
                const std::size_t first_i =
                    1 + (static_cast<int>((1 + j + k) & 1ULL) == color ? 0U
                                                                        : 1U);
                for (std::size_t i = first_i; i + 1 < shape.nx; i += 2) {
                  const std::size_t id = idx3(i, j, k, shape.nx, shape.ny);
                  const double neighbor_sum =
                      eps_x_[id] * phi_[id + 1] +
                      eps_x_[id - 1] * phi_[id - 1] +
                      eps_y_[id] * phi_[id + shape.nx] +
                      eps_y_[id - shape.nx] * phi_[id - shape.nx] +
                      eps_z_[id] * phi_[id + xy] +
                      eps_z_[id - xy] * phi_[id - xy];
                  const double eps_sum =
                      std::max(1.0e-12, denom_[id] - screen_[id]);
                  const double active_screen = screen_[id] * screen_scale;
                  double diagonal = eps_sum;
                  double nonlinear_rhs = 0.0;
                  double sinh_arg = 0.0;
                  double cosh_arg = 1.0;
                  if (active_screen > 0.0) {
                    const double raw_arg =
                        nonlinear_argument_value(phi_[id], alpha);
                    if (raw_arg > clip || raw_arg < -clip)
                      ++clipped;
                    std::tie(sinh_arg, cosh_arg) =
                        delphi_sinh_cosh_tail(raw_arg, clip);
                    const double source =
                        (delphi_field_units_active() ||
                         args_.test_nl_source == "delphi-dimensionless")
                            ? sinh_arg
                            : sinh_arg / std::max(1.0e-30, alpha);
                    diagonal += active_screen * cosh_arg;
                    nonlinear_rhs = active_screen * (cosh_arg * phi_[id] -
                                                     source);
                  }
                  const double target =
                      (neighbor_sum + rhs * rho_[id] + nonlinear_rhs) /
                      std::max(1.0e-12, diagonal);
                  const double old_phi = phi_[id];
                  const double updated = old_phi + omega * (target - old_phi);
                  const double delta = updated - old_phi;
                  phi_[id] = updated;
                  if (write_trace && iter <= 4 &&
                      std::binary_search(trace_ids.begin(), trace_ids.end(),
                                         id)) {
                    const double charge_density = rhs * rho_[id];
                    std::ostringstream row;
                    row << std::setprecision(12) << iter << ',' << color << ','
                        << i << ',' << j << ',' << k << ',' << omega << ','
                        << old_phi << ',' << neighbor_sum << ','
                        << charge_density << ',' << eps_sum << ','
                        << active_screen << ',' << sinh_arg << ','
                        << cosh_arg << ',' << nonlinear_rhs << ','
                        << std::max(1.0e-12, diagonal) << ',' << target << ','
                        << updated << ',' << delta << ",0,0,"
                        << csv_escape(delphi_field_units_active()
                                          ? "delphi-dimensionless"
                                          : args_.test_field_units)
                        << ',' << csv_escape("delphi_nwt") << ','
                        << csv_escape(args_.test_charge_grid) << '\n';
                    trace_lines[worker].push_back(row.str());
                  }
                  if (!delphi_block_stop ||
                      (iter == sweeps && color == 1)) {
                    sumsq += delta * delta;
                    max_delta = std::max(max_delta, std::abs(delta));
                    ++count;
                  }
                }
              }
            }
            local_sumsq[worker] += sumsq;
            local_max[worker] = std::max(local_max[worker], max_delta);
            local_count[worker] += count;
            local_clipped[worker] += clipped;
          });
      if (write_trace && iter <= 4) {
        for (const auto &worker_lines : trace_lines) {
          for (const std::string &row : worker_lines)
            trace_out << row;
        }
      }
    }

    const auto reduce_time_begin = std::chrono::steady_clock::now();
    double sumsq = 0.0;
    double max_delta = 0.0;
    std::size_t count = 0;
    std::size_t clipped = 0;
    for (std::size_t worker = 0; worker < workers; ++worker) {
      sumsq += local_sumsq[worker];
      max_delta = std::max(max_delta, local_max[worker]);
      count += local_count[worker];
      clipped += local_clipped[worker];
    }
    clipped_points = std::max(clipped_points, clipped);
    stats.iterations = iter;
    if (delphi_block_stop) {
      stats.rms = delphi_rms_denominator > 0
                      ? std::sqrt(sumsq /
                                  static_cast<double>(delphi_rms_denominator))
                      : 0.0;
    } else {
      stats.rms =
          count > 0 ? std::sqrt(sumsq / static_cast<double>(count)) : 0.0;
    }
    stats.max_dphi = max_delta;
    if (write_trace && iter <= 4) {
      trace_out << iter << ",summary,0,0,0," << omega
                << ",0,0,0,0,0,0,1,0,0,0,0,0," << stats.rms << ','
                << stats.max_dphi << ','
                << csv_escape(delphi_field_units_active()
                                  ? "delphi-dimensionless"
                                  : args_.test_field_units)
                << ',' << csv_escape("delphi_nwt") << ','
                << csv_escape(args_.test_charge_grid) << '\n';
    }
    add_timing_ms(TimingBucket::NwtRmsReduce,
                  std::chrono::duration<double, std::milli>(
                      std::chrono::steady_clock::now() - reduce_time_begin)
                      .count());
    if (!delphi_block_stop &&
        (stats.max_dphi <= args_.nl_dphi || stats.rms <= args_.nl_residual)) {
      stats.converged = true;
      break;
    }
  }
  if (delphi_block_stop &&
      (stats.max_dphi < args_.nl_dphi || stats.rms < args_.nl_residual))
    stats.converged = true;
  add_timing_ms(TimingBucket::NwtHalfUpdate,
                std::chrono::duration<double, std::milli>(
                    std::chrono::steady_clock::now() - update_time_begin)
                    .count());
  return stats;
}

AnalyzerImpl::SolveStats AnalyzerImpl::solve_newton_nwt_delphi_half_block(
    const GridShape &shape, int sweeps, double screen_scale,
    TitanExecutor &executor, std::size_t &clipped_points) {
  const auto update_time_begin = std::chrono::steady_clock::now();
  SolveStats stats;
  if (shape.n == 0 || sweeps <= 0)
    return stats;

  double omega =
      args_.omega == "auto"
          ? std::clamp(2.0 / (1.0 + std::sin(kPi / static_cast<double>(
                                             std::max({shape.nx, shape.ny,
                                                       shape.nz})))), 1.0,
                       1.95)
          : args_.omega_value;
  if (args_.omega == "auto" && args_.test_operator_convention == "delphi" &&
      !args_.test_pqr.empty()) {
    omega = delphi_relaxation_omega(shape);
  }
  if (args_.omega == "auto" && delphi_nwt_parity_active())
    omega = 1.0;

  const std::size_t interior_z = shape.nz > 2 ? shape.nz - 2 : 0;
  const std::size_t workers = effective_workers(executor, interior_z);
  const double rhs = nonlinear_rhs_scale(shape, nonlinear_alpha());
  prepare_nwt_half_workspace(shape, rhs);
  const double alpha = nonlinear_alpha();
  const double clip = effective_nonlinear_clip();
  const std::size_t xy = shape.nx * shape.ny;
  const bool delphi_block_stop = delphi_nwt_parity_active();
  const std::size_t delphi_rms_denominator =
      delphi_block_stop ? (shape.n + 1) / 2 : 0;
  if (nwt_half_workspace_.local_sumsq.size() != workers) {
    nwt_half_workspace_.local_sumsq.assign(workers, 0.0);
    nwt_half_workspace_.local_max.assign(workers, 0.0);
    nwt_half_workspace_.local_count.assign(workers, 0);
    nwt_half_workspace_.local_clipped.assign(workers, 0);
  }
  auto &local_sumsq = nwt_half_workspace_.local_sumsq;
  auto &local_max = nwt_half_workspace_.local_max;
  auto &local_count = nwt_half_workspace_.local_count;
  auto &local_clipped = nwt_half_workspace_.local_clipped;
  auto &even_phi = nwt_half_workspace_.even_phi;
  auto &odd_phi = nwt_half_workspace_.odd_phi;
  const auto &rhs_even = nwt_half_workspace_.rhs_even;
  const auto &rhs_odd = nwt_half_workspace_.rhs_odd;
  const auto &eps_sum_even = nwt_half_workspace_.eps_sum_even;
  const auto &eps_sum_odd = nwt_half_workspace_.eps_sum_odd;
  const auto &screen_even = nwt_half_workspace_.screen_even;
  const auto &screen_odd = nwt_half_workspace_.screen_odd;
  const auto &homo_eps_even = nwt_half_workspace_.homo_eps_even;
  const auto &homo_eps_odd = nwt_half_workspace_.homo_eps_odd;
  clipped_points = 0;
  if (active_timing_) {
    active_timing_->nwt_homo_cells = nwt_half_workspace_.homo_cells;
    active_timing_->nwt_homo_linear_cells =
        nwt_half_workspace_.homo_linear_cells;
    active_timing_->nwt_homo_salt_cells =
        nwt_half_workspace_.homo_salt_cells;
  }
  const bool argument_is_dimensionless =
      delphi_field_units_active() || nonlinear_dimensionless_mode();
  const double argument_scale = argument_is_dimensionless ? 1.0 : alpha;
  const bool use_delphi_source_units =
      delphi_field_units_active() ||
      args_.test_nl_source == "delphi-dimensionless";
  const double source_scale =
      use_delphi_source_units ? 1.0 : 1.0 / std::max(1.0e-30, alpha);
  const double effective_clip = std::max(1.0e-12, clip);
  const double sinh_clip = std::sinh(effective_clip);

  if (nwt_half_workspace_.packed_ready) {
    const auto update_cells =
        [&](const std::vector<NwtHalfWorkspace::PackedCell> &cells,
            AlignedDoubleVector &write_phi, const AlignedDoubleVector &read_phi,
            bool homogeneous, bool accumulate_rms) {
          if (cells.empty())
            return;
          executor.parallel_for(
              cells.size(), workers, 256,
              [&](std::size_t worker, std::size_t begin, std::size_t end) {
                double sumsq = 0.0;
                double max_delta = 0.0;
                std::size_t count = 0;
                std::size_t clipped = 0;
                for (std::size_t cell_index = begin; cell_index < end;
                     ++cell_index) {
                  const auto &cell = cells[cell_index];
                  const double old_phi = write_phi[cell.h];
                  double neighbor_sum = 0.0;
                  if (homogeneous) {
                    neighbor_sum =
                        cell.ex_p *
                        (read_phi[cell.xp] + read_phi[cell.xm] +
                         read_phi[cell.yp] + read_phi[cell.ym] +
                         read_phi[cell.zp] + read_phi[cell.zm]);
                  } else {
                    neighbor_sum =
                        cell.ex_p * read_phi[cell.xp] +
                        cell.ex_m * read_phi[cell.xm] +
                        cell.ey_p * read_phi[cell.yp] +
                        cell.ey_m * read_phi[cell.ym] +
                        cell.ez_p * read_phi[cell.zp] +
                        cell.ez_m * read_phi[cell.zm];
                  }
                  const double active_screen = cell.screen * screen_scale;
                  double diagonal = cell.eps_sum;
                  double nonlinear_rhs = 0.0;
                  if (active_screen > 0.0) {
                    const double raw_arg = old_phi * argument_scale;
                    if (raw_arg > effective_clip || raw_arg < -effective_clip)
                      ++clipped;
                    const auto sinh_cosh =
                        delphi_sinh_cosh_tail_precomputed(
                            raw_arg, effective_clip, sinh_clip);
                    const double source = sinh_cosh.first * source_scale;
                    diagonal += active_screen * sinh_cosh.second;
                    nonlinear_rhs =
                        active_screen *
                        (sinh_cosh.second * old_phi - source);
                  }
                  const double target =
                      (neighbor_sum + cell.rhs + nonlinear_rhs) /
                      std::max(1.0e-12, diagonal);
                  const double updated = old_phi + omega * (target - old_phi);
                  const double delta = updated - old_phi;
                  write_phi[cell.h] = updated;
                  if (accumulate_rms) {
                    sumsq += delta * delta;
                    max_delta = std::max(max_delta, std::abs(delta));
                    ++count;
                  }
                }
                local_sumsq[worker] += sumsq;
                local_max[worker] = std::max(local_max[worker], max_delta);
                local_count[worker] += count;
                local_clipped[worker] += clipped;
              });
        };

    auto run_timed_cells =
        [&](const std::vector<NwtHalfWorkspace::PackedCell> &cells,
            AlignedDoubleVector &write_phi, const AlignedDoubleVector &read_phi,
            bool homogeneous, bool accumulate_rms, TimingBucket bucket) {
          const auto begin = std::chrono::steady_clock::now();
          if (!cells.empty() && active_timing_) {
            ++active_timing_->nwt_dispatch_calls;
            if (bucket == TimingBucket::NwtHomoUpdate)
              active_timing_->nwt_homo_updates += cells.size();
          }
          update_cells(cells, write_phi, read_phi, homogeneous, accumulate_rms);
          add_timing_ms(bucket,
                        std::chrono::duration<double, std::milli>(
                            std::chrono::steady_clock::now() - begin)
                            .count());
        };

    const bool row_span_requested =
        args_.nwt_homo_kernel == "row-span" ||
        args_.nwt_homo_kernel == "avx2" ||
        args_.nwt_homo_kernel == "auto";
#if TITAN_SIMD_X86
    const bool row_span_avx2 =
        row_span_requested && args_.nwt_homo_kernel != "row-span" &&
        titan_simd_dispatch::runtime_has_avx2();
#else
    const bool row_span_avx2 = false;
#endif
    const bool salt_tiled_requested =
        args_.nwt_homo_kernel == "salt-tiled" ||
        args_.nwt_homo_kernel == "auto";
#if TITAN_SIMD_X86
    const bool salt_tiled_avx2 =
        salt_tiled_requested && titan_simd_dispatch::runtime_has_avx2();
#else
    const bool salt_tiled_avx2 = false;
#endif
    if (salt_tiled_requested && nwt_half_workspace_.salt_tile_ready) {
      const auto tile_cell_count =
          [](const std::vector<NwtHalfWorkspace::HomoTile> &tiles) {
            std::uint64_t total = 0;
            for (const auto &tile : tiles)
              total += tile.count;
            return total;
          };
      const std::uint64_t even_homo_updates =
          tile_cell_count(nwt_half_workspace_.even_homo_linear_tiles) +
          tile_cell_count(nwt_half_workspace_.even_homo_salt_tiles);
      const std::uint64_t odd_homo_updates =
          tile_cell_count(nwt_half_workspace_.odd_homo_linear_tiles) +
          tile_cell_count(nwt_half_workspace_.odd_homo_salt_tiles);
      const bool collect_nwt_detail = args_.timing == "detail";
      std::vector<double> local_stencil_ms(workers, 0.0);
      std::vector<double> local_salt_ms(workers, 0.0);
      std::vector<double> local_write_ms(workers, 0.0);
      std::vector<double> local_worker_ms(workers, 0.0);

      const auto update_homo_tiles =
          [&](const std::vector<NwtHalfWorkspace::HomoTile> &linear_tiles,
              const std::vector<NwtHalfWorkspace::HomoTile> &salt_tiles,
              AlignedDoubleVector &write_phi,
              const AlignedDoubleVector &read_phi,
              const AlignedDoubleVector &rhs_vec,
              const AlignedDoubleVector &eps_sum_vec,
              const AlignedDoubleVector &screen_vec,
              const AlignedDoubleVector &eps_vec, int color,
              bool accumulate_rms, std::uint64_t update_count) {
            if (linear_tiles.empty() && salt_tiles.empty())
              return;
            if (collect_nwt_detail) {
              std::fill(local_stencil_ms.begin(), local_stencil_ms.end(), 0.0);
              std::fill(local_salt_ms.begin(), local_salt_ms.end(), 0.0);
              std::fill(local_write_ms.begin(), local_write_ms.end(), 0.0);
              std::fill(local_worker_ms.begin(), local_worker_ms.end(), 0.0);
            }
            const auto begin_time = std::chrono::steady_clock::now();
            const NwtHomoOffsets offsets =
                nwt_homo_offsets(shape.nx, shape.ny, color);
            if (active_timing_) {
              active_timing_->nwt_homo_updates += update_count;
              ++active_timing_->nwt_dispatch_calls;
            }
            executor.parallel_do(workers, [&](std::size_t worker) {
              const auto worker_begin = collect_nwt_detail
                                            ? std::chrono::steady_clock::now()
                                            : std::chrono::steady_clock::time_point{};
              double sumsq = 0.0;
              double max_delta = 0.0;
              std::size_t count = 0;
              std::size_t clipped = 0;

              const auto process_linear_tile =
                  [&](const NwtHalfWorkspace::HomoTile &tile) {
                    const std::size_t end =
                        static_cast<std::size_t>(tile.begin) + tile.count;
                    for (std::size_t h = tile.begin; h < end; ++h) {
                      const auto base = static_cast<std::ptrdiff_t>(h);
                      const double old_phi = write_phi[h];
                      const double neighbor_sum =
                          eps_vec[h] *
                          (read_phi[base + offsets.xp] +
                           read_phi[base + offsets.xm] +
                           read_phi[base + offsets.yp] +
                           read_phi[base + offsets.ym] +
                           read_phi[base + offsets.zp] +
                           read_phi[base + offsets.zm]);
                      const double target =
                          (neighbor_sum + rhs_vec[h]) /
                          std::max(1.0e-12, eps_sum_vec[h]);
                      const double updated =
                          old_phi + omega * (target - old_phi);
                      const double delta = updated - old_phi;
                      write_phi[h] = updated;
                      if (accumulate_rms) {
                        sumsq += delta * delta;
                        max_delta = std::max(max_delta, std::abs(delta));
                        ++count;
                      }
                    }
                  };

              const auto process_salt_tile =
                  [&](const NwtHalfWorkspace::HomoTile &tile) {
                    const std::size_t end =
                        static_cast<std::size_t>(tile.begin) + tile.count;
                    if (screen_scale == 0.0) {
                      for (std::size_t h = tile.begin; h < end; ++h) {
                        const auto base = static_cast<std::ptrdiff_t>(h);
                        const double old_phi = write_phi[h];
                        const double neighbor_sum =
                            eps_vec[h] *
                            (read_phi[base + offsets.xp] +
                             read_phi[base + offsets.xm] +
                             read_phi[base + offsets.yp] +
                             read_phi[base + offsets.ym] +
                             read_phi[base + offsets.zp] +
                             read_phi[base + offsets.zm]);
                        const double target =
                            (neighbor_sum + rhs_vec[h]) /
                            std::max(1.0e-12, eps_sum_vec[h]);
                        const double updated =
                            old_phi + omega * (target - old_phi);
                        const double delta = updated - old_phi;
                        write_phi[h] = updated;
                        if (accumulate_rms) {
                          sumsq += delta * delta;
                          max_delta = std::max(max_delta, std::abs(delta));
                          ++count;
                        }
                      }
                      return;
                    }
                    for (std::size_t h = tile.begin; h < end; ++h) {
                      const auto base = static_cast<std::ptrdiff_t>(h);
                      const double old_phi = write_phi[h];
                      const double neighbor_sum =
                          eps_vec[h] *
                          (read_phi[base + offsets.xp] +
                           read_phi[base + offsets.xm] +
                           read_phi[base + offsets.yp] +
                           read_phi[base + offsets.ym] +
                           read_phi[base + offsets.zp] +
                           read_phi[base + offsets.zm]);
                      const double active_screen =
                          screen_vec[h] * screen_scale;
                      const double raw_arg = old_phi * argument_scale;
                      if (raw_arg > effective_clip ||
                          raw_arg < -effective_clip)
                        ++clipped;
                      const auto sinh_cosh =
                          delphi_sinh_cosh_tail_single_sinh(
                              raw_arg, effective_clip, sinh_clip);
                      const double source = sinh_cosh.first * source_scale;
                      const double diagonal =
                          eps_sum_vec[h] + active_screen * sinh_cosh.second;
                      const double nonlinear_rhs =
                          active_screen *
                          (sinh_cosh.second * old_phi - source);
                      const double target =
                          (neighbor_sum + rhs_vec[h] + nonlinear_rhs) /
                          std::max(1.0e-12, diagonal);
                      const double updated =
                          old_phi + omega * (target - old_phi);
                      const double delta = updated - old_phi;
                      write_phi[h] = updated;
                      if (accumulate_rms) {
                        sumsq += delta * delta;
                        max_delta = std::max(max_delta, std::abs(delta));
                        ++count;
                      }
                    }
                  };

              const std::size_t linear_begin =
                  linear_tiles.size() * worker / workers;
              const std::size_t linear_end =
                  linear_tiles.size() * (worker + 1) / workers;
              const auto linear_begin_time =
                  collect_nwt_detail
                      ? std::chrono::steady_clock::now()
                      : std::chrono::steady_clock::time_point{};
              for (std::size_t tile_index = linear_begin;
                   tile_index < linear_end; ++tile_index) {
#if TITAN_SIMD_X86
                if (salt_tiled_avx2) {
                  update_nwt_homo_linear_span_avx2(
                      write_phi.data(), read_phi.data(), rhs_vec.data(),
                      eps_sum_vec.data(), eps_vec.data(),
                      linear_tiles[tile_index].begin,
                      linear_tiles[tile_index].count, offsets, omega,
                      accumulate_rms, sumsq, max_delta, count);
                } else
#endif
                {
                  process_linear_tile(linear_tiles[tile_index]);
                }
              }
              if (collect_nwt_detail) {
                local_stencil_ms[worker] =
                    elapsed_ms(linear_begin_time,
                               std::chrono::steady_clock::now());
              }

              const std::size_t salt_begin =
                  salt_tiles.size() * worker / workers;
              const std::size_t salt_end =
                  salt_tiles.size() * (worker + 1) / workers;
              const auto salt_begin_time =
                  collect_nwt_detail
                      ? std::chrono::steady_clock::now()
                      : std::chrono::steady_clock::time_point{};
              for (std::size_t tile_index = salt_begin; tile_index < salt_end;
                   ++tile_index) {
                process_salt_tile(salt_tiles[tile_index]);
              }
              if (collect_nwt_detail) {
                local_salt_ms[worker] =
                    elapsed_ms(salt_begin_time,
                               std::chrono::steady_clock::now());
                local_worker_ms[worker] =
                    elapsed_ms(worker_begin, std::chrono::steady_clock::now());
              }

              local_sumsq[worker] += sumsq;
              local_max[worker] = std::max(local_max[worker], max_delta);
              local_count[worker] += count;
              local_clipped[worker] += clipped;
            });
            const double wall_ms =
                elapsed_ms(begin_time, std::chrono::steady_clock::now());
            add_timing_ms(TimingBucket::NwtHomoUpdate, wall_ms);
            if (collect_nwt_detail) {
              double max_worker_ms = 0.0;
              double max_stencil_ms = 0.0;
              double max_salt_ms = 0.0;
              double max_write_ms = 0.0;
              for (std::size_t worker = 0; worker < workers; ++worker) {
                max_worker_ms =
                    std::max(max_worker_ms, local_worker_ms[worker]);
                max_stencil_ms =
                    std::max(max_stencil_ms, local_stencil_ms[worker]);
                max_salt_ms = std::max(max_salt_ms, local_salt_ms[worker]);
                max_write_ms = std::max(max_write_ms, local_write_ms[worker]);
              }
              add_timing_ms(TimingBucket::NwtHomoStencil, max_stencil_ms);
              add_timing_ms(TimingBucket::NwtHomoSaltMath, max_salt_ms);
              add_timing_ms(TimingBucket::NwtHomoWrite, max_write_ms);
              add_timing_ms(TimingBucket::NwtThreadWait,
                            std::max(0.0, wall_ms - max_worker_ms));
            }
          };

      for (int iter = 1; iter <= sweeps; ++iter) {
        std::fill(local_sumsq.begin(), local_sumsq.end(), 0.0);
        std::fill(local_max.begin(), local_max.end(), 0.0);
        std::fill(local_count.begin(), local_count.end(), 0);
        std::fill(local_clipped.begin(), local_clipped.end(), 0);

        update_homo_tiles(nwt_half_workspace_.even_homo_linear_tiles,
                          nwt_half_workspace_.even_homo_salt_tiles, even_phi,
                          odd_phi, rhs_even, eps_sum_even, screen_even,
                          homo_eps_even, 0, false, even_homo_updates);
        run_timed_cells(nwt_half_workspace_.even_hetero_cells, even_phi,
                        odd_phi, false, false,
                        TimingBucket::NwtHeteroUpdate);
        const bool accumulate_odd =
            !delphi_block_stop || (iter == sweeps);
        update_homo_tiles(nwt_half_workspace_.odd_homo_linear_tiles,
                          nwt_half_workspace_.odd_homo_salt_tiles, odd_phi,
                          even_phi, rhs_odd, eps_sum_odd, screen_odd,
                          homo_eps_odd, 1, accumulate_odd, odd_homo_updates);
        run_timed_cells(nwt_half_workspace_.odd_hetero_cells, odd_phi,
                        even_phi, false, accumulate_odd,
                        TimingBucket::NwtHeteroUpdate);

        const auto reduce_time_begin = std::chrono::steady_clock::now();
        double sumsq = 0.0;
        double max_delta = 0.0;
        std::size_t count = 0;
        std::size_t clipped = 0;
        for (std::size_t worker = 0; worker < workers; ++worker) {
          sumsq += local_sumsq[worker];
          max_delta = std::max(max_delta, local_max[worker]);
          count += local_count[worker];
          clipped += local_clipped[worker];
        }
        clipped_points = std::max(clipped_points, clipped);
        stats.iterations = iter;
        if (delphi_block_stop) {
          stats.rms =
              delphi_rms_denominator > 0
                  ? std::sqrt(sumsq /
                              static_cast<double>(delphi_rms_denominator))
                  : 0.0;
        } else {
          stats.rms =
              count > 0 ? std::sqrt(sumsq / static_cast<double>(count)) : 0.0;
        }
        stats.max_dphi = max_delta;
        const double reduce_ms =
            std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - reduce_time_begin)
                .count();
        add_timing_ms(TimingBucket::NwtRmsReduce, reduce_ms);
        add_timing_ms(TimingBucket::NwtReduce, reduce_ms);
        add_timing_ms(TimingBucket::NwtBlockOverhead, reduce_ms);
        if (!delphi_block_stop &&
            (stats.max_dphi <= args_.nl_dphi ||
             stats.rms <= args_.nl_residual)) {
          stats.converged = true;
          break;
        }
      }

      if (delphi_block_stop &&
          (stats.max_dphi < args_.nl_dphi || stats.rms < args_.nl_residual))
        stats.converged = true;
      add_timing_ms(TimingBucket::NwtHalfUpdate,
                    std::chrono::duration<double, std::milli>(
                        std::chrono::steady_clock::now() - update_time_begin)
                        .count());
      return stats;
    }
    if (row_span_requested && nwt_half_workspace_.row_span_ready) {
      const auto span_cell_count =
          [](const std::vector<NwtHalfWorkspace::HomoSpan> &spans) {
            std::uint64_t total = 0;
            for (const auto &span : spans)
              total += span.count;
            return total;
          };
      const std::uint64_t even_homo_updates =
          span_cell_count(nwt_half_workspace_.even_homo_linear_spans) +
          span_cell_count(nwt_half_workspace_.even_homo_salt_spans);
      const std::uint64_t odd_homo_updates =
          span_cell_count(nwt_half_workspace_.odd_homo_linear_spans) +
          span_cell_count(nwt_half_workspace_.odd_homo_salt_spans);

      const auto update_homo_spans =
          [&](const std::vector<NwtHalfWorkspace::HomoSpan> &linear_spans,
              const std::vector<NwtHalfWorkspace::HomoSpan> &salt_spans,
              AlignedDoubleVector &write_phi,
              const AlignedDoubleVector &read_phi,
              const AlignedDoubleVector &rhs_vec,
              const AlignedDoubleVector &eps_sum_vec,
              const AlignedDoubleVector &screen_vec,
              const AlignedDoubleVector &eps_vec, int color,
              bool accumulate_rms, std::uint64_t update_count) {
            if (linear_spans.empty() && salt_spans.empty())
              return;
            const auto begin_time = std::chrono::steady_clock::now();
            const NwtHomoOffsets offsets =
                nwt_homo_offsets(shape.nx, shape.ny, color);
            if (active_timing_) {
              active_timing_->nwt_homo_updates += update_count;
              ++active_timing_->nwt_dispatch_calls;
            }
            executor.parallel_do(workers, [&](std::size_t worker) {
              double sumsq = 0.0;
              double max_delta = 0.0;
              std::size_t count = 0;
              std::size_t clipped = 0;

              const auto process_linear_span_scalar =
                  [&](const NwtHalfWorkspace::HomoSpan &span) {
                    const std::size_t end =
                        static_cast<std::size_t>(span.begin) + span.count;
                    for (std::size_t h = span.begin; h < end; ++h) {
                      const auto base = static_cast<std::ptrdiff_t>(h);
                      const double old_phi = write_phi[h];
                      const double neighbor_sum =
                          eps_vec[h] *
                          (read_phi[base + offsets.xp] +
                           read_phi[base + offsets.xm] +
                           read_phi[base + offsets.yp] +
                           read_phi[base + offsets.ym] +
                           read_phi[base + offsets.zp] +
                           read_phi[base + offsets.zm]);
                      const double target =
                          (neighbor_sum + rhs_vec[h]) /
                          std::max(1.0e-12, eps_sum_vec[h]);
                      const double updated =
                          old_phi + omega * (target - old_phi);
                      const double delta = updated - old_phi;
                      write_phi[h] = updated;
                      if (accumulate_rms) {
                        sumsq += delta * delta;
                        max_delta = std::max(max_delta, std::abs(delta));
                        ++count;
                      }
                    }
                  };

              const auto process_salt_span =
                  [&](const NwtHalfWorkspace::HomoSpan &span) {
                    const std::size_t end =
                        static_cast<std::size_t>(span.begin) + span.count;
                    for (std::size_t h = span.begin; h < end; ++h) {
                      const auto base = static_cast<std::ptrdiff_t>(h);
                      const double old_phi = write_phi[h];
                      const double neighbor_sum =
                          eps_vec[h] *
                          (read_phi[base + offsets.xp] +
                           read_phi[base + offsets.xm] +
                           read_phi[base + offsets.yp] +
                           read_phi[base + offsets.ym] +
                           read_phi[base + offsets.zp] +
                           read_phi[base + offsets.zm]);
                      const double active_screen =
                          screen_vec[h] * screen_scale;
                      double diagonal = eps_sum_vec[h];
                      double nonlinear_rhs = 0.0;
                      if (active_screen > 0.0) {
                        const double raw_arg = old_phi * argument_scale;
                        if (raw_arg > effective_clip ||
                            raw_arg < -effective_clip)
                          ++clipped;
                        const auto sinh_cosh =
                            delphi_sinh_cosh_tail_precomputed(
                                raw_arg, effective_clip, sinh_clip);
                        const double source = sinh_cosh.first * source_scale;
                        diagonal += active_screen * sinh_cosh.second;
                        nonlinear_rhs =
                            active_screen *
                            (sinh_cosh.second * old_phi - source);
                      }
                      const double target =
                          (neighbor_sum + rhs_vec[h] + nonlinear_rhs) /
                          std::max(1.0e-12, diagonal);
                      const double updated =
                          old_phi + omega * (target - old_phi);
                      const double delta = updated - old_phi;
                      write_phi[h] = updated;
                      if (accumulate_rms) {
                        sumsq += delta * delta;
                        max_delta = std::max(max_delta, std::abs(delta));
                        ++count;
                      }
                    }
                  };

              const auto process_linear_range =
                  [&](std::size_t begin, std::size_t end) {
                    for (std::size_t span_index = begin; span_index < end;
                         ++span_index) {
#if TITAN_SIMD_X86
                      if (row_span_avx2) {
                        update_nwt_homo_linear_span_avx2(
                            write_phi.data(), read_phi.data(), rhs_vec.data(),
                            eps_sum_vec.data(), eps_vec.data(),
                            linear_spans[span_index].begin,
                            linear_spans[span_index].count, offsets, omega,
                            accumulate_rms, sumsq, max_delta, count);
                      } else
#endif
                      {
                        process_linear_span_scalar(linear_spans[span_index]);
                      }
                    }
                  };

              const std::size_t linear_begin =
                  linear_spans.size() * worker / workers;
              const std::size_t linear_end =
                  linear_spans.size() * (worker + 1) / workers;
              process_linear_range(linear_begin, linear_end);

              const std::size_t salt_begin =
                  salt_spans.size() * worker / workers;
              const std::size_t salt_end =
                  salt_spans.size() * (worker + 1) / workers;
              for (std::size_t span_index = salt_begin;
                   span_index < salt_end; ++span_index) {
                process_salt_span(salt_spans[span_index]);
              }

              local_sumsq[worker] += sumsq;
              local_max[worker] = std::max(local_max[worker], max_delta);
              local_count[worker] += count;
              local_clipped[worker] += clipped;
            });
            add_timing_ms(TimingBucket::NwtHomoUpdate,
                          std::chrono::duration<double, std::milli>(
                              std::chrono::steady_clock::now() - begin_time)
                              .count());
          };

      for (int iter = 1; iter <= sweeps; ++iter) {
        std::fill(local_sumsq.begin(), local_sumsq.end(), 0.0);
        std::fill(local_max.begin(), local_max.end(), 0.0);
        std::fill(local_count.begin(), local_count.end(), 0);
        std::fill(local_clipped.begin(), local_clipped.end(), 0);

        update_homo_spans(nwt_half_workspace_.even_homo_linear_spans,
                          nwt_half_workspace_.even_homo_salt_spans, even_phi,
                          odd_phi, rhs_even, eps_sum_even, screen_even,
                          homo_eps_even, 0, false, even_homo_updates);
        run_timed_cells(nwt_half_workspace_.even_hetero_cells, even_phi,
                        odd_phi, false, false,
                        TimingBucket::NwtHeteroUpdate);
        const bool accumulate_odd =
            !delphi_block_stop || (iter == sweeps);
        update_homo_spans(nwt_half_workspace_.odd_homo_linear_spans,
                          nwt_half_workspace_.odd_homo_salt_spans, odd_phi,
                          even_phi, rhs_odd, eps_sum_odd, screen_odd,
                          homo_eps_odd, 1, accumulate_odd, odd_homo_updates);
        run_timed_cells(nwt_half_workspace_.odd_hetero_cells, odd_phi,
                        even_phi, false, accumulate_odd,
                        TimingBucket::NwtHeteroUpdate);

        const auto reduce_time_begin = std::chrono::steady_clock::now();
        double sumsq = 0.0;
        double max_delta = 0.0;
        std::size_t count = 0;
        std::size_t clipped = 0;
        for (std::size_t worker = 0; worker < workers; ++worker) {
          sumsq += local_sumsq[worker];
          max_delta = std::max(max_delta, local_max[worker]);
          count += local_count[worker];
          clipped += local_clipped[worker];
        }
        clipped_points = std::max(clipped_points, clipped);
        stats.iterations = iter;
        if (delphi_block_stop) {
          stats.rms =
              delphi_rms_denominator > 0
                  ? std::sqrt(sumsq /
                              static_cast<double>(delphi_rms_denominator))
                  : 0.0;
        } else {
          stats.rms =
              count > 0 ? std::sqrt(sumsq / static_cast<double>(count)) : 0.0;
        }
        stats.max_dphi = max_delta;
        const double reduce_ms =
            std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - reduce_time_begin)
                .count();
        add_timing_ms(TimingBucket::NwtRmsReduce, reduce_ms);
        add_timing_ms(TimingBucket::NwtReduce, reduce_ms);
        if (!delphi_block_stop &&
            (stats.max_dphi <= args_.nl_dphi ||
             stats.rms <= args_.nl_residual)) {
          stats.converged = true;
          break;
        }
      }

      if (delphi_block_stop &&
          (stats.max_dphi < args_.nl_dphi || stats.rms < args_.nl_residual))
        stats.converged = true;
      add_timing_ms(TimingBucket::NwtHalfUpdate,
                    std::chrono::duration<double, std::milli>(
                        std::chrono::steady_clock::now() - update_time_begin)
                        .count());
      return stats;
    }

    for (int iter = 1; iter <= sweeps; ++iter) {
      std::fill(local_sumsq.begin(), local_sumsq.end(), 0.0);
      std::fill(local_max.begin(), local_max.end(), 0.0);
      std::fill(local_count.begin(), local_count.end(), 0);
      std::fill(local_clipped.begin(), local_clipped.end(), 0);

      run_timed_cells(nwt_half_workspace_.even_homo_cells, even_phi, odd_phi,
                      true, false, TimingBucket::NwtHomoUpdate);
      run_timed_cells(nwt_half_workspace_.even_hetero_cells, even_phi, odd_phi,
                      false, false, TimingBucket::NwtHeteroUpdate);
      const bool accumulate_odd =
          !delphi_block_stop || (iter == sweeps);
      run_timed_cells(nwt_half_workspace_.odd_homo_cells, odd_phi, even_phi,
                      true, accumulate_odd, TimingBucket::NwtHomoUpdate);
      run_timed_cells(nwt_half_workspace_.odd_hetero_cells, odd_phi, even_phi,
                      false, accumulate_odd, TimingBucket::NwtHeteroUpdate);

      const auto reduce_time_begin = std::chrono::steady_clock::now();
      double sumsq = 0.0;
      double max_delta = 0.0;
      std::size_t count = 0;
      std::size_t clipped = 0;
      for (std::size_t worker = 0; worker < workers; ++worker) {
        sumsq += local_sumsq[worker];
        max_delta = std::max(max_delta, local_max[worker]);
        count += local_count[worker];
        clipped += local_clipped[worker];
      }
      clipped_points = std::max(clipped_points, clipped);
      stats.iterations = iter;
      if (delphi_block_stop) {
        stats.rms =
            delphi_rms_denominator > 0
                ? std::sqrt(sumsq /
                            static_cast<double>(delphi_rms_denominator))
                : 0.0;
      } else {
        stats.rms =
            count > 0 ? std::sqrt(sumsq / static_cast<double>(count)) : 0.0;
      }
      stats.max_dphi = max_delta;
      const double reduce_ms =
          std::chrono::duration<double, std::milli>(
              std::chrono::steady_clock::now() - reduce_time_begin)
              .count();
      add_timing_ms(TimingBucket::NwtRmsReduce, reduce_ms);
      add_timing_ms(TimingBucket::NwtReduce, reduce_ms);
      if (!delphi_block_stop &&
          (stats.max_dphi <= args_.nl_dphi || stats.rms <= args_.nl_residual)) {
        stats.converged = true;
        break;
      }
    }

    if (delphi_block_stop &&
        (stats.max_dphi < args_.nl_dphi || stats.rms < args_.nl_residual))
      stats.converged = true;
    add_timing_ms(TimingBucket::NwtHalfUpdate,
                  std::chrono::duration<double, std::milli>(
                      std::chrono::steady_clock::now() - update_time_begin)
                      .count());
    return stats;
  }

  for (int iter = 1; iter <= sweeps; ++iter) {
    std::fill(local_sumsq.begin(), local_sumsq.end(), 0.0);
    std::fill(local_max.begin(), local_max.end(), 0.0);
    std::fill(local_count.begin(), local_count.end(), 0);
    std::fill(local_clipped.begin(), local_clipped.end(), 0);

    for (int color = 0; color < 2; ++color) {
      executor.parallel_for(
          interior_z, workers, 1,
          [&](std::size_t worker, std::size_t begin, std::size_t end) {
            double sumsq = 0.0;
            double max_delta = 0.0;
            std::size_t count = 0;
            std::size_t clipped = 0;
            for (std::size_t kk = begin; kk < end; ++kk) {
              const std::size_t k = kk + 1;
              for (std::size_t j = 1; j + 1 < shape.ny; ++j) {
                const std::size_t first_i =
                    1 + (static_cast<int>((1 + j + k) & 1ULL) == color ? 0U
                                                                        : 1U);
                for (std::size_t i = first_i; i + 1 < shape.nx; i += 2) {
                  const std::size_t id = idx3(i, j, k, shape.nx, shape.ny);
                  const std::size_t h = id / 2;
                  double old_phi = 0.0;
                  double neighbor_sum = 0.0;
                  double eps_sum = 1.0;
                  double active_screen = 0.0;
                  double rhs_rho = 0.0;
                  if (color == 0) {
                    old_phi = even_phi[h];
                    neighbor_sum =
                        eps_x_[id] * odd_phi[(id + 1) / 2] +
                        eps_x_[id - 1] * odd_phi[(id - 1) / 2] +
                        eps_y_[id] * odd_phi[(id + shape.nx) / 2] +
                        eps_y_[id - shape.nx] *
                            odd_phi[(id - shape.nx) / 2] +
                        eps_z_[id] * odd_phi[(id + xy) / 2] +
                        eps_z_[id - xy] * odd_phi[(id - xy) / 2];
                    eps_sum = eps_sum_even[h];
                    active_screen = screen_even[h] * screen_scale;
                    rhs_rho = rhs_even[h];
                  } else {
                    old_phi = odd_phi[h];
                    neighbor_sum =
                        eps_x_[id] * even_phi[(id + 1) / 2] +
                        eps_x_[id - 1] * even_phi[(id - 1) / 2] +
                        eps_y_[id] * even_phi[(id + shape.nx) / 2] +
                        eps_y_[id - shape.nx] *
                            even_phi[(id - shape.nx) / 2] +
                        eps_z_[id] * even_phi[(id + xy) / 2] +
                        eps_z_[id - xy] * even_phi[(id - xy) / 2];
                    eps_sum = eps_sum_odd[h];
                    active_screen = screen_odd[h] * screen_scale;
                    rhs_rho = rhs_odd[h];
                  }
                  double diagonal = eps_sum;
                  double nonlinear_rhs = 0.0;
                  if (active_screen > 0.0) {
                    const double raw_arg = nonlinear_argument_value(old_phi,
                                                                    alpha);
                    if (raw_arg > clip || raw_arg < -clip)
                      ++clipped;
                    const auto sinh_cosh = delphi_sinh_cosh_tail(raw_arg, clip);
                    const double source =
                        (delphi_field_units_active() ||
                         args_.test_nl_source == "delphi-dimensionless")
                            ? sinh_cosh.first
                            : sinh_cosh.first / std::max(1.0e-30, alpha);
                    diagonal += active_screen * sinh_cosh.second;
                    nonlinear_rhs =
                        active_screen *
                        (sinh_cosh.second * old_phi - source);
                  }
                  const double target =
                      (neighbor_sum + rhs_rho + nonlinear_rhs) /
                      std::max(1.0e-12, diagonal);
                  const double updated = old_phi + omega * (target - old_phi);
                  const double delta = updated - old_phi;
                  if (color == 0)
                    even_phi[h] = updated;
                  else
                    odd_phi[h] = updated;
                  if (!delphi_block_stop || (iter == sweeps && color == 1)) {
                    sumsq += delta * delta;
                    max_delta = std::max(max_delta, std::abs(delta));
                    ++count;
                  }
                }
              }
            }
            local_sumsq[worker] += sumsq;
            local_max[worker] = std::max(local_max[worker], max_delta);
            local_count[worker] += count;
            local_clipped[worker] += clipped;
          });
    }

    const auto reduce_time_begin = std::chrono::steady_clock::now();
    double sumsq = 0.0;
    double max_delta = 0.0;
    std::size_t count = 0;
    std::size_t clipped = 0;
    for (std::size_t worker = 0; worker < workers; ++worker) {
      sumsq += local_sumsq[worker];
      max_delta = std::max(max_delta, local_max[worker]);
      count += local_count[worker];
      clipped += local_clipped[worker];
    }
    clipped_points = std::max(clipped_points, clipped);
    stats.iterations = iter;
    if (delphi_block_stop) {
      stats.rms = delphi_rms_denominator > 0
                      ? std::sqrt(sumsq /
                                  static_cast<double>(delphi_rms_denominator))
                      : 0.0;
    } else {
      stats.rms =
          count > 0 ? std::sqrt(sumsq / static_cast<double>(count)) : 0.0;
    }
    stats.max_dphi = max_delta;
    add_timing_ms(TimingBucket::NwtRmsReduce,
                  std::chrono::duration<double, std::milli>(
                      std::chrono::steady_clock::now() - reduce_time_begin)
                      .count());
    if (!delphi_block_stop &&
        (stats.max_dphi <= args_.nl_dphi || stats.rms <= args_.nl_residual)) {
      stats.converged = true;
      break;
    }
  }

  if (delphi_block_stop &&
      (stats.max_dphi < args_.nl_dphi || stats.rms < args_.nl_residual))
    stats.converged = true;
  add_timing_ms(TimingBucket::NwtHalfUpdate,
                std::chrono::duration<double, std::milli>(
                    std::chrono::steady_clock::now() - update_time_begin)
                    .count());
  return stats;
}

AnalyzerImpl::FieldSolve AnalyzerImpl::solve_newton_nonlinear_field(
    const std::vector<Atom> &atoms, double eps_in, double eps_out, double kappa,
    TitanExecutor &executor, const FrameResult *map_result) {
  FieldSolve field =
      solve_linear_field(atoms, eps_in, eps_out, kappa, executor, nullptr);
  field.linear_seed_solver = field.solver_mode;
  field.pbe_model_effective = "nonlinear";
  field.nl_solver = "newton";
  if (args_.nl_inner == "auto") {
    const bool delphi_nwt_fixture =
        !args_.test_pqr.empty() && args_.surface == "delphi-vdw";
    field.nl_inner_solver =
        delphi_nwt_fixture
            ? std::string("sor")
            : (field.shape.n > 250000 ? std::string("vcycle")
                                       : std::string("sor"));
  } else {
    field.nl_inner_solver = args_.nl_inner;
  }
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
  const auto nonlinear_time_begin = std::chrono::steady_clock::now();
  if (delphi_nwt_parity_active()) {
    for (std::size_t k = 1; k + 1 < field.shape.nz; ++k) {
      for (std::size_t j = 1; j + 1 < field.shape.ny; ++j) {
        for (std::size_t i = 1; i + 1 < field.shape.nx; ++i) {
          phi_[idx3(i, j, k, field.shape.nx, field.shape.ny)] = 0.0;
        }
      }
    }
    field.linear_seed_solver = "delphi_boundary_zero";
  }

  const double alpha = nonlinear_alpha();
  const double clip = effective_nonlinear_clip();
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
  if (delphi_nwt_parity_active())
    screen_scales = {1.0};
  else if (args_.nl_continuation == "auto" && args_.salt > 0.5)
    screen_scales = {0.1, 0.25, 0.5, 0.75, 1.0};
  else if (args_.nl_continuation == "auto")
    screen_scales = {0.25, 0.5, 0.75, 1.0};
  else
    screen_scales = {1.0};

  AlignedDoubleVector nonlinear_before;
  AlignedDoubleVector nwt_before;
  std::ofstream nwt_progress;
  const bool nwt_progress_strict =
      delphi_nwt_parity_active() &&
      (args_.nwt_progress == "strict" ||
       (args_.nwt_progress == "auto" && args_.parity_dump != "off"));
  if (nwt_progress_strict) {
    nwt_progress.open("pbe_nwt_solver_progress.csv", std::ios::out | std::ios::trunc);
    if (!nwt_progress.is_open())
      throw std::runtime_error("PBE: failed to open pbe_nwt_solver_progress.csv.");
    nwt_progress
        << "block,total_iterations,screen_scale,sweeps,rms,max_dphi,"
           "operator_residual,operator_residual_max,clipped_points,"
           "stop_reason,field_units,operator_convention\n";
    nwt_progress << std::setprecision(12);
  }

  for (double screen_scale : screen_scales) {
    final_screen_scale = screen_scale;
    ++continuation_steps;
    std::size_t clipped = 0;
    const bool fast_delphi_nwt =
        delphi_nwt_parity_active() && field.nl_inner_solver == "sor" &&
        !nwt_progress_strict;
    if (fast_delphi_nwt) {
      residual = std::numeric_limits<double>::infinity();
      last_residual_after = residual;
    } else {
      const auto source_time_begin = std::chrono::steady_clock::now();
      clipped = update_nonlinear_source(field.shape, source, alpha, clip);
      add_timing_ms(TimingBucket::NwtSourceUpdate,
                    std::chrono::duration<double, std::milli>(
                        std::chrono::steady_clock::now() - source_time_begin)
                        .count());
      max_clipped = std::max(max_clipped, clipped);
      const auto residual_time_begin = std::chrono::steady_clock::now();
      residual = compute_nonlinear_residual_vector(
          field.shape, phi_, source, screen_scale, residual_vec, executor);
      add_timing_ms(TimingBucket::NwtOperatorResidual,
                    std::chrono::duration<double, std::milli>(
                        std::chrono::steady_clock::now() - residual_time_begin)
                        .count());
      last_residual_after = residual;
    }
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
    if (delphi_nwt_parity_active() && field.nl_inner_solver == "sor") {
      const int block_size = std::max(1, args_.nl_block);
      while (total_iterations < args_.nl_max_iter && !step_done) {
        last_residual_before = residual;
        std::size_t nwt_clipped = 0;
        const int sweeps =
            std::max(1, std::min(block_size, args_.nl_max_iter - total_iterations));
        SolveStats nwt_stats = solve_newton_nwt_sor_block(
            field.shape, sweeps, screen_scale, executor, nwt_clipped);
        std::size_t nwt_source_clipped = 0;
        double operator_residual = nwt_stats.rms;
        if (nwt_progress_strict) {
          const auto source_time_begin = std::chrono::steady_clock::now();
          nwt_source_clipped =
              update_nonlinear_source(field.shape, source, alpha, clip);
          add_timing_ms(TimingBucket::NwtSourceUpdate,
                        std::chrono::duration<double, std::milli>(
                            std::chrono::steady_clock::now() -
                            source_time_begin)
                            .count());
          const auto residual_time_begin = std::chrono::steady_clock::now();
          operator_residual = compute_nonlinear_residual_vector(
              field.shape, phi_, source, screen_scale, residual_vec, executor);
          add_timing_ms(TimingBucket::NwtOperatorResidual,
                        std::chrono::duration<double, std::milli>(
                            std::chrono::steady_clock::now() -
                            residual_time_begin)
                            .count());
        }
        double operator_residual_max = 0.0;
        if (nwt_progress_strict) {
          for (double value : residual_vec)
            operator_residual_max =
                std::max(operator_residual_max, std::abs(value));
        } else {
          operator_residual_max = nwt_stats.max_dphi;
        }
        finite = std::isfinite(operator_residual) &&
                 std::isfinite(nwt_stats.rms) &&
                 std::isfinite(nwt_stats.max_dphi) &&
                 std::all_of(phi_.begin(), phi_.end(),
                             [](double v) { return std::isfinite(v); });
        const int used_inner = std::max(1, nwt_stats.iterations);
        inner_iterations += used_inner;
        total_iterations += used_inner;
        ++blocks;
        residual = std::isfinite(operator_residual) ? operator_residual : residual;
        last_residual_after = residual;
        last_step_norm = nwt_stats.rms;
        max_step = nwt_stats.max_dphi;
        max_clipped =
            std::max(max_clipped, std::max(nwt_clipped, nwt_source_clipped));
        last_stats = nwt_stats;
        last_inner_effective = "nwt_sor_block";
        last_inner_levels = 1;
        last_inner_cycles = used_inner;
        last_inner_residual_before = last_residual_before;
        last_inner_residual_after = residual;
        last_correction_residual = residual;
        last_correction_max_residual = operator_residual_max;
        line_search_status = "nwt_block_stop";
        last_accepted_lambda = 1.0;
        if (!finite) {
          stop_reason = "nonfinite";
          step_done = true;
        } else if (nwt_stats.max_dphi < args_.nl_dphi ||
            nwt_stats.rms < args_.nl_residual) {
          last_stats.converged = true;
          converged = true;
          stop_reason = nwt_stats.max_dphi < args_.nl_dphi ? "dphi_tol"
                                                           : "rms_tol";
          step_done = true;
        }
        if (nwt_progress.is_open()) {
          const auto progress_io_begin = std::chrono::steady_clock::now();
          nwt_progress << blocks << ',' << total_iterations << ','
                       << screen_scale << ',' << used_inner << ','
                       << nwt_stats.rms << ',' << nwt_stats.max_dphi << ','
                       << operator_residual << ',' << operator_residual_max
                       << ',' << max_clipped << ',' << csv_escape(stop_reason)
                       << ','
                       << csv_escape(delphi_field_units_active()
                                         ? "delphi-dimensionless"
                                         : args_.test_field_units)
                       << ',' << csv_escape("delphi_nwt") << '\n';
          add_timing_ms(TimingBucket::NwtProgressIo,
                        std::chrono::duration<double, std::milli>(
                            std::chrono::steady_clock::now() -
                            progress_io_begin)
                            .count());
        }
        if (step_done)
          break;
      }
      if (!converged && finite && total_iterations >= args_.nl_max_iter)
        stop_reason = "work_budget";
      break;
    }
    while (total_iterations < args_.nl_max_iter && !step_done) {
      nonlinear_before = phi_;
      last_residual_before = residual;
      if (screen_scale >= 1.0 && residual <= args_.nl_residual) {
        converged = true;
        last_stats.converged = true;
        stop_reason = "residual_tol";
        step_done = true;
        break;
      }

      if (field.nl_inner_solver == "sor") {
        nwt_before = phi_;
        std::size_t nwt_clipped = 0;
        SolveStats nwt_stats = solve_newton_nwt_sor_block(
            field.shape,
            std::max(1, std::min(args_.nl_inner_iter,
                                 args_.nl_max_iter - total_iterations)),
            screen_scale, executor, nwt_clipped);
        const std::size_t nwt_source_clipped =
            update_nonlinear_source(field.shape, source, alpha, clip);
        const double nwt_residual = compute_nonlinear_residual_vector(
            field.shape, phi_, source, screen_scale, residual_vec, executor);
        double nwt_residual_max = 0.0;
        for (double value : residual_vec)
          nwt_residual_max = std::max(nwt_residual_max, std::abs(value));
        const bool nwt_finite =
            std::isfinite(nwt_residual) &&
            std::all_of(phi_.begin(), phi_.end(),
                        [](double v) { return std::isfinite(v); });
        if (nwt_finite &&
            (nwt_residual < residual ||
             nwt_residual <= stage_residual_target ||
             nwt_stats.max_dphi <= args_.nl_dphi)) {
          const int used_inner = std::max(1, nwt_stats.iterations);
          inner_iterations += used_inner;
          total_iterations += used_inner;
          ++blocks;
          residual = nwt_residual;
          last_residual_after = residual;
          last_step_norm = nwt_stats.rms;
          max_step = std::max(max_step, nwt_stats.max_dphi);
          max_clipped =
              std::max(max_clipped, std::max(nwt_clipped, nwt_source_clipped));
          last_stats = nwt_stats;
          last_stats.rms = residual;
          last_inner_effective = "nwt_sor";
          last_inner_levels = 1;
          last_inner_cycles = used_inner;
          last_inner_residual_before = last_residual_before;
          last_inner_residual_after = residual;
          last_correction_residual = residual;
          last_correction_max_residual = nwt_residual_max;
          line_search_status = "nwt_sor";
          last_accepted_lambda = 1.0;
          if (residual < best_stage_residual * 0.999) {
            best_stage_residual = residual;
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
                      nwt_stats.max_dphi <= args_.nl_dphi)) {
            stop_reason = "stage_target";
            step_done = true;
          }
          continue;
        }
        phi_ = nwt_before;
        update_nonlinear_source(field.shape, source, alpha, clip);
        compute_nonlinear_residual_vector(field.shape, phi_, source,
                                          screen_scale, residual_vec,
                                          executor);
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
        phi_ = nonlinear_before;
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
              phi_[id] = nonlinear_before[id] + step;
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
        phi_ = nonlinear_before;
        update_nonlinear_source(field.shape, source, alpha, clip);
        std::size_t fallback_clipped = 0;
        SolveStats fallback_stats = solve_newton_nwt_sor_block(
            field.shape, std::max(1, std::min(args_.nl_inner_iter, 50)),
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
          line_search_status = "fallback_nwt_sor";
          last_inner_effective = "nwt_sor_fallback";
          last_inner_cycles = fallback_stats.iterations;
          inner_iterations += fallback_stats.iterations;
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
        phi_ = nonlinear_before;
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

  if (nwt_half_shape_matches(field.shape)) {
    flush_nwt_half_workspace(field.shape);
    const auto source_time_begin = std::chrono::steady_clock::now();
    const std::size_t final_source_clipped =
        update_nonlinear_source(field.shape, source, alpha, clip);
    add_timing_ms(TimingBucket::NwtSourceUpdate,
                  std::chrono::duration<double, std::milli>(
                      std::chrono::steady_clock::now() - source_time_begin)
                      .count());
    max_clipped = std::max(max_clipped, final_source_clipped);
    const auto residual_time_begin = std::chrono::steady_clock::now();
    residual = compute_nonlinear_residual_vector(
        field.shape, phi_, source, final_screen_scale, residual_vec, executor);
    add_timing_ms(TimingBucket::NwtOperatorResidual,
                  std::chrono::duration<double, std::milli>(
                      std::chrono::steady_clock::now() - residual_time_begin)
                      .count());
    last_residual_after = residual;
  }

  const double reported_nl_residual =
      delphi_nwt_parity_active() ? last_stats.rms : residual;
  field.nl_iterations = total_iterations;
  field.nl_blocks = blocks;
  field.nl_max_dphi = max_step;
  field.nl_residual =
      std::isfinite(reported_nl_residual) ? reported_nl_residual : 0.0;
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
  add_timing_ms(TimingBucket::Solver,
                elapsed_ms(nonlinear_time_begin,
                           std::chrono::steady_clock::now()));
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

AnalyzerImpl::FieldSolve AnalyzerImpl::solve_nonlinear_field(
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

  const auto nonlinear_time_begin = std::chrono::steady_clock::now();
  const double alpha = nonlinear_alpha();
  const double clip = effective_nonlinear_clip();
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

  AlignedDoubleVector picard_before;
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
      if (args_.nl_line_search == "auto")
        picard_before = phi_;
      else
        picard_before.clear();
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
          picard_before.size() == phi_.size()) {
        phi_ = picard_before;
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
  add_timing_ms(TimingBucket::Solver,
                elapsed_ms(nonlinear_time_begin,
                           std::chrono::steady_clock::now()));

  if (map_result != nullptr)
    maybe_write_maps(field.shape, *map_result, true);
  return field;
}

AnalyzerImpl::FieldSolve AnalyzerImpl::solve_physical_field(
    const std::vector<Atom> &atoms, double eps_in, double eps_out, double kappa,
    TitanExecutor &executor, const FrameResult *map_result) {
  if (args_.model == "nonlinear")
    return solve_nonlinear_field(atoms, eps_in, eps_out, kappa, executor,
                                 map_result);
  return solve_linear_field(atoms, eps_in, eps_out, kappa, executor,
                            map_result);
}

std::vector<double>
AnalyzerImpl::compute_sasa_per_atom(const std::vector<Atom> &atoms,
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

double AnalyzerImpl::compute_sasa_total(const std::vector<Atom> &atoms,
                                    TitanExecutor &executor) const {
  const std::vector<double> area_by_atom =
      compute_sasa_per_atom(atoms, executor);
  return TitanSASA::sum_sasa(
      std::span<const double>(area_by_atom.data(), area_by_atom.size()));
}

void AnalyzerImpl::write_energy_header_if_needed() {
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
           "face_refine,pbe_backend,pbe_model,pbe_model_effective,linear_seed_solver,"
           "nl_iterations,nl_blocks,nl_relax,nl_max_dphi,"
           "nl_clipped_points,nl_status,boundary_mode,"
           "test_boundary_convention,test_field_units,charge_assignment,"
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
           "vdwms_midpoint_rms,iscm_water_minus_vacuum_kT,"
           "iscm_water_minus_vacuum_kJmol,"
           "iscm_water_minus_vacuum_kcalmol,iscm_finalization_mode,"
           "iscm_position_source,iscm_position_rms_vs_oracle,"
           "iscm_position_oracle_matches,iscm_native_independent,"
           "nonlinear_stress_kT,nonlinear_stress_kJmol,"
           "nonlinear_stress_kcalmol,nonlinear_osmotic_kT,"
           "nonlinear_osmotic_kJmol,nonlinear_osmotic_kcalmol,"
           "nonlinear_energy_correction_kT,"
           "nonlinear_energy_correction_kJmol,"
           "nonlinear_energy_correction_kcalmol,"
           "nonlinear_energy_status,rpbe_coulomb_norm,"
           "rpbe_grad_coulomb_norm,rpbe_grad_epsilon_norm,"
           "rpbe_source_rms,rpbe_source_max,rpbe_source_status\n";
    energy_header_written_ = true;
  }
}

void AnalyzerImpl::append_energy_row(const FrameResult &result) {
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
              << result.pbe_backend << ',' << result.pbe_model << ','
              << result.pbe_model_effective << ',' << result.linear_seed_solver
              << ',' << result.nl_iterations
              << ',' << result.nl_blocks << ',' << result.nl_relax << ','
              << result.nl_max_dphi << ',' << result.nl_clipped_points << ','
              << result.nl_status << ',' << result.boundary_mode << ','
              << result.test_boundary_convention << ','
              << result.test_field_units << ',' << result.charge_assignment
              << ',' << result.rhs_scale << ','
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
              << result.vdwms_midpoint_rms << ','
              << result.iscm_water_minus_vacuum_kt << ','
              << result.iscm_water_minus_vacuum_kj_mol << ','
              << result.iscm_water_minus_vacuum_kcal_mol << ','
              << csv_escape(result.iscm_finalization_mode) << ','
              << csv_escape(result.iscm_position_source) << ','
              << result.iscm_position_rms_vs_oracle << ','
              << result.iscm_position_oracle_matches << ','
              << (result.iscm_native_independent ? "true" : "false") << ','
              << result.nonlinear_stress_kt << ','
              << result.nonlinear_stress_kj_mol << ','
              << result.nonlinear_stress_kcal_mol << ','
              << result.nonlinear_osmotic_kt << ','
              << result.nonlinear_osmotic_kj_mol << ','
              << result.nonlinear_osmotic_kcal_mol << ','
              << result.nonlinear_energy_correction_kt << ','
              << result.nonlinear_energy_correction_kj_mol << ','
              << result.nonlinear_energy_correction_kcal_mol << ','
              << csv_escape(result.nonlinear_energy_status) << ','
              << result.rpbe_coulomb_norm << ','
              << result.rpbe_grad_coulomb_norm << ','
              << result.rpbe_grad_epsilon_norm << ','
              << result.rpbe_source_rms << ',' << result.rpbe_source_max
              << ',' << csv_escape(result.rpbe_source_status) << '\n';
}

bool AnalyzerImpl::pbe_timing_enabled() const noexcept {
  return args_.timing == "summary" || args_.timing == "detail";
}

void AnalyzerImpl::add_timing_ms(TimingBucket bucket, double ms) const noexcept {
  if (active_timing_ == nullptr || ms <= 0.0 || !std::isfinite(ms))
    return;
  switch (bucket) {
  case TimingBucket::Setup:
    active_timing_->setup_ms += ms;
    break;
  case TimingBucket::Surface:
    active_timing_->surface_ms += ms;
    break;
  case TimingBucket::Epsilon:
    active_timing_->epsilon_ms += ms;
    break;
  case TimingBucket::Boundary:
    active_timing_->boundary_ms += ms;
    break;
  case TimingBucket::Solver:
    active_timing_->solver_ms += ms;
    break;
  case TimingBucket::Energy:
    active_timing_->energy_ms += ms;
    break;
  case TimingBucket::MapIo:
    active_timing_->map_io_ms += ms;
    break;
  case TimingBucket::RpbSource:
    active_timing_->rpb_source_ms += ms;
    break;
  case TimingBucket::FocusSetup:
    active_timing_->focus_setup_ms += ms;
    break;
  case TimingBucket::Probe:
    active_timing_->probe_ms += ms;
    break;
  case TimingBucket::FocusParentSolve:
    active_timing_->focus_parent_solve_ms += ms;
    break;
  case TimingBucket::FocusChildSetup:
    active_timing_->focus_child_setup_ms += ms;
    break;
  case TimingBucket::FocusBoundaryInterp:
    active_timing_->focus_boundary_interp_ms += ms;
    break;
  case TimingBucket::FocusChildSurface:
    active_timing_->focus_child_surface_ms += ms;
    break;
  case TimingBucket::FocusChildEpsilon:
    active_timing_->focus_child_epsilon_ms += ms;
    break;
  case TimingBucket::FocusChildSolver:
    active_timing_->focus_child_solver_ms += ms;
    break;
  case TimingBucket::ProbeEnergy:
    active_timing_->probe_energy_ms += ms;
    break;
  case TimingBucket::AuditOutput:
    active_timing_->audit_output_ms += ms;
    break;
  case TimingBucket::NwtHalfUpdate:
    active_timing_->nwt_half_update_ms += ms;
    break;
  case TimingBucket::NwtRmsReduce:
    active_timing_->nwt_rms_reduce_ms += ms;
    break;
  case TimingBucket::NwtOperatorResidual:
    active_timing_->nwt_operator_residual_ms += ms;
    break;
  case TimingBucket::NwtSourceUpdate:
    active_timing_->nwt_source_update_ms += ms;
    break;
  case TimingBucket::NwtProgressIo:
    active_timing_->nwt_progress_io_ms += ms;
    break;
  case TimingBucket::NwtCoeffPack:
    active_timing_->nwt_coeff_pack_ms += ms;
    break;
  case TimingBucket::NwtHomoUpdate:
    active_timing_->nwt_homo_update_ms += ms;
    break;
  case TimingBucket::NwtHeteroUpdate:
    active_timing_->nwt_hetero_update_ms += ms;
    break;
  case TimingBucket::NwtNonlinTerm:
    active_timing_->nwt_nonlin_term_ms += ms;
    break;
  case TimingBucket::NwtReduce:
    active_timing_->nwt_reduce_ms += ms;
    break;
  case TimingBucket::NwtHomoStencil:
    active_timing_->nwt_homo_stencil_ms += ms;
    break;
  case TimingBucket::NwtHomoSaltMath:
    active_timing_->nwt_homo_salt_math_ms += ms;
    break;
  case TimingBucket::NwtHomoWrite:
    active_timing_->nwt_homo_write_ms += ms;
    break;
  case TimingBucket::NwtThreadWait:
    active_timing_->nwt_thread_wait_ms += ms;
    break;
  case TimingBucket::NwtBlockOverhead:
    active_timing_->nwt_block_overhead_ms += ms;
    break;
  case TimingBucket::LinearHalfUpdate:
    active_timing_->linear_half_update_ms += ms;
    break;
  case TimingBucket::EpsilonCenterDensity:
    active_timing_->epsilon_center_density_ms += ms;
    break;
  case TimingBucket::EpsilonMidpointDensity:
    active_timing_->epsilon_midpoint_density_ms += ms;
    break;
  case TimingBucket::EpsilonIonExclusion:
    active_timing_->epsilon_ion_exclusion_ms += ms;
    break;
  case TimingBucket::EpsilonFinalize:
    active_timing_->epsilon_finalize_ms += ms;
    break;
  }
}

void AnalyzerImpl::write_timing_header_if_needed() {
  if (!pbe_timing_enabled())
    return;
  if (!timing_out_.is_open()) {
    timing_out_.open("pbe_timing.csv", std::ios::out | std::ios::trunc);
    if (!timing_out_.is_open())
      throw std::runtime_error("PBE: failed to open pbe_timing.csv.");
  }
  if (!timing_header_written_) {
    timing_out_
        << "trajectory_id,frame_index,global_index,time_ps,fixture,mode,"
           "pbe_backend,atom_count,focus_child_atom_count,num_threads,scope,"
           "tabi_surface_build_seconds,tabi_surface_panels,"
           "tabi_surface_point_count_requested,"
           "tabi_panel_density_per_A2,tabi_panel_budget_status,"
           "tabi_surface_exposed_points,tabi_surface_buried_points,"
           "tabi_tree_build_seconds,tabi_tree_nodes,tabi_tree_leaf_count,"
           "tabi_near_interactions,tabi_far_interactions,tabi_matvec_mode,"
           "tabi_status,tabi_convergence_reason,tabi_reaction_energy_source,"
           "tabi_gmres_restart_requested,tabi_gmres_restart_effective,"
           "tabi_best_residual,tabi_best_residual_iter,tabi_final_residual,"
           "tabi_best_preconditioned_residual,"
           "tabi_final_preconditioned_residual,"
           "tabi_true_residual_abs,tabi_true_residual_rel,"
           "tabi_true_residual_probe_mode,"
           "tabi_true_residual_probe_scope,"
           "tabi_true_residual_probe_interval,"
           "tabi_true_residual_sample_count,"
           "tabi_true_residual_probe_calls,"
           "tabi_true_residual_probe_seconds,"
           "tabi_sampled_true_residual_rel,"
           "tabi_sampled_true_residual_rows,"
           "tabi_sampled_probe_calls,"
           "tabi_sampled_probe_seconds,"
           "tabi_full_true_residual_confirm_calls,"
           "tabi_full_probe_seconds,"
           "tabi_rhs_seconds,tabi_gmres_seconds,"
           "tabi_matvec_seconds_total,tabi_matvec_seconds_avg,"
           "tabi_matvec_calls,tabi_tree_vs_direct_matvec_rel,"
           "tabi_matvec_audit_rows,tabi_matvec_audit_seconds,"
           "tabi_orthogonalization_seconds,"
           "tabi_preconditioner_seconds,tabi_surface_neighbor_seconds,"
           "tabi_surface_mesh_seconds,tabi_solid_angle_seconds,"
           "tabi_tree_moment_seconds,tabi_tree_traversal_seconds,"
           "tabi_matvec_cost_profile,tabi_leaf_panels_per_second,"
           "tabi_clean_leaf_panels_per_second,tabi_tree_traversal_fraction,"
           "tabi_near_eval_seconds,tabi_far_eval_seconds,"
           "tabi_singular_near_eval_seconds,tabi_leaf_near_eval_seconds,"
           "tabi_leaf_panel_eval_seconds,"
           "tabi_leaf_tile_eval_seconds,tabi_leaf_residual_eval_seconds,"
           "tabi_smooth_leaf_node_eval_seconds,"
           "tabi_leaf_panel_kernel,tabi_leaf_tile_status,"
           "tabi_leaf_panel_pairs_per_second,"
           "tabi_singular_near_interactions,tabi_leaf_near_interactions,"
           "tabi_far_node_interactions,"
           "tabi_leaf_panel_interactions,tabi_leaf_tile_interactions,"
           "tabi_leaf_tile_panel_interactions,"
           "tabi_leaf_residual_panel_interactions,"
           "tabi_leaf_work_profile,tabi_leaf_work_mac_terminal_panels,"
           "tabi_leaf_work_physical_gap_panels,"
           "tabi_leaf_work_fragmentation_panels,"
           "tabi_leaf_work_convergence_panels,"
           "tabi_leaf_work_convergence_estimator_panels,"
           "tabi_leaf_work_convergence_rho_panels,"
           "tabi_leaf_work_convergence_aspect_panels,"
           "tabi_leaf_work_surface_margin_panels,"
           "tabi_leaf_work_surface_margin_hard_panels,"
           "tabi_leaf_work_surface_margin_soft_panels,"
           "tabi_leaf_work_k3_rho_cap_panels,"
           "tabi_leaf_work_singular_fallback_panels,"
           "tabi_leaf_work_scalar_tail_panels,"
           "tabi_leaf_work_profile_accounted_panels,"
           "tabi_leaf_work_mac_terminal_frac,"
           "tabi_leaf_work_convergence_frac,"
           "tabi_leaf_work_surface_margin_frac,"
           "tabi_leaf_work_singular_fallback_frac,"
           "tabi_leaf_work_scalar_tail_frac,"
           "tabi_smooth_leaf_node_interactions,"
           "tabi_leaf_pruned_interactions,"
           "tabi_leaf_mac_candidates,"
           "tabi_leaf_mac_accepted,"
           "tabi_leaf_mac_rejected,"
           "tabi_leaf_mac_far_panels,"
           "tabi_leaf_mac_direct_panels_saved_estimate,"
           "tabi_adaptive_smooth_leaf_accepted,"
           "tabi_adaptive_smooth_leaf_rejected,"
           "tabi_adaptive_span_tile_accepted,"
           "tabi_adaptive_span_tile_rejected,"
           "tabi_adaptive_exact_leaf_fallbacks,"
           "tabi_adaptive_reject_one_ring_overlap,"
           "tabi_adaptive_reject_error_bound_failed,"
           "tabi_adaptive_reject_p4_check_failed,"
           "tabi_adaptive_reject_too_close,"
           "tabi_leaf_work_reduction_ratio,"
           "tabi_target_cluster_interactions,"
           "tabi_target_cluster_targets,"
           "tabi_target_cluster_panels_pruned,"
           "tabi_target_cluster_rejected,"
           "tabi_target_cluster_eval_seconds,"
           "tabi_leaf_pair_batch_interactions,"
           "tabi_leaf_pair_batch_targets,"
           "tabi_leaf_pair_batch_source_panels,"
           "tabi_leaf_pair_batch_residual_panels,"
           "tabi_leaf_pair_batch_fallbacks,"
           "tabi_leaf_pair_batch_seconds,"
           "tabi_leaf_pair_batch_kernel,"
           "tabi_leaf_pair_batch_vectorized_sources,"
           "tabi_leaf_pair_batch_vector_lanes,"
           "tabi_leaf_pair_batch_pairs_per_second,"
           "tabi_leaf_kernel_profile,"
           "tabi_jcp_direct_kernel,"
           "tabi_leaf_kernel_scalar_tasks,"
           "tabi_leaf_kernel_batch_tasks,"
           "tabi_leaf_kernel_source_panel_evals,"
           "tabi_leaf_kernel_scalar_seconds,"
           "tabi_leaf_kernel_target_avx2_seconds,"
           "tabi_leaf_kernel_tail_seconds,"
           "tabi_leaf_pair_batch_soa_sources,"
           "tabi_leaf_pair_batch_soa_seconds,"
           "tabi_leaf_pair_batch_scalar_tail_seconds,"
           "tabi_leaf_pair_batch_clean_source_panels_per_second,"
           "tabi_leaf_pair_selection_profile,"
           "tabi_leaf_pair_probe_cache,"
           "tabi_leaf_pair_target_group,"
           "tabi_leaf_pair_max_lookahead,"
           "tabi_leaf_pair_selection_candidate_scans,"
           "tabi_leaf_pair_selection_clean_rejects,"
           "tabi_leaf_pair_selection_probe_calls,"
           "tabi_leaf_pair_selection_probe_node_visits,"
           "tabi_leaf_pair_selection_probe_seconds,"
           "tabi_leaf_pair_selection_seconds,"
           "tabi_leaf_pair_selection_probe_cache_hits,"
           "tabi_leaf_pair_selection_probe_cache_misses,"
           "tabi_leaf_pair_selection_insufficient_targets,"
           "tabi_leaf_pair_selection_same_target_leaf_rejects,"
           "tabi_leaf_pair_directness_probe,"
           "tabi_leaf_pair_directness_cache,"
           "tabi_leaf_pair_directness_verify,"
           "tabi_leaf_pair_directness_ancestor_probe_calls,"
           "tabi_leaf_pair_directness_ancestor_probe_node_visits,"
           "tabi_leaf_pair_directness_ancestor_probe_seconds,"
           "tabi_leaf_pair_directness_pair_cache_hits,"
           "tabi_leaf_pair_directness_pair_cache_misses,"
           "tabi_leaf_pair_directness_all_direct_proofs,"
           "tabi_leaf_pair_directness_mixed_fallbacks,"
           "tabi_leaf_pair_directness_verify_checks,"
           "tabi_leaf_pair_directness_verify_mismatches,"
           "tabi_leaf_pair_scheduler,"
           "tabi_leaf_pair_worklist_max_active,"
           "tabi_leaf_pair_worklist_verify,"
           "tabi_leaf_pair_worklist_tasks,"
           "tabi_leaf_pair_worklist_groups,"
           "tabi_leaf_pair_worklist_batched_groups,"
           "tabi_leaf_pair_worklist_scalar_flush_tasks,"
           "tabi_leaf_pair_worklist_forced_flushes,"
           "tabi_leaf_pair_worklist_active_groups_max,"
           "tabi_leaf_pair_worklist_seconds,"
           "tabi_leaf_pair_worklist_verify_checks,"
           "tabi_leaf_pair_worklist_verify_mismatches,"
           "tabi_leaf_block_verify,"
           "tabi_leaf_block_target_leaves,"
           "tabi_leaf_block_targets,"
           "tabi_leaf_block_source_groups,"
           "tabi_leaf_block_batched_groups,"
           "tabi_leaf_block_scalar_tail_tasks,"
           "tabi_leaf_block_group_seconds,"
           "tabi_leaf_block_flush_seconds,"
           "tabi_leaf_block_writeback_seconds,"
           "tabi_leaf_block_verify_checks,"
           "tabi_leaf_block_verify_mismatches,"
           "tabi_tree_preset,tabi_adaptive_level,tabi_cluster_policy,"
           "tabi_target_batch_size,tabi_leaf_pair_min_source_panels,"
           "tabi_leaf_pair_min_targets,tabi_leaf_near_policy,"
           "tabi_interaction_mode,tabi_node_bounds_mode,"
           "tabi_tree_builder,tabi_expansion_center,"
           "tabi_fragmentation_safety,tabi_mac_safety,"
           "tabi_mac_gap_factor,"
           "tabi_tree_depth_max,"
           "tabi_mac_reject_radius_limited,"
           "tabi_mac_accept_rate,tabi_node_radius_max,"
           "tabi_node_radius_mean,tabi_node_aspect_ratio_max,"
           "tabi_node_aspect_ratio_mean,tabi_leaf_panel_work_per_matvec,"
           "tabi_far_node_work_per_matvec,"
           "tabi_far_geometric_unsafe_count,"
           "tabi_far_geometric_unsafe_panels,"
           "tabi_far_geometric_min_gap,"
           "tabi_far_geometric_gap_ratio_min,"
           "tabi_far_geometric_gap_ratio_mean,"
           "tabi_far_geometric_audit_samples,"
           "tabi_mac_physical_gap_rejects,"
           "tabi_mac_convergence_rho_rejects,"
           "tabi_mac_convergence_aspect_rejects,"
           "tabi_mac_convergence_rho_reject_panels,"
           "tabi_mac_convergence_aspect_reject_panels,"
           "tabi_mac_surface_margin_rejects,"
           "tabi_mac_surface_margin_reject_panels,"
           "tabi_mac_surface_margin_applied,"
           "tabi_mac_surface_margin_skipped,"
           "tabi_mac_surface_margin_adaptive_rejects,"
           "tabi_mac_surface_margin_adaptive_reject_panels,"
           "tabi_mac_surface_margin_skipped_accepted_far,"
           "tabi_mac_surface_margin_hard_applied,"
           "tabi_mac_surface_margin_soft_applied,"
           "tabi_mac_surface_margin_hard_rejects,"
           "tabi_mac_surface_margin_hard_reject_panels,"
           "tabi_mac_surface_margin_soft_rejects,"
           "tabi_mac_surface_margin_soft_reject_panels,"
           "tabi_fragmented_far_count,"
           "tabi_fragmented_far_panels,"
           "tabi_fragmented_split_rejects,"
           "tabi_fragmentation_center_gap_ratio_max,"
           "tabi_fragmentation_center_gap_ratio_mean,"
           "tabi_moment_audit,tabi_moment_construction,"
           "tabi_moment_sum,tabi_convergence_audit,"
           "tabi_convergence_safety,tabi_convergence_rho_cap,"
           "tabi_mac_aspect_cap,tabi_convergence_estimator,"
           "tabi_convergence_tail_threshold,"
           "tabi_convergence_estimator_error_cap,"
           "tabi_convergence_unstable_rho_cap,"
           "tabi_convergence_calibration,"
           "tabi_convergence_dielectric_amp,"
           "tabi_convergence_estimator_model,"
           "tabi_convergence_channel_calibration,"
           "tabi_convergence_split_policy,"
           "tabi_convergence_split_selectivity,"
           "tabi_convergence_split_absolute_floor,"
           "tabi_convergence_overreject_audit,"
           "tabi_convergence_overreject_detail,"
           "tabi_convergence_shape_audit,"
           "tabi_convergence_kernel_response_audit,"
           "tabi_channel_missed_semantics,"
           "tabi_accepted_miss_diagnostics,"
           "tabi_accepted_miss_absolute_floor,"
           "tabi_convergence_k3_accepted_rho_cap,"
           "tabi_mac_distance_policy,tabi_mac_surface_margin_factor,"
           "tabi_mac_surface_margin_policy,"
           "tabi_mac_surface_margin_aspect_threshold,"
           "tabi_mac_surface_margin_rho_trigger,"
           "tabi_mac_surface_margin_soft_aspect_threshold,"
           "tabi_mac_surface_margin_soft_rho_trigger,"
           "tabi_mac_surface_margin_soft_factor,"
           "tabi_convergence_estimator_samples,"
           "tabi_convergence_estimator_unstable_far_count,"
           "tabi_convergence_estimator_unstable_far_panels,"
           "tabi_convergence_estimator_split_rejects,"
           "tabi_convergence_estimator_split_reject_panels,"
           "tabi_convergence_low_abs_reaccepted_far_count,"
           "tabi_convergence_low_abs_reaccepted_far_panels,"
           "tabi_convergence_estimator_tail_ratio_max,"
           "tabi_convergence_estimator_tail_ratio_mean,"
           "tabi_convergence_estimator_error_max,"
           "tabi_convergence_estimator_error_mean,"
           "tabi_convergence_estimator_abs_error_max,"
           "tabi_convergence_estimator_abs_error_mean,"
           "tabi_convergence_tail_norm_p_max,"
           "tabi_convergence_tail_norm_p_mean,"
           "tabi_convergence_tail_norm_pm1_max,"
           "tabi_convergence_tail_norm_pm1_mean,"
           "tabi_convergence_effective_jcp_amp_max,"
           "tabi_convergence_effective_jcp_amp_mean,"
           "tabi_convergence_dielectric_amp_max,"
           "tabi_convergence_dielectric_amp_mean,"
           "tabi_convergence_k1_effective_jcp_amp_max,"
           "tabi_convergence_k2_effective_jcp_amp_max,"
           "tabi_convergence_k3_effective_jcp_amp_max,"
           "tabi_convergence_k4_effective_jcp_amp_max,"
           "tabi_convergence_k1_effective_jcp_amp_mean,"
           "tabi_convergence_k2_effective_jcp_amp_mean,"
           "tabi_convergence_k3_effective_jcp_amp_mean,"
           "tabi_convergence_k4_effective_jcp_amp_mean,"
           "tabi_convergence_k1_estimated_error_max,"
           "tabi_convergence_k2_estimated_error_max,"
           "tabi_convergence_k3_estimated_error_max,"
           "tabi_convergence_k4_estimated_error_max,"
           "tabi_convergence_k1_unstable_far_count,"
           "tabi_convergence_k2_unstable_far_count,"
           "tabi_convergence_k3_unstable_far_count,"
           "tabi_convergence_k4_unstable_far_count,"
           "tabi_convergence_k1_unstable_far_panels,"
           "tabi_convergence_k2_unstable_far_panels,"
           "tabi_convergence_k3_unstable_far_panels,"
           "tabi_convergence_k4_unstable_far_panels,"
           "tabi_convergence_k1_split_rejects,"
           "tabi_convergence_k2_split_rejects,"
           "tabi_convergence_k3_split_rejects,"
           "tabi_convergence_k4_split_rejects,"
           "tabi_convergence_k1_split_reject_panels,"
           "tabi_convergence_k2_split_reject_panels,"
           "tabi_convergence_k3_split_reject_panels,"
           "tabi_convergence_k4_split_reject_panels,"
           "tabi_convergence_k1_low_abs_reaccepted,"
           "tabi_convergence_k2_low_abs_reaccepted,"
           "tabi_convergence_k3_low_abs_reaccepted,"
           "tabi_convergence_k4_low_abs_reaccepted,"
           "tabi_convergence_k1_low_abs_reaccepted_panels,"
           "tabi_convergence_k2_low_abs_reaccepted_panels,"
           "tabi_convergence_k3_low_abs_reaccepted_panels,"
           "tabi_convergence_k4_low_abs_reaccepted_panels,"
           "tabi_convergence_k1_high_error_caught,"
           "tabi_convergence_k2_high_error_caught,"
           "tabi_convergence_k3_high_error_caught,"
           "tabi_convergence_k4_high_error_caught,"
           "tabi_convergence_k1_high_error_missed,"
           "tabi_convergence_k2_high_error_missed,"
           "tabi_convergence_k3_high_error_missed,"
           "tabi_convergence_k4_high_error_missed,"
           "tabi_convergence_k1_high_error_caught_by_split,"
           "tabi_convergence_k2_high_error_caught_by_split,"
           "tabi_convergence_k3_high_error_caught_by_split,"
           "tabi_convergence_k4_high_error_caught_by_split,"
           "tabi_convergence_estimator_high_error_caught,"
           "tabi_convergence_estimator_high_error_missed,"
           "tabi_convergence_estimator_high_error_caught_by_split,"
           "tabi_convergence_estimator_phantom_high_error_count,"
           "tabi_convergence_estimator_phantom_high_error_abs_max,"
           "tabi_convergence_overreject_count,"
           "tabi_convergence_overreject_panels,"
           "tabi_convergence_overreject_caught_count,"
           "tabi_convergence_overreject_caught_panels,"
           "tabi_convergence_overreject_abs_error_max,"
           "tabi_convergence_overreject_abs_error_mean,"
           "tabi_convergence_overreject_abs_error_p95,"
           "tabi_convergence_overreject_abs_error_p99,"
           "tabi_convergence_overreject_estimated_abs_error_max,"
           "tabi_convergence_overreject_estimated_abs_error_mean,"
           "tabi_convergence_overreject_required_scale_max,"
           "tabi_convergence_overreject_required_scale_p95,"
           "tabi_convergence_overreject_required_scale_p99,"
           "tabi_convergence_k1_overreject_count,"
           "tabi_convergence_k2_overreject_count,"
           "tabi_convergence_k3_overreject_count,"
           "tabi_convergence_k4_overreject_count,"
           "tabi_convergence_k1_overreject_caught_count,"
           "tabi_convergence_k2_overreject_caught_count,"
           "tabi_convergence_k3_overreject_caught_count,"
           "tabi_convergence_k4_overreject_caught_count,"
           "tabi_convergence_k1_overreject_panels,"
           "tabi_convergence_k2_overreject_panels,"
           "tabi_convergence_k3_overreject_panels,"
           "tabi_convergence_k4_overreject_panels,"
           "tabi_convergence_k1_overreject_abs_error_max,"
           "tabi_convergence_k2_overreject_abs_error_max,"
           "tabi_convergence_k3_overreject_abs_error_max,"
           "tabi_convergence_k4_overreject_abs_error_max,"
           "tabi_convergence_k1_overreject_abs_error_mean,"
           "tabi_convergence_k2_overreject_abs_error_mean,"
           "tabi_convergence_k3_overreject_abs_error_mean,"
           "tabi_convergence_k4_overreject_abs_error_mean,"
           "tabi_convergence_k1_overreject_required_scale_max,"
           "tabi_convergence_k2_overreject_required_scale_max,"
           "tabi_convergence_k3_overreject_required_scale_max,"
           "tabi_convergence_k4_overreject_required_scale_max,"
           "tabi_convergence_k1_overreject_required_scale_p99,"
           "tabi_convergence_k2_overreject_required_scale_p99,"
           "tabi_convergence_k3_overreject_required_scale_p99,"
           "tabi_convergence_k4_overreject_required_scale_p99,"
           "tabi_phantom_high_error_k1_count,"
           "tabi_phantom_high_error_k2_count,"
           "tabi_phantom_high_error_k3_count,"
           "tabi_phantom_high_error_k4_count,"
           "tabi_phantom_high_error_k1_abs_max,"
           "tabi_phantom_high_error_k2_abs_max,"
           "tabi_phantom_high_error_k3_abs_max,"
           "tabi_phantom_high_error_k4_abs_max,"
           "tabi_accepted_miss_worst_channel,"
           "tabi_accepted_miss_k1_count,tabi_accepted_miss_k2_count,"
           "tabi_accepted_miss_k3_count,tabi_accepted_miss_k4_count,"
           "tabi_accepted_miss_k1_panels,tabi_accepted_miss_k2_panels,"
           "tabi_accepted_miss_k3_panels,tabi_accepted_miss_k4_panels,"
           "tabi_accepted_miss_k1_rho_min,tabi_accepted_miss_k2_rho_min,"
           "tabi_accepted_miss_k3_rho_min,tabi_accepted_miss_k4_rho_min,"
           "tabi_accepted_miss_k1_rho_mean,"
           "tabi_accepted_miss_k2_rho_mean,"
           "tabi_accepted_miss_k3_rho_mean,"
           "tabi_accepted_miss_k4_rho_mean,"
           "tabi_accepted_miss_k1_rho_max,tabi_accepted_miss_k2_rho_max,"
           "tabi_accepted_miss_k3_rho_max,tabi_accepted_miss_k4_rho_max,"
           "tabi_accepted_miss_k1_tail_ratio_mean,"
           "tabi_accepted_miss_k2_tail_ratio_mean,"
           "tabi_accepted_miss_k3_tail_ratio_mean,"
           "tabi_accepted_miss_k4_tail_ratio_mean,"
           "tabi_accepted_miss_k1_tail_ratio_max,"
           "tabi_accepted_miss_k2_tail_ratio_max,"
           "tabi_accepted_miss_k3_tail_ratio_max,"
           "tabi_accepted_miss_k4_tail_ratio_max,"
           "tabi_accepted_miss_k1_estimated_error_mean,"
           "tabi_accepted_miss_k2_estimated_error_mean,"
           "tabi_accepted_miss_k3_estimated_error_mean,"
           "tabi_accepted_miss_k4_estimated_error_mean,"
           "tabi_accepted_miss_k1_estimated_error_max,"
           "tabi_accepted_miss_k2_estimated_error_max,"
           "tabi_accepted_miss_k3_estimated_error_max,"
           "tabi_accepted_miss_k4_estimated_error_max,"
           "tabi_convergence_k3_rho_cap_rejects,"
           "tabi_convergence_k3_rho_cap_reject_panels,"
           "tabi_moment_audit_nodes,"
           "tabi_moment_translation_rel_error_max,"
           "tabi_moment_translation_rel_error_mean,"
           "tabi_moment_translation_offset_max,"
           "tabi_moment_translation_offset_mean,"
           "tabi_moment_cancellation_ratio_max,"
           "tabi_moment_cancellation_ratio_mean,"
           "tabi_moment_dynamic_range_max,"
           "tabi_expansion_audit,tabi_expansion_audit_channel,"
           "tabi_far_expansion_audit_pairs,"
           "tabi_far_expansion_audit_panels,"
           "tabi_far_expansion_direct_seconds,"
           "tabi_far_expansion_rel_error_max,"
           "tabi_far_expansion_rel_error_mean,"
           "tabi_far_expansion_rel_error_p95,"
           "tabi_far_expansion_rel_error_p99,"
           "tabi_far_expansion_worst_channel,"
           "tabi_far_expansion_k1_rel_error_max,"
           "tabi_far_expansion_k2_rel_error_max,"
           "tabi_far_expansion_k3_rel_error_max,"
           "tabi_far_expansion_k4_rel_error_max,"
           "tabi_far_expansion_k1_rel_error_mean,"
           "tabi_far_expansion_k2_rel_error_mean,"
           "tabi_far_expansion_k3_rel_error_mean,"
           "tabi_far_expansion_k4_rel_error_mean,"
           "tabi_far_expansion_k1_rel_error_p99,"
           "tabi_far_expansion_k2_rel_error_p99,"
           "tabi_far_expansion_k3_rel_error_p99,"
           "tabi_far_expansion_k4_rel_error_p99,"
           "tabi_far_expansion_mean_gate_pass,"
           "tabi_far_convergence_rho_max,"
           "tabi_far_convergence_rho_mean,"
           "tabi_accepted_far_aspect_ratio_max,"
           "tabi_accepted_far_aspect_ratio_mean,"
           "tabi_quad_buffer_auto_candidate,"
           "tabi_near_integration,tabi_quad_buffer,"
           "tabi_near_self_pairs,"
           "tabi_near_one_ring_pairs,tabi_near_duffy_pairs,"
           "tabi_near_centroid_pairs,tabi_near_fallbacks,"
           "tabi_near_geometric_buffer_pairs,"
           "tabi_near_geometric_buffer_gauss_pairs,"
           "tabi_near_geometric_buffer_adaptive_pairs,"
           "tabi_near_geometric_buffer_centroid_pairs,"
           "tabi_near_geometric_buffer_min_gap,"
           "tabi_near_geometric_buffer_mean_gap_ratio,"
           "tabi_near_geometric_buffer_max_edge_mean,"
           "tabi_near_integration_seconds,tabi_duffy_subdivisions_max,"
           "tabi_duffy_subdivisions_mean,tabi_near_kernel_policy,"
           "tabi_near_k1_accepted,tabi_near_k1_rejected,"
           "tabi_near_k2_accepted,tabi_near_k2_rejected,"
           "tabi_near_k3_accepted,tabi_near_k3_rejected,"
           "tabi_near_k4_accepted,tabi_near_k4_rejected,"
           "tabi_near_one_ring_k1_accepted,"
           "tabi_near_one_ring_k1_rejected,"
           "tabi_near_one_ring_k2_accepted,"
           "tabi_near_one_ring_k2_rejected,"
           "tabi_near_one_ring_k3_accepted,"
           "tabi_near_one_ring_k3_rejected,"
           "tabi_near_one_ring_k4_accepted,"
           "tabi_near_one_ring_k4_rejected,"
           "tabi_near_k1_policy,tabi_near_k2_policy,"
           "tabi_near_k3_policy,tabi_near_k4_policy,"
           "tabi_surface_geometry,tabi_mesh_vertices,tabi_mesh_triangles,"
           "tabi_boundary_nodes,tabi_avg_valence,tabi_nonmanifold_edges,"
           "tabi_solid_angle_fallbacks,tabi_diag_min,tabi_diag_max,"
            "tabi_taylor_order,tabi_derivative_mode,"
            "tabi_interaction_build_seconds,"
           "tabi_tree_eval_seconds,"
           "tabi_formulation,tabi_energy_projection,"
           "tabi_preconditioner_requested,tabi_preconditioner_effective,"
           "tabi_preconditioner_block_near,"
           "tabi_preconditioner_block_panel_cap,"
           "tabi_preconditioner_overlap_panel_cap,"
           "tabi_preconditioner_ilu_neighbor_cap,"
           "tabi_preconditioner_ilu_drop_tol,"
           "tabi_preconditioner_memory_budget_bytes,"
           "tabi_preconditioner_estimated_bytes,"
           "tabi_preconditioner_valid_blocks,"
           "tabi_preconditioner_valid_panels,"
           "tabi_preconditioner_ill_conditioned_blocks,"
           "tabi_preconditioner_quality_rejected_blocks,"
           "tabi_preconditioner_quality_probe_max_gain,"
           "tabi_preconditioner_quality_status,"
           "tabi_preconditioner_setup_seconds,"
           "tabi_preconditioner_apply_seconds,"
           "tabi_preconditioner_apply_calls,"
           "tabi_preconditioner_sparse_nnz_blocks,"
           "tabi_preconditioner_sparse_avg_degree,"
           "tabi_preconditioner_sparse_max_degree,"
           "tabi_preconditioner_coarse_aggregate_cap,"
           "tabi_preconditioner_coarse_aggregates,"
           "tabi_preconditioner_coarse_dimension,"
           "tabi_preconditioner_coarse_base,"
           "tabi_preconditioner_coarse_matrix_bytes,"
           "tabi_preconditioner_coarse_pivot_ratio,"
           "tabi_preconditioner_coarse_status,"
           "tabi_preconditioner_coarse_setup_matvec_calls,"
           "tabi_preconditioner_coarse_setup_seconds,"
           "tabi_preconditioner_coarse_apply_seconds,"
           "tabi_preconditioner_coarse_apply_calls,"
           "tabi_preconditioner_coarse_rp_identity_error_max,"
           "tabi_preconditioner_coarse_verify_rel_error_max,"
           "tabi_surface_generator,tabi_discretization,"
           "tabi_fallback_used,"
           "tabi_rhs_scale,tabi_dielectric_contrast,tabi_operator_contrast,"
           "tabi_unknown_count,tabi_bie_system,tabi_lhs_diag_phi,"
           "tabi_lhs_diag_phin,tabi_energy_equation,"
            "tabi_kernel_convention,tabi_source_norm_s1,"
           "tabi_source_norm_s2,tabi_phi_norm,tabi_phin_norm,"
           "tabi_kappa_inv_A,tabi_coulomb_factor_kcal_mol,"
           "tabi_boundary_energy_kJmol,tabi_reference_energy_kJmol,"
           "tabi_selected_energy_scale,"
           "setup_seconds,surface_seconds,epsilon_seconds,boundary_seconds,"
           "solver_seconds,energy_seconds,map_io_seconds,rpb_source_seconds,"
           "focus_setup_seconds,probe_seconds,focus_parent_solve_seconds,"
           "focus_child_setup_seconds,focus_boundary_interp_seconds,"
           "focus_child_surface_seconds,focus_child_epsilon_seconds,"
           "focus_child_solver_seconds,probe_energy_seconds,"
           "audit_output_seconds,nwt_half_update_seconds,"
           "nwt_rms_reduce_seconds,nwt_operator_residual_seconds,"
           "nwt_source_update_seconds,nwt_progress_io_seconds,"
           "nwt_coeff_pack_seconds,nwt_homo_update_seconds,"
           "nwt_hetero_update_seconds,nwt_nonlin_term_seconds,"
           "nwt_reduce_seconds,nwt_homo_stencil_seconds,"
           "nwt_homo_salt_math_seconds,nwt_homo_write_seconds,"
           "nwt_thread_wait_seconds,nwt_block_overhead_seconds,"
           "linear_half_update_seconds,"
           "nwt_homo_cells,nwt_homo_linear_cells,nwt_homo_salt_cells,"
           "nwt_homo_updates,nwt_dispatch_calls,nwt_homo_ns_per_update,"
           "epsilon_center_density_seconds,epsilon_midpoint_density_seconds,"
           "epsilon_ion_exclusion_seconds,epsilon_finalize_seconds,"
           "total_seconds,iterations,blocks,"
           "status\n";
    timing_header_written_ = true;
  }
}

void AnalyzerImpl::append_timing_row(const FrameResult &result,
                                     const TimingBuckets &timing) {
  if (!pbe_timing_enabled())
    return;
  write_timing_header_if_needed();
  auto seconds = [](double ms) { return ms / 1000.0; };
  const std::size_t iterations =
      result.solvent_iterations + result.vacuum_iterations;
  const std::string mode =
      result.pbe_model_effective.empty()
          ? result.solver_mode
          : (result.pbe_model_effective + ":" + result.solver_mode);
  std::string scope = "single_grid";
  if (result.focus_enabled && result.focus_levels > 1)
    scope = result.probe_status == "off" ? "focus_full_workflow"
                                         : "focus_full_workflow_probe";
  else if (result.focus_enabled)
    scope = result.probe_status == "off" ? "focus_single_grid"
                                         : "focus_single_grid_probe";
  else if (result.pbe_biomodel == "rpbe")
    scope = "rpbe_full_workflow";
  else if (result.pbe_backend == "tabi")
    scope = "tabi_boundary_integral";
  timing_out_ << result.trajectory_id << ',' << result.frame_index << ','
              << result.global_index << ',' << std::setprecision(10)
              << result.time_ps << ',' << csv_escape(current_source_label_)
              << ',' << csv_escape(mode) << ','
              << csv_escape(result.pbe_backend) << ',' << result.atom_count
              << ',' << timing.focus_child_atom_count << ',' << configured_threads_
              << ',' << csv_escape(scope) << ','
              << seconds(result.tabi_surface_build_ms) << ','
              << result.tabi_surface_panels << ','
              << result.tabi_surface_point_count_requested << ','
              << result.tabi_panel_density_per_a2 << ','
              << csv_escape(result.tabi_panel_budget_status) << ','
              << result.tabi_surface_exposed_points << ','
              << result.tabi_surface_buried_points << ','
              << seconds(result.tabi_tree_build_ms) << ','
              << result.tabi_tree_nodes << ','
              << result.tabi_tree_leaf_count << ','
              << result.tabi_near_interactions << ','
              << result.tabi_far_interactions << ','
              << csv_escape(result.tabi_matvec_mode) << ','
              << csv_escape(result.tabi_status) << ','
              << csv_escape(result.tabi_convergence_reason) << ','
              << csv_escape(result.tabi_reaction_energy_source) << ','
              << result.tabi_gmres_restart_requested << ','
              << result.tabi_gmres_restart_effective << ','
              << result.tabi_best_residual << ','
              << result.tabi_best_residual_iter << ','
              << result.tabi_final_residual << ','
              << result.tabi_best_preconditioned_residual << ','
              << result.tabi_final_preconditioned_residual << ','
              << result.tabi_true_residual_abs << ','
              << result.tabi_true_residual_rel << ','
              << csv_escape(result.tabi_true_residual_probe_mode) << ','
              << csv_escape(result.tabi_true_residual_probe_scope) << ','
              << result.tabi_true_residual_probe_interval << ','
              << result.tabi_true_residual_sample_count << ','
              << result.tabi_true_residual_probe_calls << ','
              << seconds(result.tabi_true_residual_probe_ms) << ','
              << result.tabi_sampled_true_residual_rel << ','
              << result.tabi_sampled_true_residual_rows << ','
              << result.tabi_sampled_probe_calls << ','
              << seconds(result.tabi_sampled_probe_ms) << ','
              << result.tabi_full_true_residual_confirm_calls << ','
              << seconds(result.tabi_full_probe_ms) << ','
              << seconds(result.tabi_rhs_ms) << ','
              << seconds(result.tabi_gmres_ms) << ','
              << seconds(result.tabi_matvec_ms_total) << ','
              << seconds(result.tabi_matvec_ms_avg) << ','
              << result.tabi_matvec_calls << ','
              << result.tabi_tree_vs_direct_matvec_rel << ','
              << result.tabi_matvec_audit_rows << ','
              << seconds(result.tabi_matvec_audit_ms) << ','
              << seconds(result.tabi_orthogonalization_ms) << ','
              << seconds(result.tabi_preconditioner_ms) << ','
              << seconds(result.tabi_surface_neighbor_ms) << ','
              << seconds(result.tabi_surface_mesh_ms) << ','
              << seconds(result.tabi_solid_angle_ms) << ','
              << seconds(result.tabi_tree_moment_ms) << ','
              << seconds(result.tabi_tree_traversal_ms) << ','
              << csv_escape(result.tabi_matvec_cost_profile) << ','
              << result.tabi_leaf_panels_per_second << ','
              << result.tabi_clean_leaf_panels_per_second << ','
              << result.tabi_tree_traversal_fraction << ','
              << seconds(result.tabi_near_eval_ms) << ','
              << seconds(result.tabi_far_eval_ms) << ','
              << seconds(result.tabi_singular_near_eval_ms) << ','
              << seconds(result.tabi_leaf_near_eval_ms) << ','
              << seconds(result.tabi_leaf_panel_eval_ms) << ','
              << seconds(result.tabi_leaf_tile_eval_ms) << ','
              << seconds(result.tabi_leaf_residual_eval_ms) << ','
              << seconds(result.tabi_smooth_leaf_node_eval_ms) << ','
              << csv_escape(result.tabi_leaf_panel_kernel) << ','
              << csv_escape(result.tabi_leaf_tile_status) << ','
              << result.tabi_leaf_panel_pairs_per_second << ','
              << result.tabi_singular_near_interactions << ','
              << result.tabi_leaf_near_interactions << ','
              << result.tabi_far_node_interactions << ','
              << result.tabi_leaf_panel_interactions << ','
              << result.tabi_leaf_tile_interactions << ','
              << result.tabi_leaf_tile_panel_interactions << ','
              << result.tabi_leaf_residual_panel_interactions << ','
              << csv_escape(result.tabi_leaf_work_profile) << ','
              << result.tabi_leaf_work_mac_terminal_panels << ','
              << result.tabi_leaf_work_physical_gap_panels << ','
              << result.tabi_leaf_work_fragmentation_panels << ','
              << result.tabi_leaf_work_convergence_panels << ','
              << result.tabi_leaf_work_convergence_estimator_panels << ','
              << result.tabi_leaf_work_convergence_rho_panels << ','
              << result.tabi_leaf_work_convergence_aspect_panels << ','
              << result.tabi_leaf_work_surface_margin_panels << ','
              << result.tabi_leaf_work_surface_margin_hard_panels << ','
              << result.tabi_leaf_work_surface_margin_soft_panels << ','
              << result.tabi_leaf_work_k3_rho_cap_panels << ','
              << result.tabi_leaf_work_singular_fallback_panels << ','
              << result.tabi_leaf_work_scalar_tail_panels << ','
              << result.tabi_leaf_work_profile_accounted_panels << ','
              << result.tabi_leaf_work_mac_terminal_frac << ','
              << result.tabi_leaf_work_convergence_frac << ','
              << result.tabi_leaf_work_surface_margin_frac << ','
              << result.tabi_leaf_work_singular_fallback_frac << ','
              << result.tabi_leaf_work_scalar_tail_frac << ','
              << result.tabi_smooth_leaf_node_interactions << ','
              << result.tabi_leaf_pruned_interactions << ','
              << result.tabi_leaf_mac_candidates << ','
              << result.tabi_leaf_mac_accepted << ','
              << result.tabi_leaf_mac_rejected << ','
              << result.tabi_leaf_mac_far_panels << ','
              << result.tabi_leaf_mac_direct_panels_saved_estimate << ','
              << result.tabi_adaptive_smooth_leaf_accepted << ','
              << result.tabi_adaptive_smooth_leaf_rejected << ','
              << result.tabi_adaptive_span_tile_accepted << ','
              << result.tabi_adaptive_span_tile_rejected << ','
              << result.tabi_adaptive_exact_leaf_fallbacks << ','
              << result.tabi_adaptive_reject_one_ring_overlap << ','
              << result.tabi_adaptive_reject_error_bound_failed << ','
              << result.tabi_adaptive_reject_p4_check_failed << ','
              << result.tabi_adaptive_reject_too_close << ','
              << result.tabi_leaf_work_reduction_ratio << ','
              << result.tabi_target_cluster_interactions << ','
              << result.tabi_target_cluster_targets << ','
              << result.tabi_target_cluster_panels_pruned << ','
              << result.tabi_target_cluster_rejected << ','
              << seconds(result.tabi_target_cluster_eval_ms) << ','
              << result.tabi_leaf_pair_batch_interactions << ','
              << result.tabi_leaf_pair_batch_targets << ','
              << result.tabi_leaf_pair_batch_source_panels << ','
              << result.tabi_leaf_pair_batch_residual_panels << ','
              << result.tabi_leaf_pair_batch_fallbacks << ','
              << seconds(result.tabi_leaf_pair_batch_ms) << ','
              << csv_escape(result.tabi_leaf_pair_batch_kernel) << ','
              << result.tabi_leaf_pair_batch_vectorized_sources << ','
              << result.tabi_leaf_pair_batch_vector_lanes << ','
              << result.tabi_leaf_pair_batch_pairs_per_second << ','
              << csv_escape(result.tabi_leaf_kernel_profile) << ','
              << csv_escape(result.tabi_jcp_direct_kernel) << ','
              << result.tabi_leaf_kernel_scalar_tasks << ','
              << result.tabi_leaf_kernel_batch_tasks << ','
              << result.tabi_leaf_kernel_source_panel_evals << ','
              << seconds(result.tabi_leaf_kernel_scalar_ms) << ','
              << seconds(result.tabi_leaf_kernel_target_avx2_ms) << ','
              << seconds(result.tabi_leaf_kernel_tail_ms) << ','
              << result.tabi_leaf_pair_batch_soa_sources << ','
              << seconds(result.tabi_leaf_pair_batch_soa_ms) << ','
              << seconds(result.tabi_leaf_pair_batch_scalar_tail_ms) << ','
              << result.tabi_leaf_pair_batch_clean_source_panels_per_second
              << ','
              << csv_escape(result.tabi_leaf_pair_selection_profile) << ','
              << csv_escape(result.tabi_leaf_pair_probe_cache) << ','
              << csv_escape(result.tabi_leaf_pair_target_group) << ','
              << result.tabi_leaf_pair_max_lookahead << ','
              << result.tabi_leaf_pair_selection_candidate_scans << ','
              << result.tabi_leaf_pair_selection_clean_rejects << ','
              << result.tabi_leaf_pair_selection_probe_calls << ','
              << result.tabi_leaf_pair_selection_probe_node_visits << ','
              << seconds(result.tabi_leaf_pair_selection_probe_ms) << ','
              << seconds(result.tabi_leaf_pair_selection_ms) << ','
              << result.tabi_leaf_pair_selection_probe_cache_hits << ','
              << result.tabi_leaf_pair_selection_probe_cache_misses << ','
              << result.tabi_leaf_pair_selection_insufficient_targets << ','
              << result.tabi_leaf_pair_selection_same_target_leaf_rejects
              << ',' << csv_escape(result.tabi_leaf_pair_directness_probe)
              << ',' << csv_escape(result.tabi_leaf_pair_directness_cache)
              << ',' << csv_escape(result.tabi_leaf_pair_directness_verify)
              << ','
              << result.tabi_leaf_pair_directness_ancestor_probe_calls << ','
              << result.tabi_leaf_pair_directness_ancestor_probe_node_visits
              << ','
              << seconds(result.tabi_leaf_pair_directness_ancestor_probe_ms)
              << ',' << result.tabi_leaf_pair_directness_pair_cache_hits
              << ',' << result.tabi_leaf_pair_directness_pair_cache_misses
              << ',' << result.tabi_leaf_pair_directness_all_direct_proofs
              << ',' << result.tabi_leaf_pair_directness_mixed_fallbacks
              << ',' << result.tabi_leaf_pair_directness_verify_checks
              << ',' << result.tabi_leaf_pair_directness_verify_mismatches
              << ',' << csv_escape(result.tabi_leaf_pair_scheduler)
              << ',' << result.tabi_leaf_pair_worklist_max_active
              << ',' << csv_escape(result.tabi_leaf_pair_worklist_verify)
              << ',' << result.tabi_leaf_pair_worklist_tasks
              << ',' << result.tabi_leaf_pair_worklist_groups
              << ',' << result.tabi_leaf_pair_worklist_batched_groups
              << ',' << result.tabi_leaf_pair_worklist_scalar_flush_tasks
              << ',' << result.tabi_leaf_pair_worklist_forced_flushes
              << ',' << result.tabi_leaf_pair_worklist_active_groups_max
              << ',' << seconds(result.tabi_leaf_pair_worklist_ms)
              << ',' << result.tabi_leaf_pair_worklist_verify_checks
              << ',' << result.tabi_leaf_pair_worklist_verify_mismatches
              << ',' << csv_escape(result.tabi_leaf_block_verify)
              << ',' << result.tabi_leaf_block_target_leaves
              << ',' << result.tabi_leaf_block_targets
              << ',' << result.tabi_leaf_block_source_groups
              << ',' << result.tabi_leaf_block_batched_groups
              << ',' << result.tabi_leaf_block_scalar_tail_tasks
              << ',' << seconds(result.tabi_leaf_block_group_ms)
              << ',' << seconds(result.tabi_leaf_block_flush_ms)
              << ',' << seconds(result.tabi_leaf_block_writeback_ms)
              << ',' << result.tabi_leaf_block_verify_checks
              << ',' << result.tabi_leaf_block_verify_mismatches
              << ','
              << csv_escape(result.tabi_tree_preset) << ','
              << csv_escape(result.tabi_adaptive_level) << ','
              << csv_escape(result.tabi_cluster_policy) << ','
              << result.tabi_target_batch_size << ','
              << result.tabi_leaf_pair_min_source_panels << ','
              << result.tabi_leaf_pair_min_targets << ','
              << csv_escape(result.tabi_leaf_near_policy) << ','
              << csv_escape(result.tabi_interaction_mode) << ','
              << csv_escape(result.tabi_node_bounds_mode) << ','
              << csv_escape(result.tabi_tree_builder) << ','
              << csv_escape(result.tabi_expansion_center) << ','
              << csv_escape(result.tabi_fragmentation_safety) << ','
              << csv_escape(result.tabi_mac_safety) << ','
              << result.tabi_mac_gap_factor << ','
              << result.tabi_tree_depth_max << ','
              << result.tabi_mac_reject_radius_limited << ','
              << result.tabi_mac_accept_rate << ','
              << result.tabi_node_radius_max << ','
              << result.tabi_node_radius_mean << ','
              << result.tabi_node_aspect_ratio_max << ','
              << result.tabi_node_aspect_ratio_mean << ','
              << result.tabi_leaf_panel_work_per_matvec << ','
              << result.tabi_far_node_work_per_matvec << ','
              << result.tabi_far_geometric_unsafe_count << ','
              << result.tabi_far_geometric_unsafe_panels << ','
              << result.tabi_far_geometric_min_gap << ','
              << result.tabi_far_geometric_gap_ratio_min << ','
              << result.tabi_far_geometric_gap_ratio_mean << ','
              << result.tabi_far_geometric_audit_samples << ','
              << result.tabi_mac_physical_gap_rejects << ','
              << result.tabi_mac_convergence_rho_rejects << ','
              << result.tabi_mac_convergence_aspect_rejects << ','
              << result.tabi_mac_convergence_rho_reject_panels << ','
              << result.tabi_mac_convergence_aspect_reject_panels << ','
              << result.tabi_mac_surface_margin_rejects << ','
              << result.tabi_mac_surface_margin_reject_panels << ','
              << result.tabi_mac_surface_margin_applied << ','
              << result.tabi_mac_surface_margin_skipped << ','
              << result.tabi_mac_surface_margin_adaptive_rejects << ','
              << result.tabi_mac_surface_margin_adaptive_reject_panels << ','
              << result.tabi_mac_surface_margin_skipped_accepted_far << ','
              << result.tabi_mac_surface_margin_hard_applied << ','
              << result.tabi_mac_surface_margin_soft_applied << ','
              << result.tabi_mac_surface_margin_hard_rejects << ','
              << result.tabi_mac_surface_margin_hard_reject_panels << ','
              << result.tabi_mac_surface_margin_soft_rejects << ','
              << result.tabi_mac_surface_margin_soft_reject_panels << ','
              << result.tabi_fragmented_far_count << ','
              << result.tabi_fragmented_far_panels << ','
              << result.tabi_fragmented_split_rejects << ','
              << result.tabi_fragmentation_center_gap_ratio_max << ','
              << result.tabi_fragmentation_center_gap_ratio_mean << ','
              << csv_escape(result.tabi_moment_audit) << ','
              << csv_escape(result.tabi_moment_construction) << ','
              << csv_escape(result.tabi_moment_sum) << ','
              << csv_escape(result.tabi_convergence_audit) << ','
              << csv_escape(result.tabi_convergence_safety) << ','
              << result.tabi_convergence_rho_cap << ','
              << result.tabi_mac_aspect_cap << ','
              << csv_escape(result.tabi_convergence_estimator) << ','
              << result.tabi_convergence_tail_threshold << ','
              << result.tabi_convergence_estimator_error_cap << ','
              << result.tabi_convergence_unstable_rho_cap << ','
              << csv_escape(result.tabi_convergence_calibration) << ','
              << csv_escape(result.tabi_convergence_dielectric_amp) << ','
              << csv_escape(result.tabi_convergence_estimator_model) << ','
              << csv_escape(result.tabi_convergence_channel_calibration)
              << ','
              << csv_escape(result.tabi_convergence_split_policy) << ','
              << csv_escape(result.tabi_convergence_split_selectivity) << ','
               << result.tabi_convergence_split_absolute_floor << ','
               << csv_escape(result.tabi_convergence_overreject_audit) << ','
               << csv_escape(result.tabi_convergence_overreject_detail) << ','
               << csv_escape(result.tabi_convergence_shape_audit) << ','
               << csv_escape(result.tabi_convergence_kernel_response_audit) << ','
               << csv_escape(result.tabi_channel_missed_semantics) << ','
              << csv_escape(result.tabi_accepted_miss_diagnostics) << ','
              << result.tabi_accepted_miss_absolute_floor << ','
              << result.tabi_convergence_k3_accepted_rho_cap << ','
              << csv_escape(result.tabi_mac_distance_policy) << ','
              << result.tabi_mac_surface_margin_factor << ','
              << csv_escape(result.tabi_mac_surface_margin_policy) << ','
              << result.tabi_mac_surface_margin_aspect_threshold << ','
              << result.tabi_mac_surface_margin_rho_trigger << ','
              << result.tabi_mac_surface_margin_soft_aspect_threshold << ','
              << result.tabi_mac_surface_margin_soft_rho_trigger << ','
              << result.tabi_mac_surface_margin_soft_factor << ','
              << result.tabi_convergence_estimator_samples << ','
              << result.tabi_convergence_estimator_unstable_far_count << ','
              << result.tabi_convergence_estimator_unstable_far_panels << ','
              << result.tabi_convergence_estimator_split_rejects << ','
              << result.tabi_convergence_estimator_split_reject_panels << ','
              << result.tabi_convergence_low_abs_reaccepted_far_count << ','
              << result.tabi_convergence_low_abs_reaccepted_far_panels << ','
              << result.tabi_convergence_estimator_tail_ratio_max << ','
              << result.tabi_convergence_estimator_tail_ratio_mean << ','
              << result.tabi_convergence_estimator_error_max << ','
              << result.tabi_convergence_estimator_error_mean << ','
              << result.tabi_convergence_estimator_abs_error_max << ','
              << result.tabi_convergence_estimator_abs_error_mean << ','
              << result.tabi_convergence_tail_norm_p_max << ','
              << result.tabi_convergence_tail_norm_p_mean << ','
              << result.tabi_convergence_tail_norm_pm1_max << ','
              << result.tabi_convergence_tail_norm_pm1_mean << ','
              << result.tabi_convergence_effective_jcp_amp_max << ','
              << result.tabi_convergence_effective_jcp_amp_mean << ','
              << result.tabi_convergence_dielectric_amp_max << ','
              << result.tabi_convergence_dielectric_amp_mean << ','
              << result.tabi_convergence_channel_effective_jcp_amp_max[0]
              << ','
              << result.tabi_convergence_channel_effective_jcp_amp_max[1]
              << ','
              << result.tabi_convergence_channel_effective_jcp_amp_max[2]
              << ','
              << result.tabi_convergence_channel_effective_jcp_amp_max[3]
              << ','
              << result.tabi_convergence_channel_effective_jcp_amp_mean[0]
              << ','
              << result.tabi_convergence_channel_effective_jcp_amp_mean[1]
              << ','
              << result.tabi_convergence_channel_effective_jcp_amp_mean[2]
              << ','
              << result.tabi_convergence_channel_effective_jcp_amp_mean[3]
              << ','
              << result.tabi_convergence_channel_estimated_error_max[0]
              << ','
              << result.tabi_convergence_channel_estimated_error_max[1]
              << ','
              << result.tabi_convergence_channel_estimated_error_max[2]
              << ','
              << result.tabi_convergence_channel_estimated_error_max[3]
              << ','
              << result.tabi_convergence_channel_unstable_far_count[0]
              << ','
              << result.tabi_convergence_channel_unstable_far_count[1]
              << ','
              << result.tabi_convergence_channel_unstable_far_count[2]
              << ','
              << result.tabi_convergence_channel_unstable_far_count[3]
              << ','
              << result.tabi_convergence_channel_unstable_far_panels[0]
              << ','
              << result.tabi_convergence_channel_unstable_far_panels[1]
              << ','
              << result.tabi_convergence_channel_unstable_far_panels[2]
              << ','
              << result.tabi_convergence_channel_unstable_far_panels[3]
              << ','
              << result.tabi_convergence_channel_split_rejects[0] << ','
              << result.tabi_convergence_channel_split_rejects[1] << ','
              << result.tabi_convergence_channel_split_rejects[2] << ','
              << result.tabi_convergence_channel_split_rejects[3] << ','
              << result.tabi_convergence_channel_split_reject_panels[0]
              << ','
              << result.tabi_convergence_channel_split_reject_panels[1]
              << ','
              << result.tabi_convergence_channel_split_reject_panels[2]
              << ','
              << result.tabi_convergence_channel_split_reject_panels[3]
              << ','
              << result.tabi_convergence_channel_low_abs_reaccepted_count[0]
              << ','
              << result.tabi_convergence_channel_low_abs_reaccepted_count[1]
              << ','
              << result.tabi_convergence_channel_low_abs_reaccepted_count[2]
              << ','
              << result.tabi_convergence_channel_low_abs_reaccepted_count[3]
              << ','
              << result.tabi_convergence_channel_low_abs_reaccepted_panels[0]
              << ','
              << result.tabi_convergence_channel_low_abs_reaccepted_panels[1]
              << ','
              << result.tabi_convergence_channel_low_abs_reaccepted_panels[2]
              << ','
              << result.tabi_convergence_channel_low_abs_reaccepted_panels[3]
              << ','
              << result.tabi_convergence_channel_high_error_caught[0] << ','
              << result.tabi_convergence_channel_high_error_caught[1] << ','
              << result.tabi_convergence_channel_high_error_caught[2] << ','
              << result.tabi_convergence_channel_high_error_caught[3] << ','
              << result.tabi_convergence_channel_high_error_missed[0] << ','
              << result.tabi_convergence_channel_high_error_missed[1] << ','
              << result.tabi_convergence_channel_high_error_missed[2] << ','
              << result.tabi_convergence_channel_high_error_missed[3] << ','
              << result.tabi_convergence_channel_high_error_caught_by_split[0]
              << ','
              << result.tabi_convergence_channel_high_error_caught_by_split[1]
              << ','
              << result.tabi_convergence_channel_high_error_caught_by_split[2]
              << ','
              << result.tabi_convergence_channel_high_error_caught_by_split[3]
              << ','
              << result.tabi_convergence_estimator_high_error_caught << ','
              << result.tabi_convergence_estimator_high_error_missed << ','
              << result.tabi_convergence_estimator_high_error_caught_by_split
              << ','
              << result.tabi_convergence_estimator_phantom_high_error_count
              << ','
              << result.tabi_convergence_estimator_phantom_high_error_abs_max
              << ','
              << result.tabi_convergence_overreject_count << ','
              << result.tabi_convergence_overreject_panels << ','
              << result.tabi_convergence_overreject_caught_count << ','
              << result.tabi_convergence_overreject_caught_panels << ','
               << result.tabi_convergence_overreject_abs_error_max << ','
               << result.tabi_convergence_overreject_abs_error_mean << ','
               << result.tabi_convergence_overreject_abs_error_p95 << ','
               << result.tabi_convergence_overreject_abs_error_p99 << ','
               << result.tabi_convergence_overreject_estimated_abs_error_max << ','
               << result.tabi_convergence_overreject_estimated_abs_error_mean << ','
               << result.tabi_convergence_overreject_required_scale_max << ','
               << result.tabi_convergence_overreject_required_scale_p95 << ','
               << result.tabi_convergence_overreject_required_scale_p99 << ','
               << result.tabi_convergence_channel_overreject_count[0] << ','
              << result.tabi_convergence_channel_overreject_count[1] << ','
              << result.tabi_convergence_channel_overreject_count[2] << ','
              << result.tabi_convergence_channel_overreject_count[3] << ','
              << result.tabi_convergence_channel_overreject_caught_count[0]
              << ','
              << result.tabi_convergence_channel_overreject_caught_count[1]
              << ','
               << result.tabi_convergence_channel_overreject_caught_count[2]
               << ','
               << result.tabi_convergence_channel_overreject_caught_count[3]
               << ','
               << result.tabi_convergence_channel_overreject_panels[0] << ','
               << result.tabi_convergence_channel_overreject_panels[1] << ','
               << result.tabi_convergence_channel_overreject_panels[2] << ','
               << result.tabi_convergence_channel_overreject_panels[3] << ','
               << result.tabi_convergence_channel_overreject_abs_error_max[0] << ','
               << result.tabi_convergence_channel_overreject_abs_error_max[1] << ','
               << result.tabi_convergence_channel_overreject_abs_error_max[2] << ','
               << result.tabi_convergence_channel_overreject_abs_error_max[3] << ','
               << result.tabi_convergence_channel_overreject_abs_error_mean[0] << ','
               << result.tabi_convergence_channel_overreject_abs_error_mean[1] << ','
               << result.tabi_convergence_channel_overreject_abs_error_mean[2] << ','
               << result.tabi_convergence_channel_overreject_abs_error_mean[3] << ','
               << result.tabi_convergence_channel_overreject_required_scale_max[0] << ','
               << result.tabi_convergence_channel_overreject_required_scale_max[1] << ','
               << result.tabi_convergence_channel_overreject_required_scale_max[2] << ','
               << result.tabi_convergence_channel_overreject_required_scale_max[3] << ','
               << result.tabi_convergence_channel_overreject_required_scale_p99[0] << ','
               << result.tabi_convergence_channel_overreject_required_scale_p99[1] << ','
               << result.tabi_convergence_channel_overreject_required_scale_p99[2] << ','
               << result.tabi_convergence_channel_overreject_required_scale_p99[3] << ','
              << result.tabi_phantom_high_error_channel_count[0] << ','
              << result.tabi_phantom_high_error_channel_count[1] << ','
              << result.tabi_phantom_high_error_channel_count[2] << ','
              << result.tabi_phantom_high_error_channel_count[3] << ','
              << result.tabi_phantom_high_error_channel_abs_max[0] << ','
              << result.tabi_phantom_high_error_channel_abs_max[1] << ','
              << result.tabi_phantom_high_error_channel_abs_max[2] << ','
              << result.tabi_phantom_high_error_channel_abs_max[3] << ','
              << csv_escape(result.tabi_accepted_miss_worst_channel) << ','
              << result.tabi_accepted_miss_channel_count[0] << ','
              << result.tabi_accepted_miss_channel_count[1] << ','
              << result.tabi_accepted_miss_channel_count[2] << ','
              << result.tabi_accepted_miss_channel_count[3] << ','
              << result.tabi_accepted_miss_channel_panels[0] << ','
              << result.tabi_accepted_miss_channel_panels[1] << ','
              << result.tabi_accepted_miss_channel_panels[2] << ','
              << result.tabi_accepted_miss_channel_panels[3] << ','
              << result.tabi_accepted_miss_channel_rho_min[0] << ','
              << result.tabi_accepted_miss_channel_rho_min[1] << ','
              << result.tabi_accepted_miss_channel_rho_min[2] << ','
              << result.tabi_accepted_miss_channel_rho_min[3] << ','
              << result.tabi_accepted_miss_channel_rho_mean[0] << ','
              << result.tabi_accepted_miss_channel_rho_mean[1] << ','
              << result.tabi_accepted_miss_channel_rho_mean[2] << ','
              << result.tabi_accepted_miss_channel_rho_mean[3] << ','
              << result.tabi_accepted_miss_channel_rho_max[0] << ','
              << result.tabi_accepted_miss_channel_rho_max[1] << ','
              << result.tabi_accepted_miss_channel_rho_max[2] << ','
              << result.tabi_accepted_miss_channel_rho_max[3] << ','
              << result.tabi_accepted_miss_channel_tail_ratio_mean[0] << ','
              << result.tabi_accepted_miss_channel_tail_ratio_mean[1] << ','
              << result.tabi_accepted_miss_channel_tail_ratio_mean[2] << ','
              << result.tabi_accepted_miss_channel_tail_ratio_mean[3] << ','
              << result.tabi_accepted_miss_channel_tail_ratio_max[0] << ','
              << result.tabi_accepted_miss_channel_tail_ratio_max[1] << ','
              << result.tabi_accepted_miss_channel_tail_ratio_max[2] << ','
              << result.tabi_accepted_miss_channel_tail_ratio_max[3] << ','
              << result.tabi_accepted_miss_channel_estimated_error_mean[0]
              << ','
              << result.tabi_accepted_miss_channel_estimated_error_mean[1]
              << ','
              << result.tabi_accepted_miss_channel_estimated_error_mean[2]
              << ','
              << result.tabi_accepted_miss_channel_estimated_error_mean[3]
              << ','
              << result.tabi_accepted_miss_channel_estimated_error_max[0]
              << ','
              << result.tabi_accepted_miss_channel_estimated_error_max[1]
              << ','
              << result.tabi_accepted_miss_channel_estimated_error_max[2]
              << ','
              << result.tabi_accepted_miss_channel_estimated_error_max[3]
              << ','
              << result.tabi_convergence_k3_rho_cap_rejects << ','
              << result.tabi_convergence_k3_rho_cap_reject_panels << ','
              << result.tabi_moment_audit_nodes << ','
              << result.tabi_moment_translation_rel_error_max << ','
              << result.tabi_moment_translation_rel_error_mean << ','
              << result.tabi_moment_translation_offset_max << ','
              << result.tabi_moment_translation_offset_mean << ','
              << result.tabi_moment_cancellation_ratio_max << ','
              << result.tabi_moment_cancellation_ratio_mean << ','
              << result.tabi_moment_dynamic_range_max << ','
              << csv_escape(result.tabi_expansion_audit) << ','
              << csv_escape(result.tabi_expansion_audit_channel) << ','
              << result.tabi_far_expansion_audit_pairs << ','
              << result.tabi_far_expansion_audit_panels << ','
              << seconds(result.tabi_far_expansion_direct_ms) << ','
              << result.tabi_far_expansion_rel_error_max << ','
              << result.tabi_far_expansion_rel_error_mean << ','
              << result.tabi_far_expansion_rel_error_p95 << ','
              << result.tabi_far_expansion_rel_error_p99 << ','
              << csv_escape(result.tabi_far_expansion_worst_channel) << ','
              << result.tabi_far_expansion_channel_rel_error_max[0] << ','
              << result.tabi_far_expansion_channel_rel_error_max[1] << ','
              << result.tabi_far_expansion_channel_rel_error_max[2] << ','
              << result.tabi_far_expansion_channel_rel_error_max[3] << ','
              << result.tabi_far_expansion_channel_rel_error_mean[0] << ','
              << result.tabi_far_expansion_channel_rel_error_mean[1] << ','
              << result.tabi_far_expansion_channel_rel_error_mean[2] << ','
              << result.tabi_far_expansion_channel_rel_error_mean[3] << ','
              << result.tabi_far_expansion_channel_rel_error_p99[0] << ','
              << result.tabi_far_expansion_channel_rel_error_p99[1] << ','
              << result.tabi_far_expansion_channel_rel_error_p99[2] << ','
              << result.tabi_far_expansion_channel_rel_error_p99[3] << ','
              << (result.tabi_far_expansion_mean_gate_pass ? 1 : 0) << ','
              << result.tabi_far_convergence_rho_max << ','
              << result.tabi_far_convergence_rho_mean << ','
              << result.tabi_accepted_far_aspect_ratio_max << ','
              << result.tabi_accepted_far_aspect_ratio_mean << ','
              << csv_escape(result.tabi_quad_buffer_auto_candidate) << ','
              << csv_escape(result.tabi_near_integration) << ','
              << csv_escape(result.tabi_quad_buffer) << ','
              << result.tabi_near_self_pairs << ','
              << result.tabi_near_one_ring_pairs << ','
              << result.tabi_near_duffy_pairs << ','
              << result.tabi_near_centroid_pairs << ','
              << result.tabi_near_fallbacks << ','
              << result.tabi_near_geometric_buffer_pairs << ','
              << result.tabi_near_geometric_buffer_gauss_pairs << ','
              << result.tabi_near_geometric_buffer_adaptive_pairs << ','
              << result.tabi_near_geometric_buffer_centroid_pairs << ','
              << result.tabi_near_geometric_buffer_min_gap << ','
              << result.tabi_near_geometric_buffer_mean_gap_ratio << ','
              << result.tabi_near_geometric_buffer_max_edge_mean << ','
              << seconds(result.tabi_near_integration_ms) << ','
              << result.tabi_duffy_subdivisions_max << ','
              << result.tabi_duffy_subdivisions_mean << ','
              << csv_escape(result.tabi_near_kernel_policy) << ','
              << result.tabi_near_k1_accepted << ','
              << result.tabi_near_k1_rejected << ','
              << result.tabi_near_k2_accepted << ','
              << result.tabi_near_k2_rejected << ','
              << result.tabi_near_k3_accepted << ','
              << result.tabi_near_k3_rejected << ','
              << result.tabi_near_k4_accepted << ','
              << result.tabi_near_k4_rejected << ','
              << result.tabi_near_one_ring_k1_accepted << ','
              << result.tabi_near_one_ring_k1_rejected << ','
              << result.tabi_near_one_ring_k2_accepted << ','
              << result.tabi_near_one_ring_k2_rejected << ','
              << result.tabi_near_one_ring_k3_accepted << ','
              << result.tabi_near_one_ring_k3_rejected << ','
              << result.tabi_near_one_ring_k4_accepted << ','
              << result.tabi_near_one_ring_k4_rejected << ','
              << csv_escape(result.tabi_near_k1_policy) << ','
              << csv_escape(result.tabi_near_k2_policy) << ','
              << csv_escape(result.tabi_near_k3_policy) << ','
              << csv_escape(result.tabi_near_k4_policy) << ','
              << csv_escape(result.tabi_surface_geometry) << ','
              << result.tabi_mesh_vertices << ','
              << result.tabi_mesh_triangles << ','
              << result.tabi_boundary_nodes << ','
              << result.tabi_avg_valence << ','
              << result.tabi_nonmanifold_edges << ','
              << result.tabi_solid_angle_fallbacks << ','
              << result.tabi_diag_min << ','
              << result.tabi_diag_max << ','
              << result.tabi_taylor_order << ','
              << csv_escape(result.tabi_derivative_mode) << ','
              << seconds(result.tabi_interaction_build_ms) << ','
              << seconds(result.tabi_tree_eval_ms) << ','
              << csv_escape(result.tabi_formulation) << ','
              << csv_escape(result.tabi_energy_projection) << ','
              << csv_escape(result.tabi_preconditioner_requested) << ','
              << csv_escape(result.tabi_preconditioner_effective) << ','
              << csv_escape(result.tabi_preconditioner_block_near) << ','
              << result.tabi_preconditioner_block_panel_cap << ','
              << result.tabi_preconditioner_overlap_panel_cap << ','
              << result.tabi_preconditioner_ilu_neighbor_cap << ','
              << result.tabi_preconditioner_ilu_drop_tol << ','
              << result.tabi_preconditioner_memory_budget_bytes << ','
              << result.tabi_preconditioner_estimated_bytes << ','
              << result.tabi_preconditioner_valid_blocks << ','
              << result.tabi_preconditioner_valid_panels << ','
              << result.tabi_preconditioner_ill_conditioned_blocks << ','
              << result.tabi_preconditioner_quality_rejected_blocks << ','
              << result.tabi_preconditioner_quality_probe_max_gain << ','
              << csv_escape(result.tabi_preconditioner_quality_status) << ','
              << seconds(result.tabi_preconditioner_setup_ms) << ','
              << seconds(result.tabi_preconditioner_apply_ms) << ','
              << result.tabi_preconditioner_apply_calls << ','
              << result.tabi_preconditioner_sparse_nnz_blocks << ','
              << result.tabi_preconditioner_sparse_avg_degree << ','
              << result.tabi_preconditioner_sparse_max_degree << ','
              << result.tabi_preconditioner_coarse_aggregate_cap << ','
              << result.tabi_preconditioner_coarse_aggregates << ','
              << result.tabi_preconditioner_coarse_dimension << ','
              << csv_escape(result.tabi_preconditioner_coarse_base) << ','
              << result.tabi_preconditioner_coarse_matrix_bytes << ','
              << result.tabi_preconditioner_coarse_pivot_ratio << ','
              << csv_escape(result.tabi_preconditioner_coarse_status) << ','
              << result.tabi_preconditioner_coarse_setup_matvec_calls << ','
              << seconds(result.tabi_preconditioner_coarse_setup_ms) << ','
              << seconds(result.tabi_preconditioner_coarse_apply_ms) << ','
              << result.tabi_preconditioner_coarse_apply_calls << ','
              << result.tabi_preconditioner_coarse_rp_identity_error_max << ','
              << result.tabi_preconditioner_coarse_verify_rel_error_max << ','
              << csv_escape(result.tabi_surface_generator) << ','
              << csv_escape(result.tabi_discretization) << ','
              << (result.tabi_fallback_used ? "true" : "false") << ','
              << result.tabi_rhs_scale << ','
              << result.tabi_dielectric_contrast << ','
              << result.tabi_operator_contrast << ','
              << result.tabi_unknown_count << ','
              << csv_escape(result.tabi_bie_system) << ','
              << result.tabi_lhs_diag_phi << ','
              << result.tabi_lhs_diag_phin << ','
              << csv_escape(result.tabi_energy_equation) << ','
              << csv_escape(result.tabi_kernel_convention) << ','
              << result.tabi_source_norm_s1 << ','
              << result.tabi_source_norm_s2 << ','
              << result.tabi_phi_norm << ','
              << result.tabi_phin_norm << ','
              << result.tabi_kappa_inv_a << ','
              << result.tabi_coulomb_factor_kcal_mol << ','
              << result.tabi_boundary_energy_kj_mol << ','
              << result.tabi_reference_energy_kj_mol << ','
              << result.tabi_selected_energy_scale << ','
              << seconds(timing.setup_ms)
              << ',' << seconds(timing.surface_ms) << ','
              << seconds(timing.epsilon_ms) << ','
              << seconds(timing.boundary_ms) << ','
              << seconds(timing.solver_ms) << ','
              << seconds(timing.energy_ms) << ','
              << seconds(timing.map_io_ms) << ','
              << seconds(timing.rpb_source_ms) << ','
              << seconds(timing.focus_setup_ms) << ','
              << seconds(timing.probe_ms) << ','
              << seconds(timing.focus_parent_solve_ms) << ','
              << seconds(timing.focus_child_setup_ms) << ','
              << seconds(timing.focus_boundary_interp_ms) << ','
              << seconds(timing.focus_child_surface_ms) << ','
              << seconds(timing.focus_child_epsilon_ms) << ','
              << seconds(timing.focus_child_solver_ms) << ','
              << seconds(timing.probe_energy_ms) << ','
              << seconds(timing.audit_output_ms) << ','
              << seconds(timing.nwt_half_update_ms) << ','
              << seconds(timing.nwt_rms_reduce_ms) << ','
              << seconds(timing.nwt_operator_residual_ms) << ','
              << seconds(timing.nwt_source_update_ms) << ','
              << seconds(timing.nwt_progress_io_ms) << ','
              << seconds(timing.nwt_coeff_pack_ms) << ','
              << seconds(timing.nwt_homo_update_ms) << ','
              << seconds(timing.nwt_hetero_update_ms) << ','
              << seconds(timing.nwt_nonlin_term_ms) << ','
              << seconds(timing.nwt_reduce_ms) << ','
              << seconds(timing.nwt_homo_stencil_ms) << ','
              << seconds(timing.nwt_homo_salt_math_ms) << ','
              << seconds(timing.nwt_homo_write_ms) << ','
              << seconds(timing.nwt_thread_wait_ms) << ','
              << seconds(timing.nwt_block_overhead_ms) << ','
              << seconds(timing.linear_half_update_ms) << ','
              << timing.nwt_homo_cells << ','
              << timing.nwt_homo_linear_cells << ','
              << timing.nwt_homo_salt_cells << ','
              << timing.nwt_homo_updates << ','
              << timing.nwt_dispatch_calls << ','
              << (timing.nwt_homo_updates > 0
                      ? (timing.nwt_homo_update_ms * 1.0e6) /
                            static_cast<double>(timing.nwt_homo_updates)
                      : 0.0)
              << ','
              << seconds(timing.epsilon_center_density_ms) << ','
              << seconds(timing.epsilon_midpoint_density_ms) << ','
              << seconds(timing.epsilon_ion_exclusion_ms) << ','
              << seconds(timing.epsilon_finalize_ms) << ','
              << seconds(timing.total_ms) << ',' << iterations << ','
              << result.nl_blocks << ',' << csv_escape(result.status) << '\n';
}

void AnalyzerImpl::write_tabi_profile_header_if_needed() {
  if (args_.timing != "detail")
    return;
  if (!tabi_profile_out_.is_open()) {
    tabi_profile_out_.open("pbe_tabi_profile.csv",
                           std::ios::out | std::ios::trunc);
    if (!tabi_profile_out_.is_open())
      throw std::runtime_error("PBE: failed to open pbe_tabi_profile.csv.");
  }
  if (!tabi_profile_header_written_) {
    tabi_profile_out_
        << "trajectory_id,frame_index,global_index,fixture,atom_count,"
           "surface_geometry,panels,mesh_vertices,mesh_triangles,"
           "surface_point_count_requested,panel_density_per_A2,"
           "panel_budget_status,"
           "matvec_mode,fmm_m2l_mode,fmm_m2l_expansion_order,"
           "fmm_m2l_near_pairs,"
           "fmm_m2l_far_pairs,fmm_m2l_setup_seconds,"
           "fmm_m2l_apply_seconds,fmm_m2l_finite_outputs,"
           "status,convergence_reason,reaction_energy_source,"
           "iterations,residual,residual_history,gmres_restart_requested,"
           "gmres_restart_effective,best_residual,best_residual_iter,"
           "final_residual,best_preconditioned_residual,"
           "final_preconditioned_residual,true_residual_abs,"
           "true_residual_rel,true_residual_probe_mode,"
           "true_residual_probe_scope,true_residual_probe_interval,"
           "true_residual_sample_count,true_residual_probe_calls,"
           "true_residual_probe_seconds,sampled_true_residual_rel,"
           "sampled_true_residual_rows,sampled_probe_calls,"
           "sampled_probe_seconds,full_true_residual_confirm_calls,"
           "full_probe_seconds,total_seconds,"
           "surface_build_seconds,surface_neighbor_seconds,"
           "surface_mesh_seconds,solid_angle_seconds,rhs_seconds,"
           "gmres_seconds,matvec_seconds_total,matvec_seconds_avg,"
           "matvec_calls,tree_vs_direct_matvec_rel,matvec_audit_rows,"
           "matvec_audit_seconds,orthogonalization_seconds,"
           "preconditioner_seconds,"
           "tree_build_seconds,tree_moment_seconds,tree_traversal_seconds,"
           "matvec_cost_profile,leaf_panels_per_second,"
           "clean_leaf_panels_per_second,tree_traversal_fraction,"
           "tree_eval_seconds,near_eval_seconds,far_eval_seconds,"
           "singular_near_eval_seconds,leaf_near_eval_seconds,"
           "leaf_panel_eval_seconds,leaf_tile_eval_seconds,"
           "leaf_residual_eval_seconds,smooth_leaf_node_eval_seconds,"
           "leaf_panel_kernel,leaf_tile_status,leaf_panel_pairs_per_second,"
           "leaf_near_policy,cluster_policy,target_batch_size,"
           "leaf_pair_min_source_panels,leaf_pair_min_targets,"
           "leaf_pair_selection_profile,leaf_pair_probe_cache,"
           "leaf_pair_target_group,leaf_pair_max_lookahead,"
           "leaf_pair_directness_probe,leaf_pair_directness_cache,"
           "leaf_pair_directness_verify,"
           "leaf_pair_scheduler,leaf_pair_worklist_max_active,"
           "leaf_pair_worklist_verify,leaf_block_verify,"
           "tree_preset,interaction_mode,node_bounds_mode,"
           "tree_builder,expansion_center,fragmentation_safety,"
           "mac_safety,mac_gap_factor,tree_depth_max,"
           "mac_reject_radius_limited,"
           "mac_accept_rate,node_radius_max,node_radius_mean,"
           "node_aspect_ratio_max,node_aspect_ratio_mean,"
           "leaf_panel_work_per_matvec,far_node_work_per_matvec,"
           "far_geometric_unsafe_count,far_geometric_unsafe_panels,"
           "far_geometric_min_gap,far_geometric_gap_ratio_min,"
           "far_geometric_gap_ratio_mean,far_geometric_audit_samples,"
           "mac_physical_gap_rejects,mac_convergence_rho_rejects,"
           "mac_convergence_aspect_rejects,"
           "mac_convergence_rho_reject_panels,"
           "mac_convergence_aspect_reject_panels,"
           "mac_surface_margin_rejects,"
           "mac_surface_margin_reject_panels,"
           "mac_surface_margin_applied,"
           "mac_surface_margin_skipped,"
           "mac_surface_margin_adaptive_rejects,"
           "mac_surface_margin_adaptive_reject_panels,"
           "mac_surface_margin_skipped_accepted_far,"
           "mac_surface_margin_hard_applied,"
           "mac_surface_margin_soft_applied,"
           "mac_surface_margin_hard_rejects,"
           "mac_surface_margin_hard_reject_panels,"
           "mac_surface_margin_soft_rejects,"
           "mac_surface_margin_soft_reject_panels,"
           "fragmented_far_count,"
           "fragmented_far_panels,fragmented_split_rejects,"
           "fragmentation_center_gap_ratio_max,"
           "fragmentation_center_gap_ratio_mean,"
           "moment_audit,moment_construction,moment_sum,"
           "convergence_audit,convergence_safety,convergence_rho_cap,"
           "mac_aspect_cap,convergence_estimator,"
           "convergence_tail_threshold,convergence_estimator_error_cap,"
           "convergence_unstable_rho_cap,convergence_calibration,"
           "convergence_dielectric_amp,"
           "convergence_estimator_model,"
           "convergence_channel_calibration,"
           "convergence_split_policy,"
           "convergence_split_selectivity,"
           "convergence_split_absolute_floor,"
           "convergence_overreject_audit,"
           "convergence_overreject_detail,"
           "convergence_shape_audit,"
           "convergence_kernel_response_audit,"
           "channel_missed_semantics,"
           "accepted_miss_diagnostics,"
           "accepted_miss_absolute_floor,"
           "convergence_k3_accepted_rho_cap,"
           "mac_distance_policy,mac_surface_margin_factor,"
           "mac_surface_margin_policy,"
           "mac_surface_margin_aspect_threshold,"
           "mac_surface_margin_rho_trigger,"
           "mac_surface_margin_soft_aspect_threshold,"
           "mac_surface_margin_soft_rho_trigger,"
           "mac_surface_margin_soft_factor,"
           "convergence_estimator_samples,"
           "convergence_estimator_unstable_far_count,"
           "convergence_estimator_unstable_far_panels,"
           "convergence_estimator_split_rejects,"
           "convergence_estimator_split_reject_panels,"
           "convergence_low_abs_reaccepted_far_count,"
           "convergence_low_abs_reaccepted_far_panels,"
           "convergence_estimator_tail_ratio_max,"
           "convergence_estimator_tail_ratio_mean,"
           "convergence_estimator_error_max,"
           "convergence_estimator_error_mean,"
           "convergence_estimator_abs_error_max,"
           "convergence_estimator_abs_error_mean,"
           "convergence_tail_norm_p_max,convergence_tail_norm_p_mean,"
           "convergence_tail_norm_pm1_max,convergence_tail_norm_pm1_mean,"
           "convergence_effective_jcp_amp_max,"
           "convergence_effective_jcp_amp_mean,"
           "convergence_dielectric_amp_max,"
           "convergence_dielectric_amp_mean,"
           "convergence_k1_effective_jcp_amp_max,"
           "convergence_k2_effective_jcp_amp_max,"
           "convergence_k3_effective_jcp_amp_max,"
           "convergence_k4_effective_jcp_amp_max,"
           "convergence_k1_effective_jcp_amp_mean,"
           "convergence_k2_effective_jcp_amp_mean,"
           "convergence_k3_effective_jcp_amp_mean,"
           "convergence_k4_effective_jcp_amp_mean,"
           "convergence_k1_estimated_error_max,"
           "convergence_k2_estimated_error_max,"
           "convergence_k3_estimated_error_max,"
           "convergence_k4_estimated_error_max,"
           "convergence_k1_unstable_far_count,"
           "convergence_k2_unstable_far_count,"
           "convergence_k3_unstable_far_count,"
           "convergence_k4_unstable_far_count,"
           "convergence_k1_unstable_far_panels,"
           "convergence_k2_unstable_far_panels,"
           "convergence_k3_unstable_far_panels,"
           "convergence_k4_unstable_far_panels,"
           "convergence_k1_split_rejects,"
           "convergence_k2_split_rejects,"
           "convergence_k3_split_rejects,"
           "convergence_k4_split_rejects,"
           "convergence_k1_split_reject_panels,"
           "convergence_k2_split_reject_panels,"
           "convergence_k3_split_reject_panels,"
           "convergence_k4_split_reject_panels,"
           "convergence_k1_low_abs_reaccepted,"
           "convergence_k2_low_abs_reaccepted,"
           "convergence_k3_low_abs_reaccepted,"
           "convergence_k4_low_abs_reaccepted,"
           "convergence_k1_low_abs_reaccepted_panels,"
           "convergence_k2_low_abs_reaccepted_panels,"
           "convergence_k3_low_abs_reaccepted_panels,"
           "convergence_k4_low_abs_reaccepted_panels,"
           "convergence_k1_high_error_caught,"
           "convergence_k2_high_error_caught,"
           "convergence_k3_high_error_caught,"
           "convergence_k4_high_error_caught,"
           "convergence_k1_high_error_missed,"
           "convergence_k2_high_error_missed,"
           "convergence_k3_high_error_missed,"
           "convergence_k4_high_error_missed,"
           "convergence_k1_high_error_caught_by_split,"
           "convergence_k2_high_error_caught_by_split,"
           "convergence_k3_high_error_caught_by_split,"
           "convergence_k4_high_error_caught_by_split,"
           "convergence_estimator_high_error_caught,"
           "convergence_estimator_high_error_missed,"
           "convergence_estimator_high_error_caught_by_split,"
           "convergence_estimator_phantom_high_error_count,"
           "convergence_estimator_phantom_high_error_abs_max,"
           "convergence_overreject_count,"
           "convergence_overreject_panels,"
           "convergence_overreject_caught_count,"
           "convergence_overreject_caught_panels,"
           "convergence_overreject_abs_error_max,"
           "convergence_overreject_abs_error_mean,"
           "convergence_overreject_abs_error_p95,"
           "convergence_overreject_abs_error_p99,"
           "convergence_overreject_estimated_abs_error_max,"
           "convergence_overreject_estimated_abs_error_mean,"
           "convergence_overreject_required_scale_max,"
           "convergence_overreject_required_scale_p95,"
           "convergence_overreject_required_scale_p99,"
           "convergence_k1_overreject_count,"
           "convergence_k2_overreject_count,"
           "convergence_k3_overreject_count,"
           "convergence_k4_overreject_count,"
           "convergence_k1_overreject_caught_count,"
           "convergence_k2_overreject_caught_count,"
           "convergence_k3_overreject_caught_count,"
           "convergence_k4_overreject_caught_count,"
           "convergence_k1_overreject_panels,convergence_k2_overreject_panels,"
           "convergence_k3_overreject_panels,convergence_k4_overreject_panels,"
           "convergence_k1_overreject_abs_error_max,"
           "convergence_k2_overreject_abs_error_max,"
           "convergence_k3_overreject_abs_error_max,"
           "convergence_k4_overreject_abs_error_max,"
           "convergence_k1_overreject_abs_error_mean,"
           "convergence_k2_overreject_abs_error_mean,"
           "convergence_k3_overreject_abs_error_mean,"
           "convergence_k4_overreject_abs_error_mean,"
           "convergence_k1_overreject_required_scale_max,"
           "convergence_k2_overreject_required_scale_max,"
           "convergence_k3_overreject_required_scale_max,"
           "convergence_k4_overreject_required_scale_max,"
           "convergence_k1_overreject_required_scale_p99,"
           "convergence_k2_overreject_required_scale_p99,"
           "convergence_k3_overreject_required_scale_p99,"
           "convergence_k4_overreject_required_scale_p99,"
           "phantom_high_error_k1_count,phantom_high_error_k2_count,"
           "phantom_high_error_k3_count,phantom_high_error_k4_count,"
           "phantom_high_error_k1_abs_max,"
           "phantom_high_error_k2_abs_max,"
           "phantom_high_error_k3_abs_max,"
           "phantom_high_error_k4_abs_max,"
           "accepted_miss_worst_channel,"
           "accepted_miss_k1_count,accepted_miss_k2_count,"
           "accepted_miss_k3_count,accepted_miss_k4_count,"
           "accepted_miss_k1_panels,accepted_miss_k2_panels,"
           "accepted_miss_k3_panels,accepted_miss_k4_panels,"
           "accepted_miss_k1_rho_min,accepted_miss_k2_rho_min,"
           "accepted_miss_k3_rho_min,accepted_miss_k4_rho_min,"
           "accepted_miss_k1_rho_mean,accepted_miss_k2_rho_mean,"
           "accepted_miss_k3_rho_mean,accepted_miss_k4_rho_mean,"
           "accepted_miss_k1_rho_max,accepted_miss_k2_rho_max,"
           "accepted_miss_k3_rho_max,accepted_miss_k4_rho_max,"
           "accepted_miss_k1_tail_ratio_mean,"
           "accepted_miss_k2_tail_ratio_mean,"
           "accepted_miss_k3_tail_ratio_mean,"
           "accepted_miss_k4_tail_ratio_mean,"
           "accepted_miss_k1_tail_ratio_max,"
           "accepted_miss_k2_tail_ratio_max,"
           "accepted_miss_k3_tail_ratio_max,"
           "accepted_miss_k4_tail_ratio_max,"
           "accepted_miss_k1_estimated_error_mean,"
           "accepted_miss_k2_estimated_error_mean,"
           "accepted_miss_k3_estimated_error_mean,"
           "accepted_miss_k4_estimated_error_mean,"
           "accepted_miss_k1_estimated_error_max,"
           "accepted_miss_k2_estimated_error_max,"
           "accepted_miss_k3_estimated_error_max,"
           "accepted_miss_k4_estimated_error_max,"
           "convergence_k3_rho_cap_rejects,"
           "convergence_k3_rho_cap_reject_panels,"
           "moment_audit_nodes,"
           "moment_translation_rel_error_max,"
           "moment_translation_rel_error_mean,"
           "moment_translation_offset_max,moment_translation_offset_mean,"
           "moment_cancellation_ratio_max,moment_cancellation_ratio_mean,"
           "moment_dynamic_range_max,"
           "expansion_audit,expansion_audit_channel,"
           "far_expansion_audit_pairs,"
           "far_expansion_audit_panels,far_expansion_direct_seconds,"
           "far_expansion_rel_error_max,far_expansion_rel_error_mean,"
           "far_expansion_rel_error_p95,far_expansion_rel_error_p99,"
           "far_expansion_worst_channel,"
           "far_expansion_k1_rel_error_max,"
           "far_expansion_k2_rel_error_max,"
           "far_expansion_k3_rel_error_max,"
           "far_expansion_k4_rel_error_max,"
           "far_expansion_k1_rel_error_mean,"
           "far_expansion_k2_rel_error_mean,"
           "far_expansion_k3_rel_error_mean,"
           "far_expansion_k4_rel_error_mean,"
           "far_expansion_k1_rel_error_p99,"
           "far_expansion_k2_rel_error_p99,"
           "far_expansion_k3_rel_error_p99,"
           "far_expansion_k4_rel_error_p99,"
           "far_expansion_mean_gate_pass,"
           "far_convergence_rho_max,far_convergence_rho_mean,"
           "accepted_far_aspect_ratio_max,accepted_far_aspect_ratio_mean,"
           "quad_buffer_auto_candidate,"
           "near_integration,quad_buffer,near_self_pairs,near_one_ring_pairs,"
           "near_duffy_pairs,near_centroid_pairs,near_fallbacks,"
           "near_geometric_buffer_pairs,near_geometric_buffer_gauss_pairs,"
           "near_geometric_buffer_adaptive_pairs,"
           "near_geometric_buffer_centroid_pairs,"
           "near_geometric_buffer_min_gap,"
           "near_geometric_buffer_mean_gap_ratio,"
           "near_geometric_buffer_max_edge_mean,"
           "near_integration_seconds,duffy_subdivisions_max,"
           "duffy_subdivisions_mean,near_kernel_policy,"
           "near_k1_accepted,near_k1_rejected,near_k2_accepted,"
           "near_k2_rejected,near_k3_accepted,near_k3_rejected,"
           "near_k4_accepted,near_k4_rejected,"
           "near_one_ring_k1_accepted,near_one_ring_k1_rejected,"
           "near_one_ring_k2_accepted,near_one_ring_k2_rejected,"
           "near_one_ring_k3_accepted,near_one_ring_k3_rejected,"
           "near_one_ring_k4_accepted,near_one_ring_k4_rejected,"
           "near_k1_policy,"
           "near_k2_policy,near_k3_policy,near_k4_policy,"
           "near_interactions,far_interactions,singular_near_interactions,"
           "leaf_near_interactions,far_node_interactions,"
           "leaf_panel_interactions,leaf_tile_interactions,"
           "leaf_tile_panel_interactions,leaf_residual_panel_interactions,"
           "leaf_work_profile,leaf_work_mac_terminal_panels,"
           "leaf_work_physical_gap_panels,"
           "leaf_work_fragmentation_panels,"
           "leaf_work_convergence_panels,"
           "leaf_work_convergence_estimator_panels,"
           "leaf_work_convergence_rho_panels,"
           "leaf_work_convergence_aspect_panels,"
           "leaf_work_surface_margin_panels,"
           "leaf_work_surface_margin_hard_panels,"
           "leaf_work_surface_margin_soft_panels,"
           "leaf_work_k3_rho_cap_panels,"
           "leaf_work_singular_fallback_panels,"
           "leaf_work_scalar_tail_panels,"
           "leaf_work_profile_accounted_panels,"
           "leaf_work_mac_terminal_frac,"
           "leaf_work_convergence_frac,"
           "leaf_work_surface_margin_frac,"
           "leaf_work_singular_fallback_frac,"
           "leaf_work_scalar_tail_frac,"
           "smooth_leaf_node_interactions,"
           "leaf_pruned_interactions,leaf_mac_candidates,"
           "leaf_mac_accepted,leaf_mac_rejected,leaf_mac_far_panels,"
           "leaf_mac_direct_panels_saved_estimate,"
           "adaptive_smooth_leaf_accepted,"
           "adaptive_smooth_leaf_rejected,adaptive_span_tile_accepted,"
           "adaptive_span_tile_rejected,adaptive_exact_leaf_fallbacks,"
           "adaptive_reject_one_ring_overlap,"
           "adaptive_reject_error_bound_failed,"
           "adaptive_reject_p4_check_failed,adaptive_reject_too_close,"
           "leaf_work_reduction_ratio,target_cluster_interactions,"
           "target_cluster_targets,target_cluster_panels_pruned,"
           "target_cluster_rejected,target_cluster_eval_seconds,"
           "leaf_pair_batch_interactions,leaf_pair_batch_targets,"
           "leaf_pair_batch_source_panels,"
           "leaf_pair_batch_residual_panels,"
           "leaf_pair_batch_fallbacks,leaf_pair_batch_seconds,"
           "leaf_pair_batch_kernel,leaf_pair_batch_vectorized_sources,"
           "leaf_pair_batch_vector_lanes,"
           "leaf_pair_batch_pairs_per_second,"
           "leaf_kernel_profile,jcp_direct_kernel,"
           "leaf_kernel_scalar_tasks,leaf_kernel_batch_tasks,"
           "leaf_kernel_source_panel_evals,"
           "leaf_kernel_scalar_seconds,"
           "leaf_kernel_target_avx2_seconds,"
           "leaf_kernel_tail_seconds,"
           "leaf_pair_batch_soa_sources,"
           "leaf_pair_batch_soa_seconds,"
           "leaf_pair_batch_scalar_tail_seconds,"
           "leaf_pair_batch_clean_source_panels_per_second,"
           "leaf_pair_selection_candidate_scans,"
           "leaf_pair_selection_clean_rejects,"
           "leaf_pair_selection_probe_calls,"
           "leaf_pair_selection_probe_node_visits,"
           "leaf_pair_selection_probe_seconds,"
           "leaf_pair_selection_seconds,"
           "leaf_pair_selection_probe_cache_hits,"
           "leaf_pair_selection_probe_cache_misses,"
           "leaf_pair_selection_insufficient_targets,"
           "leaf_pair_selection_same_target_leaf_rejects,"
           "leaf_pair_directness_ancestor_probe_calls,"
           "leaf_pair_directness_ancestor_probe_node_visits,"
           "leaf_pair_directness_ancestor_probe_seconds,"
           "leaf_pair_directness_pair_cache_hits,"
           "leaf_pair_directness_pair_cache_misses,"
           "leaf_pair_directness_all_direct_proofs,"
           "leaf_pair_directness_mixed_fallbacks,"
           "leaf_pair_directness_verify_checks,"
           "leaf_pair_directness_verify_mismatches,"
           "leaf_pair_worklist_tasks,leaf_pair_worklist_groups,"
           "leaf_pair_worklist_batched_groups,"
           "leaf_pair_worklist_scalar_flush_tasks,"
           "leaf_pair_worklist_forced_flushes,"
           "leaf_pair_worklist_active_groups_max,"
           "leaf_pair_worklist_seconds,"
           "leaf_pair_worklist_verify_checks,"
           "leaf_pair_worklist_verify_mismatches,"
           "leaf_block_target_leaves,leaf_block_targets,"
           "leaf_block_source_groups,leaf_block_batched_groups,"
           "leaf_block_scalar_tail_tasks,leaf_block_group_seconds,"
           "leaf_block_flush_seconds,leaf_block_writeback_seconds,"
           "leaf_block_verify_checks,leaf_block_verify_mismatches,"
           "adaptive_level,"
           "tree_nodes,tree_leaf_count,"
            "taylor_order,derivative_mode,diag_min,diag_max,solid_angle_fallbacks,"
           "nonmanifold_edges,formulation,energy_projection,fallback_used,"
           "preconditioner_requested,preconditioner_effective,"
           "preconditioner_block_near,"
           "preconditioner_block_panel_cap,preconditioner_overlap_panel_cap,"
           "preconditioner_ilu_neighbor_cap,preconditioner_ilu_drop_tol,"
           "preconditioner_memory_budget_bytes,"
           "preconditioner_estimated_bytes,preconditioner_valid_blocks,"
           "preconditioner_valid_panels,preconditioner_ill_conditioned_blocks,"
           "preconditioner_quality_rejected_blocks,"
           "preconditioner_quality_probe_max_gain,"
           "preconditioner_quality_status,"
           "preconditioner_setup_seconds,preconditioner_apply_seconds,"
           "preconditioner_apply_calls,"
           "preconditioner_sparse_nnz_blocks,"
           "preconditioner_sparse_avg_degree,"
           "preconditioner_sparse_max_degree,"
           "preconditioner_coarse_aggregate_cap,preconditioner_coarse_aggregates,"
           "preconditioner_coarse_dimension,preconditioner_coarse_base,"
           "preconditioner_coarse_matrix_bytes,preconditioner_coarse_pivot_ratio,"
           "preconditioner_coarse_status,preconditioner_coarse_setup_matvec_calls,"
           "preconditioner_coarse_setup_seconds,preconditioner_coarse_apply_seconds,"
           "preconditioner_coarse_apply_calls,"
           "preconditioner_coarse_rp_identity_error_max,"
           "preconditioner_coarse_verify_rel_error_max,"
           "surface_generator,discretization,"
           "rhs_scale,dielectric_contrast,operator_contrast,unknown_count,"
           "bie_system,lhs_diag_phi,lhs_diag_phin,energy_equation,"
            "kernel_convention,source_norm_s1,source_norm_s2,phi_norm,"
           "phin_norm,kappa_inv_A,"
           "boundary_energy_kJmol,reference_energy_kJmol,"
           "selected_energy_scale,memory_bytes,energy_kJmol\n";
    tabi_profile_header_written_ = true;
  }
}

void AnalyzerImpl::append_tabi_profile_row(const FrameResult &result,
                                           const TimingBuckets &timing) {
  if (args_.timing != "detail" || result.pbe_backend != "tabi")
    return;
  write_tabi_profile_header_if_needed();
  auto seconds = [](double ms) { return ms / 1000.0; };
  tabi_profile_out_
      << result.trajectory_id << ',' << result.frame_index << ','
      << result.global_index << ',' << csv_escape(current_source_label_) << ','
      << result.atom_count << ',' << csv_escape(result.tabi_surface_geometry)
      << ',' << result.tabi_surface_panels << ',' << result.tabi_mesh_vertices
      << ',' << result.tabi_mesh_triangles << ','
      << result.tabi_surface_point_count_requested << ','
      << result.tabi_panel_density_per_a2 << ','
      << csv_escape(result.tabi_panel_budget_status) << ','
      << csv_escape(result.tabi_matvec_mode) << ','
      << csv_escape(result.tabi_fmm_m2l_mode) << ','
      << result.tabi_fmm_m2l_expansion_order << ','
      << result.tabi_fmm_m2l_near_pairs << ','
      << result.tabi_fmm_m2l_far_pairs << ','
      << seconds(result.tabi_fmm_m2l_setup_ms) << ','
      << seconds(result.tabi_fmm_m2l_apply_ms) << ','
      << (result.tabi_fmm_m2l_finite_outputs ? "true" : "false") << ','
      << csv_escape(result.tabi_status) << ','
      << csv_escape(result.tabi_convergence_reason) << ','
      << csv_escape(result.tabi_reaction_energy_source) << ','
      << result.solvent_iterations << ',' << result.linear_residual << ','
      << csv_escape(result.tabi_residual_history) << ','
      << result.tabi_gmres_restart_requested << ','
      << result.tabi_gmres_restart_effective << ','
      << result.tabi_best_residual << ','
      << result.tabi_best_residual_iter << ','
      << result.tabi_final_residual << ','
      << result.tabi_best_preconditioned_residual << ','
      << result.tabi_final_preconditioned_residual << ','
      << result.tabi_true_residual_abs << ','
      << result.tabi_true_residual_rel << ','
      << csv_escape(result.tabi_true_residual_probe_mode) << ','
      << csv_escape(result.tabi_true_residual_probe_scope) << ','
      << result.tabi_true_residual_probe_interval << ','
      << result.tabi_true_residual_sample_count << ','
      << result.tabi_true_residual_probe_calls << ','
      << seconds(result.tabi_true_residual_probe_ms) << ','
      << result.tabi_sampled_true_residual_rel << ','
      << result.tabi_sampled_true_residual_rows << ','
      << result.tabi_sampled_probe_calls << ','
      << seconds(result.tabi_sampled_probe_ms) << ','
      << result.tabi_full_true_residual_confirm_calls << ','
      << seconds(result.tabi_full_probe_ms) << ','
      << seconds(timing.total_ms) << ','
      << seconds(result.tabi_surface_build_ms) << ','
      << seconds(result.tabi_surface_neighbor_ms) << ','
      << seconds(result.tabi_surface_mesh_ms) << ','
      << seconds(result.tabi_solid_angle_ms) << ','
      << seconds(result.tabi_rhs_ms) << ','
      << seconds(result.tabi_gmres_ms) << ','
      << seconds(result.tabi_matvec_ms_total) << ','
      << seconds(result.tabi_matvec_ms_avg) << ','
      << result.tabi_matvec_calls << ','
      << result.tabi_tree_vs_direct_matvec_rel << ','
      << result.tabi_matvec_audit_rows << ','
      << seconds(result.tabi_matvec_audit_ms) << ','
      << seconds(result.tabi_orthogonalization_ms) << ','
      << seconds(result.tabi_preconditioner_ms) << ','
      << seconds(result.tabi_tree_build_ms) << ','
      << seconds(result.tabi_tree_moment_ms) << ','
      << seconds(result.tabi_tree_traversal_ms) << ','
      << csv_escape(result.tabi_matvec_cost_profile) << ','
      << result.tabi_leaf_panels_per_second << ','
      << result.tabi_clean_leaf_panels_per_second << ','
      << result.tabi_tree_traversal_fraction << ','
      << seconds(result.tabi_tree_eval_ms) << ','
      << seconds(result.tabi_near_eval_ms) << ','
      << seconds(result.tabi_far_eval_ms) << ','
      << seconds(result.tabi_singular_near_eval_ms) << ','
      << seconds(result.tabi_leaf_near_eval_ms) << ','
      << seconds(result.tabi_leaf_panel_eval_ms) << ','
      << seconds(result.tabi_leaf_tile_eval_ms) << ','
      << seconds(result.tabi_leaf_residual_eval_ms) << ','
      << seconds(result.tabi_smooth_leaf_node_eval_ms) << ','
      << csv_escape(result.tabi_leaf_panel_kernel) << ','
      << csv_escape(result.tabi_leaf_tile_status) << ','
      << result.tabi_leaf_panel_pairs_per_second << ','
      << csv_escape(result.tabi_leaf_near_policy) << ','
      << csv_escape(result.tabi_cluster_policy) << ','
      << result.tabi_target_batch_size << ','
      << result.tabi_leaf_pair_min_source_panels << ','
      << result.tabi_leaf_pair_min_targets << ','
      << csv_escape(result.tabi_leaf_pair_selection_profile) << ','
      << csv_escape(result.tabi_leaf_pair_probe_cache) << ','
      << csv_escape(result.tabi_leaf_pair_target_group) << ','
      << result.tabi_leaf_pair_max_lookahead << ','
      << csv_escape(result.tabi_leaf_pair_directness_probe) << ','
      << csv_escape(result.tabi_leaf_pair_directness_cache) << ','
      << csv_escape(result.tabi_leaf_pair_directness_verify) << ','
      << csv_escape(result.tabi_leaf_pair_scheduler) << ','
      << result.tabi_leaf_pair_worklist_max_active << ','
      << csv_escape(result.tabi_leaf_pair_worklist_verify) << ','
      << csv_escape(result.tabi_leaf_block_verify) << ','
      << csv_escape(result.tabi_tree_preset) << ','
      << csv_escape(result.tabi_interaction_mode) << ','
      << csv_escape(result.tabi_node_bounds_mode) << ','
      << csv_escape(result.tabi_tree_builder) << ','
      << csv_escape(result.tabi_expansion_center) << ','
      << csv_escape(result.tabi_fragmentation_safety) << ','
      << csv_escape(result.tabi_mac_safety) << ','
      << result.tabi_mac_gap_factor << ','
      << result.tabi_tree_depth_max << ','
      << result.tabi_mac_reject_radius_limited << ','
      << result.tabi_mac_accept_rate << ','
      << result.tabi_node_radius_max << ','
      << result.tabi_node_radius_mean << ','
      << result.tabi_node_aspect_ratio_max << ','
      << result.tabi_node_aspect_ratio_mean << ','
      << result.tabi_leaf_panel_work_per_matvec << ','
      << result.tabi_far_node_work_per_matvec << ','
      << result.tabi_far_geometric_unsafe_count << ','
      << result.tabi_far_geometric_unsafe_panels << ','
      << result.tabi_far_geometric_min_gap << ','
      << result.tabi_far_geometric_gap_ratio_min << ','
      << result.tabi_far_geometric_gap_ratio_mean << ','
      << result.tabi_far_geometric_audit_samples << ','
      << result.tabi_mac_physical_gap_rejects << ','
      << result.tabi_mac_convergence_rho_rejects << ','
      << result.tabi_mac_convergence_aspect_rejects << ','
      << result.tabi_mac_convergence_rho_reject_panels << ','
      << result.tabi_mac_convergence_aspect_reject_panels << ','
      << result.tabi_mac_surface_margin_rejects << ','
      << result.tabi_mac_surface_margin_reject_panels << ','
      << result.tabi_mac_surface_margin_applied << ','
      << result.tabi_mac_surface_margin_skipped << ','
      << result.tabi_mac_surface_margin_adaptive_rejects << ','
      << result.tabi_mac_surface_margin_adaptive_reject_panels << ','
      << result.tabi_mac_surface_margin_skipped_accepted_far << ','
      << result.tabi_mac_surface_margin_hard_applied << ','
      << result.tabi_mac_surface_margin_soft_applied << ','
      << result.tabi_mac_surface_margin_hard_rejects << ','
      << result.tabi_mac_surface_margin_hard_reject_panels << ','
      << result.tabi_mac_surface_margin_soft_rejects << ','
      << result.tabi_mac_surface_margin_soft_reject_panels << ','
      << result.tabi_fragmented_far_count << ','
      << result.tabi_fragmented_far_panels << ','
      << result.tabi_fragmented_split_rejects << ','
      << result.tabi_fragmentation_center_gap_ratio_max << ','
      << result.tabi_fragmentation_center_gap_ratio_mean << ','
      << csv_escape(result.tabi_moment_audit) << ','
      << csv_escape(result.tabi_moment_construction) << ','
      << csv_escape(result.tabi_moment_sum) << ','
      << csv_escape(result.tabi_convergence_audit) << ','
      << csv_escape(result.tabi_convergence_safety) << ','
      << result.tabi_convergence_rho_cap << ','
      << result.tabi_mac_aspect_cap << ','
      << csv_escape(result.tabi_convergence_estimator) << ','
      << result.tabi_convergence_tail_threshold << ','
      << result.tabi_convergence_estimator_error_cap << ','
      << result.tabi_convergence_unstable_rho_cap << ','
      << csv_escape(result.tabi_convergence_calibration) << ','
      << csv_escape(result.tabi_convergence_dielectric_amp) << ','
      << csv_escape(result.tabi_convergence_estimator_model) << ','
      << csv_escape(result.tabi_convergence_channel_calibration) << ','
      << csv_escape(result.tabi_convergence_split_policy) << ','
      << csv_escape(result.tabi_convergence_split_selectivity) << ','
      << result.tabi_convergence_split_absolute_floor << ','
      << csv_escape(result.tabi_convergence_overreject_audit) << ','
      << csv_escape(result.tabi_convergence_overreject_detail) << ','
      << csv_escape(result.tabi_convergence_shape_audit) << ','
      << csv_escape(result.tabi_convergence_kernel_response_audit) << ','
      << csv_escape(result.tabi_channel_missed_semantics) << ','
      << csv_escape(result.tabi_accepted_miss_diagnostics) << ','
      << result.tabi_accepted_miss_absolute_floor << ','
      << result.tabi_convergence_k3_accepted_rho_cap << ','
      << csv_escape(result.tabi_mac_distance_policy) << ','
      << result.tabi_mac_surface_margin_factor << ','
      << csv_escape(result.tabi_mac_surface_margin_policy) << ','
      << result.tabi_mac_surface_margin_aspect_threshold << ','
      << result.tabi_mac_surface_margin_rho_trigger << ','
      << result.tabi_mac_surface_margin_soft_aspect_threshold << ','
      << result.tabi_mac_surface_margin_soft_rho_trigger << ','
      << result.tabi_mac_surface_margin_soft_factor << ','
      << result.tabi_convergence_estimator_samples << ','
      << result.tabi_convergence_estimator_unstable_far_count << ','
      << result.tabi_convergence_estimator_unstable_far_panels << ','
      << result.tabi_convergence_estimator_split_rejects << ','
      << result.tabi_convergence_estimator_split_reject_panels << ','
      << result.tabi_convergence_low_abs_reaccepted_far_count << ','
      << result.tabi_convergence_low_abs_reaccepted_far_panels << ','
      << result.tabi_convergence_estimator_tail_ratio_max << ','
      << result.tabi_convergence_estimator_tail_ratio_mean << ','
      << result.tabi_convergence_estimator_error_max << ','
      << result.tabi_convergence_estimator_error_mean << ','
      << result.tabi_convergence_estimator_abs_error_max << ','
      << result.tabi_convergence_estimator_abs_error_mean << ','
      << result.tabi_convergence_tail_norm_p_max << ','
      << result.tabi_convergence_tail_norm_p_mean << ','
      << result.tabi_convergence_tail_norm_pm1_max << ','
      << result.tabi_convergence_tail_norm_pm1_mean << ','
      << result.tabi_convergence_effective_jcp_amp_max << ','
      << result.tabi_convergence_effective_jcp_amp_mean << ','
      << result.tabi_convergence_dielectric_amp_max << ','
      << result.tabi_convergence_dielectric_amp_mean << ','
      << result.tabi_convergence_channel_effective_jcp_amp_max[0] << ','
      << result.tabi_convergence_channel_effective_jcp_amp_max[1] << ','
      << result.tabi_convergence_channel_effective_jcp_amp_max[2] << ','
      << result.tabi_convergence_channel_effective_jcp_amp_max[3] << ','
      << result.tabi_convergence_channel_effective_jcp_amp_mean[0] << ','
      << result.tabi_convergence_channel_effective_jcp_amp_mean[1] << ','
      << result.tabi_convergence_channel_effective_jcp_amp_mean[2] << ','
      << result.tabi_convergence_channel_effective_jcp_amp_mean[3] << ','
      << result.tabi_convergence_channel_estimated_error_max[0] << ','
      << result.tabi_convergence_channel_estimated_error_max[1] << ','
      << result.tabi_convergence_channel_estimated_error_max[2] << ','
      << result.tabi_convergence_channel_estimated_error_max[3] << ','
      << result.tabi_convergence_channel_unstable_far_count[0] << ','
      << result.tabi_convergence_channel_unstable_far_count[1] << ','
      << result.tabi_convergence_channel_unstable_far_count[2] << ','
      << result.tabi_convergence_channel_unstable_far_count[3] << ','
      << result.tabi_convergence_channel_unstable_far_panels[0] << ','
      << result.tabi_convergence_channel_unstable_far_panels[1] << ','
      << result.tabi_convergence_channel_unstable_far_panels[2] << ','
      << result.tabi_convergence_channel_unstable_far_panels[3] << ','
      << result.tabi_convergence_channel_split_rejects[0] << ','
      << result.tabi_convergence_channel_split_rejects[1] << ','
      << result.tabi_convergence_channel_split_rejects[2] << ','
      << result.tabi_convergence_channel_split_rejects[3] << ','
      << result.tabi_convergence_channel_split_reject_panels[0] << ','
      << result.tabi_convergence_channel_split_reject_panels[1] << ','
      << result.tabi_convergence_channel_split_reject_panels[2] << ','
      << result.tabi_convergence_channel_split_reject_panels[3] << ','
      << result.tabi_convergence_channel_low_abs_reaccepted_count[0] << ','
      << result.tabi_convergence_channel_low_abs_reaccepted_count[1] << ','
      << result.tabi_convergence_channel_low_abs_reaccepted_count[2] << ','
      << result.tabi_convergence_channel_low_abs_reaccepted_count[3] << ','
      << result.tabi_convergence_channel_low_abs_reaccepted_panels[0] << ','
      << result.tabi_convergence_channel_low_abs_reaccepted_panels[1] << ','
      << result.tabi_convergence_channel_low_abs_reaccepted_panels[2] << ','
      << result.tabi_convergence_channel_low_abs_reaccepted_panels[3] << ','
      << result.tabi_convergence_channel_high_error_caught[0] << ','
      << result.tabi_convergence_channel_high_error_caught[1] << ','
      << result.tabi_convergence_channel_high_error_caught[2] << ','
      << result.tabi_convergence_channel_high_error_caught[3] << ','
      << result.tabi_convergence_channel_high_error_missed[0] << ','
      << result.tabi_convergence_channel_high_error_missed[1] << ','
      << result.tabi_convergence_channel_high_error_missed[2] << ','
      << result.tabi_convergence_channel_high_error_missed[3] << ','
      << result.tabi_convergence_channel_high_error_caught_by_split[0] << ','
      << result.tabi_convergence_channel_high_error_caught_by_split[1] << ','
      << result.tabi_convergence_channel_high_error_caught_by_split[2] << ','
      << result.tabi_convergence_channel_high_error_caught_by_split[3] << ','
      << result.tabi_convergence_estimator_high_error_caught << ','
      << result.tabi_convergence_estimator_high_error_missed << ','
      << result.tabi_convergence_estimator_high_error_caught_by_split << ','
      << result.tabi_convergence_estimator_phantom_high_error_count << ','
      << result.tabi_convergence_estimator_phantom_high_error_abs_max << ','
      << result.tabi_convergence_overreject_count << ','
      << result.tabi_convergence_overreject_panels << ','
      << result.tabi_convergence_overreject_caught_count << ','
      << result.tabi_convergence_overreject_caught_panels << ','
      << result.tabi_convergence_overreject_abs_error_max << ','
      << result.tabi_convergence_overreject_abs_error_mean << ','
      << result.tabi_convergence_overreject_abs_error_p95 << ','
      << result.tabi_convergence_overreject_abs_error_p99 << ','
      << result.tabi_convergence_overreject_estimated_abs_error_max << ','
      << result.tabi_convergence_overreject_estimated_abs_error_mean << ','
      << result.tabi_convergence_overreject_required_scale_max << ','
      << result.tabi_convergence_overreject_required_scale_p95 << ','
      << result.tabi_convergence_overreject_required_scale_p99 << ','
      << result.tabi_convergence_channel_overreject_count[0] << ','
      << result.tabi_convergence_channel_overreject_count[1] << ','
      << result.tabi_convergence_channel_overreject_count[2] << ','
      << result.tabi_convergence_channel_overreject_count[3] << ','
      << result.tabi_convergence_channel_overreject_caught_count[0] << ','
      << result.tabi_convergence_channel_overreject_caught_count[1] << ','
      << result.tabi_convergence_channel_overreject_caught_count[2] << ','
      << result.tabi_convergence_channel_overreject_caught_count[3] << ','
      << result.tabi_convergence_channel_overreject_panels[0] << ','
      << result.tabi_convergence_channel_overreject_panels[1] << ','
      << result.tabi_convergence_channel_overreject_panels[2] << ','
      << result.tabi_convergence_channel_overreject_panels[3] << ','
      << result.tabi_convergence_channel_overreject_abs_error_max[0] << ','
      << result.tabi_convergence_channel_overreject_abs_error_max[1] << ','
      << result.tabi_convergence_channel_overreject_abs_error_max[2] << ','
      << result.tabi_convergence_channel_overreject_abs_error_max[3] << ','
      << result.tabi_convergence_channel_overreject_abs_error_mean[0] << ','
      << result.tabi_convergence_channel_overreject_abs_error_mean[1] << ','
      << result.tabi_convergence_channel_overreject_abs_error_mean[2] << ','
      << result.tabi_convergence_channel_overreject_abs_error_mean[3] << ','
      << result.tabi_convergence_channel_overreject_required_scale_max[0] << ','
      << result.tabi_convergence_channel_overreject_required_scale_max[1] << ','
      << result.tabi_convergence_channel_overreject_required_scale_max[2] << ','
      << result.tabi_convergence_channel_overreject_required_scale_max[3] << ','
      << result.tabi_convergence_channel_overreject_required_scale_p99[0] << ','
      << result.tabi_convergence_channel_overreject_required_scale_p99[1] << ','
      << result.tabi_convergence_channel_overreject_required_scale_p99[2] << ','
      << result.tabi_convergence_channel_overreject_required_scale_p99[3] << ','
      << result.tabi_phantom_high_error_channel_count[0] << ','
      << result.tabi_phantom_high_error_channel_count[1] << ','
      << result.tabi_phantom_high_error_channel_count[2] << ','
      << result.tabi_phantom_high_error_channel_count[3] << ','
      << result.tabi_phantom_high_error_channel_abs_max[0] << ','
      << result.tabi_phantom_high_error_channel_abs_max[1] << ','
      << result.tabi_phantom_high_error_channel_abs_max[2] << ','
      << result.tabi_phantom_high_error_channel_abs_max[3] << ','
      << csv_escape(result.tabi_accepted_miss_worst_channel) << ','
      << result.tabi_accepted_miss_channel_count[0] << ','
      << result.tabi_accepted_miss_channel_count[1] << ','
      << result.tabi_accepted_miss_channel_count[2] << ','
      << result.tabi_accepted_miss_channel_count[3] << ','
      << result.tabi_accepted_miss_channel_panels[0] << ','
      << result.tabi_accepted_miss_channel_panels[1] << ','
      << result.tabi_accepted_miss_channel_panels[2] << ','
      << result.tabi_accepted_miss_channel_panels[3] << ','
      << result.tabi_accepted_miss_channel_rho_min[0] << ','
      << result.tabi_accepted_miss_channel_rho_min[1] << ','
      << result.tabi_accepted_miss_channel_rho_min[2] << ','
      << result.tabi_accepted_miss_channel_rho_min[3] << ','
      << result.tabi_accepted_miss_channel_rho_mean[0] << ','
      << result.tabi_accepted_miss_channel_rho_mean[1] << ','
      << result.tabi_accepted_miss_channel_rho_mean[2] << ','
      << result.tabi_accepted_miss_channel_rho_mean[3] << ','
      << result.tabi_accepted_miss_channel_rho_max[0] << ','
      << result.tabi_accepted_miss_channel_rho_max[1] << ','
      << result.tabi_accepted_miss_channel_rho_max[2] << ','
      << result.tabi_accepted_miss_channel_rho_max[3] << ','
      << result.tabi_accepted_miss_channel_tail_ratio_mean[0] << ','
      << result.tabi_accepted_miss_channel_tail_ratio_mean[1] << ','
      << result.tabi_accepted_miss_channel_tail_ratio_mean[2] << ','
      << result.tabi_accepted_miss_channel_tail_ratio_mean[3] << ','
      << result.tabi_accepted_miss_channel_tail_ratio_max[0] << ','
      << result.tabi_accepted_miss_channel_tail_ratio_max[1] << ','
      << result.tabi_accepted_miss_channel_tail_ratio_max[2] << ','
      << result.tabi_accepted_miss_channel_tail_ratio_max[3] << ','
      << result.tabi_accepted_miss_channel_estimated_error_mean[0] << ','
      << result.tabi_accepted_miss_channel_estimated_error_mean[1] << ','
      << result.tabi_accepted_miss_channel_estimated_error_mean[2] << ','
      << result.tabi_accepted_miss_channel_estimated_error_mean[3] << ','
      << result.tabi_accepted_miss_channel_estimated_error_max[0] << ','
      << result.tabi_accepted_miss_channel_estimated_error_max[1] << ','
      << result.tabi_accepted_miss_channel_estimated_error_max[2] << ','
      << result.tabi_accepted_miss_channel_estimated_error_max[3] << ','
      << result.tabi_convergence_k3_rho_cap_rejects << ','
      << result.tabi_convergence_k3_rho_cap_reject_panels << ','
      << result.tabi_moment_audit_nodes << ','
      << result.tabi_moment_translation_rel_error_max << ','
      << result.tabi_moment_translation_rel_error_mean << ','
      << result.tabi_moment_translation_offset_max << ','
      << result.tabi_moment_translation_offset_mean << ','
      << result.tabi_moment_cancellation_ratio_max << ','
      << result.tabi_moment_cancellation_ratio_mean << ','
      << result.tabi_moment_dynamic_range_max << ','
      << csv_escape(result.tabi_expansion_audit) << ','
      << csv_escape(result.tabi_expansion_audit_channel) << ','
      << result.tabi_far_expansion_audit_pairs << ','
      << result.tabi_far_expansion_audit_panels << ','
      << seconds(result.tabi_far_expansion_direct_ms) << ','
      << result.tabi_far_expansion_rel_error_max << ','
      << result.tabi_far_expansion_rel_error_mean << ','
      << result.tabi_far_expansion_rel_error_p95 << ','
      << result.tabi_far_expansion_rel_error_p99 << ','
      << csv_escape(result.tabi_far_expansion_worst_channel) << ','
      << result.tabi_far_expansion_channel_rel_error_max[0] << ','
      << result.tabi_far_expansion_channel_rel_error_max[1] << ','
      << result.tabi_far_expansion_channel_rel_error_max[2] << ','
      << result.tabi_far_expansion_channel_rel_error_max[3] << ','
      << result.tabi_far_expansion_channel_rel_error_mean[0] << ','
      << result.tabi_far_expansion_channel_rel_error_mean[1] << ','
      << result.tabi_far_expansion_channel_rel_error_mean[2] << ','
      << result.tabi_far_expansion_channel_rel_error_mean[3] << ','
      << result.tabi_far_expansion_channel_rel_error_p99[0] << ','
      << result.tabi_far_expansion_channel_rel_error_p99[1] << ','
      << result.tabi_far_expansion_channel_rel_error_p99[2] << ','
      << result.tabi_far_expansion_channel_rel_error_p99[3] << ','
      << (result.tabi_far_expansion_mean_gate_pass ? 1 : 0) << ','
      << result.tabi_far_convergence_rho_max << ','
      << result.tabi_far_convergence_rho_mean << ','
      << result.tabi_accepted_far_aspect_ratio_max << ','
      << result.tabi_accepted_far_aspect_ratio_mean << ','
      << csv_escape(result.tabi_quad_buffer_auto_candidate) << ','
      << csv_escape(result.tabi_near_integration) << ','
      << csv_escape(result.tabi_quad_buffer) << ','
      << result.tabi_near_self_pairs << ','
      << result.tabi_near_one_ring_pairs << ','
      << result.tabi_near_duffy_pairs << ','
      << result.tabi_near_centroid_pairs << ','
      << result.tabi_near_fallbacks << ','
      << result.tabi_near_geometric_buffer_pairs << ','
      << result.tabi_near_geometric_buffer_gauss_pairs << ','
      << result.tabi_near_geometric_buffer_adaptive_pairs << ','
      << result.tabi_near_geometric_buffer_centroid_pairs << ','
      << result.tabi_near_geometric_buffer_min_gap << ','
      << result.tabi_near_geometric_buffer_mean_gap_ratio << ','
      << result.tabi_near_geometric_buffer_max_edge_mean << ','
      << seconds(result.tabi_near_integration_ms) << ','
      << result.tabi_duffy_subdivisions_max << ','
      << result.tabi_duffy_subdivisions_mean << ','
      << csv_escape(result.tabi_near_kernel_policy) << ','
      << result.tabi_near_k1_accepted << ','
      << result.tabi_near_k1_rejected << ','
      << result.tabi_near_k2_accepted << ','
      << result.tabi_near_k2_rejected << ','
      << result.tabi_near_k3_accepted << ','
      << result.tabi_near_k3_rejected << ','
      << result.tabi_near_k4_accepted << ','
      << result.tabi_near_k4_rejected << ','
      << result.tabi_near_one_ring_k1_accepted << ','
      << result.tabi_near_one_ring_k1_rejected << ','
      << result.tabi_near_one_ring_k2_accepted << ','
      << result.tabi_near_one_ring_k2_rejected << ','
      << result.tabi_near_one_ring_k3_accepted << ','
      << result.tabi_near_one_ring_k3_rejected << ','
      << result.tabi_near_one_ring_k4_accepted << ','
      << result.tabi_near_one_ring_k4_rejected << ','
      << csv_escape(result.tabi_near_k1_policy) << ','
      << csv_escape(result.tabi_near_k2_policy) << ','
      << csv_escape(result.tabi_near_k3_policy) << ','
      << csv_escape(result.tabi_near_k4_policy) << ','
      << result.tabi_near_interactions << ','
      << result.tabi_far_interactions << ','
      << result.tabi_singular_near_interactions << ','
      << result.tabi_leaf_near_interactions << ','
      << result.tabi_far_node_interactions << ','
      << result.tabi_leaf_panel_interactions << ','
      << result.tabi_leaf_tile_interactions << ','
      << result.tabi_leaf_tile_panel_interactions << ','
      << result.tabi_leaf_residual_panel_interactions << ','
      << csv_escape(result.tabi_leaf_work_profile) << ','
      << result.tabi_leaf_work_mac_terminal_panels << ','
      << result.tabi_leaf_work_physical_gap_panels << ','
      << result.tabi_leaf_work_fragmentation_panels << ','
      << result.tabi_leaf_work_convergence_panels << ','
      << result.tabi_leaf_work_convergence_estimator_panels << ','
      << result.tabi_leaf_work_convergence_rho_panels << ','
      << result.tabi_leaf_work_convergence_aspect_panels << ','
      << result.tabi_leaf_work_surface_margin_panels << ','
      << result.tabi_leaf_work_surface_margin_hard_panels << ','
      << result.tabi_leaf_work_surface_margin_soft_panels << ','
      << result.tabi_leaf_work_k3_rho_cap_panels << ','
      << result.tabi_leaf_work_singular_fallback_panels << ','
      << result.tabi_leaf_work_scalar_tail_panels << ','
      << result.tabi_leaf_work_profile_accounted_panels << ','
      << result.tabi_leaf_work_mac_terminal_frac << ','
      << result.tabi_leaf_work_convergence_frac << ','
      << result.tabi_leaf_work_surface_margin_frac << ','
      << result.tabi_leaf_work_singular_fallback_frac << ','
      << result.tabi_leaf_work_scalar_tail_frac << ','
      << result.tabi_smooth_leaf_node_interactions << ','
      << result.tabi_leaf_pruned_interactions << ','
      << result.tabi_leaf_mac_candidates << ','
      << result.tabi_leaf_mac_accepted << ','
      << result.tabi_leaf_mac_rejected << ','
      << result.tabi_leaf_mac_far_panels << ','
      << result.tabi_leaf_mac_direct_panels_saved_estimate << ','
      << result.tabi_adaptive_smooth_leaf_accepted << ','
      << result.tabi_adaptive_smooth_leaf_rejected << ','
      << result.tabi_adaptive_span_tile_accepted << ','
      << result.tabi_adaptive_span_tile_rejected << ','
      << result.tabi_adaptive_exact_leaf_fallbacks << ','
      << result.tabi_adaptive_reject_one_ring_overlap << ','
      << result.tabi_adaptive_reject_error_bound_failed << ','
      << result.tabi_adaptive_reject_p4_check_failed << ','
      << result.tabi_adaptive_reject_too_close << ','
      << result.tabi_leaf_work_reduction_ratio << ','
      << result.tabi_target_cluster_interactions << ','
      << result.tabi_target_cluster_targets << ','
      << result.tabi_target_cluster_panels_pruned << ','
      << result.tabi_target_cluster_rejected << ','
      << seconds(result.tabi_target_cluster_eval_ms) << ','
      << result.tabi_leaf_pair_batch_interactions << ','
      << result.tabi_leaf_pair_batch_targets << ','
      << result.tabi_leaf_pair_batch_source_panels << ','
      << result.tabi_leaf_pair_batch_residual_panels << ','
      << result.tabi_leaf_pair_batch_fallbacks << ','
      << seconds(result.tabi_leaf_pair_batch_ms) << ','
      << csv_escape(result.tabi_leaf_pair_batch_kernel) << ','
      << result.tabi_leaf_pair_batch_vectorized_sources << ','
      << result.tabi_leaf_pair_batch_vector_lanes << ','
      << result.tabi_leaf_pair_batch_pairs_per_second << ','
      << csv_escape(result.tabi_leaf_kernel_profile) << ','
      << csv_escape(result.tabi_jcp_direct_kernel) << ','
      << result.tabi_leaf_kernel_scalar_tasks << ','
      << result.tabi_leaf_kernel_batch_tasks << ','
      << result.tabi_leaf_kernel_source_panel_evals << ','
      << seconds(result.tabi_leaf_kernel_scalar_ms) << ','
      << seconds(result.tabi_leaf_kernel_target_avx2_ms) << ','
      << seconds(result.tabi_leaf_kernel_tail_ms) << ','
      << result.tabi_leaf_pair_batch_soa_sources << ','
      << seconds(result.tabi_leaf_pair_batch_soa_ms) << ','
      << seconds(result.tabi_leaf_pair_batch_scalar_tail_ms) << ','
      << result.tabi_leaf_pair_batch_clean_source_panels_per_second << ','
      << result.tabi_leaf_pair_selection_candidate_scans << ','
      << result.tabi_leaf_pair_selection_clean_rejects << ','
      << result.tabi_leaf_pair_selection_probe_calls << ','
      << result.tabi_leaf_pair_selection_probe_node_visits << ','
      << seconds(result.tabi_leaf_pair_selection_probe_ms) << ','
      << seconds(result.tabi_leaf_pair_selection_ms) << ','
      << result.tabi_leaf_pair_selection_probe_cache_hits << ','
      << result.tabi_leaf_pair_selection_probe_cache_misses << ','
      << result.tabi_leaf_pair_selection_insufficient_targets << ','
      << result.tabi_leaf_pair_selection_same_target_leaf_rejects << ','
      << result.tabi_leaf_pair_directness_ancestor_probe_calls << ','
      << result.tabi_leaf_pair_directness_ancestor_probe_node_visits << ','
      << seconds(result.tabi_leaf_pair_directness_ancestor_probe_ms) << ','
      << result.tabi_leaf_pair_directness_pair_cache_hits << ','
      << result.tabi_leaf_pair_directness_pair_cache_misses << ','
      << result.tabi_leaf_pair_directness_all_direct_proofs << ','
      << result.tabi_leaf_pair_directness_mixed_fallbacks << ','
      << result.tabi_leaf_pair_directness_verify_checks << ','
      << result.tabi_leaf_pair_directness_verify_mismatches << ','
      << result.tabi_leaf_pair_worklist_tasks << ','
      << result.tabi_leaf_pair_worklist_groups << ','
      << result.tabi_leaf_pair_worklist_batched_groups << ','
      << result.tabi_leaf_pair_worklist_scalar_flush_tasks << ','
      << result.tabi_leaf_pair_worklist_forced_flushes << ','
      << result.tabi_leaf_pair_worklist_active_groups_max << ','
      << seconds(result.tabi_leaf_pair_worklist_ms) << ','
      << result.tabi_leaf_pair_worklist_verify_checks << ','
      << result.tabi_leaf_pair_worklist_verify_mismatches << ','
      << result.tabi_leaf_block_target_leaves << ','
      << result.tabi_leaf_block_targets << ','
      << result.tabi_leaf_block_source_groups << ','
      << result.tabi_leaf_block_batched_groups << ','
      << result.tabi_leaf_block_scalar_tail_tasks << ','
      << seconds(result.tabi_leaf_block_group_ms) << ','
      << seconds(result.tabi_leaf_block_flush_ms) << ','
      << seconds(result.tabi_leaf_block_writeback_ms) << ','
      << result.tabi_leaf_block_verify_checks << ','
      << result.tabi_leaf_block_verify_mismatches << ','
      << csv_escape(result.tabi_adaptive_level) << ','
      << result.tabi_tree_nodes << ','
      << result.tabi_tree_leaf_count << ',' << result.tabi_taylor_order << ','
      << csv_escape(result.tabi_derivative_mode) << ','
      << result.tabi_diag_min << ',' << result.tabi_diag_max << ','
      << result.tabi_solid_angle_fallbacks << ','
      << result.tabi_nonmanifold_edges << ','
      << csv_escape(result.tabi_formulation) << ','
      << csv_escape(result.tabi_energy_projection) << ','
      << (result.tabi_fallback_used ? "true" : "false") << ','
      << csv_escape(result.tabi_preconditioner_requested) << ','
      << csv_escape(result.tabi_preconditioner_effective) << ','
      << csv_escape(result.tabi_preconditioner_block_near) << ','
      << result.tabi_preconditioner_block_panel_cap << ','
      << result.tabi_preconditioner_overlap_panel_cap << ','
      << result.tabi_preconditioner_ilu_neighbor_cap << ','
      << result.tabi_preconditioner_ilu_drop_tol << ','
      << result.tabi_preconditioner_memory_budget_bytes << ','
      << result.tabi_preconditioner_estimated_bytes << ','
      << result.tabi_preconditioner_valid_blocks << ','
      << result.tabi_preconditioner_valid_panels << ','
      << result.tabi_preconditioner_ill_conditioned_blocks << ','
      << result.tabi_preconditioner_quality_rejected_blocks << ','
      << result.tabi_preconditioner_quality_probe_max_gain << ','
      << csv_escape(result.tabi_preconditioner_quality_status) << ','
      << seconds(result.tabi_preconditioner_setup_ms) << ','
      << seconds(result.tabi_preconditioner_apply_ms) << ','
      << result.tabi_preconditioner_apply_calls << ','
      << result.tabi_preconditioner_sparse_nnz_blocks << ','
      << result.tabi_preconditioner_sparse_avg_degree << ','
      << result.tabi_preconditioner_sparse_max_degree << ','
      << result.tabi_preconditioner_coarse_aggregate_cap << ','
      << result.tabi_preconditioner_coarse_aggregates << ','
      << result.tabi_preconditioner_coarse_dimension << ','
      << csv_escape(result.tabi_preconditioner_coarse_base) << ','
      << result.tabi_preconditioner_coarse_matrix_bytes << ','
      << result.tabi_preconditioner_coarse_pivot_ratio << ','
      << csv_escape(result.tabi_preconditioner_coarse_status) << ','
      << result.tabi_preconditioner_coarse_setup_matvec_calls << ','
      << seconds(result.tabi_preconditioner_coarse_setup_ms) << ','
      << seconds(result.tabi_preconditioner_coarse_apply_ms) << ','
      << result.tabi_preconditioner_coarse_apply_calls << ','
      << result.tabi_preconditioner_coarse_rp_identity_error_max << ','
      << result.tabi_preconditioner_coarse_verify_rel_error_max << ','
      << csv_escape(result.tabi_surface_generator) << ','
      << csv_escape(result.tabi_discretization) << ','
      << result.tabi_rhs_scale << ',' << result.tabi_dielectric_contrast << ','
      << result.tabi_operator_contrast << ',' << result.tabi_unknown_count
      << ',' << csv_escape(result.tabi_bie_system) << ','
      << result.tabi_lhs_diag_phi << ',' << result.tabi_lhs_diag_phin << ','
      << csv_escape(result.tabi_energy_equation) << ','
      << csv_escape(result.tabi_kernel_convention) << ','
      << result.tabi_source_norm_s1 << ',' << result.tabi_source_norm_s2
      << ',' << result.tabi_phi_norm << ',' << result.tabi_phin_norm << ','
      << result.tabi_kappa_inv_a << ','
      << result.tabi_boundary_energy_kj_mol << ','
      << result.tabi_reference_energy_kj_mol << ','
      << result.tabi_selected_energy_scale << ',' << result.memory_bytes << ','
      << result.energy_kj_mol << '\n';
}

void AnalyzerImpl::write_tabi_near_block_audit_header_if_needed() {
  if (args_.timing != "detail")
    return;
  if (!tabi_near_block_audit_out_.is_open()) {
    tabi_near_block_audit_out_.open("pbe_tabi_near_block_audit.csv",
                                    std::ios::out | std::ios::trunc);
    if (!tabi_near_block_audit_out_.is_open()) {
      throw std::runtime_error(
          "PBE: failed to open pbe_tabi_near_block_audit.csv.");
    }
  }
  if (!tabi_near_block_audit_header_written_) {
    tabi_near_block_audit_out_
        << "trajectory_id,frame_index,global_index,fixture,block,pair_type,"
           "orientation,accepted,rejected,candidate_worse_than_centroid,"
           "reference_nonconverged,orientation_sign_mismatch,"
           "strict_ref_reject,"
           "max_candidate_ref_rel,mean_candidate_ref_rel,"
           "max_centroid_ref_rel,k2_laplace_sign_ok,"
           "k2_yukawa_smooth_ref_rel,k2_orientation_pair_ref_rel\n";
    tabi_near_block_audit_header_written_ = true;
  }
}

void AnalyzerImpl::append_tabi_near_block_audit_rows(
    const FrameResult &result) {
  if (args_.timing != "detail" || result.pbe_backend != "tabi")
    return;
  write_tabi_near_block_audit_header_if_needed();
  static constexpr std::array<const char *, kNearBlockCount> kBlockLabels = {
      "K1", "K2", "K3", "K4"};
  static constexpr std::array<const char *, kNearPairTypeCount>
      kPairTypeLabels = {"self", "edge", "vertex", "smooth_near"};
  static constexpr std::array<const char *, kNearOrientationCount>
      kOrientationLabels = {"same", "flipped", "mixed"};
  for (std::size_t block = 0; block < kNearBlockCount; ++block) {
    for (std::size_t pair_type = 0; pair_type < kNearPairTypeCount;
         ++pair_type) {
      for (std::size_t orientation = 0; orientation < kNearOrientationCount;
           ++orientation) {
        const std::size_t index =
            (block * kNearPairTypeCount + pair_type) * kNearOrientationCount +
            orientation;
        if (index >= result.tabi_near_block_audit.size())
          continue;
        const auto &bucket = result.tabi_near_block_audit[index];
        const double mean_candidate_ref_rel =
            bucket.sample_count > 0
                ? bucket.sum_candidate_ref_rel /
                      static_cast<double>(bucket.sample_count)
                : 0.0;
        tabi_near_block_audit_out_
            << result.trajectory_id << ',' << result.frame_index << ','
            << result.global_index << ',' << csv_escape(current_source_label_)
            << ',' << kBlockLabels[block] << ','
            << kPairTypeLabels[pair_type] << ','
            << kOrientationLabels[orientation] << ',' << bucket.accepted
            << ',' << bucket.rejected << ','
            << bucket.candidate_worse_than_centroid << ','
            << bucket.reference_nonconverged << ','
            << bucket.orientation_sign_mismatch << ','
            << bucket.strict_ref_reject << ','
            << bucket.max_candidate_ref_rel << ',' << mean_candidate_ref_rel
            << ',' << bucket.max_centroid_ref_rel << ','
            << bucket.k2_laplace_sign_ok << ','
            << bucket.max_k2_yukawa_smooth_ref_rel << ','
            << bucket.max_k2_orientation_pair_ref_rel << '\n';
      }
    }
  }
}

void AnalyzerImpl::write_tabi_geometric_near_audit_header_if_needed() {
  if (args_.timing != "detail")
    return;
  if (!tabi_geometric_near_audit_out_.is_open()) {
    tabi_geometric_near_audit_out_.open("pbe_tabi_geometric_near_audit.csv",
                                        std::ios::out | std::ios::trunc);
    if (!tabi_geometric_near_audit_out_.is_open()) {
      throw std::runtime_error(
          "PBE: failed to open pbe_tabi_geometric_near_audit.csv.");
    }
  }
  if (!tabi_geometric_near_audit_header_written_) {
    tabi_geometric_near_audit_out_
        << "trajectory_id,frame_index,global_index,fixture,gap_bucket,"
           "center_bucket,block,kernel,pairs,min_gap,max_gap_ratio,"
           "mean_gap_ratio,mean_center_ratio\n";
    tabi_geometric_near_audit_header_written_ = true;
  }
}

void AnalyzerImpl::append_tabi_geometric_near_audit_rows(
    const FrameResult &result) {
  if (args_.timing != "detail" || result.pbe_backend != "tabi")
    return;
  write_tabi_geometric_near_audit_header_if_needed();
  static constexpr std::array<const char *, kGeometricGapBucketCount>
      kGapLabels = {"lt_0p25", "0p25_0p50", "0p50_0p75", "gt_0p75"};
  static constexpr std::array<const char *, kGeometricCenterBucketCount>
      kCenterLabels = {"lt_1p0", "1p0_2p0", "ge_2p0"};
  static constexpr std::array<const char *, kGeometricKernelCount>
      kKernelLabels = {"centroid", "gauss7", "subdivided"};
  static constexpr std::array<const char *, kNearBlockCount> kBlockLabels = {
      "K1", "K2", "K3", "K4"};
  for (std::size_t gap = 0; gap < kGeometricGapBucketCount; ++gap) {
    for (std::size_t center = 0; center < kGeometricCenterBucketCount;
         ++center) {
      for (std::size_t kernel = 0; kernel < kGeometricKernelCount; ++kernel) {
        const std::size_t index =
            (gap * kGeometricCenterBucketCount + center) *
                kGeometricKernelCount +
            kernel;
        if (index >= result.tabi_geometric_near_audit.size())
          continue;
        const auto &bucket = result.tabi_geometric_near_audit[index];
        if (bucket.pairs == 0)
          continue;
        const double denom = static_cast<double>(bucket.pairs);
        const double mean_gap_ratio = bucket.sum_gap_ratio / denom;
        const double mean_center_ratio = bucket.sum_center_ratio / denom;
        for (const char *block_label : kBlockLabels) {
          tabi_geometric_near_audit_out_
              << result.trajectory_id << ',' << result.frame_index << ','
              << result.global_index << ','
              << csv_escape(current_source_label_) << ','
              << kGapLabels[gap] << ',' << kCenterLabels[center] << ','
              << block_label << ',' << kKernelLabels[kernel] << ','
              << bucket.pairs << ',' << bucket.min_gap << ','
              << bucket.max_gap_ratio << ',' << mean_gap_ratio << ','
              << mean_center_ratio << '\n';
        }
      }
    }
  }
}

void AnalyzerImpl::write_tabi_far_geometric_audit_header_if_needed() {
  if (args_.timing != "detail")
    return;
  if (!tabi_far_geometric_audit_out_.is_open()) {
    tabi_far_geometric_audit_out_.open("pbe_tabi_far_geometric_audit.csv",
                                       std::ios::out | std::ios::trunc);
    if (!tabi_far_geometric_audit_out_.is_open()) {
      throw std::runtime_error(
          "PBE: failed to open pbe_tabi_far_geometric_audit.csv.");
    }
  }
  if (!tabi_far_geometric_audit_header_written_) {
    tabi_far_geometric_audit_out_
        << "trajectory_id,frame_index,global_index,fixture,tree_builder,"
           "node_bounds,theta,leaf,node_depth,gap_ratio,node_panel_count,"
           "mac_safety,mac_gap_factor,unsafe_count,unsafe_panels,min_gap,"
           "min_gap_ratio,mean_gap_ratio,audit_samples,"
           "physical_gap_rejects,convergence_safety,convergence_rho_cap,"
           "mac_aspect_cap,convergence_estimator,"
           "convergence_tail_threshold,convergence_estimator_error_cap,"
           "convergence_unstable_rho_cap,convergence_calibration,"
           "convergence_estimator_samples,"
           "convergence_estimator_unstable_far_count,"
           "convergence_estimator_split_rejects,"
           "convergence_estimator_tail_ratio_max,"
           "convergence_estimator_error_max,convergence_tail_norm_p_max,"
           "convergence_tail_norm_pm1_max,"
           "convergence_effective_jcp_amp_max,"
           "far_expansion_mean_gate_pass,convergence_rho_rejects,"
           "convergence_aspect_rejects,convergence_rho_reject_panels,"
           "convergence_aspect_reject_panels\n";
    tabi_far_geometric_audit_header_written_ = true;
  }
}

void AnalyzerImpl::append_tabi_far_geometric_audit_row(
    const FrameResult &result) {
  if (args_.timing != "detail" || result.pbe_backend != "tabi")
    return;
  write_tabi_far_geometric_audit_header_if_needed();
  tabi_far_geometric_audit_out_
      << result.trajectory_id << ',' << result.frame_index << ','
      << result.global_index << ',' << csv_escape(current_source_label_) << ','
      << csv_escape(result.tabi_tree_builder) << ','
      << csv_escape(result.tabi_node_bounds_mode) << ','
      << args_.tabi_tree_opening_theta << ','
      << args_.tabi_tree_leaf_size << ','
      << "all,all,all," << csv_escape(result.tabi_mac_safety) << ','
      << result.tabi_mac_gap_factor << ','
      << result.tabi_far_geometric_unsafe_count << ','
      << result.tabi_far_geometric_unsafe_panels << ','
      << result.tabi_far_geometric_min_gap << ','
      << result.tabi_far_geometric_gap_ratio_min << ','
      << result.tabi_far_geometric_gap_ratio_mean << ','
      << result.tabi_far_geometric_audit_samples << ','
      << result.tabi_mac_physical_gap_rejects << ','
      << csv_escape(result.tabi_convergence_safety) << ','
      << result.tabi_convergence_rho_cap << ','
      << result.tabi_mac_aspect_cap << ','
      << csv_escape(result.tabi_convergence_estimator) << ','
      << result.tabi_convergence_tail_threshold << ','
      << result.tabi_convergence_estimator_error_cap << ','
      << result.tabi_convergence_unstable_rho_cap << ','
      << csv_escape(result.tabi_convergence_calibration) << ','
      << result.tabi_convergence_estimator_samples << ','
      << result.tabi_convergence_estimator_unstable_far_count << ','
      << result.tabi_convergence_estimator_split_rejects << ','
      << result.tabi_convergence_estimator_tail_ratio_max << ','
      << result.tabi_convergence_estimator_error_max << ','
      << result.tabi_convergence_tail_norm_p_max << ','
      << result.tabi_convergence_tail_norm_pm1_max << ','
      << result.tabi_convergence_effective_jcp_amp_max << ','
      << (result.tabi_far_expansion_mean_gate_pass ? 1 : 0) << ','
      << result.tabi_mac_convergence_rho_rejects << ','
      << result.tabi_mac_convergence_aspect_rejects << ','
      << result.tabi_mac_convergence_rho_reject_panels << ','
      << result.tabi_mac_convergence_aspect_reject_panels << '\n';
}

void AnalyzerImpl::write_tabi_moment_audit_header_if_needed() {
  if (args_.timing != "detail")
    return;
  if (!tabi_moment_audit_out_.is_open()) {
    tabi_moment_audit_out_.open("pbe_tabi_moment_audit.csv",
                                std::ios::out | std::ios::trunc);
    if (!tabi_moment_audit_out_.is_open()) {
      throw std::runtime_error(
          "PBE: failed to open pbe_tabi_moment_audit.csv.");
    }
  }
  if (!tabi_moment_audit_header_written_) {
    tabi_moment_audit_out_
        << "trajectory_id,frame_index,global_index,fixture,tree_builder,"
           "node_bounds,moment_audit,moment_construction,moment_sum,"
           "nodes,translation_rel_error_max,translation_rel_error_mean,"
           "translation_offset_max,translation_offset_mean,"
           "cancellation_ratio_max,cancellation_ratio_mean,"
           "dynamic_range_max,status\n";
    tabi_moment_audit_header_written_ = true;
  }
}

void AnalyzerImpl::append_tabi_moment_audit_row(const FrameResult &result) {
  if (args_.timing != "detail" || result.pbe_backend != "tabi" ||
      result.tabi_moment_audit == "off")
    return;
  write_tabi_moment_audit_header_if_needed();
  tabi_moment_audit_out_
      << result.trajectory_id << ',' << result.frame_index << ','
      << result.global_index << ',' << csv_escape(current_source_label_) << ','
      << csv_escape(result.tabi_tree_builder) << ','
      << csv_escape(result.tabi_node_bounds_mode) << ','
      << csv_escape(result.tabi_moment_audit) << ','
      << csv_escape(result.tabi_moment_construction) << ','
      << csv_escape(result.tabi_moment_sum) << ','
      << result.tabi_moment_audit_nodes << ','
      << result.tabi_moment_translation_rel_error_max << ','
      << result.tabi_moment_translation_rel_error_mean << ','
      << result.tabi_moment_translation_offset_max << ','
      << result.tabi_moment_translation_offset_mean << ','
      << result.tabi_moment_cancellation_ratio_max << ','
      << result.tabi_moment_cancellation_ratio_mean << ','
      << result.tabi_moment_dynamic_range_max << ','
      << csv_escape(result.status) << '\n';
}

void AnalyzerImpl::write_tabi_surface_quality_header_if_needed() {
  if (args_.timing != "detail")
    return;
  if (!tabi_surface_quality_out_.is_open()) {
    tabi_surface_quality_out_.open("pbe_tabi_surface_quality.csv",
                                   std::ios::out | std::ios::trunc);
    if (!tabi_surface_quality_out_.is_open()) {
      throw std::runtime_error(
          "PBE: failed to open pbe_tabi_surface_quality.csv.");
    }
  }
  if (!tabi_surface_quality_header_written_) {
    tabi_surface_quality_out_
        << "trajectory_id,frame_index,global_index,fixture,pbe_backend,"
           "status,sphere_mesh_mode,panel_count,total_area,"
           "area_rel_error_vs_4piR2,triangle_area_min,triangle_area_mean,"
           "triangle_area_max,triangle_aspect_max,normal_radial_dot_min,"
           "normal_radial_dot_mean,centroid_radial_error_max,"
           "planar_area,planar_area_rel_error_vs_4piR2,enclosed_volume,"
           "enclosed_volume_rel_error_vs_4piR3_over3,triangle_area_cv\n";
    tabi_surface_quality_header_written_ = true;
  }
}

void AnalyzerImpl::append_tabi_surface_quality_row(const FrameResult &result) {
  if (args_.timing != "detail" || result.pbe_backend != "tabi")
    return;
  write_tabi_surface_quality_header_if_needed();
  tabi_surface_quality_out_
      << result.trajectory_id << ',' << result.frame_index << ','
      << result.global_index << ',' << csv_escape(current_source_label_) << ','
      << csv_escape(result.pbe_backend) << ',' << csv_escape(result.status)
      << ',' << csv_escape(result.tabi_sphere_mesh_mode) << ','
      << result.tabi_surface_panels << ','
      << result.tabi_surface_quality_total_area_a2 << ','
      << result.tabi_surface_quality_area_rel_error_vs_4piR2 << ','
      << result.tabi_surface_quality_triangle_area_min << ','
      << result.tabi_surface_quality_triangle_area_mean << ','
      << result.tabi_surface_quality_triangle_area_max << ','
      << result.tabi_surface_quality_triangle_aspect_max << ','
      << result.tabi_surface_quality_normal_radial_dot_min << ','
      << result.tabi_surface_quality_normal_radial_dot_mean << ','
      << result.tabi_surface_quality_centroid_radial_error_max << ','
      << result.tabi_surface_quality_planar_area_a2 << ','
      << result.tabi_surface_quality_planar_area_rel_error_vs_4piR2 << ','
      << result.tabi_surface_quality_enclosed_volume_a3 << ','
      << result.tabi_surface_quality_enclosed_volume_rel_error_vs_4piR3_over3
      << ',' << result.tabi_surface_quality_triangle_area_cv << '\n';
}

void AnalyzerImpl::write_memory_profile_header_if_needed() {
  if (!pbe_timing_enabled())
    return;
  if (!memory_profile_out_.is_open()) {
    memory_profile_out_.open("pbe_memory_profile.csv",
                             std::ios::out | std::ios::trunc);
    if (!memory_profile_out_.is_open())
      throw std::runtime_error("PBE: failed to open pbe_memory_profile.csv.");
  }
  if (!memory_profile_header_written_) {
    memory_profile_out_
        << "trajectory_id,frame_index,global_index,fixture,backend,"
           "atom_count,grid_points,panels,unknown_count,formulation,"
           "matvec_mode,preconditioner_effective,surface_generator,"
           "discretization,estimated_bytes_total,estimated_bytes_surface,"
           "estimated_bytes_solver,estimated_bytes_operator,"
           "estimated_bytes_tree,rss_before_bytes,rss_after_bytes,"
           "rss_delta_bytes,rss_delta_signed_bytes,peak_working_set_bytes,total_seconds,"
           "solver_seconds,status\n";
    memory_profile_header_written_ = true;
  }
}

void AnalyzerImpl::append_memory_profile_row(const FrameResult &result,
                                             const TimingBuckets &timing) {
  if (!pbe_timing_enabled())
    return;
  write_memory_profile_header_if_needed();
  auto seconds = [](double ms) { return ms / 1000.0; };
  const bool is_tabi = result.pbe_backend == "tabi";
  const std::size_t panels = is_tabi ? result.tabi_surface_panels : 0;
  const std::size_t unknowns = is_tabi ? result.tabi_unknown_count : 0;
  std::size_t surface_bytes = result.memory_estimated_surface_bytes;
  std::size_t solver_bytes = result.memory_estimated_solver_bytes;
  std::size_t operator_bytes = result.memory_estimated_operator_bytes;
  std::size_t tree_bytes = result.memory_estimated_tree_bytes;
  if (surface_bytes == 0 || solver_bytes == 0 || operator_bytes == 0) {
    if (is_tabi) {
      surface_bytes = panels * 96u + result.tabi_mesh_vertices * 56u +
                      result.tabi_mesh_triangles * 80u;
      const std::size_t basis_count = 40u;
      solver_bytes = std::max<std::size_t>(1, unknowns) * sizeof(double) *
                     basis_count;
      operator_bytes =
          std::max<std::size_t>(1, unknowns) * sizeof(double) * 8u;
      tree_bytes = result.tabi_tree_nodes * 96u +
                   result.tabi_tree_leaf_count * 64u;
      tree_bytes += panels * 14u * sizeof(double) +
                    panels * sizeof(std::size_t);
      operator_bytes += panels * 2u * sizeof(double);
      const std::size_t matvec_calls =
          std::max<std::size_t>(1, result.tabi_matvec_calls);
      const std::size_t interaction_entries =
          (result.tabi_singular_near_interactions +
           result.tabi_leaf_residual_panel_interactions +
           result.tabi_leaf_tile_interactions +
           result.tabi_far_node_interactions +
           result.tabi_smooth_leaf_node_interactions) /
          matvec_calls;
      tree_bytes += interaction_entries * sizeof(std::size_t) +
                    (panels + 1u) * 5u * sizeof(std::size_t);
    } else {
      const std::size_t n = std::max<std::size_t>(1, result.grid_points);
      surface_bytes =
          (result.solute_points + result.solvent_points +
           result.ion_accessible_points + result.stern_excluded_points) *
          sizeof(std::uint32_t);
      solver_bytes = n * sizeof(double) * 8u;
      operator_bytes = n * sizeof(double) * 6u;
      tree_bytes = 0;
    }
  }
  const std::size_t estimated_total =
      std::max(result.memory_bytes,
               surface_bytes + solver_bytes + operator_bytes + tree_bytes);
  memory_profile_out_
      << result.trajectory_id << ',' << result.frame_index << ','
      << result.global_index << ',' << csv_escape(current_source_label_) << ','
      << csv_escape(result.pbe_backend) << ',' << result.atom_count << ','
      << result.grid_points << ',' << panels << ',' << unknowns << ','
      << csv_escape(is_tabi ? result.tabi_formulation : std::string("grid"))
      << ','
      << csv_escape(is_tabi ? result.tabi_matvec_mode : result.solver_mode)
      << ','
      << csv_escape(is_tabi ? result.tabi_preconditioner_effective
                            : std::string("grid"))
      << ','
      << csv_escape(is_tabi ? result.tabi_surface_generator
                            : std::string("grid"))
      << ','
      << csv_escape(is_tabi ? result.tabi_discretization
                            : std::string("finite_difference"))
      << ',' << estimated_total << ',' << surface_bytes << ','
      << solver_bytes << ',' << operator_bytes << ',' << tree_bytes << ','
      << result.memory_rss_before_bytes << ','
      << result.memory_rss_after_bytes << ','
      << result.memory_rss_delta_bytes << ','
      << result.memory_rss_delta_signed_bytes << ','
      << result.memory_peak_working_set_bytes << ','
      << seconds(timing.total_ms) << ',' << seconds(timing.solver_ms) << ','
      << csv_escape(result.status) << '\n';
}

void AnalyzerImpl::write_tabi_formulation_audit_header_if_needed() const {
  if (!tabi_formulation_audit_out_.is_open()) {
    tabi_formulation_audit_out_.open("pbe_tabi_formulation_audit.csv",
                                     std::ios::out | std::ios::trunc);
    if (!tabi_formulation_audit_out_.is_open())
      throw std::runtime_error(
          "PBE: failed to open pbe_tabi_formulation_audit.csv.");
  }
  if (!tabi_formulation_audit_header_written_) {
    tabi_formulation_audit_out_
        << "trajectory_id,frame_index,global_index,fixture,pbe_backend,"
           "status,tabi_status,convergence_reason,reaction_energy_source,"
           "formulation,energy_projection,preconditioner_requested,"
           "preconditioner_effective,surface_generator,discretization,"
           "fallback_used,rhs_scale,"
           "dielectric_contrast,operator_contrast,unknown_count,bie_system,"
            "lhs_diag_phi,lhs_diag_phin,energy_equation,kernel_convention,"
            "derivative_mode,source_norm_s1,source_norm_s2,phi_norm,phin_norm,"
           "diagonal_mode,diag_min,"
           "diag_max,kappa_inv_A,coulomb_factor_kcal_mol,"
           "boundary_energy_kJmol,reference_energy_kJmol,"
           "selected_energy_scale,selected_energy_kJmol,residual,iterations,"
           "surface_geometry,panels,solid_angle_fallbacks,matvec_mode,"
           "near_integration,near_self_pairs,near_one_ring_pairs,"
           "near_duffy_pairs,near_centroid_pairs,near_fallbacks,"
           "near_integration_seconds,duffy_subdivisions_max,"
           "duffy_subdivisions_mean,near_kernel_policy,"
           "near_k1_accepted,near_k1_rejected,near_k2_accepted,"
           "near_k2_rejected,near_k3_accepted,near_k3_rejected,"
           "near_k4_accepted,near_k4_rejected,"
           "near_one_ring_k1_accepted,near_one_ring_k1_rejected,"
           "near_one_ring_k2_accepted,near_one_ring_k2_rejected,"
           "near_one_ring_k3_accepted,near_one_ring_k3_rejected,"
           "near_one_ring_k4_accepted,near_one_ring_k4_rejected,"
           "near_k1_policy,"
           "near_k2_policy,near_k3_policy,near_k4_policy\n";
    tabi_formulation_audit_header_written_ = true;
  }
}

void AnalyzerImpl::append_tabi_formulation_audit_row(
    const FrameResult &result) const {
  if (result.pbe_backend != "tabi")
    return;
  write_tabi_formulation_audit_header_if_needed();
  auto seconds = [](double ms) { return ms / 1000.0; };
  tabi_formulation_audit_out_
      << result.trajectory_id << ',' << result.frame_index << ','
      << result.global_index << ',' << csv_escape(current_source_label_)
      << ',' << csv_escape(result.pbe_backend) << ','
      << csv_escape(result.status) << ',' << csv_escape(result.tabi_status)
      << ',' << csv_escape(result.tabi_convergence_reason) << ','
      << csv_escape(result.tabi_reaction_energy_source) << ','
      << csv_escape(result.tabi_formulation) << ','
      << csv_escape(result.tabi_energy_projection) << ','
      << csv_escape(result.tabi_preconditioner_requested) << ','
      << csv_escape(result.tabi_preconditioner_effective) << ','
      << csv_escape(result.tabi_surface_generator) << ','
      << csv_escape(result.tabi_discretization) << ','
      << (result.tabi_fallback_used ? "true" : "false") << ','
      << result.tabi_rhs_scale << ',' << result.tabi_dielectric_contrast
      << ',' << result.tabi_operator_contrast << ','
      << result.tabi_unknown_count << ',' << csv_escape(result.tabi_bie_system)
      << ',' << result.tabi_lhs_diag_phi << ','
      << result.tabi_lhs_diag_phin << ','
      << csv_escape(result.tabi_energy_equation) << ','
      << csv_escape(result.tabi_kernel_convention) << ','
      << csv_escape(result.tabi_derivative_mode) << ','
      << result.tabi_source_norm_s1 << ',' << result.tabi_source_norm_s2
      << ',' << result.tabi_phi_norm << ',' << result.tabi_phin_norm << ','
      << csv_escape(args_.tabi_diagonal) << ',' << result.tabi_diag_min << ','
      << result.tabi_diag_max << ',' << result.tabi_kappa_inv_a << ','
      << result.tabi_coulomb_factor_kcal_mol << ','
      << result.tabi_boundary_energy_kj_mol << ','
      << result.tabi_reference_energy_kj_mol << ','
      << result.tabi_selected_energy_scale << ',' << result.energy_kj_mol
      << ',' << result.linear_residual << ',' << result.solvent_iterations
      << ',' << csv_escape(result.tabi_surface_geometry) << ','
      << result.tabi_surface_panels << ','
      << result.tabi_solid_angle_fallbacks << ','
      << csv_escape(result.tabi_matvec_mode) << ','
      << csv_escape(result.tabi_near_integration) << ','
      << result.tabi_near_self_pairs << ','
      << result.tabi_near_one_ring_pairs << ','
      << result.tabi_near_duffy_pairs << ','
      << result.tabi_near_centroid_pairs << ','
      << result.tabi_near_fallbacks << ','
      << seconds(result.tabi_near_integration_ms) << ','
      << result.tabi_duffy_subdivisions_max << ','
      << result.tabi_duffy_subdivisions_mean << ','
      << csv_escape(result.tabi_near_kernel_policy) << ','
      << result.tabi_near_k1_accepted << ','
      << result.tabi_near_k1_rejected << ','
      << result.tabi_near_k2_accepted << ','
      << result.tabi_near_k2_rejected << ','
      << result.tabi_near_k3_accepted << ','
      << result.tabi_near_k3_rejected << ','
      << result.tabi_near_k4_accepted << ','
      << result.tabi_near_k4_rejected << ','
      << result.tabi_near_one_ring_k1_accepted << ','
      << result.tabi_near_one_ring_k1_rejected << ','
      << result.tabi_near_one_ring_k2_accepted << ','
      << result.tabi_near_one_ring_k2_rejected << ','
      << result.tabi_near_one_ring_k3_accepted << ','
      << result.tabi_near_one_ring_k3_rejected << ','
      << result.tabi_near_one_ring_k4_accepted << ','
      << result.tabi_near_one_ring_k4_rejected << ','
      << csv_escape(result.tabi_near_k1_policy) << ','
      << csv_escape(result.tabi_near_k2_policy) << ','
      << csv_escape(result.tabi_near_k3_policy) << ','
      << csv_escape(result.tabi_near_k4_policy) << '\n';
}

void AnalyzerImpl::maybe_write_parity_stage(const FrameResult &result,
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
           "nl_inner_effective,nl_status,boundary_mode,"
           "test_boundary_convention,test_field_units,"
           "test_operator_convention,charge_assignment,"
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
           "iscm_vacuum_rf_kJmol,iscm_water_minus_vacuum_kT,"
           "iscm_water_minus_vacuum_kJmol,iscm_corrected_rxn_kT,"
           "iscm_corrected_rxn_kJmol,iscm_finalization_mode,"
           "iscm_position_source,iscm_position_rms_vs_oracle,"
           "iscm_position_oracle_matches,iscm_native_independent,iscm_status,"
           "parity_map_dump,"
           "test_charge_grid,test_surface_map,test_boundary_map,"
           "status,test_midpoint_epsilon,delphi_face_crossings,"
           "delphi_eps_face_matches,delphi_boundary_grid_points,"
            "delphi_charged_boundary_points,delphi_iscm_status,"
            "test_vdwms_remap,vdwms_enabled,vdwms_boundary_initial,"
            "vdwms_boundary_final,vdwms_added,vdwms_removed,"
            "vdwms_exposed_points,vdwms_midpoint_rms,"
            "midpoint_eps_sum_min,midpoint_eps_sum_max,"
            "midpoint_eps_sum_mean,midpoint_eps_sum_interior_mean,"
            "nonlinear_stress_kT,nonlinear_stress_kJmol,"
            "nonlinear_osmotic_kT,nonlinear_osmotic_kJmol,"
            "nonlinear_energy_correction_kT,"
            "nonlinear_energy_correction_kJmol,"
            "nonlinear_energy_status\n";
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
      << result.test_boundary_convention << ',' << result.test_field_units
      << ',' << result.test_operator_convention << ','
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
      << result.iscm_water_minus_vacuum_kt << ','
      << result.iscm_water_minus_vacuum_kj_mol << ','
      << result.iscm_corrected_kt << ',' << result.iscm_corrected_kj_mol
      << ',' << csv_escape(result.iscm_finalization_mode) << ','
      << csv_escape(result.iscm_position_source) << ','
      << result.iscm_position_rms_vs_oracle << ','
      << result.iscm_position_oracle_matches << ','
      << (result.iscm_native_independent ? "true" : "false") << ','
      << csv_escape(result.iscm_status) << ','
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
      << audit.midpoint_eps_sum_interior_mean << ','
      << result.nonlinear_stress_kt << ','
      << result.nonlinear_stress_kj_mol << ','
      << result.nonlinear_osmotic_kt << ','
      << result.nonlinear_osmotic_kj_mol << ','
      << result.nonlinear_energy_correction_kt << ','
      << result.nonlinear_energy_correction_kj_mol << ','
      << csv_escape(result.nonlinear_energy_status)
      << '\n';
}

bool AnalyzerImpl::barnase_forensic_enabled() const {
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

void AnalyzerImpl::maybe_write_iscm_surface_charge_samples(
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

  std::ofstream terms;
  if (!iscm_energy_terms_header_written_) {
    terms.open("pbe_iscm_energy_terms.csv", std::ios::out | std::ios::trunc);
    if (!terms.is_open())
      throw std::runtime_error(
          "PBE: failed to open pbe_iscm_energy_terms.csv.");
    terms << "trajectory_id,frame_index,global_index,rank,i,j,k,"
             "surface_x_A,surface_y_A,surface_z_A,induced_charge,"
             "atom_potential_sum,energy_contribution_kT,"
             "energy_contribution_kJmol,iscm_status\n";
    iscm_energy_terms_header_written_ = true;
  } else {
    terms.open("pbe_iscm_energy_terms.csv", std::ios::out | std::ios::app);
    if (!terms.is_open())
      throw std::runtime_error(
          "PBE: failed to append pbe_iscm_energy_terms.csv.");
  }
  terms << std::setprecision(12);
  std::size_t rank = 0;
  for (const IscmSample &sample : audit.samples) {
    if (sample.kind != "top_energy")
      continue;
    ++rank;
    terms << result.trajectory_id << ',' << result.frame_index << ','
          << result.global_index << ',' << rank << ',' << sample.i << ','
          << sample.j << ',' << sample.k << ',' << sample.x << ','
          << sample.y << ',' << sample.z << ',' << sample.induced_charge << ','
          << sample.atom_potential_sum << ',' << sample.energy_contribution_kt
          << ',' << sample.energy_contribution_kt * kt_kj << ','
          << csv_escape(audit.status) << '\n';
  }
}

void AnalyzerImpl::maybe_write_barnase_forensic_audits(
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
          "iscm_water_minus_vacuum_kJmol,iscm_corrected_rxn_kJmol,"
          "iscm_finalization_mode,iscm_position_source,"
          "iscm_position_rms_vs_oracle,iscm_position_oracle_matches,"
          "iscm_native_independent,iscm_status,"
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
     << solve.iscm_water_minus_vacuum_kj_mol << ','
     << solve.iscm_corrected_kj_mol << ','
     << csv_escape(solve.iscm_finalization_mode) << ','
     << csv_escape(solve.iscm_position_source) << ','
     << solve.iscm_position_rms_vs_oracle << ','
     << solve.iscm_position_oracle_matches << ','
     << (solve.iscm_native_independent ? "true" : "false") << ','
     << csv_escape(solve.iscm_status)
     << ',' << csv_escape(solve.vdwms_enabled) << ','
     << solve.vdwms_boundary_initial << ',' << solve.vdwms_boundary_final
     << ',' << solve.vdwms_added << ',' << solve.vdwms_removed << ','
     << solve.vdwms_exposed_points << ',' << solve.vdwms_midpoint_rms << ','
     << (solve.converged ? "ok" : "not_converged") << '\n';

  maybe_write_atom_phi_audit(result, atoms, solvent, vacuum, solve);
  maybe_write_linear_forensic_points(result, atoms, solvent, vacuum, solve);
}

void AnalyzerImpl::maybe_write_atom_phi_audit(
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

void AnalyzerImpl::maybe_write_linear_forensic_points(
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

  const double rhs = effective_rhs_scale(shape);
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

bool AnalyzerImpl::parity_map_dump_enabled() const noexcept {
  return !args_.test_pqr.empty() && args_.parity_map_dump != "off";
}

std::vector<AnalyzerImpl::MapSamplePoint>
AnalyzerImpl::collect_map_sample_points(const GridShape &shape,
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

void AnalyzerImpl::maybe_write_parity_map_dumps(
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
             "vacuum_phi,reaction_phi,atom_energy_hint_kJmol,is_boundary,"
             "field_units,solvent_phi_raw,solvent_phi_native,"
             "vacuum_phi_raw,vacuum_phi_native,reaction_phi_native\n";
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
      const double to_native = field_to_native_factor();
      const double to_raw =
          delphi_field_units_active() ? 1.0 : delphi_epkt();
      const double solvent_raw = solvent.phi_grid[id] * to_raw;
      const double vacuum_raw = vacuum.phi_grid[id] * to_raw;
      const double solvent_native = solvent.phi_grid[id] * to_native;
      const double vacuum_native = vacuum.phi_grid[id] * to_native;
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
          << (is_boundary(i, j, k, shape) ? 1 : 0) << ','
          << csv_escape(args_.test_field_units) << ',' << solvent_raw << ','
          << solvent_native << ',' << vacuum_raw << ',' << vacuum_native
          << ',' << (solvent_native - vacuum_native) << '\n';
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
             "eps_x,eps_y,eps_z,screen,field_units,solvent_phi_raw,"
             "solvent_phi_native,vacuum_phi_raw,vacuum_phi_native\n";
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
      const double to_native = field_to_native_factor();
      const double to_raw =
          delphi_field_units_active() ? 1.0 : delphi_epkt();
      out << result.trajectory_id << ',' << result.frame_index << ','
          << result.global_index << ',' << csv_escape(kind) << ',' << rank
          << ',' << i << ',' << j << ',' << k << ','
          << shape.ox + static_cast<double>(i) * shape.h << ','
          << shape.oy + static_cast<double>(j) * shape.h << ','
          << shape.oz + static_cast<double>(k) * shape.h << ','
          << solvent.phi_grid[id] << ',' << vacuum.phi_grid[id] << ','
          << rho_at(id) << ',' << inside_at(id) << ',' << eps_sum_at(id)
          << ',' << eps_x_at(id) << ',' << eps_y_at(id) << ','
          << eps_z_at(id) << ',' << screen_at(id) << ','
          << csv_escape(args_.test_field_units) << ','
          << solvent.phi_grid[id] * to_raw << ','
          << solvent.phi_grid[id] * to_native << ','
          << vacuum.phi_grid[id] * to_raw << ','
          << vacuum.phi_grid[id] * to_native << '\n';
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
    if (!operator_stage_header_written_) {
      out.open("pbe_operator_stage.csv", std::ios::out | std::ios::trunc);
      if (!out.is_open())
        throw std::runtime_error("PBE: failed to open pbe_operator_stage.csv.");
      out << "trajectory_id,frame_index,global_index,point_kind,rank,i,j,k,"
             "rhs_scale_effective,rho,denom,screen,eps_x_plus,eps_x_minus,"
             "eps_y_plus,eps_y_minus,eps_z_plus,eps_z_minus,neighbor_sum,"
             "rhs,F,residual_abs,phi,phi_native,field_units,"
             "rhs_delta_vs_oracle,screen_delta_vs_oracle,"
             "denom_delta_vs_oracle,fixed_boundary_used,"
             "operator_convention\n";
      operator_stage_header_written_ = true;
    } else {
      out.open("pbe_operator_stage.csv", std::ios::out | std::ios::app);
      if (!out.is_open())
        throw std::runtime_error("PBE: failed to append pbe_operator_stage.csv.");
    }
    out << std::setprecision(12);
    const double rhs = effective_rhs_scale(shape);
    const double to_native = field_to_native_factor();
    struct OperatorOracleRow {
      bool rhs_present{false};
      bool screen_present{false};
      bool denom_present{false};
      double rhs{0.0};
      double screen{0.0};
      double denom{0.0};
    };
    std::unordered_map<std::size_t, OperatorOracleRow> operator_oracle;
    {
      std::ifstream oracle("pbe_oracle_operator_stage.csv");
      if (oracle.is_open()) {
        std::string line;
        if (std::getline(oracle, line)) {
          const auto header = csv_header_map(split_csv_simple(line));
          const std::string fixture_filter = read_oracle_fixture_filter();
          while (std::getline(oracle, line)) {
            if (trim_copy(line).empty())
              continue;
            const auto fields = split_csv_simple(line);
            if (!oracle_fixture_matches(fields, header, fixture_filter))
              continue;
            double vacuum_flag = 0.0;
            if (csv_double(fields, header, {"vacuum"}, vacuum_flag) &&
                std::llround(vacuum_flag) != 0)
              continue;
            std::size_t i = 0, j = 0, k = 0;
            if (!csv_index(fields, header, "i", i) ||
                !csv_index(fields, header, "j", j) ||
                !csv_index(fields, header, "k", k) || i >= shape.nx ||
                j >= shape.ny || k >= shape.nz)
              continue;
            OperatorOracleRow row;
            row.rhs_present =
                csv_double(fields, header,
                           {"charge_map_scaled", "rhs", "rhs_value",
                            "charge_density"},
                           row.rhs);
            row.screen_present =
                csv_double(fields, header, {"screen", "screen_term"},
                           row.screen);
            row.denom_present =
                csv_double(fields, header,
                           {"denom", "eps_sum_plus_screen",
                            "epsilon_sum_plus_screen"},
                           row.denom);
            operator_oracle[idx3(i, j, k, shape.nx, shape.ny)] = row;
          }
        }
      }
    }
    auto write_operator = [&](const std::string &kind, std::size_t rank,
                              std::size_t id) {
      if (id >= shape.n)
        return;
      std::size_t i = 0, j = 0, k = 0;
      coords(id, i, j, k);
      if (is_boundary(i, j, k, shape) || i == 0 || j == 0 || k == 0 ||
          i + 1 >= shape.nx || j + 1 >= shape.ny || k + 1 >= shape.nz)
        return;
      const std::size_t xy = shape.nx * shape.ny;
      const double eps_x_plus = eps_x_at(id);
      const double eps_x_minus = eps_x_at(id - 1);
      const double eps_y_plus = eps_y_at(id);
      const double eps_y_minus = eps_y_at(id - shape.nx);
      const double eps_z_plus = eps_z_at(id);
      const double eps_z_minus = eps_z_at(id - xy);
      const double denom = eps_x_plus + eps_x_minus + eps_y_plus +
                           eps_y_minus + eps_z_plus + eps_z_minus +
                           screen_at(id);
      const double neighbor_sum =
          eps_x_plus * solvent.phi_grid[id + 1] +
          eps_x_minus * solvent.phi_grid[id - 1] +
          eps_y_plus * solvent.phi_grid[id + shape.nx] +
          eps_y_minus * solvent.phi_grid[id - shape.nx] +
          eps_z_plus * solvent.phi_grid[id + xy] +
          eps_z_minus * solvent.phi_grid[id - xy];
      const double rhs_value = rhs * rho_at(id);
      const double residual =
          denom * solvent.phi_grid[id] - neighbor_sum - rhs_value;
      const auto oracle_it = operator_oracle.find(id);
      const bool has_oracle = oracle_it != operator_oracle.end();
      const bool fixed_boundary_used =
          is_boundary(i + 1, j, k, shape) || is_boundary(i - 1, j, k, shape) ||
          is_boundary(i, j + 1, k, shape) || is_boundary(i, j - 1, k, shape) ||
          is_boundary(i, j, k + 1, shape) || is_boundary(i, j, k - 1, shape);
      out << result.trajectory_id << ',' << result.frame_index << ','
          << result.global_index << ',' << csv_escape(kind) << ',' << rank
          << ',' << i << ',' << j << ',' << k << ',' << rhs << ','
          << rho_at(id) << ',' << denom << ',' << screen_at(id) << ','
          << eps_x_plus << ',' << eps_x_minus << ',' << eps_y_plus << ','
          << eps_y_minus << ',' << eps_z_plus << ',' << eps_z_minus << ','
          << neighbor_sum << ',' << rhs_value << ',' << residual << ','
          << std::abs(residual) << ',' << solvent.phi_grid[id] << ','
          << solvent.phi_grid[id] * to_native << ','
          << csv_escape(args_.test_field_units);
      if (has_oracle && oracle_it->second.rhs_present)
        out << ',' << (rhs_value - oracle_it->second.rhs);
      else
        out << ",NA";
      if (has_oracle && oracle_it->second.screen_present)
        out << ',' << (screen_at(id) - oracle_it->second.screen);
      else
        out << ",NA";
      if (has_oracle && oracle_it->second.denom_present)
        out << ',' << (denom - oracle_it->second.denom);
      else
        out << ",NA";
      out << ',' << (fixed_boundary_used ? 1 : 0) << ','
          << csv_escape(args_.test_operator_convention) << '\n';
    };
    if (full) {
      std::size_t rank = 0;
      for (std::size_t id = 0; id < shape.n; ++id)
        write_operator("full_operator", ++rank, id);
    } else {
      for (const MapSamplePoint &sample : samples)
        write_operator(sample.kind, sample.rank, sample.id);
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

void AnalyzerImpl::maybe_write_parity_points(
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

AnalyzerImpl::ReactionSolve AnalyzerImpl::solve_reaction_field(
    const std::vector<Atom> &atoms, TitanExecutor &executor,
    const FrameResult *map_result) {
  ReactionSolve solve;
  TimingBuckets solvent_timing;
  TimingBuckets vacuum_timing;
  TimingBuckets reuse_timing;
  const bool capture_pbsa_timing =
      args_.pbsa && args_.pbsa_perf == "detail" &&
      !stability_probe_solve_active_;
  struct PbsaTimingScope {
    AnalyzerImpl &self;
    TimingBuckets *previous;
    bool capture;
    PbsaTimingScope(AnalyzerImpl &impl, TimingBuckets &timing, bool enabled)
        : self(impl), previous(impl.active_timing_),
          capture(enabled) {
      if (capture)
        self.active_timing_ = &timing;
    }
    ~PbsaTimingScope() {
      if (capture)
        self.active_timing_ = previous;
    }
  };
  const bool validate_reuse =
      args_.pbsa_grid_reuse == "validate" && args_.pbsa &&
      !stability_probe_solve_active_ && args_.backend == "grid" &&
      args_.model == "linear" && args_.biomodel == "pbe" &&
      args_.dielectric_model == "two-dielectric" &&
      args_.charge_assignment == "trilinear" &&
      args_.test_energy_convention == "auto" && args_.focus == "off" &&
      effective_solver_mode() == "sor";
  const double solvent_kappa = debye_kappa_inv_a(args_.salt, args_.temp);
  FieldSolve solvent;
  {
    PbsaTimingScope timing_scope(*this, solvent_timing, capture_pbsa_timing);
    solvent = solve_physical_field(atoms, args_.int_diel, args_.ext_diel,
                                   solvent_kappa, executor, map_result);
  }
  const ForensicAudit solvent_forensic =
      barnase_forensic_enabled() ? compute_forensic_audit(solvent.shape)
                                 : ForensicAudit{};
  const bool wants_iscm =
      args_.test_energy_convention == "delphi-iscm" ||
      args_.test_energy_convention == "delphi-corrected" ||
      delphi_nwt_parity_active() ||
      barnase_forensic_enabled();
  const AlignedDoubleVector solvent_screen_for_iscm =
      wants_iscm ? screen_ : AlignedDoubleVector{};
  FieldSolve vacuum;
  {
    PbsaTimingScope timing_scope(*this, vacuum_timing, capture_pbsa_timing);
    vacuum = solve_linear_field(atoms, args_.int_diel, args_.int_diel, 0.0,
                                executor, nullptr);
  }
  FieldSolve reused_vacuum;
  if (validate_reuse) {
    struct ReuseScope {
      AnalyzerImpl &self;
      ReuseScope(AnalyzerImpl &impl, const GridShape &shape) : self(impl) {
        self.pbsa_reuse_shape_ = shape;
        self.pbsa_reuse_vacuum_active_ = true;
      }
      ~ReuseScope() { self.pbsa_reuse_vacuum_active_ = false; }
    } reuse_scope(*this, solvent.shape);
    PbsaTimingScope timing_scope(*this, reuse_timing, capture_pbsa_timing);
    reused_vacuum = solve_linear_field(atoms, args_.int_diel, args_.int_diel,
                                       0.0, executor, nullptr);
  }

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
  solve.pbsa_solvent_timing = solvent_timing;
  solve.pbsa_vacuum_timing = vacuum_timing;
  solve.pbsa_reuse_timing = reuse_timing;
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
  if (validate_reuse) {
    solve.pbsa_reuse_surface_hits =
        reused_vacuum.pbsa_reuse_surface_hit ? 1U : 0U;
    solve.pbsa_reuse_operator_hits =
        reused_vacuum.pbsa_reuse_operator_hit ? 1U : 0U;
  }

  const auto energy_time_begin = std::chrono::steady_clock::now();
  double rpb_source_time_ms = 0.0;
  double probe_time_ms = 0.0;
  solve.atom_reaction_kj_mol.assign(atoms.size(), 0.0);
  const double kt_kcal =
      kKBKcalPerMolK * std::max(1.0, args_.temp);
  solve.solvent_phase_kcal_mol =
      atom_phase_energy_kcal(atoms, solvent.atom_phi);
  solve.vacuum_phase_kcal_mol =
      atom_phase_energy_kcal(atoms, vacuum.atom_phi);
  if (validate_reuse) {
    const double reused_vacuum_phase_kcal =
        atom_phase_energy_kcal(atoms, reused_vacuum.atom_phi);
    solve.pbsa_reuse_delta_kj_mol =
        (reused_vacuum_phase_kcal - solve.vacuum_phase_kcal_mol) * kKcalToKJ;
  }
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
  NonlinearEnergyAudit nonlinear_energy_audit;
  RpbAudit rpb_audit;
  AlignedDoubleVector rpb_source_map;
  AlignedDoubleVector rpb_coulomb_map;
  const AlignedDoubleVector *rpb_source_top_map = &rpb_source_map;
  const AlignedDoubleVector *rpb_coulomb_top_map = &rpb_coulomb_map;
  if (has_grid_energy) {
    const AlignedDoubleVector *solvent_energy_rho = &solvent.rho_grid;
    if (solvent.rho_energy_grid.size() == solvent.shape.n)
      solvent_energy_rho = &solvent.rho_energy_grid;
    const AlignedDoubleVector *vacuum_energy_rho = &vacuum.rho_grid;
    if (vacuum.rho_energy_grid.size() == solvent.shape.n)
      vacuum_energy_rho = &vacuum.rho_energy_grid;
    if (solvent_energy_rho->size() != solvent.shape.n)
      solvent_energy_rho = &rho_;
    if (vacuum_energy_rho->size() != solvent.shape.n)
      vacuum_energy_rho = solvent_energy_rho;
    const double grid_solvent_phase_kcal =
        grid_phase_energy_kcal(solvent.shape, solvent.phi_grid,
                               *solvent_energy_rho);
    const double grid_vacuum_phase_kcal =
        grid_phase_energy_kcal(solvent.shape, vacuum.phi_grid,
                               *vacuum_energy_rho);
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
    if (args_.model == "nonlinear") {
      nonlinear_energy_audit = compute_nonlinear_energy_audit(
          solvent.shape, solvent.phi_grid, solvent_screen_for_iscm);
      copy_nonlinear_energy_audit_to_solve(nonlinear_energy_audit, solve);
    }
    if (args_.biomodel == "rpbe") {
      const bool wants_rpb_source_top =
          map_result != nullptr && args_.parity_dump != "off";
      const bool has_solvent_rpb_audit =
          solvent.rpbe_source_status != "off" &&
          solvent.rpbe_source_status != "missing_maps" &&
          solvent.rpbe_source_status != "empty_interior";
      const bool has_solvent_rpb_top =
          solvent.rpbe_source_grid.size() == solvent.shape.n &&
          solvent.rpbe_coulomb_grid.size() == solvent.shape.n;
      if (has_solvent_rpb_audit &&
          (!wants_rpb_source_top || has_solvent_rpb_top)) {
        rpb_audit.coulomb_norm = solvent.rpbe_coulomb_norm;
        rpb_audit.grad_coulomb_norm = solvent.rpbe_grad_coulomb_norm;
        rpb_audit.grad_epsilon_norm = solvent.rpbe_grad_epsilon_norm;
        rpb_audit.source_rms = solvent.rpbe_source_rms;
        rpb_audit.source_max = solvent.rpbe_source_max;
        rpb_audit.sampled_points =
            solvent.shape.n > 0
                ? (solvent.shape.nx > 2 && solvent.shape.ny > 2 &&
                           solvent.shape.nz > 2
                       ? (solvent.shape.nx - 2) * (solvent.shape.ny - 2) *
                             (solvent.shape.nz - 2)
                       : 0)
                : 0;
        rpb_audit.status = solvent.rpbe_source_status;
        if (has_solvent_rpb_top) {
          rpb_source_top_map = &solvent.rpbe_source_grid;
          rpb_coulomb_top_map = &solvent.rpbe_coulomb_grid;
        }
      } else {
        const auto rpb_time_begin = std::chrono::steady_clock::now();
        rpb_audit = compute_rpb_source_audit(solvent.shape, atoms, executor,
                                             solvent.eps_x_grid,
                                             solvent.eps_y_grid,
                                             solvent.eps_z_grid,
                                             wants_rpb_source_top
                                                 ? &rpb_source_map
                                                 : nullptr,
                                             wants_rpb_source_top
                                                 ? &rpb_coulomb_map
                                                 : nullptr);
        const double rpb_ms =
            elapsed_ms(rpb_time_begin, std::chrono::steady_clock::now());
        rpb_source_time_ms += rpb_ms;
        add_timing_ms(TimingBucket::RpbSource, rpb_ms);
      }
      solve.rpbe_coulomb_norm = rpb_audit.coulomb_norm;
      solve.rpbe_grad_coulomb_norm = rpb_audit.grad_coulomb_norm;
      solve.rpbe_grad_epsilon_norm = rpb_audit.grad_epsilon_norm;
      solve.rpbe_source_rms = rpb_audit.source_rms;
      solve.rpbe_source_max = rpb_audit.source_max;
      solve.rpbe_source_status = rpb_audit.status;
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
                                  "delphi-corrected" ||
                              delphi_nwt_parity_active());
    if (forced_iscm) {
      solve.solvent_phase_kcal_mol = solve.iscm_water_phase_kcal_mol;
      solve.vacuum_phase_kcal_mol = solve.iscm_vacuum_phase_kcal_mol;
      solve.energy_kcal_mol = solve.iscm_corrected_kcal_mol;
      solve.reaction_energy_source =
          delphi_nwt_parity_active()
              ? "delphi-nonlinear-iscm"
              : (args_.test_energy_convention == "delphi-corrected"
                     ? "delphi-corrected-iscm"
                     : "delphi-iscm");
      if (args_.model == "nonlinear" &&
          nonlinear_energy_audit.status == "ok") {
        solve.energy_kcal_mol +=
            nonlinear_energy_audit.total_correction_kcal_mol;
        solve.reaction_energy_source += "+nonlinear";
      }
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
        0.5 * kCoulombKcalPerMol * field_to_native_factor() * atoms[i].q *
        (solvent.atom_phi[i] - vacuum.atom_phi[i]);
    solve.atom_reaction_kj_mol[i] = atom_kcal_mol * kKcalToKJ;
  }
  if (!args_.probe_file.empty()) {
    const auto probe_time_begin = std::chrono::steady_clock::now();
    solve.probe_energy_status = "missing_field";
    if (solvent.phi_grid.size() == solvent.shape.n) {
      const std::vector<Atom> probe_atoms = load_probe_atoms(args_.probe_file);
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
      double q_phi = 0.0;
      std::size_t matched = 0;
      std::string matched_ids;
      std::string matched_groups;
      auto cubic_interp = [](const double p[4], double x) {
        return x * x * (p[3] - 2.0 * p[2] + p[1]) / 2.0 +
               x * (4.0 * p[2] - 3.0 * p[1] - p[3]) / 2.0 + p[1];
      };
      auto interpolate_probe_phi = [&](double x, double y, double z,
                                       double &phi) {
        const double gx = (x - solvent.shape.ox) / solvent.shape.h;
        const double gy = (y - solvent.shape.oy) / solvent.shape.h;
        const double gz = (z - solvent.shape.oz) / solvent.shape.h;
        if (gx >= 1.0 && gy >= 1.0 && gz >= 1.0 &&
            gx < static_cast<double>(solvent.shape.nx - 2) &&
            gy < static_cast<double>(solvent.shape.ny - 2) &&
            gz < static_cast<double>(solvent.shape.nz - 2)) {
          const auto ix = static_cast<std::size_t>(std::floor(gx));
          const auto iy = static_cast<std::size_t>(std::floor(gy));
          const auto iz = static_cast<std::size_t>(std::floor(gz));
          if (std::abs(gx - static_cast<double>(ix)) <= 1.0e-6 &&
              std::abs(gy - static_cast<double>(iy)) <= 1.0e-6 &&
              std::abs(gz - static_cast<double>(iz)) <= 1.0e-6) {
            phi = solvent.phi_grid[idx3(ix, iy, iz, solvent.shape.nx,
                                        solvent.shape.ny)];
            return true;
          }
          const double fx = gx - static_cast<double>(ix);
          const double fy = gy - static_cast<double>(iy);
          const double fz = gz - static_cast<double>(iz);
          double xvals[4]{};
          for (int lx = 0; lx < 4; ++lx) {
            double yvals[4]{};
            for (int ly = 0; ly < 4; ++ly) {
              double zvals[4]{};
              for (int lz = 0; lz < 4; ++lz) {
                const std::size_t px =
                    static_cast<std::size_t>(static_cast<int>(ix) + lx - 1);
                const std::size_t py =
                    static_cast<std::size_t>(static_cast<int>(iy) + ly - 1);
                const std::size_t pz =
                    static_cast<std::size_t>(static_cast<int>(iz) + lz - 1);
                zvals[lz] =
                    solvent.phi_grid[idx3(px, py, pz, solvent.shape.nx,
                                          solvent.shape.ny)];
              }
              yvals[ly] = cubic_interp(zvals, fz);
            }
            xvals[lx] = cubic_interp(yvals, fy);
          }
          phi = cubic_interp(xvals, fx);
          return true;
        }
        return interpolate_field_value(solvent.shape, solvent.phi_grid, x, y,
                                       z, phi);
      };
      for (std::size_t i = 0; i < probe_atoms.size(); ++i) {
        const Atom &atom = probe_atoms[i];
        double phi = 0.0;
        if (!interpolate_probe_phi(atom.x, atom.y, atom.z, phi))
          continue;
        q_phi += atom.q * phi;
        if (!matched_ids.empty())
          matched_ids += ';';
        matched_ids += std::to_string(i + 1);
        if (!matched_groups.empty())
          matched_groups += ';';
        matched_groups += uppercase_letters(atom.resname) +
                          std::to_string(atom.resnr) + ':' +
                          uppercase_alnum(atom.atom_name);
        ++matched;
      }
      if (map_result != nullptr && args_.parity_dump != "off") {
        const std::string sample_filename = "pbe_probe_atom_samples.csv";
        const bool need_header =
            !std::filesystem::exists(sample_filename) ||
            std::filesystem::file_size(sample_filename) == 0;
        std::ofstream probe_out(sample_filename, std::ios::out | std::ios::app);
        if (probe_out.is_open()) {
          if (need_header) {
            probe_out
                << "trajectory_id,frame_index,global_index,time_ps,atom_index,"
                   "atom_name,resname,resnr,x_A,y_A,z_A,q,phi_native,"
                   "phi_kT_per_e,contribution_kT,status\n";
          }
          const double kt_kj = kKBKJPerMolK * std::max(1.0, args_.temp);
          const double native_to_kt_per_e =
              kt_kj > 0.0 ? (kCoulombKcalPerMol * kKcalToKJ / kt_kj) : 0.0;
          for (std::size_t i = 0; i < probe_atoms.size(); ++i) {
            const Atom &atom = probe_atoms[i];
            double phi = 0.0;
            const bool ok = interpolate_probe_phi(atom.x, atom.y, atom.z, phi);
            const double phi_native = ok ? field_to_native_factor() * phi : 0.0;
            const double phi_kt = phi_native * native_to_kt_per_e;
            probe_out << map_result->trajectory_id << ','
                      << map_result->frame_index << ','
                      << map_result->global_index << ',' << std::setprecision(10)
                      << map_result->time_ps << ',' << (i + 1) << ','
                      << csv_escape(atom.atom_name) << ','
                      << csv_escape(atom.resname) << ',' << atom.resnr << ','
                      << atom.x << ',' << atom.y << ',' << atom.z << ','
                      << atom.q << ',' << phi_native << ',' << phi_kt << ','
                      << (ok ? 0.5 * atom.q * phi_kt : 0.0) << ','
                      << (ok ? "ok" : "outside_grid") << '\n';
          }
        }
      }
      solve.probe_matched_atoms = matched;
      solve.probe_matched_ids = matched_ids;
      solve.probe_matched_groups = matched_groups;
      solve.probe_matching_semantics = "delphi_frc_probe_file";
      if (matched > 0) {
        solve.probe_energy_kcal_mol =
            0.5 * kCoulombKcalPerMol * field_to_native_factor() * q_phi;
        solve.probe_energy_kj_mol =
            solve.probe_energy_kcal_mol * kKcalToKJ;
        const double kt_kj = kKBKJPerMolK * std::max(1.0, args_.temp);
        solve.probe_energy_kt =
            kt_kj > 0.0 ? solve.probe_energy_kj_mol / kt_kj : 0.0;
        solve.probe_energy_status = "ok";
      } else {
        solve.probe_energy_status = "no_matching_atoms";
      }
    }
    const double probe_ms =
        elapsed_ms(probe_time_begin, std::chrono::steady_clock::now());
    probe_time_ms += probe_ms;
    add_timing_ms(TimingBucket::Probe, probe_ms);
    add_timing_ms(TimingBucket::ProbeEnergy, probe_ms);
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
  const double energy_ms =
      elapsed_ms(energy_time_begin, std::chrono::steady_clock::now()) -
      rpb_source_time_ms - probe_time_ms;
  add_timing_ms(TimingBucket::Energy, std::max(0.0, energy_ms));
  solve.pbsa_energy_ms = std::max(0.0, energy_ms);
  if (map_result != nullptr) {
    const auto audit_output_begin = std::chrono::steady_clock::now();
    maybe_write_iscm_surface_charge_samples(*map_result, iscm_audit);
    maybe_write_rpb_source_audit(*map_result, rpb_audit);
    maybe_write_rpb_source_top(*map_result, solve.shape, *rpb_source_top_map,
                               *rpb_coulomb_top_map);
    maybe_write_barnase_forensic_audits(*map_result, atoms, solvent, vacuum,
                                        solve, solvent_forensic);
    maybe_write_parity_map_dumps(*map_result, atoms, solvent, vacuum, solve);
    add_timing_ms(TimingBucket::AuditOutput,
                  elapsed_ms(audit_output_begin,
                             std::chrono::steady_clock::now()));
  }
  return solve;
}

AnalyzerImpl::ReactionSolve AnalyzerImpl::solve_tabi_reaction_field(
    const std::vector<Atom> &atoms, TitanExecutor &executor) const {
  ReactionSolve solve;
  TabiPrepareProgressWriter prepare_progress(args_, atoms.size(),
                                             configured_threads_);
  prepare_progress.write("tabi_prepare_start", "running");
  std::vector<Tabi::Atom> tabi_atoms;
  tabi_atoms.reserve(atoms.size());
  double min_x = std::numeric_limits<double>::infinity();
  double min_y = std::numeric_limits<double>::infinity();
  double min_z = std::numeric_limits<double>::infinity();
  double max_x = -std::numeric_limits<double>::infinity();
  double max_y = -std::numeric_limits<double>::infinity();
  double max_z = -std::numeric_limits<double>::infinity();
  for (const Atom &atom : atoms) {
    Tabi::Atom tabi_atom;
    tabi_atom.x = atom.x;
    tabi_atom.y = atom.y;
    tabi_atom.z = atom.z;
    tabi_atom.q = atom.q;
    tabi_atom.radius = effective_surface_radius(atom);
    tabi_atoms.push_back(tabi_atom);
    const double surface_radius =
        tabi_atom.radius + (args_.surface == "sas" ? args_.probe_radius : 0.0);
    min_x = std::min(min_x, atom.x - surface_radius);
    min_y = std::min(min_y, atom.y - surface_radius);
    min_z = std::min(min_z, atom.z - surface_radius);
    max_x = std::max(max_x, atom.x + surface_radius);
    max_y = std::max(max_y, atom.y + surface_radius);
    max_z = std::max(max_z, atom.z + surface_radius);
  }
  prepare_progress.write("pqr_loaded", "running");

  Tabi::Options options;
  options.eps_in = args_.int_diel;
  options.eps_out = args_.ext_diel;
  options.salt_m = args_.salt;
  options.temperature_k = args_.temp;
  options.max_iter = args_.max_iter;
  options.residual_tol = args_.rms_tol;
  options.gmres_restart = args_.tabi_gmres_restart;
  options.tree_preset = args_.tabi_tree_preset;
  options.sum_mode = args_.tabi_sum;
  options.probe_radius = args_.probe_radius;
  options.surface_mode = args_.surface;
  options.surface_geometry = args_.tabi_geometry;
  options.diagonal_mode = args_.tabi_diagonal;
  options.surface_point_count = args_.tabi_surface_points;
  options.highres_surface = args_.tabi_highres_surface;
  options.sphere_mesh_mode = args_.tabi_sphere_mesh;
  options.tree_leaf_size = args_.tabi_tree_leaf_size;
  options.tree_opening_theta = args_.tabi_tree_opening_theta;
  options.taylor_order = args_.tabi_taylor_order;
  options.derivative_mode = args_.tabi_derivatives;
  options.formulation = args_.tabi_formulation;
  options.energy_projection = args_.tabi_energy_projection;
  options.preconditioner = args_.tabi_precond;
  options.preconditioner_block_near = args_.tabi_precond_block_near;
  options.preconditioner_block_panel_cap = args_.tabi_precond_block_panel_cap;
  options.preconditioner_overlap_panel_cap =
      args_.tabi_precond_overlap_panel_cap;
  options.preconditioner_ilu_neighbor_cap =
      args_.tabi_precond_ilu_neighbor_cap;
  options.preconditioner_ilu_drop_tol = args_.tabi_precond_ilu_drop_tol;
  options.preconditioner_memory_budget_mb = args_.tabi_precond_memory_budget_mb;
  options.preconditioner_coarse_aggregate_cap =
      args_.tabi_precond_coarse_aggregate_cap;
  options.preconditioner_coarse_base = args_.tabi_precond_coarse_base;
  options.preconditioner_coarse_verify = args_.tabi_precond_coarse_verify;
  options.near_integration = args_.tabi_near_integration;
  options.quad_buffer = args_.tabi_quad_buffer;
  options.leaf_near_policy = args_.tabi_leaf_near_policy;
  options.leaf_panel_kernel = args_.tabi_leaf_panel_kernel;
  options.adaptive_level = args_.tabi_adaptive_level;
  options.cluster_policy = args_.tabi_cluster_policy;
  options.leaf_pair_kernel = args_.tabi_leaf_pair_kernel;
  options.leaf_kernel_profile = args_.tabi_leaf_kernel_profile;
  options.leaf_work_profile = args_.tabi_leaf_work_profile;
  options.matvec_cost_profile = args_.tabi_matvec_cost_profile;
  options.jcp_direct_kernel = args_.tabi_jcp_direct_kernel;
  options.leaf_pair_selection_profile =
      args_.tabi_leaf_pair_selection_profile;
  options.leaf_pair_probe_cache = args_.tabi_leaf_pair_probe_cache;
  options.leaf_pair_target_group = args_.tabi_leaf_pair_target_group;
  options.leaf_pair_directness_probe =
      args_.tabi_leaf_pair_directness_probe;
  options.leaf_pair_directness_cache =
      args_.tabi_leaf_pair_directness_cache;
  options.leaf_pair_directness_verify =
      args_.tabi_leaf_pair_directness_verify;
  options.leaf_pair_scheduler = args_.tabi_leaf_pair_scheduler;
  options.leaf_pair_worklist_max_active =
      args_.tabi_leaf_pair_worklist_max_active;
  options.leaf_pair_worklist_verify =
      args_.tabi_leaf_pair_worklist_verify;
  options.leaf_block_verify = args_.tabi_leaf_block_verify;
  options.leaf_mac = args_.tabi_leaf_mac;
  options.fmm_m2l = args_.tabi_fmm_m2l;
  options.target_batch_size = args_.tabi_target_batch_size;
  options.leaf_pair_min_source_panels =
      args_.tabi_leaf_pair_min_source_panels;
  options.leaf_pair_min_targets = args_.tabi_leaf_pair_min_targets;
  options.leaf_pair_max_lookahead =
      args_.tabi_leaf_pair_max_lookahead;
  options.interaction_mode = args_.tabi_interaction_mode;
  options.node_bounds = args_.tabi_node_bounds;
  options.tree_builder = args_.tabi_tree_builder;
  options.expansion_center = args_.tabi_expansion_center;
  options.fragmentation_safety = args_.tabi_fragmentation_safety;
  options.mac_safety = args_.tabi_mac_safety;
  options.mac_gap_factor = args_.tabi_mac_gap_factor;
  options.moment_audit = args_.tabi_moment_audit;
  options.moment_construction = args_.tabi_moment_construction;
  options.moment_sum = args_.tabi_moment_sum;
  options.convergence_audit = args_.tabi_convergence_audit;
  options.convergence_safety = args_.tabi_convergence_safety;
  options.convergence_rho_cap = args_.tabi_convergence_rho_cap;
  options.mac_aspect_cap = args_.tabi_mac_aspect_cap;
  options.convergence_estimator = args_.tabi_convergence_estimator;
  options.convergence_tail_threshold =
      args_.tabi_convergence_tail_threshold;
  options.convergence_estimator_error_cap =
      args_.tabi_convergence_estimator_error_cap;
  options.convergence_unstable_rho_cap =
      args_.tabi_convergence_unstable_rho_cap;
  options.convergence_calibration = args_.tabi_convergence_calibration;
  options.convergence_dielectric_amp =
      args_.tabi_convergence_dielectric_amp;
  options.convergence_estimator_model =
      args_.tabi_convergence_estimator_model;
  options.convergence_channel_calibration =
      args_.tabi_convergence_channel_calibration;
  options.convergence_split_policy = args_.tabi_convergence_split_policy;
  options.convergence_split_selectivity =
      args_.tabi_convergence_split_selectivity;
  options.convergence_split_absolute_floor =
      args_.tabi_convergence_split_absolute_floor;
  options.convergence_overreject_audit =
      args_.tabi_convergence_overreject_audit;
  options.convergence_overreject_detail =
      args_.tabi_convergence_overreject_detail;
  options.convergence_shape_audit = args_.tabi_convergence_shape_audit;
  options.convergence_kernel_response_audit =
      args_.tabi_convergence_kernel_response_audit;
  options.channel_missed_semantics = args_.tabi_channel_missed_semantics;
  options.accepted_miss_diagnostics =
      args_.tabi_accepted_miss_diagnostics;
  options.accepted_miss_absolute_floor =
      args_.tabi_accepted_miss_absolute_floor;
  options.convergence_k3_accepted_rho_cap =
      args_.tabi_convergence_k3_accepted_rho_cap;
  options.mac_distance_policy = args_.tabi_mac_distance_policy;
  options.mac_surface_margin_factor =
      args_.tabi_mac_surface_margin_factor;
  options.mac_surface_margin_policy =
      args_.tabi_mac_surface_margin_policy;
  options.mac_surface_margin_aspect_threshold =
      args_.tabi_mac_surface_margin_aspect_threshold;
  options.mac_surface_margin_rho_trigger =
      args_.tabi_mac_surface_margin_rho_trigger;
  options.mac_surface_margin_soft_aspect_threshold =
      args_.tabi_mac_surface_margin_soft_aspect_threshold;
  options.mac_surface_margin_soft_rho_trigger =
      args_.tabi_mac_surface_margin_soft_rho_trigger;
  options.mac_surface_margin_soft_factor =
      args_.tabi_mac_surface_margin_soft_factor;
  options.expansion_audit = args_.tabi_expansion_audit;
  options.expansion_audit_channel = args_.tabi_expansion_audit_channel;
  options.expansion_audit_sample_count =
      args_.tabi_expansion_audit_sample_count;
  options.expansion_audit_max_panels =
      args_.tabi_expansion_audit_max_panels;
  options.matvec_audit = args_.tabi_matvec_audit;
  options.matvec_audit_sample_count =
      args_.tabi_matvec_audit_sample_count;
  options.fmm_cartesian_audit = args_.tabi_fmm_cartesian_audit;
  options.fmm_local = args_.tabi_fmm_local;
  options.fmm_local_refine = args_.tabi_fmm_local_refine;
  options.fmm_local_refine_tol = args_.tabi_fmm_local_refine_tol;
  options.fmm_local_cache = args_.tabi_fmm_local_cache;
  options.fmm_local_cache_budget_mb = args_.tabi_fmm_local_cache_budget_mb;
  options.fmm_local_fallback = args_.tabi_fmm_local_fallback;
  options.fmm_local_plan_builder = args_.tabi_fmm_local_plan_builder;
  options.true_residual_probe = args_.tabi_true_residual_probe;
  options.true_residual_probe_scope = args_.tabi_true_residual_probe_scope;
  options.true_residual_probe_interval =
      args_.tabi_true_residual_probe_interval;
  options.true_residual_sample_count = args_.tabi_true_residual_sample_count;
  options.progress = args_.tabi_progress;
  if (options.progress == "auto" && args_.timing != "detail")
    options.progress = "off";
  options.progress_interval_seconds = args_.tabi_progress_interval;
  prepare_progress.write("tabi_options_built", "running");
  prepare_progress.write("solver_enter", "running");
  const Tabi::Result tabi = Tabi::solve_linear_pb(tabi_atoms, options, executor);

  solve.pbe_backend = "tabi";
  solve.tabi_surface_build_ms = tabi.surface_build_ms;
  solve.tabi_surface_panels = tabi.panel_count;
  solve.tabi_surface_point_count_requested =
      tabi.surface_point_count_requested;
  solve.tabi_panel_density_per_a2 = tabi.panel_density_per_a2;
  solve.tabi_panel_budget_status = tabi.panel_budget_status;
  solve.tabi_surface_exposed_points = tabi.exposed_points;
  solve.tabi_surface_buried_points = tabi.buried_points;
  solve.tabi_tree_build_ms = tabi.tree_build_ms;
  solve.tabi_tree_nodes = tabi.tree_nodes;
  solve.tabi_tree_leaf_count = tabi.tree_leaf_count;
  solve.tabi_near_interactions = tabi.near_interactions;
  solve.tabi_far_interactions = tabi.far_interactions;
  solve.tabi_matvec_mode = tabi.matvec_mode;
  solve.tabi_fmm_m2l_mode = tabi.fmm_m2l_mode;
  solve.tabi_fmm_m2l_expansion_order = tabi.fmm_m2l_expansion_order;
  solve.tabi_fmm_m2l_near_pairs = tabi.fmm_m2l_near_pairs;
  solve.tabi_fmm_m2l_far_pairs = tabi.fmm_m2l_far_pairs;
  solve.tabi_fmm_m2l_setup_ms = tabi.fmm_m2l_setup_ms;
  solve.tabi_fmm_m2l_apply_ms = tabi.fmm_m2l_apply_ms;
  solve.tabi_fmm_m2l_finite_outputs = tabi.fmm_m2l_finite_outputs;
  solve.tabi_status = tabi.status;
  solve.tabi_convergence_reason = tabi.convergence_reason;
  solve.tabi_residual_history = tabi.residual_history;
  solve.tabi_gmres_restart_requested = tabi.gmres_restart_requested;
  solve.tabi_gmres_restart_effective = tabi.gmres_restart_effective;
  solve.tabi_best_residual = tabi.best_residual;
  solve.tabi_best_residual_iter = tabi.best_residual_iter;
  solve.tabi_final_residual = tabi.final_residual;
  solve.tabi_best_preconditioned_residual =
      tabi.best_preconditioned_residual;
  solve.tabi_final_preconditioned_residual =
      tabi.final_preconditioned_residual;
  solve.tabi_true_residual_abs = tabi.true_residual_abs;
  solve.tabi_true_residual_rel = tabi.true_residual_rel;
  solve.tabi_true_residual_probe_mode = tabi.true_residual_probe_mode;
  solve.tabi_true_residual_probe_scope = tabi.true_residual_probe_scope;
  solve.tabi_true_residual_probe_interval =
      tabi.true_residual_probe_interval;
  solve.tabi_true_residual_sample_count = tabi.true_residual_sample_count;
  solve.tabi_true_residual_probe_calls = tabi.true_residual_probe_calls;
  solve.tabi_true_residual_probe_ms = tabi.true_residual_probe_ms;
  solve.tabi_sampled_true_residual_rel = tabi.sampled_true_residual_rel;
  solve.tabi_sampled_true_residual_rows = tabi.sampled_true_residual_rows;
  solve.tabi_sampled_probe_calls = tabi.sampled_probe_calls;
  solve.tabi_sampled_probe_ms = tabi.sampled_probe_ms;
  solve.tabi_full_true_residual_confirm_calls =
      tabi.full_true_residual_confirm_calls;
  solve.tabi_full_probe_ms = tabi.full_true_residual_confirm_ms;
  solve.tabi_reaction_energy_source = tabi.reaction_energy_source;
  solve.tabi_rhs_ms = tabi.rhs_ms;
  solve.tabi_gmres_ms = tabi.gmres_ms;
  solve.tabi_matvec_ms_total = tabi.matvec_ms_total;
  solve.tabi_matvec_ms_avg = tabi.matvec_ms_avg;
  solve.tabi_matvec_calls = tabi.matvec_calls;
  solve.tabi_tree_vs_direct_matvec_rel =
      tabi.tree_vs_direct_matvec_rel;
  solve.tabi_matvec_audit_rows = tabi.matvec_audit_rows;
  solve.tabi_matvec_audit_ms = tabi.matvec_audit_ms;
  solve.tabi_orthogonalization_ms = tabi.orthogonalization_ms;
  solve.tabi_preconditioner_ms = tabi.preconditioner_ms;
  solve.tabi_surface_neighbor_ms = tabi.surface_neighbor_ms;
  solve.tabi_surface_mesh_ms = tabi.surface_mesh_ms;
  solve.tabi_solid_angle_ms = tabi.solid_angle_ms;
  solve.tabi_tree_moment_ms = tabi.tree_moment_ms;
  solve.tabi_tree_traversal_ms = tabi.tree_traversal_ms;
  solve.tabi_matvec_cost_profile = tabi.matvec_cost_profile;
  solve.tabi_leaf_panels_per_second = tabi.leaf_panels_per_second;
  solve.tabi_clean_leaf_panels_per_second =
      tabi.clean_leaf_panels_per_second;
  solve.tabi_tree_traversal_fraction = tabi.tree_traversal_fraction;
  solve.tabi_near_eval_ms = tabi.near_eval_ms;
  solve.tabi_far_eval_ms = tabi.far_eval_ms;
  solve.tabi_singular_near_eval_ms = tabi.singular_near_eval_ms;
  solve.tabi_leaf_near_eval_ms = tabi.leaf_near_eval_ms;
  solve.tabi_leaf_panel_eval_ms = tabi.leaf_panel_eval_ms;
  solve.tabi_leaf_tile_eval_ms = tabi.leaf_tile_eval_ms;
  solve.tabi_leaf_residual_eval_ms = tabi.leaf_residual_eval_ms;
  solve.tabi_smooth_leaf_node_eval_ms = tabi.smooth_leaf_node_eval_ms;
  solve.tabi_singular_near_interactions = tabi.singular_near_interactions;
  solve.tabi_leaf_near_interactions = tabi.leaf_near_interactions;
  solve.tabi_far_node_interactions = tabi.far_node_interactions;
  solve.tabi_leaf_panel_interactions = tabi.leaf_panel_interactions;
  solve.tabi_leaf_tile_interactions = tabi.leaf_tile_interactions;
  solve.tabi_leaf_tile_panel_interactions =
      tabi.leaf_tile_panel_interactions;
  solve.tabi_leaf_residual_panel_interactions =
      tabi.leaf_residual_panel_interactions;
  solve.tabi_smooth_leaf_node_interactions =
      tabi.smooth_leaf_node_interactions;
  solve.tabi_leaf_pruned_interactions = tabi.leaf_pruned_interactions;
  solve.tabi_leaf_mac_candidates = tabi.leaf_mac_candidates;
  solve.tabi_leaf_mac_accepted = tabi.leaf_mac_accepted;
  solve.tabi_leaf_mac_rejected = tabi.leaf_mac_rejected;
  solve.tabi_leaf_mac_far_panels = tabi.leaf_mac_far_panels;
  solve.tabi_leaf_mac_direct_panels_saved_estimate =
      tabi.leaf_mac_direct_panels_saved_estimate;
  solve.tabi_leaf_work_profile = tabi.leaf_work_profile;
  solve.tabi_leaf_work_mac_terminal_panels =
      tabi.leaf_work_mac_terminal_panels;
  solve.tabi_leaf_work_physical_gap_panels =
      tabi.leaf_work_physical_gap_panels;
  solve.tabi_leaf_work_fragmentation_panels =
      tabi.leaf_work_fragmentation_panels;
  solve.tabi_leaf_work_convergence_panels =
      tabi.leaf_work_convergence_panels;
  solve.tabi_leaf_work_convergence_estimator_panels =
      tabi.leaf_work_convergence_estimator_panels;
  solve.tabi_leaf_work_convergence_rho_panels =
      tabi.leaf_work_convergence_rho_panels;
  solve.tabi_leaf_work_convergence_aspect_panels =
      tabi.leaf_work_convergence_aspect_panels;
  solve.tabi_leaf_work_surface_margin_panels =
      tabi.leaf_work_surface_margin_panels;
  solve.tabi_leaf_work_surface_margin_hard_panels =
      tabi.leaf_work_surface_margin_hard_panels;
  solve.tabi_leaf_work_surface_margin_soft_panels =
      tabi.leaf_work_surface_margin_soft_panels;
  solve.tabi_leaf_work_k3_rho_cap_panels =
      tabi.leaf_work_k3_rho_cap_panels;
  solve.tabi_leaf_work_singular_fallback_panels =
      tabi.leaf_work_singular_fallback_panels;
  solve.tabi_leaf_work_scalar_tail_panels =
      tabi.leaf_work_scalar_tail_panels;
  solve.tabi_leaf_work_profile_accounted_panels =
      tabi.leaf_work_profile_accounted_panels;
  solve.tabi_leaf_work_mac_terminal_frac =
      tabi.leaf_work_mac_terminal_frac;
  solve.tabi_leaf_work_convergence_frac =
      tabi.leaf_work_convergence_frac;
  solve.tabi_leaf_work_surface_margin_frac =
      tabi.leaf_work_surface_margin_frac;
  solve.tabi_leaf_work_singular_fallback_frac =
      tabi.leaf_work_singular_fallback_frac;
  solve.tabi_leaf_work_scalar_tail_frac =
      tabi.leaf_work_scalar_tail_frac;
  solve.tabi_adaptive_smooth_leaf_accepted =
      tabi.adaptive_smooth_leaf_accepted;
  solve.tabi_adaptive_smooth_leaf_rejected =
      tabi.adaptive_smooth_leaf_rejected;
  solve.tabi_adaptive_span_tile_accepted =
      tabi.adaptive_span_tile_accepted;
  solve.tabi_adaptive_span_tile_rejected =
      tabi.adaptive_span_tile_rejected;
  solve.tabi_adaptive_exact_leaf_fallbacks =
      tabi.adaptive_exact_leaf_fallbacks;
  solve.tabi_adaptive_reject_one_ring_overlap =
      tabi.adaptive_reject_one_ring_overlap;
  solve.tabi_adaptive_reject_error_bound_failed =
      tabi.adaptive_reject_error_bound_failed;
  solve.tabi_adaptive_reject_p4_check_failed =
      tabi.adaptive_reject_p4_check_failed;
  solve.tabi_adaptive_reject_too_close =
      tabi.adaptive_reject_too_close;
  solve.tabi_leaf_work_reduction_ratio =
      tabi.leaf_work_reduction_ratio;
  solve.tabi_target_cluster_interactions =
      tabi.target_cluster_interactions;
  solve.tabi_target_cluster_targets = tabi.target_cluster_targets;
  solve.tabi_target_cluster_panels_pruned =
      tabi.target_cluster_panels_pruned;
  solve.tabi_target_cluster_rejected = tabi.target_cluster_rejected;
  solve.tabi_target_cluster_eval_ms = tabi.target_cluster_eval_ms;
  solve.tabi_leaf_pair_batch_interactions =
      tabi.leaf_pair_batch_interactions;
  solve.tabi_leaf_pair_batch_targets = tabi.leaf_pair_batch_targets;
  solve.tabi_leaf_pair_batch_source_panels =
      tabi.leaf_pair_batch_source_panels;
  solve.tabi_leaf_pair_batch_residual_panels =
      tabi.leaf_pair_batch_residual_panels;
  solve.tabi_leaf_pair_batch_fallbacks = tabi.leaf_pair_batch_fallbacks;
  solve.tabi_leaf_pair_batch_ms = tabi.leaf_pair_batch_ms;
  solve.tabi_leaf_pair_batch_kernel = tabi.leaf_pair_batch_kernel;
  solve.tabi_leaf_pair_batch_vectorized_sources =
      tabi.leaf_pair_batch_vectorized_sources;
  solve.tabi_leaf_pair_batch_vector_lanes =
      tabi.leaf_pair_batch_vector_lanes;
  solve.tabi_leaf_pair_batch_pairs_per_second =
      tabi.leaf_pair_batch_pairs_per_second;
  solve.tabi_leaf_kernel_profile = tabi.leaf_kernel_profile;
  solve.tabi_leaf_work_profile = tabi.leaf_work_profile;
  solve.tabi_jcp_direct_kernel = tabi.jcp_direct_kernel;
  solve.tabi_leaf_kernel_scalar_tasks = tabi.leaf_kernel_scalar_tasks;
  solve.tabi_leaf_kernel_batch_tasks = tabi.leaf_kernel_batch_tasks;
  solve.tabi_leaf_kernel_source_panel_evals =
      tabi.leaf_kernel_source_panel_evals;
  solve.tabi_leaf_kernel_scalar_ms = tabi.leaf_kernel_scalar_ms;
  solve.tabi_leaf_kernel_target_avx2_ms =
      tabi.leaf_kernel_target_avx2_ms;
  solve.tabi_leaf_kernel_tail_ms = tabi.leaf_kernel_tail_ms;
  solve.tabi_leaf_pair_batch_soa_sources =
      tabi.leaf_pair_batch_soa_sources;
  solve.tabi_leaf_pair_batch_soa_ms = tabi.leaf_pair_batch_soa_ms;
  solve.tabi_leaf_pair_batch_scalar_tail_ms =
      tabi.leaf_pair_batch_scalar_tail_ms;
  solve.tabi_leaf_pair_batch_clean_source_panels_per_second =
      tabi.leaf_pair_batch_clean_source_panels_per_second;
  solve.tabi_leaf_pair_selection_profile =
      tabi.leaf_pair_selection_profile;
  solve.tabi_leaf_pair_probe_cache = tabi.leaf_pair_probe_cache;
  solve.tabi_leaf_pair_target_group = tabi.leaf_pair_target_group;
  solve.tabi_leaf_pair_directness_probe = tabi.leaf_pair_directness_probe;
  solve.tabi_leaf_pair_directness_cache = tabi.leaf_pair_directness_cache;
  solve.tabi_leaf_pair_directness_verify = tabi.leaf_pair_directness_verify;
  solve.tabi_leaf_pair_max_lookahead = tabi.leaf_pair_max_lookahead;
  solve.tabi_leaf_pair_selection_candidate_scans =
      tabi.leaf_pair_selection_candidate_scans;
  solve.tabi_leaf_pair_selection_clean_rejects =
      tabi.leaf_pair_selection_clean_rejects;
  solve.tabi_leaf_pair_selection_probe_calls =
      tabi.leaf_pair_selection_probe_calls;
  solve.tabi_leaf_pair_selection_probe_node_visits =
      tabi.leaf_pair_selection_probe_node_visits;
  solve.tabi_leaf_pair_selection_probe_ms =
      tabi.leaf_pair_selection_probe_ms;
  solve.tabi_leaf_pair_selection_ms = tabi.leaf_pair_selection_ms;
  solve.tabi_leaf_pair_selection_probe_cache_hits =
      tabi.leaf_pair_selection_probe_cache_hits;
  solve.tabi_leaf_pair_selection_probe_cache_misses =
      tabi.leaf_pair_selection_probe_cache_misses;
  solve.tabi_leaf_pair_selection_insufficient_targets =
      tabi.leaf_pair_selection_insufficient_targets;
  solve.tabi_leaf_pair_selection_same_target_leaf_rejects =
      tabi.leaf_pair_selection_same_target_leaf_rejects;
  solve.tabi_leaf_pair_directness_ancestor_probe_calls =
      tabi.leaf_pair_directness_ancestor_probe_calls;
  solve.tabi_leaf_pair_directness_ancestor_probe_node_visits =
      tabi.leaf_pair_directness_ancestor_probe_node_visits;
  solve.tabi_leaf_pair_directness_ancestor_probe_ms =
      tabi.leaf_pair_directness_ancestor_probe_ms;
  solve.tabi_leaf_pair_directness_pair_cache_hits =
      tabi.leaf_pair_directness_pair_cache_hits;
  solve.tabi_leaf_pair_directness_pair_cache_misses =
      tabi.leaf_pair_directness_pair_cache_misses;
  solve.tabi_leaf_pair_directness_all_direct_proofs =
      tabi.leaf_pair_directness_all_direct_proofs;
  solve.tabi_leaf_pair_directness_mixed_fallbacks =
      tabi.leaf_pair_directness_mixed_fallbacks;
  solve.tabi_leaf_pair_directness_verify_checks =
      tabi.leaf_pair_directness_verify_checks;
  solve.tabi_leaf_pair_directness_verify_mismatches =
      tabi.leaf_pair_directness_verify_mismatches;
  solve.tabi_leaf_pair_scheduler = tabi.leaf_pair_scheduler;
  solve.tabi_leaf_pair_worklist_max_active =
      tabi.leaf_pair_worklist_max_active;
  solve.tabi_leaf_pair_worklist_verify = tabi.leaf_pair_worklist_verify;
  solve.tabi_leaf_pair_worklist_tasks = tabi.leaf_pair_worklist_tasks;
  solve.tabi_leaf_pair_worklist_groups = tabi.leaf_pair_worklist_groups;
  solve.tabi_leaf_pair_worklist_batched_groups =
      tabi.leaf_pair_worklist_batched_groups;
  solve.tabi_leaf_pair_worklist_scalar_flush_tasks =
      tabi.leaf_pair_worklist_scalar_flush_tasks;
  solve.tabi_leaf_pair_worklist_forced_flushes =
      tabi.leaf_pair_worklist_forced_flushes;
  solve.tabi_leaf_pair_worklist_active_groups_max =
      tabi.leaf_pair_worklist_active_groups_max;
  solve.tabi_leaf_pair_worklist_ms = tabi.leaf_pair_worklist_ms;
  solve.tabi_leaf_pair_worklist_verify_checks =
      tabi.leaf_pair_worklist_verify_checks;
  solve.tabi_leaf_pair_worklist_verify_mismatches =
      tabi.leaf_pair_worklist_verify_mismatches;
  solve.tabi_leaf_block_verify = tabi.leaf_block_verify;
  solve.tabi_leaf_block_target_leaves = tabi.leaf_block_target_leaves;
  solve.tabi_leaf_block_targets = tabi.leaf_block_targets;
  solve.tabi_leaf_block_source_groups = tabi.leaf_block_source_groups;
  solve.tabi_leaf_block_batched_groups = tabi.leaf_block_batched_groups;
  solve.tabi_leaf_block_scalar_tail_tasks =
      tabi.leaf_block_scalar_tail_tasks;
  solve.tabi_leaf_block_group_ms = tabi.leaf_block_group_ms;
  solve.tabi_leaf_block_flush_ms = tabi.leaf_block_flush_ms;
  solve.tabi_leaf_block_writeback_ms = tabi.leaf_block_writeback_ms;
  solve.tabi_leaf_block_verify_checks = tabi.leaf_block_verify_checks;
  solve.tabi_leaf_block_verify_mismatches =
      tabi.leaf_block_verify_mismatches;
  solve.tabi_tree_preset = tabi.tree_preset;
  solve.tabi_adaptive_level = tabi.adaptive_level;
  solve.tabi_cluster_policy = tabi.cluster_policy;
  solve.tabi_target_batch_size = tabi.target_batch_size;
  solve.tabi_leaf_pair_min_source_panels =
      tabi.leaf_pair_min_source_panels;
  solve.tabi_leaf_pair_min_targets = tabi.leaf_pair_min_targets;
  solve.tabi_leaf_near_policy = tabi.leaf_near_policy;
  solve.tabi_leaf_panel_kernel = tabi.leaf_panel_kernel;
  solve.tabi_interaction_mode = tabi.interaction_mode;
  solve.tabi_node_bounds_mode = tabi.node_bounds_mode;
  solve.tabi_tree_builder = tabi.tree_builder;
  solve.tabi_expansion_center = tabi.expansion_center;
  solve.tabi_fragmentation_safety = tabi.fragmentation_safety;
  solve.tabi_mac_safety = tabi.mac_safety;
  solve.tabi_mac_gap_factor = tabi.mac_gap_factor;
  solve.tabi_tree_depth_max = tabi.tree_depth_max;
  solve.tabi_mac_reject_radius_limited = tabi.mac_reject_radius_limited;
  solve.tabi_mac_accept_rate = tabi.mac_accept_rate;
  solve.tabi_node_radius_max = tabi.node_radius_max;
  solve.tabi_node_radius_mean = tabi.node_radius_mean;
  solve.tabi_node_aspect_ratio_max = tabi.node_aspect_ratio_max;
  solve.tabi_node_aspect_ratio_mean = tabi.node_aspect_ratio_mean;
  solve.tabi_leaf_panel_work_per_matvec =
      tabi.leaf_panel_work_per_matvec;
  solve.tabi_far_node_work_per_matvec =
      tabi.far_node_work_per_matvec;
  solve.tabi_far_geometric_unsafe_count =
      tabi.far_geometric_unsafe_count;
  solve.tabi_far_geometric_unsafe_panels =
      tabi.far_geometric_unsafe_panels;
  solve.tabi_far_geometric_min_gap = tabi.far_geometric_min_gap;
  solve.tabi_far_geometric_gap_ratio_min =
      tabi.far_geometric_gap_ratio_min;
  solve.tabi_far_geometric_gap_ratio_mean =
      tabi.far_geometric_gap_ratio_mean;
  solve.tabi_far_geometric_audit_samples =
      tabi.far_geometric_audit_samples;
  solve.tabi_mac_physical_gap_rejects =
      tabi.mac_physical_gap_rejects;
  solve.tabi_mac_convergence_rho_rejects =
      tabi.mac_convergence_rho_rejects;
  solve.tabi_mac_convergence_aspect_rejects =
      tabi.mac_convergence_aspect_rejects;
  solve.tabi_mac_convergence_rho_reject_panels =
      tabi.mac_convergence_rho_reject_panels;
  solve.tabi_mac_convergence_aspect_reject_panels =
      tabi.mac_convergence_aspect_reject_panels;
  solve.tabi_mac_surface_margin_rejects =
      tabi.mac_surface_margin_rejects;
  solve.tabi_mac_surface_margin_reject_panels =
      tabi.mac_surface_margin_reject_panels;
  solve.tabi_mac_surface_margin_applied =
      tabi.mac_surface_margin_applied;
  solve.tabi_mac_surface_margin_skipped =
      tabi.mac_surface_margin_skipped;
  solve.tabi_mac_surface_margin_adaptive_rejects =
      tabi.mac_surface_margin_adaptive_rejects;
  solve.tabi_mac_surface_margin_adaptive_reject_panels =
      tabi.mac_surface_margin_adaptive_reject_panels;
  solve.tabi_mac_surface_margin_skipped_accepted_far =
      tabi.mac_surface_margin_skipped_accepted_far;
  solve.tabi_mac_surface_margin_hard_applied =
      tabi.mac_surface_margin_hard_applied;
  solve.tabi_mac_surface_margin_soft_applied =
      tabi.mac_surface_margin_soft_applied;
  solve.tabi_mac_surface_margin_hard_rejects =
      tabi.mac_surface_margin_hard_rejects;
  solve.tabi_mac_surface_margin_hard_reject_panels =
      tabi.mac_surface_margin_hard_reject_panels;
  solve.tabi_mac_surface_margin_soft_rejects =
      tabi.mac_surface_margin_soft_rejects;
  solve.tabi_mac_surface_margin_soft_reject_panels =
      tabi.mac_surface_margin_soft_reject_panels;
  solve.tabi_fragmented_far_count = tabi.fragmented_far_count;
  solve.tabi_fragmented_far_panels = tabi.fragmented_far_panels;
  solve.tabi_fragmented_split_rejects = tabi.fragmented_split_rejects;
  solve.tabi_fragmentation_center_gap_ratio_max =
      tabi.fragmentation_center_gap_ratio_max;
  solve.tabi_fragmentation_center_gap_ratio_mean =
      tabi.fragmentation_center_gap_ratio_mean;
  solve.tabi_moment_audit = tabi.moment_audit;
  solve.tabi_moment_construction = tabi.moment_construction;
  solve.tabi_moment_sum = tabi.moment_sum;
  solve.tabi_convergence_audit = tabi.convergence_audit;
  solve.tabi_convergence_safety = tabi.convergence_safety;
  solve.tabi_convergence_rho_cap = tabi.convergence_rho_cap;
  solve.tabi_mac_aspect_cap = tabi.mac_aspect_cap;
  solve.tabi_convergence_estimator = tabi.convergence_estimator;
  solve.tabi_convergence_tail_threshold =
      tabi.convergence_tail_threshold;
  solve.tabi_convergence_estimator_error_cap =
      tabi.convergence_estimator_error_cap;
  solve.tabi_convergence_unstable_rho_cap =
      tabi.convergence_unstable_rho_cap;
  solve.tabi_convergence_calibration = tabi.convergence_calibration;
  solve.tabi_convergence_dielectric_amp =
      tabi.convergence_dielectric_amp;
  solve.tabi_convergence_estimator_model =
      tabi.convergence_estimator_model;
  solve.tabi_convergence_channel_calibration =
      tabi.convergence_channel_calibration;
  solve.tabi_convergence_split_policy = tabi.convergence_split_policy;
  solve.tabi_convergence_split_selectivity =
      tabi.convergence_split_selectivity;
  solve.tabi_convergence_split_absolute_floor =
      tabi.convergence_split_absolute_floor;
  solve.tabi_convergence_overreject_audit =
      tabi.convergence_overreject_audit;
  solve.tabi_convergence_overreject_detail =
      tabi.convergence_overreject_detail;
  solve.tabi_convergence_shape_audit = tabi.convergence_shape_audit;
  solve.tabi_convergence_kernel_response_audit =
      tabi.convergence_kernel_response_audit;
  solve.tabi_channel_missed_semantics = tabi.channel_missed_semantics;
  solve.tabi_accepted_miss_diagnostics =
      tabi.accepted_miss_diagnostics;
  solve.tabi_accepted_miss_absolute_floor =
      tabi.accepted_miss_absolute_floor;
  solve.tabi_convergence_k3_accepted_rho_cap =
      tabi.convergence_k3_accepted_rho_cap;
  solve.tabi_mac_distance_policy = tabi.mac_distance_policy;
  solve.tabi_mac_surface_margin_factor =
      tabi.mac_surface_margin_factor;
  solve.tabi_mac_surface_margin_policy =
      tabi.mac_surface_margin_policy;
  solve.tabi_mac_surface_margin_aspect_threshold =
      tabi.mac_surface_margin_aspect_threshold;
  solve.tabi_mac_surface_margin_rho_trigger =
      tabi.mac_surface_margin_rho_trigger;
  solve.tabi_mac_surface_margin_soft_aspect_threshold =
      tabi.mac_surface_margin_soft_aspect_threshold;
  solve.tabi_mac_surface_margin_soft_rho_trigger =
      tabi.mac_surface_margin_soft_rho_trigger;
  solve.tabi_mac_surface_margin_soft_factor =
      tabi.mac_surface_margin_soft_factor;
  solve.tabi_convergence_estimator_samples =
      tabi.convergence_estimator_samples;
  solve.tabi_convergence_estimator_unstable_far_count =
      tabi.convergence_estimator_unstable_far_count;
  solve.tabi_convergence_estimator_unstable_far_panels =
      tabi.convergence_estimator_unstable_far_panels;
  solve.tabi_convergence_estimator_split_rejects =
      tabi.convergence_estimator_split_rejects;
  solve.tabi_convergence_estimator_split_reject_panels =
      tabi.convergence_estimator_split_reject_panels;
  solve.tabi_convergence_low_abs_reaccepted_far_count =
      tabi.convergence_low_abs_reaccepted_far_count;
  solve.tabi_convergence_low_abs_reaccepted_far_panels =
      tabi.convergence_low_abs_reaccepted_far_panels;
  solve.tabi_convergence_estimator_tail_ratio_max =
      tabi.convergence_estimator_tail_ratio_max;
  solve.tabi_convergence_estimator_tail_ratio_mean =
      tabi.convergence_estimator_tail_ratio_mean;
  solve.tabi_convergence_estimator_error_max =
      tabi.convergence_estimator_error_max;
  solve.tabi_convergence_estimator_error_mean =
      tabi.convergence_estimator_error_mean;
  solve.tabi_convergence_estimator_abs_error_max =
      tabi.convergence_estimator_abs_error_max;
  solve.tabi_convergence_estimator_abs_error_mean =
      tabi.convergence_estimator_abs_error_mean;
  solve.tabi_convergence_tail_norm_p_max =
      tabi.convergence_tail_norm_p_max;
  solve.tabi_convergence_tail_norm_p_mean =
      tabi.convergence_tail_norm_p_mean;
  solve.tabi_convergence_tail_norm_pm1_max =
      tabi.convergence_tail_norm_pm1_max;
  solve.tabi_convergence_tail_norm_pm1_mean =
      tabi.convergence_tail_norm_pm1_mean;
  solve.tabi_convergence_effective_jcp_amp_max =
      tabi.convergence_effective_jcp_amp_max;
  solve.tabi_convergence_effective_jcp_amp_mean =
      tabi.convergence_effective_jcp_amp_mean;
  solve.tabi_convergence_dielectric_amp_max =
      tabi.convergence_dielectric_amp_max;
  solve.tabi_convergence_dielectric_amp_mean =
      tabi.convergence_dielectric_amp_mean;
  solve.tabi_convergence_channel_estimated_error_max =
      tabi.convergence_channel_estimated_error_max;
  solve.tabi_convergence_channel_estimated_error_mean =
      tabi.convergence_channel_estimated_error_mean;
  solve.tabi_convergence_channel_effective_jcp_amp_max =
      tabi.convergence_channel_effective_jcp_amp_max;
  solve.tabi_convergence_channel_effective_jcp_amp_mean =
      tabi.convergence_channel_effective_jcp_amp_mean;
  solve.tabi_convergence_channel_unstable_far_count =
      tabi.convergence_channel_unstable_far_count;
  solve.tabi_convergence_channel_unstable_far_panels =
      tabi.convergence_channel_unstable_far_panels;
  solve.tabi_convergence_channel_split_rejects =
      tabi.convergence_channel_split_rejects;
  solve.tabi_convergence_channel_split_reject_panels =
      tabi.convergence_channel_split_reject_panels;
  solve.tabi_convergence_channel_low_abs_reaccepted_count =
      tabi.convergence_channel_low_abs_reaccepted_count;
  solve.tabi_convergence_channel_low_abs_reaccepted_panels =
      tabi.convergence_channel_low_abs_reaccepted_panels;
  solve.tabi_convergence_channel_high_error_caught =
      tabi.convergence_channel_high_error_caught;
  solve.tabi_convergence_channel_high_error_missed =
      tabi.convergence_channel_high_error_missed;
  solve.tabi_convergence_channel_high_error_caught_by_split =
      tabi.convergence_channel_high_error_caught_by_split;
  solve.tabi_convergence_estimator_high_error_caught =
      tabi.convergence_estimator_high_error_caught;
  solve.tabi_convergence_estimator_high_error_missed =
      tabi.convergence_estimator_high_error_missed;
  solve.tabi_convergence_estimator_high_error_caught_by_split =
      tabi.convergence_estimator_high_error_caught_by_split;
  solve.tabi_convergence_estimator_phantom_high_error_count =
      tabi.convergence_estimator_phantom_high_error_count;
  solve.tabi_convergence_estimator_phantom_high_error_abs_max =
      tabi.convergence_estimator_phantom_high_error_abs_max;
  solve.tabi_convergence_overreject_count =
      tabi.convergence_overreject_count;
  solve.tabi_convergence_overreject_panels =
      tabi.convergence_overreject_panels;
  solve.tabi_convergence_overreject_caught_count =
      tabi.convergence_overreject_caught_count;
  solve.tabi_convergence_overreject_caught_panels =
      tabi.convergence_overreject_caught_panels;
  solve.tabi_convergence_overreject_abs_error_max =
      tabi.convergence_overreject_abs_error_max;
  solve.tabi_convergence_overreject_abs_error_mean =
      tabi.convergence_overreject_abs_error_mean;
  solve.tabi_convergence_overreject_abs_error_p95 =
      tabi.convergence_overreject_abs_error_p95;
  solve.tabi_convergence_overreject_abs_error_p99 =
      tabi.convergence_overreject_abs_error_p99;
  solve.tabi_convergence_overreject_estimated_abs_error_max =
      tabi.convergence_overreject_estimated_abs_error_max;
  solve.tabi_convergence_overreject_estimated_abs_error_mean =
      tabi.convergence_overreject_estimated_abs_error_mean;
  solve.tabi_convergence_overreject_required_scale_max =
      tabi.convergence_overreject_required_scale_max;
  solve.tabi_convergence_overreject_required_scale_p95 =
      tabi.convergence_overreject_required_scale_p95;
  solve.tabi_convergence_overreject_required_scale_p99 =
      tabi.convergence_overreject_required_scale_p99;
  solve.tabi_convergence_channel_overreject_count =
      tabi.convergence_channel_overreject_count;
  solve.tabi_convergence_channel_overreject_caught_count =
      tabi.convergence_channel_overreject_caught_count;
  solve.tabi_convergence_channel_overreject_panels =
      tabi.convergence_channel_overreject_panels;
  solve.tabi_convergence_channel_overreject_abs_error_max =
      tabi.convergence_channel_overreject_abs_error_max;
  solve.tabi_convergence_channel_overreject_abs_error_mean =
      tabi.convergence_channel_overreject_abs_error_mean;
  solve.tabi_convergence_channel_overreject_required_scale_max =
      tabi.convergence_channel_overreject_required_scale_max;
  solve.tabi_convergence_channel_overreject_required_scale_p99 =
      tabi.convergence_channel_overreject_required_scale_p99;
  solve.tabi_phantom_high_error_channel_count =
      tabi.phantom_high_error_channel_count;
  solve.tabi_phantom_high_error_channel_abs_max =
      tabi.phantom_high_error_channel_abs_max;
  solve.tabi_accepted_miss_channel_count =
      tabi.accepted_miss_channel_count;
  solve.tabi_accepted_miss_channel_panels =
      tabi.accepted_miss_channel_panels;
  solve.tabi_accepted_miss_channel_rho_min =
      tabi.accepted_miss_channel_rho_min;
  solve.tabi_accepted_miss_channel_rho_mean =
      tabi.accepted_miss_channel_rho_mean;
  solve.tabi_accepted_miss_channel_rho_max =
      tabi.accepted_miss_channel_rho_max;
  solve.tabi_accepted_miss_channel_tail_ratio_mean =
      tabi.accepted_miss_channel_tail_ratio_mean;
  solve.tabi_accepted_miss_channel_tail_ratio_max =
      tabi.accepted_miss_channel_tail_ratio_max;
  solve.tabi_accepted_miss_channel_estimated_error_mean =
      tabi.accepted_miss_channel_estimated_error_mean;
  solve.tabi_accepted_miss_channel_estimated_error_max =
      tabi.accepted_miss_channel_estimated_error_max;
  solve.tabi_accepted_miss_worst_channel =
      tabi.accepted_miss_worst_channel;
  solve.tabi_convergence_k3_rho_cap_rejects =
      tabi.convergence_k3_rho_cap_rejects;
  solve.tabi_convergence_k3_rho_cap_reject_panels =
      tabi.convergence_k3_rho_cap_reject_panels;
  solve.tabi_mac_surface_margin_rejects =
      tabi.mac_surface_margin_rejects;
  solve.tabi_mac_surface_margin_reject_panels =
      tabi.mac_surface_margin_reject_panels;
  solve.tabi_moment_audit_nodes = tabi.moment_audit_nodes;
  solve.tabi_moment_translation_rel_error_max =
      tabi.moment_translation_rel_error_max;
  solve.tabi_moment_translation_rel_error_mean =
      tabi.moment_translation_rel_error_mean;
  solve.tabi_moment_translation_offset_max =
      tabi.moment_translation_offset_max;
  solve.tabi_moment_translation_offset_mean =
      tabi.moment_translation_offset_mean;
  solve.tabi_moment_cancellation_ratio_max =
      tabi.moment_cancellation_ratio_max;
  solve.tabi_moment_cancellation_ratio_mean =
      tabi.moment_cancellation_ratio_mean;
  solve.tabi_moment_dynamic_range_max =
      tabi.moment_dynamic_range_max;
  solve.tabi_expansion_audit = tabi.expansion_audit;
  solve.tabi_expansion_audit_channel = tabi.expansion_audit_channel;
  solve.tabi_far_expansion_audit_pairs =
      tabi.far_expansion_audit_pairs;
  solve.tabi_far_expansion_audit_panels =
      tabi.far_expansion_audit_panels;
  solve.tabi_far_expansion_direct_ms =
      tabi.far_expansion_direct_ms;
  solve.tabi_far_expansion_rel_error_max =
      tabi.far_expansion_rel_error_max;
  solve.tabi_far_expansion_rel_error_mean =
      tabi.far_expansion_rel_error_mean;
  solve.tabi_far_expansion_rel_error_p95 =
      tabi.far_expansion_rel_error_p95;
  solve.tabi_far_expansion_rel_error_p99 =
      tabi.far_expansion_rel_error_p99;
  solve.tabi_far_expansion_channel_rel_error_max =
      tabi.far_expansion_channel_rel_error_max;
  solve.tabi_far_expansion_channel_rel_error_mean =
      tabi.far_expansion_channel_rel_error_mean;
  solve.tabi_far_expansion_channel_rel_error_p99 =
      tabi.far_expansion_channel_rel_error_p99;
  solve.tabi_far_expansion_worst_channel =
      tabi.far_expansion_worst_channel;
  solve.tabi_far_expansion_mean_gate_pass =
      tabi.far_expansion_mean_gate_pass;
  solve.tabi_far_convergence_rho_max =
      tabi.far_convergence_rho_max;
  solve.tabi_far_convergence_rho_mean =
      tabi.far_convergence_rho_mean;
  solve.tabi_accepted_far_aspect_ratio_max =
      tabi.accepted_far_aspect_ratio_max;
  solve.tabi_accepted_far_aspect_ratio_mean =
      tabi.accepted_far_aspect_ratio_mean;
  solve.tabi_quad_buffer_auto_candidate = tabi.quad_buffer_auto_candidate;
  solve.tabi_leaf_tile_status = tabi.leaf_tile_status;
  solve.tabi_leaf_panel_pairs_per_second =
      tabi.leaf_panel_pairs_per_second;
  solve.tabi_near_integration = tabi.near_integration;
  solve.tabi_quad_buffer = tabi.quad_buffer;
  solve.tabi_near_self_pairs = tabi.near_self_pairs;
  solve.tabi_near_one_ring_pairs = tabi.near_one_ring_pairs;
  solve.tabi_near_duffy_pairs = tabi.near_duffy_pairs;
  solve.tabi_near_centroid_pairs = tabi.near_centroid_pairs;
  solve.tabi_near_fallbacks = tabi.near_fallbacks;
  solve.tabi_near_geometric_buffer_pairs =
      tabi.near_geometric_buffer_pairs;
  solve.tabi_near_geometric_buffer_gauss_pairs =
      tabi.near_geometric_buffer_gauss_pairs;
  solve.tabi_near_geometric_buffer_adaptive_pairs =
      tabi.near_geometric_buffer_adaptive_pairs;
  solve.tabi_near_geometric_buffer_centroid_pairs =
      tabi.near_geometric_buffer_centroid_pairs;
  solve.tabi_near_geometric_buffer_min_gap =
      tabi.near_geometric_buffer_min_gap;
  solve.tabi_near_geometric_buffer_mean_gap_ratio =
      tabi.near_geometric_buffer_mean_gap_ratio;
  solve.tabi_near_geometric_buffer_max_edge_mean =
      tabi.near_geometric_buffer_max_edge_mean;
  solve.tabi_near_integration_ms = tabi.near_integration_ms;
  solve.tabi_duffy_subdivisions_max = tabi.duffy_subdivisions_max;
  solve.tabi_duffy_subdivisions_mean = tabi.duffy_subdivisions_mean;
  solve.tabi_near_kernel_policy = tabi.near_kernel_policy;
  solve.tabi_near_k1_accepted = tabi.near_k1_accepted;
  solve.tabi_near_k1_rejected = tabi.near_k1_rejected;
  solve.tabi_near_k2_accepted = tabi.near_k2_accepted;
  solve.tabi_near_k2_rejected = tabi.near_k2_rejected;
  solve.tabi_near_k3_accepted = tabi.near_k3_accepted;
  solve.tabi_near_k3_rejected = tabi.near_k3_rejected;
  solve.tabi_near_k4_accepted = tabi.near_k4_accepted;
  solve.tabi_near_k4_rejected = tabi.near_k4_rejected;
  solve.tabi_near_one_ring_k1_accepted =
      tabi.near_one_ring_k1_accepted;
  solve.tabi_near_one_ring_k1_rejected =
      tabi.near_one_ring_k1_rejected;
  solve.tabi_near_one_ring_k2_accepted =
      tabi.near_one_ring_k2_accepted;
  solve.tabi_near_one_ring_k2_rejected =
      tabi.near_one_ring_k2_rejected;
  solve.tabi_near_one_ring_k3_accepted =
      tabi.near_one_ring_k3_accepted;
  solve.tabi_near_one_ring_k3_rejected =
      tabi.near_one_ring_k3_rejected;
  solve.tabi_near_one_ring_k4_accepted =
      tabi.near_one_ring_k4_accepted;
  solve.tabi_near_one_ring_k4_rejected =
      tabi.near_one_ring_k4_rejected;
  solve.tabi_near_k1_policy = tabi.near_k1_policy;
  solve.tabi_near_k2_policy = tabi.near_k2_policy;
  solve.tabi_near_k3_policy = tabi.near_k3_policy;
  solve.tabi_near_k4_policy = tabi.near_k4_policy;
  for (std::size_t i = 0; i < solve.tabi_near_block_audit.size(); ++i) {
    const auto &src = tabi.near_block_audit[i];
    auto &dst = solve.tabi_near_block_audit[i];
    dst.accepted = src.accepted;
    dst.rejected = src.rejected;
    dst.candidate_worse_than_centroid = src.candidate_worse_than_centroid;
    dst.reference_nonconverged = src.reference_nonconverged;
    dst.orientation_sign_mismatch = src.orientation_sign_mismatch;
    dst.strict_ref_reject = src.strict_ref_reject;
    dst.max_candidate_ref_rel = src.max_candidate_ref_rel;
    dst.sum_candidate_ref_rel = src.sum_candidate_ref_rel;
    dst.max_centroid_ref_rel = src.max_centroid_ref_rel;
    dst.k2_laplace_sign_ok = src.k2_laplace_sign_ok;
    dst.max_k2_yukawa_smooth_ref_rel =
        src.max_k2_yukawa_smooth_ref_rel;
    dst.max_k2_orientation_pair_ref_rel =
        src.max_k2_orientation_pair_ref_rel;
    dst.sample_count = src.sample_count;
  }
  for (std::size_t i = 0; i < solve.tabi_geometric_near_audit.size(); ++i) {
    const auto &src = tabi.geometric_near_audit[i];
    auto &dst = solve.tabi_geometric_near_audit[i];
    dst.pairs = src.pairs;
    dst.min_gap = src.min_gap;
    dst.max_gap_ratio = src.max_gap_ratio;
    dst.sum_gap_ratio = src.sum_gap_ratio;
    dst.sum_center_ratio = src.sum_center_ratio;
  }
  solve.tabi_surface_geometry = tabi.surface_geometry;
  solve.tabi_sphere_mesh_mode = tabi.sphere_mesh_mode;
  solve.tabi_surface_quality_total_area_a2 =
      tabi.surface_quality_total_area_a2;
  solve.tabi_surface_quality_area_rel_error_vs_4piR2 =
      tabi.surface_quality_area_rel_error_vs_4piR2;
  solve.tabi_surface_quality_triangle_area_min =
      tabi.surface_quality_triangle_area_min;
  solve.tabi_surface_quality_triangle_area_mean =
      tabi.surface_quality_triangle_area_mean;
  solve.tabi_surface_quality_triangle_area_max =
      tabi.surface_quality_triangle_area_max;
  solve.tabi_surface_quality_triangle_aspect_max =
      tabi.surface_quality_triangle_aspect_max;
  solve.tabi_surface_quality_normal_radial_dot_min =
      tabi.surface_quality_normal_radial_dot_min;
  solve.tabi_surface_quality_normal_radial_dot_mean =
      tabi.surface_quality_normal_radial_dot_mean;
  solve.tabi_surface_quality_centroid_radial_error_max =
      tabi.surface_quality_centroid_radial_error_max;
  solve.tabi_surface_quality_planar_area_a2 =
      tabi.surface_quality_planar_area_a2;
  solve.tabi_surface_quality_planar_area_rel_error_vs_4piR2 =
      tabi.surface_quality_planar_area_rel_error_vs_4piR2;
  solve.tabi_surface_quality_enclosed_volume_a3 =
      tabi.surface_quality_enclosed_volume_a3;
  solve.tabi_surface_quality_enclosed_volume_rel_error_vs_4piR3_over3 =
      tabi.surface_quality_enclosed_volume_rel_error_vs_4piR3_over3;
  solve.tabi_surface_quality_triangle_area_cv =
      tabi.surface_quality_triangle_area_cv;
  solve.tabi_mesh_vertices = tabi.mesh_vertices;
  solve.tabi_mesh_triangles = tabi.mesh_triangles;
  solve.tabi_boundary_nodes = tabi.boundary_nodes;
  solve.tabi_avg_valence = tabi.avg_valence;
  solve.tabi_nonmanifold_edges = tabi.nonmanifold_edges;
  solve.tabi_solid_angle_fallbacks = tabi.solid_angle_fallbacks;
  solve.tabi_diag_min = tabi.diag_min;
  solve.tabi_diag_max = tabi.diag_max;
  solve.tabi_taylor_order = tabi.taylor_order;
  solve.tabi_derivative_mode = tabi.derivative_mode;
  solve.tabi_interaction_build_ms = tabi.interaction_build_ms;
  solve.tabi_tree_eval_ms = tabi.tree_eval_ms;
  solve.tabi_formulation = tabi.formulation;
  solve.tabi_energy_projection = tabi.energy_projection;
  solve.tabi_preconditioner_requested = tabi.preconditioner_requested;
  solve.tabi_preconditioner_effective = tabi.preconditioner_effective;
  solve.tabi_preconditioner_block_near = tabi.preconditioner_block_near;
  solve.tabi_preconditioner_block_panel_cap =
      tabi.preconditioner_block_panel_cap;
  solve.tabi_preconditioner_overlap_panel_cap =
      tabi.preconditioner_overlap_panel_cap;
  solve.tabi_preconditioner_ilu_neighbor_cap =
      tabi.preconditioner_ilu_neighbor_cap;
  solve.tabi_preconditioner_ilu_drop_tol = tabi.preconditioner_ilu_drop_tol;
  solve.tabi_preconditioner_memory_budget_bytes =
      tabi.preconditioner_memory_budget_bytes;
  solve.tabi_preconditioner_estimated_bytes =
      tabi.preconditioner_estimated_bytes;
  solve.tabi_preconditioner_valid_blocks = tabi.preconditioner_valid_blocks;
  solve.tabi_preconditioner_valid_panels = tabi.preconditioner_valid_panels;
  solve.tabi_preconditioner_ill_conditioned_blocks =
      tabi.preconditioner_ill_conditioned_blocks;
  solve.tabi_preconditioner_quality_rejected_blocks =
      tabi.preconditioner_quality_rejected_blocks;
  solve.tabi_preconditioner_quality_probe_max_gain =
      tabi.preconditioner_quality_probe_max_gain;
  solve.tabi_preconditioner_quality_status =
      tabi.preconditioner_quality_status;
  solve.tabi_preconditioner_setup_ms = tabi.preconditioner_setup_ms;
  solve.tabi_preconditioner_apply_ms = tabi.preconditioner_apply_ms;
  solve.tabi_preconditioner_apply_calls = tabi.preconditioner_apply_calls;
  solve.tabi_preconditioner_sparse_nnz_blocks =
      tabi.preconditioner_sparse_nnz_blocks;
  solve.tabi_preconditioner_sparse_avg_degree =
      tabi.preconditioner_sparse_avg_degree;
  solve.tabi_preconditioner_sparse_max_degree =
      tabi.preconditioner_sparse_max_degree;
  solve.tabi_preconditioner_coarse_aggregate_cap =
      tabi.preconditioner_coarse_aggregate_cap;
  solve.tabi_preconditioner_coarse_aggregates =
      tabi.preconditioner_coarse_aggregates;
  solve.tabi_preconditioner_coarse_dimension =
      tabi.preconditioner_coarse_dimension;
  solve.tabi_preconditioner_coarse_base = tabi.preconditioner_coarse_base;
  solve.tabi_preconditioner_coarse_matrix_bytes =
      tabi.preconditioner_coarse_matrix_bytes;
  solve.tabi_preconditioner_coarse_pivot_ratio =
      tabi.preconditioner_coarse_pivot_ratio;
  solve.tabi_preconditioner_coarse_status = tabi.preconditioner_coarse_status;
  solve.tabi_preconditioner_coarse_setup_matvec_calls =
      tabi.preconditioner_coarse_setup_matvec_calls;
  solve.tabi_preconditioner_coarse_setup_ms =
      tabi.preconditioner_coarse_setup_ms;
  solve.tabi_preconditioner_coarse_apply_ms =
      tabi.preconditioner_coarse_apply_ms;
  solve.tabi_preconditioner_coarse_apply_calls =
      tabi.preconditioner_coarse_apply_calls;
  solve.tabi_preconditioner_coarse_rp_identity_error_max =
      tabi.preconditioner_coarse_rp_identity_error_max;
  solve.tabi_preconditioner_coarse_verify_rel_error_max =
      tabi.preconditioner_coarse_verify_rel_error_max;
  solve.tabi_surface_generator = tabi.surface_generator;
  solve.tabi_discretization = tabi.discretization;
  solve.tabi_rhs_scale = tabi.rhs_scale;
  solve.tabi_dielectric_contrast = tabi.dielectric_contrast;
  solve.tabi_operator_contrast = tabi.operator_contrast;
  solve.tabi_unknown_count = tabi.unknown_count;
  solve.tabi_bie_system = tabi.bie_system;
  solve.tabi_lhs_diag_phi = tabi.lhs_diag_phi;
  solve.tabi_lhs_diag_phin = tabi.lhs_diag_phin;
  solve.tabi_energy_equation = tabi.energy_equation;
  solve.tabi_kernel_convention = tabi.kernel_convention;
  solve.tabi_source_norm_s1 = tabi.source_norm_s1;
  solve.tabi_source_norm_s2 = tabi.source_norm_s2;
  solve.tabi_phi_norm = tabi.phi_norm;
  solve.tabi_phin_norm = tabi.phin_norm;
  solve.tabi_kappa_inv_a = tabi.kappa_inv_a;
  solve.tabi_coulomb_factor_kcal_mol = tabi.coulomb_factor_kcal_mol;
  solve.tabi_boundary_energy_kj_mol = tabi.boundary_energy_kcal_mol * kKcalToKJ;
  solve.tabi_reference_energy_kj_mol =
      tabi.reference_energy_kcal_mol * kKcalToKJ;
  solve.tabi_selected_energy_scale = tabi.selected_energy_scale;
  solve.tabi_fallback_used = tabi.fallback_used;
  solve.shape.n = tabi.panel_count;
  solve.shape.nx = std::max<std::size_t>(1, tabi.panel_count);
  solve.shape.ny = 1;
  solve.shape.nz = 1;
  solve.shape.h = 0.0;
  if (!atoms.empty()) {
    solve.shape.ox = min_x;
    solve.shape.oy = min_y;
    solve.shape.oz = min_z;
    solve.focus_parent_grid =
        format_fixed(0.5 * (min_x + max_x), 6) + "," +
        format_fixed(0.5 * (min_y + max_y), 6) + "," +
        format_fixed(0.5 * (min_z + max_z), 6);
  }
  solve.shape.memory_bytes =
      tabi.peak_memory_bytes > 0 ? tabi.peak_memory_bytes : tabi.memory_bytes;
  solve.peak_memory_bytes = solve.shape.memory_bytes;
  solve.solute_points = tabi.panel_count;
  solve.solvent_points = 0;
  solve.solvent_total_iterations = tabi.iterations;
  solve.vacuum_total_iterations = 0;
  solve.solvent.iterations = tabi.iterations;
  solve.solvent.rms = tabi.residual;
  solve.solvent.max_dphi = tabi.residual;
  solve.solvent.converged = tabi.converged;
  solve.vacuum.converged = true;
  solve.pbe_model_effective = "linear";
  solve.linear_seed_solver = "none";
  solve.solver_mode = "tabi_" + tabi.sum_mode_effective;
  solve.mg_levels = 1;
  solve.focus_levels = 1;
  solve.boundary_rms = tabi.residual;
  solve.boundary_max = tabi.residual;
  solve.linear_residual = tabi.residual;
  solve.charge_assigned = tabi.total_charge;
  solve.energy_kcal_mol = tabi.energy_kcal_mol;
  solve.energy_kj_mol = tabi.energy_kj_mol;
  solve.energy_kt = tabi.energy_kt;
  solve.atom_corrected_kcal_mol = tabi.energy_kcal_mol;
  solve.atom_corrected_kj_mol = tabi.energy_kj_mol;
  solve.atom_corrected_kt = tabi.energy_kt;
  solve.grid_corrected_kcal_mol = tabi.energy_kcal_mol;
  solve.grid_corrected_kj_mol = tabi.energy_kj_mol;
  solve.grid_corrected_kt = tabi.energy_kt;
  solve.solvent_phase_kcal_mol = tabi.energy_kcal_mol;
  solve.solvent_phase_kj_mol = tabi.energy_kj_mol;
  solve.solvent_phase_kt = tabi.energy_kt;
  solve.coulomb_energy_kcal_mol =
      analytic_coulomb_energy_kcal(atoms, args_.int_diel);
  solve.coulomb_energy_kj_mol = solve.coulomb_energy_kcal_mol * kKcalToKJ;
  const double kt_kcal = kKBKcalPerMolK * std::max(1.0, args_.temp);
  solve.coulomb_energy_kt =
      kt_kcal > 0.0 ? solve.coulomb_energy_kcal_mol / kt_kcal : 0.0;
  solve.reaction_energy_source = tabi.reaction_energy_source;
  solve.converged = tabi.status == "ok";
  return solve;
}

void AnalyzerImpl::compute_stability_ensemble(const std::vector<Atom> &atoms,
                                              TitanExecutor &executor,
                                              const FrameResult &result) {
  last_stability_ensemble_ = StabilityEnsembleState{};
  if (args_.stability_mode != "ensemble" || atoms.empty() ||
      stability_probe_solve_active_)
    return;

  struct ProbeSolveScope {
    AnalyzerImpl &self;
    bool previous;
    explicit ProbeSolveScope(AnalyzerImpl &impl)
        : self(impl), previous(impl.stability_probe_solve_active_) {
      self.stability_probe_solve_active_ = true;
    }
    ~ProbeSolveScope() { self.stability_probe_solve_active_ = previous; }
  } probe_scope(*this);

  std::vector<double> energies;
  energies.reserve(7);
  energies.push_back(result.energy_kj_mol);
  constexpr double kJitterA = 1.0e-7;
  const std::array<std::array<double, 3>, 6> offsets{{
      {{kJitterA, 0.0, 0.0}},
      {{-kJitterA, 0.0, 0.0}},
      {{0.0, kJitterA, 0.0}},
      {{0.0, -kJitterA, 0.0}},
      {{0.0, 0.0, kJitterA}},
      {{0.0, 0.0, -kJitterA}},
  }};
  for (const auto &offset : offsets) {
    std::vector<Atom> shifted = atoms;
    for (Atom &atom : shifted) {
      atom.x += offset[0];
      atom.y += offset[1];
      atom.z += offset[2];
    }
    const ReactionSolve probe = solve_reaction_field(shifted, executor, nullptr);
    if (probe.converged && std::isfinite(probe.energy_kj_mol))
      energies.push_back(probe.energy_kj_mol);
  }
  if (energies.empty())
    return;
  const double sum =
      std::accumulate(energies.begin(), energies.end(), 0.0);
  const double mean = sum / static_cast<double>(energies.size());
  double m2 = 0.0;
  for (double value : energies) {
    const double delta = value - mean;
    m2 += delta * delta;
  }
  StabilityEnsembleState state;
  state.valid = true;
  state.samples = energies.size();
  state.mean_kj_mol = mean;
  state.stddev_kj_mol =
      energies.size() > 1
          ? std::sqrt(m2 / static_cast<double>(energies.size() - 1))
          : 0.0;
  const auto mm = std::minmax_element(energies.begin(), energies.end());
  state.min_kj_mol = *mm.first;
  state.max_kj_mol = *mm.second;
  last_stability_ensemble_ = state;
}

FrameResult AnalyzerImpl::consume_frame(const MOL &mol1,
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
  const ProcessMemorySample frame_memory_before = sample_process_memory();
  TimingBuckets frame_timing;
  struct ActiveTimingScope {
    AnalyzerImpl &self;
    TimingBuckets *previous;
    ActiveTimingScope(AnalyzerImpl &impl, TimingBuckets *next)
        : self(impl), previous(impl.active_timing_) {
      self.active_timing_ = next;
    }
    ~ActiveTimingScope() { self.active_timing_ = previous; }
  } timing_scope(*this, pbe_timing_enabled() ? &frame_timing : nullptr);
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
  result.pbe_backend = args_.backend;
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
  result.test_boundary_convention = effective_boundary_convention();
  result.test_field_units = args_.test_field_units;
  result.test_operator_convention =
      (args_.test_operator_convention == "native" && args_.focus == "delphi")
          ? std::string("delphi")
          : args_.test_operator_convention;
  result.charge_assignment = args_.charge_assignment;
  for (const Atom &atom : atoms_)
    result.total_charge += atom.q;
  last_stability_ensemble_ = StabilityEnsembleState{};
  stability_hysteresis_cells_kept_frame_ = 0;
  stability_hysteresis_cells_released_frame_ = 0;
  add_timing_ms(TimingBucket::Setup,
                elapsed_ms(frame_time_begin,
                           std::chrono::steady_clock::now()));

  if (args_.backend == "tabi") {
    const auto tabi_time_begin = std::chrono::steady_clock::now();
    const ReactionSolve tabi_solve =
        solve_tabi_reaction_field(atoms_, executor);
    const auto tabi_time_end = std::chrono::steady_clock::now();
    const double tabi_total_ms = elapsed_ms(tabi_time_begin, tabi_time_end);
    add_timing_ms(TimingBucket::Surface, tabi_solve.tabi_surface_build_ms);
    add_timing_ms(TimingBucket::Solver,
                  std::max(0.0,
                           tabi_total_ms - tabi_solve.tabi_surface_build_ms));

    result.nx = tabi_solve.shape.nx;
    result.ny = tabi_solve.shape.ny;
    result.nz = tabi_solve.shape.nz;
    result.grid_points = tabi_solve.shape.n;
    result.memory_bytes =
        std::max(tabi_solve.shape.memory_bytes, tabi_solve.peak_memory_bytes);
    result.grid_origin_x_A = tabi_solve.shape.ox;
    result.grid_origin_y_A = tabi_solve.shape.oy;
    result.grid_origin_z_A = tabi_solve.shape.oz;
    result.grid_spacing_A = tabi_solve.shape.h;
    result.grid_center_x_A = tabi_solve.shape.ox;
    result.grid_center_y_A = tabi_solve.shape.oy;
    result.grid_center_z_A = tabi_solve.shape.oz;
    result.solute_points = tabi_solve.solute_points;
    result.solvent_points = tabi_solve.solvent_points;
    result.ion_accessible_points = tabi_solve.ion_accessible_points;
    result.stern_excluded_points = tabi_solve.stern_excluded_points;
    result.tabi_surface_build_ms = tabi_solve.tabi_surface_build_ms;
    result.tabi_surface_panels = tabi_solve.tabi_surface_panels;
    result.tabi_surface_point_count_requested =
        tabi_solve.tabi_surface_point_count_requested;
    result.tabi_panel_density_per_a2 =
        tabi_solve.tabi_panel_density_per_a2;
    result.tabi_panel_budget_status =
        tabi_solve.tabi_panel_budget_status;
    result.tabi_surface_exposed_points =
        tabi_solve.tabi_surface_exposed_points;
    result.tabi_surface_buried_points =
        tabi_solve.tabi_surface_buried_points;
    result.tabi_tree_build_ms = tabi_solve.tabi_tree_build_ms;
    result.tabi_tree_nodes = tabi_solve.tabi_tree_nodes;
    result.tabi_tree_leaf_count = tabi_solve.tabi_tree_leaf_count;
    result.tabi_near_interactions = tabi_solve.tabi_near_interactions;
    result.tabi_far_interactions = tabi_solve.tabi_far_interactions;
    result.tabi_matvec_mode = tabi_solve.tabi_matvec_mode;
    result.tabi_fmm_m2l_mode = tabi_solve.tabi_fmm_m2l_mode;
    result.tabi_fmm_m2l_expansion_order =
        tabi_solve.tabi_fmm_m2l_expansion_order;
    result.tabi_fmm_m2l_near_pairs =
        tabi_solve.tabi_fmm_m2l_near_pairs;
    result.tabi_fmm_m2l_far_pairs = tabi_solve.tabi_fmm_m2l_far_pairs;
    result.tabi_fmm_m2l_setup_ms = tabi_solve.tabi_fmm_m2l_setup_ms;
    result.tabi_fmm_m2l_apply_ms = tabi_solve.tabi_fmm_m2l_apply_ms;
    result.tabi_fmm_m2l_finite_outputs =
        tabi_solve.tabi_fmm_m2l_finite_outputs;
    result.tabi_status = tabi_solve.tabi_status;
    result.tabi_convergence_reason = tabi_solve.tabi_convergence_reason;
    result.tabi_residual_history = tabi_solve.tabi_residual_history;
    result.tabi_gmres_restart_requested =
        tabi_solve.tabi_gmres_restart_requested;
    result.tabi_gmres_restart_effective =
        tabi_solve.tabi_gmres_restart_effective;
    result.tabi_best_residual = tabi_solve.tabi_best_residual;
    result.tabi_best_residual_iter = tabi_solve.tabi_best_residual_iter;
    result.tabi_final_residual = tabi_solve.tabi_final_residual;
    result.tabi_best_preconditioned_residual =
        tabi_solve.tabi_best_preconditioned_residual;
    result.tabi_final_preconditioned_residual =
        tabi_solve.tabi_final_preconditioned_residual;
    result.tabi_true_residual_abs = tabi_solve.tabi_true_residual_abs;
    result.tabi_true_residual_rel = tabi_solve.tabi_true_residual_rel;
    result.tabi_true_residual_probe_mode =
        tabi_solve.tabi_true_residual_probe_mode;
    result.tabi_true_residual_probe_scope =
        tabi_solve.tabi_true_residual_probe_scope;
    result.tabi_true_residual_probe_interval =
        tabi_solve.tabi_true_residual_probe_interval;
    result.tabi_true_residual_sample_count =
        tabi_solve.tabi_true_residual_sample_count;
    result.tabi_true_residual_probe_calls =
        tabi_solve.tabi_true_residual_probe_calls;
    result.tabi_true_residual_probe_ms =
        tabi_solve.tabi_true_residual_probe_ms;
    result.tabi_sampled_true_residual_rel =
        tabi_solve.tabi_sampled_true_residual_rel;
    result.tabi_sampled_true_residual_rows =
        tabi_solve.tabi_sampled_true_residual_rows;
    result.tabi_sampled_probe_calls = tabi_solve.tabi_sampled_probe_calls;
    result.tabi_sampled_probe_ms = tabi_solve.tabi_sampled_probe_ms;
    result.tabi_full_true_residual_confirm_calls =
        tabi_solve.tabi_full_true_residual_confirm_calls;
    result.tabi_full_probe_ms = tabi_solve.tabi_full_probe_ms;
    result.tabi_reaction_energy_source =
        tabi_solve.tabi_reaction_energy_source;
    result.tabi_rhs_ms = tabi_solve.tabi_rhs_ms;
    result.tabi_gmres_ms = tabi_solve.tabi_gmres_ms;
    result.tabi_matvec_ms_total = tabi_solve.tabi_matvec_ms_total;
    result.tabi_matvec_ms_avg = tabi_solve.tabi_matvec_ms_avg;
    result.tabi_matvec_calls = tabi_solve.tabi_matvec_calls;
    result.tabi_tree_vs_direct_matvec_rel =
        tabi_solve.tabi_tree_vs_direct_matvec_rel;
    result.tabi_matvec_audit_rows = tabi_solve.tabi_matvec_audit_rows;
    result.tabi_matvec_audit_ms = tabi_solve.tabi_matvec_audit_ms;
    result.tabi_orthogonalization_ms =
        tabi_solve.tabi_orthogonalization_ms;
    result.tabi_preconditioner_ms = tabi_solve.tabi_preconditioner_ms;
    result.tabi_surface_neighbor_ms = tabi_solve.tabi_surface_neighbor_ms;
    result.tabi_surface_mesh_ms = tabi_solve.tabi_surface_mesh_ms;
    result.tabi_solid_angle_ms = tabi_solve.tabi_solid_angle_ms;
    result.tabi_tree_moment_ms = tabi_solve.tabi_tree_moment_ms;
    result.tabi_tree_traversal_ms = tabi_solve.tabi_tree_traversal_ms;
    result.tabi_matvec_cost_profile =
        tabi_solve.tabi_matvec_cost_profile;
    result.tabi_leaf_panels_per_second =
        tabi_solve.tabi_leaf_panels_per_second;
    result.tabi_clean_leaf_panels_per_second =
        tabi_solve.tabi_clean_leaf_panels_per_second;
    result.tabi_tree_traversal_fraction =
        tabi_solve.tabi_tree_traversal_fraction;
    result.tabi_near_eval_ms = tabi_solve.tabi_near_eval_ms;
    result.tabi_far_eval_ms = tabi_solve.tabi_far_eval_ms;
    result.tabi_singular_near_eval_ms =
        tabi_solve.tabi_singular_near_eval_ms;
    result.tabi_leaf_near_eval_ms = tabi_solve.tabi_leaf_near_eval_ms;
    result.tabi_leaf_panel_eval_ms = tabi_solve.tabi_leaf_panel_eval_ms;
    result.tabi_leaf_tile_eval_ms = tabi_solve.tabi_leaf_tile_eval_ms;
    result.tabi_leaf_residual_eval_ms =
        tabi_solve.tabi_leaf_residual_eval_ms;
    result.tabi_smooth_leaf_node_eval_ms =
        tabi_solve.tabi_smooth_leaf_node_eval_ms;
    result.tabi_singular_near_interactions =
        tabi_solve.tabi_singular_near_interactions;
    result.tabi_leaf_near_interactions =
        tabi_solve.tabi_leaf_near_interactions;
    result.tabi_far_node_interactions =
        tabi_solve.tabi_far_node_interactions;
    result.tabi_leaf_panel_interactions =
        tabi_solve.tabi_leaf_panel_interactions;
    result.tabi_leaf_tile_interactions =
        tabi_solve.tabi_leaf_tile_interactions;
    result.tabi_leaf_tile_panel_interactions =
        tabi_solve.tabi_leaf_tile_panel_interactions;
    result.tabi_leaf_residual_panel_interactions =
        tabi_solve.tabi_leaf_residual_panel_interactions;
    result.tabi_smooth_leaf_node_interactions =
        tabi_solve.tabi_smooth_leaf_node_interactions;
    result.tabi_leaf_pruned_interactions =
        tabi_solve.tabi_leaf_pruned_interactions;
    result.tabi_leaf_mac_candidates = tabi_solve.tabi_leaf_mac_candidates;
    result.tabi_leaf_mac_accepted = tabi_solve.tabi_leaf_mac_accepted;
    result.tabi_leaf_mac_rejected = tabi_solve.tabi_leaf_mac_rejected;
    result.tabi_leaf_mac_far_panels = tabi_solve.tabi_leaf_mac_far_panels;
    result.tabi_leaf_mac_direct_panels_saved_estimate =
        tabi_solve.tabi_leaf_mac_direct_panels_saved_estimate;
    result.tabi_leaf_work_profile = tabi_solve.tabi_leaf_work_profile;
    result.tabi_leaf_work_mac_terminal_panels =
        tabi_solve.tabi_leaf_work_mac_terminal_panels;
    result.tabi_leaf_work_physical_gap_panels =
        tabi_solve.tabi_leaf_work_physical_gap_panels;
    result.tabi_leaf_work_fragmentation_panels =
        tabi_solve.tabi_leaf_work_fragmentation_panels;
    result.tabi_leaf_work_convergence_panels =
        tabi_solve.tabi_leaf_work_convergence_panels;
    result.tabi_leaf_work_convergence_estimator_panels =
        tabi_solve.tabi_leaf_work_convergence_estimator_panels;
    result.tabi_leaf_work_convergence_rho_panels =
        tabi_solve.tabi_leaf_work_convergence_rho_panels;
    result.tabi_leaf_work_convergence_aspect_panels =
        tabi_solve.tabi_leaf_work_convergence_aspect_panels;
    result.tabi_leaf_work_surface_margin_panels =
        tabi_solve.tabi_leaf_work_surface_margin_panels;
    result.tabi_leaf_work_surface_margin_hard_panels =
        tabi_solve.tabi_leaf_work_surface_margin_hard_panels;
    result.tabi_leaf_work_surface_margin_soft_panels =
        tabi_solve.tabi_leaf_work_surface_margin_soft_panels;
    result.tabi_leaf_work_k3_rho_cap_panels =
        tabi_solve.tabi_leaf_work_k3_rho_cap_panels;
    result.tabi_leaf_work_singular_fallback_panels =
        tabi_solve.tabi_leaf_work_singular_fallback_panels;
    result.tabi_leaf_work_scalar_tail_panels =
        tabi_solve.tabi_leaf_work_scalar_tail_panels;
    result.tabi_leaf_work_profile_accounted_panels =
        tabi_solve.tabi_leaf_work_profile_accounted_panels;
    result.tabi_leaf_work_mac_terminal_frac =
        tabi_solve.tabi_leaf_work_mac_terminal_frac;
    result.tabi_leaf_work_convergence_frac =
        tabi_solve.tabi_leaf_work_convergence_frac;
    result.tabi_leaf_work_surface_margin_frac =
        tabi_solve.tabi_leaf_work_surface_margin_frac;
    result.tabi_leaf_work_singular_fallback_frac =
        tabi_solve.tabi_leaf_work_singular_fallback_frac;
    result.tabi_leaf_work_scalar_tail_frac =
        tabi_solve.tabi_leaf_work_scalar_tail_frac;
    result.tabi_adaptive_smooth_leaf_accepted =
        tabi_solve.tabi_adaptive_smooth_leaf_accepted;
    result.tabi_adaptive_smooth_leaf_rejected =
        tabi_solve.tabi_adaptive_smooth_leaf_rejected;
    result.tabi_adaptive_span_tile_accepted =
        tabi_solve.tabi_adaptive_span_tile_accepted;
    result.tabi_adaptive_span_tile_rejected =
        tabi_solve.tabi_adaptive_span_tile_rejected;
    result.tabi_adaptive_exact_leaf_fallbacks =
        tabi_solve.tabi_adaptive_exact_leaf_fallbacks;
    result.tabi_adaptive_reject_one_ring_overlap =
        tabi_solve.tabi_adaptive_reject_one_ring_overlap;
    result.tabi_adaptive_reject_error_bound_failed =
        tabi_solve.tabi_adaptive_reject_error_bound_failed;
    result.tabi_adaptive_reject_p4_check_failed =
        tabi_solve.tabi_adaptive_reject_p4_check_failed;
    result.tabi_adaptive_reject_too_close =
        tabi_solve.tabi_adaptive_reject_too_close;
    result.tabi_leaf_work_reduction_ratio =
        tabi_solve.tabi_leaf_work_reduction_ratio;
    result.tabi_target_cluster_interactions =
        tabi_solve.tabi_target_cluster_interactions;
    result.tabi_target_cluster_targets =
        tabi_solve.tabi_target_cluster_targets;
    result.tabi_target_cluster_panels_pruned =
        tabi_solve.tabi_target_cluster_panels_pruned;
    result.tabi_target_cluster_rejected =
        tabi_solve.tabi_target_cluster_rejected;
    result.tabi_target_cluster_eval_ms =
        tabi_solve.tabi_target_cluster_eval_ms;
    result.tabi_leaf_pair_batch_interactions =
        tabi_solve.tabi_leaf_pair_batch_interactions;
    result.tabi_leaf_pair_batch_targets =
        tabi_solve.tabi_leaf_pair_batch_targets;
    result.tabi_leaf_pair_batch_source_panels =
        tabi_solve.tabi_leaf_pair_batch_source_panels;
    result.tabi_leaf_pair_batch_residual_panels =
        tabi_solve.tabi_leaf_pair_batch_residual_panels;
    result.tabi_leaf_pair_batch_fallbacks =
        tabi_solve.tabi_leaf_pair_batch_fallbacks;
    result.tabi_leaf_pair_batch_ms =
        tabi_solve.tabi_leaf_pair_batch_ms;
    result.tabi_leaf_pair_batch_kernel =
        tabi_solve.tabi_leaf_pair_batch_kernel;
    result.tabi_leaf_pair_batch_vectorized_sources =
        tabi_solve.tabi_leaf_pair_batch_vectorized_sources;
    result.tabi_leaf_pair_batch_vector_lanes =
        tabi_solve.tabi_leaf_pair_batch_vector_lanes;
    result.tabi_leaf_pair_batch_pairs_per_second =
        tabi_solve.tabi_leaf_pair_batch_pairs_per_second;
    result.tabi_leaf_kernel_profile = tabi_solve.tabi_leaf_kernel_profile;
    result.tabi_leaf_work_profile = tabi_solve.tabi_leaf_work_profile;
    result.tabi_jcp_direct_kernel = tabi_solve.tabi_jcp_direct_kernel;
    result.tabi_leaf_kernel_scalar_tasks =
        tabi_solve.tabi_leaf_kernel_scalar_tasks;
    result.tabi_leaf_kernel_batch_tasks =
        tabi_solve.tabi_leaf_kernel_batch_tasks;
    result.tabi_leaf_kernel_source_panel_evals =
        tabi_solve.tabi_leaf_kernel_source_panel_evals;
    result.tabi_leaf_kernel_scalar_ms =
        tabi_solve.tabi_leaf_kernel_scalar_ms;
    result.tabi_leaf_kernel_target_avx2_ms =
        tabi_solve.tabi_leaf_kernel_target_avx2_ms;
    result.tabi_leaf_kernel_tail_ms =
        tabi_solve.tabi_leaf_kernel_tail_ms;
    result.tabi_leaf_pair_batch_soa_sources =
        tabi_solve.tabi_leaf_pair_batch_soa_sources;
    result.tabi_leaf_pair_batch_soa_ms =
        tabi_solve.tabi_leaf_pair_batch_soa_ms;
    result.tabi_leaf_pair_batch_scalar_tail_ms =
        tabi_solve.tabi_leaf_pair_batch_scalar_tail_ms;
    result.tabi_leaf_pair_batch_clean_source_panels_per_second =
        tabi_solve.tabi_leaf_pair_batch_clean_source_panels_per_second;
    result.tabi_leaf_pair_selection_profile =
        tabi_solve.tabi_leaf_pair_selection_profile;
    result.tabi_leaf_pair_probe_cache =
        tabi_solve.tabi_leaf_pair_probe_cache;
    result.tabi_leaf_pair_target_group =
        tabi_solve.tabi_leaf_pair_target_group;
    result.tabi_leaf_pair_directness_probe =
        tabi_solve.tabi_leaf_pair_directness_probe;
    result.tabi_leaf_pair_directness_cache =
        tabi_solve.tabi_leaf_pair_directness_cache;
    result.tabi_leaf_pair_directness_verify =
        tabi_solve.tabi_leaf_pair_directness_verify;
    result.tabi_leaf_pair_max_lookahead =
        tabi_solve.tabi_leaf_pair_max_lookahead;
    result.tabi_leaf_pair_selection_candidate_scans =
        tabi_solve.tabi_leaf_pair_selection_candidate_scans;
    result.tabi_leaf_pair_selection_clean_rejects =
        tabi_solve.tabi_leaf_pair_selection_clean_rejects;
    result.tabi_leaf_pair_selection_probe_calls =
        tabi_solve.tabi_leaf_pair_selection_probe_calls;
    result.tabi_leaf_pair_selection_probe_node_visits =
        tabi_solve.tabi_leaf_pair_selection_probe_node_visits;
    result.tabi_leaf_pair_selection_probe_ms =
        tabi_solve.tabi_leaf_pair_selection_probe_ms;
    result.tabi_leaf_pair_selection_ms =
        tabi_solve.tabi_leaf_pair_selection_ms;
    result.tabi_leaf_pair_selection_probe_cache_hits =
        tabi_solve.tabi_leaf_pair_selection_probe_cache_hits;
    result.tabi_leaf_pair_selection_probe_cache_misses =
        tabi_solve.tabi_leaf_pair_selection_probe_cache_misses;
    result.tabi_leaf_pair_selection_insufficient_targets =
        tabi_solve.tabi_leaf_pair_selection_insufficient_targets;
    result.tabi_leaf_pair_selection_same_target_leaf_rejects =
        tabi_solve.tabi_leaf_pair_selection_same_target_leaf_rejects;
    result.tabi_leaf_pair_directness_ancestor_probe_calls =
        tabi_solve.tabi_leaf_pair_directness_ancestor_probe_calls;
    result.tabi_leaf_pair_directness_ancestor_probe_node_visits =
        tabi_solve.tabi_leaf_pair_directness_ancestor_probe_node_visits;
    result.tabi_leaf_pair_directness_ancestor_probe_ms =
        tabi_solve.tabi_leaf_pair_directness_ancestor_probe_ms;
    result.tabi_leaf_pair_directness_pair_cache_hits =
        tabi_solve.tabi_leaf_pair_directness_pair_cache_hits;
    result.tabi_leaf_pair_directness_pair_cache_misses =
        tabi_solve.tabi_leaf_pair_directness_pair_cache_misses;
    result.tabi_leaf_pair_directness_all_direct_proofs =
        tabi_solve.tabi_leaf_pair_directness_all_direct_proofs;
    result.tabi_leaf_pair_directness_mixed_fallbacks =
        tabi_solve.tabi_leaf_pair_directness_mixed_fallbacks;
    result.tabi_leaf_pair_directness_verify_checks =
        tabi_solve.tabi_leaf_pair_directness_verify_checks;
    result.tabi_leaf_pair_directness_verify_mismatches =
        tabi_solve.tabi_leaf_pair_directness_verify_mismatches;
    result.tabi_leaf_pair_scheduler =
        tabi_solve.tabi_leaf_pair_scheduler;
    result.tabi_leaf_pair_worklist_max_active =
        tabi_solve.tabi_leaf_pair_worklist_max_active;
    result.tabi_leaf_pair_worklist_verify =
        tabi_solve.tabi_leaf_pair_worklist_verify;
    result.tabi_leaf_pair_worklist_tasks =
        tabi_solve.tabi_leaf_pair_worklist_tasks;
    result.tabi_leaf_pair_worklist_groups =
        tabi_solve.tabi_leaf_pair_worklist_groups;
    result.tabi_leaf_pair_worklist_batched_groups =
        tabi_solve.tabi_leaf_pair_worklist_batched_groups;
    result.tabi_leaf_pair_worklist_scalar_flush_tasks =
        tabi_solve.tabi_leaf_pair_worklist_scalar_flush_tasks;
    result.tabi_leaf_pair_worklist_forced_flushes =
        tabi_solve.tabi_leaf_pair_worklist_forced_flushes;
    result.tabi_leaf_pair_worklist_active_groups_max =
        tabi_solve.tabi_leaf_pair_worklist_active_groups_max;
    result.tabi_leaf_pair_worklist_ms =
        tabi_solve.tabi_leaf_pair_worklist_ms;
    result.tabi_leaf_pair_worklist_verify_checks =
        tabi_solve.tabi_leaf_pair_worklist_verify_checks;
    result.tabi_leaf_pair_worklist_verify_mismatches =
        tabi_solve.tabi_leaf_pair_worklist_verify_mismatches;
    result.tabi_leaf_block_verify = tabi_solve.tabi_leaf_block_verify;
    result.tabi_leaf_block_target_leaves =
        tabi_solve.tabi_leaf_block_target_leaves;
    result.tabi_leaf_block_targets = tabi_solve.tabi_leaf_block_targets;
    result.tabi_leaf_block_source_groups =
        tabi_solve.tabi_leaf_block_source_groups;
    result.tabi_leaf_block_batched_groups =
        tabi_solve.tabi_leaf_block_batched_groups;
    result.tabi_leaf_block_scalar_tail_tasks =
        tabi_solve.tabi_leaf_block_scalar_tail_tasks;
    result.tabi_leaf_block_group_ms = tabi_solve.tabi_leaf_block_group_ms;
    result.tabi_leaf_block_flush_ms = tabi_solve.tabi_leaf_block_flush_ms;
    result.tabi_leaf_block_writeback_ms =
        tabi_solve.tabi_leaf_block_writeback_ms;
    result.tabi_leaf_block_verify_checks =
        tabi_solve.tabi_leaf_block_verify_checks;
    result.tabi_leaf_block_verify_mismatches =
        tabi_solve.tabi_leaf_block_verify_mismatches;
    result.tabi_tree_preset = tabi_solve.tabi_tree_preset;
    result.tabi_adaptive_level = tabi_solve.tabi_adaptive_level;
    result.tabi_cluster_policy = tabi_solve.tabi_cluster_policy;
    result.tabi_target_batch_size = tabi_solve.tabi_target_batch_size;
    result.tabi_leaf_pair_min_source_panels =
        tabi_solve.tabi_leaf_pair_min_source_panels;
    result.tabi_leaf_pair_min_targets =
        tabi_solve.tabi_leaf_pair_min_targets;
    result.tabi_leaf_near_policy = tabi_solve.tabi_leaf_near_policy;
    result.tabi_leaf_panel_kernel = tabi_solve.tabi_leaf_panel_kernel;
    result.tabi_interaction_mode = tabi_solve.tabi_interaction_mode;
    result.tabi_node_bounds_mode = tabi_solve.tabi_node_bounds_mode;
    result.tabi_tree_builder = tabi_solve.tabi_tree_builder;
    result.tabi_expansion_center = tabi_solve.tabi_expansion_center;
    result.tabi_fragmentation_safety = tabi_solve.tabi_fragmentation_safety;
    result.tabi_mac_safety = tabi_solve.tabi_mac_safety;
    result.tabi_mac_gap_factor = tabi_solve.tabi_mac_gap_factor;
    result.tabi_tree_depth_max = tabi_solve.tabi_tree_depth_max;
    result.tabi_mac_reject_radius_limited =
        tabi_solve.tabi_mac_reject_radius_limited;
    result.tabi_mac_accept_rate = tabi_solve.tabi_mac_accept_rate;
    result.tabi_node_radius_max = tabi_solve.tabi_node_radius_max;
    result.tabi_node_radius_mean = tabi_solve.tabi_node_radius_mean;
    result.tabi_node_aspect_ratio_max =
        tabi_solve.tabi_node_aspect_ratio_max;
    result.tabi_node_aspect_ratio_mean =
        tabi_solve.tabi_node_aspect_ratio_mean;
    result.tabi_leaf_panel_work_per_matvec =
        tabi_solve.tabi_leaf_panel_work_per_matvec;
    result.tabi_far_node_work_per_matvec =
        tabi_solve.tabi_far_node_work_per_matvec;
    result.tabi_far_geometric_unsafe_count =
        tabi_solve.tabi_far_geometric_unsafe_count;
    result.tabi_far_geometric_unsafe_panels =
        tabi_solve.tabi_far_geometric_unsafe_panels;
    result.tabi_far_geometric_min_gap =
        tabi_solve.tabi_far_geometric_min_gap;
    result.tabi_far_geometric_gap_ratio_min =
        tabi_solve.tabi_far_geometric_gap_ratio_min;
    result.tabi_far_geometric_gap_ratio_mean =
        tabi_solve.tabi_far_geometric_gap_ratio_mean;
    result.tabi_far_geometric_audit_samples =
        tabi_solve.tabi_far_geometric_audit_samples;
    result.tabi_mac_physical_gap_rejects =
        tabi_solve.tabi_mac_physical_gap_rejects;
    result.tabi_mac_convergence_rho_rejects =
        tabi_solve.tabi_mac_convergence_rho_rejects;
    result.tabi_mac_convergence_aspect_rejects =
        tabi_solve.tabi_mac_convergence_aspect_rejects;
    result.tabi_mac_convergence_rho_reject_panels =
        tabi_solve.tabi_mac_convergence_rho_reject_panels;
    result.tabi_mac_convergence_aspect_reject_panels =
        tabi_solve.tabi_mac_convergence_aspect_reject_panels;
    result.tabi_mac_surface_margin_rejects =
        tabi_solve.tabi_mac_surface_margin_rejects;
    result.tabi_mac_surface_margin_reject_panels =
        tabi_solve.tabi_mac_surface_margin_reject_panels;
    result.tabi_mac_surface_margin_applied =
        tabi_solve.tabi_mac_surface_margin_applied;
    result.tabi_mac_surface_margin_skipped =
        tabi_solve.tabi_mac_surface_margin_skipped;
    result.tabi_mac_surface_margin_adaptive_rejects =
        tabi_solve.tabi_mac_surface_margin_adaptive_rejects;
    result.tabi_mac_surface_margin_adaptive_reject_panels =
        tabi_solve.tabi_mac_surface_margin_adaptive_reject_panels;
    result.tabi_mac_surface_margin_skipped_accepted_far =
        tabi_solve.tabi_mac_surface_margin_skipped_accepted_far;
    result.tabi_mac_surface_margin_hard_applied =
        tabi_solve.tabi_mac_surface_margin_hard_applied;
    result.tabi_mac_surface_margin_soft_applied =
        tabi_solve.tabi_mac_surface_margin_soft_applied;
    result.tabi_mac_surface_margin_hard_rejects =
        tabi_solve.tabi_mac_surface_margin_hard_rejects;
    result.tabi_mac_surface_margin_hard_reject_panels =
        tabi_solve.tabi_mac_surface_margin_hard_reject_panels;
    result.tabi_mac_surface_margin_soft_rejects =
        tabi_solve.tabi_mac_surface_margin_soft_rejects;
    result.tabi_mac_surface_margin_soft_reject_panels =
        tabi_solve.tabi_mac_surface_margin_soft_reject_panels;
    result.tabi_fragmented_far_count =
        tabi_solve.tabi_fragmented_far_count;
    result.tabi_fragmented_far_panels =
        tabi_solve.tabi_fragmented_far_panels;
    result.tabi_fragmented_split_rejects =
        tabi_solve.tabi_fragmented_split_rejects;
    result.tabi_fragmentation_center_gap_ratio_max =
        tabi_solve.tabi_fragmentation_center_gap_ratio_max;
    result.tabi_fragmentation_center_gap_ratio_mean =
        tabi_solve.tabi_fragmentation_center_gap_ratio_mean;
    result.tabi_moment_audit = tabi_solve.tabi_moment_audit;
    result.tabi_moment_construction =
        tabi_solve.tabi_moment_construction;
    result.tabi_moment_sum = tabi_solve.tabi_moment_sum;
    result.tabi_convergence_audit = tabi_solve.tabi_convergence_audit;
    result.tabi_convergence_safety = tabi_solve.tabi_convergence_safety;
    result.tabi_convergence_rho_cap = tabi_solve.tabi_convergence_rho_cap;
    result.tabi_mac_aspect_cap = tabi_solve.tabi_mac_aspect_cap;
    result.tabi_convergence_estimator =
        tabi_solve.tabi_convergence_estimator;
    result.tabi_convergence_tail_threshold =
        tabi_solve.tabi_convergence_tail_threshold;
    result.tabi_convergence_estimator_error_cap =
        tabi_solve.tabi_convergence_estimator_error_cap;
    result.tabi_convergence_unstable_rho_cap =
        tabi_solve.tabi_convergence_unstable_rho_cap;
    result.tabi_convergence_calibration =
        tabi_solve.tabi_convergence_calibration;
    result.tabi_convergence_dielectric_amp =
        tabi_solve.tabi_convergence_dielectric_amp;
    result.tabi_convergence_estimator_model =
        tabi_solve.tabi_convergence_estimator_model;
    result.tabi_convergence_channel_calibration =
        tabi_solve.tabi_convergence_channel_calibration;
    result.tabi_convergence_split_policy =
        tabi_solve.tabi_convergence_split_policy;
    result.tabi_convergence_split_selectivity =
        tabi_solve.tabi_convergence_split_selectivity;
    result.tabi_convergence_split_absolute_floor =
        tabi_solve.tabi_convergence_split_absolute_floor;
    result.tabi_convergence_overreject_audit =
        tabi_solve.tabi_convergence_overreject_audit;
    result.tabi_convergence_overreject_detail =
        tabi_solve.tabi_convergence_overreject_detail;
    result.tabi_convergence_shape_audit =
        tabi_solve.tabi_convergence_shape_audit;
    result.tabi_convergence_kernel_response_audit =
        tabi_solve.tabi_convergence_kernel_response_audit;
    result.tabi_channel_missed_semantics =
        tabi_solve.tabi_channel_missed_semantics;
    result.tabi_accepted_miss_diagnostics =
        tabi_solve.tabi_accepted_miss_diagnostics;
    result.tabi_accepted_miss_absolute_floor =
        tabi_solve.tabi_accepted_miss_absolute_floor;
    result.tabi_convergence_k3_accepted_rho_cap =
        tabi_solve.tabi_convergence_k3_accepted_rho_cap;
    result.tabi_mac_distance_policy =
        tabi_solve.tabi_mac_distance_policy;
    result.tabi_mac_surface_margin_factor =
        tabi_solve.tabi_mac_surface_margin_factor;
    result.tabi_mac_surface_margin_policy =
        tabi_solve.tabi_mac_surface_margin_policy;
    result.tabi_mac_surface_margin_aspect_threshold =
        tabi_solve.tabi_mac_surface_margin_aspect_threshold;
    result.tabi_mac_surface_margin_rho_trigger =
        tabi_solve.tabi_mac_surface_margin_rho_trigger;
    result.tabi_mac_surface_margin_soft_aspect_threshold =
        tabi_solve.tabi_mac_surface_margin_soft_aspect_threshold;
    result.tabi_mac_surface_margin_soft_rho_trigger =
        tabi_solve.tabi_mac_surface_margin_soft_rho_trigger;
    result.tabi_mac_surface_margin_soft_factor =
        tabi_solve.tabi_mac_surface_margin_soft_factor;
    result.tabi_convergence_estimator_samples =
        tabi_solve.tabi_convergence_estimator_samples;
    result.tabi_convergence_estimator_unstable_far_count =
        tabi_solve.tabi_convergence_estimator_unstable_far_count;
    result.tabi_convergence_estimator_unstable_far_panels =
        tabi_solve.tabi_convergence_estimator_unstable_far_panels;
    result.tabi_convergence_estimator_split_rejects =
        tabi_solve.tabi_convergence_estimator_split_rejects;
    result.tabi_convergence_estimator_split_reject_panels =
        tabi_solve.tabi_convergence_estimator_split_reject_panels;
    result.tabi_convergence_low_abs_reaccepted_far_count =
        tabi_solve.tabi_convergence_low_abs_reaccepted_far_count;
    result.tabi_convergence_low_abs_reaccepted_far_panels =
        tabi_solve.tabi_convergence_low_abs_reaccepted_far_panels;
    result.tabi_convergence_estimator_tail_ratio_max =
        tabi_solve.tabi_convergence_estimator_tail_ratio_max;
    result.tabi_convergence_estimator_tail_ratio_mean =
        tabi_solve.tabi_convergence_estimator_tail_ratio_mean;
    result.tabi_convergence_estimator_error_max =
        tabi_solve.tabi_convergence_estimator_error_max;
    result.tabi_convergence_estimator_error_mean =
        tabi_solve.tabi_convergence_estimator_error_mean;
    result.tabi_convergence_estimator_abs_error_max =
        tabi_solve.tabi_convergence_estimator_abs_error_max;
    result.tabi_convergence_estimator_abs_error_mean =
        tabi_solve.tabi_convergence_estimator_abs_error_mean;
    result.tabi_convergence_tail_norm_p_max =
        tabi_solve.tabi_convergence_tail_norm_p_max;
    result.tabi_convergence_tail_norm_p_mean =
        tabi_solve.tabi_convergence_tail_norm_p_mean;
    result.tabi_convergence_tail_norm_pm1_max =
        tabi_solve.tabi_convergence_tail_norm_pm1_max;
    result.tabi_convergence_tail_norm_pm1_mean =
        tabi_solve.tabi_convergence_tail_norm_pm1_mean;
    result.tabi_convergence_effective_jcp_amp_max =
        tabi_solve.tabi_convergence_effective_jcp_amp_max;
    result.tabi_convergence_effective_jcp_amp_mean =
        tabi_solve.tabi_convergence_effective_jcp_amp_mean;
    result.tabi_convergence_dielectric_amp_max =
        tabi_solve.tabi_convergence_dielectric_amp_max;
    result.tabi_convergence_dielectric_amp_mean =
        tabi_solve.tabi_convergence_dielectric_amp_mean;
    result.tabi_convergence_channel_estimated_error_max =
        tabi_solve.tabi_convergence_channel_estimated_error_max;
    result.tabi_convergence_channel_estimated_error_mean =
        tabi_solve.tabi_convergence_channel_estimated_error_mean;
    result.tabi_convergence_channel_effective_jcp_amp_max =
        tabi_solve.tabi_convergence_channel_effective_jcp_amp_max;
    result.tabi_convergence_channel_effective_jcp_amp_mean =
        tabi_solve.tabi_convergence_channel_effective_jcp_amp_mean;
    result.tabi_convergence_channel_unstable_far_count =
        tabi_solve.tabi_convergence_channel_unstable_far_count;
    result.tabi_convergence_channel_unstable_far_panels =
        tabi_solve.tabi_convergence_channel_unstable_far_panels;
    result.tabi_convergence_channel_split_rejects =
        tabi_solve.tabi_convergence_channel_split_rejects;
    result.tabi_convergence_channel_split_reject_panels =
        tabi_solve.tabi_convergence_channel_split_reject_panels;
    result.tabi_convergence_channel_low_abs_reaccepted_count =
        tabi_solve.tabi_convergence_channel_low_abs_reaccepted_count;
    result.tabi_convergence_channel_low_abs_reaccepted_panels =
        tabi_solve.tabi_convergence_channel_low_abs_reaccepted_panels;
    result.tabi_convergence_channel_high_error_caught =
        tabi_solve.tabi_convergence_channel_high_error_caught;
    result.tabi_convergence_channel_high_error_missed =
        tabi_solve.tabi_convergence_channel_high_error_missed;
    result.tabi_convergence_channel_high_error_caught_by_split =
        tabi_solve.tabi_convergence_channel_high_error_caught_by_split;
    result.tabi_convergence_estimator_high_error_caught =
        tabi_solve.tabi_convergence_estimator_high_error_caught;
    result.tabi_convergence_estimator_high_error_missed =
        tabi_solve.tabi_convergence_estimator_high_error_missed;
    result.tabi_convergence_estimator_high_error_caught_by_split =
        tabi_solve.tabi_convergence_estimator_high_error_caught_by_split;
    result.tabi_convergence_estimator_phantom_high_error_count =
        tabi_solve.tabi_convergence_estimator_phantom_high_error_count;
    result.tabi_convergence_estimator_phantom_high_error_abs_max =
        tabi_solve.tabi_convergence_estimator_phantom_high_error_abs_max;
    result.tabi_convergence_overreject_count =
        tabi_solve.tabi_convergence_overreject_count;
    result.tabi_convergence_overreject_panels =
        tabi_solve.tabi_convergence_overreject_panels;
    result.tabi_convergence_overreject_caught_count =
        tabi_solve.tabi_convergence_overreject_caught_count;
    result.tabi_convergence_overreject_caught_panels =
        tabi_solve.tabi_convergence_overreject_caught_panels;
    result.tabi_convergence_overreject_abs_error_max =
        tabi_solve.tabi_convergence_overreject_abs_error_max;
    result.tabi_convergence_overreject_abs_error_mean =
        tabi_solve.tabi_convergence_overreject_abs_error_mean;
    result.tabi_convergence_overreject_abs_error_p95 =
        tabi_solve.tabi_convergence_overreject_abs_error_p95;
    result.tabi_convergence_overreject_abs_error_p99 =
        tabi_solve.tabi_convergence_overreject_abs_error_p99;
    result.tabi_convergence_overreject_estimated_abs_error_max =
        tabi_solve.tabi_convergence_overreject_estimated_abs_error_max;
    result.tabi_convergence_overreject_estimated_abs_error_mean =
        tabi_solve.tabi_convergence_overreject_estimated_abs_error_mean;
    result.tabi_convergence_overreject_required_scale_max =
        tabi_solve.tabi_convergence_overreject_required_scale_max;
    result.tabi_convergence_overreject_required_scale_p95 =
        tabi_solve.tabi_convergence_overreject_required_scale_p95;
    result.tabi_convergence_overreject_required_scale_p99 =
        tabi_solve.tabi_convergence_overreject_required_scale_p99;
    result.tabi_convergence_channel_overreject_count =
        tabi_solve.tabi_convergence_channel_overreject_count;
    result.tabi_convergence_channel_overreject_caught_count =
        tabi_solve.tabi_convergence_channel_overreject_caught_count;
    result.tabi_convergence_channel_overreject_panels =
        tabi_solve.tabi_convergence_channel_overreject_panels;
    result.tabi_convergence_channel_overreject_abs_error_max =
        tabi_solve.tabi_convergence_channel_overreject_abs_error_max;
    result.tabi_convergence_channel_overreject_abs_error_mean =
        tabi_solve.tabi_convergence_channel_overreject_abs_error_mean;
    result.tabi_convergence_channel_overreject_required_scale_max =
        tabi_solve.tabi_convergence_channel_overreject_required_scale_max;
    result.tabi_convergence_channel_overreject_required_scale_p99 =
        tabi_solve.tabi_convergence_channel_overreject_required_scale_p99;
    result.tabi_phantom_high_error_channel_count =
        tabi_solve.tabi_phantom_high_error_channel_count;
    result.tabi_phantom_high_error_channel_abs_max =
        tabi_solve.tabi_phantom_high_error_channel_abs_max;
    result.tabi_accepted_miss_channel_count =
        tabi_solve.tabi_accepted_miss_channel_count;
    result.tabi_accepted_miss_channel_panels =
        tabi_solve.tabi_accepted_miss_channel_panels;
    result.tabi_accepted_miss_channel_rho_min =
        tabi_solve.tabi_accepted_miss_channel_rho_min;
    result.tabi_accepted_miss_channel_rho_mean =
        tabi_solve.tabi_accepted_miss_channel_rho_mean;
    result.tabi_accepted_miss_channel_rho_max =
        tabi_solve.tabi_accepted_miss_channel_rho_max;
    result.tabi_accepted_miss_channel_tail_ratio_mean =
        tabi_solve.tabi_accepted_miss_channel_tail_ratio_mean;
    result.tabi_accepted_miss_channel_tail_ratio_max =
        tabi_solve.tabi_accepted_miss_channel_tail_ratio_max;
    result.tabi_accepted_miss_channel_estimated_error_mean =
        tabi_solve.tabi_accepted_miss_channel_estimated_error_mean;
    result.tabi_accepted_miss_channel_estimated_error_max =
        tabi_solve.tabi_accepted_miss_channel_estimated_error_max;
    result.tabi_accepted_miss_worst_channel =
        tabi_solve.tabi_accepted_miss_worst_channel;
    result.tabi_convergence_k3_rho_cap_rejects =
        tabi_solve.tabi_convergence_k3_rho_cap_rejects;
    result.tabi_convergence_k3_rho_cap_reject_panels =
        tabi_solve.tabi_convergence_k3_rho_cap_reject_panels;
    result.tabi_mac_surface_margin_rejects =
        tabi_solve.tabi_mac_surface_margin_rejects;
    result.tabi_mac_surface_margin_reject_panels =
        tabi_solve.tabi_mac_surface_margin_reject_panels;
    result.tabi_mac_surface_margin_applied =
        tabi_solve.tabi_mac_surface_margin_applied;
    result.tabi_mac_surface_margin_skipped =
        tabi_solve.tabi_mac_surface_margin_skipped;
    result.tabi_mac_surface_margin_adaptive_rejects =
        tabi_solve.tabi_mac_surface_margin_adaptive_rejects;
    result.tabi_mac_surface_margin_adaptive_reject_panels =
        tabi_solve.tabi_mac_surface_margin_adaptive_reject_panels;
    result.tabi_mac_surface_margin_skipped_accepted_far =
        tabi_solve.tabi_mac_surface_margin_skipped_accepted_far;
    result.tabi_mac_surface_margin_hard_applied =
        tabi_solve.tabi_mac_surface_margin_hard_applied;
    result.tabi_mac_surface_margin_soft_applied =
        tabi_solve.tabi_mac_surface_margin_soft_applied;
    result.tabi_mac_surface_margin_hard_rejects =
        tabi_solve.tabi_mac_surface_margin_hard_rejects;
    result.tabi_mac_surface_margin_hard_reject_panels =
        tabi_solve.tabi_mac_surface_margin_hard_reject_panels;
    result.tabi_mac_surface_margin_soft_rejects =
        tabi_solve.tabi_mac_surface_margin_soft_rejects;
    result.tabi_mac_surface_margin_soft_reject_panels =
        tabi_solve.tabi_mac_surface_margin_soft_reject_panels;
    result.tabi_moment_audit_nodes = tabi_solve.tabi_moment_audit_nodes;
    result.tabi_moment_translation_rel_error_max =
        tabi_solve.tabi_moment_translation_rel_error_max;
    result.tabi_moment_translation_rel_error_mean =
        tabi_solve.tabi_moment_translation_rel_error_mean;
    result.tabi_moment_translation_offset_max =
        tabi_solve.tabi_moment_translation_offset_max;
    result.tabi_moment_translation_offset_mean =
        tabi_solve.tabi_moment_translation_offset_mean;
    result.tabi_moment_cancellation_ratio_max =
        tabi_solve.tabi_moment_cancellation_ratio_max;
    result.tabi_moment_cancellation_ratio_mean =
        tabi_solve.tabi_moment_cancellation_ratio_mean;
    result.tabi_moment_dynamic_range_max =
        tabi_solve.tabi_moment_dynamic_range_max;
    result.tabi_expansion_audit = tabi_solve.tabi_expansion_audit;
    result.tabi_expansion_audit_channel =
        tabi_solve.tabi_expansion_audit_channel;
    result.tabi_far_expansion_audit_pairs =
        tabi_solve.tabi_far_expansion_audit_pairs;
    result.tabi_far_expansion_audit_panels =
        tabi_solve.tabi_far_expansion_audit_panels;
    result.tabi_far_expansion_direct_ms =
        tabi_solve.tabi_far_expansion_direct_ms;
    result.tabi_far_expansion_rel_error_max =
        tabi_solve.tabi_far_expansion_rel_error_max;
    result.tabi_far_expansion_rel_error_mean =
        tabi_solve.tabi_far_expansion_rel_error_mean;
    result.tabi_far_expansion_rel_error_p95 =
        tabi_solve.tabi_far_expansion_rel_error_p95;
    result.tabi_far_expansion_rel_error_p99 =
        tabi_solve.tabi_far_expansion_rel_error_p99;
    result.tabi_far_expansion_channel_rel_error_max =
        tabi_solve.tabi_far_expansion_channel_rel_error_max;
    result.tabi_far_expansion_channel_rel_error_mean =
        tabi_solve.tabi_far_expansion_channel_rel_error_mean;
    result.tabi_far_expansion_channel_rel_error_p99 =
        tabi_solve.tabi_far_expansion_channel_rel_error_p99;
    result.tabi_far_expansion_worst_channel =
        tabi_solve.tabi_far_expansion_worst_channel;
    result.tabi_far_expansion_mean_gate_pass =
        tabi_solve.tabi_far_expansion_mean_gate_pass;
    result.tabi_far_convergence_rho_max =
        tabi_solve.tabi_far_convergence_rho_max;
    result.tabi_far_convergence_rho_mean =
        tabi_solve.tabi_far_convergence_rho_mean;
    result.tabi_accepted_far_aspect_ratio_max =
        tabi_solve.tabi_accepted_far_aspect_ratio_max;
    result.tabi_accepted_far_aspect_ratio_mean =
        tabi_solve.tabi_accepted_far_aspect_ratio_mean;
    result.tabi_quad_buffer_auto_candidate =
        tabi_solve.tabi_quad_buffer_auto_candidate;
    result.tabi_leaf_tile_status = tabi_solve.tabi_leaf_tile_status;
    result.tabi_leaf_panel_pairs_per_second =
        tabi_solve.tabi_leaf_panel_pairs_per_second;
    result.tabi_near_integration = tabi_solve.tabi_near_integration;
    result.tabi_quad_buffer = tabi_solve.tabi_quad_buffer;
    result.tabi_near_self_pairs = tabi_solve.tabi_near_self_pairs;
    result.tabi_near_one_ring_pairs = tabi_solve.tabi_near_one_ring_pairs;
    result.tabi_near_duffy_pairs = tabi_solve.tabi_near_duffy_pairs;
    result.tabi_near_centroid_pairs = tabi_solve.tabi_near_centroid_pairs;
    result.tabi_near_fallbacks = tabi_solve.tabi_near_fallbacks;
    result.tabi_near_geometric_buffer_pairs =
        tabi_solve.tabi_near_geometric_buffer_pairs;
    result.tabi_near_geometric_buffer_gauss_pairs =
        tabi_solve.tabi_near_geometric_buffer_gauss_pairs;
    result.tabi_near_geometric_buffer_adaptive_pairs =
        tabi_solve.tabi_near_geometric_buffer_adaptive_pairs;
    result.tabi_near_geometric_buffer_centroid_pairs =
        tabi_solve.tabi_near_geometric_buffer_centroid_pairs;
    result.tabi_near_geometric_buffer_min_gap =
        tabi_solve.tabi_near_geometric_buffer_min_gap;
    result.tabi_near_geometric_buffer_mean_gap_ratio =
        tabi_solve.tabi_near_geometric_buffer_mean_gap_ratio;
    result.tabi_near_geometric_buffer_max_edge_mean =
        tabi_solve.tabi_near_geometric_buffer_max_edge_mean;
    result.tabi_near_integration_ms = tabi_solve.tabi_near_integration_ms;
    result.tabi_duffy_subdivisions_max =
        tabi_solve.tabi_duffy_subdivisions_max;
    result.tabi_duffy_subdivisions_mean =
        tabi_solve.tabi_duffy_subdivisions_mean;
    result.tabi_near_kernel_policy = tabi_solve.tabi_near_kernel_policy;
    result.tabi_near_k1_accepted = tabi_solve.tabi_near_k1_accepted;
    result.tabi_near_k1_rejected = tabi_solve.tabi_near_k1_rejected;
    result.tabi_near_k2_accepted = tabi_solve.tabi_near_k2_accepted;
    result.tabi_near_k2_rejected = tabi_solve.tabi_near_k2_rejected;
    result.tabi_near_k3_accepted = tabi_solve.tabi_near_k3_accepted;
    result.tabi_near_k3_rejected = tabi_solve.tabi_near_k3_rejected;
    result.tabi_near_k4_accepted = tabi_solve.tabi_near_k4_accepted;
    result.tabi_near_k4_rejected = tabi_solve.tabi_near_k4_rejected;
    result.tabi_near_one_ring_k1_accepted =
        tabi_solve.tabi_near_one_ring_k1_accepted;
    result.tabi_near_one_ring_k1_rejected =
        tabi_solve.tabi_near_one_ring_k1_rejected;
    result.tabi_near_one_ring_k2_accepted =
        tabi_solve.tabi_near_one_ring_k2_accepted;
    result.tabi_near_one_ring_k2_rejected =
        tabi_solve.tabi_near_one_ring_k2_rejected;
    result.tabi_near_one_ring_k3_accepted =
        tabi_solve.tabi_near_one_ring_k3_accepted;
    result.tabi_near_one_ring_k3_rejected =
        tabi_solve.tabi_near_one_ring_k3_rejected;
    result.tabi_near_one_ring_k4_accepted =
        tabi_solve.tabi_near_one_ring_k4_accepted;
    result.tabi_near_one_ring_k4_rejected =
        tabi_solve.tabi_near_one_ring_k4_rejected;
    result.tabi_near_k1_policy = tabi_solve.tabi_near_k1_policy;
    result.tabi_near_k2_policy = tabi_solve.tabi_near_k2_policy;
    result.tabi_near_k3_policy = tabi_solve.tabi_near_k3_policy;
    result.tabi_near_k4_policy = tabi_solve.tabi_near_k4_policy;
    result.tabi_near_block_audit = tabi_solve.tabi_near_block_audit;
    result.tabi_geometric_near_audit =
        tabi_solve.tabi_geometric_near_audit;
    result.tabi_surface_geometry = tabi_solve.tabi_surface_geometry;
    result.tabi_sphere_mesh_mode = tabi_solve.tabi_sphere_mesh_mode;
    result.tabi_surface_quality_total_area_a2 =
        tabi_solve.tabi_surface_quality_total_area_a2;
    result.tabi_surface_quality_area_rel_error_vs_4piR2 =
        tabi_solve.tabi_surface_quality_area_rel_error_vs_4piR2;
    result.tabi_surface_quality_triangle_area_min =
        tabi_solve.tabi_surface_quality_triangle_area_min;
    result.tabi_surface_quality_triangle_area_mean =
        tabi_solve.tabi_surface_quality_triangle_area_mean;
    result.tabi_surface_quality_triangle_area_max =
        tabi_solve.tabi_surface_quality_triangle_area_max;
    result.tabi_surface_quality_triangle_aspect_max =
        tabi_solve.tabi_surface_quality_triangle_aspect_max;
    result.tabi_surface_quality_normal_radial_dot_min =
        tabi_solve.tabi_surface_quality_normal_radial_dot_min;
    result.tabi_surface_quality_normal_radial_dot_mean =
        tabi_solve.tabi_surface_quality_normal_radial_dot_mean;
    result.tabi_surface_quality_centroid_radial_error_max =
        tabi_solve.tabi_surface_quality_centroid_radial_error_max;
    result.tabi_surface_quality_planar_area_a2 =
        tabi_solve.tabi_surface_quality_planar_area_a2;
    result.tabi_surface_quality_planar_area_rel_error_vs_4piR2 =
        tabi_solve.tabi_surface_quality_planar_area_rel_error_vs_4piR2;
    result.tabi_surface_quality_enclosed_volume_a3 =
        tabi_solve.tabi_surface_quality_enclosed_volume_a3;
    result.tabi_surface_quality_enclosed_volume_rel_error_vs_4piR3_over3 =
        tabi_solve.tabi_surface_quality_enclosed_volume_rel_error_vs_4piR3_over3;
    result.tabi_surface_quality_triangle_area_cv =
        tabi_solve.tabi_surface_quality_triangle_area_cv;
    result.tabi_mesh_vertices = tabi_solve.tabi_mesh_vertices;
    result.tabi_mesh_triangles = tabi_solve.tabi_mesh_triangles;
    result.tabi_boundary_nodes = tabi_solve.tabi_boundary_nodes;
    result.tabi_avg_valence = tabi_solve.tabi_avg_valence;
    result.tabi_nonmanifold_edges = tabi_solve.tabi_nonmanifold_edges;
    result.tabi_solid_angle_fallbacks =
        tabi_solve.tabi_solid_angle_fallbacks;
    result.tabi_diag_min = tabi_solve.tabi_diag_min;
    result.tabi_diag_max = tabi_solve.tabi_diag_max;
    result.tabi_taylor_order = tabi_solve.tabi_taylor_order;
    result.tabi_derivative_mode = tabi_solve.tabi_derivative_mode;
    result.tabi_interaction_build_ms = tabi_solve.tabi_interaction_build_ms;
    result.tabi_tree_eval_ms = tabi_solve.tabi_tree_eval_ms;
    result.tabi_formulation = tabi_solve.tabi_formulation;
    result.tabi_energy_projection = tabi_solve.tabi_energy_projection;
    result.tabi_preconditioner_requested =
        tabi_solve.tabi_preconditioner_requested;
    result.tabi_preconditioner_effective =
        tabi_solve.tabi_preconditioner_effective;
    result.tabi_preconditioner_block_near =
        tabi_solve.tabi_preconditioner_block_near;
    result.tabi_preconditioner_block_panel_cap =
        tabi_solve.tabi_preconditioner_block_panel_cap;
    result.tabi_preconditioner_overlap_panel_cap =
        tabi_solve.tabi_preconditioner_overlap_panel_cap;
    result.tabi_preconditioner_ilu_neighbor_cap =
        tabi_solve.tabi_preconditioner_ilu_neighbor_cap;
    result.tabi_preconditioner_ilu_drop_tol =
        tabi_solve.tabi_preconditioner_ilu_drop_tol;
    result.tabi_preconditioner_memory_budget_bytes =
        tabi_solve.tabi_preconditioner_memory_budget_bytes;
    result.tabi_preconditioner_estimated_bytes =
        tabi_solve.tabi_preconditioner_estimated_bytes;
    result.tabi_preconditioner_valid_blocks =
        tabi_solve.tabi_preconditioner_valid_blocks;
    result.tabi_preconditioner_valid_panels =
        tabi_solve.tabi_preconditioner_valid_panels;
    result.tabi_preconditioner_ill_conditioned_blocks =
        tabi_solve.tabi_preconditioner_ill_conditioned_blocks;
    result.tabi_preconditioner_quality_rejected_blocks =
        tabi_solve.tabi_preconditioner_quality_rejected_blocks;
    result.tabi_preconditioner_quality_probe_max_gain =
        tabi_solve.tabi_preconditioner_quality_probe_max_gain;
    result.tabi_preconditioner_quality_status =
        tabi_solve.tabi_preconditioner_quality_status;
    result.tabi_preconditioner_setup_ms =
        tabi_solve.tabi_preconditioner_setup_ms;
    result.tabi_preconditioner_apply_ms =
        tabi_solve.tabi_preconditioner_apply_ms;
    result.tabi_preconditioner_apply_calls =
        tabi_solve.tabi_preconditioner_apply_calls;
    result.tabi_preconditioner_sparse_nnz_blocks =
        tabi_solve.tabi_preconditioner_sparse_nnz_blocks;
    result.tabi_preconditioner_sparse_avg_degree =
        tabi_solve.tabi_preconditioner_sparse_avg_degree;
    result.tabi_preconditioner_sparse_max_degree =
        tabi_solve.tabi_preconditioner_sparse_max_degree;
    result.tabi_preconditioner_coarse_aggregate_cap =
        tabi_solve.tabi_preconditioner_coarse_aggregate_cap;
    result.tabi_preconditioner_coarse_aggregates =
        tabi_solve.tabi_preconditioner_coarse_aggregates;
    result.tabi_preconditioner_coarse_dimension =
        tabi_solve.tabi_preconditioner_coarse_dimension;
    result.tabi_preconditioner_coarse_base =
        tabi_solve.tabi_preconditioner_coarse_base;
    result.tabi_preconditioner_coarse_matrix_bytes =
        tabi_solve.tabi_preconditioner_coarse_matrix_bytes;
    result.tabi_preconditioner_coarse_pivot_ratio =
        tabi_solve.tabi_preconditioner_coarse_pivot_ratio;
    result.tabi_preconditioner_coarse_status =
        tabi_solve.tabi_preconditioner_coarse_status;
    result.tabi_preconditioner_coarse_setup_matvec_calls =
        tabi_solve.tabi_preconditioner_coarse_setup_matvec_calls;
    result.tabi_preconditioner_coarse_setup_ms =
        tabi_solve.tabi_preconditioner_coarse_setup_ms;
    result.tabi_preconditioner_coarse_apply_ms =
        tabi_solve.tabi_preconditioner_coarse_apply_ms;
    result.tabi_preconditioner_coarse_apply_calls =
        tabi_solve.tabi_preconditioner_coarse_apply_calls;
    result.tabi_preconditioner_coarse_rp_identity_error_max =
        tabi_solve.tabi_preconditioner_coarse_rp_identity_error_max;
    result.tabi_preconditioner_coarse_verify_rel_error_max =
        tabi_solve.tabi_preconditioner_coarse_verify_rel_error_max;
    result.tabi_surface_generator = tabi_solve.tabi_surface_generator;
    result.tabi_discretization = tabi_solve.tabi_discretization;
    result.tabi_rhs_scale = tabi_solve.tabi_rhs_scale;
    result.tabi_dielectric_contrast = tabi_solve.tabi_dielectric_contrast;
    result.tabi_operator_contrast = tabi_solve.tabi_operator_contrast;
    result.tabi_unknown_count = tabi_solve.tabi_unknown_count;
    result.tabi_bie_system = tabi_solve.tabi_bie_system;
    result.tabi_lhs_diag_phi = tabi_solve.tabi_lhs_diag_phi;
    result.tabi_lhs_diag_phin = tabi_solve.tabi_lhs_diag_phin;
    result.tabi_energy_equation = tabi_solve.tabi_energy_equation;
    result.tabi_kernel_convention = tabi_solve.tabi_kernel_convention;
    result.tabi_source_norm_s1 = tabi_solve.tabi_source_norm_s1;
    result.tabi_source_norm_s2 = tabi_solve.tabi_source_norm_s2;
    result.tabi_phi_norm = tabi_solve.tabi_phi_norm;
    result.tabi_phin_norm = tabi_solve.tabi_phin_norm;
    result.tabi_kappa_inv_a = tabi_solve.tabi_kappa_inv_a;
    result.tabi_coulomb_factor_kcal_mol =
        tabi_solve.tabi_coulomb_factor_kcal_mol;
    result.tabi_boundary_energy_kj_mol =
        tabi_solve.tabi_boundary_energy_kj_mol;
    result.tabi_reference_energy_kj_mol =
        tabi_solve.tabi_reference_energy_kj_mol;
    result.tabi_selected_energy_scale =
        tabi_solve.tabi_selected_energy_scale;
    result.tabi_fallback_used = tabi_solve.tabi_fallback_used;
    result.pbe_model_effective = tabi_solve.pbe_model_effective;
    result.linear_seed_solver = tabi_solve.linear_seed_solver;
    result.focus_parent_grid = tabi_solve.focus_parent_grid;
    result.probe_status = "off";
    result.rhs_scale = tabi_solve.rhs_scale;
    result.charge_assigned = tabi_solve.charge_assigned;
    result.charge_moment_error = tabi_solve.charge_moment_error;
    result.boundary_rms = tabi_solve.boundary_rms;
    result.boundary_max = tabi_solve.boundary_max;
    result.linear_residual = tabi_solve.linear_residual;
    result.nl_status = "linear";
    result.focus_enabled = false;
    result.focus_levels = 1;
    result.solver_mode = tabi_solve.solver_mode;
    result.mg_levels = 1;
    result.solvent_iterations = tabi_solve.solvent_total_iterations;
    result.vacuum_iterations = tabi_solve.vacuum_total_iterations;
    result.solvent_iterations_final = tabi_solve.solvent.iterations;
    result.vacuum_iterations_final = tabi_solve.vacuum.iterations;
    result.solvent_rms = tabi_solve.solvent.rms;
    result.vacuum_rms = tabi_solve.vacuum.rms;
    result.solvent_max_dphi = tabi_solve.solvent.max_dphi;
    result.vacuum_max_dphi = tabi_solve.vacuum.max_dphi;
    result.solvent_converged = tabi_solve.solvent.converged;
    result.vacuum_converged = tabi_solve.vacuum.converged;
    result.energy_kcal_mol = tabi_solve.energy_kcal_mol;
    result.energy_kj_mol = tabi_solve.energy_kj_mol;
    result.energy_kt = tabi_solve.energy_kt;
    result.solvent_phase_kcal_mol = tabi_solve.solvent_phase_kcal_mol;
    result.solvent_phase_kj_mol = tabi_solve.solvent_phase_kj_mol;
    result.solvent_phase_kt = tabi_solve.solvent_phase_kt;
    result.atom_corrected_kcal_mol = tabi_solve.atom_corrected_kcal_mol;
    result.atom_corrected_kj_mol = tabi_solve.atom_corrected_kj_mol;
    result.atom_corrected_kt = tabi_solve.atom_corrected_kt;
    result.grid_corrected_kcal_mol = tabi_solve.grid_corrected_kcal_mol;
    result.grid_corrected_kj_mol = tabi_solve.grid_corrected_kj_mol;
    result.grid_corrected_kt = tabi_solve.grid_corrected_kt;
    result.coulomb_energy_kcal_mol = tabi_solve.coulomb_energy_kcal_mol;
    result.coulomb_energy_kj_mol = tabi_solve.coulomb_energy_kj_mol;
    result.coulomb_energy_kt = tabi_solve.coulomb_energy_kt;
    result.reaction_energy_source = tabi_solve.reaction_energy_source;
    result.complex_energy_kcal_mol = tabi_solve.energy_kcal_mol;
    result.complex_energy_kj_mol = tabi_solve.energy_kj_mol;
    result.status = result.tabi_status == "off" || result.tabi_status.empty()
                        ? ((result.solvent_converged &&
                            result.vacuum_converged)
                               ? "ok"
                               : "failed")
                        : result.tabi_status;
    result.pbsa_time_frame_ms =
        elapsed_ms(frame_time_begin, std::chrono::steady_clock::now());
    const ProcessMemorySample frame_memory_after = sample_process_memory();
    result.memory_rss_before_bytes = frame_memory_before.working_set_bytes;
    result.memory_rss_after_bytes = frame_memory_after.working_set_bytes;
    result.memory_rss_delta_bytes =
        frame_memory_after.working_set_bytes >= frame_memory_before.working_set_bytes
            ? frame_memory_after.working_set_bytes -
                  frame_memory_before.working_set_bytes
            : 0;
    result.memory_rss_delta_signed_bytes =
        static_cast<std::int64_t>(frame_memory_after.working_set_bytes) -
        static_cast<std::int64_t>(frame_memory_before.working_set_bytes);
    result.memory_peak_working_set_bytes =
        frame_memory_after.peak_working_set_bytes;
    append_energy_row(result);
    frame_timing.total_ms =
        elapsed_ms(frame_time_begin, std::chrono::steady_clock::now());
    append_timing_row(result, frame_timing);
    append_memory_profile_row(result, frame_timing);
    append_tabi_profile_row(result, frame_timing);
    append_tabi_formulation_audit_row(result);
    append_tabi_near_block_audit_rows(result);
    append_tabi_geometric_near_audit_rows(result);
    append_tabi_far_geometric_audit_row(result);
    append_tabi_moment_audit_row(result);
    append_tabi_surface_quality_row(result);

    ++frames_total_;
    if (result.status == "ok")
      ++frames_converged_;
    else
      ++frames_failed_;
    max_grid_points_seen_ = std::max(max_grid_points_seen_, result.grid_points);
    max_memory_bytes_seen_ =
        std::max(max_memory_bytes_seen_, result.memory_bytes);
    max_tabi_surface_panels_seen_ =
        std::max(max_tabi_surface_panels_seen_, result.tabi_surface_panels);
    max_tabi_tree_nodes_seen_ =
        std::max(max_tabi_tree_nodes_seen_, result.tabi_tree_nodes);
    max_tabi_tree_leaf_count_seen_ = std::max(
        max_tabi_tree_leaf_count_seen_, result.tabi_tree_leaf_count);
    max_tabi_near_interactions_seen_ = std::max(
        max_tabi_near_interactions_seen_, result.tabi_near_interactions);
    max_tabi_far_interactions_seen_ = std::max(
        max_tabi_far_interactions_seen_, result.tabi_far_interactions);
    max_tabi_mesh_vertices_seen_ =
        std::max(max_tabi_mesh_vertices_seen_, result.tabi_mesh_vertices);
    max_tabi_mesh_triangles_seen_ =
        std::max(max_tabi_mesh_triangles_seen_, result.tabi_mesh_triangles);
    max_tabi_boundary_nodes_seen_ =
        std::max(max_tabi_boundary_nodes_seen_, result.tabi_boundary_nodes);
    max_tabi_nonmanifold_edges_seen_ = std::max(
        max_tabi_nonmanifold_edges_seen_, result.tabi_nonmanifold_edges);
    max_tabi_solid_angle_fallbacks_seen_ =
        std::max(max_tabi_solid_angle_fallbacks_seen_,
                 result.tabi_solid_angle_fallbacks);
    min_tabi_diag_seen_ = std::min(min_tabi_diag_seen_, result.tabi_diag_min);
    max_tabi_diag_seen_ = std::max(max_tabi_diag_seen_, result.tabi_diag_max);
    max_solute_points_seen_ =
        std::max(max_solute_points_seen_, result.solute_points);
    max_solvent_points_seen_ =
        std::max(max_solvent_points_seen_, result.solvent_points);
    max_ion_accessible_points_seen_ = std::max(
        max_ion_accessible_points_seen_, result.ion_accessible_points);
    max_stern_excluded_points_seen_ = std::max(
        max_stern_excluded_points_seen_, result.stern_excluded_points);
    max_focus_levels_seen_ =
        std::max(max_focus_levels_seen_, result.focus_levels);
    max_mg_levels_seen_ = std::max(max_mg_levels_seen_, result.mg_levels);
    max_mg_cycles_seen_ = std::max(max_mg_cycles_seen_, result.mg_cycles);
    max_linear_residual_seen_ =
        std::max(max_linear_residual_seen_, result.linear_residual);
    max_boundary_rms_seen_ =
        std::max(max_boundary_rms_seen_, result.boundary_rms);
    max_boundary_max_seen_ =
        std::max(max_boundary_max_seen_, result.boundary_max);

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
        (static_cast<double>(result.solvent_iterations) - iter_mean_) /
        frame_n;
    mg_cycle_mean_ +=
        (static_cast<double>(result.mg_cycles) - mg_cycle_mean_) / frame_n;
    if (log.is_open() && result.status != "ok") {
      TitanUI::log_kv(log, "PBE Warning",
                      "frame " + std::to_string(result.global_index) +
                          " did not meet convergence thresholds");
    }
    return result;
  }

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
  result.rpbe_coulomb_norm = complex_solve.rpbe_coulomb_norm;
  result.rpbe_grad_coulomb_norm = complex_solve.rpbe_grad_coulomb_norm;
  result.rpbe_grad_epsilon_norm = complex_solve.rpbe_grad_epsilon_norm;
  result.rpbe_source_rms = complex_solve.rpbe_source_rms;
  result.rpbe_source_max = complex_solve.rpbe_source_max;
  result.rpbe_source_status = complex_solve.rpbe_source_status;
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
  result.iscm_water_minus_vacuum_kcal_mol =
      complex_solve.iscm_water_minus_vacuum_kcal_mol;
  result.iscm_water_minus_vacuum_kj_mol =
      complex_solve.iscm_water_minus_vacuum_kj_mol;
  result.iscm_water_minus_vacuum_kt =
      complex_solve.iscm_water_minus_vacuum_kt;
  result.nonlinear_stress_kt = complex_solve.nonlinear_stress_kt;
  result.nonlinear_stress_kj_mol = complex_solve.nonlinear_stress_kj_mol;
  result.nonlinear_stress_kcal_mol = complex_solve.nonlinear_stress_kcal_mol;
  result.nonlinear_osmotic_kt = complex_solve.nonlinear_osmotic_kt;
  result.nonlinear_osmotic_kj_mol = complex_solve.nonlinear_osmotic_kj_mol;
  result.nonlinear_osmotic_kcal_mol =
      complex_solve.nonlinear_osmotic_kcal_mol;
  result.nonlinear_energy_correction_kt =
      complex_solve.nonlinear_energy_correction_kt;
  result.nonlinear_energy_correction_kj_mol =
      complex_solve.nonlinear_energy_correction_kj_mol;
  result.nonlinear_energy_correction_kcal_mol =
      complex_solve.nonlinear_energy_correction_kcal_mol;
  result.nonlinear_energy_status = complex_solve.nonlinear_energy_status;
  result.iscm_finalization_mode = complex_solve.iscm_finalization_mode;
  result.iscm_position_source = complex_solve.iscm_position_source;
  result.iscm_position_rms_vs_oracle =
      complex_solve.iscm_position_rms_vs_oracle;
  result.iscm_position_oracle_matches =
      complex_solve.iscm_position_oracle_matches;
  result.iscm_native_independent = complex_solve.iscm_native_independent;
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
  compute_stability_ensemble(atoms_, executor, result);

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
      result.pbsa_complex_reuse_delta_kj_mol =
          complex_solve.pbsa_reuse_delta_kj_mol;
      result.pbsa_receptor_reuse_delta_kj_mol =
          mol1_solve.pbsa_reuse_delta_kj_mol;
      result.pbsa_ligand_reuse_delta_kj_mol =
          mol2_solve.pbsa_reuse_delta_kj_mol;
      result.pbsa_complex_reuse_surface_hits =
          complex_solve.pbsa_reuse_surface_hits;
      result.pbsa_complex_reuse_operator_hits =
          complex_solve.pbsa_reuse_operator_hits;
      result.pbsa_receptor_reuse_surface_hits =
          mol1_solve.pbsa_reuse_surface_hits;
      result.pbsa_receptor_reuse_operator_hits =
          mol1_solve.pbsa_reuse_operator_hits;
      result.pbsa_ligand_reuse_surface_hits =
          mol2_solve.pbsa_reuse_surface_hits;
      result.pbsa_ligand_reuse_operator_hits =
          mol2_solve.pbsa_reuse_operator_hits;
      const std::size_t reuse_hits =
          result.pbsa_complex_reuse_surface_hits +
          result.pbsa_receptor_reuse_surface_hits +
          result.pbsa_ligand_reuse_surface_hits;
      if (args_.pbsa_grid_reuse == "validate")
        result.pbsa_grid_reuse_effective =
            reuse_hits == 3 ? "validate" : "ineligible";
      else if (args_.pbsa_grid_reuse == "auto")
        result.pbsa_grid_reuse_effective = "pending_gate";
      else
        result.pbsa_grid_reuse_effective = "off";
      if (args_.pbsa_perf == "detail") {
        auto append_timing =
            [&](const std::string &endpoint, const std::string &phase,
                const std::string &variant, const TimingBuckets &timing,
                std::size_t surface_hits, std::size_t operator_hits,
                const std::string &status) {
              PbsaEndpointTimingRow row;
              row.trajectory_id = result.trajectory_id;
              row.frame_index = result.frame_index;
              row.global_index = result.global_index;
              row.time_ps = result.time_ps;
              row.endpoint = endpoint;
              row.phase = phase;
              row.variant = variant;
              row.timing = timing;
              row.reuse_surface_hits = surface_hits;
              row.reuse_operator_hits = operator_hits;
              row.status = status;
              pbsa_endpoint_timing_rows_.push_back(std::move(row));
            };
        auto append_endpoint = [&](const std::string &endpoint,
                                   const ReactionSolve &endpoint_solve,
                                   const std::string &status) {
          append_timing(endpoint, "solvent", "baseline",
                        endpoint_solve.pbsa_solvent_timing, 0U, 0U, status);
          const bool reused = endpoint_solve.pbsa_reuse_surface_hits > 0U ||
                              endpoint_solve.pbsa_reuse_operator_hits > 0U;
          append_timing(endpoint, "vacuum", "baseline",
                        endpoint_solve.pbsa_vacuum_timing, 0U, 0U, status);
          if (reused) {
            append_timing(endpoint, "vacuum", "reuse",
                          endpoint_solve.pbsa_reuse_timing,
                          endpoint_solve.pbsa_reuse_surface_hits,
                          endpoint_solve.pbsa_reuse_operator_hits, status);
          }
          TimingBuckets energy_timing;
          energy_timing.energy_ms = endpoint_solve.pbsa_energy_ms;
          append_timing(endpoint, "reaction_energy", "baseline",
                        energy_timing, 0U, 0U, status);
        };
        append_endpoint("complex", complex_solve, result.pbsa_complex_status);
        append_endpoint("receptor", mol1_solve, result.pbsa_receptor_status);
        append_endpoint("ligand", mol2_solve, result.pbsa_ligand_status);
      }
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
    const std::string pbsa_decomp = effective_pbsa_decomp();
    const bool legacy_pbsa_decomp =
        pbsa_decomp == "residue" || pbsa_decomp == "hybrid";
    if (legacy_pbsa_decomp || edecomp_sink_) {
      const auto residue_time_begin = std::chrono::steady_clock::now();
      const std::vector<ResidueFrameRow> residue_rows =
          compute_pbsa_residue_rows(mol1, *mol2, mol1_atoms, mol2_atoms,
                                    atoms_, complex_solve, mol1_solve,
                                    mol2_solve, result, executor,
                                    energy_provider);
      const auto residue_time_end = std::chrono::steady_clock::now();
      result.pbsa_time_residue_ms =
          elapsed_ms(residue_time_begin, residue_time_end);
      if (legacy_pbsa_decomp)
        accumulate_pbsa_residue_rows(result, residue_rows);
      if (edecomp_sink_) {
        std::vector<EDecompRow> sink_rows;
        sink_rows.reserve(residue_rows.size());
        for (const ResidueFrameRow &row : residue_rows) {
          EDecompRow sink_row;
          sink_row.molecule = row.key.molecule;
          sink_row.unit_type = row.key.unit_type;
          sink_row.chain = row.key.chain;
          sink_row.resnr = row.key.resnr;
          sink_row.respos = row.key.respos;
          sink_row.resname = row.key.resname;
          sink_row.atom_index = row.key.atom_index;
          sink_row.atom_name = row.key.atom_name;
          sink_row.elec_kj_mol = row.elec_kj_mol;
          sink_row.vdw_kj_mol = row.vdw_kj_mol;
          sink_row.polar_kj_mol = row.polar_kj_mol;
          sink_row.nonpolar_kj_mol = row.nonpolar_kj_mol;
          sink_rows.push_back(std::move(sink_row));
        }
        edecomp_sink_(sink_rows);
      }
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
  const auto probe_time_begin = std::chrono::steady_clock::now();
  result.probe_status = maybe_write_probe_sites(result, atoms_, complex_solve);
  add_timing_ms(TimingBucket::Probe,
                elapsed_ms(probe_time_begin,
                           std::chrono::steady_clock::now()));

  maybe_write_parity_stage(result, complex_solve);
  maybe_write_stability_audit(result, atoms_, complex_solve);
  const ProcessMemorySample frame_memory_after = sample_process_memory();
  result.memory_rss_before_bytes = frame_memory_before.working_set_bytes;
  result.memory_rss_after_bytes = frame_memory_after.working_set_bytes;
  result.memory_rss_delta_bytes =
      frame_memory_after.working_set_bytes >= frame_memory_before.working_set_bytes
          ? frame_memory_after.working_set_bytes -
                frame_memory_before.working_set_bytes
          : 0;
  result.memory_rss_delta_signed_bytes =
      static_cast<std::int64_t>(frame_memory_after.working_set_bytes) -
      static_cast<std::int64_t>(frame_memory_before.working_set_bytes);
  result.memory_peak_working_set_bytes =
      frame_memory_after.peak_working_set_bytes;
  append_energy_row(result);
  frame_timing.total_ms =
      elapsed_ms(frame_time_begin, std::chrono::steady_clock::now());
  append_timing_row(result, frame_timing);
  append_memory_profile_row(result, frame_timing);

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

void AnalyzerImpl::finalize(std::fstream &log) {
  if (!enabled())
    return;
  if (energy_out_.is_open())
    energy_out_.flush();
  if (timing_out_.is_open())
    timing_out_.flush();
  if (tabi_profile_out_.is_open())
    tabi_profile_out_.flush();
  if (memory_profile_out_.is_open())
    memory_profile_out_.flush();
  if (tabi_formulation_audit_out_.is_open())
    tabi_formulation_audit_out_.flush();
  if (tabi_near_block_audit_out_.is_open())
    tabi_near_block_audit_out_.flush();
  if (tabi_geometric_near_audit_out_.is_open())
    tabi_geometric_near_audit_out_.flush();
  if (tabi_surface_quality_out_.is_open())
    tabi_surface_quality_out_.flush();
  if (stability_audit_out_.is_open())
    stability_audit_out_.flush();
  if (stability_samples_out_.is_open())
    stability_samples_out_.flush();

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
    summary << "pbe_backend " << args_.backend << "\n";
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
    summary << "tabi_sum " << args_.tabi_sum << "\n";
    summary << "tabi_tree_preset " << args_.tabi_tree_preset << "\n";
    summary << "tabi_surface_geometry " << args_.tabi_geometry << "\n";
    summary << "tabi_sphere_mesh " << args_.tabi_sphere_mesh << "\n";
    summary << "tabi_diagonal " << args_.tabi_diagonal << "\n";
    summary << "tabi_taylor_order " << args_.tabi_taylor_order << "\n";
    summary << "tabi_formulation " << args_.tabi_formulation << "\n";
    summary << "tabi_energy_projection " << args_.tabi_energy_projection
            << "\n";
    summary << "tabi_near_integration " << args_.tabi_near_integration << "\n";
    summary << "tabi_quad_buffer " << args_.tabi_quad_buffer << "\n";
    summary << "tabi_interaction_mode " << args_.tabi_interaction_mode << "\n";
    summary << "tabi_precond_block_near " << args_.tabi_precond_block_near
            << "\n";
    summary << "tabi_node_bounds " << args_.tabi_node_bounds << "\n";
    summary << "tabi_tree_builder " << args_.tabi_tree_builder << "\n";
    summary << "tabi_expansion_center " << args_.tabi_expansion_center
            << "\n";
    summary << "tabi_fragmentation_safety "
            << args_.tabi_fragmentation_safety << "\n";
    summary << "tabi_mac_safety " << args_.tabi_mac_safety << "\n";
    summary << "tabi_mac_gap_factor "
            << format_fixed(args_.tabi_mac_gap_factor, 6) << "\n";
    summary << "tabi_moment_audit " << args_.tabi_moment_audit << "\n";
    summary << "tabi_moment_construction "
            << args_.tabi_moment_construction << "\n";
    summary << "tabi_moment_sum " << args_.tabi_moment_sum << "\n";
    summary << "tabi_convergence_audit "
            << args_.tabi_convergence_audit << "\n";
    summary << "tabi_convergence_safety "
            << args_.tabi_convergence_safety << "\n";
    summary << "tabi_convergence_rho_cap "
            << format_fixed(args_.tabi_convergence_rho_cap, 6) << "\n";
    summary << "tabi_mac_aspect_cap "
            << format_fixed(args_.tabi_mac_aspect_cap, 6) << "\n";
    summary << "tabi_convergence_estimator "
            << args_.tabi_convergence_estimator << "\n";
    summary << "tabi_convergence_tail_threshold "
            << format_fixed(args_.tabi_convergence_tail_threshold, 6)
            << "\n";
    summary << "tabi_convergence_estimator_error_cap "
            << format_fixed(args_.tabi_convergence_estimator_error_cap, 8)
            << "\n";
    summary << "tabi_convergence_unstable_rho_cap "
            << format_fixed(args_.tabi_convergence_unstable_rho_cap, 6)
            << "\n";
    summary << "tabi_convergence_calibration "
            << args_.tabi_convergence_calibration << "\n";
    summary << "tabi_expansion_audit " << args_.tabi_expansion_audit
            << "\n";
    summary << "tabi_expansion_audit_sample_count "
            << args_.tabi_expansion_audit_sample_count << "\n";
    summary << "tabi_expansion_audit_max_panels "
            << args_.tabi_expansion_audit_max_panels << "\n";
    summary << "tabi_matvec_audit " << args_.tabi_matvec_audit << "\n";
    summary << "tabi_matvec_audit_sample_count "
            << args_.tabi_matvec_audit_sample_count << "\n";
    summary << "tabi_true_residual_probe "
            << args_.tabi_true_residual_probe << "\n";
    summary << "tabi_true_residual_probe_scope "
            << args_.tabi_true_residual_probe_scope << "\n";
    summary << "tabi_true_residual_probe_interval "
            << args_.tabi_true_residual_probe_interval << "\n";
    summary << "tabi_true_residual_sample_count "
            << args_.tabi_true_residual_sample_count << "\n";
    summary << "tabi_formulation_audit_csv "
            << "pbe_tabi_formulation_audit.csv\n";
    summary << "tabi_profile_csv "
            << (args_.timing == "detail" ? "pbe_tabi_profile.csv" : "off")
            << "\n";
    summary << "tabi_near_block_audit_csv "
            << (args_.timing == "detail" ? "pbe_tabi_near_block_audit.csv"
                                          : "off")
            << "\n";
    summary << "tabi_geometric_near_audit_csv "
            << (args_.timing == "detail"
                    ? "pbe_tabi_geometric_near_audit.csv"
                    : "off")
            << "\n";
    summary << "tabi_far_geometric_audit_csv "
            << (args_.timing == "detail"
                    ? "pbe_tabi_far_geometric_audit.csv"
                    : "off")
            << "\n";
    summary << "tabi_moment_audit_csv "
            << (args_.timing == "detail"
                    ? "pbe_tabi_moment_audit.csv"
                    : "off")
            << "\n";
    summary << "tabi_surface_quality_csv "
            << (args_.timing == "detail" ? "pbe_tabi_surface_quality.csv"
                                          : "off")
            << "\n";
    summary << "memory_profile_csv "
            << (pbe_timing_enabled() ? "pbe_memory_profile.csv" : "off")
            << "\n";
    summary << "tabi_surface_panels_max " << max_tabi_surface_panels_seen_
            << "\n";
    summary << "tabi_mesh_vertices_max " << max_tabi_mesh_vertices_seen_
            << "\n";
    summary << "tabi_mesh_triangles_max " << max_tabi_mesh_triangles_seen_
            << "\n";
    summary << "tabi_boundary_nodes_max " << max_tabi_boundary_nodes_seen_
            << "\n";
    summary << "tabi_diag_min " << format_fixed(min_tabi_diag_seen_, 6)
            << "\n";
    summary << "tabi_diag_max " << format_fixed(max_tabi_diag_seen_, 6)
            << "\n";
    summary << "tabi_nonmanifold_edges_max "
            << max_tabi_nonmanifold_edges_seen_ << "\n";
    summary << "tabi_solid_angle_fallbacks_max "
            << max_tabi_solid_angle_fallbacks_seen_ << "\n";
    summary << "tabi_tree_nodes_max " << max_tabi_tree_nodes_seen_ << "\n";
    summary << "tabi_tree_leaf_count_max " << max_tabi_tree_leaf_count_seen_
            << "\n";
    summary << "tabi_near_interactions_max "
            << max_tabi_near_interactions_seen_ << "\n";
    summary << "tabi_far_interactions_max "
            << max_tabi_far_interactions_seen_ << "\n";
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
      summary << "test_boundary_convention "
              << args_.test_boundary_convention << "\n";
      summary << "test_field_units " << args_.test_field_units << "\n";
      summary << "test_operator_convention "
              << args_.test_operator_convention << "\n";
      summary << "test_rhs_map " << args_.test_rhs_map << "\n";
      summary << "test_screen_map " << args_.test_screen_map << "\n";
      summary << "test_denom_map " << args_.test_denom_map << "\n";
      if (args_.test_midpoint_diff != "off") {
        summary << "test_midpoint_diff " << args_.test_midpoint_diff << "\n";
        summary << "midpoint_face_diff_csv pbe_midpoint_face_diff.csv\n";
      }
      if (args_.test_charge_grid_diff != "off") {
        summary << "test_charge_grid_diff " << args_.test_charge_grid_diff
                << "\n";
        summary << "charge_grid_diff_csv pbe_charge_grid_diff.csv\n";
        summary << "charge_grid_diff_summary_csv "
                   "pbe_charge_grid_diff_summary.csv\n";
      }
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
      summary << "test_charge_grid_diff " << args_.test_charge_grid_diff
              << "\n";
      summary << "test_surface_map " << args_.test_surface_map << "\n";
      summary << "test_boundary_map " << args_.test_boundary_map << "\n";
      summary << "test_boundary_convention "
              << args_.test_boundary_convention << "\n";
      summary << "test_field_units " << args_.test_field_units << "\n";
      summary << "test_operator_convention "
              << args_.test_operator_convention << "\n";
      summary << "test_rhs_map " << args_.test_rhs_map << "\n";
      summary << "test_screen_map " << args_.test_screen_map << "\n";
      summary << "test_denom_map " << args_.test_denom_map << "\n";
      summary << "field_samples_csv pbe_field_samples.csv\n";
      summary << "boundary_samples_csv pbe_boundary_samples.csv\n";
      summary << "operator_stage_csv pbe_operator_stage.csv\n";
      if (args_.parity_dump == "all")
        summary << "sor_trace_csv pbe_sor_trace.csv\n";
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
    TitanUI::log_kv(log, "PBE Backend", args_.backend);
    if (args_.backend == "tabi") {
      TitanUI::log_kv(log, "TABI Geometry", args_.tabi_geometry);
      TitanUI::log_kv(log, "TABI Diagonal", args_.tabi_diagonal);
      TitanUI::log_kv(log, "TABI Sum", args_.tabi_sum);
      TitanUI::log_kv(log, "TABI Tree Preset", args_.tabi_tree_preset);
      TitanUI::log_kv(log, "TABI Formulation", args_.tabi_formulation);
      TitanUI::log_kv(log, "TABI Energy Projection",
                      args_.tabi_energy_projection);
      TitanUI::log_kv(log, "TABI Near Integration",
                      args_.tabi_near_integration);
      TitanUI::log_kv(log, "TABI Quad Buffer", args_.tabi_quad_buffer);
      TitanUI::log_kv(log, "TABI Interaction Mode",
                      args_.tabi_interaction_mode);
      TitanUI::log_kv(log, "TABI Preconditioner Block Near",
                      args_.tabi_precond_block_near);
      TitanUI::log_kv(log, "TABI Node Bounds", args_.tabi_node_bounds);
      TitanUI::log_kv(log, "TABI Tree Builder", args_.tabi_tree_builder);
      TitanUI::log_kv(log, "TABI Expansion Center",
                      args_.tabi_expansion_center);
      TitanUI::log_kv(log, "TABI Fragmentation Safety",
                      args_.tabi_fragmentation_safety);
      TitanUI::log_kv(log, "TABI MAC Safety", args_.tabi_mac_safety);
      TitanUI::log_kv(log, "TABI MAC Gap Factor",
                      format_fixed(args_.tabi_mac_gap_factor, 6));
      TitanUI::log_kv(log, "TABI Moment Audit", args_.tabi_moment_audit);
      TitanUI::log_kv(log, "TABI Moment Construction",
                      args_.tabi_moment_construction);
      TitanUI::log_kv(log, "TABI Moment Sum", args_.tabi_moment_sum);
      TitanUI::log_kv(log, "TABI Convergence Audit",
                      args_.tabi_convergence_audit);
      TitanUI::log_kv(log, "TABI Convergence Safety",
                      args_.tabi_convergence_safety);
      TitanUI::log_kv(log, "TABI Convergence Rho Cap",
                      format_fixed(args_.tabi_convergence_rho_cap, 6));
      TitanUI::log_kv(log, "TABI MAC Aspect Cap",
                      format_fixed(args_.tabi_mac_aspect_cap, 6));
      TitanUI::log_kv(log, "TABI Convergence Estimator",
                      args_.tabi_convergence_estimator);
      TitanUI::log_kv(
          log, "TABI Convergence Tail Threshold",
          format_fixed(args_.tabi_convergence_tail_threshold, 6));
      TitanUI::log_kv(
          log, "TABI Convergence Estimator Error Cap",
          format_fixed(args_.tabi_convergence_estimator_error_cap, 8));
      TitanUI::log_kv(
          log, "TABI Convergence Unstable Rho Cap",
          format_fixed(args_.tabi_convergence_unstable_rho_cap, 6));
      TitanUI::log_kv(log, "TABI Convergence Calibration",
                      args_.tabi_convergence_calibration);
      TitanUI::log_kv(log, "TABI Expansion Audit",
                      args_.tabi_expansion_audit);
      TitanUI::log_kv(log, "TABI True Residual Probe",
                      args_.tabi_true_residual_probe);
      TitanUI::log_kv(log, "TABI True Residual Probe Scope",
                      args_.tabi_true_residual_probe_scope);
      TitanUI::log_kv(log, "TABI True Residual Probe Interval",
                      std::to_string(args_.tabi_true_residual_probe_interval));
      TitanUI::log_kv(log, "TABI True Residual Sample Count",
                      std::to_string(args_.tabi_true_residual_sample_count));
      TitanUI::log_kv(log, "TABI Formulation Audit",
                      "pbe_tabi_formulation_audit.csv");
      TitanUI::log_kv(log, "TABI Profile CSV",
                      args_.timing == "detail" ? "pbe_tabi_profile.csv"
                                                : "off");
      TitanUI::log_kv(log, "TABI Near Block Audit",
                      args_.timing == "detail"
                          ? "pbe_tabi_near_block_audit.csv"
                          : "off");
      TitanUI::log_kv(log, "TABI Geometric Near Audit",
                      args_.timing == "detail"
                          ? "pbe_tabi_geometric_near_audit.csv"
                          : "off");
      TitanUI::log_kv(log, "TABI Far Geometric Audit",
                      args_.timing == "detail"
                          ? "pbe_tabi_far_geometric_audit.csv"
                          : "off");
      TitanUI::log_kv(log, "TABI Moment Audit",
                      args_.timing == "detail"
                          ? "pbe_tabi_moment_audit.csv"
                          : "off");
      TitanUI::log_kv(log, "TABI Surface Quality",
                      args_.timing == "detail"
                          ? "pbe_tabi_surface_quality.csv"
                          : "off");
    }
    TitanUI::log_kv(log, "PBE Memory Profile",
                    pbe_timing_enabled() ? "pbe_memory_profile.csv" : "off");
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
      TitanUI::log_kv(log, "Charge Grid Diff", args_.test_charge_grid_diff);
      TitanUI::log_kv(log, "Surface Map Replay", args_.test_surface_map);
      TitanUI::log_kv(log, "Midpoint Epsilon Replay",
                      args_.test_midpoint_epsilon);
      TitanUI::log_kv(log, "Midpoint Face Diff", args_.test_midpoint_diff);
      TitanUI::log_kv(log, "Boundary Map Replay", args_.test_boundary_map);
      TitanUI::log_kv(log, "Boundary Convention",
                      args_.test_boundary_convention);
      TitanUI::log_kv(log, "Field Units", args_.test_field_units);
      TitanUI::log_kv(log, "Operator Convention",
                      args_.test_operator_convention);
      TitanUI::log_kv(log, "RHS Replay", args_.test_rhs_map);
      TitanUI::log_kv(log, "Screen Replay", args_.test_screen_map);
      TitanUI::log_kv(log, "Denom Replay", args_.test_denom_map);
      if (args_.test_midpoint_diff != "off")
        TitanUI::log_kv(log, "Midpoint Face Diff CSV",
                        "pbe_midpoint_face_diff.csv");
      if (args_.test_charge_grid_diff != "off") {
        TitanUI::log_kv(log, "Charge Grid Diff CSV",
                        "pbe_charge_grid_diff.csv");
        TitanUI::log_kv(log, "Charge Grid Diff Summary CSV",
                        "pbe_charge_grid_diff_summary.csv");
      }
      if (args_.parity_map_dump != "off") {
        TitanUI::log_kv(log, "Field Samples CSV", "pbe_field_samples.csv");
        TitanUI::log_kv(log, "Boundary Samples CSV",
                        "pbe_boundary_samples.csv");
        TitanUI::log_kv(log, "Operator Stage CSV",
                        "pbe_operator_stage.csv");
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
