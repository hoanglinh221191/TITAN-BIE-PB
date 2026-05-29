#include "titan_dihedral.h"

#include "BoxMatrix.h"
#include "info.h"
#include "titan_executor.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <limits>
#include <optional>
#include <sstream>
#include <stdexcept>

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kDeg = 180.0 / kPi;
constexpr double kRad = kPi / 180.0;
constexpr double kEps = 1.0e-12;

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

std::string normalized_atom_name(std::string value) {
  return uppercase_copy(trim_copy(std::move(value)));
}

std::string join_output_path(const std::string &base_dir,
                             const std::string &filename) {
  if (base_dir.empty() || base_dir == ".")
    return filename;
  return (std::filesystem::path(base_dir) / filename).string();
}

double dot3(const TitanDihedralAnalyzer::Point &a,
            const TitanDihedralAnalyzer::Point &b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

TitanDihedralAnalyzer::Point
cross3(const TitanDihedralAnalyzer::Point &a,
       const TitanDihedralAnalyzer::Point &b) {
  return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z,
          a.x * b.y - a.y * b.x};
}

double norm3(const TitanDihedralAnalyzer::Point &a) {
  return std::sqrt(dot3(a, a));
}

TitanDihedralAnalyzer::Point
sub3(const TitanDihedralAnalyzer::Point &a,
     const TitanDihedralAnalyzer::Point &b) {
  return {a.x - b.x, a.y - b.y, a.z - b.z};
}

TitanDihedralAnalyzer::Point
add3(const TitanDihedralAnalyzer::Point &a,
     const TitanDihedralAnalyzer::Point &b) {
  return {a.x + b.x, a.y + b.y, a.z + b.z};
}

TitanDihedralAnalyzer::Point scale3(const TitanDihedralAnalyzer::Point &a,
                                    double scale) {
  return {a.x * scale, a.y * scale, a.z * scale};
}

} // namespace

double TitanDihedralAnalyzer::compute_dihedral_deg(const Point &p0,
                                                   const Point &p1,
                                                   const Point &p2,
                                                   const Point &p3) {
  const Point b0 = scale3(sub3(p0, p1), -1.0);
  const Point b1 = sub3(p2, p1);
  const Point b2 = sub3(p3, p2);
  const double b1_norm = norm3(b1);
  if (b1_norm <= kEps)
    return std::numeric_limits<double>::quiet_NaN();

  const Point b1_unit = scale3(b1, 1.0 / b1_norm);
  const Point v = sub3(b0, scale3(b1_unit, dot3(b0, b1_unit)));
  const Point w = sub3(b2, scale3(b1_unit, dot3(b2, b1_unit)));
  if (norm3(v) <= kEps || norm3(w) <= kEps)
    return std::numeric_limits<double>::quiet_NaN();

  const double x = dot3(v, w);
  const double y = dot3(cross3(b1_unit, v), w);
  double angle = std::atan2(y, x) * kDeg;
  while (angle > 180.0)
    angle -= 360.0;
  while (angle <= -180.0)
    angle += 360.0;
  return angle;
}

void TitanDihedralAnalyzer::configure(const DihedralArgs &args) {
  requested_ = args.flag_dihedral;
  kind_ = args.kind;
  residue_selection_ = args.residues;
  stride_ = args.stride;
  export_tica_ = args.export_tica;
  if (kind_ != "backbone")
    throw std::invalid_argument("Dihedral: v1 supports backbone only.");
  if (stride_ < 1)
    throw std::invalid_argument("Dihedral: -dihedral-stride must be >= 1.");
  enabled_ = requested_ && initialized_;
}

void TitanDihedralAnalyzer::init(const Multi_MOL &mol, const arguments &,
                                 std::fstream &log) {
  if (!requested_)
    return;
  if (!initialized_) {
    build_residue_map(mol.info);
    parse_residue_selection();
    for (std::size_t i = 0; i < residues_.size(); ++i) {
      if (residue_selected(residues_[i]))
        selected_indices_.push_back(i);
    }
    rama_valid_samples_.assign(residues_.size(), 0);
    initialized_ = !selected_indices_.empty();
    enabled_ = requested_ && initialized_;
    if (log.is_open()) {
      log << "Dihedral kind: " << kind_ << '\n';
      log << "Dihedral selected residues: " << selected_indices_.size()
          << '\n';
      log << "Dihedral PBC policy: " << pbc_policy_ << '\n';
    }
  }
  reset_trajectory_state();
}

void TitanDihedralAnalyzer::build_residue_map(
    const std::vector<mol_info> &mol_info_vec) {
  residues_.clear();
  if (mol_info_vec.empty())
    return;

  std::size_t prev_chain = std::numeric_limits<std::size_t>::max();
  std::size_t prev_resnr = std::numeric_limits<std::size_t>::max();
  std::size_t local_index = 0;

  for (std::size_t atom_idx = 0; atom_idx < mol_info_vec.size(); ++atom_idx) {
    const auto &atom_info = mol_info_vec[atom_idx];
    const std::size_t chainnum = atom_info.get_chainnum();
    const std::size_t resnr = atom_info.get_resnr();
    const bool new_residue =
        residues_.empty() || chainnum != prev_chain || resnr != prev_resnr;
    if (new_residue) {
      if (residues_.empty() || chainnum != prev_chain)
        local_index = 0;
      else
        ++local_index;
      residues_.push_back(ResidueBackbone{
          residues_.size() + 1, resnr, chainnum, local_index,
          atom_info.get_chain(), atom_info.get_res_name(), -1, -1, -1});
      prev_chain = chainnum;
      prev_resnr = resnr;
    }

    ResidueBackbone &residue = residues_.back();
    const std::string atom_name = normalized_atom_name(atom_info.get_atom());
    if (atom_name == "N" && residue.atom_n < 0) {
      residue.atom_n = static_cast<int>(atom_idx);
    } else if (atom_name == "CA" && residue.atom_ca < 0) {
      residue.atom_ca = static_cast<int>(atom_idx);
    } else if (atom_name == "C" && residue.atom_c < 0) {
      residue.atom_c = static_cast<int>(atom_idx);
    }
  }
}

void TitanDihedralAnalyzer::parse_residue_selection() {
  ranges_.clear();
  const std::string raw = trim_copy(residue_selection_);
  if (raw.empty() || raw == "all")
    return;

  std::stringstream ss(raw);
  std::string token;
  while (std::getline(ss, token, ',')) {
    token = trim_copy(token);
    if (token.empty())
      continue;

    ResidueRange range;
    const auto colon = token.find(':');
    std::string span = token;
    if (colon != std::string::npos) {
      range.has_chain = true;
      range.chain = token.substr(0, colon);
      span = token.substr(colon + 1);
      if (range.chain.empty() || span.empty())
        throw std::invalid_argument(
            "Dihedral: invalid -dihedral-residues selector.");
    }

    const auto dash = span.find('-');
    try {
      if (dash == std::string::npos) {
        range.begin = static_cast<std::size_t>(std::stoul(span));
        range.end = range.begin;
      } else {
        range.begin = static_cast<std::size_t>(std::stoul(span.substr(0, dash)));
        range.end = static_cast<std::size_t>(std::stoul(span.substr(dash + 1)));
      }
    } catch (...) {
      throw std::invalid_argument(
          "Dihedral: invalid -dihedral-residues selector.");
    }
    if (range.begin == 0 || range.end == 0 || range.begin > range.end)
      throw std::invalid_argument(
          "Dihedral: invalid -dihedral-residues selector.");
    ranges_.push_back(std::move(range));
  }
}

bool TitanDihedralAnalyzer::residue_selected(
    const ResidueBackbone &residue) const {
  if (ranges_.empty())
    return true;
  for (const ResidueRange &range : ranges_) {
    if (range.has_chain && range.chain != residue.chain)
      continue;
    if (residue.resnr >= range.begin && residue.resnr <= range.end)
      return true;
  }
  return false;
}

void TitanDihedralAnalyzer::reset_trajectory_state() {
  frames_seen_ = 0;
  frames_used_ = 0;
  finite_counts_ = {};
  missing_counts_ = {};
  frames_.clear();
  rama_counts_.clear();
  rama_valid_samples_.assign(residues_.size(), 0);
}

bool TitanDihedralAnalyzer::adjacent_residues(std::size_t lhs,
                                              std::size_t rhs) const noexcept {
  if (lhs >= residues_.size() || rhs >= residues_.size())
    return false;
  const auto &a = residues_[lhs];
  const auto &b = residues_[rhs];
  return a.chainnum == b.chainnum && b.local_index == a.local_index + 1;
}

bool TitanDihedralAnalyzer::has_valid_box(const box &frame_box) const noexcept {
  return frame_box.get_a() > 0.0 && frame_box.get_b() > 0.0 &&
         frame_box.get_c() > 0.0;
}

TitanDihedralAnalyzer::Point
TitanDihedralAnalyzer::atom_point(const MOL &frame, int atom_idx) const {
  if (atom_idx < 0)
    return {std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN()};
  const std::size_t ai = static_cast<std::size_t>(atom_idx);
  if (ai >= frame.pos.x.size() || ai >= frame.pos.y.size() ||
      ai >= frame.pos.z.size()) {
    return {std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN(),
            std::numeric_limits<double>::quiet_NaN()};
  }
  return {frame.pos.x[ai], frame.pos.y[ai], frame.pos.z[ai]};
}

double TitanDihedralAnalyzer::torsion(const MOL &frame, int atom0, int atom1,
                                      int atom2, int atom3, bool use_mic,
                                      const BoxMatrix *box_matrix) const {
  if (atom0 < 0 || atom1 < 0 || atom2 < 0 || atom3 < 0)
    return std::numeric_limits<double>::quiet_NaN();

  const Point raw0 = atom_point(frame, atom0);
  const Point raw1 = atom_point(frame, atom1);
  const Point raw2 = atom_point(frame, atom2);
  const Point raw3 = atom_point(frame, atom3);
  if (!std::isfinite(raw0.x) || !std::isfinite(raw1.x) ||
      !std::isfinite(raw2.x) || !std::isfinite(raw3.x))
    return std::numeric_limits<double>::quiet_NaN();

  if (!use_mic || box_matrix == nullptr)
    return compute_dihedral_deg(raw0, raw1, raw2, raw3);

  Point v10 = sub3(raw0, raw1);
  Point v12 = sub3(raw2, raw1);
  Point v23 = sub3(raw3, raw2);
  box_matrix->apply_mic(v10.x, v10.y, v10.z);
  box_matrix->apply_mic(v12.x, v12.y, v12.z);
  box_matrix->apply_mic(v23.x, v23.y, v23.z);
  const Point p1 = raw1;
  const Point p0 = add3(p1, v10);
  const Point p2 = add3(p1, v12);
  const Point p3 = add3(p2, v23);
  return compute_dihedral_deg(p0, p1, p2, p3);
}

void TitanDihedralAnalyzer::process_frame(const Multi_MOL &mol,
                                          std::size_t snap_idx,
                                          TitanExecutor &executor) {
  if (!enabled_)
    return;
  const std::size_t trajectory_frame = static_cast<std::size_t>(frames_seen_);
  ++frames_seen_;
  if (stride_ > 1 &&
      (trajectory_frame % static_cast<std::size_t>(stride_)) != 0)
    return;
  if (snap_idx >= mol.m_mol.size())
    return;

  const MOL &frame = mol.m_mol[snap_idx];
  const bool use_mic = has_valid_box(frame.Box);
  std::optional<BoxMatrix> box_matrix;
  if (use_mic)
    box_matrix.emplace(frame.Box);
  ++frames_used_;

  FrameRecord record;
  record.frame = trajectory_frame;
  record.time_ps = frame.get_time();
  record.samples.resize(selected_indices_.size());

  const std::size_t worker_count =
      std::max<std::size_t>(1, executor.capacity());
  std::vector<std::array<std::uint64_t, 3>> worker_finite(worker_count);
  std::vector<std::array<std::uint64_t, 3>> worker_missing(worker_count);
  std::vector<std::map<std::tuple<std::size_t, int, int>, std::uint64_t>>
      worker_rama_counts(worker_count);
  std::vector<std::map<std::size_t, std::uint64_t>> worker_rama_valid(
      worker_count);

  executor.parallel_for(selected_indices_.size(), worker_count, 16,
                        [&](std::size_t worker_id, std::size_t begin,
                            std::size_t end) {
    for (std::size_t selected_pos = begin; selected_pos < end; ++selected_pos) {
      const std::size_t residue_index = selected_indices_[selected_pos];
      const ResidueBackbone &residue = residues_[residue_index];
      const ResidueBackbone *prev = nullptr;
      const ResidueBackbone *next = nullptr;
      if (residue_index > 0 &&
          adjacent_residues(residue_index - 1, residue_index))
        prev = &residues_[residue_index - 1];
      if (residue_index + 1 < residues_.size() &&
          adjacent_residues(residue_index, residue_index + 1))
        next = &residues_[residue_index + 1];

      double phi = std::numeric_limits<double>::quiet_NaN();
      double psi = std::numeric_limits<double>::quiet_NaN();
      double omega = std::numeric_limits<double>::quiet_NaN();

      if (prev != nullptr) {
        phi = torsion(frame, prev->atom_c, residue.atom_n, residue.atom_ca,
                      residue.atom_c, use_mic,
                      box_matrix ? &(*box_matrix) : nullptr);
      }
      if (next != nullptr) {
        psi = torsion(frame, residue.atom_n, residue.atom_ca, residue.atom_c,
                      next->atom_n, use_mic,
                      box_matrix ? &(*box_matrix) : nullptr);
        omega = torsion(frame, residue.atom_ca, residue.atom_c, next->atom_n,
                        next->atom_ca, use_mic,
                        box_matrix ? &(*box_matrix) : nullptr);
      }

      const std::array<double, 3> angles{phi, psi, omega};
      for (std::size_t i = 0; i < angles.size(); ++i) {
        if (std::isfinite(angles[i]))
          ++worker_finite[worker_id][i];
        else
          ++worker_missing[worker_id][i];
      }
      if (std::isfinite(phi) && std::isfinite(psi)) {
        ++worker_rama_valid[worker_id][residue_index];
        ++worker_rama_counts[worker_id][std::make_tuple(
            residue_index, ramachandran_bin(phi), ramachandran_bin(psi))];
      }
      record.samples[selected_pos] =
          AngleSample{residue_index, phi, psi, omega};
    }
  });

  for (std::size_t worker_id = 0; worker_id < worker_count; ++worker_id) {
    for (std::size_t i = 0; i < finite_counts_.size(); ++i) {
      finite_counts_[i] += worker_finite[worker_id][i];
      missing_counts_[i] += worker_missing[worker_id][i];
    }
    for (const auto &[residue_index, count] : worker_rama_valid[worker_id])
      rama_valid_samples_[residue_index] += count;
    for (const auto &[key, count] : worker_rama_counts[worker_id])
      rama_counts_[key] += count;
  }

  frames_.push_back(std::move(record));
}

std::string
TitanDihedralAnalyzer::residue_label(const ResidueBackbone &residue) const {
  std::ostringstream out;
  out << residue.chain << ':' << residue.resnr << ':' << residue.resname;
  return out.str();
}

std::string TitanDihedralAnalyzer::format_angle(double value) const {
  if (!std::isfinite(value))
    return "nan";
  std::ostringstream out;
  out << std::fixed << std::setprecision(6) << value;
  return out.str();
}

int TitanDihedralAnalyzer::ramachandran_bin(double angle) const {
  if (angle >= 180.0)
    return 180;
  return static_cast<int>(std::floor(angle / 10.0) * 10.0);
}

void TitanDihedralAnalyzer::export_timeseries(const std::string &filename,
                                              std::fstream &log) const {
  std::ofstream out(filename);
  out << "frame,time_ps,residue,chain,resnr,resname,phi_deg,psi_deg,omega_deg\n";
  out << std::fixed << std::setprecision(6);
  for (const FrameRecord &frame : frames_) {
    for (const AngleSample &sample : frame.samples) {
      const ResidueBackbone &residue = residues_[sample.residue_index];
      out << frame.frame << ',' << frame.time_ps << ','
          << residue_label(residue) << ',' << residue.chain << ','
          << residue.resnr << ',' << residue.resname << ','
          << format_angle(sample.phi) << ',' << format_angle(sample.psi) << ','
          << format_angle(sample.omega) << '\n';
    }
  }
  if (log.is_open())
    log << "Dihedral timeseries: " << filename << '\n';
}

void TitanDihedralAnalyzer::export_ramachandran(const std::string &filename,
                                                std::fstream &log) const {
  std::ofstream out(filename);
  out << "residue,chain,resnr,resname,phi_bin_deg,psi_bin_deg,count,occupancy,valid_samples\n";
  out << std::fixed << std::setprecision(6);
  for (const auto &entry : rama_counts_) {
    const std::size_t residue_index = std::get<0>(entry.first);
    const int phi_bin = std::get<1>(entry.first);
    const int psi_bin = std::get<2>(entry.first);
    const std::uint64_t count = entry.second;
    const std::uint64_t valid = residue_index < rama_valid_samples_.size()
                                    ? rama_valid_samples_[residue_index]
                                    : 0;
    const double occupancy =
        valid > 0 ? static_cast<double>(count) / static_cast<double>(valid)
                  : 0.0;
    const ResidueBackbone &residue = residues_[residue_index];
    out << residue_label(residue) << ',' << residue.chain << ','
        << residue.resnr << ',' << residue.resname << ',' << phi_bin << ','
        << psi_bin << ',' << count << ',' << occupancy << ',' << valid
        << '\n';
  }
  if (log.is_open())
    log << "Dihedral Ramachandran: " << filename << '\n';
}

void TitanDihedralAnalyzer::export_summary(const std::string &filename,
                                           std::fstream &log) const {
  std::ofstream out(filename);
  out << "kind: " << kind_ << '\n';
  out << "residue_selection: " << residue_selection_ << '\n';
  out << "stride: " << stride_ << '\n';
  out << "export_tica: " << (export_tica_ ? "true" : "false") << '\n';
  out << "pbc_policy: " << pbc_policy_ << '\n';
  out << "frames_seen: " << frames_seen_ << '\n';
  out << "frames_used: " << frames_used_ << '\n';
  out << "residues_total: " << residues_.size() << '\n';
  out << "residues_selected: " << selected_indices_.size() << '\n';
  out << "finite_phi: " << finite_counts_[0] << '\n';
  out << "finite_psi: " << finite_counts_[1] << '\n';
  out << "finite_omega: " << finite_counts_[2] << '\n';
  out << "missing_phi: " << missing_counts_[0] << '\n';
  out << "missing_psi: " << missing_counts_[1] << '\n';
  out << "missing_omega: " << missing_counts_[2] << '\n';
  if (log.is_open())
    log << "Dihedral summary: " << filename << '\n';
}

void TitanDihedralAnalyzer::export_feature_matrix(const std::string &filename,
                                                  std::fstream &log) const {
  std::ofstream out(filename);
  out << "frame,time_ps";
  for (const std::size_t residue_index : selected_indices_) {
    const std::string label = residue_label(residues_[residue_index]);
    out << ',' << label << "_phi_sin," << label << "_phi_cos";
    out << ',' << label << "_psi_sin," << label << "_psi_cos";
    out << ',' << label << "_omega_sin," << label << "_omega_cos";
  }
  out << '\n';
  out << std::fixed << std::setprecision(6);

  for (const FrameRecord &frame : frames_) {
    out << frame.frame << ',' << frame.time_ps;
    for (const AngleSample &sample : frame.samples) {
      const std::array<double, 3> angles{sample.phi, sample.psi, sample.omega};
      for (double angle : angles) {
        if (!std::isfinite(angle)) {
          out << ",nan,nan";
        } else {
          out << ',' << std::sin(angle * kRad) << ','
              << std::cos(angle * kRad);
        }
      }
    }
    out << '\n';
  }
  if (log.is_open())
    log << "Dihedral feature matrix: " << filename << '\n';
}

void TitanDihedralAnalyzer::finalize_trajectory(
    const std::string &base_dir, const info &info_output, const arguments &args,
    const std::vector<std::string> &input_line, int input_name,
    std::fstream &log) {
  if (!enabled_)
    return;

  export_timeseries(
      join_output_path(base_dir, info_output.get_dihedral_timeseries_filename(
                                     args, input_line, input_name)),
      log);
  export_ramachandran(
      join_output_path(base_dir, info_output.get_dihedral_ramachandran_filename(
                                     args, input_line, input_name)),
      log);
  export_summary(
      join_output_path(base_dir, info_output.get_dihedral_summary_filename(
                                     args, input_line, input_name)),
      log);
  if (export_tica_) {
    export_feature_matrix(
        join_output_path(base_dir,
                         info_output.get_dihedral_feature_matrix_filename(
                             args, input_line, input_name)),
        log);
  }
  reset_trajectory_state();
}
