#include "titan_contact_map.h"

#include "info.h"
#include "titan_sq_core.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <sstream>
#include <unordered_map>

namespace {

std::string lowercase_copy(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
    return static_cast<char>(std::tolower(c));
  });
  return value;
}

std::string normalize_domain(const std::string &domain) {
  const std::string normalized =
      lowercase_copy(titan_contact::trim_copy(domain));
  if (normalized == "intra" || normalized == "inter" || normalized == "both")
    return normalized;
  throw std::invalid_argument(
      "CMap: -cmap-domain must be one of intra, inter, both.");
}

std::string join_output_path(const std::string &base_dir,
                             const std::string &filename) {
  if (base_dir.empty() || base_dir == ".")
    return filename;
  return (std::filesystem::path(base_dir) / filename).string();
}

bool residue_less(const titan_contact::ResidueKey &lhs,
                  const titan_contact::ResidueKey &rhs) {
  if (lhs.chainnum != rhs.chainnum)
    return lhs.chainnum < rhs.chainnum;
  if (lhs.resnr != rhs.resnr)
    return lhs.resnr < rhs.resnr;
  return lhs.ordinal < rhs.ordinal;
}

void log_cmap_value(std::fstream &log, const std::string &key,
                    const std::string &value) {
  if (log.is_open())
    log << key << ": " << value << '\n';
}

} // namespace

void TitanContactMapAnalyzer::configure(const CMapArgs &args) {
  requested_ = args.flag_cmap;
  mode_ = titan_contact::parse_contact_mode(args.mode);
  domain_ = normalize_domain(args.domain);
  cutoff_a_ = args.cutoff_ang;
  if (!std::isfinite(cutoff_a_) || cutoff_a_ <= 0.0) {
    throw std::invalid_argument("CMap: -cmap-cut must be finite and > 0.");
  }
  min_seqsep_ = args.min_seqsep;
  enabled_ = requested_ && initialized_ && !contact_atoms_.empty() &&
             cutoff_a_ > 0.0;
}

void TitanContactMapAnalyzer::reset_trajectory_state() {
  warned_invalid_box_ = false;
  total_frames_seen_ = 0;
  valid_frames_ = 0;
  skipped_invalid_box_frames_ = 0;
  skipped_incomplete_frames_ = 0;
  records_.clear();
}

void TitanContactMapAnalyzer::init(const Multi_MOL &mol1, const arguments &,
                                   std::fstream &log) {
  contact_atoms_.clear();
  residues_.clear();
  active_atom_mask_.clear();
  atom_index_to_contact_pos_.clear();

  const std::size_t info_atoms = mol1.info.size();
  const std::size_t frame_atoms =
      mol1.m_mol.empty() ? 0 : mol1.m_mol.front().size();
  atom_count_ = (info_atoms > 0 && frame_atoms > 0)
                    ? std::min(info_atoms, frame_atoms)
                    : std::max(info_atoms, frame_atoms);

  active_atom_mask_.assign(atom_count_, 0);
  atom_index_to_contact_pos_.assign(atom_count_, -1);
  std::unordered_map<std::size_t, std::size_t> chain_residue_ordinals;

  for (std::size_t atom_idx = 0; atom_idx < atom_count_ &&
                                 atom_idx < mol1.info.size();
       ++atom_idx) {
    const mol_info &atom_info = mol1.info[atom_idx];
    if (!titan_contact::is_selected_contact_atom(atom_info, mode_))
      continue;

    titan_contact::ResidueKey residue =
        titan_contact::residue_key_from_atom(atom_info, 0);
    std::size_t residue_ordinal = residues_.size();
    auto existing =
        std::find_if(residues_.begin(), residues_.end(), [&](const auto &item) {
          return item == residue;
        });
    if (existing == residues_.end()) {
      residue.ordinal = ++chain_residue_ordinals[residue.chainnum];
      residues_.push_back(residue);
    } else {
      residue_ordinal =
          static_cast<std::size_t>(std::distance(residues_.begin(), existing));
    }

    ContactAtom contact_atom{};
    contact_atom.atom_index = atom_idx;
    contact_atom.residue = residues_[residue_ordinal];
    atom_index_to_contact_pos_[atom_idx] =
        static_cast<int>(contact_atoms_.size());
    contact_atoms_.push_back(std::move(contact_atom));
    active_atom_mask_[atom_idx] = 1;
  }

  initialized_ = true;
  enabled_ = requested_ && !contact_atoms_.empty() && cutoff_a_ > 0.0;
  reset_trajectory_state();

  log_cmap_value(log, "CMap Mode", titan_contact::contact_mode_name(mode_));
  log_cmap_value(log, "CMap Cutoff", std::to_string(cutoff_a_) + " A");
  log_cmap_value(log, "CMap Domain", domain_);
  log_cmap_value(log, "CMap Contact Atoms",
                 std::to_string(contact_atoms_.size()));
  log_cmap_value(log, "CMap Residues", std::to_string(residues_.size()));
}

bool TitanContactMapAnalyzer::has_valid_box(const box &frame_box) const
    noexcept {
  return titan_sq_detail::has_positive_box_lengths(frame_box);
}

bool TitanContactMapAnalyzer::pair_allowed(
    const titan_contact::ResidueKey &lhs,
    const titan_contact::ResidueKey &rhs) const {
  if (titan_contact::same_residue(lhs, rhs))
    return false;

  const bool same_chain = titan_contact::same_chain(lhs, rhs);
  if (domain_ == "intra" && !same_chain)
    return false;
  if (domain_ == "inter" && same_chain)
    return false;

  if (same_chain && min_seqsep_ > 0) {
    const std::size_t diff = lhs.ordinal > rhs.ordinal
                                 ? lhs.ordinal - rhs.ordinal
                                 : rhs.ordinal - lhs.ordinal;
    if (diff < static_cast<std::size_t>(min_seqsep_))
      return false;
  }

  return true;
}

void TitanContactMapAnalyzer::process_frame(const Multi_MOL &mol1,
                                            std::size_t snap_idx,
                                            TitanExecutor & /*executor*/) {
  if (!enabled_ || snap_idx >= mol1.m_mol.size())
    return;

  ++total_frames_seen_;
  const MOL &frame = mol1.m_mol[snap_idx];
  if (!has_valid_box(frame.Box)) {
    ++skipped_invalid_box_frames_;
    warned_invalid_box_ = true;
    return;
  }

  if (contact_atoms_.size() < 2)
    return;

  const auto xs = frame.get_pos_x();
  const auto ys = frame.get_pos_y();
  const auto zs = frame.get_pos_z();
  const std::size_t atom_limit =
      std::min({xs.size(), ys.size(), zs.size()});
  const bool all_contact_atoms_present = std::all_of(
      contact_atoms_.begin(), contact_atoms_.end(),
      [&](const ContactAtom &atom) { return atom.atom_index < atom_limit; });
  if (!all_contact_atoms_present) {
    ++skipped_incomplete_frames_;
    return;
  }
  ++valid_frames_;
  if (atom_limit < 2)
    return;

  std::unordered_map<titan_contact::ResiduePairKey, double,
                     titan_contact::ResiduePairKeyHash>
      frame_nearest;

  titan_sq_detail::for_each_active_pair_within_cutoff(
      xs, ys, zs, frame.Box, active_atom_mask_, cutoff_a_, neighbor_grid_,
      [&](std::size_t lhs, std::size_t rhs, double distance_a) {
        const int lhs_pos =
            lhs < atom_index_to_contact_pos_.size()
                ? atom_index_to_contact_pos_[lhs]
                : -1;
        const int rhs_pos =
            rhs < atom_index_to_contact_pos_.size()
                ? atom_index_to_contact_pos_[rhs]
                : -1;
        if (lhs_pos < 0 || rhs_pos < 0)
          return;

        const auto &lhs_residue =
            contact_atoms_[static_cast<std::size_t>(lhs_pos)].residue;
        const auto &rhs_residue =
            contact_atoms_[static_cast<std::size_t>(rhs_pos)].residue;
        if (!pair_allowed(lhs_residue, rhs_residue))
          return;

        const titan_contact::ResiduePairKey key =
            titan_contact::make_ordered_pair(lhs_residue, rhs_residue);
        auto found = frame_nearest.find(key);
        if (found == frame_nearest.end() || distance_a < found->second)
          frame_nearest[key] = distance_a;
      });

  for (const auto &[key, distance_a] : frame_nearest) {
    auto [it, inserted] = records_.try_emplace(key);
    ContactRecord &record = it->second;
    if (inserted) {
      record.key = key;
      record.min_distance_a = distance_a;
    } else {
      record.min_distance_a = std::min(record.min_distance_a, distance_a);
    }
    ++record.contact_frames;
    record.distance_sum_a += distance_a;
    ++record.distance_samples;
  }
}

std::string TitanContactMapAnalyzer::residue_label(
    const titan_contact::ResidueKey &residue) const {
  std::ostringstream out;
  out << residue.chain << ':' << residue.resnr << ':' << residue.resname;
  return out.str();
}

void TitanContactMapAnalyzer::export_summary(const std::string &filename,
                                             std::fstream &log) const {
  std::ofstream out(filename);
  if (!out) {
    log_cmap_value(log, "CMap Summary", "failed to open " + filename);
    return;
  }
  out << std::fixed << std::setprecision(6);
  out << "module: titan_cmap\n";
  out << "mode: " << titan_contact::contact_mode_name(mode_) << '\n';
  out << "cutoff_A: " << cutoff_a_ << '\n';
  out << "domain: " << domain_ << '\n';
  out << "min_seqsep: " << min_seqsep_ << '\n';
  out << "frames_seen: " << total_frames_seen_ << '\n';
  out << "valid_frames: " << valid_frames_ << '\n';
  out << "skipped_invalid_box_frames: " << skipped_invalid_box_frames_ << '\n';
  out << "skipped_incomplete_frames: " << skipped_incomplete_frames_ << '\n';
  out << "contact_atoms: " << contact_atoms_.size() << '\n';
  out << "residues: " << residues_.size() << '\n';
  out << "contact_edges: " << records_.size() << '\n';
  log_cmap_value(log, "CMap Summary", filename);
}

void TitanContactMapAnalyzer::export_edges(const std::string &filename,
                                           std::fstream &log) const {
  std::ofstream out(filename);
  if (!out) {
    log_cmap_value(log, "CMap Edges", "failed to open " + filename);
    return;
  }
  out << "chain_i,resnr_i,resname_i,chain_j,resnr_j,resname_j,probability,"
         "contact_frames,mean_distance_A,min_distance_A\n";

  std::vector<ContactRecord> records;
  records.reserve(records_.size());
  for (const auto &[_, record] : records_)
    records.push_back(record);
  std::sort(records.begin(), records.end(),
            [](const ContactRecord &lhs, const ContactRecord &rhs) {
              if (residue_less(lhs.key.first, rhs.key.first))
                return true;
              if (residue_less(rhs.key.first, lhs.key.first))
                return false;
              return residue_less(lhs.key.second, rhs.key.second);
            });

  out << std::fixed << std::setprecision(6);
  for (const ContactRecord &record : records) {
    const double probability =
        valid_frames_ > 0
            ? static_cast<double>(record.contact_frames) /
                  static_cast<double>(valid_frames_)
            : 0.0;
    const double mean_distance =
        record.distance_samples > 0
            ? record.distance_sum_a / static_cast<double>(record.distance_samples)
            : 0.0;
    out << record.key.first.chain << ',' << record.key.first.resnr << ','
        << record.key.first.resname << ',' << record.key.second.chain << ','
        << record.key.second.resnr << ',' << record.key.second.resname << ','
        << probability << ',' << record.contact_frames << ',' << mean_distance
        << ',' << record.min_distance_a << '\n';
  }
  log_cmap_value(log, "CMap Edges", filename);
}

void TitanContactMapAnalyzer::export_matrix(const std::string &filename,
                                            std::fstream &log) const {
  std::ofstream out(filename);
  if (!out) {
    log_cmap_value(log, "CMap Matrix", "failed to open " + filename);
    return;
  }

  std::unordered_map<titan_contact::ResiduePairKey, double,
                     titan_contact::ResiduePairKeyHash>
      probabilities;
  for (const auto &[key, record] : records_) {
    probabilities[key] =
        valid_frames_ > 0
            ? static_cast<double>(record.contact_frames) /
                  static_cast<double>(valid_frames_)
            : 0.0;
  }

  out << "residue";
  for (const auto &residue : residues_)
    out << ',' << residue_label(residue);
  out << '\n';

  out << std::fixed << std::setprecision(6);
  for (const auto &row_residue : residues_) {
    out << residue_label(row_residue);
    for (const auto &col_residue : residues_) {
      double probability = 0.0;
      if (!titan_contact::same_residue(row_residue, col_residue)) {
        const auto key = titan_contact::make_ordered_pair(row_residue, col_residue);
        auto found = probabilities.find(key);
        if (found != probabilities.end())
          probability = found->second;
      }
      out << ',' << probability;
    }
    out << '\n';
  }
  log_cmap_value(log, "CMap Matrix", filename);
}

void TitanContactMapAnalyzer::finalize_trajectory(
    const std::string &base_dir, const info &info_output, const arguments &args,
    const std::vector<std::string> &input_line, int input_name,
    std::fstream &log) {
  if (!initialized_ || !requested_)
    return;

  const std::string summary_filename = join_output_path(
      base_dir, info_output.get_cmap_summary_filename(args, input_line, input_name));
  const std::string edges_filename = join_output_path(
      base_dir, info_output.get_cmap_edges_filename(args, input_line, input_name));
  const std::string matrix_filename = join_output_path(
      base_dir, info_output.get_cmap_matrix_filename(args, input_line, input_name));

  export_summary(summary_filename, log);
  export_edges(edges_filename, log);
  export_matrix(matrix_filename, log);
  reset_trajectory_state();
}
