#include "titan_rg.h"

#include "titan_executor.h"
#include "titan_gyration.h"
#include "titan_ui.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>

namespace {

const char *trace_kind_name(TitanRgAnalyzer::TraceKind kind) noexcept {
  switch (kind) {
  case TitanRgAnalyzer::TraceKind::CA:
    return "CA";
  case TitanRgAnalyzer::TraceKind::C4Prime:
    return "C4'";
  case TitanRgAnalyzer::TraceKind::MartiniBB:
    return "BB";
  case TitanRgAnalyzer::TraceKind::None:
  default:
    return "none";
  }
}

int trace_priority(TitanRgAnalyzer::TraceKind kind,
                   const std::string &atom_name) noexcept {
  switch (kind) {
  case TitanRgAnalyzer::TraceKind::CA:
    return atom_name == "CA" ? 2 : (atom_name == "CA " ? 1 : -1);
  case TitanRgAnalyzer::TraceKind::C4Prime:
    if (atom_name == "C4'")
      return 2;
    if (atom_name == "C4*")
      return 1;
    return -1;
  case TitanRgAnalyzer::TraceKind::MartiniBB:
    if (atom_name == "BB")
      return 2;
    if (atom_name == "BB2")
      return 1;
    return -1;
  case TitanRgAnalyzer::TraceKind::None:
  default:
    return -1;
  }
}

std::size_t effective_worker_count(std::size_t requested_threads,
                                   std::size_t work_items) noexcept {
  if (requested_threads == 0 || work_items < 512)
    return 1;
  return std::max<std::size_t>(
      1, std::min<std::size_t>(requested_threads, work_items / 256));
}

} // namespace

std::string TitanRgAnalyzer::uppercase_copy(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(),
                 [](unsigned char c) {
                   return static_cast<char>(std::toupper(c));
                 });
  return value;
}

std::string TitanRgAnalyzer::trim_copy(std::string value) {
  const auto first = value.find_first_not_of(" \t");
  if (first == std::string::npos)
    return {};
  const auto last = value.find_last_not_of(" \t");
  return value.substr(first, last - first + 1);
}

bool TitanRgAnalyzer::is_traceable_atom(const std::string &atom_name,
                                        TraceKind trace_kind) noexcept {
  return trace_priority(trace_kind, atom_name) >= 0;
}

void TitanRgAnalyzer::configure(const RgArgs &args) {
  requested_ = args.flag_rg2;
  hist_bins_ = std::max(1, args.hist_bins);
  enabled_ = requested_ && initialized_ && total_trace_count_ > 0;
}

TitanRgAnalyzer::TraceKind
TitanRgAnalyzer::select_trace_kind(const std::string &mol_type,
                                   const std::string &force_field) const {
  const std::string type = uppercase_copy(mol_type);
  const std::string ff = uppercase_copy(force_field);

  if (ff == "MARTINI") {
    if (type == "LIG" || type == "LIGAND")
      return TraceKind::None;
    return TraceKind::MartiniBB;
  }

  if (type == "PROTEIN" || type == "PROTEIN-USER")
    return TraceKind::CA;
  if (type == "RNA" || type == "DNA")
    return TraceKind::C4Prime;
  return TraceKind::None;
}

void TitanRgAnalyzer::build_trace_atoms(const Multi_MOL &mol,
                                        const std::string &mol_type,
                                        const std::string &force_field,
                                        std::vector<TraceAtom> &trace_atoms) {
  trace_atoms.clear();
  const TraceKind trace_kind = select_trace_kind(mol_type, force_field);
  if (trace_kind == TraceKind::None || mol.info.empty())
    return;

  std::size_t prev_chain = std::numeric_limits<std::size_t>::max();
  std::size_t prev_resnr = std::numeric_limits<std::size_t>::max();
  std::size_t ordinal = 0;
  int best_priority = -1;
  TraceAtom current{};
  bool has_residue = false;

  auto flush_current = [&]() {
    if (has_residue && current.atom_index >= 0) {
      current.ordinal = trace_atoms.size() + 1;
      trace_atoms.push_back(current);
    }
  };

  for (std::size_t atom_idx = 0; atom_idx < mol.info.size(); ++atom_idx) {
    const auto &atom = mol.info[atom_idx];
    const std::size_t chainnum = atom.get_chainnum();
    const std::size_t resnr = atom.get_resnr();
    const bool new_residue =
        !has_residue || chainnum != prev_chain || resnr != prev_resnr;

    if (new_residue) {
      flush_current();
      prev_chain = chainnum;
      prev_resnr = resnr;
      ++ordinal;
      has_residue = true;
      best_priority = -1;
      current = TraceAtom{};
      current.ordinal = ordinal;
      current.resnr = resnr;
      current.chainnum = chainnum;
      current.chain = atom.get_chain();
      current.resname = atom.get_res_name();
    }

    const std::string atom_name = uppercase_copy(trim_copy(atom.get_atom()));
    if (!is_traceable_atom(atom_name, trace_kind))
      continue;

    const int priority = trace_priority(trace_kind, atom_name);
    if (priority > best_priority) {
      best_priority = priority;
      current.atom_index = static_cast<int>(atom_idx);
      current.atom_name = trim_copy(atom.get_atom());
    }
  }

  flush_current();
}

void TitanRgAnalyzer::reset_runtime_buffers() {
  raw_x_.assign(total_trace_count_, 0.0);
  raw_y_.assign(total_trace_count_, 0.0);
  raw_z_.assign(total_trace_count_, 0.0);
  reset_trajectory_state();
}

void TitanRgAnalyzer::reset_trajectory_state() {
  rg2_series_.clear();
  rg_series_.clear();
  time_series_ps_.clear();
}

void TitanRgAnalyzer::init(const Multi_MOL &mol1, const Multi_MOL &mol2,
                           const arguments &args, std::fstream &log) {
  nocal_mode_ = args.flag_nocal;
  force_field_ = args.force_f;

  if (initialized_) {
    enabled_ = requested_ && total_trace_count_ > 0;
    return;
  }

  build_trace_atoms(mol1, args.mol_1_type, args.force_f, mol1_trace_);
  if (!nocal_mode_)
    build_trace_atoms(mol2, args.mol_2_type, args.force_f, mol2_trace_);

  total_trace_count_ = mol1_trace_.size() + mol2_trace_.size();
  initialized_ = true;
  enabled_ = requested_ && total_trace_count_ > 0;

  if (total_trace_count_ == 0) {
    TitanUI::print_warning(
        "RG2: no supported trace atoms found for the selected molecule types.");
    if (log.is_open()) {
      TitanUI::log_kv(log, "RG2 Trace mol1", "0");
      TitanUI::log_kv(log, "RG2 Trace mol2", "0");
    }
    return;
  }

  reset_runtime_buffers();

  if (log.is_open()) {
    const TraceKind mol1_kind = select_trace_kind(args.mol_1_type, args.force_f);
    TitanUI::log_kv(log, "RG2 Trace mol1",
                    std::to_string(mol1_trace_.size()) + " (" +
                        trace_kind_name(mol1_kind) + ")");
    if (!nocal_mode_) {
      const TraceKind mol2_kind =
          select_trace_kind(args.mol_2_type, args.force_f);
      TitanUI::log_kv(log, "RG2 Trace mol2",
                      std::to_string(mol2_trace_.size()) + " (" +
                          trace_kind_name(mol2_kind) + ")");
    }
  }
}

void TitanRgAnalyzer::gather_trace_coordinates(const Multi_MOL &mol1,
                                               const Multi_MOL &mol2,
                                               std::size_t snap_idx,
                                               TitanExecutor &executor) {
  const std::size_t mol1_count = mol1_trace_.size();
  const std::size_t total_count = total_trace_count_;
  if (total_count == 0)
    return;

  const MOL &mol1_frame = mol1.m_mol[snap_idx];
  const std::size_t mol2_snap =
      mol2.m_mol.empty()
          ? 0
          : (mol2.m_mol.size() > 1
                 ? std::min<std::size_t>(snap_idx, mol2.m_mol.size() - 1)
                 : 0);
  const MOL *mol2_frame = mol2.m_mol.empty() ? nullptr : &mol2.m_mol[mol2_snap];

  auto gather_range = [&](std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      if (i < mol1_count) {
        const int atom_idx = mol1_trace_[i].atom_index;
        if (atom_idx >= 0 &&
            static_cast<std::size_t>(atom_idx) < mol1_frame.pos.x.size() &&
            static_cast<std::size_t>(atom_idx) < mol1_frame.pos.y.size() &&
            static_cast<std::size_t>(atom_idx) < mol1_frame.pos.z.size()) {
          raw_x_[i] = mol1_frame.pos.x[static_cast<std::size_t>(atom_idx)];
          raw_y_[i] = mol1_frame.pos.y[static_cast<std::size_t>(atom_idx)];
          raw_z_[i] = mol1_frame.pos.z[static_cast<std::size_t>(atom_idx)];
        } else {
          raw_x_[i] = raw_y_[i] = raw_z_[i] = 0.0;
        }
        continue;
      }

      const std::size_t local = i - mol1_count;
      const int atom_idx = mol2_trace_[local].atom_index;
      if (mol2_frame != nullptr && atom_idx >= 0 &&
          static_cast<std::size_t>(atom_idx) < mol2_frame->pos.x.size() &&
          static_cast<std::size_t>(atom_idx) < mol2_frame->pos.y.size() &&
          static_cast<std::size_t>(atom_idx) < mol2_frame->pos.z.size()) {
        raw_x_[i] = mol2_frame->pos.x[static_cast<std::size_t>(atom_idx)];
        raw_y_[i] = mol2_frame->pos.y[static_cast<std::size_t>(atom_idx)];
        raw_z_[i] = mol2_frame->pos.z[static_cast<std::size_t>(atom_idx)];
      } else {
        raw_x_[i] = raw_y_[i] = raw_z_[i] = 0.0;
      }
    }
  };

  const std::size_t worker_count =
      effective_worker_count(executor.capacity(), total_count);
  if (worker_count <= 1) {
    gather_range(0, total_count);
    return;
  }

  executor.parallel_for(total_count, worker_count, 1,
                        [&](std::size_t, std::size_t begin, std::size_t end) {
                          gather_range(begin, end);
                        });
}

void TitanRgAnalyzer::process_frame(const Multi_MOL &mol1, const Multi_MOL &mol2,
                                    std::size_t snap_idx,
                                    TitanExecutor &executor) {
  if (!enabled_ || total_trace_count_ == 0)
    return;
  if (snap_idx >= mol1.m_mol.size())
    return;

  gather_trace_coordinates(mol1, mol2, snap_idx, executor);
  const auto stats = titan_gyration::build_point_stats_soa(
      raw_x_.data(), raw_y_.data(), raw_z_.data(), total_trace_count_);
  const double rg2 = titan_gyration::compute_rg2_from_stats(stats);
  const double rg = titan_gyration::compute_rg_from_stats(stats);

  rg2_series_.push_back(rg2);
  rg_series_.push_back(rg);
  time_series_ps_.push_back(mol1.m_mol[snap_idx].time_of_snap);
}

TitanRgAnalyzer::OutputSummary TitanRgAnalyzer::compute_summary() const {
  return compute_summary_for_series(rg2_series_, rg_series_);
}

TitanRgAnalyzer::OutputSummary TitanRgAnalyzer::compute_summary_for_series(
    const std::vector<double> &rg2_series,
    const std::vector<double> &rg_series) {
  OutputSummary summary;
  summary.frame_count = rg2_series.size();
  if (rg2_series.empty() || rg_series.empty())
    return summary;

  summary.mean_rg2 =
      std::accumulate(rg2_series.begin(), rg2_series.end(), 0.0) /
      static_cast<double>(rg2_series.size());
  summary.mean_rg =
      std::accumulate(rg_series.begin(), rg_series.end(), 0.0) /
      static_cast<double>(rg_series.size());
  summary.min_rg2 = *std::min_element(rg2_series.begin(), rg2_series.end());
  summary.max_rg2 = *std::max_element(rg2_series.begin(), rg2_series.end());
  summary.min_rg = *std::min_element(rg_series.begin(), rg_series.end());
  summary.max_rg = *std::max_element(rg_series.begin(), rg_series.end());

  double var_rg2 = 0.0;
  for (const double value : rg2_series) {
    const double diff = value - summary.mean_rg2;
    var_rg2 += diff * diff;
  }
  double var_rg = 0.0;
  for (const double value : rg_series) {
    const double diff = value - summary.mean_rg;
    var_rg += diff * diff;
  }
  var_rg2 /= static_cast<double>(rg2_series.size());
  var_rg /= static_cast<double>(rg_series.size());
  summary.std_rg2 = std::sqrt(std::max(0.0, var_rg2));
  summary.std_rg = std::sqrt(std::max(0.0, var_rg));
  return summary;
}

void TitanRgAnalyzer::export_timeseries(const std::string &filename,
                                        const OutputSummary &summary,
                                        std::fstream &log) const {
  std::ofstream out(filename);
  if (!out.is_open()) {
    TitanUI::print_warning("RG2: failed to open " + filename);
    return;
  }

  out << "# TITAN Rg2 time series\n";
  out << "# trace_atoms_total " << total_trace_count_ << "\n";
  out << "# frames " << summary.frame_count << "\n";
  out << "# mean_rg2_A2 " << std::fixed << std::setprecision(6)
      << summary.mean_rg2 << "\n";
  out << "# mean_rg_A " << std::fixed << std::setprecision(6)
      << summary.mean_rg << "\n";
  out << "# columns = frame time_ps rg2_A2 rg_A\n";
  for (std::size_t i = 0; i < rg2_series_.size(); ++i) {
    out << std::setw(10) << std::right << (i + 1) << " "
        << std::setw(14) << std::fixed << std::setprecision(4)
        << (i < time_series_ps_.size() ? time_series_ps_[i] : 0.0) << " "
        << std::setw(14) << std::fixed << std::setprecision(6)
        << rg2_series_[i] << " " << std::setw(12) << std::fixed
        << std::setprecision(6) << (i < rg_series_.size() ? rg_series_[i] : 0.0)
        << "\n";
  }

  if (log.is_open())
    TitanUI::log_kv(log, "RG2 Output", filename);
}

void TitanRgAnalyzer::export_summary(const std::string &filename,
                                     const OutputSummary &summary,
                                     std::fstream &log) const {
  std::ofstream out(filename);
  if (!out.is_open()) {
    TitanUI::print_warning("RG2: failed to open " + filename);
    return;
  }

  out << "# columns = n_frames mean_rg2_A2 std_rg2_A2 min_rg2_A2 max_rg2_A2 "
         "mean_rg_A std_rg_A min_rg_A max_rg_A\n";
  out << summary.frame_count << " " << std::fixed << std::setprecision(6)
      << summary.mean_rg2 << " " << summary.std_rg2 << " " << summary.min_rg2
      << " " << summary.max_rg2 << " " << summary.mean_rg << " "
      << summary.std_rg << " " << summary.min_rg << " " << summary.max_rg
      << "\n";

  if (log.is_open())
    TitanUI::log_kv(log, "RG2 Summary Output", filename);
}

void TitanRgAnalyzer::export_histogram(const std::string &filename,
                                       const OutputSummary &summary,
                                       std::fstream &log) const {
  std::ofstream out(filename);
  if (!out.is_open()) {
    TitanUI::print_warning("RG2: failed to open " + filename);
    return;
  }

  const int bins = std::max(1, hist_bins_);
  const double min_rg2 = summary.min_rg2;
  const double max_rg2 = summary.max_rg2;
  const double width =
      (bins > 0 && max_rg2 > min_rg2) ? (max_rg2 - min_rg2) / bins : 1.0;
  std::vector<std::size_t> counts(static_cast<std::size_t>(bins), 0);
  for (const double value : rg2_series_) {
    std::size_t idx = 0;
    if (max_rg2 > min_rg2) {
      const double scaled = (value - min_rg2) / width;
      idx = std::min<std::size_t>(static_cast<std::size_t>(scaled),
                                  counts.size() - 1);
    }
    ++counts[idx];
  }

  out << "# columns = bin_left_A2 bin_center_A2 bin_right_A2 count probability\n";
  for (std::size_t i = 0; i < counts.size(); ++i) {
    const double left = min_rg2 + static_cast<double>(i) * width;
    const double right =
        (max_rg2 > min_rg2) ? (left + width) : min_rg2;
    const double center = (left + right) * 0.5;
    const double probability =
        rg2_series_.empty()
            ? 0.0
            : static_cast<double>(counts[i]) / static_cast<double>(rg2_series_.size());
    out << std::fixed << std::setprecision(6) << left << " " << center << " "
        << right << " " << counts[i] << " " << probability << "\n";
  }

  if (log.is_open())
    TitanUI::log_kv(log, "RG2 Histogram Output", filename);
}

void TitanRgAnalyzer::export_aggregate_timeseries(const std::string &filename,
                                                  std::fstream &log) const {
  std::ofstream out(filename);
  if (!out.is_open()) {
    TitanUI::print_warning("RG2: failed to open " + filename);
    return;
  }

  out << "# TITAN aggregated Rg2 time series\n";
  out << "# trajectories " << aggregate_summaries_.size() << "\n";
  out << "# frames " << aggregate_frames_.size() << "\n";
  out << "# columns = traj_index frame time_ps rg2_A2 rg_A\n";
  for (const auto &frame : aggregate_frames_) {
    out << std::setw(10) << std::right << frame.traj_index << " "
        << std::setw(10) << std::right << frame.frame_index << " "
        << std::setw(14) << std::fixed << std::setprecision(4) << frame.time_ps
        << " " << std::setw(14) << std::fixed << std::setprecision(6)
        << frame.rg2 << " " << std::setw(12) << std::fixed
        << std::setprecision(6) << frame.rg << "\n";
  }

  if (log.is_open())
    TitanUI::log_kv(log, "RG2 Aggregate Output", filename);
}

void TitanRgAnalyzer::export_aggregate_summary(const std::string &filename,
                                               std::fstream &log) const {
  std::ofstream out(filename);
  if (!out.is_open()) {
    TitanUI::print_warning("RG2: failed to open " + filename);
    return;
  }

  out << "# columns = traj_index n_frames mean_rg2_A2 std_rg2_A2 min_rg2_A2 "
         "max_rg2_A2 mean_rg_A std_rg_A min_rg_A max_rg_A label\n";
  for (const auto &record : aggregate_summaries_) {
    const auto &summary = record.summary;
    out << record.traj_index << " " << summary.frame_count << " " << std::fixed
        << std::setprecision(6) << summary.mean_rg2 << " " << summary.std_rg2
        << " " << summary.min_rg2 << " " << summary.max_rg2 << " "
        << summary.mean_rg << " " << summary.std_rg << " " << summary.min_rg
        << " " << summary.max_rg << " " << record.label << "\n";
  }

  if (!aggregate_rg2_series_.empty() && !aggregate_rg_series_.empty()) {
    const OutputSummary global_summary =
        compute_summary_for_series(aggregate_rg2_series_, aggregate_rg_series_);
    out << -1 << " " << global_summary.frame_count << " " << std::fixed
        << std::setprecision(6) << global_summary.mean_rg2 << " "
        << global_summary.std_rg2 << " " << global_summary.min_rg2 << " "
        << global_summary.max_rg2 << " " << global_summary.mean_rg << " "
        << global_summary.std_rg << " " << global_summary.min_rg << " "
        << global_summary.max_rg << " ALL\n";
  }

  if (log.is_open())
    TitanUI::log_kv(log, "RG2 Aggregate Summary", filename);
}

void TitanRgAnalyzer::export_aggregate_histogram(const std::string &filename,
                                                 std::fstream &log) const {
  std::ofstream out(filename);
  if (!out.is_open()) {
    TitanUI::print_warning("RG2: failed to open " + filename);
    return;
  }

  if (aggregate_rg2_series_.empty()) {
    out << "# columns = bin_left_A2 bin_center_A2 bin_right_A2 count probability\n";
    if (log.is_open())
      TitanUI::log_kv(log, "RG2 Aggregate Histogram", filename);
    return;
  }

  const OutputSummary summary =
      compute_summary_for_series(aggregate_rg2_series_, aggregate_rg_series_);
  const int bins = std::max(1, hist_bins_);
  const double min_rg2 = summary.min_rg2;
  const double max_rg2 = summary.max_rg2;
  const double width =
      (bins > 0 && max_rg2 > min_rg2) ? (max_rg2 - min_rg2) / bins : 1.0;
  std::vector<std::size_t> counts(static_cast<std::size_t>(bins), 0);
  for (const double value : aggregate_rg2_series_) {
    std::size_t idx = 0;
    if (max_rg2 > min_rg2) {
      const double scaled = (value - min_rg2) / width;
      idx = std::min<std::size_t>(static_cast<std::size_t>(scaled),
                                  counts.size() - 1);
    }
    ++counts[idx];
  }

  out << "# columns = bin_left_A2 bin_center_A2 bin_right_A2 count probability\n";
  for (std::size_t i = 0; i < counts.size(); ++i) {
    const double left = min_rg2 + static_cast<double>(i) * width;
    const double right = (max_rg2 > min_rg2) ? (left + width) : min_rg2;
    const double center = (left + right) * 0.5;
    const double probability =
        static_cast<double>(counts[i]) /
        static_cast<double>(aggregate_rg2_series_.size());
    out << std::fixed << std::setprecision(6) << left << " " << center << " "
        << right << " " << counts[i] << " " << probability << "\n";
  }

  if (log.is_open())
    TitanUI::log_kv(log, "RG2 Aggregate Histogram", filename);
}

void TitanRgAnalyzer::finalize_trajectory(const std::string &rg2_filename,
                                          const std::string &summary_filename,
                                          const std::string &hist_filename,
                                          std::size_t traj_index,
                                          const std::string &traj_label,
                                          std::fstream &log) {
  if (!enabled_ || rg2_series_.empty() || total_trace_count_ == 0) {
    reset_trajectory_state();
    return;
  }

  const OutputSummary summary = compute_summary();
  for (std::size_t i = 0; i < rg2_series_.size(); ++i) {
    aggregate_frames_.push_back({traj_index, i + 1,
                                 i < time_series_ps_.size() ? time_series_ps_[i]
                                                             : 0.0,
                                 rg2_series_[i],
                                 i < rg_series_.size() ? rg_series_[i] : 0.0});
  }
  aggregate_summaries_.push_back({traj_index, traj_label, summary});
  aggregate_rg2_series_.insert(aggregate_rg2_series_.end(), rg2_series_.begin(),
                               rg2_series_.end());
  aggregate_rg_series_.insert(aggregate_rg_series_.end(), rg_series_.begin(),
                              rg_series_.end());

  export_timeseries(rg2_filename, summary, log);
  export_summary(summary_filename, summary, log);
  export_histogram(hist_filename, summary, log);

  TitanUI::print_info("[RG2] Exported Rg2 to " + rg2_filename);
  TitanUI::print_info("[RG2] Exported summary to " + summary_filename);
  TitanUI::print_info("[RG2] Exported histogram to " + hist_filename);

  reset_trajectory_state();
}

void TitanRgAnalyzer::finalize_all(const std::string &rg2_all_filename,
                                   const std::string &summary_all_filename,
                                   const std::string &hist_all_filename,
                                   std::fstream &log) {
  if (!enabled_ || aggregate_summaries_.size() <= 1 ||
      aggregate_rg2_series_.empty() || aggregate_rg_series_.empty()) {
    return;
  }

  export_aggregate_timeseries(rg2_all_filename, log);
  export_aggregate_summary(summary_all_filename, log);
  export_aggregate_histogram(hist_all_filename, log);

  TitanUI::print_info("[RG2] Exported aggregate Rg2 to " + rg2_all_filename);
  TitanUI::print_info("[RG2] Exported aggregate summary to " +
                      summary_all_filename);
  TitanUI::print_info("[RG2] Exported aggregate histogram to " +
                      hist_all_filename);
}
