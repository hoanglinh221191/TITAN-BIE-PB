#include "titan_llps_cluster.h"
#include "BoxMatrix.h"
#include "box_type.h"
#include "info.h"
#include "titan_ui.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

namespace {

constexpr double kPi = 3.14159265358979323846264338327950288;
constexpr double kEps = 1.0e-12;

double deg_to_rad(double v) { return v * kPi / 180.0; }

double box_volume_nm3(const box &Box) {
  const double a = Box.get_a();
  const double b = Box.get_b();
  const double c = Box.get_c();
  if (a <= 0.0 || b <= 0.0 || c <= 0.0)
    return 0.0;
  const double ca = std::cos(deg_to_rad(Box.get_alpha()));
  const double cb = std::cos(deg_to_rad(Box.get_beta()));
  const double cg = std::cos(deg_to_rad(Box.get_gamma()));
  const double metric =
      1.0 + 2.0 * ca * cb * cg - ca * ca - cb * cb - cg * cg;
  return a * b * c * std::sqrt(std::max(0.0, metric)) / 1000.0;
}

uint64_t choose2(uint64_t n) { return (n > 1) ? (n * (n - 1) / 2) : 0; }

double safe_partition(double dense, double dilute) {
  if (dense <= 0.0 && dilute <= 0.0)
    return 0.0;
  return (dense + kEps) / (dilute + kEps);
}

double transfer_proxy(double partition) {
  if (partition <= 0.0)
    return 0.0;
  return -std::log(partition);
}

std::string region_label(double center_nm, double core_limit_nm,
                         double interface_limit_nm) {
  if (center_nm <= core_limit_nm)
    return "core";
  if (center_nm <= interface_limit_nm)
    return "interface";
  return "dilute";
}

struct BlockStats {
  int blocks{0};
  double mean{0.0};
  double sem{0.0};
};

template <typename Fn>
BlockStats block_reduce(const std::vector<LLPSPaperAccumulator::FrameRow> &rows,
                        int requested_blocks, Fn value_of) {
  BlockStats out;
  if (rows.empty())
    return out;
  int blocks = requested_blocks > 0 ? requested_blocks : 5;
  blocks = std::min<int>(blocks, static_cast<int>(rows.size()));
  if (rows.size() < 20 && requested_blocks <= 0)
    blocks = std::min<int>(blocks, std::max<int>(1, rows.size() / 2));
  blocks = std::max(1, blocks);
  out.blocks = blocks;

  std::vector<double> block_means;
  block_means.reserve(static_cast<std::size_t>(blocks));
  for (int b = 0; b < blocks; ++b) {
    const std::size_t begin =
        static_cast<std::size_t>((static_cast<uint64_t>(b) * rows.size()) /
                                 static_cast<uint64_t>(blocks));
    const std::size_t end =
        static_cast<std::size_t>((static_cast<uint64_t>(b + 1) * rows.size()) /
                                 static_cast<uint64_t>(blocks));
    if (end <= begin)
      continue;
    double sum = 0.0;
    for (std::size_t i = begin; i < end; ++i)
      sum += value_of(rows[i]);
    block_means.push_back(sum / static_cast<double>(end - begin));
  }
  if (block_means.empty())
    return out;
  out.mean = std::accumulate(block_means.begin(), block_means.end(), 0.0) /
             static_cast<double>(block_means.size());
  if (block_means.size() > 1) {
    double var = 0.0;
    for (double v : block_means) {
      const double d = v - out.mean;
      var += d * d;
    }
    var /= static_cast<double>(block_means.size() - 1);
    out.sem =
        std::sqrt(var) / std::sqrt(static_cast<double>(block_means.size()));
  }
  return out;
}

double drift_slope_per_ns(
    const std::vector<LLPSPaperAccumulator::FrameRow> &rows) {
  if (rows.size() < 2)
    return 0.0;
  double sx = 0.0;
  double sy = 0.0;
  double sxx = 0.0;
  double sxy = 0.0;
  for (const auto &row : rows) {
    sx += row.time_ns;
    sy += row.dense_fraction;
    sxx += row.time_ns * row.time_ns;
    sxy += row.time_ns * row.dense_fraction;
  }
  const double n = static_cast<double>(rows.size());
  const double denom = n * sxx - sx * sx;
  if (std::abs(denom) < kEps)
    return 0.0;
  return (n * sxy - sx * sy) / denom;
}

} // namespace

void LLPSPaperAccumulator::RunningScalar::add(double v) {
  ++n;
  sum += v;
  sum_sq += v * v;
}

double LLPSPaperAccumulator::RunningScalar::mean() const {
  return n > 0 ? sum / static_cast<double>(n) : 0.0;
}

double LLPSPaperAccumulator::RunningScalar::sem() const {
  if (n < 2)
    return 0.0;
  const double m = mean();
  const double var =
      std::max(0.0, sum_sq / static_cast<double>(n) - m * m);
  return std::sqrt(var / static_cast<double>(n));
}

void LLPSPaperAccumulator::configure(const LLPSPaperConfig &config) {
  config_ = config;
  config_.profile_bin_nm = std::max(0.01, config_.profile_bin_nm);
  if (config_.geometry != "droplet" && config_.geometry != "slab")
    config_.geometry = "auto";
  frames_.clear();
  contact_rows_.clear();
  profile_bins_.clear();
  skipped_no_dense_ = 0;
  skipped_bad_box_ = 0;
  observed_geometry_ = "none";
}

void LLPSPaperAccumulator::process_frame(
    double time_ps, int frame_idx, const std::vector<LLPSNode> &nodes,
    const std::vector<bool> &is_dense, int giant_cluster_id, const box &Box) {
  if (!config_.enabled)
    return;
  const double volume_nm3 = box_volume_nm3(Box);
  if (volume_nm3 <= 0.0 || nodes.empty() || is_dense.size() != nodes.size()) {
    ++skipped_bad_box_;
    return;
  }

  const double box_nm[3] = {Box.get_a() / 10.0, Box.get_b() / 10.0,
                            Box.get_c() / 10.0};
  int axis = 0;
  if (box_nm[1] > box_nm[axis])
    axis = 1;
  if (box_nm[2] > box_nm[axis])
    axis = 2;

  FrameRow row;
  row.time_ns = time_ps / 1000.0;
  row.frame = frame_idx;
  row.box_volume_nm3 = volume_nm3;

  int first_dense = -1;
  for (std::size_t i = 0; i < nodes.size(); ++i) {
    const bool protein = nodes[i].type == 1;
    row.total_protein_nodes += protein ? 1 : 0;
    row.total_rna_nodes += protein ? 0 : 1;
    if (is_dense[i]) {
      if (first_dense < 0)
        first_dense = static_cast<int>(i);
      row.dense_protein_nodes += protein ? 1 : 0;
      row.dense_rna_nodes += protein ? 0 : 1;
    }
  }
  row.dilute_protein_nodes = row.total_protein_nodes - row.dense_protein_nodes;
  row.dilute_rna_nodes = row.total_rna_nodes - row.dense_rna_nodes;
  row.dense_fraction =
      static_cast<double>(row.dense_protein_nodes + row.dense_rna_nodes) /
      static_cast<double>(nodes.size());
  row.dilute_fraction = 1.0 - row.dense_fraction;

  if (first_dense < 0 || giant_cluster_id <= 0) {
    ++skipped_no_dense_;
    row.geometry = "none";
    row.dilute_volume_nm3 = volume_nm3;
    if (row.dilute_volume_nm3 > 0.0) {
      row.dilute_protein_conc =
          static_cast<double>(row.dilute_protein_nodes) / row.dilute_volume_nm3;
      row.dilute_rna_conc =
          static_cast<double>(row.dilute_rna_nodes) / row.dilute_volume_nm3;
    }
    frames_.push_back(row);
    return;
  }

  const BoxMatrix bmat(Box);
  const LLPSNode &ref = nodes[static_cast<std::size_t>(first_dense)];
  std::vector<std::array<double, 3>> dense_xyz;
  dense_xyz.reserve(nodes.size());
  std::array<double, 3> min_nm = {std::numeric_limits<double>::max(),
                                  std::numeric_limits<double>::max(),
                                  std::numeric_limits<double>::max()};
  std::array<double, 3> max_nm = {std::numeric_limits<double>::lowest(),
                                  std::numeric_limits<double>::lowest(),
                                  std::numeric_limits<double>::lowest()};
  double cx = 0.0;
  double cy = 0.0;
  double cz = 0.0;
  for (std::size_t i = 0; i < nodes.size(); ++i) {
    if (!is_dense[i])
      continue;
    double dx = nodes[i].x - ref.x;
    double dy = nodes[i].y - ref.y;
    double dz = nodes[i].z - ref.z;
    bmat.apply_mic(dx, dy, dz);
    const double ux = ref.x + dx;
    const double uy = ref.y + dy;
    const double uz = ref.z + dz;
    dense_xyz.push_back({ux, uy, uz});
    cx += ux;
    cy += uy;
    cz += uz;
    const double nm[3] = {ux / 10.0, uy / 10.0, uz / 10.0};
    for (int d = 0; d < 3; ++d) {
      min_nm[d] = std::min(min_nm[d], nm[d]);
      max_nm[d] = std::max(max_nm[d], nm[d]);
    }
  }
  const double inv_dense = 1.0 / static_cast<double>(dense_xyz.size());
  cx *= inv_dense;
  cy *= inv_dense;
  cz *= inv_dense;

  std::array<double, 3> extent_nm = {max_nm[0] - min_nm[0],
                                     max_nm[1] - min_nm[1],
                                     max_nm[2] - min_nm[2]};
  const double max_extent =
      std::max(extent_nm[0], std::max(extent_nm[1], extent_nm[2]));
  const double min_extent =
      std::max(0.001, std::min(extent_nm[0], std::min(extent_nm[1], extent_nm[2])));
  std::string geometry = config_.geometry;
  if (geometry == "auto") {
    const bool spans_long_axis =
        box_nm[axis] > 0.0 && extent_nm[axis] > 0.55 * box_nm[axis];
    geometry = (spans_long_axis || max_extent / min_extent > 2.2) ? "slab"
                                                                  : "droplet";
  }
  row.geometry = geometry;
  if (observed_geometry_ == "none")
    observed_geometry_ = geometry;

  double radius_nm = 0.0;
  double max_profile_distance_nm = 0.0;
  std::vector<double> dense_distance_nm(nodes.size(), 0.0);
  for (std::size_t i = 0; i < nodes.size(); ++i) {
    double dx = nodes[i].x - cx;
    double dy = nodes[i].y - cy;
    double dz = nodes[i].z - cz;
    bmat.apply_mic(dx, dy, dz);
    const double r = std::sqrt(dx * dx + dy * dy + dz * dz) / 10.0;
    const double slab_d =
        std::abs((axis == 0 ? dx : (axis == 1 ? dy : dz))) / 10.0;
    const double dist = (geometry == "slab") ? slab_d : r;
    dense_distance_nm[i] = dist;
    max_profile_distance_nm = std::max(max_profile_distance_nm, dist);
    if (is_dense[i])
      radius_nm = std::max(radius_nm, dist);
  }
  radius_nm = std::max(radius_nm, config_.profile_bin_nm);
  row.radius_nm = radius_nm;

  if (geometry == "slab") {
    const double thickness = std::max(config_.profile_bin_nm, extent_nm[axis]);
    double cross_section = 1.0;
    for (int d = 0; d < 3; ++d) {
      if (d != axis)
        cross_section *= std::max(0.0, box_nm[d]);
    }
    row.dense_volume_nm3 =
        std::min(volume_nm3, std::max(config_.profile_bin_nm * cross_section,
                                      thickness * cross_section));
    row.interface_width_nm = std::max(config_.profile_bin_nm, 0.25 * thickness);
  } else {
    row.dense_volume_nm3 =
        std::min(volume_nm3, (4.0 / 3.0) * kPi * radius_nm * radius_nm *
                                 radius_nm);
    row.interface_width_nm = std::max(config_.profile_bin_nm, 0.20 * radius_nm);
  }
  row.dilute_volume_nm3 = std::max(kEps, volume_nm3 - row.dense_volume_nm3);
  row.dense_volume_nm3 = std::max(kEps, row.dense_volume_nm3);
  row.dense_protein_conc =
      static_cast<double>(row.dense_protein_nodes) / row.dense_volume_nm3;
  row.dilute_protein_conc =
      static_cast<double>(row.dilute_protein_nodes) / row.dilute_volume_nm3;
  row.dense_rna_conc =
      static_cast<double>(row.dense_rna_nodes) / row.dense_volume_nm3;
  row.dilute_rna_conc =
      static_cast<double>(row.dilute_rna_nodes) / row.dilute_volume_nm3;
  row.protein_partition =
      safe_partition(row.dense_protein_conc, row.dilute_protein_conc);
  row.rna_partition = safe_partition(row.dense_rna_conc, row.dilute_rna_conc);
  row.protein_transfer_proxy = transfer_proxy(row.protein_partition);
  row.rna_transfer_proxy = transfer_proxy(row.rna_partition);
  row.dense_protein_rna_ratio =
      static_cast<double>(row.dense_protein_nodes) /
      std::max(1.0, static_cast<double>(row.dense_rna_nodes));
  const double global_rna_fraction =
      static_cast<double>(row.total_rna_nodes) /
      std::max(1.0, static_cast<double>(row.total_protein_nodes +
                                        row.total_rna_nodes));
  const double dense_rna_fraction =
      static_cast<double>(row.dense_rna_nodes) /
      std::max(1.0, static_cast<double>(row.dense_protein_nodes +
                                        row.dense_rna_nodes));
  row.rna_enrichment = safe_partition(dense_rna_fraction, global_rna_fraction);

  uint64_t core_p = 0;
  uint64_t core_r = 0;
  uint64_t interface_p = 0;
  uint64_t interface_r = 0;
  const double core_limit = 0.60 * radius_nm;
  for (std::size_t i = 0; i < nodes.size(); ++i) {
    if (!is_dense[i])
      continue;
    const bool is_core = dense_distance_nm[i] <= core_limit;
    if (is_core) {
      if (nodes[i].type == 1)
        ++core_p;
      else
        ++core_r;
    } else {
      if (nodes[i].type == 1)
        ++interface_p;
      else
        ++interface_r;
    }
  }
  const double dense_nodes =
      static_cast<double>(row.dense_protein_nodes + row.dense_rna_nodes);
  row.core_fraction = static_cast<double>(core_p + core_r) /
                      std::max(1.0, dense_nodes);
  row.interface_fraction =
      static_cast<double>(interface_p + interface_r) / std::max(1.0, dense_nodes);

  ContactFrameRow contact;
  contact.time_ns = row.time_ns;
  contact.frame = frame_idx;
  const uint64_t dense_p = static_cast<uint64_t>(row.dense_protein_nodes);
  const uint64_t dense_r = static_cast<uint64_t>(row.dense_rna_nodes);
  contact.pp = choose2(dense_p);
  contact.prna = dense_p * dense_r;
  contact.rr = choose2(dense_r);
  contact.pp_core = choose2(core_p);
  contact.prna_core = core_p * core_r;
  contact.rr_core = choose2(core_r);
  contact.pp_interface = choose2(interface_p);
  contact.prna_interface = interface_p * interface_r;
  contact.rr_interface = choose2(interface_r);
  contact_rows_.push_back(contact);

  const std::size_t needed_bins =
      static_cast<std::size_t>(std::floor(
          std::max(radius_nm, max_profile_distance_nm) /
          config_.profile_bin_nm)) +
      2;
  if (profile_bins_.size() < needed_bins)
    profile_bins_.resize(needed_bins);
  for (std::size_t i = 0; i < nodes.size(); ++i) {
    const std::size_t bin = std::min<std::size_t>(
        profile_bins_.size() - 1,
        static_cast<std::size_t>(dense_distance_nm[i] /
                                 config_.profile_bin_nm));
    auto &dst = profile_bins_[bin];
    if (nodes[i].type == 1) {
      dst.protein_count += 1.0;
      dst.sticker_count += 1.0;
      if (is_dense[i])
        dst.dense_protein_count += 1.0;
    } else {
      dst.rna_count += 1.0;
      if (is_dense[i])
        dst.dense_rna_count += 1.0;
    }
  }
  for (std::size_t bin = 0; bin < profile_bins_.size(); ++bin) {
    auto &dst = profile_bins_[bin];
    const double r0 = static_cast<double>(bin) * config_.profile_bin_nm;
    const double r1 = r0 + config_.profile_bin_nm;
    double shell_volume = 0.0;
    if (geometry == "slab") {
      double cross_section = 1.0;
      for (int d = 0; d < 3; ++d) {
        if (d != axis)
          cross_section *= std::max(0.0, box_nm[d]);
      }
      shell_volume = 2.0 * config_.profile_bin_nm * cross_section;
      if (r0 >= 0.5 * box_nm[axis])
        shell_volume = 0.0;
    } else {
      shell_volume = (4.0 / 3.0) * kPi * (r1 * r1 * r1 - r0 * r0 * r0);
    }
    if (shell_volume > 0.0) {
      dst.volume_nm3_sum += shell_volume;
      ++dst.samples;
    }
  }

  frames_.push_back(row);
}

void LLPSPaperAccumulator::export_outputs(const info &info_output,
                                          std::fstream &log) const {
  if (!config_.enabled)
    return;
  const std::string geometry =
      (observed_geometry_ == "slab") ? "slab" : "droplet";
  const std::string profile_file =
      (geometry == "slab") ? info_output.get_llps_slab_profiles_filename()
                            : info_output.get_llps_radial_profiles_filename();
  const int requested_blocks = config_.requested_blocks;

  const auto protein_dense =
      block_reduce(frames_, requested_blocks,
                   [](const FrameRow &r) { return r.dense_protein_conc; });
  const auto protein_dilute =
      block_reduce(frames_, requested_blocks,
                   [](const FrameRow &r) { return r.dilute_protein_conc; });
  const auto rna_dense =
      block_reduce(frames_, requested_blocks,
                   [](const FrameRow &r) { return r.dense_rna_conc; });
  const auto rna_dilute =
      block_reduce(frames_, requested_blocks,
                   [](const FrameRow &r) { return r.dilute_rna_conc; });
  const auto protein_partition =
      block_reduce(frames_, requested_blocks,
                   [](const FrameRow &r) { return r.protein_partition; });
  const auto rna_partition =
      block_reduce(frames_, requested_blocks,
                   [](const FrameRow &r) { return r.rna_partition; });
  const auto dense_ratio =
      block_reduce(frames_, requested_blocks,
                   [](const FrameRow &r) { return r.dense_protein_rna_ratio; });
  const auto rna_enrichment =
      block_reduce(frames_, requested_blocks,
                   [](const FrameRow &r) { return r.rna_enrichment; });
  const auto radius =
      block_reduce(frames_, requested_blocks,
                   [](const FrameRow &r) { return r.radius_nm; });
  const auto dense_fraction =
      block_reduce(frames_, requested_blocks,
                   [](const FrameRow &r) { return r.dense_fraction; });

  {
    std::ofstream out(info_output.get_llps_phase_partition_summary_filename());
    out << std::fixed << std::setprecision(8);
    out << "# TITAN LLPS paper phase/partition summary\n";
    out << "# Concentrations are node-count densities in nm^-3; transfer "
           "proxies are -ln(K) in kBT-like units.\n";
    out << "Metric Mean SEM Blocks Unit\n";
    auto emit = [&](const std::string &name, const BlockStats &s,
                    const std::string &unit) {
      out << name << " " << s.mean << " " << s.sem << " " << s.blocks << " "
          << unit << "\n";
    };
    emit("protein_dense_concentration", protein_dense, "nodes_nm^-3");
    emit("protein_dilute_concentration", protein_dilute, "nodes_nm^-3");
    emit("rna_dense_concentration", rna_dense, "nodes_nm^-3");
    emit("rna_dilute_concentration", rna_dilute, "nodes_nm^-3");
    emit("protein_partition_coefficient", protein_partition, "unitless");
    emit("rna_partition_coefficient", rna_partition, "unitless");
    emit("protein_transfer_free_energy_proxy",
         block_reduce(frames_, requested_blocks,
                      [](const FrameRow &r) {
                        return r.protein_transfer_proxy;
                      }),
         "kBT_proxy");
    emit("rna_transfer_free_energy_proxy",
         block_reduce(frames_, requested_blocks,
                      [](const FrameRow &r) { return r.rna_transfer_proxy; }),
         "kBT_proxy");
    emit("dense_protein_rna_ratio", dense_ratio, "node_ratio");
    emit("rna_enrichment", rna_enrichment, "unitless");
    emit("dense_fraction", dense_fraction, "fraction");
    emit("condensate_radius_or_halfwidth", radius, "nm");
  }

  {
    std::ofstream out(profile_file);
    out << std::fixed << std::setprecision(8);
    out << "BinStart_nm,BinEnd_nm,BinCenter_nm,Geometry,Frames,"
           "ShellVolumeMean_nm3,ProteinDensity,RNADensity,StickerDensity,"
           "WaterDensity,IonDensity,DenseProteinDensity,DenseRNADensity,"
           "RegionLabel\n";
    const double core_limit = 0.60 * radius.mean;
    const double interface_limit = std::max(core_limit, radius.mean);
    for (std::size_t i = 0; i < profile_bins_.size(); ++i) {
      const auto &bin = profile_bins_[i];
      const double r0 = static_cast<double>(i) * config_.profile_bin_nm;
      const double r1 = r0 + config_.profile_bin_nm;
      const double center = 0.5 * (r0 + r1);
      const double volume =
          bin.samples > 0 ? bin.volume_nm3_sum / static_cast<double>(bin.samples)
                          : 0.0;
      const double norm = std::max(kEps, bin.volume_nm3_sum);
      out << r0 << "," << r1 << "," << center << "," << geometry << ","
          << bin.samples << "," << volume << ","
          << bin.protein_count / norm << "," << bin.rna_count / norm << ","
          << bin.sticker_count / norm << "," << bin.water_count / norm << ","
          << bin.ion_count / norm << "," << bin.dense_protein_count / norm
          << "," << bin.dense_rna_count / norm << ","
          << region_label(center, core_limit, interface_limit) << "\n";
    }
  }

  {
    std::ofstream out(info_output.get_llps_composition_timeseries_filename());
    out << std::fixed << std::setprecision(8);
    out << "Time_ns,Frame,Geometry,TotalProteinNodes,TotalRNANodes,"
           "DenseProteinNodes,DenseRNANodes,DiluteProteinNodes,"
           "DiluteRNANodes,DenseFraction,CoreFraction,InterfaceFraction,"
           "DiluteFraction,DenseProteinConc,DenseRNAConc,DiluteProteinConc,"
           "DiluteRNAConc,ProteinPartition,RNAPartition,"
           "DenseProteinRNARatio,RNAEnrichment,RadiusOrHalfWidth_nm,"
           "InterfaceWidth_nm\n";
    for (const auto &r : frames_) {
      out << r.time_ns << "," << r.frame << "," << r.geometry << ","
          << r.total_protein_nodes << "," << r.total_rna_nodes << ","
          << r.dense_protein_nodes << "," << r.dense_rna_nodes << ","
          << r.dilute_protein_nodes << "," << r.dilute_rna_nodes << ","
          << r.dense_fraction << "," << r.core_fraction << ","
          << r.interface_fraction << "," << r.dilute_fraction << ","
          << r.dense_protein_conc << "," << r.dense_rna_conc << ","
          << r.dilute_protein_conc << "," << r.dilute_rna_conc << ","
          << r.protein_partition << "," << r.rna_partition << ","
          << r.dense_protein_rna_ratio << "," << r.rna_enrichment << ","
          << r.radius_nm << "," << r.interface_width_nm << "\n";
    }
  }

  {
    std::ofstream out(info_output.get_llps_contact_taxonomy_filename());
    out << "Time_ns,Frame,Metric,PP,PRNA,RR,PP_Core,PRNA_Core,RR_Core,"
           "PP_Interface,PRNA_Interface,RR_Interface\n";
    out << std::fixed << std::setprecision(8);
    for (const auto &r : contact_rows_) {
      out << r.time_ns << "," << r.frame << "," << r.metric << "," << r.pp
          << "," << r.prna << "," << r.rr << "," << r.pp_core << ","
          << r.prna_core << "," << r.rr_core << "," << r.pp_interface << ","
          << r.prna_interface << "," << r.rr_interface << "\n";
    }
  }

  {
    RunningScalar pp;
    RunningScalar prna;
    RunningScalar rr;
    for (const auto &r : contact_rows_) {
      pp.add(static_cast<double>(r.pp));
      prna.add(static_cast<double>(r.prna));
      rr.add(static_cast<double>(r.rr));
    }
    std::ofstream out(info_output.get_llps_contact_taxonomy_summary_filename());
    out << std::fixed << std::setprecision(8);
    out << "# Contact taxonomy currently reports dense-cluster co-membership "
           "pair proxies from DBSCAN labels.\n";
    out << "Metric Mean SEM Frames\n";
    out << "PP_cluster_comembership_proxy " << pp.mean() << " " << pp.sem()
        << " " << pp.n << "\n";
    out << "PRNA_cluster_comembership_proxy " << prna.mean() << " "
        << prna.sem() << " " << prna.n << "\n";
    out << "RR_cluster_comembership_proxy " << rr.mean() << " " << rr.sem()
        << " " << rr.n << "\n";
  }

  const double slope = drift_slope_per_ns(frames_);
  const double span_ns = frames_.empty()
                             ? 0.0
                             : (frames_.back().time_ns - frames_.front().time_ns);
  const bool insufficient_blocks = dense_fraction.blocks < 2;
  const bool drift_warning = std::abs(slope * span_ns) > 0.10;
  const bool short_warning = frames_.size() < 10;
  {
    std::ofstream out(info_output.get_llps_convergence_report_filename());
    out << std::fixed << std::setprecision(8);
    out << "frames_used " << frames_.size() << "\n";
    out << "frames_skipped_no_stable_dense_phase " << skipped_no_dense_ << "\n";
    out << "frames_skipped_bad_box " << skipped_bad_box_ << "\n";
    out << "blocks " << dense_fraction.blocks << "\n";
    out << "dense_fraction_mean " << dense_fraction.mean << "\n";
    out << "dense_fraction_sem " << dense_fraction.sem << "\n";
    out << "dense_fraction_drift_slope_per_ns " << slope << "\n";
    out << "time_span_ns " << span_ns << "\n";
    out << "insufficient_blocks " << (insufficient_blocks ? "yes" : "no")
        << "\n";
    out << "short_trajectory_warning " << (short_warning ? "yes" : "no")
        << "\n";
    out << "equilibration_drift_warning " << (drift_warning ? "yes" : "no")
        << "\n";
    out << "stable_dense_phase "
        << (skipped_no_dense_ < frames_.size() ? "yes" : "no") << "\n";
  }

  {
    std::ofstream out(info_output.get_llps_paper_summary_md_filename());
    out << "# TITAN LLPS Paper Summary\n\n";
    out << "- Frames used: " << frames_.size() << "\n";
    out << "- Geometry: " << geometry << "\n";
    out << "- Protein dense concentration: " << protein_dense.mean << " +/- "
        << protein_dense.sem << " nodes/nm^3\n";
    out << "- RNA dense concentration: " << rna_dense.mean << " +/- "
        << rna_dense.sem << " nodes/nm^3\n";
    out << "- Protein partition coefficient: " << protein_partition.mean
        << " +/- " << protein_partition.sem << "\n";
    out << "- RNA partition coefficient: " << rna_partition.mean << " +/- "
        << rna_partition.sem << "\n";
    out << "- Dense protein:RNA ratio: " << dense_ratio.mean << " +/- "
        << dense_ratio.sem << "\n";
    out << "- RNA enrichment: " << rna_enrichment.mean << " +/- "
        << rna_enrichment.sem << "\n";
    out << "- Dense fraction drift slope: " << slope << " per ns\n";
    if (skipped_no_dense_ >= frames_.size())
      out << "- Status: no_stable_dense_phase\n";
    else if (insufficient_blocks || drift_warning || short_warning)
      out << "- Status: convergence_warning\n";
    else
      out << "- Status: ok\n";
  }

  {
    std::ofstream out(info_output.get_llps_paper_summary_json_filename());
    out << std::fixed << std::setprecision(8);
    out << "{\n";
    out << "  \"frames_used\": " << frames_.size() << ",\n";
    out << "  \"frames_skipped_no_stable_dense_phase\": " << skipped_no_dense_
        << ",\n";
    out << "  \"geometry\": \"" << geometry << "\",\n";
    out << "  \"protein_dense_concentration_mean\": " << protein_dense.mean
        << ",\n";
    out << "  \"protein_dense_concentration_sem\": " << protein_dense.sem
        << ",\n";
    out << "  \"rna_dense_concentration_mean\": " << rna_dense.mean << ",\n";
    out << "  \"rna_dense_concentration_sem\": " << rna_dense.sem << ",\n";
    out << "  \"protein_partition_mean\": " << protein_partition.mean << ",\n";
    out << "  \"protein_partition_sem\": " << protein_partition.sem << ",\n";
    out << "  \"rna_partition_mean\": " << rna_partition.mean << ",\n";
    out << "  \"rna_partition_sem\": " << rna_partition.sem << ",\n";
    out << "  \"dense_protein_rna_ratio_mean\": " << dense_ratio.mean << ",\n";
    out << "  \"rna_enrichment_mean\": " << rna_enrichment.mean << ",\n";
    out << "  \"dense_fraction_mean\": " << dense_fraction.mean << ",\n";
    out << "  \"dense_fraction_sem\": " << dense_fraction.sem << ",\n";
    out << "  \"dense_fraction_drift_slope_per_ns\": " << slope << ",\n";
    out << "  \"insufficient_blocks\": "
        << (insufficient_blocks ? "true" : "false") << ",\n";
    out << "  \"equilibration_drift_warning\": "
        << (drift_warning ? "true" : "false") << ",\n";
    out << "  \"short_trajectory_warning\": "
        << (short_warning ? "true" : "false") << "\n";
    out << "}\n";
  }

  TitanUI::log_section(log, "LLPS Paper Output Audit");
  TitanUI::log_kv(log, info_output.get_llps_phase_partition_summary_filename(),
                  "emitted");
  TitanUI::log_kv(log, profile_file, "emitted");
  TitanUI::log_kv(log, info_output.get_llps_composition_timeseries_filename(),
                  "emitted");
  TitanUI::log_kv(log, info_output.get_llps_contact_taxonomy_filename(),
                  "emitted");
  TitanUI::log_kv(log, info_output.get_llps_contact_taxonomy_summary_filename(),
                  "emitted");
  TitanUI::log_kv(log, info_output.get_llps_convergence_report_filename(),
                  "emitted");
  TitanUI::log_kv(log, info_output.get_llps_paper_summary_md_filename(),
                  "emitted");
  TitanUI::log_kv(log, info_output.get_llps_paper_summary_json_filename(),
                  "emitted");
}

void TitanLLPSCluster::configure_paper(const LLPSPaperConfig &config) {
  paper_accumulator_.configure(config);
}

void TitanLLPSCluster::export_paper_outputs(const info &info_output,
                                            std::fstream &log) const {
  paper_accumulator_.export_outputs(info_output, log);
}
