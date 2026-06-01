#ifndef TITAN_MEMBRANE_H
#define TITAN_MEMBRANE_H

#include "Multi_MOL.h"
#include "aligned_allocator.h"
#include "arguments.h"

#include <array>
#include <cstdint>
#include <fstream>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

class TitanExecutor;
class info;

class TitanMembraneAnalyzer {
public:
  struct FrameSummary {
    std::size_t frame{0};
    double time_ps{0.0};
    double mean_apl_a2{0.0};
    double mean_thickness_a{0.0};
    double interdigitation{0.0};
    double penetration_depth_a{0.0};
    std::size_t upper_count{0};
    std::size_t lower_count{0};
  };

  TitanMembraneAnalyzer() = default;

  void configure(const MembraneArgs &args);
  void init(const Multi_MOL &lipids, std::fstream &log);
  void process_frame(const Multi_MOL &lipids, std::size_t snap_idx,
                     TitanExecutor &executor);
  void finalize_trajectory(const std::string &base_dir, const info &info_output,
                           const arguments &args,
                           const std::vector<std::string> &input_line,
                           int input_name, std::fstream &log);

  [[nodiscard]] bool enabled() const noexcept { return enabled_; }
  [[nodiscard]] bool initialized() const noexcept { return initialized_; }
  [[nodiscard]] const std::vector<FrameSummary> &frames() const noexcept {
    return frame_summaries_;
  }
  [[nodiscard]] std::size_t lipid_count() const noexcept {
    return lipid_plans_.size();
  }

  [[nodiscard]] static double order_parameter(double dx, double dy, double dz,
                                              double nx, double ny,
                                              double nz) noexcept;

private:
  struct LipidPlan {
    std::string id;
    std::string chain;
    std::string resname;
    std::size_t resnr{0};
    std::size_t head_atom{0};
    std::vector<std::vector<std::size_t>> tails;
  };

  struct ThicknessCell {
    double sum{0.0};
    double sum_sq{0.0};
    std::size_t count{0};
  };

  struct OrderAccumulator {
    double sum{0.0};
    double sum_sq{0.0};
    std::size_t count{0};
  };

  struct FlipState {
    int stable{0};
    int pending{0};
    int pending_count{0};
  };

  struct FlipEvent {
    std::size_t frame{0};
    double time_ps{0.0};
    std::string lipid_id;
    std::string from_leaflet;
    std::string to_leaflet;
  };

  struct AplRow {
    std::size_t frame{0};
    double time_ps{0.0};
    std::string lipid_id;
    std::string leaflet;
    double apl_a2{0.0};
  };

  struct LeafletRow {
    std::size_t frame{0};
    double time_ps{0.0};
    std::string lipid_id;
    std::string leaflet;
    double signed_distance_a{0.0};
  };

  void reset_trajectory_state();
  void build_topology_cache(const Multi_MOL &lipids, std::fstream &log);
  [[nodiscard]] std::array<double, 3>
  resolve_normal(const MOL &frame) const;
  void build_tangent_basis(const std::array<double, 3> &normal,
                           std::array<double, 3> &u,
                           std::array<double, 3> &v) const noexcept;
  void classify_leaflets(const MOL &frame, const std::array<double, 3> &normal,
                         std::vector<int> &leaflets,
                         AlignedDoubleVector &signed_distance) const;
  void compute_apl(const MOL &frame, const std::array<double, 3> &normal,
                   const std::vector<int> &leaflets, TitanExecutor &executor,
                   std::vector<double> &apl, double &surface_area) const;
  void update_order_parameters(const MOL &frame,
                               const std::array<double, 3> &normal);
  void update_flipflop(const std::vector<int> &leaflets, std::size_t frame,
                       double time_ps);
  void compute_interdigitation(const MOL &frame,
                               const std::array<double, 3> &normal,
                               const std::vector<int> &leaflets,
                               double &overlap, double &penetration) const;
  void update_thickness_map(const MOL &frame,
                            const std::array<double, 3> &normal,
                            const std::vector<int> &leaflets,
                            double &mean_thickness);
  [[nodiscard]] std::pair<int, int> resolved_grid(const MOL *frame) const;
  [[nodiscard]] std::string leaflet_name(int leaflet) const;
  void load_user_leaflets();

  void export_leaflets(const std::string &filename) const;
  void export_apl_timeseries(const std::string &filename) const;
  void export_apl_per_lipid(const std::string &filename) const;
  void export_thickness_timeseries(const std::string &filename) const;
  void export_thickness_map(const std::string &filename) const;
  void export_order_parameters(const std::string &filename) const;
  void export_interdigitation(const std::string &filename) const;
  void export_flipflop_events(const std::string &filename) const;
  void export_summary(const std::string &filename) const;
  void accumulate_trajectory();
  void export_average(const std::string &base_dir,
                      const info &info_output) const;

  MembraneArgs args_{};
  bool requested_{false};
  bool enabled_{false};
  bool initialized_{false};
  std::vector<LipidPlan> lipid_plans_;
  std::unordered_map<std::string, int> user_leaflets_;

  mutable std::array<double, 3> previous_normal_{0.0, 0.0, 1.0};
  mutable bool have_previous_normal_{false};
  AlignedDoubleVector head_x_;
  AlignedDoubleVector head_y_;
  AlignedDoubleVector head_z_;

  std::vector<FrameSummary> frame_summaries_;
  std::vector<AplRow> apl_rows_;
  std::vector<LeafletRow> leaflet_rows_;
  std::map<std::string, OrderAccumulator> order_accumulators_;
  std::vector<FlipState> flip_states_;
  std::vector<FlipEvent> flip_events_;
  std::vector<ThicknessCell> thickness_cells_;
  int grid_x_{0};
  int grid_y_{0};
  std::size_t skipped_frames_{0};
  std::size_t unknown_lipids_{0};

  ThicknessCell aggregate_apl_;
  ThicknessCell aggregate_thickness_;
  ThicknessCell aggregate_interdigitation_;
  std::vector<ThicknessCell> aggregate_thickness_cells_;
  std::map<std::string, ThicknessCell> aggregate_order_parameters_;
  std::size_t aggregate_trajectories_{0};
};

#endif // TITAN_MEMBRANE_H
