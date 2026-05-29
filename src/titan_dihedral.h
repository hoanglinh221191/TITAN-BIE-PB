#ifndef TITAN_DIHEDRAL_H
#define TITAN_DIHEDRAL_H

#include "Multi_MOL.h"
#include "arguments.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <map>
#include <string>
#include <tuple>
#include <vector>

class TitanExecutor;
struct BoxMatrix;
class info;

class TitanDihedralAnalyzer {
public:
  struct Point {
    double x{};
    double y{};
    double z{};
  };

  TitanDihedralAnalyzer() = default;

  static double compute_dihedral_deg(const Point &p0, const Point &p1,
                                     const Point &p2, const Point &p3);

  void configure(const DihedralArgs &args);
  void init(const Multi_MOL &mol, const arguments &args, std::fstream &log);
  void process_frame(const Multi_MOL &mol, std::size_t snap_idx,
                     TitanExecutor &executor);
  void finalize_trajectory(const std::string &base_dir, const info &info_output,
                           const arguments &args,
                           const std::vector<std::string> &input_line,
                           int input_name, std::fstream &log);

  [[nodiscard]] bool enabled() const noexcept { return enabled_; }
  [[nodiscard]] bool initialized() const noexcept { return initialized_; }

private:
  struct ResidueBackbone {
    std::size_t ordinal{0};
    std::size_t resnr{0};
    std::size_t chainnum{0};
    std::size_t local_index{0};
    std::string chain{" "};
    std::string resname{"UNK"};
    int atom_n{-1};
    int atom_ca{-1};
    int atom_c{-1};
  };

  struct ResidueRange {
    bool has_chain{false};
    std::string chain{};
    std::size_t begin{0};
    std::size_t end{0};
  };

  struct AngleSample {
    std::size_t residue_index{0};
    double phi{0.0};
    double psi{0.0};
    double omega{0.0};
  };

  struct FrameRecord {
    std::size_t frame{0};
    double time_ps{0.0};
    std::vector<AngleSample> samples;
  };

  void reset_trajectory_state();
  void build_residue_map(const std::vector<mol_info> &mol_info);
  void parse_residue_selection();
  [[nodiscard]] bool residue_selected(const ResidueBackbone &residue) const;
  [[nodiscard]] bool adjacent_residues(std::size_t lhs,
                                       std::size_t rhs) const noexcept;
  [[nodiscard]] bool has_valid_box(const box &frame_box) const noexcept;
  [[nodiscard]] Point atom_point(const MOL &frame, int atom_idx) const;
  [[nodiscard]] double torsion(const MOL &frame, int atom0, int atom1,
                               int atom2, int atom3, bool use_mic,
                               const BoxMatrix *box_matrix) const;
  [[nodiscard]] std::string residue_label(const ResidueBackbone &residue) const;
  [[nodiscard]] std::string format_angle(double value) const;
  [[nodiscard]] int ramachandran_bin(double angle) const;
  void export_timeseries(const std::string &filename, std::fstream &log) const;
  void export_ramachandran(const std::string &filename,
                           std::fstream &log) const;
  void export_summary(const std::string &filename, std::fstream &log) const;
  void export_feature_matrix(const std::string &filename,
                             std::fstream &log) const;

  bool requested_{false};
  bool enabled_{false};
  bool initialized_{false};
  bool export_tica_{false};
  std::string kind_{"backbone"};
  std::string residue_selection_{"all"};
  int stride_{1};
  std::string pbc_policy_{"mic_bond_vectors_when_box_valid"};

  std::vector<ResidueBackbone> residues_;
  std::vector<ResidueRange> ranges_;
  std::vector<std::size_t> selected_indices_;
  std::vector<FrameRecord> frames_;
  std::map<std::tuple<std::size_t, int, int>, std::uint64_t> rama_counts_;
  std::vector<std::uint64_t> rama_valid_samples_;

  std::uint64_t frames_seen_{0};
  std::uint64_t frames_used_{0};
  std::array<std::uint64_t, 3> finite_counts_{};
  std::array<std::uint64_t, 3> missing_counts_{};
};

#endif // TITAN_DIHEDRAL_H
