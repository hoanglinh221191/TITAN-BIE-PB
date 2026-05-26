#ifndef TITAN_CONTACT_MAP_H
#define TITAN_CONTACT_MAP_H

#include "Multi_MOL.h"
#include "arguments.h"
#include "titan_contact_primitives.h"
#include "titan_neighbor_grid.h"

#include <cstddef>
#include <cstdint>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>

class TitanExecutor;
class info;

class TitanContactMapAnalyzer {
public:
  TitanContactMapAnalyzer() = default;

  void configure(const CMapArgs &args);
  void init(const Multi_MOL &mol1, const arguments &args, std::fstream &log);
  void process_frame(const Multi_MOL &mol1, std::size_t snap_idx,
                     TitanExecutor &executor);
  void finalize_trajectory(const std::string &base_dir, const info &info_output,
                           const arguments &args,
                           const std::vector<std::string> &input_line,
                           int input_name, std::fstream &log);

  [[nodiscard]] bool enabled() const noexcept { return enabled_; }
  [[nodiscard]] bool initialized() const noexcept { return initialized_; }

private:
  struct ContactAtom {
    std::size_t atom_index{};
    titan_contact::ResidueKey residue;
  };

  struct ContactRecord {
    titan_contact::ResiduePairKey key;
    std::uint64_t contact_frames{0};
    double min_distance_a{0.0};
    double distance_sum_a{0.0};
    std::uint64_t distance_samples{0};
  };

  void reset_trajectory_state();
  [[nodiscard]] bool pair_allowed(const titan_contact::ResidueKey &lhs,
                                  const titan_contact::ResidueKey &rhs) const;
  [[nodiscard]] bool has_valid_box(const box &frame_box) const noexcept;
  [[nodiscard]] std::string residue_label(
      const titan_contact::ResidueKey &residue) const;
  void export_summary(const std::string &filename, std::fstream &log) const;
  void export_edges(const std::string &filename, std::fstream &log) const;
  void export_matrix(const std::string &filename, std::fstream &log) const;

  bool requested_{false};
  bool enabled_{false};
  bool initialized_{false};
  bool warned_invalid_box_{false};

  titan_contact::ContactMode mode_{titan_contact::ContactMode::SidechainHeavy};
  std::string domain_{"both"};
  double cutoff_a_{4.5};
  int min_seqsep_{3};
  std::size_t atom_count_{0};

  std::uint64_t total_frames_seen_{0};
  std::uint64_t valid_frames_{0};
  std::uint64_t skipped_invalid_box_frames_{0};
  std::uint64_t skipped_incomplete_frames_{0};

  std::vector<ContactAtom> contact_atoms_;
  std::vector<titan_contact::ResidueKey> residues_;
  std::vector<unsigned char> active_atom_mask_;
  std::vector<int> atom_index_to_contact_pos_;
  std::unordered_map<titan_contact::ResiduePairKey, ContactRecord,
                     titan_contact::ResiduePairKeyHash>
      records_;
  titan_neighbor::HalfShellGrid neighbor_grid_;
};

#endif // TITAN_CONTACT_MAP_H
