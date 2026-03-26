// titan_llps_cluster.h — Phase Separation Morphology via DBSCAN
// Fix 1: BoxMatrix MIC (triclinic support)
// Fix 2: Cell dedup + std::floor
// Fix 3: Zero-allocation member buffers (head_, nxt_, adj_, bfs_q_)
// Fix 4: min_pts includes self + snprintf PDB + streaming CSV
// Fix 5: Hybrid-36 resSeq encoding for >9999 chains
#ifndef TITAN_LLPS_CLUSTER_H
#define TITAN_LLPS_CLUSTER_H

#include "MOL.h"
#include "MOL_info.h"
#include "box_type.h"

#include <fstream>
#include <string>
#include <vector>

struct LLPSNode {
  int id;
  int type;       // 1: Protein, 2: RNA/Ligand
  int chain_idx;
  double x, y, z;
  int cluster_id{0};   // 0: unclassified, -1: noise, >0: cluster
  bool is_core{false}; // DBSCAN core point
};

struct ClusterMetrics {
  double time_ps;
  int num_nodes;
  int num_clusters;               // condensates with ≥2 unique pro chains
  int largest_cluster_pro_chains;
  int largest_cluster_rna_chains;
  int noise_count;

  // Per-phase Radius of Gyration (Å)
  double rg_pro_dilute{0.0};
  double rg_pro_dense{0.0};
  double rg_rna_dilute{0.0};
  double rg_rna_dense{0.0};
};

class TitanLLPSCluster {
public:
  TitanLLPSCluster() = default;
  ~TitanLLPSCluster();

  // Call once before the frame loop.
  // out_csv_name: if non-empty, CSV is opened immediately for streaming writes
  //               (data survives crashes — written per-frame, not at end).
  void init(double eps_angstrom, int min_pts,
            const std::string &out_pdb_name,
            const std::string &out_csv_name = "");

  // Enable Fuzzy EM refinement after DBSCAN (Stage 2)
  void enable_fuzzy(bool v) { fuzzy_enabled_ = v; }

  // Call for every trajectory frame.
  void process_frame(const MOL &pro,
                     const std::vector<mol_info> &pro_info,
                     const MOL &lig,
                     const std::vector<mol_info> &lig_info,
                     double current_time_ps, int frame_idx);

  // Write time-series CSV after the frame loop.
  // No-op if streaming CSV was opened in init(); otherwise writes all at once.
  void export_timeseries(const std::string &out_csv) const;

private:
  double eps_{0.0};
  double eps_sq_{0.0};
  int min_pts_{2};
  int frame_count_{0};  // counts frames for periodic CSV flush
  bool fuzzy_enabled_{false}; // Stage 2: Fuzzy EM refinement
  int write_pdb_stride_{1};   // only write PDB every N frames

  // [Issue 4] Streaming CSV — opened in init(), written per-frame
  std::ofstream ts_out_;
  bool streaming_csv_{false};

  // Fallback in-memory buffer (used only when no CSV name given to init)
  std::vector<ClusterMetrics> timeseries_;

  std::ofstream vmd_out_;

  // [Fix 3] Persistent buffers reused every frame — zero new/delete overhead
  std::vector<int> head_;
  std::vector<int> nxt_;
  std::vector<std::vector<int>> adj_;
  std::vector<int> bfs_q_;  // [Issue 3] BFS queue as contiguous vector

  std::vector<LLPSNode>
  extract_nodes(const MOL &pro, const std::vector<mol_info> &pro_info,
                const MOL &lig,
                const std::vector<mol_info> &lig_info) const;

  // O(N) Cell-List DBSCAN with BoxMatrix MIC (triclinic safe)
  void run_dbscan(std::vector<LLPSNode> &nodes, const box &Box);

  // Fuzzy EM refinement (Stage 2) — classifies beads as dense/dilute,
  // outputs per-bead B-factor (membership) and occupancy (phase tag).
  void fuzzy_refine(const std::vector<LLPSNode>& nodes, const box& Box,
                    int giant_cluster_id,
                    std::vector<bool>& is_dense,
                    std::vector<double>& bfactor,
                    std::vector<double>& occupancy);

  // Streaming PDB with per-bead B-factor, occupancy, and chain ID
  void write_frame_to_pdb(const std::vector<LLPSNode> &nodes,
                          int frame_idx,
                          const std::vector<double>& bfactor,
                          const std::vector<double>& occupancy);

  // [Issue 5] Hybrid-36 encode n → 4-char string for PDB resSeq
  // Supports up to ~2.7M unique chain IDs without silent wrap.
  static std::string encode_hybrid36(int n);
};

#endif // TITAN_LLPS_CLUSTER_H
