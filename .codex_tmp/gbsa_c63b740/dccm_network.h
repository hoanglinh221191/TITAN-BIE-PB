#ifndef DCCM_NETWORK_H
#define DCCM_NETWORK_H

#include "Multi_MOL.h"
#include "arguments.h"
#include <Eigen/Dense>
#include <array>
#include <string>
#include <tuple>
#include <vector>

class DCCMNetwork {
public:
  DCCMNetwork() = default;
  ~DCCMNetwork() = default;

  // Khởi tạo các Node (CA atoms)
  void init(const Multi_MOL &pro_in, const Multi_MOL &lig,
            const arguments &args);

  // Lưu trữ tọa độ cho mỗi Frame
  void process_frame(const Multi_MOL &pro_in, const Multi_MOL &lig,
                     size_t chunk_snap_idx);

  // Tính toán siêu tốc và xuất file
  void calculate_and_export(size_t num_threads);

  // Information-Theory Allosteric Analysis (GC + Transfer Entropy)
  void calculate_allostery_information(size_t num_threads, int tau_lag,
                                       int num_bins);

  // ── LSP-MD: Local Spatial Pattern (Kornev, JCTC 2025) ──────────────────────
  struct LSPMDConfig {
    double cutoff_ang    = 12.0;  // Cα–Cα cutoff (Ångström)
    bool   use_avg_cutoff = true; // true: avg distance cutoff; false: frame 0
  };

  struct LSPMDResult {
    std::vector<double> dc;  // [num_ca] Degree Centrality
    std::vector<double> bc;  // [num_ca] Betweenness Centrality
    std::vector<std::tuple<int,int,double>> edges; // (i, j, W_ij)
  };

  void calculate_lspmd(size_t num_threads, const LSPMDConfig& cfg,
                       LSPMDResult& result);
  void export_lspmd(const LSPMDResult& result, double weight_threshold = 0.0);

  // (Tùy chọn) Cho phép PCA lấy data
  const std::vector<std::string> &get_resnames() const { return ca_resname; }
  const std::vector<int> &get_resnrs() const { return ca_resnrs; }
  const std::vector<std::string> &get_chains() const { return ca_chain; }
  const std::vector<int> &get_is_lig() const { return ca_is_lig; }
  const std::vector<int> &get_ca_indices() const { return ca_indices; }
  // Định dạng lại cho PCA nếu cần (giữ tương thích ngược)
  std::vector<std::vector<std::array<double, 3>>> get_ca_traj_aos() const;

  // Direct SoA access for TICA (zero-copy per-frame)
  size_t get_total_frames() const { return total_frames; }
  const std::vector<std::vector<double>> &get_traj_x() const { return traj_x; }
  const std::vector<std::vector<double>> &get_traj_y() const { return traj_y; }
  const std::vector<std::vector<double>> &get_traj_z() const { return traj_z; }

  // Cβ (or N for Gly) SoA access
  const std::vector<int> &get_cb_indices() const { return cb_indices; }
  const std::vector<std::vector<double>> &get_traj_cbx() const { return traj_cbx; }
  const std::vector<std::vector<double>> &get_traj_cby() const { return traj_cby; }
  const std::vector<std::vector<double>> &get_traj_cbz() const { return traj_cbz; }

private:
  std::vector<int> ca_indices;
  std::vector<int> ca_resnrs;
  std::vector<int> ca_is_lig;
  std::vector<std::string> ca_chain;
  std::vector<std::string> ca_resname;

  // Cβ (or N for Gly) global atom index per CA node (-1 = not found, fallback to CA)
  std::vector<int> cb_indices;

  // Lưu trữ theo trục (Structure of Arrays - SoA) để tối ưu hóa bộ nhớ Cache
  // cho AVX2
  std::vector<std::vector<double>> traj_x; // [frame][atom] — CA
  std::vector<std::vector<double>> traj_y;
  std::vector<std::vector<double>> traj_z;

  // Cβ/N trajectory (same layout as CA)
  std::vector<std::vector<double>> traj_cbx; // [frame][atom]
  std::vector<std::vector<double>> traj_cby;
  std::vector<std::vector<double>> traj_cbz;

  size_t total_frames = 0;
  bool kabsch_aligned_ = false; // guard against double Kabsch alignment

  // Helper: find Cβ (or N for Gly) atom index within a residue
  static int find_cb_for_residue(const Multi_MOL &mol, int resnr,
                                 const std::string &chain, bool is_gly);
};

#endif // DCCM_NETWORK_H
