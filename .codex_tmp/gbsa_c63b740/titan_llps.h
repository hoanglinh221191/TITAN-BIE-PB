#pragma once
#ifndef TITAN_LLPS_H
#define TITAN_LLPS_H

#include "MOL.h"
#include "MOL_info.h"
#include "info.h"
#include <atomic>
#include <condition_variable>
#include <cstdint>
#include <fstream>
#include <mutex>
#include <optional>
#include <queue>
#include <string>
#include <thread>
#include <vector>

// ============================================================================
// Structs
// ============================================================================

/// Per-frame snapshot of LLPS state (macro-node level)
struct LLPSFrameResult {
  int frame_idx;
  float time_ns;
  float gcc_fraction;   // N_gcc_macronodes / N_macronodes
  int n_clusters;       // number of connected components
  int gcc_size;         // macro-nodes in giant connected component
  float gcc_rg_nm;      // Rg of GCC (nm), PBC-corrected
  float mean_degree;    // <k> = 2*edges / n_macronodes
  float density_dense;  // macro-nodes/nm³ inside GCC (sphere approx)
  float density_dilute; // macro-nodes/nm³ outside GCC
};

/// Full-trajectory LLPS results
struct LLPSResult {
  std::vector<LLPSFrameResult> frames;
  std::vector<double> ct_lags_ns;
  std::vector<double> ct_values;
  double tau_ns = 0.0;
  bool is_percolated = false;
  double percolation_time_ns = -1.0;
  std::vector<double> pk_mean; // P(k) mean over trajectory
  std::vector<double> pk_std;
  // Residue-level contact frequency (optional, sticker beads)
  std::vector<std::vector<float>> contact_freq; // [n_sticker][n_sticker]
  std::vector<int> sticker_residue_indices;
};

// ============================================================================
// TitanLLPS
// ============================================================================
class TitanLLPS {
public:
  // ---- Construction --------------------------------------------------------
  // n_chains       : number of CG chains from impl_martini chainnum assignment
  // chain_start    : chain_start[i] = first atom index of chain i
  // chain_len      : chain_len[i]   = number of atoms in chain i
  // r_cut_nm       : contact cutoff (nm), default 0.7 for Martini
  // perc_thr       : GCC fraction threshold to declare percolation
  // dimer_size     : group every N consecutive chains into 1 macro-node
  // (default 1) contacts_on    : enable residue contact frequency map
  TitanLLPS(int n_chains, const std::vector<int> &chain_start,
            const std::vector<int> &chain_len, float r_cut_nm = 0.7f,
            float perc_thr = 0.5f, int dimer_size = 1, bool contacts_on = false,
            float voxel_size_nm = 0.5f);

  // ---- Setup ---------------------------------------------------------------
  // Parse sticker bead names (comma-sep) and auto-detect backbone bead per
  // macro-node
  void set_sticker_beads(const std::vector<mol_info> &info,
                         const std::string &sticker_str);

  // ---- Hot loop ------------------------------------------------------------
  // Called AFTER Titan's PBC solver for each frame.
  // mol_frame : pro_in.m_mol[chunk_snap_idx] — provides pos.x/y/z in Å
  // (post-PBC) info      : pro_in.info — provides atom names for sticker mask
  // box_nm    : box side lengths in nm (from header or 0 if non-periodic PDB)
  void process_frame(const MOL &mol_frame, const std::vector<mol_info> &info,
                     float box_x_nm, float box_y_nm, float box_z_nm,
                     float time_ns, int frame_idx,
                     const MOL *rna_frame = nullptr,
                     const std::vector<mol_info> *rna_info = nullptr);

  // ---- Async pipeline (jthread) -------------------------------------------
  // Enable pipeline: spins up a dedicated worker thread.
  void start_async();

  // Async version: packages data and pushes to queue. Returns immediately.
  void process_frame_async(const MOL &mol_frame,
                           const std::vector<mol_info> &info, float box_x_nm,
                           float box_y_nm, float box_z_nm, float time_ns,
                           int frame_idx, const MOL *rna_frame = nullptr,
                           const std::vector<mol_info> *rna_info = nullptr);

  // Block until queue is drained (call before finalize()).
  void flush_async();

  // ---- Post-processing -----------------------------------------------------
  void finalize(int n_threads, int max_lag_frames);

  // ---- Export --------------------------------------------------------------
  void export_results(const info &info_output, std::fstream &log) const;

  // ---- Accessors -----------------------------------------------------------
  const LLPSResult &result() const { return result_; }
  bool is_percolated() const { return result_.is_percolated; }
  int n_frames_processed() const {
    return static_cast<int>(result_.frames.size());
  }
  int n_macronodes() const { return n_macro_; }

private:
  // ---- Async pipeline (jthread worker) -------------------------------------
  struct FrameWorkItem {
    std::vector<float> pos_nm;                // n_atoms*3, already in nm
    std::optional<std::vector<float>> rna_nm; // present if rna_frame
    float bx, by, bz;
    float time_ns;
    int frame_idx;
    int n_rna = 0;
  };
  std::queue<FrameWorkItem> work_queue_;
  std::mutex queue_mtx_;
  std::condition_variable queue_cv_;
  bool worker_stop_ = false;
  std::jthread worker_thread_;
  bool async_mode_ = false;
  static constexpr int MAX_QUEUE_SIZE_ = 4;
  std::condition_variable queue_full_cv_; // backpressure

  void worker_loop(std::stop_token stoken);

  // Internal implementation accepting pre-converted nm coordinates.
  void process_frame_impl(const std::vector<float> &pos_nm,
                          const float *rna_nm_data, int n_rna, float bx,
                          float by, float bz, float time_ns, int frame_idx);

  // ---- Configuration -------------------------------------------------------
  int n_chains_;   // raw chains from Martini
  int dimer_size_; // chains per macro-node
  int n_macro_;    // n_chains_ / dimer_size_
  float r_cut_;    // contact cutoff nm
  float r_cut2_;   // r_cut²
  float perc_thr_;
  bool contacts_on_;

  // ---- Chain layout --------------------------------------------------------
  std::vector<int> chain_start_; // atom offset per chain
  std::vector<int> chain_len_;   // atom count per chain
  // macro-node → list of raw chain indices
  std::vector<std::vector<int>> macro_chains_;
  // macro-node → representative backbone atom index (for Rg center)
  std::vector<int> bb_atom_;

  // ---- Sticker mask --------------------------------------------------------
  std::vector<bool> sticker_mask_;    // [n_atoms_total]
  std::vector<int> sticker_atom_idx_; // sorted list of sticker atom indices

  // ---- Frame storage for C(t) — SLIDING WINDOW ----------------------------
  int n_pairs_;              // n_macro_*(n_macro_-1)/2
  int n_words_;              // ceil(n_pairs_/64)
  int max_lag_frames_ = 500; // sliding window size
  // Ring buffer: adj_ring_[frame_idx % max_lag_frames_][word]
  std::vector<std::vector<uint64_t>> adj_ring_;
  int ring_head_ = 0;    // next write position in ring
  int ring_count_ = 0;   // current number of valid entries
  int max_ring_cap_ = 0; // actual ring capacity (memory-budgeted)
  // Online C(t) accumulators: ct_num_[lag], ct_den_[lag]
  std::vector<double> ct_num_; // numerator accumulator per lag
  std::vector<double> ct_den_; // denominator accumulator per lag
  std::vector<int> ct_count_;  // number of pairs counted per lag
  int ct_total_frames_ = 0;    // total frames seen (for ring indexing)
  // mean_adj built incrementally (pair contact count over all frames)
  std::vector<int> pair_contact_cnt_; // [n_pairs_] — count how many frames each
                                      // pair is in contact

  // ---- Degree distribution (ONLINE) ----------------------------------------
  int pk_max_k_ = -1;
  std::vector<double> pk_sum_; // [k] — sum of P(k) across frames
  std::vector<double> pk_sq_;  // [k] — sum of P(k)^2 across frames

  // ---- Contact frequency accumulator (sticker beads) ----------------------
  int n_stickers_;
  std::vector<std::vector<float>> contact_accum_; // [n_stickers][n_stickers]

  // ---- RNA detection (via mol_type) ----------------------------------------
  bool has_rna_ = false;

  // ---- Precomputed atom → macro-node map (fix #1: avoid per-frame rebuild) --
  std::vector<int> atom_macro_;

  // ---- Generation counter for seen[] (fix #2: avoid per-frame alloc+zero) ---
  std::vector<int> seen_gen_;
  int cur_gen_ = 0;

  // ---- Reusable buffers (fix #3: avoid per-frame heap alloc) ----------------
  mutable std::vector<float> shifted_pos_buf_;
  mutable std::vector<float> rna_nm_buf_;
  mutable std::vector<float> pos_nm_buf_;                 // B1
  mutable std::vector<int> cell_head_buf_;                // B2
  mutable std::vector<int> next_atom_buf_;                // B2
  mutable std::vector<std::pair<int, int>> contacts_buf_; // B3
  mutable std::vector<int> degree_buf_;                   // B4
  mutable std::vector<int> gcc_macros_buf_;               // B5
  std::vector<int> sticker_local_idx_; // B7: atom → sticker index (-1 if not)

  // ---- Atom-to-chain precomputed map ---------------------------------------
  std::vector<int> atom_chain_;

  // ---- Feature 1: 3D Density Map -------------------------------------------
  float voxel_size_nm_ = 0.5f;
  int grid_nx_ = 0, grid_ny_ = 0, grid_nz_ = 0;
  mutable std::vector<uint64_t> vox_pro_;
  mutable std::vector<uint64_t> vox_rna_;
  int density_frame_count_ = 0;
  float grid_bx_ = 0, grid_by_ = 0, grid_bz_ = 0;

  // ---- Feature 2: Polymer Scaling (RNA only) --------------------------------
  int max_rna_chain_len_ = 0;
  std::vector<int> rna_chain_starts_; // auto-detected chain boundaries
  std::vector<double> sum_R2_s_;
  std::vector<uint64_t> count_R2_s_;

  // ---- Results -------------------------------------------------------------
  LLPSResult result_;
  int percolation_frame_idx_ = -1;

  // ---- Union-Find ----------------------------------------------------------
  struct UF {
    std::vector<int> parent, rank_, size_;
    void init(int n);  // allocate once
    void reset(int n); // reset each frame (no realloc)
    int find(int x);
    void merge(int a, int b);
    int largest_component() const;
  } uf_;

  // ---- Helpers -------------------------------------------------------------
  // Build macro-node contact list via cell-list O(N).
  // pos_nm: flat float[n_atoms*3] in nm (post Å→nm conversion + PBC)
  // box: nm
  void build_contacts(const std::vector<float> &pos_nm, float bx, float by,
                      float bz,
                      std::vector<std::pair<int, int>> &macro_contacts);

  // Compute GCC Rg using macro-node backbone atoms (nm, PBC unwrapped)
  float compute_gcc_rg(const std::vector<float> &pos_nm, float bx, float by,
                       float bz, const std::vector<int> &gcc_macros) const;

  static float pbc_dist2(float dx, float dy, float dz, float bx, float by,
                         float bz);

  inline int pair_index(int i, int j) const;
  inline void set_adj_bit(std::vector<uint64_t> &words, int i, int j) const;

  double fit_tau(const std::vector<double> &ct,
                 const std::vector<double> &lags_ns) const;
};

#endif // TITAN_LLPS_H
