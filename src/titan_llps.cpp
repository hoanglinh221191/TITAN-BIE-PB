// titan_llps.cpp - LLPS Analysis (Martini CG, macro-node level)
#include "titan_llps.h"
#include "MOL.h"
#include "MOL_info.h"
#include "info.h"
#include "memory_utils.h"
#include "titan_executor.h"
#include "titan_gyration.h"
#include "titan_simd_dispatch.h"
#include "titan_ui.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <sstream>
#include <thread>

// ============================================================================
// Union-Find
// ============================================================================





// ============================================================================
// Helpers
// ============================================================================



// ============================================================================
// Constructor
// ============================================================================
TitanLLPS::TitanLLPS(int n_chains, const std::vector<int> &chain_start,
                     const std::vector<int> &chain_len, float r_cut_nm,
                     float perc_thr, int dimer_size, bool contacts_on,
                     float voxel_size_nm, bool dynamics_enabled,
                     bool dynamics_detailed, float prna_cutoff_nm,
                     const std::string &kernel_mode,
                     const std::string &density_mode)
    : n_chains_(n_chains), dimer_size_(std::max(1, dimer_size)),
      n_macro_(0), r_cut_(r_cut_nm), r_cut2_(r_cut_nm * r_cut_nm),
      perc_thr_(perc_thr), contacts_on_(contacts_on),
      dynamics_enabled_(dynamics_enabled),
      dynamics_detailed_(dynamics_detailed),
      prna_cutoff_nm_(prna_cutoff_nm > 0.0f ? prna_cutoff_nm : r_cut_nm),
      prna_cutoff2_nm_(prna_cutoff_nm_ * prna_cutoff_nm_),
      kernel_mode_(parse_kernel_mode(kernel_mode)), density_mode_(density_mode),
      chain_start_(chain_start),
      chain_len_(chain_len), n_stickers_(0), voxel_size_nm_(voxel_size_nm) {
  n_macro_ = (n_chains_ + dimer_size_ - 1) / dimer_size_;
  macro_chains_.resize(n_macro_);
  for (int c = 0; c < n_chains_; ++c)
    macro_chains_[c / dimer_size_].push_back(c);

  bb_atom_.assign(n_macro_, -1);
  n_pairs_ = n_macro_ * (n_macro_ - 1) / 2;
  n_words_ = (n_pairs_ + 63) / 64;
  uf_.init(n_macro_);

  int total_atoms = 0;
  for (int c = 0; c < n_chains_; ++c) {
    if (chain_start_[c] >= 0)
      total_atoms = std::max(total_atoms, chain_start_[c] + chain_len_[c]);
  }
  atom_chain_.assign(total_atoms, -1);
  atom_residue_.assign(total_atoms, -1);
  for (int c = 0; c < n_chains_; ++c) {
    if (chain_start_[c] < 0)
      continue;
    for (int a = chain_start_[c]; a < chain_start_[c] + chain_len_[c]; ++a)
      atom_chain_[a] = c;
  }

  atom_macro_.assign(total_atoms, -1);
  for (int m = 0; m < n_macro_; ++m) {
    for (int c : macro_chains_[m]) {
      if (chain_start_[c] < 0)
        continue;
      for (int a = chain_start_[c]; a < chain_start_[c] + chain_len_[c]; ++a)
        atom_macro_[a] = m;
    }
  }
  seen_gen_.assign(static_cast<std::size_t>(n_macro_) * n_macro_, 0);
  prev_neighbor_offsets_.assign(static_cast<std::size_t>(n_chains_) + 1, 0);
}

// ============================================================================
// Setup
// ============================================================================
void TitanLLPS::set_sticker_beads(const std::vector<mol_info> &info,
                                  const std::string &sticker_str) {
  std::vector<std::string> sticker_names;
  std::istringstream ss(sticker_str);
  std::string tok;
  while (std::getline(ss, tok, ',')) {
    tok.erase(0, tok.find_first_not_of(" \t"));
    tok.erase(tok.find_last_not_of(" \t") + 1);
    if (!tok.empty())
      sticker_names.push_back(tok);
  }

  sticker_mask_.assign(info.size(), false);
  sticker_local_idx_.assign(info.size(), -1);
  for (std::size_t i = 0; i < info.size(); ++i) {
    if (i < atom_residue_.size())
      atom_residue_[i] = static_cast<int>(info[i].get_resnr());
    for (const auto &name : sticker_names) {
      if (info[i].get_atom() == name) {
        sticker_mask_[i] = true;
        break;
      }
    }
  }

  static const std::vector<std::string> bb_names = {"BB", "BB1", "BB2", "CA"};
  for (int m = 0; m < n_macro_; ++m) {
    for (int c : macro_chains_[m]) {
      for (int a = chain_start_[c]; a < chain_start_[c] + chain_len_[c]; ++a) {
        if (a >= static_cast<int>(info.size()))
          break;
        for (const auto &bb : bb_names) {
          if (info[a].get_atom() == bb) {
            bb_atom_[m] = a;
            goto next_macro;
          }
        }
      }
    }
    if (!macro_chains_[m].empty() && chain_len_[macro_chains_[m][0]] > 0)
      bb_atom_[m] = chain_start_[macro_chains_[m][0]];
  next_macro:;
  }

  sticker_atom_idx_.clear();
  result_.sticker_residue_indices.clear();
  std::unordered_map<int, int> residue_to_local;
  residue_to_local.reserve(info.size());
  for (std::size_t i = 0; i < sticker_mask_.size(); ++i) {
    if (!sticker_mask_[i])
      continue;
    sticker_atom_idx_.push_back(static_cast<int>(i));
    const int resnr = static_cast<int>(info[i].get_resnr());
    auto [it, inserted] =
        residue_to_local.emplace(resnr,
                                 static_cast<int>(result_.sticker_residue_indices.size()));
    if (inserted)
      result_.sticker_residue_indices.push_back(resnr);
    sticker_local_idx_[i] = it->second;
  }
  n_stickers_ = static_cast<int>(result_.sticker_residue_indices.size());

  if (contacts_on_ && n_stickers_ > 0) {
    const std::size_t upper_pairs =
        static_cast<std::size_t>(n_stickers_) * (n_stickers_ - 1) / 2;
    contact_accum_upper_.assign(upper_pairs, 0.0f);
    sticker_seen_gen_.assign(upper_pairs, 0);
    sticker_cur_gen_ = 0;
  } else {
    contact_accum_upper_.clear();
    sticker_seen_gen_.clear();
    sticker_cur_gen_ = 0;
  }
}

void TitanLLPS::set_rna_metadata(const std::vector<mol_info> &rna_info) {
  rna_atom_chain_.assign(rna_info.size(), -1);
  rna_atom_residue_.assign(rna_info.size(), -1);
  rna_chain_raw_ids_.clear();
  rna_chain_raw_to_local_.clear();
  for (std::size_t i = 0; i < rna_info.size(); ++i) {
    const int raw_chain = static_cast<int>(rna_info[i].get_chainnum());
    auto it = rna_chain_raw_to_local_.find(raw_chain);
    int local = -1;
    if (it == rna_chain_raw_to_local_.end()) {
      local = static_cast<int>(rna_chain_raw_ids_.size());
      rna_chain_raw_to_local_.emplace(raw_chain, local);
      rna_chain_raw_ids_.push_back(raw_chain);
    } else {
      local = it->second;
    }
    rna_atom_chain_[i] = local;
    rna_atom_residue_[i] = static_cast<int>(rna_info[i].get_resnr());
  }
  prev_neighbor_offsets_.assign(static_cast<std::size_t>(n_chains_) +
                                    rna_chain_raw_ids_.size() + 1,
                                0);
  prev_neighbor_flat_.clear();
}

void TitanLLPS::configure_dynamics(bool enabled, bool detailed,
                                   float prna_cutoff_nm) {
  dynamics_enabled_ = enabled;
  dynamics_detailed_ = detailed;
  prna_cutoff_nm_ = (prna_cutoff_nm > 0.0f) ? prna_cutoff_nm : r_cut_;
  prna_cutoff2_nm_ = prna_cutoff_nm_ * prna_cutoff_nm_;
}

void TitanLLPS::configure_trace_outputs(const info &info_output) {
  contact_trace_streaming_ = false;
  if (contact_trace_out_.is_open())
    contact_trace_out_.close();
  if (!dynamics_enabled_ || !dynamics_detailed_)
    return;

  contact_trace_out_.open(info_output.get_llps_chain_contact_trace_filename());
  if (!contact_trace_out_.is_open()) {
    TitanUI::print_warning("[LLPS] Cannot open " +
                           info_output.get_llps_chain_contact_trace_filename());
    return;
  }
  contact_trace_out_ << "Time_ns,ContactType,IdA,IdB,State,Duration_ns\n";
  contact_trace_out_ << std::fixed << std::setprecision(6);
  contact_trace_streaming_ = true;
  contact_trace_rows_.clear();
}

// ============================================================================
// Distance helpers
// ============================================================================





LlpsKernelMode TitanLLPS::parse_kernel_mode(const std::string &mode) {
  std::string lowered = mode;
  std::transform(lowered.begin(), lowered.end(), lowered.begin(),
                 [](unsigned char c) {
                   return static_cast<char>(std::tolower(c));
                 });
  if (lowered == "scalar")
    return LlpsKernelMode::Scalar;
  if (lowered == "avx2")
    return LlpsKernelMode::Avx2;
  return LlpsKernelMode::Auto;
}

void TitanLLPS::expand_contact_map_result(int n_frames) {
  result_.contact_freq.assign(
      static_cast<std::size_t>(n_stickers_),
      std::vector<float>(static_cast<std::size_t>(n_stickers_), 0.0f));
  if (n_stickers_ <= 1 || n_frames <= 0 || contact_accum_upper_.empty())
    return;
  const float inv_frames = 1.0f / static_cast<float>(n_frames);
  for (int i = 0; i < n_stickers_; ++i) {
    for (int j = i + 1; j < n_stickers_; ++j) {
      const float value =
          contact_accum_upper_[contact_upper_index(i, j)] * inv_frames;
      result_.contact_freq[static_cast<std::size_t>(i)]
                          [static_cast<std::size_t>(j)] = value;
      result_.contact_freq[static_cast<std::size_t>(j)]
                          [static_cast<std::size_t>(i)] = value;
    }
  }
}

// ============================================================================
// Contact builders
// ============================================================================


// ============================================================================
// GCC Rg
// ============================================================================

// ============================================================================
// Process frame
// ============================================================================
void TitanLLPS::process_frame(const MOL &mol_frame,
                              const std::vector<mol_info> & /*info*/,
                              float box_x_nm, float box_y_nm, float box_z_nm,
                              float time_ns, int frame_idx,
                              const MOL *rna_frame,
                              const std::vector<mol_info> *rna_info) {
  const auto pack_begin = std::chrono::steady_clock::now();
  const auto &ax = mol_frame.get_pos_x();
  const auto &ay = mol_frame.get_pos_y();
  const auto &az = mol_frame.get_pos_z();
  const int n_atoms = static_cast<int>(ax.size());

  constexpr float kAngToNm = 0.1f;
  pos_nm_buf_.resize(static_cast<std::size_t>(n_atoms) * 3);
  for (int i = 0; i < n_atoms; ++i) {
    pos_nm_buf_[3 * i + 0] = static_cast<float>(ax[i]) * kAngToNm;
    pos_nm_buf_[3 * i + 1] = static_cast<float>(ay[i]) * kAngToNm;
    pos_nm_buf_[3 * i + 2] = static_cast<float>(az[i]) * kAngToNm;
  }

  const float *rna_data = nullptr;
  int n_rna = 0;
  if (rna_frame) {
    if (rna_info && rna_atom_chain_.empty())
      set_rna_metadata(*rna_info);
    const auto &rx = rna_frame->get_pos_x();
    const auto &ry = rna_frame->get_pos_y();
    const auto &rz = rna_frame->get_pos_z();
    n_rna = static_cast<int>(rx.size());
    rna_nm_buf_.resize(static_cast<std::size_t>(n_rna) * 3);
    for (int i = 0; i < n_rna; ++i) {
      rna_nm_buf_[3 * i + 0] = static_cast<float>(rx[i]) * kAngToNm;
      rna_nm_buf_[3 * i + 1] = static_cast<float>(ry[i]) * kAngToNm;
      rna_nm_buf_[3 * i + 2] = static_cast<float>(rz[i]) * kAngToNm;
    }
    rna_data = rna_nm_buf_.data();
  }
  process_frame_impl(pos_nm_buf_, rna_data, n_rna, box_x_nm, box_y_nm, box_z_nm,
                     time_ns, frame_idx);
  perf_.pack_seconds += std::chrono::duration<double>(
                            std::chrono::steady_clock::now() - pack_begin)
                            .count();
}

void TitanLLPS::process_frame_impl(const std::vector<float> &pos_nm,
                                   const float *rna_nm_data, int n_rna, float bx,
                                   float by, float bz, float time_ns,
                                   int frame_idx) {
  const auto frame_begin = std::chrono::steady_clock::now();
  const int n_atoms = static_cast<int>(pos_nm.size()) / 3;

  macro_contact_keys_buf_.clear();
  auto &contacts = macro_contact_keys_buf_;
  std::vector<uint64_t> *sticker_pairs =
      dynamics_enabled_ ? &sticker_pair_keys_buf_ : nullptr;
  std::vector<uint64_t> *protein_pairs =
      dynamics_enabled_ ? &protein_chain_pair_keys_buf_ : nullptr;
  const auto contact_begin = std::chrono::steady_clock::now();
  build_contacts(pos_nm, bx, by, bz, contacts, sticker_pairs, protein_pairs);
  perf_.contact_seconds += std::chrono::duration<double>(
                               std::chrono::steady_clock::now() - contact_begin)
                               .count();
  perf_.macro_contact_pairs += static_cast<uint64_t>(contacts.size());
  if (sticker_pairs)
    perf_.sticker_contact_pairs += static_cast<uint64_t>(sticker_pairs->size());
  if (protein_pairs)
    perf_.protein_chain_pairs += static_cast<uint64_t>(protein_pairs->size());

  if (dynamics_enabled_) {
    const auto prna_begin = std::chrono::steady_clock::now();
    build_protein_rna_contacts(pos_nm, rna_nm_data, n_rna, bx, by, bz,
                               protein_rna_pair_keys_buf_);
    perf_.protein_rna_seconds +=
        std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                      prna_begin)
            .count();
    perf_.protein_rna_pairs +=
        static_cast<uint64_t>(protein_rna_pair_keys_buf_.size());
    const auto dyn_begin = std::chrono::steady_clock::now();
    update_contact_dynamics(time_ns, frame_idx, sticker_pair_keys_buf_,
                            protein_chain_pair_keys_buf_,
                            protein_rna_pair_keys_buf_);
    perf_.dynamics_seconds += std::chrono::duration<double>(
                                  std::chrono::steady_clock::now() - dyn_begin)
                                  .count();
  }

  if (adj_ring_.empty()) {
    const std::size_t avail = get_available_ram();
    const std::size_t ring_budget =
        (avail > 0) ? (avail / 4) : (256ULL * 1024ULL * 1024ULL);
    const std::size_t bytes_per_slot =
        static_cast<std::size_t>(n_words_) * sizeof(uint64_t);
    const int ideal_cap = max_lag_frames_ * 2;
    int ring_cap_local = ideal_cap;
    if (bytes_per_slot > 0) {
      int budget_cap = static_cast<int>(ring_budget / bytes_per_slot);
      budget_cap = std::max(budget_cap, max_lag_frames_ + 1);
      ring_cap_local = std::min(ring_cap_local, budget_cap);
    }
    ring_cap_local = std::max(ring_cap_local, 2);
    max_ring_cap_ = ring_cap_local;
    adj_ring_.resize(static_cast<std::size_t>(max_ring_cap_),
                     std::vector<uint64_t>(static_cast<std::size_t>(n_words_),
                                           0));
    ct_num_.assign(static_cast<std::size_t>(max_lag_frames_), 0.0);
    ct_den_.assign(static_cast<std::size_t>(max_lag_frames_), 0.0);
    ct_count_.assign(static_cast<std::size_t>(max_lag_frames_), 0);
    pair_contact_cnt_.assign(static_cast<std::size_t>(n_pairs_), 0);
  }

  const int slot = ct_total_frames_ % max_ring_cap_;
  auto &bits = adj_ring_[static_cast<std::size_t>(slot)];
  std::fill(bits.begin(), bits.end(), uint64_t(0));
  for (uint64_t key : contacts) {
    const auto [ci, cj] = unpack_pair(key);
    set_adj_bit(bits, ci, cj);
  }

  for (int word = 0; word < n_words_; ++word) {
    uint64_t mask = bits[static_cast<std::size_t>(word)];
    while (mask) {
      const int bit = __builtin_ctzll(mask);
      const int idx = word * 64 + bit;
      if (idx < n_pairs_)
        ++pair_contact_cnt_[static_cast<std::size_t>(idx)];
      mask &= mask - 1;
    }
  }
  ring_count_ = std::min(ring_count_ + 1, max_ring_cap_);
  ++ct_total_frames_;

  uf_.reset(n_macro_);
  degree_buf_.assign(static_cast<std::size_t>(n_macro_), 0);
  auto &degree = degree_buf_;
  const auto graph_begin = std::chrono::steady_clock::now();
  for (uint64_t key : contacts) {
    const auto [ci, cj] = unpack_pair(key);
    uf_.merge(ci, cj);
    ++degree[static_cast<std::size_t>(ci)];
    ++degree[static_cast<std::size_t>(cj)];
  }
  perf_.graph_seconds += std::chrono::duration<double>(
                             std::chrono::steady_clock::now() - graph_begin)
                             .count();

  int gcc_size = uf_.largest_component();
  int gcc_root = -1;
  for (int macro = 0; macro < n_macro_; ++macro) {
    if (uf_.find(macro) == macro && uf_.size_[macro] == gcc_size) {
      gcc_root = macro;
      break;
    }
  }

  gcc_macros_buf_.clear();
  auto &gcc_macros = gcc_macros_buf_;
  int n_clusters = 0;
  for (int macro = 0; macro < n_macro_; ++macro) {
    const int root = uf_.find(macro);
    if (root == macro)
      ++n_clusters;
    if (gcc_root >= 0 && root == gcc_root)
      gcc_macros.push_back(macro);
  }

  const float gcc_fraction =
      (n_macro_ > 0) ? static_cast<float>(gcc_size) / n_macro_ : 0.0f;
  const float gcc_rg = compute_gcc_rg(pos_nm, bx, by, bz, gcc_macros);
  const float v_box = (bx > 0 && by > 0 && bz > 0) ? bx * by * bz : 1.0f;
  constexpr float kPi = 3.14159265f;
  const float v_gcc = (4.0f / 3.0f) * kPi * gcc_rg * gcc_rg * gcc_rg;
  const float v_dil = v_box - v_gcc;
  const float rho_dense = (v_gcc > 1e-9f) ? gcc_size / v_gcc : 0.0f;
  const float rho_dilute =
      (v_dil > 1e-9f) ? (n_macro_ - gcc_size) / v_dil : 0.0f;

  float mean_k = 0.0f;
  for (int d : degree)
    mean_k += static_cast<float>(d);
  if (n_macro_ > 0)
    mean_k /= static_cast<float>(n_macro_);

  result_.frames.push_back({frame_idx, time_ns, gcc_fraction, n_clusters,
                            gcc_size, gcc_rg, mean_k, rho_dense, rho_dilute});

  if (!result_.is_percolated && gcc_fraction >= perc_thr_) {
    result_.is_percolated = true;
    result_.percolation_time_ns = time_ns;
    percolation_frame_idx_ = static_cast<int>(result_.frames.size()) - 1;
  }

  int local_max_k = 0;
  for (int d : degree)
    local_max_k = std::max(local_max_k, d);
  if (local_max_k > pk_max_k_) {
    pk_sum_.resize(static_cast<std::size_t>(local_max_k + 1), 0.0);
    pk_sq_.resize(static_cast<std::size_t>(local_max_k + 1), 0.0);
    pk_max_k_ = local_max_k;
  }
  for (int d : degree) {
    const double pk = (n_macro_ > 0) ? (1.0 / n_macro_) : 0.0;
    pk_sum_[static_cast<std::size_t>(d)] += pk;
    pk_sq_[static_cast<std::size_t>(d)] += pk * pk;
  }

  if (density_mode_ != "off" && bx > 0 && by > 0 && bz > 0 &&
      !gcc_macros.empty()) {
    const auto density_begin = std::chrono::steady_clock::now();
    double com_x = 0.0;
    double com_y = 0.0;
    double com_z = 0.0;
    int gcc_atom_count = 0;
    for (int macro : gcc_macros) {
      for (int chain : macro_chains_[macro]) {
        for (int atom = chain_start_[chain];
             atom < chain_start_[chain] + chain_len_[chain]; ++atom) {
          com_x += pos_nm[3 * atom + 0];
          com_y += pos_nm[3 * atom + 1];
          com_z += pos_nm[3 * atom + 2];
          ++gcc_atom_count;
        }
      }
    }
    if (gcc_atom_count > 0) {
      com_x /= gcc_atom_count;
      com_y /= gcc_atom_count;
      com_z /= gcc_atom_count;
    }
    const float sx = bx * 0.5f - static_cast<float>(com_x);
    const float sy = by * 0.5f - static_cast<float>(com_y);
    const float sz = bz * 0.5f - static_cast<float>(com_z);
    shifted_pos_buf_.resize(pos_nm.size());
    auto &shifted_pos = shifted_pos_buf_;
    for (int atom = 0; atom < n_atoms; ++atom) {
      float x = pos_nm[3 * atom + 0] + sx;
      float y = pos_nm[3 * atom + 1] + sy;
      float z = pos_nm[3 * atom + 2] + sz;
      x -= bx * std::floor(x / bx);
      y -= by * std::floor(y / by);
      z -= bz * std::floor(z / bz);
      shifted_pos[3 * atom + 0] = x;
      shifted_pos[3 * atom + 1] = y;
      shifted_pos[3 * atom + 2] = z;
    }

    const int nx = std::max(1, static_cast<int>(bx / voxel_size_nm_));
    const int ny = std::max(1, static_cast<int>(by / voxel_size_nm_));
    const int nz = std::max(1, static_cast<int>(bz / voxel_size_nm_));
    const std::size_t needed = static_cast<std::size_t>(nx) * ny * nz;
    const std::size_t density_bytes = needed * 2 * sizeof(uint64_t);
    if (density_mode_ == "projection" && !density_disabled_by_budget_) {
      density_disabled_by_budget_ = true;
      TitanUI::print_warning(
          "LLPS density: -llps-density=projection disables full 3D DX density "
          "accumulation to reduce RAM.");
    } else if (density_mode_ == "auto" && !density_disabled_by_budget_) {
      const size_t available = get_available_ram();
      if (available > 0 &&
          static_cast<long double>(density_bytes) >
              static_cast<long double>(available) * 0.20L) {
        density_disabled_by_budget_ = true;
        TitanUI::print_warning(
            "LLPS density: auto policy skipped full 3D DX accumulation because "
            "the voxel grid exceeds 20% of currently available RAM. Use "
            "-llps-density=full to force it.");
      }
    }

    if (!density_disabled_by_budget_) {
    if (needed > vox_pro_.size()) {
      vox_pro_.assign(needed, 0);
      vox_rna_.assign(needed, 0);
      grid_nx_ = nx;
      grid_ny_ = ny;
      grid_nz_ = nz;
      grid_bx_ = bx;
      grid_by_ = by;
      grid_bz_ = bz;
      density_frame_count_ = 0;
    }

    for (int atom = 0; atom < n_atoms; ++atom) {
      int ix = static_cast<int>(shifted_pos[3 * atom + 0] / voxel_size_nm_);
      int iy = static_cast<int>(shifted_pos[3 * atom + 1] / voxel_size_nm_);
      int iz = static_cast<int>(shifted_pos[3 * atom + 2] / voxel_size_nm_);
      ix = std::clamp(ix, 0, grid_nx_ - 1);
      iy = std::clamp(iy, 0, grid_ny_ - 1);
      iz = std::clamp(iz, 0, grid_nz_ - 1);
      const std::size_t flat =
          static_cast<std::size_t>(ix) * grid_ny_ * grid_nz_ +
          static_cast<std::size_t>(iy) * grid_nz_ + static_cast<std::size_t>(iz);
      ++vox_pro_[flat];
    }
    if (rna_nm_data && n_rna > 0) {
      for (int atom = 0; atom < n_rna; ++atom) {
        float x = rna_nm_data[3 * atom + 0] + sx;
        float y = rna_nm_data[3 * atom + 1] + sy;
        float z = rna_nm_data[3 * atom + 2] + sz;
        x -= bx * std::floor(x / bx);
        y -= by * std::floor(y / by);
        z -= bz * std::floor(z / bz);
        const int ix =
            std::clamp(static_cast<int>(x / voxel_size_nm_), 0, grid_nx_ - 1);
        const int iy =
            std::clamp(static_cast<int>(y / voxel_size_nm_), 0, grid_ny_ - 1);
        const int iz =
            std::clamp(static_cast<int>(z / voxel_size_nm_), 0, grid_nz_ - 1);
        const std::size_t flat =
            static_cast<std::size_t>(ix) * grid_ny_ * grid_nz_ +
            static_cast<std::size_t>(iy) * grid_nz_ +
            static_cast<std::size_t>(iz);
        ++vox_rna_[flat];
      }
      has_rna_ = true;
    }
    ++density_frame_count_;
    }
    perf_.density_seconds += std::chrono::duration<double>(
                                 std::chrono::steady_clock::now() -
                                 density_begin)
                                 .count();
  }

  if (rna_nm_data && n_rna > 0 && bx > 0 && by > 0 && bz > 0) {
    const auto polymer_begin = std::chrono::steady_clock::now();
    const int N = n_rna;
    if (rna_chain_starts_.empty()) {
      rna_chain_starts_.push_back(0);
      const float half_bx = bx * 0.5f;
      const float half_by = by * 0.5f;
      const float half_bz = bz * 0.5f;
      for (int k = 0; k < N - 1; ++k) {
        float dx = rna_nm_data[3 * k + 0] - rna_nm_data[3 * (k + 1) + 0];
        float dy = rna_nm_data[3 * k + 1] - rna_nm_data[3 * (k + 1) + 1];
        float dz = rna_nm_data[3 * k + 2] - rna_nm_data[3 * (k + 1) + 2];
        if (dx > half_bx)
          dx -= bx;
        else if (dx < -half_bx)
          dx += bx;
        if (dy > half_by)
          dy -= by;
        else if (dy < -half_by)
          dy += by;
        if (dz > half_bz)
          dz -= bz;
        else if (dz < -half_bz)
          dz += bz;
        if (dx * dx + dy * dy + dz * dz > 2.25f)
          rna_chain_starts_.push_back(k + 1);
      }
      rna_chain_starts_.push_back(N);
    }

    int local_max_len = 0;
    for (std::size_t chain = 0; chain + 1 < rna_chain_starts_.size(); ++chain) {
      local_max_len =
          std::max(local_max_len,
                   rna_chain_starts_[chain + 1] - rna_chain_starts_[chain]);
    }
    if (max_rna_chain_len_ == 0 && local_max_len > 1) {
      max_rna_chain_len_ = local_max_len;
      sum_R2_s_.assign(static_cast<std::size_t>(max_rna_chain_len_), 0.0);
      count_R2_s_.assign(static_cast<std::size_t>(max_rna_chain_len_), 0);
    }

    for (std::size_t chain = 0; chain + 1 < rna_chain_starts_.size(); ++chain) {
      const int start = rna_chain_starts_[chain];
      const int end = rna_chain_starts_[chain + 1];
      const int len = end - start;
      for (int s = 1; s < len && s < max_rna_chain_len_; ++s) {
        double r2_sum = 0.0;
        for (int k = start; k < end - s; ++k) {
          const float dx = rna_nm_data[3 * k + 0] - rna_nm_data[3 * (k + s) + 0];
          const float dy = rna_nm_data[3 * k + 1] - rna_nm_data[3 * (k + s) + 1];
          const float dz = rna_nm_data[3 * k + 2] - rna_nm_data[3 * (k + s) + 2];
          r2_sum += static_cast<double>(dx * dx + dy * dy + dz * dz);
        }
        sum_R2_s_[static_cast<std::size_t>(s)] += r2_sum;
        count_R2_s_[static_cast<std::size_t>(s)] +=
            static_cast<uint64_t>(len - s);
      }
    }
    perf_.polymer_scaling_seconds +=
        std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                      polymer_begin)
            .count();
  }
  ++perf_.frames_processed;
  (void)frame_begin;
}

// ============================================================================
// Async pipeline
// ============================================================================




// ============================================================================
// Finalize
// ============================================================================
double TitanLLPS::fit_tau(const std::vector<double> &ct,
                          const std::vector<double> &lags) const {
  if (ct.size() < 2 || ct.size() != lags.size())
    return 0.0;
  constexpr double inv_e = 1.0 / 2.71828182845904523536;
  for (std::size_t i = 1; i < ct.size(); ++i) {
    if (ct[i] <= inv_e && ct[i - 1] > ct[i]) {
      const double frac = (ct[i - 1] - inv_e) / (ct[i - 1] - ct[i]);
      return lags[i - 1] + frac * (lags[i] - lags[i - 1]);
    }
  }
  double tau = 0.0;
  for (std::size_t i = 1; i < ct.size(); ++i)
    tau += 0.5 * (ct[i - 1] + ct[i]) * (lags[i] - lags[i - 1]);
  return tau;
}

void TitanLLPS::finalize(TitanExecutor &executor, int max_lag_frames) {
  const auto finalize_begin = std::chrono::steady_clock::now();
  const int N = static_cast<int>(result_.frames.size());
  if (N == 0)
    return;

  result_.pk_mean.resize(static_cast<std::size_t>(pk_max_k_ + 1), 0.0);
  result_.pk_std.resize(static_cast<std::size_t>(pk_max_k_ + 1), 0.0);
  for (int k = 0; k <= pk_max_k_; ++k) {
    result_.pk_mean[static_cast<std::size_t>(k)] = pk_sum_[static_cast<std::size_t>(k)] / N;
    result_.pk_std[static_cast<std::size_t>(k)] = std::sqrt(std::max(
        0.0, pk_sq_[static_cast<std::size_t>(k)] / N -
                 result_.pk_mean[static_cast<std::size_t>(k)] *
                     result_.pk_mean[static_cast<std::size_t>(k)]));
  }

  const int perc_frame = (percolation_frame_idx_ >= 0) ? percolation_frame_idx_ : 0;
  const int post_perc_frames = N - perc_frame;
  if (post_perc_frames >= 2 && n_words_ > 0 && ring_count_ >= 2) {
    std::vector<uint64_t> mean_adj(static_cast<std::size_t>(n_words_), 0);
    for (int idx = 0; idx < n_pairs_; ++idx) {
      if (pair_contact_cnt_[static_cast<std::size_t>(idx)] * 2 >= post_perc_frames)
        mean_adj[static_cast<std::size_t>(idx >> 6)] |=
            (uint64_t(1) << (idx & 63));
    }

    const int usable_ring = ring_count_;
    const int n_lags = std::min(max_lag_frames, usable_ring / 2);
    result_.ct_lags_ns.assign(static_cast<std::size_t>(n_lags), 0.0);
    result_.ct_values.assign(static_cast<std::size_t>(n_lags), 0.0);
    const double dt_ns =
        (N >= 2) ? (result_.frames[1].time_ns - result_.frames[0].time_ns) : 1.0;
    const int ring_cap = static_cast<int>(adj_ring_.size());
    const std::size_t workers =
        std::max<std::size_t>(1, std::min<std::size_t>(executor.capacity(),
                                                       static_cast<std::size_t>(n_lags)));
    perf_.finalize_executor_workers = workers;
    const bool use_avx2_ct = kernel_allows_avx2();
    const auto ct_begin = std::chrono::steady_clock::now();
    executor.parallel_do(workers, [&](std::size_t worker) {
      const int begin = static_cast<int>(worker * static_cast<std::size_t>(n_lags) /
                                         workers);
      const int end = static_cast<int>((worker + 1) *
                                       static_cast<std::size_t>(n_lags) /
                                       workers);
      for (int lag = begin; lag < end; ++lag) {
        result_.ct_lags_ns[static_cast<std::size_t>(lag)] = lag * dt_ns;
        double num = 0.0;
        double denom = 0.0;
        const int oldest_global = ct_total_frames_ - usable_ring;
        for (int gi = oldest_global; gi + lag < ct_total_frames_; ++gi) {
          const int slot0 = gi % ring_cap;
          const int slotL = (gi + lag) % ring_cap;
          const auto counts =
              use_avx2_ct
                  ? titan_simd_dispatch::xor_overlap_popcount(
                        adj_ring_[static_cast<std::size_t>(slot0)].data(),
                        adj_ring_[static_cast<std::size_t>(slotL)].data(),
                        mean_adj.data(), static_cast<std::size_t>(n_words_))
                  : titan_simd_dispatch::xor_overlap_popcount_scalar(
                        adj_ring_[static_cast<std::size_t>(slot0)].data(),
                        adj_ring_[static_cast<std::size_t>(slotL)].data(),
                        mean_adj.data(), static_cast<std::size_t>(n_words_));
          num += static_cast<double>(counts.overlap);
          denom += static_cast<double>(counts.reference);
        }
        result_.ct_values[static_cast<std::size_t>(lag)] =
            (denom > 0.0) ? (num / denom) : 0.0;
      }
    });

    if (!result_.ct_values.empty() && result_.ct_values.front() > 1e-12) {
      for (double &value : result_.ct_values)
        value /= result_.ct_values.front();
    }
    result_.tau_ns = fit_tau(result_.ct_values, result_.ct_lags_ns);
    perf_.finalize_ct_seconds +=
        std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                      ct_begin)
            .count();
  } else {
    TitanUI::print_info("LLPS: not enough post-percolation frames for C(t).");
  }

  adj_ring_.clear();
  adj_ring_.shrink_to_fit();
  pair_contact_cnt_.clear();
  pair_contact_cnt_.shrink_to_fit();

  if (contacts_on_ && n_stickers_ > 0 && N > 0) {
    const auto contact_map_begin = std::chrono::steady_clock::now();
    expand_contact_map_result(N);
    perf_.finalize_contact_map_seconds +=
        std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                      contact_map_begin)
            .count();
  }

  finalize_contact_dynamics();
  perf_.finalize_seconds += std::chrono::duration<double>(
                                std::chrono::steady_clock::now() -
                                finalize_begin)
                                .count();
}

// ============================================================================
// Contact dynamics
// ============================================================================



// ============================================================================
// Export
// ============================================================================
