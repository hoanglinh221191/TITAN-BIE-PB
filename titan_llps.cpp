// titan_llps.cpp — LLPS Analysis (Martini CG, macro-node level)
// Corrected: MOL coordinate API, dimer macro-nodes, PBC-post-solve hook
#include "titan_llps.h"
#include "memory_utils.h"
#include "titan_gyration.h"
#include "titan_simd_dispatch.h"
#include "titan_ui.h"

#include <algorithm>
#include <bit> // std::popcount (C++20)
#include <cmath>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <thread>
#include <unordered_map>

// ============================================================================
// Union-Find
// ============================================================================
void TitanLLPS::UF::init(int n) {
  parent.resize(n);
  rank_.resize(n, 0);
  size_.resize(n, 1);
  std::iota(parent.begin(), parent.end(), 0);
}
void TitanLLPS::UF::reset(int n) {
  // Fix #4: no realloc, just fill existing buffers
  std::iota(parent.begin(), parent.begin() + n, 0);
  std::fill(rank_.begin(), rank_.begin() + n, 0);
  std::fill(size_.begin(), size_.begin() + n, 1);
}
int TitanLLPS::UF::find(int x) {
  while (parent[x] != x) {
    parent[x] = parent[parent[x]];
    x = parent[x];
  }
  return x;
}
void TitanLLPS::UF::merge(int a, int b) {
  a = find(a);
  b = find(b);
  if (a == b)
    return;
  if (rank_[a] < rank_[b])
    std::swap(a, b);
  parent[b] = a;
  size_[a] += size_[b];
  if (rank_[a] == rank_[b])
    ++rank_[a];
}
int TitanLLPS::UF::largest_component() const {
  return static_cast<int>(*std::max_element(size_.begin(), size_.end()));
}

// ============================================================================
// Constructor — builds dimer macro-node map
// ============================================================================
TitanLLPS::TitanLLPS(int n_chains, const std::vector<int> &chain_start,
                     const std::vector<int> &chain_len, float r_cut_nm,
                     float perc_thr, int dimer_size, bool contacts_on,
                     float voxel_size_nm)
    : n_chains_(n_chains), dimer_size_(std::max(1, dimer_size)),
      r_cut_(r_cut_nm), r_cut2_(r_cut_nm * r_cut_nm), perc_thr_(perc_thr),
      contacts_on_(contacts_on), chain_start_(chain_start),
      chain_len_(chain_len), n_stickers_(0), voxel_size_nm_(voxel_size_nm) {
  // Build macro-node → raw chain indices
  n_macro_ = (n_chains_ + dimer_size_ - 1) / dimer_size_;
  macro_chains_.resize(n_macro_);
  for (int c = 0; c < n_chains_; ++c)
    macro_chains_[c / dimer_size_].push_back(c);

  bb_atom_.assign(n_macro_, -1);

  n_pairs_ = n_macro_ * (n_macro_ - 1) / 2;
  n_words_ = (n_pairs_ + 63) / 64;

  uf_.init(n_macro_);

  // Build atom-to-chain map
  int total_atoms = 0;
  for (int c = 0; c < n_chains_; ++c)
    if (chain_start_[c] >= 0)
      total_atoms = std::max(total_atoms, chain_start_[c] + chain_len_[c]);
  atom_chain_.assign(total_atoms, -1);
  for (int c = 0; c < n_chains_; ++c)
    if (chain_start_[c] >= 0)
      for (int a = chain_start_[c]; a < chain_start_[c] + chain_len_[c]; ++a)
        atom_chain_[a] = c;

  // Fix #1: Build atom_macro_ once (atom → macro-node)
  atom_macro_.assign(total_atoms, -1);
  for (int m = 0; m < n_macro_; ++m)
    for (int c : macro_chains_[m])
      if (chain_start_[c] >= 0)
        for (int a = chain_start_[c]; a < chain_start_[c] + chain_len_[c]; ++a)
          if (a < total_atoms)
            atom_macro_[a] = m;

  // Fix #2: Pre-allocate seen generation counter
  seen_gen_.assign(n_macro_ * n_macro_, 0);
}

// ============================================================================
// set_sticker_beads — parse comma-sep names, detect backbone per macro-node
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
  for (size_t i = 0; i < info.size(); ++i)
    for (const auto &s : sticker_names)
      if (info[i].get_atom() == s) {
        sticker_mask_[i] = true;
        break;
      }

  // Auto-detect backbone (first BB/BB1/BB2 atom) per macro-node
  static const std::vector<std::string> BB = {"BB", "BB1", "BB2", "CA"};
  for (int m = 0; m < n_macro_; ++m) {
    for (int c : macro_chains_[m]) {
      for (int a = chain_start_[c]; a < chain_start_[c] + chain_len_[c]; ++a) {
        if (a >= static_cast<int>(info.size()))
          break;
        for (const auto &b : BB)
          if (info[a].get_atom() == b) {
            bb_atom_[m] = a;
            goto next_macro;
          }
      }
    }
    // fallback: first atom of first chain in macro-node
    if (!macro_chains_[m].empty() && chain_len_[macro_chains_[m][0]] > 0)
      bb_atom_[m] = chain_start_[macro_chains_[m][0]];
  next_macro:;
  }

  sticker_atom_idx_.clear();
  sticker_local_idx_.assign(info.size(), -1);
  result_.sticker_residue_indices.clear();

  // Collapse sticker beads onto residue numbers so llps_contact_map is
  // sequence/residue indexed instead of raw bead indexed.
  std::unordered_map<size_t, int> residue_to_local;
  residue_to_local.reserve(info.size());

  for (size_t i = 0; i < sticker_mask_.size(); ++i) {
    if (!sticker_mask_[i])
      continue;

    sticker_atom_idx_.push_back(static_cast<int>(i));
    const size_t resnr = info[i].get_resnr();
    auto [it, inserted] =
        residue_to_local.emplace(resnr,
                                 static_cast<int>(result_.sticker_residue_indices.size()));
    if (inserted)
      result_.sticker_residue_indices.push_back(static_cast<int>(resnr));
    sticker_local_idx_[i] = it->second;
  }
  n_stickers_ = static_cast<int>(result_.sticker_residue_indices.size());

  if (contacts_on_ && n_stickers_ > 0) {
    contact_accum_.assign(n_stickers_, std::vector<float>(n_stickers_, 0.f));
    sticker_seen_gen_.assign(static_cast<size_t>(n_stickers_) * n_stickers_, 0);
    sticker_cur_gen_ = 0;
  } else {
    contact_accum_.clear();
    sticker_seen_gen_.clear();
    sticker_cur_gen_ = 0;
  }
}

// ============================================================================
// PBC min-image squared distance (nm)
// ============================================================================
float TitanLLPS::pbc_dist2(float dx, float dy, float dz, float bx, float by,
                           float bz) {
  if (bx > 0) {
    float half_bx = bx * 0.5f;
    if (dx > half_bx)
      dx -= bx;
    else if (dx < -half_bx)
      dx += bx;
  }
  if (by > 0) {
    float half_by = by * 0.5f;
    if (dy > half_by)
      dy -= by;
    else if (dy < -half_by)
      dy += by;
  }
  if (bz > 0) {
    float half_bz = bz * 0.5f;
    if (dz > half_bz)
      dz -= bz;
    else if (dz < -half_bz)
      dz += bz;
  }
  return dx * dx + dy * dy + dz * dz;
}

// ============================================================================
// Bit helpers
// ============================================================================
inline int TitanLLPS::pair_index(int i, int j) const {
  return i * (2 * n_macro_ - i - 1) / 2 + (j - i - 1);
}
inline void TitanLLPS::set_adj_bit(std::vector<uint64_t> &words, int i,
                                   int j) const {
  if (i > j)
    std::swap(i, j);
  int idx = pair_index(i, j);
  words[idx >> 6] |= (uint64_t(1) << (idx & 63));
}

// ============================================================================
// Cell-list contact detection (pos in nm, box in nm)
// Returns macro-node pairs in contact
// ============================================================================
void TitanLLPS::build_contacts(
    const std::vector<float> &pos_nm, float bx, float by, float bz,
    std::vector<std::pair<int, int>> &macro_contacts) {
  macro_contacts.clear();
  const int n_atoms = static_cast<int>(pos_nm.size()) / 3;

  // Fix #1: Use precomputed atom_macro_ instead of rebuilding
  // (atom_macro_ built once in constructor)

  // Use loose bounding box based on actual atom positions if box is non-PBC
  float Lx = bx, Ly = by, Lz = bz;
  float shift_x = 0.0f;
  float shift_y = 0.0f;
  float shift_z = 0.0f;
  const bool periodic = (Lx > 0 && Ly > 0 && Lz > 0);
  if (!periodic) {
    // Compute bounding box from coords (non-PBC case)
    float xmin = 1e9, xmax = -1e9, ymin = 1e9, ymax = -1e9, zmin = 1e9,
          zmax = -1e9;
    for (int i = 0; i < n_atoms; ++i) {
      xmin = std::min(xmin, pos_nm[3 * i]);
      xmax = std::max(xmax, pos_nm[3 * i]);
      ymin = std::min(ymin, pos_nm[3 * i + 1]);
      ymax = std::max(ymax, pos_nm[3 * i + 1]);
      zmin = std::min(zmin, pos_nm[3 * i + 2]);
      zmax = std::max(zmax, pos_nm[3 * i + 2]);
    }
    // Add cutoff padding, use as non-periodic cell size
    Lx = xmax - xmin + 2 * r_cut_;
    Ly = ymax - ymin + 2 * r_cut_;
    Lz = zmax - zmin + 2 * r_cut_;
    shift_x = -xmin + r_cut_;
    shift_y = -ymin + r_cut_;
    shift_z = -zmin + r_cut_;
  }

  titan_neighbor::OrthorhombicGridSpec spec;
  spec.cutoff = r_cut_;
  spec.len_x = Lx;
  spec.len_y = Ly;
  spec.len_z = Lz;
  spec.shift_x = shift_x;
  spec.shift_y = shift_y;
  spec.shift_z = shift_z;
  spec.periodic = periodic;

  const auto active_fn = [&](std::size_t a) noexcept {
    return a < atom_macro_.size() && atom_macro_[a] >= 0;
  };
  const auto coord_fn = [&](std::size_t a, double &x, double &y,
                            double &z) noexcept {
    x = static_cast<double>(pos_nm[3 * a + 0]);
    y = static_cast<double>(pos_nm[3 * a + 1]);
    z = static_cast<double>(pos_nm[3 * a + 2]);
  };

  // Fix #2: Use generation counter instead of alloc+zero
  ++cur_gen_;
  ++sticker_cur_gen_;

  contact_grid_.for_each_pair(
      static_cast<std::size_t>(n_atoms), spec, coord_fn, active_fn,
      [&](std::size_t lhs, std::size_t rhs) {
        const int a = static_cast<int>(lhs);
        const int b = static_cast<int>(rhs);
        const int mi = atom_macro_[lhs];
        const int mj = atom_macro_[rhs];
        if (mj < 0 || mj == mi)
          return;

        const float d2 = pbc_dist2(pos_nm[3 * a] - pos_nm[3 * b],
                                   pos_nm[3 * a + 1] - pos_nm[3 * b + 1],
                                   pos_nm[3 * a + 2] - pos_nm[3 * b + 2],
                                   periodic ? bx : 0.0f, periodic ? by : 0.0f,
                                   periodic ? bz : 0.0f);
        if (d2 >= r_cut2_)
          return;

        const int key = std::min(mi, mj) * n_macro_ + std::max(mi, mj);
        if (seen_gen_[static_cast<std::size_t>(key)] != cur_gen_) {
          seen_gen_[static_cast<std::size_t>(key)] = cur_gen_;
          macro_contacts.emplace_back(std::min(mi, mj), std::max(mi, mj));
        }

        if (!contacts_on_ || lhs >= sticker_local_idx_.size() ||
            rhs >= sticker_local_idx_.size()) {
          return;
        }

        int si = sticker_local_idx_[lhs];
        int sj = sticker_local_idx_[rhs];
        if (si < 0 || sj < 0 || si == sj)
          return;
        if (si > sj)
          std::swap(si, sj);
        const int sticker_key = si * n_stickers_ + sj;
        if (sticker_key < static_cast<int>(sticker_seen_gen_.size()) &&
            sticker_seen_gen_[static_cast<std::size_t>(sticker_key)] ==
                sticker_cur_gen_) {
          return;
        }
        if (sticker_key < static_cast<int>(sticker_seen_gen_.size()))
          sticker_seen_gen_[static_cast<std::size_t>(sticker_key)] =
              sticker_cur_gen_;
        contact_accum_[static_cast<std::size_t>(si)]
                     [static_cast<std::size_t>(sj)] += 1.f;
        contact_accum_[static_cast<std::size_t>(sj)]
                     [static_cast<std::size_t>(si)] += 1.f;
      });
}

// ============================================================================
// GCC Rg (nm), PBC-unwrapped around reference backbone atom
// ============================================================================
float TitanLLPS::compute_gcc_rg(const std::vector<float> &pos_nm, float bx,
                                float by, float bz,
                                const std::vector<int> &gcc_macros) const {
  if (gcc_macros.empty())
    return 0.f;
  const int ref = bb_atom_[gcc_macros[0]];
  if (ref < 0)
    return 0.f;
  const float rx = pos_nm[3 * ref], ry = pos_nm[3 * ref + 1],
              rz = pos_nm[3 * ref + 2];

  gcc_rg_x_buf_.resize(gcc_macros.size());
  gcc_rg_y_buf_.resize(gcc_macros.size());
  gcc_rg_z_buf_.resize(gcc_macros.size());

  std::size_t cnt = 0;
  for (int m : gcc_macros) {
    int b = bb_atom_[m];
    if (b < 0)
      continue;
    float dx = pos_nm[3 * b] - rx, dy = pos_nm[3 * b + 1] - ry,
          dz = pos_nm[3 * b + 2] - rz;
    if (bx > 0)
      dx -= bx * std::round(dx / bx);
    if (by > 0)
      dy -= by * std::round(dy / by);
    if (bz > 0)
      dz -= bz * std::round(dz / bz);
    gcc_rg_x_buf_[cnt] = static_cast<double>(rx + dx);
    gcc_rg_y_buf_[cnt] = static_cast<double>(ry + dy);
    gcc_rg_z_buf_[cnt] = static_cast<double>(rz + dz);
    ++cnt;
  }
  if (cnt == 0)
    return 0.f;
  const titan_gyration::PointStats stats = titan_gyration::build_point_stats_soa(
      gcc_rg_x_buf_.data(), gcc_rg_y_buf_.data(), gcc_rg_z_buf_.data(), cnt);
  return static_cast<float>(titan_gyration::compute_rg_from_stats(stats));
}

// ============================================================================
// process_frame — thin wrapper: convert Å → nm, then call impl
// ============================================================================
void TitanLLPS::process_frame(const MOL &mol_frame,
                              const std::vector<mol_info> &info, float box_x_nm,
                              float box_y_nm, float box_z_nm, float time_ns,
                              int frame_idx, const MOL *rna_frame,
                              const std::vector<mol_info> *rna_info) {
  const auto &ax = mol_frame.get_pos_x();
  const auto &ay = mol_frame.get_pos_y();
  const auto &az = mol_frame.get_pos_z();
  const int n_atoms = static_cast<int>(ax.size());

  constexpr float ANG_TO_NM = 0.1f;
  pos_nm_buf_.resize(n_atoms * 3);
  for (int i = 0; i < n_atoms; ++i) {
    pos_nm_buf_[3 * i + 0] = static_cast<float>(ax[i]) * ANG_TO_NM;
    pos_nm_buf_[3 * i + 1] = static_cast<float>(ay[i]) * ANG_TO_NM;
    pos_nm_buf_[3 * i + 2] = static_cast<float>(az[i]) * ANG_TO_NM;
  }

  const float *rna_data = nullptr;
  int n_rna = 0;
  if (rna_frame) {
    const auto &rx = rna_frame->get_pos_x();
    const auto &ry = rna_frame->get_pos_y();
    const auto &rz = rna_frame->get_pos_z();
    n_rna = static_cast<int>(rx.size());
    rna_nm_buf_.resize(n_rna * 3);
    for (int i = 0; i < n_rna; ++i) {
      rna_nm_buf_[3 * i + 0] = static_cast<float>(rx[i]) * ANG_TO_NM;
      rna_nm_buf_[3 * i + 1] = static_cast<float>(ry[i]) * ANG_TO_NM;
      rna_nm_buf_[3 * i + 2] = static_cast<float>(rz[i]) * ANG_TO_NM;
    }
    rna_data = rna_nm_buf_.data();
  }

  process_frame_impl(pos_nm_buf_, rna_data, n_rna, box_x_nm, box_y_nm, box_z_nm,
                     time_ns, frame_idx);
}

// ============================================================================
// process_frame_impl — core logic, operates on pre-converted nm arrays
// ============================================================================
void TitanLLPS::process_frame_impl(const std::vector<float> &pos_nm,
                                   const float *rna_nm_data, int n_rna,
                                   float bx, float by, float bz, float time_ns,
                                   int frame_idx) {
  const int n_atoms = static_cast<int>(pos_nm.size()) / 3;

  // B3: Reuse contacts buffer instead of per-frame heap alloc
  contacts_buf_.clear();
  auto &contacts = contacts_buf_;
  build_contacts(pos_nm, bx, by, bz, contacts);

  // ---- Ring-buffer adjacency (MEMORY-BUDGETED) ----
  // Lazy init: cap ring buffer to 1/4 of currently available RAM
  if (adj_ring_.empty()) {
    // Budget = 1/4 of free RAM (adaptive to system load)
    size_t avail = get_available_ram();
    size_t RING_BUDGET_BYTES =
        avail > 0 ? avail / 4 : 256ULL * 1024 * 1024; // fallback 256 MB
    const size_t bytes_per_slot =
        static_cast<size_t>(n_words_) * sizeof(uint64_t);
    const int ideal_cap = max_lag_frames_ * 2;
    int ring_cap_local = ideal_cap;
    if (bytes_per_slot > 0) {
      int budget_cap = static_cast<int>(RING_BUDGET_BYTES / bytes_per_slot);
      budget_cap =
          std::max(budget_cap, max_lag_frames_ + 1); // at least max_lag+1
      ring_cap_local = std::min(ring_cap_local, budget_cap);
    }
    ring_cap_local = std::max(ring_cap_local, 2); // absolute minimum
    max_ring_cap_ = ring_cap_local;
    adj_ring_.resize(max_ring_cap_, std::vector<uint64_t>(n_words_, 0));
    double total_ring_mb =
        static_cast<double>(max_ring_cap_) * bytes_per_slot / (1024.0 * 1024.0);
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(1);
    oss << "[LLPS] Ring buffer: " << max_ring_cap_ << " slots x "
        << bytes_per_slot << " B/slot = ";
    if (total_ring_mb < 1.0)
      oss << std::setprecision(1)
          << (static_cast<double>(max_ring_cap_) * bytes_per_slot / 1024.0)
          << " KB";
    else
      oss << total_ring_mb << " MB";
    TitanUI::print_info(oss.str());
    ct_num_.assign(max_lag_frames_, 0.0);
    ct_den_.assign(max_lag_frames_, 0.0);
    ct_count_.assign(max_lag_frames_, 0);
    pair_contact_cnt_.assign(n_pairs_, 0);
  }
  // Write into ring slot (overwrite oldest if full)
  const int slot = ct_total_frames_ % max_ring_cap_;
  auto &bits = adj_ring_[slot];
  std::fill(bits.begin(), bits.end(), uint64_t(0));
  for (const auto &[ci, cj] : contacts)
    set_adj_bit(bits, ci, cj);

  // Update pair_contact_cnt_ for mean_adj (incremental)
  for (int w = 0; w < n_words_; ++w) {
    uint64_t mask = bits[w];
    while (mask) {
      int bit = __builtin_ctzll(mask);
      int pair_idx = w * 64 + bit;
      if (pair_idx < n_pairs_)
        pair_contact_cnt_[pair_idx]++;
      mask &= mask - 1;
    }
  }
  ring_count_ = std::min(ring_count_ + 1, max_ring_cap_);
  ++ct_total_frames_;

  // Fix #4: Use reset (no realloc) instead of init
  uf_.reset(n_macro_);
  // B4: Reuse degree buffer instead of per-frame heap alloc
  degree_buf_.assign(n_macro_, 0);
  auto &degree = degree_buf_;
  for (const auto &[ci, cj] : contacts) {
    uf_.merge(ci, cj);
    degree[ci]++;
    degree[cj]++;
  }

  int gcc_size = uf_.largest_component();
  int gcc_root = -1;
  for (int m = 0; m < n_macro_; ++m)
    if (uf_.find(m) == m && uf_.size_[m] == gcc_size) {
      gcc_root = m;
      break;
    }

  // B5: Reuse gcc_macros buffer + merge two scan loops into one
  gcc_macros_buf_.clear();
  auto &gcc_macros = gcc_macros_buf_;
  int n_clusters = 0;
  for (int m = 0; m < n_macro_; ++m) {
    int root = uf_.find(m);
    if (root == m)
      ++n_clusters;
    if (gcc_root >= 0 && root == gcc_root)
      gcc_macros.push_back(m);
  }

  const float gcc_fraction = static_cast<float>(gcc_size) / n_macro_;
  float gcc_rg = compute_gcc_rg(pos_nm, bx, by, bz, gcc_macros);

  const float V_box = (bx > 0 && by > 0 && bz > 0) ? bx * by * bz : 1.f;
  constexpr float PI = 3.14159265f;
  const float V_gcc = (4.f / 3.f) * PI * gcc_rg * gcc_rg * gcc_rg;
  const float V_dil = V_box - V_gcc;
  float rho_dense = (V_gcc > 1e-9f) ? gcc_size / V_gcc : 0.f;
  float rho_dilute = (V_dil > 1e-9f) ? (n_macro_ - gcc_size) / V_dil : 0.f;

  float mean_k = 0.f;
  for (int d : degree)
    mean_k += d;
  mean_k /= n_macro_;

  result_.frames.push_back({frame_idx, time_ns, gcc_fraction, n_clusters,
                            gcc_size, gcc_rg, mean_k, rho_dense, rho_dilute});

  if (!result_.is_percolated && gcc_fraction >= perc_thr_) {
    result_.is_percolated = true;
    result_.percolation_time_ns = time_ns;
    percolation_frame_idx_ = static_cast<int>(result_.frames.size()) - 1;
  }

  // ---- Online P(k) accumulation (no per-frame storage) ----
  int local_max_k = 0;
  for (int d : degree)
    local_max_k = std::max(local_max_k, d);
  if (local_max_k > pk_max_k_) {
    pk_sum_.resize(local_max_k + 1, 0.0);
    pk_sq_.resize(local_max_k + 1, 0.0);
    pk_max_k_ = local_max_k;
  }
  for (int d : degree) {
    double pk = 1.0 / n_macro_; // this atom contributes 1/N to bin d
    pk_sum_[d] += pk;
    pk_sq_[d] += pk * pk;
  }

  // B7: Sticker contacts now accumulated inside build_contacts via cell-list
  // (no more O(N²) brute-force loop here)

  // ========================================================================
  // Feature 1: 3D Volumetric Density Map (GCC-centered)
  // ========================================================================
  if (bx > 0 && by > 0 && bz > 0 && !gcc_macros.empty()) {
    // Step 1 — GCC COM centering (shifted copy)
    double com_x = 0, com_y = 0, com_z = 0;
    int gcc_atom_count = 0;
    for (int m : gcc_macros) {
      for (int c : macro_chains_[m]) {
        for (int a = chain_start_[c]; a < chain_start_[c] + chain_len_[c];
             ++a) {
          com_x += pos_nm[3 * a];
          com_y += pos_nm[3 * a + 1];
          com_z += pos_nm[3 * a + 2];
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

    // Fix #3: Reuse shifted_pos_buf_ instead of per-frame alloc
    shifted_pos_buf_.resize(pos_nm.size());
    auto &shifted_pos = shifted_pos_buf_;
    for (int a = 0; a < n_atoms; ++a) {
      float x = pos_nm[3 * a] + sx;
      float y = pos_nm[3 * a + 1] + sy;
      float z = pos_nm[3 * a + 2] + sz;
      // PBC wrap: x -= bx * floor(x / bx)
      x -= bx * std::floor(x / bx);
      y -= by * std::floor(y / by);
      z -= bz * std::floor(z / bz);
      shifted_pos[3 * a] = x;
      shifted_pos[3 * a + 1] = y;
      shifted_pos[3 * a + 2] = z;
    }

    // Step 2 — Grid allocation (lazy)
    int nx = std::max(1, static_cast<int>(bx / voxel_size_nm_));
    int ny = std::max(1, static_cast<int>(by / voxel_size_nm_));
    int nz = std::max(1, static_cast<int>(bz / voxel_size_nm_));
    size_t needed = static_cast<size_t>(nx) * ny * nz;
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
      TitanUI::print_warning(
          "[LLPS] Density grid (re)allocated: " + std::to_string(nx) + "x" +
          std::to_string(ny) + "x" + std::to_string(nz));
    }

    // Step 3 — Voxel accumulation (protein atoms from pro_in)
    for (int a = 0; a < n_atoms; ++a) {
      int ix = static_cast<int>(shifted_pos[3 * a] / voxel_size_nm_);
      int iy = static_cast<int>(shifted_pos[3 * a + 1] / voxel_size_nm_);
      int iz = static_cast<int>(shifted_pos[3 * a + 2] / voxel_size_nm_);
      ix = std::clamp(ix, 0, grid_nx_ - 1);
      iy = std::clamp(iy, 0, grid_ny_ - 1);
      iz = std::clamp(iz, 0, grid_nz_ - 1);
      size_t flat_idx =
          static_cast<size_t>(ix) * grid_ny_ * grid_nz_ + iy * grid_nz_ + iz;
      vox_pro_[flat_idx]++;
    }
    // Step 3b — RNA density from pre-converted rna_nm_data
    if (rna_nm_data && n_rna > 0) {
      for (int a = 0; a < n_rna; ++a) {
        float x = rna_nm_data[3 * a] + sx;
        float y = rna_nm_data[3 * a + 1] + sy;
        float z = rna_nm_data[3 * a + 2] + sz;
        x -= bx * std::floor(x / bx);
        y -= by * std::floor(y / by);
        z -= bz * std::floor(z / bz);
        int ix =
            std::clamp(static_cast<int>(x / voxel_size_nm_), 0, grid_nx_ - 1);
        int iy =
            std::clamp(static_cast<int>(y / voxel_size_nm_), 0, grid_ny_ - 1);
        int iz =
            std::clamp(static_cast<int>(z / voxel_size_nm_), 0, grid_nz_ - 1);
        vox_rna_[static_cast<size_t>(ix) * grid_ny_ * grid_nz_ + iy * grid_nz_ +
                 iz]++;
      }
      has_rna_ = true;
    }
    ++density_frame_count_;
  }

  // ========================================================================
  // Feature 2: Intra-chain Polymer Scaling (RNA only)
  // ========================================================================
  if (rna_nm_data && n_rna > 0 && bx > 0 && by > 0 && bz > 0) {
    const int N = n_rna;

    // 1. Auto-detect RNA chain boundaries via bond-length threshold (once)
    // Martini bond length < 1.0 nm; if gap > 1.5 nm → new chain
    if (rna_chain_starts_.empty()) {
      rna_chain_starts_.push_back(0);
      const float half_bx = bx * 0.5f, half_by = by * 0.5f, half_bz = bz * 0.5f;
      for (int k = 0; k < N - 1; ++k) {
        float dx = rna_nm_data[3 * k] - rna_nm_data[3 * (k + 1)];
        float dy = rna_nm_data[3 * k + 1] - rna_nm_data[3 * (k + 1) + 1];
        float dz = rna_nm_data[3 * k + 2] - rna_nm_data[3 * (k + 1) + 2];
        // PBC for adjacent-bead bond check (bead may straddle box boundary)
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
        if (dx * dx + dy * dy + dz * dz > 2.25f) // d > 1.5 nm
          rna_chain_starts_.push_back(k + 1);
      }
      rna_chain_starts_.push_back(N);
    }

    // 2. Init scaling arrays based on longest actual RNA chain
    int local_max_len = 0;
    for (size_t c = 0; c < rna_chain_starts_.size() - 1; ++c)
      local_max_len = std::max(local_max_len,
                               rna_chain_starts_[c + 1] - rna_chain_starts_[c]);

    if (max_rna_chain_len_ == 0 && local_max_len > 1) {
      max_rna_chain_len_ = local_max_len;
      sum_R2_s_.assign(max_rna_chain_len_, 0.0);
      count_R2_s_.assign(max_rna_chain_len_, 0);
    }

    // 3. Compute R²(s) PER CHAIN only (no inter-chain mixing, no PBC wrapping)
    for (size_t c = 0; c < rna_chain_starts_.size() - 1; ++c) {
      int start = rna_chain_starts_[c];
      int end = rna_chain_starts_[c + 1];
      int len = end - start;

      for (int s = 1; s < len && s < max_rna_chain_len_; ++s) {
        double r2_sum_local = 0.0;
        for (int k = start; k < end - s; ++k) {
          float dx = rna_nm_data[3 * k] - rna_nm_data[3 * (k + s)];
          float dy = rna_nm_data[3 * k + 1] - rna_nm_data[3 * (k + s) + 1];
          float dz = rna_nm_data[3 * k + 2] - rna_nm_data[3 * (k + s) + 2];
          // NO PBC here: intra-chain distance must be physical, not wrapped
          r2_sum_local += static_cast<double>(dx * dx + dy * dy + dz * dz);
        }
        sum_R2_s_[s] += r2_sum_local;
        count_R2_s_[s] += static_cast<uint64_t>(len - s);
      }
    }
  }
}

// ============================================================================
// finalize — C(t) + P(k), parallelized over lags
// ============================================================================
void TitanLLPS::finalize(int n_threads, int max_lag_frames) {
  const int N = static_cast<int>(result_.frames.size());
  if (N == 0)
    return;

  // ---- P(k) from online accumulators ----
  result_.pk_mean.resize(pk_max_k_ + 1, 0.0);
  result_.pk_std.resize(pk_max_k_ + 1, 0.0);
  for (int k = 0; k <= pk_max_k_; ++k) {
    result_.pk_mean[k] = pk_sum_[k] / N;
    result_.pk_std[k] = std::sqrt(
        std::max(0.0, pk_sq_[k] / N - result_.pk_mean[k] * result_.pk_mean[k]));
  }

  // ---- C(t) using ring buffer ----
  // mean_adj from pair_contact_cnt_ (accumulated over ALL frames)
  const int perc_frame =
      (percolation_frame_idx_ >= 0) ? percolation_frame_idx_ : 0;
  const int post_perc_frames = N - perc_frame;

  if (post_perc_frames < 2 || n_words_ == 0 || ring_count_ < 2) {
    TitanUI::print_info("LLPS: not enough post-percolation frames for C(t).");
    // Free ring buffer memory
    adj_ring_.clear();
    adj_ring_.shrink_to_fit();
    pair_contact_cnt_.clear();
    pair_contact_cnt_.shrink_to_fit();
    goto finalize_contacts;
  }

  {
    // Build mean_adj: bit=1 if pair in contact >50% of post-percolation frames
    std::vector<uint64_t> mean_adj(n_words_, 0);
    for (int idx = 0; idx < n_pairs_; ++idx) {
      if (pair_contact_cnt_[idx] * 2 >= post_perc_frames)
        mean_adj[idx >> 6] |= (uint64_t(1) << (idx & 63));
    }

    // C(t) over ring buffer entries
    // ring buffer holds the last ring_count_ frames;
    // we compute C(t) for lag = 0..min(max_lag_frames, ring_count_/2)
    const int usable_ring = ring_count_;
    const int n_lags = std::min(max_lag_frames, usable_ring / 2);
    result_.ct_lags_ns.resize(n_lags);
    result_.ct_values.resize(n_lags, 0.0);

    // Determine frame time step from result_.frames
    const double dt_ns =
        (N >= 2) ? (result_.frames[1].time_ns - result_.frames[0].time_ns)
                 : 1.0;

    const int ring_cap = static_cast<int>(adj_ring_.size());

    // Map ring indices: the most recent frame is at slot (ct_total_frames_-1) %
    // ring_cap The oldest valid frame is at slot (ct_total_frames_ -
    // ring_count_) % ring_cap We want to iterate over sequential frame pairs
    // within the ring.

    const int n_t = std::max(1, n_threads);
    const int chunk = (n_lags + n_t - 1) / n_t;
    std::vector<std::thread> workers;
    workers.reserve(n_t);
    for (int t = 0; t < n_t; ++t) {
      int s = t * chunk, e = std::min(s + chunk, n_lags);
      if (s >= e)
        continue;
      workers.push_back(std::thread([this, s, e, ring_cap, usable_ring,
                                     &mean_adj, dt_ns]() {
        for (int lag = s; lag < e; ++lag) {
          result_.ct_lags_ns[lag] = lag * dt_ns;
          double num = 0, denom = 0;
          // Iterate over all valid (t0, t0+lag) pairs in ring
          const int oldest_global = ct_total_frames_ - usable_ring;
          for (int gi = oldest_global; gi + lag < ct_total_frames_; ++gi) {
            const int slot0 = gi % ring_cap;
            const int slotL = (gi + lag) % ring_cap;
            const auto &A0 = adj_ring_[slot0];
            const auto &AL = adj_ring_[slotL];
            const auto counts = titan_simd_dispatch::xor_overlap_popcount(
                A0.data(), AL.data(), mean_adj.data(),
                static_cast<std::size_t>(n_words_));
            num += static_cast<double>(counts.overlap);
            denom += static_cast<double>(counts.reference);
          }
          result_.ct_values[lag] = (denom > 0) ? num / denom : 0.0;
        }
      }));
    }
    for (auto &w : workers)
      w.join();

    if (!result_.ct_values.empty() && result_.ct_values[0] > 1e-12)
      for (auto &v : result_.ct_values)
        v /= result_.ct_values[0];

    result_.tau_ns = fit_tau(result_.ct_values, result_.ct_lags_ns);
  }

  // Free ring buffer memory after C(t) is computed
  adj_ring_.clear();
  adj_ring_.shrink_to_fit();
  pair_contact_cnt_.clear();
  pair_contact_cnt_.shrink_to_fit();

finalize_contacts:
  if (contacts_on_ && n_stickers_ > 0 && N > 0) {
    result_.contact_freq = contact_accum_;
    for (auto &row : result_.contact_freq)
      for (auto &v : row)
        v /= static_cast<float>(N);
  }
}

// ============================================================================
// fit_tau
// ============================================================================
double TitanLLPS::fit_tau(const std::vector<double> &ct,
                          const std::vector<double> &lags) const {
  if (ct.size() < 2)
    return 0.0;
  constexpr double inv_e = 1.0 / 2.71828182845904523536;
  for (size_t i = 1; i < ct.size(); ++i) {
    if (ct[i] <= inv_e) {
      double f = (ct[i - 1] - inv_e) / (ct[i - 1] - ct[i]);
      return lags[i - 1] + f * (lags[i] - lags[i - 1]);
    }
  }
  double tau = 0;
  for (size_t i = 1; i < ct.size(); ++i)
    tau += 0.5 * (ct[i - 1] + ct[i]) * (lags[i] - lags[i - 1]);
  return tau;
}

// ============================================================================
// Async pipeline — worker_loop, start_async, process_frame_async, flush_async
// ============================================================================

void TitanLLPS::worker_loop(titan::stop_token stoken) {
  while (true) {
    FrameWorkItem item;
    {
      std::unique_lock<std::mutex> lk(queue_mtx_);
      queue_cv_.wait(lk, [&] {
        return !work_queue_.empty() || stoken.stop_requested() || worker_stop_;
      });
      if (work_queue_.empty()) {
        if (stoken.stop_requested() || worker_stop_)
          return;
        continue;
      }
      item = std::move(work_queue_.front());
      work_queue_.pop();
    }
    // Notify producer that queue has space (backpressure release)
    queue_full_cv_.notify_one();

    // Process using the internal impl path
    process_frame_impl(item.pos_nm, item.rna_nm ? item.rna_nm->data() : nullptr,
                       item.n_rna, item.bx, item.by, item.bz, item.time_ns,
                       item.frame_idx);
  }
}

void TitanLLPS::start_async() {
  async_mode_ = true;
  worker_stop_ = false;
  worker_thread_ =
      titan::jthread([this](titan::stop_token st) { worker_loop(st); });
}

void TitanLLPS::process_frame_async(
    const MOL &mol_frame, const std::vector<mol_info> & /*info*/, float bx,
    float by, float bz, float time_ns, int frame_idx, const MOL *rna_frame,
    const std::vector<mol_info> * /*rna_info*/) {
  // Package: convert Å → nm on the calling thread (cheap O(N) copy)
  const auto &ax = mol_frame.get_pos_x();
  const auto &ay = mol_frame.get_pos_y();
  const auto &az = mol_frame.get_pos_z();
  const int n_atoms = static_cast<int>(ax.size());

  FrameWorkItem item;
  item.pos_nm.resize(n_atoms * 3);
  constexpr float ANG_TO_NM = 0.1f;
  for (int i = 0; i < n_atoms; ++i) {
    item.pos_nm[3 * i + 0] = static_cast<float>(ax[i]) * ANG_TO_NM;
    item.pos_nm[3 * i + 1] = static_cast<float>(ay[i]) * ANG_TO_NM;
    item.pos_nm[3 * i + 2] = static_cast<float>(az[i]) * ANG_TO_NM;
  }
  item.bx = bx;
  item.by = by;
  item.bz = bz;
  item.time_ns = time_ns;
  item.frame_idx = frame_idx;
  item.n_rna = 0;

  if (rna_frame) {
    const auto &rx = rna_frame->get_pos_x();
    const auto &ry = rna_frame->get_pos_y();
    const auto &rz = rna_frame->get_pos_z();
    const int n_rna = static_cast<int>(rx.size());
    item.n_rna = n_rna;
    std::vector<float> rna_buf(n_rna * 3);
    for (int i = 0; i < n_rna; ++i) {
      rna_buf[3 * i + 0] = static_cast<float>(rx[i]) * ANG_TO_NM;
      rna_buf[3 * i + 1] = static_cast<float>(ry[i]) * ANG_TO_NM;
      rna_buf[3 * i + 2] = static_cast<float>(rz[i]) * ANG_TO_NM;
    }
    item.rna_nm = std::move(rna_buf);
  }

  {
    std::unique_lock<std::mutex> lk(queue_mtx_);
    // Backpressure: wait if queue is full
    queue_full_cv_.wait(lk, [&] {
      return static_cast<int>(work_queue_.size()) < MAX_QUEUE_SIZE_;
    });
    work_queue_.push(std::move(item));
  }
  queue_cv_.notify_one();
}

void TitanLLPS::flush_async() {
  if (!async_mode_)
    return;
  // Wait for queue to drain, then stop worker
  {
    std::unique_lock<std::mutex> lk(queue_mtx_);
    queue_cv_.wait(lk, [&] { return work_queue_.empty(); });
    worker_stop_ = true;
  }
  queue_cv_.notify_all();
  if (worker_thread_.joinable())
    worker_thread_.join();
  async_mode_ = false;
}

// ============================================================================
// export_results
// ============================================================================
void TitanLLPS::export_results(const info &info_output,
                               std::fstream &log) const {
  TitanUI::print_module_header("LLPS Analysis — Export");
  const int N = static_cast<int>(result_.frames.size());
  TitanUI::log_kv(log, "LLPS Frames", std::to_string(N));
  TitanUI::log_kv(log, "LLPS Macro-Nodes",
                  std::to_string(n_macro_) +
                      " (dimer_size=" + std::to_string(dimer_size_) + ")");
  TitanUI::log_kv(log, "LLPS Percolated", result_.is_percolated ? "Yes" : "No");
  if (result_.is_percolated)
    TitanUI::log_kv(log, "Percolation Time (ns)",
                    std::to_string(result_.percolation_time_ns));
  TitanUI::log_kv(log, "tau (ns)", std::to_string(result_.tau_ns));

  // 1. Timeseries
  {
    std::ofstream f(info_output.get_llps_timeseries_filename());
    if (f.is_open()) {
      f << "# frame  time_ns  gcc_fraction  gcc_size  gcc_rg_nm  mean_degree  "
           "rho_dense  rho_dilute\n";
      f << std::fixed << std::setprecision(6);
      for (const auto &fr : result_.frames)
        f << fr.frame_idx << "  " << fr.time_ns << "  " << fr.gcc_fraction
          << "  " << fr.gcc_size << "  " << fr.gcc_rg_nm << "  "
          << fr.mean_degree << "  " << fr.density_dense << "  "
          << fr.density_dilute << "\n";
      TitanUI::print_info("  Written: " +
                          info_output.get_llps_timeseries_filename());
    }
  }
  // 2. C(t)
  if (!result_.ct_lags_ns.empty()) {
    std::ofstream f(info_output.get_llps_ct_filename());
    if (f.is_open()) {
      f << "# lag_ns  C(t)\n" << std::fixed << std::setprecision(8);
      for (size_t i = 0; i < result_.ct_lags_ns.size(); ++i)
        f << result_.ct_lags_ns[i] << "  " << result_.ct_values[i] << "\n";
      TitanUI::print_info("  Written: " + info_output.get_llps_ct_filename());
    }
  }
  // 3. Degree distribution
  {
    std::ofstream f(info_output.get_llps_degree_dist_filename());
    if (f.is_open()) {
      f << "# k  P(k)_mean  P(k)_std\n" << std::fixed << std::setprecision(6);
      for (size_t k = 0; k < result_.pk_mean.size(); ++k)
        f << k << "  " << result_.pk_mean[k] << "  " << result_.pk_std[k]
          << "\n";
      TitanUI::print_info("  Written: " +
                          info_output.get_llps_degree_dist_filename());
    }
  }
  // 4. Contact map
  if (contacts_on_ && !result_.contact_freq.empty()) {
    std::ofstream f(info_output.get_llps_contact_map_filename());
    if (f.is_open()) {
      f << "# residue_i  residue_j  contact_frequency\n"
        << "# n_residues=" << n_stickers_ << "\n"
        << "# contact map is aggregated by residue number, not bead index\n"
        << std::fixed << std::setprecision(6);
      for (int i = 0; i < n_stickers_; ++i) {
        for (int j = i + 1; j < n_stickers_; ++j) {
          if (result_.contact_freq[i][j] > 1e-6f) {
            f << result_.sticker_residue_indices[i] << "  "
              << result_.sticker_residue_indices[j] << "  "
              << result_.contact_freq[i][j] << "\n";
          }
        }
      }
      TitanUI::print_info("  Written: " +
                          info_output.get_llps_contact_map_filename());
    }
  }
  // 5. Summary
  {
    std::ofstream f(info_output.get_llps_summary_filename());
    if (f.is_open()) {
      f << "percolated: " << (result_.is_percolated ? "yes" : "no") << "\n";
      f << "percolation_time_ns: " << result_.percolation_time_ns << "\n";
      f << "tau_ns: " << result_.tau_ns << "\n";
      f << "n_macronodes: " << n_macro_ << "\n";
      f << "dimer_size: " << dimer_size_ << "\n";
      double mr = 0, mf = 0;
      for (const auto &fr : result_.frames) {
        mr += fr.gcc_rg_nm;
        mf += fr.gcc_fraction;
      }
      if (N > 0) {
        mr /= N;
        mf /= N;
      }
      f << "mean_gcc_rg_nm: " << mr << "\n";
      f << "mean_gcc_fraction: " << mf << "\n";
      TitanUI::print_info("  Written: " +
                          info_output.get_llps_summary_filename());
    }
  }
  // 6. Density maps (OpenDX format)
  if (density_frame_count_ > 0 && grid_nx_ > 0) {
    const double voxel_vol =
        static_cast<double>(voxel_size_nm_) * voxel_size_nm_ * voxel_size_nm_;
    const size_t ntotal = static_cast<size_t>(grid_nx_) * grid_ny_ * grid_nz_;

    auto write_dx = [&](const std::string &filename,
                        const std::vector<uint64_t> &vox,
                        const std::string &label) {
      std::ofstream f(filename);
      if (!f.is_open())
        return;
      f << "object 1 class gridpositions counts " << grid_nx_ << " " << grid_ny_
        << " " << grid_nz_ << "\n";
      f << "origin 0 0 0\n";
      f << "delta " << voxel_size_nm_ << " 0 0\n";
      f << "delta 0 " << voxel_size_nm_ << " 0\n";
      f << "delta 0 0 " << voxel_size_nm_ << "\n";
      f << "object 2 class gridconnections counts " << grid_nx_ << " "
        << grid_ny_ << " " << grid_nz_ << "\n";
      f << "object 3 class array type double rank 0 items " << ntotal
        << " data follows\n";
      f << std::fixed << std::setprecision(6);
      int col = 0;
      for (size_t i = 0; i < ntotal; ++i) {
        double density =
            static_cast<double>(vox[i]) / (density_frame_count_ * voxel_vol);
        if (col > 0)
          f << ' ';
        f << density;
        if (++col == 3) {
          f << '\n';
          col = 0;
        }
      }
      if (col > 0)
        f << '\n';
      f << "object \"" << label << "\" class field\n";
      TitanUI::print_info("  Written: " + filename);
    };

    write_dx(info_output.get_llps_density_pro_filename(), vox_pro_,
             "density_pro");
    write_dx(info_output.get_llps_density_rna_filename(), vox_rna_,
             "density_rna");
  }
  // 7. Polymer scaling law (RNA)
  if (max_rna_chain_len_ > 1) {
    std::ofstream f(info_output.get_llps_polymer_scaling_filename());
    if (f.is_open()) {
      f << "# s  mean_R2_nm2  R_rms_nm  count\n";
      f << std::fixed << std::setprecision(6);
      for (int s = 1; s < max_rna_chain_len_; ++s) {
        if (count_R2_s_[s] == 0)
          continue;
        double mean_R2 = sum_R2_s_[s] / static_cast<double>(count_R2_s_[s]);
        double R_rms = std::sqrt(mean_R2);
        f << s << "  " << mean_R2 << "  " << R_rms << "  " << count_R2_s_[s]
          << "\n";
      }
      TitanUI::print_info("  Written: " +
                          info_output.get_llps_polymer_scaling_filename());
    }
  }
}
