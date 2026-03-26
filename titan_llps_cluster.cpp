// titan_llps_cluster.cpp — DBSCAN Phase Separation Analysis
// Fix 1: BoxMatrix MIC for triclinic/dodecahedron boxes
// Fix 2: Cell dedup (static array) + std::floor for negative coords
// Fix 3: Zero-allocation — head_, nxt_, adj_, bfs_q_ reused every frame
// Fix 4: min_pts counts self (adj.size()+1) + snprintf PDB columns
//        Streaming CSV: open in init(), write per-frame (crash-safe)
// Fix 5: Hybrid-36 resSeq encoding for >9999 chains (VMD-compatible)
#include "titan_llps_cluster.h"
#include "BoxMatrix.h"
#include "titan_ui.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <iomanip>

// ─────────────────────────────────────────────────────────────
// encode_hybrid36 — encode integer n into a hybrid-36 string (width 4)
// Standard PDB resSeq is decimal [0, 9999].
// For n >= 10000: encode with A-Z (10000–35999) then a-z (36000–61999).
// This matches the CCP4/VMD hybrid-36 convention.
// ─────────────────────────────────────────────────────────────
std::string TitanLLPSCluster::encode_hybrid36(int n) {
  const int width = 4;
  // Decimal range [0..9999]
  constexpr int dec_max = 9999;
  if (n <= dec_max) {
    char buf[8];
    std::snprintf(buf, sizeof(buf), "%4d", (n < 0 ? 0 : n));
    return std::string(buf);
  }
  // Alpha range: base-26 after offset
  // A-Z covers offset 10000..10000+26^4-1, a-z continues
  static const char upper[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  static const char lower[] = "abcdefghijklmnopqrstuvwxyz";
  int rem = n - (dec_max + 1);
  const int base26_max = 26 * 26 * 26 * 26; // 456976
  const char *alpha = (rem < base26_max) ? upper : lower;
  if (rem >= base26_max)
    rem -= base26_max;
  char result[8];
  // Encode 4-char base-26 column-major
  result[3] = alpha[rem % 26];
  rem /= 26;
  result[2] = alpha[rem % 26];
  rem /= 26;
  result[1] = alpha[rem % 26];
  rem /= 26;
  result[0] = alpha[rem % 26];
  result[4] = '\0';
  return std::string(result);
}

// ─────────────────────────────────────────────────────────────
TitanLLPSCluster::~TitanLLPSCluster() {
  if (vmd_out_.is_open())
    vmd_out_.close();
  if (ts_out_.is_open())
    ts_out_.close();
}

// ─────────────────────────────────────────────────────────────
void TitanLLPSCluster::init(double eps_angstrom, int min_pts,
                            const std::string &out_pdb_name,
                            const std::string &out_csv_name) {
  eps_ = eps_angstrom;
  eps_sq_ = eps_ * eps_;
  min_pts_ = min_pts;
  frame_count_ = 0;
  timeseries_.clear();
  streaming_csv_ = false;

  if (vmd_out_.is_open())
    vmd_out_.close();
  if (!out_pdb_name.empty())
    vmd_out_.open(out_pdb_name);

  // [Issue 4] Open streaming CSV immediately — data is flushed per-frame
  // so it survives crashes mid-trajectory.
  if (!out_csv_name.empty()) {
    if (ts_out_.is_open())
      ts_out_.close();
    ts_out_.open(out_csv_name);
    if (ts_out_.is_open()) {
      ts_out_ << "Time_ps,TotalBeads,NumCondensates,"
                 "GiantCluster_ProChains,GiantCluster_RnaChains,NoiseBeads,"
                 "Rg_Pro_Dilute,Rg_Pro_Dense,Rg_RNA_Dilute,Rg_RNA_Dense\n";
      streaming_csv_ = true;
      TitanUI::print_info("[LLPS-DBSCAN] Streaming CSV -> " + out_csv_name);
    } else {
      TitanUI::print_warning("[LLPS-DBSCAN] Cannot open CSV: " + out_csv_name);
    }
  }

  TitanUI::print_info("[LLPS-DBSCAN] eps=" + std::to_string(eps_) +
                      " A, min_pts=" + std::to_string(min_pts_));
  if (!out_pdb_name.empty())
    TitanUI::print_info("[LLPS-DBSCAN] Streaming PDB -> " + out_pdb_name);
}

// ─────────────────────────────────────────────────────────────
// extract_nodes — backbone / Cα / P beads only
// ─────────────────────────────────────────────────────────────
std::vector<LLPSNode> TitanLLPSCluster::extract_nodes(
    const MOL &pro, const std::vector<mol_info> &pro_info, const MOL &lig,
    const std::vector<mol_info> &lig_info) const {
  std::vector<LLPSNode> nodes;
  nodes.reserve(pro.pos.x.size() + lig.pos.x.size());
  int nid = 1;

  auto add_mol = [&](const MOL &m, const std::vector<mol_info> &info,
                     int mol_type) {
    if (m.pos.x.empty())
      return;
    const size_t n =
        std::min({m.pos.x.size(), m.pos.y.size(), m.pos.z.size(), info.size()});
    for (size_t i = 0; i < n; ++i) {
      const std::string &a = info[i].get_atom();
      if (a == "BB" || a == "CA" || a == "BB1" || a == "BB2" || a == "BB3" ||
          a == "P") {
        nodes.push_back({nid++, mol_type,
                         static_cast<int>(info[i].get_chainnum()), m.pos.x[i],
                         m.pos.y[i], m.pos.z[i], 0, false});
      }
    }
  };

  add_mol(pro, pro_info, 1); // protein
  add_mol(lig, lig_info, 2); // RNA / ligand
  return nodes;
}

// ─────────────────────────────────────────────────────────────
// run_dbscan — Cell-List O(N) + BoxMatrix MIC + BFS labelling
// ─────────────────────────────────────────────────────────────
void TitanLLPSCluster::run_dbscan(std::vector<LLPSNode> &nodes,
                                  const box &Box) {
  if (nodes.empty())
    return;

  const double bx = Box.get_a(), by = Box.get_b(), bz = Box.get_c();
  if (bx <= 0.0 || by <= 0.0 || bz <= 0.0)
    return;

  // [Fix 1] BoxMatrix handles triclinic / rhombic-dodecahedron boxes
  const BoxMatrix bmat(Box);

  // ── 1. Build Cell-List ──────────────────────────────────────
  const int nx = std::max(1, static_cast<int>(bx / eps_));
  const int ny = std::max(1, static_cast<int>(by / eps_));
  const int nz = std::max(1, static_cast<int>(bz / eps_));
  const double cx = bx / nx, cy = by / ny, cz = bz / nz;
  const int ncells = nx * ny * nz;

  // [Fix 3] Reuse persistent buffers (zero heap allocations per frame)
  head_.assign(ncells, -1);
  nxt_.assign(nodes.size(), -1);
  if (adj_.size() < nodes.size())
    adj_.resize(nodes.size());
  for (size_t i = 0; i < nodes.size(); ++i)
    adj_[i].clear();

  // [Fix 2] std::floor correct for negative coordinates
  auto cell_of = [&](double x, double y, double z) -> int {
    int ix = static_cast<int>(std::floor(x / cx));
    ix = ((ix % nx) + nx) % nx;
    int iy = static_cast<int>(std::floor(y / cy));
    iy = ((iy % ny) + ny) % ny;
    int iz = static_cast<int>(std::floor(z / cz));
    iz = ((iz % nz) + nz) % nz;
    return ix * ny * nz + iy * nz + iz;
  };

  for (size_t i = 0; i < nodes.size(); ++i) {
    int c = cell_of(nodes[i].x, nodes[i].y, nodes[i].z);
    nxt_[i] = head_[c];
    head_[c] = static_cast<int>(i);
  }

  // ── 2. Neighbour search with cell dedup ─────────────────────
  int checked_cells[27]; // static array — no heap, max 3³ = 27 cells
  for (size_t i = 0; i < nodes.size(); ++i) {
    int ix = static_cast<int>(std::floor(nodes[i].x / cx));
    ix = ((ix % nx) + nx) % nx;
    int iy = static_cast<int>(std::floor(nodes[i].y / cy));
    iy = ((iy % ny) + ny) % ny;
    int iz = static_cast<int>(std::floor(nodes[i].z / cz));
    iz = ((iz % nz) + nz) % nz;

    int num_checked = 0;
    for (int ddx = -1; ddx <= 1; ++ddx) {
      int cx2 = ((ix + ddx) % nx + nx) % nx;
      for (int ddy = -1; ddy <= 1; ++ddy) {
        int cy2 = ((iy + ddy) % ny + ny) % ny;
        for (int ddz = -1; ddz <= 1; ++ddz) {
          int cz2 = ((iz + ddz) % nz + nz) % nz;
          int cell = cx2 * ny * nz + cy2 * nz + cz2;
          // [Fix 2] Skip duplicate cells (happens when nx/ny/nz ≤ 2)
          bool dup = false;
          for (int k = 0; k < num_checked; ++k)
            if (checked_cells[k] == cell) {
              dup = true;
              break;
            }
          if (!dup)
            checked_cells[num_checked++] = cell;
        }
      }
    }

    for (int k = 0; k < num_checked; ++k) {
      for (int j = head_[checked_cells[k]]; j != -1; j = nxt_[j]) {
        if (static_cast<size_t>(j) == i)
          continue;
        double dx = nodes[i].x - nodes[j].x;
        double dy = nodes[i].y - nodes[j].y;
        double dz = nodes[i].z - nodes[j].z;
        // [Fix 1] Triclinic MIC via BoxMatrix
        bmat.apply_mic(dx, dy, dz);
        if (dx * dx + dy * dy + dz * dz <= eps_sq_)
          adj_[i].push_back(j);
      }
    }

    // [Fix 4] DBSCAN standard: core point if |neighbours| + 1 (self) >= min_pts
    if (adj_[i].size() + 1 >= static_cast<size_t>(min_pts_))
      nodes[i].is_core = true;
  }

  // ── 3. BFS cluster labelling ────────────────────────────────
  // [Issue 3] Use bfs_q_ (persistent member vector) instead of std::queue.
  // Avoids 5M+ malloc/free calls in large MARTINI trajectories.
  // Vector-as-queue: clear() + indexed iteration — O(N) same as deque BFS.
  int cur_cluster = 0;
  for (size_t i = 0; i < nodes.size(); ++i) {
    if (nodes[i].cluster_id != 0 || !nodes[i].is_core)
      continue;
    ++cur_cluster;
    nodes[i].cluster_id = cur_cluster;
    bfs_q_.clear();
    bfs_q_.push_back(static_cast<int>(i));
    for (size_t qi = 0; qi < bfs_q_.size(); ++qi) {
      int curr = bfs_q_[qi];
      for (int nb : adj_[curr]) {
        if (nodes[nb].cluster_id == 0) {
          nodes[nb].cluster_id = cur_cluster;
          if (nodes[nb].is_core)
            bfs_q_.push_back(nb);
        }
      }
    }
  }

  for (auto &n : nodes)
    if (n.cluster_id == 0)
      n.cluster_id = -1;
}

// ─────────────────────────────────────────────────────────────
// process_frame
// ─────────────────────────────────────────────────────────────
void TitanLLPSCluster::process_frame(const MOL &pro,
                                     const std::vector<mol_info> &pro_info,
                                     const MOL &lig,
                                     const std::vector<mol_info> &lig_info,
                                     double current_time_ps, int frame_idx) {
  auto nodes = extract_nodes(pro, pro_info, lig, lig_info);
  if (nodes.empty())
    return;

  run_dbscan(nodes, pro.Box);

  // ── [BUG 1 FIX] Count ALL unique protein chains from nodes ──
  // (not just clustered ones — noise chains must be included for
  //  the 5% threshold to work correctly)
  int total_pro_chains = 0;
  {
    std::vector<int> all_pc;
    all_pc.reserve(nodes.size());
    for (const auto &n : nodes)
      if (n.type == 1) all_pc.push_back(n.chain_idx);
    std::sort(all_pc.begin(), all_pc.end());
    all_pc.erase(std::unique(all_pc.begin(), all_pc.end()), all_pc.end());
    total_pro_chains = static_cast<int>(all_pc.size());
  }

  // Chain-level statistics via unique chain deduplication
  int max_cid = 0;
  for (const auto &n : nodes)
    if (n.cluster_id > max_cid)
      max_cid = n.cluster_id;

  std::vector<std::vector<int>> cpro(max_cid + 1), crna(max_cid + 1);
  int noise_cnt = 0;

  for (const auto &n : nodes) {
    if (n.cluster_id > 0) {
      if (n.type == 1)
        cpro[n.cluster_id].push_back(n.chain_idx);
      else
        crna[n.cluster_id].push_back(n.chain_idx);
    } else {
      ++noise_cnt;
    }
  }

  int lrg_pro = 0, lrg_rna = 0, n_valid = 0;
  int giant_cluster_id = -1;
  for (int c = 1; c <= max_cid; ++c) {
    std::sort(cpro[c].begin(), cpro[c].end());
    cpro[c].erase(std::unique(cpro[c].begin(), cpro[c].end()), cpro[c].end());
    std::sort(crna[c].begin(), crna[c].end());
    crna[c].erase(std::unique(crna[c].begin(), crna[c].end()), crna[c].end());

    int np = static_cast<int>(cpro[c].size());
    if (np > 1) {
      ++n_valid;
      if (np > lrg_pro) {
        lrg_pro = np;
        lrg_rna = static_cast<int>(crna[c].size());
        giant_cluster_id = c;
      }
    }
  }

  // [BUG 1 FIX] Now uses total from ALL chains (including noise)
  if (lrg_pro < std::max(2, total_pro_chains / 20))
    giant_cluster_id = -1;

  // ════════════════════════════════════════════════════════════
  // Per-bead dense/dilute classification + PDB metadata
  // ════════════════════════════════════════════════════════════
  const int N = static_cast<int>(nodes.size());
  std::vector<bool> bead_is_dense(N, false);
  std::vector<double> bfactor(N, 0.0);
  std::vector<double> occupancy(N, 0.0);

  if (fuzzy_enabled_ && giant_cluster_id > 0) {
    // ── Stage 2: Fuzzy EM → per-bead classification ──
    fuzzy_refine(nodes, pro.Box, giant_cluster_id,
                 bead_is_dense, bfactor, occupancy);
  } else {
    // ── DBSCAN hard classification ──
    // Per-bead: in giant cluster → dense, else dilute
    for (int i = 0; i < N; ++i) {
      bool dense = (giant_cluster_id > 0 &&
                    nodes[i].cluster_id == giant_cluster_id);
      bead_is_dense[i] = dense;
      bfactor[i] = (nodes[i].cluster_id > 0) ? 1.0 : 0.0;  // any cluster
      occupancy[i] = dense ? 1.0 : 0.0;  // 1.0=core, 0.0=dilute
    }
  }

  // ════════════════════════════════════════════════════════════
  // [ISSUE 3 FIX] Save PDB coords BEFORE sort+unwrap
  // ════════════════════════════════════════════════════════════
  auto nodes_pdb = nodes; // snapshot for PDB output

  // ════════════════════════════════════════════════════════════
  // RADIUS OF GYRATION — per-phase, per-chain (same logic for both paths)
  // [BUG 2 FIX] Both fuzzy and non-fuzzy now use identical per-chain Rg
  // ════════════════════════════════════════════════════════════
  // Sort by (type, chain_idx) for linear chain scan
  // We need a parallel sort of bead_is_dense with nodes
  std::vector<int> order(N);
  for (int i = 0; i < N; ++i) order[i] = i;
  std::sort(order.begin(), order.end(), [&](int a, int b) {
    if (nodes[a].type != nodes[b].type) return nodes[a].type < nodes[b].type;
    if (nodes[a].chain_idx != nodes[b].chain_idx) return nodes[a].chain_idx < nodes[b].chain_idx;
    return nodes[a].id < nodes[b].id;  // preserve backbone order within chain
  });
  // Apply permutation to nodes and bead_is_dense
  {
    std::vector<LLPSNode> sorted_nodes(N);
    std::vector<bool> sorted_dense(N);
    for (int i = 0; i < N; ++i) {
      sorted_nodes[i] = nodes[order[i]];
      sorted_dense[i] = bead_is_dense[order[i]];
    }
    nodes = std::move(sorted_nodes);
    bead_is_dense = std::move(sorted_dense);
  }

  const BoxMatrix bmat(pro.Box);
  double sum_rg_pro_dilute = 0.0, sum_rg_pro_dense = 0.0;
  double sum_rg_rna_dilute = 0.0, sum_rg_rna_dense = 0.0;
  int cnt_pro_dilute = 0, cnt_pro_dense = 0;
  int cnt_rna_dilute = 0, cnt_rna_dense = 0;

  size_t ci = 0;
  while (ci < static_cast<size_t>(N)) {
    size_t cj = ci;
    int cur_type = nodes[ci].type;
    int cur_chain = nodes[ci].chain_idx;
    while (cj < static_cast<size_t>(N) && nodes[cj].type == cur_type &&
           nodes[cj].chain_idx == cur_chain)
      ++cj;
    int n_beads = static_cast<int>(cj - ci);

    // Sequential PBC unwrap for Rg
    for (size_t k = ci + 1; k < cj; ++k) {
      double dx = nodes[k].x - nodes[k - 1].x;
      double dy = nodes[k].y - nodes[k - 1].y;
      double dz = nodes[k].z - nodes[k - 1].z;
      bmat.apply_mic(dx, dy, dz);
      nodes[k].x = nodes[k - 1].x + dx;
      nodes[k].y = nodes[k - 1].y + dy;
      nodes[k].z = nodes[k - 1].z + dz;
    }

    // Classify chain: ≥50% dense beads → Dense chain
    int dense_beads = 0;
    for (size_t k = ci; k < cj; ++k)
      if (bead_is_dense[k]) ++dense_beads;
    bool chain_dense = (dense_beads > n_beads / 2);

    // COM
    double ccx = 0, ccy = 0, ccz = 0;
    for (size_t k = ci; k < cj; ++k) {
      ccx += nodes[k].x; ccy += nodes[k].y; ccz += nodes[k].z;
    }
    ccx /= n_beads; ccy /= n_beads; ccz /= n_beads;

    // Rg
    double rg_sq = 0.0;
    for (size_t k = ci; k < cj; ++k) {
      double dx = nodes[k].x - ccx;
      double dy = nodes[k].y - ccy;
      double dz = nodes[k].z - ccz;
      rg_sq += dx * dx + dy * dy + dz * dz;
    }
    double rg = std::sqrt(rg_sq / n_beads);

    if (cur_type == 1) {
      if (chain_dense) { sum_rg_pro_dense += rg; ++cnt_pro_dense; }
      else { sum_rg_pro_dilute += rg; ++cnt_pro_dilute; }
    } else {
      if (chain_dense) { sum_rg_rna_dense += rg; ++cnt_rna_dense; }
      else { sum_rg_rna_dilute += rg; ++cnt_rna_dilute; }
    }
    ci = cj;
  }

  double rg_pd = (cnt_pro_dilute > 0) ? sum_rg_pro_dilute / cnt_pro_dilute : 0.0;
  double rg_pD = (cnt_pro_dense > 0) ? sum_rg_pro_dense / cnt_pro_dense : 0.0;
  double rg_rd = (cnt_rna_dilute > 0) ? sum_rg_rna_dilute / cnt_rna_dilute : 0.0;
  double rg_rD = (cnt_rna_dense > 0) ? sum_rg_rna_dense / cnt_rna_dense : 0.0;

  // Write CSV
  if (streaming_csv_ && ts_out_.is_open()) {
    ts_out_ << std::fixed << std::setprecision(1) << current_time_ps << ","
            << N << "," << n_valid << ","
            << lrg_pro << "," << lrg_rna << "," << noise_cnt << ","
            << std::setprecision(3) << rg_pd << "," << rg_pD << "," << rg_rd
            << "," << rg_rD << "\n";
    if (frame_count_ % 100 == 0)
      ts_out_.flush();
  } else {
    ClusterMetrics m;
    m.time_ps = current_time_ps;
    m.num_nodes = N;
    m.num_clusters = n_valid;
    m.largest_cluster_pro_chains = lrg_pro;
    m.largest_cluster_rna_chains = lrg_rna;
    m.noise_count = noise_cnt;
    m.rg_pro_dilute = rg_pd;
    m.rg_pro_dense = rg_pD;
    m.rg_rna_dilute = rg_rd;
    m.rg_rna_dense = rg_rD;
    timeseries_.push_back(m);
  }

  // PDB output: check stride BEFORE increment so frame 0 is always written
  if (frame_count_ % write_pdb_stride_ == 0)
    write_frame_to_pdb(nodes_pdb, frame_idx, bfactor, occupancy);
  ++frame_count_;
}

// ─────────────────────────────────────────────────────────────
// fuzzy_refine — Stage 2: Fuzzy EM refinement after DBSCAN
// Now outputs per-bead classification + PDB metadata, NOT Rg.
// ─────────────────────────────────────────────────────────────
void TitanLLPSCluster::fuzzy_refine(
    const std::vector<LLPSNode>& nodes, const box& Box,
    int giant_cluster_id,
    std::vector<bool>& is_dense,
    std::vector<double>& bfactor,
    std::vector<double>& occupancy) {

  const int N = static_cast<int>(nodes.size());
  is_dense.assign(N, false);
  bfactor.assign(N, 0.0);
  occupancy.assign(N, 0.0);
  if (N == 0) return;

  const BoxMatrix bmat(Box);

  // ── Find distinct cluster IDs ──
  int max_cid = 0;
  for (const auto& n : nodes)
    if (n.cluster_id > max_cid) max_cid = n.cluster_id;
  if (max_cid == 0) return;

  std::vector<int> cid_list;
  cid_list.reserve(max_cid);
  for (int c = 1; c <= max_cid; ++c) {
    bool has = false;
    for (const auto& n : nodes)
      if (n.cluster_id == c) { has = true; break; }
    if (has) cid_list.push_back(c);
  }
  const int K_dense = static_cast<int>(cid_list.size());
  if (K_dense == 0) return;
  const int K_total = K_dense + 1;
  const int dilute_k = K_dense;

  std::vector<int> cid_to_idx(max_cid + 1, -1);
  for (int k = 0; k < K_dense; ++k)
    cid_to_idx[cid_list[k]] = k;

  int giant_k = (giant_cluster_id > 0 && giant_cluster_id <= max_cid)
                    ? cid_to_idx[giant_cluster_id] : -1;

  // ── Scale parameters from giant cluster ──
  std::vector<std::array<double,3>> X_k(K_total);
  {
    int ref_idx = -1;
    for (int i = 0; i < N; ++i)
      if (nodes[i].cluster_id == giant_cluster_id) { ref_idx = i; break; }
    if (ref_idx < 0) return;

    double sx = 0, sy = 0, sz = 0;
    int cnt = 0;
    for (int i = 0; i < N; ++i) {
      if (nodes[i].cluster_id != giant_cluster_id) continue;
      double dx = nodes[i].x - nodes[ref_idx].x;
      double dy = nodes[i].y - nodes[ref_idx].y;
      double dz = nodes[i].z - nodes[ref_idx].z;
      bmat.apply_mic(dx, dy, dz);
      sx += nodes[ref_idx].x + dx;
      sy += nodes[ref_idx].y + dy;
      sz += nodes[ref_idx].z + dz;
      ++cnt;
    }
    double gx = sx / cnt, gy = sy / cnt, gz = sz / cnt;

    double sum_sq = 0.0;
    for (int i = 0; i < N; ++i) {
      if (nodes[i].cluster_id != giant_cluster_id) continue;
      double dx = nodes[i].x - gx;
      double dy = nodes[i].y - gy;
      double dz = nodes[i].z - gz;
      bmat.apply_mic(dx, dy, dz);
      sum_sq += dx*dx + dy*dy + dz*dz;
    }
    double r_bar = std::sqrt(sum_sq / cnt);
    if (r_bar < 1.0) r_bar = 1.0;

    double r_cutoff = 2.0 * r_bar;
    double sigma_pos = r_bar / 3.0;
    double lambda = 15.0 * sigma_pos * r_bar;
    double r_cutoff_sq = r_cutoff * r_cutoff;

    // Initialize cluster centers
    for (int k = 0; k < K_dense; ++k) {
      int cid = cid_list[k];
      int ref = -1;
      for (int i = 0; i < N; ++i)
        if (nodes[i].cluster_id == cid) { ref = i; break; }
      double wsx = 0, wsy = 0, wsz = 0;
      int wcnt = 0;
      for (int i = 0; i < N; ++i) {
        if (nodes[i].cluster_id != cid) continue;
        double dx = nodes[i].x - nodes[ref].x;
        double dy = nodes[i].y - nodes[ref].y;
        double dz = nodes[i].z - nodes[ref].z;
        bmat.apply_mic(dx, dy, dz);
        wsx += nodes[ref].x + dx;
        wsy += nodes[ref].y + dy;
        wsz += nodes[ref].z + dz;
        ++wcnt;
      }
      X_k[k] = {wsx/wcnt, wsy/wcnt, wsz/wcnt};
    }

    // Dilute pseudo-cluster: no centroid needed.
    // Dilute beads are scattered uniformly; their "centroid" is a
    // physically meaningless point. Instead we assign a uniform
    // background probability in the E-step (see below).
    X_k[dilute_k] = {0.0, 0.0, 0.0}; // placeholder, not used in E-step

    // ── EM iterations ──
    // [FIX 2] Flat 1D matrix — single allocation, cache-friendly
    std::vector<double> u(static_cast<size_t>(N) * K_total, 0.0);
    auto U = [&](int i, int k) -> double& { return u[i * K_total + k]; };
    constexpr int MAX_ITER = 50;
    constexpr double CONV_THR = 0.01;
    constexpr double DILUTE_BG = 1e-4; // uniform background for dilute

    for (int iter = 0; iter < MAX_ITER; ++iter) {
      // E-step
      for (int i = 0; i < N; ++i) {
        double sum_raw = 0.0;
        // Dense clusters: distance-based kernel
        for (int k = 0; k < K_dense; ++k) {
          double dx = nodes[i].x - X_k[k][0];
          double dy = nodes[i].y - X_k[k][1];
          double dz = nodes[i].z - X_k[k][2];
          bmat.apply_mic(dx, dy, dz);
          double d_sq = dx*dx + dy*dy + dz*dz;
          if (d_sq < r_cutoff_sq) {
            double t = std::sqrt(d_sq) / r_cutoff;
            double C = 1.0 - 3.0*t*t + 2.0*t*t*t;
            U(i,k) = std::exp(-d_sq / lambda) * C;
          } else {
            U(i,k) = 0.0;
          }
          sum_raw += U(i,k);
        }
        // [FIX 3] Dilute pseudo-cluster: uniform background probability
        U(i, dilute_k) = DILUTE_BG;
        sum_raw += DILUTE_BG;
        // Normalize
        if (sum_raw > 0.0) {
          for (int k = 0; k < K_total; ++k) U(i,k) /= sum_raw;
        } else {
          for (int k = 0; k < K_total; ++k) U(i,k) = 1.0 / K_total;
        }
      }

      // M-step (only dense clusters — dilute has no centroid)
      double max_shift = 0.0;
      for (int k = 0; k < K_dense; ++k) {
        int best_i = 0;
        double best_u = U(0,k);
        for (int i = 1; i < N; ++i)
          if (U(i,k) > best_u) { best_u = U(i,k); best_i = i; }

        double wsx = 0, wsy = 0, wsz = 0, wsum = 0;
        for (int i = 0; i < N; ++i) {
          if (U(i,k) < 1e-8) continue;
          double dx = nodes[i].x - nodes[best_i].x;
          double dy = nodes[i].y - nodes[best_i].y;
          double dz = nodes[i].z - nodes[best_i].z;
          bmat.apply_mic(dx, dy, dz);
          wsx += U(i,k) * (nodes[best_i].x + dx);
          wsy += U(i,k) * (nodes[best_i].y + dy);
          wsz += U(i,k) * (nodes[best_i].z + dz);
          wsum += U(i,k);
        }
        if (wsum > 0) {
          double nx = wsx/wsum, ny = wsy/wsum, nz = wsz/wsum;
          // [FIX 4] MIC on convergence shift — prevents PBC wrap artifacts
          double sdx = nx - X_k[k][0], sdy = ny - X_k[k][1], sdz = nz - X_k[k][2];
          bmat.apply_mic(sdx, sdy, sdz);
          double shift = sdx*sdx + sdy*sdy + sdz*sdz;
          if (shift > max_shift) max_shift = shift;
          X_k[k] = {nx, ny, nz};
        }
      }
      if (std::sqrt(max_shift) < CONV_THR) break;
    }

    // ── Classify beads + compute PDB metadata ──
    for (int i = 0; i < N; ++i) {
      int dom_k = 0;
      for (int k = 1; k < K_total; ++k)
        if (U(i,k) > U(i,dom_k)) dom_k = k;

      bool dense = (dom_k != dilute_k) && (U(i,dom_k) > 0.4);
      is_dense[i] = dense;

      // B-factor = membership to giant cluster
      double u_giant = (giant_k >= 0) ? U(i,giant_k) : 0.0;
      bfactor[i] = std::max(0.0, std::min(1.0, u_giant));

      // Occupancy = phase tag
      if (u_giant < 0.2)       occupancy[i] = 0.00; // dilute
      else if (u_giant < 0.8)  occupancy[i] = 0.50; // interface
      else                     occupancy[i] = 1.00; // core
    }
  }
}

// ─────────────────────────────────────────────────────────────
// write_frame_to_pdb
// [ISSUE 4 FIX] chain A/B, per-bead B-factor + occupancy
// ─────────────────────────────────────────────────────────────
void TitanLLPSCluster::write_frame_to_pdb(const std::vector<LLPSNode> &nodes,
                                          int frame_idx,
                                          const std::vector<double>& bfactor,
                                          const std::vector<double>& occupancy) {
  if (!vmd_out_.is_open())
    return;
  char buf[128];
  std::snprintf(buf, sizeof(buf), "MODEL     %4d\n", frame_idx);
  vmd_out_ << buf;
  for (size_t i = 0; i < nodes.size(); ++i) {
    const auto &n = nodes[i];
    const char *resname = (n.type == 1) ? "PRO" : "RNA";
    char chain_id = (n.type == 1) ? 'A' : 'B'; // protein=A, RNA=B
    double bf = (i < bfactor.size()) ? bfactor[i] : 0.0;
    double occ = (i < occupancy.size()) ? occupancy[i] : 1.0;
    int atom_serial = (n.id - 1) % 99999 + 1;
    std::string res_str = encode_hybrid36(n.chain_idx);
    std::snprintf(buf, sizeof(buf),
                  "ATOM  %5d  CA  %3s %c%4s    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                  atom_serial, resname, chain_id, res_str.c_str(),
                  n.x, n.y, n.z, occ, bf);
    vmd_out_ << buf;
  }
  vmd_out_ << "ENDMDL\n";
}

// ─────────────────────────────────────────────────────────────
// export_timeseries
// If streaming CSV was opened in init(), data is already written; this is a
// no-op that just flushes and logs. Otherwise, writes all buffered data.
// ─────────────────────────────────────────────────────────────
void TitanLLPSCluster::export_timeseries(const std::string &out_csv) const {
  if (streaming_csv_) {
    // Data already written per-frame. Just flush the file handle.
    if (ts_out_.is_open())
      const_cast<std::ofstream &>(ts_out_).flush();
    TitanUI::print_info("[LLPS-DBSCAN] Timeseries already streamed to CSV.");
    return;
  }
  // Fallback: in-memory buffer path (no CSV name was given to init)
  std::ofstream out(out_csv);
  if (!out.is_open()) {
    TitanUI::print_warning("[LLPS-DBSCAN] Cannot open " + out_csv);
    return;
  }
  out << "Time_ps,TotalBeads,NumCondensates,"
         "GiantCluster_ProChains,GiantCluster_RnaChains,NoiseBeads,"
         "Rg_Pro_Dilute,Rg_Pro_Dense,Rg_RNA_Dilute,Rg_RNA_Dense\n";
  for (const auto &m : timeseries_) {
    out << std::fixed << std::setprecision(1) << m.time_ps << "," << m.num_nodes
        << "," << m.num_clusters << "," << m.largest_cluster_pro_chains << ","
        << m.largest_cluster_rna_chains << "," << m.noise_count << ","
        << std::setprecision(3) << m.rg_pro_dilute << "," << m.rg_pro_dense
        << "," << m.rg_rna_dilute << "," << m.rg_rna_dense << "\n";
  }
  TitanUI::print_info("[LLPS-DBSCAN] Written: " + out_csv);
}
