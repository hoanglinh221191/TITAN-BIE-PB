#include "dccm_network.h"
#include "titan_alignment.h"
#include "titan_ui.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <queue>
#include <thread>


#ifdef __AVX__
#include <immintrin.h>
#endif

// ── Helper: find Cβ (or N for Gly) atom index within a residue ────────────────
int DCCMNetwork::find_cb_for_residue(const Multi_MOL &mol, int resnr,
                                     const std::string &chain, bool is_gly) {
  const std::string target = is_gly ? "N" : "CB";
  for (size_t i = 0; i < mol.info.size(); ++i) {
    std::string at = mol.info[i].get_atom();
    // trim whitespace
    auto p1 = at.find_first_not_of(" \t");
    auto p2 = at.find_last_not_of(" \t");
    if (p1 != std::string::npos)
      at = at.substr(p1, p2 - p1 + 1);
    if (at == target &&
        mol.info[i].get_resnr() == resnr &&
        mol.info[i].get_chain() == chain)
      return static_cast<int>(i);
  }
  return -1; // not found
}

// Hàm dot_product tối ưu bằng tập lệnh SIMD AVX2
inline double fast_dot_product(const double *a, const double *b, size_t n) {
  double total = 0.0;
  size_t i = 0;
#ifdef __AVX__
  __m256d sum_vec = _mm256_setzero_pd();
  for (; i + 3 < n; i += 4) {
    // Dùng loadu (unaligned) an toàn cho std::vector
    __m256d va = _mm256_loadu_pd(&a[i]);
    __m256d vb = _mm256_loadu_pd(&b[i]);
    sum_vec = _mm256_fmadd_pd(va, vb, sum_vec); // Fused Multiply-Add
  }
  double buffer[4];
  _mm256_storeu_pd(buffer, sum_vec);
  total = buffer[0] + buffer[1] + buffer[2] + buffer[3];
#endif
  // Phần dư còn lại nếu n không chia hết cho 4
  for (; i < n; ++i) {
    total += a[i] * b[i];
  }
  return total;
}

void DCCMNetwork::init(const Multi_MOL &pro_in, const Multi_MOL &lig,
                       const arguments &args) {
  // CRITICAL FIX: Clear all state from previous trajectories (Batch Mode
  // Safety)
  ca_indices.clear();
  ca_resnrs.clear();
  ca_is_lig.clear();
  ca_chain.clear();
  ca_resname.clear();
  cb_indices.clear();
  traj_x.clear();
  traj_y.clear();
  traj_z.clear();
  traj_cbx.clear();
  traj_cby.clear();
  traj_cbz.clear();
  total_frames = 0;
  kabsch_aligned_ = false;

  bool is_lig_rec = (args.mol_1_type == "ligand" || args.mol_1_type == "lig");
  bool is_lig_lg = (args.mol_2_type == "ligand" || args.mol_2_type == "lig");

  // Collect Receptor Nodes
  for (size_t i = 0; i < pro_in.info.size(); ++i) {
    std::string atom_name = pro_in.info[i].get_atom();
    if (is_lig_rec || atom_name == "CA" || atom_name == "CA " ||
        atom_name == "BB" || atom_name == "BB2" || atom_name == "C4'" ||
        atom_name == "C4*") {
      ca_indices.emplace_back((int)i);
      ca_resnrs.emplace_back(is_lig_rec ? (i + 1) : pro_in.info[i].get_resnr());
      ca_is_lig.emplace_back(0); // 0 = Receptor
      ca_chain.emplace_back(pro_in.info[i].get_chain());
      ca_resname.emplace_back(is_lig_rec ? pro_in.info[i].get_atom()
                                         : pro_in.info[i].get_res_name());
    }
  }
  // Collect Ligand Nodes
  for (size_t i = 0; i < lig.info.size(); ++i) {
    std::string atom_name = lig.info[i].get_atom();
    if (is_lig_lg || atom_name == "CA" || atom_name == "CA " ||
        atom_name == "BB" || atom_name == "BB2" || atom_name == "C4'" ||
        atom_name == "C4*") {
      ca_indices.emplace_back((int)i);
      ca_resnrs.emplace_back(is_lig_lg ? (i + 1) : lig.info[i].get_resnr());
      ca_is_lig.emplace_back(1); // 1 = Ligand
      ca_chain.emplace_back(lig.info[i].get_chain());
      ca_resname.emplace_back(is_lig_lg ? lig.info[i].get_atom()
                                        : lig.info[i].get_res_name());
    }
  }

  // ── Build cb_indices: for each CA node, find matching Cβ (or N for GLY) ────
  cb_indices.resize(ca_indices.size());
  for (size_t i = 0; i < ca_indices.size(); ++i) {
    // Determine residue name (trimmed)
    std::string rn = ca_resname[i];
    auto p1 = rn.find_first_not_of(" \t");
    auto p2 = rn.find_last_not_of(" \t");
    if (p1 != std::string::npos)
      rn = rn.substr(p1, p2 - p1 + 1);
    bool is_gly = (rn == "GLY");

    const Multi_MOL &mol_ref = (ca_is_lig[i] == 0) ? pro_in : lig;
    cb_indices[i] = find_cb_for_residue(mol_ref, ca_resnrs[i],
                                        ca_chain[i], is_gly);
  }

  total_frames = 0;
}

void DCCMNetwork::process_frame(const Multi_MOL &pro_in, const Multi_MOL &lig,
                                size_t chunk_snap_idx) {
  if (ca_indices.empty())
    return;

  std::vector<double> cx(ca_indices.size()), cy(ca_indices.size()),
      cz(ca_indices.size());
  std::vector<double> cbx(ca_indices.size()), cby(ca_indices.size()),
      cbz(ca_indices.size());
  size_t l_snap = (lig.m_mol.size() > 1) ? chunk_snap_idx : 0;

  for (size_t i = 0; i < ca_indices.size(); ++i) {
    int idx = ca_indices[i];
    if (ca_is_lig[i] == 0) {
      cx[i] = pro_in.m_mol[chunk_snap_idx].pos.x[idx];
      cy[i] = pro_in.m_mol[chunk_snap_idx].pos.y[idx];
      cz[i] = pro_in.m_mol[chunk_snap_idx].pos.z[idx];
    } else {
      cx[i] = lig.m_mol[l_snap].pos.x[idx];
      cy[i] = lig.m_mol[l_snap].pos.y[idx];
      cz[i] = lig.m_mol[l_snap].pos.z[idx];
    }

    // Cβ/N coord: use cb_indices if valid, else fallback to CA
    int cb_idx = cb_indices[i];
    if (cb_idx >= 0) {
      if (ca_is_lig[i] == 0) {
        cbx[i] = pro_in.m_mol[chunk_snap_idx].pos.x[cb_idx];
        cby[i] = pro_in.m_mol[chunk_snap_idx].pos.y[cb_idx];
        cbz[i] = pro_in.m_mol[chunk_snap_idx].pos.z[cb_idx];
      } else {
        cbx[i] = lig.m_mol[l_snap].pos.x[cb_idx];
        cby[i] = lig.m_mol[l_snap].pos.y[cb_idx];
        cbz[i] = lig.m_mol[l_snap].pos.z[cb_idx];
      }
    } else {
      // Fallback: use CA coords
      cbx[i] = cx[i];
      cby[i] = cy[i];
      cbz[i] = cz[i];
    }
  }
  traj_x.emplace_back(std::move(cx));
  traj_y.emplace_back(std::move(cy));
  traj_z.emplace_back(std::move(cz));
  traj_cbx.emplace_back(std::move(cbx));
  traj_cby.emplace_back(std::move(cby));
  traj_cbz.emplace_back(std::move(cbz));
  total_frames++;
}

std::vector<std::vector<std::array<double, 3>>>
DCCMNetwork::get_ca_traj_aos() const {
  std::vector<std::vector<std::array<double, 3>>> aos(
      total_frames, std::vector<std::array<double, 3>>(ca_indices.size()));
  for (size_t f = 0; f < total_frames; ++f) {
    for (size_t a = 0; a < ca_indices.size(); ++a) {
      aos[f][a][0] = traj_x[f][a];
      aos[f][a][1] = traj_y[f][a];
      aos[f][a][2] = traj_z[f][a];
    }
  }
  return aos;
}

void DCCMNetwork::calculate_and_export(size_t num_threads) {
  if (total_frames == 0 || ca_indices.empty())
    return;
  size_t num_ca = ca_indices.size();

  // ── STEP 0: Kabsch align (guarded — skip if already aligned) ───────────────
  if (!kabsch_aligned_) {
    Eigen::MatrixXd ref(num_ca, 3);
    for (size_t a = 0; a < num_ca; ++a) {
      ref(a, 0) = traj_x[0][a];
      ref(a, 1) = traj_y[0][a];
      ref(a, 2) = traj_z[0][a];
    }
    ref.rowwise() -= ref.colwise().mean().eval();

    Eigen::MatrixXd mobile(num_ca, 3);
    Eigen::Matrix3d R;
    for (size_t f = 0; f < total_frames; ++f) {
      for (size_t a = 0; a < num_ca; ++a) {
        mobile(a, 0) = traj_x[f][a];
        mobile(a, 1) = traj_y[f][a];
        mobile(a, 2) = traj_z[f][a];
      }
      Eigen::RowVector3d ca_centroid = mobile.colwise().mean();
      TitanAlignment::align_frame_inplace(ref, mobile, R);
      for (size_t a = 0; a < num_ca; ++a) {
        traj_x[f][a] = mobile(a, 0);
        traj_y[f][a] = mobile(a, 1);
        traj_z[f][a] = mobile(a, 2);
      }
      if (!traj_cbx.empty()) {
        for (size_t a = 0; a < num_ca; ++a) {
          Eigen::RowVector3d cb(traj_cbx[f][a], traj_cby[f][a],
                                traj_cbz[f][a]);
          cb -= ca_centroid;
          cb = (R * cb.transpose()).transpose();
          traj_cbx[f][a] = cb(0);
          traj_cby[f][a] = cb(1);
          traj_cbz[f][a] = cb(2);
        }
      }
    }
    kabsch_aligned_ = true;
  }

  // 1. Pivot Data [atom][frame] for cache locality
  std::vector<std::vector<double>> dev_x(num_ca,
                                         std::vector<double>(total_frames));
  std::vector<std::vector<double>> dev_y(num_ca,
                                         std::vector<double>(total_frames));
  std::vector<std::vector<double>> dev_z(num_ca,
                                         std::vector<double>(total_frames));

  for (size_t a = 0; a < num_ca; ++a) {
    double sum_x = 0, sum_y = 0, sum_z = 0;
    for (size_t f = 0; f < total_frames; ++f) {
      sum_x += traj_x[f][a];
      sum_y += traj_y[f][a];
      sum_z += traj_z[f][a];
    }
    double avg_x = sum_x / static_cast<double>(total_frames);
    double avg_y = sum_y / static_cast<double>(total_frames);
    double avg_z = sum_z / static_cast<double>(total_frames);

    for (size_t f = 0; f < total_frames; ++f) {
      dev_x[a][f] = traj_x[f][a] - avg_x;
      dev_y[a][f] = traj_y[f][a] - avg_y;
      dev_z[a][f] = traj_z[f][a] - avg_z;
    }
  }

  // 2. Tính Ma trận Hiệp phương sai (Multithreaded + AVX2)
  //    FIX: pair-based distribution for balanced workload,
  //         no mirror-write in workers to avoid data race.
  std::vector<std::vector<double>> cov(num_ca,
                                       std::vector<double>(num_ca, 0.0));

  // Build upper-triangle pair list (including diagonal)
  std::vector<std::pair<size_t, size_t>> pairs;
  pairs.reserve(num_ca * (num_ca + 1) / 2);
  for (size_t i = 0; i < num_ca; ++i)
    for (size_t j = i; j < num_ca; ++j)
      pairs.emplace_back(i, j);

  // Worker: each thread processes a contiguous slice of pairs
  auto calc_cov_worker = [&](size_t p_start, size_t p_end) {
    for (size_t p = p_start; p < p_end; ++p) {
      auto [i, j] = pairs[p];
      double dot_x =
          fast_dot_product(dev_x[i].data(), dev_x[j].data(), total_frames);
      double dot_y =
          fast_dot_product(dev_y[i].data(), dev_y[j].data(), total_frames);
      double dot_z =
          fast_dot_product(dev_z[i].data(), dev_z[j].data(), total_frames);

      cov[i][j] = (dot_x + dot_y + dot_z) / static_cast<double>(total_frames);
    }
  };

  // Distribute pairs equally among threads
  std::vector<std::thread> workers;
  size_t total_pairs = pairs.size();
  size_t chunk = (total_pairs + num_threads - 1) / num_threads;
  for (size_t t = 0; t < num_threads; ++t) {
    size_t start = t * chunk;
    size_t end = std::min(start + chunk, total_pairs);
    if (start < end) {
      workers.emplace_back(calc_cov_worker, start, end);
    }
  }
  for (auto &w : workers) {
    if (w.joinable())
      w.join();
  }

  // Single-threaded symmetry fill (safe, no race)
  for (size_t i = 0; i < num_ca; ++i)
    for (size_t j = i + 1; j < num_ca; ++j)
      cov[j][i] = cov[i][j];

  // 3. Tính Correlation và Xuất File
  std::ofstream edges_out("network_edges.csv");
  edges_out << "Source,Target,Correlation\n";
  for (size_t i = 0; i < num_ca; ++i) {
    for (size_t j = i + 1; j < num_ca; ++j) {
      double denom = std::sqrt(cov[i][i] * cov[j][j]);
      if (denom > 1e-12) {
        double corr = cov[i][j] / denom;
        if (std::abs(corr) > 0.4) {
          edges_out << ca_resname[i] << ca_resnrs[i] << "_" << ca_chain[i]
                    << "," << ca_resname[j] << ca_resnrs[j] << "_"
                    << ca_chain[j] << "," << corr << "\n";
        }
      }
    }
  }
  edges_out.close();

  std::ofstream nodes_out("network_nodes.csv");
  nodes_out << "Id,Label,Type,AtomName\n";
  for (size_t i = 0; i < num_ca; ++i) {
    nodes_out << ca_resname[i] << ca_resnrs[i] << "_" << ca_chain[i] << ","
              << ca_resnrs[i] << "," << (ca_is_lig[i] ? "Ligand" : "Receptor")
              << "," << ca_resname[i] << "\n";
  }
  nodes_out.close();
}

// =============================================================================
// Information-Theory Allosteric Analysis: GC + Transfer Entropy
// =============================================================================
void DCCMNetwork::calculate_allostery_information(size_t num_threads,
                                                  int tau_lag, int num_bins) {
  if (total_frames == 0 || ca_indices.empty())
    return;
  size_t num_ca = ca_indices.size();
  size_t tau = static_cast<size_t>(std::max(1, tau_lag));
  if (tau >= total_frames) {
    std::cerr << "[Allostery] tau_lag (" << tau << ") >= total_frames ("
              << total_frames << "). Skipping.\n";
    return;
  }

  // ── STEP 1: Signed Projection ──────────────────────────────────────────────
  // For each atom, compute lightweight 3×3 covariance matrix from trajectory,
  // extract dominant eigenvector, and project displacements onto it.
  // This preserves sign (anti-correlation visible) and reduces to 1D.
  std::vector<std::vector<double>> proj(num_ca,
                                        std::vector<double>(total_frames));
  {
    // Compute per-atom mean (traj already Kabsch-aligned)
    std::vector<double> mean_x(num_ca, 0.0), mean_y(num_ca, 0.0),
        mean_z(num_ca, 0.0);
    for (size_t f = 0; f < total_frames; ++f) {
      for (size_t a = 0; a < num_ca; ++a) {
        mean_x[a] += traj_x[f][a];
        mean_y[a] += traj_y[f][a];
        mean_z[a] += traj_z[f][a];
      }
    }
    double inv_f = 1.0 / static_cast<double>(total_frames);
    for (size_t a = 0; a < num_ca; ++a) {
      mean_x[a] *= inv_f;
      mean_y[a] *= inv_f;
      mean_z[a] *= inv_f;
    }

    // Per-atom 3×3 covariance → dominant eigenvector → signed projection
    for (size_t a = 0; a < num_ca; ++a) {
      Eigen::Matrix3d cov3 = Eigen::Matrix3d::Zero();
      for (size_t f = 0; f < total_frames; ++f) {
        Eigen::Vector3d d(traj_x[f][a] - mean_x[a], traj_y[f][a] - mean_y[a],
                          traj_z[f][a] - mean_z[a]);
        cov3.noalias() += d * d.transpose();
      }
      cov3 /= static_cast<double>(total_frames);

      Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es(cov3);
      // Largest eigenvalue is last (ascending order)
      Eigen::Vector3d e1 = es.eigenvectors().col(2);

      for (size_t f = 0; f < total_frames; ++f) {
        Eigen::Vector3d d(traj_x[f][a] - mean_x[a], traj_y[f][a] - mean_y[a],
                          traj_z[f][a] - mean_z[a]);
        proj[a][f] = d.dot(e1); // signed 1D projection
      }
    }
  }

  // ── STEP 2: Min-Max Scaling to [0, 1] ─────────────────────────────────────
  std::vector<double> p_min(num_ca), p_range(num_ca);
  for (size_t a = 0; a < num_ca; ++a) {
    auto [mn, mx] = std::minmax_element(proj[a].begin(), proj[a].end());
    p_min[a] = *mn;
    p_range[a] = (*mx - *mn);
    if (p_range[a] < 1e-15)
      p_range[a] = 1e-15; // avoid division by zero for frozen atoms
    for (size_t f = 0; f < total_frames; ++f)
      proj[a][f] = (proj[a][f] - p_min[a]) / p_range[a];
  }

  // ── STEP 3: Adaptive Binning ──────────────────────────────────────────────
  int B = num_bins;
  if (B <= 0)
    B = std::max(8, static_cast<int>(std::cbrt(total_frames / 5.0)));

  // ── STEP 4: Precompute per-atom 1D Entropy (shared, read-only) ────────────
  // Each atom's entropy is needed by many pairs, compute once.
  std::vector<double> H1(num_ca, 0.0); // H(X_i)
  {
    std::vector<double> hist(B, 0.0);
    for (size_t a = 0; a < num_ca; ++a) {
      std::fill(hist.begin(), hist.end(), 0.5); // Laplace pseudocount
      for (size_t f = 0; f < total_frames; ++f) {
        int bin = static_cast<int>(proj[a][f] * B);
        if (bin >= B)
          bin = B - 1;
        hist[bin] += 1.0;
      }
      double norm = std::accumulate(hist.begin(), hist.end(), 0.0);
      double H = 0.0;
      for (int b = 0; b < B; ++b) {
        double p = hist[b] / norm;
        if (p > 0.0)
          H -= p * std::log(p);
      }
      H1[a] = H;
    }
  }

  // ── STEP 5: Pair-based multithreaded GC + TE ──────────────────────────────
  // Build upper-triangle pair list (excluding diagonal — self-MI is trivial)
  std::vector<std::pair<size_t, size_t>> pairs;
  pairs.reserve(num_ca * (num_ca - 1) / 2);
  for (size_t i = 0; i < num_ca; ++i)
    for (size_t j = i + 1; j < num_ca; ++j)
      pairs.emplace_back(i, j);

  size_t n_pairs = pairs.size();
  // Result arrays
  std::vector<double> gc_vals(n_pairs, 0.0);
  std::vector<double> net_te_vals(n_pairs, 0.0);

  auto worker = [&](size_t p_start, size_t p_end) {
    // Thread-local histogram buffers — allocated once, reused per pair
    std::vector<double> hist2d(static_cast<size_t>(B) * B, 0.0);
    std::vector<double> hist2d_b(static_cast<size_t>(B) * B, 0.0); // for TE
    std::vector<double> hist3d(static_cast<size_t>(B) * B * B, 0.0);
    std::vector<double> hist1d_te(B, 0.0);

    for (size_t p = p_start; p < p_end; ++p) {
      size_t i = pairs[p].first;
      size_t j = pairs[p].second;

      // ── GC: Mutual Information ────────────────────────────────────────
      // Joint 2D histogram H(X_i, X_j)
      std::fill(hist2d.begin(), hist2d.end(), 0.5); // Laplace
      for (size_t f = 0; f < total_frames; ++f) {
        int bi = static_cast<int>(proj[i][f] * B);
        int bj = static_cast<int>(proj[j][f] * B);
        if (bi >= B)
          bi = B - 1;
        if (bj >= B)
          bj = B - 1;
        hist2d[bi * B + bj] += 1.0;
      }
      double norm2 = std::accumulate(hist2d.begin(), hist2d.end(), 0.0);
      double H_ij = 0.0;
      for (auto c : hist2d) {
        double p = c / norm2;
        if (p > 0.0)
          H_ij -= p * std::log(p);
      }

      double MI = H1[i] + H1[j] - H_ij;
      if (MI < 0.0)
        MI = 0.0; // numerical noise guard
      double GC = std::sqrt(1.0 - std::exp(-2.0 * MI));
      gc_vals[p] = GC;

      // ── Transfer Entropy i→j and j→i ──────────────────────────────────
      // TE(i→j) = H(Yₜ,Yₜ₋τ) + H(Yₜ₋τ,Xₜ₋τ) - H(Yₜ₋τ) - H(Yₜ,Yₜ₋τ,Xₜ₋τ)
      auto calc_te = [&](size_t src, size_t dst) -> double {
        // H(Yₜ₋τ) — 1D marginal
        std::fill(hist1d_te.begin(), hist1d_te.end(), 0.5);
        for (size_t f = tau; f < total_frames; ++f) {
          int b = static_cast<int>(proj[dst][f - tau] * B);
          if (b >= B)
            b = B - 1;
          hist1d_te[b] += 1.0;
        }
        double norm1 = std::accumulate(hist1d_te.begin(), hist1d_te.end(), 0.0);
        double H_Yt_tau = 0.0;
        for (int b = 0; b < B; ++b) {
          double pp = hist1d_te[b] / norm1;
          if (pp > 0.0)
            H_Yt_tau -= pp * std::log(pp);
        }

        // H(Yₜ, Yₜ₋τ) — 2D
        std::fill(hist2d.begin(), hist2d.end(), 0.5);
        for (size_t f = tau; f < total_frames; ++f) {
          int b1 = static_cast<int>(proj[dst][f] * B);
          int b2 = static_cast<int>(proj[dst][f - tau] * B);
          if (b1 >= B)
            b1 = B - 1;
          if (b2 >= B)
            b2 = B - 1;
          hist2d[b1 * B + b2] += 1.0;
        }
        double n2 = std::accumulate(hist2d.begin(), hist2d.end(), 0.0);
        double H_Yt_Yt_tau = 0.0;
        for (auto c : hist2d) {
          double pp = c / n2;
          if (pp > 0.0)
            H_Yt_Yt_tau -= pp * std::log(pp);
        }

        // H(Yₜ₋τ, Xₜ₋τ) — 2D
        std::fill(hist2d_b.begin(), hist2d_b.end(), 0.5);
        for (size_t f = tau; f < total_frames; ++f) {
          int b1 = static_cast<int>(proj[dst][f - tau] * B);
          int b2 = static_cast<int>(proj[src][f - tau] * B);
          if (b1 >= B)
            b1 = B - 1;
          if (b2 >= B)
            b2 = B - 1;
          hist2d_b[b1 * B + b2] += 1.0;
        }
        double n2b = std::accumulate(hist2d_b.begin(), hist2d_b.end(), 0.0);
        double H_Yt_tau_Xt_tau = 0.0;
        for (auto c : hist2d_b) {
          double pp = c / n2b;
          if (pp > 0.0)
            H_Yt_tau_Xt_tau -= pp * std::log(pp);
        }

        // H(Yₜ, Yₜ₋τ, Xₜ₋τ) — 3D
        std::fill(hist3d.begin(), hist3d.end(), 0.5);
        for (size_t f = tau; f < total_frames; ++f) {
          int b1 = static_cast<int>(proj[dst][f] * B);
          int b2 = static_cast<int>(proj[dst][f - tau] * B);
          int b3 = static_cast<int>(proj[src][f - tau] * B);
          if (b1 >= B)
            b1 = B - 1;
          if (b2 >= B)
            b2 = B - 1;
          if (b3 >= B)
            b3 = B - 1;
          hist3d[(b1 * B + b2) * B + b3] += 1.0;
        }
        double n3 = std::accumulate(hist3d.begin(), hist3d.end(), 0.0);
        double H_Yt_Yt_tau_Xt_tau = 0.0;
        for (auto c : hist3d) {
          double pp = c / n3;
          if (pp > 0.0)
            H_Yt_Yt_tau_Xt_tau -= pp * std::log(pp);
        }

        return H_Yt_Yt_tau + H_Yt_tau_Xt_tau - H_Yt_tau - H_Yt_Yt_tau_Xt_tau;
      };

      double te_ij = calc_te(i, j);
      double te_ji = calc_te(j, i);
      net_te_vals[p] = te_ij - te_ji;
    }
  };

  // Launch threads with balanced pair distribution
  std::vector<std::thread> workers;
  size_t chunk = (n_pairs + num_threads - 1) / num_threads;
  for (size_t t = 0; t < num_threads; ++t) {
    size_t start = t * chunk;
    size_t end = std::min(start + chunk, n_pairs);
    if (start < end)
      workers.emplace_back(worker, start, end);
  }
  for (auto &w : workers) {
    if (w.joinable())
      w.join();
  }

  // ── STEP 6: Export Results ─────────────────────────────────────────────────
  // GC edges
  std::ofstream gc_out("network_gc_edges.csv");
  gc_out << "Source,Target,GC\n";
  for (size_t p = 0; p < n_pairs; ++p) {
    if (gc_vals[p] > 0.35) {
      size_t i = pairs[p].first;
      size_t j = pairs[p].second;
      gc_out << ca_resname[i] << ca_resnrs[i] << "_" << ca_chain[i] << ","
             << ca_resname[j] << ca_resnrs[j] << "_" << ca_chain[j] << ","
             << gc_vals[p] << "\n";
    }
  }
  gc_out.close();

  // TE edges (directed — source is the driver)
  std::ofstream te_out("network_te_edges.csv");
  te_out << "Source,Target,Net_TE\n";
  for (size_t p = 0; p < n_pairs; ++p) {
    if (std::abs(net_te_vals[p]) > 0.05) {
      size_t drv, tgt;
      if (net_te_vals[p] > 0.0) {
        drv = pairs[p].first;
        tgt = pairs[p].second;
      } else {
        drv = pairs[p].second;
        tgt = pairs[p].first;
      }
      te_out << ca_resname[drv] << ca_resnrs[drv] << "_" << ca_chain[drv] << ","
             << ca_resname[tgt] << ca_resnrs[tgt] << "_" << ca_chain[tgt] << ","
             << std::abs(net_te_vals[p]) << "\n";
    }
  }
  te_out.close();
}

// =============================================================================
// LSP-MD: Local Spatial Pattern (Kornev, JCTC 2025, 21, 8699)
// =============================================================================

// ── Brandes weighted betweenness centrality (Dijkstra variant) ────────────────
namespace {
void brandes_weighted_bc(
    size_t n,
    const std::vector<std::vector<std::pair<size_t, double>>> &adj,
    std::vector<double> &bc_out,
    size_t num_threads = 1) {
  bc_out.assign(n, 0.0);

  // Per-thread partial sums (avoids false sharing / contention)
  std::vector<std::vector<double>> partial(num_threads,
                                           std::vector<double>(n, 0.0));

  auto worker = [&](size_t t, size_t s_start, size_t s_end) {
    auto &bc_local = partial[t];
    // Allocate once per thread, reuse across source iterations
    std::vector<double> dist(n), sigma(n), delta(n);
    std::vector<std::vector<size_t>> pred(n);
    std::vector<size_t> order;
    order.reserve(n);

    for (size_t s = s_start; s < s_end; ++s) {
      std::fill(dist.begin(), dist.end(), 1e30);
      std::fill(sigma.begin(), sigma.end(), 0.0);
      std::fill(delta.begin(), delta.end(), 0.0);
      for (auto &p : pred) p.clear(); // clear keeps capacity
      order.clear();

      using PQEntry = std::pair<double, size_t>;
      std::priority_queue<PQEntry, std::vector<PQEntry>,
                          std::greater<PQEntry>> pq;

      dist[s] = 0.0;
      sigma[s] = 1.0;
      pq.push({0.0, s});

      while (!pq.empty()) {
        auto [d_u, u] = pq.top();
        pq.pop();
        if (d_u > dist[u])
          continue;
        order.push_back(u);
        for (auto &[v, w] : adj[u]) {
          double alt = dist[u] + w;
          if (alt < dist[v] - 1e-12) {
            dist[v] = alt;
            sigma[v] = sigma[u];
            pred[v].clear();
            pred[v].push_back(u);
            pq.push({alt, v});
          } else if (std::abs(alt - dist[v]) < 1e-12) {
            sigma[v] += sigma[u];
            pred[v].push_back(u);
          }
        }
      }

      // Back-propagation of dependencies
      for (auto it = order.rbegin(); it != order.rend(); ++it) {
        size_t w = *it;
        for (size_t p : pred[w]) {
          double coeff = (sigma[p] / sigma[w]) * (1.0 + delta[w]);
          delta[p] += coeff;
        }
        if (w != s) {
          bc_local[w] += delta[w];
        }
      }
    }
  };

  // Distribute source nodes across threads
  std::vector<std::thread> threads;
  size_t chunk = (n + num_threads - 1) / num_threads;
  for (size_t t = 0; t < num_threads; ++t) {
    size_t s = t * chunk;
    size_t e = std::min(s + chunk, n);
    if (s < e)
      threads.emplace_back(worker, t, s, e);
  }
  for (auto &th : threads)
    if (th.joinable())
      th.join();

  // Reduce partial sums
  for (size_t t = 0; t < num_threads; ++t)
    for (size_t i = 0; i < n; ++i)
      bc_out[i] += partial[t][i];

  // Normalize for undirected graph
  if (n > 2) {
    double norm = static_cast<double>((n - 1) * (n - 2)) / 2.0;
    for (auto &v : bc_out)
      v /= norm;
  }
}
} // anonymous namespace

// ── AVX2: compute 4 distances in parallel ────────────────────────────────────
namespace {
#ifdef __AVX__
inline void compute_4dist_avx(
    double cax_i, double cay_i, double caz_i,
    double cbx_i, double cby_i, double cbz_i,
    double cax_j, double cay_j, double caz_j,
    double cbx_j, double cby_j, double cbz_j,
    double &d1, double &d2, double &d3, double &d4) {
  // _mm256_set_pd order: arg3→lane3, arg2→lane2, arg1→lane1, arg0→lane0
  // buf[0]=d1(CA-CA), buf[1]=d2(CA-CB), buf[2]=d3(CB-CA), buf[3]=d4(CB-CB)
  __m256d dx = _mm256_set_pd(cbx_i - cbx_j, cbx_i - cax_j,
                              cax_i - cbx_j, cax_i - cax_j);
  __m256d dy = _mm256_set_pd(cby_i - cby_j, cby_i - cay_j,
                              cay_i - cby_j, cay_i - cay_j);
  __m256d dz = _mm256_set_pd(cbz_i - cbz_j, cbz_i - caz_j,
                              caz_i - cbz_j, caz_i - caz_j);
  __m256d r2 = _mm256_fmadd_pd(
      dz, dz, _mm256_fmadd_pd(dy, dy, _mm256_mul_pd(dx, dx)));
  __m256d r = _mm256_sqrt_pd(r2);
  double buf[4];
  _mm256_storeu_pd(buf, r); // unaligned store — safe on stack
  d1 = buf[0]; d2 = buf[1]; d3 = buf[2]; d4 = buf[3];
}
#endif

inline void compute_4dist_scalar(
    double cax_i, double cay_i, double caz_i,
    double cbx_i, double cby_i, double cbz_i,
    double cax_j, double cay_j, double caz_j,
    double cbx_j, double cby_j, double cbz_j,
    double &d1, double &d2, double &d3, double &d4) {
  double dx, dy, dz;
  dx = cax_i - cax_j; dy = cay_i - cay_j; dz = caz_i - caz_j;
  d1 = std::sqrt(dx*dx + dy*dy + dz*dz);
  dx = cax_i - cbx_j; dy = cay_i - cby_j; dz = caz_i - cbz_j;
  d2 = std::sqrt(dx*dx + dy*dy + dz*dz);
  dx = cbx_i - cax_j; dy = cby_i - cay_j; dz = cbz_i - caz_j;
  d3 = std::sqrt(dx*dx + dy*dy + dz*dz);
  dx = cbx_i - cbx_j; dy = cby_i - cby_j; dz = cbz_i - cbz_j;
  d4 = std::sqrt(dx*dx + dy*dy + dz*dz);
}
} // anonymous namespace

void DCCMNetwork::calculate_lspmd(size_t num_threads,
                                   const LSPMDConfig &cfg,
                                   LSPMDResult &result) {
  if (total_frames == 0 || ca_indices.empty())
    return;

  if (total_frames < 10) {
    std::cerr << "[LSP-MD] WARNING: Too few frames (" << total_frames
              << " < 10), results may be unreliable.\n";
  }

  size_t num_ca = ca_indices.size();
  double cutoff = cfg.cutoff_ang;

  // ── STEP 0: Kabsch alignment (guarded — skip if already aligned) ───────────
  if (!kabsch_aligned_) {
    Eigen::MatrixXd ref(num_ca, 3);
    for (size_t a = 0; a < num_ca; ++a) {
      ref(a, 0) = traj_x[0][a];
      ref(a, 1) = traj_y[0][a];
      ref(a, 2) = traj_z[0][a];
    }
    ref.rowwise() -= ref.colwise().mean().eval();

    Eigen::MatrixXd mobile(num_ca, 3);
    Eigen::Matrix3d R;
    for (size_t f = 0; f < total_frames; ++f) {
      for (size_t a = 0; a < num_ca; ++a) {
        mobile(a, 0) = traj_x[f][a];
        mobile(a, 1) = traj_y[f][a];
        mobile(a, 2) = traj_z[f][a];
      }
      Eigen::RowVector3d ca_centroid = mobile.colwise().mean();
      TitanAlignment::align_frame_inplace(ref, mobile, R);
      for (size_t a = 0; a < num_ca; ++a) {
        traj_x[f][a] = mobile(a, 0);
        traj_y[f][a] = mobile(a, 1);
        traj_z[f][a] = mobile(a, 2);
      }
      if (!traj_cbx.empty()) {
        for (size_t a = 0; a < num_ca; ++a) {
          Eigen::RowVector3d cb(traj_cbx[f][a], traj_cby[f][a],
                                traj_cbz[f][a]);
          cb -= ca_centroid;
          cb = (R * cb.transpose()).transpose();
          traj_cbx[f][a] = cb(0);
          traj_cby[f][a] = cb(1);
          traj_cbz[f][a] = cb(2);
        }
      }
    }
    kabsch_aligned_ = true;
  }

  // ── Pre-transpose to atom-major layout [atom][frame] for cache locality ────
  size_t F = total_frames;
  std::vector<std::vector<double>> ax(num_ca, std::vector<double>(F));
  std::vector<std::vector<double>> ay(num_ca, std::vector<double>(F));
  std::vector<std::vector<double>> az(num_ca, std::vector<double>(F));
  std::vector<std::vector<double>> bx(num_ca, std::vector<double>(F));
  std::vector<std::vector<double>> by(num_ca, std::vector<double>(F));
  std::vector<std::vector<double>> bz(num_ca, std::vector<double>(F));

  for (size_t f = 0; f < F; ++f)
    for (size_t a = 0; a < num_ca; ++a) {
      ax[a][f] = traj_x[f][a]; ay[a][f] = traj_y[f][a]; az[a][f] = traj_z[f][a];
      bx[a][f] = traj_cbx[f][a]; by[a][f] = traj_cby[f][a]; bz[a][f] = traj_cbz[f][a];
    }

  // ── STEP 1: Build upper-triangle pair list and compute W_ij ────────────────
  std::vector<std::pair<size_t, size_t>> pairs;
  pairs.reserve(num_ca * (num_ca - 1) / 2);
  for (size_t i = 0; i < num_ca; ++i)
    for (size_t j = i + 1; j < num_ca; ++j)
      pairs.emplace_back(i, j);

  size_t n_pairs = pairs.size();
  std::vector<double> W_vals(n_pairs,
      std::numeric_limits<double>::quiet_NaN());

  auto worker = [&](size_t p_start, size_t p_end) {
    for (size_t p = p_start; p < p_end; ++p) {
      size_t i = pairs[p].first;
      size_t j = pairs[p].second;

      // Pointers into transposed arrays — sequential access in inner loop
      const double *axi = ax[i].data(), *ayi = ay[i].data(), *azi = az[i].data();
      const double *bxi = bx[i].data(), *byi = by[i].data(), *bzi = bz[i].data();
      const double *axj = ax[j].data(), *ayj = ay[j].data(), *azj = az[j].data();
      const double *bxj = bx[j].data(), *byj = by[j].data(), *bzj = bz[j].data();

      // ── Single pass: compute distances + cache for deviation ──
      std::vector<double> d1v(F), d2v(F), d3v(F), d4v(F);
      double sum_d1 = 0.0, sum_d2 = 0.0, sum_d3 = 0.0, sum_d4 = 0.0;
      for (size_t f = 0; f < F; ++f) {
#ifdef __AVX__
        compute_4dist_avx(axi[f], ayi[f], azi[f], bxi[f], byi[f], bzi[f],
                          axj[f], ayj[f], azj[f], bxj[f], byj[f], bzj[f],
                          d1v[f], d2v[f], d3v[f], d4v[f]);
#else
        compute_4dist_scalar(axi[f], ayi[f], azi[f], bxi[f], byi[f], bzi[f],
                             axj[f], ayj[f], azj[f], bxj[f], byj[f], bzj[f],
                             d1v[f], d2v[f], d3v[f], d4v[f]);
#endif
        sum_d1 += d1v[f]; sum_d2 += d2v[f];
        sum_d3 += d3v[f]; sum_d4 += d4v[f];
      }

      double inv_f = 1.0 / static_cast<double>(F);
      double avg_d1 = sum_d1 * inv_f;
      double avg_d2 = sum_d2 * inv_f;
      double avg_d3 = sum_d3 * inv_f;
      double avg_d4 = sum_d4 * inv_f;

      // ── Cutoff check ──
      if (cfg.use_avg_cutoff) {
        if (avg_d1 > cutoff)
          continue;
      } else {
        double dx0 = ax[i][0] - ax[j][0];
        double dy0 = ay[i][0] - ay[j][0];
        double dz0 = az[i][0] - az[j][0];
        if (std::sqrt(dx0 * dx0 + dy0 * dy0 + dz0 * dz0) > cutoff)
          continue;
      }

      // ── Reuse cached distances for deviation ──
      double sum_Z2 = 0.0;
      for (size_t f = 0; f < F; ++f) {
        double dd1 = d1v[f] - avg_d1, dd2 = d2v[f] - avg_d2;
        double dd3 = d3v[f] - avg_d3, dd4 = d4v[f] - avg_d4;
        sum_Z2 += dd1*dd1 + dd2*dd2 + dd3*dd3 + dd4*dd4;
      }

      double Z_avg = std::sqrt(sum_Z2 * inv_f);
      W_vals[p] = std::exp(-Z_avg);
    }
  };

  // ── Launch multithreaded pair computation ──
  {
    std::vector<std::thread> workers;
    size_t chunk = (n_pairs + num_threads - 1) / num_threads;
    for (size_t t = 0; t < num_threads; ++t) {
      size_t start = t * chunk;
      size_t end = std::min(start + chunk, n_pairs);
      if (start < end)
        workers.emplace_back(worker, start, end);
    }
    for (auto &w : workers)
      if (w.joinable())
        w.join();
  }

  // ── Collect valid edges (NaN = skipped) ──
  result.edges.clear();
  for (size_t p = 0; p < n_pairs; ++p) {
    if (!std::isnan(W_vals[p])) {
      result.edges.emplace_back(static_cast<int>(pairs[p].first),
                                static_cast<int>(pairs[p].second),
                                W_vals[p]);
    }
  }

  // ── STEP 2: Degree Centrality ──────────────────────────────────────────────
  result.dc.assign(num_ca, 0.0);
  for (auto &[i, j, w] : result.edges) {
    result.dc[i] += w;
    result.dc[j] += w;
  }

  // ── STEP 3: Betweenness Centrality (Brandes + Dijkstra) ────────────────────
  // Build adjacency list with D_ij = -ln(W_ij) as distance
  std::vector<std::vector<std::pair<size_t, double>>> adj(num_ca);
  for (auto &[i, j, w] : result.edges) {
    double d_weight = (w > 1e-15) ? -std::log(w) : 35.0;
    adj[i].emplace_back(static_cast<size_t>(j), d_weight);
    adj[j].emplace_back(static_cast<size_t>(i), d_weight);
  }

  brandes_weighted_bc(num_ca, adj, result.bc, num_threads);
}

void DCCMNetwork::export_lspmd(const LSPMDResult &result,
                                double weight_threshold) {
  const auto &resname = ca_resname;
  const auto &resnrs = ca_resnrs;
  const auto &chains = ca_chain;
  size_t N = resname.size();

  // --- 1. Edges ---
  {
    std::ofstream ofs("lspmd_edges.csv");
    ofs << "Source,Target,Z_avg,Weight,D_ij\n";
    ofs << std::fixed << std::setprecision(6);
    for (const auto &[i, j, W] : result.edges) {
      if (W < weight_threshold)
        continue;
      double Z = -std::log(W > 1e-15 ? W : 1e-15);
      double D = Z; // D_ij = -ln(W_ij) = <Z>
      ofs << resname[i] << resnrs[i] << "_" << chains[i] << ","
          << resname[j] << resnrs[j] << "_" << chains[j] << "," << Z << ","
          << W << "," << D << "\n";
    }
  }

  // --- 2. Degree Centrality ---
  {
    double dc_max = *std::max_element(result.dc.begin(), result.dc.end());
    if (dc_max < 1e-15)
      dc_max = 1.0;
    std::ofstream ofs("lspmd_degree_centrality.csv");
    ofs << "Residue,ResNr,Chain,IsLig,DC,DC_norm\n";
    ofs << std::fixed << std::setprecision(6);
    for (size_t i = 0; i < N; ++i) {
      ofs << resname[i] << "," << resnrs[i] << "," << chains[i] << ","
          << ca_is_lig[i] << "," << result.dc[i] << ","
          << result.dc[i] / dc_max << "\n";
    }
  }

  // --- 3. Betweenness Centrality ---
  {
    double bc_max = *std::max_element(result.bc.begin(), result.bc.end());
    if (bc_max < 1e-15)
      bc_max = 1.0;
    std::ofstream ofs("lspmd_betweenness_centrality.csv");
    ofs << "Residue,ResNr,Chain,IsLig,BC,BC_norm\n";
    ofs << std::fixed << std::setprecision(6);
    for (size_t i = 0; i < N; ++i) {
      ofs << resname[i] << "," << resnrs[i] << "," << chains[i] << ","
          << ca_is_lig[i] << "," << result.bc[i] << ","
          << result.bc[i] / bc_max << "\n";
    }
  }

  TitanUI::print_info("[LSP-MD] Exported: lspmd_edges.csv, "
                      "lspmd_degree_centrality.csv, "
                      "lspmd_betweenness_centrality.csv");
}
