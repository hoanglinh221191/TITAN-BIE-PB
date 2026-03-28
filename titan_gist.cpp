// titan_gist.cpp — GIST Water Thermodynamic Analysis
// Grid Inhomogeneous Solvation Theory — density, enthalpy, entropy maps
// Overhaul: merged solute arrays, AVX2 E_sw kernel, jthread parallelism
#include "titan_gist.h"
#include "titan_ui.h"

#include <algorithm>
#include <cstdint>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <thread>

#ifdef __AVX__
#include "functions.h" // hsum_double_avx
#include "simd_compat.h"
#include <immintrin.h>
#endif

namespace {

struct LocalGistVoxelAccum {
  uint32_t count{0};
  double energy_sw{0.0};
  double dipole_x{0.0};
  double dipole_y{0.0};
  double dipole_z{0.0};
};

struct LocalGistAccum {
  std::unordered_map<size_t, size_t> voxel_to_slot;
  std::vector<size_t> touched_voxels;
  std::vector<LocalGistVoxelAccum> values;

  void reserve(size_t expected) {
    touched_voxels.reserve(expected);
    values.reserve(expected);
    voxel_to_slot.reserve(expected * 2 + 1);
  }

  void add(size_t voxel_idx, double dipole_x, double dipole_y,
           double dipole_z, double energy_sw) {
    auto [it, inserted] =
        voxel_to_slot.try_emplace(voxel_idx, touched_voxels.size());
    if (inserted) {
      touched_voxels.push_back(voxel_idx);
      values.emplace_back();
    }

    auto &slot = values[it->second];
    ++slot.count;
    slot.energy_sw += energy_sw;
    slot.dipole_x += dipole_x;
    slot.dipole_y += dipole_y;
    slot.dipole_z += dipole_z;
  }
};

void merge_local_gist_accum(std::vector<GistVoxel> &grid,
                            const LocalGistAccum &accum) {
  for (size_t i = 0; i < accum.touched_voxels.size(); ++i) {
    const size_t voxel_idx = accum.touched_voxels[i];
    const auto &slot = accum.values[i];
    auto &dst = grid[voxel_idx];
    dst.count += slot.count;
    dst.energy_sw += slot.energy_sw;
    dst.dipole_x += slot.dipole_x;
    dst.dipole_y += slot.dipole_y;
    dst.dipole_z += slot.dipole_z;
  }
}

} // namespace

// ============================================================================
// Water model constants (TIP3P / SPC/E)
// ============================================================================
struct WaterModelParams {
  double q_O, q_H, sigma_O, eps_O, rho_bulk, E_ww_bulk;
};

static WaterModelParams get_water_params(const std::string &model) {
  if (model == "TIP3P" || model == "tip3p")
    return {-0.8340, +0.4170, 3.1507, 0.6360, 0.03340, -11.60};
  if (model == "SPC" || model == "SPCE" || model == "spc" || model == "spce")
    return {-0.8476, +0.4238, 3.1660, 0.6500, 0.03316, -12.10};
  TitanUI::print_warning("[GIST] Unknown water model '" + model +
                         "', defaulting to TIP3P");
  return {-0.8340, +0.4170, 3.1507, 0.6360, 0.03340, -11.60};
}

// ============================================================================
// Kabsch alignment (minimal — used when TICA is off)
// ============================================================================
static Eigen::Matrix3d kabsch_rotation(const Eigen::MatrixXd &ref_centered,
                                       const Eigen::MatrixXd &mobile_centered) {
  Eigen::Matrix3d H = mobile_centered.transpose() * ref_centered;
  Eigen::JacobiSVD<Eigen::Matrix3d> svd(H, Eigen::ComputeFullU |
                                               Eigen::ComputeFullV);
  Eigen::Matrix3d V = svd.matrixV();
  Eigen::Matrix3d U = svd.matrixU();
  double d = (V * U.transpose()).determinant();
  Eigen::Matrix3d D = Eigen::Matrix3d::Identity();
  D(2, 2) = (d < 0) ? -1.0 : 1.0;
  return V * D * U.transpose();
}

// ============================================================================
// init — build grid from ligand AABB + padding, cache merged solute params
// ============================================================================
void TitanGIST::init(const MOL &lig_ref, double padding, double resolution,
                     double temperature, double cutoff, size_t n_threads,
                     const std::string &water_model,
                     const Multi_MOL_parameters &rec_params,
                     const std::vector<mol_info> &rec_info,
                     const Multi_MOL_parameters &lig_params,
                     const std::vector<mol_info> &lig_info) {
  voxel_size_ = resolution;
  temperature_ = temperature;
  cut2_ = cutoff * cutoff;
  nt_ = (n_threads > 0) ? n_threads : 1;

  // Water model
  auto wm = get_water_params(water_model);
  q_O_ = wm.q_O;
  q_H_ = wm.q_H;
  sigma_O_ = wm.sigma_O;
  eps_O_ = wm.eps_O;
  rho_bulk_ = wm.rho_bulk;
  E_ww_bulk_ = wm.E_ww_bulk;

  // Bounding box from ligand (frame 0)
  const auto &lx = lig_ref.get_pos_x();
  const auto &ly = lig_ref.get_pos_y();
  const auto &lz = lig_ref.get_pos_z();
  if (lx.empty()) {
    TitanUI::print_warning("[GIST] No ligand atoms — cannot build grid");
    return;
  }

  double xmin = 1e9, ymin = 1e9, zmin = 1e9;
  double xmax = -1e9, ymax = -1e9, zmax = -1e9;
  for (size_t i = 0; i < lx.size(); ++i) {
    xmin = std::min(xmin, lx[i]);
    xmax = std::max(xmax, lx[i]);
    ymin = std::min(ymin, ly[i]);
    ymax = std::max(ymax, ly[i]);
    zmin = std::min(zmin, lz[i]);
    zmax = std::max(zmax, lz[i]);
  }

  origin_ = {xmin - padding, ymin - padding, zmin - padding};
  nx_ = std::max(
      1, static_cast<int>(std::ceil((xmax - xmin + 2 * padding) / resolution)));
  ny_ = std::max(
      1, static_cast<int>(std::ceil((ymax - ymin + 2 * padding) / resolution)));
  nz_ = std::max(
      1, static_cast<int>(std::ceil((zmax - zmin + 2 * padding) / resolution)));

  size_t total = static_cast<size_t>(nx_) * ny_ * nz_;
  if (total > 5000000) {
    TitanUI::print_warning("[GIST] Grid too large (" + std::to_string(total) +
                           " voxels). Consider increasing --gist-resolution.");
  }
  grid_.resize(total);
  TitanUI::print_info("[GIST] Grid: " + std::to_string(nx_) + "x" +
                      std::to_string(ny_) + "x" + std::to_string(nz_) + " = " +
                      std::to_string(total) + " voxels (" +
                      std::to_string(resolution) +
                      " Å), cutoff=" + std::to_string(cutoff) + " Å");

  // --- Collect heavy atom indices for rec and lig ---
  rec_heavy_idx_.clear();
  for (size_t i = 0; i < rec_info.size(); ++i) {
    if (rec_info[i].get_atom()[0] != 'H')
      rec_heavy_idx_.push_back(i);
  }
  lig_heavy_idx_.clear();
  for (size_t i = 0; i < lig_info.size(); ++i) {
    if (lig_info[i].get_atom()[0] != 'H')
      lig_heavy_idx_.push_back(i);
  }

  n_rec_heavy_ = rec_heavy_idx_.size();
  n_solute_real_ = n_rec_heavy_ + lig_heavy_idx_.size();
  // Pad to multiple of 4 for AVX2
  n_solute_ = ((n_solute_real_ + 3) / 4) * 4;

  // Allocate merged arrays (padded with zeros)
  solute_x_.assign(n_solute_, 0.0);
  solute_y_.assign(n_solute_, 0.0);
  solute_z_.assign(n_solute_, 0.0);
  solute_charge_.assign(n_solute_, 0.0);
  solute_sigma_half_.assign(n_solute_, 0.0);
  solute_eps_mixed_.assign(n_solute_, 0.0);

  // Fill receptor portion [0, n_rec_heavy_)
  for (size_t j = 0; j < n_rec_heavy_; ++j) {
    size_t i = rec_heavy_idx_[j];
    if (i < rec_params.charge.size())
      solute_charge_[j] = rec_params.charge[i];
    double sigma_i = 0.0, eps_i = 0.0;
    if (i < rec_params.sigma.size())
      sigma_i = rec_params.sigma[i];
    if (i < rec_params.epsilon.size())
      eps_i = rec_params.epsilon[i];
    // Precompute mixing rules
    solute_sigma_half_[j] = 0.5 * (sigma_O_ + sigma_i);
    solute_eps_mixed_[j] = std::sqrt(eps_O_ * eps_i);
  }

  // Fill ligand portion [n_rec_heavy_, n_rec_heavy_ + lig_heavy_idx_.size())
  for (size_t j = 0; j < lig_heavy_idx_.size(); ++j) {
    size_t idx = n_rec_heavy_ + j;
    size_t i = lig_heavy_idx_[j];
    if (i < lig_params.charge.size())
      solute_charge_[idx] = lig_params.charge[i];
    double sigma_i = 0.0, eps_i = 0.0;
    if (i < lig_params.sigma.size())
      sigma_i = lig_params.sigma[i];
    if (i < lig_params.epsilon.size())
      eps_i = lig_params.epsilon[i];
    solute_sigma_half_[idx] = 0.5 * (sigma_O_ + sigma_i);
    solute_eps_mixed_[idx] = std::sqrt(eps_O_ * eps_i);
  }

  // Detect Cα indices for self-alignment fallback
  for (size_t i = 0; i < rec_info.size(); ++i) {
    const std::string &aname = rec_info[i].get_atom();
    if (aname == "CA" || aname == "BB" || aname == "BB1")
      ca_indices_.push_back(i);
  }

  TitanUI::print_info("[GIST] Solute heavy: " + std::to_string(n_solute_real_) +
                      " (rec=" + std::to_string(n_rec_heavy_) +
                      " lig=" + std::to_string(lig_heavy_idx_.size()) +
                      "), padded to " + std::to_string(n_solute_) +
                      ", Cα: " + std::to_string(ca_indices_.size()) +
                      ", threads: " + std::to_string(nt_));
}

// ============================================================================
// voxel_idx — returns flat index or -1 if outside grid
// ============================================================================
inline int TitanGIST::voxel_idx(double x, double y, double z) const noexcept {
  int ix = static_cast<int>((x - origin_[0]) / voxel_size_);
  int iy = static_cast<int>((y - origin_[1]) / voxel_size_);
  int iz = static_cast<int>((z - origin_[2]) / voxel_size_);
  if (ix < 0 || ix >= nx_ || iy < 0 || iy >= ny_ || iz < 0 || iz >= nz_)
    return -1;
  return ix * ny_ * nz_ + iy * nz_ + iz;
}

// ============================================================================
// compute_esw — Solute-Water energy for one water molecule (O, H1, H2)
// Coulomb (O + H1 + H2) + LJ (O only) on merged solute array
// Uses AVX2 with cutoff mask + NaN guard, scalar tail for remainder
// ============================================================================
double TitanGIST::compute_esw(double ox, double oy, double oz, double h1x,
                              double h1y, double h1z, double h2x, double h2y,
                              double h2z) const {
  constexpr double K_E = 138.93546; // kJ·Å/(mol·e²), Coulomb constant
  constexpr double MIN_D2 = 0.01;   // 0.1 Å², NaN guard
  double e = 0.0;

#ifdef __AVX__
  // Broadcast water atom coords
  const __m256d v_ox = _mm256_set1_pd(ox);
  const __m256d v_oy = _mm256_set1_pd(oy);
  const __m256d v_oz = _mm256_set1_pd(oz);
  const __m256d v_h1x = _mm256_set1_pd(h1x);
  const __m256d v_h1y = _mm256_set1_pd(h1y);
  const __m256d v_h1z = _mm256_set1_pd(h1z);
  const __m256d v_h2x = _mm256_set1_pd(h2x);
  const __m256d v_h2y = _mm256_set1_pd(h2y);
  const __m256d v_h2z = _mm256_set1_pd(h2z);

  // Broadcast constants
  const __m256d v_qO_KE = _mm256_set1_pd(K_E * q_O_);
  const __m256d v_qH_KE = _mm256_set1_pd(K_E * q_H_);
  const __m256d v_cut2 = _mm256_set1_pd(cut2_);
  const __m256d v_min_d2 = _mm256_set1_pd(MIN_D2);
  const __m256d v_four = _mm256_set1_pd(4.0);

  __m256d v_e_acc = _mm256_setzero_pd();

  // AVX2 loop over solute atoms, 4 at a time (array is padded)
  const size_t n4 = (n_solute_ / 4) * 4;
  for (size_t pi = 0; pi < n4; pi += 4) {
    // Load solute coords (used for O, H1, H2 distances)
    const __m256d vsx = _mm256_loadu_pd(&solute_x_[pi]);
    const __m256d vsy = _mm256_loadu_pd(&solute_y_[pi]);
    const __m256d vsz = _mm256_loadu_pd(&solute_z_[pi]);

    // d²(O, solute)
    const __m256d dx_O = _mm256_sub_pd(v_ox, vsx);
    const __m256d dy_O = _mm256_sub_pd(v_oy, vsy);
    const __m256d dz_O = _mm256_sub_pd(v_oz, vsz);
    const __m256d d2_O = _mm256_add_pd(
        _mm256_add_pd(_mm256_mul_pd(dx_O, dx_O), _mm256_mul_pd(dy_O, dy_O)),
        _mm256_mul_pd(dz_O, dz_O));

    // Cutoff + NaN mask on O-solute distance
    const __m256d cut_mask = _mm256_cmp_pd(d2_O, v_cut2, _CMP_LT_OQ);
    const __m256d valid_mask = _mm256_cmp_pd(d2_O, v_min_d2, _CMP_GT_OQ);
    const __m256d mask = _mm256_and_pd(cut_mask, valid_mask);

    // Early skip if all 4 lanes are outside cutoff
    if (_mm256_movemask_pd(mask) == 0)
      continue;

    const __m256d d_O = _mm256_sqrt_pd(d2_O);

    // Load solute charge, precomputed sigma_half, eps_mixed
    const __m256d vq = _mm256_loadu_pd(&solute_charge_[pi]);
    const __m256d v_sighf = _mm256_loadu_pd(&solute_sigma_half_[pi]);
    const __m256d v_epsmx = _mm256_loadu_pd(&solute_eps_mixed_[pi]);

    // --- Coulomb O ---
    // Use safe reciprocal: mask ensures d_O > sqrt(MIN_D2)
    __m256d safe_d_O = _mm256_blendv_pd(_mm256_set1_pd(1.0), d_O, mask);
    __m256d e_coul_O = _mm256_div_pd(_mm256_mul_pd(v_qO_KE, vq), safe_d_O);
    e_coul_O = _mm256_and_pd(e_coul_O, mask);

    // --- Coulomb H1 ---
    // H1/H2 Coulomb gated by O-solute cutoff (H near O ⇒ also near solute).
    // This is a deliberate approximation consistent with cpptraj/GIST.
    const __m256d dx_H1 = _mm256_sub_pd(v_h1x, vsx);
    const __m256d dy_H1 = _mm256_sub_pd(v_h1y, vsy);
    const __m256d dz_H1 = _mm256_sub_pd(v_h1z, vsz);
    const __m256d d2_H1 = _mm256_add_pd(
        _mm256_add_pd(_mm256_mul_pd(dx_H1, dx_H1), _mm256_mul_pd(dy_H1, dy_H1)),
        _mm256_mul_pd(dz_H1, dz_H1));
    const __m256d d_H1 = _mm256_sqrt_pd(
        _mm256_blendv_pd(_mm256_set1_pd(1.0), d2_H1,
                         _mm256_cmp_pd(d2_H1, v_min_d2, _CMP_GT_OQ)));
    __m256d e_coul_H1 = _mm256_div_pd(_mm256_mul_pd(v_qH_KE, vq), d_H1);
    e_coul_H1 = _mm256_and_pd(e_coul_H1, mask);

    // --- Coulomb H2 ---
    const __m256d dx_H2 = _mm256_sub_pd(v_h2x, vsx);
    const __m256d dy_H2 = _mm256_sub_pd(v_h2y, vsy);
    const __m256d dz_H2 = _mm256_sub_pd(v_h2z, vsz);
    const __m256d d2_H2 = _mm256_add_pd(
        _mm256_add_pd(_mm256_mul_pd(dx_H2, dx_H2), _mm256_mul_pd(dy_H2, dy_H2)),
        _mm256_mul_pd(dz_H2, dz_H2));
    const __m256d d_H2 = _mm256_sqrt_pd(
        _mm256_blendv_pd(_mm256_set1_pd(1.0), d2_H2,
                         _mm256_cmp_pd(d2_H2, v_min_d2, _CMP_GT_OQ)));
    __m256d e_coul_H2 = _mm256_div_pd(_mm256_mul_pd(v_qH_KE, vq), d_H2);
    e_coul_H2 = _mm256_and_pd(e_coul_H2, mask);

    // --- LJ O-solute (Lorentz-Berthelot) ---
    // eps_ij = eps_mixed² (since eps_mixed = sqrt(eps_O * eps_i))
    const __m256d eps_ij = _mm256_mul_pd(v_epsmx, v_epsmx);
    // sr = sigma_half / d_O
    const __m256d sr = _mm256_div_pd(v_sighf, safe_d_O);
    const __m256d sr2 = _mm256_mul_pd(sr, sr);
    const __m256d sr4 = _mm256_mul_pd(sr2, sr2);
    const __m256d sr6 = _mm256_mul_pd(sr4, sr2);
    // e_LJ = 4 * eps_ij * (sr6² - sr6)  [FMA: fused mul-sub]
    const __m256d sr_term = simd_compat::fmsub_pd(sr6, sr6, sr6);
    __m256d e_lj = _mm256_mul_pd(v_four, _mm256_mul_pd(eps_ij, sr_term));
    e_lj = _mm256_and_pd(e_lj, mask);

    // Accumulate
    v_e_acc =
        _mm256_add_pd(v_e_acc, _mm256_add_pd(_mm256_add_pd(e_coul_O, e_coul_H1),
                                             _mm256_add_pd(e_coul_H2, e_lj)));
  }

  e += hsum_double_avx(v_e_acc);

  // Scalar tail: should not be needed since n_solute_ is padded to multiple of
  // 4, but guard for safety. Pad atoms have charge=0 and eps=0, contributing 0
  // energy.

#else
  // Pure scalar fallback (no AVX)
  for (size_t pi = 0; pi < n_solute_real_; ++pi) {
    double dx = ox - solute_x_[pi], dy = oy - solute_y_[pi],
           dz = oz - solute_z_[pi];
    double d2 = dx * dx + dy * dy + dz * dz;
    if (d2 > cut2_ || d2 < MIN_D2)
      continue;
    double d = std::sqrt(d2);
    double qp = solute_charge_[pi];

    // Coulomb: O-solute
    e += K_E * q_O_ * qp / d;
    // Coulomb: H1-solute
    {
      double ddx = h1x - solute_x_[pi], ddy = h1y - solute_y_[pi],
             ddz = h1z - solute_z_[pi];
      double dd = std::sqrt(ddx * ddx + ddy * ddy + ddz * ddz);
      if (dd > 1e-6)
        e += K_E * q_H_ * qp / dd;
    }
    // Coulomb: H2-solute
    {
      double ddx = h2x - solute_x_[pi], ddy = h2y - solute_y_[pi],
             ddz = h2z - solute_z_[pi];
      double dd = std::sqrt(ddx * ddx + ddy * ddy + ddz * ddz);
      if (dd > 1e-6)
        e += K_E * q_H_ * qp / dd;
    }
    // LJ: O-solute (Lorentz-Berthelot mixing, precomputed)
    double eps_ij = solute_eps_mixed_[pi] * solute_eps_mixed_[pi];
    if (eps_ij > 1e-12) {
      double sr = solute_sigma_half_[pi] / d;
      double sr6 = sr * sr * sr * sr * sr * sr;
      e += 4.0 * eps_ij * (sr6 * sr6 - sr6);
    }
  }
#endif

  return e;
}

// ============================================================================
// process_frame — transform water, bin into grid, accumulate E_sw + dipole
// ============================================================================
void TitanGIST::process_frame(const MOL &pro, const MOL &lig, const MOL &water,
                              const std::vector<mol_info> &pro_info,
                              const Eigen::MatrixXd &ref_ca_centered,
                              const Eigen::Matrix3d &R, bool is_frame_zero,
                              bool tica_available) {
  if (grid_.empty())
    return;

  const auto &wx = water.get_pos_x();
  const auto &wy = water.get_pos_y();
  const auto &wz = water.get_pos_z();
  const size_t n_wat = wx.size();
  if (n_wat < 3)
    return;

  const auto &px = pro.get_pos_x();
  const auto &py = pro.get_pos_y();
  const auto &pz = pro.get_pos_z();

  // --- Step 1: Build or reuse Kabsch alignment ---
  // BUG A fix: Initialize with safe defaults (prevents UB if both blocks skip)
  Eigen::Matrix3d R_actual = Eigen::Matrix3d::Identity();
  Eigen::Vector3d t_mobile = ref_com_;

  if (is_frame_zero && !ref_frame_locked_) {
    // Cache receptor coords in reference frame (frame 0 = reference)
    for (size_t j = 0; j < n_rec_heavy_; ++j) {
      size_t i = rec_heavy_idx_[j];
      solute_x_[j] = px[i];
      solute_y_[j] = py[i];
      solute_z_[j] = pz[i];
    }

    // Cache ligand coords
    const auto &lgx = lig.get_pos_x();
    const auto &lgy = lig.get_pos_y();
    const auto &lgz = lig.get_pos_z();
    for (size_t j = 0; j < lig_heavy_idx_.size(); ++j) {
      size_t idx = n_rec_heavy_ + j;
      size_t i = lig_heavy_idx_[j];
      if (i < lgx.size()) {
        solute_x_[idx] = lgx[i];
        solute_y_[idx] = lgy[i];
        solute_z_[idx] = lgz[i];
      }
    }

    // Build self-alignment reference from Cα (BUG2: prefer tica_ca_indices_)
    const auto &ca_set =
        tica_ca_indices_.empty() ? ca_indices_ : tica_ca_indices_;
    if (!ca_set.empty()) {
      int n_ca = static_cast<int>(ca_set.size());
      self_ref_centered_.resize(n_ca, 3);
      for (int a = 0; a < n_ca; ++a) {
        size_t idx = ca_set[a];
        self_ref_centered_(a, 0) = px[idx];
        self_ref_centered_(a, 1) = py[idx];
        self_ref_centered_(a, 2) = pz[idx];
      }
      ref_com_ = self_ref_centered_.colwise().mean();
      self_ref_centered_.rowwise() -= ref_com_.transpose();
      ref_built_ = true;
    }

    ref_frame_locked_ = true;

    // BUG4: Frame 0 uses identity transform — continue to binning
    R_actual = Eigen::Matrix3d::Identity();
    t_mobile = ref_com_;
  }

  // --- Step 2: Compute rotation for subsequent frames (or multi-file frame 0)
  if (ref_frame_locked_ && !(is_frame_zero && ref_built_)) {
    // BUG2: use tica_ca_indices_ when TICA is on (matches DCCM's Cα subset)
    const auto &ca_rt =
        tica_ca_indices_.empty() ? ca_indices_ : tica_ca_indices_;
    if (tica_available && ref_ca_centered.rows() > 0) {
      R_actual = R;
      t_mobile = Eigen::Vector3d::Zero();
      for (size_t idx : ca_rt) {
        t_mobile[0] += px[idx];
        t_mobile[1] += py[idx];
        t_mobile[2] += pz[idx];
      }
      if (!ca_rt.empty())
        t_mobile /= static_cast<double>(ca_rt.size());
    } else if (ref_built_ && !ca_rt.empty()) {
      int n_ca = static_cast<int>(ca_rt.size());
      Eigen::MatrixXd mobile(n_ca, 3);
      for (int a = 0; a < n_ca; ++a) {
        size_t idx = ca_rt[a];
        mobile(a, 0) = px[idx];
        mobile(a, 1) = py[idx];
        mobile(a, 2) = pz[idx];
      }
      t_mobile = mobile.colwise().mean();
      Eigen::MatrixXd mobile_centered = mobile;
      mobile_centered.rowwise() -= t_mobile.transpose();
      R_actual = kabsch_rotation(self_ref_centered_, mobile_centered);
    } else {
      // BUG B fix: use ref_com_ as identity fallback (not Zero)
      R_actual = Eigen::Matrix3d::Identity();
      t_mobile = ref_com_;
    }
  }

  // Update ligand portion of solute coords in reference frame
  const auto &lgx = lig.get_pos_x();
  const auto &lgy = lig.get_pos_y();
  const auto &lgz = lig.get_pos_z();
  for (size_t j = 0; j < lig_heavy_idx_.size(); ++j) {
    size_t idx = n_rec_heavy_ + j;
    size_t i = lig_heavy_idx_[j];
    if (i < lgx.size()) {
      Eigen::Vector3d p{lgx[i] - t_mobile[0], lgy[i] - t_mobile[1],
                        lgz[i] - t_mobile[2]};
      Eigen::Vector3d q = R_actual.transpose() * p + ref_com_;
      solute_x_[idx] = q[0];
      solute_y_[idx] = q[1];
      solute_z_[idx] = q[2];
    }
  }

  // --- Step 3: Transform water coordinates into reference frame ---
  wat_buf_.resize(n_wat);
  for (size_t i = 0; i < n_wat; ++i) {
    Eigen::Vector3d p{wx[i] - t_mobile[0], wy[i] - t_mobile[1],
                      wz[i] - t_mobile[2]};
    Eigen::Vector3d q = R_actual.transpose() * p + ref_com_;
    wat_buf_[i] = {q[0], q[1], q[2]};
  }

  // --- Step 4: Grid binning + dipole + E_sw (multi-threaded) ---
  const size_t n_water_mols = n_wat / 3; // number of complete water molecules

  if (nt_ <= 1 || n_water_mols < 32) {
    // Single-threaded path
    for (size_t w = 0; w + 2 < n_wat; w += 3) {
      const auto &O = wat_buf_[w];
      const auto &H1 = wat_buf_[w + 1];
      const auto &H2 = wat_buf_[w + 2];

      int vi = voxel_idx(O[0], O[1], O[2]);
      if (vi < 0)
        continue;

      grid_[vi].count++;

      // Dipole: O → midpoint(H1, H2), normalize
      double mx = 0.5 * (H1[0] + H2[0]) - O[0];
      double my = 0.5 * (H1[1] + H2[1]) - O[1];
      double mz = 0.5 * (H1[2] + H2[2]) - O[2];
      double mlen = std::sqrt(mx * mx + my * my + mz * mz);
      if (mlen > 1e-9) {
        mx /= mlen;
        my /= mlen;
        mz /= mlen;
      }
      grid_[vi].dipole_x += mx;
      grid_[vi].dipole_y += my;
      grid_[vi].dipole_z += mz;

      // E_sw
      grid_[vi].energy_sw += compute_esw(O[0], O[1], O[2], H1[0], H1[1], H1[2],
                                         H2[0], H2[1], H2[2]);
    }
  } else {
    // Multi-threaded path: each thread accumulates into a sparse local buffer.
    std::vector<std::jthread> threads;
    threads.reserve(nt_);
    std::vector<LocalGistAccum> local_accums(nt_);

    auto worker = [&](size_t tid) {
      // Divide water molecules among threads
      size_t w_per_thread = n_water_mols / nt_;
      size_t w_start = tid * w_per_thread;
      size_t w_end = (tid == nt_ - 1) ? n_water_mols : (tid + 1) * w_per_thread;
      auto &accum = local_accums[tid];
      accum.reserve(std::max<size_t>(1, w_end - w_start));

      for (size_t wm = w_start; wm < w_end; ++wm) {
        size_t w = wm * 3;
        if (w + 2 >= n_wat)
          break;
        const auto &O = wat_buf_[w];
        const auto &H1 = wat_buf_[w + 1];
        const auto &H2 = wat_buf_[w + 2];

        int vi = voxel_idx(O[0], O[1], O[2]);
        if (vi < 0)
          continue;

        double mx = 0.5 * (H1[0] + H2[0]) - O[0];
        double my = 0.5 * (H1[1] + H2[1]) - O[1];
        double mz = 0.5 * (H1[2] + H2[2]) - O[2];
        double mlen = std::sqrt(mx * mx + my * my + mz * mz);
        if (mlen > 1e-9) {
          mx /= mlen;
          my /= mlen;
          mz /= mlen;
        }
        const double esw = compute_esw(O[0], O[1], O[2], H1[0], H1[1],
                                       H1[2], H2[0], H2[1], H2[2]);
        accum.add(static_cast<size_t>(vi), mx, my, mz, esw);
      }
    };

    for (size_t t = 0; t < nt_; ++t) {
      threads.emplace_back(worker, t);
    }

    // jthread auto-joins at destruction, but we must wait here before reducing
    // thread-local grids.
    for (auto &th : threads) {
      if (th.joinable()) {
        th.join();
      }
    }

    // Serial reduce only the voxels that were actually touched.
    for (const auto &accum : local_accums) {
      merge_local_gist_accum(grid_, accum);
    }
  }
}

// ============================================================================
// write_dx — write one OpenDX file (VMD/PyMOL compatible)
// ============================================================================
void TitanGIST::write_dx(const std::string &path,
                         const std::vector<double> &data,
                         const std::string &quantity_name) const {
  std::ofstream out(path);
  if (!out.is_open()) {
    TitanUI::print_warning("[GIST] Cannot open " + path + " for writing");
    return;
  }

  int ntotal = nx_ * ny_ * nz_;
  out << std::fixed << std::setprecision(6);
  out << "object 1 class gridpositions counts " << nx_ << " " << ny_ << " "
      << nz_ << "\n";
  out << "origin " << origin_[0] << " " << origin_[1] << " " << origin_[2]
      << "\n";
  out << "delta " << voxel_size_ << " 0.000000 0.000000\n";
  out << "delta 0.000000 " << voxel_size_ << " 0.000000\n";
  out << "delta 0.000000 0.000000 " << voxel_size_ << "\n";
  out << "object 2 class gridconnections counts " << nx_ << " " << ny_ << " "
      << nz_ << "\n";
  out << "object 3 class array type double rank 0 items " << ntotal
      << " data follows\n";

  int count = 0;
  for (int ix = 0; ix < nx_; ++ix) {
    for (int iy = 0; iy < ny_; ++iy) {
      for (int iz = 0; iz < nz_; ++iz) {
        int idx = ix * ny_ * nz_ + iy * nz_ + iz;
        out << data[idx];
        ++count;
        out << ((count % 3 == 0) ? '\n' : ' ');
      }
    }
  }
  if (count % 3 != 0)
    out << '\n';

  out << "attribute \"dep\" string \"positions\"\n";
  out << "object \"" << quantity_name << "\" class field\n";
  out.close();
}

// ============================================================================
// export_dx_maps — compute thermodynamics and write 7 .dx files
// ============================================================================
void TitanGIST::export_dx_maps(const std::string &prefix, int n_frames,
                               int min_count) const {
  if (grid_.empty() || n_frames < 1)
    return;

  const int ntotal = nx_ * ny_ * nz_;
  const double V_voxel = voxel_size_ * voxel_size_ * voxel_size_; // ų
  constexpr double R_gas = 8.314462e-3;                           // kJ/(mol·K)

  std::vector<double> density(ntotal, 0.0);
  std::vector<double> dH_sw(ntotal, 0.0);
  std::vector<double> dH_ww(ntotal, 0.0);
  std::vector<double> dH_total(ntotal, 0.0);
  std::vector<double> TdS_trans(ntotal, 0.0);
  std::vector<double> TdS_orient(ntotal, 0.0);
  std::vector<double> dG(ntotal, 0.0);

  int occupied_voxels = 0;
  for (int vi = 0; vi < ntotal; ++vi) {
    uint32_t cnt = grid_[vi].count;
    if (cnt < static_cast<uint32_t>(min_count))
      continue; // All outputs stay 0.0

    ++occupied_voxels;
    double c = static_cast<double>(cnt);

    // g(r) = count / (n_frames * V_voxel * rho_bulk)
    double g = c / (n_frames * V_voxel * rho_bulk_);
    density[vi] = g;

    // ΔH_sw = <E_sw> per visit
    dH_sw[vi] = grid_[vi].energy_sw / c;

    // ΔH_ww = <E_ww>/2 − E_ww_bulk/2  [Phase 1: skip when E_ww=0]
    if (grid_[vi].energy_ww != 0.0) {
      dH_ww[vi] = grid_[vi].energy_ww / c / 2.0 - E_ww_bulk_ / 2.0;
    } else {
      dH_ww[vi] = 0.0; // Phase 1: E_ww not computed
    }

    dH_total[vi] = dH_sw[vi] + dH_ww[vi];

    // −TΔS_trans = R·T·ln(g)
    if (g > 1e-12)
      TdS_trans[vi] = R_gas * temperature_ * std::log(g);

    // −TΔS_orient = R·T·(1 − |<d̂>|²)
    // |<d̂>|² = (Σdx)² + (Σdy)² + (Σdz)² / count²
    // Approximation: -T*dS_orient ≈ RT*(1 - |<d̂>|²)
    // Full formula: Lazaridis 1998 / Nguyen 2012 JCTC 8:2189
    double avg_d2 = (grid_[vi].dipole_x * grid_[vi].dipole_x +
                     grid_[vi].dipole_y * grid_[vi].dipole_y +
                     grid_[vi].dipole_z * grid_[vi].dipole_z) /
                    (c * c);
    TdS_orient[vi] = R_gas * temperature_ * (1.0 - avg_d2);

    dG[vi] = dH_total[vi] + TdS_trans[vi] + TdS_orient[vi];
  }

  TitanUI::print_info(
      "[GIST] Occupied voxels (count >= " + std::to_string(min_count) +
      "): " + std::to_string(occupied_voxels) + " / " + std::to_string(ntotal));

  // Write 7 .dx files
  write_dx(prefix + "_density.dx", density, "density");
  write_dx(prefix + "_dH_sw.dx", dH_sw, "dH_sw");
  write_dx(prefix + "_dH_ww.dx", dH_ww, "dH_ww");
  write_dx(prefix + "_dH_total.dx", dH_total, "dH_total");
  write_dx(prefix + "_TdS_trans.dx", TdS_trans, "TdS_trans");
  write_dx(prefix + "_TdS_orient.dx", TdS_orient, "TdS_orient");
  write_dx(prefix + "_dG.dx", dG, "dG");

  TitanUI::print_info("[GIST] Exported 7 .dx maps → " + prefix + "_*.dx");
}
