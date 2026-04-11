#include "gbsa_intrin.h"
#include "charmm_nonbonded.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>

#ifdef __AVX__
#include "avx_math.h"
#include <immintrin.h>

// Helper to horizontal sum AVX register
static inline double hsum_avx(__m256d v) {
  __m128d vlow = _mm256_castpd256_pd128(v);
  __m128d vhigh = _mm256_extractf128_pd(v, 1);
  vlow = _mm_add_pd(vlow, vhigh);
  __m128d shuf = _mm_shuffle_pd(vlow, vlow, _MM_SHUFFLE2(0, 1));
  vlow = _mm_add_pd(vlow, shuf);
  return _mm_cvtsd_f64(vlow);
}

static inline __m256d charmm_lj_energy_gbsa_avx(__m256d sigma_vec,
                                                __m256d epsilon_vec,
                                                __m256d distance_vec) {
  const __m256d avx_zero = _mm256_setzero_pd();
  const __m256d avx_four = _mm256_set1_pd(4.0);
  const __m256d sigma_valid = _mm256_cmp_pd(sigma_vec, avx_zero, _CMP_GT_OQ);
  const __m256d epsilon_valid =
      _mm256_cmp_pd(epsilon_vec, avx_zero, _CMP_NEQ_OQ);
  const __m256d valid_mask = _mm256_and_pd(sigma_valid, epsilon_valid);

  const __m256d sr = _mm256_div_pd(sigma_vec, distance_vec);
  const __m256d sr2 = _mm256_mul_pd(sr, sr);
  const __m256d sr6 = _mm256_mul_pd(_mm256_mul_pd(sr2, sr2), sr2);
  const __m256d sr12 = _mm256_mul_pd(sr6, sr6);
  __m256d energy = _mm256_mul_pd(
      _mm256_mul_pd(avx_four, epsilon_vec), _mm256_sub_pd(sr12, sr6));
  return _mm256_and_pd(energy, valid_mask);
}

void cal_effective_born_radii_avx2(const MOL &mol,
                                   const AlignedDoubleVector &vdw_radii,
                                   AlignedDoubleVector &radii, int igb,
                                   size_t begin, size_t end, const MOL *env_mol,
                                   const AlignedDoubleVector *env_vdw_radii) {

  size_t size = mol.size();

  double alpha = 0.0, beta = 0.0, gamma = 0.0;
  if (igb == 2) {
    alpha = 0.8;
    beta = 0.0;
    gamma = 2.91;
  } else if (igb == 5) {
    alpha = 1.0;
    beta = 0.8;
    gamma = 4.85;
  }

  __m256d avx_one = _mm256_set1_pd(1.0);

  for (size_t i = begin; i < end; ++i) {
    double pos_x = mol.pos.x[i], pos_y = mol.pos.y[i], pos_z = mol.pos.z[i];
    double rho_i = vdw_radii[i]; // In Angstroms!
    if (rho_i < 1.0)
      rho_i = 1.5; // fallback
    __m256d v_pos_x = _mm256_set1_pd(pos_x);
    __m256d v_pos_y = _mm256_set1_pd(pos_y);
    __m256d v_pos_z = _mm256_set1_pd(pos_z);

    __m256d psi_sum = _mm256_setzero_pd();
    size_t size_vec = size / 4;
    for (size_t j = 0; j < size_vec * 4; j += 4) {
      __m256d v_x = _mm256_load_pd(&mol.pos.x[j]);
      __m256d v_y = _mm256_load_pd(&mol.pos.y[j]);
      __m256d v_z = _mm256_load_pd(&mol.pos.z[j]);
      __m256d v_rho_j = _mm256_load_pd(&vdw_radii[j]);

      __m256d dx = _mm256_sub_pd(v_x, v_pos_x);
      __m256d dy = _mm256_sub_pd(v_y, v_pos_y);
      __m256d dz = _mm256_sub_pd(v_z, v_pos_z);

      __m256d dx2 = _mm256_mul_pd(dx, dx);
      __m256d dy2 = _mm256_mul_pd(dy, dy);
      __m256d dz2 = _mm256_mul_pd(dz, dz);

      __m256d r2 = _mm256_add_pd(_mm256_add_pd(dx2, dy2), dz2);
      __m256d mask = _mm256_cmp_pd(r2, _mm256_set1_pd(1e-6), _CMP_GT_OQ);

      // Psi integral: (rho_j^3 / 3.0) / r^4
      __m256d rho_j_3 = _mm256_mul_pd(v_rho_j, _mm256_mul_pd(v_rho_j, v_rho_j));
      __m256d vol_j = _mm256_div_pd(rho_j_3, _mm256_set1_pd(3.0));
      __m256d inv_r2 = _mm256_div_pd(avx_one, r2);
      __m256d term = _mm256_mul_pd(vol_j, _mm256_mul_pd(inv_r2, inv_r2));
      term = _mm256_and_pd(term, mask); // zero out self-interaction
      psi_sum = _mm256_add_pd(psi_sum, term);
    }

    double psi = hsum_avx(psi_sum);
    // Scalar remainder
    for (size_t j = size_vec * 4; j < size; ++j) {
      if (i == j)
        continue;
      double dx = mol.pos.x[j] - pos_x;
      double dy = mol.pos.y[j] - pos_y;
      double dz = mol.pos.z[j] - pos_z;
      double r2 = dx * dx + dy * dy + dz * dz;
      if (r2 > 1e-6) {
        double vol_j = vdw_radii[j] * vdw_radii[j] * vdw_radii[j] / 3.0;
        psi += vol_j / (r2 * r2);
      }
    }

    if (env_mol != nullptr && env_vdw_radii != nullptr) {
      size_t env_size = env_mol->size();
      size_t env_size_vec = env_size / 4;
      __m256d env_psi_sum = _mm256_setzero_pd();
      for (size_t j = 0; j < env_size_vec * 4; j += 4) {
        __m256d v_x = _mm256_load_pd(&env_mol->pos.x[j]);
        __m256d v_y = _mm256_load_pd(&env_mol->pos.y[j]);
        __m256d v_z = _mm256_load_pd(&env_mol->pos.z[j]);
        __m256d v_rho_j = _mm256_load_pd(&(*env_vdw_radii)[j]);

        __m256d dx = _mm256_sub_pd(v_x, v_pos_x);
        __m256d dy = _mm256_sub_pd(v_y, v_pos_y);
        __m256d dz = _mm256_sub_pd(v_z, v_pos_z);
        __m256d dx2 = _mm256_mul_pd(dx, dx);
        __m256d dy2 = _mm256_mul_pd(dy, dy);
        __m256d dz2 = _mm256_mul_pd(dz, dz);
        __m256d r2 = _mm256_add_pd(_mm256_add_pd(dx2, dy2), dz2);
        __m256d mask = _mm256_cmp_pd(r2, _mm256_set1_pd(1e-6), _CMP_GT_OQ);

        __m256d rho_j_3 =
            _mm256_mul_pd(v_rho_j, _mm256_mul_pd(v_rho_j, v_rho_j));
        __m256d vol_j = _mm256_div_pd(rho_j_3, _mm256_set1_pd(3.0));
        __m256d inv_r2 = _mm256_div_pd(avx_one, r2);
        __m256d term = _mm256_mul_pd(vol_j, _mm256_mul_pd(inv_r2, inv_r2));
        term = _mm256_and_pd(term, mask);
        env_psi_sum = _mm256_add_pd(env_psi_sum, term);
      }
      psi += hsum_avx(env_psi_sum);
      for (size_t j = env_size_vec * 4; j < env_size; ++j) {
        double dx = env_mol->pos.x[j] - pos_x;
        double dy = env_mol->pos.y[j] - pos_y;
        double dz = env_mol->pos.z[j] - pos_z;
        double r2 = dx * dx + dy * dy + dz * dz;
        if (r2 > 1e-6) {
          double vol_j = (*env_vdw_radii)[j] * (*env_vdw_radii)[j] *
                         (*env_vdw_radii)[j] / 3.0;
          psi += vol_j / (r2 * r2);
        }
      }
    }

    psi *= 0.1; // Convert integral from Angstrom^-4 to nm^-4

    double inv_Ri = 1.0 / rho_i;
    if (igb == 1) {
      inv_Ri -= psi;
    } else {
      double tanh_arg =
          alpha * psi - beta * psi * psi + gamma * psi * psi * psi;
      inv_Ri -= (1.0 / 25.0) *
                std::tanh(tanh_arg); // Rmax = 25.0 A Born radius ceiling
    }

    if (inv_Ri < 0.001)
      inv_Ri = 0.001; // ceiling
    radii[i] = (1.0 / inv_Ri) *
               0.1; // Convert Born radius from A to nm for use in energy kernel
  }
}

// -----------------------------------------------------------------------------
// Shrake-Rupley SASA Calculation
// -----------------------------------------------------------------------------
struct Point3D {
  double x, y, z;
};

static std::vector<Point3D> generate_golden_spiral(int num_points) {
  std::vector<Point3D> points(num_points);
  double inc = 3.14159265358979323846 * (3.0 - std::sqrt(5.0));
  double off = 2.0 / num_points;
  for (int k = 0; k < num_points; ++k) {
    double y = k * off - 1.0 + (off / 2.0);
    double r = std::sqrt(1.0 - y * y);
    double phi = k * inc;
    points[k] = {std::cos(phi) * r, y, std::sin(phi) * r};
  }
  return points;
}

static const int SASA_NUM_POINTS = 96;
static const std::vector<Point3D> SASA_SPHERE_POINTS =
    generate_golden_spiral(SASA_NUM_POINTS);

struct SASA_Grid {
  std::unordered_map<int, std::vector<size_t>> cells;
  double cell_size;
};

inline SASA_Grid build_sasa_grid(const MOL &mol, std::span<const double> radii,
                                 double probe_radius) {
  double max_radius = 0.0;
  for (double r : radii)
    max_radius = std::max(max_radius, r);
  double cell_size = (max_radius + probe_radius) * 2.0;
  if (cell_size < 1.0)
    cell_size = 3.0;

  SASA_Grid grid;
  grid.cell_size = cell_size;

  auto get_cell_key = [cell_size](double x, double y, double z) -> int {
    int cx = static_cast<int>(std::floor(x / cell_size));
    int cy = static_cast<int>(std::floor(y / cell_size));
    int cz = static_cast<int>(std::floor(z / cell_size));
    size_t h = static_cast<size_t>(cx);
    h = h * 1000003 + static_cast<size_t>(cy);
    h = h * 1000003 + static_cast<size_t>(cz);
    return static_cast<int>(h);
  };
  size_t num_atoms = mol.size();
  for (size_t i = 0; i < num_atoms; ++i) {
    grid.cells[get_cell_key(mol.pos.x[i], mol.pos.y[i], mol.pos.z[i])]
        .push_back(i);
  }
  return grid;
}

std::vector<double>
compute_sasa_shrake_rupley(const MOL &mol, std::span<const double> radii,
                           const std::vector<size_t> &targets,
                           const SASA_Grid &grid, double probe_radius = 1.4) {
  size_t num_atoms = mol.size();
  std::vector<double> sasa(num_atoms, 0.0);
  double area_per_point = (4.0 * 3.14159265358979323846) / SASA_NUM_POINTS;

  double cell_size = grid.cell_size;
  auto get_cell_key = [cell_size](double x, double y, double z) -> int {
    int cx = static_cast<int>(std::floor(x / cell_size));
    int cy = static_cast<int>(std::floor(y / cell_size));
    int cz = static_cast<int>(std::floor(z / cell_size));
    size_t h = static_cast<size_t>(cx);
    h = h * 1000003 + static_cast<size_t>(cy);
    h = h * 1000003 + static_cast<size_t>(cz);
    return static_cast<int>(h);
  };

  // Evaluate ONLY targeted atoms for this thread
  for (size_t i : targets) {
    double x_i = mol.pos.x[i], y_i = mol.pos.y[i], z_i = mol.pos.z[i];
    double R_i = radii[i] + probe_radius;
    double R_i_sq = R_i * R_i;

    std::vector<size_t> neighbors;
    neighbors.reserve(64);
    for (int dx = -1; dx <= 1; ++dx) {
      for (int dy = -1; dy <= 1; ++dy) {
        for (int dz = -1; dz <= 1; ++dz) {
          int key = get_cell_key(x_i + dx * cell_size, y_i + dy * cell_size,
                                 z_i + dz * cell_size);
          if (grid.cells.count(key)) {
            for (size_t j : grid.cells.at(key)) {
              if (i != j) {
                double R_j = radii[j] + probe_radius;
                double dist_sq = (x_i - mol.pos.x[j]) * (x_i - mol.pos.x[j]) +
                                 (y_i - mol.pos.y[j]) * (y_i - mol.pos.y[j]) +
                                 (z_i - mol.pos.z[j]) * (z_i - mol.pos.z[j]);
                if (dist_sq < (R_i + R_j) * (R_i + R_j)) {
                  neighbors.push_back(j);
                }
              }
            }
          }
        }
      }
    }

#ifdef __AVX__
    // AVX2 96-Point Intersection Check (24 loops of 4 points)
    __m256d px[24], py[24], pz[24];
    __m256d buried_mask[24];
    for (int k = 0; k < 24; ++k) {
      buried_mask[k] = _mm256_setzero_pd();
      px[k] = _mm256_set_pd(x_i + SASA_SPHERE_POINTS[k * 4 + 3].x * R_i,
                            x_i + SASA_SPHERE_POINTS[k * 4 + 2].x * R_i,
                            x_i + SASA_SPHERE_POINTS[k * 4 + 1].x * R_i,
                            x_i + SASA_SPHERE_POINTS[k * 4 + 0].x * R_i);
      py[k] = _mm256_set_pd(y_i + SASA_SPHERE_POINTS[k * 4 + 3].y * R_i,
                            y_i + SASA_SPHERE_POINTS[k * 4 + 2].y * R_i,
                            y_i + SASA_SPHERE_POINTS[k * 4 + 1].y * R_i,
                            y_i + SASA_SPHERE_POINTS[k * 4 + 0].y * R_i);
      pz[k] = _mm256_set_pd(z_i + SASA_SPHERE_POINTS[k * 4 + 3].z * R_i,
                            z_i + SASA_SPHERE_POINTS[k * 4 + 2].z * R_i,
                            z_i + SASA_SPHERE_POINTS[k * 4 + 1].z * R_i,
                            z_i + SASA_SPHERE_POINTS[k * 4 + 0].z * R_i);
    }

    for (size_t j : neighbors) {
      __m256d jx = _mm256_set1_pd(mol.pos.x[j]);
      __m256d jy = _mm256_set1_pd(mol.pos.y[j]);
      __m256d jz = _mm256_set1_pd(mol.pos.z[j]);
      __m256d Rj2 =
          _mm256_set1_pd((radii[j] + probe_radius) * (radii[j] + probe_radius));

      for (int k = 0; k < 24; ++k) {
        __m256d dx = _mm256_sub_pd(px[k], jx);
        __m256d dy = _mm256_sub_pd(py[k], jy);
        __m256d dz = _mm256_sub_pd(pz[k], jz);
        __m256d dist_sq = _mm256_add_pd(
            _mm256_add_pd(_mm256_mul_pd(dx, dx), _mm256_mul_pd(dy, dy)),
            _mm256_mul_pd(dz, dz));
        __m256d mask = _mm256_cmp_pd(dist_sq, Rj2, _CMP_LT_OQ);
        buried_mask[k] = _mm256_or_pd(buried_mask[k], mask);
      }
    }

    int exposed_points = 0;
    for (int k = 0; k < 24; ++k) {
      int mask_int = _mm256_movemask_pd(buried_mask[k]);
      if ((mask_int & 1) == 0)
        exposed_points++;
      if ((mask_int & 2) == 0)
        exposed_points++;
      if ((mask_int & 4) == 0)
        exposed_points++;
      if ((mask_int & 8) == 0)
        exposed_points++;
    }
#else
    // Scalar fallback if targeted atoms > 0 but not vectorizable fully
    // Count how many points are accessible
    int exposed_points = 0;
    for (const auto &pt : SASA_SPHERE_POINTS) {
      double px_val = x_i + pt.x * R_i, py_val = y_i + pt.y * R_i,
             pz_val = z_i + pt.z * R_i;
      bool buried = false;
      for (size_t j : neighbors) {
        double R_j = radii[j] + probe_radius;
        double dist_sq = (px_val - mol.pos.x[j]) * (px_val - mol.pos.x[j]) +
                         (py_val - mol.pos.y[j]) * (py_val - mol.pos.y[j]) +
                         (pz_val - mol.pos.z[j]) * (pz_val - mol.pos.z[j]);
        if (dist_sq < R_j * R_j) {
          buried = true;
          break;
        }
      }
      if (!buried)
        exposed_points++;
    }
#endif
    sasa[i] = exposed_points * area_per_point * R_i_sq;
  }
  return sasa;
}

static void build_complex_geometry(const MOL &rec, const MOL &lg,
                                   const Multi_MOL_parameters &rec_parameters,
                                   const Multi_MOL_parameters &lg_parameters,
                                   MOL &complex_mol,
                                   std::vector<double> &complex_radii) {
  const size_t rec_size = rec.size();
  const size_t lg_size = lg.size();
  complex_mol.pos.resize(rec_size + lg_size);
  complex_radii.resize(rec_size + lg_size);

  for (size_t i = 0; i < rec_size; ++i) {
    complex_mol.pos.x[i] = rec.pos.x[i];
    complex_mol.pos.y[i] = rec.pos.y[i];
    complex_mol.pos.z[i] = rec.pos.z[i];
    complex_radii[i] = rec_parameters.radius[i];
  }

  for (size_t j = 0; j < lg_size; ++j) {
    complex_mol.pos.x[rec_size + j] = lg.pos.x[j];
    complex_mol.pos.y[rec_size + j] = lg.pos.y[j];
    complex_mol.pos.z[rec_size + j] = lg.pos.z[j];
    complex_radii[rec_size + j] = lg_parameters.radius[j];
  }
}

void build_gbsa_sasa_cache(const MOL &rec, const MOL &lg,
                           const Multi_MOL_parameters &rec_parameters,
                           const Multi_MOL_parameters &lg_parameters,
                           GBSASasaCache &cache_out) {
  const size_t rec_size = rec.size();
  const size_t lg_size = lg.size();

  std::vector<size_t> rec_targets(rec_size);
  std::iota(rec_targets.begin(), rec_targets.end(), size_t{0});

  std::vector<size_t> lig_targets(lg_size);
  std::iota(lig_targets.begin(), lig_targets.end(), size_t{0});

  std::vector<size_t> com_targets(rec_size + lg_size);
  std::iota(com_targets.begin(), com_targets.end(), size_t{0});

  MOL complex_mol;
  std::vector<double> complex_radii;
  build_complex_geometry(rec, lg, rec_parameters, lg_parameters, complex_mol,
                         complex_radii);

  const SASA_Grid rec_grid = build_sasa_grid(
      rec, std::span<const double>(rec_parameters.radius.data(), rec_size),
      1.4);
  cache_out.sasa_rec_free = compute_sasa_shrake_rupley(
      rec, std::span<const double>(rec_parameters.radius.data(), rec_size),
      rec_targets, rec_grid);

  const SASA_Grid lig_grid = build_sasa_grid(
      lg, std::span<const double>(lg_parameters.radius.data(), lg_size), 1.4);
  cache_out.sasa_lig_free = compute_sasa_shrake_rupley(
      lg, std::span<const double>(lg_parameters.radius.data(), lg_size),
      lig_targets, lig_grid);

  const SASA_Grid com_grid =
      build_sasa_grid(complex_mol,
                      std::span<const double>(complex_radii.data(),
                                              rec_size + lg_size),
                      1.4);
  cache_out.sasa_complex = compute_sasa_shrake_rupley(
      complex_mol,
      std::span<const double>(complex_radii.data(), rec_size + lg_size),
      com_targets, com_grid);
}

int cal_energy_gbsa_avx2(
    const MOL &rec, const MOL &lg, const std::vector<mol_info> &info_rec,
    const std::vector<mol_info> &info_lg, cache_sigma_epsilon &CachedData,
    size_t rec_begin, size_t rec_end, size_t lg_begin, size_t lg_end,
    GBSA_energy_snap &Energy, size_t num_snap, const arguments &Arguments,
    const Multi_MOL_parameters &rec_parameters,
    const Multi_MOL_parameters &lg_parameters,
    const std::vector<double> &rec_lcpo_p1,
    const std::vector<double> &rec_lcpo_p2,
    const std::vector<double> &lg_lcpo_p1,
    const std::vector<double> &lg_lcpo_p2,
    const AlignedDoubleVector &rec_free_radii,
    const AlignedDoubleVector &lg_free_radii,
    const AlignedDoubleVector &rec_complex_radii,
    const AlignedDoubleVector &lg_complex_radii, GBSA_Aggregator *aggregator,
    info_ener_res *infoener_res, const GBSASasaCache *sasa_cache,
    bool write_pair_cache) {
  (void)num_snap;
  (void)rec_lcpo_p1;
  (void)rec_lcpo_p2;
  (void)lg_lcpo_p1;
  (void)lg_lcpo_p2;

  bool is_lig_rec =
      (Arguments.mol_1_type == "ligand" || Arguments.mol_1_type == "lig") &&
      !Arguments.flag_use_common_topo;
  bool is_lig_lg =
      (Arguments.mol_2_type == "ligand" || Arguments.mol_2_type == "lig") &&
      !Arguments.flag_use_common_topo;
  __m256d avx_one = _mm256_set1_pd(1.0);
  __m256d avx_four = _mm256_set1_pd(4.0);
  __m256d kappa = _mm256_set1_pd(
      0.73 * std::sqrt(Arguments.saltcon)); // Debye length approx

  double ext_diel = Arguments.ext_diel;
  double int_diel = Arguments.int_diel;
  __m256d avx_half = _mm256_set1_pd(0.5);
  bool is_oplsaa = (Arguments.force_f == "OPLSAA");

  bool use_alpb = Arguments.flag_alpb;
  double alpb_alpha = 2.5; // Empirical effective electrostatic size
  __m256d avx_int_diel = _mm256_set1_pd(int_diel);
  __m256d avx_ext_diel = _mm256_set1_pd(ext_diel);
  __m256d avx_alpb_alpha_term = _mm256_set1_pd(int_diel * alpb_alpha);

  double self_gpol_sum = 0.0;
  double kappa_val = 0.73 * std::sqrt(Arguments.saltcon);
  const bool is_charmm =
      Arguments.force_f.find("CHARMM") != std::string::npos;
  const bool can_write_pair_cache =
      write_pair_cache && !CachedData.elec_ener.empty();

  auto add_pair_map = [&](size_t rec_idx, size_t lig_idx, double pair_elec,
                          double pair_vdw) {
    if (infoener_res == nullptr)
      return;

    size_t rec_cn = info_rec[rec_idx].get_chainnum();
    size_t lg_cn = info_lg[lig_idx].get_chainnum();
    size_t rec_in = is_lig_rec ? (rec_idx + 1)
                               : info_rec[rec_idx].get_internal_resnr();
    size_t lg_in = is_lig_lg ? (lig_idx + 1)
                             : info_lg[lig_idx].get_internal_resnr();

    if (is_lig_rec && !is_lig_lg) {
      infoener_res->add_energy(std::make_tuple(lg_cn, rec_cn, lg_in, rec_in),
                               pair_elec, pair_vdw);
    } else {
      infoener_res->add_energy(std::make_tuple(rec_cn, lg_cn, rec_in, lg_in),
                               pair_elec, pair_vdw);
    }
  };

  const size_t rec_size = rec.size();
  const size_t lg_size = lg.size();
  const bool collect_decomposition = (aggregator != nullptr);
  struct LocalDecompEntry {
    bool is_ligand;
    int chainnum;
    int resnr;
    bool detail_atom;
    const mol_info *m_info;
    double ev = 0.0;
    double ee = 0.0;
    double gpol = 0.0;
    double gnp = 0.0;
  };

  auto make_decomp_key = [](int chainnum, int resnr) -> std::uint64_t {
    return (static_cast<std::uint64_t>(static_cast<std::uint32_t>(chainnum))
            << 32) |
           static_cast<std::uint32_t>(resnr);
  };

  std::vector<LocalDecompEntry> local_decomp_entries;
  if (collect_decomposition) {
    local_decomp_entries.reserve(rec_size + (lg_end - lg_begin));
  }

  auto register_decomp_slot =
      [&](std::unordered_map<std::uint64_t, size_t> &slot_map, bool is_ligand,
          int chainnum, int resnr, bool detail_atom,
          const mol_info *m_info) -> size_t {
    const std::uint64_t key = make_decomp_key(chainnum, resnr);
    const auto it = slot_map.find(key);
    if (it != slot_map.end()) {
      return it->second;
    }

    const size_t slot = local_decomp_entries.size();
    slot_map.emplace(key, slot);
    local_decomp_entries.push_back(
        {is_ligand, chainnum, resnr, detail_atom, m_info});
    return slot;
  };

  std::unordered_map<std::uint64_t, size_t> rec_decomp_slots_by_key;
  std::vector<size_t> rec_decomp_slot(rec_size, size_t{0});
  if (collect_decomposition) {
    rec_decomp_slots_by_key.reserve(rec_size);
    for (size_t i = 0; i < rec_size; ++i) {
      const int rec_id = is_lig_rec ? static_cast<int>(i + 1)
                                    : info_rec[i].get_resnr();
      rec_decomp_slot[i] = register_decomp_slot(
          rec_decomp_slots_by_key, false, info_rec[i].get_chainnum(), rec_id,
          is_lig_rec, &info_rec[i]);
    }
  }

  std::unordered_map<std::uint64_t, size_t> lig_decomp_slots_by_key;
  std::vector<size_t> lig_decomp_slot(lg_end - lg_begin, size_t{0});
  if (collect_decomposition) {
    lig_decomp_slots_by_key.reserve(lg_end - lg_begin);
    for (size_t j = lg_begin; j < lg_end; ++j) {
      const int lig_id = is_lig_lg ? static_cast<int>(j + 1)
                                   : info_lg[j].get_resnr();
      lig_decomp_slot[j - lg_begin] = register_decomp_slot(
          lig_decomp_slots_by_key, true, info_lg[j].get_chainnum(), lig_id,
          is_lig_lg, &info_lg[j]);
    }
  }

  std::vector<double> local_sasa_rec_free;
  std::vector<double> local_sasa_lig_free;
  std::vector<double> local_sasa_complex;
  const std::vector<double> *sasa_rec_free_ptr = nullptr;
  const std::vector<double> *sasa_lig_free_ptr = nullptr;
  const std::vector<double> *sasa_complex_ptr = nullptr;

#if 0
  for (size_t i = 0; i < rec_size; ++i) {
    complex_mol.pos.x[i] = rec.pos.x[i];
    complex_mol.pos.y[i] = rec.pos.y[i];
    complex_mol.pos.z[i] = rec.pos.z[i];
    complex_radii[i] = rec_parameters.radius[i];
  }
  for (size_t j = 0; j < lg_size; ++j) {
    complex_mol.pos.x[rec_size + j] = lg.pos.x[j];
    complex_mol.pos.y[rec_size + j] = lg.pos.y[j];
    complex_mol.pos.z[rec_size + j] = lg.pos.z[j];
    // THÊM DÒNG NÀY ĐỂ TRẢ LẠI BÁN KÍNH VDW CHO LIGAND
    complex_radii[rec_size + j] = lg_parameters.radius[j];
  }

  // --- Thread-Specific Shrake-Rupley SASA Calculation ---
  std::vector<size_t> rec_targets;
  rec_targets.reserve(rec_end - rec_begin);
  for (size_t i = rec_begin; i < rec_end; ++i)
    rec_targets.push_back(i);

  std::vector<size_t> lig_targets;
  lig_targets.reserve(lg_end - lg_begin);
  for (size_t j = lg_begin; j < lg_end; ++j)
    lig_targets.push_back(j);

  std::vector<size_t> com_targets;
  com_targets.reserve((rec_end - rec_begin) + (lg_end - lg_begin));
  for (size_t i = rec_begin; i < rec_end; ++i)
    com_targets.push_back(i);
  for (size_t j = lg_begin; j < lg_end; ++j)
    com_targets.push_back(rec_size + j);

  SASA_Grid rec_grid = build_sasa_grid(
      rec, std::span<const double>(rec_parameters.radius.data(), rec_size),
      1.4);
  std::vector<double> sasa_rec_free = compute_sasa_shrake_rupley(
      rec, std::span<const double>(rec_parameters.radius.data(), rec_size),
      rec_targets, rec_grid);

  SASA_Grid lig_grid = build_sasa_grid(
      lg, std::span<const double>(lg_parameters.radius.data(), lg_size), 1.4);
  std::vector<double> sasa_lig_free = compute_sasa_shrake_rupley(
      lg, std::span<const double>(lg_parameters.radius.data(), lg_size),
      lig_targets, lig_grid);

  SASA_Grid com_grid = build_sasa_grid(
      complex_mol,
      std::span<const double>(complex_radii.data(), rec_size + lg_size), 1.4);
  std::vector<double> sasa_complex = compute_sasa_shrake_rupley(
      complex_mol,
      std::span<const double>(complex_radii.data(), rec_size + lg_size),
      com_targets, com_grid);
#endif

  const bool use_precomputed_sasa =
      sasa_cache != nullptr &&
      sasa_cache->sasa_rec_free.size() == rec_size &&
      sasa_cache->sasa_lig_free.size() == lg_size &&
      sasa_cache->sasa_complex.size() == rec_size + lg_size;

  if (use_precomputed_sasa) {
    sasa_rec_free_ptr = &sasa_cache->sasa_rec_free;
    sasa_lig_free_ptr = &sasa_cache->sasa_lig_free;
    sasa_complex_ptr = &sasa_cache->sasa_complex;
  } else {
    std::vector<size_t> rec_targets;
    rec_targets.reserve(rec_end - rec_begin);
    for (size_t i = rec_begin; i < rec_end; ++i)
      rec_targets.push_back(i);

    std::vector<size_t> lig_targets;
    lig_targets.reserve(lg_end - lg_begin);
    for (size_t j = lg_begin; j < lg_end; ++j)
      lig_targets.push_back(j);

    std::vector<size_t> com_targets;
    com_targets.reserve((rec_end - rec_begin) + (lg_end - lg_begin));
    for (size_t i = rec_begin; i < rec_end; ++i)
      com_targets.push_back(i);
    for (size_t j = lg_begin; j < lg_end; ++j)
      com_targets.push_back(rec_size + j);

    MOL complex_mol;
    std::vector<double> complex_radii;
    build_complex_geometry(rec, lg, rec_parameters, lg_parameters, complex_mol,
                           complex_radii);

    const SASA_Grid rec_grid = build_sasa_grid(
        rec, std::span<const double>(rec_parameters.radius.data(), rec_size),
        1.4);
    local_sasa_rec_free = compute_sasa_shrake_rupley(
        rec, std::span<const double>(rec_parameters.radius.data(), rec_size),
        rec_targets, rec_grid);

    const SASA_Grid lig_grid = build_sasa_grid(
        lg, std::span<const double>(lg_parameters.radius.data(), lg_size), 1.4);
    local_sasa_lig_free = compute_sasa_shrake_rupley(
        lg, std::span<const double>(lg_parameters.radius.data(), lg_size),
        lig_targets, lig_grid);

    const SASA_Grid com_grid =
        build_sasa_grid(complex_mol,
                        std::span<const double>(complex_radii.data(),
                                                rec_size + lg_size),
                        1.4);
    local_sasa_complex = compute_sasa_shrake_rupley(
        complex_mol,
        std::span<const double>(complex_radii.data(), rec_size + lg_size),
        com_targets, com_grid);

    sasa_rec_free_ptr = &local_sasa_rec_free;
    sasa_lig_free_ptr = &local_sasa_lig_free;
    sasa_complex_ptr = &local_sasa_complex;
  }

  double gamma_kj = 0.0072 * 4.184;
  std::vector<double> dG_np_rec(rec_size, 0.0);
  std::vector<double> dG_np_lig(lg_size, 0.0);

  for (size_t i = rec_begin; i < rec_end; ++i) {
    dG_np_rec[i] =
        gamma_kj * ((*sasa_complex_ptr)[i] - (*sasa_rec_free_ptr)[i]);
    Energy.G_nonpolar += dG_np_rec[i];
  }
  for (size_t j = lg_begin; j < lg_end; ++j) {
    dG_np_lig[j] = gamma_kj *
                   ((*sasa_complex_ptr)[rec_size + j] - (*sasa_lig_free_ptr)[j]);
    Energy.G_nonpolar += dG_np_lig[j];
  }

  for (size_t i = rec_begin; i < rec_end; ++i) {
    double q = rec_parameters.charge[i];
    double dG_self;
    if (use_alpb) {
      double term_com_in = 1.0 / (int_diel * rec_complex_radii[i]);
      double term_com_out =
          std::exp(-kappa_val * rec_complex_radii[i]) /
          (ext_diel * rec_complex_radii[i] + int_diel * alpb_alpha);
      double term_free_in = 1.0 / (int_diel * rec_free_radii[i]);
      double term_free_out =
          std::exp(-kappa_val * rec_free_radii[i]) /
          (ext_diel * rec_free_radii[i] + int_diel * alpb_alpha);
      dG_self = -138.935 * 0.5 * q *
                q * // 138.935 = 1/(4pe0) in units of kJ A/(mol e2)
                ((term_com_in - term_com_out) - (term_free_in - term_free_out));
    } else {
      double term_com_in = 1.0 / (int_diel * rec_complex_radii[i]);
      double term_com_out = std::exp(-kappa_val * rec_complex_radii[i]) /
                            (ext_diel * rec_complex_radii[i]);
      double term_free_in = 1.0 / (int_diel * rec_free_radii[i]);
      double term_free_out = std::exp(-kappa_val * rec_free_radii[i]) /
                             (ext_diel * rec_free_radii[i]);
      dG_self = -138.935 * 0.5 * q *
                q * // 138.935 = 1/(4pe0) in units of kJ A/(mol e2)
                ((term_com_in - term_com_out) - (term_free_in - term_free_out));
    }
    self_gpol_sum += dG_self;
  }
  for (size_t j = lg_begin; j < lg_end; ++j) {
    double q = lg_parameters.charge[j];
    double dG_self;
    if (use_alpb) {
      double term_com_in = 1.0 / (int_diel * lg_complex_radii[j]);
      double term_com_out =
          std::exp(-kappa_val * lg_complex_radii[j]) /
          (ext_diel * lg_complex_radii[j] + int_diel * alpb_alpha);
      double term_free_in = 1.0 / (int_diel * lg_free_radii[j]);
      double term_free_out =
          std::exp(-kappa_val * lg_free_radii[j]) /
          (ext_diel * lg_free_radii[j] + int_diel * alpb_alpha);
      dG_self = -138.935 * 0.5 * q * q *
                ((term_com_in - term_com_out) - (term_free_in - term_free_out));
    } else {
      double term_com_in = 1.0 / (int_diel * lg_complex_radii[j]);
      double term_com_out = std::exp(-kappa_val * lg_complex_radii[j]) /
                            (ext_diel * lg_complex_radii[j]);
      double term_free_in = 1.0 / (int_diel * lg_free_radii[j]);
      double term_free_out = std::exp(-kappa_val * lg_free_radii[j]) /
                             (ext_diel * lg_free_radii[j]);
      dG_self = -138.935 * 0.5 * q * q *
                ((term_com_in - term_com_out) - (term_free_in - term_free_out));
    }
    self_gpol_sum += dG_self;
  }
  // FIX: DO NOT divide by scalar_eu_par here. Keep it in kJ/mol.
  Energy.G_polar += self_gpol_sum;

  __m256d v_scale = _mm256_set1_pd(0.01);


  for (size_t j = lg_begin; j < lg_end; ++j) {
    double pos_x = lg.pos.x[j], pos_y = lg.pos.y[j], pos_z = lg.pos.z[j];
    double R_j = lg_complex_radii[j];
    __m256d v_pos_x = _mm256_set1_pd(pos_x);
    __m256d v_pos_y = _mm256_set1_pd(pos_y);
    __m256d v_pos_z = _mm256_set1_pd(pos_z);
    __m256d v_R_j = _mm256_set1_pd(R_j);
    __m256d v_q_j = _mm256_set1_pd(lg_parameters.charge[j]);

    size_t size_vector_rec = rec.size() / 4;
    for (size_t i = 0; i < size_vector_rec * 4; i += 4) {
      __m256d v_x = _mm256_load_pd(&rec.pos.x[i]);
      __m256d v_y = _mm256_load_pd(&rec.pos.y[i]);
      __m256d v_z = _mm256_load_pd(&rec.pos.z[i]);
      __m256d v_q_i = _mm256_load_pd(&rec_parameters.charge[i]);
      __m256d v_R_i = _mm256_load_pd(&rec_complex_radii[i]);

      __m256d dx = _mm256_sub_pd(v_x, v_pos_x);
      __m256d dy = _mm256_sub_pd(v_y, v_pos_y);
      __m256d dz = _mm256_sub_pd(v_z, v_pos_z);
      __m256d r2_A = _mm256_add_pd(
          _mm256_add_pd(_mm256_mul_pd(dx, dx), _mm256_mul_pd(dy, dy)),
          _mm256_mul_pd(dz, dz));
      __m256d r2_nm = _mm256_mul_pd(r2_A, v_scale);
      __m256d r_nm = _mm256_sqrt_pd(r2_nm);

      // VDW Mixing Rules
      __m256d avx_rec_epsilon = _mm256_load_pd(&rec_parameters.epsilon[i]);
      __m256d avx_rec_sigma = _mm256_load_pd(&rec_parameters.sigma[i]);
      __m256d avx_temp_epsilon_lg = _mm256_set1_pd(lg_parameters.epsilon[j]);
      __m256d avx_temp_sigma_lg = _mm256_set1_pd(lg_parameters.sigma[j]);

      __m256d eps_ij = _mm256_sqrt_pd(_mm256_mul_pd(
          avx_temp_epsilon_lg, avx_rec_epsilon)); // Geometric mean
      // FIX: Guard against r ~ 0 causing LJ explosion.
      // Cap r6_nm to a minimum to prevent (sigma/r)^12 -> infinity.
      __m256d r6_nm_raw = _mm256_mul_pd(_mm256_mul_pd(r2_nm, r2_nm), r2_nm);
      __m256d r6_min = _mm256_set1_pd(1e-12); // floor for r^6 in nm^6
      __m256d r6_nm = _mm256_max_pd(r6_nm_raw, r6_min);
      // Mask: zero VdW for pairs with r_nm < 0.05 nm (0.5 A) — bonded/overlap
      __m256d r_min_mask =
          _mm256_cmp_pd(r_nm, _mm256_set1_pd(0.05), _CMP_GE_OQ);

      __m256d vdw_res;
      if (is_charmm) {
        alignas(32) double sigma_lane[4]{};
        alignas(32) double epsilon_lane[4]{};
        charmm::pair_params_block(CachedData, rec_parameters, lg_parameters, i,
                                  4, j, sigma_lane, epsilon_lane);
        const __m256d sigma_vec = _mm256_loadu_pd(sigma_lane);
        const __m256d epsilon_vec = _mm256_loadu_pd(epsilon_lane);
        vdw_res = charmm_lj_energy_gbsa_avx(sigma_vec, epsilon_vec, r_nm);
      } else if (is_oplsaa) {
        __m256d four_eps = _mm256_mul_pd(_mm256_set1_pd(4.0), eps_ij);
        __m256d sig_ij =
            _mm256_sqrt_pd(_mm256_mul_pd(avx_temp_sigma_lg, avx_rec_sigma));
        __m256d sig2 = _mm256_mul_pd(sig_ij, sig_ij);
        __m256d sig6 = _mm256_mul_pd(_mm256_mul_pd(sig2, sig2), sig2);
        __m256d sr6 = _mm256_mul_pd(sig6, _mm256_div_pd(avx_one, r6_nm));
        __m256d sr12 = _mm256_mul_pd(sr6, sr6);
        vdw_res = _mm256_mul_pd(four_eps, _mm256_sub_pd(sr12, sr6));
      } else {
        __m256d four_eps = _mm256_mul_pd(_mm256_set1_pd(4.0), eps_ij);
        __m256d sig_ij = _mm256_mul_pd(
            _mm256_add_pd(avx_temp_sigma_lg, avx_rec_sigma), avx_half);
        __m256d sig2 = _mm256_mul_pd(sig_ij, sig_ij);
        __m256d sig6 = _mm256_mul_pd(_mm256_mul_pd(sig2, sig2), sig2);
        __m256d sr6 = _mm256_mul_pd(sig6, _mm256_div_pd(avx_one, r6_nm));
        __m256d sr12 = _mm256_mul_pd(sr6, sr6);
        vdw_res = _mm256_mul_pd(four_eps, _mm256_sub_pd(sr12, sr6));
      }
      // Apply distance mask: zero out VdW for too-close pairs
      vdw_res = _mm256_and_pd(vdw_res, r_min_mask);

      // Elec (Vacuum Coulomb MUST NOT be divided by int_diel)
      __m256d elec_res = _mm256_div_pd(
          _mm256_mul_pd(_mm256_set1_pd(138.935), _mm256_mul_pd(v_q_i, v_q_j)),
          r_nm);

      // G_polar
      __m256d RiRj = _mm256_mul_pd(v_R_i, v_R_j);
      __m256d exp_term = f_exp_avx(
          _mm256_mul_pd(_mm256_set1_pd(-1.0),
                        _mm256_div_pd(r2_nm, _mm256_mul_pd(avx_four, RiRj))));
      __m256d fGB =
          _mm256_sqrt_pd(_mm256_add_pd(r2_nm, _mm256_mul_pd(RiRj, exp_term)));

      // Proper DH Screening: exp(-kappa * fGB) applied ONLY to external
      // dielectric
      __m256d screen_exp = f_exp_avx(
          _mm256_mul_pd(_mm256_set1_pd(-1.0), _mm256_mul_pd(kappa, fGB)));
      __m256d v_qij = _mm256_mul_pd(v_q_i, v_q_j);
      __m256d gpol_res;

      if (use_alpb) {
        __m256d term_com_in =
            _mm256_div_pd(avx_one, _mm256_mul_pd(avx_int_diel, fGB));
        __m256d term_com_out = _mm256_div_pd(
            screen_exp, _mm256_add_pd(_mm256_mul_pd(avx_ext_diel, fGB),
                                      avx_alpb_alpha_term));
        gpol_res = _mm256_mul_pd(
            _mm256_set1_pd(-138.935),
            _mm256_mul_pd(v_qij, _mm256_sub_pd(term_com_in, term_com_out)));
      } else {
        __m256d term_com_in =
            _mm256_div_pd(avx_one, _mm256_mul_pd(avx_int_diel, fGB));
        __m256d term_com_out =
            _mm256_div_pd(screen_exp, _mm256_mul_pd(avx_ext_diel, fGB));
        gpol_res = _mm256_mul_pd(
            _mm256_set1_pd(-138.935),
            _mm256_mul_pd(v_qij, _mm256_sub_pd(term_com_in, term_com_out)));
      }

      __m256d E_vdw_scaled = vdw_res;
      __m256d E_elec_scaled = elec_res;
      __m256d G_polar_scaled = gpol_res;
      E_elec_scaled = _mm256_and_pd(E_elec_scaled, r_min_mask);
      G_polar_scaled = _mm256_and_pd(G_polar_scaled, r_min_mask);

      double vdw_a[4], elec_a[4], gpol_a[4];
      _mm256_storeu_pd(vdw_a, E_vdw_scaled);
      _mm256_storeu_pd(elec_a, E_elec_scaled);
      _mm256_storeu_pd(gpol_a, G_polar_scaled);

      Energy.E_vdw += vdw_a[0] + vdw_a[1] + vdw_a[2] + vdw_a[3];
      Energy.E_elec += elec_a[0] + elec_a[1] + elec_a[2] + elec_a[3];
      Energy.G_polar += gpol_a[0] + gpol_a[1] + gpol_a[2] + gpol_a[3];

      // FIX: Clean decomposition assignment
      for (int k = 0; k < 4; ++k) {
        if (i + k >= rec_size)
          continue;
        size_t pos = rec_size * j + (i + k);
        if (can_write_pair_cache) {
          CachedData.set_elec_vdw(pos, elec_a[k], vdw_a[k]);
        }
        if (collect_decomposition) {
          auto &rec_decomp = local_decomp_entries[rec_decomp_slot[i + k]];
          rec_decomp.ev += vdw_a[k];
          rec_decomp.ee += elec_a[k];
          rec_decomp.gpol += gpol_a[k];

          auto &lig_decomp =
              local_decomp_entries[lig_decomp_slot[j - lg_begin]];
          lig_decomp.ev += vdw_a[k];
          lig_decomp.ee += elec_a[k];
          lig_decomp.gpol += gpol_a[k];
        }
        add_pair_map(i + k, j, elec_a[k] + gpol_a[k], vdw_a[k]);
      }
    }

    // Scalar Remainder Loop
    for (size_t i = size_vector_rec * 4; i < rec.size(); ++i) {
      double dx = rec.pos.x[i] - pos_x;
      double dy = rec.pos.y[i] - pos_y;
      double dz = rec.pos.z[i] - pos_z;
      double r2_A = dx * dx + dy * dy + dz * dz;
      double r2_nm = r2_A * 0.01;
      double r_nm = std::sqrt(r2_nm);
      double r6_nm = r2_nm * r2_nm * r2_nm;

      // FIX: Guard against r ~ 0 causing LJ explosion
      if (r_nm < 0.05) { // 0.5 Angstrom — overlap/bonded
        if (can_write_pair_cache) {
          size_t pos = rec.size() * j + i;
          CachedData.set_elec_vdw(pos, 0., 0.);
        }
        continue;
      }
      if (r6_nm < 1e-12)
        r6_nm = 1e-12; // floor

      double vdw;
      if (Arguments.force_f.find("CHARMM") != std::string::npos) {
        double sigma_ij{}, epsilon_ij{};
        charmm::pair_params(CachedData, rec_parameters, lg_parameters, i, j,
                            sigma_ij, epsilon_ij);
        vdw = charmm::lj_energy(sigma_ij, epsilon_ij, r_nm);
      } else {
        double eps_ij =
            std::sqrt(lg_parameters.epsilon[j] * rec_parameters.epsilon[i]);
        if (is_oplsaa) {
          double sig_ij =
              std::sqrt(lg_parameters.sigma[j] * rec_parameters.sigma[i]);
          double sig2 = sig_ij * sig_ij;
          double sig6 = sig2 * sig2 * sig2;
          double sr6 = sig6 / r6_nm;
          vdw = 4.0 * eps_ij * (sr6 * sr6 - sr6);
        } else {
          double sig_ij =
              0.5 * (lg_parameters.sigma[j] + rec_parameters.sigma[i]);
          double sig2 = sig_ij * sig_ij;
          double sig6 = sig2 * sig2 * sig2;
          double sr6 = sig6 / r6_nm;
          vdw = 4.0 * eps_ij * (sr6 * sr6 - sr6);
        }
      }

      // Elec (Remove int_diel division)
      double elec =
          138.935 * lg_parameters.charge[j] * rec_parameters.charge[i] / r_nm;

      // G_polar
      double RiRj = rec_complex_radii[i] * R_j;
      double exp_arg = -0.25 * r2_nm / (RiRj);
      double fGB = std::sqrt(r2_nm + RiRj * std::exp(exp_arg));
      double exp_screen = std::exp(-kappa_val * fGB);
      double gpol;

      if (use_alpb) {
        double term_com_in = 1.0 / (int_diel * fGB);
        double term_com_out =
            exp_screen / (ext_diel * fGB + int_diel * alpb_alpha);
        gpol = -138.935 * lg_parameters.charge[j] * rec_parameters.charge[i] *
               (term_com_in - term_com_out);
      } else {
        double term_com_in = 1.0 / (int_diel * fGB);
        double term_com_out = exp_screen / (ext_diel * fGB);
        gpol = -138.935 * lg_parameters.charge[j] * rec_parameters.charge[i] *
               (term_com_in - term_com_out);
      }

      // CRITICAL: Remove all / scalar_eu_par divisions below this point.

      if (can_write_pair_cache) {
        size_t pos = rec.size() * j + i;
        CachedData.set_elec_vdw(pos, elec, vdw);
      }

      Energy.E_vdw += vdw;
      Energy.E_elec += elec;
      Energy.G_polar += gpol;

      if (collect_decomposition) {
        auto &rec_decomp = local_decomp_entries[rec_decomp_slot[i]];
        rec_decomp.ev += vdw;
        rec_decomp.ee += elec;
        rec_decomp.gpol += gpol;

        auto &lig_decomp = local_decomp_entries[lig_decomp_slot[j - lg_begin]];
        lig_decomp.ev += vdw;
        lig_decomp.ee += elec;
        lig_decomp.gpol += gpol;
      }
      add_pair_map(i, j, elec + gpol, vdw);
    }
  }

  // CRITICAL: Final injection of SASA nonpolar energy into decomposition
  // structure Because Gnonpol does not distribute purely pairwise, we assign it
  // post-loop directly into the relevant aggregates
  for (size_t i = rec_begin; i < rec_end; ++i) {
    if (collect_decomposition && std::abs(dG_np_rec[i]) > 1e-8) {
      auto &rec_decomp = local_decomp_entries[rec_decomp_slot[i]];
      rec_decomp.gnp += dG_np_rec[i];
    }
  }
  for (size_t j = lg_begin; j < lg_end; ++j) {
    if (collect_decomposition && std::abs(dG_np_lig[j]) > 1e-8) {
      auto &lig_decomp = local_decomp_entries[lig_decomp_slot[j - lg_begin]];
      lig_decomp.gnp += dG_np_lig[j];
    }
  }

  if (collect_decomposition) {
    for (const auto &entry : local_decomp_entries) {
      if (entry.ev == 0.0 && entry.ee == 0.0 && entry.gpol == 0.0 &&
          entry.gnp == 0.0) {
        continue;
      }
      aggregator->add_decomposition(entry.is_ligand, entry.chainnum,
                                    entry.resnr, entry.detail_atom,
                                    entry.m_info, entry.ev, entry.ee,
                                    entry.gpol, entry.gnp);
    }
  }

  Energy.total_energy =
      Energy.E_vdw + Energy.E_elec + Energy.G_polar + Energy.G_nonpolar;
  return 0;
}
#endif // __AVX__

void GBSA_Aggregator::initialize_trajectory(size_t traj_index) {
  if (trajectories.size() <= traj_index) {
    trajectories.resize(traj_index + 1);
  }
}

void GBSA_Aggregator::add_snapshot(size_t traj_index, GBSA_energy_snap snap) {
  initialize_trajectory(traj_index);
  trajectories[traj_index].emplace_back(std::move(snap));
}

void GBSA_Aggregator::add_decomposition(bool is_ligand, int chainnum, int resnr,
                                        bool detail_atom,
                                        const mol_info *m_info, double ev,
                                        double ee, double gpol, double gnp) {
  DecompKey key{is_ligand, chainnum, resnr};
  auto it = decomp_index.find(key);
  if (it != decomp_index.end()) {
    auto &d = decomp[it->second];
    d.E_vdw += ev;
    d.E_elec += ee;
    d.G_polar += gpol;
    d.G_nonpolar += gnp;
    d.total += (ev + ee + gpol + gnp);
    d.count++;
  } else {
    GBSA_decomposition d;
    d.is_ligand = is_ligand;
    d.chainnum = chainnum;
    d.resnr = resnr;
    d.unit_type_str = is_ligand ? "lig" : "protein";
    d.chain_str = m_info->get_chain();
    d.name_str = detail_atom ? m_info->get_atom() : m_info->get_res_name();
    d.E_vdw = ev;
    d.E_elec = ee;
    d.G_polar = gpol;
    d.G_nonpolar = gnp;
    d.total = ev + ee + gpol + gnp;
    d.count = 1;
    decomp_index[key] = decomp.size();
    decomp.emplace_back(std::move(d));
  }
}

void GBSA_Aggregator::merge(const GBSA_Aggregator &other) {
  for (const auto &od : other.decomp) {
    if (od.count == 0)
      continue;
    DecompKey key{od.is_ligand, od.chainnum, od.resnr};
    auto it = decomp_index.find(key);
    if (it != decomp_index.end()) {
      auto &md = decomp[it->second];
      md.E_vdw += od.E_vdw;
      md.E_elec += od.E_elec;
      md.G_polar += od.G_polar;
      md.G_nonpolar += od.G_nonpolar;
      md.total += od.total;
      md.count += od.count;
    } else {
      decomp_index[key] = decomp.size();
      decomp.emplace_back(od);
    }
  }
}

static double
calc_residue_interaction_entropy(const std::vector<double> &energies,
                                 double eu_par) {
  if (energies.empty())
    return 0.0;
  double avg_e = 0.0;
  for (double e : energies)
    avg_e += e;
  avg_e /= energies.size();

  double sigma = 0.0;
  for (double e : energies)
    sigma += (e - avg_e) * (e - avg_e);
  sigma = std::sqrt(sigma / energies.size());

  // RT MUST match the internal unit of 'energies' (kJ/mol) to make the exponent
  // dimensionless.
  const double kB_kJ = 1.38064e-23 * 6.02214e23 / 1000.0; // kJ/(mol·K)
  const double T = 300.0;                                 // K
  const double RT = kB_kJ * T;

  std::vector<double> valid_e;
  for (double e : energies) {
    if (std::abs(e - avg_e) <= 2.0 * sigma) {
      valid_e.push_back((e - avg_e) / RT); // Now rigorously unitless
    }
  }

  if (valid_e.empty())
    return 0.0;
  double max_e = valid_e[0];
  for (double e : valid_e)
    if (e > max_e)
      max_e = e;

  double log_e = 0.0;
  for (double e : valid_e)
    log_e += std::exp(e - max_e);
  log_e = max_e + std::log(log_e);

  // Convert the final TS energy back to the user's requested output unit via
  // eu_par
  return (-std::log(valid_e.size()) + log_e) * (RT / eu_par);
}

void GBSA_Aggregator::commit_frame_decompositions(
    const GBSA_Aggregator &frame_data, size_t current_frame_idx) {
  // 1. Ensure all newly discovered residues from frame_data exist in global
  // decomp
  for (const auto &fd : frame_data.decomp) {
    DecompKey key{fd.is_ligand, fd.chainnum, fd.resnr};
    if (decomp_index.find(key) == decomp_index.end()) {
      GBSA_decomposition new_d = fd;
      new_d.total = 0;
      new_d.count = 0;
      new_d.E_vdw = 0;
      new_d.E_elec = 0;
      new_d.G_polar = 0;
      new_d.G_nonpolar = 0;
      new_d.frame_energies.resize(current_frame_idx,
                                  0.0); // Pad previous frames with 0.0
      decomp_index[key] = decomp.size();
      decomp.emplace_back(std::move(new_d));
    }
  }

  // 2. Append current frame energy to ALL global residues
  for (auto &gd : decomp) {
    double current_energy = 0.0;
    DecompKey key{gd.is_ligand, gd.chainnum, gd.resnr};
    auto it = frame_data.decomp_index.find(key);
    if (it != frame_data.decomp_index.end()) {
      const auto &fd = frame_data.decomp[it->second];
      current_energy = fd.total;
      gd.E_vdw += fd.E_vdw;
      gd.E_elec += fd.E_elec;
      gd.G_polar += fd.G_polar;
      gd.G_nonpolar += fd.G_nonpolar;
      gd.total += fd.total;
      gd.count++;
    }
    gd.frame_energies.push_back(current_energy);
  }
}

void GBSA_Aggregator::print_all_stats(const std::string &dir_prefix) {
  std::string eu_str = "(" + Arguments.eu + ")";
  // 1. Raw Snapshots
  for (size_t t = 0; t < trajectories.size(); ++t) {
    std::string fname_snap = dir_prefix + "/traj_" + std::to_string(t) + "_" +
                             info_output.get_gbsa_snapshots_filename();
    std::fstream fsnap(fname_snap, std::ios::out);
    if (!fsnap.is_open()) {
      std::cerr << "[GBSA] ERROR: Cannot open output file: " << fname_snap
                << "\n";
      return;
    }
    fsnap << "Time(ps) | E_vdw" << eu_str << " | E_elec" << eu_str
          << " | G_polar" << eu_str << " | G_nonpolar" << eu_str << " | Total"
          << eu_str << "\n";
    for (const auto &s : trajectories[t]) {
      fsnap << std::fixed << std::setprecision(4) << s.time << " | "
            << s.E_vdw / Arguments.eu_par << " | "
            << s.E_elec / Arguments.eu_par << " | "
            << s.G_polar / Arguments.eu_par << " | "
            << s.G_nonpolar / Arguments.eu_par << " | "
            << s.total_energy / Arguments.eu_par << "\n";
    }
    fsnap.close();
  }

  // 2. Time-series Average
  std::string fname_ts =
      dir_prefix + "/" + info_output.get_gbsa_timeseries_avg_filename();
  std::fstream fts(fname_ts, std::ios::out);
  if (!fts.is_open()) {
    std::cerr << "[GBSA] ERROR: Cannot open output file: " << fname_ts << "\n";
    return;
  }
  fts << "Time(ps) | Avg_Total" << eu_str << " | SD_Total" << eu_str
      << " | Avg_VDW" << eu_str << " | Avg_ELEC" << eu_str << " | Avg_GPOL"
      << eu_str << " | Avg_GNONPOL" << eu_str << "\n";

  if (!trajectories.empty()) {
    size_t num_frames = 0;
    for (const auto &traj : trajectories) {
      num_frames = std::max(num_frames, traj.size());
    }
    for (size_t f = 0; f < num_frames; ++f) {
      double avg_tot = 0., avg_vdw = 0., avg_elec = 0., avg_gpol = 0.,
             avg_gnpol = 0.;
      double avg_time_ps = 0.0;
      int count = 0;
      for (size_t t = 0; t < trajectories.size(); ++t) {
        if (f < trajectories[t].size()) {
          avg_time_ps += trajectories[t][f].time;
          avg_tot += trajectories[t][f].total_energy;
          avg_vdw += trajectories[t][f].E_vdw;
          avg_elec += trajectories[t][f].E_elec;
          avg_gpol += trajectories[t][f].G_polar;
          avg_gnpol += trajectories[t][f].G_nonpolar;
          count++;
        }
      }
      if (count > 0) {
        avg_time_ps /= count;
        avg_tot /= count;
        avg_vdw /= count;
        avg_elec /= count;
        avg_gpol /= count;
        avg_gnpol /= count;
        double var_tot = 0.;
        for (size_t t = 0; t < trajectories.size(); ++t) {
          if (f < trajectories[t].size()) {
            var_tot += (trajectories[t][f].total_energy - avg_tot) *
                       (trajectories[t][f].total_energy - avg_tot);
          }
        }
        double sd_tot = std::sqrt(var_tot / count);
        fts << std::fixed << std::setprecision(4) << avg_time_ps << " | "
            << avg_tot / Arguments.eu_par << " | " << sd_tot / Arguments.eu_par
            << " | " << avg_vdw / Arguments.eu_par << " | "
            << avg_elec / Arguments.eu_par << " | "
            << avg_gpol / Arguments.eu_par << " | "
            << avg_gnpol / Arguments.eu_par << "\n";
      }
    }
  }
  fts.close();

  // 3. Summary File & Calculate Total Frames
  std::string fname_sum =
      dir_prefix + "/" + info_output.get_gbsa_summary_filename();
  std::fstream fsum(fname_sum, std::ios::out);
  if (!fsum.is_open()) {
    std::cerr << "[GBSA] ERROR: Cannot open output file: " << fname_sum << "\n";
    return;
  }
  fsum << "GRAND AVERAGE OVERALL TRAJECTORIES AND FRAMES " << eu_str << ":\n";
  double g_tot = 0., g_vdw = 0., g_elec = 0., g_pol = 0., g_nonpol = 0.;

  int total_frames = 0;
  std::vector<double> global_total_energies;

  for (size_t t = 0; t < trajectories.size(); ++t) {
    total_frames += trajectories[t].size();
    for (size_t f = 0; f < trajectories[t].size(); ++f) {
      g_tot += trajectories[t][f].total_energy;
      g_vdw += trajectories[t][f].E_vdw;
      g_elec += trajectories[t][f].E_elec;
      g_pol += trajectories[t][f].G_polar;
      g_nonpol += trajectories[t][f].G_nonpolar;

      global_total_energies.push_back(trajectories[t][f].total_energy);
    }
  }

  if (total_frames > 0) {
    double avg_tot = (g_tot / total_frames) / Arguments.eu_par;
    double avg_vdw = (g_vdw / total_frames) / Arguments.eu_par;
    double avg_elec = (g_elec / total_frames) / Arguments.eu_par;
    double avg_gpol = (g_pol / total_frames) / Arguments.eu_par;
    double avg_gnonpol = (g_nonpol / total_frames) / Arguments.eu_par;

    double global_minus_TdS = calc_residue_interaction_entropy(
        global_total_energies, Arguments.eu_par);
    double global_dG_binding = avg_tot + global_minus_TdS;

    fsum << "Total(dH): " << avg_tot << " | VDW: " << avg_vdw
         << " | Elec: " << avg_elec << " | Gpol: " << avg_gpol
         << " | Gnonpol: " << avg_gnonpol << "\n";
    fsum << "-TdS_IE  : " << global_minus_TdS << "\n";
    fsum << "Final_dG : " << global_dG_binding << "\n\n";
  }

  fsum << "AVERAGE PER TRAJECTORY " << eu_str << ":\n";
  for (size_t t = 0; t < trajectories.size(); ++t) {
    double t_tot = 0.;
    std::vector<double> traj_energies;
    for (const auto &s : trajectories[t]) {
      t_tot += s.total_energy;
      traj_energies.push_back(s.total_energy);
    }
    if (!trajectories[t].empty()) {
      double t_avg = (t_tot / trajectories[t].size()) / Arguments.eu_par;
      double t_minus_TdS =
          calc_residue_interaction_entropy(traj_energies, Arguments.eu_par);
      double t_dG = t_avg + t_minus_TdS;

      fsum << "Trajectory " << t << " Avg Total(dH) = " << t_avg
           << " | -TdS_IE = " << t_minus_TdS << " | Final dG = " << t_dG
           << "\n";
    }
  }
  fsum.close();

  // 4. Decomposition (Split into Receptor and Ligand)
  std::string base_dec =
      dir_prefix + "/" + info_output.get_gbsa_decomposition_filename();

  size_t dot_pos = base_dec.find_last_of(".");
  std::string prefix_dec =
      (dot_pos != std::string::npos) ? base_dec.substr(0, dot_pos) : base_dec;

  std::string fname_dec_rec = prefix_dec + "_receptor.dat";
  std::string fname_dec_lig = prefix_dec + "_ligand.dat";

  std::fstream fdec_rec(fname_dec_rec, std::ios::out);
  if (!fdec_rec.is_open()) {
    std::cerr << "[GBSA] ERROR: Cannot open output file: " << fname_dec_rec
              << "\n";
    return;
  }
  std::fstream fdec_lig(fname_dec_lig, std::ios::out);
  if (!fdec_lig.is_open()) {
    std::cerr << "[GBSA] ERROR: Cannot open output file: " << fname_dec_lig
              << "\n";
    return;
  }

  std::string header = "Chain | ID | ResName/AtomName | E_vdw" + eu_str +
                       " | E_elec" + eu_str + " | G_polar" + eu_str +
                       " | G_nonpolar" + eu_str + " | Total" + eu_str +
                       " | -TdS_IE" + eu_str + " | dG_Binding" + eu_str + "\n";
  fdec_rec << header;
  fdec_lig << header;

  // FIX: Chia năng lượng phân rã cho total_frames thay vì d.count
  double div = (total_frames > 0) ? (double)total_frames : 1.0;

  for (const auto &d : decomp) {
    // Chỉ in ra những unit thực sự có dữ liệu
    if (d.count > 0) {
      // Enforce user-requested unit conversion
      double avg_vdw = (d.E_vdw / div) / Arguments.eu_par;
      double avg_elec = (d.E_elec / div) / Arguments.eu_par;
      double avg_gpol = (d.G_polar / div) / Arguments.eu_par;
      double avg_gnp = (d.G_nonpolar / div) / Arguments.eu_par;
      double avg_total = (d.total / div) / Arguments.eu_par;

      double minus_TdS =
          calc_residue_interaction_entropy(d.frame_energies, Arguments.eu_par);
      double dG_binding = avg_total + minus_TdS;

      std::fstream &target_file = d.is_ligand ? fdec_lig : fdec_rec;
      target_file << std::fixed << d.chain_str << " | " << std::setprecision(4)
                  << d.resnr << " | " << d.name_str << " | " << avg_vdw << " | "
                  << avg_elec << " | " << avg_gpol << " | " << avg_gnp << " | "
                  << avg_total << " | " << minus_TdS << " | " << dG_binding
                  << "\n";
    }
  }

  fdec_rec.close();
  fdec_lig.close();
}
