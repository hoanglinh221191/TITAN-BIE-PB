#include "gbsa_intrin.h"
#include "charmm_nonbonded.h"
#include "titan_simd_dispatch.h"
#include "titan_sasa.h"
#include "titan_ui.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <span>

#ifdef __AVX__
#include "avx_math.h"
#include <immintrin.h>

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

  for (size_t i = begin; i < end; ++i) {
    double pos_x = mol.pos.x[i], pos_y = mol.pos.y[i], pos_z = mol.pos.z[i];
    double rho_i = vdw_radii[i]; // In Angstroms!
    if (rho_i < 1.0)
      rho_i = 1.5; // fallback
    double psi = titan_simd_dispatch::gbsa_inverse_r4_volume_sum(
        mol.pos.x.data(), mol.pos.y.data(), mol.pos.z.data(), vdw_radii.data(),
        size, pos_x, pos_y, pos_z, 1e-6, true);

    if (env_mol != nullptr && env_vdw_radii != nullptr) {
      psi += titan_simd_dispatch::gbsa_inverse_r4_volume_sum(
          env_mol->pos.x.data(), env_mol->pos.y.data(), env_mol->pos.z.data(),
          env_vdw_radii->data(), env_mol->size(), pos_x, pos_y, pos_z, 1e-6,
          true);
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

static std::vector<TitanSASA::Atom>
make_sasa_atoms(const MOL &mol, std::span<const double> radii) {
  const size_t n = mol.size();
  std::vector<TitanSASA::Atom> atoms;
  atoms.reserve(n);
  for (size_t i = 0; i < n; ++i) {
    const double radius = (i < radii.size()) ? radii[i] : 0.0;
    atoms.push_back({mol.pos.x[i], mol.pos.y[i], mol.pos.z[i], radius});
  }
  return atoms;
}

static std::vector<double>
compute_gbsa_sasa(const MOL &mol, std::span<const double> radii,
                  std::span<const size_t> targets = {}) {
  std::vector<TitanSASA::Atom> atoms = make_sasa_atoms(mol, radii);
  TitanSASA::Options options;
  options.probe_radius = TitanSASA::kDefaultProbeRadiusA;
  options.min_cell_size = 3.0;
  options.point_count = TitanSASA::kDefaultPointCount;
  return TitanSASA::compute_per_atom_sasa(
      std::span<const TitanSASA::Atom>(atoms.data(), atoms.size()), options,
      targets);
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

  cache_out.sasa_rec_free = compute_gbsa_sasa(
      rec, std::span<const double>(rec_parameters.radius.data(), rec_size),
      std::span<const size_t>(rec_targets.data(), rec_targets.size()));

  cache_out.sasa_lig_free = compute_gbsa_sasa(
      lg, std::span<const double>(lg_parameters.radius.data(), lg_size),
      std::span<const size_t>(lig_targets.data(), lig_targets.size()));

  cache_out.sasa_complex = compute_gbsa_sasa(
      complex_mol,
      std::span<const double>(complex_radii.data(), rec_size + lg_size),
      std::span<const size_t>(com_targets.data(), com_targets.size()));
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

    local_sasa_rec_free = compute_gbsa_sasa(
        rec, std::span<const double>(rec_parameters.radius.data(), rec_size),
        std::span<const size_t>(rec_targets.data(), rec_targets.size()));

    local_sasa_lig_free = compute_gbsa_sasa(
        lg, std::span<const double>(lg_parameters.radius.data(), lg_size),
        std::span<const size_t>(lig_targets.data(), lig_targets.size()));

    local_sasa_complex = compute_gbsa_sasa(
        complex_mol,
        std::span<const double>(complex_radii.data(), rec_size + lg_size),
        std::span<const size_t>(com_targets.data(), com_targets.size()));

    sasa_rec_free_ptr = &local_sasa_rec_free;
    sasa_lig_free_ptr = &local_sasa_lig_free;
    sasa_complex_ptr = &local_sasa_complex;
  }

  double gamma_kj = TitanSASA::kNonpolarGammaKJPerMolA2;
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
      } else {
        alignas(32) double sigma_lane[4]{};
        alignas(32) double epsilon_lane[4]{};
        alignas(32) double r6_lane[4]{};
        alignas(32) double vdw_lane[4]{};
        const __m256d sig_ij = is_oplsaa
                                   ? _mm256_sqrt_pd(_mm256_mul_pd(
                                         avx_temp_sigma_lg, avx_rec_sigma))
                                   : _mm256_mul_pd(
                                         _mm256_add_pd(avx_temp_sigma_lg,
                                                       avx_rec_sigma),
                                         avx_half);
        _mm256_storeu_pd(sigma_lane, sig_ij);
        _mm256_storeu_pd(epsilon_lane, eps_ij);
        _mm256_storeu_pd(r6_lane, r6_nm);
        titan_simd_dispatch::gbsa_lj_mixed_sigma_epsilon_r6_fourlanes(
            sigma_lane, epsilon_lane, r6_lane, vdw_lane);
        vdw_res = _mm256_loadu_pd(vdw_lane);
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
      new_d.frame_mm_energies.resize(current_frame_idx,
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
      current_energy = fd.E_vdw + fd.E_elec;
      gd.E_vdw += fd.E_vdw;
      gd.E_elec += fd.E_elec;
      gd.G_polar += fd.G_polar;
      gd.G_nonpolar += fd.G_nonpolar;
      gd.total += fd.total;
      gd.count++;
    }
    gd.frame_mm_energies.push_back(current_energy);
  }
}

void GBSA_Aggregator::print_all_stats(const std::string &dir_prefix,
                                      std::fstream *log) {
  const bool capture_performance = Arguments.gbsa_perf != "off";
  const auto entropy_output_begin =
      capture_performance ? std::chrono::steady_clock::now()
                          : std::chrono::steady_clock::time_point{};
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
  fsum << "interaction_entropy_source E_MM_elec_plus_vdw\n";
  fsum << "GRAND AVERAGE OVERALL TRAJECTORIES AND FRAMES " << eu_str << ":\n";
  double g_tot = 0., g_vdw = 0., g_elec = 0., g_pol = 0., g_nonpol = 0.;

  int total_frames = 0;
  std::vector<double> global_mm_energies;

  for (size_t t = 0; t < trajectories.size(); ++t) {
    total_frames += trajectories[t].size();
    for (size_t f = 0; f < trajectories[t].size(); ++f) {
      g_tot += trajectories[t][f].total_energy;
      g_vdw += trajectories[t][f].E_vdw;
      g_elec += trajectories[t][f].E_elec;
      g_pol += trajectories[t][f].G_polar;
      g_nonpol += trajectories[t][f].G_nonpolar;

      global_mm_energies.push_back(trajectories[t][f].E_vdw +
                                   trajectories[t][f].E_elec);
    }
  }

  if (total_frames > 0) {
    double avg_tot = (g_tot / total_frames) / Arguments.eu_par;
    double avg_vdw = (g_vdw / total_frames) / Arguments.eu_par;
    double avg_elec = (g_elec / total_frames) / Arguments.eu_par;
    double avg_gpol = (g_pol / total_frames) / Arguments.eu_par;
    double avg_gnonpol = (g_nonpol / total_frames) / Arguments.eu_par;

    double global_minus_TdS = calc_residue_interaction_entropy(
        global_mm_energies, Arguments.eu_par);
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
    std::vector<double> traj_mm_energies;
    for (const auto &s : trajectories[t]) {
      t_tot += s.total_energy;
      traj_mm_energies.push_back(s.E_vdw + s.E_elec);
    }
    if (!trajectories[t].empty()) {
      double t_avg = (t_tot / trajectories[t].size()) / Arguments.eu_par;
      double t_minus_TdS =
          calc_residue_interaction_entropy(traj_mm_energies, Arguments.eu_par);
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
  fdec_rec << "# interaction_entropy_source E_MM_elec_plus_vdw\n" << header;
  fdec_lig << "# interaction_entropy_source E_MM_elec_plus_vdw\n" << header;

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
          calc_residue_interaction_entropy(d.frame_mm_energies,
                                           Arguments.eu_par);
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

  if (!capture_performance)
    return;
  const double entropy_output_ms =
      std::chrono::duration<double, std::milli>(
          std::chrono::steady_clock::now() - entropy_output_begin)
          .count();

  double precompute_ms = 0.0;
  double energy_kernel_ms = 0.0;
  double decomp_merge_ms = 0.0;
  double total_ms = 0.0;
  size_t frame_count = 0;
  for (const auto &trajectory : trajectories) {
    for (const auto &snapshot : trajectory) {
      precompute_ms += snapshot.time_precompute_ms;
      energy_kernel_ms += snapshot.time_energy_kernel_ms;
      decomp_merge_ms += snapshot.time_decomposition_merge_ms;
      total_ms += snapshot.time_total_ms;
      ++frame_count;
    }
  }

  std::ofstream perf(dir_prefix + "/gbsa_performance.dat",
                     std::ios::out | std::ios::trunc);
  if (perf.is_open()) {
    perf << "# TITAN MM/GBSA performance summary\n";
    perf << "timing_granularity natural_boundaries_fused_energy_kernel\n";
    perf << "interaction_entropy_source E_MM_elec_plus_vdw\n";
    perf << "frames " << frame_count << "\n";
    perf << std::fixed << std::setprecision(6);
    perf << "gbsa_precompute_ms_total " << precompute_ms << "\n";
    perf << "gbsa_energy_kernel_ms_total " << energy_kernel_ms << "\n";
    perf << "decomposition_merge_ms_total " << decomp_merge_ms << "\n";
    perf << "entropy_output_ms " << entropy_output_ms << "\n";
    perf << "frame_total_ms_total " << total_ms << "\n";
    perf << "total_ms_total " << total_ms << "\n";
  } else {
    TitanUI::print_warning("GBSA: failed to write gbsa_performance.dat.");
  }

  if (Arguments.gbsa_perf == "detail") {
    std::ofstream frames(dir_prefix + "/gbsa_performance_frames.csv",
                         std::ios::out | std::ios::trunc);
    if (frames.is_open()) {
      frames << "trajectory_id,frame_index,time_ps,gbsa_precompute_ms,"
                "gbsa_energy_kernel_ms,decomposition_merge_ms,"
                "frame_total_ms,status\n";
      frames << std::setprecision(10);
      for (size_t trajectory_id = 0; trajectory_id < trajectories.size();
           ++trajectory_id) {
        for (size_t frame_index = 0;
             frame_index < trajectories[trajectory_id].size(); ++frame_index) {
          const auto &snapshot = trajectories[trajectory_id][frame_index];
          frames << trajectory_id << ',' << frame_index << ',' << snapshot.time
                 << ',' << snapshot.time_precompute_ms << ','
                 << snapshot.time_energy_kernel_ms << ','
                 << snapshot.time_decomposition_merge_ms << ','
                 << snapshot.time_total_ms << ",ok\n";
        }
      }
    } else {
      TitanUI::print_warning(
          "GBSA: failed to write gbsa_performance_frames.csv.");
    }
  }

  if (log != nullptr && log->is_open()) {
    TitanUI::log_section(*log, "MM/GBSA COMPONENT TIMING");
    TitanUI::log_kv(*log, "Timing Granularity",
                    "natural_boundaries_fused_energy_kernel");
    TitanUI::log_kv(*log, "IE Source", "E_MM_elec_plus_vdw");
    TitanUI::log_time(*log, "GBSA Precompute", precompute_ms / 1000.0);
    TitanUI::log_time(*log, "GBSA Energy Kernel (MM + GB pair fused)",
                      energy_kernel_ms / 1000.0);
    TitanUI::log_time(*log, "Decomposition Merge",
                      decomp_merge_ms / 1000.0);
    TitanUI::log_time(*log, "Entropy + Output", entropy_output_ms / 1000.0);
    TitanUI::log_time(*log, "Frame Total", total_ms / 1000.0);
  }
}
