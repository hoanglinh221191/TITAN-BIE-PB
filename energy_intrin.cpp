#include "cache_sigma_epsilon.h"
#include "force.h"
#include "functions.h"
#include "martini_nonbonded.h"
#include "miscellaneous.h"
#include "new_emap.h"
#include "stat_ener_res_res.h"
#include "uti_func.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <set>
#include <string>
#include <thread>
#include <unordered_set>
#include <utility>
#include <vector>

#ifdef __AVX__
#include "avx_math.h"
#include <immintrin.h>

int cal_one_ele(const MOL &rec, const MOL &lg, cache_sigma_epsilon &CachedData,
                const arguments &Arguments, double pos_x, double pos_y,
                double pos_z, double temp_sigma, double temp_charge,
                double temp_epsilon_lg,
                const Multi_MOL_parameters &rec_parameters,
                const Multi_MOL_parameters &lg_parameters,
                const vector<mol_info> &info_rec,
                const vector<mol_info> &info_lg, cache_energy &Energy, size_t i,
                size_t j);

int cal_energy_amber_intrin(const MOL &rec, const MOL &lg,
                            const vector<mol_info> &info_rec,
                            const vector<mol_info> &info_lg,
                            cache_sigma_epsilon &CachedData, size_t begin,
                            size_t end, cache_energy &Energy, size_t num_snap,
                            const arguments &Arguments,
                            const Multi_MOL_parameters &rec_parameters,
                            const Multi_MOL_parameters &lg_parameters) {
  Energy.set_energy(0., 0.);
  __m256d avx_half = _mm256_set1_pd(0.5);
  __m256d v_scale = _mm256_set1_pd(0.01);
  const size_t rec_size{rec.size()};
  const size_t size_vector_rec{rec.size() / 4};
  double epsilon_in = Arguments.epsilon_in;

  if (CachedData.elec_ener.empty()) {
    for (size_t j{begin}; j < end; ++j) {
      const double pos_x{lg.pos.x[j]}, pos_y{lg.pos.y[j]}, pos_z{lg.pos.z[j]};

      __m256d avx_temp_epsilon_lg =
          _mm256_set1_pd(16. * lg_parameters.epsilon[j]);
      __m256d avx_temp_charge_lg =
          _mm256_set1_pd(138.935457839 * lg_parameters.charge[j] / epsilon_in);
      __m256d avx_temp_sigma = _mm256_set1_pd(lg_parameters.sigma[j]);

      __m256d v_pos_x = _mm256_set1_pd(pos_x);
      __m256d v_pos_y = _mm256_set1_pd(pos_y);
      __m256d v_pos_z = _mm256_set1_pd(pos_z);

      for (size_t i = 0; i < size_vector_rec * 4; i += 4) {
        __m256d v_x = _mm256_loadu_pd(&rec.pos.x[i]);
        __m256d v_y = _mm256_loadu_pd(&rec.pos.y[i]);
        __m256d v_z = _mm256_loadu_pd(&rec.pos.z[i]);

        __m256d dx = _mm256_sub_pd(v_x, v_pos_x);
        __m256d dy = _mm256_sub_pd(v_y, v_pos_y);
        __m256d dz = _mm256_sub_pd(v_z, v_pos_z);

        __m256d dx2 = _mm256_mul_pd(dx, dx);
        __m256d dy2 = _mm256_mul_pd(dy, dy);
        __m256d dz2 = _mm256_mul_pd(dz, dz);

        __m256d r2 = _mm256_add_pd(_mm256_add_pd(dx2, dy2), dz2);
        r2 = _mm256_mul_pd(r2, v_scale);
        __m256d avx_array_d_rec = _mm256_sqrt_pd(r2);

        __m256d avx_array_charge_rec =
            _mm256_loadu_pd(&rec_parameters.charge[i]);
        __m256d avx_elec_result =
            _mm256_mul_pd(avx_temp_charge_lg, avx_array_charge_rec);
        avx_elec_result = _mm256_div_pd(avx_elec_result, avx_array_d_rec);

        __m256d vdw_temp = _mm256_loadu_pd(&rec_parameters.epsilon[i]);
        vdw_temp = _mm256_mul_pd(avx_temp_epsilon_lg, vdw_temp);
        vdw_temp = _mm256_sqrt_pd(vdw_temp);

        __m256d vdw_temp2 = _mm256_loadu_pd(&rec_parameters.sigma[i]);
        vdw_temp2 = _mm256_add_pd(avx_temp_sigma, vdw_temp2);
        vdw_temp2 = _mm256_mul_pd(vdw_temp2, avx_half);

        vdw_temp2 = _mm256_div_pd(vdw_temp2, avx_array_d_rec);
        vdw_temp2 = _mm256_mul_pd(vdw_temp2, vdw_temp2);
        __m256d vdw_temp3 = vdw_temp2;
        vdw_temp2 = _mm256_mul_pd(vdw_temp2, vdw_temp2);
        vdw_temp2 = _mm256_mul_pd(vdw_temp2, vdw_temp3);
        vdw_temp3 =
            _mm256_sub_pd(_mm256_mul_pd(vdw_temp2, vdw_temp2), vdw_temp2);
        vdw_temp3 = _mm256_mul_pd(vdw_temp, vdw_temp3);

        Energy.elec_ener += hsum_double_avx(avx_elec_result);
        Energy.vdw_ener += hsum_double_avx(vdw_temp3);
      }

      if (size_vector_rec * 4 < rec_size) {
        for (size_t i{size_vector_rec * 4}; i < rec_size; ++i) {
          double d{0.1 * sqrt(square(rec.pos.x[i] - pos_x) +
                              square(rec.pos.y[i] - pos_y) +
                              square(rec.pos.z[i] - pos_z))},
              c_6{c6(0.5 * (lg_parameters.sigma[j] + rec_parameters.sigma[i]) /
                     d)};
          if (Arguments.distance_energy_cutoff > 0. &&
              d > Arguments.distance_energy_cutoff) {
            continue;
          }
          if (Arguments.self_int == true &&
              (d < 0.5 * (lg_parameters.sigma[j] + rec_parameters.sigma[i]) ||
               (Arguments.flag_interchain_energy == true &&
                info_lg[j].get_chain() == info_rec[i].get_chain()))) {
            continue;
          }
          Energy.elec_ener +=
              (138.935457839 * lg_parameters.charge[j] / epsilon_in) *
              rec_parameters.charge[i] / d;
          if (lg_parameters.sigma[j] < 0. || rec_parameters.sigma[i] < 0.)
            Energy.vdw_ener += sqrt(16. * lg_parameters.epsilon[j] *
                                    rec_parameters.epsilon[i]) *
                               (square(c_6) - 0.000);
          else
            Energy.vdw_ener += sqrt(16. * lg_parameters.epsilon[j] *
                                    rec_parameters.epsilon[i]) *
                               (square(c_6) - c_6);
        }
      }
    }
  } else {
    for (size_t j{begin}; j < end; ++j) {
      const double pos_x{lg.pos.x[j]}, pos_y{lg.pos.y[j]}, pos_z{lg.pos.z[j]};
      double elec_temp[4]{}, vdw_temp_d[4]{};

      __m256d avx_temp_epsilon_lg =
          _mm256_set1_pd(16. * lg_parameters.epsilon[j]);
      __m256d avx_temp_charge_lg =
          _mm256_set1_pd(138.935457839 * lg_parameters.charge[j] / epsilon_in);
      __m256d avx_temp_sigma = _mm256_set1_pd(lg_parameters.sigma[j]);

      __m256d v_pos_x = _mm256_set1_pd(pos_x);
      __m256d v_pos_y = _mm256_set1_pd(pos_y);
      __m256d v_pos_z = _mm256_set1_pd(pos_z);

      for (size_t i = 0; i < size_vector_rec * 4; i += 4) {
        __m256d v_x = _mm256_loadu_pd(&rec.pos.x[i]);
        __m256d v_y = _mm256_loadu_pd(&rec.pos.y[i]);
        __m256d v_z = _mm256_loadu_pd(&rec.pos.z[i]);

        __m256d dx = _mm256_sub_pd(v_x, v_pos_x);
        __m256d dy = _mm256_sub_pd(v_y, v_pos_y);
        __m256d dz = _mm256_sub_pd(v_z, v_pos_z);

        __m256d dx2 = _mm256_mul_pd(dx, dx);
        __m256d dy2 = _mm256_mul_pd(dy, dy);
        __m256d dz2 = _mm256_mul_pd(dz, dz);

        __m256d r2 = _mm256_add_pd(_mm256_add_pd(dx2, dy2), dz2);
        r2 = _mm256_mul_pd(r2, v_scale);
        __m256d avx_array_d_rec = _mm256_sqrt_pd(r2);

        __m256d avx_array_charge_rec =
            _mm256_loadu_pd(&rec_parameters.charge[i]);
        __m256d avx_elec_result =
            _mm256_mul_pd(avx_temp_charge_lg, avx_array_charge_rec);
        avx_elec_result = _mm256_div_pd(avx_elec_result, avx_array_d_rec);

        __m256d vdw_temp = _mm256_loadu_pd(&rec_parameters.epsilon[i]);
        vdw_temp = _mm256_mul_pd(avx_temp_epsilon_lg, vdw_temp);
        vdw_temp = _mm256_sqrt_pd(vdw_temp);

        __m256d vdw_temp2 = _mm256_loadu_pd(&rec_parameters.sigma[i]);
        vdw_temp2 = _mm256_add_pd(avx_temp_sigma, vdw_temp2);
        vdw_temp2 = _mm256_mul_pd(vdw_temp2, avx_half);

        vdw_temp2 = _mm256_div_pd(vdw_temp2, avx_array_d_rec);
        vdw_temp2 = _mm256_mul_pd(vdw_temp2, vdw_temp2);
        __m256d vdw_temp3 = vdw_temp2;
        vdw_temp2 = _mm256_mul_pd(vdw_temp2, vdw_temp2);
        vdw_temp2 = _mm256_mul_pd(vdw_temp2, vdw_temp3);
        vdw_temp3 =
            _mm256_sub_pd(_mm256_mul_pd(vdw_temp2, vdw_temp2), vdw_temp2);
        vdw_temp3 = _mm256_mul_pd(vdw_temp, vdw_temp3);

        _mm256_storeu_pd(elec_temp, avx_elec_result);
        _mm256_storeu_pd(vdw_temp_d, vdw_temp3);

        if (Arguments.distance_energy_cutoff < 0.) {
          Energy.elec_ener += hsum_double_avx(avx_elec_result);
          Energy.vdw_ener += hsum_double_avx(vdw_temp3);
          size_t pos_base = rec_size * j + i;
          for (int k = 0; k < 4; ++k) {
            CachedData.set_elec_vdw(pos_base + k, elec_temp[k], vdw_temp_d[k]);
          }
        } else {
          double d_arr[4];
          _mm256_storeu_pd(d_arr, avx_array_d_rec);
          size_t pos_base = rec_size * j + i;
          for (int k = 0; k < 4; ++k) {
            size_t pos = pos_base + k;
            if (d_arr[k] > Arguments.distance_energy_cutoff) {
              CachedData.set_elec_vdw(pos, 0., 0.);
              continue;
            }
            CachedData.set_elec_vdw(pos, elec_temp[k], vdw_temp_d[k]);
            Energy.elec_ener += elec_temp[k];
            Energy.vdw_ener += vdw_temp_d[k];
          }
        }
      }

      if (size_vector_rec * 4 < rec_size) {
        for (size_t i{size_vector_rec * 4}; i < rec_size; ++i) {
          double d{0.1 * sqrt(square(rec.pos.x[i] - pos_x) +
                              square(rec.pos.y[i] - pos_y) +
                              square(rec.pos.z[i] - pos_z))},
              c_6{c6(0.5 * (lg_parameters.sigma[j] + rec_parameters.sigma[i]) /
                     d)};
          size_t pos{rec_size * j + i};

          if (Arguments.distance_energy_cutoff > 0. &&
              d > Arguments.distance_energy_cutoff) {
            CachedData.set_elec_vdw(pos, 0., 0.);
            continue;
          }
          if (Arguments.self_int == true &&
              (d < 0.5 * (lg_parameters.sigma[j] + rec_parameters.sigma[i]) ||
               (Arguments.flag_interchain_energy == true &&
                info_lg[j].get_chain() == info_rec[i].get_chain()))) {
            CachedData.set_elec_vdw(pos, 0., 0.);
            continue;
          }

          double elec_val =
              (138.935457839 * lg_parameters.charge[j] / epsilon_in) *
              rec_parameters.charge[i] / d;
          double vdw_val = 0.0;
          if (lg_parameters.sigma[j] < 0. || rec_parameters.sigma[i] < 0.)
            vdw_val = sqrt(16. * lg_parameters.epsilon[j] *
                           rec_parameters.epsilon[i]) *
                      (square(c_6) - 0.000);
          else
            vdw_val = sqrt(16. * lg_parameters.epsilon[j] *
                           rec_parameters.epsilon[i]) *
                      (square(c_6) - c_6);

          CachedData.set_elec_vdw(pos, elec_val, vdw_val);
          Energy += CachedData.get_energy(pos);
        }
      }
    }
  }
  return 0;
}

int cal_one_ele(const MOL &rec, const MOL &lg, cache_sigma_epsilon &CachedData,
                const arguments &Arguments, double pos_x, double pos_y,
                double pos_z, double temp_sigma, double temp_charge,
                double temp_epsilon_lg,
                const Multi_MOL_parameters &rec_parameters,
                const Multi_MOL_parameters &lg_parameters,
                const vector<mol_info> &info_rec,
                const vector<mol_info> &info_lg, cache_energy &Energy, size_t i,
                size_t j) {
  double d{0.1 *
           sqrt(square(rec.pos.x[i] - pos_x) + square(rec.pos.y[i] - pos_y) +
                square(rec.pos.z[i] - pos_z))};
  if (d < 1e-10)
    return 0; // Guard: beads at same position
  double c_6{c6(0.5 * (temp_sigma + rec_parameters[i].get_sigma()) / d)};

  if (Arguments.distance_energy_cutoff > 0. &&
      d > Arguments.distance_energy_cutoff) {
    size_t pos{rec.size() * j + i};
    CachedData.set_elec_vdw(pos, 0., 0.);
    return 0;
  }

  if (Arguments.self_int == true &&
      (d < 0.5 * (temp_sigma + rec_parameters[i].get_sigma()) ||
       (Arguments.flag_interchain_energy == true &&
        info_lg[j].get_chain() == info_rec[i].get_chain()))) {
    return 0;
  }

  if (CachedData.elec_ener.empty()) {
    Energy.elec_ener += temp_charge * rec_parameters[i].get_charge() / d;
    if (temp_sigma < 0. || rec_parameters[i].get_sigma() < 0.)
      Energy.vdw_ener +=
          sqrt(temp_epsilon_lg * rec_parameters[i].get_epsilon()) *
          (square(c_6) - 0.000);
    else
      Energy.vdw_ener +=
          sqrt(temp_epsilon_lg * rec_parameters[i].get_epsilon()) *
          (square(c_6) - c_6);
  } else {
    size_t pos{rec.size() * j + i};
    if (temp_sigma < 0. || rec_parameters[i].get_sigma() < 0.)
      CachedData.set_elec_vdw(
          pos, temp_charge * rec_parameters[i].get_charge() / d,
          sqrt(rec_parameters[i].get_epsilon() * temp_epsilon_lg) *
              (square(c_6) - 0.000));
    else
      CachedData.set_elec_vdw(
          pos, temp_charge * rec_parameters[i].get_charge() / d,
          sqrt(rec_parameters[i].get_epsilon() * temp_epsilon_lg) *
              (square(c_6) - c_6));
    Energy += CachedData.get_energy(pos);
  }
  return 0;
}

int cal_energy_amber_screening_intrin(
    const MOL &rec, const MOL &lg, const vector<mol_info> &info_rec,
    const vector<mol_info> &info_lg, cache_sigma_epsilon &CachedData,
    size_t begin, size_t end, cache_energy &Energy, size_t num_snap,
    const arguments &Arguments, const Multi_MOL_parameters &rec_parameters,
    const Multi_MOL_parameters &lg_parameters) {
  Energy.set_energy(0., 0.);
  double epsilon_in{Arguments.epsilon_in},
      lambda_debye{3.07 / sqrt(Arguments.salt_concentration)};

  __m256d avx_epsilon_in = _mm256_set1_pd(epsilon_in);
  __m256d avx_neg_inv_lambda = _mm256_set1_pd(-1.0 / lambda_debye);
  __m256d avx_one = _mm256_set1_pd(1.0);
  __m256d avx_half = _mm256_set1_pd(0.5);

  for (size_t j{begin}; j < end; ++j) {
    const double pos_x{lg.pos.x[j]}, pos_y{lg.pos.y[j]}, pos_z{lg.pos.z[j]};
    const size_t rec_size{rec.size()}, size_vector_rec{rec.size() / 4};
    double elec_temp[4]{}, vdw_temp_d[4]{}, d[4]{};

    __m256d avx_temp_epsilon_lg =
        _mm256_set1_pd(16. * lg_parameters[j].get_epsilon());
    __m256d avx_temp_charge_lg =
        _mm256_set1_pd(138.935457839 * lg_parameters[j].get_charge());
    __m256d avx_temp_sigma = _mm256_set1_pd(lg_parameters[j].get_sigma());

    __m256d v_pos_x = _mm256_set1_pd(pos_x);
    __m256d v_pos_y = _mm256_set1_pd(pos_y);
    __m256d v_pos_z = _mm256_set1_pd(pos_z);
    __m256d v_scale = _mm256_set1_pd(0.01);

    for (size_t i = 0; i < size_vector_rec * 4; i += 4) {
      __m256d v_x = _mm256_loadu_pd(&rec.pos.x[i]);
      __m256d v_y = _mm256_loadu_pd(&rec.pos.y[i]);
      __m256d v_z = _mm256_loadu_pd(&rec.pos.z[i]);

      __m256d dx = _mm256_sub_pd(v_x, v_pos_x);
      __m256d dy = _mm256_sub_pd(v_y, v_pos_y);
      __m256d dz = _mm256_sub_pd(v_z, v_pos_z);

      __m256d dx2 = _mm256_mul_pd(dx, dx);
      __m256d dy2 = _mm256_mul_pd(dy, dy);
      __m256d dz2 = _mm256_mul_pd(dz, dz);

      __m256d r2 = _mm256_add_pd(_mm256_add_pd(dx2, dy2), dz2);
      r2 = _mm256_mul_pd(r2, v_scale);
      __m256d avx_array_d_rec = _mm256_sqrt_pd(r2);                // d
      __m256d avx_inv_d = _mm256_div_pd(avx_one, avx_array_d_rec); // 1/d

      // Screening term: exp(-d / lambda)
      __m256d avx_exp_arg = _mm256_mul_pd(avx_array_d_rec, avx_neg_inv_lambda);
      __m256d avx_exp_val = f_exp_avx(avx_exp_arg);

      // Elec = C * q1 * q2 * exp / d / eps_in
      __m256d avx_array_charge_rec = _mm256_loadu_pd(&rec_parameters.charge[i]);

      __m256d avx_elec_result =
          _mm256_mul_pd(avx_temp_charge_lg, avx_array_charge_rec);
      avx_elec_result = _mm256_mul_pd(avx_elec_result, avx_exp_val);
      avx_elec_result = _mm256_mul_pd(avx_elec_result, avx_inv_d);
      avx_elec_result = _mm256_div_pd(avx_elec_result, avx_epsilon_in);

      // VdW Calculation
      __m256d vdw_temp = _mm256_loadu_pd(&rec_parameters.epsilon[i]);
      vdw_temp = _mm256_mul_pd(avx_temp_epsilon_lg, vdw_temp);
      vdw_temp = _mm256_sqrt_pd(vdw_temp);

      __m256d vdw_temp2 = _mm256_loadu_pd(&rec_parameters.sigma[i]);
      vdw_temp2 = _mm256_add_pd(avx_temp_sigma, vdw_temp2);
      vdw_temp2 = _mm256_mul_pd(vdw_temp2, avx_half);

      vdw_temp2 = _mm256_div_pd(vdw_temp2, avx_array_d_rec);
      vdw_temp2 = _mm256_mul_pd(vdw_temp2, vdw_temp2);
      __m256d vdw_temp3 = vdw_temp2;
      vdw_temp2 = _mm256_mul_pd(vdw_temp2, vdw_temp2);
      vdw_temp2 = _mm256_mul_pd(vdw_temp2, vdw_temp3); // c6

      vdw_temp3 = _mm256_sub_pd(_mm256_mul_pd(vdw_temp2, vdw_temp2),
                                vdw_temp2); // c12 - c6
      vdw_temp3 = _mm256_mul_pd(vdw_temp, vdw_temp3);

      if (CachedData.elec_ener.empty()) {
        Energy.elec_ener += hsum_double_avx(avx_elec_result);
        Energy.vdw_ener += hsum_double_avx(vdw_temp3);
      } else {
        _mm256_storeu_pd(elec_temp, avx_elec_result);
        _mm256_storeu_pd(vdw_temp_d, vdw_temp3);

        if (Arguments.distance_energy_cutoff < 0.) {
          Energy.elec_ener += hsum_double_avx(avx_elec_result);
          Energy.vdw_ener += hsum_double_avx(vdw_temp3);
          CachedData.set_elec_vdw(rec.size() * j + i, elec_temp, vdw_temp_d, 4);
        } else {
          _mm256_storeu_pd(d, avx_array_d_rec);
          for (size_t k{}; k < 4; ++k) {
            size_t pos{rec.size() * j + i + k};
            if (d[k] > Arguments.distance_energy_cutoff) {
              CachedData.set_elec_vdw(pos, 0., 0.);
              continue;
            }
            CachedData.set_elec_vdw(pos, elec_temp[k], vdw_temp_d[k]);
            Energy.elec_ener += elec_temp[k];
            Energy.vdw_ener += vdw_temp_d[k];
          }
        }
      }
    }

    if (size_vector_rec * 4 < rec_size) {
      for (size_t i{size_vector_rec * 4}; i < rec_size; ++i) {
        double d{0.1 * sqrt(square(rec.pos.x[i] - pos_x) +
                            square(rec.pos.y[i] - pos_y) +
                            square(rec.pos.z[i] - pos_z))};
        double c_6{
            c6(0.5 *
               (lg_parameters[j].get_sigma() + rec_parameters[i].get_sigma()) /
               d)};

        if (Arguments.distance_energy_cutoff > 0. &&
            d > Arguments.distance_energy_cutoff) {
          size_t pos{rec.size() * j + i};
          CachedData.set_elec_vdw(pos, 0., 0.);
          continue;
        }

        double term_exp = exp(-d / lambda_debye);
        double e_val = (138.935457839 * lg_parameters[j].get_charge() *
                        rec_parameters[i].get_charge() * term_exp) /
                       (d * epsilon_in);
        double v_val = sqrt(lg_parameters[j].get_epsilon() * 16. *
                            rec_parameters[i].get_epsilon()) *
                       (square(c_6) - c_6);

        if (CachedData.elec_ener.empty()) {
          Energy.elec_ener += e_val;
          Energy.vdw_ener += v_val;
        } else {
          size_t pos{rec.size() * j + i};
          CachedData.set_elec_vdw(pos, e_val, v_val);
          Energy += CachedData.get_energy(pos);
        }
      }
    }
  }
  return 0;
}

int cal_energy_gromos_intrin(const MOL &rec, const MOL &lg,
                             const vector<mol_info> &info_rec,
                             const vector<mol_info> &info_lg,
                             cache_sigma_epsilon &CachedData, size_t begin,
                             size_t end, cache_energy &Energy, size_t num_snap,
                             const arguments &Arguments,
                             const Multi_MOL_parameters &rec_parameters,
                             const Multi_MOL_parameters &lg_parameters) {
  Energy.set_energy(0., 0.);
  __m256d avx_epsilon_in = _mm256_set1_pd(Arguments.epsilon_in);
  __m256d avx_one = _mm256_set1_pd(1.0);

  for (size_t j{begin}; j < end; ++j) {
    const double pos_x{lg.pos.x[j]}, pos_y{lg.pos.y[j]}, pos_z{lg.pos.z[j]};
    const size_t rec_size{rec.size()}, size_vector_rec{rec.size() / 4};
    double elec_temp[4]{}, vdw_temp_d[4]{}, d[4]{};

    __m256d avx_temp_charge_lg =
        _mm256_set1_pd(138.935457839 * lg_parameters[j].get_charge());
    __m256d avx_temp_epsilon_lg =
        _mm256_set1_pd(lg_parameters[j].get_epsilon());
    __m256d avx_temp_sigma_lg = _mm256_set1_pd(lg_parameters[j].get_sigma());

    __m256d v_pos_x = _mm256_set1_pd(pos_x);
    __m256d v_pos_y = _mm256_set1_pd(pos_y);
    __m256d v_pos_z = _mm256_set1_pd(pos_z);
    __m256d v_scale = _mm256_set1_pd(0.01);

    for (size_t i = 0; i < size_vector_rec * 4; i += 4) {
      __m256d v_x = _mm256_loadu_pd(&rec.pos.x[i]);
      __m256d v_y = _mm256_loadu_pd(&rec.pos.y[i]);
      __m256d v_z = _mm256_loadu_pd(&rec.pos.z[i]);

      __m256d dx = _mm256_sub_pd(v_x, v_pos_x);
      __m256d dy = _mm256_sub_pd(v_y, v_pos_y);
      __m256d dz = _mm256_sub_pd(v_z, v_pos_z);

      __m256d dx2 = _mm256_mul_pd(dx, dx);
      __m256d dy2 = _mm256_mul_pd(dy, dy);
      __m256d dz2 = _mm256_mul_pd(dz, dz);

      __m256d r2 = _mm256_add_pd(_mm256_add_pd(dx2, dy2), dz2);
      r2 = _mm256_mul_pd(r2, v_scale);
      __m256d avx_array_d_rec = _mm256_sqrt_pd(r2);                // d
      __m256d avx_inv_d = _mm256_div_pd(avx_one, avx_array_d_rec); // 1/d

      // Elec = C * q1 * q2 / d / eps_in
      __m256d avx_array_charge_rec = _mm256_loadu_pd(&rec_parameters.charge[i]);
      __m256d avx_elec_result =
          _mm256_mul_pd(avx_temp_charge_lg, avx_array_charge_rec);
      avx_elec_result = _mm256_mul_pd(avx_elec_result, avx_inv_d);
      avx_elec_result = _mm256_div_pd(avx_elec_result, avx_epsilon_in);

      // Gromos VdW: sqrt(eps1*eps2)/r^12 - sqrt(sig1*sig2)/r^6
      __m256d avx_rec_epsilon = _mm256_loadu_pd(&rec_parameters.epsilon[i]);
      __m256d avx_rec_sigma = _mm256_loadu_pd(&rec_parameters.sigma[i]);

      __m256d term_eps =
          _mm256_sqrt_pd(_mm256_mul_pd(avx_temp_epsilon_lg, avx_rec_epsilon));
      __m256d term_sig =
          _mm256_sqrt_pd(_mm256_mul_pd(avx_temp_sigma_lg, avx_rec_sigma));

      __m256d r2_vdw = _mm256_mul_pd(avx_array_d_rec, avx_array_d_rec);
      __m256d r6 = _mm256_mul_pd(r2_vdw, r2_vdw); // r^4
      r6 = _mm256_mul_pd(r6, r2_vdw);             // r^6

      __m256d inv_r6 = _mm256_div_pd(avx_one, r6);
      __m256d inv_r12 = _mm256_mul_pd(inv_r6, inv_r6);

      __m256d vdw_res = _mm256_sub_pd(_mm256_mul_pd(term_eps, inv_r12),
                                      _mm256_mul_pd(term_sig, inv_r6));

      if (CachedData.elec_ener.empty()) {
        Energy.elec_ener += hsum_double_avx(avx_elec_result);
        Energy.vdw_ener += hsum_double_avx(vdw_res);
      } else {
        _mm256_storeu_pd(elec_temp, avx_elec_result);
        _mm256_storeu_pd(vdw_temp_d, vdw_res);
        if (Arguments.distance_energy_cutoff < 0.) {
          Energy.elec_ener += hsum_double_avx(avx_elec_result);
          Energy.vdw_ener += hsum_double_avx(vdw_res);
          CachedData.set_elec_vdw(rec.size() * j + i, elec_temp, vdw_temp_d, 4);
        } else {
          _mm256_storeu_pd(d, avx_array_d_rec);
          for (size_t k{}; k < 4; ++k) {
            size_t pos{rec.size() * j + i + k};
            if (d[k] > Arguments.distance_energy_cutoff) {
              CachedData.set_elec_vdw(pos, 0., 0.);
              continue;
            }
            CachedData.set_elec_vdw(pos, elec_temp[k], vdw_temp_d[k]);
            Energy.elec_ener += elec_temp[k];
            Energy.vdw_ener += vdw_temp_d[k];
          }
        }
      }
    }

    // Remainder
    if (size_vector_rec * 4 < rec_size) {
      for (size_t i{size_vector_rec * 4}; i < rec_size; ++i) {
        double d{0.1 * sqrt(square(rec.pos.x[i] - pos_x) +
                            square(rec.pos.y[i] - pos_y) +
                            square(rec.pos.z[i] - pos_z))};
        if (Arguments.distance_energy_cutoff > 0. &&
            d > Arguments.distance_energy_cutoff) {
          size_t pos{rec.size() * j + i};
          CachedData.set_elec_vdw(pos, 0., 0.);
          continue;
        }

        double temp_ch = 138.935457839 * lg_parameters[j].get_charge() *
                         rec_parameters[i].get_charge() /
                         (d * Arguments.epsilon_in);
        double c_6 = c6(d); // c6(d) = d^6
        double eps_term = sqrt(lg_parameters[j].get_epsilon() *
                               rec_parameters[i].get_epsilon());
        double sig_term =
            sqrt(lg_parameters[j].get_sigma() * rec_parameters[i].get_sigma());
        double vdw_val = eps_term / square(c_6) - sig_term / c_6;

        if (CachedData.elec_ener.empty()) {
          Energy.elec_ener += temp_ch;
          Energy.vdw_ener += vdw_val;
        } else {
          size_t pos{rec.size() * j + i};
          CachedData.set_elec_vdw(pos, temp_ch, vdw_val);
          Energy += CachedData.get_energy(pos);
        }
      }
    }
  }
  return 0;
}

void cal_energy_opls_intrin(const MOL &rec, const MOL &lg,
                            const vector<mol_info> &info_rec,
                            const vector<mol_info> &info_lg,
                            cache_sigma_epsilon &CachedData, size_t begin,
                            size_t end, cache_energy &Energy, size_t num_snap,
                            const arguments &Arguments,
                            const Multi_MOL_parameters &rec_parameters,
                            const Multi_MOL_parameters &lg_parameters) {
  Energy.set_energy(0., 0.);
  __m256d avx_one = _mm256_set1_pd(1.0);

  for (size_t j{begin}; j < end; ++j) {
    const double pos_x{lg.pos.x[j]}, pos_y{lg.pos.y[j]}, pos_z{lg.pos.z[j]};
    const size_t rec_size{rec.size()}, size_vector_rec{rec.size() / 4};
    double elec_temp[4]{}, vdw_temp_d[4]{}, d[4]{};

    __m256d avx_temp_epsilon_lg =
        _mm256_set1_pd(16. * lg_parameters[j].get_epsilon());
    __m256d avx_temp_charge_lg = _mm256_set1_pd(
        138.935457839 * lg_parameters[j].get_charge() / Arguments.epsilon_in);
    __m256d avx_temp_sigma = _mm256_set1_pd(lg_parameters[j].get_sigma());

    __m256d v_pos_x = _mm256_set1_pd(pos_x);
    __m256d v_pos_y = _mm256_set1_pd(pos_y);
    __m256d v_pos_z = _mm256_set1_pd(pos_z);
    __m256d v_scale = _mm256_set1_pd(0.01);

    for (size_t i = 0; i < size_vector_rec * 4; i += 4) {
      __m256d v_x = _mm256_loadu_pd(&rec.pos.x[i]);
      __m256d v_y = _mm256_loadu_pd(&rec.pos.y[i]);
      __m256d v_z = _mm256_loadu_pd(&rec.pos.z[i]);

      __m256d dx = _mm256_sub_pd(v_x, v_pos_x);
      __m256d dy = _mm256_sub_pd(v_y, v_pos_y);
      __m256d dz = _mm256_sub_pd(v_z, v_pos_z);

      __m256d dx2 = _mm256_mul_pd(dx, dx);
      __m256d dy2 = _mm256_mul_pd(dy, dy);
      __m256d dz2 = _mm256_mul_pd(dz, dz);

      __m256d r2 = _mm256_add_pd(_mm256_add_pd(dx2, dy2), dz2);
      r2 = _mm256_mul_pd(r2, v_scale);
      __m256d avx_array_d_rec = _mm256_sqrt_pd(r2);
      __m256d avx_inv_d = _mm256_div_pd(avx_one, avx_array_d_rec);
      __m256d avx_inv_d2 = _mm256_mul_pd(avx_inv_d, avx_inv_d);

      __m256d avx_array_charge_rec = _mm256_loadu_pd(&rec_parameters.charge[i]);
      __m256d avx_elec_result =
          _mm256_mul_pd(avx_temp_charge_lg, avx_array_charge_rec);
      avx_elec_result = _mm256_mul_pd(avx_elec_result, avx_inv_d);

      __m256d vdw_temp = _mm256_loadu_pd(&rec_parameters.epsilon[i]);
      vdw_temp = _mm256_sqrt_pd(_mm256_mul_pd(avx_temp_epsilon_lg, vdw_temp));

      __m256d vdw_temp2 = _mm256_loadu_pd(&rec_parameters.sigma[i]);
      vdw_temp2 = _mm256_mul_pd(avx_temp_sigma, vdw_temp2);
      vdw_temp2 = _mm256_sqrt_pd(vdw_temp2); // sqrt(sigma1 * sigma2)

      // Calculate (sigma/d)^6 using d^2
      vdw_temp2 = _mm256_mul_pd(vdw_temp2, vdw_temp2);  // sigma_comb^2
      vdw_temp2 = _mm256_mul_pd(vdw_temp2, avx_inv_d2); // (sigma/d)^2

      __m256d vdw_temp3 = vdw_temp2;
      vdw_temp2 = _mm256_mul_pd(vdw_temp2, vdw_temp2); // ^4
      vdw_temp2 = _mm256_mul_pd(vdw_temp2, vdw_temp3); // ^6 (c6)

      vdw_temp3 = _mm256_sub_pd(_mm256_mul_pd(vdw_temp2, vdw_temp2), vdw_temp2);
      vdw_temp3 = _mm256_mul_pd(vdw_temp, vdw_temp3); // eps * (c12 - c6)

      if (CachedData.elec_ener.empty()) {
        Energy.elec_ener += hsum_double_avx(avx_elec_result);
        Energy.vdw_ener += hsum_double_avx(vdw_temp3);
      } else {
        _mm256_storeu_pd(elec_temp, avx_elec_result);
        _mm256_storeu_pd(vdw_temp_d, vdw_temp3);
        if (Arguments.distance_energy_cutoff < 0.) {
          Energy.elec_ener += hsum_double_avx(avx_elec_result);
          Energy.vdw_ener += hsum_double_avx(vdw_temp3);
          CachedData.set_elec_vdw(rec.size() * j + i, elec_temp, vdw_temp_d, 4);
        } else {
          _mm256_storeu_pd(d, avx_array_d_rec);
          for (size_t k{}; k < 4; ++k) {
            size_t pos{rec.size() * j + i + k};
            if (d[k] > Arguments.distance_energy_cutoff) {
              CachedData.set_elec_vdw(pos, 0., 0.);
              continue;
            }
            const double temp_sigma{lg_parameters[j].get_sigma()};
            if (Arguments.self_int == true &&
                (c6(sqrt(temp_sigma * rec_parameters[i + k].get_sigma()) /
                    d[k]) > 1. ||
                 (Arguments.flag_interchain_energy == true &&
                  info_lg[j].get_chain() == info_rec[i + k].get_chain()))) {
              CachedData.set_elec_vdw(pos, 0., 0.);
              continue;
            }

            CachedData.set_elec_vdw(pos, elec_temp[k], vdw_temp_d[k]);

            Energy.elec_ener += elec_temp[k];
            Energy.vdw_ener += vdw_temp_d[k];
          }
        }
      }
    }

    if (size_vector_rec * 4 < rec_size)
      for (size_t i{size_vector_rec * 4}; i < rec_size; ++i) {
        double d{0.1 * sqrt(square(rec.pos.x[i] - pos_x) +
                            square(rec.pos.y[i] - pos_y) +
                            square(rec.pos.z[i] - pos_z))},
            c_6{c6(sqrt(lg_parameters[j].get_sigma() *
                        rec_parameters[i].get_sigma()) /
                   d)};

        if (Arguments.self_int == true &&
            (c_6 > 1. || (Arguments.flag_interchain_energy == true &&
                          info_lg[j].get_chain() == info_rec[i].get_chain()))) {
          continue;
        }

        if (Arguments.distance_energy_cutoff > 0. &&
            d > Arguments.distance_energy_cutoff) {
          size_t pos{rec.size() * j + i};
          CachedData.set_elec_vdw(pos, 0., 0.);
          continue;
        }

        double e_val = 138.935457839 * lg_parameters[j].get_charge() *
                       rec_parameters[i].get_charge() /
                       (d * Arguments.epsilon_in);
        double v_val = sqrt(lg_parameters[j].get_epsilon() * 16. *
                            rec_parameters[i].get_epsilon()) *
                       (square(c_6) - c_6);

        if (CachedData.elec_ener.empty()) {
          Energy.elec_ener += e_val;
          Energy.vdw_ener += v_val;
        } else {
          size_t pos{rec.size() * j + i};
          CachedData.set_elec_vdw(pos, e_val, v_val);
          Energy += CachedData.get_energy(pos);
        }
      }
  }
}

// MARTINI: AVX2 energy kernel using the topology-derived LJ pair matrix and
// Martini reaction-field electrostatics. This stays MARTINI-only and does not
// alter all-atom kernels.
int cal_energy_martini_intrin(const MOL &rec, const MOL &lg,
                              const vector<mol_info> &info_rec,
                              const vector<mol_info> &info_lg,
                              cache_sigma_epsilon &CachedData, size_t begin,
                              size_t end, cache_energy &Energy, size_t num_snap,
                              const arguments &Arguments,
                              const Multi_MOL_parameters &rec_parameters,
                              const Multi_MOL_parameters &lg_parameters) {
  Energy.set_energy(0., 0.);
  const auto elec_cfg = martini::electrostatic_config(Arguments);
  const __m256d v_scale = _mm256_set1_pd(0.01);
  const __m256d avx_one = _mm256_set1_pd(1.0);
  const __m256d avx_zero = _mm256_setzero_pd();
  const __m256d avx_neg_one = _mm256_set1_pd(-1.0);
  const __m256d avx_half = _mm256_set1_pd(0.5);
  const __m256d avx_four = _mm256_set1_pd(4.0);
  const __m256d avx_charge_scale =
      _mm256_set1_pd(martini::kCoulomb / elec_cfg.epsilon_r);
  const __m256d avx_k_rf = _mm256_set1_pd(elec_cfg.k_rf);
  const __m256d avx_c_rf = _mm256_set1_pd(elec_cfg.c_rf);
  const __m256d avx_min_d = _mm256_set1_pd(1e-10);
  const size_t rec_size = rec.size();
  const size_t size_vector_rec = rec_size / 4;
  const bool write_cache = !CachedData.elec_ener.empty();
  const bool has_pair_table = CachedData.has_martini_pair_table();
  const int martini_types = CachedData.martini_types();
  const double *pair_sigma = CachedData.martini_sigma_ptr();
  const double *pair_epsilon = CachedData.martini_epsilon_ptr();
  const double *rec_x = rec.pos.x.data();
  const double *rec_y = rec.pos.y.data();
  const double *rec_z = rec.pos.z.data();
  const double *rec_charge = rec_parameters.charge.data();
  const double *rec_sigma = rec_parameters.sigma.data();
  const double *rec_epsilon = rec_parameters.epsilon.data();
  const int32_t *rec_type = rec_parameters.type_id_ptr();

  for (size_t j{begin}; j < end; ++j) {
    const double pos_x{lg.pos.x[j]}, pos_y{lg.pos.y[j]}, pos_z{lg.pos.z[j]};
    const int lg_type = lg_parameters.get_type_id(j);
    const __m256d avx_temp_charge_lg = _mm256_mul_pd(
        avx_charge_scale, _mm256_set1_pd(lg_parameters.charge[j]));
    const __m256d avx_lg_sigma = _mm256_set1_pd(lg_parameters.sigma[j]);
    const __m256d avx_lg_eps = _mm256_set1_pd(lg_parameters.epsilon[j]);
    const __m256d v_pos_x = _mm256_set1_pd(pos_x);
    const __m256d v_pos_y = _mm256_set1_pd(pos_y);
    const __m256d v_pos_z = _mm256_set1_pd(pos_z);
    const int lg_row_offset = lg_type >= 0 ? lg_type * martini_types : 0;

    for (size_t i = 0; i < size_vector_rec * 4; i += 4) {
      const __m256d v_x = _mm256_load_pd(&rec_x[i]);
      const __m256d v_y = _mm256_load_pd(&rec_y[i]);
      const __m256d v_z = _mm256_load_pd(&rec_z[i]);

      const __m256d dx = _mm256_sub_pd(v_x, v_pos_x);
      const __m256d dy = _mm256_sub_pd(v_y, v_pos_y);
      const __m256d dz = _mm256_sub_pd(v_z, v_pos_z);

      const __m256d dx2 = _mm256_mul_pd(dx, dx);
      const __m256d dy2 = _mm256_mul_pd(dy, dy);
      const __m256d dz2 = _mm256_mul_pd(dz, dz);

      __m256d r2 = _mm256_add_pd(_mm256_add_pd(dx2, dy2), dz2);
      r2 = _mm256_mul_pd(r2, v_scale);
      __m256d avx_d = _mm256_sqrt_pd(r2);
      const __m256d avx_d_mask = _mm256_cmp_pd(avx_d, avx_min_d, _CMP_GT_OQ);
      __m256d avx_d_safe =
          _mm256_blendv_pd(avx_one, avx_d, avx_d_mask);

      __m256d sigma_vec = avx_neg_one;
      __m256d epsilon_vec = avx_neg_one;
      if (has_pair_table && lg_type >= 0) {
        alignas(16) int pair_idx[4];
        alignas(32) double lookup_mask_vals[4];
        for (int k = 0; k < 4; ++k) {
          const int rec_type_lane = rec_type ? rec_type[i + k] : -1;
          if (rec_type_lane >= 0) {
            pair_idx[k] = lg_row_offset + rec_type_lane;
            lookup_mask_vals[k] = -1.0;
          } else {
            pair_idx[k] = 0;
            lookup_mask_vals[k] = 0.0;
          }
        }
#ifdef __AVX2__
        const __m128i gather_idx =
            _mm_loadu_si128(reinterpret_cast<const __m128i *>(pair_idx));
        sigma_vec = _mm256_i32gather_pd(pair_sigma, gather_idx, 8);
        epsilon_vec = _mm256_i32gather_pd(pair_epsilon, gather_idx, 8);
#else
        alignas(32) double sigma_lane[4], epsilon_lane[4];
        for (int k = 0; k < 4; ++k) {
          sigma_lane[k] = pair_sigma[pair_idx[k]];
          epsilon_lane[k] = pair_epsilon[pair_idx[k]];
        }
        sigma_vec = _mm256_loadu_pd(sigma_lane);
        epsilon_vec = _mm256_loadu_pd(epsilon_lane);
#endif
        const __m256d lookup_mask = _mm256_loadu_pd(lookup_mask_vals);
        sigma_vec = _mm256_blendv_pd(avx_neg_one, sigma_vec, lookup_mask);
        epsilon_vec = _mm256_blendv_pd(avx_neg_one, epsilon_vec, lookup_mask);
      }

      __m256d fallback_sigma =
          _mm256_mul_pd(_mm256_add_pd(avx_lg_sigma, _mm256_load_pd(&rec_sigma[i])),
                        avx_half);
      __m256d fallback_epsilon = _mm256_sqrt_pd(
          _mm256_mul_pd(avx_lg_eps, _mm256_load_pd(&rec_epsilon[i])));
      __m256d sigma_missing =
          _mm256_cmp_pd(sigma_vec, avx_zero, _CMP_LT_OQ);
      __m256d epsilon_missing =
          _mm256_cmp_pd(epsilon_vec, avx_zero, _CMP_LT_OQ);
      sigma_vec = _mm256_blendv_pd(sigma_vec, fallback_sigma, sigma_missing);
      epsilon_vec =
          _mm256_blendv_pd(epsilon_vec, fallback_epsilon, epsilon_missing);

      __m256d inv_d = _mm256_div_pd(avx_one, avx_d_safe);
      __m256d avx_elec_result =
          _mm256_mul_pd(avx_temp_charge_lg, _mm256_load_pd(&rec_charge[i]));
      __m256d rf_term = _mm256_sub_pd(
          _mm256_add_pd(inv_d, _mm256_mul_pd(avx_k_rf, r2)), avx_c_rf);
      avx_elec_result = _mm256_mul_pd(avx_elec_result, rf_term);

      __m256d sr = _mm256_div_pd(sigma_vec, avx_d_safe);
      __m256d sr2 = _mm256_mul_pd(sr, sr);
      __m256d sr4 = _mm256_mul_pd(sr2, sr2);
      __m256d sr6 = _mm256_mul_pd(sr4, sr2);
      __m256d sr12 = _mm256_mul_pd(sr6, sr6);
      __m256d vdw_result =
          _mm256_mul_pd(_mm256_mul_pd(avx_four, epsilon_vec),
                        _mm256_sub_pd(sr12, sr6));

      alignas(32) double d_arr[4], sigma_arr[4], lane_mask_vals[4];
      _mm256_storeu_pd(d_arr, avx_d);
      _mm256_storeu_pd(sigma_arr, sigma_vec);

      for (int k = 0; k < 4; ++k) {
        bool keep = d_arr[k] > 1e-10;
        if (keep)
          keep = d_arr[k] <= elec_cfg.cutoff_nm;
        if (keep && Arguments.self_int == true) {
          keep = d_arr[k] >= sigma_arr[k];
          if (keep && Arguments.flag_interchain_energy == true)
            keep = info_lg[j].get_chain() != info_rec[i + k].get_chain();
        }
        lane_mask_vals[k] = keep ? -1.0 : 0.0;
      }

      const __m256d final_mask = _mm256_loadu_pd(lane_mask_vals);
      avx_elec_result = _mm256_and_pd(avx_elec_result, final_mask);
      vdw_result = _mm256_and_pd(vdw_result, final_mask);

      Energy.elec_ener += hsum_double_avx(avx_elec_result);
      Energy.vdw_ener += hsum_double_avx(vdw_result);

      if (write_cache) {
        alignas(32) double elec_temp[4], vdw_temp_d[4];
        _mm256_storeu_pd(elec_temp, avx_elec_result);
        _mm256_storeu_pd(vdw_temp_d, vdw_result);
        const size_t pos_base = rec_size * j + i;
        for (int k = 0; k < 4; ++k)
          CachedData.set_elec_vdw(pos_base + k, elec_temp[k], vdw_temp_d[k]);
      }
    }

    if (size_vector_rec * 4 < rec_size) {
      for (size_t i{size_vector_rec * 4}; i < rec_size; ++i) {
        double d{0.1 * sqrt(square(rec.pos.x[i] - pos_x) +
                            square(rec.pos.y[i] - pos_y) +
                            square(rec.pos.z[i] - pos_z))};
        size_t pos{rec_size * j + i};
        if (d < 1e-10) {
          if (write_cache)
            CachedData.set_elec_vdw(pos, 0., 0.);
          continue;
        }

        double sigma_ij{}, epsilon_ij{};
        martini::pair_params(CachedData, rec_parameters, lg_parameters, i, j,
                             sigma_ij, epsilon_ij);
        if (!martini::within_cutoff(elec_cfg, d)) {
          if (write_cache)
            CachedData.set_elec_vdw(pos, 0., 0.);
          continue;
        }
        if (martini::excluded_pair(Arguments, info_rec, info_lg, i, j, d,
                                   sigma_ij)) {
          if (write_cache)
            CachedData.set_elec_vdw(pos, 0., 0.);
          continue;
        }

        const double elec_val = martini::electro_energy(
            elec_cfg, rec_parameters.charge[i], lg_parameters.charge[j], d);
        const double vdw_val = martini::lj_energy(sigma_ij, epsilon_ij, d);

        if (write_cache)
          CachedData.set_elec_vdw(pos, elec_val, vdw_val);
        Energy.elec_ener += elec_val;
        Energy.vdw_ener += vdw_val;
      }
    }
  }
  return 0;
}
#endif
