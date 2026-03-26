#include "cache_sigma_epsilon.h"
#include "force.h"
#include "functions.h"
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

// MARTINI: Lorentz-Berthelot arithmetic mean for sigma, geometric mean for
// epsilon via sqrt(16*e1 * e2) = 4*sqrt(e1*e2). Sigma in nm, epsilon in
// kJ/mol. epsilon_in should be set to ~15.0 (Martini relative permittivity).
// This function is intentionally identical to cal_energy_amber_intrin.
int cal_energy_martini_intrin(const MOL &rec, const MOL &lg,
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

        // Guard: mask lanes where d < 1e-10 to avoid 0/0 = NaN
        __m256d avx_d_eps = _mm256_set1_pd(1e-10);
        __m256d avx_d_mask =
            _mm256_cmp_pd(avx_array_d_rec, avx_d_eps, _CMP_GT_OQ);
        __m256d avx_d_safe =
            _mm256_blendv_pd(_mm256_set1_pd(1.0), avx_array_d_rec, avx_d_mask);

        __m256d avx_array_charge_rec =
            _mm256_loadu_pd(&rec_parameters.charge[i]);
        __m256d avx_elec_result =
            _mm256_mul_pd(avx_temp_charge_lg, avx_array_charge_rec);
        avx_elec_result = _mm256_div_pd(avx_elec_result, avx_d_safe);
        avx_elec_result = _mm256_and_pd(avx_elec_result, avx_d_mask);

        __m256d vdw_temp = _mm256_loadu_pd(&rec_parameters.epsilon[i]);
        vdw_temp = _mm256_mul_pd(avx_temp_epsilon_lg, vdw_temp);
        vdw_temp = _mm256_sqrt_pd(vdw_temp);

        __m256d vdw_temp2 = _mm256_loadu_pd(&rec_parameters.sigma[i]);
        vdw_temp2 = _mm256_add_pd(avx_temp_sigma, vdw_temp2);
        vdw_temp2 = _mm256_mul_pd(vdw_temp2, avx_half);

        vdw_temp2 = _mm256_div_pd(vdw_temp2, avx_d_safe);
        vdw_temp2 = _mm256_mul_pd(vdw_temp2, vdw_temp2);
        __m256d vdw_temp3 = vdw_temp2;
        vdw_temp2 = _mm256_mul_pd(vdw_temp2, vdw_temp2);
        vdw_temp2 = _mm256_mul_pd(vdw_temp2, vdw_temp3);
        vdw_temp3 =
            _mm256_sub_pd(_mm256_mul_pd(vdw_temp2, vdw_temp2), vdw_temp2);
        vdw_temp3 = _mm256_mul_pd(vdw_temp, vdw_temp3);
        vdw_temp3 = _mm256_and_pd(vdw_temp3, avx_d_mask); // zero out d=0 lanes

        Energy.elec_ener += hsum_double_avx(avx_elec_result);
        Energy.vdw_ener += hsum_double_avx(vdw_temp3);
      }

      if (size_vector_rec * 4 < rec_size) {
        for (size_t i{size_vector_rec * 4}; i < rec_size; ++i) {
          double d{0.1 * sqrt(square(rec.pos.x[i] - pos_x) +
                              square(rec.pos.y[i] - pos_y) +
                              square(rec.pos.z[i] - pos_z))};
          if (d < 1e-10)
            continue; // Guard: beads at same position
          double c_6{
              c6(0.5 * (lg_parameters.sigma[j] + rec_parameters.sigma[i]) / d)};
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

        // Guard: mask lanes where d < 1e-10 to avoid 0/0 = NaN
        __m256d avx_d_eps = _mm256_set1_pd(1e-10);
        __m256d avx_d_mask =
            _mm256_cmp_pd(avx_array_d_rec, avx_d_eps, _CMP_GT_OQ);
        __m256d avx_d_safe =
            _mm256_blendv_pd(_mm256_set1_pd(1.0), avx_array_d_rec, avx_d_mask);

        __m256d avx_array_charge_rec =
            _mm256_loadu_pd(&rec_parameters.charge[i]);
        __m256d avx_elec_result =
            _mm256_mul_pd(avx_temp_charge_lg, avx_array_charge_rec);
        avx_elec_result = _mm256_div_pd(avx_elec_result, avx_d_safe);
        avx_elec_result = _mm256_and_pd(avx_elec_result, avx_d_mask);

        __m256d vdw_temp = _mm256_loadu_pd(&rec_parameters.epsilon[i]);
        vdw_temp = _mm256_mul_pd(avx_temp_epsilon_lg, vdw_temp);
        vdw_temp = _mm256_sqrt_pd(vdw_temp);

        __m256d vdw_temp2 = _mm256_loadu_pd(&rec_parameters.sigma[i]);
        vdw_temp2 = _mm256_add_pd(avx_temp_sigma, vdw_temp2);
        vdw_temp2 = _mm256_mul_pd(vdw_temp2, avx_half);

        vdw_temp2 = _mm256_div_pd(vdw_temp2, avx_d_safe);
        vdw_temp2 = _mm256_mul_pd(vdw_temp2, vdw_temp2);
        __m256d vdw_temp3 = vdw_temp2;
        vdw_temp2 = _mm256_mul_pd(vdw_temp2, vdw_temp2);
        vdw_temp2 = _mm256_mul_pd(vdw_temp2, vdw_temp3);
        vdw_temp3 =
            _mm256_sub_pd(_mm256_mul_pd(vdw_temp2, vdw_temp2), vdw_temp2);
        vdw_temp3 = _mm256_mul_pd(vdw_temp, vdw_temp3);
        vdw_temp3 = _mm256_and_pd(vdw_temp3, avx_d_mask); // zero out d=0 lanes

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
                              square(rec.pos.z[i] - pos_z))};
          if (d < 1e-10) {
            CachedData.set_elec_vdw(rec_size * j + i, 0., 0.);
            continue;
          }
          double c_6{
              c6(0.5 * (lg_parameters.sigma[j] + rec_parameters.sigma[i]) / d)};
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
#endif
