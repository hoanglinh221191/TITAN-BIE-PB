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
                const vector<Multi_MOL_parameters> &rec_parameters,
                const vector<Multi_MOL_parameters> &lg_parameters,
                const vector<mol_info> &info_rec,
                const vector<mol_info> &info_lg, cache_energy &Energy, size_t i,
                size_t j);

int cal_energy_amber_intrin(const MOL &rec, const MOL &lg,
                            const vector<mol_info> &info_rec,
                            const vector<mol_info> &info_lg,
                            cache_sigma_epsilon &CachedData, size_t begin,
                            size_t end, cache_energy &Energy, size_t num_snap,
                            const arguments &Arguments,
                            const vector<Multi_MOL_parameters> &rec_parameters,
                            const vector<Multi_MOL_parameters> &lg_parameters) {
  Energy.set_energy(0., 0.);
  for (size_t j{begin}; j < end; ++j) {
    const double pos_x{lg.pos.x[j]}, pos_y{lg.pos.y[j]}, pos_z{lg.pos.z[j]};
    const size_t rec_size{rec.size()}, size_vector_rec{rec.size() / 4};
    double elec_temp[4]{}, vdw_temp_d[4]{}, d[4]{};

    __m256d avx_temp_epsilon_lg =
        _mm256_set1_pd(16. * lg_parameters[j].get_epsilon());
    __m256d avx_temp_charge_lg = _mm256_set1_pd(
        138.935457839 * lg_parameters[j].get_charge() / Arguments.epsilon_in);
    __m256d avx_temp_sigma = _mm256_set1_pd(lg_parameters[j].get_sigma());
    __m256d avx_half = _mm256_set1_pd(0.5);

    for (size_t i{}; i <= (size_vector_rec - 1) * 4; i += 4) {
      __m256d avx_array_d_rec = _mm256_set_pd(
          0.01 * (square(rec.pos.x[i] - pos_x) + square(rec.pos.y[i] - pos_y) +
                  square(rec.pos.z[i] - pos_z)),
          0.01 * (square(rec.pos.x[i + 1] - pos_x) + square(rec.pos.y[i + 1] - pos_y) +
                  square(rec.pos.z[i + 1] - pos_z)),
          0.01 * (square(rec.pos.x[i + 2] - pos_x) + square(rec.pos.y[i + 2] - pos_y) +
                  square(rec.pos.z[i + 2] - pos_z)),
          0.01 * (square(rec.pos.x[i + 3] - pos_x) + square(rec.pos.y[i + 3] - pos_y) +
                  square(rec.pos.z[i + 3] - pos_z)));

      avx_array_d_rec = _mm256_sqrt_pd(avx_array_d_rec);

      __m256d avx_array_charge_rec = _mm256_set_pd(
          rec_parameters[i].get_charge(), rec_parameters[i + 1].get_charge(),
          rec_parameters[i + 2].get_charge(),
          rec_parameters[i + 3].get_charge());
      __m256d avx_elec_result =
          _mm256_mul_pd(avx_temp_charge_lg, avx_array_charge_rec);
      avx_elec_result = _mm256_div_pd(avx_elec_result, avx_array_d_rec);

      __m256d vdw_temp = _mm256_set_pd(rec_parameters[i].get_epsilon(),
                                       rec_parameters[i + 1].get_epsilon(),
                                       rec_parameters[i + 2].get_epsilon(),
                                       rec_parameters[i + 3].get_epsilon());
      vdw_temp = _mm256_mul_pd(avx_temp_epsilon_lg, vdw_temp);
      vdw_temp = _mm256_sqrt_pd(vdw_temp);

      __m256d vdw_temp2 = _mm256_set_pd(
          rec_parameters[i].get_sigma(), rec_parameters[i + 1].get_sigma(),
          rec_parameters[i + 2].get_sigma(), rec_parameters[i + 3].get_sigma());
      vdw_temp2 = _mm256_add_pd(avx_temp_sigma, vdw_temp2);
      vdw_temp2 = _mm256_mul_pd(vdw_temp2, avx_half);

      vdw_temp2 = _mm256_div_pd(vdw_temp2, avx_array_d_rec);
      vdw_temp2 = _mm256_mul_pd(vdw_temp2, vdw_temp2);
      __m256d vdw_temp3 = vdw_temp2;
      vdw_temp2 = _mm256_mul_pd(vdw_temp2, vdw_temp2);
      vdw_temp2 = _mm256_mul_pd(vdw_temp2, vdw_temp3);
      // vdw_temp3 = _mm256_fmsub_pd(vdw_temp2, vdw_temp2, vdw_temp2);
      vdw_temp3 = _mm256_sub_pd(_mm256_mul_pd(vdw_temp2, vdw_temp2), vdw_temp2);
      vdw_temp3 = _mm256_mul_pd(vdw_temp, vdw_temp3);

      if (CachedData.elec_ener.empty()) {
        Energy.elec_ener += hsum_double_avx(avx_elec_result);
        Energy.vdw_ener += hsum_double_avx(vdw_temp3);
      } else {
        double elec_arr[4], vdw_arr[4];
        _mm256_storeu_pd(elec_arr, avx_elec_result);
        _mm256_storeu_pd(vdw_arr, vdw_temp3);
        size_t pos = rec_size * j + i;
        for (int k = 0; k < 4; ++k) {
          CachedData.set_elec_vdw(pos + 3 - k, elec_arr[k], vdw_arr[k]);
        }
        Energy.elec_ener += hsum_double_avx(avx_elec_result);
        Energy.vdw_ener += hsum_double_avx(vdw_temp3);
      }
    }

    if (size_vector_rec * 4 < rec_size)
      for (size_t i{size_vector_rec * 4}; i < rec_size; ++i) {
        const double temp_charge{138.935457839 * lg_parameters[j].get_charge() /
                                 Arguments.epsilon_in},
            temp_epsilon_lg{16. * lg_parameters[j].get_epsilon()},
            temp_sigma{lg_parameters[j].get_sigma()};
        cal_one_ele(rec, lg, CachedData, Arguments, pos_x, pos_y, pos_z,
                    temp_sigma, temp_charge, temp_epsilon_lg, rec_parameters,
                    lg_parameters, info_rec, info_lg, Energy, i, j);
      }
  }
  return 0;
}

int cal_one_ele(const MOL &rec, const MOL &lg, cache_sigma_epsilon &CachedData,
                const arguments &Arguments, double pos_x, double pos_y,
                double pos_z, double temp_sigma, double temp_charge,
                double temp_epsilon_lg,
                const vector<Multi_MOL_parameters> &rec_parameters,
                const vector<Multi_MOL_parameters> &lg_parameters,
                const vector<mol_info> &info_rec,
                const vector<mol_info> &info_lg, cache_energy &Energy, size_t i,
                size_t j) {
  double d{0.1 * sqrt(square(rec.pos.x[i] - pos_x) + square(rec.pos.y[i] - pos_y) +
                      square(rec.pos.z[i] - pos_z))},
      c_6{c6(0.5 * (temp_sigma + rec_parameters[i].get_sigma()) / d)};

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
#endif
