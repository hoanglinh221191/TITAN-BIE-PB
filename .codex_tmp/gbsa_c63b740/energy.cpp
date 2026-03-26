#include "cache_sigma_epsilon.h"
#include "force.h"
#include "functions.h"
#include "miscellaneous.h"
#include "new_emap.h"
#include "stat_ener_res_res.h"
#include "uti_func.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <set>
#include <string>
#include <string_view>
using std::string_view;
#include <thread>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <vector>

#ifdef __AVX__
#include "avx_math.h"
#include <immintrin.h>
#else
#warning AVX is not available.Code will not compile!
#endif
#include <functional>

// AVX is not available.

using namespace std;
int imp_cal_energy_lg(vector<data_out> &lg_ener_inter, const Multi_MOL &rec,
                      const Multi_MOL &lg, size_t num_snap, size_t begin,
                      size_t end, const vector<d_min> &dmin,
                      const cache_sigma_epsilon &CachedData, bool flag_all_mode,
                      const Multi_MOL_parameters &rec_parameters,
                      const Multi_MOL_parameters &lg_parameters,
                      unordered_map<pair<size_t, size_t>, pair<size_t, size_t>,
                                    pair_hash> &map_list_res,
                      vector<pair<size_t, size_t>> &list_res,
                      info_ener_res &infoener_res, const arguments &Arguments);
int imp_cal_energy_rec(
    vector<data_out> &rec_ener_inter, info_ener_res &infoener_res,
    const Multi_MOL &rec, const Multi_MOL &lg, size_t num_snap, size_t begin,
    size_t end, const vector<d_min> &dmin,
    const cache_sigma_epsilon &CachedData, bool flag_all_mode,
    const Multi_MOL_parameters &rec_parameters,
    const Multi_MOL_parameters &lg_parameters,
    unordered_map<pair<size_t, size_t>, pair<size_t, size_t>, pair_hash>
        &map_list_res,
    vector<pair<size_t, size_t>> &list_res, const arguments &Arguments,
    const vector<bool> &is_lg_dmin_res);
void print_energy_map(size_t nt, vector<vector<data_out>> &data,
                      fstream &energy_file, bool flag_all_mode,
                      const arguments &Arguments,
                      vector<energy_res> &Energy_res, size_t num_snap,
                      info &info_output);
// bool check_in_energy_list(vector<data_out>& lg_ener_inter, string chain,
// size_t resnr);

inline int add_energy(const cache_sigma_epsilon &CachedData, size_t i,
                      double &electro_energy, double &van_der_waals_energy,
                      const Multi_MOL_parameters &rec_parameters,
                      const Multi_MOL_parameters &lg_parameters, const MOL &rec,
                      const MOL &lg, size_t snap, const arguments &Arguments,
                      const vector<mol_info> &info_rec,
                      const vector<mol_info> &info_lg, size_t rec_pos = 0,
                      size_t lg_pos = 0) {
  if (i < CachedData.elec_ener.size()) {
    electro_energy += CachedData.elec_ener[i];
    van_der_waals_energy += CachedData.vdw_ener[i];

    if (CachedData.vdw_ener[i] > 500.0) {
      // Intentionally removed I/O inside hot loop
    }
  } else {
    if (rec_pos >= rec.size() || lg_pos >= lg.size())
      return 1;
    double dx{rec.pos.x[rec_pos] - lg.pos.x[lg_pos]},
        dy{rec.pos.y[rec_pos] - lg.pos.y[lg_pos]},
        dz{rec.pos.z[rec_pos] - lg.pos.z[lg_pos]},
        charge{138.935457839 * rec_parameters[rec_pos].get_charge() *
               lg_parameters[lg_pos].get_charge() / Arguments.epsilon_in},
        d{0.1 * sqrt(square(dx) + square(dy) + square(dz))};

    if (CachedData.cutoff_dis > 0. && d > CachedData.cutoff_dis)
      return 1;

    electro_energy += charge / d;

    double vdw_initial = van_der_waals_energy;
    if (CachedData.force_field_name.find("AMBER") != string::npos ||
        CachedData.force_field_name.find("CHARMM") != string::npos) {
      double temp_vdw{sqrt(rec_parameters[rec_pos].get_epsilon() *
                           lg_parameters[lg_pos].get_epsilon())},
          c_6{c6(0.5 *
                 (lg_parameters[lg_pos].get_sigma() +
                  rec_parameters[rec_pos].get_sigma()) /
                 d)},
          c_12{square(c_6)}, sd{square(1. / d)},
          temp_par_vdw{temp_vdw * (12. * c_12 * sd - 6. * c_6 * sd)};

      if (lg_parameters[lg_pos].get_sigma() < 0. ||
          rec_parameters[rec_pos].get_sigma() < 0.)
        c_6 = 0.000;

      if (Arguments.self_int == true &&
          (c_6 > 1. ||
           (Arguments.flag_interchain_energy == true &&
            info_lg[lg_pos].get_chain() == info_rec[rec_pos].get_chain())))
        return 1;

      van_der_waals_energy += temp_vdw * (c_12 - c_6);
    } else if (CachedData.force_field_name.find("OPLS") != string::npos) {
      double temp_vdw{sqrt(rec_parameters[rec_pos].get_epsilon() *
                           lg_parameters[lg_pos].get_epsilon())},
          c_6{c6(sqrt(lg_parameters[lg_pos].get_sigma() *
                      rec_parameters[rec_pos].get_sigma()) /
                 d)},
          c_12{square(c_6)}, sd{square(1. / d)},
          temp_par_vdw{temp_vdw * (12. * c_12 * sd - 6. * c_6 * sd)};

      if (lg_parameters[lg_pos].get_sigma() < 0. ||
          rec_parameters[rec_pos].get_sigma() < 0.)
        c_6 = 0.000;

      if (Arguments.self_int == true &&
          (c_6 > 1. ||
           (Arguments.flag_interchain_energy == true &&
            info_lg[lg_pos].get_chain() == info_rec[rec_pos].get_chain())))
        return 1;

      van_der_waals_energy += temp_vdw * (c_12 - c_6);
    } else if (CachedData.force_field_name.find("MARTINI") != string::npos) {
      // MARTINI: Lorentz-Berthelot arithmetic mean for sigma (same as AMBER)
      double sigma_ij = 0.5 * (lg_parameters[lg_pos].get_sigma() +
                               rec_parameters[rec_pos].get_sigma());
      double c_6 = c6(sigma_ij / d);
      double c_12 = square(c_6);

      if (lg_parameters[lg_pos].get_sigma() < 0. ||
          rec_parameters[rec_pos].get_sigma() < 0.)
        c_6 = 0.0;

      if (Arguments.self_int == true &&
          (d < sigma_ij ||
           (Arguments.flag_interchain_energy == true &&
            info_lg[lg_pos].get_chain() == info_rec[rec_pos].get_chain())))
        return 1;

      double temp_vdw = sqrt(rec_parameters[rec_pos].get_epsilon() *
                             lg_parameters[lg_pos].get_epsilon());
      van_der_waals_energy += temp_vdw * (c_12 - c_6);
      // NOTE: electro_energy already accumulated above — do NOT add it here.
    } else {
      double sd{1. / square(d)}, temp_elec{sd * charge / Arguments.epsilon_in},
          c_6{c6(d)}, c_12{square(c_6)};
      double epsilon{sqrt(rec_parameters[rec_pos].get_epsilon() *
                          lg_parameters[lg_pos].get_epsilon())},
          sigma{sqrt(rec_parameters[rec_pos].get_sigma() *
                     lg_parameters[lg_pos].get_sigma())},
          temp_vdw{12. * sd * epsilon / c_12 - 6. * sd * sigma / c_6};

      van_der_waals_energy += epsilon / c_12 - sigma / c_6;
    }

    if (van_der_waals_energy - vdw_initial > 500.0) {
      // Intentionally removed I/O inside hot loop
    }
  }

  return 0;
}

int helper_print_energy(size_t nt, vector<vector<data_out>> &data_lg,
                        vector<vector<data_out>> &data_rec,
                        fstream &lg_ener_inter, fstream &rec_ener_inter,
                        bool flag_all_mode, const arguments &Arguments,
                        vector<energy_res> &Energy_res_mol1,
                        vector<energy_res> &Energy_res_mol2, size_t num_snap,
                        info &info_output) {
  print_energy_map(nt, data_lg, lg_ener_inter, flag_all_mode, Arguments,
                   Energy_res_mol2, num_snap, info_output);
  print_energy_map(nt, data_rec, rec_ener_inter, flag_all_mode, Arguments,
                   Energy_res_mol1, num_snap, info_output);
  return 0;
}

int wrapper_print_energy(size_t nt, vector<vector<data_out>> &data_lg,
                         vector<vector<data_out>> &data_rec, bool flag_all_mode,
                         const arguments &Arguments,
                         vector<energy_res> &Energy_res_mol1,
                         vector<energy_res> &Energy_res_mol2, size_t num_snap,
                         info &info_output) {
  fstream lg_ener_inter("ener_lg.txt", ios::out | ios::app);
  // lg_ener_inter << "Atom" << setw(12) << "Electro." << setw(12) << "Van der
  // W." << '\n';
  fstream rec_ener_inter("ener_rec.txt", ios::out | ios::app);
  // rec_ener_inter << "Res" << setw(12) << "Electro." << setw(12) << "Van der
  // W." << '\n';

  if (nt >= 2) {
    // --- P2: unordered_map replaces std::map for O(1) aggregation ---
    // pair_hash from miscellaneous.h supports pair<string, size_t>

    // Aggregate Ligand Data
    unordered_map<pair<string, size_t>, data_out, pair_hash> aggregated_lg;
    for (const auto &thread_vec : data_lg) {
      for (const auto &item : thread_vec) {
        auto [it, inserted] =
            aggregated_lg.try_emplace({item.chain, item.resnr}, item);
        if (!inserted) {
          it->second.E_E += item.E_E;
          it->second.E_V += item.E_V;
          it->second.F_elec += item.F_elec;
          it->second.F_vdw += item.F_vdw;
        }
      }
    }

    vector<vector<data_out>> aggregated_vec_lg(1);
    aggregated_vec_lg[0].reserve(aggregated_lg.size());
    for (auto const &[key, val] : aggregated_lg) {
      aggregated_vec_lg[0].emplace_back(std::move(val));
    }

    // Aggregate Receptor Data
    unordered_map<pair<string, size_t>, data_out, pair_hash> aggregated_rec;
    for (const auto &thread_vec : data_rec) {
      for (const auto &item : thread_vec) {
        auto [it, inserted] =
            aggregated_rec.try_emplace({item.chain, item.resnr}, item);
        if (!inserted) {
          it->second.E_E += item.E_E;
          it->second.E_V += item.E_V;
          it->second.F_elec += item.F_elec;
          it->second.F_vdw += item.F_vdw;
        }
      }
    }

    vector<vector<data_out>> aggregated_vec_rec(1);
    aggregated_vec_rec[0].reserve(aggregated_rec.size());
    for (auto const &[key, val] : aggregated_rec) {
      aggregated_vec_rec[0].emplace_back(std::move(val));
    }

    thread print(print_energy_map, 1, ref(aggregated_vec_lg),
                 ref(lg_ener_inter), flag_all_mode, ref(Arguments),
                 ref(Energy_res_mol2), num_snap, ref(info_output));
    thread print2(print_energy_map, 1, ref(aggregated_vec_rec),
                  ref(rec_ener_inter), flag_all_mode, ref(Arguments),
                  ref(Energy_res_mol1), num_snap, ref(info_output));
    print.join();
    print2.join();
  } else {
    helper_print_energy(nt, data_lg, data_rec, lg_ener_inter, rec_ener_inter,
                        flag_all_mode, Arguments, Energy_res_mol1,
                        Energy_res_mol2, num_snap, info_output);
  }
  return 0;
}

int cal_energy_amber_screening(const MOL &rec, const MOL &lg,
                               const vector<mol_info> &info_rec,
                               const vector<mol_info> &info_lg,
                               cache_sigma_epsilon &CachedData, size_t begin,
                               size_t end, cache_energy &Energy,
                               size_t num_snap, const arguments &Arguments,
                               const Multi_MOL_parameters &rec_parameters,
                               const Multi_MOL_parameters &lg_parameters) {
  double epsilon_in{Arguments.epsilon_in},
      lambda_debye{3.07 / sqrt(Arguments.salt_concentration)};
  Energy.set_energy(0., 0.);

  for (size_t j{begin}; j < end; ++j) {
    double temp_epsilon_lg{16. * lg_parameters[j].get_epsilon()},
        temp_charge{138.935457839 * lg_parameters[j].get_charge()},
        temp_sigma{lg_parameters[j].get_sigma()};
    size_t rec_size{rec.size()};
    for (size_t i{}; i < rec_size; ++i) {
      double d{0.1 * sqrt(square(rec.pos.x[i] - lg.pos.x[j]) +
                          square(rec.pos.y[i] - lg.pos.y[j]) +
                          square(rec.pos.z[i] - lg.pos.z[j]))},
          c_6{c6(0.5 * (temp_sigma + rec_parameters[i].get_sigma()) / d)};

      if (Arguments.distance_energy_cutoff > 0. &&
          d > Arguments.distance_energy_cutoff) [[unlikely]] {
        size_t pos{rec.size() * j + i};
        CachedData.set_elec_vdw(pos, 0., 0.);
        continue;
      }

      if (Arguments.self_int == true &&
          (d < 0.5 * (temp_sigma + rec_parameters[i].get_sigma()) ||
           (Arguments.flag_interchain_energy == true &&
            info_lg[j].get_chain() == info_rec[i].get_chain()))) {
        continue;
      }

      if (CachedData.elec_ener.empty()) {
        Energy.elec_ener += temp_charge * rec_parameters[i].get_charge() *
                            exp(-d / lambda_debye) / d / epsilon_in;
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
              pos,
              temp_charge * rec_parameters[i].get_charge() *
                  exp(-d / lambda_debye) / d / epsilon_in,
              sqrt(rec_parameters[i].get_epsilon() * temp_epsilon_lg) *
                  (square(c_6) - 0.000));
        else
          CachedData.set_elec_vdw(
              pos,
              temp_charge * rec_parameters[i].get_charge() *
                  exp(-d / lambda_debye) / d / epsilon_in,
              sqrt(rec_parameters[i].get_epsilon() * temp_epsilon_lg) *
                  (square(c_6) - c_6));
        Energy += CachedData.get_energy(pos);
      }
    }
  }
  return 0;
}

void cal_energy_amber(const MOL &rec, const MOL &lg,
                      const vector<mol_info> &info_rec,
                      const vector<mol_info> &info_lg,
                      cache_sigma_epsilon &CachedData, size_t begin, size_t end,
                      cache_energy &Energy, size_t num_snap,
                      const arguments &Arguments,
                      const Multi_MOL_parameters &rec_parameters,
                      const Multi_MOL_parameters &lg_parameters) {
  Energy.set_energy(0., 0.);
  double epsilon_in{Arguments.epsilon_in};
  for (size_t j{begin}; j < end; ++j) {
    const double pos_x{lg.pos.x[j]}, pos_y{lg.pos.y[j]}, pos_z{lg.pos.z[j]};
    const size_t rec_size{rec.size()};
    for (size_t i{}; i < rec_size; ++i) {
      double d{0.1 * sqrt(square(rec.pos.x[i] - pos_x) +
                          square(rec.pos.y[i] - pos_y) +
                          square(rec.pos.z[i] - pos_z))},
          c_6{c6(
              0.5 *
              (lg_parameters[j].get_sigma() + rec_parameters[i].get_sigma()) /
              d)};

      if (Arguments.distance_energy_cutoff > 0. &&
          d > Arguments.distance_energy_cutoff) [[unlikely]] {
        size_t pos{rec.size() * j + i};
        CachedData.set_elec_vdw(pos, 0., 0.);
        continue;
      }

      if (Arguments.self_int == true &&
          (d < 0.5 * (lg_parameters[j].get_sigma() +
                      rec_parameters[i].get_sigma()) ||
           (Arguments.flag_interchain_energy == true &&
            info_lg[j].get_chain() == info_rec[i].get_chain()))) {
        continue;
      }

      double temp_charge{138.935457839 * lg_parameters[j].get_charge()};
      double e_val = temp_charge * rec_parameters[i].get_charge() / d /
                     Arguments.epsilon_in;
      double temp_epsilon_lg{16. * lg_parameters[j].get_epsilon()};
      double v_val;
      if (lg_parameters[j].get_sigma() < 0. ||
          rec_parameters[i].get_sigma() < 0.)
        v_val = sqrt(temp_epsilon_lg * rec_parameters[i].get_epsilon()) *
                (square(c_6) - 0.000);
      else
        v_val = sqrt(temp_epsilon_lg * rec_parameters[i].get_epsilon()) *
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
void cal_energy_opls(const MOL &rec, const MOL &lg,
                     const vector<mol_info> &info_rec,
                     const vector<mol_info> &info_lg,
                     cache_sigma_epsilon &CachedData, size_t begin, size_t end,
                     cache_energy &Energy, size_t num_snap,
                     const arguments &Arguments,
                     const Multi_MOL_parameters &rec_parameters,
                     const Multi_MOL_parameters &lg_parameters) {
  Energy.set_energy(0., 0.);
  for (size_t j{begin}; j < end; ++j) {
    double temp_epsilon_lg{16. * lg_parameters[j].get_epsilon()},
        temp_charge{138.935457839 * lg_parameters[j].get_charge() /
                    Arguments.epsilon_in},
        temp_sigma{lg_parameters[j].get_sigma()}, pos_x{lg.pos.x[j]},
        pos_y{lg.pos.y[j]}, pos_z{lg.pos.z[j]};
    size_t rec_size{rec.size()};

    for (size_t i{}; i < rec_size; ++i) {
      double d{0.1 * sqrt(square(rec.pos.x[i] - pos_x) +
                          square(rec.pos.y[i] - pos_y) +
                          square(rec.pos.z[i] - pos_z))},
          c_6{c6(sqrt(temp_sigma * rec_parameters[i].get_sigma()) / d)};

      if (Arguments.self_int == true &&
          (c_6 > 1. || (Arguments.flag_interchain_energy == true &&
                        info_lg[j].get_chain() == info_rec[i].get_chain()))) {
        continue;
      }

      if (Arguments.distance_energy_cutoff > 0. &&
          d > Arguments.distance_energy_cutoff) [[unlikely]] {
        size_t pos{rec.size() * j + i};
        CachedData.set_elec_vdw(pos, 0., 0.);
        continue;
      }

      if (CachedData.elec_ener.empty()) {
        Energy.elec_ener += temp_charge * rec_parameters[i].get_charge() / d;
        Energy.vdw_ener +=
            sqrt(temp_epsilon_lg * rec_parameters[i].get_epsilon()) *
            (square(c_6) - c_6);
      } else {
        size_t pos{rec.size() * j + i};
        CachedData.set_elec_vdw(
            pos, temp_charge * rec_parameters[i].get_charge() / d,
            sqrt(rec_parameters[i].get_epsilon() * temp_epsilon_lg) *
                (square(c_6) - c_6));
        Energy += CachedData.get_energy(pos);
      }
    }
  }
}

#ifdef __AVX__

#include "avx_math.h"

#endif

void cal_energy_gromos(const MOL &rec, const MOL &lg,
                       const vector<mol_info> &info_rec,
                       const vector<mol_info> &info_lg,
                       cache_sigma_epsilon &CachedData, size_t begin,
                       size_t end, cache_energy &Energy, size_t num_snap,
                       const arguments &Arguments,
                       const Multi_MOL_parameters &rec_parameters,
                       const Multi_MOL_parameters &lg_parameters) {
  Energy.set_energy(0., 0.);
  for (size_t j{begin}; j < end; ++j) {
    double temp_charge{138.935457839 * lg_parameters[j].get_charge() /
                       Arguments.epsilon_in},
        temp_epsilon{lg_parameters[j].get_epsilon()},
        temp_sigma{lg_parameters[j].get_sigma()}, pos_x{lg.pos.x[j]},
        pos_y{lg.pos.y[j]}, pos_z{lg.pos.z[j]};
    size_t rec_size{rec.size()};
    for (size_t i{}; i < rec_size; ++i) {
      double d{0.1 * sqrt(square(rec.pos.x[i] - pos_x) +
                          square(rec.pos.y[i] - pos_y) +
                          square(rec.pos.z[i] - pos_z))},
          c_6{c6(d)};

      if (num_snap == 0 && i == 1719 && j == 2584) {
        ofstream err("global_debug.log", ios::app);
        err << "GlobalPair: i=" << i << " j=" << j << " D=" << d
            << " c6=" << c_6 << " Cutoff=" << Arguments.distance_energy_cutoff
            << " Self=" << Arguments.self_int << endl;
      }

      if (Arguments.distance_energy_cutoff > 0. &&
          d > Arguments.distance_energy_cutoff) [[unlikely]] {
        size_t pos{rec.size() * j + i};
        CachedData.set_elec_vdw(pos, 0., 0.);
        continue;
      }

      if (Arguments.self_int == true &&
          (i == j || (Arguments.flag_interchain_energy == true &&
                      info_lg[j].get_chain() == info_rec[i].get_chain()))) {
        continue;
      }

      if (CachedData.elec_ener.empty()) {
        Energy.elec_ener += temp_charge * rec_parameters[i].get_charge() / d;
        Energy.vdw_ener +=
            sqrt(rec_parameters[i].get_epsilon() * temp_epsilon) / square(c_6) -
            sqrt(rec_parameters[i].get_sigma() * temp_sigma) / c_6;
      } else {
        size_t pos{rec.size() * j + i};
        CachedData.set_elec_vdw(
            pos, temp_charge * rec_parameters[i].get_charge() / d,
            sqrt(rec_parameters[i].get_epsilon() * temp_epsilon) / square(c_6) -
                sqrt(rec_parameters[i].get_sigma() * temp_sigma) / c_6);
        Energy += CachedData.get_energy(pos);
      }
    }
  }
}

// MARTINI: Lorentz-Berthelot arithmetic mean for sigma, geometric mean for
// epsilon via sqrt(16*e1 * e2) = 4*sqrt(e1*e2). Sigma in nm, epsilon in
// kJ/mol. This function is intentionally identical to cal_energy_amber.
void cal_energy_martini(const MOL &rec, const MOL &lg,
                        const vector<mol_info> &info_rec,
                        const vector<mol_info> &info_lg,
                        cache_sigma_epsilon &CachedData, size_t begin,
                        size_t end, cache_energy &Energy, size_t num_snap,
                        const arguments &Arguments,
                        const Multi_MOL_parameters &rec_parameters,
                        const Multi_MOL_parameters &lg_parameters) {
  Energy.set_energy(0., 0.);
  double epsilon_in{Arguments.epsilon_in};
  for (size_t j{begin}; j < end; ++j) {
    const double pos_x{lg.pos.x[j]}, pos_y{lg.pos.y[j]}, pos_z{lg.pos.z[j]};
    const size_t rec_size{rec.size()};
    for (size_t i{}; i < rec_size; ++i) {
      double d{0.1 * sqrt(square(rec.pos.x[i] - pos_x) +
                          square(rec.pos.y[i] - pos_y) +
                          square(rec.pos.z[i] - pos_z))};
      if (d < 1e-10) {
        if (!CachedData.elec_ener.empty()) {
          size_t pos{rec.size() * j + i};
          CachedData.set_elec_vdw(pos, 0., 0.);
        }
        continue; // Guard: beads at same position, avoid 0/0 NaN
      }
      double c_6{c6(
          0.5 * (lg_parameters[j].get_sigma() + rec_parameters[i].get_sigma()) /
          d)};

      if (Arguments.distance_energy_cutoff > 0. &&
          d > Arguments.distance_energy_cutoff) [[unlikely]] {
        if (!CachedData.elec_ener.empty()) {
          size_t pos{rec.size() * j + i};
          CachedData.set_elec_vdw(pos, 0., 0.);
        }
        continue;
      }

      if (Arguments.self_int == true &&
          (d < 0.5 * (lg_parameters[j].get_sigma() +
                      rec_parameters[i].get_sigma()) ||
           (Arguments.flag_interchain_energy == true &&
            info_lg[j].get_chain() == info_rec[i].get_chain()))) {
        continue;
      }

      double temp_charge{138.935457839 * lg_parameters[j].get_charge()};
      double e_val = temp_charge * rec_parameters[i].get_charge() / d /
                     Arguments.epsilon_in;
      double temp_epsilon_lg{16. * lg_parameters[j].get_epsilon()};
      double v_val;
      if (lg_parameters[j].get_sigma() < 0. ||
          rec_parameters[i].get_sigma() < 0.)
        v_val = sqrt(temp_epsilon_lg * rec_parameters[i].get_epsilon()) *
                (square(c_6) - 0.000);
      else
        v_val = sqrt(temp_epsilon_lg * rec_parameters[i].get_epsilon()) *
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

int imp_cal_energy_lg_atom(vector<data_out> &lg_ener_inter,
                           const Multi_MOL &rec, const Multi_MOL &lg,
                           size_t num_snap, size_t begin, size_t end,
                           const vector<d_min> &dmin,
                           const cache_sigma_epsilon &CachedData,
                           const Multi_MOL_parameters &rec_parameters,
                           const Multi_MOL_parameters &lg_parameters,
                           const arguments &Arguments, bool is_swapped = false,
                           info_ener_res *infoener_res = nullptr) {

  // Hoist invariant bools outside all loops
  const bool need_emap =
      infoener_res &&
      (Arguments.energy_map == 1 || Arguments.energy_map == 4);
  const bool is_lig_mol1 =
      (Arguments.mol_1_type == "ligand" || Arguments.mol_1_type == "lig");
  const bool is_lig_mol2 =
      (Arguments.mol_2_type == "ligand" || Arguments.mol_2_type == "lig");

  // --- P1: Build compact lookup tables matching original key convention ---
  // j always indexes rec param, i always indexes lg param.
  // Original code: when !swapped, r_info=rec[j], l_info=lg[i]
  //                when  swapped, r_info=lg[i],  l_info=rec[j]
  // Tuple always: (r_chain, l_chain, r_resnr, l_resnr)
  // So j-side lookup from rec.info, i-side from lg.info — no swap undo needed.
  size_t n_j_res = 0, n_i_res = 0;
  vector<size_t> j_res_lookup, i_res_lookup;
  vector<size_t> j_chain_keys, j_resnr_keys;
  vector<size_t> i_chain_keys, i_resnr_keys;

  if (need_emap) {
    // is_lig flags follow the PARAMETER, not the argument label.
    // When !swapped: rec=mol1→is_lig_mol1, lg=mol2→is_lig_mol2
    // When  swapped: rec=mol2→is_lig_mol2, lg=mol1→is_lig_mol1
    const bool j_is_lig = is_swapped ? is_lig_mol2 : is_lig_mol1;
    const bool i_is_lig = is_swapped ? is_lig_mol1 : is_lig_mol2;

    // Build j-side lookup from rec.info (j iterates over rec)
    j_res_lookup.resize(rec.info.size());
    {
      unordered_map<pair<size_t, size_t>, size_t, pair_hash> seen;
      for (size_t a = 0; a < rec.info.size(); ++a) {
        size_t chin = rec.info[a].get_chainnum();
        size_t resnr =
            j_is_lig ? (a + 1) : rec.info[a].get_internal_resnr();
        auto [it, inserted] = seen.try_emplace({chin, resnr}, n_j_res);
        if (inserted) {
          j_chain_keys.push_back(chin);
          j_resnr_keys.push_back(resnr);
          ++n_j_res;
        }
        j_res_lookup[a] = it->second;
      }
    }

    // Build i-side lookup from lg.info (i iterates over lg)
    i_res_lookup.resize(lg.info.size());
    {
      unordered_map<pair<size_t, size_t>, size_t, pair_hash> seen;
      for (size_t a = 0; a < lg.info.size(); ++a) {
        size_t chin = lg.info[a].get_chainnum();
        size_t resnr =
            i_is_lig ? (a + 1) : lg.info[a].get_internal_resnr();
        auto [it, inserted] = seen.try_emplace({chin, resnr}, n_i_res);
        if (inserted) {
          i_chain_keys.push_back(chin);
          i_resnr_keys.push_back(resnr);
          ++n_i_res;
        }
        i_res_lookup[a] = it->second;
      }
    }
  }

  // Per-thread local flat accumulation arrays (no sharing between threads)
  const size_t flat_size = n_j_res * n_i_res;
  vector<double> local_ee, local_ev;
  if (need_emap && flat_size > 0) {
    local_ee.assign(flat_size, 0.0);
    local_ev.assign(flat_size, 0.0);
  }

  // --- Main atom-pair loop ---
  for (size_t i = begin; i < end; ++i) {
    double E_E{}, E_V{};
    for (size_t j{}; j < rec.m_mol[num_snap].size(); ++j) {
      size_t stride =
          is_swapped ? lg.m_mol[num_snap].size() : rec.m_mol[num_snap].size();
      size_t cache_idx = is_swapped ? stride * j + i : stride * i + j;

      double e_temp{}, v_temp{};
      add_energy(CachedData, cache_idx, e_temp, v_temp, rec_parameters,
                 lg_parameters, rec.m_mol[num_snap], lg.m_mol[num_snap],
                 num_snap, Arguments, rec.info, lg.info, j, i);

      E_E += e_temp;
      E_V += v_temp;

      // --- P1: O(1) flat array write (j-side × i-side) ---
      if (need_emap) {
        if (i < i_res_lookup.size() && j < j_res_lookup.size()) {
          size_t flat_idx = j_res_lookup[j] * n_i_res + i_res_lookup[i];
          local_ee[flat_idx] += e_temp;
          local_ev[flat_idx] += v_temp;
        }
      }
    }
    bool is_current_lig = is_swapped ? is_lig_mol1 : is_lig_mol2;
    size_t out_resnr = is_current_lig ? (i + 1) : lg.info[i].get_resnr();

    lg_ener_inter.emplace_back(out_resnr, lg.info[i].get_chain(),
                               lg.info[i].get_original_res_name(), E_E, E_V);
  }

  // --- P1: Single-pass flush to infoener_res ---
  // Original convention: tuple = (r_chain, l_chain, r_resnr, l_resnr)
  //   !swapped: r=rec[j]=j_side, l=lg[i]=i_side → (j, i) order
  //    swapped: r=lg[i]=i_side, l=rec[j]=j_side → (i, j) order
  if (need_emap && flat_size > 0) {
    for (size_t rj = 0; rj < n_j_res; ++rj) {
      for (size_t ri = 0; ri < n_i_res; ++ri) {
        size_t flat_idx = rj * n_i_res + ri;
        if (local_ee[flat_idx] != 0.0 || local_ev[flat_idx] != 0.0) {
          if (!is_swapped) {
            infoener_res->add_energy(
                make_tuple(j_chain_keys[rj], i_chain_keys[ri],
                           j_resnr_keys[rj], i_resnr_keys[ri]),
                local_ee[flat_idx], local_ev[flat_idx]);
          } else {
            infoener_res->add_energy(
                make_tuple(i_chain_keys[ri], j_chain_keys[rj],
                           i_resnr_keys[ri], j_resnr_keys[rj]),
                local_ee[flat_idx], local_ev[flat_idx]);
          }
        }
      }
    }
  }

  return 0;
}



int imp_cal_energy_lg(vector<data_out> &lg_ener_inter, const Multi_MOL &rec,
                      const Multi_MOL &lg, size_t num_snap, size_t begin,
                      size_t end, const vector<d_min> &dmin,
                      const cache_sigma_epsilon &CachedData, bool flag_all_mode,
                      const Multi_MOL_parameters &rec_parameters,
                      const Multi_MOL_parameters &lg_parameters,
                      unordered_map<pair<size_t, size_t>, pair<size_t, size_t>,
                                    pair_hash> &map_list_res,
                      vector<pair<size_t, size_t>> &list_res,
                      info_ener_res &infoener_res, const arguments &Arguments) {
  int last_rec_chain = -1, last_lg_chain = -1, last_rec_res = -1,
      last_lg_res = -1;
  double sum_ee_emap{}, sum_ev_emap{};
  auto add_emap_local = [&](int r_c, int l_c, int r_r, int l_r, double ee,
                            double ev) {
    if (ee == 0.0 && ev == 0.0)
      return;
    if (last_rec_chain != r_c || last_lg_chain != l_c || last_rec_res != r_r ||
        last_lg_res != l_r) {
      if (last_rec_chain != -1) {
        infoener_res.add_energy(make_tuple(last_rec_chain, last_lg_chain,
                                           last_rec_res, last_lg_res),
                                sum_ee_emap, sum_ev_emap);
      }
      last_rec_chain = r_c;
      last_lg_chain = l_c;
      last_rec_res = r_r;
      last_lg_res = l_r;
      sum_ee_emap = 0;
      sum_ev_emap = 0;
    }
    sum_ee_emap += ee;
    sum_ev_emap += ev;
  };

  for (size_t ite{begin}; ite < end; ++ite) {
    double E_E{}, E_V{};
    pair<size_t, size_t> ptemp{map_list_res[list_res[ite]]};
    size_t j_start = dmin[ptemp.first].base[ptemp.second].atom_num - 1;
    while (j_start > 0 &&
           lg.info[j_start - 1].get_internal_resnr() ==
               lg.info[j_start].get_internal_resnr() &&
           lg.info[j_start - 1].get_chain() == lg.info[j_start].get_chain())
      j_start--;

    size_t num = j_start + 300 > lg.m_mol[num_snap].size()
                     ? lg.m_mol[num_snap].size()
                     : j_start + 300;

    for (size_t j = j_start; j < num; ++j) {
      if (lg.info[j].get_internal_resnr() ==
              lg.info[dmin[ptemp.first].base[ptemp.second].atom_num - 1]
                  .get_internal_resnr() &&
          lg.info[j].get_chain() ==
              lg.info[dmin[ptemp.first].base[ptemp.second].atom_num - 1]
                  .get_chain()) {
        if (flag_all_mode == true) {
          for (size_t k{}; k < rec.m_mol[num_snap].size(); ++k) {
            double ee{}, ev{};
            add_energy(CachedData, rec.m_mol[0].size() * j + k, ee, ev,
                       rec_parameters, lg_parameters, rec.m_mol[num_snap],
                       lg.m_mol[num_snap], num_snap, Arguments, rec.info,
                       lg.info, k, j);
            E_E += ee;
            E_V += ev;

            if (Arguments.energy_map == 1 || Arguments.energy_map == 4) {
              add_emap_local(rec.info[k].get_chainnum(),
                             lg.info[j].get_chainnum(),
                             rec.info[k].get_internal_resnr(),
                             lg.info[j].get_internal_resnr(), ee, ev);
            }
          }
        } else {
          size_t beginpos{
              dmin[ptemp.first].base[ptemp.second].atom_num_rec - 1 >= 0
                  ? dmin[ptemp.first].base[ptemp.second].atom_num_rec - 1
                  : 0};
          size_t k_start = beginpos;
          while (k_start > 0 &&
                 rec.info[k_start - 1].get_internal_resnr() ==
                     rec.info[k_start].get_internal_resnr() &&
                 rec.info[k_start - 1].get_chain() ==
                     rec.info[k_start].get_chain())
            k_start--;

          for (size_t k{k_start}; k < rec.m_mol[num_snap].size(); ++k) {
            if (rec.info[k].get_internal_resnr() ==
                    rec.info[dmin[ptemp.first].base[ptemp.second].atom_num_rec -
                             1]
                        .get_internal_resnr() &&
                rec.info[k].get_chain() ==
                    rec.info[dmin[ptemp.first].base[ptemp.second].atom_num_rec -
                             1]
                        .get_chain()) {
              double ee{}, ev{};
              add_energy(CachedData, rec.m_mol[0].size() * j + k, ee, ev,
                         rec_parameters, lg_parameters, rec.m_mol[num_snap],
                         lg.m_mol[num_snap], num_snap, Arguments, rec.info,
                         lg.info, k, j);
              E_E += ee;
              E_V += ev;
              if (Arguments.energy_map == 1 || Arguments.energy_map == 4) {
                bool is_lig_rec = (Arguments.mol_1_type == "ligand" ||
                                   Arguments.mol_1_type == "lig");
                bool is_lig_lg = (Arguments.mol_2_type == "ligand" ||
                                  Arguments.mol_2_type == "lig");
                size_t rec_in =
                    is_lig_rec ? (k + 1) : rec.info[k].get_internal_resnr();
                size_t lg_in =
                    is_lig_lg ? (j + 1) : lg.info[j].get_internal_resnr();
                add_emap_local(rec.info[k].get_chainnum(),
                               lg.info[j].get_chainnum(), rec_in, lg_in, ee,
                               ev);
              }
            } else
              break;
          }
        }
      } else
        break;
    }

    lg_ener_inter.emplace_back(dmin[ptemp.first].base[ptemp.second].resnr_ab,
                               dmin[ptemp.first].base[ptemp.second].chain_ab,
                               dmin[ptemp.first].base[ptemp.second].res_name_ab,
                               E_E, E_V);
  }

  if (last_rec_chain != -1 && (sum_ee_emap != 0.0 || sum_ev_emap != 0.0) &&
      (Arguments.energy_map == 1 || Arguments.energy_map == 4)) {
    infoener_res.add_energy(
        make_tuple(last_rec_chain, last_lg_chain, last_rec_res, last_lg_res),
        sum_ee_emap, sum_ev_emap);
  }

  /*unordered_set<pair<size_t, size_t>, pair_hash>
  energy_list;

  for (size_t ii{ begin }; ii < end; ++ii)
  {
          for (size_t i{}; i < dmin[ii].base.size(); ++i)
          {
                  if
  (energy_list.emplace(dmin[ii].base[i].chainid_ab,
  dmin[ii].base[i].resnr_ab).second == false) continue;

                  double E_E{}, E_V{};
                  size_t num = dmin[ii].base[i].atom_num -
  1 + 300 > lg.m_mol[num_snap].size() ?
  lg.m_mol[num_snap].size() :
  dmin[ii].base[i].atom_num - 1 + 300;

                  for (size_t j = dmin[ii].base[i].atom_num
  - 1; j < num; ++j)
                  {
                          if (lg.info[j].get_resnr() ==
  lg.info[dmin[ii].base[i].atom_num - 1].get_resnr()
                                  && lg.info[j].get_chain()
  == lg.info[dmin[ii].base[i].atom_num - 1].get_chain())
                          {

                                  if (flag_all_mode ==
  true)
                                  {
                                          for (size_t k{};
  k < rec.m_mol[num_snap].size(); ++k)
  add_energy(CachedData, rec.m_mol[0].size() * j +
  k, E_E, E_V, rec_parameters, lg_parameters,
  rec.m_mol[num_snap], lg.m_mol[num_snap], num_snap,
  Arguments, rec.info, lg.info, k, j);
                                  }
                                  else
                                  {
                                          size_t beginpos{
  dmin[ii].base[i].atom_num_rec - 1 >= 0 ?
  dmin[ii].base[i].atom_num_rec - 1 : 0
  }; for (size_t k{ beginpos }; k <
  rec.m_mol[num_snap].size(); ++k)
                                          {
                                                  if
  (rec.info[k].get_resnr() ==
  rec.info[dmin[ii].base[i].atom_num_rec - 1].get_resnr()
                                                          &&
  rec.info[k].get_chain() ==
  rec.info[dmin[ii].base[i].atom_num_rec - 1].get_chain())
                                                  {
                                                          add_energy(CachedData,
  rec.m_mol[0].size() * j + k, E_E, E_V,
  rec_parameters, lg_parameters, rec.m_mol[num_snap],
  lg.m_mol[num_snap], num_snap, Arguments, rec.info,
  lg.info, k, j);
                                                  }
                                                  else
  break;
                                          }
                                  }
                          }
                          else
                                  break;
                  }

                  lg_ener_inter.emplace_back(dmin[ii].base[i].resnr_ab,
  dmin[ii].base[i].chain_ab, dmin[ii].base[i].res_name_ab,
  E_E, E_V);
          }
  }
  */
  return 0;
}

int imp_cal_energy_rec(
    vector<data_out> &rec_ener_inter, info_ener_res &infoener_res,
    const Multi_MOL &rec, const Multi_MOL &lg, size_t num_snap, size_t begin,
    size_t end, const vector<d_min> &dmin,
    const cache_sigma_epsilon &CachedData, bool flag_all_mode,
    const Multi_MOL_parameters &rec_parameters,
    const Multi_MOL_parameters &lg_parameters,
    unordered_map<pair<size_t, size_t>, pair<size_t, size_t>, pair_hash>
        &map_list_res,
    vector<pair<size_t, size_t>> &list_res, const arguments &Arguments,
    const vector<bool> &is_lg_dmin_res) {
  int last_rec_chain = -1, last_lg_chain = -1, last_rec_res = -1,
      last_lg_res = -1;
  double sum_ee_emap{}, sum_ev_emap{};
  auto add_emap_local = [&](int r_c, int l_c, int r_r, int l_r, double ee,
                            double ev) {
    if (ee == 0.0 && ev == 0.0)
      return;
    if (last_rec_chain != r_c || last_lg_chain != l_c || last_rec_res != r_r ||
        last_lg_res != l_r) {
      if (last_rec_chain != -1) {
        infoener_res.add_energy(make_tuple(last_rec_chain, last_lg_chain,
                                           last_rec_res, last_lg_res),
                                sum_ee_emap, sum_ev_emap);
      }
      last_rec_chain = r_c;
      last_lg_chain = l_c;
      last_rec_res = r_r;
      last_lg_res = l_r;
      sum_ee_emap = 0;
      sum_ev_emap = 0;
    }
    sum_ee_emap += ee;
    sum_ev_emap += ev;
  };

  for (size_t ite{begin}; ite < end; ++ite) {
    double E_E{}, E_V{};
    pair<size_t, size_t> ptemp{map_list_res[list_res[ite]]};
    if (flag_all_mode == true) {
      size_t j_limit = rec.m_mol[num_snap].size();
      size_t lg_limit = lg.m_mol[num_snap].size();
      auto rec_info_dmin{
          rec.info[dmin[ptemp.first].base[ptemp.second].atom_num_rec - 1]};

      size_t j_start = dmin[ptemp.first].base[ptemp.second].atom_num_rec - 1;
      while (j_start > 0 &&
             rec.info[j_start - 1].get_internal_resnr() ==
                 rec.info[j_start].get_internal_resnr() &&
             rec.info[j_start - 1].get_chain() == rec.info[j_start].get_chain())
        j_start--;

      size_t num = j_start + 300 > j_limit ? j_limit : j_start + 300;

      for (size_t j{j_start}; j < num; ++j) {
        if (rec.info[j].get_internal_resnr() ==
                rec_info_dmin.get_internal_resnr() &&
            rec.info[j].get_chainnum() == rec_info_dmin.get_chainnum()) {
          for (size_t k{}; k < lg_limit; ++k) {
            double ee{}, ev{};
            add_energy(CachedData, rec.m_mol[0].size() * k + j, ee, ev,
                       rec_parameters, lg_parameters, rec.m_mol[num_snap],
                       lg.m_mol[num_snap], num_snap, Arguments, rec.info,
                       lg.info, j, k);

            E_E += ee;
            E_V += ev;
            if (Arguments.energy_map == 1 || Arguments.energy_map == 4) {
              if (flag_all_mode == false || k >= lg.m_mol[num_snap].size() ||
                  is_lg_dmin_res.empty() ||
                  !is_lg_dmin_res[lg.info[k].get_internal_resnr()]) {
                add_emap_local(rec.info[j].get_chainnum(),
                               lg.info[k].get_chainnum(),
                               rec.info[j].get_internal_resnr(),
                               lg.info[k].get_internal_resnr(), ee, ev);
              }
            }
          }
        } else
          break;
      }
    } else {
      size_t j_start = dmin[ptemp.first].base[ptemp.second].atom_num_rec - 1;
      while (j_start > 0 &&
             rec.info[j_start - 1].get_internal_resnr() ==
                 rec.info[j_start].get_internal_resnr() &&
             rec.info[j_start - 1].get_chain() == rec.info[j_start].get_chain())
        j_start--;

      size_t num =
          j_start + 300 > rec.info.size() ? rec.info.size() : j_start + 300;
      for (size_t j = j_start; j < num; ++j) {
        auto &rec_info_dmin =
            rec.info[dmin[ptemp.first].base[ptemp.second].atom_num_rec - 1];
        auto &rec_info = rec.info[j];

        if (rec_info.get_internal_resnr() ==
                rec_info_dmin.get_internal_resnr() &&
            rec_info.get_chain() == rec_info_dmin.get_chain()) {
          size_t beginpos{dmin[ptemp.first].base[ptemp.second].atom_num - 1 >= 0
                              ? dmin[ptemp.first].base[ptemp.second].atom_num -
                                    1
                              : 0};
          size_t k_start = beginpos;
          while (k_start > 0 &&
                 lg.info[k_start - 1].get_internal_resnr() ==
                     lg.info[k_start].get_internal_resnr() &&
                 lg.info[k_start - 1].get_chain() ==
                     lg.info[k_start].get_chain())
            k_start--;

          for (size_t k{k_start}; k < lg.m_mol[num_snap].size(); ++k)
            if (lg.info[k].get_internal_resnr() ==
                    lg.info[dmin[ptemp.first].base[ptemp.second].atom_num - 1]
                        .get_internal_resnr() &&
                lg.info[k].get_chain() ==
                    lg.info[dmin[ptemp.first].base[ptemp.second].atom_num - 1]
                        .get_chain()) {
              double ee{}, ev{};
              add_energy(CachedData, rec.m_mol[0].size() * k + j, ee, ev,
                         rec_parameters, lg_parameters, rec.m_mol[num_snap],
                         lg.m_mol[num_snap], num_snap, Arguments, rec.info,
                         lg.info, j, k);
              E_E += ee;
              E_V += ev;

              if (Arguments.energy_map == 1 || Arguments.energy_map == 4) {
                if (is_lg_dmin_res.empty() ||
                    !is_lg_dmin_res[lg.info[k].get_internal_resnr()]) {
                  bool is_lig_rec = (Arguments.mol_1_type == "ligand" ||
                                     Arguments.mol_1_type == "lig");
                  bool is_lig_lg = (Arguments.mol_2_type == "ligand" ||
                                    Arguments.mol_2_type == "lig");
                  size_t rec_in =
                      is_lig_rec ? (j + 1) : rec_info.get_internal_resnr();
                  size_t lg_in =
                      is_lig_lg ? (k + 1) : lg.info[k].get_internal_resnr();
                  add_emap_local(rec_info.get_chainnum(),
                                 lg.info[k].get_chainnum(), rec_in, lg_in, ee,
                                 ev);
                }
              }
            } else
              break;
        }
      }
    }
    rec_ener_inter.emplace_back(
        dmin[ptemp.first].base[ptemp.second].resnr_rec,
        dmin[ptemp.first].base[ptemp.second].chain_rec,
        dmin[ptemp.first].base[ptemp.second].res_name_rec, E_E, E_V);
  }

  if (last_rec_chain != -1 && (sum_ee_emap != 0.0 || sum_ev_emap != 0.0) &&
      (Arguments.energy_map == 1 || Arguments.energy_map == 4)) {
    infoener_res.add_energy(
        make_tuple(last_rec_chain, last_lg_chain, last_rec_res, last_lg_res),
        sum_ee_emap, sum_ev_emap);
  }

  /*unordered_set<pair<size_t, size_t>, pair_hash>
  energy_list;

  for (size_t ii{begin}; ii < end; ++ii)
  {
          for (size_t i{}; i < dmin[ii].base.size(); ++i)
          {
                  if
  (energy_list.emplace(dmin[ii].base[i].chainid_rec,
  dmin[ii].base[i].resnr_rec).second == false) continue;

                  double E_E{}, E_V{};

                  if (flag_all_mode == true)
                  {
                          for (size_t k{}; k <
  lg.m_mol[num_snap].size(); ++k)
                          {
                                  size_t num =
  dmin[ii].base[i].atom_num_rec - 1
  + 300 > rec.info.size() ? rec.info.size() :
  dmin[ii].base[i].atom_num_rec - 1
  + 300; auto& rec_info_dmin =
  rec.info[dmin[ii].base[i].atom_num_rec - 1]; for (size_t
  j = dmin[ii].base[i].atom_num_rec - 1; j < num; ++j)
                                  {
                                          auto& rec_info =
  rec.info[j]; if (rec_info.get_resnr() ==
  rec_info_dmin.get_resnr()
                                                  &&
  rec_info.get_chain() == rec_info_dmin.get_chain())
                                          {
                                                  double
  ee{}, ev{}; add_energy(CachedData,
  rec.m_mol[0].size() * k + j, ee, ev,
  rec_parameters, lg_parameters, rec.m_mol[num_snap],
  lg.m_mol[num_snap], num_snap, Arguments, rec.info,
  lg.info, j, k); E_E += ee; E_V += ev;

                                                  if
  (Arguments.energy_map == 4)
                                                  {
                                                          infoener_res.add_energy(tuple(rec_info.get_chainnum(),
  lg.info[k].get_chainnum(), rec_info.get_resnr(),
  lg.info[k].get_resnr()), ee, ev);
                                                  }
                                          }
                                          else break;
                                  }
                          }
                  }
                  else
                  {
                          size_t num =
  dmin[ii].base[i].atom_num_rec - 1 + 300 > rec.info.size()
  ? rec.info.size() : dmin[ii].base[i].atom_num_rec - 1 +
  300; for (size_t j = dmin[ii].base[i].atom_num_rec - 1; j
  < num; ++j)
                          {
                                  auto& rec_info_dmin =
  rec.info[dmin[ii].base[i].atom_num_rec - 1]; auto&
  rec_info = rec.info[j];

                                  if (rec_info.get_resnr()
  == rec_info_dmin.get_resnr()
                                          &&
  rec_info.get_chain() == rec_info_dmin.get_chain())
                                  {
                                          size_t beginpos{
  dmin[ii].base[i].atom_num - 1 >= 0 ?
  dmin[ii].base[i].atom_num - 1 : 0 }; for (size_t k{
  beginpos }; k < lg.m_mol[num_snap].size(); ++k)
  if (lg.info[k].get_resnr() ==
  lg.info[dmin[ii].base[i].atom_num - 1].get_resnr()
                                                          &&
  lg.info[k].get_chain() ==
  lg.info[dmin[ii].base[i].atom_num - 1].get_chain())
                                                  {
                                                          double ee{}, ev{};
                                                          add_energy(CachedData,
  rec.m_mol[0].size() * k + j, ee, ev,
  rec_parameters, lg_parameters, rec.m_mol[num_snap],
  lg.m_mol[num_snap], num_snap, Arguments, rec.info,
  lg.info, j, k); E_E += ee; E_V += ev;

                                                          if
  (Arguments.energy_map == 4)
                                                          {
                                                                  infoener_res.add_energy(tuple(rec_info.get_chainnum(),
  lg.info[k].get_chainnum(), rec_info.get_resnr(),
  lg.info[k].get_resnr()), ee, ev);
                                                          }
                                                  }
                                                  else
                                                          break;
                                  }
                          }
                  }

                  rec_ener_inter.emplace_back(dmin[ii].base[i].resnr_rec,
  dmin[ii].base[i].chain_rec,
  dmin[ii].base[i].res_name_rec, E_E, E_V);
          }
  }
  //cout << "infoener res size: " <<
  infoener_res.emap_ee.size() <<
  '\n';*/
  return 0;
}

int imp_cal_energy_rec_res(multi_snap_emap &rec_ener_inter,
                           const Multi_MOL &rec, const Multi_MOL &lg,
                           size_t num_snap, size_t begin, size_t end,
                           const vector<d_min> &dmin,
                           const cache_sigma_epsilon &CachedData,
                           const arguments &Arguments) {
  unordered_set<tuple<size_t, size_t, size_t, size_t>, tuple_hash>
      dmin_noduplicate;
  for (size_t ii{0}; ii < dmin.size(); ++ii) {
    multi_emap_per_res emap_list_temp{num_snap};
    for (size_t i{}; i < dmin[ii].base.size(); ++i) {
      if (dmin_noduplicate
              .emplace(dmin[ii].base[i].chainid_ab,
                       dmin[ii].base[i].chainid_rec, dmin[ii].base[i].resnr_ab,
                       dmin[ii].base[i].resnr_rec)
              .second == false)
        continue;
      double E_E{}, E_V{};
      size_t num = dmin[ii].base[i].atom_num_rec - 1 + 250 > rec.info.size()
                       ? rec.info.size()
                       : dmin[ii].base[i].atom_num_rec - 1 + 250;
      // if (Arguments.flag_topo_2 == true)
      //	num = dmin[ii].base[i].atom_num_rec;

      for (size_t j = dmin[ii].base[i].atom_num_rec - 1; j < num; ++j) {
        size_t num2{dmin[ii].base[i].atom_num - 1 + 250 > lg.info.size()
                        ? lg.info.size()
                        : dmin[ii].base[i].atom_num - 1 + 250};
        // if (Arguments.flag_topo_1 == true)
        // num2 = dmin[ii].base[i].atom_num;

        for (size_t k = dmin[ii].base[i].atom_num - 1; k < num2; ++k)
          if (lg.info[k].get_resnr() ==
                  lg.info[dmin[ii].base[i].atom_num - 1].get_resnr() &&
              lg.info[k].get_chain() ==
                  lg.info[dmin[ii].base[i].atom_num - 1].get_chain() &&
              rec.info[j].get_resnr() ==
                  rec.info[dmin[ii].base[i].atom_num_rec - 1].get_resnr() &&
              rec.info[j].get_chain() ==
                  rec.info[dmin[ii].base[i].atom_num_rec - 1].get_chain()) {
            add_energy(CachedData, rec.m_mol[0].size() * k + j, E_E, E_V,
                       rec.parameters, lg.parameters, rec.m_mol[num_snap],
                       lg.m_mol[num_snap], num_snap, Arguments, rec.info,
                       lg.info, j, k);
          }
      }
      size_t resnr_lg{Arguments.flag_topo_1 == true
                          ? dmin[ii].base[i].atom_num
                          : lg.info[dmin[ii].base[i].atom_num - 1].get_resnr()};
      size_t resnr_rec{
          Arguments.flag_topo_2 == true
              ? dmin[ii].base[i].atom_num_rec
              : rec.info[dmin[ii].base[i].atom_num_rec - 1].get_resnr()};

      emap_list_temp.ResList.emplace_back(
          resnr_rec,
          rec.info[dmin[ii].base[i].atom_num_rec - 1].get_original_res_name(),
          rec.info[dmin[ii].base[i].atom_num_rec - 1].get_chain(),
          data_out(
              resnr_lg, lg.info[dmin[ii].base[i].atom_num - 1].get_chain(),
              lg.info[dmin[ii].base[i].atom_num - 1].get_original_res_name(),
              E_E, E_V));
    }
    rec_ener_inter.emap_pdb.emplace_back(emap_list_temp);
  }
  return 0;
}

void energyMapOf2Pro(vector<thread> &threads, info_ener_res &infoener_res,
                     const Multi_MOL &rec, const Multi_MOL &lg, size_t num_snap,
                     const vector<d_min> &dmin,
                     const cache_sigma_epsilon &CachedData, size_t nt,
                     bool flag_all_mode, const arguments &Arguments,
                     vector<energy_res> &Energy_res_mol1,
                     vector<energy_res> &Energy_res_mol2, info &info_output,
                     const Multi_MOL_parameters &rec_parameters,
                     const Multi_MOL_parameters &lg_parameters) {
  auto enermaptimestart = std::chrono::steady_clock::now();
  vector<vector<data_out>> data_lg(nt), data_rec(nt);
  vector<info_ener_res> infoener_res_v(nt);

  unordered_map<pair<size_t, size_t>, pair<size_t, size_t>, pair_hash>
      map_list_res, map_list_res_mol1;
  vector<pair<size_t, size_t>> list_res_mol2, list_res_mol1;
  for (size_t ii{}; ii < dmin.size(); ++ii) {
    for (size_t i{}; i < dmin[ii].base.size(); ++i) {
      // cout << dmin[ii].base[i].chainid_ab << " " <<
      // dmin[ii].base[i].resnr_ab
      // << '\n';
      if (true ==
          map_list_res
              .try_emplace(pair<size_t, size_t>(dmin[ii].base[i].chainid_ab,
                                                dmin[ii].base[i].resnr_ab),
                           pair<size_t, size_t>(ii, i))
              .second)
        list_res_mol2.emplace_back(dmin[ii].base[i].chainid_ab,
                                   dmin[ii].base[i].resnr_ab);
      if (true ==
          map_list_res_mol1
              .try_emplace(pair<size_t, size_t>(dmin[ii].base[i].chainid_rec,
                                                dmin[ii].base[i].resnr_rec),
                           pair<size_t, size_t>(ii, i))
              .second)
        list_res_mol1.emplace_back(dmin[ii].base[i].chainid_rec,
                                   dmin[ii].base[i].resnr_rec);
    }
  }
  // cout << "size_of map mol2: " << list_res_mol2.size()
  // << '\n'; cout << "size of map mol1: " <<
  // list_res_mol1.size() << '\n';

  for (size_t i{}; i < nt; ++i)
    threads[i] =
        thread(imp_cal_energy_lg, ref(data_lg[i]), ref(rec), ref(lg), num_snap,
               static_cast<size_t>(i * list_res_mol2.size() / nt),
               static_cast<size_t>((i + 1) * list_res_mol2.size() / nt),
               ref(dmin), ref(CachedData), flag_all_mode, ref(rec_parameters),
               ref(lg_parameters), ref(map_list_res), ref(list_res_mol2),
               ref(infoener_res_v[i]), ref(Arguments));
  for (auto &t : threads)
    t.join();

  vector<bool> is_lg_dmin_res;
  if (flag_all_mode) {
    size_t lg_max_resnr = 0;
    for (const auto &a : lg.info) {
      if (a.get_internal_resnr() > lg_max_resnr)
        lg_max_resnr = a.get_internal_resnr();
    }
    is_lg_dmin_res.assign(lg_max_resnr + 1, false);
    for (size_t ii{}; ii < dmin.size(); ++ii) {
      for (size_t i{}; i < dmin[ii].base.size(); ++i) {
        size_t atom_idx = dmin[ii].base[i].atom_num - 1;
        if (atom_idx < lg.info.size()) {
          is_lg_dmin_res[lg.info[atom_idx].get_internal_resnr()] = true;
        }
      }
    }
  }

  cache_sigma_epsilon empty_cache;
  for (size_t i{}; i < nt; ++i)
    threads[i] = thread(
        imp_cal_energy_rec, ref(data_rec[i]), ref(infoener_res_v[i]), ref(rec),
        ref(lg), num_snap, static_cast<size_t>(i * list_res_mol1.size() / nt),
        static_cast<size_t>((i + 1) * list_res_mol1.size() / nt), ref(dmin),
        ref(CachedData), flag_all_mode, ref(rec_parameters), ref(lg_parameters),
        ref(map_list_res_mol1), ref(list_res_mol1), ref(Arguments),
        std::cref(is_lg_dmin_res));
  for (auto &t : threads)
    t.join();
  // cout << "time of cal map snap " << num_snap << " (ms):
  // "
  //	<<
  // std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now()
  //- enermaptimestart).count() << " ";
  if (Arguments.energy_map == 1 || Arguments.energy_map == 4)
    for (auto &ii : infoener_res_v) {
      for (auto &ele : ii.position_map) {
        infoener_res.add_energy(ele.first,
                                ii.emap_ee[ele.second],
                                ii.emap_ev[ele.second]);
      }
    }

  wrapper_print_energy(nt, data_lg, data_rec, flag_all_mode, Arguments,
                       Energy_res_mol1, Energy_res_mol2, num_snap, info_output);
}

void energyMapOfProLig(vector<thread> &threads, info_ener_res &infoener_res,
                       const Multi_MOL &rec, const Multi_MOL &lg,
                       size_t num_snap, const vector<d_min> &dmin,
                       const cache_sigma_epsilon &CachedData, size_t nt,
                       bool flag_all_mode, const arguments &Arguments,
                       vector<energy_res> &Energy_res_mol1,
                       vector<energy_res> &Energy_res_mol2, info &info_output,
                       const Multi_MOL_parameters &rec_parameters,
                       const Multi_MOL_parameters &lg_parameters) {
  vector<vector<data_out>> data_lg(nt), data_rec(nt);
  vector<info_ener_res> infoener_res_v(nt);

  unordered_map<pair<size_t, size_t>, pair<size_t, size_t>, pair_hash>
      map_list_res, map_list_res_mol1;
  vector<pair<size_t, size_t>> list_res_mol2, list_res_mol1;
  for (size_t ii{}; ii < dmin.size(); ++ii) {
    for (size_t i{}; i < dmin[ii].base.size(); ++i) {
      if (true ==
          map_list_res
              .try_emplace(pair<size_t, size_t>(dmin[ii].base[i].chainid_ab,
                                                dmin[ii].base[i].resnr_ab),
                           pair<size_t, size_t>(ii, i))
              .second)
        list_res_mol2.emplace_back(dmin[ii].base[i].chainid_ab,
                                   dmin[ii].base[i].resnr_ab);
      if (true ==
          map_list_res_mol1
              .try_emplace(pair<size_t, size_t>(dmin[ii].base[i].chainid_rec,
                                                dmin[ii].base[i].resnr_rec),
                           pair<size_t, size_t>(ii, i))
              .second)
        list_res_mol1.emplace_back(dmin[ii].base[i].chainid_rec,
                                   dmin[ii].base[i].resnr_rec);
    }
  }

  for (size_t i{}; i < nt; ++i)
    threads[i] = thread(
        imp_cal_energy_lg_atom, ref(data_lg[i]), ref(rec), ref(lg), num_snap,
        static_cast<size_t>(i * lg.m_mol[num_snap].size() / nt),
        std::min(static_cast<size_t>((i + 1) * lg.m_mol[num_snap].size() / nt),
                 lg.info.size()),
        ref(dmin), ref(CachedData), ref(rec_parameters), ref(lg_parameters),
        ref(Arguments), false, &infoener_res_v[i]);
  for (auto &t : threads)
    t.join();

  for (size_t i{}; i < nt; ++i)
    threads[i] = thread(
        imp_cal_energy_lg_atom, ref(data_rec[i]), ref(lg), ref(rec), num_snap,
        static_cast<size_t>(i * rec.m_mol[num_snap].size() / nt),
        std::min(static_cast<size_t>((i + 1) * rec.m_mol[num_snap].size() / nt),
                 rec.info.size()),
        ref(dmin), ref(CachedData), ref(lg_parameters), ref(rec_parameters),
        ref(Arguments), true, &infoener_res_v[i]);
  for (auto &t : threads)
    t.join();

  for (auto &ii : infoener_res_v) {
    for (auto &ele : ii.position_map) {
      infoener_res.add_energy(ele.first,
                              ii.emap_ee[ele.second],
                              ii.emap_ev[ele.second]);
    }
  }

  wrapper_print_energy(nt, data_lg, data_rec, flag_all_mode, Arguments,
                       Energy_res_mol1, Energy_res_mol2, num_snap, info_output);
}

int energyMapOf2Ligs(const Multi_MOL &rec, const Multi_MOL &lg, size_t num_snap,
                     const vector<d_min> &dmin,
                     const cache_sigma_epsilon &CachedData, size_t nt,
                     const arguments &Arguments,
                     vector<energy_res> &Energy_res_mol1,
                     vector<energy_res> &Energy_res_mol2, info &info_output,
                     const Multi_MOL_parameters &rec_parameters,
                     const Multi_MOL_parameters &lg_parameters,
                     info_ener_res &infoener_res) {
  vector<thread> ener(nt);
  vector<vector<data_out>> data_lg(nt), data_rec(nt);
  vector<info_ener_res> infoener_res_v(nt);

  for (size_t i{}; i < nt; ++i)
    ener[i] = thread(
        imp_cal_energy_lg_atom, ref(data_lg[i]), ref(rec), ref(lg), num_snap,
        static_cast<size_t>(i * lg.m_mol[num_snap].size() / nt),
        std::min(static_cast<size_t>((i + 1) * lg.m_mol[num_snap].size() / nt),
                 lg.info.size()),
        ref(dmin), ref(CachedData), ref(rec_parameters), ref(lg_parameters),
        ref(Arguments), false, &infoener_res_v[i]);
  for (auto &t : ener)
    t.join();

  for (size_t i{}; i < nt; ++i)
    ener[i] = thread(
        imp_cal_energy_lg_atom, ref(data_rec[i]), ref(lg), ref(rec), num_snap,
        static_cast<size_t>(i * rec.m_mol[num_snap].size() / nt),
        std::min(static_cast<size_t>((i + 1) * rec.m_mol[num_snap].size() / nt),
                 rec.info.size()),
        ref(dmin), ref(CachedData), ref(lg_parameters), ref(rec_parameters),
        ref(Arguments), true, &infoener_res_v[i]);
  for (auto &t : ener)
    t.join();

  for (auto &ii : infoener_res_v) {
    for (auto &ele : ii.position_map) {
      infoener_res.add_energy(ele.first,
                              ii.emap_ee[ele.second],
                              ii.emap_ev[ele.second]);
    }
  }

  wrapper_print_energy(1, data_lg, data_rec, true, Arguments, Energy_res_mol1,
                       Energy_res_mol2, num_snap, info_output);

  return 0;
}

int wrap_cal_energy_map(vector<thread> &threads, Multi_MOL &mol1,
                        Multi_MOL &mol2, arguments &Arguments, size_t num_snap,
                        const vector<d_min> &dmin,
                        const cache_sigma_epsilon &CachedData,
                        vector<energy_res> &Energy_res_mol1,
                        vector<energy_res> &Energy_res_mol2, info &info_output,
                        info_ener_res &infoener_res) {
  // cout << Arguments.mol_1_type << "-" <<
  // Arguments.mol_2_type
  //	<<" "<< Arguments.energy_map << " " <<
  // Arguments.emap_type << '\n';
  const bool is_lig1 = (Arguments.mol_1_type == "lig");
  const bool is_lig2 = (Arguments.mol_2_type == "lig");

  if ((is_lig1 && is_lig2) ||
      (Arguments.energy_map == 3 && Arguments.emap_type == "atom-atom")) {
    energyMapOf2Ligs(mol1, mol2, num_snap, dmin, CachedData, Arguments.nt,
                     Arguments, Energy_res_mol1, Energy_res_mol2, info_output,
                     mol1.parameters, mol2.parameters, infoener_res);
  } else if ((!is_lig1 && is_lig2) ||
             (Arguments.energy_map == 3 &&
              Arguments.emap_type == "residue-atom")) {
    energyMapOfProLig(threads, infoener_res, mol1, mol2, num_snap, dmin,
                      CachedData, Arguments.nt, Arguments.emap_level == "all",
                      Arguments, Energy_res_mol1, Energy_res_mol2, info_output,
                      mol1.parameters, mol2.parameters);
  } else if ((is_lig1 && !is_lig2) ||
             (Arguments.energy_map == 3 &&
              Arguments.emap_type == "atom-residue")) {
    energyMapOfProLig(threads, infoener_res, mol2, mol1, num_snap, dmin,
                      CachedData, Arguments.nt, Arguments.emap_level == "all",
                      Arguments, Energy_res_mol1, Energy_res_mol2, info_output,
                      mol2.parameters, mol1.parameters);
  } else {
    energyMapOf2Pro(threads, infoener_res, mol1, mol2, num_snap, dmin,
                    CachedData, Arguments.nt, true, Arguments, Energy_res_mol1,
                    Energy_res_mol2, info_output, mol1.parameters,
                    mol2.parameters);
  }

  return 0;
}
