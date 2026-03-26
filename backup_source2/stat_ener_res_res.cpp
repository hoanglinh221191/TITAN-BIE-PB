#include "stat_ener_res_res.h"
#include "Multi_MOL.h"
#include "info.h"
#include <chrono>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;

/*bool info_ener_res::check_in(pair<info_ener, tuple<size_t, double, double>>&
rhs)
{
        for (auto& ele : info_energy_respair)
        {
                if (ele.first.resnum_1 == rhs.first.resnum_1 &&
ele.first.resnum_2 == rhs.first.resnum_2
                        && ele.first.chain_1 == rhs.first.chain_1 &&
ele.first.chain_2 == rhs.first.chain_2
                        && ele.first.chain_1id == rhs.first.chain_1id &&
ele.first.chain_2id == rhs.first.chain_2id) return true; else continue;
        }
        return false;
}

bool info_ener_res::check_in(info_ener& rhs)
{
        for (auto& ele : info_energy_respair)
        {
                if (ele.first.resnum_1 == rhs.resnum_1 && ele.first.resnum_2 ==
rhs.resnum_2
                        && ele.first.chain_1 == rhs.chain_1 && ele.first.chain_2
== rhs.chain_2
                        && ele.first.chain_1id == rhs.chain_1id &&
ele.first.chain_2id == rhs.chain_2id) return true; else continue;
        }
        return false;
}*/

size_t info_ener_res::find(tuple<size_t, size_t, size_t, size_t> &&postemp,
                           size_t begin, size_t end, size_t &result) {
  if (position_map.find(postemp) != position_map.end()) {
    result = position_map[postemp];
    return position_map[postemp];
  } else {
    result = info_energy_respair.size();
    return result;
  }
}

size_t info_ener_res::find(const info_ener &rhs, size_t begin, size_t end,
                           size_t &result) {
  tuple<size_t, size_t, size_t, size_t> postemp(rhs.chain_1id, rhs.chain_2id,
                                                rhs.resnum_1, rhs.resnum_2);
  if (position_map.find(postemp) != position_map.end()) {
    result = position_map[postemp];
    return position_map[postemp];
  } else {
    result = info_energy_respair.size();
    return info_energy_respair.size();
  }
}

int info_ener_res::add_energy(
    pair<info_ener, tuple<size_t, double, double>> &rhs) {
  size_t pos{info_energy_respair.size()};
  tuple<size_t, size_t, size_t, size_t> temp(
      rhs.first.chain_1id, rhs.first.chain_2id, rhs.first.resnum_1,
      rhs.first.resnum_2);
  if (position_map.try_emplace(temp, pos).second == true) {
    info_energy_respair.emplace_back(rhs.second);
  } else {
    pos = position_map[temp];
    get<0>(info_energy_respair[pos]) += get<0>(rhs.second);
    get<1>(info_energy_respair[pos]) += get<1>(rhs.second);
    get<2>(info_energy_respair[pos]) += get<2>(rhs.second);
  }

  return 0;
}

int info_ener_res::add_energy(tuple<size_t, size_t, size_t, size_t> &&ref,
                              double ee, double ev) {
  size_t pos{info_energy_respair.size()};
  if (position_map.try_emplace(ref, pos).second == false) {
    pos = position_map[ref];
    ++get<0>(info_energy_respair[pos]);
    get<1>(info_energy_respair[pos]) += ee;
    get<2>(info_energy_respair[pos]) += ev;
  } else {
    info_energy_respair.emplace_back(tuple<size_t, double, double>(1u, ee, ev));
  }
  return 0;
}

int info_ener_res::add_energy(const tuple<size_t, size_t, size_t, size_t> &ref,
                              double ee, double ev) {
  size_t pos{info_energy_respair.size()};
  if (position_map.try_emplace(ref, pos).second == false) {
    pos = position_map[ref];
    ++get<0>(info_energy_respair[pos]);
    get<1>(info_energy_respair[pos]) += ee;
    get<2>(info_energy_respair[pos]) += ev;
  } else {
    info_energy_respair.emplace_back(tuple<size_t, double, double>(1u, ee, ev));
  }
  return 0;
}

int info_ener_res::wrapper_impl_ener_maptype(Multi_MOL &mol1, Multi_MOL &mol2,
                                             info &output_info,
                                             const arguments &Arguments,
                                             size_t num_snap) {
  vector<vector<double>> emap_elec, emap_vdw;
  info_ener_1mol_res mol1_ener_res, mol2_ener_res;

  this->cal_energy_maptype(mol1, mol2, emap_elec, emap_vdw, mol1_ener_res,
                           mol2_ener_res, Arguments, num_snap);

  this->print_energy_maptype(mol1, mol2, emap_elec, emap_vdw, num_snap,
                             output_info, Arguments);

  this->print_energy_pdbtype(mol1, mol2, mol1_ener_res, mol2_ener_res, num_snap,
                             output_info, Arguments);

  return 0;
}

int info_ener_res::build_pos_map(Multi_MOL &mol1, Multi_MOL &mol2,
                                 vector<vector<size_t>> &pos_map,
                                 const arguments &Arguments) {
  size_t pos1{};
  for (size_t i{}; i < mol1.Index.size() - 1; ++i)
    for (size_t j{}; j < mol1.Index[i].get_size_v(); ++j) {
      size_t pos2{};
      for (size_t k{}; k < mol2.Index.size() - 1; ++k)
        for (size_t l{}; l < mol2.Index[k].get_size_v(); ++l) {
          if (pos1 >= pos_map.size() || pos2 >= pos_map[pos1].size()) {
            break;
          }

          ;

          size_t int_resnr1 = mol1.info[mol1.Index[i].get_res_pos_e(j) - 1]
                                  .get_internal_resnr();
          size_t int_resnr2 = mol2.info[mol2.Index[k].get_res_pos_e(l) - 1]
                                  .get_internal_resnr();
          vector<size_t> find_result(Arguments.nt);
          find(tuple<size_t, size_t, size_t, size_t>(
                   mol1.Index[i].get_chainid(), mol2.Index[k].get_chainid(),
                   int_resnr1, int_resnr2),
               0, info_energy_respair.size(), ref(find_result[i]));

          size_t pos{info_energy_respair.size()};
          for (auto &elec : find_result)
            if (elec < info_energy_respair.size())
              pos = elec;

          pos_map[pos1][pos2] = pos;

          ++pos2;
        }
      ++pos1;
    }
  return 0;
}

size_t info_ener_res::determine_pos(Multi_MOL &mol1, Multi_MOL &mol2,
                                    const arguments &Arguments, size_t i,
                                    size_t j, size_t l, size_t k) {
  size_t pos{};
  size_t int_resnr1 =
      mol1.info[mol1.Index[i].get_res_pos_e(j) - 1].get_internal_resnr();
  size_t int_resnr2 =
      mol2.info[mol2.Index[k].get_res_pos_e(l) - 1].get_internal_resnr();
  find(make_tuple(mol1.Index[i].get_chainid(), mol2.Index[k].get_chainid(),
                  int_resnr1, int_resnr2),
       0, info_energy_respair.size(), pos);
  return pos;
}

int info_ener_res::cal_energy_maptype(Multi_MOL &mol1, Multi_MOL &mol2,
                                      vector<vector<double>> &emap_elec,
                                      vector<vector<double>> &emap_vdw,
                                      info_ener_1mol_res &mol1_ener_res,
                                      info_ener_1mol_res &mol2_ener_res,
                                      const arguments &Arguments,
                                      size_t num_snap) {
  if (info_energy_respair.empty()) {
    cout << "info energy is empty" << std::endl;
    return 1;
  }

  size_t num_res_1{}, num_res_2{};

  for (size_t i{}; i < mol1.Index.size() - 1; ++i)
    num_res_1 += mol1.Index[i].get_size_v();

  for (size_t i{}; i < mol2.Index.size() - 1; ++i)
    num_res_2 += mol2.Index[i].get_size_v();

  if (num_res_1 == 0 || num_res_2 == 0) {
    return 1;
  }

  emap_elec = vector<vector<double>>(num_res_1);
  emap_vdw = vector<vector<double>>(num_res_1);
  for (size_t i{}; i < emap_elec.size(); ++i) {
    emap_elec[i] = vector<double>(num_res_2, 0.);
    emap_vdw[i] = vector<double>(num_res_2, 0.);
  }

  // vector<thread> threads(Arguments.nt);

  struct fast_info {
    string chain;
    size_t chainid;
    size_t res_num;
  };
  vector<fast_info> fast1, fast2;
  unordered_map<pair<size_t, size_t>, size_t, pair_hash> lookup1, lookup2;

  size_t pos1{};
  for (size_t i{}; i < mol1.Index.size() - 1; ++i) {
    for (size_t j{}; j < mol1.Index[i].get_size_v(); ++j) {
      size_t int_resnr =
          mol1.info[mol1.Index[i].get_res_pos_e(j) - 1].get_internal_resnr();
      fast1.push_back({mol1.Index[i].get_chain(), mol1.Index[i].get_chainid(),
                       mol1.Index[i].get_res_num_e(j)});
      lookup1[make_pair(mol1.Index[i].get_chainid(), int_resnr)] = pos1++;
    }
  }

  size_t pos2{};
  for (size_t k{}; k < mol2.Index.size() - 1; ++k) {
    for (size_t l{}; l < mol2.Index[k].get_size_v(); ++l) {
      size_t int_resnr =
          mol2.info[mol2.Index[k].get_res_pos_e(l) - 1].get_internal_resnr();
      fast2.push_back({mol2.Index[k].get_chain(), mol2.Index[k].get_chainid(),
                       mol2.Index[k].get_res_num_e(l)});
      lookup2[make_pair(mol2.Index[k].get_chainid(), int_resnr)] = pos2++;
    }
  }

  size_t match_count = 0;
  size_t miss_count = 0;
  for (const auto &[key, pos] : position_map) {
    if (pos >= info_energy_respair.size())
      continue;

    size_t chain1 = get<0>(key);
    size_t chain2 = get<1>(key);
    size_t res1 = get<2>(key);
    size_t res2 = get<3>(key);

    auto it1 = lookup1.find(make_pair(chain1, res1));
    auto it2 = lookup2.find(make_pair(chain2, res2));

    if (it1 != lookup1.end() && it2 != lookup2.end()) {
      match_count++;
      size_t p1 = it1->second;
      size_t p2 = it2->second;

      emap_elec[p1][p2] += get<1>(info_energy_respair[pos]);
      emap_vdw[p1][p2] += get<2>(info_energy_respair[pos]);

      info_ener_1mol ener_1mol_temp(fast1[p1].chainid, fast1[p1].res_num);
      ++get<0>(mol1_ener_res.info_energy_1mol_respair[ener_1mol_temp]);
      get<1>(mol1_ener_res.info_energy_1mol_respair[ener_1mol_temp]) +=
          get<1>(info_energy_respair[pos]);
      get<2>(mol1_ener_res.info_energy_1mol_respair[ener_1mol_temp]) +=
          get<2>(info_energy_respair[pos]);

      info_ener_1mol ener_1mol_temp2(fast2[p2].chainid, fast2[p2].res_num);
      ++get<0>(mol2_ener_res.info_energy_1mol_respair[ener_1mol_temp2]);
      get<1>(mol2_ener_res.info_energy_1mol_respair[ener_1mol_temp2]) +=
          get<1>(info_energy_respair[pos]);
      get<2>(mol2_ener_res.info_energy_1mol_respair[ener_1mol_temp2]) +=
          get<2>(info_energy_respair[pos]);
    } else {
      miss_count++;
      if (miss_count < 20) {
        std::cout << " (it1 matched? " << (it1 != lookup1.end()) << ") "
                  << "chain2=" << chain2 << " res2=" << res2
                  << " (it2 matched? " << (it2 != lookup2.end()) << ")"
                  << std::endl;
      }
    }
  }

  return 0;
}

int info_ener_res::print_energy_maptype(Multi_MOL &mol1, Multi_MOL &mol2,
                                        vector<vector<double>> &emap_elec,
                                        vector<vector<double>> &emap_vdw,
                                        size_t num_snap, info &output_info,
                                        const arguments &Arguments) {
  fstream output_elec("electrostatic" + output_info.getOuputEnerMaptypeAve(),
                      ios::out),
      output_vdw("vanderWaals" + output_info.getOuputEnerMaptypeAve(),
                 ios::out),
      output_total("total" + output_info.getOuputEnerMaptypeAve(), ios::out);

  output_elec.precision(4);
  output_vdw.precision(4);
  output_total.precision(4);
  for (size_t i{}; i < emap_elec.size(); ++i) {
    for (size_t j{}; j < emap_elec[i].size(); ++j) {
      output_elec << fixed << left << setw(9)
                  << emap_elec[i][j] / num_snap / Arguments.eu_par << " ";
      output_vdw << fixed << left << setw(9)
                 << emap_vdw[i][j] / num_snap / Arguments.eu_par << " ";
      output_total << fixed << left << setw(9)
                   << (emap_elec[i][j] + emap_vdw[i][j]) / num_snap /
                          Arguments.eu_par
                   << " ";
    }
    output_elec << '\n';
    output_vdw << '\n';
    output_total << '\n';
  }

  fstream output_matrix(output_info.getOutputAverageEmapMatrix(), ios::out);
  if (output_matrix.is_open()) {
    size_t max_internal_1 = 0;
    for (auto &a : mol1.info) {
      if (a.get_internal_resnr() > max_internal_1)
        max_internal_1 = a.get_internal_resnr();
    }
    vector<size_t> int_to_resnr1(max_internal_1 + 1, 0);
    for (auto &a : mol1.info)
      int_to_resnr1[a.get_internal_resnr()] = a.get_resnr();

    size_t max_internal_2 = 0;
    for (auto &a : mol2.info) {
      if (a.get_internal_resnr() > max_internal_2)
        max_internal_2 = a.get_internal_resnr();
    }
    vector<size_t> int_to_resnr2(max_internal_2 + 1, 0);
    for (auto &a : mol2.info)
      int_to_resnr2[a.get_internal_resnr()] = a.get_resnr();

    vector<size_t> res_nums1;
    for (size_t i{}; i < mol1.Index.size() - 1; ++i) {
      for (size_t j{}; j < mol1.Index[i].get_size_v(); ++j) {
        res_nums1.push_back(mol1.Index[i].get_res_num_e(j));
      }
    }

    bool is_lig = (mol2.get_mol_type() == "ligand");
    vector<size_t> res_nums2;
    for (size_t k{}; k < mol2.Index.size() - 1; ++k) {
      for (size_t l{}; l < mol2.Index[k].get_size_v(); ++l) {
        if (is_lig) {
          res_nums2.push_back(mol2.Index[k].get_res_pos_e(l));
        } else {
          res_nums2.push_back(mol2.Index[k].get_res_num_e(l));
        }
      }
    }

    output_matrix.precision(4);
    output_matrix << left << setw(9) << " ";
    if (emap_elec.size() > 0 && emap_elec[0].size() > 0 &&
        res_nums2.size() >= emap_elec[0].size()) {
      for (size_t j{}; j < emap_elec[0].size(); ++j) {
        output_matrix << left << setw(9) << ("#" + to_string(res_nums2[j]))
                      << " ";
      }
    }
    output_matrix << '\n';

    for (size_t i{}; i < emap_elec.size(); ++i) {
      if (i < res_nums1.size()) {
        output_matrix << left << setw(9) << ("#" + to_string(res_nums1[i]))
                      << " ";
      } else {
        output_matrix << left << setw(9) << ("#?") << " ";
      }
      for (size_t j{}; j < emap_elec[i].size(); ++j) {
        output_matrix << fixed << left << setw(9)
                      << (emap_elec[i][j] + emap_vdw[i][j]) / num_snap /
                             Arguments.eu_par
                      << " ";
      }
      output_matrix << '\n';
    }
    output_matrix.close();
  }

  return 0;
}

int info_ener_res::print_energy_pdbtype(Multi_MOL &mol1, Multi_MOL &mol2,
                                        info_ener_1mol_res &mol1_ener_res,
                                        info_ener_1mol_res &mol2_ener_res,
                                        size_t num_snap, info &output_info,
                                        const arguments &Arguments) {
  fstream output("total" + output_info.getOutputEnerPdbtypeAve(), ios::out);
  output.precision(3);

  size_t num{};
  for (size_t l{}; l < mol1.m_mol[0].Molecule.size(); ++l) {
    l > 99998 ? num = 99999 : num = l;
    size_t factor{get<0>(mol1_ener_res.info_energy_1mol_respair[info_ener_1mol(
        mol1.info[l].get_chainnum(), mol1.info[l].get_resnr())])};
    if (factor == 0 || Arguments.mol_1_type.find("protein") != string::npos)
      factor = 1;
    string chain{" "};
    if (!isdigit(mol1.info[l].get_chain()[0]))
      chain = mol1.info[l].get_chain();

    output << setw(4) << "ATOM" << setw(7) << num + 1 << setw(5)
           << mol1.info[l].get_atom() << " " << setw(4) << left
           << mol1.info[l].get_original_res_name() << right << setw(1) << chain
           << setw(4) << mol1.info[l].get_resnr() << setw(12) << fixed
           << mol1.m_mol[0].pos.x[l] << setw(8) << fixed
           << mol1.m_mol[0].pos.y[l] << setw(8)
           << mol1.m_mol[0].pos.z[l] << setw(6) << "1.00" << setw(7)
           << (get<1>(mol1_ener_res.info_energy_1mol_respair[info_ener_1mol(
                   mol1.info[l].get_chainnum(), mol1.info[l].get_resnr())]) /
                   num_snap +
               get<2>(mol1_ener_res.info_energy_1mol_respair[info_ener_1mol(
                   mol1.info[l].get_chainnum(), mol1.info[l].get_resnr())]) /
                   num_snap) /
                  Arguments.eu_par / factor
           << '\n';
  }

  for (size_t l{}; l < mol2.m_mol[0].Molecule.size(); ++l) {
    size_t factor{1};
    if (Arguments.mol_2_type == "lig")
      factor = get<0>(mol2_ener_res.info_energy_1mol_respair[info_ener_1mol(
                   mol2.info[l].get_chainnum(), mol2.info[l].get_resnr())]) /
               get<0>(mol1_ener_res.info_energy_1mol_respair[info_ener_1mol(
                   mol1.info[0].get_chainnum(), mol1.info[0].get_resnr())]);

    factor *= get<0>(mol2_ener_res.info_energy_1mol_respair[info_ener_1mol(
        mol2.info[l].get_chainnum(), mol2.info[l].get_resnr())]);

    if (factor == 0 || Arguments.mol_2_type.find("protein") != string::npos ||
        get<0>(mol2_ener_res.info_energy_1mol_respair[info_ener_1mol(
            mol2.info[l].get_chainnum(), mol2.info[l].get_resnr())]) *
                factor <=
            0)
      factor = 1;

    l > 99998 ? num = 99999 : num = l;
    string chain{" "};
    if (!isdigit(mol2.info[l].get_chain()[0]))
      chain = mol2.info[l].get_chain();
    output << setw(4) << "ATOM" << setw(7) << num + 1 << setw(5)
           << mol2.info[l].get_atom() << " " << setw(4) << left
           << mol2.info[l].get_original_res_name() << right << setw(1) << chain
           << setw(4) << mol2.info[l].get_resnr() << setw(12) << fixed
           << mol2.m_mol[0].pos.x[l] << setw(8) << fixed
           << mol2.m_mol[0].pos.y[l] << setw(8)
           << mol2.m_mol[0].pos.z[l] << setw(6) << "1.00" << setw(7)
           << (get<1>(mol2_ener_res.info_energy_1mol_respair[info_ener_1mol(
                   mol2.info[l].get_chainnum(), mol2.info[l].get_resnr())]) /
                   num_snap +
               get<2>(mol2_ener_res.info_energy_1mol_respair[info_ener_1mol(
                   mol2.info[l].get_chainnum(), mol2.info[l].get_resnr())]) /
                   num_snap) /
                  Arguments.eu_par / factor
           << '\n';
  }
  output.close();
  return 0;
}

int info_ener_res::print_average_emap(Multi_MOL &mol1, Multi_MOL &mol2,
                                      size_t num_snap, info &output_info,
                                      const arguments &Arguments) {
  fstream output(output_info.getOutputAverageEmap(), ios::out);
  if (!output.is_open())
    return 1;

  size_t max_internal_1 = 0;
  for (auto &a : mol1.info) {
    if (a.get_internal_resnr() > max_internal_1)
      max_internal_1 = a.get_internal_resnr();
  }
  vector<size_t> int_to_resnr1(max_internal_1 + 1, 0);
  for (auto &a : mol1.info)
    int_to_resnr1[a.get_internal_resnr()] = a.get_resnr();

  size_t max_internal_2 = 0;
  for (auto &a : mol2.info) {
    if (a.get_internal_resnr() > max_internal_2)
      max_internal_2 = a.get_internal_resnr();
  }
  vector<size_t> int_to_resnr2(max_internal_2 + 1, 0);
  for (auto &a : mol2.info)
    int_to_resnr2[a.get_internal_resnr()] = a.get_resnr();

  output << "#Res1(Rec) | Res2(Lig) | Avg Elec | Avg VdW | Avg Total\n";
  output.precision(4);

  for (const auto &pair : position_map) {
    size_t mol1_chainid = get<0>(pair.first);
    size_t mol2_chainid = get<1>(pair.first);

    size_t int_res1 = get<2>(pair.first);
    size_t int_res2 = get<3>(pair.first);

    size_t mol1_resnr =
        int_res1 < int_to_resnr1.size() ? int_to_resnr1[int_res1] : int_res1;
    size_t mol2_resnr =
        int_res2 < int_to_resnr2.size() ? int_to_resnr2[int_res2] : int_res2;
    size_t idx = pair.second;

    double ee = get<1>(info_energy_respair[idx]) / num_snap / Arguments.eu_par;
    double ev = get<2>(info_energy_respair[idx]) / num_snap / Arguments.eu_par;
    double total = ee + ev;

    output << left << setw(8) << mol1_resnr << " " << setw(4) << mol1_chainid
           << " | " << left << setw(8) << mol2_resnr << " " << setw(4)
           << mol2_chainid << " | " << right << fixed << setw(10) << ee << " | "
           << right << fixed << setw(10) << ev << " | " << right << fixed
           << setw(10) << total << '\n';
  }
  output.close();
  return 0;
}
