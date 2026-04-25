#include "stat_ener_res_res.h"
#include "Multi_MOL.h"
#include "info.h"
#include "titan_energy_stats.h"
#include <algorithm>
#include <chrono>
#include <cctype>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <tuple>
#include <vector>

using namespace std;

namespace {

double average_component(double value, size_t num_snap,
                         const arguments &Arguments) {
  return value / (num_snap > 0 ? num_snap : 1) / Arguments.eu_par;
}

void write_matrix_tsv(const string &filename,
                      const titan_energy_stats::EnergyAxis &axis1,
                      const titan_energy_stats::EnergyAxis &axis2,
                      const vector<vector<double>> &matrix, size_t num_snap,
                      const arguments &Arguments) {
  fstream out(filename, ios::out);
  if (!out.is_open())
    return;

  out << "entity";
  for (const auto &entry : axis2.entries)
    out << '\t' << titan_energy_stats::axis_label(entry);
  out << '\n';

  out << fixed << setprecision(6);
  for (size_t i{}; i < matrix.size(); ++i) {
    const string label =
        i < axis1.entries.size() ? titan_energy_stats::axis_label(axis1.entries[i])
                                 : ("row_" + to_string(i));
    out << label;
    for (size_t j{}; j < matrix[i].size(); ++j)
      out << '\t' << average_component(matrix[i][j], num_snap, Arguments);
    out << '\n';
  }
}

void write_total_matrix_tsv(const string &filename,
                            const titan_energy_stats::EnergyAxis &axis1,
                            const titan_energy_stats::EnergyAxis &axis2,
                            const vector<vector<double>> &emap_elec,
                            const vector<vector<double>> &emap_vdw,
                            size_t num_snap, const arguments &Arguments) {
  fstream out(filename, ios::out);
  if (!out.is_open())
    return;

  out << "entity";
  for (const auto &entry : axis2.entries)
    out << '\t' << titan_energy_stats::axis_label(entry);
  out << '\n';

  out << fixed << setprecision(6);
  for (size_t i{}; i < emap_elec.size(); ++i) {
    const string label =
        i < axis1.entries.size() ? titan_energy_stats::axis_label(axis1.entries[i])
                                 : ("row_" + to_string(i));
    out << label;
    for (size_t j{}; j < emap_elec[i].size(); ++j)
      out << '\t'
          << average_component(emap_elec[i][j] + emap_vdw[i][j], num_snap,
                               Arguments);
    out << '\n';
  }
}

tuple<size_t, double, double>
find_entity_energy(const info_ener_1mol_res &stats, size_t chain_id,
                   size_t key_id) {
  const auto it = stats.info_energy_1mol_respair.find(
      info_ener_1mol(chain_id, key_id));
  if (it == stats.info_energy_1mol_respair.end())
    return {0u, 0.0, 0.0};
  return it->second;
}

double pdb_bfactor_for_atom(const Multi_MOL &mol,
                            const info_ener_1mol_res &stats,
                            size_t atom_index0, size_t num_snap,
                            const arguments &Arguments) {
  const auto key = titan_energy_stats::axis_lookup_key_for_atom(mol, atom_index0);
  const auto [count, elec, vdw] =
      find_entity_energy(stats, key.first, key.second);
  (void)count;
  return average_component(elec + vdw, num_snap, Arguments);
}

string pdb_chain_label(const mol_info &info) {
  const auto &chain = info.get_chain();
  if (chain.empty() || std::isdigit(static_cast<unsigned char>(chain[0])))
    return " ";
  return chain;
}

void write_pdb_atoms(fstream &output, const Multi_MOL &mol,
                     const info_ener_1mol_res &stats, size_t num_snap,
                     const arguments &Arguments) {
  if (mol.m_mol.empty())
    return;
  size_t num{};
  for (size_t atom{}; atom < mol.m_mol[0].size() && atom < mol.info.size();
       ++atom) {
    atom > 99998 ? num = 99999 : num = atom;
    const string chain = pdb_chain_label(mol.info[atom]);
    output << setw(4) << "ATOM" << setw(7) << num + 1 << setw(5)
           << mol.info[atom].get_atom() << " " << setw(4) << left
           << mol.info[atom].get_original_res_name() << right << setw(1)
           << chain << setw(4) << mol.info[atom].get_resnr() << setw(12)
           << fixed << mol.m_mol[0].pos.x[atom] << setw(8) << fixed
           << mol.m_mol[0].pos.y[atom] << setw(8) << mol.m_mol[0].pos.z[atom]
           << setw(6) << "1.00" << setw(7)
           << pdb_bfactor_for_atom(mol, stats, atom, num_snap, Arguments)
           << '\n';
  }
}

struct AverageEmapRow {
  size_t left{};
  size_t right{};
  size_t idx{};
};

} // namespace

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
    result = emap_ee.size();
    return emap_ee.size();
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
    result = emap_ee.size();
    return emap_ee.size();
  }
}

int info_ener_res::add_energy(
    pair<info_ener, tuple<size_t, double, double>> &rhs) {
  size_t pos = emap_ee.size();
  tuple<size_t, size_t, size_t, size_t> temp(
      rhs.first.chain_1id, rhs.first.chain_2id, rhs.first.resnum_1,
      rhs.first.resnum_2);
  auto [it, inserted] = position_map.try_emplace(temp, pos);
  if (inserted) {
    emap_count.push_back(get<0>(rhs.second));
    emap_ee.push_back(get<1>(rhs.second));
    emap_ev.push_back(get<2>(rhs.second));
  } else {
    pos = it->second;
    emap_count[pos] += get<0>(rhs.second);
    emap_ee[pos] += get<1>(rhs.second);
    emap_ev[pos] += get<2>(rhs.second);
  }

  return 0;
}

int info_ener_res::add_energy(tuple<size_t, size_t, size_t, size_t> &&ref,
                              double ee, double ev) {
  size_t pos = emap_ee.size();
  auto [it, inserted] = position_map.try_emplace(std::move(ref), pos);
  if (inserted) {
    emap_count.push_back(1u);
    emap_ee.push_back(ee);
    emap_ev.push_back(ev);
  } else {
    pos = it->second;
    ++emap_count[pos];
    emap_ee[pos] += ee;
    emap_ev[pos] += ev;
  }
  return 0;
}

int info_ener_res::add_energy(const tuple<size_t, size_t, size_t, size_t> &ref,
                              double ee, double ev) {
  size_t pos = emap_ee.size();
  auto [it, inserted] = position_map.try_emplace(ref, pos);
  if (inserted) {
    emap_count.push_back(1u);
    emap_ee.push_back(ee);
    emap_ev.push_back(ev);
  } else {
    pos = it->second;
    ++emap_count[pos];
    emap_ee[pos] += ee;
    emap_ev[pos] += ev;
  }
  return 0;
}

int info_ener_res::wrapper_impl_ener_maptype(Multi_MOL &mol1, Multi_MOL &mol2,
                                             info &output_info,
                                             const arguments &Arguments,
                                             size_t num_snap,
                                             const std::string &prefix) {
  vector<vector<double>> emap_elec, emap_vdw;
  info_ener_1mol_res mol1_ener_res, mol2_ener_res;

  this->cal_energy_maptype(mol1, mol2, emap_elec, emap_vdw, mol1_ener_res,
                           mol2_ener_res, Arguments, num_snap);

  this->print_energy_maptype(mol1, mol2, emap_elec, emap_vdw, num_snap,
                             output_info, Arguments, prefix);

  this->print_energy_pdbtype(mol1, mol2, mol1_ener_res, mol2_ener_res, num_snap,
                             output_info, Arguments, prefix);

  return 0;
}

int info_ener_res::build_pos_map(Multi_MOL &mol1, Multi_MOL &mol2,
                                 vector<vector<size_t>> &pos_map,
                                 const arguments &Arguments) {
  (void)Arguments;
  const auto axis1 = titan_energy_stats::build_energy_axis(mol1);
  const auto axis2 = titan_energy_stats::build_energy_axis(mol2);
  for (size_t pos1{}; pos1 < axis1.entries.size() && pos1 < pos_map.size();
       ++pos1) {
    for (size_t pos2{}; pos2 < axis2.entries.size() &&
                       pos2 < pos_map[pos1].size();
         ++pos2) {
      const auto &left = axis1.entries[pos1];
      const auto &right = axis2.entries[pos2];
      size_t pos{emap_ee.size()};
      find(tuple<size_t, size_t, size_t, size_t>(
               left.chain_id, right.chain_id, left.key_id, right.key_id),
           0, emap_ee.size(), pos);
      pos_map[pos1][pos2] = pos;
    }
  }
  return 0;
}

size_t info_ener_res::determine_pos(Multi_MOL &mol1, Multi_MOL &mol2,
                                    const arguments &Arguments, size_t i,
                                    size_t j, size_t l, size_t k) {
  (void)Arguments;
  size_t pos{};
  if (i >= mol1.Index.size() || k >= mol2.Index.size() ||
      j >= mol1.Index[i].get_size_v() || l >= mol2.Index[k].get_size_v())
    return emap_ee.size();

  const size_t atom1 = mol1.Index[i].get_res_pos_e(j);
  const size_t atom2 = mol2.Index[k].get_res_pos_e(l);
  if (atom1 == 0 || atom2 == 0)
    return emap_ee.size();

  const auto key1 = titan_energy_stats::axis_lookup_key_for_atom(mol1, atom1 - 1);
  const auto key2 = titan_energy_stats::axis_lookup_key_for_atom(mol2, atom2 - 1);
  find(tuple<size_t, size_t, size_t, size_t>(key1.first, key2.first,
                                             key1.second, key2.second),
       0, emap_ee.size(), pos);
  return pos;
}

int info_ener_res::cal_energy_maptype(Multi_MOL &mol1, Multi_MOL &mol2,
                                      vector<vector<double>> &emap_elec,
                                      vector<vector<double>> &emap_vdw,
                                      info_ener_1mol_res &mol1_ener_res,
                                      info_ener_1mol_res &mol2_ener_res,
                                      const arguments &Arguments,
                                      size_t num_snap) {
  (void)Arguments;
  (void)num_snap;
  if (emap_ee.empty()) {
    cout << "info energy is empty" << std::endl;
    return 1;
  }

  const auto axis1 = titan_energy_stats::build_energy_axis(mol1);
  const auto axis2 = titan_energy_stats::build_energy_axis(mol2);
  const size_t num_res_1 = axis1.entries.size();
  const size_t num_res_2 = axis2.entries.size();

  if (num_res_1 == 0 || num_res_2 == 0) {
    return 1;
  }

  emap_elec.assign(num_res_1, vector<double>(num_res_2, 0.0));
  emap_vdw.assign(num_res_1, vector<double>(num_res_2, 0.0));

  size_t match_count = 0;
  size_t miss_count = 0;
  for (const auto &[key, pos] : position_map) {
    if (pos >= emap_ee.size())
      continue;

    const size_t chain1 = get<0>(key);
    const size_t chain2 = get<1>(key);
    const size_t entity1 = get<2>(key);
    const size_t entity2 = get<3>(key);

    const auto it1 = axis1.lookup.find(make_pair(chain1, entity1));
    const auto it2 = axis2.lookup.find(make_pair(chain2, entity2));

    if (it1 == axis1.lookup.end() || it2 == axis2.lookup.end()) {
      ++miss_count;
      continue;
    }

    ++match_count;
    const size_t p1 = it1->second;
    const size_t p2 = it2->second;
    if (p1 >= emap_elec.size() || p2 >= emap_elec[p1].size())
      continue;

    emap_elec[p1][p2] += emap_ee[pos];
    emap_vdw[p1][p2] += emap_ev[pos];

    const auto &entry1 = axis1.entries[p1];
    const auto &entry2 = axis2.entries[p2];
    auto &left =
        mol1_ener_res.info_energy_1mol_respair[info_ener_1mol(
            entry1.chain_id, entry1.key_id)];
    ++get<0>(left);
    get<1>(left) += emap_ee[pos];
    get<2>(left) += emap_ev[pos];

    auto &right =
        mol2_ener_res.info_energy_1mol_respair[info_ener_1mol(
            entry2.chain_id, entry2.key_id)];
    ++get<0>(right);
    get<1>(right) += emap_ee[pos];
    get<2>(right) += emap_ev[pos];
  }
  if (miss_count > 0) {
    cout << "WARNING: skipped " << miss_count
         << " energy-map entries that did not match current molecule axes; "
         << "matched " << match_count << ".\n";
  }
  return 0;
}

int info_ener_res::print_energy_maptype(Multi_MOL &mol1, Multi_MOL &mol2,
                                        vector<vector<double>> &emap_elec,
                                        vector<vector<double>> &emap_vdw,
                                        size_t num_snap, info &output_info,
                                        const arguments &Arguments,
                                        const std::string &prefix) {
  const auto axis1 = titan_energy_stats::build_energy_axis(mol1);
  const auto axis2 = titan_energy_stats::build_energy_axis(mol2);
  write_matrix_tsv(prefix + "electrostatic" +
                       output_info.getOuputEnerMaptypeAve(),
                   axis1, axis2, emap_elec, num_snap, Arguments);
  write_matrix_tsv(prefix + "vanderWaals" + output_info.getOuputEnerMaptypeAve(),
                   axis1, axis2, emap_vdw, num_snap, Arguments);
  write_total_matrix_tsv(prefix + "total" + output_info.getOuputEnerMaptypeAve(),
                         axis1, axis2, emap_elec, emap_vdw, num_snap,
                         Arguments);
  write_total_matrix_tsv(prefix + output_info.getOutputAverageEmapMatrix(),
                         axis1, axis2, emap_elec, emap_vdw, num_snap,
                         Arguments);
  return 0;
}

int info_ener_res::print_energy_pdbtype(Multi_MOL &mol1, Multi_MOL &mol2,
                                        info_ener_1mol_res &mol1_ener_res,
                                        info_ener_1mol_res &mol2_ener_res,
                                        size_t num_snap, info &output_info,
                                        const arguments &Arguments,
                                        const std::string &prefix) {
  fstream output(prefix + "total" + output_info.getOutputEnerPdbtypeAve(),
                 ios::out);
  output.precision(3);
  write_pdb_atoms(output, mol1, mol1_ener_res, num_snap, Arguments);
  write_pdb_atoms(output, mol2, mol2_ener_res, num_snap, Arguments);
  output.close();
  return 0;
}

int info_ener_res::print_average_emap(Multi_MOL &mol1, Multi_MOL &mol2,
                                      size_t num_snap, info &output_info,
                                      const arguments &Arguments,
                                      const std::string &prefix) {
  fstream output(prefix + output_info.getOutputAverageEmap(), ios::out);
  if (!output.is_open())
    return 1;

  const auto axis1 = titan_energy_stats::build_energy_axis(mol1);
  const auto axis2 = titan_energy_stats::build_energy_axis(mol2);
  vector<AverageEmapRow> rows;
  rows.reserve(position_map.size());
  for (const auto &pair : position_map) {
    const auto key = pair.first;
    const auto it1 = axis1.lookup.find(make_pair(get<0>(key), get<2>(key)));
    const auto it2 = axis2.lookup.find(make_pair(get<1>(key), get<3>(key)));
    if (it1 == axis1.lookup.end() || it2 == axis2.lookup.end())
      continue;
    if (pair.second >= emap_ee.size() || pair.second >= emap_ev.size())
      continue;
    rows.push_back({it1->second, it2->second, pair.second});
  }
  sort(rows.begin(), rows.end(), [](const AverageEmapRow &a,
                                    const AverageEmapRow &b) {
    return tie(a.left, a.right, a.idx) < tie(b.left, b.right, b.idx);
  });

  output << "left_label\tleft_kind\tleft_chain_id\tleft_chain\tleft_key_id"
         << "\tleft_resnr\tleft_internal_resnr\tleft_atom_index"
         << "\tright_label\tright_kind\tright_chain_id\tright_chain"
         << "\tright_key_id\tright_resnr\tright_internal_resnr"
         << "\tright_atom_index\tavg_elec_" << Arguments.eu << "\tavg_vdw_"
         << Arguments.eu << "\tavg_total_" << Arguments.eu
         << "\tn_frames\tn_pairs\n";
  output << fixed << setprecision(6);

  for (const auto &row : rows) {
    const auto &left = axis1.entries[row.left];
    const auto &right = axis2.entries[row.right];
    const double ee = average_component(emap_ee[row.idx], num_snap, Arguments);
    const double ev = average_component(emap_ev[row.idx], num_snap, Arguments);
    const size_t count = row.idx < emap_count.size() ? emap_count[row.idx] : 0u;
    output << titan_energy_stats::axis_label(left) << '\t'
           << (left.atom_level ? "atom" : "residue") << '\t' << left.chain_id
           << '\t' << titan_energy_stats::tsv_escape(left.chain) << '\t'
           << left.key_id << '\t' << left.resnr << '\t'
           << left.internal_resnr << '\t' << left.atom_index1 << '\t'
           << titan_energy_stats::axis_label(right) << '\t'
           << (right.atom_level ? "atom" : "residue") << '\t'
           << right.chain_id << '\t'
           << titan_energy_stats::tsv_escape(right.chain) << '\t'
           << right.key_id << '\t' << right.resnr << '\t'
           << right.internal_resnr << '\t' << right.atom_index1 << '\t' << ee
           << '\t' << ev << '\t' << ee + ev << '\t'
           << (num_snap > 0 ? num_snap : 1) << '\t' << count << '\n';
  }

  output.close();
  return 0;
}
