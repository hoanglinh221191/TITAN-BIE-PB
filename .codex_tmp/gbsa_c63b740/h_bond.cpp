#include "h_bond.h"
#include <cmath>
#include <thread>
// #include<mutex>
#include "MOL.h"
#include "functions.h"
#include "miscellaneous.h"
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace std;
double PI = 3.1415926535897;
// std::mutex emplace_mutex_1;
// std::mutex emplace_mutex_2;
double cal_angle_3_atoms(double x1, double y1, double z1, double x2, double y2,
                         double z2, double x3, double y3, double z3);
inline double cal_cos_angle_3_atoms(double x1, double y1, double z1, double x2,
                                    double y2, double z2, double x3, double y3,
                                    double z3) {
  double vector_1[3]{x1 - x2, y1 - y2, z1 - z2};
  double vector_2[3]{x3 - x2, y3 - y2, z3 - z2};
  double d1 = std::sqrt(vector_1[0] * vector_1[0] + vector_1[1] * vector_1[1] +
                        vector_1[2] * vector_1[2]);
  double d2 = std::sqrt(vector_2[0] * vector_2[0] + vector_2[1] * vector_2[1] +
                        vector_2[2] * vector_2[2]);
  return (vector_1[0] * vector_2[0] + vector_1[1] * vector_2[1] +
          vector_1[2] * vector_2[2]) /
         (d1 * d2);
} // Returns cos(angle)
int print_hbond_per_snap(MOL &pro_in, MOL &lig, vector<mol_info> &info_rec,
                         vector<mol_info> &info_lg, const arguments &Arguments);
inline bool check_sidechain_hbond(const arguments &Arguments, string atom) {
  return (!Arguments.flag_sidechain_hbond ||
          (Arguments.flag_sidechain_hbond && atom != "N" && atom != "O"));
}
int h_group::add_groups(const MOL &mol_in, const MOL &ref, size_t i,
                        const size_t remain, char _name, double bond_length,
                        vector<mol_info> &info) {
  name = _name;
  index = i;
  size_t next{i + remain > mol_in.size() ? mol_in.size() : i + remain};
  for (size_t j{i}; j < next; ++j) {
    if (info[j].get_resnr() == info[i].get_resnr()) {
      if ((info[j].get_atom_type()[0] == 'H' ||
           (info[j].get_atom().length() >= 4 && info[j].get_atom()[1] == 'H') ||
           info[j].get_atom()[0] == 'H') &&
          std::sqrt(std::pow(mol_in.pos.x[i] - mol_in.pos.x[j], 2) +
                    std::pow(mol_in.pos.y[i] - mol_in.pos.y[j], 2) +
                    std::pow(mol_in.pos.z[i] - mol_in.pos.z[j], 2)) < bond_length)
        h_index.emplace_back(j);
    } else
      break;
  }
  return 0;
}
int h_groups::find_groups(const MOL &mol_in, const MOL &ref,
                          vector<mol_info> &info, const arguments &Arguments) {
  size_t remain{120};
  for (size_t i{}; i < mol_in.size(); ++i) {
    h_group group_temp;
    if ((info[i].get_atom_type()[0] == 'O' || info[i].get_atom()[0] == 'O') &&
        check_sidechain_hbond(Arguments, info[i].get_atom())) {
      group_temp.add_groups(mol_in, ref, i, remain, 'O', 1.1, info);
      groups.emplace_back(group_temp);
    } else if ((info[i].get_atom_type()[0] == 'N' ||
                info[i].get_atom()[0] == 'N') &&
               check_sidechain_hbond(Arguments, info[i].get_atom())) {
      group_temp.add_groups(mol_in, ref, i, remain, 'N', 1.5, info);
      groups.emplace_back(group_temp);
    } else if (info[i].get_atom_type()[0] == 'F' ||
               info[i].get_atom()[0] == 'F') {
      group_temp.add_groups(mol_in, ref, i, remain, 'F', 1.1, info);
      groups.emplace_back(group_temp);
    } else
      continue;
  }
  return 0;
}
void cal_h_bonds(MOL &pro_in, MOL &lig, size_t begin, size_t end,
                 const arguments &Arguments, h_groups &groups_1,
                 h_groups &groups_2, h_bonds_wrapper &dat_pro_in,
                 h_bonds_wrapper &dat_lig) {
  for (size_t i{begin}; i < end; ++i) {
    for (size_t j{}; j < groups_2.groups.size(); ++j) {
      double d{std::sqrt(std::pow(pro_in.pos.x[groups_1.groups[i].index] -
                                      lig.pos.x[groups_2.groups[j].index],
                                  2) +
                         std::pow(pro_in.pos.y[groups_1.groups[i].index] -
                                      lig.pos.y[groups_2.groups[j].index],
                                  2) +
                         std::pow(pro_in.pos.z[groups_1.groups[i].index] -
                                      lig.pos.z[groups_2.groups[j].index],
                                  2))};

      if ((!groups_1.groups[i].h_index.empty() ||
           !groups_2.groups[j].h_index.empty()) &&
          d <= Arguments.dis_h && d > 0.01) {
        const double cos_cutoff = cos(Arguments.angle_h * PI / 180.0);
        for (size_t k{}; k < groups_1.groups[i].h_index.size(); ++k) {
          // Optimized angle check: cos(angle) <= cos(cutoff)  (since angle >=
          // cutoff and cos is decreasing in 0-180) angle = acos(dot / (d1*d2))
          // dot / (d1*d2) <= cos_cutoff
          // Note: cal_angle_3_atoms returned fabs(acos(...)), implying inputs
          // in range [-1,1] -> [0, PI]. Angle is always positive 0..180.

          // Recalculate dot product and dists inline or call helper
          // Helper:
          // double vector_1[3]{ mol1.x - mol2.x, ... }
          // double vector_2[3]{ mol3.x - mol2.x, ... }
          // dot = v1.v2
          // d1 = mol1.dis(mol2) (already computed? no d is pro[i] to lig[j].
          // Wait. Angle is D-H...A.
          // groups_1[i].index is Donor? h_index[k] is H?
          // pro_in.Molecule[groups_1.groups[i].index] is the heavy atom (D).
          // pro_in.Molecule[groups_1.groups[i].h_index[k]] is the H atom.
          // lig.Molecule[groups_2.groups[j].index] is the Acceptor (A).
          // The angle is usually D-H...A or N-H...O.
          // Angle at H? Or angle at D?
          // cal_angle_3_atoms(D, H, A) -> angle at H (arg2 is center).
          // Standard H-bond definition is usually angle(D-H-A) > 120.

          // Calling original function signature: cal_angle_3_atoms(mol1, mol2,
          // mol3) Inside: center is mol2. Call 1: (D, H, A). Center is H.

          double angle_cos_val = cal_cos_angle_3_atoms(
              pro_in.pos.x[groups_1.groups[i].index],
              pro_in.pos.y[groups_1.groups[i].index],
              pro_in.pos.z[groups_1.groups[i].index],
              pro_in.pos.x[groups_1.groups[i].h_index[k]],
              pro_in.pos.y[groups_1.groups[i].h_index[k]],
              pro_in.pos.z[groups_1.groups[i].h_index[k]],
              lig.pos.x[groups_2.groups[j].index], lig.pos.y[groups_2.groups[j].index],
              lig.pos.z[groups_2.groups[j].index]);
          if (angle_cos_val <= cos_cutoff) {
            // std::lock_guard<std::mutex> lock(emplace_mutex_1);
            //  We need angle in degrees for output? Yes, emplace_back takes
            //  angle. So we must compute acos ONLY if it passes the check.
            double angle_deg = acos(angle_cos_val) * 180.0 / PI;
            dat_lig.dat.emplace_back(
                groups_1.groups[i].index, groups_2.groups[j].index,
                groups_1.groups[i].h_index[k], angle_deg, d);
          }
        }
        for (size_t k{}; k < groups_2.groups[j].h_index.size(); ++k) {
          double angle_cos_val = cal_cos_angle_3_atoms(
              lig.pos.x[groups_2.groups[j].index], lig.pos.y[groups_2.groups[j].index],
              lig.pos.z[groups_2.groups[j].index],
              lig.pos.x[groups_2.groups[j].h_index[k]],
              lig.pos.y[groups_2.groups[j].h_index[k]],
              lig.pos.z[groups_2.groups[j].h_index[k]],
              pro_in.pos.x[groups_1.groups[i].index],
              pro_in.pos.y[groups_1.groups[i].index],
              pro_in.pos.z[groups_1.groups[i].index]);
          if (angle_cos_val <= cos_cutoff) {
            double angle_deg = acos(angle_cos_val) * 180.0 / PI;
            // std::lock_guard<std::mutex> lock(emplace_mutex_2);
            dat_pro_in.dat.emplace_back(
                groups_2.groups[j].index, groups_1.groups[i].index,
                groups_2.groups[j].h_index[k], angle_deg, d);
          }
        }
      }
    }
  }
}
int find_h_bonds(MOL &pro_in, MOL &pro_in_ref, MOL &lig, MOL &lig_ref,
                 h_groups &groups_1, h_groups &groups_2,
                 const arguments &Arguments, size_t num_snap,
                 vector<mol_info> &info_rec, vector<mol_info> &info_lg,
                 vector<vector<size_t>> &Hbonds_vector) {
  if (num_snap == 0) {
    if (Arguments.nt >= 2) {
      thread t1(&h_groups::find_groups, &groups_1, ref(pro_in), ref(pro_in_ref),
                ref(info_rec), ref(Arguments));
      thread t2(&h_groups::find_groups, &groups_2, ref(lig), ref(lig_ref),
                ref(info_lg), ref(Arguments));
      t1.join();
      t2.join();
    } else {
      groups_1.find_groups(pro_in, pro_in_ref, info_rec, Arguments);
      groups_2.find_groups(lig, lig_ref, info_lg, Arguments);
    }
  }

  vector<h_bonds_wrapper> h_bonds_pro{Arguments.nt}, h_bonds_lig{Arguments.nt};

  vector<thread> threads_h(Arguments.nt);
  for (size_t i{}; i < Arguments.nt; ++i) {
    threads_h[i] = thread(cal_h_bonds, ref(pro_in), ref(lig),
                          i * groups_1.groups.size() / Arguments.nt,
                          (i + 1) * groups_1.groups.size() / Arguments.nt,
                          ref(Arguments), ref(groups_1), ref(groups_2),
                          ref(h_bonds_pro[i]), ref(h_bonds_lig[i]));
  }
  for (auto &t : threads_h) {
    t.join();
  }

  for (size_t i{}; i < Arguments.nt; ++i) {
    for (size_t j{}; j < h_bonds_pro[i].dat.size(); ++j)
      pro_in.list_h_bonds.dat.emplace_back(h_bonds_pro[i].dat[j]);
    for (size_t j{}; j < h_bonds_lig[i].dat.size(); ++j)
      lig.list_h_bonds.dat.emplace_back(h_bonds_lig[i].dat[j]);
  }

  // cout << "pos proc..";
  lig.post_proc_hbonds(lig, info_lg, info_rec);
  // cout << " lig done ";
  pro_in.post_proc_hbonds(pro_in, info_rec, info_lg);
  // cout << "beg.. "<<num_snap;
  print_hbond_per_snap(pro_in, lig, info_rec, info_lg, Arguments);
  // cout << " printed  ";
  if (num_snap >= Hbonds_vector.size()) {
    Hbonds_vector.emplace_back(vector<size_t>(
        1, pro_in.list_h_bonds.dat.size() + lig.list_h_bonds.dat.size()));
  } else {
    Hbonds_vector[num_snap].emplace_back(pro_in.list_h_bonds.dat.size() +
                                         lig.list_h_bonds.dat.size());
  }
  // cout << "exported!\n";
  return 0;
}
inline double cal_angle_3_atoms(double x1, double y1, double z1, double x2,
                                double y2, double z2, double x3, double y3,
                                double z3) {
  double vector_1[3]{x1 - x2, y1 - y2, z1 - z2};
  double vector_2[3]{x3 - x2, y3 - y2, z3 - z2};
  double d1 = std::sqrt(vector_1[0] * vector_1[0] + vector_1[1] * vector_1[1] +
                        vector_1[2] * vector_1[2]);
  double d2 = std::sqrt(vector_2[0] * vector_2[0] + vector_2[1] * vector_2[1] +
                        vector_2[2] * vector_2[2]);
  return fabs(acos((vector_1[0] * vector_2[0] + vector_1[1] * vector_2[1] +
                    vector_1[2] * vector_2[2]) /
                   (d1 * d2)) *
              180. / PI);
}
int MOL::post_proc_hbonds(const MOL &mol, vector<mol_info> &info,
                          vector<mol_info> &counter_info) {
  vector<pair<char, int>> max_hbonds;
  max_hbonds.emplace_back('O', 2);
  max_hbonds.emplace_back('N', 3);
  max_hbonds.emplace_back('H', 1);
  max_hbonds.emplace_back('F', 1);
  for (size_t i{}; i < max_hbonds.size(); ++i)
    imp_post_proc(max_hbonds[i], mol, info, counter_info);
  return 0;
}
int MOL::imp_post_proc(const pair<char, size_t> &atom_dat, const MOL &mol,
                       vector<mol_info> &info, vector<mol_info> &counter_info) {
  vector<hbonds_pos> pos_summary;
  for (size_t i{}; i < list_h_bonds.dat.size(); ++i) {
    // if (info.size() <= list_h_bonds.dat[i].get_index_acceptor()) { cout <<
    // info.size() << " smaller than acceptor index " <<
    // list_h_bonds.dat[i].get_index_acceptor() << '\n'; continue; }
    string name{info[list_h_bonds.dat[i].get_index_acceptor()].get_atom_type()};
    int count_appear{};
    size_t index{list_h_bonds.dat[i].get_index_acceptor()};
    if (atom_dat.first == 'H') {
      // if (counter_info.size() <= list_h_bonds.dat[i].get_index_h()) { cout <<
      // counter_info.size() << " smaller than h index " <<
      // list_h_bonds.dat[i].get_index_h() << '\n'; continue; }
      name = counter_info[list_h_bonds.dat[i].get_index_h()].get_atom_type();
      index = list_h_bonds.dat[i].get_index_h();
    }
    if (name.find(atom_dat.first) != string::npos) {
      for (size_t j{}; j < pos_summary.size(); ++j) {
        if (index == pos_summary[j].index_mol) {
          pos_summary[j].index_rep.emplace_back(i);
          ++count_appear;
        } else
          continue;
      }
      if (count_appear > 0)
        count_appear = 0;
      else
        pos_summary.emplace_back(index, i);
    }
  }
  for (size_t i{}; i < pos_summary.size(); ++i) {
    if (pos_summary[i].index_rep.size() > atom_dat.second) {
      while (pos_summary[i].index_rep.size() != atom_dat.second) {
        double max_dis{0.};
        size_t pos_max_dis{}, pos_summary_del{};
        for (size_t j{}; j < pos_summary[i].index_rep.size(); ++j) {
          if (max_dis <
              list_h_bonds.dat[pos_summary[i].index_rep[j]].get_dis()) {
            max_dis = list_h_bonds.dat[pos_summary[i].index_rep[j]].get_dis();
            pos_max_dis = pos_summary[i].index_rep[j];
            pos_summary_del = j;
          }
        }
        // if (pos_summary_del >= pos_summary[i].index_rep.size())cout << "
        // !!!pos_summary ";
        pos_summary[i].index_rep.erase(pos_summary[i].index_rep.begin() +
                                       pos_summary_del);
        // if (pos_max_dis >= list_h_bonds.dat.size())cout << " !!!list h
        // bond.dat ";
        list_h_bonds.dat.erase(list_h_bonds.dat.begin() + pos_max_dis);

        for (size_t k{}; k < pos_summary.size(); ++k)
          for (size_t l{}; l < pos_summary[k].index_rep.size(); ++l)
            if (pos_summary[k].index_rep[l] > pos_max_dis)
              --pos_summary[k].index_rep[l];
      }
    }
  }
  return 0;
}
void h_groups::print_group(int m) {
  fstream groups_out("out_hbonds_groups_" + to_string(m) + ".dat", ios::out);
  for (size_t i{}; i < groups.size(); ++i) {
    groups_out << groups[i].name << " " << groups[i].index << '\n';
    for (size_t j{}; j < groups[i].h_index.size(); ++j)
      groups_out << groups[i].h_index[j] << "   ";
    groups_out << '\n';
  }
}

int stat_hbonds(vector<vector<size_t>> &Hbonds_vector, const info &infos,
                arguments &Arguments, vector<double> &time_of_snaps,
                vector<string> &list_name) {
  if (Arguments.flag_h_bond == 0)
    return 1;
  fstream statHB(infos.get_stathb_name(), ios::out);
  statHB << "#Time (" << Arguments.tu << ")   Number of H-bonds   SD\n";
  for (size_t i{}; i < Hbonds_vector.size(); ++i) {
    double mean{}, sd{};
    for (size_t j{}; j < Hbonds_vector[i].size(); ++j) {
      mean += Hbonds_vector[i][j];
    }
    mean /= Hbonds_vector[i].size();
    for (size_t j{}; j < Hbonds_vector[i].size(); ++j)
      sd += square(Hbonds_vector[i][j] - mean);
    sd = sqrt(sd / Hbonds_vector[i].size());
    statHB << setw(8) << time_of_snaps[i] << setw(12) << mean << setw(12) << sd
           << '\n';
  }
  statHB.close();
  fstream stat_aveHB(infos.get_stat_ave_hb_name(), ios::out);
  if (!Hbonds_vector.empty()) {
    size_t num_traj = 0;
    for (size_t k = 0; k < Hbonds_vector.size(); ++k) {
      if (Hbonds_vector[k].size() > num_traj) {
        num_traj = Hbonds_vector[k].size();
      }
    }
    for (size_t i{}; i < num_traj; ++i) {
      double mean{}, sd{};
      size_t Size{}; // Count valid trajectories that have this column `i`

      for (size_t j{}; j < Hbonds_vector.size(); ++j) {
        if (i < Hbonds_vector[j].size()) {
          mean += Hbonds_vector[j][i];
          Size++;
        }
      }

      if (Size > 0) {
        mean /= Size;
        double variance_sum = 0.0;
        for (size_t j{}; j < Hbonds_vector.size(); ++j) {
          if (i < Hbonds_vector[j].size()) {
            variance_sum += square(Hbonds_vector[j][i] - mean);
          }
        }
        string traj_name = (list_name.size() > i)
                               ? list_name[i]
                               : ("Traj_" + to_string(i + 1));
        stat_aveHB << setw(15) << traj_name << setw(12) << mean << setw(12)
                   << sqrt(variance_sum / Size) << '\n';
      }
    }
  }
  stat_aveHB.close();
  return 0;
}

int print_hbond_per_snap(MOL &pro_in, MOL &lig, vector<mol_info> &info_rec,
                         vector<mol_info> &info_lg,
                         const arguments &Arguments) {
  fstream out_h_bonds("out_hydrogen_bonds.dat", ios::app);
  for (size_t i{}; i < pro_in.list_h_bonds.dat.size(); ++i) {
    out_h_bonds
        << setw(4) << left
        << info_rec[pro_in.list_h_bonds.dat[i].get_index_acceptor()].get_atom()
        << setw(3) << right
        << info_rec[pro_in.list_h_bonds.dat[i].get_index_acceptor()].get_chain()
        << setw(5) << right
        << info_rec[pro_in.list_h_bonds.dat[i].get_index_acceptor()].get_resnr()
        << setw(6) << right
        << info_lg[pro_in.list_h_bonds.dat[i].get_index_h()].get_atom()
        << setw(5) << right
        << info_lg[pro_in.list_h_bonds.dat[i].get_index_h()].get_resnr()
        << '\n';
  }
  for (size_t i{}; i < lig.list_h_bonds.dat.size(); ++i) {
    if (lig.list_h_bonds.dat[i].get_index_h() >= info_rec.size())
      cout << "info_rec.size smaller than index h\n";
    if (lig.list_h_bonds.dat[i].get_index_acceptor() > info_lg.size())
      cout << "info_lg size smaller than index acc\n";

    out_h_bonds
        << setw(4) << left
        << info_lg[lig.list_h_bonds.dat[i].get_index_acceptor()].get_atom()
        << setw(3) << right
        << info_lg[lig.list_h_bonds.dat[i].get_index_acceptor()].get_chain()
        << setw(5) << right
        << info_lg[lig.list_h_bonds.dat[i].get_index_acceptor()].get_resnr()
        << setw(6) << right
        << info_rec[lig.list_h_bonds.dat[i].get_index_h()].get_atom() << setw(5)
        << right << info_rec[lig.list_h_bonds.dat[i].get_index_h()].get_resnr()
        << '\n';
  }
  out_h_bonds.close();
  fstream out_h_bonds_rec("out_hydrogen_bonds_receptor.dat", ios::app);
  for (size_t i{}; i < pro_in.list_h_bonds.dat.size(); ++i) {
    out_h_bonds_rec
        << setw(4) << left
        << info_rec[pro_in.list_h_bonds.dat[i].get_index_acceptor()].get_atom()
        << setw(5) << right
        << info_rec[pro_in.list_h_bonds.dat[i].get_index_acceptor()].get_resnr()
        << setw(6) << right
        << info_lg[pro_in.list_h_bonds.dat[i].get_index_h()].get_atom()
        << setw(5) << right
        << info_lg[pro_in.list_h_bonds.dat[i].get_index_h()].get_resnr()
        << '\n';
  }
  out_h_bonds_rec.close();
  fstream out_h_bonds_lig("out_hydrogen_bonds_ligand.dat", ios::app);
  for (size_t i{}; i < lig.list_h_bonds.dat.size(); ++i) {
    out_h_bonds_lig
        << setw(4) << left
        << info_lg[lig.list_h_bonds.dat[i].get_index_acceptor()].get_atom()
        << setw(5) << right
        << info_lg[lig.list_h_bonds.dat[i].get_index_acceptor()].get_resnr()
        << setw(6) << right
        << info_rec[lig.list_h_bonds.dat[i].get_index_h()].get_atom() << setw(5)
        << right << info_rec[lig.list_h_bonds.dat[i].get_index_h()].get_resnr()
        << '\n';
  }
  out_h_bonds_lig.close();
  return 0;
}