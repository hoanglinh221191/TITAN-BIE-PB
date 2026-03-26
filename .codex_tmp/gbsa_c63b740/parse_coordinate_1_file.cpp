#include "MOL.h"
#include "functions.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

using namespace std;

int helper_parse_coor_not_copy(Multi_MOL &pro_in, Multi_MOL &lig,
                               Multi_MOL &water_traj, groups_cal &groups,
                               arguments &Arguments, box &Box_temp,
                               bool &flag_set_time, double &time_of_snap_temp,
                               fstream &input_pdb, size_t &I,
                               bool &flag_set_box, ::info &output_info,
                               string protein_1, size_t &atom_num,
                               string input_name) {

  if (flag_set_time == false) {
    time_of_snap_temp = static_cast<double>(I);
    int check_time{
        check_and_settime(time_of_snap_temp, protein_1, input_pdb, Arguments)};
    if (check_time == 0)
      flag_set_time = true;
  }

  if (protein_1.length() >= 55 && protein_1.substr(0, 6) == "CRYST1") {
    Box_temp =
        box(stod(protein_1.substr(6, 10)), stod(protein_1.substr(16, 9)),
            stod(protein_1.substr(25, 9)), stod(protein_1.substr(34, 7)),
            stod(protein_1.substr(41, 7)), stod(protein_1.substr(48, 7)));
    if (flag_set_box == false && I != 0)
      flag_set_box = true;
  }
  if (protein_1.length() >= 55 && protein_1.substr(0, 4) == "ATOM") {
    ++atom_num;
    if (Arguments.flag_end_set && time_of_snap_temp > Arguments.end_time)
      return 1;
    if (Arguments.flag_begin_set == false && I == 0)
      output_info.set_current_begin_time(time_of_snap_temp);
    if (check_in_time(Arguments, time_of_snap_temp) == true) {
      string atom{protein_1.substr(12, 4)};
      atom.erase(remove(atom.begin(), atom.end(), ' '), atom.end());
      string resname{protein_1.substr(17, 4)};
      resname.erase(remove(resname.begin(), resname.end(), ' '), resname.end());
      if (!groups.empty()) {
        if (groups.get_group(groups.group_1_index())->check_in(atom_num))
          pro_in.emplace(I, Box_temp, time_of_snap_temp, Arguments.force_f,
                         mol_info(stoull(protein_1.substr(22, 4)),
                                  string("NULL"), atom, protein_1.substr(21, 1),
                                  resname),
                         [&] {
                           MOL m;
                           m.pos.x.emplace_back(stod(protein_1.substr(31, 8)));
                           m.pos.y.emplace_back(stod(protein_1.substr(39, 8)));
                           m.pos.z.emplace_back(stod(protein_1.substr(47, 8)));
                           return m;
                         }());
        if (groups.get_group(groups.group_2_index())->check_in(atom_num))
          lig.emplace(I, Box_temp, time_of_snap_temp, Arguments.force_f,
                      mol_info(stoull(protein_1.substr(22, 4)), string("NULL"),
                               atom, protein_1.substr(21, 1), resname),
                      [&] {
                        MOL m;
                        m.pos.x.emplace_back(stod(protein_1.substr(31, 8)));
                        m.pos.y.emplace_back(stod(protein_1.substr(39, 8)));
                        m.pos.z.emplace_back(stod(protein_1.substr(47, 8)));
                        return m;
                      }());
        if (Arguments.n_wat >= 0 && Arguments.n_wat < groups.group_num() &&
            groups.get_group(Arguments.n_wat)->check_in(atom_num))
          water_traj.emplace(
              I, Box_temp, time_of_snap_temp, Arguments.force_f,
              mol_info(stoull(protein_1.substr(22, 4)), string("NULL"), atom,
                       protein_1.substr(21, 1), resname),
              [&] {
                MOL m;
                m.pos.x.emplace_back(stod(protein_1.substr(31, 8)));
                m.pos.y.emplace_back(stod(protein_1.substr(39, 8)));
                m.pos.z.emplace_back(stod(protein_1.substr(47, 8)));
                return m;
              }(),
              "water");
      } else {
        if (atom_num < Arguments.begin_mol_2)
          pro_in.emplace(I, Box_temp, time_of_snap_temp, Arguments.force_f,
                         mol_info(stoull(protein_1.substr(22, 4)),
                                  string("NULL"), atom, protein_1.substr(21, 1),
                                  resname),
                         [&] {
                           MOL m;
                           m.pos.x.emplace_back(stod(protein_1.substr(31, 8)));
                           m.pos.y.emplace_back(stod(protein_1.substr(39, 8)));
                           m.pos.z.emplace_back(stod(protein_1.substr(47, 8)));
                           return m;
                         }());
        else
          lig.emplace(I, Box_temp, time_of_snap_temp, Arguments.force_f,
                      mol_info(stoull(protein_1.substr(22, 4)), string("NULL"),
                               atom, protein_1.substr(21, 1), resname),
                      [&] {
                        MOL m;
                        m.pos.x.emplace_back(stod(protein_1.substr(31, 8)));
                        m.pos.y.emplace_back(stod(protein_1.substr(39, 8)));
                        m.pos.z.emplace_back(stod(protein_1.substr(47, 8)));
                        return m;
                      }());
      }
    } else {
      flag_set_time = false;
      flag_set_box = false;
    }
  }
  if ((protein_1 == "ENDMDL" || input_pdb.eof() || protein_1 == "END") &&
      check_in_time(Arguments, time_of_snap_temp) == true) {
    if (pro_in.m_mol.empty() == true || lig.m_mol.empty() == true) {
      cout << "Error occured in " << input_name
           << " file. Cannot parse atom coordinates of MOL(s)\n";
      cout << pro_in.m_mol.size() << " " << lig.m_mol.size();
      EXIT_FAILURE;
    }

    pro_in.check_and_pop_size_2last_snaps(I, " of molecule 1 ");
    lig.check_and_pop_size_2last_snaps(I, " of molecule 2 ");

    if (flag_set_box == false && I != 0 && I < pro_in.m_mol.size() &&
        I < lig.m_mol.size()) {
      pro_in.m_mol[I].Box = pro_in.m_mol[I - 1].Box;
      lig.m_mol[I].Box = lig.m_mol[I - 1].Box;
    } else
      flag_set_box = false;

    ++I;
    atom_num = 0;
    flag_set_time = false;
  }

  return 0;
}

void take_coordinate(Multi_MOL &pro_in, Multi_MOL &lig, Multi_MOL &water_traj,
                     groups_cal &groups, string input_name, string force_f,
                     arguments &Arguments, const amino_acid &amine,
                     ::info &output_info, int flag_raw, size_t max_frames) {
  box Box_temp;
  double time_of_snap_temp{};
  bool flag_set_time{false}, flag_set_box{false};
  fstream input_pdb;
  input_pdb.open(input_name, ios::in);
  if (!input_pdb.is_open()) {
    cout << "Can not open " << input_name << " file!" << '\n';
    exit(EXIT_FAILURE);
  }

  size_t I{}, atom_num{};
  bool flag_end{true};
  while (flag_end == true) {
    if (max_frames > 0 && pro_in.m_mol.size() >= max_frames) {
      flag_end = false;
      break;
    }
    if (max_frames > 0 && lig.m_mol.size() >= max_frames) {
      flag_end = false;
      break;
    }

    vector<string> buffer_string;
    size_t count_snap{};
    for (string protein_1{}; getline(input_pdb, protein_1);) {
      buffer_string.emplace_back(protein_1);
      if (protein_1 == "ENDMDL" || protein_1 == "END") {
        ++count_snap;
        if (count_snap == 10) {
          for (size_t i{}; i < buffer_string.size(); ++i)
            if (1 == helper_parse_coor_not_copy(
                         pro_in, lig, water_traj, groups, Arguments, Box_temp,
                         flag_set_time, time_of_snap_temp, input_pdb, I,
                         flag_set_box, output_info, buffer_string[i], atom_num,
                         input_name)) {
              flag_end = false;
              break;
            }
          buffer_string.clear();
          count_snap = 0;
        }
      }
      if (input_pdb.eof()) {
        flag_end = false;
        break;
      }
      if (flag_end == false)
        break;
    }
    if (!buffer_string.empty()) {
      if (buffer_string[buffer_string.size() - 1].length() >= 1 &&
          buffer_string[buffer_string.size() - 1] != "ENDMDL" &&
          buffer_string[buffer_string.size() - 1] != "END")
        buffer_string.emplace_back("ENDMDL");
      for (size_t i{}; i < buffer_string.size(); ++i)
        if (1 == helper_parse_coor_not_copy(
                     pro_in, lig, water_traj, groups, Arguments, Box_temp,
                     flag_set_time, time_of_snap_temp, input_pdb, I,
                     flag_set_box, output_info, buffer_string[i], atom_num,
                     input_name)) {
          flag_end = false;
          buffer_string.clear();
          break;
        }
    }
    flag_end = false;
  }
  input_pdb.close();

  if (pro_in.m_mol.size() >= 2)
    output_info.set_current_dt(
        pro_in.m_mol[pro_in.m_mol.size() - 1].time_of_snap -
        pro_in.m_mol[pro_in.m_mol.size() - 2].time_of_snap);

  if (Arguments.flag_end_set == false && !pro_in.m_mol.empty())
    output_info.set_current_end_time(
        pro_in.m_mol[pro_in.m_mol.size() - 1].time_of_snap);

  if (Arguments.flag_begin_set == true && pro_in.m_mol.empty() &&
      time_of_snap_temp < Arguments.begin_time) {
    cout << "WARNING!!! Your begin time is " << Arguments.begin_time
         << " ps while time of the last snapshot is smaller than that!\n";
  }

  pro_in.gen_vector_parameters();
  lig.gen_vector_parameters();
  if (Arguments.n_wat >= 0)
    water_traj.gen_vector_parameters();
}

int check_and_settime(double &time_of_snap, string protein_1,
                      fstream &input_pdb, arguments &Arguments) {
  if (protein_1.length() >= 5 && protein_1.find("TITLE") != string::npos) {
    size_t found = protein_1.find("=");
    size_t found_2 = protein_1.find("step"),
           time_length{protein_1.length() < found_2 ? protein_1.length()
                                                    : found_2};

    if (found != string::npos && protein_1.length() > found + 2) {
      time_of_snap =
          (stod(protein_1.substr(found + 1, time_length - found - 1)));
      return 0;
    } else {
      getline(input_pdb, protein_1);
      found = protein_1.find("t=");
      found_2 = protein_1.find("step");
      time_length = protein_1.length() < found_2 ? protein_1.length() : found_2;
      if (found != string::npos) {
        time_of_snap =
            (stod(protein_1.substr(found + 2, time_length - found - 2)));
      }
      return 0;
    }
  } else
    return 1;
}