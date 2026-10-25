// [ moleculetype ], [ atomtypes ] and [ atoms ] readers of a GROMACS
// topology, shared by the protein and ligand paths.
#include "Multi_MOL.h"
#include "arguments.h"
#include "functions.h"
#include "info.h"
#include "miscellaneous.h"
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <string_view>
using std::string_view;
#include <vector>

using namespace std;
using index_type = class index;


int process_moleculetype(fstream &input_topo_protein, string &lig_name,
                         string &protein_1) {
  if (protein_1.find("moleculetype") != string::npos) {
    while (getline(input_topo_protein, protein_1)) {
      if (protein_1.length() > 1 && protein_1.find("[ ") == string::npos) {
        if (protein_1[0] != ';') {
          istringstream is(protein_1);
          int nrexcl{};
          is >> lig_name >> nrexcl;
          return 0;
        } else
          continue;
      } else
        return 1;
    }
    return 0;
  } else
    return 1;
}
int process_atomtypes(fstream &input_topo_protein, string &protein_1,
                      para &parameters, string force_f, string &lig_name) {
  if (protein_1.find("atomtypes") != string::npos) {
    while (getline(input_topo_protein, protein_1)) {
      if (protein_1.length() <= 1 || protein_1.find(";") == 0)
        continue;
      if (protein_1.find(";") != string::npos) {
        protein_1 = protein_1.substr(0, protein_1.find(";"));
      }

      if (protein_1.length() > 10 && protein_1.find("[ ") == string::npos &&
          protein_1.find(" ]") == string::npos) {
        string s_t1{}, s_t2{}, s_t3{}, s_t4{}, s_t5{}, s_t6{}, s_t7{};
        istringstream i_temp(protein_1);
        i_temp >> s_t1 >> s_t2 >> s_t3 >> s_t4 >> s_t5 >> s_t6 >> s_t7;

        if (s_t3.empty() || s_t6.empty() || s_t7.empty())
          continue;

        try {
          double mass{1.0};
          if (stod(s_t3) > 0.1)
            mass = stod(s_t3);
          if (force_f.find("AMBER") != string::npos ||
              force_f.find("CHARMM") != string::npos) {
            if (std::fabs(stod(s_t7)) <= 0.0001 &&
                std::fabs(stod(s_t6)) <= 0.0001)
              parameters.base.emplace_back(stod(s_t6), stod(s_t7), mass, 0.0,
                                           s_t1);
            else
              parameters.base.emplace_back(
                  stod(s_t6), stod(s_t7), mass,
                  5.0 * pow(2.0, 1.0 / 6.0) * stod(s_t6), s_t1);
            parameters.register_atom_type(s_t1);
          } else {
            if (std::fabs(stod(s_t7)) <= 0.0001 &&
                std::fabs(stod(s_t6)) <= 0.0001)
              parameters.base.emplace_back(stod(s_t6), stod(s_t7), mass, 0.0,
                                           s_t1);
            else
              parameters.base.emplace_back(
                  stod(s_t6), stod(s_t7), mass,
                  10 * 0.5 *
                      pow((stod(s_t7) * 2 * pow(stod(s_t6), -1.0)),
                          (1.0 / 6.0)),
                  s_t1);
            parameters.register_atom_type(s_t1);
          }
        } catch (const std::exception &) {
          continue;
        }
      } else if (protein_1.find("moleculetype") != string::npos) {
        process_moleculetype(input_topo_protein, lig_name, protein_1);
        break;
      } else
        break;
    }
    return 0;
  }
  return 1;
}
int Multi_MOL::process_atom_from_topo(fstream &input_topo_protein,
                                      string &protein_1, string &lig_name,
                                      int &i,
                                      vector<topo_info_temp> &info_topo) {
  if (protein_1.find("[ atoms ]") != string::npos) {
    while (getline(input_topo_protein, protein_1)) {
      string line{protein_1};
      const size_t first_non_ws{line.find_first_not_of(" \t\r")};
      if (first_non_ws == string::npos)
        continue;

      line.erase(0, first_non_ws);

      if (line[0] == '[')
        break;

      if (line[0] == ';' || line[0] == '#')
        continue;

      if (const size_t comment_pos = line.find(';'); comment_pos != string::npos)
        line = line.substr(0, comment_pos);

      if (line.find_first_not_of(" \t\r") == string::npos)
        continue;

      // Every well-formed row inside [ atoms ] belongs to the active
      // moleculetype. Matching against lig_name breaks RNA/DNA because those
      // rows contain residue names (GUA/URA/...) instead of the parent name.
      istringstream topo_is(line);
      int i_temp{}, chrg{}, resnr{};
      string a_type{}, a_name{}, res_name{};
      double charge{}, mass{};
      if (!(topo_is >> i_temp >> a_type >> resnr >> res_name >> a_name >>
            chrg >> charge >> mass))
        continue;

      info_topo.emplace_back(a_name, a_type, res_name, charge, mass, resnr);
      ++i;
    }
  } else
    return 1;
  return 0;
}

int Multi_MOL::post_process_atom_from_topo(vector<topo_info_temp> &info_topo) {
  if (info_topo.empty()) {
    std::cerr << "[INFO] Standard topology parser: 0 atoms parsed from "
              << "topo file. If using MARTINI, this is expected.\n";
    return 1;
  }

  if (this->parameters.size() < info.size()) {
    this->parameters.resize(info.size());
  }

  // ASSIGN INTERNAL RESNUM to PDB atoms to allow unique identification
  // This separates "Display ResNr" (which can be duplicate/confusing) from
  // "Calculation ID". We do this purely based on PDB connectivity/changes.
  int internal_cnt = 0;
  if (!info.empty()) {
    info[0].set_internal_resnr(0);
  }
  for (size_t k = 1; k < info.size(); ++k) {
    if (info[k].get_resnr() != info[k - 1].get_resnr() ||
        info[k].get_chain() != info[k - 1].get_chain()) {
      internal_cnt++;
    }
    info[k].set_internal_resnr(internal_cnt);
  }

  // AUTO-SUBSET FOR COMMON TOPOLOGY SYSTEMS
  // If the parsed topology is larger than the extracted PDB group, we use the
  // original PDB atom IDs to pluck the exact corresponding topology atom,
  // bypassing string-matching.
  if (info_topo.size() > info.size() && !info.empty()) {
    const bool ids_were_set =
        std::all_of(info.begin(), info.end(), [](const mol_info &atom_info) {
          return atom_info.has_original_atom_id();
        });

    if (ids_were_set) {
      bool can_subset = true;
      for (size_t i = 0; i < info.size(); ++i) {
        if (info[i].get_original_atom_id() >= info_topo.size()) {
          can_subset = false;
          break;
        }
      }
      if (can_subset) {
        vector<topo_info_temp> subset_topo;
        subset_topo.reserve(info.size());
        for (size_t i = 0; i < info.size(); ++i) {
          subset_topo.push_back(info_topo[info[i].get_original_atom_id()]);
        }
        info_topo =
            std::move(subset_topo); // Now info.size() == info_topo.size()
      }
    }
  }

  // 1-1 MAPPING LOGIC (Common Topo)
  if (info.size() == info_topo.size()) {
    //<< "). Using sequential 1-1 parameter assignment.\n";
    int mismatch_count = 0;
    for (size_t i = 0; i < info.size(); ++i) {
      // Warning for sanity check
      if (info[i].get_res_name() != info_topo[i].res_name ||
          info[i].get_atom() != info_topo[i].atom_name) {
        if (mismatch_count < 10) { // Limit spam
          //<< " PDB: " << info[i].get_res_name() << ":"
          //<< info[i].get_atom() << " Topo: " << info_topo[i].res_name
          //<< ":" << info_topo[i].atom_name << endl;
        }
        mismatch_count++;
      }

      // Assign Parameters sequentially
      info[i].set_atom_type(info_topo[i].atom_type);
      this->parameters[i].set_charge(info_topo[i].atom_charge);
      this->parameters[i].set_mass(info_topo[i].atom_mass);
      // Do NOT overwrite info[i].resnr - keep PDB numbering for output.
    }
    if (mismatch_count > 0) {
      //<< endl;
    }
  } else {
    //<< info_topo.size() << "). Falling back to Strict Match Search.\n";
    // Fallback Search Logic
    int matched_count = 0;
    vector<bool> consumed(info_topo.size(), false);
    for (size_t i{}; i < info.size(); ++i) {
      int best_match_index = -1;
      for (size_t j = 0; j < info_topo.size(); ++j) {
        if (consumed[j])
          continue;

        bool names_match = (info[i].get_res_name() == info_topo[j].res_name);
        // Simple alias check
        if (!names_match) {
          string pdb_name = info[i].get_res_name();
          string topo_name = info_topo[j].res_name;
          if (pdb_name.length() == topo_name.length() + 1 &&
              (pdb_name[0] == 'N' || pdb_name[0] == 'C')) {
            if (pdb_name.substr(1) == topo_name)
              names_match = true;
          }
        }

        if (names_match && info[i].get_atom() == info_topo[j].atom_name) {
          best_match_index = j;
          break;
        }
      }

      if (best_match_index != -1) {
        info[i].set_atom_type(info_topo[best_match_index].atom_type);
        this->parameters[i].set_charge(info_topo[best_match_index].atom_charge);
        this->parameters[i].set_mass(info_topo[best_match_index].atom_mass);
        consumed[best_match_index] = true;
        matched_count++;
      } else {
        //<< ":" << info[i].get_atom() << ") not found in Topo!\n";
        this->parameters[i].set_charge(0.0); // Reset to 0
      }
    }
    //<< " atoms.\n";
  }

  return 0;
}
