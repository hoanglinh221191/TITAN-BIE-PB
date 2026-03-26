#include "Multi_MOL.h"
#include "functions.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <string_view>
using std::string_view;
#include <vector>

int Multi_MOL::topo_protein(std::string_view input_topo,
                            std::string_view force_f, para &parameters,
                            const std::vector<topo_info_temp> &extra_bonds) {
  if (input_topo.length() > 5 &&
      input_topo.substr(input_topo.length() - 4) != ".top" &&
      input_topo.substr(input_topo.length() - 4) != ".itp") {
    cout << "The file type of topo " << input_topo
         << " is not supported. This software supports only .top and .itp\n";
    return 1;
  }
  //<< info.size() << "\n";

  fstream input(std::string(input_topo), ios::in);
  if (!input.is_open()) {
    cout << "Cannot open topo file " << input_topo << "\n";
    return 1;
  } else {
    cout << "Parsing topo from file " << input_topo << "\n";
  }

  vector<string> content_input;
  std::vector<std::pair<std::string, int>> system_molecules;
  string main_content{};

  for (string temp{}; getline(input, temp);) {
    if (temp.find("[ molecules ]") != string::npos) {
      main_content = "molecules";
    } else if (main_content == "molecules" && temp.length() > 2 &&
               temp[0] != ';' && temp[0] != '#') {
      std::istringstream iss(temp);
      std::string m_name;
      int m_count;
      if (iss >> m_name >> m_count) {
        system_molecules.emplace_back(m_name, m_count);
      }
    }

    if (temp.find(";") != string::npos)
      temp = temp.substr(0, temp.find(";"));
    if (temp.length() < 1)
      continue;
    if (temp.find("#include") != string::npos) {
      size_t pos_1{}, pos_2{};
      string temp_2{};
      pos_1 = temp.find("\"");
      if (pos_1 != string::npos)
        temp_2 = temp.substr(pos_1 + 1);
      pos_2 = temp_2.find("\"");
      if (pos_2 != string::npos) {
        temp_2 = temp.substr(pos_1 + 1, pos_2);
        temp_2.erase(remove(temp_2.begin(), temp_2.end(), ' '), temp_2.end());
        content_input.emplace_back(temp_2);
      }
    }
  }

  if (content_input.empty()) {
    cout << "cannot parse " << input_topo << " file!!!\n";
    return 1;
  }

  bool flag_in{false};
  vector<topo_info_temp> info_topo_total_temp;
  vector<topo_info_temp> info_topo;

  for (auto &filename : content_input) {
    // cout << filename << "\n";
    string content{}, molname{};
    fstream input(filename, ios::in);
    vector<string> data;
    for (string temp{}; getline(input, temp);) {
      /*if (temp.find(";") != string::npos)
              temp = temp.substr(0, temp.find(";"));
      if (temp.length() < 1)
              continue;*/
      data.emplace_back(temp);

      if (temp.find("[ moleculetype ]") != string::npos) {
        content = "moleculetype";
      }

      if (content == "moleculetype" && temp.find("[ atoms ]") != string::npos) {
        for (int j = (int)data.size() - 1; j >= 0; --j) {
          if (data[j].find("; name	nrexcl") != string::npos &&
              j + 1 < (int)data.size()) {
            std::istringstream iss(data[j + 1]);
            iss >> molname;
            break;
          }
        }
      }

      if (temp.find("[ atomtypes ]") != string::npos) {
        content = "atomtypes";
      }

      if (temp.find("[ molecules ]") != string::npos) {
        content = "molecules";
      } else if (content == "molecules" && temp.length() > 2 &&
                 temp[0] != ';' && temp[0] != '#') {
        std::istringstream iss(temp);
        std::string m_name;
        int m_count;
        // BUG 4 FIX: Removed 'if (filename == input_topo)' since
        // system_molecules was already parsed correctly in the first pass and
        // this guard was always false.
      }
    }

    input.close();
    input.open(filename, ios::in);

    for (string temp{}; getline(input, temp);) {
      process_atomtypes(input, temp, parameters, std::string(force_f), molname);
    }
    input.close();
    int atomnum{};
    input.open(filename, ios::in);
    // BUG 3 FIX: Removed hardcoded molname = "Protein";
    // process_atom_from_topo needs actual names. We let process_moleculetype
    // set it later. If we passed molname as it was found in the first include
    // pass, we're good.
    vector<topo_info_temp> info_topo_temp;
    for (string temp{}; getline(input, temp);) {
      if (temp.find("moleculetype") != string::npos) {
        process_moleculetype(input, molname, temp);
      }
      this->process_atom_from_topo(input, temp, molname, atomnum,
                                   info_topo_temp);
    }
    if (info_topo_temp.size() == info.size()) {
      //<< info_topo_temp.size() << " atoms)\n";
      info_topo = info_topo_temp;
      flag_in = true;

      // BUG 5 FIX: When the parser bypasses pre_proc_pdb and assigns arrays
      // directly, we MUST assign original_atom_id to each PDB info atom so that
      // subsetting works.
      for (size_t a = 0; a < info.size(); ++a) {
        info[a].set_original_atom_id(a);
      }
    } else {
      //<< ". Parsed: " << info_topo_temp.size()
      //<< ", Expected: " << info.size() << ". Accumulating.\n";
      for (auto &i : info_topo_temp)
        info_topo_total_temp.emplace_back(i);
    }

    // Check if this molecule is in system_molecules and requires multiple
    // copies
    int num_copies = 1;
    for (const auto &smol : system_molecules) {
      if (smol.first == molname) {
        num_copies = smol.second;
        break;
      }
    }

    if (num_copies > 1 && !info_topo_temp.empty()) { // Append copies if defined
      for (int c = 1; c < num_copies; ++c) {
        for (auto &i : info_topo_temp) {
          info_topo_total_temp.emplace_back(i);
        }
      }
    }
  }

  if (flag_in == false)
    info_topo = info_topo_total_temp;

  // BUG 5 FIX: If we have multiple copies added to info_topo_total_temp,
  // we must also initialize original_atom_ids if the sizes match exactly
  // to avoid subset array bounds crashes.
  if (info_topo.size() == info.size()) {
    for (size_t a = 0; a < info.size(); ++a) {
      info[a].set_original_atom_id(a);
    }
  }

  this->post_process_atom_from_topo(info_topo);
  // cout << "info topo size: " << info_topo.size() << '\n';
  // cout << "parameter size: " << parameters.base.size() << '\n';
  return 0;
}
