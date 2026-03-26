#include "Multi_MOL.h"
#include "titan_ui.h"
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;

struct LJParam {
  double sigma;
  double epsilon;
};
struct AtomDef {
  std::string type;
  int resnr;
  std::string resname;
  std::string atomname;
  double charge;
};
struct MolTemplate {
  std::vector<AtomDef> atoms;
};

class MartiniTopology {
public:
  std::vector<std::pair<std::string, int>> system_molecules;

private:
  std::unordered_map<std::string, LJParam> nonbond_params;
  std::unordered_map<std::string, MolTemplate> molecule_types;
  std::vector<AtomDef> flat_system;

  std::string get_directive(const std::string &line) {
    size_t p1 = line.find('[');
    size_t p2 = line.find(']');
    if (p1 != std::string::npos && p2 != std::string::npos && p2 > p1) {
      std::string dir = line.substr(p1 + 1, p2 - p1 - 1);
      size_t s = dir.find_first_not_of(" \t");
      size_t e = dir.find_last_not_of(" \t");
      if (s != std::string::npos)
        return dir.substr(s, e - s + 1);
    }
    return "";
  }

public:
  int files_read = 0;

  void parse_files(const std::string &top_file, std::fstream &log_file) {
    std::string base_dir = "";
    size_t last_slash = top_file.find_last_of("/\\");
    if (last_slash != std::string::npos) {
      base_dir = top_file.substr(0, last_slash + 1);
    }

    std::vector<std::string> lines;

    std::function<void(const std::string &)> read_file =
        [&](const std::string &filename) {
          std::fstream f(filename, std::ios::in);
          if (!f.is_open()) {
            TitanUI::print_info(
                "\033[1;33m[WARNING]\033[0m MARTINI Parser: Cannot open "
                "included file: " +
                filename + " — check if file exists.");
            TitanUI::log_message(
                log_file, "[WARNING] Cannot open included file: " + filename);
            return;
          }
          files_read++;
          std::string line;
          while (std::getline(f, line)) {
            size_t comment_pos = line.find(';');
            if (comment_pos != std::string::npos)
              line = line.substr(0, comment_pos); // Xóa comment

            size_t start = line.find_first_not_of(" \t\r\n");
            if (start == std::string::npos)
              continue; // Bỏ qua dòng rỗng
            size_t end = line.find_last_not_of(" \t\r\n");
            line = line.substr(start, end - start + 1);

            // Phân giải #include
            if (line.find("#include") != std::string::npos) {
              size_t q1 = line.find('"');
              size_t q2 = line.find('"', q1 + 1);
              if (q1 != std::string::npos && q2 != std::string::npos) {
                std::string inc_file = line.substr(q1 + 1, q2 - q1 - 1);
                read_file(base_dir + inc_file);
              }
              continue;
            }
            lines.emplace_back(line);
          }
        };

    read_file(top_file);

    std::string current_dir = "";
    std::string current_mol = "";

    for (const auto &line : lines) {
      std::string dir = get_directive(line);
      if (!dir.empty()) {
        current_dir = dir;
        if (dir == "moleculetype")
          current_mol = "";
        continue;
      }

      std::istringstream iss(line);
      std::vector<std::string> tokens;
      std::string token;
      while (iss >> token)
        tokens.emplace_back(token);
      if (tokens.empty())
        continue;

      if (current_dir == "nonbond_params") {
        if (tokens.size() >= 5) {
          if (tokens[0] == tokens[1]) {
            try {
              double sigma = std::stod(tokens[3]);
              double epsilon = std::stod(tokens[4]);
              nonbond_params[tokens[0]] = {sigma, epsilon};
            } catch (...) {
            }
          }
        }
      } else if (current_dir == "atomtypes") {
        if (tokens.size() >= 6) {
          try {
            double sigma = std::stod(tokens[tokens.size() - 2]);
            double epsilon = std::stod(tokens[tokens.size() - 1]);
            nonbond_params[tokens[0]] = {sigma, epsilon};
          } catch (...) {
          }
        }
      } else if (current_dir == "moleculetype") {
        if (current_mol.empty())
          current_mol = tokens[0];
      } else if (current_dir == "atoms") {
        if (!current_mol.empty() && tokens.size() >= 7) {
          try {
            std::string type = tokens[1];
            int resnr = std::stoi(tokens[2]);
            std::string resname = tokens[3];
            std::string atomname = tokens[4];
            double charge = std::stod(tokens[6]);
            molecule_types[current_mol].atoms.emplace_back(type, resnr, resname,
                                                           atomname, charge);
          } catch (...) {
          }
        }
      } else if (current_dir == "molecules") {
        if (tokens.size() >= 2) {
          try {
            system_molecules.emplace_back(tokens[0], std::stoi(tokens[1]));
          } catch (...) {
          }
        }
      }
    }

    for (const auto &comp : system_molecules) {
      if (molecule_types.find(comp.first) == molecule_types.end()) {
        TitanUI::log_message(
            log_file,
            "[ERROR FATAL] Molecule name '" + comp.first +
                "' in [ molecules ] DOES NOT MATCH any [ moleculetype ]!");
      }
    }
  }

  void build_flat_system() {
    for (const auto &comp : system_molecules) {
      std::string m_name = comp.first;
      int count = comp.second;
      if (molecule_types.count(m_name)) {
        for (int c = 0; c < count; ++c) {
          for (const auto &atom : molecule_types[m_name].atoms) {
            flat_system.emplace_back(atom);
          }
        }
      }
    }
  }

  void map_to_titan(std::vector<mol_info> &info,
                    Multi_MOL_parameters &parameters, std::fstream &log_file) {
    if (info.size() != flat_system.size()) {
      TitanUI::print_info("\033[1;33m[WARNING]\033[0m MARTINI Parser: Topology "
                          "flattened atoms (" +
                          std::to_string(flat_system.size()) +
                          ") differs from PDB atoms (" +
                          std::to_string(info.size()) +
                          "). Engaging robust Name-Match fallback mapping.");
    }

    int successfully_mapped = 0;

    std::unordered_map<std::string, std::pair<std::string, double>>
        fallback_dict;
    for (const auto &atom : flat_system) {
      std::string key = atom.resname + ":" + atom.atomname;
      fallback_dict[key] = {atom.type, atom.charge};
    }

    for (size_t i = 0; i < info.size(); ++i) {
      std::string atype = "";
      double acharge = 0.0;

      if (i < flat_system.size() && info.size() == flat_system.size()) {
        atype = flat_system[i].type;
        acharge = flat_system[i].charge;
      } else {
        std::string key = info[i].get_res_name() + ":" + info[i].get_atom();
        if (fallback_dict.count(key)) {
          atype = fallback_dict[key].first;
          acharge = fallback_dict[key].second;
        } else {
          continue; // Unmapped
        }
      }

      parameters[i].set_charge(acharge);
      info[i].set_atom_type(atype);

      auto it = nonbond_params.find(atype);
      if (it != nonbond_params.end()) {
        double sigma = it->second.sigma;
        double epsilon = it->second.epsilon;
        double radius = sigma * 1.122462048309372981 / 2.0; // r_min/2 approx
        parameters[i].set_epsilon_sigma_radius(epsilon, sigma, radius);
        successfully_mapped++;
      }

      if (info[i].get_atom().find("BB") != std::string::npos ||
          info[i].get_atom().find("PO4") != std::string::npos) {
        info[i].set_side_chain(false);
      } else {
        info[i].set_side_chain(true);
      }
    }

    TitanUI::log_kv(log_file, "VdW Types Loaded",
                    std::to_string(nonbond_params.size()));
    TitanUI::log_kv(log_file, "Templates Loaded",
                    std::to_string(molecule_types.size()));
    TitanUI::log_kv(log_file, "Flat Atoms Expanded",
                    std::to_string(flat_system.size()));
    TitanUI::log_kv(log_file, "Beads Mapped",
                    std::to_string(successfully_mapped) + " / " +
                        std::to_string(info.size()));

    // ---- Assign chainnum per bead ----------------------------------------
    // Walk system_molecules in order; for each copy of each molecule type,
    // stamp chainnum on the corresponding mol_info slice.
    // This allows downstream modules (e.g. LLPS) to group beads by chain.
    {
      size_t bead_offset = 0;
      size_t chain_id = 0;
      for (const auto &[mol_name, mol_count] : system_molecules) {
        auto it = molecule_types.find(mol_name);
        if (it == molecule_types.end()) {
          // Unknown molecule type (e.g. solvent W, ions NA/CL):
          // these are NOT in Titan's filtered PDB, so do NOT advance
          // bead_offset. Just skip entirely.
          TitanUI::log_message(log_file, "[INFO] Skipping non-macromolecule '" +
                                             mol_name + "' (" +
                                             std::to_string(mol_count) +
                                             " copies) for chain assignment.");
          continue;
        }
        const size_t beads_per_mol = it->second.atoms.size();
        for (int c = 0; c < mol_count; ++c) {
          for (size_t b = 0; b < beads_per_mol; ++b) {
            size_t idx = bead_offset + b;
            if (idx < info.size())
              info[idx].set_chainnum(chain_id);
          }
          bead_offset += beads_per_mol;
          ++chain_id;
        }
      }
      TitanUI::log_kv(log_file, "LLPS Chains Detected",
                      std::to_string(chain_id));
    }

    TitanUI::print_info("MARTINI Parser: Extracted " +
                        std::to_string(nonbond_params.size()) +
                        " VdW rules, mapped " +
                        std::to_string(successfully_mapped) + " atoms.");
  }

  size_t get_flat_size() const { return flat_system.size(); }

  int get_molecule_count(const std::string &target_mol) const {
    for (const auto &comp : system_molecules) {
      if (comp.first == target_mol) {
        return comp.second;
      }
    }
    return 1;
  }
};

int Multi_MOL::parse_martini_topology(const arguments &Arguments,
                                      std::fstream &log_file) {
  if (info.empty())
    return 0;

  TitanUI::print_info(
      "Initializing Object-Oriented MARTINI Topology Parser...");
  TitanUI::log_section(log_file, "MARTINI Topology Parser");

  MartiniTopology parser;

  // Step 2: Parse Topology Structure
  parser.parse_files(Arguments.topo_name, log_file);
  TitanUI::log_kv(log_file, "Files Parsed", std::to_string(parser.files_read));

  // Step 3: Flatten Templates to Linear Sequence
  parser.build_flat_system();

  // Step 4: Map Properties and Assign Respos/Chains
  parser.map_to_titan(info, parameters, log_file);

  // Auto-detect the primary molecule copies (e.g., molecule_0 or generic
  // Protein) Assuming the primary molecule type is the first one in the [
  // molecules ] section if not explicit
  if (!parser.system_molecules.empty()) {
    int auto_copies = parser.system_molecules.front().second;
    TitanUI::print_info("MARTINI Parser: Auto-detected " +
                        std::to_string(auto_copies) +
                        " copies of primary molecule '" +
                        parser.system_molecules.front().first + "'.");

    // We can't modify const arguments, but we can store this if needed
    // elsewhere, or we handle the copies generation right here based on
    // auto_copies if the design requires it. For now, Titan relies on
    // `Arguments.num_martini_copies` in `generate_martini_copies`. We will need
    // to update the caller or remove the const constraint if we must overwrite
    // it.
  }

  // Re-use `generate_index()` to build the system index using the PDB's native
  // chain/resnr
  // this->generate_index(); // REMOVED because PDB parser already does it!
  return 0;
}
