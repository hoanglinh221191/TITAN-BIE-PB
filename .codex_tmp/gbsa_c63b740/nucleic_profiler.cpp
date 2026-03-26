#include "nucleic_profiler.h"
#include <fstream>
#include <iostream>

void NucleicProfiler::reset() {
  interaction_timeline.clear();
  global_recognition_profile.clear();
  resnum_to_name.clear();
  total_frames_processed = 0;
}

std::vector<double> NucleicProfiler::calc_residue_center(
    const MOL &frame, const std::vector<mol_info> &info,
    const class index &chain_idx, size_t res_idx_in_chain,
    const std::string &center_type) {
  std::vector<double> center = {0.0, 0.0, 0.0};
  int count = 0;

  size_t start_idx = (res_idx_in_chain == 0)
                         ? chain_idx.get_res_pos_e(0) - 1
                         : chain_idx.get_res_pos_e(res_idx_in_chain - 1);
  size_t end_idx = chain_idx.get_res_pos_e(res_idx_in_chain) - 1;

  for (size_t i = start_idx; i <= end_idx && i < info.size(); ++i) {
    std::string atom_name = info[i].get_atom();
    bool match = false;

    if (center_type == "AromaticRing") {
      if (atom_name == "CG" || atom_name == "CD1" || atom_name == "CD2" ||
          atom_name == "CE1" || atom_name == "CE2" || atom_name == "CE3" ||
          atom_name == "CZ" || atom_name == "CZ2" || atom_name == "CZ3" ||
          atom_name == "CH2") {
        match = true;
      }
    } else if (center_type == "CationicGroup") {
      if (atom_name == "NH1" || atom_name == "NH2" || atom_name == "NZ" ||
          atom_name == "ND1" || atom_name == "NE2") {
        match = true;
      }
    } else if (center_type == "Phosphate") {
      if (atom_name == "P" || atom_name == "OP1" || atom_name == "OP2" ||
          atom_name == "O1P" || atom_name == "O2P") {
        match = true;
      }
    } else if (center_type == "Nucleobase") {
      if (atom_name == "N1" || atom_name == "C2" || atom_name == "N3" ||
          atom_name == "C4" || atom_name == "C5" || atom_name == "C6" ||
          atom_name == "N7" || atom_name == "C8" || atom_name == "N9") {
        match = true;
      }
    }

    if (match) {
      center[0] += frame.pos.x[i];
      center[1] += frame.pos.y[i];
      center[2] += frame.pos.z[i];
      count++;
    }
  }

  if (count > 0) {
    center[0] /= count;
    center[1] /= count;
    center[2] /= count;
  }
  return center;
}

void NucleicProfiler::process_frame(const MOL &pro_frame, const MOL &na_frame,
                                    const std::vector<mol_info> &pro_info,
                                    const std::vector<mol_info> &na_info,
                                    const std::vector<class index> &pro_index,
                                    const std::vector<class index> &na_index,
                                    size_t num_snap, double time_ps) {

  total_frames_processed++;

  // Step 1: Gather all protein atoms
  struct AtomInfo {
    int resnr;
    std::string resname;
    double x, y, z;
    AtomInfo(int r, std::string n, double x, double y, double z)
        : resnr(r), resname(std::move(n)), x(x), y(y), z(z) {}
  };
  std::vector<AtomInfo> pro_atoms;

  for (const auto &chn : pro_index) {
    if (const_cast<class index &>(chn).get_type() != "protein")
      continue;
    for (size_t r = 1; r < const_cast<class index &>(chn).get_size_v(); ++r) {
      size_t start_idx =
          (r == 1) ? chn.get_res_pos_e(0) - 1 : chn.get_res_pos_e(r - 1);
      size_t end_idx = chn.get_res_pos_e(r) - 1;

      for (size_t i = start_idx;
           i <= end_idx && i < pro_info.size() && i < pro_frame.pos.x.size(); ++i) {
        pro_atoms.emplace_back(static_cast<int>(pro_info[i].get_resnr()),
                               pro_info[i].get_res_name(), pro_frame.pos.x[i],
                               pro_frame.pos.y[i], pro_frame.pos.z[i]);
      }
    }
  }

  // To prevent double counting in the same frame
  std::set<int> base_contact_resnrs;
  std::set<int> backbone_contact_resnrs;

  // Step 2: Iterate over Nucleic Acid Atoms and checking squared distance
  // < 12.25
  for (const auto &chn : na_index) {
    std::string type = const_cast<class index &>(chn).get_type();
    if (type != "rna" && type != "dna" && type != "complex")
      continue;

    for (size_t r = 1; r < const_cast<class index &>(chn).get_size_v(); ++r) {
      size_t start_idx =
          (r == 1) ? chn.get_res_pos_e(0) - 1 : chn.get_res_pos_e(r - 1);
      size_t end_idx = chn.get_res_pos_e(r) - 1;

      for (size_t i = start_idx;
           i <= end_idx && i < na_info.size() && i < na_frame.pos.x.size(); ++i) {
        double na_x = na_frame.pos.x[i];
        double na_y = na_frame.pos.y[i];
        double na_z = na_frame.pos.z[i];
        std::string atom_name = na_info[i].get_atom();

        bool is_bb = is_backbone(atom_name);
        bool is_b = is_base(atom_name);
        if (!is_bb && !is_b)
          continue;

        // Compare against protein atoms
        for (const auto &p_atom : pro_atoms) {
          double dx = na_x - p_atom.x;
          double dy = na_y - p_atom.y;
          double dz = na_z - p_atom.z;
          double dist_sq = dx * dx + dy * dy + dz * dz;

          if (dist_sq < 12.25) { // 3.5^2
            if (is_bb && backbone_contact_resnrs.find(p_atom.resnr) ==
                             backbone_contact_resnrs.end()) {
              backbone_contact_resnrs.insert(p_atom.resnr);
              resnum_to_name.insert_or_assign(p_atom.resnr, p_atom.resname);
            }
            if (is_b && base_contact_resnrs.find(p_atom.resnr) ==
                            base_contact_resnrs.end()) {
              base_contact_resnrs.insert(p_atom.resnr);
              resnum_to_name.insert_or_assign(p_atom.resnr, p_atom.resname);
            }
          }
        }
      }
    }
  }

  // Step 3: Global update
  for (int resnr : base_contact_resnrs) {
    global_recognition_profile[resnr].first++;
  }
  for (int resnr : backbone_contact_resnrs) {
    global_recognition_profile[resnr].second++;
  }
}

void NucleicProfiler::export_results(const std::string &csv_filename) {
  if (total_frames_processed == 0)
    return;

  std::ofstream out(csv_filename);
  if (!out)
    return;

  out << "Protein_ResNum,ResName,Base_Contact_Frames,Backbone_Contact_Frames,"
         "Base_Occupancy_%,Backbone_Occupancy_%,Recognition_Type\n";

  for (const auto &[resnr, counts] : global_recognition_profile) {
    int base_contacts = counts.first;
    int backbone_contacts = counts.second;

    double base_occ =
        (static_cast<double>(base_contacts) / total_frames_processed) * 100.0;
    double back_occ =
        (static_cast<double>(backbone_contacts) / total_frames_processed) *
        100.0;

    if (base_occ > 5.0 || back_occ > 5.0) {
      std::string type = "Mixed_Binding";
      if (base_contacts > backbone_contacts * 1.5) {
        type = "Highly_Specific";
      } else if (backbone_contacts > base_contacts * 1.5) {
        type = "Non_Specific_Backbone";
      }

      out << resnr << "," << resnum_to_name[resnr] << "," << base_contacts
          << "," << backbone_contacts << "," << base_occ << "," << back_occ
          << "," << type << "\n";
    }
  }
}
