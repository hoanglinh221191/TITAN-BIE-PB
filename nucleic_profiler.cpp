#include "nucleic_profiler.h"
#include "polymer_classification.h"
#include <fstream>
#include <iostream>

namespace PC = PolymerClassification;

namespace {

std::pair<size_t, size_t>
residue_atom_range(const std::vector<class index> &indices, size_t chain_pos,
                   size_t res_idx, size_t total_atoms) {
  if (chain_pos >= indices.size())
    return {0, 0};
  const auto &chn = indices[chain_pos];
  if (res_idx >= chn.get_size_v())
    return {0, 0};

  const size_t start_idx = chn.get_res_pos_e(res_idx) - 1;
  size_t end_idx = total_atoms;
  if (res_idx + 1 < chn.get_size_v()) {
    end_idx = chn.get_res_pos_e(res_idx + 1) - 1;
  } else if (chain_pos + 1 < indices.size()) {
    end_idx = indices[chain_pos + 1].get_res_pos_e(0) - 1;
  }
  return {start_idx, std::min(end_idx, total_atoms)};
}

bool atom_has_nucleic_signature(const std::string &atom_name) {
  const std::string atom = PC::normalize_nucleic_atom_name(atom_name);
  return atom == "P" || atom == "OP1" || atom == "OP2" || atom == "OP3" ||
         atom == "O1P" || atom == "O2P" || atom == "O3P" ||
         atom == "C1'" || atom == "C2'" || atom == "O2'" || atom == "C3'" ||
         atom == "O3'" || atom == "C4'" || atom == "O4'" || atom == "C5'" ||
         atom == "O5'" || atom == "N9" || atom == "C8" || atom == "N7" ||
         atom == "O6" || atom == "N6" || atom == "N2" || atom == "O2" ||
         atom == "O4" || atom == "N4";
}

bool residue_is_protein_candidate(const class index &chn, size_t res_idx) {
  const auto kind = PC::classify_residue_polymer(chn.get_res_name_e(res_idx));
  if (kind == PC::PolymerKind::Protein)
    return true;

  const std::string type = PC::uppercase_copy(chn.get_type());
  return kind == PC::PolymerKind::Other &&
         type.find("PROTEIN") != std::string::npos;
}

bool residue_is_nucleic_candidate(const std::vector<class index> &indices,
                                  size_t chain_pos, size_t res_idx,
                                  const std::vector<mol_info> &info) {
  if (chain_pos >= indices.size())
    return false;
  const auto &chn = indices[chain_pos];
  const auto kind = PC::classify_residue_polymer(chn.get_res_name_e(res_idx));
  if (kind == PC::PolymerKind::Rna || kind == PC::PolymerKind::Dna)
    return true;
  if (kind == PC::PolymerKind::Protein)
    return false;

  const auto [start_idx, end_idx] =
      residue_atom_range(indices, chain_pos, res_idx, info.size());
  if (start_idx >= end_idx || start_idx >= info.size())
    return false;

  const std::string chain_type = PC::uppercase_copy(chn.get_type());
  const bool chain_hints_nucleic =
      PC::is_nucleic_type(chain_type) || chain_type == "COMPLEX";
  for (size_t i = start_idx; i < end_idx && i < info.size(); ++i) {
    if (atom_has_nucleic_signature(info[i].get_atom()))
      return true;
  }
  return chain_hints_nucleic;
}

} // namespace

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

  for (size_t chain_pos = 0; chain_pos < pro_index.size(); ++chain_pos) {
    const auto &chn = pro_index[chain_pos];
    for (size_t r = 0; r < chn.get_size_v(); ++r) {
      if (!residue_is_protein_candidate(chn, r))
        continue;

      const auto [start_idx, end_idx] =
          residue_atom_range(pro_index, chain_pos, r, pro_info.size());
      for (size_t i = start_idx;
           i < end_idx && i < pro_info.size() && i < pro_frame.pos.x.size();
           ++i) {
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
  for (size_t chain_pos = 0; chain_pos < na_index.size(); ++chain_pos) {
    const auto &chn = na_index[chain_pos];
    for (size_t r = 0; r < chn.get_size_v(); ++r) {
      if (!residue_is_nucleic_candidate(na_index, chain_pos, r, na_info))
        continue;

      const auto [start_idx, end_idx] =
          residue_atom_range(na_index, chain_pos, r, na_info.size());
      for (size_t i = start_idx;
           i < end_idx && i < na_info.size() && i < na_frame.pos.x.size(); ++i) {
        double na_x = na_frame.pos.x[i];
        double na_y = na_frame.pos.y[i];
        double na_z = na_frame.pos.z[i];
        std::string atom_name =
            PC::normalize_nucleic_atom_name(na_info[i].get_atom());

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
