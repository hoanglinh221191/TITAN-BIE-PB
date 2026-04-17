#ifndef NUCLEIC_PROFILER_H
#define NUCLEIC_PROFILER_H

#include "MOL.h"
#include "MOL_info.h"
#include "arguments.h"
#include "index_v2.h"
#include <map>
#include <set>
#include <string>
#include <vector>


struct NA_Interaction {
  int frame;
  double time_ps;
  std::string type;
  std::string pro_res;
  std::string na_res;
  double distance;
};

class NucleicProfiler {
private:
  std::vector<NA_Interaction> interaction_timeline;

  std::map<int, std::pair<int, int>>
      global_recognition_profile; // Key: Protein ResNum, Value: <Base_Contacts,
                                  // Backbone_Contacts>
  std::map<int, std::string> resnum_to_name;
  size_t total_frames_processed = 0;

  std::vector<double> calc_residue_center(const MOL &frame,
                                          const std::vector<mol_info> &info,
                                          const class index &chain_idx,
                                          size_t res_idx_in_chain,
                                          const std::string &center_type);

  inline bool is_backbone(const std::string &atom_name) const {
    return (atom_name == "P" || atom_name == "OP1" || atom_name == "OP2" ||
            atom_name == "OP3" || atom_name == "O1P" || atom_name == "O2P" ||
            atom_name == "O3P" ||
            atom_name == "O5'" || atom_name == "C5'" || atom_name == "C4'" ||
            atom_name == "O4'" || atom_name == "C3'" || atom_name == "O3'" ||
            atom_name == "C2'" || atom_name == "C1'" || atom_name == "O2'");
  }

  inline bool is_base(const std::string &atom_name) const {
    if (atom_name.empty())
      return false;
    char first = atom_name[0];
    if (first != 'N' && first != 'C' && first != 'O')
      return false;
    if (is_backbone(atom_name))
      return false;
    return true;
  }

public:
  void reset();
  void process_frame(const MOL &pro_frame, const MOL &na_frame,
                     const std::vector<mol_info> &pro_info,
                     const std::vector<mol_info> &na_info,
                     const std::vector<class index> &pro_index,
                     const std::vector<class index> &na_index, size_t num_snap,
                     double time_ps);

  void export_results(const std::string &csv_filename);
};

#endif // NUCLEIC_PROFILER_H
