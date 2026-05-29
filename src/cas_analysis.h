#ifndef CAS_ANALYSIS_H
#define CAS_ANALYSIS_H

#include "MOL.h"
#include "arguments.h"
#include "functions.h"
#include "gbsa_intrin.h"
#include "groups.h"
#include <mutex>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

struct CAS_Result {
  int res_id;
  std::string res_name;
  double ddG_vdw;
  double ddG_elec;
  double ddG_polar;
  double ddG_nonpolar;
  double ddG_total;
  double distance_sum;
  int contact_frames;

  CAS_Result()
      : res_id(0), res_name(""), ddG_vdw(0.0), ddG_elec(0.0), ddG_polar(0.0),
        ddG_nonpolar(0.0), ddG_total(0.0), distance_sum(0.0),
        contact_frames(0) {}
};

struct CASMutationAtom {
  size_t atom_index{};
  double charge{};
  double epsilon{};
  double sigma{};
  double radius{};
};

class CASAnalyzer {
private:
  std::unordered_map<int, CAS_Result> cas_receptor_data;
  std::unordered_map<int, CAS_Result> cas_ligand_data;
  std::unordered_map<int, std::vector<CASMutationAtom>> receptor_mutation_plans;
  std::unordered_map<int, std::vector<CASMutationAtom>> ligand_mutation_plans;
  std::vector<int> unique_res_rec;
  std::vector<int> unique_res_lig;
  size_t global_total_frames;
  bool has_reference;
  MOL ref_pro_frame;
  MOL ref_lig_frame;
  std::vector<mol_info> ref_pro_info;
  std::vector<mol_info> ref_lig_info;

  // (mut_free_radii / mut_complex_radii moved to per-worker locals)

public:
  CASAnalyzer();

  void init(Multi_MOL &pro_in, Multi_MOL &lig, const std::string &mol_2_type);

  void process_frame(Multi_MOL &pro_in, Multi_MOL &lig, size_t num_snap,
                     size_t num_threads, const std::vector<d_min> &dis_min,
                     cache_sigma_epsilon &CachedData,
                     const arguments &Arguments, info &info_output,
                     const AlignedDoubleVector &rec_free_radii,
                     const AlignedDoubleVector &lg_free_radii,
                     const AlignedDoubleVector &rec_complex_radii,
                     const AlignedDoubleVector &lg_complex_radii,
                     const GBSA_energy_snap &frame_gbsa_total);

  void export_results(const std::string &rec_file, const std::string &lig_file,
                      double eu_par);

  void export_hotspot_pdb(const std::string &out_pdb_filename, double eu_par);
};

#endif // CAS_ANALYSIS_H
