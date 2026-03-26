#ifndef GBSA_INTRIN_H
#define GBSA_INTRIN_H

#include "MOL.h"
#include "arguments.h"
#include "cache_sigma_epsilon.h"
#include "info.h"
#include "stat_ener_res_res.h"
#include <map>
#include <string>
#include <unordered_map>

#include <vector>

// Structure to hold GBSA outputs per snap
struct GBSA_energy_snap {
  double total_energy;
  double E_vdw;
  double E_elec;
  double G_polar;
  double G_nonpolar;
  double time;
  GBSA_energy_snap()
      : total_energy(0.), E_vdw(0.), E_elec(0.), G_polar(0.), G_nonpolar(0.),
        time(0.) {}
};

struct GBSA_decomposition {
  bool is_ligand;
  int chainnum;
  int resnr;
  std::string unit_type_str;
  std::string chain_str;
  std::string name_str;
  double E_vdw;
  double E_elec;
  double G_polar;
  double G_nonpolar;
  double total;
  size_t count;
  std::vector<double> frame_energies;
  GBSA_decomposition()
      : is_ligand(false), chainnum(0), resnr(0), E_vdw(0.), E_elec(0.),
        G_polar(0.), G_nonpolar(0.), total(0.), count(0) {}
};

class GBSA_Aggregator {
private:
  std::vector<std::vector<GBSA_energy_snap>> trajectories;
  std::vector<GBSA_decomposition> decomp;

  using DecompKey = std::tuple<bool, int, int>; // (is_ligand, chainnum, resnr)
  struct DecompKeyHash {
    size_t operator()(const DecompKey &k) const {
      size_t h = std::hash<bool>()(std::get<0>(k));
      h ^= std::hash<int>()(std::get<1>(k)) + 0x9e3779b9 + (h << 6) + (h >> 2);
      h ^= std::hash<int>()(std::get<2>(k)) + 0x9e3779b9 + (h << 6) + (h >> 2);
      return h;
    }
  };
  std::unordered_map<DecompKey, size_t, DecompKeyHash> decomp_index;

public:
  struct pair_hash {
    template <class T1, class T2>
    std::size_t operator()(const std::pair<T1, T2> &p) const {
      return std::hash<T1>()(p.first) ^ (std::hash<T2>()(p.second) << 1);
    }
  };
  std::unordered_map<std::pair<int, int>, double, pair_hash> pairwise_energies;

private:
  const arguments &Arguments;
  const info &info_output;
  double saltcon, ext_diel, int_diel;
  int igb;
  std::string prefix;

public:
  GBSA_Aggregator(const arguments &arg, const info &info)
      : Arguments(arg), info_output(info) {
    saltcon = arg.saltcon;
    ext_diel = arg.ext_diel;
    int_diel = arg.int_diel;
    igb = arg.igb;
    prefix = "";
  }

  void add_snapshot(size_t traj_index, GBSA_energy_snap snap);
  void add_decomposition(bool is_ligand, int chainnum, int resnr,
                         bool detail_atom, const mol_info *m_info, double ev,
                         double ee, double gpol, double gnp);
  void add_pairwise(int rec_id, int lig_id, double energy) {
    pairwise_energies[{rec_id, lig_id}] += energy;
  }
  void merge(const GBSA_Aggregator &other);
  void commit_frame_decompositions(const GBSA_Aggregator &frame_data,
                                   size_t current_frame_idx);
  void print_all_stats(const std::string &dir_prefix);
  void initialize_trajectory(size_t traj_index);

  const std::vector<GBSA_decomposition> &get_decomp() const { return decomp; }
};

#ifdef __AVX__

// Calculate Effective Born Radii AVX2 Kernel
void cal_effective_born_radii_avx2(
    const MOL &mol, const AlignedDoubleVector &vdw_radii,
    AlignedDoubleVector &radii, int igb, size_t begin, size_t end,
    const MOL *env_mol = nullptr,
    const AlignedDoubleVector *env_vdw_radii = nullptr);

int cal_energy_gbsa_avx2(
    const MOL &rec, const MOL &lg, const std::vector<mol_info> &info_rec,
    const std::vector<mol_info> &info_lg, cache_sigma_epsilon &CachedData,
    size_t rec_begin, size_t rec_end, size_t lg_begin, size_t lg_end,
    GBSA_energy_snap &Energy, size_t num_snap, const arguments &Arguments,
    const Multi_MOL_parameters &rec_parameters,
    const Multi_MOL_parameters &lg_parameters,
    const std::vector<double> &rec_lcpo_p1,
    const std::vector<double> &rec_lcpo_p2,
    const std::vector<double> &lg_lcpo_p1,
    const std::vector<double> &lg_lcpo_p2,
    const AlignedDoubleVector &rec_free_radii,
    const AlignedDoubleVector &lg_free_radii,
    const AlignedDoubleVector &rec_complex_radii,
    const AlignedDoubleVector &lg_complex_radii, GBSA_Aggregator &aggregator,
    info_ener_res *infoener_res = nullptr);

#endif // __AVX__

#endif // GBSA_INTRIN_H
