#ifndef INFO_H
#define INFO_H
#include "arguments.h"
#include "groups.h"
#include <string>

using std::string;
class info {
  string Build, PdbFileName, TotalEnergyFileName, TotalEnergyTimeFileName,
      SideChainName, LogFile, EmapLgFileName, EmapRecFileName, EntropyFile,
      EntropyQuasiFile, group_1, group_2, stat_hbonds, stat_average_HB,
      stat_maxforce_res, emap, dis_energy_2mol, res_list_hb_rec, res_list_hb_lg,
      side_chain_out, side_chain_time, out_energy_res_motif, dis_com_2mol,
      dipole_filename;
  double current_begin_time, current_end_time;
  double current_dt;

public:
  info(const arguments &Arguments)
      : Build{Arguments.get_build_number()}, PdbFileName{Arguments.pdb_name},
        TotalEnergyFileName{"cal.txt"}, TotalEnergyTimeFileName{"cal.xvg"},
        SideChainName{"out_sc.txt"}, LogFile{"log.txt"},
        EmapLgFileName{"ener_lg.txt"}, EmapRecFileName{"ener_rec.txt"},
        EntropyFile{"out_entropy.dat"},
        EntropyQuasiFile{"out_entropy_quasiHarmonic.dat"}, group_1{"mol1"},
        group_2{"mol2"}, stat_hbonds{"stat_hydrogenbonds.dat"},
        stat_average_HB{"stat_average_HB.dat"},
        stat_maxforce_res{"_out_stat_maxforce_res.dat"}, emap{"out_emap_"},
        dis_energy_2mol{"out_dis_com_energy.dat"},
        dis_com_2mol{"out_dis_com.dat"}, res_list_hb_rec{}, res_list_hb_lg{},
        out_energy_res_motif{"_stat_energy_res.dat"},
        current_begin_time{1000000000000.}, current_end_time{}, current_dt{},
        side_chain_out{"out_sc.txt"}, side_chain_time{"out_sc_time.dat"} {}
  string getBuild() const { return Build; }
  int get_info_group(const groups_cal &groups) {
    if (groups.get_parse_index_flag()) {
      group_1 = groups.group_name(groups.group_1_index());
      group_2 = groups.group_name(groups.group_2_index());
    } else {
    }
    return 1;
  }
  void setPdbName(const arguments &Arguments) {
    PdbFileName = Arguments.pdb_name;
  }
  string get_group_1_name() const { return group_1; }
  string get_group_2_name() const { return group_2; }
  string get_stathb_name() const { return stat_hbonds; }
  string get_stat_ave_hb_name() const { return stat_average_HB; }
  string get_stat_maxforce_name() const { return stat_maxforce_res; }
  string get_emap() const { return emap; }
  string getPdbFileName() const { return PdbFileName; }
  string get_pdboutname(const arguments &Arguments) const {
    return Arguments.get_pdbOutput();
  }
  string getEntropyQuasi() const { return EntropyQuasiFile; }
  string get_res_list_hb_rec_filename() {
    return "out_list_res_HB_" + group_1 + ".dat";
  }
  string get_res_list_hb_lg_filename() {
    return "out_list_res_HB_" + group_2 + ".dat";
  }
  string get_res_list_sc_contact_rec_filename() {
    return "out_list_res_sidechain_contact_" + group_1 + ".dat";
  }
  string get_res_list_sc_contact_lg_filename() {
    return "out_list_res_sidechain_contact_" + group_2 + ".dat";
  }
  string getDisCom2Mol() const { return dis_com_2mol; }
  string getDisEner2Mol() const { return dis_energy_2mol; }
  string getDipoleFilename() const { return "_dipole_moment.dat"; }
  string getDistributionDistanceCharged() const {
    return "_distributionChargedAtoms.dat";
  }
  string getDistributionAngleDipoleEfield() const {
    return "_distributionDipoleEfield.dat";
  }
  string getOutputDistanceRes2COM() const { return "_distances_res_2COM.dat"; }
  string getOuputAverageContactMap() const { return "average_contact_map.dat"; }
  string getOuputEnerMaptypeAve() const { return "_energy_maptype.dat"; }
  string getOutputAverageEmap() const { return "average_emap.dat"; }
  string getOutputAverageEmapMatrix() const {
    return "average_emap_matrix.dat";
  }
  string getOutputEnerPdbtypeAve() const { return "_energy_pdbtype.pdb"; }

  string get_DE_output_filename() { return "DE_out_vector.dat"; }
  string get_DE_log_filename() { return "DE_log.dat"; }
  string get_DE_pdb_vector_filename() { return "DE_vector.pdb"; }

  string getSCout() const { return side_chain_out; }
  string getContmap() const { return "contact_map_traj_"; }
  string getSCtimeout() const { return side_chain_time; }
  void set_current_begin_time(double t) { current_begin_time = t; }
  double get_current_begin_time() const { return current_begin_time; }
  void set_current_end_time(double t) { current_end_time = t; }
  double get_current_end_time() const { return current_end_time; }
  void set_current_dt(double dt) { current_dt = dt; }
  double get_current_dt() const { return current_dt; }
  void set_out_energy_res_motif(string motif) { out_energy_res_motif = motif; }
  string get_out_energy_res_motif() { return out_energy_res_motif; }
  string get_out_mindis_subgroups() { return "mininum_distance_"; }
  string get_out_histogram_mindis() { return "histo_mindistance_"; }
};
#endif
