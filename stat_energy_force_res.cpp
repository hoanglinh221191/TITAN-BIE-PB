#include "arguments.h"
#include "folder.h"
#include "info.h"
#include "miscellaneous.h"
#include "titan_energy_stats.h"

#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

using namespace std;

namespace {

double energy_out(double value, const arguments &Arguments) {
  return value / Arguments.eu_par;
}

double force_out(double value, const arguments &Arguments) {
  return value * Arguments.fu_par;
}

void print_energy_component_table(
    const vector<titan_energy_stats::ResidueEnergySummary> &summaries,
    const arguments &Arguments, const string &filename,
    double titan_energy_stats::ResidueEnergySummary::*mean_member,
    double titan_energy_stats::ResidueEnergySummary::*sd_member) {
  fstream out(filename, ios::out);
  out << "# ChainResidResname  Average energy (" << Arguments.eu << ")    SD ("
      << Arguments.eu << ")\n";
  for (const auto &summary : summaries) {
    out << summary.key.chain << setw(5) << left << summary.key.resnr << setw(4)
        << summary.key.resname << right << setw(12) << setprecision(6)
        << energy_out(summary.*mean_member, Arguments) << setw(12)
        << setprecision(6) << energy_out(summary.*sd_member, Arguments)
        << '\n';
  }
}

void print_force_component_table(
    const vector<titan_energy_stats::ResidueEnergySummary> &summaries,
    const arguments &Arguments, const string &filename,
    double titan_energy_stats::ResidueEnergySummary::*mean_member,
    double titan_energy_stats::ResidueEnergySummary::*sd_member) {
  fstream out(filename, ios::out);
  out << "# ChainResidResname  Average force (" << Arguments.fu << ")    SD ("
      << Arguments.fu << ")\n";
  for (const auto &summary : summaries) {
    out << summary.key.chain << setw(5) << left << summary.key.resnr << setw(4)
        << summary.key.resname << right << setw(12) << setprecision(6)
        << force_out(summary.*mean_member, Arguments) << setw(12)
        << setprecision(6) << force_out(summary.*sd_member, Arguments) << '\n';
  }
}

void print_energy_stat_tsv(
    const vector<titan_energy_stats::ResidueEnergySummary> &summaries,
    info &info_output, const arguments &Arguments, const string &name) {
  (void)info_output;
  fstream out("out_" + name + "_energy_stat_res.tsv", ios::out);
  out << "chain\tresnr\tresname\tsamples\tmean_elec_" << Arguments.eu
      << "\tsd_elec_" << Arguments.eu << "\tmean_vdw_" << Arguments.eu
      << "\tsd_vdw_" << Arguments.eu << "\tmean_total_" << Arguments.eu
      << "\tsd_total_" << Arguments.eu << '\n';
  out << setprecision(10);
  for (const auto &summary : summaries) {
    out << titan_energy_stats::tsv_escape(summary.key.chain) << '\t'
        << summary.key.resnr << '\t'
        << titan_energy_stats::tsv_escape(summary.key.resname) << '\t'
        << summary.samples << '\t'
        << energy_out(summary.mean_elec, Arguments) << '\t'
        << energy_out(summary.sd_elec, Arguments) << '\t'
        << energy_out(summary.mean_vdw, Arguments) << '\t'
        << energy_out(summary.sd_vdw, Arguments) << '\t'
        << energy_out(summary.mean_total, Arguments) << '\t'
        << energy_out(summary.sd_total, Arguments) << '\n';
  }
}

void print_force_stat_tsv(
    const vector<titan_energy_stats::ResidueEnergySummary> &summaries,
    const arguments &Arguments, const string &name) {
  fstream out("out_" + name + "_force_stat_res.tsv", ios::out);
  out << "chain\tresnr\tresname\tsamples\tmean_elec_mag_" << Arguments.fu
      << "\tsd_elec_mag_" << Arguments.fu << "\tmean_vdw_mag_" << Arguments.fu
      << "\tsd_vdw_mag_" << Arguments.fu << "\tmean_total_mag_"
      << Arguments.fu << "\tsd_total_mag_" << Arguments.fu
      << "\tmean_elec_fx\tmean_elec_fy\tmean_elec_fz"
      << "\tsd_elec_fx\tsd_elec_fy\tsd_elec_fz"
      << "\tmean_vdw_fx\tmean_vdw_fy\tmean_vdw_fz"
      << "\tsd_vdw_fx\tsd_vdw_fy\tsd_vdw_fz"
      << "\tmean_total_fx\tmean_total_fy\tmean_total_fz"
      << "\tsd_total_fx\tsd_total_fy\tsd_total_fz\n";
  out << setprecision(10);
  for (const auto &summary : summaries) {
    out << titan_energy_stats::tsv_escape(summary.key.chain) << '\t'
        << summary.key.resnr << '\t'
        << titan_energy_stats::tsv_escape(summary.key.resname) << '\t'
        << summary.samples << '\t'
        << force_out(summary.mean_force_elec_mag, Arguments) << '\t'
        << force_out(summary.sd_force_elec_mag, Arguments) << '\t'
        << force_out(summary.mean_force_vdw_mag, Arguments) << '\t'
        << force_out(summary.sd_force_vdw_mag, Arguments) << '\t'
        << force_out(summary.mean_force_total_mag, Arguments) << '\t'
        << force_out(summary.sd_force_total_mag, Arguments) << '\t'
        << force_out(summary.mean_force_elec.Fx, Arguments) << '\t'
        << force_out(summary.mean_force_elec.Fy, Arguments) << '\t'
        << force_out(summary.mean_force_elec.Fz, Arguments) << '\t'
        << force_out(summary.sd_force_elec.Fx, Arguments) << '\t'
        << force_out(summary.sd_force_elec.Fy, Arguments) << '\t'
        << force_out(summary.sd_force_elec.Fz, Arguments) << '\t'
        << force_out(summary.mean_force_vdw.Fx, Arguments) << '\t'
        << force_out(summary.mean_force_vdw.Fy, Arguments) << '\t'
        << force_out(summary.mean_force_vdw.Fz, Arguments) << '\t'
        << force_out(summary.sd_force_vdw.Fx, Arguments) << '\t'
        << force_out(summary.sd_force_vdw.Fy, Arguments) << '\t'
        << force_out(summary.sd_force_vdw.Fz, Arguments) << '\t'
        << force_out(summary.mean_force_total.Fx, Arguments) << '\t'
        << force_out(summary.mean_force_total.Fy, Arguments) << '\t'
        << force_out(summary.mean_force_total.Fz, Arguments) << '\t'
        << force_out(summary.sd_force_total.Fx, Arguments) << '\t'
        << force_out(summary.sd_force_total.Fy, Arguments) << '\t'
        << force_out(summary.sd_force_total.Fz, Arguments) << '\n';
  }
}

} // namespace

int print_stat_energy_res(
    vector<titan_energy_stats::ResidueEnergySummary> &summaries,
    info &info_output, const arguments &Arguments, string name) {
  print_energy_component_table(
      summaries, Arguments,
      "out_" + name + "_electrostatic_energy" +
          info_output.get_out_energy_res_motif(),
      &titan_energy_stats::ResidueEnergySummary::mean_elec,
      &titan_energy_stats::ResidueEnergySummary::sd_elec);
  print_energy_component_table(
      summaries, Arguments,
      "out_" + name + "_vanderwaals_energy" +
          info_output.get_out_energy_res_motif(),
      &titan_energy_stats::ResidueEnergySummary::mean_vdw,
      &titan_energy_stats::ResidueEnergySummary::sd_vdw);
  print_energy_component_table(
      summaries, Arguments,
      "out_" + name + "_total_energy" + info_output.get_out_energy_res_motif(),
      &titan_energy_stats::ResidueEnergySummary::mean_total,
      &titan_energy_stats::ResidueEnergySummary::sd_total);
  print_energy_stat_tsv(summaries, info_output, Arguments, name);

  if (Arguments.print_force_res == true) {
    print_force_component_table(
        summaries, Arguments,
        "out_" + name + "_electrostatic_force" +
            info_output.get_out_energy_res_motif(),
        &titan_energy_stats::ResidueEnergySummary::mean_force_elec_mag,
        &titan_energy_stats::ResidueEnergySummary::sd_force_elec_mag);
    print_force_component_table(
        summaries, Arguments,
        "out_" + name + "_vanderwaals_force" +
            info_output.get_out_energy_res_motif(),
        &titan_energy_stats::ResidueEnergySummary::mean_force_vdw_mag,
        &titan_energy_stats::ResidueEnergySummary::sd_force_vdw_mag);
    print_force_component_table(
        summaries, Arguments,
        "out_" + name + "_total_force" + info_output.get_out_energy_res_motif(),
        &titan_energy_stats::ResidueEnergySummary::mean_force_total_mag,
        &titan_energy_stats::ResidueEnergySummary::sd_force_total_mag);
    print_force_stat_tsv(summaries, Arguments, name);
  }
  return 0;
}

int stat_energy_res(vector<vector<energy_res>> &Energy_stat_res, int file_num,
                    const arguments &Arguments, info &info_output,
                    string mol_name) {
  (void)file_num;
  auto summaries = titan_energy_stats::summarize_residue_energy(Energy_stat_res);
  if (summaries.empty())
    return 0;
  print_stat_energy_res(summaries, info_output, Arguments, mol_name);
  return 0;
}
