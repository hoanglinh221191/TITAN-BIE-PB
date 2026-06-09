#include "titan_hotspot_report.h"

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

namespace fs = std::filesystem;

namespace {

void require(bool condition, const std::string &message) {
  if (!condition) {
    std::cerr << "FAIL titan_hotspot_report_smoke: " << message << '\n';
    std::exit(1);
  }
}

std::string read_file(const fs::path &path) {
  std::ifstream in(path);
  return std::string(std::istreambuf_iterator<char>(in),
                     std::istreambuf_iterator<char>());
}

fs::path fresh_dir(const std::string &name) {
  fs::path path = fs::temp_directory_path() / ("titan_hotspot_" + name);
  fs::remove_all(path);
  fs::create_directories(path);
  return path;
}

void write_common_inputs(const fs::path &dir) {
  {
    std::ofstream out(dir / "cas_receptor.dat");
    out << "role,chainnum,chain,resnr,resname,mutation_status,"
           "occupancy_fraction,ddG_total_mean,ddG_total_bound_mean,"
           "ddG_elec_mean,ddG_vdw_mean,ddG_polar_mean,ddG_nonpolar_mean,"
           "avg_distance_A,contact_frames,total_frames\n";
    out << "mol1,0,A,7,LYS,valid_mutation,1.0,2.0,4.0,1.0,1.0,1.0,1.0,"
           "4.5,2,2\n";
    out << "mol1,1,B,7,ASP,valid_mutation,1.0,-8.0,-10.0,-2.0,-2.0,-3.0,"
           "-3.0,5.0,2,2\n";
    out << "mol1,1,C,9,ALA,noop_alanine,1.0,0.0,0.0,0.0,0.0,0.0,0.0,"
           "5.0,2,2\n";
  }
  {
    std::ofstream out(dir / "cas_ligand.dat");
    out << "role,chainnum,chain,resnr,resname,mutation_status,"
           "occupancy_fraction,ddG_total_mean,ddG_total_bound_mean,"
           "ddG_elec_mean,ddG_vdw_mean,ddG_polar_mean,ddG_nonpolar_mean,"
           "avg_distance_A,contact_frames,total_frames\n";
    out << "mol2,0,L,3,ARG,valid_mutation,0.5,3.0,6.0,1.0,1.0,1.0,1.0,"
           "6.0,1,2\n";
  }
  {
    std::ofstream out(dir / "edecomp_mm_residue_ranking.csv");
    out << "rank,molecule_role,entity_index,entity_label,elec_signed_mean,"
           "elec_abs_mean,elec_mean_abs,vdw_signed_mean,vdw_abs_mean,"
           "vdw_mean_abs,signed_mean,abs_mean,mean_abs\n";
    out << "1,mol1,0,mol1|B:7:ASP,-5,5,5,-1,1,1,-6,6,6\n";
    out << "2,mol2,0,mol2|L:3:ARG,2,2,2,1,1,1,3,3,3\n";
    out << "3,mol1,1,mol1|A:7:LYS,1,1,1,1,1,1,2,2,2\n";
    out << "4,mol1,2,mol1|D:8:TYR,1,1,1,0,0,0,1,1,1\n";
  }
  {
    std::ofstream out(dir / "edecomp_mmgbsa_residue_ranking.csv");
    out << "rank,molecule_role,entity_index,entity_label,E_elec,E_vdw,E_MM,"
           "G_polar_GB,G_nonpolar_SASA,DeltaG_no_entropy,abs_mean,mean_abs,"
           "pooled_frame_mean,samples\n";
    out << "1,mol1,0,mol1|A:7:LYS,1,1,2,3,4,9,9,9,9,2\n";
  }
  {
    std::ofstream out(dir / "edecomp_mmpbsa_residue_ranking.csv");
    out << "rank,molecule_role,entity_index,entity_label,E_elec,E_vdw,E_MM,"
           "DeltaDeltaG_solv_PB,DeltaG_nonpolar_SASA,DeltaG_no_entropy,"
           "abs_mean,mean_abs,pooled_frame_mean,samples\n";
    out << "1,mol1,0,mol1|A:7:LYS,1,1,2,3,4,11,11,11,11,2\n";
  }
  {
    std::ofstream out(dir / "cas_hotspots_b_factor.pdb");
    out << "REMARK existing CAS template\n";
    out << "ATOM      1  CA  LYS A   7       0.000   0.000   0.000  1.00  0.00           C\n";
    out << "ATOM      2  CA  ASP B   7       1.000   0.000   0.000  1.00  0.00           C\n";
    out << "END\n";
  }
}

void check_combined_report() {
  const fs::path dir = fresh_dir("combined");
  write_common_inputs(dir);

  titan_hotspot::ReportConfig config;
  config.source = "both";
  config.rank_by = "combined";
  config.edecomp_source = "mm";
  const auto result = titan_hotspot::export_report(dir, config);

  require(result.rows_written == 5, "writes joined CAS/EDecomp residue rows");
  require(result.matched_rows == 3, "matches CAS and EDecomp by role/chain/resnr/resname");
  require(result.cas_only_rows == 1, "tracks CAS-only residues");
  require(result.edecomp_only_rows == 1, "tracks EDecomp-only residues");
  require(result.selected_edecomp_source == "mm",
          "explicit MM EDecomp source is used for merge");
  require(result.bfactor_written, "rewrites B-factor PDB when template exists");

  const std::string ranking = read_file(dir / "hotspot_residue_ranking.csv");
  require(ranking.find("source_coverage") != std::string::npos,
          "ranking exports source coverage column");
  require(ranking.find("mol1,1,B,7,ASP") != std::string::npos,
          "duplicate residue numbers stay distinct by chain");
  require(ranking.find("dimensionless") != std::string::npos,
          "ranking notes dimensionless score semantics");
  require(ranking.find("cas_only") != std::string::npos,
          "ranking includes CAS-only coverage");
  require(ranking.find("edecomp_only") != std::string::npos,
          "ranking includes EDecomp-only coverage");

  const std::string summary = read_file(dir / "hotspot_summary.dat");
  require(summary.find("CAS method: GBSA alanine-like side-chain neutralization") !=
              std::string::npos,
          "summary documents CAS method");
  require(summary.find("EDecomp note: endpoint entity decomposition") !=
              std::string::npos,
          "summary documents EDecomp endpoint semantics");
  require(summary.find("coverage_warning:") != std::string::npos,
          "summary reports source coverage mismatch");

  const std::string pdb = read_file(dir / "hotspot_bfactor.pdb");
  require(pdb.find("REMARK TITAN hotspot B-factor stores dimensionless") !=
              std::string::npos,
          "PDB remark documents score semantics");

  titan_hotspot::ReportConfig all;
  all.source = "edecomp";
  all.rank_by = "edecomp-abs-mean";
  all.edecomp_source = "all";
  const auto all_result = titan_hotspot::export_report(dir, all);
  require(all_result.selected_edecomp_source == "mmpbsa",
          "edecomp-source=all prefers MMPBSA when available");
}

void check_source_specific_reports() {
  const fs::path dir = fresh_dir("source_specific");
  write_common_inputs(dir);

  titan_hotspot::ReportConfig cas;
  cas.source = "cas";
  cas.rank_by = "cas-ddg";
  auto cas_result = titan_hotspot::export_report(dir, cas);
  require(cas_result.edecomp_rows == 0, "CAS-only report skips EDecomp rows");
  require(read_file(dir / "hotspot_residue_ranking.csv").find("cas_only") !=
              std::string::npos,
          "CAS-only report labels coverage");

  titan_hotspot::ReportConfig edecomp;
  edecomp.source = "edecomp";
  edecomp.rank_by = "edecomp-abs-mean";
  edecomp.edecomp_source = "mm";
  auto edecomp_result = titan_hotspot::export_report(dir, edecomp);
  require(edecomp_result.cas_rows == 0, "EDecomp-only report skips CAS rows");
  require(edecomp_result.selected_edecomp_source == "mm",
          "explicit MM EDecomp source is used");
}

} // namespace

int main() {
  check_combined_report();
  check_source_specific_reports();
  std::cout << "PASS titan_hotspot_report_smoke\n";
  return 0;
}
