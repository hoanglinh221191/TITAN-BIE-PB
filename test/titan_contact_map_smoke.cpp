#include "titan_contact_map.h"

#include "info.h"
#include "titan_executor.h"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

std::string arguments::set_Buildnumber() { return "titan_contact_map_smoke"; }
std::string arguments::get_build_number() const {
  return "titan_contact_map_smoke";
}

namespace {

namespace fs = std::filesystem;

mol_info make_atom(std::size_t resnr, const std::string &resname,
                   const std::string &atom, std::size_t chainnum = 0,
                   const std::string &chain = "A") {
  mol_info info(resnr, atom, atom, chain, resname);
  info.set_chainnum(chainnum);
  info.set_respos(resnr);
  info.set_internal_resnr(resnr);
  info.set_side_chain(true);
  return info;
}

box make_box() {
  box frame_box;
  frame_box.set_a(60.0);
  frame_box.set_b(60.0);
  frame_box.set_c(60.0);
  frame_box.set_alpha(90.0);
  frame_box.set_beta(90.0);
  frame_box.set_gamma(90.0);
  return frame_box;
}

Multi_MOL make_molecule(const std::vector<mol_info> &atoms,
                        const std::vector<double> &xs,
                        const std::vector<double> &masses,
                        double time_ps = 0.0) {
  Multi_MOL mol;
  mol.info = atoms;
  mol.parameters.resize(atoms.size());
  for (std::size_t i = 0; i < masses.size(); ++i)
    mol.parameters[i].set_mass(masses[i]);

  MOL frame(make_box(), time_ps, "AMBER99SB");
  const std::vector<double> ys(xs.size(), 1.0);
  const std::vector<double> zs(xs.size(), 1.0);
  frame.pos.x.assign(xs.begin(), xs.end());
  frame.pos.y.assign(ys.begin(), ys.end());
  frame.pos.z.assign(zs.begin(), zs.end());
  mol.m_mol.push_back(std::move(frame));
  return mol;
}

std::string read_file(const fs::path &path) {
  std::ifstream in(path);
  std::ostringstream buffer;
  buffer << in.rdbuf();
  return buffer.str();
}

bool contains(const std::string &haystack, const std::string &needle) {
  return haystack.find(needle) != std::string::npos;
}

void require(bool condition, const std::string &message) {
  if (!condition) {
    std::cerr << "FAIL titan_contact_map_smoke: " << message << '\n';
    std::exit(1);
  }
}

fs::path new_output_dir(const std::string &case_name) {
  const fs::path out_dir =
      fs::temp_directory_path() / ("titan_contact_map_smoke_" + case_name);
  fs::remove_all(out_dir);
  fs::create_directories(out_dir);
  return out_dir;
}

void run_sidechain_heavy_regression(TitanExecutor &executor) {
  arguments args;
  args.cmap.flag_cmap = true;
  args.cmap.mode = "sidechain-heavy";
  args.cmap.cutoff_ang = 4.5;
  args.cmap.domain = "intra";
  args.cmap.min_seqsep = 3;

  Multi_MOL mol = make_molecule(
      {make_atom(1, "ALA", "CB"), make_atom(2, "GLY", "CA"),
       make_atom(3, "VAL", "CB"), make_atom(4, "THR", "CB"),
       make_atom(5, "LEU", "CB"), make_atom(5, "LEU", "CG"),
       make_atom(9, "SER", "CB")},
      {1.0, 20.0, 20.0, 20.0, 3.0, 3.2, 20.0},
      {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0});

  const fs::path out_dir = new_output_dir("sidechain_heavy");
  std::fstream log((out_dir / "run.log").string(),
                   std::ios::out | std::ios::trunc);

  TitanContactMapAnalyzer analyzer;
  analyzer.configure(args.cmap);
  analyzer.init(mol, nullptr, args, log);
  require(analyzer.enabled(), "sidechain-heavy analyzer should be enabled");
  analyzer.process_frame(mol, nullptr, 0, executor);

  info output_info(args);
  const std::vector<std::string> input_line{"sidechain-heavy.pdb"};
  analyzer.finalize_trajectory(out_dir.string(), output_info, args, input_line,
                               0, log);

  const std::string summary =
      read_file(out_dir / "contact_probability_summary.dat");
  const std::string edges =
      read_file(out_dir / "contact_probability_edges.csv");
  const std::string matrix =
      read_file(out_dir / "contact_probability_matrix.csv");

  require(contains(summary, "mode: sidechain-heavy"),
          "explicit sidechain-heavy mode remains supported");
  require(contains(edges, "ALA") && contains(edges, "LEU"),
          "sidechain-heavy contact edge remains exported");
  require(contains(edges, "1.000000"),
          "sidechain-heavy contact probability remains one");
  require(contains(matrix, "mol1|A:1:ALA") &&
              contains(matrix, "mol1|A:5:LEU"),
          "mol1 matrix labels identify molecule role");
}

void run_default_sidechain_com_contract(TitanExecutor &executor) {
  arguments args;
  args.cmap.flag_cmap = true;
  args.cmap.min_seqsep = 0;

  // ALA mass-weighted side-chain COM is x=15.0; its geometric center is x=10.
  // LEU is at x=21, so only the mass-weighted COM lies within the default
  // 6.5 A threshold. SER at x=0 shares an atom-level contact with ALA but is
  // outside the COM cutoff, guarding against heavy-atom minimum behavior.
  Multi_MOL mol = make_molecule(
      {make_atom(1, "ALA", "CB"), make_atom(1, "ALA", "CG"),
       make_atom(4, "LEU", "CB"), make_atom(8, "SER", "CB")},
      {0.0, 20.0, 21.0, 0.0}, {1.0, 3.0, 1.0, 1.0}, 2.5);

  const fs::path out_dir = new_output_dir("sidechain_com");
  std::fstream log((out_dir / "run.log").string(),
                   std::ios::out | std::ios::trunc);

  TitanContactMapAnalyzer analyzer;
  analyzer.configure(args.cmap);
  analyzer.init(mol, nullptr, args, log);
  analyzer.process_frame(mol, nullptr, 0, executor);

  info output_info(args);
  const std::vector<std::string> input_line{"sidechain-com.pdb"};
  analyzer.finalize_trajectory(out_dir.string(), output_info, args, input_line,
                               0, log);

  const std::string summary =
      read_file(out_dir / "contact_probability_summary.dat");
  const std::string edges =
      read_file(out_dir / "contact_probability_edges.csv");
  const std::string timeseries =
      read_file(out_dir / "contact_count_timeseries.csv");
  const std::string occupancy =
      read_file(out_dir / "contact_residue_occupancy.csv");

  require(contains(summary, "mode: sidechain-com"),
          "default CMap mode is sidechain-com");
  require(contains(summary, "cutoff_A: 6.500000"),
          "default sidechain-com cutoff is 6.5 A");
  require(contains(edges, "ALA") && contains(edges, "LEU") &&
              contains(edges, "1.000000"),
          "mass-weighted COM contact within 6.5 A is exported");
  require(!contains(edges, "SER"),
          "atom-level proximity does not count outside COM cutoff");
  require(contains(timeseries, "frame,time_ps,contact_pairs") &&
              contains(timeseries, ",1"),
          "contact count timeseries output is exported");
  require(contains(occupancy, "occupancy_probability") &&
              contains(occupancy, "ALA") && contains(occupancy, "LEU"),
          "per-residue occupancy output is exported");
}

void run_interface_rectangular_map_contract(TitanExecutor &executor) {
  arguments args;
  args.cmap.flag_cmap = true;
  args.cmap.mode = "sidechain-com";
  args.cmap.cutoff_ang = 6.5;
  args.cmap.scope = "interface";

  Multi_MOL mol1 =
      make_molecule({make_atom(1, "ALA", "CB", 0, "A")},
                    {10.0}, {1.0}, 3.0);
  Multi_MOL mol2 =
      make_molecule({make_atom(1, "LEU", "CB", 0, "B"),
                     make_atom(2, "SER", "CB", 0, "B")},
                    {15.0, 30.0}, {1.0, 1.0}, 3.0);

  const fs::path out_dir = new_output_dir("interface");
  std::fstream log((out_dir / "run.log").string(),
                   std::ios::out | std::ios::trunc);

  TitanContactMapAnalyzer analyzer;
  analyzer.configure(args.cmap);
  analyzer.init(mol1, &mol2, args, log);
  analyzer.process_frame(mol1, &mol2, 0, executor);

  info output_info(args);
  const std::vector<std::string> input_line{"interface.pdb"};
  analyzer.finalize_trajectory(out_dir.string(), output_info, args, input_line,
                               0, log);

  const std::string summary =
      read_file(out_dir / "contact_probability_summary.dat");
  const std::string edges =
      read_file(out_dir / "contact_probability_edges.csv");
  const std::string matrix =
      read_file(out_dir / "contact_probability_matrix.csv");
  const std::string timeseries =
      read_file(out_dir / "contact_count_timeseries.csv");
  const std::string occupancy =
      read_file(out_dir / "contact_residue_occupancy.csv");

  require(contains(summary, "scope: interface") &&
              contains(summary, "matrix_layout: rectangular"),
          "interface summary describes a rectangular matrix");
  require(contains(matrix, "mol1|A:1:ALA") &&
              contains(matrix, "mol2|B:1:LEU") &&
              contains(matrix, "mol2|B:2:SER"),
          "interface matrix distinguishes mol1 rows and mol2 columns");
  require(contains(edges, "molecule_i,molecule_j") &&
              contains(edges, "mol1") && contains(edges, "mol2"),
          "interface edges expose molecule roles");
  require(contains(timeseries, "frame,time_ps,contact_pairs"),
          "interface timeseries output is exported");
  require(contains(occupancy, "molecule") && contains(occupancy, "mol1") &&
              contains(occupancy, "mol2"),
          "interface residue occupancy identifies molecule roles");
}

void run_pooled_average_contract(TitanExecutor &executor) {
  arguments args;
  args.cmap.flag_cmap = true;
  args.cmap.min_seqsep = 0;
  const fs::path out_dir = new_output_dir("pooled");
  std::fstream log((out_dir / "run.log").string(),
                   std::ios::out | std::ios::trunc);
  const std::vector<std::string> inputs{"traj-a.pdb", "traj-b.pdb"};
  info output_info(args);
  TitanContactMapAnalyzer analyzer;
  analyzer.configure(args.cmap);

  for (int input = 0; input < 2; ++input) {
    Multi_MOL mol = make_molecule(
        {make_atom(1, "ALA", "CB"), make_atom(4, "LEU", "CB")},
        {10.0, 14.0}, {1.0, 1.0}, static_cast<double>(input));
    analyzer.init(mol, nullptr, args, log);
    analyzer.process_frame(mol, nullptr, 0, executor);
    analyzer.finalize_trajectory(out_dir.string(), output_info, args, inputs,
                                 input, log);
  }
  analyzer.finalize_all(out_dir.string(), output_info, log);
  const std::string summary =
      read_file(out_dir / "contact_probability_summary_average.dat");
  const std::string matrix =
      read_file(out_dir / "contact_probability_matrix_average.csv");
  require(contains(summary, "valid_frames: 2"),
          "pooled summary accumulates frames from all trajectories");
  require(contains(matrix, "1.000000"),
          "pooled matrix exports pooled contact probability");
}

} // namespace

int main() {
  TitanExecutor executor(4);
  run_sidechain_heavy_regression(executor);
  run_default_sidechain_com_contract(executor);
  run_interface_rectangular_map_contract(executor);
  run_pooled_average_contract(executor);
  std::cout << "PASS titan_contact_map_smoke\n";
  return 0;
}
