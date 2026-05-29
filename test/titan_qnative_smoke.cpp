#include "titan_qnative.h"

#include "info.h"
#include "titan_executor.h"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

std::string arguments::set_Buildnumber() { return "titan_qnative_smoke"; }
std::string arguments::get_build_number() const {
  return "titan_qnative_smoke";
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
  frame_box.set_a(80.0);
  frame_box.set_b(80.0);
  frame_box.set_c(80.0);
  frame_box.set_alpha(90.0);
  frame_box.set_beta(90.0);
  frame_box.set_gamma(90.0);
  return frame_box;
}

Multi_MOL make_molecule(const std::vector<mol_info> &atoms,
                        const std::vector<std::vector<double>> &frame_xs,
                        const std::vector<double> &masses) {
  Multi_MOL mol;
  mol.info = atoms;
  mol.parameters.resize(atoms.size());
  for (std::size_t i = 0; i < masses.size(); ++i)
    mol.parameters[i].set_mass(masses[i]);

  for (std::size_t frame_idx = 0; frame_idx < frame_xs.size(); ++frame_idx) {
    MOL frame(make_box(), static_cast<double>(frame_idx), "AMBER99SB");
    const std::vector<double> ys(frame_xs[frame_idx].size(), 1.0);
    const std::vector<double> zs(frame_xs[frame_idx].size(), 1.0);
    frame.pos.x.assign(frame_xs[frame_idx].begin(), frame_xs[frame_idx].end());
    frame.pos.y.assign(ys.begin(), ys.end());
    frame.pos.z.assign(zs.begin(), zs.end());
    mol.m_mol.push_back(std::move(frame));
  }
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
    std::cerr << "FAIL titan_qnative_smoke: " << message << '\n';
    std::exit(1);
  }
}

fs::path new_output_dir(const std::string &case_name) {
  const fs::path out_dir =
      fs::temp_directory_path() / ("titan_qnative_smoke_" + case_name);
  fs::remove_all(out_dir);
  fs::create_directories(out_dir);
  return out_dir;
}

void run_hard_q_contract(TitanExecutor &executor) {
  arguments args;
  args.qnative.flag_qnative = true;
  args.qnative.min_seqsep = 0;

  const std::vector<mol_info> atoms{
      make_atom(1, "ALA", "CB"),  make_atom(4, "LEU", "CB"),
      make_atom(8, "SER", "CB"),  make_atom(12, "THR", "CB"),
      make_atom(16, "VAL", "CB"), make_atom(20, "LYS", "CB")};
  Multi_MOL mol = make_molecule(
      atoms,
      {{0.0, 5.0, 20.0, 25.0, 40.0, 45.0},
       {0.0, 5.0, 20.0, 25.0, 40.0, 60.0},
       {0.0, 20.0, 20.0, 40.0, 40.0, 60.0}},
      {1.0, 1.0, 1.0, 1.0, 1.0, 1.0});

  const fs::path out_dir = new_output_dir("hard");
  std::fstream log((out_dir / "run.log").string(),
                   std::ios::out | std::ios::trunc);

  TitanQNativeAnalyzer analyzer;
  analyzer.configure(args.qnative);
  analyzer.init(mol, nullptr, args, log);
  require(analyzer.enabled(), "qnative analyzer should be enabled");
  for (std::size_t frame = 0; frame < mol.m_mol.size(); ++frame)
    analyzer.process_frame(mol, nullptr, frame, executor);

  info output_info(args);
  const std::vector<std::string> input_line{"qnative-hard.pdb"};
  analyzer.finalize_trajectory(out_dir.string(), output_info, args, input_line,
                               0, log);

  const std::string summary = read_file(out_dir / "qnative_summary.dat");
  const std::string timeseries = read_file(out_dir / "qnative_timeseries.csv");
  const std::string native_contacts =
      read_file(out_dir / "qnative_native_contacts.csv");
  const std::string occupancy =
      read_file(out_dir / "qnative_contact_occupancy.csv");

  require(contains(summary, "mode: sidechain-com"),
          "default qnative mode is sidechain-com");
  require(contains(summary, "cutoff_A: 6.500000"),
          "sidechain-com default cutoff is 6.5 A");
  require(contains(summary, "total_native_contacts: 3"),
          "reference native set has three contacts");
  require(contains(timeseries, "frame,time_ps,Q,formed_native_contacts,total_native_contacts"),
          "qnative timeseries header is exported");
  require(contains(timeseries, "0,0.000000,1.000000,3.000000,3"),
          "frame 0 has Q=1.0");
  require(contains(timeseries, "1,1.000000,0.666667,2.000000,3"),
          "frame 1 has partial Q");
  require(contains(timeseries, "2,2.000000,0.000000,0.000000,3"),
          "frame 2 has Q=0.0");
  require(contains(native_contacts, "mol1|A:1:ALA") &&
              contains(native_contacts, "mol1|A:20:LYS"),
          "native contact list uses CMap residue labels");
  require(contains(occupancy, "occupancy_probability") &&
              contains(occupancy, "0.666667"),
          "native-contact occupancy is averaged over frames");
}

void run_interface_contract(TitanExecutor &executor) {
  arguments args;
  args.qnative.flag_qnative = true;
  args.qnative.scope = "interface";

  Multi_MOL mol1 =
      make_molecule({make_atom(1, "ALA", "CB", 0, "A")},
                    {{10.0}, {10.0}}, {1.0});
  Multi_MOL mol2 =
      make_molecule({make_atom(1, "LEU", "CB", 0, "B"),
                     make_atom(2, "SER", "CB", 0, "B")},
                    {{15.0, 35.0}, {30.0, 35.0}}, {1.0, 1.0});

  const fs::path out_dir = new_output_dir("interface");
  std::fstream log((out_dir / "run.log").string(),
                   std::ios::out | std::ios::trunc);

  TitanQNativeAnalyzer analyzer;
  analyzer.configure(args.qnative);
  analyzer.init(mol1, &mol2, args, log);
  for (std::size_t frame = 0; frame < mol1.m_mol.size(); ++frame)
    analyzer.process_frame(mol1, &mol2, frame, executor);

  info output_info(args);
  const std::vector<std::string> input_line{"qnative-interface.pdb"};
  analyzer.finalize_trajectory(out_dir.string(), output_info, args, input_line,
                               0, log);

  const std::string summary = read_file(out_dir / "qnative_summary.dat");
  const std::string native_contacts =
      read_file(out_dir / "qnative_native_contacts.csv");
  const std::string occupancy =
      read_file(out_dir / "qnative_contact_occupancy.csv");
  require(contains(summary, "scope: interface"),
          "interface scope is recorded in summary");
  require(contains(native_contacts, "mol1|A:1:ALA") &&
              contains(native_contacts, "mol2|B:1:LEU"),
          "interface native list uses mol1/mol2 labels");
  require(contains(occupancy, "0.500000"),
          "interface occupancy reflects broken native contact");
}

void run_soft_q_contract(TitanExecutor &executor) {
  arguments args;
  args.qnative.flag_qnative = true;
  args.qnative.min_seqsep = 0;
  args.qnative.soft = true;
  args.qnative.beta = 1.0;

  Multi_MOL mol = make_molecule(
      {make_atom(1, "ALA", "CB"), make_atom(4, "LEU", "CB")},
      {{0.0, 5.0}, {0.0, 6.5}}, {1.0, 1.0});

  const fs::path out_dir = new_output_dir("soft");
  std::fstream log((out_dir / "run.log").string(),
                   std::ios::out | std::ios::trunc);

  TitanQNativeAnalyzer analyzer;
  analyzer.configure(args.qnative);
  analyzer.init(mol, nullptr, args, log);
  for (std::size_t frame = 0; frame < mol.m_mol.size(); ++frame)
    analyzer.process_frame(mol, nullptr, frame, executor);

  info output_info(args);
  const std::vector<std::string> input_line{"qnative-soft.pdb"};
  analyzer.finalize_trajectory(out_dir.string(), output_info, args, input_line,
                               0, log);

  const std::string summary = read_file(out_dir / "qnative_summary.dat");
  const std::string timeseries = read_file(out_dir / "qnative_timeseries.csv");
  require(contains(summary, "soft: true") && contains(summary, "beta: 1.000000"),
          "soft mode and beta are recorded");
  require(contains(timeseries, "1,1.000000,0.500000,0.500000,1"),
          "soft Q is continuous at the cutoff");
}

} // namespace

int main() {
  TitanExecutor executor(4);
  run_hard_q_contract(executor);
  run_interface_contract(executor);
  run_soft_q_contract(executor);
  std::cout << "PASS titan_qnative_smoke\n";
  return 0;
}
