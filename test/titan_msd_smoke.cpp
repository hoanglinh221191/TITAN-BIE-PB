#include "titan_msd.h"

#include "info.h"
#include "titan_executor.h"

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

std::string arguments::set_Buildnumber() { return "titan_msd_smoke"; }
std::string arguments::get_build_number() const { return "titan_msd_smoke"; }

namespace {

namespace fs = std::filesystem;

mol_info make_atom(std::size_t resnr, const std::string &resname,
                   const std::string &atom, const std::string &chain = "A") {
  mol_info info(resnr, atom, atom, chain, resname);
  info.set_chainnum(chain == "A" ? 0 : 1);
  info.set_respos(resnr);
  info.set_internal_resnr(resnr);
  return info;
}

box make_box(double side = 100.0) {
  box frame_box;
  frame_box.set_a(side);
  frame_box.set_b(side);
  frame_box.set_c(side);
  frame_box.set_alpha(90.0);
  frame_box.set_beta(90.0);
  frame_box.set_gamma(90.0);
  return frame_box;
}

Multi_MOL make_two_atom_mol(const std::vector<double> &xs0,
                            const std::vector<double> &xs1,
                            double box_side = 100.0,
                            bool zero_mass = true) {
  Multi_MOL mol;
  mol.info = {make_atom(1, "ALA", "CA"), make_atom(2, "GLY", "CA")};
  mol.parameters.resize(2);
  mol.parameters.mass[0] = zero_mass ? 0.0 : 12.0;
  mol.parameters.mass[1] = zero_mass ? 0.0 : 24.0;
  for (std::size_t frame_idx = 0; frame_idx < xs0.size(); ++frame_idx) {
    MOL frame(make_box(box_side), static_cast<double>(frame_idx), "AMBER99SB");
    frame.pos.x = {xs0[frame_idx], xs1[frame_idx]};
    frame.pos.y = {0.0, 0.0};
    frame.pos.z = {0.0, 0.0};
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
    std::cerr << "FAIL titan_msd_smoke: " << message << '\n';
    std::exit(1);
  }
}

void require_close(double actual, double expected, const std::string &message,
                   double tol = 1.0e-6) {
  if (std::fabs(actual - expected) > tol) {
    std::cerr << "FAIL titan_msd_smoke: " << message << " actual=" << actual
              << " expected=" << expected << '\n';
    std::exit(1);
  }
}

std::vector<double> read_msd_values(const fs::path &path) {
  std::ifstream in(path);
  std::string line;
  std::getline(in, line);
  std::vector<double> values;
  while (std::getline(in, line)) {
    std::istringstream row(line);
    std::size_t lag = 0;
    double time = 0.0;
    double msd = 0.0;
    std::uint64_t samples = 0;
    row >> lag >> time >> msd >> samples;
    values.push_back(msd);
  }
  return values;
}

fs::path new_output_dir(const std::string &case_name) {
  const fs::path out_dir =
      fs::temp_directory_path() / ("titan_msd_smoke_" + case_name);
  fs::remove_all(out_dir);
  fs::create_directories(out_dir);
  return out_dir;
}

void run_analyzer_case(const std::string &case_name, Multi_MOL &mol,
                       const MsdArgs &msd_args, TitanExecutor &executor) {
  arguments args;
  args.msd = msd_args;
  const fs::path out_dir = new_output_dir(case_name);
  std::fstream log((out_dir / "run.log").string(),
                   std::ios::out | std::ios::trunc);
  TitanMSDAnalyzer analyzer;
  analyzer.configure(args.msd);
  analyzer.init(mol, nullptr, args, log);
  require(analyzer.enabled(), "MSD analyzer should be enabled");
  for (std::size_t frame = 0; frame < mol.m_mol.size(); ++frame)
    analyzer.process_frame(mol, nullptr, frame, executor);
  info output_info(args);
  const std::vector<std::string> input_line{"msd-smoke.pdb"};
  analyzer.finalize_trajectory(out_dir.string(), output_info, args, input_line,
                               0, executor, log);
}

void run_constant_velocity_contract(TitanExecutor &executor) {
  Multi_MOL mol = make_two_atom_mol({0.0, 1.0, 2.0, 3.0},
                                    {2.0, 3.0, 4.0, 5.0}, 100.0, true);
  MsdArgs msd;
  msd.flag_msd = true;
  msd.selection = "mol1";
  msd.mode = "com";
  msd.maxlag = "3";
  msd.fit_window = "auto";
  msd.pbc = "auto";
  run_analyzer_case("constant_velocity", mol, msd, executor);

  const fs::path out_dir = fs::temp_directory_path() /
                           "titan_msd_smoke_constant_velocity";
  const std::vector<double> msd_values =
      read_msd_values(out_dir / "msd_timeseries.dat");
  require(msd_values.size() == 3, "maxlag produces three rows");
  require_close(msd_values[0], 1.0, "lag 1 COM MSD");
  require_close(msd_values[1], 4.0, "lag 2 COM MSD");
  require_close(msd_values[2], 9.0, "lag 3 COM MSD");
  const std::string summary =
      read_file(out_dir / "msd_diffusion_summary.dat");
  const std::string alpha = read_file(out_dir / "msd_alpha_fit.dat");
  require(contains(summary, "fit_status: ok") &&
              contains(summary, "D_A2_per_ps: 0.66666667"),
          "diffusion summary records D from Einstein relation");
  require(contains(alpha, "alpha: 2.00000000"),
          "alpha fit detects ballistic exponent");
}

void run_atom_mode_contract(TitanExecutor &executor) {
  Multi_MOL mol = make_two_atom_mol({0.0, 1.0, 2.0}, {0.0, 2.0, 4.0});
  MsdArgs msd;
  msd.flag_msd = true;
  msd.selection = "atom:CA";
  msd.mode = "atom";
  msd.maxlag = "2";
  run_analyzer_case("atom_mode", mol, msd, executor);
  const std::vector<double> values =
      read_msd_values(fs::temp_directory_path() /
                      "titan_msd_smoke_atom_mode" / "msd_timeseries.dat");
  require_close(values[0], 2.5, "atom mode averages two atom tracks at lag 1");
  require_close(values[1], 10.0, "atom mode averages two atom tracks at lag 2");
}

void run_pbc_policy_contract(TitanExecutor &executor) {
  Multi_MOL mol = make_two_atom_mol({9.0, 1.0, 3.0, 5.0},
                                    {9.0, 1.0, 3.0, 5.0}, 10.0, true);
  MsdArgs unwrap;
  unwrap.flag_msd = true;
  unwrap.selection = "mol1";
  unwrap.mode = "com";
  unwrap.maxlag = "3";
  unwrap.pbc = "auto";
  run_analyzer_case("pbc_auto", mol, unwrap, executor);
  const auto auto_values =
      read_msd_values(fs::temp_directory_path() / "titan_msd_smoke_pbc_auto" /
                      "msd_timeseries.dat");
  require_close(auto_values[2], 36.0, "auto unwrap keeps long displacement");

  MsdArgs mic = unwrap;
  mic.pbc = "mic";
  run_analyzer_case("pbc_mic", mol, mic, executor);
  const auto mic_values =
      read_msd_values(fs::temp_directory_path() / "titan_msd_smoke_pbc_mic" /
                      "msd_timeseries.dat");
  require_close(mic_values[2], 16.0, "mic uses minimum-image long lag");
}

void run_com_straddling_pbc_contract(TitanExecutor &executor) {
  Multi_MOL mol =
      make_two_atom_mol({9.9, 0.9}, {0.1, 1.1}, 10.0, true);
  MsdArgs msd;
  msd.flag_msd = true;
  msd.selection = "mol1";
  msd.mode = "com";
  msd.maxlag = "1";
  msd.pbc = "auto";
  run_analyzer_case("com_straddling_pbc", mol, msd, executor);
  const auto values =
      read_msd_values(fs::temp_directory_path() /
                      "titan_msd_smoke_com_straddling_pbc" /
                      "msd_timeseries.dat");
  require_close(values[0], 1.0,
                "COM reconstruction handles a body straddling the box edge");
}

void run_nonincreasing_timestamp_reject_contract(TitanExecutor &executor) {
  Multi_MOL mol = make_two_atom_mol({0.0, 1.0}, {2.0, 3.0});
  mol.m_mol[1].set_time(0.0);
  MsdArgs msd;
  msd.flag_msd = true;
  bool rejected = false;
  try {
    run_analyzer_case("nonincreasing_timestamp", mol, msd, executor);
  } catch (const std::runtime_error &) {
    rejected = true;
  }
  require(rejected, "non-increasing frame timestamps fail closed");
}

} // namespace

int main() {
  TitanExecutor executor(4);
  run_constant_velocity_contract(executor);
  run_atom_mode_contract(executor);
  run_pbc_policy_contract(executor);
  run_com_straddling_pbc_contract(executor);
  run_nonincreasing_timestamp_reject_contract(executor);
  std::cout << "PASS titan_msd_smoke\n";
  return 0;
}
