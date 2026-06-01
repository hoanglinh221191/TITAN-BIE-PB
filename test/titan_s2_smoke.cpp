#include "titan_s2.h"

#include "info.h"
#include "titan_executor.h"

#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

std::string arguments::set_Buildnumber() { return "titan_s2_smoke"; }
std::string arguments::get_build_number() const { return "titan_s2_smoke"; }

bool Multi_MOL::topology_matches_frame(std::size_t frame_index) const noexcept {
  return frame_index < m_mol.size() &&
         m_mol[frame_index].topology_consistent(info.size());
}

namespace {

namespace fs = std::filesystem;

void require(bool condition, const std::string &message) {
  if (!condition) {
    std::cerr << "FAIL titan_s2_smoke: " << message << '\n';
    std::exit(1);
  }
}

void require_close(double actual, double expected, const std::string &message,
                   double tolerance = 1.0e-6) {
  if (!std::isfinite(actual) || std::fabs(actual - expected) > tolerance) {
    std::cerr << "FAIL titan_s2_smoke: " << message << " actual=" << actual
              << " expected=" << expected << '\n';
    std::exit(1);
  }
}

template <typename Fn>
void require_runtime_error(Fn &&fn, const std::string &message) {
  bool rejected = false;
  try {
    fn();
  } catch (const std::runtime_error &) {
    rejected = true;
  }
  require(rejected, message);
}

box make_box(double side = 20.0) {
  box frame_box;
  frame_box.set_a(side);
  frame_box.set_b(side);
  frame_box.set_c(side);
  frame_box.set_alpha(90.0);
  frame_box.set_beta(90.0);
  frame_box.set_gamma(90.0);
  return frame_box;
}

mol_info make_atom(std::size_t resnr, const std::string &resname,
                   const std::string &atom) {
  mol_info result(resnr, atom, atom, "A", resname);
  result.set_chainnum(0);
  result.set_respos(resnr);
  result.set_internal_resnr(resnr);
  return result;
}

void append_xyz(MOL &frame, double x, double y, double z) {
  frame.pos.x.push_back(x);
  frame.pos.y.push_back(y);
  frame.pos.z.push_back(z);
}

Multi_MOL make_resolver_fixture() {
  Multi_MOL mol;
  mol.info = {
      make_atom(1, "ALA", "N"),  make_atom(1, "ALA", "H"),
      make_atom(1, "ALA", "HN"), make_atom(1, "ALA", "CA"),
      make_atom(1, "ALA", "C"),  make_atom(2, "GLY", "N"),
      make_atom(2, "GLY", "HN"), make_atom(2, "GLY", "CA"),
      make_atom(2, "GLY", "C"),  make_atom(3, "PRO", "N"),
      make_atom(3, "PRO", "H"),  make_atom(3, "PRO", "CA"),
      make_atom(3, "PRO", "C"),  make_atom(4, "SER", "N"),
      make_atom(4, "SER", "H1"), make_atom(4, "SER", "H2"),
      make_atom(4, "SER", "H3"), make_atom(4, "SER", "CA"),
      make_atom(4, "SER", "C")};
  MOL frame(make_box(), 0.0, "AMBER99SB");
  for (std::size_t i = 0; i < mol.info.size(); ++i)
    append_xyz(frame, static_cast<double>(i), 0.0, 0.0);
  mol.parameters.resize(mol.info.size());
  mol.m_mol.push_back(std::move(frame));
  return mol;
}

Multi_MOL make_dynamic_fixture(bool rotated_body = false) {
  Multi_MOL mol;
  mol.info = {
      make_atom(1, "ALA", "N"),  make_atom(1, "ALA", "H"),
      make_atom(1, "ALA", "CA"), make_atom(1, "ALA", "C"),
      make_atom(2, "GLY", "N"),  make_atom(2, "GLY", "HN"),
      make_atom(2, "GLY", "CA"), make_atom(2, "GLY", "C"),
      make_atom(3, "SER", "N"),  make_atom(3, "SER", "H"),
      make_atom(3, "SER", "CA"), make_atom(3, "SER", "C")};
  const std::vector<std::array<double, 3>> base{
      {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 1.0, 0.0},
      {4.0, 0.0, 0.0}, {4.0, 1.0, 0.0}, {4.0, 0.0, 0.0}, {4.0, 0.0, 1.0},
      {0.0, 4.0, 0.0}, {0.0, 4.0, 1.0}, {0.0, 4.0, 0.0}, {1.0, 4.0, 0.0}};
  for (std::size_t frame_idx = 0; frame_idx < 4; ++frame_idx) {
    MOL frame(make_box(), static_cast<double>(frame_idx), "AMBER99SB");
    for (const auto &point : base) {
      double x = point[0];
      double y = point[1];
      double z = point[2];
      if (rotated_body && frame_idx > 0) {
        const double old_x = x;
        x = -y + 7.0;
        y = old_x + 3.0;
      }
      append_xyz(frame, x, y, z);
    }
    mol.m_mol.push_back(std::move(frame));
  }
  mol.parameters.resize(mol.info.size());
  return mol;
}

std::string read_file(const fs::path &path) {
  std::ifstream input(path);
  std::ostringstream buffer;
  buffer << input.rdbuf();
  return buffer.str();
}

bool contains(const std::string &haystack, const std::string &needle) {
  return haystack.find(needle) != std::string::npos;
}

fs::path new_output_dir(const std::string &case_name) {
  const fs::path path =
      fs::temp_directory_path() / ("titan_s2_smoke_" + case_name);
  fs::remove_all(path);
  fs::create_directories(path);
  return path;
}

TitanS2OutputFiles trajectory_files() {
  return {};
}

void run_resolver_contract() {
  Multi_MOL mol = make_resolver_fixture();
  TitanS2Analyzer analyzer;
  TitanS2Config config;
  config.flag_s2 = true;
  config.align = "none";
  analyzer.configure(config);
  const fs::path out_dir = new_output_dir("resolver");
  std::fstream log((out_dir / "run.log").string(),
                   std::ios::out | std::ios::trunc);
  analyzer.init(mol, config, log);
  const auto &plans = analyzer.plans();
  require(plans.size() == 6, "default both emits independent NH and CA-C rows");
  require(plans[0].kind == "NH" && plans[0].from_atom == 0 &&
              plans[0].to_atom == 1,
          "NH resolver prefers H over HN");
  require(plans[2].kind == "NH" && plans[2].from_atom == 5 &&
              plans[2].to_atom == 6,
          "NH resolver falls back to HN");
  require(plans[4].resnr == 3 && plans[4].kind == "CA-C",
          "PRO keeps CA-C but excludes NH");
  require(plans[5].resnr == 4 && plans[5].kind == "CA-C",
          "terminal H1/H2/H3 do not create NH");
  TitanExecutor executor(1);
  analyzer.process_frame(mol, 0, executor);
  analyzer.finalize_trajectory(out_dir.string(), trajectory_files(), false,
                               executor, log);
  const std::string summary = read_file(out_dir / "s2_summary.dat");
  require(contains(summary, "missing_nh_skipped: 2") &&
              contains(summary, "pro_nh_skipped: 1") &&
              contains(summary, "warning: skipped 2 NH vector(s)"),
          "summary counts and warns for unavailable NH vectors");
}

void run_selector_contract() {
  Multi_MOL mol = make_dynamic_fixture();
  TitanS2Analyzer analyzer;
  TitanS2Config config;
  config.flag_s2 = true;
  config.residues = "A:2-3";
  config.vectors = "nh";
  analyzer.configure(config);
  std::fstream log;
  analyzer.init(mol, config, log);
  require(analyzer.plans().size() == 2, "selector filters vector plans");
  require(analyzer.plans()[0].resnr == 2 && analyzer.plans()[1].resnr == 3,
          "selector supports the dihedral A:begin-end policy");
}

void run_alignment_contract(TitanExecutor &executor) {
  Multi_MOL mol = make_dynamic_fixture(true);
  std::fstream log;

  TitanS2Config none_config;
  none_config.flag_s2 = true;
  none_config.vectors = "nh";
  none_config.align = "none";
  TitanS2Analyzer none;
  none.configure(none_config);
  none.init(mol, none_config, log);
  for (std::size_t frame = 0; frame < mol.m_mol.size(); ++frame)
    none.process_frame(mol, frame, executor);
  const auto none_rows = none.compute_order_parameters(executor);

  TitanS2Config rms_config = none_config;
  rms_config.align = "rms";
  TitanS2Analyzer rms;
  rms.configure(rms_config);
  rms.init(mol, rms_config, log);
  for (std::size_t frame = 0; frame < mol.m_mol.size(); ++frame)
    rms.process_frame(mol, frame, executor);
  const auto rms_rows = rms.compute_order_parameters(executor);

  require_close(rms_rows.front().s2, 1.0,
                "RMS alignment removes rigid-body rotation");
  require(std::isnan(rms_rows.front().fit_r2),
          "short aligned trajectory reports nan R2 for tail fallback");
  require(none_rows.front().s2 < 0.9,
          "none alignment preserves rigid-body rotation signal");
}

void run_fit_window_and_missing_box_contract(TitanExecutor &executor) {
  Multi_MOL mol = make_dynamic_fixture();
  mol.m_mol[2].Box = box{};
  const fs::path out_dir = new_output_dir("policy");
  std::fstream log((out_dir / "run.log").string(),
                   std::ios::out | std::ios::trunc);
  TitanS2Config config;
  config.flag_s2 = true;
  config.vectors = "nh";
  config.align = "none";
  config.maxlag = "2";
  config.fit_window = "1,2";
  TitanS2Analyzer analyzer;
  analyzer.configure(config);
  analyzer.init(mol, config, log);
  for (std::size_t frame = 0; frame < mol.m_mol.size(); ++frame)
    analyzer.process_frame(mol, frame, executor);
  const auto rows = analyzer.compute_order_parameters(executor);
  require(rows.front().correlation_rows == 2 && rows.front().fit_rows == 2,
          "explicit maxlag and inclusive fit window filter before kernel fit");
  analyzer.finalize_trajectory(out_dir.string(), trajectory_files(), false,
                               executor, log);
  const std::string summary = read_file(out_dir / "s2_summary.dat");
  require(contains(summary, "raw_box_fallback_frames: 1") &&
              contains(summary, "warning: raw bond-vector fallback"),
          "invalid box uses raw vectors and increments summary counter");
}

void run_production_adapter_and_single_sem_contract(TitanExecutor &executor) {
  const fs::path out_dir = new_output_dir("adapter");
  std::fstream log((out_dir / "run.log").string(),
                   std::ios::out | std::ios::trunc);
  Multi_MOL mol = make_dynamic_fixture();
  arguments args;
  args.s2.flag_s2 = true;
  args.s2.align = "none";
  args.s2.bfactor_pdb = true;
  args.flag_separation_results = true;
  args.flag_listfile = true;
  info output_info(args);
  TitanS2Analyzer analyzer;
  analyzer.configure(args.s2);
  analyzer.init(mol, args, log);
  for (std::size_t frame = 0; frame < mol.m_mol.size(); ++frame)
    analyzer.process_frame(mol, frame, executor);
  analyzer.finalize_trajectory(out_dir.string(), output_info, args,
                               {"synthetic.pdb"}, 0, executor, log);
  require(fs::exists(out_dir / "s2_correlation_functions.dat"),
          "default correlation filename follows production contract");
  require(fs::exists(out_dir / "s2_bfactor.pdb"),
          "S2Args bfactor_pdb opts into trajectory PDB");
  require(fs::exists(out_dir / "s2_bfactor_average.pdb"),
          "S2Args bfactor_pdb opts into aggregate PDB");
  const std::string average =
      read_file(out_dir / "s2_order_parameters_average.csv");
  require(contains(average, "valid_replicates") &&
              contains(average, ",nan,") &&
              contains(average, ",nan,1"),
          "single valid replicate exports nan SEM for raw plateau and clipped S2");
}

void run_structural_readiness_contract(TitanExecutor &executor) {
  TitanS2Config config;
  config.flag_s2 = true;
  config.align = "none";
  Multi_MOL parameterless = make_dynamic_fixture();
  parameterless.parameters.clear();
  TitanS2Analyzer analyzer;
  analyzer.configure(config);
  std::fstream log;
  analyzer.init(parameterless, config, log);
  analyzer.process_frame(parameterless, 0, executor);

  Multi_MOL mismatched = make_dynamic_fixture();
  mismatched.m_mol[0].pos.x.pop_back();
  TitanS2Analyzer mismatch_analyzer;
  mismatch_analyzer.configure(config);
  mismatch_analyzer.init(mismatched, config, log);
  require_runtime_error(
      [&]() { mismatch_analyzer.process_frame(mismatched, 0, executor); },
      "analyzer rejects topology/frame-size mismatch");
}

void run_output_and_aggregate_contract(TitanExecutor &executor) {
  const fs::path out_dir = new_output_dir("output");
  std::fstream log((out_dir / "run.log").string(),
                   std::ios::out | std::ios::trunc);
  TitanS2Config config;
  config.flag_s2 = true;
  config.pdb = true;
  config.average_pdb = true;
  config.align = "none";
  TitanS2Analyzer analyzer;
  analyzer.configure(config);

  for (int replicate = 0; replicate < 2; ++replicate) {
    Multi_MOL mol = make_dynamic_fixture();
    if (replicate == 1) {
      for (std::size_t frame = 1; frame < mol.m_mol.size(); ++frame)
        mol.m_mol[frame].pos.y[1] += static_cast<double>(frame);
    }
    analyzer.init(mol, config, log);
    for (std::size_t frame = 0; frame < mol.m_mol.size(); ++frame)
      analyzer.process_frame(mol, frame, executor);
    analyzer.finalize_trajectory(out_dir.string(), trajectory_files(), true,
                                 executor, log);
  }

  const std::string csv = read_file(out_dir / "s2_order_parameters.csv");
  const std::string summary = read_file(out_dir / "s2_summary.dat");
  const std::string pdb = read_file(out_dir / "s2_bfactor.pdb");
  const std::string average =
      read_file(out_dir / "s2_order_parameters_average.csv");
  const std::string average_summary =
      read_file(out_dir / "s2_summary_average.dat");
  const std::string average_pdb =
      read_file(out_dir / "s2_bfactor_average.pdb");
  const std::string correlation =
      read_file(out_dir / "s2_correlation_functions.dat");

  require(contains(csv,
                   "residue,chain,resnr,resname,vector_type,plateau_C_raw,"
                   "S2_clipped,tau_ps,fit_R2,fit_window_start_ps,"
                   "fit_window_end_ps,samples,fit_status,weighted_residual"),
          "per-trajectory order CSV exports the promised fit schema");
  require(contains(summary, "reference_source: first_analyzed_frame") &&
              contains(summary,
                       "reference_note: Bond vectors are evaluated after Kabsch "
                       "alignment to the first analyzed mol1 frame when "
                       "s2_align=rms. Interpret S2 with fit diagnostics and "
                       "trajectory length.") &&
              contains(summary,
                       "bfactor_note: s2_bfactor.pdb stores clipped S2 in the "
                       "B-factor field for visualization; NH is preferred and "
                       "CA-C is used only when NH is unavailable.") &&
              contains(summary, "missing_nh_skipped: 0") &&
              contains(summary, "clipped_s2_count:"),
          "summary exports exact reference lines and science diagnostics");
  require(contains(correlation,
                   "residue chain resnr resname vector_type lag_index "
                   "lag_time_ps correlation samples"),
          "correlation DAT exports separate residue metadata columns");
  require(contains(pdb, "ATOM") && contains(pdb, "  1.00"),
          "PDB maps NH-preferred residue S2 onto first snapshot atoms");
  require(contains(pdb, "ATOM      1    N ALA A   1") &&
              contains(pdb, "  0.37"),
          "trajectory PDB prefers NH S2 over the residue CA-C S2");
  require(contains(average,
                   "mean_plateau_C_raw,sem_plateau_C_raw,mean_S2_clipped,"
                   "sem_S2_clipped,valid_replicates") &&
              contains(average,
                       "A:1:ALA,A,1,ALA,NH,0.68500000,0.31500000,"
                       "0.68500000,0.31500000,2") &&
              contains(average_summary, "replicates: 2"),
          "aggregate CSV and summary retain replicate mean and SEM across reset");
  require(contains(average_pdb, "ATOM") && contains(average_pdb, "  0.68"),
          "average PDB is emitted when explicitly opted in");

  const fs::path root_dir = new_output_dir("root_average");
  analyzer.finalize_all(root_dir.string(), trajectory_files());
  require(fs::exists(root_dir / "s2_order_parameters_average.csv") &&
              contains(read_file(root_dir / "s2_summary_average.dat"),
                       "replicates: 2"),
          "final root aggregate export is available after separated results");
}

void run_tail_fallback_diagnostics_contract(TitanExecutor &executor) {
  Multi_MOL mol = make_dynamic_fixture();
  TitanS2Config config;
  config.flag_s2 = true;
  config.align = "none";
  config.vectors = "nh";
  TitanS2Analyzer analyzer;
  analyzer.configure(config);
  std::fstream log;
  analyzer.init(mol, config, log);
  for (std::size_t frame = 0; frame < mol.m_mol.size(); ++frame)
    analyzer.process_frame(mol, frame, executor);
  const auto rows = analyzer.compute_order_parameters(executor);
  require(!rows.empty() && rows.front().fit_ok &&
              rows.front().used_tail_fallback,
          "short trajectory uses explicit tail fallback");
  require(std::isnan(rows.front().tau_ps) && std::isnan(rows.front().fit_r2) &&
              std::isnan(rows.front().weighted_residual),
          "tail fallback does not export fake model diagnostics");
}

void run_empty_and_mismatched_aggregate_contract(TitanExecutor &executor) {
  const fs::path out_dir = new_output_dir("aggregate_guards");
  TitanS2Config config;
  config.flag_s2 = true;
  config.align = "none";
  config.vectors = "nh";
  TitanS2Analyzer analyzer;
  analyzer.configure(config);
  std::fstream log;

  Multi_MOL empty = make_dynamic_fixture();
  analyzer.init(empty, config, log);
  analyzer.finalize_trajectory(out_dir.string(), trajectory_files(), true,
                               executor, log);
  require(!fs::exists(out_dir / "s2_summary_average.dat"),
          "empty trajectory does not create a replicate aggregate");

  Multi_MOL first = make_dynamic_fixture();
  analyzer.init(first, config, log);
  for (std::size_t frame = 0; frame < first.m_mol.size(); ++frame)
    analyzer.process_frame(first, frame, executor);
  analyzer.finalize_trajectory(out_dir.string(), trajectory_files(), true,
                               executor, log);
  require(contains(read_file(out_dir / "s2_summary_average.dat"),
                   "replicates: 1"),
          "first valid trajectory creates one aggregate replicate");

  Multi_MOL mismatched = make_dynamic_fixture();
  mismatched.info.erase(mismatched.info.begin() + 1);
  for (MOL &frame : mismatched.m_mol) {
    frame.pos.x.erase(frame.pos.x.begin() + 1);
    frame.pos.y.erase(frame.pos.y.begin() + 1);
    frame.pos.z.erase(frame.pos.z.begin() + 1);
  }
  analyzer.init(mismatched, config, log);
  for (std::size_t frame = 0; frame < mismatched.m_mol.size(); ++frame)
    analyzer.process_frame(mismatched, frame, executor);
  require_runtime_error(
      [&]() {
        analyzer.finalize_trajectory(out_dir.string(), trajectory_files(),
                                     true, executor, log);
      },
      "aggregate rejects non-identical vector plans across trajectories");
}

} // namespace

int main() {
  TitanExecutor executor(4);
  run_resolver_contract();
  run_selector_contract();
  run_alignment_contract(executor);
  run_fit_window_and_missing_box_contract(executor);
  run_production_adapter_and_single_sem_contract(executor);
  run_structural_readiness_contract(executor);
  run_output_and_aggregate_contract(executor);
  run_tail_fallback_diagnostics_contract(executor);
  run_empty_and_mismatched_aggregate_contract(executor);
  std::cout << "PASS titan_s2_smoke\n";
  return 0;
}
