#include "titan_dihedral.h"

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

std::string arguments::set_Buildnumber() { return "titan_dihedral_smoke"; }
std::string arguments::get_build_number() const {
  return "titan_dihedral_smoke";
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

Multi_MOL make_peptide(const std::vector<std::vector<double>> &xs,
                       const std::vector<std::vector<double>> &ys,
                       const std::vector<std::vector<double>> &zs) {
  Multi_MOL mol;
  mol.info = {
      make_atom(1, "ALA", "N"),  make_atom(1, "ALA", "CA"),
      make_atom(1, "ALA", "C"),  make_atom(2, "GLY", "N"),
      make_atom(2, "GLY", "CA"), make_atom(2, "GLY", "C"),
      make_atom(3, "SER", "N"),  make_atom(3, "SER", "CA"),
      make_atom(3, "SER", "C"),  make_atom(4, "THR", "N"),
      make_atom(4, "THR", "CA"), make_atom(4, "THR", "C")};

  for (std::size_t frame_idx = 0; frame_idx < xs.size(); ++frame_idx) {
    MOL frame(make_box(), static_cast<double>(frame_idx), "AMBER99SB");
    frame.pos.x.assign(xs[frame_idx].begin(), xs[frame_idx].end());
    frame.pos.y.assign(ys[frame_idx].begin(), ys[frame_idx].end());
    frame.pos.z.assign(zs[frame_idx].begin(), zs[frame_idx].end());
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
    std::cerr << "FAIL titan_dihedral_smoke: " << message << '\n';
    std::exit(1);
  }
}

void require_close(double actual, double expected, const std::string &message) {
  if (std::fabs(actual - expected) > 1.0e-6) {
    std::cerr << "FAIL titan_dihedral_smoke: " << message << " actual="
              << actual << " expected=" << expected << '\n';
    std::exit(1);
  }
}

fs::path new_output_dir(const std::string &case_name) {
  const fs::path out_dir =
      fs::temp_directory_path() / ("titan_dihedral_smoke_" + case_name);
  fs::remove_all(out_dir);
  fs::create_directories(out_dir);
  return out_dir;
}

void run_geometry_contract() {
  require_close(TitanDihedralAnalyzer::compute_dihedral_deg(
                    {1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 1.0, 0.0},
                    {0.0, 1.0, 1.0}),
                90.0, "canonical +90 degree torsion");
  require_close(TitanDihedralAnalyzer::compute_dihedral_deg(
                    {1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 1.0, 0.0},
                    {0.0, 1.0, -1.0}),
                -90.0, "canonical -90 degree torsion");
  require_close(TitanDihedralAnalyzer::compute_dihedral_deg(
                    {1.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 1.0, 0.0},
                    {1.0, 1.0, 0.0}),
                180.0, "canonical trans torsion");
}

void run_analyzer_contract(TitanExecutor &executor) {
  arguments args;
  args.dihedral.flag_dihedral = true;
  args.dihedral.kind = "backbone";
  args.dihedral.residues = "A:2-3";
  args.dihedral.stride = 2;
  args.dihedral.export_tica = true;

  const std::vector<double> x0{0.0, 1.0, 1.0, 0.0, 0.0, 0.0,
                              1.0, 1.0, 2.0, 2.0, 2.0, 3.0};
  const std::vector<double> y0{0.0, 0.0, 1.0, 0.0, 1.0, 1.0,
                              1.0, 2.0, 2.0, 2.0, 3.0, 3.0};
  const std::vector<double> z0{0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
                              1.0, 1.0, 2.0, 2.0, 2.0, 3.0};
  std::vector<double> x1 = x0;
  x1[5] += 0.5;
  Multi_MOL mol = make_peptide({x0, x1}, {y0, y0}, {z0, z0});

  const fs::path out_dir = new_output_dir("analyzer");
  std::fstream log((out_dir / "run.log").string(),
                   std::ios::out | std::ios::trunc);

  TitanDihedralAnalyzer analyzer;
  analyzer.configure(args.dihedral);
  analyzer.init(mol, args, log);
  require(analyzer.enabled(), "dihedral analyzer should be enabled");
  for (std::size_t frame = 0; frame < mol.m_mol.size(); ++frame)
    analyzer.process_frame(mol, frame, executor);

  info output_info(args);
  const std::vector<std::string> input_line{"dihedral-smoke.pdb"};
  analyzer.finalize_trajectory(out_dir.string(), output_info, args, input_line,
                               0, log);

  const std::string timeseries =
      read_file(out_dir / "dihedral_backbone_timeseries.csv");
  const std::string rama = read_file(out_dir / "dihedral_ramachandran.csv");
  const std::string summary = read_file(out_dir / "dihedral_summary.dat");
  const std::string features =
      read_file(out_dir / "dihedral_feature_matrix.csv");

  require(contains(timeseries,
                   "frame,time_ps,residue,chain,resnr,resname,phi_deg,psi_deg,omega_deg"),
          "timeseries header is exported");
  require(contains(timeseries, "0,0.000000,A:2:GLY"),
          "selected residue 2 is exported");
  require(!contains(timeseries, "1,1.000000"),
          "stride=2 skips frame 1");
  require(contains(rama,
                   "residue,chain,resnr,resname,phi_bin_deg,psi_bin_deg,count,occupancy,valid_samples"),
          "Ramachandran header is exported");
  require(contains(summary, "kind: backbone") &&
              contains(summary, "stride: 2") &&
              contains(summary, "residue_selection: A:2-3"),
          "summary records backbone config");
  require(contains(features, "A:2:GLY_phi_sin") &&
              contains(features, "A:3:SER_omega_cos"),
          "TICA feature matrix exports sin/cos backbone columns");
}

void run_stride_uses_trajectory_ordinal(TitanExecutor &executor) {
  arguments args;
  args.dihedral.flag_dihedral = true;
  args.dihedral.stride = 2;
  args.dihedral.export_tica = true;

  const std::vector<double> x0{0.0, 1.0, 1.0, 0.0, 0.0, 0.0,
                              1.0, 1.0, 2.0, 2.0, 2.0, 3.0};
  const std::vector<double> y0{0.0, 0.0, 1.0, 0.0, 1.0, 1.0,
                              1.0, 2.0, 2.0, 2.0, 3.0, 3.0};
  const std::vector<double> z0{0.0, 0.0, 0.0, 0.0, 0.0, 1.0,
                              1.0, 1.0, 2.0, 2.0, 2.0, 3.0};
  Multi_MOL mol = make_peptide({x0}, {y0}, {z0});

  const fs::path out_dir = new_output_dir("stride_ordinal");
  std::fstream log((out_dir / "run.log").string(),
                   std::ios::out | std::ios::trunc);

  TitanDihedralAnalyzer analyzer;
  analyzer.configure(args.dihedral);
  analyzer.init(mol, args, log);
  analyzer.process_frame(mol, 0, executor);
  analyzer.process_frame(mol, 0, executor);
  analyzer.process_frame(mol, 0, executor);

  info output_info(args);
  const std::vector<std::string> input_line{"dihedral-stride.pdb"};
  analyzer.finalize_trajectory(out_dir.string(), output_info, args, input_line,
                               0, log);

  const std::string summary = read_file(out_dir / "dihedral_summary.dat");
  const std::string timeseries =
      read_file(out_dir / "dihedral_backbone_timeseries.csv");
  require(contains(summary, "frames_seen: 3") &&
              contains(summary, "frames_used: 2"),
          "stride uses trajectory ordinal even when snap_idx repeats");
  require(contains(timeseries, "0,0.000000") &&
              contains(timeseries, "2,0.000000") &&
              !contains(timeseries, "1,0.000000"),
          "exported frame ids reflect trajectory ordinal under chunk-local snap_idx");
}

} // namespace

int main() {
  TitanExecutor executor(4);
  run_geometry_contract();
  run_analyzer_contract(executor);
  run_stride_uses_trajectory_ordinal(executor);
  std::cout << "PASS titan_dihedral_smoke\n";
  return 0;
}
