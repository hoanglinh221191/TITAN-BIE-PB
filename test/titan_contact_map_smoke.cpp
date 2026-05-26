#include "titan_contact_map.h"

#include "info.h"
#include "titan_executor.h"

#include <cstdlib>
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

std::string read_file(const std::string &path) {
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

} // namespace

int main() {
  arguments args;
  args.cmap.flag_cmap = true;
  args.cmap.mode = "sidechain-heavy";
  args.cmap.cutoff_ang = 4.5;
  args.cmap.domain = "intra";
  args.cmap.min_seqsep = 3;

  Multi_MOL mol;
  mol.info.push_back(make_atom(1, "ALA", "CB"));
  mol.info.push_back(make_atom(2, "GLY", "CA"));
  mol.info.push_back(make_atom(3, "VAL", "CB"));
  mol.info.push_back(make_atom(4, "THR", "CB"));
  mol.info.push_back(make_atom(5, "LEU", "CB"));
  mol.info.push_back(make_atom(5, "LEU", "CG"));
  mol.info.push_back(make_atom(9, "SER", "CB"));

  box frame_box;
  frame_box.set_a(30.0);
  frame_box.set_b(30.0);
  frame_box.set_c(30.0);
  frame_box.set_alpha(90.0);
  frame_box.set_beta(90.0);
  frame_box.set_gamma(90.0);

  MOL frame(frame_box, 0.0, "AMBER99SB");
  frame.pos.x = {1.0, 20.0, 20.0, 20.0, 3.0, 3.2, 20.0};
  frame.pos.y = {1.0, 2.0, 4.0, 6.0, 1.0, 1.0, 20.0};
  frame.pos.z = {1.0, 2.0, 4.0, 6.0, 1.0, 1.0, 20.0};
  mol.m_mol.push_back(frame);

  TitanContactMapAnalyzer analyzer;
  analyzer.configure(args.cmap);

  std::fstream log("titan_contact_map_smoke.log",
                   std::ios::out | std::ios::trunc);
  analyzer.init(mol, args, log);
  require(analyzer.enabled(), "analyzer should be enabled");
  require(analyzer.initialized(), "analyzer should initialize");

  TitanExecutor executor(1);
  analyzer.process_frame(mol, 0, executor);

  info output_info(args);
  const std::vector<std::string> input_line{"smoke.pdb"};
  analyzer.finalize_trajectory(".", output_info, args, input_line, 0, log);
  log.close();

  const std::string summary_path =
      output_info.get_cmap_summary_filename(args, input_line, 0);
  const std::string edges_path =
      output_info.get_cmap_edges_filename(args, input_line, 0);
  const std::string matrix_path =
      output_info.get_cmap_matrix_filename(args, input_line, 0);

  const std::string summary = read_file(summary_path);
  const std::string edges = read_file(edges_path);
  const std::string matrix = read_file(matrix_path);

  require(contains(summary, "module: titan_cmap"), "summary module marker");
  require(contains(summary, "mode: sidechain-heavy"), "summary mode");
  require(contains(edges, "ALA"), "edge ALA residue");
  require(contains(edges, "LEU"), "edge LEU residue");
  require(contains(edges, "1.000000"), "edge probability");
  require(contains(matrix, "A:1:ALA"), "matrix ALA label");
  require(contains(matrix, "A:5:LEU"), "matrix LEU label");
  require(contains(matrix, "1.000000"), "matrix probability");

  std::cout << "PASS titan_contact_map_smoke\n";
  return 0;
}
