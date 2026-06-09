#include "titan_membrane.h"

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

std::string arguments::set_Buildnumber() { return "titan_membrane_smoke"; }
std::string arguments::get_build_number() const { return "titan_membrane_smoke"; }

namespace {

namespace fs = std::filesystem;

void require(bool condition, const std::string &message) {
  if (!condition) {
    std::cerr << "FAIL titan_membrane_smoke: " << message << '\n';
    std::exit(1);
  }
}

void require_close(double actual, double expected, const std::string &message,
                   double tolerance = 1.0e-6) {
  if (std::fabs(actual - expected) > tolerance) {
    std::cerr << "FAIL titan_membrane_smoke: " << message
              << " actual=" << actual << " expected=" << expected << '\n';
    std::exit(1);
  }
}

box make_box(double a = 20.0, double b = 20.0, double c = 40.0) {
  box frame_box;
  frame_box.set_a(a);
  frame_box.set_b(b);
  frame_box.set_c(c);
  frame_box.set_alpha(90.0);
  frame_box.set_beta(90.0);
  frame_box.set_gamma(90.0);
  return frame_box;
}

mol_info make_atom(std::size_t resnr, const std::string &name) {
  mol_info atom(resnr, name, name, "L", "DOPC");
  atom.set_respos(resnr);
  atom.set_internal_resnr(resnr);
  return atom;
}

mol_info make_protein_atom(std::size_t resnr, const std::string &name,
                           const std::string &chain = "A",
                           const std::string &resname = "LEU") {
  mol_info atom(resnr, name, name, chain, resname);
  atom.set_respos(resnr);
  atom.set_internal_resnr(resnr);
  return atom;
}

void append_lipid(Multi_MOL &lipids, MOL &frame, std::size_t resnr, double x,
                  double y, double head_z, bool upper) {
  lipids.info.push_back(make_atom(resnr, "PO4"));
  lipids.info.push_back(make_atom(resnr, "C1A"));
  lipids.info.push_back(make_atom(resnr, "C2A"));
  const double direction = upper ? -1.0 : 1.0;
  frame.pos.x.insert(frame.pos.x.end(), {x, x, x});
  frame.pos.y.insert(frame.pos.y.end(), {y, y, y});
  frame.pos.z.insert(frame.pos.z.end(),
                     {head_z, head_z + direction, head_z + 2.0 * direction});
}

Multi_MOL make_bilayer() {
  Multi_MOL lipids;
  const std::vector<std::pair<double, double>> xy{
      {5.0, 5.0}, {15.0, 5.0}, {5.0, 15.0}, {15.0, 15.0}};
  for (std::size_t frame_idx = 0; frame_idx < 3; ++frame_idx) {
    MOL frame(make_box(), static_cast<double>(frame_idx), "MARTINI");
    for (std::size_t i = 0; i < xy.size(); ++i) {
      double upper_z = 5.0;
      bool upper = true;
      if (i == 0 && frame_idx > 0) {
        upper_z = -5.0;
        upper = false;
      }
      append_lipid(lipids, frame, i + 1, xy[i].first, xy[i].second, upper_z,
                   upper);
    }
    for (std::size_t i = 0; i < xy.size(); ++i)
      append_lipid(lipids, frame, i + 5, xy[i].first, xy[i].second, -5.0,
                   false);
    lipids.m_mol.push_back(std::move(frame));
    if (frame_idx == 0)
      lipids.parameters.resize(lipids.info.size());
    else
      lipids.info.resize(lipids.parameters.size());
  }
  return lipids;
}

Multi_MOL make_nearby_protein() {
  Multi_MOL protein;
  MOL frame(make_box(), 0.0, "MARTINI");
  for (const std::string &name : {"N", "CA", "CB", "CG"})
    protein.info.push_back(make_protein_atom(1, name));
  frame.pos.x.insert(frame.pos.x.end(), {5.0, 5.2, 5.4, 5.4});
  frame.pos.y.insert(frame.pos.y.end(), {5.0, 5.2, 5.4, 5.4});
  frame.pos.z.insert(frame.pos.z.end(), {6.0, 5.8, 5.6, 5.4});
  protein.parameters.resize(protein.info.size());
  protein.m_mol.push_back(std::move(frame));
  return protein;
}

Multi_MOL make_undulating_bilayer() {
  Multi_MOL lipids;
  MOL frame(make_box(48.0, 48.0, 40.0), 0.0, "MARTINI");
  std::size_t resnr = 1;
  for (double y : {6.0, 18.0, 30.0, 42.0}) {
    for (double x : {6.0, 18.0, 30.0, 42.0}) {
      const double midplane = 0.04 * (x - 24.0) + 0.02 * (y - 24.0);
      append_lipid(lipids, frame, resnr++, x, y, midplane + 6.0, true);
    }
  }
  for (double y : {6.0, 18.0, 30.0, 42.0}) {
    for (double x : {6.0, 18.0, 30.0, 42.0}) {
      const double midplane = 0.04 * (x - 24.0) + 0.02 * (y - 24.0);
      append_lipid(lipids, frame, resnr++, x, y, midplane - 6.0, false);
    }
  }
  lipids.parameters.resize(lipids.info.size());
  lipids.m_mol.push_back(std::move(frame));
  return lipids;
}

std::string read_file(const fs::path &path) {
  std::ifstream input(path);
  std::ostringstream buffer;
  buffer << input.rdbuf();
  return buffer.str();
}

void run_core_contract() {
  const fs::path out_dir = fs::temp_directory_path() / "titan_membrane_smoke";
  fs::remove_all(out_dir);
  fs::create_directories(out_dir);
  std::fstream log((out_dir / "run.log").string(),
                   std::ios::out | std::ios::trunc);
  Multi_MOL lipids = make_bilayer();
  arguments args;
  args.membrane.flag_membrane = true;
  args.membrane.normal = "z";
  args.membrane.grid = "2,2";
  args.membrane.flipflop_persistence = 2;
  TitanExecutor executor(4);
  TitanMembraneAnalyzer analyzer;
  analyzer.configure(args.membrane);
  analyzer.init(lipids, log);
  require(analyzer.enabled(), "analyzer is enabled after init");
  require(analyzer.lipid_count() == 8, "all DOPC residues are cached");
  for (std::size_t frame = 0; frame < lipids.m_mol.size(); ++frame)
    analyzer.process_frame(lipids, frame, executor);
  require(analyzer.frames().size() == 3, "all synthetic frames are analyzed");
  require_close(analyzer.frames()[0].mean_apl_a2, 100.0,
                "periodic Voronoi APL conserves leaflet area");
  require_close(analyzer.frames()[0].mean_thickness_a, 10.0,
                "flat bilayer thickness follows headgroup separation");
  require_close(TitanMembraneAnalyzer::order_parameter(0.0, 0.0, 1.0, 0.0,
                                                       0.0, 1.0),
                1.0, "parallel tail vector has S=1");
  require_close(TitanMembraneAnalyzer::order_parameter(1.0, 0.0, 0.0, 0.0,
                                                       0.0, 1.0),
                -0.5, "perpendicular tail vector has S=-0.5");
  info output_info(args);
  analyzer.finalize_trajectory(out_dir.string(), output_info, args,
                               {"synthetic.pdb"}, 0, log);
  require(fs::exists(out_dir / "membrane_leaflets.csv"),
          "leaflet output is written");
  require(fs::exists(out_dir / "membrane_apl_timeseries.csv"),
          "APL timeseries output is written");
  require(fs::exists(out_dir / "membrane_apl_per_lipid.csv"),
          "per-lipid APL output is written");
  require(fs::exists(out_dir / "membrane_thickness_timeseries.csv"),
          "thickness timeseries output is written");
  require(fs::exists(out_dir / "membrane_thickness_map.csv"),
          "thickness map output is written");
  require(fs::exists(out_dir / "membrane_order_parameters.csv"),
          "order parameter output is written");
  require(fs::exists(out_dir / "membrane_interdigitation.dat"),
          "interdigitation output is written");
  require(fs::exists(out_dir / "membrane_flipflop_events.csv"),
          "flip-flop output is written");
  const std::string summary = read_file(out_dir / "membrane_summary.dat");
  require(summary.find("reference_note:") != std::string::npos,
          "summary tells users how geometry is referenced");
  const std::string flipflop =
      read_file(out_dir / "membrane_flipflop_events.csv");
  require(flipflop.find("L:1:DOPC,upper,lower") != std::string::npos,
          "persistent leaflet change is emitted as flip-flop");
}

void run_auto_local_midplane_contract() {
  const fs::path out_dir =
      fs::temp_directory_path() / "titan_membrane_auto_smoke";
  fs::remove_all(out_dir);
  fs::create_directories(out_dir);
  std::fstream log((out_dir / "run.log").string(),
                   std::ios::out | std::ios::trunc);
  Multi_MOL lipids = make_undulating_bilayer();
  arguments args;
  args.membrane.flag_membrane = true;
  args.membrane.normal = "auto";
  args.membrane.grid = "2,2";
  TitanExecutor executor(4);
  TitanMembraneAnalyzer analyzer;
  analyzer.configure(args.membrane);
  analyzer.init(lipids, log);
  analyzer.process_frame(lipids, 0, executor);
  require(analyzer.frames().size() == 1,
          "auto local-midplane branch analyzes the tilted fixture");
  require(analyzer.frames()[0].upper_count == 16 &&
              analyzer.frames()[0].lower_count == 16,
          "auto local-midplane branch resolves both undulating leaflets");
}

void run_user_leaflet_contract() {
  const fs::path out_dir =
      fs::temp_directory_path() / "titan_membrane_user_leaflet_smoke";
  fs::remove_all(out_dir);
  fs::create_directories(out_dir);
  const fs::path leaflet_file = out_dir / "leaflets.csv";
  {
    std::ofstream out(leaflet_file);
    out << "lipid_id,leaflet\n";
    for (std::size_t resnr = 1; resnr <= 8; ++resnr)
      out << "L:" << resnr << ":DOPC,"
          << (resnr <= 4 ? "upper" : "lower") << '\n';
  }
  std::fstream log((out_dir / "run.log").string(),
                   std::ios::out | std::ios::trunc);
  Multi_MOL lipids = make_bilayer();
  arguments args;
  args.membrane.flag_membrane = true;
  args.membrane.normal = "z";
  args.membrane.grid = "2,2";
  args.membrane.leaflet = "user";
  args.membrane.leaflet_file = leaflet_file.string();
  TitanExecutor executor(4);
  TitanMembraneAnalyzer analyzer;
  analyzer.configure(args.membrane);
  analyzer.init(lipids, log);
  analyzer.process_frame(lipids, 0, executor);
  require(analyzer.frames()[0].upper_count == 4 &&
              analyzer.frames()[0].lower_count == 4,
          "user leaflet CSV controls runtime assignment");
}

void run_proximity_coupling_contract() {
  const fs::path out_dir =
      fs::temp_directory_path() / "titan_membrane_proximity_smoke";
  fs::remove_all(out_dir);
  fs::create_directories(out_dir);
  std::fstream log((out_dir / "run.log").string(),
                   std::ios::out | std::ios::trunc);
  Multi_MOL lipids = make_bilayer();
  Multi_MOL protein = make_nearby_protein();
  arguments args;
  args.membrane.flag_membrane = true;
  args.membrane.normal = "z";
  args.membrane.grid = "2,2";
  args.membrane.protein_contacts = true;
  args.membrane.enrichment = true;
  args.membrane.curvature = "grid";
  TitanExecutor executor(4);
  TitanMembraneAnalyzer analyzer;
  analyzer.configure(args.membrane);
  analyzer.init(lipids, log);
  analyzer.process_frame(lipids, &protein, 0, executor);
  info output_info(args);
  analyzer.finalize_trajectory(out_dir.string(), output_info, args,
                               {"synthetic.pdb"}, 0, log);
  const fs::path coupling_path = out_dir / "membrane_proximity_coupling.csv";
  require(fs::exists(coupling_path),
          "protein-proximity coupling output is written");
  const std::string coupling = read_file(coupling_path);
  require(coupling.find("apl_A2") != std::string::npos,
          "proximity coupling includes APL rows");
  require(coupling.find("thickness_A") != std::string::npos,
          "proximity coupling includes thickness rows");
  require(coupling.find("order_S") != std::string::npos,
          "proximity coupling includes order rows");
  require(coupling.find("curvature_proxy_Ainv") != std::string::npos,
          "proximity coupling includes curvature rows");
  require(coupling.find("lipid_enrichment") != std::string::npos,
          "proximity coupling includes leaflet composition rows");
  const std::string summary =
      read_file(out_dir / "membrane_coupling_summary.dat");
  require(summary.find("proximity_coupling_note: Near/far") !=
              std::string::npos,
          "coupling summary documents near/far policy");
}

} // namespace

int main() {
  run_core_contract();
  run_auto_local_midplane_contract();
  run_user_leaflet_contract();
  run_proximity_coupling_contract();
  std::cout << "PASS titan_membrane_smoke\n";
  return 0;
}
