#include "titan_energy_decomp_engine.h"

#include "arguments.h"
#include "cache_sigma_epsilon.h"
#include "titan_energy_decomposition.h"
#include "titan_executor.h"

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

std::string arguments::set_Buildnumber() {
  return "titan_energy_decomp_engine_smoke";
}
std::string arguments::get_build_number() const {
  return "titan_energy_decomp_engine_smoke";
}

namespace {

namespace edecomp = titan_energy_decomp;
namespace fs = std::filesystem;

void require(bool condition, const std::string &message) {
  if (!condition) {
    std::cerr << "FAIL titan_energy_decomp_engine_smoke: " << message << '\n';
    std::exit(1);
  }
}

void require_close(double actual, double expected, const std::string &message,
                   double tolerance = 1.0e-12) {
  if (!std::isfinite(actual) || std::fabs(actual - expected) > tolerance) {
    std::cerr << "FAIL titan_energy_decomp_engine_smoke: " << message
              << " actual=" << actual << " expected=" << expected << '\n';
    std::exit(1);
  }
}

mol_info make_atom(std::size_t resnr, const std::string &resname,
                   const std::string &atom, const std::string &chain) {
  mol_info out(resnr, atom, atom, chain, resname);
  out.set_chainnum(chain == "A" ? 0 : 1);
  out.set_respos(resnr);
  out.set_internal_resnr(resnr);
  return out;
}

box make_box() {
  box out;
  out.set_a(80.0);
  out.set_b(80.0);
  out.set_c(80.0);
  out.set_alpha(90.0);
  out.set_beta(90.0);
  out.set_gamma(90.0);
  return out;
}

Multi_MOL make_molecule(const std::vector<mol_info> &atoms,
                        const std::vector<double> &xs) {
  Multi_MOL out;
  out.info = atoms;
  out.parameters.resize(atoms.size());
  for (std::size_t i = 0; i < atoms.size(); ++i) {
    out.parameters[i].set_data(0.1 * static_cast<double>(i + 1), 0.1, 0.2,
                               0.4, 1.0);
  }
  MOL frame(make_box(), 0.0, "AMBER99SB");
  frame.pos.x.assign(xs.begin(), xs.end());
  frame.pos.y.assign(xs.size(), 0.0);
  frame.pos.z.assign(xs.size(), 0.0);
  out.m_mol.push_back(std::move(frame));
  return out;
}

cache_sigma_epsilon make_cache(const std::vector<double> &totals) {
  cache_sigma_epsilon out;
  out.elec_ener.resize(totals.size());
  out.vdw_ener.resize(totals.size());
  for (std::size_t i = 0; i < totals.size(); ++i) {
    out.elec_ener[i] = 0.25 * totals[i];
    out.vdw_ener[i] = 0.75 * totals[i];
  }
  return out;
}

edecomp::EngineConfig config(std::size_t top_n = 16) {
  edecomp::EngineConfig out;
  out.top_edge_count = top_n;
  out.edge_metric = edecomp::RankingMetric::MeanAbs;
  out.ranking_metric = edecomp::RankingMetric::SignedMean;
  out.sort_direction = edecomp::SortDirection::Descending;
  return out;
}

const edecomp::PairStatistics &at(const edecomp::DenseMatrix &matrix,
                                  std::size_t row, std::size_t col) {
  return matrix.values.at(row * matrix.cols + col);
}

edecomp::TrajectoryResult process_one(
    const Multi_MOL &mol1, const Multi_MOL &mol2,
    const std::vector<std::vector<double>> &frames, TitanExecutor &executor,
    std::size_t top_n = 16) {
  edecomp::EnergyDecompositionEngine engine(
      edecomp::make_axes(mol1, "mol1"), edecomp::make_axes(mol2, "mol2"),
      config(top_n));
  engine.begin_trajectory("single");
  for (const auto &frame : frames) {
    cache_sigma_epsilon cache = make_cache(frame);
    engine.process_frame({&cache, nullptr, nullptr, nullptr}, executor);
  }
  return engine.finalize_trajectory();
}

void run_cache_residue_atom_ranking_contract() {
  const Multi_MOL mol1 =
      make_molecule({make_atom(1, "ALA", "CA", "A"),
                     make_atom(1, "ALA", "CB", "A"),
                     make_atom(2, "GLY", "CA", "A")},
                    {0.0, 1.0, 2.0});
  const Multi_MOL mol2 =
      make_molecule({make_atom(5, "LEU", "CA", "B"),
                     make_atom(6, "SER", "CA", "B")},
                    {10.0, 20.0});
  TitanExecutor executor(4);
  const edecomp::TrajectoryResult result =
      process_one(mol1, mol2, {{1.0, -2.0, 3.0, 4.0, -5.0, 6.0},
                               {3.0, -2.0, -1.0, 2.0, -1.0, 4.0}},
                  executor, 2);

  require(result.frames == 2, "trajectory frame count");
  require(result.residues.rows == 2 && result.residues.cols == 2,
          "residue matrix dimensions");
  require(result.atoms.rows == 2 && result.atoms.cols == 3,
          "atom matrix dimensions follow mol2 rows and mol1 columns");
  require_close(at(result.residues, 0, 0).signed_mean, 0.0,
                "residue closure row0 col0 signed mean");
  require_close(at(result.residues, 0, 0).mean_abs, 1.0,
                "residue closure row0 col0 mean abs");
  require_close(at(result.residues, 1, 0).signed_mean, 0.0,
                "residue closure row1 col0 signed mean");
  require_close(at(result.residues, 1, 0).abs_mean, 0.0,
                "abs_mean is abs of signed mean");
  require_close(at(result.residues, 1, 0).mean_abs, 1.0,
                "mean_abs averages absolute frame energies");
  require_close(at(result.residues, 0, 0).signed_mean,
                at(result.atoms, 0, 0).signed_mean +
                    at(result.atoms, 0, 1).signed_mean,
                "residue signed mean closes over member atom edges");
  require(result.atom_edges.size() == 2, "atom top N edge selection");
  require(result.atom_edges[0].row == 1 && result.atom_edges[0].col == 2,
          "top atom edge sorted by mean_abs");
  require_close(result.atom_edges[0].stats.mean_abs, 5.0,
                "top atom edge mean_abs");
  require(!result.mol1_atom_rankings.empty(), "mol1 atom rankings emitted");
  require(result.mol1_atom_rankings[0].entity_index == 2,
          "rankings sorted by configured signed mean");
  require_close(result.mol1_atom_rankings[0].stats.signed_mean, 6.0,
                "ranking signed mean");
}

void run_cache_vs_fallback_contract() {
  Multi_MOL mol1 =
      make_molecule({make_atom(1, "ALA", "CA", "A"),
                     make_atom(2, "GLY", "CA", "A")},
                    {0.0, 3.0});
  Multi_MOL mol2 =
      make_molecule({make_atom(5, "LEU", "CA", "B")}, {12.0});
  arguments args;
  args.epsilon_in = 1.0;
  cache_sigma_epsilon provider_cache;
  provider_cache.force_field_name = "AMBER99SB";
  EnergyDecompositionProvider provider;
  provider.configure(mol1, mol2, provider_cache, args, -1.0, false);

  std::vector<double> totals;
  for (std::size_t mol1_atom = 0; mol1_atom < mol1.info.size(); ++mol1_atom) {
    const EnergyPairDecomposition pair =
        provider.pair_energy(mol1.m_mol[0], mol2.m_mol[0], mol1_atom, 0);
    totals.push_back(pair.total());
  }

  TitanExecutor executor(2);
  edecomp::EnergyDecompositionEngine from_cache(
      edecomp::make_axes(mol1, "mol1"), edecomp::make_axes(mol2, "mol2"),
      config());
  from_cache.begin_trajectory("cache");
  cache_sigma_epsilon cache = make_cache(totals);
  from_cache.process_frame({&cache, nullptr, nullptr, nullptr}, executor);
  const auto cache_result = from_cache.finalize_trajectory();

  edecomp::EnergyDecompositionEngine from_fallback(
      edecomp::make_axes(mol1, "mol1"), edecomp::make_axes(mol2, "mol2"),
      config());
  from_fallback.begin_trajectory("fallback");
  from_fallback.process_frame(
      {nullptr, &provider, &mol1.m_mol[0], &mol2.m_mol[0]}, executor);
  const auto fallback_result = from_fallback.finalize_trajectory();

  for (std::size_t i = 0; i < totals.size(); ++i) {
    require_close(at(cache_result.atoms, 0, i).signed_mean,
                  at(fallback_result.atoms, 0, i).signed_mean,
                  "cache and fallback total energy parity");
  }
}

void run_trajectory_aggregate_contract() {
  const Multi_MOL mol1 = make_molecule({make_atom(1, "ALA", "CA", "A")},
                                       {0.0});
  const Multi_MOL mol2 = make_molecule({make_atom(5, "LEU", "CA", "B")},
                                       {10.0});
  TitanExecutor executor(2);
  edecomp::EnergyDecompositionEngine engine(
      edecomp::make_axes(mol1, "mol1"), edecomp::make_axes(mol2, "mol2"),
      config());

  engine.begin_trajectory("short");
  cache_sigma_epsilon short_cache = make_cache({10.0});
  engine.process_frame({&short_cache, nullptr, nullptr, nullptr}, executor);
  const auto short_result = engine.finalize_trajectory();
  require(short_result.frames == 1, "short trajectory finalization");

  engine.begin_trajectory("long");
  for (int frame = 0; frame < 3; ++frame) {
    cache_sigma_epsilon long_cache = make_cache({2.0});
    engine.process_frame({&long_cache, nullptr, nullptr, nullptr}, executor);
  }
  const auto long_result = engine.finalize_trajectory();
  require(long_result.frames == 3, "long trajectory finalization");

  const edecomp::AllTrajectoryResult all = engine.finalize_all();
  require(all.trajectories == 2 && all.pooled_frames.frames == 4,
          "all-trajectory counters");
  require_close(at(all.trajectory_equal.atoms, 0, 0).signed_mean, 6.0,
                "trajectory equal aggregation");
  require_close(at(all.pooled_frames.atoms, 0, 0).signed_mean, 4.0,
                "pooled frame aggregation");
  const fs::path out_dir =
      fs::temp_directory_path() / "titan_energy_decomp_engine_all_smoke";
  fs::remove_all(out_dir);
  fs::create_directories(out_dir);
  edecomp::export_all_results(all, out_dir, "_average");
  require(fs::exists(out_dir / "edecomp_summary_average.dat") &&
              fs::exists(out_dir / "edecomp_summary_pooled_frames.dat") &&
              fs::exists(out_dir /
                         "edecomp_mm_residue_matrix_total_average.csv"),
          "aggregate export filenames follow edecomp output contract");
}

void run_nt_determinism_and_export_contract() {
  const Multi_MOL mol1 =
      make_molecule({make_atom(1, "ALA", "CA", "A"),
                     make_atom(2, "GLY", "CA", "A"),
                     make_atom(3, "SER", "CA", "A"),
                     make_atom(4, "TYR", "CA", "A")},
                    {0.0, 1.0, 2.0, 3.0});
  const Multi_MOL mol2 =
      make_molecule({make_atom(5, "LEU", "CA", "B"),
                     make_atom(6, "SER", "CA", "B"),
                     make_atom(7, "VAL", "CA", "B"),
                     make_atom(8, "ASN", "CA", "B")},
                    {10.0, 20.0, 30.0, 40.0});
  const std::vector<std::vector<double>> frames{
      {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16},
      {-1, 2, -3, 4, -5, 6, -7, 8, -9, 10, -11, 12, -13, 14, -15, 16}};
  TitanExecutor one(1);
  TitanExecutor four(4);
  const auto serial = process_one(mol1, mol2, frames, one);
  const auto parallel = process_one(mol1, mol2, frames, four);
  require(serial.atoms.values.size() == parallel.atoms.values.size(),
          "deterministic atom matrix size");
  for (std::size_t i = 0; i < serial.atoms.values.size(); ++i) {
    require(serial.atoms.values[i].signed_mean ==
                parallel.atoms.values[i].signed_mean &&
            serial.atoms.values[i].mean_abs ==
                parallel.atoms.values[i].mean_abs,
            "nt=1 and nt=4 deterministic reductions");
  }

  const fs::path out_dir =
      fs::temp_directory_path() / "titan_energy_decomp_engine_smoke";
  fs::remove_all(out_dir);
  fs::create_directories(out_dir);
  edecomp::export_result(serial, out_dir);
  const fs::path summary = out_dir / "edecomp_summary.dat";
  const fs::path residue = out_dir / "edecomp_mm_residue_matrix_total.csv";
  const fs::path residue_elec =
      out_dir / "edecomp_mm_residue_matrix_elec.csv";
  const fs::path residue_vdw = out_dir / "edecomp_mm_residue_matrix_vdw.csv";
  const fs::path edges = out_dir / "edecomp_mm_atom_edges.csv";
  require(fs::exists(summary) && fs::exists(residue) &&
              fs::exists(residue_elec) && fs::exists(residue_vdw) &&
              fs::exists(edges),
          "export filenames follow edecomp output contract");
  std::ifstream input(summary);
  std::ostringstream text;
  text << input.rdbuf();
  require(text.str().find("abs_mean = abs(signed_mean)") != std::string::npos &&
              text.str().find("mean_abs = mean(abs(frame_energy))") !=
                  std::string::npos,
          "summary documents signed aggregation semantics");
}

void run_endpoint_entity_contract() {
  const Multi_MOL mol1 =
      make_molecule({make_atom(1, "ALA", "CA", "A"),
                     make_atom(2, "GLY", "CA", "A")},
                    {0.0, 1.0});
  const Multi_MOL mol2 =
      make_molecule({make_atom(5, "LEU", "CA", "B")}, {10.0});
  edecomp::EndpointEntityDecompositionEngine engine(
      edecomp::make_axes(mol1, "mol1"), edecomp::make_axes(mol2, "mol2"),
      config(16));
  engine.begin_trajectory("endpoint");
  engine.process_frame({{edecomp::EntityKind::Residue, "mol1", 0,
                         -2.0, 1.0, -0.5, 0.25},
                        {edecomp::EntityKind::Residue, "mol2", 0,
                         1.0, 2.0, -1.0, 0.5},
                        {edecomp::EntityKind::Atom, "mol1", 0,
                         -2.0, 1.0, -0.5, 0.25}});
  engine.process_frame({{edecomp::EntityKind::Residue, "mol1", 0,
                         4.0, -1.0, -0.5, 0.25},
                        {edecomp::EntityKind::Residue, "mol2", 0,
                         1.0, 2.0, -1.0, 0.5},
                        {edecomp::EntityKind::Atom, "mol1", 0,
                         4.0, -1.0, -0.5, 0.25}});
  const auto result = engine.finalize_trajectory();
  require(result.frames == 2, "endpoint frame count");
  require(result.residue_map.size() == 3, "endpoint residue map includes all entities");
  const auto &mol1_res0 = result.residue_map[0];
  require_close(mol1_res0.stats.e_elec_mean, 1.0,
                "endpoint elec signed mean");
  require_close(mol1_res0.stats.e_mm_mean, 1.0,
                "endpoint MM signed mean");
  require_close(mol1_res0.stats.delta_g_no_entropy_mean, 0.75,
                "endpoint deltaG signed mean");
  require_close(mol1_res0.stats.abs_mean, 0.75,
                "endpoint abs_mean is abs of signed deltaG");
  require_close(mol1_res0.stats.mean_abs, 2.0,
                "endpoint mean_abs averages abs frame deltaG");
  require(!result.residue_ranking.empty() &&
              result.residue_ranking.front().molecule_role == "mol2",
          "endpoint ranking sorts by abs_mean");

  const fs::path out_dir =
      fs::temp_directory_path() / "titan_energy_decomp_endpoint_smoke";
  fs::remove_all(out_dir);
  fs::create_directories(out_dir);
  edecomp::export_endpoint_result(result, out_dir, "mmgbsa");
  require(fs::exists(out_dir / "edecomp_mmgbsa_residue_map.csv") &&
              fs::exists(out_dir / "edecomp_mmgbsa_residue_ranking.csv") &&
              fs::exists(out_dir / "edecomp_mmgbsa_atom_map.csv") &&
              fs::exists(out_dir / "edecomp_mmgbsa_atom_ranking.csv") &&
              !fs::exists(out_dir / "edecomp_mmgbsa_residue_matrix.csv"),
          "endpoint export uses map/ranking files only");
  {
    std::ifstream input(out_dir / "edecomp_mmgbsa_residue_map.csv");
    std::string header;
    std::getline(input, header);
    require(header.find("G_polar_GB") != std::string::npos,
            "MM/GBSA endpoint map uses GBSA component labels");
  }

  fs::remove_all(out_dir);
  fs::create_directories(out_dir);
  edecomp::export_endpoint_result(result, out_dir, "mmpbsa");
  require(fs::exists(out_dir / "edecomp_mmpbsa_residue_map.csv") &&
              fs::exists(out_dir / "edecomp_mmpbsa_residue_ranking.csv") &&
              fs::exists(out_dir / "edecomp_mmpbsa_atom_map.csv") &&
              fs::exists(out_dir / "edecomp_mmpbsa_atom_ranking.csv") &&
              !fs::exists(out_dir / "edecomp_mmpbsa_residue_matrix.csv"),
          "MM/PBSA endpoint export uses map/ranking files only");
  {
    std::ifstream input(out_dir / "edecomp_mmpbsa_residue_map.csv");
    std::string header;
    std::getline(input, header);
    require(header.find("DeltaDeltaG_solv_PB") != std::string::npos &&
                header.find("DeltaG_nonpolar_SASA") != std::string::npos,
            "MM/PBSA endpoint map uses PBSA component labels");
  }
  {
    std::ifstream input(out_dir / "edecomp_summary.dat");
    std::ostringstream text;
    text << input.rdbuf();
    require(text.str().find(
                "MM/PBSA maps are endpoint entity decompositions") !=
                std::string::npos,
            "MM/PBSA summary documents endpoint semantics");
  }
}

} // namespace

int main() {
  run_cache_residue_atom_ranking_contract();
  run_cache_vs_fallback_contract();
  run_trajectory_aggregate_contract();
  run_nt_determinism_and_export_contract();
  run_endpoint_entity_contract();
  std::cout << "PASS titan_energy_decomp_engine_smoke\n";
  return 0;
}
