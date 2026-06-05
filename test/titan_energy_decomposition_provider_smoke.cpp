#include "titan_energy_decomposition.h"

#include "arguments.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

std::string arguments::set_Buildnumber() {
  return "titan_energy_decomposition_provider_smoke";
}

namespace {

void require(bool condition, const std::string &message) {
  if (!condition) {
    std::cerr << "FAIL titan_energy_decomposition_provider_smoke: " << message
              << '\n';
    std::exit(1);
  }
}

void require_close(double actual, double expected, const std::string &message,
                   double tolerance = 1.0e-12) {
  if (!std::isfinite(actual) || std::fabs(actual - expected) > tolerance) {
    std::cerr << "FAIL titan_energy_decomposition_provider_smoke: " << message
              << " actual=" << actual << " expected=" << expected << '\n';
    std::exit(1);
  }
}

box make_box(double length_a) {
  box out;
  out.set_a(length_a);
  out.set_b(length_a);
  out.set_c(length_a);
  out.set_alpha(90.0);
  out.set_beta(90.0);
  out.set_gamma(90.0);
  return out;
}

Multi_MOL make_molecule(const std::vector<double> &xs, double charge_scale,
                        double sigma_scale, double epsilon_scale,
                        double box_length_a = 80.0) {
  Multi_MOL out;
  out.info.resize(xs.size());
  out.parameters.resize(xs.size());
  for (std::size_t i = 0; i < xs.size(); ++i) {
    out.parameters[i].set_data(charge_scale * static_cast<double>(i + 1), 0.1,
                               sigma_scale * static_cast<double>(i + 1),
                               epsilon_scale * static_cast<double>(i + 1),
                               1.0);
  }
  MOL frame(make_box(box_length_a), 0.0, "AMBER99SB");
  frame.pos.x.assign(xs.begin(), xs.end());
  frame.pos.y.assign(xs.size(), 0.0);
  frame.pos.z.assign(xs.size(), 0.0);
  out.m_mol.push_back(std::move(frame));
  return out;
}

EnergyDecompositionProvider make_provider(const Multi_MOL &rec,
                                          const Multi_MOL &lig,
                                          cache_sigma_epsilon &cache,
                                          arguments &args, double cutoff_nm,
                                          bool cutoff_active) {
  cache.force_field_name = "AMBER99SB";
  args.force_f = "AMBER99SB";
  args.epsilon_in = 1.0;
  EnergyDecompositionProvider provider;
  provider.configure(rec, lig, cache, args, cutoff_nm, cutoff_active);
  return provider;
}

void require_batch_matches_scalar(const EnergyDecompositionProvider &provider,
                                  const MOL &rec_frame, const MOL &lig_frame,
                                  std::size_t rec_begin, std::size_t rec_end,
                                  std::size_t lig_begin, std::size_t lig_end,
                                  const std::string &label) {
  const std::size_t rec_count = rec_end - rec_begin;
  const std::size_t lig_count = lig_end - lig_begin;
  std::vector<double> elec(rec_count * lig_count, -1.0);
  std::vector<double> vdw(rec_count * lig_count, -1.0);
  require(provider.inter_pair_energies_batch(rec_frame, lig_frame, rec_begin,
                                             rec_end, lig_begin, lig_end, elec,
                                             vdw),
          label + " batch accepted");

  for (std::size_t lig_pos = lig_begin; lig_pos < lig_end; ++lig_pos) {
    for (std::size_t rec_pos = rec_begin; rec_pos < rec_end; ++rec_pos) {
      const auto scalar =
          provider.pair_energy(rec_frame, lig_frame, rec_pos, lig_pos);
      const std::size_t row = lig_pos - lig_begin;
      const std::size_t col = rec_pos - rec_begin;
      const std::size_t index =
          titan_energy_cache_index(row, rec_count, col);
      require_close(elec[index], scalar.elec, label + " elec parity");
      require_close(vdw[index], scalar.vdw, label + " vdw parity");
    }
  }
}

void run_cache_index_and_minimal_amber_contract() {
  require(titan_energy_cache_index(2, 5, 3) == 13,
          "cache index is mol2 row-major over mol1 columns");
  Multi_MOL rec =
      make_molecule({0.0, 3.0, 7.0}, 0.1, 0.12, 0.20);
  Multi_MOL lig = make_molecule({11.0, 17.0}, -0.2, 0.15, 0.25);
  arguments args;
  cache_sigma_epsilon cache;
  const auto provider = make_provider(rec, lig, cache, args, -1.0, false);
  require_batch_matches_scalar(provider, rec.m_mol[0], lig.m_mol[0], 1, 3, 0,
                               2, "minimal AMBER subset");
}

void run_cutoff_and_mic_contract() {
  Multi_MOL rec = make_molecule({1.0}, 0.1, 0.12, 0.20, 20.0);
  Multi_MOL lig = make_molecule({19.0, 9.0}, -0.2, 0.15, 0.25, 20.0);
  arguments args;
  cache_sigma_epsilon cache;
  const auto provider = make_provider(rec, lig, cache, args, 0.3, true);
  std::vector<double> elec(2, -1.0);
  std::vector<double> vdw(2, -1.0);
  require(provider.inter_pair_energies_batch(rec.m_mol[0], lig.m_mol[0], 0, 1,
                                             0, 2, elec, vdw),
          "MIC/cutoff batch accepted");
  require_close(elec[0], provider.pair_energy(rec.m_mol[0], lig.m_mol[0], 0, 0)
                             .elec,
                "MIC wrapped pair elec matches scalar");
  require(std::fabs(elec[0]) > 0.0, "MIC wrapped pair remains inside cutoff");
  require_close(elec[1], 0.0, "pair beyond cutoff has zero elec");
  require_close(vdw[1], 0.0, "pair beyond cutoff has zero vdw");
}

void run_tail_and_determinism_contract() {
  Multi_MOL rec =
      make_molecule({0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0}, 0.1, 0.12, 0.20);
  Multi_MOL lig = make_molecule({10.0, 14.0}, -0.2, 0.15, 0.25);
  arguments args;
  cache_sigma_epsilon cache;
  const auto provider = make_provider(rec, lig, cache, args, -1.0, false);
  for (std::size_t rec_count : {5u, 6u, 7u}) {
    require_batch_matches_scalar(provider, rec.m_mol[0], lig.m_mol[0], 0,
                                 rec_count, 0, 2,
                                 "AMBER vector block plus tail");
  }

  std::vector<double> first_elec(14, -1.0);
  std::vector<double> first_vdw(14, -1.0);
  std::vector<double> second_elec(14, -1.0);
  std::vector<double> second_vdw(14, -1.0);
  require(provider.inter_pair_energies_batch(rec.m_mol[0], lig.m_mol[0], 0, 7,
                                             0, 2, first_elec, first_vdw),
          "first deterministic batch accepted");
  require(provider.inter_pair_energies_batch(rec.m_mol[0], lig.m_mol[0], 0, 7,
                                             0, 2, second_elec, second_vdw),
          "second deterministic batch accepted");
  require(first_elec == second_elec && first_vdw == second_vdw,
          "repeated batch output is deterministic");
}

} // namespace

int main() {
  run_cache_index_and_minimal_amber_contract();
  run_cutoff_and_mic_contract();
  run_tail_and_determinism_contract();
  std::cout << "PASS titan_energy_decomposition_provider_smoke\n";
  return 0;
}
