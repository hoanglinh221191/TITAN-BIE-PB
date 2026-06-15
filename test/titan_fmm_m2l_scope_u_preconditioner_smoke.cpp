#include "tabi/tabi_fmm_m2l_preconditioner.h"
#include "tabi/tabi_fmm_m2l_prototype.h"
#include "tabi/tabi_types.h"
#include "titan_executor.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using titan::tabi::fmm_m2l::FmmM2lAnalyticContext;
using titan::tabi::fmm_m2l::KernelKind;
using titan::tabi::fmm_m2l::PrototypeConfig;
using titan::tabi::fmm_m2l::PrototypePanel;
using TitanPBE::Tabi::Panel;

void fail(const std::string &message) {
  std::cerr << "FAIL titan_fmm_m2l_scope_u_preconditioner_smoke: "
            << message << '\n';
  std::exit(1);
}

void expect(bool condition, const std::string &message) {
  if (!condition)
    fail(message);
}

std::vector<Panel> make_panels(std::size_t count) {
  std::vector<Panel> panels(count);
  for (std::size_t i = 0; i < count; ++i) {
    const double t = static_cast<double>(i);
    panels[i].x = 0.35 * std::cos(0.41 * t) + 0.015 * t;
    panels[i].y = 0.35 * std::sin(0.37 * t);
    panels[i].z = 0.06 * static_cast<double>(i % 7);
    panels[i].nx = 0.0;
    panels[i].ny = 0.0;
    panels[i].nz = 1.0;
    panels[i].area = 0.45 + 0.01 * static_cast<double>(i % 5);
    panels[i].diagonal = 1.0;
  }
  return panels;
}

std::vector<PrototypePanel> to_prototype(const std::vector<Panel> &panels) {
  std::vector<PrototypePanel> converted;
  converted.reserve(panels.size());
  for (const Panel &panel : panels) {
    PrototypePanel out;
    out.x = panel.x;
    out.y = panel.y;
    out.z = panel.z;
    out.nx = panel.nx;
    out.ny = panel.ny;
    out.nz = panel.nz;
    out.area = panel.area;
    converted.push_back(out);
  }
  return converted;
}

void test_fmm_block_plan_covers_panels_once() {
  PrototypeConfig config;
  config.kernel = KernelKind::Yukawa;
  config.kappa = 0.12;
  config.max_leaf_size = 8;
  config.expansion_order = 4;
  config.separation_theta = 0.6;

  const std::vector<Panel> panels = make_panels(37);
  FmmM2lAnalyticContext context(config, to_prototype(panels), 1.0 / 80.0,
                                config.kappa);

  const auto blocks = context.block_panel_ids(6);
  expect(!blocks.empty(), "FMM block plan should contain at least one block");

  std::vector<int> seen(panels.size(), 0);
  for (const auto &block : blocks) {
    expect(!block.empty(), "FMM block plan should not emit empty blocks");
    expect(block.size() <= 6, "FMM block plan should respect requested cap");
    for (const std::size_t panel : block) {
      expect(panel < panels.size(), "FMM block plan panel id out of range");
      ++seen[panel];
    }
  }
  for (std::size_t i = 0; i < seen.size(); ++i) {
    expect(seen[i] == 1, "FMM block plan should cover each panel once");
  }
}

void test_fmm_block_preconditioner_builds_and_applies() {
  PrototypeConfig config;
  config.kernel = KernelKind::Yukawa;
  config.kappa = 0.08;
  config.max_leaf_size = 8;
  config.expansion_order = 4;
  config.separation_theta = 0.6;

  const std::vector<Panel> panels = make_panels(32);
  FmmM2lAnalyticContext context(config, to_prototype(panels), 1.0 / 80.0,
                                config.kappa);
  TitanExecutor executor(4);

  std::vector<double> diagonal(2 * panels.size(), 2.0);
  const auto preconditioner =
      TitanPBE::Tabi::build_fmm_m2l_block_jacobi_preconditioner(
          panels, context, true, 1.0 / 80.0, config.kappa, 0.5, 0.5,
          diagonal, 8, executor);

  expect(preconditioner != nullptr, "FMM block preconditioner should build");
  expect(preconditioner->panel_count == panels.size(),
         "FMM block preconditioner panel count mismatch");
  expect(preconditioner->coupled,
         "FMM block preconditioner should preserve coupled mode");
  expect(preconditioner->valid_blocks > 0,
         "FMM block preconditioner should have valid blocks");
  expect(preconditioner->valid_panel_count >=
             static_cast<std::size_t>(0.70 * static_cast<double>(panels.size())),
         "FMM block preconditioner should cover most panels");

  std::vector<double> rhs(2 * panels.size(), 0.0);
  for (std::size_t i = 0; i < rhs.size(); ++i)
    rhs[i] = 0.25 + 0.01 * static_cast<double>(i % 11);
  std::vector<double> out;
  preconditioner->apply_parallel(rhs, out, executor);
  expect(out.size() == rhs.size(), "FMM block preconditioner output size");
  for (const double value : out)
    expect(std::isfinite(value), "FMM block preconditioner non-finite output");
}

void test_fmm_block_preconditioner_accepts_large_caps() {
  PrototypeConfig config;
  config.kernel = KernelKind::Yukawa;
  config.kappa = 0.05;
  config.max_leaf_size = 256;
  config.expansion_order = 4;
  config.separation_theta = 0.6;

  const std::vector<Panel> panels = make_panels(150);
  FmmM2lAnalyticContext context(config, to_prototype(panels), 1.0 / 80.0,
                                config.kappa);
  TitanExecutor executor(4);

  const std::vector<double> diagonal(2 * panels.size(), 2.0);
  const std::size_t bytes32 =
      TitanPBE::Tabi::estimate_fmm_m2l_block_jacobi_bytes(
          context, diagonal.size(), true, 32);
  const std::size_t bytes64 =
      TitanPBE::Tabi::estimate_fmm_m2l_block_jacobi_bytes(
          context, diagonal.size(), true, 64);
  const std::size_t bytes128 =
      TitanPBE::Tabi::estimate_fmm_m2l_block_jacobi_bytes(
          context, diagonal.size(), true, 128);
  expect(bytes64 > bytes32,
         "FMM block preconditioner cap 64 should estimate more bytes than 32");
  expect(bytes128 > bytes64,
         "FMM block preconditioner cap 128 should estimate more bytes than 64");

  const auto preconditioner =
      TitanPBE::Tabi::build_fmm_m2l_block_jacobi_preconditioner(
          panels, context, true, 1.0 / 80.0, config.kappa, 0.5, 0.5,
          diagonal, 128, executor);
  expect(preconditioner != nullptr, "FMM cap 128 preconditioner should build");
  std::size_t max_block_panels = 0;
  for (const auto &block : preconditioner->blocks) {
    max_block_panels = std::max(max_block_panels, block.panels.size());
    expect(block.panels.size() <= 128,
           "FMM block preconditioner should respect cap 128");
  }
  expect(max_block_panels > 32,
         "FMM block preconditioner cap 128 should create blocks larger than 32");
}

void test_fmm_overlap_groups_are_restricted_and_deterministic() {
  PrototypeConfig config;
  config.kernel = KernelKind::Yukawa;
  config.kappa = 0.07;
  config.max_leaf_size = 8;
  config.expansion_order = 4;
  config.separation_theta = 0.6;

  const std::vector<Panel> panels = make_panels(48);
  FmmM2lAnalyticContext context(config, to_prototype(panels), 1.0 / 80.0,
                                config.kappa);

  const auto block_groups = context.block_panel_ids(6);
  const auto no_overlap_groups = context.overlap_panel_groups(6, 0);
  expect(no_overlap_groups.size() == block_groups.size(),
         "overlap cap 0 should preserve FMM block group count");
  for (std::size_t i = 0; i < block_groups.size(); ++i) {
    expect(no_overlap_groups[i].owned_panels == block_groups[i],
           "overlap cap 0 should preserve owned panels");
    expect(no_overlap_groups[i].factor_panels == block_groups[i],
           "overlap cap 0 should factor only owned panels");
    expect(no_overlap_groups[i].halo_count == 0,
           "overlap cap 0 should not emit halo panels");
  }

  const auto overlap_groups = context.overlap_panel_groups(6, 12);
  expect(overlap_groups.size() == block_groups.size(),
         "overlap groups should keep the owned block partition");

  std::vector<int> owned_seen(panels.size(), 0);
  bool saw_halo = false;
  for (const auto &group : overlap_groups) {
    expect(!group.owned_panels.empty(),
           "overlap group should not have empty owned panels");
    expect(group.owned_panels.size() <= 6,
           "overlap group should respect owned block cap");
    expect(group.factor_panels.size() >= group.owned_panels.size(),
           "factor panels should include owned panels");
    expect(group.halo_count + group.owned_panels.size() ==
               group.factor_panels.size(),
           "halo count should match factor minus owned panels");
    saw_halo = saw_halo || group.halo_count > 0;

    std::vector<int> local_seen(panels.size(), 0);
    for (std::size_t i = 0; i < group.owned_panels.size(); ++i) {
      const std::size_t panel = group.owned_panels[i];
      expect(panel < panels.size(), "owned overlap panel out of range");
      expect(group.factor_panels[i] == panel,
             "factor panels should store owned panels first");
      ++owned_seen[panel];
    }
    for (const std::size_t panel : group.factor_panels) {
      expect(panel < panels.size(), "factor overlap panel out of range");
      ++local_seen[panel];
      expect(local_seen[panel] <= 1,
             "factor panels should not duplicate panels within a group");
    }
  }
  for (std::size_t i = 0; i < owned_seen.size(); ++i)
    expect(owned_seen[i] == 1,
           "overlap owned panels should cover each panel exactly once");
  expect(saw_halo, "positive overlap cap should add at least one halo panel");
}

void test_fmm_overlap_schwarz_preconditioner_builds_and_applies() {
  PrototypeConfig config;
  config.kernel = KernelKind::Yukawa;
  config.kappa = 0.08;
  config.max_leaf_size = 8;
  config.expansion_order = 4;
  config.separation_theta = 0.6;

  const std::vector<Panel> panels = make_panels(40);
  FmmM2lAnalyticContext context(config, to_prototype(panels), 1.0 / 80.0,
                                config.kappa);
  TitanExecutor executor(4);

  std::vector<double> diagonal(2 * panels.size(), 2.0);
  const std::size_t bytes0 =
      TitanPBE::Tabi::estimate_fmm_m2l_overlap_schwarz_bytes(
          context, diagonal.size(), true, 8, 0);
  const std::size_t bytes16 =
      TitanPBE::Tabi::estimate_fmm_m2l_overlap_schwarz_bytes(
          context, diagonal.size(), true, 8, 16);
  expect(bytes16 > bytes0,
         "overlap Schwarz estimated bytes should grow with halo cap");

  const auto preconditioner =
      TitanPBE::Tabi::build_fmm_m2l_overlap_schwarz_preconditioner(
          panels, context, true, 1.0 / 80.0, config.kappa, 0.5, 0.5,
          diagonal, 8, 16, executor);

  expect(preconditioner != nullptr,
         "FMM overlap-Schwarz preconditioner should build");
  expect(preconditioner->valid_blocks > 0,
         "FMM overlap-Schwarz should have valid blocks");
  expect(preconditioner->valid_panel_count >=
             static_cast<std::size_t>(0.70 * static_cast<double>(panels.size())),
         "FMM overlap-Schwarz should cover most owned panels");

  bool saw_restricted_block = false;
  for (const auto &block : preconditioner->blocks) {
    expect(block.owned_count > 0, "overlap block should record owned count");
    expect(block.owned_count <= block.panels.size(),
           "owned count should not exceed factor panel count");
    saw_restricted_block = saw_restricted_block ||
                           block.owned_count < block.panels.size();
  }
  expect(saw_restricted_block,
         "positive overlap should create at least one restricted block");

  std::vector<double> rhs(2 * panels.size(), 0.0);
  for (std::size_t i = 0; i < rhs.size(); ++i)
    rhs[i] = 0.15 + 0.02 * static_cast<double>(i % 13);
  std::vector<double> out;
  preconditioner->apply_parallel(rhs, out, executor);
  expect(out.size() == rhs.size(), "FMM overlap-Schwarz output size");
  for (const double value : out)
    expect(std::isfinite(value), "FMM overlap-Schwarz non-finite output");
}

void test_fmm_nearfield_graph_is_sorted_and_deterministic() {
  PrototypeConfig config;
  config.kernel = KernelKind::Yukawa;
  config.kappa = 0.06;
  config.max_leaf_size = 8;
  config.expansion_order = 4;
  config.separation_theta = 0.6;

  const std::vector<Panel> panels = make_panels(44);
  FmmM2lAnalyticContext context(config, to_prototype(panels), 1.0 / 80.0,
                                config.kappa);

  const auto graph64 = context.near_panel_graph(64);
  const auto graph128 = context.near_panel_graph(128);
  expect(graph64.size() == panels.size(),
         "nearfield graph should have one row per panel");
  expect(graph128.size() == panels.size(),
         "nearfield graph cap 128 should have one row per panel");
  std::size_t edges64 = 0;
  std::size_t edges128 = 0;
  for (std::size_t row = 0; row < panels.size(); ++row) {
    const auto &neighbors64 = graph64[row];
    const auto &neighbors128 = graph128[row];
    expect(!neighbors64.empty(), "nearfield graph row should include self");
    expect(std::find(neighbors64.begin(), neighbors64.end(), row) !=
               neighbors64.end(),
           "nearfield graph row should include self panel");
    expect(std::is_sorted(neighbors64.begin(), neighbors64.end()),
           "nearfield graph neighbors should be sorted");
    expect(std::adjacent_find(neighbors64.begin(), neighbors64.end()) ==
               neighbors64.end(),
           "nearfield graph neighbors should be unique");
    expect(neighbors64.size() <= 64,
           "nearfield graph should respect neighbor cap");
    expect(neighbors128.size() >= neighbors64.size(),
           "nearfield graph should grow monotonically with cap");
    edges64 += neighbors64.size();
    edges128 += neighbors128.size();
  }
  expect(edges128 >= edges64,
         "nearfield graph total edges should grow monotonically with cap");
}

void test_fmm_nearfield_ilu0_builds_and_applies() {
  PrototypeConfig config;
  config.kernel = KernelKind::Yukawa;
  config.kappa = 0.08;
  config.max_leaf_size = 8;
  config.expansion_order = 4;
  config.separation_theta = 0.6;

  const std::vector<Panel> panels = make_panels(36);
  FmmM2lAnalyticContext context(config, to_prototype(panels), 1.0 / 80.0,
                                config.kappa);
  TitanExecutor executor(4);
  std::vector<double> diagonal(2 * panels.size(), 2.0);

  const std::size_t bytes64 =
      TitanPBE::Tabi::estimate_fmm_m2l_nearfield_ilu0_bytes(
          context, diagonal.size(), true, 64);
  const std::size_t bytes128 =
      TitanPBE::Tabi::estimate_fmm_m2l_nearfield_ilu0_bytes(
          context, diagonal.size(), true, 128);
  expect(bytes64 > 0, "nearfield ILU estimated bytes should be positive");
  expect(bytes128 >= bytes64,
         "nearfield ILU estimated bytes should grow with neighbor cap");

  const auto ilu =
      TitanPBE::Tabi::build_fmm_m2l_nearfield_ilu0_preconditioner(
          panels, context, true, 1.0 / 80.0, config.kappa, 0.5, 0.5,
          diagonal, 96, 0.0, executor);
  expect(ilu != nullptr, "nearfield ILU should build");
  expect(ilu->panel_count == panels.size(),
         "nearfield ILU panel count mismatch");
  expect(ilu->valid_panel_count == panels.size(),
         "nearfield ILU should cover every panel");
  expect(ilu->nnz_blocks >= panels.size(),
         "nearfield ILU should have at least diagonal blocks");
  expect(ilu->max_degree >= ilu->avg_degree,
         "nearfield ILU max degree should dominate average degree");

  std::vector<double> rhs(2 * panels.size(), 0.0);
  for (std::size_t i = 0; i < rhs.size(); ++i)
    rhs[i] = 0.2 + 0.015 * static_cast<double>(i % 17);
  std::vector<double> out;
  ilu->apply(rhs, out);
  expect(out.size() == rhs.size(), "nearfield ILU output size");
  for (const double value : out)
    expect(std::isfinite(value), "nearfield ILU non-finite output");
}

} // namespace

int main() {
  test_fmm_block_plan_covers_panels_once();
  test_fmm_block_preconditioner_builds_and_applies();
  test_fmm_block_preconditioner_accepts_large_caps();
  test_fmm_overlap_groups_are_restricted_and_deterministic();
  test_fmm_overlap_schwarz_preconditioner_builds_and_applies();
  test_fmm_nearfield_graph_is_sorted_and_deterministic();
  test_fmm_nearfield_ilu0_builds_and_applies();
  std::cout << "titan_fmm_m2l_scope_u_preconditioner_smoke ok\n";
  return 0;
}
