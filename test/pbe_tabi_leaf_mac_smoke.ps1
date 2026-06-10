param(
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_leaf_mac_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot 'pbe_tabi_treecode_sources.ps1')
$Compiler = 'C:\msys64\ucrt64\bin\clang++.exe'
$SourcePath = Join-Path $WorkRoot 'pbe_tabi_leaf_mac_smoke.cpp'
$ExePath = Join-Path $WorkRoot 'pbe_tabi_leaf_mac_smoke.exe'

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside leaf-MAC smoke root: $full"
    }
}

function New-CleanDir([string]$Path) {
    Assert-UnderPath -Path $Path -Root $WorkRoot
    if (Test-Path -LiteralPath $Path) {
        Remove-Item -LiteralPath $Path -Recurse -Force
    }
    New-Item -ItemType Directory -Path $Path -Force | Out-Null
}

if (-not (Test-Path -LiteralPath $Compiler)) {
    throw "Missing clang++ compiler: $Compiler"
}

New-CleanDir $WorkRoot

$source = @'
#include "titan_executor.h"
#include "tabi/tabi_types.h"
#include "tabi/tabi_treecode_internal.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void require_true(const std::string &label, bool value) {
  if (!value)
    throw std::runtime_error(label);
}

double rms(const std::vector<double> &values) {
  double sum = 0.0;
  for (double value : values)
    sum += value * value;
  return values.empty() ? 0.0
                        : std::sqrt(sum / static_cast<double>(values.size()));
}

double relative_rms(const std::vector<double> &reference,
                    const std::vector<double> &candidate) {
  require_true("same vector size", reference.size() == candidate.size());
  std::vector<double> diff(reference.size(), 0.0);
  for (std::size_t i = 0; i < reference.size(); ++i)
    diff[i] = candidate[i] - reference[i];
  return rms(diff) / std::max(1.0e-12, rms(reference));
}

void require_close(const std::string &label, double a, double b,
                   double tolerance) {
  const double scale = std::max(1.0, std::max(std::abs(a), std::abs(b)));
  if (std::abs(a - b) > tolerance * scale)
    throw std::runtime_error(label);
}

void require_same_estimator(
    const TitanPBE::Tabi::detail::ConvergenceEstimatorDecision &cached,
    const TitanPBE::Tabi::detail::ConvergenceEstimatorDecision &fallback) {
  require_true("estimator active matches", cached.active == fallback.active);
  require_true("estimator unstable matches",
               cached.unstable == fallback.unstable);
  require_true("estimator split_reject matches",
               cached.split_reject == fallback.split_reject);
  require_close("estimator source_norm matches", cached.source_norm,
                fallback.source_norm, 1.0e-14);
  require_close("estimator tail_norm_p matches", cached.tail_norm_p,
                fallback.tail_norm_p, 1.0e-14);
  require_close("estimator tail_norm_pm1 matches", cached.tail_norm_pm1,
                fallback.tail_norm_pm1, 1.0e-14);
  require_close("estimator estimated_error matches",
                cached.estimated_error, fallback.estimated_error, 1.0e-14);
  for (std::size_t channel = 0;
       channel < TitanPBE::Tabi::kJcpExpansionChannelCount; ++channel) {
    require_close("estimator channel source_norm matches",
                  cached.channel_source_norm[channel],
                  fallback.channel_source_norm[channel], 1.0e-14);
    require_close("estimator channel tail_norm_p matches",
                  cached.channel_tail_norm_p[channel],
                  fallback.channel_tail_norm_p[channel], 1.0e-14);
    require_close("estimator channel tail_norm_pm1 matches",
                  cached.channel_tail_norm_pm1[channel],
                  fallback.channel_tail_norm_pm1[channel], 1.0e-14);
    require_close("estimator channel error matches",
                  cached.channel_estimated_error[channel],
                  fallback.channel_estimated_error[channel], 1.0e-14);
    require_true("estimator channel unstable matches",
                 cached.channel_unstable[channel] ==
                     fallback.channel_unstable[channel]);
  }
}

TitanPBE::Tabi::Panel make_panel(double x, double y, double z,
                                 std::size_t index) {
  TitanPBE::Tabi::Panel p;
  p.x = x;
  p.y = y;
  p.z = z;
  const double norm = std::sqrt(x * x + y * y + z * z);
  p.nx = norm > 1.0e-12 ? x / norm : 1.0;
  p.ny = norm > 1.0e-12 ? y / norm : 0.0;
  p.nz = norm > 1.0e-12 ? z / norm : 0.0;
  p.area = 0.2 + 0.005 * static_cast<double>(index % 11);
  p.atom_index = index;
  return p;
}

std::vector<TitanPBE::Tabi::Panel> build_two_cluster_fixture() {
  std::vector<TitanPBE::Tabi::Panel> panels;
  panels.reserve(96);
  std::size_t index = 0;
  for (int cluster = 0; cluster < 3; ++cluster) {
    const double cx = cluster == 0 ? -14.0 : (cluster == 1 ? 14.0 : 0.0);
    const double cy = cluster == 2 ? 18.0 : 0.0;
    for (int i = 0; i < 32; ++i) {
      const double a = 0.19634954084936207 * static_cast<double>(i);
      const double r = 0.7 + 0.02 * static_cast<double>(i % 5);
      panels.push_back(make_panel(cx + r * std::cos(a), cy + r * std::sin(a),
                                  0.1 * std::cos(2.0 * a), index++));
    }
  }
  return panels;
}

std::vector<double> make_density(std::size_t n) {
  std::vector<double> density(n, 0.0);
  for (std::size_t i = 0; i < n; ++i)
    density[i] = 0.5 + 0.2 * std::sin(0.13 * static_cast<double>(i));
  return density;
}

std::vector<double> make_unknowns(std::size_t n) {
  std::vector<double> unknowns(2 * n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    unknowns[i] = 0.7 + 0.20 * std::sin(0.19 * static_cast<double>(i));
    unknowns[n + i] = -1.1 + 0.15 * std::cos(0.27 * static_cast<double>(i));
  }
  return unknowns;
}

TitanPBE::Tabi::TreeWorkspace make_tree(
    const std::vector<TitanPBE::Tabi::Panel> &panels, TitanExecutor &executor,
    const std::string &leaf_mac) {
  TitanPBE::Tabi::TreeWorkspace tree =
      TitanPBE::Tabi::build_tree_workspace(
          panels, executor, 16, 0.55, "exact", "tile-avx2", "strict", "off",
          "onthefly", "tight-aabb", "topdown-aabb", "off", 1.5, "off",
          "area-centroid", "off", "sampled", 256, 2048, "translated",
          "naive", "off", "off", 0.0, 0.0, "off", 1.0e4, 1.0e-3, 0.25,
          "stage3j", "off", "off", 4, 16, 2);
  tree.taylor_order = 3;
  tree.derivative_mode = "analytic";
  tree.leaf_mac = leaf_mac;
  return tree;
}

std::vector<std::size_t> make_rank(const TitanPBE::Tabi::TreeWorkspace &tree,
                                   std::size_t n) {
  std::vector<std::size_t> rank(n, 0);
  for (std::size_t pos = 0; pos < tree.order.size(); ++pos)
    if (tree.order[pos] < rank.size())
      rank[tree.order[pos]] = pos;
  return rank;
}

std::size_t find_direct_candidate_leaf(
    const std::vector<TitanPBE::Tabi::Panel> &panels,
    const TitanPBE::Tabi::TreeWorkspace &tree, std::size_t target_index) {
  const std::vector<std::size_t> rank = make_rank(tree, panels.size());
  const std::size_t target_rank = rank[target_index];
  for (std::size_t node_index = 0; node_index < tree.nodes.size();
       ++node_index) {
    const TitanPBE::Tabi::TreeNode &node = tree.nodes[node_index];
    if (node.child_count != 0 || node.panel_count == 0)
      continue;
    const bool contains_target =
        target_rank >= node.first_panel &&
        target_rank < node.first_panel + node.panel_count;
    if (contains_target)
      continue;
    const double dx = panels[target_index].x - node.cx;
    const double dy = panels[target_index].y - node.cy;
    const double dz = panels[target_index].z - node.cz;
    const double distance =
        std::sqrt(std::max(1.0e-16, dx * dx + dy * dy + dz * dz));
    const TitanPBE::Tabi::detail::MacDecision mac =
        TitanPBE::Tabi::detail::evaluate_mac(panels[target_index], node, tree,
                                             distance);
    if (mac.accepted)
      return node_index;
  }
  throw std::runtime_error("no directness leaf-MAC candidate");
}

TitanPBE::Tabi::TreecodeStats apply_jcp(
    const std::vector<TitanPBE::Tabi::Panel> &panels,
    const std::vector<double> &unknowns, TitanExecutor &executor,
    const std::string &leaf_mac, std::vector<double> &out) {
  TitanPBE::Tabi::TreeWorkspace tree = make_tree(panels, executor, leaf_mac);
  TitanPBE::Tabi::TreecodeStats stats;
  TitanPBE::Tabi::apply_jcp_treecode_operator(
      panels, unknowns, 80.0, 0.10, 0.0, 0.0, true, tree, "centroid", "off",
      executor, out, &stats);
  return stats;
}

TitanPBE::Tabi::TreecodeStats apply_jcp_forced_leaf_split(
    const std::vector<TitanPBE::Tabi::Panel> &panels,
    const std::vector<double> &unknowns, TitanExecutor &executor,
    std::vector<double> &out) {
  TitanPBE::Tabi::TreeWorkspace tree = make_tree(panels, executor, "on");
  tree.convergence_estimator = "hybrid";
  tree.convergence_estimator_model = "channel-tail-eps";
  tree.convergence_safety = "split";
  tree.convergence_estimator_error_cap = 0.0;
  tree.convergence_unstable_rho_cap = 0.0;
  tree.convergence_split_selectivity = "legacy";
  tree.convergence_k3_accepted_rho_cap = 0.0;
  tree.expansion_audit = "sampled";
  tree.expansion_audit_channel = "sampled-detail";
  TitanPBE::Tabi::TreecodeStats stats;
  TitanPBE::Tabi::apply_jcp_treecode_operator(
      panels, unknowns, 80.0, 0.10, 0.0, 0.0, true, tree, "centroid", "off",
      executor, out, &stats);
  return stats;
}

TitanPBE::Tabi::TreecodeStats apply_non_jcp(
    const std::vector<TitanPBE::Tabi::Panel> &panels,
    const std::vector<double> &density, TitanExecutor &executor,
    const std::string &leaf_mac, std::vector<double> &out) {
  TitanPBE::Tabi::TreeWorkspace tree = make_tree(panels, executor, leaf_mac);
  TitanPBE::Tabi::TreecodeStats stats;
  TitanPBE::Tabi::apply_treecode_operator(panels, density, 0.5, 0.10, tree,
                                          executor, out, &stats);
  return stats;
}

} // namespace

int main() {
  TitanExecutor executor(1);
  const std::vector<TitanPBE::Tabi::Panel> panels = build_two_cluster_fixture();
  const std::vector<double> unknowns = make_unknowns(panels.size());
  const std::vector<double> density = make_density(panels.size());

  std::vector<double> jcp_off;
  std::vector<double> jcp_audit;
  std::vector<double> jcp_on;
  const auto jcp_off_stats = apply_jcp(panels, unknowns, executor, "off", jcp_off);
  const auto jcp_audit_stats =
      apply_jcp(panels, unknowns, executor, "audit", jcp_audit);
  const auto jcp_on_stats = apply_jcp(panels, unknowns, executor, "on", jcp_on);

  require_true("audit keeps JCP output unchanged",
               relative_rms(jcp_off, jcp_audit) <= 1.0e-14);
  require_true("leaf-MAC JCP remains close to off/direct traversal",
               relative_rms(jcp_off, jcp_on) <= 5.0e-3);
  require_true("audit sees leaf-MAC candidates",
               jcp_audit_stats.leaf_mac_candidates > 0);
  require_true("audit sees leaf-MAC accepted leaves",
               jcp_audit_stats.leaf_mac_accepted > 0);
  require_true("audit does not convert leaves to far nodes",
               jcp_audit_stats.leaf_mac_far_panels == 0);
  require_true("on sees leaf-MAC candidates",
               jcp_on_stats.leaf_mac_candidates > 0);
  require_true("on converts accepted leaves to far panels",
               jcp_on_stats.leaf_mac_far_panels > 0);
  require_true("on reduces JCP leaf work",
               jcp_on_stats.leaf_panel_interactions <
                   jcp_off_stats.leaf_panel_interactions);
  require_true("on estimates saved direct panels",
               jcp_on_stats.leaf_mac_direct_panels_saved_estimate > 0);

  std::vector<double> jcp_forced_split;
  const auto jcp_forced_split_stats =
      apply_jcp_forced_leaf_split(panels, unknowns, executor,
                                  jcp_forced_split);
  require_true("forced split sees convergence split rejects",
               jcp_forced_split_stats.convergence_estimator_split_rejects > 0);
  require_true("leaf split fallback avoids leaf far accounting",
               jcp_forced_split_stats.leaf_mac_far_panels == 0);
  require_true("leaf split fallback evaluates direct leaf work",
               jcp_forced_split_stats.leaf_panel_interactions > 0);
  require_true("leaf split fallback avoids far expansion audit",
               jcp_forced_split_stats.far_expansion_audit_pairs == 0);

  std::vector<double> non_jcp_off;
  std::vector<double> non_jcp_audit;
  std::vector<double> non_jcp_on;
  const auto non_jcp_off_stats =
      apply_non_jcp(panels, density, executor, "off", non_jcp_off);
  const auto non_jcp_audit_stats =
      apply_non_jcp(panels, density, executor, "audit", non_jcp_audit);
  const auto non_jcp_on_stats =
      apply_non_jcp(panels, density, executor, "on", non_jcp_on);
  require_true("audit keeps non-JCP output unchanged",
               relative_rms(non_jcp_off, non_jcp_audit) <= 1.0e-14);
  require_true("leaf-MAC non-JCP remains close to off/direct traversal",
               relative_rms(non_jcp_off, non_jcp_on) <= 5.0e-3);
  require_true("non-JCP audit sees leaf-MAC candidates",
               non_jcp_audit_stats.leaf_mac_candidates > 0);
  require_true("non-JCP on converts accepted leaves",
               non_jcp_on_stats.leaf_mac_far_panels > 0);
  require_true("non-JCP on reduces near interactions",
               non_jcp_on_stats.near_interactions <
                   non_jcp_off_stats.near_interactions);

  TitanPBE::Tabi::TreeWorkspace probe_tree = make_tree(panels, executor, "on");
  const std::vector<std::size_t> rank = make_rank(probe_tree, panels.size());
  const TitanPBE::Tabi::detail::PackedJcpMoments moments =
      TitanPBE::Tabi::detail::compute_jcp_node_moments_packed(
          panels, unknowns, true, probe_tree, probe_tree.taylor_order,
          executor);
  require_true("JCP convergence norm cache built",
               moments.convergence_norms.size() == probe_tree.nodes.size());
  const std::size_t probe_target = 0;
  const std::size_t source_leaf =
      find_direct_candidate_leaf(panels, probe_tree, probe_target);
  const TitanPBE::Tabi::TreeNode &source_node = probe_tree.nodes[source_leaf];
  const double dx = panels[probe_target].x - source_node.cx;
  const double dy = panels[probe_target].y - source_node.cy;
  const double dz = panels[probe_target].z - source_node.cz;
  const double distance =
      std::sqrt(std::max(1.0e-16, dx * dx + dy * dy + dz * dz));
  const TitanPBE::Tabi::detail::MacDecision source_mac =
      TitanPBE::Tabi::detail::evaluate_mac(panels[probe_target], source_node,
                                           probe_tree, distance);
  TitanPBE::Tabi::detail::PackedJcpMoments fallback_moments = moments;
  fallback_moments.convergence_norms.clear();
  require_same_estimator(
      TitanPBE::Tabi::detail::evaluate_convergence_estimator(
          probe_tree, source_node, moments, source_leaf,
          probe_tree.taylor_order, true, source_mac.rho, 80.0),
      TitanPBE::Tabi::detail::evaluate_convergence_estimator(
          probe_tree, source_node, fallback_moments, source_leaf,
          probe_tree.taylor_order, true, source_mac.rho, 80.0));
  std::size_t visits = 0;
  require_true("ancestor directness sees Leaf-MAC far leaf",
               !TitanPBE::Tabi::detail::
                   would_visit_leaf_direct_for_target_ancestor_path(
                       panels, probe_tree, rank, moments, probe_target,
                       source_leaf, probe_tree.taylor_order, true, 80.0,
                       &visits));
  probe_tree.leaf_mac = "audit";
  require_true("ancestor directness keeps audit leaf direct",
               TitanPBE::Tabi::detail::
                   would_visit_leaf_direct_for_target_ancestor_path(
                       panels, probe_tree, rank, moments, probe_target,
                       source_leaf, probe_tree.taylor_order, true, 80.0,
                       &visits));

  std::cout << "PBE TABI leaf-MAC smoke PASS panels=" << panels.size()
            << " jcp_rel=" << relative_rms(jcp_off, jcp_on)
            << " jcp_saved="
            << jcp_on_stats.leaf_mac_direct_panels_saved_estimate
            << " non_jcp_rel=" << relative_rms(non_jcp_off, non_jcp_on)
            << " non_jcp_saved="
            << non_jcp_on_stats.leaf_mac_direct_panels_saved_estimate
            << "\n";
  return 0;
}
'@

Set-Content -LiteralPath $SourcePath -Value $source -Encoding ASCII

$compileArgs = @(
    '-std=c++23',
    '-O2',
    '-mavx2',
    '-mfma',
    '-DTITAN_RUNTIME_DISPATCH_HOTSPOTS=1',
    ("-I" + (Join-Path $RepoRoot "src")),
    (Join-Path $RepoRoot 'src\tabi\tabi_near.cpp'),
    (Get-PbeTabiTreecodeSources -RepoRoot $RepoRoot),
    (Join-Path $RepoRoot 'src\titan_executor.cpp'),
    $SourcePath,
    '-o',
    $ExePath
)
& $Compiler @compileArgs
if ($LASTEXITCODE -ne 0) {
    throw "Failed to compile TABI leaf-MAC smoke fixture."
}

& $ExePath
if ($LASTEXITCODE -ne 0) {
    throw "TABI leaf-MAC smoke fixture failed."
}
