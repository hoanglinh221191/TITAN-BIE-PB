param(
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_treecode_accuracy')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot 'pbe_tabi_treecode_sources.ps1')
$Compiler = 'C:\msys64\ucrt64\bin\clang++.exe'
$SourcePath = Join-Path $WorkRoot 'pbe_tabi_treecode_accuracy.cpp'
$ExePath = Join-Path $WorkRoot 'pbe_tabi_treecode_accuracy.exe'

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside TABI treecode accuracy root: $full"
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
  return values.empty() ? 0.0 : std::sqrt(sum / static_cast<double>(values.size()));
}

double relative_rms(const std::vector<double> &reference,
                    const std::vector<double> &candidate) {
  std::vector<double> diff(reference.size(), 0.0);
  for (std::size_t i = 0; i < reference.size(); ++i)
    diff[i] = candidate[i] - reference[i];
  return rms(diff) / std::max(1.0e-12, rms(reference));
}

struct Vec3 {
  double x;
  double y;
  double z;
};

Vec3 sub(Vec3 a, Vec3 b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }

Vec3 cross(Vec3 a, Vec3 b) {
  return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z,
          a.x * b.y - a.y * b.x};
}

double norm(Vec3 a) { return std::sqrt(a.x * a.x + a.y * a.y + a.z * a.z); }

TitanPBE::Tabi::Panel make_panel(Vec3 a, Vec3 b, Vec3 c, std::size_t ia,
                                 std::size_t ib, std::size_t ic) {
  TitanPBE::Tabi::Panel p;
  p.v0x = a.x;
  p.v0y = a.y;
  p.v0z = a.z;
  p.v1x = b.x;
  p.v1y = b.y;
  p.v1z = b.z;
  p.v2x = c.x;
  p.v2y = c.y;
  p.v2z = c.z;
  p.v0_index = ia;
  p.v1_index = ib;
  p.v2_index = ic;
  p.x = (a.x + b.x + c.x) / 3.0;
  p.y = (a.y + b.y + c.y) / 3.0;
  p.z = (a.z + b.z + c.z) / 3.0;
  const Vec3 normal = cross(sub(b, a), sub(c, a));
  const double nrm = norm(normal);
  require_true("manual panel normal", nrm > 0.0);
  p.area = 0.5 * nrm;
  p.nx = normal.x / nrm;
  p.ny = normal.y / nrm;
  p.nz = normal.z / nrm;
  p.has_triangle_geometry = true;
  return p;
}

std::vector<TitanPBE::Tabi::Panel> build_panels() {
  TitanPBE::Tabi::Options opt;
  opt.surface_mode = "vdw";
  opt.surface_geometry = "triangle";
  opt.diagonal_mode = "solid-angle";
  TitanPBE::Tabi::Result result;
  std::vector<TitanPBE::Tabi::Panel> panels =
      TitanPBE::Tabi::build_surface(
          {{-3.0, -0.2, 0.0, 0.4, 1.55},
           {-1.5, 0.3, 0.1, -0.2, 1.70},
           {0.0, -0.1, -0.3, 0.3, 1.65},
           {1.5, 0.2, 0.2, -0.4, 1.70},
           {3.0, -0.2, 0.0, 0.2, 1.55}},
          opt, &result);
  require_true("triangle panels", panels.size() > 900);
  require_true("solid-angle diagonal", result.diag_min > 0.0);
  return panels;
}

} // namespace

int main() {
  std::vector<TitanPBE::Tabi::Panel> panels = build_panels();
  std::vector<double> density(panels.size(), 0.0);
  for (std::size_t i = 0; i < density.size(); ++i)
    density[i] = std::sin(0.017 * static_cast<double>(i)) +
                 0.25 * std::cos(0.031 * static_cast<double>(i));

  TitanExecutor executor(2);
  std::vector<double> direct;
  TitanPBE::Tabi::TreecodeStats direct_stats;
  TitanPBE::Tabi::apply_direct_operator(panels, density, 0.42, 0.12,
                                        executor, direct, &direct_stats);
  require_true("direct matvec", direct_stats.matvec_mode == "direct");

  TitanPBE::Tabi::TreeWorkspace tree =
      TitanPBE::Tabi::build_tree_workspace(panels, executor);
  require_true("tree nodes", tree.nodes.size() > 1);

  std::vector<double> tree3;
  tree.taylor_order = 3;
  TitanPBE::Tabi::TreecodeStats stats3;
  TitanPBE::Tabi::apply_treecode_operator(panels, density, 0.42, 0.12, tree,
                                          executor, tree3, &stats3);
  const double rel3 = relative_rms(direct, tree3);
  require_true("treecode p3 mode", stats3.matvec_mode == "treecode");
  require_true("treecode p3 far interactions", stats3.far_interactions > 0);
  require_true("treecode p3 relative error", rel3 < 2.5e-3);

  std::vector<double> tree4;
  tree.taylor_order = 4;
  TitanPBE::Tabi::TreecodeStats stats4;
  TitanPBE::Tabi::apply_treecode_operator(panels, density, 0.42, 0.12, tree,
                                          executor, tree4, &stats4);
  const double rel4 = relative_rms(direct, tree4);
  require_true("treecode p4 order", stats4.taylor_order == 4);
  require_true("treecode p4 not worse", rel4 <= rel3 * 1.50 + 5.0e-5);

  TitanPBE::Tabi::TreeWorkspace topdown_tree =
      TitanPBE::Tabi::build_tree_workspace(
          panels, executor, 32, 0.45, "auto", "auto", "strict", "off",
          "cached", "moment", "topdown-aabb");
  require_true("topdown builder label",
               topdown_tree.tree_builder == "topdown-aabb");
  require_true("topdown child indirection",
               topdown_tree.child_indices.size() >=
                   topdown_tree.nodes.size() - topdown_tree.leaf_count);
  std::vector<double> topdown_tree3;
  topdown_tree.taylor_order = 3;
  TitanPBE::Tabi::TreecodeStats topdown_stats3;
  TitanPBE::Tabi::apply_treecode_operator(
      panels, density, 0.42, 0.12, topdown_tree, executor, topdown_tree3,
      &topdown_stats3);
  const double topdown_rel3 = relative_rms(direct, topdown_tree3);
  require_true("topdown treecode mode",
               topdown_stats3.matvec_mode == "treecode");
  require_true("topdown stats builder",
               topdown_stats3.tree_builder == "topdown-aabb");
  require_true("topdown relative error", topdown_rel3 < 2.5e-3);

  std::vector<double> unknowns(2 * panels.size(), 0.0);
  for (std::size_t i = 0; i < panels.size(); ++i) {
    unknowns[i] = 0.1 * std::sin(0.013 * static_cast<double>(i));
    unknowns[panels.size() + i] =
        0.05 * std::cos(0.019 * static_cast<double>(i));
  }
  std::vector<double> jcp_tree;
  TitanPBE::Tabi::TreecodeStats jcp_stats;
  TitanPBE::Tabi::apply_jcp_treecode_operator(
      panels, unknowns, 0.0125, 0.12, 0.50625, 40.5, true, tree,
      "duffy_adaptive", "off", executor, jcp_tree, &jcp_stats);
  require_true("jcp duffy near policy",
               jcp_stats.near_kernel_policy == "one_ring_cpv");
  require_true("jcp duffy one-ring pairs",
               jcp_stats.near_duffy_pairs > jcp_stats.near_self_pairs);
  require_true("jcp duffy no fallback", jcp_stats.near_fallbacks == 0);
  require_true("jcp duffy k1 decisions",
               jcp_stats.near_k1_accepted + jcp_stats.near_k1_rejected > 0);
  require_true("jcp duffy k2 decisions",
               jcp_stats.near_k2_accepted + jcp_stats.near_k2_rejected > 0);
  require_true("jcp duffy k3 decisions",
               jcp_stats.near_k3_accepted + jcp_stats.near_k3_rejected > 0);
  require_true("jcp duffy k4 decisions",
               jcp_stats.near_k4_accepted + jcp_stats.near_k4_rejected > 0);
  require_true("jcp duffy one-ring k1 decisions",
               jcp_stats.near_one_ring_k1_accepted +
                       jcp_stats.near_one_ring_k1_rejected >
                   0);
  require_true("jcp duffy one-ring k2 decisions",
               jcp_stats.near_one_ring_k2_accepted +
                       jcp_stats.near_one_ring_k2_rejected >
                   0);
  require_true("jcp duffy one-ring k3 decisions",
               jcp_stats.near_one_ring_k3_accepted +
                       jcp_stats.near_one_ring_k3_rejected >
                   0);
  require_true("jcp duffy one-ring k4 decisions",
               jcp_stats.near_one_ring_k4_accepted +
                       jcp_stats.near_one_ring_k4_rejected >
                   0);
  require_true("jcp tree output size", jcp_tree.size() == unknowns.size());

  std::vector<TitanPBE::Tabi::Panel> close_panels;
  close_panels.push_back(make_panel({0.0, 0.0, 0.0}, {1.0, 0.0, 0.0},
                                    {0.0, 1.0, 0.0}, 0, 1, 2));
  close_panels.push_back(make_panel({0.0, 0.0, 0.30}, {1.0, 0.0, 0.30},
                                    {0.0, 1.0, 0.30}, 10, 11, 12));
  std::vector<double> close_unknowns{0.25, -0.35, 0.50, -0.15};
  TitanPBE::Tabi::TreeWorkspace buffer_tree =
      TitanPBE::Tabi::build_tree_workspace(
          close_panels, executor, 8, 0.5, "exact", "tile", "strict",
          "adaptive");
  std::vector<double> buffer_out;
  TitanPBE::Tabi::TreecodeStats buffer_stats;
  TitanPBE::Tabi::apply_jcp_treecode_operator(
      close_panels, close_unknowns, 0.0125, 0.12, 0.50625, 40.5, true,
      buffer_tree, "duffy_adaptive", "adaptive", executor, buffer_out,
      &buffer_stats);
  require_true("treecode geometric buffer pair",
               buffer_stats.near_geometric_buffer_pairs > 0);
  require_true("treecode geometric pair is singular path",
               buffer_stats.singular_near_interactions >=
                   buffer_stats.near_geometric_buffer_pairs);

  std::cout << "PBE TABI treecode accuracy PASS panels=" << panels.size()
            << " rel3=" << rel3 << " rel4=" << rel4
            << " topdown_rel3=" << topdown_rel3
            << " far=" << stats3.far_interactions
            << " jcp_duffy=" << jcp_stats.near_duffy_pairs << "\n";
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
    (Join-Path $RepoRoot 'src\tabi\tabi_surface.cpp'),
    (Join-Path $RepoRoot 'src\tabi\tabi_near.cpp'),
    (Get-PbeTabiTreecodeSources -RepoRoot $RepoRoot),
    (Join-Path $RepoRoot 'src\titan_surface.cpp'),
    (Join-Path $RepoRoot 'src\titan_executor.cpp'),
    $SourcePath,
    '-o',
    $ExePath
)
& $Compiler @compileArgs
if ($LASTEXITCODE -ne 0) {
    throw "Failed to compile TABI treecode accuracy fixture."
}

& $ExePath
if ($LASTEXITCODE -ne 0) {
    throw "TABI treecode accuracy fixture failed."
}
