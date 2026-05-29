param(
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_leaf_work_attribution_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot 'pbe_tabi_treecode_sources.ps1')
$Compiler = 'C:\msys64\ucrt64\bin\clang++.exe'
$SourcePath = Join-Path $WorkRoot 'pbe_tabi_leaf_work_attribution_smoke.cpp'
$ExePath = Join-Path $WorkRoot 'pbe_tabi_leaf_work_attribution_smoke.exe'

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside leaf-work attribution smoke root: $full"
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

TitanPBE::Tabi::Panel panel(double x, double y, double z, std::size_t index) {
  TitanPBE::Tabi::Panel p;
  p.x = x;
  p.y = y;
  p.z = z;
  const double nx = 0.18 + 0.01 * static_cast<double>(index % 17);
  const double ny = (index % 2 == 0) ? 0.60 : -0.52;
  const double nz = 0.72;
  const double norm = std::sqrt(nx * nx + ny * ny + nz * nz);
  p.nx = nx / norm;
  p.ny = ny / norm;
  p.nz = nz / norm;
  p.area = 0.13 + 0.01 * static_cast<double>(index % 9);
  p.atom_index = index;
  return p;
}

std::vector<TitanPBE::Tabi::Panel> build_fixture() {
  std::vector<TitanPBE::Tabi::Panel> panels;
  panels.reserve(192);
  std::size_t index = 0;
  const double centers[6][3] = {{-18.0, 0.0, 0.0}, {18.0, 0.0, 0.0},
                                {0.0, -18.0, 0.0}, {0.0, 18.0, 0.0},
                                {-12.0, 12.0, 7.0}, {12.0, -12.0, -7.0}};
  for (const auto &center : centers) {
    for (int i = 0; i < 32; ++i) {
      const double a = 0.19634954084936207 * static_cast<double>(i);
      const double r = 1.0 + 0.05 * static_cast<double>(i % 5);
      panels.push_back(panel(center[0] + r * std::cos(a),
                             center[1] + r * std::sin(a),
                             center[2] + 0.12 * std::cos(2.0 * a),
                             index++));
    }
  }
  return panels;
}

std::vector<double> unknowns(std::size_t n) {
  std::vector<double> values(2 * n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    values[i] = 0.7 + 0.20 * std::sin(0.19 * static_cast<double>(i));
    values[n + i] = -0.9 + 0.17 * std::cos(0.23 * static_cast<double>(i));
  }
  return values;
}

TitanPBE::Tabi::TreecodeStats run_case(
    const std::vector<TitanPBE::Tabi::Panel> &panels,
    const std::vector<double> &values, TitanExecutor &executor,
    const std::string &profile, std::vector<double> &out) {
  TitanPBE::Tabi::TreeWorkspace tree =
      TitanPBE::Tabi::build_tree_workspace(
          panels, executor, 32, 0.55, "exact", "tile-avx2", "strict", "off",
          "onthefly", "tight-aabb", "topdown-aabb", "physical-gap", 1.5,
          "off", "area-centroid", "split", "off", 128, 2048, "translated",
          "naive", "off", "split", 0.0, 0.0, "hybrid", 1.0e4, 1.0e-3,
          0.25, "consecutive-aspect", "off", "on", 4, 64, 4,
          "channel-tail-eps", "off", "polynomial-eps", "accepted-only",
          "global", "off", 0.0, "surface-margin", 0.25, "risk-scaled", 3.0,
          0.55, 1.0e-7, "auto", "summary", profile, "baseline", "off",
          "off", "same-target-leaf", 32, "ancestor-path",
          "target-leaf-source-leaf", "off", "target-leaf-block", 2048,
          "off", "off", 2.0, 0.45, 0.10);
  tree.taylor_order = 4;
  tree.derivative_mode = "analytic";
  TitanPBE::Tabi::TreecodeStats stats;
  TitanPBE::Tabi::apply_jcp_treecode_operator(
      panels, values, 80.0, 0.10, 0.0, 0.0, true, tree, "centroid", "off",
      executor, out, &stats);
  return stats;
}

std::size_t aggregate_bucket_total(const TitanPBE::Tabi::TreecodeStats &s) {
  return s.leaf_work_mac_terminal_panels + s.leaf_work_physical_gap_panels +
         s.leaf_work_fragmentation_panels + s.leaf_work_convergence_panels +
         s.leaf_work_surface_margin_panels + s.leaf_work_k3_rho_cap_panels +
         s.leaf_work_singular_fallback_panels + s.leaf_work_scalar_tail_panels;
}

} // namespace

int main() {
  TitanExecutor executor(2);
  const std::vector<TitanPBE::Tabi::Panel> panels = build_fixture();
  const std::vector<double> values = unknowns(panels.size());

  std::vector<double> off_out;
  std::vector<double> summary_out;
  const TitanPBE::Tabi::TreecodeStats off =
      run_case(panels, values, executor, "off", off_out);
  const TitanPBE::Tabi::TreecodeStats summary =
      run_case(panels, values, executor, "summary", summary_out);

  require_true("profile off remains off", off.leaf_work_profile == "off");
  require_true("summary profile reported", summary.leaf_work_profile == "summary");
  require_true("profile does not perturb matvec",
               relative_rms(off_out, summary_out) <= 1.0e-12);
  require_true("leaf work is nonzero", summary.leaf_panel_interactions > 0);
  require_true("aggregate buckets equal accounted panels",
               aggregate_bucket_total(summary) ==
                   summary.leaf_work_profile_accounted_panels);
  require_true("accounted panels equal leaf interactions",
               summary.leaf_work_profile_accounted_panels ==
                   summary.leaf_panel_interactions);
  require_true("convergence sub-buckets do not exceed aggregate",
               summary.leaf_work_convergence_estimator_panels +
                       summary.leaf_work_convergence_rho_panels +
                       summary.leaf_work_convergence_aspect_panels <=
                   summary.leaf_work_convergence_panels);
  require_true("surface-margin sub-buckets do not exceed aggregate",
               summary.leaf_work_surface_margin_hard_panels +
                       summary.leaf_work_surface_margin_soft_panels <=
                   summary.leaf_work_surface_margin_panels);
  const double fraction_sum =
      summary.leaf_work_mac_terminal_frac + summary.leaf_work_convergence_frac +
      summary.leaf_work_surface_margin_frac +
      summary.leaf_work_singular_fallback_frac +
      summary.leaf_work_scalar_tail_frac;
  require_true("reported fractions are bounded",
               fraction_sum >= 0.0 && fraction_sum <= 1.0 + 1.0e-12);

  std::cout << "PBE TABI leaf-work attribution smoke PASS leaf="
            << summary.leaf_panel_interactions
            << " convergence=" << summary.leaf_work_convergence_panels
            << " margin=" << summary.leaf_work_surface_margin_panels
            << " frac_sum=" << fraction_sum << "\n";
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
    throw "Failed to compile TABI leaf-work attribution smoke fixture."
}

& $ExePath
if ($LASTEXITCODE -ne 0) {
    throw "TABI leaf-work attribution smoke fixture failed."
}
