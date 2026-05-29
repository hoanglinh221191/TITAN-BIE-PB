param(
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_convergence_shape_audit_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot 'pbe_tabi_treecode_sources.ps1')
$Compiler = 'C:\msys64\ucrt64\bin\clang++.exe'
$SourcePath = Join-Path $WorkRoot 'pbe_tabi_convergence_shape_audit_smoke.cpp'
$ExePath = Join-Path $WorkRoot 'pbe_tabi_convergence_shape_audit_smoke.exe'

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside shape audit smoke root: $full"
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
#include "tabi/tabi_treecode_internal.h"

#include <iostream>
#include <stdexcept>
#include <string>

namespace {

using namespace TitanPBE::Tabi;
using namespace TitanPBE::Tabi::detail;

void require_true(const std::string &label, bool value) {
  if (!value)
    throw std::runtime_error(label);
}

TreeNode make_node() {
  TreeNode node;
  node.min_x = -2.0;
  node.min_y = -0.5;
  node.min_z = -0.25;
  node.max_x = 2.0;
  node.max_y = 0.5;
  node.max_z = 0.25;
  node.radius = 2.0766559657295187;
  node.panel_count = 32;
  return node;
}

ConvergenceEstimatorDecision make_estimator() {
  ConvergenceEstimatorDecision estimator;
  estimator.active = true;
  estimator.unstable = true;
  estimator.split_reject = true;
  estimator.source_norm = 4.0;
  estimator.tail_norm_p = 0.010;
  estimator.tail_norm_pm1 = 0.020;
  estimator.tail_ratio = 0.005;
  estimator.effective_jcp_amp = 120.0;
  estimator.estimated_error = 0.02;
  estimator.estimated_abs_error = 0.5;
  for (std::size_t channel = 0; channel < kJcpExpansionChannelCount;
       ++channel) {
    estimator.channel_source_norm[channel] = 1.0 + static_cast<double>(channel);
    estimator.channel_tail_norm_p[channel] = 0.002 * (1.0 + channel);
    estimator.channel_tail_norm_pm1[channel] = 0.003 * (1.0 + channel);
    estimator.channel_tail_ratio[channel] = 0.01 * (1.0 + channel);
    estimator.channel_effective_jcp_amp[channel] = 10.0 * (1.0 + channel);
    estimator.channel_estimated_error[channel] = 0.02 * (1.0 + channel);
    estimator.channel_estimated_abs_error[channel] = 0.1 * (1.0 + channel);
    estimator.channel_unstable[channel] = true;
  }
  return estimator;
}

JcpChannelAccum repeat_channel(const JcpAccum &value) {
  JcpChannelAccum channels;
  for (std::size_t channel = 0; channel < kJcpExpansionChannelCount;
       ++channel)
    channels.channel[channel] = value;
  return channels;
}

void record_pair(ExpansionAuditTally &tally, const JcpAccum &taylor,
                 const JcpAccum &direct, const std::string &shape_mode) {
  const TreeNode node = make_node();
  const ConvergenceEstimatorDecision estimator = make_estimator();
  const JcpChannelAccum taylor_channels = repeat_channel(taylor);
  const JcpChannelAccum direct_channels = repeat_channel(direct);
  record_convergence_overreject_audit(
      tally, taylor, direct, taylor_channels, direct_channels, estimator, node,
      false, 1.0e-3, 1.0e-7, "off", shape_mode, "off", 5, 9, 0.50, 7.5,
      4, JcpChannelAccum{}, JcpChannelAccum{});
}

} // namespace

int main() {
  const JcpAccum tiny_direct{1.0e-6, 0.0};
  const JcpAccum tiny_taylor{1.002e-6, 0.0};
  const JcpAccum real_direct{1.0e-3, 0.0};
  const JcpAccum real_taylor{1.02e-3, 0.0};

  ExpansionAuditTally off_mode;
  ExpansionAuditTally detail_mode;
  record_pair(off_mode, tiny_taylor, tiny_direct, "off");
  record_pair(detail_mode, tiny_taylor, tiny_direct, "sampled-detail");
  record_pair(detail_mode, real_taylor, real_direct, "sampled-detail");

  require_true("shape audit off preserves base overreject count",
               off_mode.convergence_overreject_count == 1);
  require_true("shape detail records caught and overreject rows",
               detail_mode.convergence_shape_records.size() >= 2);
  require_true("first shape row keeps estimator features",
               detail_mode.convergence_shape_records[0].source_norm > 0.0 &&
                   detail_mode.convergence_shape_records[0].tail_norm_p >
                       0.0 &&
                   detail_mode.convergence_shape_records[0].rho_power > 0.0);
  bool saw_caught = false;
  bool saw_overreject = false;
  for (const auto &record : detail_mode.convergence_shape_records) {
    saw_caught = saw_caught || record.caught;
    saw_overreject = saw_overreject || !record.caught;
  }
  require_true("shape records contain both classes", saw_caught && saw_overreject);

  write_convergence_shape_calibration_csv(detail_mode);
  std::cout << "PBE TABI convergence shape audit smoke PASS records="
            << detail_mode.convergence_shape_records.size() << "\n";
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
    throw "Failed to compile TABI convergence shape audit smoke fixture."
}

Push-Location $WorkRoot
try {
    & $ExePath
    if ($LASTEXITCODE -ne 0) {
        throw "TABI convergence shape audit smoke fixture failed."
    }
} finally {
    Pop-Location
}

$ShapeCsv = Join-Path $WorkRoot 'pbe_tabi_convergence_shape_calibration.csv'
if (-not (Test-Path -LiteralPath $ShapeCsv)) {
    throw "Shape calibration CSV was not emitted: $ShapeCsv"
}
$header = (Get-Content -LiteralPath $ShapeCsv -TotalCount 1)
foreach ($column in @('source_norm', 'tail_norm_p', 'tail_norm_pm1', 'oracle_abs_over_tail_rho', 'classification')) {
    if ($header -notmatch "(^|,)$column(,|$)") {
        throw "Shape calibration CSV is missing column $column"
    }
}
$rows = @(Import-Csv -LiteralPath $ShapeCsv)
if (-not ($rows | Where-Object { $_.classification -eq 'caught' }) -or
    -not ($rows | Where-Object { $_.classification -eq 'overreject' })) {
    throw "Shape calibration CSV does not contain both caught and overreject rows."
}
