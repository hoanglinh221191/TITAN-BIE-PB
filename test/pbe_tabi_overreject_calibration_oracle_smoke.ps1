param(
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_overreject_calibration_oracle_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
. (Join-Path $PSScriptRoot 'pbe_tabi_treecode_sources.ps1')
$Compiler = 'C:\msys64\ucrt64\bin\clang++.exe'
$SourcePath = Join-Path $WorkRoot 'pbe_tabi_overreject_calibration_oracle_smoke.cpp'
$ExePath = Join-Path $WorkRoot 'pbe_tabi_overreject_calibration_oracle_smoke.exe'

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside calibration smoke root: $full"
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
#include "titan_pbe/tabi_treecode_internal.h"

#include <cmath>
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
  node.min_x = -3.0;
  node.min_y = -0.5;
  node.min_z = -0.25;
  node.max_x = 3.0;
  node.max_y = 0.5;
  node.max_z = 0.25;
  node.radius = 3.0516389039334255;
  node.panel_count = 24;
  return node;
}

ConvergenceEstimatorDecision make_estimator() {
  ConvergenceEstimatorDecision estimator;
  estimator.active = true;
  estimator.unstable = true;
  estimator.split_reject = true;
  estimator.estimated_error = 2.5e-2;
  estimator.estimated_abs_error = 1.0e-2;
  for (std::size_t channel = 0; channel < kJcpExpansionChannelCount;
       ++channel) {
    estimator.channel_estimated_error[channel] = 2.5e-2;
    estimator.channel_estimated_abs_error[channel] = 1.0e-2;
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
                 const JcpAccum &direct, const std::string &detail_mode) {
  const TreeNode node = make_node();
  const ConvergenceEstimatorDecision estimator = make_estimator();
  const JcpChannelAccum taylor_channels = repeat_channel(taylor);
  const JcpChannelAccum direct_channels = repeat_channel(direct);
  record_convergence_overreject_audit(
      tally, taylor, direct, taylor_channels, direct_channels, estimator, node,
      false, 1.0e-3, 1.0e-7, detail_mode, detail_mode, "off", 7, 11,
      0.55, 8.0, 4, JcpChannelAccum{}, JcpChannelAccum{});
}

} // namespace

int main() {
  const JcpAccum tiny_direct{1.0e-6, 0.0};
  const JcpAccum tiny_taylor{1.002e-6, 0.0};
  ExpansionAuditTally overreject;
  record_pair(overreject, tiny_taylor, tiny_direct, "sampled-detail");
  require_true("low absolute split rejection classified overreject",
               overreject.convergence_overreject_count == 1 &&
                   overreject.convergence_overreject_caught_count == 0);
  require_true("overreject stores oracle and scale samples",
               overreject.convergence_overreject_abs_error_samples.size() ==
                       1 &&
                   overreject.convergence_overreject_required_scale_samples
                           .size() == 1);
  require_true("required scale is small",
               overreject.convergence_overreject_required_scale_samples[0] <
                   1.0e-5);
  require_true("detail row marks overreject",
               !overreject.overreject_calibration_records.empty() &&
                   !overreject.overreject_calibration_records[0].caught);

  const JcpAccum real_direct{1.0e-3, 0.0};
  const JcpAccum real_taylor{1.02e-3, 0.0};
  ExpansionAuditTally caught;
  record_pair(caught, real_taylor, real_direct, "sampled-detail");
  require_true("high absolute split rejection classified caught",
               caught.convergence_overreject_count == 0 &&
                   caught.convergence_overreject_caught_count == 1);
  require_true("detail row marks caught",
               !caught.overreject_calibration_records.empty() &&
                   caught.overreject_calibration_records[0].caught);

  ExpansionAuditTally off_mode;
  ExpansionAuditTally summary_mode;
  record_pair(off_mode, tiny_taylor, tiny_direct, "off");
  record_pair(summary_mode, tiny_taylor, tiny_direct, "summary");
  require_true("detail mode preserves base classification counters",
               off_mode.convergence_overreject_count ==
                       summary_mode.convergence_overreject_count &&
                   off_mode.convergence_overreject_caught_count ==
                       summary_mode.convergence_overreject_caught_count);
  require_true("detail off stores no new samples",
               off_mode.convergence_overreject_abs_error_samples.empty());
  require_true("summary stores samples without raw detail rows",
               summary_mode.convergence_overreject_abs_error_samples.size() ==
                       1 &&
                   summary_mode.overreject_calibration_records.empty());

  write_overreject_calibration_detail_csv(overreject);
  std::cout << "PBE TABI overreject calibration oracle smoke PASS"
            << " required_scale="
            << overreject.convergence_overreject_required_scale_samples[0]
            << " caught=" << caught.convergence_overreject_caught_count
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
    "-I$RepoRoot",
    (Join-Path $RepoRoot 'titan_pbe\tabi_near.cpp'),
    (Get-PbeTabiTreecodeSources -RepoRoot $RepoRoot),
    (Join-Path $RepoRoot 'titan_executor.cpp'),
    $SourcePath,
    '-o',
    $ExePath
)
& $Compiler @compileArgs
if ($LASTEXITCODE -ne 0) {
    throw "Failed to compile TABI overreject calibration oracle smoke fixture."
}

Push-Location $WorkRoot
try {
    & $ExePath
    if ($LASTEXITCODE -ne 0) {
        throw "TABI overreject calibration oracle smoke fixture failed."
    }
} finally {
    Pop-Location
}

$DetailCsv = Join-Path $WorkRoot 'pbe_tabi_convergence_overreject_calibration.csv'
if (-not (Test-Path -LiteralPath $DetailCsv)) {
    throw "Calibration detail CSV was not emitted: $DetailCsv"
}
$rows = @(Import-Csv -LiteralPath $DetailCsv)
if ($rows.Count -lt 1 -or $rows[0].classification -ne 'overreject') {
    throw "Calibration detail CSV does not contain the expected overreject row."
}
