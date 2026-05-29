param(
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_absolute_aware_convergence_split_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
$Compiler = 'C:\msys64\ucrt64\bin\clang++.exe'
$SourcePath = Join-Path $WorkRoot 'pbe_tabi_absolute_aware_convergence_split_smoke.cpp'
$ExePath = Join-Path $WorkRoot 'pbe_tabi_absolute_aware_convergence_split_smoke.exe'

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside absolute-aware split smoke root: $full"
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

namespace {

void require_true(const char *label, bool value) {
  if (!value)
    throw std::runtime_error(label);
}

TitanPBE::Tabi::detail::ConvergenceSplitSelectivityInput base_input() {
  TitanPBE::Tabi::detail::ConvergenceSplitSelectivityInput input;
  input.unstable = true;
  input.safety_split = true;
  input.rho = 0.60;
  input.unstable_rho_cap = 0.25;
  input.estimator_model = "channel-tail-eps";
  input.channel_unstable[TitanPBE::Tabi::detail::kJcpK3PhinChannel] = true;
  return input;
}

} // namespace

int main() {
  using TitanPBE::Tabi::detail::evaluate_convergence_split_selectivity;
  using TitanPBE::Tabi::detail::kJcpK3PhinChannel;
  using TitanPBE::Tabi::detail::kJcpK4MixedChannel;

  auto legacy = base_input();
  legacy.split_selectivity = "legacy";
  legacy.split_absolute_floor = 1.0e-7;
  legacy.estimated_abs_error = 1.0e-12;
  require_true("legacy keeps split regardless of absolute floor",
               evaluate_convergence_split_selectivity(legacy).split_reject);

  auto low_abs = base_input();
  low_abs.split_selectivity = "absolute-aware";
  low_abs.split_absolute_floor = 1.0e-7;
  low_abs.estimated_abs_error = 2.0e-8;
  const auto low_result = evaluate_convergence_split_selectivity(low_abs);
  require_true("absolute-aware reaccepts low absolute error",
               !low_result.split_reject && low_result.low_abs_reaccepted);

  auto high_abs = base_input();
  high_abs.split_selectivity = "absolute-aware";
  high_abs.split_absolute_floor = 1.0e-7;
  high_abs.estimated_abs_error = 2.0e-6;
  require_true("absolute-aware still splits real high absolute error",
               evaluate_convergence_split_selectivity(high_abs).split_reject);

  auto channel = base_input();
  channel.split_selectivity = "channel-absolute-aware";
  channel.split_absolute_floor = 1.0e-7;
  channel.channel_unstable[kJcpK4MixedChannel] = true;
  channel.channel_estimated_abs_error[kJcpK3PhinChannel] = 1.0e-8;
  channel.channel_estimated_abs_error[kJcpK4MixedChannel] = 2.0e-6;
  const auto channel_result =
      evaluate_convergence_split_selectivity(channel);
  require_true("channel absolute-aware splits because K4 is real high error",
               channel_result.split_reject);
  require_true("channel absolute-aware marks K3 low-abs reaccepted",
               channel_result.channel_low_abs_reaccepted[kJcpK3PhinChannel]);

  auto k3_hot = base_input();
  k3_hot.split_policy = "k3-hot";
  k3_hot.split_selectivity = "channel-absolute-aware";
  k3_hot.split_absolute_floor = 1.0e-7;
  k3_hot.channel_unstable[kJcpK3PhinChannel] = false;
  k3_hot.channel_unstable[kJcpK4MixedChannel] = true;
  k3_hot.channel_estimated_abs_error[kJcpK4MixedChannel] = 1.0e-4;
  require_true("k3-hot policy does not split non-K3-only instability",
               !evaluate_convergence_split_selectivity(k3_hot).split_reject);

  std::cout << "absolute-aware convergence split smoke passed\n";
  return 0;
}
'@

Set-Content -LiteralPath $SourcePath -Value $source -Encoding ASCII

& $Compiler -std=c++23 -O2 $("-I" + (Join-Path $RepoRoot "src")) `
    $SourcePath (Join-Path $RepoRoot 'src\tabi\tabi_treecode_internal.cpp') `
    -o $ExePath
if ($LASTEXITCODE -ne 0) {
    throw "Failed to compile absolute-aware convergence split smoke."
}

& $ExePath
if ($LASTEXITCODE -ne 0) {
    throw "Absolute-aware convergence split smoke failed."
}
