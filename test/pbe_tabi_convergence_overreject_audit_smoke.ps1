param(
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_convergence_overreject_audit_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
$Compiler = 'C:\msys64\ucrt64\bin\clang++.exe'
$SourcePath = Join-Path $WorkRoot 'pbe_tabi_convergence_overreject_audit_smoke.cpp'
$ExePath = Join-Path $WorkRoot 'pbe_tabi_convergence_overreject_audit_smoke.exe'

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside overreject audit smoke root: $full"
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

#include <iostream>
#include <stdexcept>

void require_true(const char *label, bool value) {
  if (!value)
    throw std::runtime_error(label);
}

int main() {
  using TitanPBE::Tabi::detail::classify_convergence_overreject_audit;

  const auto low_error =
      classify_convergence_overreject_audit(true, 2.0e-3, 1.0e-9, 1.0e-3,
                                            1.0e-7);
  require_true("split-rejected low absolute error is overreject",
               low_error.low_error_overreject &&
                   !low_error.caught_high_error);

  const auto high_error =
      classify_convergence_overreject_audit(true, 2.0e-3, 1.0e-5, 1.0e-3,
                                            1.0e-7);
  require_true("split-rejected real high error is caught",
               high_error.caught_high_error &&
                   !high_error.low_error_overreject);

  const auto accepted =
      classify_convergence_overreject_audit(false, 2.0e-3, 1.0e-5, 1.0e-3,
                                            1.0e-7);
  require_true("accepted far candidate is not an overreject audit sample",
               !accepted.caught_high_error &&
                   !accepted.low_error_overreject);

  std::cout << "convergence overreject audit smoke passed\n";
  return 0;
}
'@

Set-Content -LiteralPath $SourcePath -Value $source -Encoding ASCII

& $Compiler -std=c++23 -O2 -I"$RepoRoot" `
    $SourcePath (Join-Path $RepoRoot 'titan_pbe\tabi_treecode_internal.cpp') `
    -o $ExePath
if ($LASTEXITCODE -ne 0) {
    throw "Failed to compile convergence overreject audit smoke."
}

& $ExePath
if ($LASTEXITCODE -ne 0) {
    throw "Convergence overreject audit smoke failed."
}
