param(
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_surface_mesh_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
$Compiler = 'C:\msys64\ucrt64\bin\clang++.exe'
$SourcePath = Join-Path $WorkRoot 'pbe_tabi_surface_mesh_smoke.cpp'
$ExePath = Join-Path $WorkRoot 'pbe_tabi_surface_mesh_smoke.exe'

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside TABI surface mesh smoke root: $full"
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
#include "tabi/tabi_types.h"

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

void require_true(const std::string &label, bool value) {
  if (!value)
    throw std::runtime_error(label);
}

void require_close_rel(const std::string &label, double a, double b,
                       double rel_tol) {
  const double scale = std::max(1.0, std::abs(b));
  if (std::abs(a - b) > rel_tol * scale)
    throw std::runtime_error(label + " mismatch: " + std::to_string(a) +
                             " vs " + std::to_string(b));
}

TitanPBE::Tabi::Options options() {
  TitanPBE::Tabi::Options opt;
  opt.surface_mode = "vdw";
  opt.surface_geometry = "triangle";
  opt.diagonal_mode = "solid-angle";
  opt.probe_radius = 0.0;
  return opt;
}

TitanPBE::Tabi::Result build_case(
    const std::string &label,
    const std::vector<TitanPBE::Tabi::Atom> &atoms) {
  TitanPBE::Tabi::Result result;
  std::vector<TitanPBE::Tabi::Panel> panels =
      TitanPBE::Tabi::build_surface(atoms, options(), &result);
  require_true(label + " geometry", result.surface_geometry == "triangle");
  require_true(label + " panel count", !panels.empty());
  require_true(label + " mesh vertices", result.mesh_vertices > 0);
  require_true(label + " mesh triangles",
               result.mesh_triangles == panels.size());
  require_true(label + " boundary nodes",
               result.boundary_nodes == panels.size());
  require_true(label + " diagonal finite",
               std::isfinite(result.diag_min) &&
                   std::isfinite(result.diag_max));
  require_true(label + " diagonal positive", result.diag_min > 0.0);
  require_true(label + " fallback sane",
               result.solid_angle_fallbacks < panels.size() / 20 + 1);
  for (const auto &panel : panels) {
    const auto &atom = atoms[panel.atom_index];
    const double rx = panel.x - atom.x;
    const double ry = panel.y - atom.y;
    const double rz = panel.z - atom.z;
    require_true(label + " outward normal",
                 panel.nx * rx + panel.ny * ry + panel.nz * rz > 0.0);
    require_true(label + " positive area", panel.area > 0.0);
  }
  return result;
}

} // namespace

int main() {
  const double radius = 2.0;
  const auto single =
      build_case("single_sphere", {{0.0, 0.0, 0.0, 1.0, radius}});
  require_close_rel("single sphere area", single.surface_area_a2,
                    4.0 * kPi * radius * radius, 0.08);
  require_true("single sphere nonmanifold", single.nonmanifold_edges == 0);

  const auto two =
      build_case("two_overlapping_atoms",
                 {{0.0, 0.0, 0.0, 1.0, 1.7},
                  {2.0, 0.0, 0.0, -1.0, 1.7}});
  require_true("two atom buried triangles", two.buried_points > 0);
  require_true("two atom area reduced",
               two.surface_area_a2 < 2.0 * 4.0 * kPi * 1.7 * 1.7);

  const auto mini =
      build_case("mini_pqr_like",
                 {{-3.2, 0.0, 0.0, -0.3, 1.55},
                  {-1.8, 0.1, 0.0, 0.2, 1.70},
                  {-1.16, 0.5, -1.32, 0.5, 1.70},
                  {-1.72, 1.26, -2.10, -0.5, 1.52},
                  {0.34, 0.2, -1.40, 0.1, 1.60}});
  require_true("mini mesh size", mini.mesh_triangles > 500);

  TitanPBE::Tabi::Result repeat;
  std::vector<TitanPBE::Tabi::Panel> repeat_panels =
      TitanPBE::Tabi::build_surface(
          {{-3.2, 0.0, 0.0, -0.3, 1.55},
           {-1.8, 0.1, 0.0, 0.2, 1.70},
           {-1.16, 0.5, -1.32, 0.5, 1.70},
           {-1.72, 1.26, -2.10, -0.5, 1.52},
           {0.34, 0.2, -1.40, 0.1, 1.60}},
          options(), &repeat);
  require_true("repeat triangle count",
               repeat.mesh_triangles == mini.mesh_triangles);
  require_close_rel("repeat area", repeat.surface_area_a2,
                    mini.surface_area_a2, 1.0e-12);
  require_true("repeat panels", repeat_panels.size() == mini.boundary_nodes);

  std::cout << "PBE TABI surface mesh smoke PASS panels="
            << mini.boundary_nodes << " vertices=" << mini.mesh_vertices
            << " diag=[" << mini.diag_min << "," << mini.diag_max << "]\n";
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
    (Join-Path $RepoRoot 'src\titan_surface.cpp'),
    (Join-Path $RepoRoot 'src\titan_executor.cpp'),
    $SourcePath,
    '-o',
    $ExePath
)
& $Compiler @compileArgs
if ($LASTEXITCODE -ne 0) {
    throw "Failed to compile TABI surface mesh smoke fixture."
}

& $ExePath
if ($LASTEXITCODE -ne 0) {
    throw "TABI surface mesh smoke fixture failed."
}
