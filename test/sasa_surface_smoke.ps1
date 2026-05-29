param(
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\sasa_surface_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot
$Compiler = 'C:\msys64\ucrt64\bin\clang++.exe'
$SourcePath = Join-Path $WorkRoot 'sasa_surface_smoke.cpp'
$ExePath = Join-Path $WorkRoot 'sasa_surface_smoke.exe'

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside SASA surface smoke root: $full"
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
#include "titan_sasa.h"
#include "titan_simd_dispatch.h"
#include "titan_surface.h"

#include <bit>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

void require_close(const std::string &label, double a, double b, double tol) {
  if (std::abs(a - b) > tol) {
    throw std::runtime_error(label + " mismatch: " + std::to_string(a) +
                             " vs " + std::to_string(b));
  }
}

void require_true(const std::string &label, bool value) {
  if (!value)
    throw std::runtime_error(label);
}

TitanSurface::Options options_for(const std::string &mode, bool panels = true,
                                  bool neighbors = true) {
  TitanSurface::Options options;
  options.probe_radius = TitanSASA::kDefaultProbeRadiusA;
  options.point_count = TitanSASA::kDefaultPointCount;
  options.surface_mode = "sas";
  options.acceleration = mode;
  options.emit_panels = panels;
  options.emit_neighbor_workspace = neighbors;
  return options;
}

TitanSurface::Workspace
build_ws(const std::vector<TitanSurface::Atom> &atoms,
         const std::string &mode, TitanExecutor *executor = nullptr,
         bool panels = true, bool neighbors = true) {
  return TitanSurface::build_workspace(atoms, options_for(mode, panels, neighbors),
                                       {}, executor);
}

std::vector<TitanSASA::Atom>
to_sasa_atoms(const std::vector<TitanSurface::Atom> &atoms) {
  std::vector<TitanSASA::Atom> out;
  out.reserve(atoms.size());
  for (const TitanSurface::Atom &atom : atoms)
    out.push_back({atom.x, atom.y, atom.z, atom.radius});
  return out;
}

void compare_workspaces(const std::string &label,
                        const TitanSurface::Workspace &legacy,
                        const TitanSurface::Workspace &optimized,
                        bool compare_panels = true) {
  require_true(label + " area size",
               legacy.area_by_atom.size() == optimized.area_by_atom.size());
  for (std::size_t i = 0; i < legacy.area_by_atom.size(); ++i)
    require_close(label + " area atom " + std::to_string(i),
                  legacy.area_by_atom[i], optimized.area_by_atom[i], 1.0e-9);
  require_close(label + " total area", legacy.total_area, optimized.total_area,
                1.0e-9);
  require_true(label + " exposed points",
               legacy.diagnostics.exposed_points ==
                   optimized.diagnostics.exposed_points);
  require_true(label + " buried points",
               legacy.diagnostics.buried_points ==
                   optimized.diagnostics.buried_points);
  require_true(label + " neighbor offsets",
               legacy.neighbor_offsets == optimized.neighbor_offsets);
  require_true(label + " neighbors", legacy.neighbors == optimized.neighbors);

  if (!compare_panels)
    return;
  require_true(label + " panel count",
               legacy.panels.size() == optimized.panels.size());
  for (std::size_t i = 0; i < legacy.panels.size(); ++i) {
    const TitanSurface::Panel &a = legacy.panels[i];
    const TitanSurface::Panel &b = optimized.panels[i];
    require_true(label + " panel atom " + std::to_string(i),
                 a.atom_index == b.atom_index);
    require_close(label + " panel x " + std::to_string(i), a.x, b.x, 1.0e-12);
    require_close(label + " panel y " + std::to_string(i), a.y, b.y, 1.0e-12);
    require_close(label + " panel z " + std::to_string(i), a.z, b.z, 1.0e-12);
    require_close(label + " panel area " + std::to_string(i), a.area, b.area,
                  1.0e-12);
  }
}

void check_case(const std::string &label,
                const std::vector<TitanSurface::Atom> &surface_atoms) {
  TitanSurface::Workspace legacy = build_ws(surface_atoms, "legacy");
  TitanSurface::Workspace optimized =
      build_ws(surface_atoms, "reordered_padded");
  compare_workspaces(label, legacy, optimized);
  require_true(label + " legacy mode",
               legacy.diagnostics.acceleration_effective == "legacy");
  require_true(label + " optimized mode",
               optimized.diagnostics.acceleration_effective ==
                   "reordered_padded");
  require_true(label + " padded slots multiple of four",
               optimized.diagnostics.padded_neighbor_slots % 4 == 0);

  const std::vector<TitanSASA::Atom> sasa_atoms = to_sasa_atoms(surface_atoms);
  TitanSASA::Options sasa_options;
  const std::vector<double> wrapper =
      TitanSASA::compute_per_atom_sasa(sasa_atoms, sasa_options);
  for (std::size_t i = 0; i < wrapper.size(); ++i)
    require_close(label + " wrapper atom " + std::to_string(i), wrapper[i],
                  legacy.area_by_atom[i], 1.0e-9);
}

std::vector<TitanSurface::Atom> deterministic_atoms(std::size_t count) {
  std::vector<TitanSurface::Atom> atoms;
  atoms.reserve(count);
  std::uint64_t state = 0x9e3779b97f4a7c15ull;
  auto next_unit = [&]() {
    state = state * 6364136223846793005ull + 1442695040888963407ull;
    return static_cast<double>((state >> 11) & ((1ull << 53) - 1ull)) /
           static_cast<double>(1ull << 53);
  };
  for (std::size_t i = 0; i < count; ++i) {
    atoms.push_back({30.0 * next_unit() - 15.0,
                     30.0 * next_unit() - 15.0,
                     30.0 * next_unit() - 15.0,
                     1.2 + 0.7 * next_unit()});
  }
  return atoms;
}

std::uint64_t hash_workspace(const TitanSurface::Workspace &workspace) {
  std::uint64_t hash = 1469598103934665603ull;
  auto mix = [&](std::uint64_t value) {
    hash ^= value;
    hash *= 1099511628211ull;
  };
  for (double value : workspace.area_by_atom)
    mix(std::bit_cast<std::uint64_t>(value));
  mix(workspace.diagnostics.exposed_points);
  mix(workspace.diagnostics.buried_points);
  for (const TitanSurface::Panel &panel : workspace.panels) {
    mix(panel.atom_index);
    mix(std::bit_cast<std::uint64_t>(panel.area));
  }
  return hash;
}

void check_determinism() {
  const std::vector<TitanSurface::Atom> atoms = deterministic_atoms(160);
  TitanExecutor exec1(1);
  TitanExecutor exec4(4);
  TitanSurface::Workspace legacy = build_ws(atoms, "legacy", &exec1, true, true);
  TitanSurface::Workspace optimized1 =
      build_ws(atoms, "reordered_padded", &exec1, true, true);
  TitanSurface::Workspace optimized4 =
      build_ws(atoms, "reordered_padded", &exec4, true, true);
  TitanSurface::Workspace optimized4_repeat =
      build_ws(atoms, "reordered_padded", &exec4, true, true);
  compare_workspaces("deterministic nt1", legacy, optimized1);
  compare_workspaces("deterministic nt4", legacy, optimized4);
  require_true("deterministic repeat hash",
               hash_workspace(optimized4) == hash_workspace(optimized4_repeat));

  TitanSurface::Workspace auto_ws = build_ws(atoms, "auto", &exec4, false, false);
  require_true("auto threshold optimized",
               auto_ws.diagnostics.acceleration_effective ==
                   "reordered_padded");
}

void check_simd_variable_cutoff() {
  alignas(32) double xs[4] = {0.1, 2.0, 0.0, 8.0};
  alignas(32) double ys[4] = {0.0, 0.0, 3.0, 8.0};
  alignas(32) double zs[4] = {0.0, 0.0, 0.0, 8.0};
  alignas(32) double cutoff_sq[4] = {0.25, 1.0, 16.0, -1.0};
  std::int64_t indices[4] = {0, 1, 2, 0};
  const unsigned scalar =
      titan_simd_dispatch::point_vs_gathered4_variable_cutoff_mask_scalar(
          0.0, 0.0, 0.0, xs, ys, zs, indices, cutoff_sq, -1.0, false);
  const unsigned dispatched =
      titan_simd_dispatch::point_vs_gathered4_variable_cutoff_mask(
          0.0, 0.0, 0.0, xs, ys, zs, indices, cutoff_sq, -1.0, false);
  require_true("variable cutoff SIMD agreement", scalar == dispatched);
  require_true("variable cutoff expected mask", scalar == 0x5u);
}

} // namespace

int main() {
  check_simd_variable_cutoff();
  check_case("single_sphere", {{0.0, 0.0, 0.0, 1.5}});
  check_case("two_overlapping_atoms",
             {{0.0, 0.0, 0.0, 1.7}, {2.0, 0.0, 0.0, 1.7}});
  check_case("mini_protein",
             {{-3.2, 0.0, 0.0, 1.55},
              {-1.8, 0.1, 0.0, 1.70},
              {-1.16, 0.5, -1.32, 1.70},
              {-1.72, 1.26, -2.10, 1.52},
              {0.34, 0.2, -1.40, 1.60},
              {1.12, -0.2, -0.28, 1.70}});
  check_case("receptor_ligand_complex",
             {{-4.0, -0.4, 0.0, 1.70},
              {-2.7, 0.2, 0.2, 1.55},
              {-1.4, -0.1, -0.3, 1.80},
              {-0.2, 0.3, 0.1, 1.52},
              {1.2, -0.2, 0.0, 1.70},
              {2.2, 0.4, 0.2, 1.55},
              {3.4, -0.1, -0.1, 1.80},
              {4.4, 0.2, 0.3, 1.52}});

  TitanSurface::Options vdw_options;
  vdw_options.surface_mode = "vdw";
  vdw_options.probe_radius = 9.9;
  vdw_options.acceleration = "reordered_padded";
  vdw_options.emit_panels = true;
  std::vector<TitanSurface::Atom> vdw_atoms{{0.0, 0.0, 0.0, 2.0}};
  TitanSurface::Workspace vdw =
      TitanSurface::build_workspace(vdw_atoms, vdw_options);
  require_close("vdw probe ignored", vdw.total_area, 4.0 * kPi * 4.0,
                1.0e-9);

  check_determinism();
  std::cout << "SASA surface smoke PASS acceleration="
            << vdw.diagnostics.acceleration_effective
            << " panels=" << vdw.panels.size()
            << " area=" << vdw.total_area << "\n";
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
    (Join-Path $RepoRoot 'src\titan_surface.cpp'),
    (Join-Path $RepoRoot 'src\titan_sasa.cpp'),
    (Join-Path $RepoRoot 'src\titan_executor.cpp'),
    $SourcePath,
    '-o',
    $ExePath
)
& $Compiler @compileArgs
if ($LASTEXITCODE -ne 0) {
    throw "Failed to compile SASA surface smoke fixture."
}

& $ExePath
if ($LASTEXITCODE -ne 0) {
    throw "SASA surface smoke fixture failed."
}
