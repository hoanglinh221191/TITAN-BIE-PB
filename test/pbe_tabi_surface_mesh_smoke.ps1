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

TitanPBE::Tabi::Options options(const std::string &mesh_mode = "latlong",
                                int surface_points = 96,
                                bool highres_surface = false) {
  TitanPBE::Tabi::Options opt;
  opt.surface_mode = "vdw";
  opt.surface_geometry = "triangle";
  opt.diagonal_mode = "solid-angle";
  opt.probe_radius = 0.0;
  opt.sphere_mesh_mode = mesh_mode;
  opt.surface_point_count = surface_points;
  opt.highres_surface = highres_surface;
  return opt;
}

TitanPBE::Tabi::Result build_case(
    const std::string &label,
    const std::vector<TitanPBE::Tabi::Atom> &atoms,
    const std::string &mesh_mode = "latlong", int surface_points = 96) {
  TitanPBE::Tabi::Result result;
  std::vector<TitanPBE::Tabi::Panel> panels =
      TitanPBE::Tabi::build_surface(atoms, options(mesh_mode, surface_points),
                                    &result);
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

double panel_face_radial_dot(const TitanPBE::Tabi::Panel &panel,
                             const TitanPBE::Tabi::Atom &atom) {
  const double ux = panel.v1x - panel.v0x;
  const double uy = panel.v1y - panel.v0y;
  const double uz = panel.v1z - panel.v0z;
  const double vx = panel.v2x - panel.v0x;
  const double vy = panel.v2y - panel.v0y;
  const double vz = panel.v2z - panel.v0z;
  double nx = uy * vz - uz * vy;
  double ny = uz * vx - ux * vz;
  double nz = ux * vy - uy * vx;
  const double nlen = std::sqrt(nx * nx + ny * ny + nz * nz);
  const double cx = (panel.v0x + panel.v1x + panel.v2x) / 3.0 - atom.x;
  const double cy = (panel.v0y + panel.v1y + panel.v2y) / 3.0 - atom.y;
  const double cz = (panel.v0z + panel.v1z + panel.v2z) / 3.0 - atom.z;
  const double clen = std::sqrt(cx * cx + cy * cy + cz * cz);
  double dot = (nx * cx + ny * cy + nz * cz) / (nlen * clen);
  if (dot < 0.0)
    dot = -dot;
  return dot;
}

double panel_flat_centroid_error(const TitanPBE::Tabi::Panel &panel,
                                 const TitanPBE::Tabi::Atom &atom) {
  const double cx = (panel.v0x + panel.v1x + panel.v2x) / 3.0 - atom.x;
  const double cy = (panel.v0y + panel.v1y + panel.v2y) / 3.0 - atom.y;
  const double cz = (panel.v0z + panel.v1z + panel.v2z) / 3.0 - atom.z;
  return std::abs(std::sqrt(cx * cx + cy * cy + cz * cz) - atom.radius);
}

void require_geometry_telemetry_matches_panels(
    const std::string &label, const TitanPBE::Tabi::Atom &atom,
    const std::string &mesh_mode, int surface_points) {
  TitanPBE::Tabi::Result result;
  const std::vector<TitanPBE::Tabi::Panel> panels =
      TitanPBE::Tabi::build_surface({atom}, options(mesh_mode, surface_points),
                                    &result);
  double dot_min = 1.0;
  double centroid_error_max = 0.0;
  for (const auto &panel : panels) {
    dot_min = std::min(dot_min, panel_face_radial_dot(panel, atom));
    centroid_error_max = std::max(
        centroid_error_max, panel_flat_centroid_error(panel, atom));
  }
  require_close_rel(label + " face normal dot",
                    result.surface_quality_normal_radial_dot_min, dot_min,
                    1.0e-12);
  require_close_rel(label + " flat centroid error",
                    result.surface_quality_centroid_radial_error_max,
                    centroid_error_max, 1.0e-12);
  require_true(label + " planar area finite",
               std::isfinite(result.surface_quality_planar_area_a2) &&
                   result.surface_quality_planar_area_a2 > 0.0);
  require_true(label + " enclosed volume finite",
               std::isfinite(result.surface_quality_enclosed_volume_a3) &&
                   result.surface_quality_enclosed_volume_a3 > 0.0);
  require_true(label + " triangle area cv finite",
               std::isfinite(result.surface_quality_triangle_area_cv) &&
                   result.surface_quality_triangle_area_cv >= 0.0);
}

void require_highres_icosphere_counts(const TitanPBE::Tabi::Atom &atom) {
  struct Case {
    int requested;
    std::size_t expected_panels;
  };
  const std::vector<Case> cases{{1281, 5120}, {5121, 20480},
                                {20481, 81920}};
  double previous_area_error = std::numeric_limits<double>::infinity();
  double previous_volume_error = std::numeric_limits<double>::infinity();
  for (const auto &test_case : cases) {
    TitanPBE::Tabi::Result result;
    const std::vector<TitanPBE::Tabi::Panel> panels =
        TitanPBE::Tabi::build_surface(
            {atom}, options("icosphere", test_case.requested, true), &result);
    require_true("highres panel count",
                 panels.size() == test_case.expected_panels);
    require_true("highres result triangles",
                 result.mesh_triangles == test_case.expected_panels);
    require_true("highres planar area finite",
                 std::isfinite(result.surface_quality_planar_area_a2));
    require_true("highres volume finite",
                 std::isfinite(result.surface_quality_enclosed_volume_a3));
    require_true("highres planar area convergence",
                 result.surface_quality_planar_area_rel_error_vs_4piR2 <
                     previous_area_error);
    require_true("highres volume convergence",
                 result.surface_quality_enclosed_volume_rel_error_vs_4piR3_over3 <
                     previous_volume_error);
    previous_area_error =
        result.surface_quality_planar_area_rel_error_vs_4piR2;
    previous_volume_error =
        result.surface_quality_enclosed_volume_rel_error_vs_4piR3_over3;
  }
}

} // namespace

int main() {
  const double radius = 2.0;
  const auto single =
      build_case("single_sphere", {{0.0, 0.0, 0.0, 1.0, radius}});
  require_close_rel("single sphere area", single.surface_area_a2,
                    4.0 * kPi * radius * radius, 0.08);
  require_true("single sphere nonmanifold", single.nonmanifold_edges == 0);

  const TitanPBE::Tabi::Atom sphere{0.0, 0.0, 0.0, 1.0, radius};
  require_geometry_telemetry_matches_panels("icosphere telemetry", sphere,
                                            "icosphere", 320);
  require_geometry_telemetry_matches_panels("latlong telemetry", sphere,
                                            "latlong", 81);
  const auto ico_coarse = build_case("ico_coarse", {sphere}, "icosphere", 80);
  const auto ico_fine = build_case("ico_fine", {sphere}, "icosphere", 1280);
  require_true("icosphere planar area convergence",
               ico_fine.surface_quality_planar_area_rel_error_vs_4piR2 <
                   ico_coarse.surface_quality_planar_area_rel_error_vs_4piR2);
  require_true("icosphere volume convergence",
               ico_fine.surface_quality_enclosed_volume_rel_error_vs_4piR3_over3 <
                   ico_coarse.surface_quality_enclosed_volume_rel_error_vs_4piR3_over3);
  require_highres_icosphere_counts(sphere);

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
