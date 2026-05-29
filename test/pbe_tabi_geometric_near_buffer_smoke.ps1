param(
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_geometric_near_buffer_smoke'),
    [string]$Cxx = ''
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$RepoRoot = Split-Path -Parent $PSScriptRoot

function Assert-UnderPath([string]$Path, [string]$Root) {
    $full = [System.IO.Path]::GetFullPath($Path)
    $rootBase = [System.IO.Path]::GetFullPath($Root).TrimEnd('\')
    $rootFull = $rootBase + '\'
    if (-not ($full.Equals($rootBase, [System.StringComparison]::OrdinalIgnoreCase) -or
              $full.StartsWith($rootFull, [System.StringComparison]::OrdinalIgnoreCase))) {
        throw "Refusing to clean path outside geometric-near smoke root: $full"
    }
}

function New-CleanDir([string]$Path) {
    Assert-UnderPath -Path $Path -Root $WorkRoot
    if (Test-Path -LiteralPath $Path) { Remove-Item -LiteralPath $Path -Recurse -Force }
    New-Item -ItemType Directory -Path $Path -Force | Out-Null
}

function Resolve-Cxx([string]$Requested) {
    if (-not [string]::IsNullOrWhiteSpace($Requested)) {
        if (-not (Test-Path -LiteralPath $Requested)) { throw "Requested C++ compiler not found: $Requested" }
        return $Requested
    }
    foreach ($candidate in @(
            'C:\msys64\ucrt64\bin\clang++.exe',
            'C:\msys64\mingw64\bin\clang++.exe',
            'C:\msys64\ucrt64\bin\g++.exe',
            'C:\msys64\mingw64\bin\g++.exe')) {
        if (Test-Path -LiteralPath $candidate) { return $candidate }
    }
    $cmd = Get-Command clang++ -ErrorAction SilentlyContinue
    if ($null -ne $cmd) { return $cmd.Source }
    $cmd = Get-Command g++ -ErrorAction SilentlyContinue
    if ($null -ne $cmd) { return $cmd.Source }
    throw 'No C++ compiler found for geometric-near buffer smoke.'
}

New-CleanDir $WorkRoot
$sourcePath = Join-Path $WorkRoot 'pbe_tabi_geometric_near_buffer_smoke.cpp'
$exePath = Join-Path $WorkRoot 'pbe_tabi_geometric_near_buffer_smoke.exe'

Set-Content -LiteralPath $sourcePath -Encoding ASCII -Value @'
#include <cmath>
#include <cstddef>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "tabi/tabi_near.cpp"

namespace {

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
  const Vec3 n = cross(sub(b, a), sub(c, a));
  const double nrm = norm(n);
  if (nrm <= 0.0)
    throw std::runtime_error("degenerate test panel");
  p.area = 0.5 * nrm;
  p.nx = n.x / nrm;
  p.ny = n.y / nrm;
  p.nz = n.z / nrm;
  p.has_triangle_geometry = true;
  return p;
}

void require(bool condition, const std::string &message) {
  if (!condition)
    throw std::runtime_error(message);
}

bool finite_acc(const TitanPBE::Tabi::JcpNearAccum &acc) {
  return std::isfinite(acc.phi) && std::isfinite(acc.phin);
}

TitanPBE::Tabi::JcpNearAccum run_pair(
    const std::vector<TitanPBE::Tabi::Panel> &panels,
    const std::vector<double> &unknowns, const std::string &quad_buffer,
    TitanPBE::Tabi::TreecodeStats &stats) {
  const std::vector<std::size_t> order{1};
  return TitanPBE::Tabi::accumulate_jcp_near_ordered(
      panels[0], 0, panels, unknowns, 80.0, 0.125, true, order, 0,
      order.size(), "duffy_adaptive", quad_buffer, &stats);
}

} // namespace

int main() {
  try {
    const auto target = make_panel({0.0, 0.0, 0.0}, {1.0, 0.0, 0.0},
                                   {0.0, 1.0, 0.0}, 0, 1, 2);
    const auto close_source =
        make_panel({0.0, 0.0, 0.30}, {1.0, 0.0, 0.30},
                   {0.0, 1.0, 0.30}, 10, 11, 12);
    const auto far_source =
        make_panel({0.0, 0.0, 4.0}, {1.0, 0.0, 4.0}, {0.0, 1.0, 4.0},
                   20, 21, 22);

    const auto off_decision =
        TitanPBE::Tabi::classify_geometric_near_buffer_pair(target,
                                                            close_source,
                                                            "off");
    require(!off_decision.active, "quad_buffer=off should not activate");
    const auto gauss_decision =
        TitanPBE::Tabi::classify_geometric_near_buffer_pair(target,
                                                            close_source,
                                                            "gauss7");
    require(gauss_decision.active, "close non-one-ring pair not detected");
    require(!gauss_decision.adaptive, "gauss7 should not request adaptive");
    const auto adaptive_decision =
        TitanPBE::Tabi::classify_geometric_near_buffer_pair(target,
                                                            close_source,
                                                            "adaptive");
    require(adaptive_decision.active && adaptive_decision.adaptive,
            "adaptive buffer did not request subdivided quadrature");
    const auto far_decision =
        TitanPBE::Tabi::classify_geometric_near_buffer_pair(target,
                                                            far_source,
                                                            "adaptive");
    require(!far_decision.active, "far pair should not enter buffer");

    std::vector<TitanPBE::Tabi::Panel> panels{target, close_source};
    std::vector<double> unknowns{0.25, -0.35, 0.50, -0.15};

    TitanPBE::Tabi::TreecodeStats centroid_stats;
    const auto centroid = run_pair(panels, unknowns, "off", centroid_stats);
    require(finite_acc(centroid), "centroid accumulation is non-finite");
    require(centroid_stats.near_centroid_pairs == 1u,
            "off path should use centroid for close non-one-ring pair");

    TitanPBE::Tabi::TreecodeStats gauss_stats;
    const auto gauss = run_pair(panels, unknowns, "gauss7", gauss_stats);
    require(finite_acc(gauss), "gauss7 accumulation is non-finite");
    require(gauss_stats.near_geometric_buffer_pairs == 1u,
            "gauss7 did not record geometric buffer pair");
    require(gauss_stats.near_geometric_buffer_gauss_pairs == 1u,
            "gauss7 counter mismatch");
    require(gauss_stats.near_geometric_buffer_min_gap > 0.0,
            "geometric buffer gap was not recorded");
    require(std::fabs(gauss.phi - centroid.phi) > 1.0e-10 ||
                std::fabs(gauss.phin - centroid.phin) > 1.0e-10,
            "gauss7 quadrature unexpectedly matched centroid exactly");

    TitanPBE::Tabi::TreecodeStats adaptive_stats;
    const auto adaptive = run_pair(panels, unknowns, "adaptive",
                                   adaptive_stats);
    require(finite_acc(adaptive), "adaptive accumulation is non-finite");
    require(adaptive_stats.near_geometric_buffer_pairs == 1u,
            "adaptive did not record geometric buffer pair");
    require(adaptive_stats.near_geometric_buffer_adaptive_pairs == 1u,
            "adaptive counter mismatch");
    require(adaptive_stats.geometric_near_audit[2].pairs == 1u,
            "adaptive audit bucket did not record subdivided close pair");

    panels[1] = far_source;
    TitanPBE::Tabi::TreecodeStats far_stats;
    const auto far = run_pair(panels, unknowns, "adaptive", far_stats);
    require(finite_acc(far), "far accumulation is non-finite");
    require(far_stats.near_geometric_buffer_pairs == 0u,
            "far pair should not use geometric buffer");
    require(far_stats.near_centroid_pairs == 1u,
            "far pair should fall back to centroid");

    std::cout << "PBE TABI geometric-near buffer smoke PASS"
              << " gap_ratio=" << gauss_decision.gap_ratio
              << " center_ratio=" << gauss_decision.center_ratio << "\n";
    return 0;
  } catch (const std::exception &ex) {
    std::cerr << "PBE TABI geometric-near buffer smoke FAIL: " << ex.what()
              << "\n";
    return 1;
  }
}
'@

$compiler = Resolve-Cxx $Cxx
$includeRoot = $RepoRoot -replace '\\', '/'
& $compiler -std=c++23 -O2 "-I$includeRoot" $sourcePath -o $exePath
if ($LASTEXITCODE -ne 0) {
    throw "C++ compile failed for geometric-near buffer smoke using $compiler"
}
& $exePath
if ($LASTEXITCODE -ne 0) {
    throw "geometric-near buffer smoke executable failed"
}
