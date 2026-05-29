param(
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\pbe_tabi_near_kernel_blocks_smoke'),
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
        throw "Refusing to clean path outside near-kernel block smoke root: $full"
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
    throw 'No C++ compiler found for near-kernel block smoke.'
}

New-CleanDir $WorkRoot
$sourcePath = Join-Path $WorkRoot 'pbe_tabi_near_kernel_blocks_smoke.cpp'
$exePath = Join-Path $WorkRoot 'pbe_tabi_near_kernel_blocks_smoke.exe'

Set-Content -LiteralPath $sourcePath -Encoding ASCII -Value @'
#include <algorithm>
#include <array>
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

TitanPBE::Tabi::TreecodeStats run_pair(
    const std::string &name, const std::vector<TitanPBE::Tabi::Panel> &panels,
    const std::vector<double> &unknowns, std::size_t source_index,
    bool expect_self, bool expect_one_ring, bool expect_duffy,
    bool expect_centroid) {
  TitanPBE::Tabi::TreecodeStats stats;
  const std::vector<std::size_t> order{source_index};
  const auto acc = TitanPBE::Tabi::accumulate_jcp_near_ordered(
      panels[0], 0, panels, unknowns, 80.0, 0.125, true, order, 0,
      order.size(), "duffy_adaptive", "off", &stats);
  require(finite_acc(acc), name + ": non-finite accumulation");
  require(stats.near_fallbacks == 0, name + ": unexpected near fallback");
  require(stats.near_self_pairs == (expect_self ? 1u : 0u),
          name + ": self pair counter mismatch");
  require(stats.near_one_ring_pairs == (expect_one_ring ? 1u : 0u),
          name + ": one-ring pair counter mismatch");
  require(stats.near_duffy_pairs == (expect_duffy ? 1u : 0u),
          name + ": Duffy pair counter mismatch");
  require(stats.near_centroid_pairs == (expect_centroid ? 1u : 0u),
          name + ": centroid pair counter mismatch");
  if (expect_one_ring) {
    require(stats.near_one_ring_k1_accepted + stats.near_one_ring_k1_rejected ==
                1u,
            name + ": K1 one-ring decision missing");
    require(stats.near_one_ring_k2_accepted + stats.near_one_ring_k2_rejected ==
                1u,
            name + ": K2 one-ring decision missing");
    require(stats.near_one_ring_k3_accepted + stats.near_one_ring_k3_rejected ==
                1u,
            name + ": K3 one-ring decision missing");
    require(stats.near_one_ring_k4_accepted + stats.near_one_ring_k4_rejected ==
                1u,
            name + ": K4 one-ring decision missing");
  } else {
    require(stats.near_one_ring_k1_accepted + stats.near_one_ring_k1_rejected ==
                0u,
            name + ": unexpected one-ring K1 decision");
  }
  return stats;
}

} // namespace

int main() {
  try {
    std::vector<TitanPBE::Tabi::Panel> panels;
    panels.push_back(
        make_panel({0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, 0, 1, 2));
    panels.push_back(make_panel({1.0, 0.0, 0.0}, {0.0, 0.0, 0.0},
                                {1.0, -1.0, 0.0}, 1, 0, 3));
    panels.push_back(make_panel({0.0, 0.0, 0.0}, {-1.0, 0.0, 0.0},
                                {0.0, -1.0, 0.0}, 0, 4, 5));
    panels.push_back(make_panel({0.15, 0.15, 0.05}, {0.65, 0.15, 0.05},
                                {0.15, 0.65, 0.05}, 10, 11, 12));
    panels.push_back(make_panel({1.0, 0.0, 0.0}, {1.0, -1.0, 0.0},
                                {0.0, 0.0, 0.0}, 1, 3, 0));

    const std::size_t n = panels.size();
    std::vector<double> unknowns(2 * n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
      unknowns[i] = 0.25 + 0.1 * static_cast<double>(i);
      unknowns[n + i] = -0.4 + 0.05 * static_cast<double>(i);
    }

    const auto self = run_pair("self", panels, unknowns, 0, true, false, true,
                               false);
    require(self.near_k1_accepted == 1u && self.near_k2_accepted == 1u &&
                self.near_k3_accepted == 1u && self.near_k4_accepted == 1u,
            "self ASI block decisions were not accepted");

    const auto edge = run_pair("edge-neighbor", panels, unknowns, 1, false,
                               true, true, false);
    const auto vertex = run_pair("vertex-neighbor", panels, unknowns, 2, false,
                                 true, true, false);
    const auto smooth = run_pair("smooth-near", panels, unknowns, 3, false,
                                 false, false, true);
    (void)smooth;
    const auto flipped = run_pair("edge-neighbor-flipped", panels, unknowns, 4,
                                  false, true, true, false);

    const std::size_t k1_accept = edge.near_one_ring_k1_accepted +
                                  vertex.near_one_ring_k1_accepted +
                                  flipped.near_one_ring_k1_accepted;
    const std::size_t k2_accept = edge.near_one_ring_k2_accepted +
                                  vertex.near_one_ring_k2_accepted +
                                  flipped.near_one_ring_k2_accepted;
    const std::size_t k3_accept = edge.near_one_ring_k3_accepted +
                                  vertex.near_one_ring_k3_accepted +
                                  flipped.near_one_ring_k3_accepted;
    const std::size_t k4_accept = edge.near_one_ring_k4_accepted +
                                  vertex.near_one_ring_k4_accepted +
                                  flipped.near_one_ring_k4_accepted;
    const std::size_t accepted_total =
        k1_accept + k2_accept + k3_accept + k4_accept;
    require(accepted_total > 0u,
            "one-ring CPV guard rejected every canonical candidate");

    std::cout << "PBE TABI near-kernel block smoke PASS"
              << " k1_accept=" << k1_accept << " k2_accept=" << k2_accept
              << " k3_accept=" << k3_accept << " k4_accept=" << k4_accept
              << "\n";
    return 0;
  } catch (const std::exception &ex) {
    std::cerr << "PBE TABI near-kernel block smoke FAIL: " << ex.what()
              << "\n";
    return 1;
  }
}
'@

$compiler = Resolve-Cxx $Cxx
$includeRoot = $RepoRoot -replace '\\', '/'
& $compiler -std=c++23 -O2 "-I$includeRoot" $sourcePath -o $exePath
if ($LASTEXITCODE -ne 0) {
    throw "C++ compile failed for near-kernel block smoke using $compiler"
}
& $exePath
if ($LASTEXITCODE -ne 0) {
    throw "near-kernel block smoke executable failed"
}
