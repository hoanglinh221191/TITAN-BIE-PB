param(
    [string]$Compiler = 'g++',
    [string]$WorkRoot = (Join-Path $PSScriptRoot 'tmp\mol_topology_smoke')
)

$ErrorActionPreference = 'Stop'
Set-StrictMode -Version Latest

$repo = Split-Path -Parent $PSScriptRoot
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

$src = Join-Path $WorkRoot 'mol_topology_smoke.cpp'
$exe = Join-Path $WorkRoot 'mol_topology_smoke.exe'
$top = Join-Path $WorkRoot 'mini.top'
$itp = Join-Path $WorkRoot 'mini.itp'
$standaloneTop = Join-Path $WorkRoot 'standalone_ligand.top'

@'
[ moleculetype ]
MOLA 3

[ atoms ]
1 CT 1 RES A1 1 0.0 12.0
2 CT 1 RES A2 1 0.0 12.0
3 CT 1 RES A3 1 0.0 12.0

[ bonds ]
1 2 1
2 3 1

[ moleculetype ]
MOLB 3

[ atoms ]
1 CT 1 RES B1 1 0.0 12.0
2 CT 1 RES B2 1 0.0 12.0

[ bonds ]
1 2 1
'@ | Set-Content -LiteralPath $itp -Encoding ASCII

@"
#include `"mini.itp`"

[ system ]
mini

[ molecules ]
MOLA 2
MOLB 1
"@ | Set-Content -LiteralPath $top -Encoding ASCII

@'
[ moleculetype ]
LIG 3

[ atoms ]
1 CT 1 LIG L1 1 0.0 12.0
2 CT 1 LIG L2 1 0.0 12.0

[ bonds ]
1 2 1

[ molecules ]
LIG 1
'@ | Set-Content -LiteralPath $standaloneTop -Encoding ASCII

@'
#include "mol_topology.h"
#include "MOL_info.h"

#include <iostream>
#include <vector>

namespace {

void expect(bool condition, const char *message) {
  if (!condition) {
    std::cerr << "FAIL mol_topology_smoke: " << message << "\n";
    std::exit(1);
  }
}

void expect_eq_size(size_t actual, size_t expected, const char *message) {
  if (actual != expected) {
    std::cerr << "FAIL mol_topology_smoke: " << message << ": got "
              << actual << " expected " << expected << "\n";
    std::exit(1);
  }
}

std::vector<mol_info> make_selection(std::initializer_list<size_t> source_ids) {
  std::vector<mol_info> info(source_ids.size());
  size_t i = 0;
  for (const size_t source_id : source_ids) {
    info[i].set_original_atom_id(source_id);
    ++i;
  }
  return info;
}

} // namespace

int main(int argc, char **argv) {
  expect(argc == 3, "expected topology path arguments");
  const std::string top_path = argv[1];

  auto full_info = make_selection({0, 1, 2, 3, 4, 5, 6, 7});
  auto full_result = build_gromacs_mol_topology(top_path, full_info);
  expect(full_result.has_value(), "full topology parse/project failed");
  const MolTopology &full = *full_result;
  expect(full.available(), "full topology not marked available");
  expect_eq_size(full.atom_count(), 8, "full atom_count mismatch");
  expect_eq_size(full.bonds().size(), 5, "full expanded bond count mismatch");
  expect_eq_size(full.neighbors(1).size(), 2,
                 "CSR neighbor count for atom 1 mismatch");
  expect_eq_size(full.neighbors(3).size(), 1,
                 "CSR neighbor count for atom 3 mismatch");
  expect_eq_size(full.molecule_ranges().size(), 3,
                 "molecule range expansion mismatch");

  auto subset_info = make_selection({3, 4});
  auto subset_result = build_gromacs_mol_topology(top_path, subset_info);
  expect(subset_result.has_value(), "subset topology project failed");
  expect_eq_size(subset_result->atom_count(), 2, "subset atom_count mismatch");
  expect_eq_size(subset_result->bonds().size(), 1,
                 "subset projected bond count mismatch");
  expect_eq_size(subset_result->neighbors(0).size(), 1,
                 "subset CSR neighbor mismatch");

  std::vector<mol_info> sequential_info(8);
  auto sequential_result = build_gromacs_mol_topology(top_path, sequential_info);
  expect(sequential_result.has_value(), "sequential fallback failed");
  expect_eq_size(sequential_result->bonds().size(), 5,
                 "sequential fallback bond count mismatch");

  auto standalone_ligand_info = make_selection({6, 7});
  auto standalone_result =
      build_gromacs_mol_topology(argv[2], standalone_ligand_info);
  expect(standalone_result.has_value(),
         "standalone topology local fallback failed");
  expect_eq_size(standalone_result->bonds().size(), 1,
                 "standalone topology local fallback bond count mismatch");

  auto missing_result = build_gromacs_mol_topology("missing_topology.top",
                                                  full_info);
  expect(!missing_result.has_value(), "missing topology unexpectedly parsed");

  std::cout << "PASS mol_topology_smoke\n";
  return 0;
}
'@ | Set-Content -LiteralPath $src -Encoding ASCII

& $Compiler -std=c++23 -O2 -Isrc $src src\mol_topology.cpp -o $exe
if ($LASTEXITCODE -ne 0) {
    throw 'Failed to compile mol_topology smoke'
}

& $exe $top $standaloneTop
if ($LASTEXITCODE -ne 0) {
    throw "mol_topology smoke failed with exit code $LASTEXITCODE"
}
