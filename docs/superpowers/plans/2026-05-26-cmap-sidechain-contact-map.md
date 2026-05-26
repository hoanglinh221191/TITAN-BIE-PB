# CMap Side-chain Contact Map Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add Slice A from `TITAN_STANDARD_MD_ANALYSIS_ROADMAP.md`: an opt-in `-cmap` module that exports side-chain residue contact probability matrices and sparse edge summaries.

**Architecture:** Implement a small shared contact primitive first, then a streaming `TitanContactMapAnalyzer` modeled after `TitanSqAnalyzer`. Wire it through `arguments`, `impl_engine.cpp`, `info.h`, `help.cpp`, and docs without changing existing default analysis paths.

**Tech Stack:** C++23, existing `Multi_MOL` / `MOL` / `mol_info` data model, `titan_neighbor::HalfShellGrid`, `TitanExecutor`, PowerShell smoke scripts, `compile.bat` native/portable profiles.

---

## Current Context

- Existing dirty files in this checkout include `arguments.cpp`, `arguments.h`, `engine_postproc.cpp`, `help.cpp`, `impl_engine.cpp`, and `parse_arguments.cpp`. Before implementation, read each target hunk and preserve unrelated edits.
- `mol_info::check_side_chain()` already treats glycine as side-chain-capable and marks atom side-chain state. Do not rely on this alone for contact-map policy because `-cmap-mode` needs explicit heavy-atom/backbone/CB/CA behavior.
- `TitanSqAnalyzer` is the closest module pattern: it configures from args, initializes from `Multi_MOL`, processes frames by streaming pairs through `HalfShellGrid`, then exports per-trajectory files in `finalize_trajectory`.
- `info.h` centralizes output filenames. Add CMap filenames there instead of hardcoding strings in solver code.
- Root `Makefile` picks up new root `*.cpp` automatically. No source-list edit is needed for new root `.cpp` files.

## File Structure

- Create `titan_contact_primitives.h`: residue/atom classification and contact-key structs shared by `-cmap` and future `-qnative`.
- Create `titan_contact_primitives.cpp`: side-chain heavy atom, CB, CA, and backbone selection logic.
- Create `titan_contact_map.h`: `TitanContactMapAnalyzer` public API and accumulator state.
- Create `titan_contact_map.cpp`: streaming contact detection, occupancy accumulation, and exports.
- Create `test/titan_contact_primitives_smoke.cpp`: pure C++ unit smoke for atom classification and residue-key normalization.
- Create `test/titan_contact_map_smoke.cpp`: pure C++ analyzer smoke using an in-memory `Multi_MOL` with known side-chain contacts.
- Create `test/cmap_smoke.ps1`: build/run smoke that checks output headers and expected contact probabilities.
- Modify `arguments.h`: add `CMapArgs` and `arguments::cmap`.
- Modify `parse_arguments.cpp`: parse `-cmap`, `-cmap-mode`, `-cmap-cut`, `-cmap-domain`, `-cmap-min-seqsep`.
- Modify `arguments.cpp`: validate and normalize CMap options.
- Modify `info.h`: add contact-map filename helpers.
- Modify `impl_engine.cpp`: instantiate/configure/init/process/finalize `TitanContactMapAnalyzer`, with timing log.
- Modify `help.cpp`: document CLI and outputs.
- Modify `manual.txt`: add a short public section after RMS/Rg/SQ or near contact/IFP.
- Modify `AGENT_ONBOARDING.md`: add the new module to analysis surface notes.
- Modify `TITAN_STANDARD_MD_ANALYSIS_ROADMAP.md`: mark Slice A implementation status after verification.

---

### Task 1: Contact Primitive Unit Smoke

**Files:**
- Create: `titan_contact_primitives.h`
- Create: `titan_contact_primitives.cpp`
- Create: `test/titan_contact_primitives_smoke.cpp`

- [ ] **Step 1: Write the failing primitive smoke**

Create `test/titan_contact_primitives_smoke.cpp`:

```cpp
#include "MOL_info.h"
#include "titan_contact_primitives.h"

#include <cstdlib>
#include <iostream>

namespace {

mol_info atom(std::size_t resnr, const char *res, const char *name) {
  mol_info info(resnr, "", name, "A", res);
  info.set_chainnum(0);
  return info;
}

void require(bool ok, const char *message) {
  if (!ok) {
    std::cerr << "FAIL: " << message << "\n";
    std::exit(1);
  }
}

} // namespace

int main() {
  using titan_contact::ContactMode;

  require(titan_contact::parse_contact_mode("sidechain-heavy") ==
              ContactMode::SidechainHeavy,
          "sidechain-heavy parses");
  require(titan_contact::parse_contact_mode("cb") == ContactMode::CB,
          "cb parses");
  require(titan_contact::parse_contact_mode("ca") == ContactMode::CA,
          "ca parses");
  require(titan_contact::parse_contact_mode("bb") == ContactMode::Backbone,
          "bb parses");

  require(titan_contact::is_selected_contact_atom(atom(1, "ALA", "CB"),
                                                 ContactMode::SidechainHeavy),
          "ALA CB is side-chain heavy");
  require(!titan_contact::is_selected_contact_atom(atom(1, "ALA", "CA"),
                                                  ContactMode::SidechainHeavy),
          "ALA CA is excluded from side-chain-heavy");
  require(!titan_contact::is_selected_contact_atom(atom(1, "ALA", "HB1"),
                                                  ContactMode::SidechainHeavy),
          "hydrogen is excluded from side-chain-heavy");
  require(titan_contact::is_selected_contact_atom(atom(2, "GLY", "CA"),
                                                 ContactMode::SidechainHeavy),
          "GLY CA is side-chain-heavy fallback");
  require(titan_contact::is_selected_contact_atom(atom(3, "LEU", "CB"),
                                                 ContactMode::CB),
          "CB mode keeps CB");
  require(!titan_contact::is_selected_contact_atom(atom(3, "LEU", "CG"),
                                                  ContactMode::CB),
          "CB mode excludes CG");
  require(titan_contact::is_selected_contact_atom(atom(4, "SER", "CA"),
                                                 ContactMode::CA),
          "CA mode keeps CA");
  require(titan_contact::is_selected_contact_atom(atom(4, "SER", "N"),
                                                 ContactMode::Backbone),
          "BB mode keeps N");
  require(titan_contact::is_selected_contact_atom(atom(4, "SER", "O"),
                                                 ContactMode::Backbone),
          "BB mode keeps O");

  titan_contact::ResidueKey a{0, 10, "A", "ALA", 1};
  titan_contact::ResidueKey b{0, 12, "A", "LEU", 2};
  titan_contact::ResiduePairKey pair = titan_contact::make_ordered_pair(a, b);
  require(pair.first.resnr == 10 && pair.second.resnr == 12,
          "ordered pair preserves ascending residue order");
  pair = titan_contact::make_ordered_pair(b, a);
  require(pair.first.resnr == 10 && pair.second.resnr == 12,
          "ordered pair sorts descending input");

  std::cout << "PASS titan_contact_primitives_smoke\n";
  return 0;
}
```

- [ ] **Step 2: Run the smoke and verify it fails**

Run:

```powershell
C:\msys64\ucrt64\bin\clang++.exe -std=c++23 -I. test\titan_contact_primitives_smoke.cpp titan_contact_primitives.cpp -o test\titan_contact_primitives_smoke.exe
```

Expected: compile fails because `titan_contact_primitives.h/.cpp` do not exist yet.

- [ ] **Step 3: Add primitive header**

Create `titan_contact_primitives.h`:

```cpp
#ifndef TITAN_CONTACT_PRIMITIVES_H
#define TITAN_CONTACT_PRIMITIVES_H

#include "MOL_info.h"

#include <cstddef>
#include <cstdint>
#include <string>

namespace titan_contact {

enum class ContactMode : std::uint8_t {
  SidechainHeavy,
  CB,
  CA,
  Backbone
};

struct ResidueKey {
  std::size_t chainnum{};
  std::size_t resnr{};
  std::string chain{" "};
  std::string resname{"UNK"};
  std::size_t ordinal{};
};

struct ResiduePairKey {
  ResidueKey first;
  ResidueKey second;

  bool operator==(const ResiduePairKey &other) const noexcept {
    return first.chainnum == other.first.chainnum &&
           first.resnr == other.first.resnr &&
           second.chainnum == other.second.chainnum &&
           second.resnr == other.second.resnr;
  }
};

struct ResiduePairKeyHash {
  std::size_t operator()(const ResiduePairKey &key) const noexcept;
};

std::string uppercase_copy(std::string value);
std::string trim_copy(std::string value);
ContactMode parse_contact_mode(const std::string &raw);
const char *contact_mode_name(ContactMode mode) noexcept;
bool is_backbone_atom_name(const std::string &atom_name);
bool is_hydrogen_atom_name(const std::string &atom_name);
bool is_selected_contact_atom(const mol_info &atom, ContactMode mode);
ResidueKey residue_key_from_atom(const mol_info &atom, std::size_t ordinal);
ResiduePairKey make_ordered_pair(ResidueKey lhs, ResidueKey rhs);
bool same_residue(const ResidueKey &lhs, const ResidueKey &rhs) noexcept;
bool same_chain(const ResidueKey &lhs, const ResidueKey &rhs) noexcept;

} // namespace titan_contact

#endif // TITAN_CONTACT_PRIMITIVES_H
```

- [ ] **Step 4: Add primitive implementation**

Create `titan_contact_primitives.cpp`:

```cpp
#include "titan_contact_primitives.h"

#include <algorithm>
#include <cctype>
#include <stdexcept>

namespace titan_contact {

std::string uppercase_copy(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(),
                 [](unsigned char c) {
                   return static_cast<char>(std::toupper(c));
                 });
  return value;
}

std::string trim_copy(std::string value) {
  const auto first = value.find_first_not_of(" \t");
  if (first == std::string::npos)
    return {};
  const auto last = value.find_last_not_of(" \t");
  return value.substr(first, last - first + 1);
}

ContactMode parse_contact_mode(const std::string &raw) {
  std::string mode = uppercase_copy(trim_copy(raw));
  std::replace(mode.begin(), mode.end(), '_', '-');
  if (mode == "SIDECHAIN-HEAVY" || mode == "SIDECHAIN")
    return ContactMode::SidechainHeavy;
  if (mode == "CB")
    return ContactMode::CB;
  if (mode == "CA")
    return ContactMode::CA;
  if (mode == "BB" || mode == "BACKBONE")
    return ContactMode::Backbone;
  throw std::invalid_argument(
      "CMap: -cmap-mode must be one of sidechain-heavy, cb, ca, bb.");
}

const char *contact_mode_name(ContactMode mode) noexcept {
  switch (mode) {
  case ContactMode::SidechainHeavy:
    return "sidechain-heavy";
  case ContactMode::CB:
    return "cb";
  case ContactMode::CA:
    return "ca";
  case ContactMode::Backbone:
    return "bb";
  }
  return "sidechain-heavy";
}

bool is_backbone_atom_name(const std::string &atom_name) {
  const std::string atom = uppercase_copy(trim_copy(atom_name));
  return atom == "N" || atom == "CA" || atom == "C" || atom == "O" ||
         atom == "OXT" || atom == "H" || atom == "H1" || atom == "H2" ||
         atom == "H3" || atom == "HA";
}

bool is_hydrogen_atom_name(const std::string &atom_name) {
  const std::string atom = uppercase_copy(trim_copy(atom_name));
  return !atom.empty() && atom.front() == "H"[0];
}

bool is_selected_contact_atom(const mol_info &atom, ContactMode mode) {
  const std::string atom_name = uppercase_copy(trim_copy(atom.get_atom()));
  const std::string res_name = uppercase_copy(trim_copy(atom.get_res_name()));
  if (atom_name.empty())
    return false;
  switch (mode) {
  case ContactMode::SidechainHeavy:
    if (is_hydrogen_atom_name(atom_name))
      return false;
    if (res_name.rfind("GLY", 0) == 0)
      return atom_name == "CA";
    return !is_backbone_atom_name(atom_name);
  case ContactMode::CB:
    if (res_name.rfind("GLY", 0) == 0)
      return atom_name == "CA";
    return atom_name == "CB";
  case ContactMode::CA:
    return atom_name == "CA";
  case ContactMode::Backbone:
    return is_backbone_atom_name(atom_name) && !is_hydrogen_atom_name(atom_name);
  }
  return false;
}

ResidueKey residue_key_from_atom(const mol_info &atom, std::size_t ordinal) {
  ResidueKey key{};
  key.chainnum = atom.get_chainnum();
  key.resnr = atom.get_resnr();
  key.chain = atom.get_chain();
  key.resname = atom.get_res_name();
  key.ordinal = ordinal;
  return key;
}

bool same_residue(const ResidueKey &lhs, const ResidueKey &rhs) noexcept {
  return lhs.chainnum == rhs.chainnum && lhs.resnr == rhs.resnr;
}

bool same_chain(const ResidueKey &lhs, const ResidueKey &rhs) noexcept {
  return lhs.chainnum == rhs.chainnum;
}

ResiduePairKey make_ordered_pair(ResidueKey lhs, ResidueKey rhs) {
  const bool swap_needed =
      rhs.chainnum < lhs.chainnum ||
      (rhs.chainnum == lhs.chainnum && rhs.resnr < lhs.resnr);
  if (swap_needed)
    std::swap(lhs, rhs);
  return ResiduePairKey{std::move(lhs), std::move(rhs)};
}

std::size_t ResiduePairKeyHash::operator()(
    const ResiduePairKey &key) const noexcept {
  const auto mix = [](std::size_t seed, std::size_t value) {
    seed ^= value + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2);
    return seed;
  };
  std::size_t seed = 0;
  seed = mix(seed, key.first.chainnum);
  seed = mix(seed, key.first.resnr);
  seed = mix(seed, key.second.chainnum);
  seed = mix(seed, key.second.resnr);
  return seed;
}

} // namespace titan_contact
```

- [ ] **Step 5: Run primitive smoke**

Run:

```powershell
C:\msys64\ucrt64\bin\clang++.exe -std=c++23 -I. test\titan_contact_primitives_smoke.cpp titan_contact_primitives.cpp -o test\titan_contact_primitives_smoke.exe
.\test\titan_contact_primitives_smoke.exe
```

Expected: `PASS titan_contact_primitives_smoke`.

- [ ] **Step 6: Commit Task 1**

Run:

```powershell
git add titan_contact_primitives.h titan_contact_primitives.cpp test/titan_contact_primitives_smoke.cpp
git commit -m "Add TITAN contact primitive helpers"
```

If the working tree still contains unrelated dirty files, stage only these three paths.

---

### Task 2: CMap Arguments, Validation, and Output Names

**Files:**
- Modify: `arguments.h`
- Modify: `parse_arguments.cpp`
- Modify: `arguments.cpp`
- Modify: `info.h`
- Modify: `help.cpp`

- [ ] **Step 1: Add `CMapArgs` to `arguments.h`**

Add near `SqArgs` / `HBondKineticsArgs`:

```cpp
struct CMapArgs {
  bool flag_cmap{false};
  std::string mode{"sidechain-heavy"};
  double cutoff_ang{4.5};
  std::string domain{"both"};
  int min_seqsep{3};
};
```

Add the member in `struct arguments` near other module args:

```cpp
CMapArgs cmap;
```

- [ ] **Step 2: Parse CLI and arg-file syntax**

In `parse_arguments.cpp`, add handlers in the same region as `-sq` / `-hbk`:

```cpp
if (arg == "-cmap" || arg == "--cmap") {
  cmap.flag_cmap = true;
  continue;
}
if (arg.rfind("-cmap-mode=", 0) == 0 ||
    arg.rfind("--cmap-mode=", 0) == 0) {
  cmap.flag_cmap = true;
  cmap.mode = arg.substr(arg.find('=') + 1);
  continue;
}
if ((arg == "-cmap-mode" || arg == "--cmap-mode") && i + 1 < argc) {
  cmap.flag_cmap = true;
  cmap.mode = string(argv[++i]);
  continue;
}
if (arg.rfind("-cmap-cut=", 0) == 0 ||
    arg.rfind("--cmap-cut=", 0) == 0) {
  cmap.flag_cmap = true;
  cmap.cutoff_ang = stod(arg.substr(arg.find('=') + 1));
  continue;
}
if ((arg == "-cmap-cut" || arg == "--cmap-cut") && i + 1 < argc) {
  cmap.flag_cmap = true;
  cmap.cutoff_ang = stod(string(argv[++i]));
  continue;
}
if (arg.rfind("-cmap-domain=", 0) == 0 ||
    arg.rfind("--cmap-domain=", 0) == 0) {
  cmap.flag_cmap = true;
  cmap.domain = arg.substr(arg.find('=') + 1);
  continue;
}
if ((arg == "-cmap-domain" || arg == "--cmap-domain") && i + 1 < argc) {
  cmap.flag_cmap = true;
  cmap.domain = string(argv[++i]);
  continue;
}
if (arg.rfind("-cmap-min-seqsep=", 0) == 0 ||
    arg.rfind("--cmap-min-seqsep=", 0) == 0) {
  cmap.flag_cmap = true;
  cmap.min_seqsep = stoi(arg.substr(arg.find('=') + 1));
  continue;
}
if ((arg == "-cmap-min-seqsep" || arg == "--cmap-min-seqsep") &&
    i + 1 < argc) {
  cmap.flag_cmap = true;
  cmap.min_seqsep = stoi(string(argv[++i]));
  continue;
}
```

- [ ] **Step 3: Validate CMap options**

In `arguments.cpp::post_process()`, after nearby module validation, add:

```cpp
if (cmap.flag_cmap) {
  std::transform(cmap.mode.begin(), cmap.mode.end(), cmap.mode.begin(),
                 [](unsigned char c) {
                   return static_cast<char>(std::tolower(c));
                 });
  std::replace(cmap.mode.begin(), cmap.mode.end(), '_', '-');
  if (cmap.mode != "sidechain-heavy" && cmap.mode != "sidechain" &&
      cmap.mode != "cb" && cmap.mode != "ca" && cmap.mode != "bb" &&
      cmap.mode != "backbone") {
    TitanUI::print_error(
        "CMap: -cmap-mode must be one of sidechain-heavy, cb, ca, bb.");
    exit(EXIT_FAILURE);
  }
  if (cmap.mode == "sidechain")
    cmap.mode = "sidechain-heavy";
  if (cmap.mode == "backbone")
    cmap.mode = "bb";
  if (cmap.cutoff_ang <= 0.0) {
    TitanUI::print_error("CMap: -cmap-cut must be > 0.");
    exit(EXIT_FAILURE);
  }
  std::transform(cmap.domain.begin(), cmap.domain.end(), cmap.domain.begin(),
                 [](unsigned char c) {
                   return static_cast<char>(std::tolower(c));
                 });
  if (cmap.domain != "intra" && cmap.domain != "inter" &&
      cmap.domain != "both") {
    TitanUI::print_error(
        "CMap: -cmap-domain must be one of intra, inter, both.");
    exit(EXIT_FAILURE);
  }
  if (cmap.min_seqsep < 0) {
    TitanUI::print_error("CMap: -cmap-min-seqsep must be >= 0.");
    exit(EXIT_FAILURE);
  }
}
```

- [ ] **Step 4: Add filename helpers**

In `info.h`, add near other analysis output helpers:

```cpp
string get_cmap_matrix_filename(const arguments &, const std::vector<string> &,
                                int input_name) const {
  if (input_name > 0)
    return "contact_probability_matrix_" + std::to_string(input_name) + ".csv";
  return "contact_probability_matrix.csv";
}
string get_cmap_edges_filename(const arguments &, const std::vector<string> &,
                               int input_name) const {
  if (input_name > 0)
    return "contact_probability_edges_" + std::to_string(input_name) + ".csv";
  return "contact_probability_edges.csv";
}
string get_cmap_summary_filename(const arguments &, const std::vector<string> &,
                                 int input_name) const {
  if (input_name > 0)
    return "contact_probability_summary_" + std::to_string(input_name) + ".dat";
  return "contact_probability_summary.dat";
}
```

- [ ] **Step 5: Add help text**

In `help.cpp`, add a section near IFP or SQ:

```cpp
  2.X Side-chain Contact Probability Map

    Purpose:
      Compute an untyped residue-residue contact probability matrix. The
      scientific default is side-chain heavy-atom contact, with glycine CA
      fallback.

    Main inputs:
      -cmap
      -cmap-mode <sidechain-heavy|cb|ca|bb>
      -cmap-cut <A>
      -cmap-domain <intra|inter|both>
      -cmap-min-seqsep <N>

    Typical outputs:
      contact_probability_matrix.csv
      contact_probability_edges.csv
      contact_probability_summary.dat

    Notes:
      - `sidechain-heavy` excludes backbone atoms and hydrogens.
      - `cb` and `ca` are fast/legacy anchor modes, not the default scientific
        definition.
      - `intra` applies `-cmap-min-seqsep`; `inter` ignores sequence separation.
```

Also add concise CLI rows in the compact option listing near other analysis flags.

- [ ] **Step 6: Build parser smoke**

Run:

```powershell
$env:TITAN_NO_PAUSE='1'
.\compile.bat native all
.\titan.exe -h | Select-String -Pattern '-cmap','contact_probability_matrix'
```

Expected: build succeeds and help output contains `-cmap` and `contact_probability_matrix.csv`.

- [ ] **Step 7: Commit Task 2**

Run:

```powershell
git add arguments.h parse_arguments.cpp arguments.cpp info.h help.cpp
git commit -m "Add CMap CLI and output names"
```

Stage only these paths plus any required conflict-resolution hunks.

---

### Task 3: Streaming CMap Analyzer

**Files:**
- Create: `titan_contact_map.h`
- Create: `titan_contact_map.cpp`
- Create: `test/titan_contact_map_smoke.cpp`

- [ ] **Step 1: Write analyzer smoke**

Create `test/titan_contact_map_smoke.cpp`:

```cpp
#include "Multi_MOL.h"
#include "arguments.h"
#include "info.h"
#include "titan_contact_map.h"
#include "titan_executor.h"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

namespace {

void require(bool ok, const char *message) {
  if (!ok) {
    std::cerr << "FAIL: " << message << "\n";
    std::exit(1);
  }
}

mol_info make_info(std::size_t resnr, const char *res, const char *atom,
                   std::size_t chainnum, const char *chain) {
  mol_info info(resnr, "", atom, chain, res);
  info.set_chainnum(chainnum);
  info.set_respos(resnr);
  info.set_internal_resnr(resnr);
  return info;
}

void append_atom(Multi_MOL &mol, MOL &frame, std::size_t resnr,
                 const char *res, const char *atom, double x, double y,
                 double z) {
  mol.info.push_back(make_info(resnr, res, atom, 0, "A"));
  frame.pos.x.push_back(x);
  frame.pos.y.push_back(y);
  frame.pos.z.push_back(z);
}

Multi_MOL build_fixture() {
  Multi_MOL mol("AMBER99SB");
  mol.set_mol_type("protein");

  MOL frame;
  frame.Box.set_a(40.0);
  frame.Box.set_b(40.0);
  frame.Box.set_c(40.0);
  frame.Box.set_box_type("orthorhombic");
  append_atom(mol, frame, 1, "ALA", "CA", 1.0, 1.0, 1.0);
  append_atom(mol, frame, 1, "ALA", "CB", 2.0, 1.0, 1.0);
  append_atom(mol, frame, 5, "LEU", "CA", 10.0, 1.0, 1.0);
  append_atom(mol, frame, 5, "LEU", "CB", 4.0, 1.0, 1.0);
  append_atom(mol, frame, 5, "LEU", "CG", 4.3, 1.0, 1.0);
  append_atom(mol, frame, 9, "SER", "CA", 20.0, 1.0, 1.0);
  append_atom(mol, frame, 9, "SER", "CB", 22.0, 1.0, 1.0);
  mol.m_mol.push_back(frame);
  return mol;
}

bool file_contains(const std::string &path, const std::string &needle) {
  std::ifstream in(path);
  std::ostringstream buffer;
  buffer << in.rdbuf();
  return buffer.str().find(needle) != std::string::npos;
}

} // namespace

int main() {
  arguments args;
  args.cmap.flag_cmap = true;
  args.cmap.mode = "sidechain-heavy";
  args.cmap.cutoff_ang = 4.5;
  args.cmap.domain = "intra";
  args.cmap.min_seqsep = 3;

  Multi_MOL mol = build_fixture();
  TitanContactMapAnalyzer analyzer;
  std::fstream log("test_cmap_smoke.log", std::ios::out);
  analyzer.configure(args.cmap);
  analyzer.init(mol, args, log);
  require(analyzer.enabled(), "analyzer enabled");

  TitanExecutor executor(1);
  analyzer.process_frame(mol, 0, executor);
  info output_info(args);
  std::vector<std::string> input_line{"fixture.pdb"};
  analyzer.finalize_trajectory(".", output_info, args, input_line, 0, log);

  require(file_contains("contact_probability_edges.csv", "ALA"),
          "edges include ALA");
  require(file_contains("contact_probability_edges.csv", "LEU"),
          "edges include LEU");
  require(file_contains("contact_probability_edges.csv", "1.000000"),
          "edges include full occupancy");
  require(file_contains("contact_probability_summary.dat",
                        "mode: sidechain-heavy"),
          "summary records mode");
  std::cout << "PASS titan_contact_map_smoke\n";
  return 0;
}
```

- [ ] **Step 2: Run analyzer smoke and verify it fails**

Run:

```powershell
C:\msys64\ucrt64\bin\clang++.exe -std=c++23 -I. test\titan_contact_map_smoke.cpp titan_contact_map.cpp titan_contact_primitives.cpp titan_executor.cpp -o test\titan_contact_map_smoke.exe
```

Expected: compile fails because `TitanContactMapAnalyzer` does not exist yet.

- [ ] **Step 3: Add analyzer header**

Create `titan_contact_map.h`:

```cpp
#ifndef TITAN_CONTACT_MAP_H
#define TITAN_CONTACT_MAP_H

#include "Multi_MOL.h"
#include "arguments.h"
#include "titan_contact_primitives.h"
#include "titan_neighbor_grid.h"

#include <cstddef>
#include <cstdint>
#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>

class TitanExecutor;
class info;

class TitanContactMapAnalyzer {
public:
  TitanContactMapAnalyzer() = default;

  void configure(const CMapArgs &args);
  void init(const Multi_MOL &mol1, const arguments &args, std::fstream &log);
  void process_frame(const Multi_MOL &mol1, std::size_t snap_idx,
                     TitanExecutor &executor);
  void finalize_trajectory(const std::string &base_dir, const info &info_output,
                           const arguments &args,
                           const std::vector<std::string> &input_line,
                           int input_name, std::fstream &log);

  [[nodiscard]] bool enabled() const noexcept { return enabled_; }
  [[nodiscard]] bool initialized() const noexcept { return initialized_; }

private:
  struct ContactAtom {
    std::size_t atom_index{};
    titan_contact::ResidueKey residue;
  };

  struct ContactRecord {
    titan_contact::ResiduePairKey key;
    std::uint64_t contact_frames{0};
    double min_distance_a{0.0};
    double distance_sum_a{0.0};
    std::uint64_t distance_samples{0};
  };

  void reset_trajectory_state();
  void prepare_atoms(const Multi_MOL &mol1, std::fstream &log);
  bool pair_allowed(const titan_contact::ResidueKey &lhs,
                    const titan_contact::ResidueKey &rhs) const;
  void export_summary(const std::string &filename, std::fstream &log) const;
  void export_edges(const std::string &filename, std::fstream &log) const;
  void export_matrix(const std::string &filename, std::fstream &log) const;

  bool requested_{false};
  bool enabled_{false};
  bool initialized_{false};
  bool warned_invalid_box_{false};
  titan_contact::ContactMode mode_{titan_contact::ContactMode::SidechainHeavy};
  std::string domain_{"both"};
  double cutoff_a_{4.5};
  int min_seqsep_{3};
  std::size_t atom_count_{0};
  std::uint64_t total_frames_seen_{0};
  std::uint64_t valid_frames_{0};
  std::uint64_t skipped_invalid_box_frames_{0};
  std::vector<ContactAtom> contact_atoms_;
  std::vector<unsigned char> active_atom_mask_;
  std::vector<titan_contact::ResidueKey> residues_;
  std::unordered_map<titan_contact::ResiduePairKey, ContactRecord,
                     titan_contact::ResiduePairKeyHash>
      records_;
  titan_neighbor::HalfShellGrid neighbor_grid_;
};

#endif // TITAN_CONTACT_MAP_H
```

- [ ] **Step 4: Implement analyzer**

Create `titan_contact_map.cpp` with these required behaviors:

```cpp
#include "titan_contact_map.h"

#include "info.h"
#include "titan_sq_core.h"
#include "titan_ui.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <set>

void TitanContactMapAnalyzer::configure(const CMapArgs &args) {
  requested_ = args.flag_cmap;
  mode_ = titan_contact::parse_contact_mode(args.mode);
  domain_ = args.domain;
  cutoff_a_ = args.cutoff_ang;
  min_seqsep_ = args.min_seqsep;
  enabled_ = requested_ && initialized_ && !contact_atoms_.empty();
}

void TitanContactMapAnalyzer::reset_trajectory_state() {
  warned_invalid_box_ = false;
  total_frames_seen_ = 0;
  valid_frames_ = 0;
  skipped_invalid_box_frames_ = 0;
  records_.clear();
}

void TitanContactMapAnalyzer::prepare_atoms(const Multi_MOL &mol1,
                                            std::fstream &log) {
  contact_atoms_.clear();
  residues_.clear();
  active_atom_mask_.assign(atom_count_, 0);
  std::set<std::pair<std::size_t, std::size_t>> seen_residues;

  for (std::size_t atom_idx = 0; atom_idx < mol1.info.size() &&
                                 atom_idx < atom_count_;
       ++atom_idx) {
    const mol_info &atom = mol1.info[atom_idx];
    if (!titan_contact::is_selected_contact_atom(atom, mode_))
      continue;
    titan_contact::ResidueKey residue =
        titan_contact::residue_key_from_atom(atom, residues_.size() + 1);
    const auto residue_id = std::make_pair(residue.chainnum, residue.resnr);
    if (seen_residues.insert(residue_id).second)
      residues_.push_back(residue);
    contact_atoms_.push_back(ContactAtom{atom_idx, residue});
    active_atom_mask_[atom_idx] = 1;
  }

  if (log.is_open()) {
    TitanUI::log_kv(log, "CMap Mode",
                    titan_contact::contact_mode_name(mode_));
    TitanUI::log_kv(log, "CMap Cutoff", std::to_string(cutoff_a_) + " A");
    TitanUI::log_kv(log, "CMap Domain", domain_);
    TitanUI::log_kv(log, "CMap Contact Atoms",
                    std::to_string(contact_atoms_.size()));
    TitanUI::log_kv(log, "CMap Residues", std::to_string(residues_.size()));
  }
}

void TitanContactMapAnalyzer::init(const Multi_MOL &mol1, const arguments &,
                                   std::fstream &log) {
  const std::size_t info_atoms = mol1.info.size();
  const std::size_t frame_atoms =
      mol1.m_mol.empty() ? 0 : mol1.m_mol.front().size();
  atom_count_ = (info_atoms > 0 && frame_atoms > 0)
                    ? std::min(info_atoms, frame_atoms)
                    : std::max(info_atoms, frame_atoms);
  initialized_ = true;
  prepare_atoms(mol1, log);
  reset_trajectory_state();
  enabled_ = requested_ && atom_count_ > 0 && !contact_atoms_.empty();
  if (requested_ && !enabled_) {
    TitanUI::print_warning(
        "CMap: no contact atoms matched the requested mode; module disabled.");
  }
}

bool TitanContactMapAnalyzer::pair_allowed(
    const titan_contact::ResidueKey &lhs,
    const titan_contact::ResidueKey &rhs) const {
  if (titan_contact::same_residue(lhs, rhs))
    return false;
  const bool intra = titan_contact::same_chain(lhs, rhs);
  if (domain_ == "intra" && !intra)
    return false;
  if (domain_ == "inter" && intra)
    return false;
  if (intra && min_seqsep_ > 0) {
    const auto diff =
        lhs.resnr > rhs.resnr ? lhs.resnr - rhs.resnr : rhs.resnr - lhs.resnr;
    if (diff < static_cast<std::size_t>(min_seqsep_))
      return false;
  }
  return true;
}

void TitanContactMapAnalyzer::process_frame(const Multi_MOL &mol1,
                                            std::size_t snap_idx,
                                            TitanExecutor &) {
  if (!enabled_ || snap_idx >= mol1.m_mol.size())
    return;
  ++total_frames_seen_;
  const MOL &frame = mol1.m_mol[snap_idx];
  if (!titan_sq_detail::has_positive_box_lengths(frame.Box)) {
    ++skipped_invalid_box_frames_;
    if (!warned_invalid_box_) {
      TitanUI::print_warning(
          "CMap: valid orthorhombic-style box lengths are required; invalid "
          "frames will be skipped.");
      warned_invalid_box_ = true;
    }
    return;
  }
  ++valid_frames_;

  std::unordered_map<titan_contact::ResiduePairKey, double,
                     titan_contact::ResiduePairKeyHash>
      frame_min_distances;
  const auto xs = frame.get_pos_x();
  const auto ys = frame.get_pos_y();
  const auto zs = frame.get_pos_z();
  titan_sq_detail::for_each_active_pair_within_cutoff(
      xs, ys, zs, frame.Box, active_atom_mask_, cutoff_a_, neighbor_grid_,
      [&](std::size_t lhs_atom, std::size_t rhs_atom, double distance_a) {
        const auto find_contact_atom = [](const std::vector<ContactAtom> &atoms,
                                          std::size_t atom_idx)
            -> const ContactAtom * {
          for (const auto &atom : atoms) {
            if (atom.atom_index == atom_idx)
              return &atom;
          }
          return nullptr;
        };
        const ContactAtom *lhs = find_contact_atom(contact_atoms_, lhs_atom);
        const ContactAtom *rhs = find_contact_atom(contact_atoms_, rhs_atom);
        if (!lhs || !rhs || !pair_allowed(lhs->residue, rhs->residue))
          return;
        titan_contact::ResiduePairKey key =
            titan_contact::make_ordered_pair(lhs->residue, rhs->residue);
        auto it = frame_min_distances.find(key);
        if (it == frame_min_distances.end() || distance_a < it->second)
          frame_min_distances[key] = distance_a;
      });

  for (const auto &[key, min_distance_a] : frame_min_distances) {
    ContactRecord &record = records_[key];
    if (record.contact_frames == 0)
      record.key = key;
    ++record.contact_frames;
    ++record.distance_samples;
    record.distance_sum_a += min_distance_a;
    if (record.min_distance_a == 0.0 || min_distance_a < record.min_distance_a)
      record.min_distance_a = min_distance_a;
  }
}
```

Then implement `export_summary`, `export_edges`, `export_matrix`, and `finalize_trajectory` with these exact schemas:

```text
contact_probability_summary.dat:
module: titan_cmap
mode: sidechain-heavy
cutoff_A: 4.500000
domain: intra
min_seqsep: 3
frames_seen: 1
valid_frames: 1
skipped_invalid_box_frames: 0
contact_atoms: 3
residues: 3
contact_edges: 1
```

```text
contact_probability_edges.csv:
chain_i,resnr_i,resname_i,chain_j,resnr_j,resname_j,probability,contact_frames,mean_distance_A,min_distance_A
```

```text
contact_probability_matrix.csv:
residue,chain:resnr:resname,...
chain:resnr:resname,0.000000,...
```

Important implementation detail: count each residue-pair at most once per frame.
Use the nearest atom-pair distance in that frame for `mean_distance_A` and
`min_distance_A`. Export `probability = contact_frames / valid_frames`.

- [ ] **Step 5: Run analyzer smoke**

Run:

```powershell
C:\msys64\ucrt64\bin\clang++.exe -std=c++23 -I. test\titan_contact_map_smoke.cpp titan_contact_map.cpp titan_contact_primitives.cpp titan_executor.cpp titan_sq_core.cpp -o test\titan_contact_map_smoke.exe
.\test\titan_contact_map_smoke.exe
```

Expected: `PASS titan_contact_map_smoke`.

- [ ] **Step 6: Commit Task 3**

Run:

```powershell
git add titan_contact_map.h titan_contact_map.cpp test/titan_contact_map_smoke.cpp
git commit -m "Add streaming CMap analyzer"
```

---

### Task 4: Engine Integration

**Files:**
- Modify: `impl_engine.cpp`

- [ ] **Step 1: Include and declare analyzer**

Add:

```cpp
#include "titan_contact_map.h"
```

Near other analyzer declarations:

```cpp
TitanContactMapAnalyzer cmap_analyzer;
std::chrono::milliseconds time_cmap{0};
```

- [ ] **Step 2: Configure and initialize**

Near `sq_analyzer.configure(Arguments.sq);`:

```cpp
cmap_analyzer.configure(Arguments.cmap);
```

Near `sq_analyzer.init(pro_in, Arguments, log);`:

```cpp
if (Arguments.cmap.flag_cmap) {
  cmap_analyzer.init(pro_in, Arguments, log);
}
```

- [ ] **Step 3: Process frames in both streaming branches**

In every branch where `sq_analyzer.process_frame(...)` is called, add the same pattern for CMap:

```cpp
if (cmap_analyzer.enabled()) {
  auto start_cmap = std::chrono::steady_clock::now();
  cmap_analyzer.process_frame(pro_in, chunk_snap_idx, executor);
  time_cmap += std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::steady_clock::now() - start_cmap);
}
```

For the single-frame fallback branch, use frame index `0`.

- [ ] **Step 4: Finalize per trajectory**

Near SQ finalize:

```cpp
if (cmap_analyzer.enabled()) {
  cmap_analyzer.finalize_trajectory(base_dir, info_output, Arguments,
                                    input_line, input_name, log);
}
```

- [ ] **Step 5: Add timing output**

Where timing summaries list module times, add:

```cpp
TitanUI::log_time(log, "CMap Analysis",
                  (time_cmap - time_cmap_base).count() / 1000.0);
```

If the surrounding code uses per-input base variables, define `const auto time_cmap_base = time_cmap;` beside the other base timers.

- [ ] **Step 6: Build native**

Run:

```powershell
$env:TITAN_NO_PAUSE='1'
.\compile.bat native all
```

Expected: native build succeeds and `test\bin\native\titan.exe` is updated.

- [ ] **Step 7: Commit Task 4**

Run:

```powershell
git add impl_engine.cpp
git commit -m "Wire CMap analyzer into TITAN runtime"
```

---

### Task 5: End-to-End Smoke Script

**Files:**
- Create: `test/cmap_smoke.ps1`

- [ ] **Step 1: Add smoke script**

Create `test/cmap_smoke.ps1`:

```powershell
param(
  [string]$TitanExe = ".\test\bin\native\titan.exe"
)

$ErrorActionPreference = "Stop"
$repo = Resolve-Path (Join-Path $PSScriptRoot "..")
$work = Join-Path $repo "test\cmap_inputs_local"
New-Item -ItemType Directory -Force -Path $work | Out-Null

$arg = Join-Path $work "arg_cmap.in"
@"
-f=example\2\md1_short.pdb
-1pdb
-nocal
-cmap
-cmap-mode=sidechain-heavy
-cmap-cut=4.5
-cmap-domain=intra
-cmap-min-seqsep=3
"@ | Set-Content -Encoding ASCII $arg

Push-Location $repo
try {
  & $TitanExe -a $arg
  if ($LASTEXITCODE -ne 0) {
    throw "TITAN exited with code $LASTEXITCODE"
  }
  foreach ($file in @(
      "contact_probability_matrix.csv",
      "contact_probability_edges.csv",
      "contact_probability_summary.dat")) {
    if (!(Test-Path -LiteralPath $file)) {
      throw "Missing expected output $file"
    }
  }
  Select-String -Path "contact_probability_summary.dat" -Pattern "module: titan_cmap" | Out-Null
  Select-String -Path "contact_probability_summary.dat" -Pattern "mode: sidechain-heavy" | Out-Null
  Select-String -Path "contact_probability_edges.csv" -Pattern "probability" | Out-Null
  Write-Host "PASS cmap_smoke"
}
finally {
  Pop-Location
}
```

The fixture `example\2\md1_short.pdb` exists in this checkout and keeps the
smoke fast by using single-PDB trajectory mode.

- [ ] **Step 2: Run smoke**

Run:

```powershell
.\test\cmap_smoke.ps1 -TitanExe .\test\bin\native\titan.exe
```

Expected: `PASS cmap_smoke`.

- [ ] **Step 3: Commit Task 5**

Run:

```powershell
git add test/cmap_smoke.ps1
git commit -m "Add CMap end-to-end smoke"
```

---

### Task 6: Documentation and Roadmap Sync

**Files:**
- Modify: `manual.txt`
- Modify: `AGENT_ONBOARDING.md`
- Modify: `TITAN_STANDARD_MD_ANALYSIS_ROADMAP.md`

- [ ] **Step 1: Add manual section**

Add a concise public section:

```text
2.X Side-chain Contact Probability Map

  Purpose:
    Compute residue-residue contact probabilities from a trajectory. The default
    definition is side-chain heavy-atom contact with glycine CA fallback.

  Main inputs:
    -cmap
    -cmap-mode <sidechain-heavy|cb|ca|bb>
    -cmap-cut <A>          default 4.5
    -cmap-domain <intra|inter|both>
    -cmap-min-seqsep <N>   default 3

  Typical outputs:
    contact_probability_matrix.csv
    contact_probability_edges.csv
    contact_probability_summary.dat

  Notes:
    - `sidechain-heavy` excludes backbone atoms and hydrogens.
    - `cb` and `ca` are available for fast/legacy anchor maps.
    - `bb` is a diagnostic backbone mode.
    - This module is independent of IFP; it reports untyped contact
      probability rather than H-bond/salt/hydrophobic categories.
```

- [ ] **Step 2: Update onboarding**

Add one bullet to the analysis capability list:

```text
- `-cmap` exports side-chain residue contact probability matrices and sparse
  contact edge lists; use it when the workflow needs untyped residue-residue
  contact frequency rather than DCCM correlation or typed IFP occupancy.
```

- [ ] **Step 3: Mark roadmap status**

In `TITAN_STANDARD_MD_ANALYSIS_ROADMAP.md`, under Slice A, add:

```text
Status:

- [done] Slice A v1 shipped `-cmap` side-chain contact probability output:
  parser/options, shared contact primitive, streaming analyzer, matrix/edge/
  summary outputs, native build, and local smoke verification.
```

Only add the `[done]` line after Tasks 1-5 have actually passed.

- [ ] **Step 4: Commit Task 6**

Run:

```powershell
git add manual.txt AGENT_ONBOARDING.md TITAN_STANDARD_MD_ANALYSIS_ROADMAP.md
git commit -m "Document CMap analysis module"
```

---

### Task 7: Final Verification Matrix

**Files:**
- No source edits expected

- [ ] **Step 1: Run primitive smoke**

Run:

```powershell
.\test\titan_contact_primitives_smoke.exe
```

Expected: `PASS titan_contact_primitives_smoke`.

- [ ] **Step 2: Run analyzer smoke**

Run:

```powershell
.\test\titan_contact_map_smoke.exe
```

Expected: `PASS titan_contact_map_smoke`.

- [ ] **Step 3: Run native build**

Run:

```powershell
$env:TITAN_NO_PAUSE='1'
.\compile.bat native all
```

Expected: build finishes and snapshots `test\bin\native\titan.exe`.

- [ ] **Step 4: Run portable build**

Run:

```powershell
$env:TITAN_NO_PAUSE='1'
.\compile.bat portable all
```

Expected: build finishes and snapshots `test\bin\portable\titan.exe`.

- [ ] **Step 5: Run native and portable CMap smoke**

Run:

```powershell
.\test\cmap_smoke.ps1 -TitanExe .\test\bin\native\titan.exe
.\test\cmap_smoke.ps1 -TitanExe .\test\bin\portable\titan.exe
```

Expected: both runs print `PASS cmap_smoke`.

- [ ] **Step 6: Verify default behavior remains inert**

Run a small existing benchmark case without `-cmap`:

```powershell
.\test\benchmark_harness.ps1 -Cases small -Profile native
```

Expected: harness completes. If the harness requires an array for `-Cases`, run:

```powershell
.\test\benchmark_harness.ps1 -Cases @("small") -Profile native
```

Expected: no new contact probability files are produced unless the case explicitly enables `-cmap`.

- [ ] **Step 7: Commit final verification note if needed**

If verification changes only docs with evidence, commit those docs. If no files changed, do not create an empty commit.

---

## Self-Review

- Spec coverage: Slice A requirements are covered: side-chain classifier, glycine fallback, shared primitive, parser/options/logging, neighbor-grid occupancy, matrix/edge/summary exports, smoke tests, docs/help.
- Scope control: `-qnative`, TICA feature export, graph-source `contact`, MBAR weights, and per-frame contact bits remain future slices.
- Compatibility: default behavior is inert because all runtime work is gated behind `args.cmap.flag_cmap`.
- Risk to revisit during implementation: the first production version should
  profile the O(contact-atoms) lookup inside the pair lambda and replace it with
  an atom-index lookup table if the smoke or benchmark shows it is measurable.
