# LLPS Paper Fractal Dense Phase Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add paper-ready fractal morphology outputs for Martini protein/RNA LLPS trajectories, measured from bead/residue dense-phase occupancy rather than only from dimer macro-nodes.

**Architecture:** Extend the existing `-llps-paper` bundle with a focused fractal accumulator inside the LLPS DBSCAN/paper path. The accumulator receives dense-phase nodes already produced by `TitanLLPSCluster`, computes scale-dependent bead/residue occupancy and mass-radius fits, and exports citeable files named `llps_paper_fractal_*`. The existing macro-node LLPS analyzer remains unchanged except for documentation that distinguishes macro percolation from dense-phase fractal morphology.

**Tech Stack:** C++17, existing TITAN LLPS modules, existing PowerShell smoke scripts, existing `mingw32-make PROFILE=native TARGET=titan.exe all -j1` native build gate.

---

## Scientific Contract

This feature measures morphology, not thermodynamic phase boundaries. A valid paper interpretation must combine these outputs with existing `llps_phase_partition_summary.dat`, `llps_convergence_report.dat`, contact dynamics, RNA compaction/topology outputs, and any thermodynamic convergence evidence available for the run.

Primary analysis target: dense-phase bead/residue occupancy from `TitanLLPSCluster` nodes and dense membership. This is the right target for `10 dimer N protein + 1 RNA` because `-llps-dimer-size=2` produces only about 10 protein macro-nodes, which is too small for stable per-frame fractal fitting.

Default behavior:

- `-llps-paper=true` enables fractal outputs automatically.
- `-llps-paper-fractal=<auto|true|false>` controls the feature explicitly.
- `auto` follows `-llps-paper`.
- `true` enables `-llps` and the DBSCAN paper path if needed.
- `false` suppresses fractal outputs even when `-llps-paper=true`.

Output files:

- `llps_paper_fractal_timeseries.csv`
- `llps_paper_fractal_scales.csv`
- `llps_paper_fractal_summary.dat`
- `llps_paper_fractal_readme.md`

Core columns:

- `frame`
- `time_ns`
- `dense_nodes`
- `dense_protein_nodes`
- `dense_rna_nodes`
- `df_mass_radius`
- `df_box_count`
- `kf_mass_radius`
- `fit_r2_mass_radius`
- `fit_r2_box_count`
- `scale_count_mass_radius`
- `scale_count_box_count`
- `dense_radius_nm`
- `dense_rg_nm`
- `mean_local_valency`
- `protein_rna_neighbor_fraction`
- `valid`
- `invalid_reason`

Interpretation guardrails:

- If `dense_nodes < 8`, mark the frame invalid with `invalid_reason=dense_nodes_lt_8`.
- If fewer than 3 usable scales exist for either fit, mark that fit as invalid with `nan` and record `scale_count_*`.
- If fitted `df` is outside `[0.5, 3.2]`, keep the raw value in the scales file but mark the timeseries frame invalid with `invalid_reason=df_out_of_range`.
- Summary reports block mean and SEM only across valid frames.

## File Structure

- Modify: `src/arguments.h`
  - Add `llps_paper_fractal` mode string to `LlpsArgs`.
- Modify: `src/parse_arguments.cpp`
  - Parse `-llps-paper-fractal` and `--llps-paper-fractal`.
  - Make `true` enable `flag_llps` and `llps_paper`.
- Modify: `src/info.h`
  - Add output filename helpers for the four citeable fractal outputs.
- Modify: `src/titan_llps_cluster.h`
  - Add `LLPSPaperFractalConfig`, `LLPSPaperFractalFrame`, `LLPSPaperFractalScaleRow`, `LLPSPaperFractalAccumulator`.
  - Add `LLPSPaperConfig::fractal_enabled`.
- Modify: `src/titan_llps_paper.cpp`
  - Implement the fractal accumulator and call it from `LLPSPaperAccumulator::process_frame`.
  - Export the new files from `LLPSPaperAccumulator::export_outputs`.
- Modify: `src/impl_engine.cpp`
  - Pass the parsed fractal mode into the LLPS paper config.
  - Log the effective fractal mode under `LLPS Run Context`.
- Modify: `src/help.cpp` and `manual.txt`
  - Document `-llps-paper-fractal` and the `llps_paper_fractal_*` outputs.
- Modify: `plot_llps_paper.py`
  - Load and include the new fractal summary in its generated report when present.
- Create: `test/llps_paper_fractal_smoke.cpp`
  - Unit smoke for scale fitting and output contracts without requiring a full trajectory.
- Create: `test/llps_paper_fractal_smoke.ps1`
  - Compile and run the smoke with the same local pattern as existing LLPS tests.

---

## Phase 1: CLI And Output Contract

### Task 1: Add Arguments And Filename Helpers

**Files:**
- Modify: `src/arguments.h`
- Modify: `src/parse_arguments.cpp`
- Modify: `src/info.h`
- Test: `test/llps_paper_fractal_args_smoke.cpp`
- Test: `test/llps_paper_fractal_args_smoke.ps1`

- [ ] **Step 1: Add the argument field**

In `src/arguments.h`, add this field immediately after `llps_paper_geometry`:

```cpp
  std::string llps_paper_fractal{"auto"};
```

- [ ] **Step 2: Add filename helpers**

In `src/info.h`, add these helpers after `get_llps_paper_summary_json_filename()`:

```cpp
  string get_llps_paper_fractal_timeseries_filename() const {
    return "llps_paper_fractal_timeseries.csv";
  }
  string get_llps_paper_fractal_scales_filename() const {
    return "llps_paper_fractal_scales.csv";
  }
  string get_llps_paper_fractal_summary_filename() const {
    return "llps_paper_fractal_summary.dat";
  }
  string get_llps_paper_fractal_readme_filename() const {
    return "llps_paper_fractal_readme.md";
  }
```

- [ ] **Step 3: Parse the option**

In `src/parse_arguments.cpp`, add this block after parsing `-llps-paper-geometry`:

```cpp
    if (consume_option_value(arg,
                             {"-llps-paper-fractal",
                              "--llps-paper-fractal"},
                             llps_value, i)) {
      llps.llps_paper_fractal =
          require_choice(llps_value, {"auto", "true", "false", "yes", "no",
                                      "on", "off", "1", "0"},
                         "-llps-paper-fractal");
      if (llps.llps_paper_fractal == "true" ||
          llps.llps_paper_fractal == "yes" ||
          llps.llps_paper_fractal == "on" ||
          llps.llps_paper_fractal == "1") {
        llps.flag_llps = true;
        llps.llps_paper = true;
        llps.llps_contacts = true;
        llps.llps_dynamics = true;
      }
    }
```

- [ ] **Step 4: Write the failing parser smoke**

Create `test/llps_paper_fractal_args_smoke.cpp`:

```cpp
#include "arguments.h"

#include <iostream>
#include <string>
#include <vector>

std::string arguments::set_Buildnumber() { return "llps_paper_fractal_args_smoke"; }
std::string arguments::get_build_number() const { return "llps_paper_fractal_args_smoke"; }

static void fail(const std::string &message) {
  std::cerr << "FAIL llps_paper_fractal_args_smoke: " << message << "\n";
  std::exit(1);
}

int main() {
  arguments args;
  std::vector<std::string> argv = {
      "titan.exe", "-llps-paper-fractal=true", "-llps-paper-profile-bin", "0.25"};
  args.parse_argv(argv);
  if (!args.llps.flag_llps)
    fail("fractal=true should enable LLPS");
  if (!args.llps.llps_paper)
    fail("fractal=true should enable paper bundle");
  if (!args.llps.llps_contacts)
    fail("fractal=true should enable contacts");
  if (!args.llps.llps_dynamics)
    fail("fractal=true should enable dynamics");
  if (args.llps.llps_paper_fractal != "true")
    fail("fractal mode mismatch");
  if (args.llps.llps_paper_profile_bin_nm != 0.25)
    fail("paper profile bin parse regression");
  std::cout << "PASS llps_paper_fractal_args_smoke\n";
  return 0;
}
```

Create `test/llps_paper_fractal_args_smoke.ps1`:

```powershell
param(
  [string]$WorkRoot = (Join-Path $PSScriptRoot "tmp\llps_paper_fractal_args_smoke")
)

$ErrorActionPreference = "Stop"
$repo = Split-Path $PSScriptRoot -Parent
$full = [System.IO.Path]::GetFullPath($WorkRoot)
$allowed = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot "tmp"))
if (-not $full.StartsWith($allowed, [System.StringComparison]::OrdinalIgnoreCase)) {
  throw "Refusing path outside test tmp root: $full"
}
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null
$exe = Join-Path $WorkRoot "llps_paper_fractal_args_smoke.exe"
g++ -std=c++17 -O2 -I "$repo\src" `
  "$repo\test\llps_paper_fractal_args_smoke.cpp" `
  "$repo\src\arguments.cpp" "$repo\src\parse_arguments.cpp" `
  "$repo\src\folder.cpp" "$repo\src\titan_ui.cpp" `
  -o $exe
if ($LASTEXITCODE -ne 0) { throw "Failed to compile fractal args smoke" }
& $exe
if ($LASTEXITCODE -ne 0) { throw "Fractal args smoke failed" }
Write-Host "PASS llps_paper_fractal_args_smoke"
```

- [ ] **Step 5: Run the failing smoke**

Run:

```powershell
.\test\llps_paper_fractal_args_smoke.ps1
```

Expected before implementation: compile failure or option parse failure referencing `llps_paper_fractal`.

- [ ] **Step 6: Run the passing smoke**

Run:

```powershell
.\test\llps_paper_fractal_args_smoke.ps1
```

Expected after implementation:

```text
PASS llps_paper_fractal_args_smoke
```

- [ ] **Step 7: Commit Phase 1**

Run:

```powershell
git add src\arguments.h src\parse_arguments.cpp src\info.h test\llps_paper_fractal_args_smoke.cpp test\llps_paper_fractal_args_smoke.ps1
git commit -m "feat: add LLPS paper fractal CLI contract"
```

---

## Phase 2: Dense-Phase Fractal Accumulator

### Task 2: Add Data Structures And Fit Helpers

**Files:**
- Modify: `src/titan_llps_cluster.h`
- Modify: `src/titan_llps_paper.cpp`
- Test: `test/llps_paper_fractal_smoke.cpp`
- Test: `test/llps_paper_fractal_smoke.ps1`

- [ ] **Step 1: Add public result structures**

In `src/titan_llps_cluster.h`, add these declarations before `class LLPSPaperAccumulator`:

```cpp
struct LLPSPaperFractalConfig {
  bool enabled{false};
  double min_scale_nm{0.5};
  double max_scale_fraction{0.75};
};

struct LLPSPaperFractalFrame {
  int frame{0};
  double time_ns{0.0};
  int dense_nodes{0};
  int dense_protein_nodes{0};
  int dense_rna_nodes{0};
  double df_mass_radius{std::numeric_limits<double>::quiet_NaN()};
  double df_box_count{std::numeric_limits<double>::quiet_NaN()};
  double kf_mass_radius{std::numeric_limits<double>::quiet_NaN()};
  double fit_r2_mass_radius{0.0};
  double fit_r2_box_count{0.0};
  int scale_count_mass_radius{0};
  int scale_count_box_count{0};
  double dense_radius_nm{0.0};
  double dense_rg_nm{0.0};
  double mean_local_valency{0.0};
  double protein_rna_neighbor_fraction{0.0};
  bool valid{false};
  std::string invalid_reason{"not_evaluated"};
};

struct LLPSPaperFractalScaleRow {
  int frame{0};
  double time_ns{0.0};
  std::string method{"mass_radius"};
  double scale_nm{0.0};
  double value{0.0};
  double log_scale{0.0};
  double log_value{0.0};
};
```

Add this member to `LLPSPaperConfig`:

```cpp
  bool fractal_enabled{false};
```

- [ ] **Step 2: Add the accumulator class declaration**

In `src/titan_llps_cluster.h`, add this class before `class LLPSPaperAccumulator`:

```cpp
class LLPSPaperFractalAccumulator {
public:
  void configure(const LLPSPaperFractalConfig &config);
  bool enabled() const { return config_.enabled; }
  void process_frame(double time_ns, int frame,
                     const std::vector<LLPSNode> &nodes,
                     const std::vector<bool> &is_dense,
                     int giant_cluster_id, const box &Box);
  void export_outputs(const info &info_output, std::fstream &log,
                      int requested_blocks) const;
  const std::vector<LLPSPaperFractalFrame> &frames() const { return frames_; }
  const std::vector<LLPSPaperFractalScaleRow> &scales() const {
    return scales_;
  }

#ifdef TITAN_LLPS_TESTING
  static LLPSPaperFractalFrame debug_fit_points(
      const std::vector<std::array<double, 3>> &points_nm);
#endif

private:
  LLPSPaperFractalConfig config_;
  std::vector<LLPSPaperFractalFrame> frames_;
  std::vector<LLPSPaperFractalScaleRow> scales_;
};
```

Add this private member to `LLPSPaperAccumulator`:

```cpp
  LLPSPaperFractalAccumulator fractal_;
```

- [ ] **Step 3: Implement least-squares helpers**

In `src/titan_llps_paper.cpp`, add this helper in the anonymous namespace:

```cpp
struct LinearFitResult {
  double slope{std::numeric_limits<double>::quiet_NaN()};
  double intercept{std::numeric_limits<double>::quiet_NaN()};
  double r2{0.0};
  int n{0};
};

LinearFitResult fit_line(const std::vector<std::pair<double, double>> &xy) {
  LinearFitResult out;
  out.n = static_cast<int>(xy.size());
  if (xy.size() < 3)
    return out;
  double sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0, syy = 0.0;
  for (const auto &p : xy) {
    sx += p.first;
    sy += p.second;
    sxx += p.first * p.first;
    sxy += p.first * p.second;
    syy += p.second * p.second;
  }
  const double n = static_cast<double>(xy.size());
  const double denom = n * sxx - sx * sx;
  if (std::abs(denom) < kEps)
    return out;
  out.slope = (n * sxy - sx * sy) / denom;
  out.intercept = (sy - out.slope * sx) / n;
  const double ss_tot = syy - sy * sy / n;
  double ss_res = 0.0;
  for (const auto &p : xy) {
    const double yhat = out.slope * p.first + out.intercept;
    const double d = p.second - yhat;
    ss_res += d * d;
  }
  out.r2 = ss_tot > kEps ? std::max(0.0, 1.0 - ss_res / ss_tot) : 0.0;
  return out;
}
```

- [ ] **Step 4: Implement `LLPSPaperFractalAccumulator::configure`**

In `src/titan_llps_paper.cpp`, add:

```cpp
void LLPSPaperFractalAccumulator::configure(
    const LLPSPaperFractalConfig &config) {
  config_ = config;
  config_.min_scale_nm = std::max(0.1, config_.min_scale_nm);
  config_.max_scale_fraction =
      std::clamp(config_.max_scale_fraction, 0.25, 1.0);
  frames_.clear();
  scales_.clear();
}
```

- [ ] **Step 5: Implement the dense point extraction**

In `LLPSPaperFractalAccumulator::process_frame`, start with this body:

```cpp
void LLPSPaperFractalAccumulator::process_frame(
    double time_ns, int frame, const std::vector<LLPSNode> &nodes,
    const std::vector<bool> &is_dense, int giant_cluster_id, const box &Box) {
  if (!config_.enabled)
    return;

  LLPSPaperFractalFrame row;
  row.frame = frame;
  row.time_ns = time_ns;

  if (nodes.empty() || is_dense.size() != nodes.size() || giant_cluster_id <= 0) {
    row.invalid_reason = "no_dense_cluster";
    frames_.push_back(row);
    return;
  }

  std::vector<std::array<double, 3>> points_nm;
  points_nm.reserve(nodes.size());
  for (std::size_t i = 0; i < nodes.size(); ++i) {
    if (!is_dense[i])
      continue;
    const bool protein = nodes[i].type == 1;
    row.dense_protein_nodes += protein ? 1 : 0;
    row.dense_rna_nodes += protein ? 0 : 1;
    points_nm.push_back({nodes[i].x / 10.0, nodes[i].y / 10.0,
                         nodes[i].z / 10.0});
  }
  row.dense_nodes = static_cast<int>(points_nm.size());
  if (row.dense_nodes < 8) {
    row.invalid_reason = "dense_nodes_lt_8";
    frames_.push_back(row);
    return;
  }

  row = compute_fractal_frame_from_points(points_nm, row, frame, time_ns);
  frames_.push_back(row);
}
```

Add `compute_fractal_frame_from_points` as a private anonymous-namespace helper in the same file. It must:

- Compute center of mass from `points_nm`.
- Compute `dense_rg_nm = sqrt(mean(r^2))`.
- Compute `dense_radius_nm = max(r)`.
- Build mass-radius scales from radii at `25%, 35%, 50%, 70%, 90%, 100%` of `dense_radius_nm`.
- For each radius, count points within that radius.
- Fit `log(count)` vs `log(radius)` with `fit_line`.
- Set `df_mass_radius = slope`, `kf_mass_radius = exp(intercept)`, and `fit_r2_mass_radius`.
- Build box-counting scales from powers of two divisions of the dense bounding box.
- For each grid scale, count occupied boxes and fit `log(Nbox)` vs `log(1/scale_nm)`.
- Set `df_box_count = slope` and `fit_r2_box_count`.
- Mark valid when both fits have at least 3 scales and `df_mass_radius` is in `[0.5, 3.2]`.

- [ ] **Step 6: Add the unit smoke**

Create `test/llps_paper_fractal_smoke.cpp`:

```cpp
#define TITAN_LLPS_TESTING 1

#include "titan_llps_cluster.h"

#include <array>
#include <cmath>
#include <iostream>
#include <vector>

static void fail(const char *message) {
  std::cerr << "FAIL llps_paper_fractal_smoke: " << message << "\n";
  std::exit(1);
}

int main() {
  std::vector<std::array<double, 3>> points;
  for (int x = 0; x < 4; ++x) {
    for (int y = 0; y < 4; ++y) {
      for (int z = 0; z < 4; ++z) {
        points.push_back({0.25 * x, 0.25 * y, 0.25 * z});
      }
    }
  }
  const auto result = LLPSPaperFractalAccumulator::debug_fit_points(points);
  if (!result.valid)
    fail("compact cubic dense phase should be valid");
  if (result.scale_count_mass_radius < 3)
    fail("mass-radius scale count too small");
  if (result.scale_count_box_count < 3)
    fail("box-count scale count too small");
  if (result.df_mass_radius < 2.0 || result.df_mass_radius > 3.3)
    fail("mass-radius df outside compact 3D range");
  if (result.df_box_count < 2.0 || result.df_box_count > 3.3)
    fail("box-count df outside compact 3D range");
  if (result.fit_r2_mass_radius < 0.85)
    fail("mass-radius fit r2 too low");
  if (result.fit_r2_box_count < 0.85)
    fail("box-count fit r2 too low");
  std::cout << "PASS llps_paper_fractal_smoke\n";
  return 0;
}
```

Create `test/llps_paper_fractal_smoke.ps1`:

```powershell
param(
  [string]$WorkRoot = (Join-Path $PSScriptRoot "tmp\llps_paper_fractal_smoke")
)

$ErrorActionPreference = "Stop"
$repo = Split-Path $PSScriptRoot -Parent
$full = [System.IO.Path]::GetFullPath($WorkRoot)
$allowed = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot "tmp"))
if (-not $full.StartsWith($allowed, [System.StringComparison]::OrdinalIgnoreCase)) {
  throw "Refusing path outside test tmp root: $full"
}
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null
$exe = Join-Path $WorkRoot "llps_paper_fractal_smoke.exe"
g++ -std=c++17 -O2 -DTITAN_LLPS_TESTING -I "$repo\src" `
  "$repo\test\llps_paper_fractal_smoke.cpp" `
  "$repo\src\titan_llps_paper.cpp" "$repo\src\titan_ui.cpp" `
  -o $exe
if ($LASTEXITCODE -ne 0) { throw "Failed to compile fractal smoke" }
& $exe
if ($LASTEXITCODE -ne 0) { throw "Fractal smoke failed" }
Write-Host "PASS llps_paper_fractal_smoke"
```

- [ ] **Step 7: Run the failing smoke**

Run:

```powershell
.\test\llps_paper_fractal_smoke.ps1
```

Expected before implementation: compile failure referencing `LLPSPaperFractalAccumulator`.

- [ ] **Step 8: Run the passing smoke**

Run:

```powershell
.\test\llps_paper_fractal_smoke.ps1
```

Expected after implementation:

```text
PASS llps_paper_fractal_smoke
```

- [ ] **Step 9: Commit Phase 2**

Run:

```powershell
git add src\titan_llps_cluster.h src\titan_llps_paper.cpp test\llps_paper_fractal_smoke.cpp test\llps_paper_fractal_smoke.ps1
git commit -m "feat: compute LLPS paper fractal dense-phase descriptors"
```

---

## Phase 3: Wire Into LLPS Paper Bundle

### Task 3: Call The Fractal Accumulator From The Existing Paper Path

**Files:**
- Modify: `src/titan_llps_paper.cpp`
- Modify: `src/impl_engine.cpp`
- Test: `test/llps_paper_fractal_integration_smoke.ps1`

- [ ] **Step 1: Configure the nested accumulator**

In `LLPSPaperAccumulator::configure`, after normalizing `config_`, add:

```cpp
  fractal_.configure(
      {config_.enabled && config_.fractal_enabled, config_.profile_bin_nm, 0.75});
```

- [ ] **Step 2: Process frames**

At the end of `LLPSPaperAccumulator::process_frame`, before returning, add:

```cpp
  fractal_.process_frame(row.time_ns, frame_idx, nodes, is_dense,
                         giant_cluster_id, Box);
```

For the early `first_dense < 0 || giant_cluster_id <= 0` branch, add the same call before `frames_.push_back(row); return;` so invalid frames are represented in `llps_paper_fractal_timeseries.csv`.

- [ ] **Step 3: Export outputs**

In `LLPSPaperAccumulator::export_outputs`, after existing paper outputs are written, add:

```cpp
  fractal_.export_outputs(info_output, log, config_.requested_blocks);
```

- [ ] **Step 4: Pass effective config from engine**

In `src/impl_engine.cpp`, replace the existing `configure_paper` call:

```cpp
      dbscan_llps->configure_paper(
          {Arguments.llps.llps_paper,
           Arguments.llps.llps_paper_profile_bin_nm,
           Arguments.llps.llps_paper_blocks,
           Arguments.llps.llps_paper_geometry});
```

with:

```cpp
      const bool llps_fractal_enabled =
          Arguments.llps.llps_paper_fractal == "auto"
              ? Arguments.llps.llps_paper
              : (Arguments.llps.llps_paper_fractal == "true" ||
                 Arguments.llps.llps_paper_fractal == "yes" ||
                 Arguments.llps.llps_paper_fractal == "on" ||
                 Arguments.llps.llps_paper_fractal == "1");
      dbscan_llps->configure_paper(
          {Arguments.llps.llps_paper,
           Arguments.llps.llps_paper_profile_bin_nm,
           Arguments.llps.llps_paper_blocks,
           Arguments.llps.llps_paper_geometry,
           llps_fractal_enabled});
```

Also log:

```cpp
      TitanUI::log_kv(log, "LLPS Paper Fractal",
                      llps_fractal_enabled ? "enabled" : "disabled");
```

- [ ] **Step 5: Create the integration smoke**

Create `test/llps_paper_fractal_integration_smoke.ps1`:

```powershell
param(
  [string]$TitanExe = (Join-Path (Split-Path $PSScriptRoot -Parent) "titan.exe"),
  [string]$WorkRoot = (Join-Path $PSScriptRoot "tmp\llps_paper_fractal_integration_smoke")
)

$ErrorActionPreference = "Stop"
$repo = Split-Path $PSScriptRoot -Parent
$full = [System.IO.Path]::GetFullPath($WorkRoot)
$allowed = [System.IO.Path]::GetFullPath((Join-Path $PSScriptRoot "tmp"))
if (-not $full.StartsWith($allowed, [System.StringComparison]::OrdinalIgnoreCase)) {
  throw "Refusing path outside test tmp root: $full"
}
New-Item -ItemType Directory -Force -Path $WorkRoot | Out-Null

$arg = Join-Path $WorkRoot "arg.in"
@"
-llps-paper=true
-llps-paper-fractal=true
-llps-dimer-size=2
-llps-rcut=0.7
-llps-dbscan-eps=15
-llps-dbscan-minpts=3
-llps-sticker=BB,SC1,SC2,SC3,SC4,BB1,BB2,BB3
"@ | Set-Content -Path $arg -Encoding ASCII

Push-Location $WorkRoot
try {
  & $TitanExe -a $arg
  if ($LASTEXITCODE -ne 0) { throw "TITAN exited with code $LASTEXITCODE" }
  foreach ($name in @(
    "llps_paper_fractal_timeseries.csv",
    "llps_paper_fractal_scales.csv",
    "llps_paper_fractal_summary.dat",
    "llps_paper_fractal_readme.md")) {
    if (-not (Test-Path $name)) { throw "Missing $name" }
  }
  $header = Get-Content "llps_paper_fractal_timeseries.csv" -First 1
  if ($header -notmatch "df_mass_radius" -or $header -notmatch "df_box_count") {
    throw "Timeseries header missing fractal columns: $header"
  }
}
finally {
  Pop-Location
}
Write-Host "PASS llps_paper_fractal_integration_smoke"
```

If the repository lacks a tiny LLPS trajectory fixture usable by `titan.exe`, convert this script during implementation into a C++ integration smoke that seeds `LLPSPaperAccumulator` directly with synthetic `LLPSNode` rows. Keep the output filenames and header checks identical.

- [ ] **Step 6: Run the integration smoke**

Run:

```powershell
.\test\llps_paper_fractal_integration_smoke.ps1
```

Expected:

```text
PASS llps_paper_fractal_integration_smoke
```

- [ ] **Step 7: Commit Phase 3**

Run:

```powershell
git add src\titan_llps_paper.cpp src\impl_engine.cpp test\llps_paper_fractal_integration_smoke.ps1
git commit -m "feat: wire fractal descriptors into LLPS paper bundle"
```

---

## Phase 4: Export Summary And Citeable README

### Task 4: Write Paper-Ready Files

**Files:**
- Modify: `src/titan_llps_paper.cpp`
- Modify: `plot_llps_paper.py`
- Test: `test/llps_paper_fractal_smoke.ps1`

- [ ] **Step 1: Export timeseries**

Implement `LLPSPaperFractalAccumulator::export_outputs` with this timeseries header:

```cpp
timeseries << "frame,time_ns,dense_nodes,dense_protein_nodes,dense_rna_nodes,"
              "df_mass_radius,df_box_count,kf_mass_radius,"
              "fit_r2_mass_radius,fit_r2_box_count,"
              "scale_count_mass_radius,scale_count_box_count,"
              "dense_radius_nm,dense_rg_nm,mean_local_valency,"
              "protein_rna_neighbor_fraction,valid,invalid_reason\n";
```

Write one row per processed frame. Write non-finite numeric values as `nan` using:

```cpp
auto write_number = [](std::ostream &os, double value) {
  if (std::isfinite(value))
    os << value;
  else
    os << "nan";
};
```

- [ ] **Step 2: Export scales**

Write `llps_paper_fractal_scales.csv` with:

```text
frame,time_ns,method,scale_nm,value,log_scale,log_value
```

This file is the audit trail for reviewer-facing plots. It contains all scale rows, including frames later excluded from summary means.

- [ ] **Step 3: Export summary**

Write `llps_paper_fractal_summary.dat` with exact keys:

```text
enabled: yes
frames_total: <N>
frames_valid: <Nvalid>
frames_invalid: <Ninvalid>
mean_df_mass_radius: <value>
sem_df_mass_radius: <value>
mean_df_box_count: <value>
sem_df_box_count: <value>
mean_kf_mass_radius: <value>
sem_kf_mass_radius: <value>
mean_dense_radius_nm: <value>
sem_dense_radius_nm: <value>
mean_dense_rg_nm: <value>
sem_dense_rg_nm: <value>
mean_fit_r2_mass_radius: <value>
mean_fit_r2_box_count: <value>
interpretation: morphology_descriptor_not_phase_boundary
```

Use the existing `block_reduce` pattern for SEM, with `config_.requested_blocks`.

- [ ] **Step 4: Export README**

Write `llps_paper_fractal_readme.md` with this content:

```markdown
# TITAN LLPS Paper Fractal Outputs

These files report dense-phase fractal morphology descriptors from bead/residue occupancy in the LLPS DBSCAN dense phase.

Primary citation files:

- `llps_paper_fractal_timeseries.csv`: per-frame fractal descriptors.
- `llps_paper_fractal_scales.csv`: scale-by-scale audit data for log-log fits.
- `llps_paper_fractal_summary.dat`: block-averaged summary statistics.

The descriptors quantify dense-phase morphology. They are not standalone evidence for LLPS, binodal boundaries, or spinodal behavior.
```

- [ ] **Step 5: Teach the plot aggregator to read the summary**

In `plot_llps_paper.py`, add a tolerant loader:

```python
def read_key_value_dat(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    if not path.exists():
        return out
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if ":" not in line or line.lstrip().startswith("#"):
            continue
        key, value = line.split(":", 1)
        out[key.strip()] = value.strip()
    return out
```

Add `llps_paper_fractal_summary.dat` to the replicate loader and include these fields in the generated quality/report markdown when present:

```text
mean_df_mass_radius
sem_df_mass_radius
mean_df_box_count
sem_df_box_count
interpretation
```

- [ ] **Step 6: Run the smoke**

Run:

```powershell
.\test\llps_paper_fractal_smoke.ps1
```

Expected:

```text
PASS llps_paper_fractal_smoke
```

- [ ] **Step 7: Commit Phase 4**

Run:

```powershell
git add src\titan_llps_paper.cpp plot_llps_paper.py
git commit -m "feat: export citeable LLPS paper fractal summaries"
```

---

## Phase 5: Documentation And Roadmap Sync

### Task 5: Update User-Facing Documentation

**Files:**
- Modify: `src/help.cpp`
- Modify: `manual.txt`
- Modify: `TITAN_STANDARD_MD_ANALYSIS_ROADMAP.md`

- [ ] **Step 1: Document the option in help**

In the LLPS option list in `src/help.cpp`, add:

```text
  -llps-paper-fractal=<auto|true|false>
                        Dense-phase bead/residue fractal descriptors for
                        citeable LLPS paper outputs. Default: auto
```

- [ ] **Step 2: Document the outputs in help**

In the `with -llps-paper:` output list in `src/help.cpp`, add:

```text
      llps_paper_fractal_timeseries.csv
      llps_paper_fractal_scales.csv
      llps_paper_fractal_summary.dat
      llps_paper_fractal_readme.md
```

- [ ] **Step 3: Mirror the same docs in `manual.txt`**

Apply the same option and output text to `manual.txt` so help/manual stay synchronized.

- [ ] **Step 4: Update the roadmap ledger**

In `TITAN_STANDARD_MD_ANALYSIS_ROADMAP.md`, add a planned/shipped entry under the LLPS analysis section:

```markdown
### LLPS Paper Fractal Dense-Phase Morphology

Status: planned for implementation.

Scope:
- Dense-phase bead/residue occupancy descriptors for Martini protein/RNA LLPS.
- Citeable `llps_paper_fractal_*` output files under `-llps-paper`.
- Guardrails that prevent small macro-node systems from over-reporting fractal dimension.

Non-scope:
- FracVAL structure generation.
- Claiming LLPS phase boundaries from fractal descriptors alone.
```

When the implementation ships, change `Status: planned for implementation.` to `Status: shipped.` in the same commit as the passing smoke evidence.

- [ ] **Step 5: Commit Phase 5**

Run:

```powershell
git add src\help.cpp manual.txt TITAN_STANDARD_MD_ANALYSIS_ROADMAP.md
git commit -m "docs: describe LLPS paper fractal outputs"
```

---

## Phase 6: Verification Gate

### Task 6: Build, Smoke, And Cleanup

**Files:**
- No source edits unless a verification failure identifies a specific defect.

- [ ] **Step 1: Run focused smokes**

Run:

```powershell
.\test\llps_paper_fractal_args_smoke.ps1
.\test\llps_paper_fractal_smoke.ps1
.\test\llps_paper_fractal_integration_smoke.ps1
```

Expected:

```text
PASS llps_paper_fractal_args_smoke
PASS llps_paper_fractal_smoke
PASS llps_paper_fractal_integration_smoke
```

- [ ] **Step 2: Run existing LLPS smokes**

Run:

```powershell
g++ -std=c++17 -O2 -DTITAN_LLPS_TESTING -I src test\llps_neighbor_grid_smoke.cpp src\titan_llps.cpp src\titan_llps_contacts.cpp src\titan_ui.cpp -o test\tmp\llps_neighbor_grid_smoke.exe
.\test\tmp\llps_neighbor_grid_smoke.exe
.\test\llps_dense_msd_smoke.ps1
```

Expected:

```text
llps_neighbor_grid_smoke PASS
llps_dense_msd_smoke PASS
```

If `test\llps_dense_msd_smoke.ps1` does not exist in this checkout, compile and run `test\llps_dense_msd_smoke.cpp` using the same command shape as `test\llps_paper_fractal_smoke.ps1`, adding `src\titan_llps_cluster_dynamics.cpp` and any directly required split LLPS cluster source files reported by the linker.

- [ ] **Step 3: Run native build**

Run:

```powershell
mingw32-make PROFILE=native TARGET=titan.exe all -j1
```

Expected: build completes and writes `titan.exe`.

- [ ] **Step 4: Run whitespace cleanup gate**

Run:

```powershell
git diff --check
```

Expected: no output.

- [ ] **Step 5: Inspect final diff**

Run:

```powershell
git diff --stat
git diff -- src\arguments.h src\parse_arguments.cpp src\info.h src\titan_llps_cluster.h src\titan_llps_paper.cpp src\impl_engine.cpp src\help.cpp manual.txt TITAN_STANDARD_MD_ANALYSIS_ROADMAP.md
```

Expected: only LLPS paper fractal, docs, tests, and aggregator changes are present.

- [ ] **Step 6: Commit final verification notes if roadmap status changed**

Run:

```powershell
git add TITAN_STANDARD_MD_ANALYSIS_ROADMAP.md
git commit -m "docs: mark LLPS fractal dense-phase outputs shipped"
```

Only run this commit if Phase 5 was committed before implementation and the roadmap status changed after successful verification.

---

## Self-Review

Spec coverage:

- Bead/residue dense phase is the primary analysis target.
- Output filenames include `paper` and are citeable.
- The plan uses existing `-llps-paper` rather than adding a parallel paper bundle.
- The plan includes parser, help, manual, roadmap, tests, native build, and `git diff --check`.

Placeholder scan:

- No unresolved placeholder tokens are present.
- Every new file has an exact path.
- Every output file has an exact name.

Type consistency:

- `LlpsArgs::llps_paper_fractal` is parsed in `src/parse_arguments.cpp` and consumed in `src/impl_engine.cpp`.
- `LLPSPaperConfig::fractal_enabled` is passed to `LLPSPaperFractalAccumulator`.
- `info` filename helpers match the documented output filenames.
