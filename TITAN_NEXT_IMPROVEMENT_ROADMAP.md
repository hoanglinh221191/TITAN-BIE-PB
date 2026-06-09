# TITAN Next Improvement Roadmap Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Turn the current TITAN feature set into a more integrated, measurable, and scientifically defensible MD analysis platform.

**Architecture:** Keep existing public modules stable, then add shared infrastructure that lets observables, kinetics, validation, reporting, and performance telemetry talk to each other. Prioritize opt-in paths so default `-tica`, `-msm`, `-tpt`, PBE/grid, GBSA/PBSA, and force-map behavior remain unchanged unless a new option is explicitly selected.

**Tech Stack:** C++23, TitanExecutor, aligned SoA buffers, AVX2/FMA with scalar fallback, PowerShell smoke tests, native/portable build profiles, existing `manual.txt`, `src/help.cpp`, `arguments.*`, `parse_arguments.cpp`, and roadmap files.

---

## Summary

TITAN already has strong individual modules: CMap, QNative, Dihedral, MSD,
Membrane, S2, CAS, EDecomp, validation, MBAR overlays, graph analytics, and
optimized PBE/TABI diagnostics. The next improvements should focus less on
adding isolated analyzers and more on:

1. shared feature/CV export for TICA/MSM/TPT,
2. unified performance and output observability,
3. uncertainty/convergence diagnostics across standard observables,
4. domain extensions where TITAN already has a scientific foothold,
5. benchmark-gated solver work only where evidence supports it.

## Priority Ledger

| Rank | Area | Recommendation | Why |
| ---: | --- | --- | --- |
| 1 | Feature + science | Shared `-feature-export` and mixed-source TICA input | Highest leverage: turns `-cmap`, `-qnative`, `-dihedral`, `-s2`, `-membrane`, and LJ features into real CV sources for MSM/TPT. |
| 2 | Performance | Unified per-module performance telemetry | Makes future optimization measurable instead of intuition-driven. |
| 3 | Stability + UX | Output manifest and `-dry-run` argument/schema validation | Reduces confusion from many outputs and catches bad runs before expensive analysis. |
| 4 | Science | Cross-module uncertainty and convergence framework | Gives users confidence in contact/Q/MSD/S2/membrane/EDecomp rankings. |
| 5 | Domain features | Membrane protein-lipid extension and CAS/EDecomp cross-ranking | Builds on implemented Phase 5 and unified energy decomposition. |
| 6 | Solver performance | PBE/TABI only through benchmark-driven experimental gates | Current evidence favors grid defaults; TABI should remain opt-in unless new medium/large benchmarks justify promotion. |

## Non-Goals

- Do not change default `-tica`, `-msm`, `-tpt`, `-pbe=grid`, `-edecomp`, or
  `-fmap` semantics.
- Do not add another isolated observable before wiring existing observables
  into the feature/CV layer.
- Do not promote TABI/FMM-style paths based on small-fixture speedups alone.
- Do not silently pool independent trajectories into one kinetic trajectory.

---

## Phase 1: Shared Feature / CV Export

**Category:** Feature + science

**Status:** Phase 5 v1 implemented. Protein-lipid contacts, leaflet enrichment,
curvature grid proxy, sidecar outputs, multi-input aggregate summaries, and
`test/membrane_contacts_smoke.ps1` are in place. Contact code is integrated
inside `TitanMembraneAnalyzer` to reuse its lipid topology/leaflet cache.

**Status:** Recommended next phase.

**Goal:** Create a shared feature matrix layer that can feed TICA/MSM/TPT from
new observables without copying trajectory data again.

**Files:**
- Create: `src/titan_observable_features.h`
- Create: `src/titan_observable_features.cpp`
- Modify: `src/arguments.h`
- Modify: `src/arguments.cpp`
- Modify: `src/parse_arguments.cpp`
- Modify: `src/help.cpp`
- Modify: `src/info.h`
- Modify: `src/impl_engine.cpp`
- Modify: `manual.txt`
- Modify: `TITAN_STANDARD_MD_ANALYSIS_ROADMAP.md`
- Test: `test/titan_observable_features_smoke.cpp`
- Test: `test/feature_export_smoke.ps1`

**Public Interface:**

```text
-feature-export=true|false
-feature-source=<auto|dccm|lj|dihedral|qnative|cmap|mixed>
```

**Output Contract:**

```text
features_by_frame[_N].csv
features_metadata[_N].csv
feature_export_summary[_N].dat
```

For TICA opt-in:

```text
-tica=true -feature-source=dihedral
-tica=true -feature-source=mixed
```

Default protein-like `-tica` continues to use the current network/DCCM path
unless the user explicitly sets a feature source.

**Implementation Tasks:**

- [x] Add `FeatureExportArgs` with defaults: disabled, source `auto`.
- [x] Add parser validation for allowed sources and auto-enable explicit
  `dihedral`, `qnative`, and `cmap` producers.
- [x] Implement feature metadata: source, name, unit, periodic flag, residue
  label, and role.
- [x] Register QNative scalar Q.
- [x] Register Dihedral sin/cos columns for enabled phi/psi/omega/chi1/chi2.
- [x] Register CMap scalar `contact_count` and normalized count.
- [x] Keep S2 and membrane as future feature producers until their
  frame-feature semantics are designed.
- [x] Add TICA input switch for explicit `-feature-source` while preserving
  default DCCM/LJ behavior in `auto`.
- [x] Add feature export outputs: `features_by_frame[_N].csv`,
  `features_metadata[_N].csv`, and `feature_export_summary[_N].dat`.
- [ ] Add broader long-run smoke proving `-feature-source=mixed` with full
  TICA/MSM/TPT on a production-like fixture.

**Scientific Notes:**

- Circular torsions must use `sin(theta),cos(theta)` rather than raw degrees.
- QNative can be a scalar CV and native-contact bits can be sparse feature
  channels.
- Independent trajectory chunks must be indexed and not concatenated into one
  kinetic trajectory unless the existing MSM workflow explicitly supports that
  input mode.

---

## Phase 2: Unified Performance Observatory

**Category:** Performance

**Status:** Phase 2A implemented. The shared observer, CLI, outputs, standard-MD timers, and CAS/EDecomp/GBSA/PBSA summary bridges are in place. Exact peak-memory accounting remains future benchmark-driven tuning.

**Goal:** Add one cross-module performance report so future optimization work
has comparable evidence.

**Files:**
- Create: `src/titan_performance_observer.h`
- Create: `src/titan_performance_observer.cpp`
- Modify: `src/impl_engine.cpp`
- Modify: `src/info.h`
- Modify: analyzers that already expose perf counters: CAS, EDecomp, GBSA, PBSA, MSD, S2, Membrane, CMap, QNative, Dihedral
- Modify: `manual.txt`
- Modify: `src/help.cpp`
- Test: `test/titan_performance_observer_smoke.cpp`
- Test: `test/performance_observer_smoke.ps1`

**Public Interface:**

```text
-perf-report=<off|summary|detail>
```

Default: `off` to avoid output churn in existing runs.

**Outputs:**

```text
titan_performance_summary.dat
titan_performance_frames.csv      # only with detail
titan_performance_modules.csv
```

**Metrics:**

- frames seen and frames analyzed,
- module wall time,
- per-frame average time,
- export time,
- peak scratch capacity,
- estimated memory footprint,
- worker count and grain policy,
- warning count and skipped-frame count.

**Implementation Tasks:**

- [x] Add a lightweight `PerformanceObserver` with `begin(module)`, `end(module)`, and frame-scoped timers.
- [x] Route existing CAS timing into the observer while keeping `cas_performance.dat` unchanged.
- [x] Route EDecomp timing into the observer without changing EDecomp files.
- [x] Bridge GBSA/PBSA frame-level timing at summary granularity without refactoring solver internals.
- [x] Add CMap/QNative/Dihedral/MSD/S2/Membrane module timers around compute and export sections.
- [x] Add deterministic CSV ordering by module name and frame index.
- [x] Add smoke that verifies output exists and contains all active modules.
- [x] Add smoke that verifies turning the report off creates no shared perf files.
- [ ] Add exact peak-memory counters for every module once scratch ownership is standardized.

**Performance Notes:**

- Do not place high-frequency timers inside innermost AVX loops.
- Use one timer per frame/module section; keep detail mode opt-in if overhead is measurable.

---

## Phase 3: Output Manifest And Dry-Run Validation

**Category:** Stability + feature UX

**Status:** Phase 3 v1 implemented. Primary output manifest, dry-run early exit,
run-config summary, public schema JSON, and equals-syntax smokes are in place.
Future tuning can expand the schema to every hidden/debug option and add more
fine-grained trajectory-index prediction for complex batch layouts.

**Goal:** Make TITAN runs inspectable before and after execution.

**Files:**
- Create: `src/titan_output_manifest.h`
- Create: `src/titan_output_manifest.cpp`
- Create: `src/titan_config_schema.h`
- Create: `src/titan_config_schema.cpp`
- Modify: `src/arguments.h`
- Modify: `src/arguments.cpp`
- Modify: `src/parse_arguments.cpp`
- Modify: `src/info.h`
- Modify: `src/help.cpp`
- Modify: `manual.txt`
- Test: `test/output_manifest_smoke.ps1`
- Test: `test/dry_run_validation_smoke.ps1`

**Public Interface:**

```text
-dry-run=true|false
-output-manifest=true|false
-schema-export=<off|json>
```

**Outputs:**

```text
titan_outputs_manifest.csv
titan_run_config_summary.dat
titan_schema.json                 # only when requested
```

**Implementation Tasks:**

- [x] Emit manifest rows for primary expected outputs: module, filename,
  trajectory index/pattern, units, averaging policy, and semantic note.
- [x] Add `-dry-run=true` path that parses inputs, validates options, resolves
  modules, writes expected outputs/schema/config sidecars, then exits before
  coordinate loading and trajectory processing.
- [x] Add schema export for public option names, types, defaults, accepted
  values, and incompatible-option notes.
- [x] Keep validation for common expensive-run mistakes in parser/post-process:
  invalid membrane group, `-cas=true -nocal=true`,
  missing `-gbsa` for `-edecomp-source=gbsa`, and missing `-pbsa` for
  `-edecomp-source=pbsa`.
- [x] Add smoke with `arg.in` equals syntax.
- [x] Add smoke proving normal run writes manifest and expected module outputs.
- [ ] Expand manifest/schema coverage to every low-level debug option and rare
  sidecar once the public output registry is centralized.

**UX Notes:**

- Prefer explicit migration messages over silent aliases.
- Manifest rows should include scientific warnings such as "MM/GBSA maps are endpoint entity decompositions, not pair matrices."

---

## Phase 4: Cross-Module Uncertainty And Convergence

**Category:** Science

**Goal:** Add block/bootstrap diagnostics for standard observables and rankings,
without changing default point estimates.

**Files:**
- Create: `src/titan_uncertainty.h`
- Create: `src/titan_uncertainty.cpp`
- Modify: CMap, QNative, Dihedral, MSD, S2, Membrane, EDecomp ranking exporters
- Modify: `src/arguments.h`
- Modify: `src/arguments.cpp`
- Modify: `src/help.cpp`
- Modify: `manual.txt`
- Test: `test/titan_uncertainty_smoke.cpp`
- Test: `test/uncertainty_standard_observables_smoke.ps1`

**Public Interface:**

```text
-uncertainty=<off|block|bootstrap>
-uncertainty-blocks=<auto|N>
-uncertainty-bootstrap=<N>
-uncertainty-seed=<N>
```

**Outputs:**

```text
uncertainty_summary.dat
*_uncertainty.csv                 # module-specific sidecars
```

**Implementation Tasks:**

- [x] Add common block partition helper with deterministic frame assignment.
- [x] Add mean, SEM, confidence interval, and effective sample count helper.
- [x] Add CMap/QNative uncertainty sidecars for probabilities and Q means.
- [x] Add Dihedral uncertainty for rotamer occupancy.
- [x] Add MSD V1 sidecar for MSD lag-series diagnostics; slope/D/alpha
      refit uncertainty remains a future in-analyzer refinement.
- [x] Add S2 uncertainty sidecar for plateau and clipped S2 entity rows.
- [x] Add Membrane uncertainty sidecar for APL, thickness, and order parameter.
- [x] Add EDecomp uncertainty sidecar for rankings with both `abs_mean` and `mean_abs`.
- [x] Keep MSM/TPT validation bootstrap unchanged; this phase extends observable-level diagnostics only.

**Implementation Status:** Phase 4 V2 implemented as opt-in sidecar diagnostics
from finalized module outputs. Existing point-estimate schemas are unchanged.
MSD now emits fit-aware rows for slope, diffusion constants, alpha, and
prefactor; S2 and Membrane sidecars include richer entity/observable labels and
Membrane interdigitation diagnostics. Ranking-only outputs remain marked as
point-estimate/entity-distribution diagnostics instead of fake time uncertainty.
Future work, if needed, is in-analyzer per-frame block state for CMap,
QNative, and EDecomp.

**Scientific Notes:**

- Independent replicas should aggregate by trajectory-equal means and report
  pooled-frame diagnostics separately.
- `SEM=nan` is acceptable for one valid replicate.

---

## Phase 5: Membrane Protein-Lipid Extension

**Category:** Feature + science

**Goal:** Extend Phase 5 membrane analysis from lipid-only observables to
protein-lipid coupling and membrane remodeling.

**Files:**
- Implemented in: `src/titan_membrane.h`
- Implemented in: `src/titan_membrane.cpp`
- Modify: `src/titan_membrane.h`
- Modify: `src/titan_membrane.cpp`
- Modify: `src/arguments.h`
- Modify: `src/arguments.cpp`
- Modify: `src/help.cpp`
- Modify: `manual.txt`
- Test: `test/membrane_contacts_smoke.ps1`

**Public Interface:**

```text
-membrane-protein-contacts=true|false
-membrane-contact-cut=<A>
-membrane-enrichment=<true|false>
-membrane-curvature=<off|grid>
```

**Outputs:**

```text
membrane_protein_lipid_contacts.csv
membrane_lipid_enrichment.csv
membrane_curvature_map.csv
membrane_proximity_coupling.csv
membrane_coupling_summary.dat
membrane_lipid_enrichment_average.csv
membrane_curvature_map_average.csv
membrane_proximity_coupling_average.csv
membrane_coupling_summary_average.dat
```

**Implementation Tasks:**

- [x] Implement mol1 protein residue side-chain COM to lipid head/tail bead
  contacts with MIC and `TitanExecutor`.
- [x] Export lipid resname enrichment/depletion around protein by leaflet.
- [x] Compute local curvature proxy from headgroup surface grid.
- [x] Correlate local thickness/APL/order parameter with protein proximity.
- [x] Add smoke with synthetic MARTINI protein-bilayer fixture.

**Scientific Notes:**

- Keep MARTINI-first mapping.
- Protein-proximity coupling uses per-frame near/far lipid masks from structural
  protein-lipid contacts and reports APL, thickness, order, curvature, contact
  probability, and leaflet-composition deltas.
- Atomistic lipid mapping remains a future slice unless a fixture and mapping table are available.

---

## Phase 6: CAS + EDecomp Scientific Integration

**Category:** Science + feature UX

**Status:** Implemented. `-hotspot-report` now exports a residue-level
CAS/EDecomp consensus ranking, summary and B-factor visualization when a CAS
PDB template is available.

**Goal:** Let users compare alanine-like CAS hotspots with MM/MMGBSA/MMPBSA
energy decomposition rankings in one consistent hotspot report.

**Files:**
- Create: `src/titan_hotspot_report.h`
- Create: `src/titan_hotspot_report.cpp`
- Modify: parser/validation, `src/impl_engine.cpp`, `src/info.h`
- Modify: `src/help.cpp`
- Modify: `manual.txt`
- Test: `test/hotspot_report_smoke.ps1`

**Public Interface:**

```text
-hotspot-report=true|false
-hotspot-source=<cas|edecomp|both>
-hotspot-rank-by=<cas-ddg|edecomp-abs-mean|combined>
```

**Outputs:**

```text
hotspot_residue_ranking.csv
hotspot_summary.dat
hotspot_bfactor.pdb
```

**Implementation Tasks:**

- [x] Join CAS rows and EDecomp ranking rows by stable molecule role + chain + residue key.
- [x] Add combined ranking score that reports components rather than hiding them.
- [x] Keep CAS method note: alanine-like side-chain neutralization, no relaxation/repacking.
- [x] Add summary warning when CAS and EDecomp sources cover different residue sets.
- [x] Add smoke with receptor/ligand duplicate residue numbers across chains.

**Scientific Notes:**

- Do not call CAS a full alanine mutation unless structural relaxation is implemented.
- EDecomp solvation terms remain endpoint entity maps, not pair matrices.

---

## Phase 7: Benchmark-Gated PBE/TABI Direction

**Category:** Performance + science

**Status:** Implemented. `test/pbe_tabi_phase7_gate.ps1` now wraps the
existing PBE/TABI perf matrix into a repeatable decision gate and exports
normalized matrix, decision and summary files. The current production decision
remains conservative: keep `grid` as the default and keep TABI opt-in /
experimental unless a full medium/large run emits `promote_candidate`. The
gate uses research PQR fixtures from `example/1/rec_027.pqr` and
`example/4/rec_0.pqr`. Latest local run with 60-second per-case timeout kept
TABI experimental: grid solved both fixtures in about `0.6 s` and `1.5 s`,
while TABI direct and treecode timed out on both.

**Goal:** Decide whether further TABI investment is justified by medium/large
fixtures before adding more solver machinery.

**Files:**
- Modify: `TITAN_PBE_TABI_ROADMAP.md`
- Add: `test/pbe_tabi_phase7_gate.ps1`
- Test: Phase 7 gate quick smoke plus existing PBE/TABI matrix and compare scripts

**Outputs:**

```text
test/tmp/pbe_tabi_phase7_gate/pbe_tabi_phase7_matrix.csv
test/tmp/pbe_tabi_phase7_gate/pbe_tabi_phase7_decision.csv
test/tmp/pbe_tabi_phase7_gate/pbe_tabi_phase7_summary.dat
```

**Decision Gate:**

Continue TABI implementation only if medium/large evidence shows one of:

- a real memory advantage that matters for target workloads,
- a runtime path competitive with grid at accepted accuracy,
- a scientific accuracy path that grid cannot provide for the target use case.

Otherwise:

- keep grid as production default,
- keep TABI hidden/experimental or opt-in,
- focus performance work on shared runtime, output, and observable pipelines.

**Implementation Tasks:**

- [x] Add a repeatable gate that compares grid, TABI direct and TABI treecode.
- [x] Record time, RAM, accuracy status, GMRES iterations and failure modes in a normalized CSV.
- [x] Export an explicit decision row: `promote_candidate`,
  `keep_experimental`, `pause_new_solver_work`, or quick-smoke
  non-authoritative status.
- [x] Keep production default unchanged: `grid` remains default and TABI
  remains hidden/experimental unless the gate passes.
- [x] Do not add FMM, new scheduler layers, or another preconditioner until
  the decision gate passes.

---

## Verification Matrix

Every implementation slice should run the narrow smoke for its module plus:

```powershell
mingw32-make PROFILE=native TARGET=titan.exe all -j1
git diff --check
```

Recommended regression set after major integration work:

```powershell
test\cmap_smoke.ps1
test\qnative_smoke.ps1
test\dihedral_smoke.ps1
test\msd_smoke.ps1
test\membrane_smoke.ps1
test\s2_smoke.ps1
test\edecomp_smoke.ps1
test\cas_smoke.ps1
test\gbsa_perf_entropy_smoke.ps1
test\pbsa_smoke.ps1
```

For TICA/MSM/TPT-facing changes, also run:

```powershell
test\feature_export_smoke.ps1
```

## Roadmap Completion Criteria

This roadmap is considered complete when:

- `-feature-export` can feed explicit mixed features into TICA without changing
  default TICA behavior.
- Shared performance reports identify per-module compute/export cost.
- TITAN can emit an output manifest and dry-run config summary.
- Standard observables have uncertainty sidecars where scientifically useful.
- Membrane and hotspot extensions are available as opt-in features.
- PBE/TABI has an evidence-backed keep/pause/promote decision.
