# TITAN Energy Decomposition Refactor Roadmap

## Goal

Replace the legacy energy-map CLI surface with one explicit energy
decomposition contract while preserving the existing force-map family.
The implementation now includes the new CLI, MM elec+vdW provider/engine
path, the public MM output contract, and endpoint adapters for MM/GBSA and
MM/PBSA.

## CLI Contract

Use `-option=value` syntax in `arg.in` files.

| Option | Values | Default | Notes |
| --- | --- | --- | --- |
| `-edecomp=true` | `true`, `false` | `false` | Enables unified energy decomposition. |
| `-edecomp-source=mm` | `mm`, `gbsa`, `pbsa`, `all` | `mm` | `gbsa` requires `-gbsa`; `pbsa` requires `-pbsa=true`. |
| `-edecomp-level=residue` | `residue`, `atom`, `both` | `residue` | Selects output granularity. |
| `-edecomp-edge-top=200` | positive integer, `all` | `200` | Limits ranked edge output. |
| `-edecomp-rank-by=abs-mean` | `abs-mean`, `mean-abs` | `abs-mean` | Selects edge ranking metric. |
| `-edecomp-frame-output=off` | `off`, `all` | `off` | Controls per-frame edge output. |
| `-edecomp-perf=summary` | `off`, `summary`, `detail` | `summary` | Controls performance diagnostics. |

`-edecomp-source=all` does not auto-enable GBSA or PBSA. It only requests
the explicitly enabled providers.

## Migration Contract

The parser rejects these legacy energy-only options with a migration error:

- `-emap`
- `-emap-level`
- `-energy-res`
- `-new-emap`
- `-useForceOptionforEnergy`

Use `-edecomp=true` and the `-edecomp-*` options instead. The `-fmap`,
`-fmap-level`, and `-force-res` family remains unchanged.

## Legacy Removal Status

Implemented. Legacy fields such as `energy_map`, `emap_type`, `emap_level`,
`print_energy_res`, `flag_new_emap`, and `flag_energy_same_force` have been
removed from `arguments` and their engine call sites. The old `new_emap.*`,
`stat_ener_res_res.*`, and residue-energy statistic writer path have been
removed from the build. `-useForceOptionforEnergy` now fails with guidance to
use `-edecomp=true` for energy decomposition and `-fmap` for force maps.

## Implementation Status

| Area | Status | Notes |
| --- | --- | --- |
| CLI contract | Implemented | `EDecompArgs` parses `-edecomp*`; legacy `-emap`, `-emap-level`, `-energy-res`, `-new-emap`, and `-useForceOptionforEnergy` fail with migration guidance. |
| Arg-file equals syntax | Implemented for EDecomp smoke path | `-f=`, `-s=`, `-n=`, `-g=`, and `-ff=` are accepted so new `arg.in` smoke files can stay on `-option=value`. |
| MM pair provider | Implemented | `EnergyDecompositionProvider` exposes cache-indexed and batch fallback pair energy access with scalar/AVX2-compatible callers. |
| Shared MM engine | Implemented | Uses `TitanExecutor`, stable residue/atom axes, dense residue matrices, sparse atom edges, rankings, deterministic trajectory aggregation. |
| Streaming MM wiring | Implemented | `impl_engine` runs MM decomposition from `CachedData` when available and provider fallback otherwise. |
| MM outputs | Implemented | Emits `edecomp_mm_residue_matrix_{elec,vdw,total}.csv`, residue/atom edges, residue/atom rankings, and `edecomp_summary.dat`. |
| GBSA adapter | Implemented | Existing GBSA solver outputs remain; unified `edecomp_mmgbsa_*` entity maps/rankings are wired for `-edecomp-source=gbsa|all`. |
| PBSA adapter | Implemented | Existing PBSA solver/decomposition outputs remain; unified `edecomp_mmpbsa_*` endpoint maps/rankings are wired for `-edecomp-source=pbsa|all`. |
| Legacy code removal | Implemented | Unified `-edecomp` is the only energy decomposition/reporting framework; `-fmap` remains independent. |

## Integration Phases

### Phase 1: CLI and Documentation

- Add `EDecompArgs` defaults and parser validation.
- Reject legacy energy-only flags with a clear migration message.
- Keep the force-map family untouched.
- Add focused parser/CLI smoke coverage.

### Phase 2: Provider Integration

- Route MM pair decomposition through the unified engine adapter.
- Route GBSA decomposition only when `-gbsa` is explicitly enabled.
- Route PBSA decomposition only when `-pbsa=true` is explicitly enabled.
- For `source=all`, aggregate only providers explicitly enabled by the user.

### Phase 3: Output Contract

- Emit stable residue, atom, or combined edge outputs.
- Apply `edge-top`, `rank-by`, and frame-output policy consistently.
- Add summary and detail performance diagnostics without changing numerical
  results.

### Phase 4: Legacy Removal

- Implemented: transitional energy-only fields and call sites are gone.
- Implemented: `-useForceOptionforEnergy` rejects with a clear migration error.
- Implemented: `fmap` behavior remains stable and independent.

## Roadmap Completion

Complete. Remaining work is future tuning only: performance probes on larger
production trajectories, optional uncertainty diagnostics, and any new source
adapters that should plug into the unified endpoint engine.

## Verification

- Run `test/edecomp_args_smoke.ps1`.
- Run `test/edecomp_smoke.ps1`.
- Run `test/titan_energy_decomposition_provider_smoke.cpp`.
- Run `test/titan_energy_decomp_engine_smoke.cpp`.
- Run existing GBSA regression smoke after GBSA adapter changes.
- Compile the focused parser harness generated by that script.
- Run `git diff --check`.
- Run engine-owned integration smokes when the provider slice is merged.
