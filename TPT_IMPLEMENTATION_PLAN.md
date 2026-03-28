# TITAN TPT Implementation Plan

This document is the detailed implementation roadmap for adding
Transition Path Theory (TPT) to TITAN.

Primary rollout order:

1. LJ REMD first
2. protein REMD after LJ is stable
3. kinetic-TPT only after ensemble-TPT is working and validated

================================================================================
1. Goals
================================================================================

Add a TPT module that can:

- compute forward and backward committors
- compute reactive flux and net flux on MSM state networks
- identify dominant reactive pathways and bottleneck edges
- report basin populations and reactive populations
- export files that can be plotted directly for LJ and protein workflows

The first scientifically safe target is:

- ensemble-TPT on MSMs built from REMD data

The first deployment target is:

- LJ T-REMD in NVT

Then extend to:

- LJ T-REMD in NPT
- protein REMD in NPT

Absolute kinetic claims must be kept separate from mechanism/flux claims.

================================================================================
2. TPT modes
================================================================================

Two distinct operating modes are planned.

2.1 Ensemble-TPT

Use case:

- REMD trajectories
- DHAM/TRAM-like reweighting targets
- mechanism-focused analysis

Interpretation:

- strong for pathway topology, bottlenecks, committor ordering, and reactive
  flux structure
- not the first choice for absolute kinetics when replica exchange breaks
  continuous physical dynamics

2.2 Kinetic-TPT

Use case:

- plain MD trajectories at a fixed target ensemble
- or demultiplexed trajectories that are continuous in time at one target
  temperature/pressure

Interpretation:

- can support stronger claims about rates and MFPT
- requires better validation of MSM Markovianity and sufficient back-and-forth
  transitions between basins

================================================================================
3. Reusable TITAN infrastructure already available
================================================================================

The current codebase already provides most of the hard prerequisites.

Existing components to reuse:

- `TitanMSM`
  - microstate TPM
  - count matrix
  - stationary distribution
  - discrete trajectory
  - PCCA macrostate assignment
  - macrostate trajectory export
- `TitanDHAM`
  - reweighted stationary distribution and TPM at target conditions
  - frame weights for downstream reweighted analyses
- `TitanKinetics`
  - TICA feature pipeline
- LJ feature backend
  - cluster, density, CN, q4/q6, q4_bar/q6_bar, solid-like fractions
- protein network/TICA/MSM pipeline
  - DCCM-derived TICA features

Main files expected to be touched:

- `arguments.h`
- `parse_arguments.cpp`
- `arguments.cpp`
- `info.h`
- `impl_engine.cpp`
- `engine_postproc.cpp`
- new module:
  - `titan_tpt.h`
  - `titan_tpt.cpp`

================================================================================
4. Proposed CLI and outputs
================================================================================

4.1 New argument block

Add a new `TptArgs` struct to `arguments.h`.

Phase-1 flags:

- `-tpt`
- `-tpt-mode ensemble|kinetic`
- `-tpt-source <spec>`
- `-tpt-sink <spec>`
- `-tpt-top-paths <N>`
- `-tpt-flux-cut <val>`
- `-tpt-use-dham yes|no`
- `-tpt-macro yes|no`

Optional later flags:

- `-tpt-coarse micro|macro`
- `-tpt-path-fraction <val>`
- `-tpt-basin-file <file>`
- `-tpt-state-weights msm|dham`
- `-tpt-allow-nonreversible yes|no`

4.2 Basin specification

Phase 0 should support only the simplest and most robust basin definitions:

- microstate list
  - example: `-tpt-source micro:0,1,2`
- macrostate list
  - example: `-tpt-source macro:0`
- TICA rectangle
  - example: `-tpt-source tica_rect:-3.0,-1.0,0.5,2.0`

LJ-specific basin definitions can be added in phase 1:

- threshold on `largest_cluster_fraction`
- threshold on `solid_like_fraction`
- threshold on `q6_bar`

Protein-specific basin definitions can be added later:

- macrostate IDs
- TICA/PCA rectangles
- representative-basin metadata from PCA/TICA FES

4.3 Output files

New `info.h` getters should produce:

- `tpt_summary.dat`
- `tpt_committor.dat`
- `tpt_flux_edges.dat`
- `tpt_paths.dat`
- `tpt_state_labels.dat`
- `tpt_reactive_fes_2d.dat`
- `tpt_flux_map_2d.dat`

Optional later:

- `tpt_reactive_projection.dat`
- `tpt_macro_flux_edges.dat`

================================================================================
5. Phase 0: core TPT engine
================================================================================

Goal:

Build a general TPT solver that works on a Markov model regardless of whether
the states came from LJ or protein workflows.

5.1 Implementation tasks

Add `TitanTPT` in `titan_tpt.h/.cpp`.

Core inputs:

- transition matrix `P`
- stationary distribution `pi`
- source basin `A`
- sink basin `B`

Core computations:

- validate `A` and `B`
  - non-empty
  - disjoint
  - inside state index range
- compute forward committor `q+`
  - boundary conditions:
    - `q+(i)=0` for `i in A`
    - `q+(i)=1` for `i in B`
- compute backward committor `q-`
  - for reversible chain this can be simplified; keep code general enough to
    support future nonreversible handling if needed
- compute reactive flux:
  - `f_ij = pi_i * q-_i * P_ij * q+_j`
- compute net flux:
  - `f^+_ij = max(0, f_ij - f_ji)`
- compute total reaction rate and optional MFPT-like summaries
- extract dominant pathways from net-flux graph
  - start with greedy widest-path / highest-capacity path extraction
  - subtract path flux iteratively

5.2 Numerical policy

Phase 0 should assume:

- row-stochastic TPM
- non-negative probabilities
- finite stationary distribution

Phase 0 should reject:

- NaN/Inf in `P` or `pi`
- negative probabilities beyond tolerance
- disconnected `A` or `B` with zero outgoing/incoming support

5.3 Integration point

Add `run_tpt_postproc(...)` in `engine_postproc.cpp`.

Initial order:

1. TICA
2. MSM
3. DHAM
4. TPT
5. reweighted FES/PCA/UMAP if relevant

TPT should choose the transition model source as follows:

- if `-tpt-use-dham=yes` and DHAM converged:
  - use DHAM reweighted TPM / stationary distribution
- else:
  - use MSM TPM / stationary distribution

================================================================================
6. Phase 1: LJ REMD NVT ensemble-TPT
================================================================================

This is the first real deployment target.

Scientific question class:

- gas -> liquid
- liquid -> dense droplet
- liquid -> solid-like ordering

6.1 State space

Use current LJ phase-2 features:

- `largest_cluster_fraction`
- `largest_cluster_rg`
- `largest_cluster_rmax`
- `mean_cn`
- `std_cn`
- `dense_fraction`
- `mean_dense_q4/q6`
- `mean_dense_q4_bar/q6_bar`
- `solid_like_fraction`
- `fcc/hcp/bcc/ico-like fractions`
- shell enrichments

Workflow:

- LJ features -> TICA -> MSM -> TPT

6.2 Basin definition priority

Phase 1 should support:

- macrostate-based basins
  - recommended default
- LJ-threshold basins
  - `cluster_frac >= x`
  - `solid_like_frac >= y`
  - `q6_bar >= z`

Recommended first UI:

- if PCCA exists, allow:
  - `-tpt-source macro:0`
  - `-tpt-sink macro:1`

6.3 Outputs and interpretation

Expected user value:

- identify dominant condensation pathways
- identify pre-critical region in MSM state space
- tell whether ordering occurs before or after densification

Claims allowed in documentation:

- pathway structure
- bottlenecks
- committor ordering
- reactive flux topology

Claims to avoid:

- absolute physical rates unless data are demultiplexed/continuous

================================================================================
7. Phase 2: LJ REMD NPT ensemble-TPT
================================================================================

Scientific question class:

- densification vs ordering under fluctuating volume
- gas/liquid/supercritical transitions
- liquid -> solid-like ordering in a variable-density box

7.1 State space adjustment

NPT should explicitly keep:

- `rho_bulk`
- or `volume`

inside the TICA/MSM input space.

Recommended workflow:

- LJ features + `rho_bulk` -> TICA -> MSM -> DHAM(NPT) -> TPT

7.2 Model source

Preferred source for TPT in NPT:

- DHAM-reweighted TPM and stationary distribution at target `(T,P)`

Fallback:

- plain MSM TPM if DHAM does not converge

7.3 Interpretation

Phase-2 TPT should answer:

- does the system densify before it becomes ordered?
- what pathway fraction goes through dense-liquid intermediates?
- is there a distinct pre-ordered dense cluster bottleneck?

================================================================================
8. Phase 3: LJ kinetic-TPT
================================================================================

Only start this phase after phases 0-2 are working.

Data requirements:

- plain NVT/NPT MD at a fixed target ensemble
  or
- demultiplexed REMD trajectories continuous at the chosen target state point

Validation requirements:

- enough back transitions between basins
- ITS reasonably flat near chosen lag
- CK-like checks acceptable
- basin occupancy not vanishingly small

Deliverables:

- kinetic-TPT mode using MSM TPM built from continuous trajectories
- rate and MFPT reporting only in this mode

================================================================================
9. Phase 4: protein REMD NPT ensemble-TPT
================================================================================

Use current protein pathway after LJ is stable.

Workflow:

- protein DCCM/TICA features -> MSM -> optional DHAM(NPT) -> TPT

Recommended first basin definitions:

- macrostate IDs from PCCA
- TICA rectangle

Protein-specific use cases:

- folded/open vs compact/expanded transitions
- binding-competent vs binding-incompetent macrostates
- RNA-contact-rich vs RNA-contact-poor basins

Interpretation:

- strong on mechanism and bottlenecks
- weak on absolute kinetics for raw REMD trajectories

================================================================================
10. Phase 5: protein kinetic-TPT
================================================================================

Only after ensemble-TPT is validated.

Data requirements:

- conventional MD or demux trajectories continuous at one target `(T,P)`
- enough transitions between metastable basins

This phase can later serve:

- folding/unfolding
- open/closed transitions
- binding/unbinding-like coarse observables

================================================================================
11. Special note: Martini LLPS trajectories and kinetic-TPT
================================================================================

For Martini protein/RNA LLPS, kinetic-TPT is potentially useful but only when
the trajectories satisfy kinetic assumptions.

Three trajectories can be enough to attempt kinetic-TPT if all of the following
are true:

- they are plain continuous MD trajectories, not REMD
- all runs share the same Hamiltonian, composition, box size policy, and target
  `(T,P)`
- they contain multiple transitions between candidate basins, not only a single
  one-way condensation event
- the chosen basin definition is coarse enough to be revisited repeatedly
- MSM lag validation is acceptable

Three trajectories are usually not enough for strong kinetic claims if:

- each trajectory only condenses once and never redissolves
- the system remains trapped in one condensed state after nucleation
- there are too few A<->B recrossings
- basin states are extremely imbalanced

Recommended practical rule:

- if your three Martini trajectories each show repeated transitions between
  dilute, pre-condensed, and condensed basins, kinetic-TPT is worth trying
- if they mainly show one-way LLPS formation, use ensemble/mechanism-focused
  TPT first and treat rates as exploratory only

Recommended first Martini LLPS basins:

- `A`: dilute / disconnected state
- `I`: pre-condensed or percolating intermediate
- `B`: condensed/percolated state

Recommended first Martini LLPS features:

- largest-cluster fraction
- number of intermolecular contacts
- protein-protein and protein-RNA contact fractions
- dense-phase occupancy/percolation measure
- Rg or compactness of the largest condensate

================================================================================
12. CLI phase rollout
================================================================================

Phase 0:

- parse `-tpt`
- parse source/sink
- postproc only, using existing MSM outputs

Phase 1:

- add LJ threshold basin parsing
- add macrostate basin helpers

Phase 2:

- add `-tpt-use-dham`
- route to DHAM-reweighted model if available

Phase 3+:

- add `-tpt-mode kinetic`
- enforce stricter validation/warnings for kinetic mode

================================================================================
13. Validation plan
================================================================================

13.1 Unit-style validation

Create a small synthetic MSM with known A/B flux:

- 3-state and 4-state toy TPMs
- verify committor monotonicity
- verify `q+(A)=0`, `q+(B)=1`
- verify non-negative net flux

13.2 Integration validation on LJ

Use existing `rho03` examples:

- NVT LJ smoke run with TICA/MSM
- NPT LJ reweighted run with DHAM

Check:

- TPT files are written
- flux graph is non-empty
- dominant path extraction is stable under small flux cut changes

13.3 Integration validation on protein

Use the REMD NPT dataset already used for DHAM smoke tests.

Check:

- macrostate-based basin selection works
- TPT runs even if DHAM falls back

13.4 Visualization validation

Check:

- `tpt_reactive_fes_2d.dat` is written whenever TICA/MSM centroids have at
  least two dimensions
- the reactive FES maximum aligns with high-reactive-density microstates
- `tpt_flux_map_2d.dat` arrows terminate on valid TICA centroids and conserve
  the same net-flux ordering as `tpt_flux_edges.dat`

================================================================================
14. Performance and implementation constraints
================================================================================

- keep TPT entirely in post-processing
- no per-frame heap-heavy TPT work in the streaming loop
- operate on dense `Eigen::MatrixXd` initially because MSM state counts are
  usually moderate
- only optimize to sparse graph kernels if state counts become large enough to
  justify it
- use `TitanUI` consistently for terminal/log output

================================================================================
15. Recommended immediate implementation order
================================================================================

Immediate sprint:

1. add `TptArgs`
2. create `TitanTPT`
3. add micro/macro basin parsing
4. wire `run_tpt_postproc(...)` after MSM/DHAM
5. export `summary/committor/flux/paths`
6. validate on LJ NVT first

Next sprint:

1. add DHAM-backed ensemble-TPT for LJ NPT
2. add LJ threshold basins
3. add protein macrostate TPT

Later sprint:

1. kinetic-TPT mode with stricter warnings and validation

================================================================================
16. Post-phase enhancement backlog
================================================================================

These items are worth adding after the current core/ensemble workflow is stable.

16.1 Visualization and interpretation

- reactive FES and 2D flux-map export on the TICA plane
- optional macrostate-level flux map
- direct plotting helpers for TITAN output

16.2 Model quality and uncertainty

- bootstrap uncertainty for committors, fluxes, dominant paths, rates, and
  MFPT [implemented and smoke-tested for MSM-backed / demux-backed models;
  DHAM-backed bootstrap pending]
- stronger kinetic-TPT guardrails
  - round-trip sufficiency
  - occupancy imbalance
  - optional CK / Markov sanity checks [first-pass implemented]

16.3 Data preparation helpers

- DHAM window-matching helper to align trajectory and energy windows before
  reweighting [implemented]
- demultiplexed trajectory support for kinetic-TPT [implemented]
- automatic checks for missing or truncated energy windows [implemented with
  `dham_window_audit.dat`]

16.4 Basin-definition presets

- LJ presets
  - gas / dense-liquid / solid-like basin helpers from cluster/order metrics
    [implemented via `ljpreset:*` or `lj:*` preset names]
- protein presets
  - macrostate presets [implemented via `protein:major_macro` and
    `protein:minor_macro`]
  - TICA/PCA rectangles
  - RMSD/DSSP-informed basin helpers [implemented; e.g. `protein:rmsd_low`,
    `protein:helical`, `protein:coil`]
  - contact-informed basin helpers [pending]

16.5 Richer kinetics summaries

- all-pairs macrostate MFPT / rate summaries [implemented]
- pathway clustering / filtering to merge nearly equivalent dominant paths
  [implemented]
- bottleneck-state and bottleneck-edge ranking tables [implemented]
