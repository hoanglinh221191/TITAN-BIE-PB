# TITAN Standard MD Analysis Roadmap

Tai lieu nay bo sung roadmap tinh nang phan tich MD chuan cho TITAN, tap trung
vao cac observable con thieu so voi pipeline protein / membrane MD hien dai.

Muc tieu khong phai viet them solver lon moi, ma la mo rong lop observable va
collective variables de cac module da co san nhu TICA, MSM, TPT, graph,
validation, MBAR-context, va AWH-context co dau vao khoa hoc tot hon.

Ngay tao: 2026-05-26

================================================================================
1. Nguyen tac thiet ke
================================================================================

- Giu module moi o che do opt-in, khong lam doi output cu neu user khong bat
  flag moi.
- Tai su dung trace selection, PBC unwrap/MIC, neighborhood grid, sparse graph,
  TitanExecutor, va validation harness hien co.
- Moi observable moi phai co:
  - timeseries per-frame khi no co y nghia nhu CV
  - summary theo residue/contact/lipid/selection
  - output dang CSV/dat on dinh de feed vao Python/R/plot tool
  - smoke fixture nho va native/portable regression neu co the
- Cac feature co the lam CV phai export duoc vao shared feature-matrix layer de
  TICA/MSM/TPT consume ma khong copy lai trajectory.
- Contact map residue-residue phai dung side-chain contact lam mac dinh, khong
  dung CA/CB distance lam dinh nghia scientific default.

================================================================================
2. Phase 0: Shared Observable / CV Foundation
================================================================================

Muc tieu:

- Tao lop dung chung de cac module moi xuat feature timeseries vao TICA/MSM/TPT.
- Chuan hoa selection, reference, frame weights, va PBC policy cho standard MD
  analysis.

Cong viec:

- Them `titan_observable_features.{h,cpp}` hoac helper tuong duong:
  - dang ky feature channel: `dccm`, `lj`, `dihedral`, `qnative`, `cmap`,
    `msd`, `membrane`, `s2`
  - export matrix `features_by_frame.csv` khi user bat explicit flag
  - ho tro feature metadata: name, residue id, atom ids, unit, periodic angle
- Them CLI chung:
  - `-feature-export`
  - `-feature-source <auto|dccm|lj|dihedral|qnative|cmap|mixed>`
  - `-feature-stride <N>`
  - `-feature-standardize <none|zscore|circular>`
- TICA integration:
  - giu default cu: protein-like jobs van di qua network/DCCM neu user khong
    chon feature source moi
  - khi user chon `-feature-source dihedral|qnative|mixed`, TICA doc feature
    matrix moi thay vi chi dung DCCM trace
  - circular features cua dihedral nen export thanh `sin(theta), cos(theta)`
    thay vi angle raw cho TICA

Exit criteria:

- Default `-tica` output khong doi neu khong bat feature-source moi.
- `-feature-export` co fixture nho tao duoc feature matrix va metadata.
- TICA smoke voi feature source moi sinh `tica_*` binh thuong.

================================================================================
3. Phase 1: Side-chain Residue Contact Probability Map (`-cmap`)
================================================================================

Muc tieu:

- Bo sung contact frequency map standalone, khac voi DCCM correlation va khac
  voi IFP typed interactions.
- Tao ma tran xac suat contact residue-residue theo dinh nghia side-chain COM;
  matrix la vuong cho `mol1` va chu nhat cho `interface`.

Scientific default:

- Dinh nghia contact mac dinh: khoang cach giua hai mass-weighted side-chain
  COM <= 6.5 A.
- `sidechain-com` dung dung semantic legacy cua `MOL::take_com`: residue khong
  co atom side-chain hop le fallback ve COM toan residue; glycine theo
  side-chain marking hien co.
- Ho tro tuy chon anchor khac chi nhu compatibility / speed mode:
  - `sidechain-com` default
  - `sidechain-heavy` atom-pair packing mode
  - `cb` quick mode
  - `ca` legacy / coarse mode
  - `bb` backbone diagnostic only

CLI de xuat:

- `-cmap`
- `-cmap-scope <mol1|interface>`
- `-cmap-mode <sidechain-com|sidechain-heavy|cb|ca|bb>`
- `-cmap-cut <A>` default 6.5 cho `sidechain-com`, 4.5 cho atom-based mode
- `-cmap-domain <intra|inter|both>`
- `-cmap-min-seqsep <N>` default 3 cho intra-chain protein
- `-cmap-weight-source <uniform|mbar|auto>` chi cho equilibrium contact
  probability, khong doi kinetics
- `-cmap-export-frame-bits <true|false>` de xuat sparse per-frame contact
  timeseries khi can feed TICA / clustering

Outputs:

- `contact_probability_matrix.csv`
  - matrix residue x residue, value = occupancy probability
- `contact_probability_edges.csv`
  - sparse edge list: residue_i, residue_j, probability, mean_distance,
    min_distance, contact_frames
- `contact_probability_summary.dat`
  - cutoff, mode, domain, frames, residues, contact count, dense/sparse stats
- `contact_count_timeseries.csv`
  - so contact residue-pair theo frame
- `contact_residue_occupancy.csv`
  - occupancy theo residue va role molecule
- `contact_probability_{matrix,edges,summary,residue_occupancy}_average.*`
  - pooled outputs cho multi-trajectory
- `contact_probability_bfactor.pdb` optional
  - per-residue contact degree/strength mapped to B-factor

Integration:

- `-graph-source contact` can consume `contact_probability_edges.csv`.
- `-community=louvain` can run on side-chain contact network without DCCM.
- `-allostery-knn-prior=contact` should reuse this contact prior when available
  instead of recomputing a separate contact heuristic.
- `-awh-ifp` can optionally overlay untyped contact occupancy per RC bin later.

Implementation notes:

- Use shared `HalfShellGrid` / neighbor backend for side-chain atom pairs.
- Store per-residue side-chain atom lists once from topology.
- For multi-trajectory aggregation:
  - report per-trajectory occupancy
  - report global occupancy
  - keep SEM across trajectories when inputs are independent replicates

Verification:

- Unit fixture with known residue side-chain contacts.
- Intra-chain `min_seqsep` guard test.
- Native/portable hash or numeric tolerance check for matrix/edges.
- Confirm `-cmap` does not alter IFP/DCCM outputs when combined.

================================================================================
4. Phase 2: Native Contacts Q (`-qnative`)
================================================================================

Muc tieu:

- Them fraction of native contacts Q nhu reaction coordinate chuan cho folding,
  unfolding, binding, and conformational state definition.
- Reuse contact-map core de tranh hai dinh nghia contact khac nhau.

Scientific default:

- Native contacts duoc xac dinh tu reference structure bang mass-weighted
  side-chain COM contact <= 6.5 A, cung default voi `-cmap`.
- `Q(t) = formed_native_contacts(t) / total_native_contacts`.
- Contact formed co the dung hard cutoff hoac soft switching function.

CLI de xuat:

- `-qnative`
- `-qnative-ref <pdb|gro>`
- `-qnative-scope <mol1|interface>`
- `-qnative-mode <sidechain-com|sidechain-heavy|cb|ca|bb>`
- `-qnative-cut <A>` default 6.5 cho `sidechain-com`, 4.5 cho atom modes
- `-qnative-soft <true|false>`
- `-qnative-beta <val>` cho switching function
- `-qnative-domain <intra|inter|both>`
- `-qnative-min-seqsep <N>`

Outputs:

- `qnative_timeseries.csv`
  - frame, time_ps, Q, formed_native_contacts, total_native_contacts
- `qnative_native_contacts.csv`
  - native contact list from reference
- `qnative_contact_occupancy.csv`
  - each native contact occupancy over trajectory
- `qnative_summary.dat`
  - mean/std/min/max Q, folded/unfolded occupancy by threshold if requested

Integration:

- TICA:
  - export Q as one scalar CV, and optionally native-contact bits as a sparse
    feature block.
- MSM:
  - allow clustering in mixed space: TICA components + Q.
  - optional state labels by Q range for interpretation.
- TPT:
  - add source/sink selectors:
    - `qnative_range:min,max`
    - `qnative_lt:x`
    - `qnative_gt:x`
  - this complements current `tica_rect` selector and makes folding paths much
    easier to define.
- AWH:
  - Q can be overlayed as exploratory context per RC bin.

Verification:

- Synthetic 3-frame trajectory: Q = 1.0, partial, 0.0.
- Reference/contact-list parser handles residue numbering and chain IDs.
- TPT selector test for `qnative_range`.
- Default TPT behavior unchanged without Q selector.

Slice B implementation ledger:

- `QNativeArgs` mirrors CMap defaults and validation.
- `titan_qnative.{h,cpp}` streams hard or soft Q from the same contact
  primitive semantics as CMap.
- Outputs are `qnative_timeseries.csv`, `qnative_native_contacts.csv`,
  `qnative_contact_occupancy.csv`, and `qnative_summary.dat`, with indexed and
  `_average` forms for multi-input runs.
- TPT basin parser accepts `qnative_range:min,max`, `qnative_lt:x`, and
  `qnative_gt:x`; selectors resolve on per-microstate mean Q.
- Remaining future item: mixed-source feature export into TICA stays after
  qnative + dihedral.

================================================================================
5. Phase 3: Dihedral Analysis (`-dihedral`)
================================================================================

Muc tieu:

- Bo sung backbone/side-chain dihedral analysis:
  - phi/psi Ramachandran per residue
  - omega cis/trans diagnostics
  - chi1/chi2 rotamer populations
  - rotamer transition counts
  - dihedral timeseries for TICA additional CVs

CLI de xuat:

- `-dihedral`
- `-dihedral-kind <backbone|sidechain|both>`
- `-dihedral-residues <all|range/list>`
- `-dihedral-export-tica <true|false>`
- `-dihedral-rotamer-lib <default|none>`
- `-dihedral-stride <N>`

Outputs:

- `dihedral_backbone_timeseries.csv`
  - frame, time_ps, residue, phi, psi, omega
- `dihedral_sidechain_timeseries.csv`
  - frame, time_ps, residue, chi1, chi2 in v1; chi3/chi4 optional future
- `dihedral_ramachandran.csv`
  - residue, phi_bin, psi_bin, occupancy
- `dihedral_rotamer_populations.csv`
  - residue, chi state, occupancy
- `dihedral_rotamer_transitions.csv`
  - residue, chi, from_state, to_state, count
- `dihedral_feature_matrix.csv` optional
  - sin/cos encoded angles for TICA
- `dihedral_summary.dat`

Implementation notes:

- Build residue templates for protein atom names:
  - phi: C(i-1)-N(i)-CA(i)-C(i)
  - psi: N(i)-CA(i)-C(i)-N(i+1)
  - omega: CA(i)-C(i)-N(i+1)-CA(i+1)
  - chi definitions by residue type
- Missing atoms produce `nan` angle samples, not fatal error.
- Shared SoA torsion kernel uses TitanExecutor, aligned vectors, and AVX2
  MIC/bond-vector preparation when available; final atan2 remains scalar per
  torsion for accuracy and portability.

Integration:

- TICA feature source:
  - `-feature-source dihedral`
  - `-feature-source mixed` can combine DCCM trace, Q, and dihedral sin/cos
- MSM/TPT:
  - use dihedral-informed TICA to improve slow-mode extraction for
    conformational changes.
- Validation:
  - VAMP-2 candidate sweeps can compare DCCM-only vs dihedral vs mixed feature
    spaces later.

Verification:

- Static peptide fixture with known phi/psi.
- Rotamer transition fixture for chi1 state changes.
- Angle periodicity test around -180/180 degrees.
- TICA smoke with sin/cos feature matrix.

================================================================================
6. Phase 4: Standalone MSD / Diffusion (`-msd`)
================================================================================

Muc tieu:

- Tach MSD/diffusion khoi LLPS dense-phase-only output.
- Slice 4A ho tro translational diffusion va anomalous exponent alpha.
- Slice 4B ho tro rotational diffusion bang Kabsch orientation quaternion va
  Favro quaternion covariance tensor.

CLI de xuat:

- `-msd`
- `-msd-selection <mol1|mol2|trace|backbone|atom:N,CA|residue:10,12-20,A:30>`
- `-msd-mode <com|atom>`
- `-msd-maxlag <N|auto>` default `auto = min(frames - 1, 1000)`
- `-msd-fit-window <auto|tmin,tmax>`
- `-msd-pbc <auto|unwrap|mic>`
- `-msd-rot <true|false>`; `true` duoc ho tro khi co `-msd=true`
- `-msd-rot-fit-window <auto|tmin,tmax>`

Outputs:

- `msd_timeseries.dat`
  - lag_index, lag_time_ps, MSD_A2, samples
- `msd_diffusion_summary.dat`
  - frames, points, selection, PBC policy, fit window, slope, D, R2
- `msd_alpha_fit.dat`
  - anomalous exponent alpha from MSD ~ t^alpha
- `rotational_msd.dat` optional
  - lag index, lag time, rotational MSD, samples
- `rotational_quaternion_timeseries.csv` optional
  - orientation quaternion theo frame
- `rotational_covariance_fit.csv` optional
  - Favro quaternion covariance tensor va fit diagnostics
- `rotational_diffusion_summary.dat` optional
  - rotational diffusion fit, `reference_note`, va finite-size note

Implementation notes:

- Translational MSD needs unwrapped coordinates or robust image correction.
- Rotational MSD dung orientation quaternion tu Kabsch fit voi first analyzed
  frame lam reference; fail closed neu selection khong dinh nghia duoc stable
  body frame.
- Rotational diffusion dung Favro quaternion covariance tensor. Summary phai
  giu `reference_note` va finite-size note de interpretation khong mat caveat.
- Multi-time-origin algorithm should avoid O(F^2) blow-up for long trajectories.

Verification:

- Brownian synthetic trajectory with known D.
- Ballistic/constant-velocity sanity fixture.
- PBC wrap/unwrap test.
- LLPS outputs remain unchanged; standalone `-msd` is separate.
- Rotational synthetic trajectory validates quaternion continuity, tensor fit,
  and first-analyzed-frame reference note.

Slice 4A implementation ledger:

- `MsdArgs` adds opt-in standalone translational MSD without changing LLPS
  dense-phase MSD filenames or behavior.
- `titan_msd.{h,cpp}` stores selected COM or atom tracks in aligned flattened
  SoA arrays and computes bounded multi-time-origin lags through
  `TitanExecutor`.
- COM mode follows the hardened mass policy: mass-weighted when possible,
  unit weights when selected total mass is zero.
- `-msd-pbc=auto` unwraps by consecutive-frame MIC when all used frames have
  valid boxes; otherwise it documents raw-coordinate policy in the summary.
- Outputs are `msd_timeseries.dat`, `msd_diffusion_summary.dat`, and
  `msd_alpha_fit.dat`, with `_N` suffixes for multi-input.
- Slice 4B extends `-msd-rot=true` as an opt-in rotational path when `-msd=true`.
- `-msd-rot-fit-window=<auto|tmin,tmax>` controls the rotational covariance fit
  independently from the translational fit.
- Rotational orientation is defined by Kabsch fit to the first analyzed frame.
- Rotational diffusion is derived from the Favro quaternion covariance tensor.
- Outputs are `rotational_msd.dat`, `rotational_quaternion_timeseries.csv`,
  `rotational_covariance_fit.csv`, and `rotational_diffusion_summary.dat`.
- The rotational summary records `reference_note` and a finite-size note.

================================================================================
7. Phase 5: Lipid / Membrane Analysis (`-membrane`)
================================================================================

Muc tieu:

- Ho tro MARTINI membrane/bilayer workflows ma TITAN da co force-field support
  nhung chua co analysis surface.

Core features:

- Leaflet assignment
- Area per lipid (APL)
- Membrane thickness
- Tail order proxy for MARTINI consecutive tail-bead vectors
- Interdigitation and flip-flop diagnostics

CLI:

- `-membrane`
- `-membrane-group <auto|N>`
- `-membrane-lipids <auto|resname:list>`
- `-membrane-head <auto|atomname:list>`
- `-membrane-tail <auto|atomname:list>`
- `-membrane-normal <z|x|y|auto>`
- `-membrane-grid <auto|Nx,Ny>`
- `-membrane-leaflet <auto|upper-lower|user>`
- `-membrane-leaflet-file <csv>`
- `-membrane-flipflop-persistence <N>`

Outputs:

- `membrane_leaflets.csv`
  - lipid id, leaflet assignment over time or consensus leaflet
- `membrane_apl_timeseries.csv`
- `membrane_apl_per_lipid.csv`
- `membrane_thickness_timeseries.csv`
- `membrane_thickness_map.csv`
- `membrane_order_parameters.csv`
  - lipid type, tail, bead/bond index, S_CD-like value
- `membrane_interdigitation.dat`
- `membrane_flipflop_events.csv`
- `membrane_summary.dat`
- multi-input aggregates: `membrane_summary_average.dat`,
  `membrane_thickness_map_average.csv`, and
  `membrane_order_parameters_average.csv`

Implementation notes:

- MARTINI first:
  - use phosphate/head beads for leaflet and thickness
  - use tail bead vectors for coarse-grained order proxies
- Atomistic path can be added after MARTINI defaults are stable.
- Leaflet assignment uses local periodic PCA midplanes for auto mode, not only
  a global z split.
- Lipids use a dedicated coordinate stream and do not consume `mol2`.

Integration:

- Thickness/APL/order timeseries can feed feature export later.
- Contact map can support protein-lipid side-chain/lipid contact extension in a
  later phase, but that is not in v1.

Verification:

- Small MARTINI bilayer fixture with known two-leaflet topology.
- Tilted/undulating bilayer sanity case.
- Unknown lipid templates warn and skip. A requested analysis with fewer than
  four supported lipids fails clearly instead of emitting misleading output.

================================================================================
8. Phase 6: NMR Order Parameters S2 (`-s2`)
================================================================================

Muc tieu:

- Tinh Lipari-Szabo model-free order parameters S2 tu trajectory de validate MD
  voi NMR relaxation data va mo ta rigidity/flexibility theo residue.

Scientific core:

- For each bond vector mu:
  - C(t) = <P2(mu(0) dot mu(t))>
  - S2 = plateau C(t -> infinity) or long-lag average
- Primary vectors:
  - exact backbone N-H using H or HN atoms
  - exact C-alpha to C bond as an independent row

Public arg.in contract:

- `-s2=true|false` or bare `-s2`
- `-s2-vector=<nh|ca-c|both>`; default `both`
- `-s2-align=<rms|none>`; default `rms`
- `-s2-maxlag=<auto|N>`; default `auto`, explicit `N >= 1`
- `-s2-fit-window=<auto|tmin,tmax>`; default `auto`, require `tmax > tmin`
- `-s2-residues=<all|list/range>`
- `-s2-bfactor-pdb=<true|false>`; default `false`

V1 compatibility notes:

- `-s2-align=provided` is rejected clearly because v1 supports first-frame RMS
  alignment or no alignment only.
- Legacy `-s2-plateau-window` is rejected clearly; use `-s2-fit-window`.

Outputs:

- `s2_order_parameters.csv`
  - chain, residue_number, residue_name, vector_type, plateau_C_raw,
    clipped S2, tau_ps, fit_R2, fit_window_start_ps, fit_window_end_ps,
    samples, fit_status
- `s2_correlation_functions.dat`
  - chain, residue_number, residue_name, vector_type, lag_index, lag_time_ps,
    C(t), samples
- `s2_summary.dat`
- `s2_bfactor.pdb` optional
- multi-input aggregates: `s2_order_parameters_average.csv`,
  `s2_summary_average.dat`, and optional `s2_bfactor_average.pdb`

Implementation notes:

- N-H rows require exact H/HN atoms. PRO and terminal residues without a
  backbone hydrogen are skipped with a warning.
- CA-C uses the exact vector and is emitted as an independent row.
- `-s2-align=rms` uses Kabsch alignment to the first analyzed frame.
- TitanExecutor operates on aligned SoA buffers with AVX2/FMA kernels and a
  scalar fallback.
- Fit each trajectory independently. Multi-input aggregates report mean and
  SEM; do not pool independent trajectories into one C(t).
- Use bounded lag windows to avoid noisy tails being interpreted as plateau.

Verification:

- Rigid-body trajectory should produce S2 near 1 after alignment.
- Randomized vector trajectory should decay toward low S2.
- Missing-H fixture should warn predictably.

================================================================================
9. Recommended Execution Order
================================================================================

1. `-cmap` side-chain residue contact probability map
   - Highest immediate value.
   - Reuses neighbor backend and sparse graph.
   - Also becomes foundation for Q native and contact-prior allostery.

2. `-qnative`
   - Small extra solver cost after contact core exists.
   - Unlocks folding/binding CVs and `qnative_range` TPT selectors.

3. `-dihedral`
   - Important for protein conformational changes.
   - Unlocks better TICA features when DCCM trace is insufficient.

4. Shared `-feature-export` / TICA mixed feature integration
   - Do after `cmap/qnative/dihedral` have real feature producers.
   - Avoid designing an abstract feature system with no users.

5. Standalone `-msd`
   - Useful and self-contained, but less coupled to MSM/TPT.

6. `-membrane`
   - Larger domain-specific surface; start with MARTINI-only v1.

7. `-s2`
   - Valuable validation observable, but depends on alignment/vector handling
     quality and careful plateau diagnostics.

================================================================================
10. Compatibility / Risk Register
================================================================================

- Contact definition drift:
  - `-cmap` and `-qnative` must share the same contact primitive by default.
  - Mass-weighted side-chain COM contact at `<= 6.5 A` is the scientific
    default.
  - Side-chain heavy atom contact is an explicit packing/detail mode.
  - CA/CB modes are explicit fast/legacy options, not the main definition.

- TICA feature ambiguity:
  - Default `-tica` must keep current behavior.
  - New features require explicit `-feature-source` or module-specific
    `-dihedral-export-tica`.

- PBC / unwrap:
  - Dihedral, MSD, membrane, and contact analysis need different PBC handling.
  - Each module must log its PBC policy and fail closed when input is ambiguous.

- Multi-trajectory semantics:
  - For probability outputs, report both pooled occupancy and replicate SEM
    where meaningful.
  - Do not silently treat independent replicates as one long kinetic trajectory
    for MSM/TPT.

- Runtime cost:
  - Side-chain contact maps can be expensive if done naively.
  - Use neighbor grid and sparse outputs by default; dense matrix is output
    artifact, not internal working representation for every frame.

- MARTINI atom naming:
  - Membrane v1 should define supported lipid templates and clear skip/warn
    behavior for unknown lipids.

================================================================================
11. Roadmap Acceptance Gates
================================================================================

For each shipped module:

- Parser/help/manual entries added together.
- `arg.in` syntax and CLI syntax both verified if module has new options.
- Smoke fixture in `test/` with deterministic outputs.
- Native and portable build pass.
- Existing benchmark cases for DCCM/TICA/MSM/TPT/IFP remain unchanged unless
  the user explicitly opted into new feature source.
- New output files documented in `manual.txt`, `help.cpp`, and onboarding notes.
- If a module feeds TICA/MSM/TPT, add at least one end-to-end smoke proving the
  feature source can be consumed without changing default behavior.

================================================================================
12. First Implementation Slice Proposal
================================================================================

Slice A: `-cmap` side-chain COM contact map v2

Status: implemented in this branch.

- Replace legacy map output with side-chain COM contact probability output.
- Add `ContactPrimitive` shared between future `-cmap` and `-qnative`.
- Implement `-cmap` parser/options/logging.
- Compute mass-weighted side-chain COM occupancy; retain atom modes explicitly.
- Export:
  - `contact_probability_matrix.csv`
  - `contact_probability_edges.csv`
  - `contact_probability_summary.dat`
  - `contact_count_timeseries.csv`
  - `contact_residue_occupancy.csv`
- Add smoke test with known contacts.
- Update docs/help.

Implemented notes:

- Core contact helpers live in `titan_contact_primitives.{h,cpp}`.
- Streaming analyzer lives in `titan_contact_map.{h,cpp}` and is wired through
  `impl_engine.cpp`.
- CLI options are `-cmap`, `-cmap-scope`, `-cmap-mode`, `-cmap-cut`,
  `-cmap-domain`, and `-cmap-min-seqsep`.
- Default is `sidechain-com` at `6.5 A`; `sidechain-heavy` is explicit.
- Deprecated `-contact-map=yes` migrates to interface sidechain-COM CMap.
- End-to-end smoke is `test/cmap_smoke.ps1`.
- XDR and multi-input runs may use indexed output names such as
  `contact_probability_edges_0.csv` to avoid overwriting per-trajectory files.

Slice B: `-qnative` v1

Status: implemented in this branch.

- Reuse `ContactPrimitive`.
- Build native contact list from reference.
- Export Q timeseries and native contact occupancy.
- Add TPT selector parser for `qnative_range`.
- Add synthetic folding/unfolding smoke.

Slice C: `-dihedral` backbone + sidechain chi1/chi2 v1

Status: implemented in this branch.

- Implement phi/psi/omega backbone torsions.
- Implement sidechain chi1/chi2 rotamer populations and transition counts.
- Export backbone/sidechain timeseries, Ramachandran, summary, and sin/cos
  feature matrix for TICA prep.
- Use TitanExecutor + aligned SoA torsion kernel with AVX2 MIC prep.
- Chi3/chi4 and direct mixed-source TICA integration remain future optional
  expansions.

Slice 4A: `-msd` translational MSD v1

Status: implemented in this branch.

- Add standalone translational MSD separate from LLPS dense-phase MSD.
- Support `mol1`, `mol2`, `trace`, `backbone`, `atom:<names>`, and
  `residue:<ranges>` selections.
- Support COM and atom-averaged modes with aligned SoA storage.
- Compute bounded multi-time-origin lags through `TitanExecutor`.
- Export `msd_timeseries.dat`, `msd_diffusion_summary.dat`, and
  `msd_alpha_fit.dat` with indexed multi-input filenames.
- Keep `-msd-rot=false` as the default.

Slice 4B: `-msd-rot=true` rotational diffusion v1

Status: implemented in this branch. Phase 4 complete.

- Support `-msd-rot=true` together with `-msd=true`.
- Add independent `-msd-rot-fit-window=<auto|tmin,tmax>`.
- Use the first analyzed frame as the Kabsch reference orientation.
- Derive rotational diffusion from the Favro quaternion covariance tensor.
- Export `rotational_msd.dat`, `rotational_quaternion_timeseries.csv`,
  `rotational_covariance_fit.csv`, and `rotational_diffusion_summary.dat`.
- Preserve `reference_note` and finite-size note in the rotational summary.

Completed target: Phase 5 membrane analysis.

Phase 4, Phase 5, and the standalone Phase 6 S2 implementation are complete in
this branch.

Phase 5: membrane analysis

Status: implemented in this branch. Phase 5 complete.

Membrane analysis uses a dedicated lipid coordinate stream, so
`-nocal=true -membrane=true` does not enable the energy pipeline or consume
`mol2`. Auto normal uses PCA with sign continuity. Auto leaflet assignment
uses periodic local-neighborhood PCA midplanes to remain useful for undulating
bilayers.

Public arg.in contract:

- `-membrane=true|false` or bare `-membrane`
- `-membrane-group=<auto|N>`
- `-membrane-lipids=<auto|resname:list>`
- `-membrane-head=<auto|atomname:list>`
- `-membrane-tail=<auto|atomname:list>`
- `-membrane-normal=<auto|x|y|z>`
- `-membrane-grid=<auto|Nx,Ny>`
- `-membrane-leaflet=<auto|upper-lower|user>`
- `-membrane-leaflet-file=<csv>`
- `-membrane-flipflop-persistence=<N>`

Defaults are inert: membrane analysis is off, selector fields use `auto`,
`-membrane-leaflet-file` is empty, and flip-flop persistence is 3 frames.
Validation requires a non-negative explicit group, positive grid dimensions,
positive persistence, and a leaflet CSV exactly when leaflet mode is `user`.

Slice 5A: membrane APL + thickness v1

Status: implemented in this branch.

- Resolve lipid/head/tail selectors and membrane normal.
- Assign leaflets automatically, as upper/lower, or from a user CSV.
- Export `membrane_leaflets.csv`, `membrane_apl_timeseries.csv`,
  `membrane_apl_per_lipid.csv`, `membrane_thickness_timeseries.csv`,
  `membrane_thickness_map.csv`, and `membrane_summary.dat`.
- Add aggregate `membrane_summary_average.dat` and
  `membrane_thickness_map_average.csv`.

Slice 5B: order parameters + interdigitation + flip-flop v1

Status: implemented in this branch.

- Export `membrane_order_parameters.csv`, `membrane_interdigitation.dat`, and
  `membrane_flipflop_events.csv`.
- Require persistent leaflet reassignment for the configured number of frames
  before reporting a flip-flop event.
- Add aggregate `membrane_order_parameters_average.csv`.

Implementation notes:

- APL uses periodic tiled-neighbor Voronoi clipping per leaflet.
- Auto grid sizing is aspect-aware and targets about 8 lipids per leaflet per
  cell, clamped to 4..64 bins per axis.
- Headgroup projection has an AVX2 path with scalar fallback. Voronoi polygon
  clipping remains scalar for numerical stability.
- `membrane_summary.dat` records a `reference_note` explaining the head-bead
  geometry, PCA normal, and local-midplane policy.

Phase 6: NMR order parameters S2

Status: implemented in this branch. Phase 6 standalone S2 runtime integration
is complete.

- Add inert `S2Args` defaults and parser support for bare or equals booleans.
- Preserve `-s2=false` as one token when parsing `arg.in`.
- Validate vector kind, alignment policy, maximum lag, fit window, and
  non-empty residue selection.
- Emit exact H/HN N-H rows, warn and skip PRO or terminal residues without a
  backbone hydrogen, and emit exact CA-C vectors as independent rows.
- Use first-analyzed-frame Kabsch alignment and TitanExecutor aligned SoA
  buffers with AVX2/FMA kernels plus scalar fallback.
- Fit each trajectory independently and aggregate mean and SEM without pooling
  independent C(t) curves.
- Export indexed per-trajectory results and stable aggregate outputs.
- Cover parser behavior and standalone end-to-end runtime with S2 smoke tests.
