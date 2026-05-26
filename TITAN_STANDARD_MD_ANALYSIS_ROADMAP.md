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
- Tao ma tran xac suat contact residue-residue N x N cho moi residue pair theo
  dinh nghia side-chain contact.

Scientific default:

- Dinh nghia contact mac dinh: bat ky cap side-chain heavy atom nao giua hai
  residue co khoang cach <= cutoff.
- Glycine khong co side-chain: dung CA hoac pseudo-CB fallback va danh dau ro
  trong metadata.
- Loai bo backbone atom khoi default de contact map phan anh packing / binding
  interface hon la proximity cua backbone.
- Ho tro tuy chon anchor khac chi nhu compatibility / speed mode:
  - `sidechain-heavy` default
  - `cb` quick mode
  - `ca` legacy / coarse mode
  - `bb` backbone diagnostic only

CLI de xuat:

- `-cmap`
- `-cmap-mode <sidechain-heavy|cb|ca|bb>`
- `-cmap-cut <A>` default 4.5 cho side-chain heavy atom
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
- `contact_frame_contacts.csv` optional
  - sparse per-frame contact bits for downstream CV/TICA
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

- Native contacts duoc xac dinh tu reference structure bang side-chain heavy atom
  contact, cung default voi `-cmap`.
- `Q(t) = formed_native_contacts(t) / total_native_contacts`.
- Contact formed co the dung hard cutoff hoac soft switching function.

CLI de xuat:

- `-qnative`
- `-qnative-ref <pdb|gro>`
- `-qnative-mode <sidechain-heavy|cb|ca>`
- `-qnative-cut <A>`
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
  - frame, time_ps, residue, chi1, chi2, chi3, chi4 where available
- `dihedral_ramachandran.csv`
  - residue, phi_bin, psi_bin, occupancy
- `dihedral_rotamer_populations.csv`
  - residue, chi state, occupancy
- `dihedral_rotamer_transitions.csv`
  - residue, from_state, to_state, count, rate_per_ns
- `dihedral_feature_matrix.csv` optional
  - sin/cos encoded angles for TICA
- `dihedral_summary.dat`

Implementation notes:

- Build residue templates for protein atom names:
  - phi: C(i-1)-N(i)-CA(i)-C(i)
  - psi: N(i)-CA(i)-C(i)-N(i+1)
  - omega: CA(i)-C(i)-N(i+1)-CA(i+1)
  - chi definitions by residue type
- Missing atoms produce explicit skipped-residue diagnostics, not fatal error.
- PBC unwrap should be applied before angle calculation when needed.

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
- Ho tro selection tuy y cho translational diffusion, rotational diffusion, va
  anomalous diffusion exponent alpha.

CLI de xuat:

- `-msd`
- `-msd-selection <mol1|mol2|protein|rna|residue-range|atom-selection>`
- `-msd-mode <com|atom|residue|molecule>`
- `-msd-maxlag <N>`
- `-msd-fit-window <tmin,tmax>`
- `-msd-pbc <unwrap|mic|auto>`
- `-msd-rot <true|false>`

Outputs:

- `msd_timeseries.dat`
  - lag_time_ps, MSD_A2, samples
- `msd_diffusion_summary.dat`
  - D_trans from Einstein relation, fit window, R2, uncertainty
- `msd_alpha_fit.dat`
  - anomalous exponent alpha from MSD ~ t^alpha
- `rotational_msd.dat` optional
- `rotational_diffusion_summary.dat` optional

Implementation notes:

- Translational MSD needs unwrapped coordinates or robust image correction.
- Rotational MSD should use quaternion orientation autocorrelation for rigid
  body selections; fail closed if selection cannot define a stable body frame.
- Multi-time-origin algorithm should avoid O(F^2) blow-up for long trajectories.

Verification:

- Brownian synthetic trajectory with known D.
- Ballistic/constant-velocity sanity fixture.
- PBC wrap/unwrap test.
- LLPS outputs remain unchanged; standalone `-msd` is separate.

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
- Tail order parameters for MARTINI / atomistic mapping
- Interdigitation and flip-flop diagnostics

CLI de xuat:

- `-membrane`
- `-membrane-lipids <auto|resname:list>`
- `-membrane-head <auto|atomname:list>`
- `-membrane-tail <auto|atomname:list>`
- `-membrane-normal <z|x|y|auto>`
- `-membrane-grid <Nx,Ny>`
- `-membrane-leaflet <auto|upper-lower|user>`

Outputs:

- `membrane_leaflets.csv`
  - lipid id, leaflet assignment over time or consensus leaflet
- `membrane_apl_timeseries.csv`
- `membrane_thickness_timeseries.csv`
- `membrane_thickness_map.csv`
- `membrane_order_parameters.csv`
  - lipid type, tail, bead/bond index, S_CD-like value
- `membrane_interdigitation.dat`
- `membrane_summary.dat`

Implementation notes:

- MARTINI first:
  - use phosphate/head beads for leaflet and thickness
  - use tail bead vectors for coarse-grained order proxies
- Atomistic path can be added after MARTINI defaults are stable.
- Leaflet assignment should be robust to undulations via midplane/grid model,
  not only global z split.

Integration:

- Thickness/APL/order timeseries can feed feature export later.
- Contact map can support protein-lipid side-chain/lipid contact extension in a
  later phase, but that is not in v1.

Verification:

- Small MARTINI bilayer fixture with known two-leaflet topology.
- Tilted/undulating bilayer sanity case.
- No-op behavior on non-membrane systems: warn and skip, not crash.

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
  - backbone N-H when H exists
  - C-alpha to C bond / surrogate vector when hydrogens are missing

CLI de xuat:

- `-s2`
- `-s2-vector <NH|CA-C|both>`
- `-s2-align <rms|none|provided>`
- `-s2-maxlag <N>`
- `-s2-plateau-window <tmin,tmax>`
- `-s2-residues <all|range/list>`

Outputs:

- `s2_order_parameters.csv`
  - residue, vector_type, S2, plateau_start, plateau_end, samples
- `s2_correlation_functions.dat`
  - residue/vector, lag_time_ps, C(t)
- `s2_summary.dat`
- `s2_bfactor.pdb` optional

Implementation notes:

- Reuse RMSD alignment machinery when `-s2-align rms`.
- If hydrogens are absent, warn and fall back only when requested or when a
  documented surrogate is available.
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
  - Side-chain heavy atom contact is the scientific default.
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

Slice A: `-cmap` side-chain contact map v1

Status: implemented in this branch.

- Add side-chain atom classifier and glycine fallback.
- Add `ContactPrimitive` shared between future `-cmap` and `-qnative`.
- Implement `-cmap` parser/options/logging.
- Compute side-chain contact occupancy using neighbor grid.
- Export:
  - `contact_probability_matrix.csv`
  - `contact_probability_edges.csv`
  - `contact_probability_summary.dat`
- Add smoke test with known contacts.
- Update docs/help.

Implemented notes:

- Core contact helpers live in `titan_contact_primitives.{h,cpp}`.
- Streaming analyzer lives in `titan_contact_map.{h,cpp}` and is wired through
  `impl_engine.cpp`.
- CLI options are `-cmap`, `-cmap-mode`, `-cmap-cut`, `-cmap-domain`, and
  `-cmap-min-seqsep`.
- End-to-end smoke is `test/cmap_smoke.ps1`.
- XDR and multi-input runs may use indexed output names such as
  `contact_probability_edges_0.csv` to avoid overwriting per-trajectory files.

Slice B: `-qnative` v1

- Reuse `ContactPrimitive`.
- Build native contact list from reference.
- Export Q timeseries and native contact occupancy.
- Add TPT selector parser for `qnative_range`.
- Add synthetic folding/unfolding smoke.

Slice C: `-dihedral` backbone v1

- Implement phi/psi/omega only first.
- Export timeseries and Ramachandran occupancy.
- Add sin/cos feature export for TICA.
- Add peptide fixture and angle periodicity tests.

Recommended first target: Slice A.

Ly do: no dung ngay comment cua user ve residue-residue contact probability
matrix, mo khoa graph/allostery/Q-native, va it phu thuoc nhat vao refactor TICA.
