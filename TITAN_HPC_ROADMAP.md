# TITAN HPC / Compute Roadmap

Tai lieu nay chuyen static review compute codebase thanh roadmap thuc thi.
Muc tieu la nang dong thoi:

- hieu nang tren bai toan lon
- do on dinh so hoc va runtime
- kha nang mo rong cho cac module Titan moi

Roadmap uu tien nhung thay doi co tac dong lon, giam debt kien truc, va mo duong
cho cac workflow PCA / MSM / TPT / DCCM / LJ / GIST ve sau.

================================================================================
1. Muc tieu chien luoc
================================================================================

Trong 1-2 vong phat trien tiep theo, Titan nen dat duoc:

- giam memory blow-up o cac path dense nhat
- giam copy / allocation trong hot path streaming va analysis
- giu duoc AVX2/FMA path nhanh nhung co fallback ro rang
- them benchmark va regression gate de moi toi uu hoa deu do duoc
- xay dung backend dung chung cho sparse graph, neighbor search, va flat buffers

Nguyen tac uu tien:

1. Do dac truoc khi toi uu hoa
2. Sua data layout truoc micro-opt
3. Uu tien sparse / low-rank truoc khi them thread
4. Khong mo rong module moi tren nen tang dense de vo han

1.1 Trang thai hien tai (2026-03-30)
================================================================================

Da hoan thanh:

- [done] Phase 0 `Benchmark harness + regression guard` duoi `test/`
  - benchmark / regression deu chay tu thu muc con duoi `test/`
  - da co baseline cho chunk streamer, screened electrostatics, PCA, va GIST
  - output benchmark duoc tong hop thanh markdown summary, doi chieu `nan` /
    `inf`, va ghi peak RAM

- [done] Phase 0 `Portable build profile + local stress suites`
  - `Makefile` va `compile.bat` da co `native` / `portable` profile ro rang
  - `portable` hien tai bo `-march=native` nhung giu explicit
    `AVX2 + FMA + generic tuning`; runtime dispatch da duoc chot o Phase 4
  - `test/benchmark_harness.ps1` da co suite `small` / `medium` / `large`,
    ghi build profile vao summary, va uu tien binary snapshot trong `test/bin`
  - da verify `small portable`, `medium native`, va `large native` bang repo-
    local cases

- [done] Phase 1 / 3.2 `PCA scratch reuse` trong `titan_pca.cpp`
  - `Eigen::MatrixXd frame_mat(num_ca, 3)` duoc tai su dung qua vong lap
    alignment
  - `std::vector<double> dev(dimensions)` duoc tai su dung qua vong lap
    covariance
  - regression `example/3` da duoc doi chieu hash-identical cho output PCA

- [done] Phase 1 / 3.3 `GIST reduction` trong `titan_gist.cpp`
  - bo full-grid clone cho moi thread
  - thay bang sparse per-thread accumulator + touched-voxel serial merge
  - doi chieu `nt=1` va `nt=4` cho output `gist_*.dx` khop nhau

- [done] Phase 1 / 3.4 `Defensive numerical guards`
  - them validation input cho screened electrostatics trong `arguments.cpp`
  - them near-zero / invalid distance guard trong `energy.cpp` va
    `energy_intrin.cpp`
  - fail-fast ro rang cho `salt_concentration <= 0` va `epsilon_in <= 0`

- [done] Phase 1 / 3.1 `Chunk streamer copy removal`
  - PDB split path khong con deep-copy `init_pro.m_mol[0]` /
    `init_lig.m_mol[0]` cho moi frame
  - giu nguyen output path va semantics cua XDR / topology

- [done] Phase 2 / 4.1 `Flat trajectory backend` slice 1 trong `dccm_network.cpp`
  - `cov`, `dev_x/dev_y/dev_z`, va LSP-MD transposed caches da chuyen sang
    flat/aligned row-major helper
  - LSP-MD da bo allocation scratch theo tung pair
  - regression `screening_valid` va `pca_example3` da pass, output hash-identical
    voi baseline cho `lspmd_*.csv`, `network_edges.csv`, `pca_fel.csv`, va
    `pca_fes_2d.dat`

- [done] Phase 2 / 4.1 `Flat trajectory backend` slice 2 trong `dccm_network.cpp`
  - `calculate_allostery_information_optimized()` da chuyen `proj_bins` sang
    flat/aligned row-major
  - entropy / TE kernels doc row pointers truc tiep, bo nested-vector access
  - regression `dccm_allostery`, `pca_post`, va `lspmd_post` da pass,
    hash-identical voi baseline

- [done] Phase 2 / 4.1 `Flat trajectory backend` slice 3 trong `dccm_network.cpp`
  - legacy `calculate_allostery_information()` da gom vao shared flat kernel
  - optimized allostery dispatch chon `uint8_t / uint16_t / uint32_t` theo so
    bin, khong con fallback sang legacy path
  - regression `dccm_allostery`, `pca_post`, va `lspmd_post` van hash-identical
    voi baseline

- [done] Phase 2 / 4.1 `Flat trajectory backend` slice 4 trong
  `engine_postproc.cpp` / `titan_pca.cpp` / `dccm_network.*`
  - PCA khong con goi `get_ca_traj_aos()` va khong con tao AoS trajectory copy
    tu DCCM
  - `TitanPCA` nhan truc tiep `traj_x/y/z` va materialize flat row-major aligned
    buffer cho alignment / covariance / projection
  - regression `dccm_allostery`, `pca_post`, va `lspmd_post` van
    hash-identical voi baseline

- [done] Phase 3 / 5.2 `Sparse TPT` slice 1 trong `titan_tpt.cpp` /
  `titan_tpt.h`
  - bo persistent dense `net_flux_`; flux graph duoc emit truc tiep vao
    `edges_` sau flux reduction
  - dominant path extraction da chuyen sang sparse residual edge list thay vi
    clone dense matrix
  - regression core `test/tpt_sparse_flux_unit_driver.cpp` pass cho `nt=1` va
    `nt=4`; smoke end-to-end `test/tpt_sparse_phase3_cli_post_escalated` giu
    hash-identical cho `tpt_committor.dat`, `tpt_flux_edges.dat`,
    `tpt_flux_map_2d.dat`, `tpt_paths.dat`, `tpt_reactive_fes_2d.dat`, va
    `tpt_state_labels.dat`

- [done] Phase 3 / 5.2 `Sparse TPT` slice 2 trong `titan_tpt.cpp`
  - `solve_committor_system()` da chuyen sang `Eigen::SparseMatrix` +
    `SparseLU`, kem dense fallback de giu do ben cho case kho
  - regression core da doi chieu `q_forward` / `q_backward` voi dense reference
    trong `test/tpt_sparse_flux_unit_driver.cpp`
  - smoke `test/tpt_sparse_phase3_committor_post2` giu hash-identical cho
    `tpt_committor.dat`, `tpt_flux_edges.dat`, `tpt_flux_map_2d.dat`,
    `tpt_paths.dat`, `tpt_reactive_fes_2d.dat`, va `tpt_state_labels.dat`;
    cac field khoa hoc chung trong `tpt_summary.dat` va `tpt_kinetics.dat`
    van khop baseline

- [done] Phase 3 / 5.2 `Sparse TPT` slice 3 trong `titan_tpt.cpp` /
  `titan_tpt.h`
  - bo persistent dense `reverse_tpm_`; backward committor duoc build truc tiep
    tu `tpm_` + `stationary_dist_`
  - `compute_flux_matrices()` giu `partial_edges` ngay trong worker, bo luot
    quet dense pair thu hai khi sinh `edges_`
  - smoke `test/tpt_sparse_phase3_post_final` giu hash-identical cho
    `tpt_committor.dat`, `tpt_flux_edges.dat`, `tpt_flux_map_2d.dat`,
    `tpt_paths.dat`, `tpt_reactive_fes_2d.dat`, va `tpt_state_labels.dat`;
    cac field khoa hoc chung trong `tpt_summary.dat` va `tpt_kinetics.dat`
    van khop baseline

- [done] Phase 3 / 5.2 `Sparse TPT` slice 4 trong `titan_tpt.cpp` /
  `titan_tpt.h`
  - bo persistent dense `tpm_`; transition graph da duoc normalize va luu bang
    CSR forward / reverse
  - `solve_committor_system()` va `compute_flux_matrices()` da doc adjacency
    sparse thay vi scan dense row, giu semantics khoa hoc qua `noise_floor`
    cuoi cung
  - unit `test/tpt_sparse_flux_unit_driver.cpp` pass lai sau khi doi chieu voi
    dense reference; smoke `test/tpt_sparse_phase4_post2` giu hash-identical
    cho `tpt_committor.dat`, `tpt_flux_edges.dat`, `tpt_flux_map_2d.dat`,
    `tpt_paths.dat`, `tpt_reactive_fes_2d.dat`, `tpt_state_labels.dat`,
    `tpt_summary.dat`, va `tpt_kinetics.dat`

- [done] Phase 2 / 4.2 `Shared sparse graph backend` slice 1 trong
  `titan_sparse_graph.h`, `titan_tpt.*`, va `titan_lj.*`
  - them backend header-only gom `CSRTopology` va `WeightedCSRGraph`
  - TPT va LJ da dung chung sparse backend thay vi tu mang CSR/adjacency rieng
  - TPT regression `test/tpt_sparse_flux_unit_driver.cpp` va smoke
    `test/shared_sparse_graph_tpt_post2` deu khop hash voi baseline
  - LJ smoke `test/shared_sparse_graph_lj_post2` khop cho
    `lj_features_0.dat`, `lj_summary_0.dat`, `lj_timeseries_0.dat`,
    `msm_*.dat`, va `tica_*`; `lj_cluster_rho_r_0.dat` van co sai khac giua
    cac rerun tuong duong nen chua duoc xem la hash gate on dinh

- [done] Phase 2 / 4.2 `Shared sparse graph backend` slice 2 trong
  `titan_sparse_graph.h` va `dccm_network.cpp`
  - DCCM LSP-MD da migrate weighted adjacency / Brandes path sang
    `WeightedCSRGraph`, khong con phu thuoc vao
    `vector<vector<pair<size_t,double>>>`
  - regression all-atom `test/shared_sparse_graph_dccm_lspmd_post1` khop hash
    voi baseline `test/dccm_lspmd_post5` cho `lspmd_edges.csv`,
    `lspmd_degree_centrality.csv`, va `lspmd_betweenness_centrality.csv`
  - regression LJ `test/shared_sparse_graph_lj_post3` giu nguyen
    `lj_features_0.dat`, `lj_summary_0.dat`, `lj_timeseries_0.dat`,
    `msm_*.dat`, va `tica_*`, xac nhan duong LJ va all-atom khong anh huong
    nhau qua backend chung; `lj_cluster_rho_r_0.dat` van la file rerun-variant

- [done] Phase 2 / 4.2 `Shared sparse graph backend` slice 3 trong
  `titan_lj.cpp`
  - da fix race trong `accumulate_radial_density()` bang cach join
    `std::jthread` truoc khi reduce `local_cluster`
  - hai rerun LJ tuong duong `test/shared_sparse_graph_lj_post4a` va
    `test/shared_sparse_graph_lj_post4b` da khop hash cho
    `lj_cluster_rho_r_0.dat`
  - regression all-atom `test/shared_sparse_graph_dccm_lspmd_post2` van khop
    hash voi `test/shared_sparse_graph_dccm_lspmd_post1` cho `lspmd_*.csv`
  - residual rerun-variance cua `lj_cluster_rho_r_0.dat` da duoc dong

- [done] Phase 2 / 4.2 `Shared sparse graph backend` slice 4 trong
  `dccm_network.cpp`
  - Brandes weighted BC da chuyen `partial` sang flat/aligned row-major
    scratch thay cho `vector<vector<double>>`
  - predecessor scratch da reserve theo degree cua sparse graph de giam
    heap churn trong loop theo source
  - regression all-atom `test/shared_sparse_graph_dccm_lspmd_post3` van khop
    hash voi `test/shared_sparse_graph_dccm_lspmd_post2` cho `lspmd_edges.csv`,
    `lspmd_degree_centrality.csv`, va `lspmd_betweenness_centrality.csv`

- [done] Phase 2 / 4.3 `Shared neighborhood backend` slice 1 trong
  `titan_neighbor_grid.h`, `titan_lj.cpp`, `titan_llps.cpp`, va
  `dccm_network.cpp`
  - da them `HalfShellGrid` dung chung cho half-shell cell-list traversal
  - LJ `build_contact_graph()`, LLPS `build_contacts()`, va DCCM
    `build_lspmd_candidate_pairs()` da dung chung backend nay
  - helper da phan biet ro periodic / non-periodic traversal, tranh wrapped-cell
    false positives ngoai PBC
  - repo-local smoke `test/shared_neighborhood_lj_smoke.exe` va
    `test/llps_neighbor_grid_smoke.exe` da pass
  - rerun LLPS `test/shared_neighbor_llps_post2` va
    `test/shared_neighbor_llps_post3` deterministic; hai file
    `example/5/llps_contact_map.dat` va `example/5/llps_density_rna.dx`
    khong con phu hop de dung lam strict hash gate cho slice nay

Con lam duoc ngay:

- [done] Phase 3 / 5.3 `DCCM / allostery pruning` slice 1
  - allostery da co residue-distance pruning opt-in qua
    `-allostery-cut` + `-allostery-cut-mode`
  - transfer entropy da co the tat qua `-allostery-no-te`, giu file
    `network_te_edges.csv` o dang header-only de khong vo workflow cu
  - regression `test/dccm_allostery_pruning_smoke` giu hash-identical voi
    `test/dccm_allostery_smoke`, run `test/dccm_allostery_pruned_cut12`
    xac nhan pruning giam pairs xuong `9342 / 81810`, va
    `test/allostery_pruning_lspmd_guard_post2` giu hash-identical cho
    `lspmd_*.csv`
  - invalid input `test/allostery_pruning_invalid_cut` da fail-fast dung thong
    diep `Allostery: -allostery-cut must be >= 0.`
- [done] Phase 3 / 5.3 `DCCM / allostery pruning` slice 2
  - da them ranking-aware `frame-voted kNN` pruning va union an toan voi
    cutoff-pruning trong `dccm_network.cpp`
  - `test/allostery_pruning_knn_only_post1..post3` cho thay pruning
    deterministic (`913 / 81810` candidate pairs, `network_*.csv`
    hash-identical qua rerun)
  - `test/allostery_pruning_cut_knn_union_post1` giu
    `9342 / 81810` candidate pairs, khong nho hon cutoff-only baseline
  - `test/allostery_pruning_invalid_knn` va
    `test/allostery_pruning_invalid_cut_mode` deu fail-fast dung thong diep
- [done] Phase 3 / 5.3 `DCCM / allostery pruning` slice 3
  - da them `-allostery-knn-prior=<frame|contact>` va
    `-allostery-knn-contact-cut=<A>` de ho tro true contact-frequency prior
    tren anchor `Cbeta/N`
  - da them `-allostery-prune-window=<avg|first>`; `-allostery-cut-mode`
    duoc giu lam alias tuong thich nguoc
  - `test/allostery_knn_contact_slice3_post1..post2` cho thay
    `999 / 81810` candidate pairs va output allostery deterministic qua rerun
  - `test/allostery_invalid_knn_prior_slice3` va
    `test/allostery_invalid_contact_cut_slice3` deu fail-fast dung thong diep
- [done] Phase 2 / 4.2 `Shared sparse graph backend` slice 3
- da quyet toan hash gate cho `lj_cluster_rho_r_0.dat`
  - da ra not Brandes scratch trong `dccm_network.cpp`; phan con lai neu lam
    tiep se nghieng sang pruning / neighborhood reuse o muc thuat toan cao hon

================================================================================
2. Phase 0: Measurement + Guard Rails
================================================================================

Muc tieu:

- tao baseline hieu nang / RAM / numerical stability de cac phase sau co the do

Cong viec:

- them benchmark nho cho `energy_intrin.cpp`, `chunk_streamer.cpp`,
  `titan_pca.cpp`, `titan_tpt.cpp`, `dccm_network.cpp`, `titan_gist.cpp`
- chuan hoa 3 kich thuoc dataset test:
  - small: debug / CI
  - medium: workstation thong thuong
  - large: stress test cho memory va scaling
- them regression check cho:
  - NaN / Inf
  - zero distance
  - zero row sum
  - ill-conditioned solve
  - salt concentration <= 0 trong path screened electrostatics
- tach build profile:
  - portable: khong phu thuoc `-march=native`, giu explicit `AVX2 + FMA`
  - native: giu toi uu local
- quy uoc regression / benchmark:
  - chay tu thu muc `test/<run_name>` rieng
  - khong ghi output vao thu muc goc chua source code

Deliverables:

- benchmark harness co output wall time + peak RAM
- bo input regression nho de test lai sau moi refactor
- tai lieu build profile ro rang

Status:

- [done] Da co baseline benchmark / regression ban dau duoi `test/`
- [done] Co markdown summary cho chunk streamer, screening, PCA, GIST, va suite
  local `small` / `medium` / `large`
- [done] Peak RAM da duoc dua vao benchmark harness PowerShell
- [done] `Makefile` / `compile.bat` da tach ro `native` / `portable`
- [done] `BUILD_PROFILES.md` da tai lieu hoa build profile va benchmark suites
- [note] `portable` hien tai la profile `AVX2 + FMA + generic tune`; runtime
  CPU dispatch / scalar-portable van thuoc Phase 4

Files / modules lien quan:

- `Makefile`
- `energy_intrin.cpp`
- `chunk_streamer.cpp`
- `titan_pca.cpp`
- `titan_tpt.cpp`
- `dccm_network.cpp`
- `titan_gist.cpp`

Success criteria:

- co baseline truoc khi sua code nang
- moi thay doi o phase sau deu co the doi chieu bang so

================================================================================
3. Phase 1: Quick Wins (1-2 tuan)
================================================================================

Muc tieu:

- giam allocation, copy, va overhead de lay loi ich nhanh ma it rui ro

3.1 Chunk streamer

- bo pattern copy `init_pro.m_mol[0]` / `init_lig.m_mol[0]` cho moi frame roi moi
  overwrite coordinates
- tach frame template immutable khoi coordinate storage
- neu khoi luong cong viec moi chunk nho, tranh respawn worker qua thuong xuyen

Expected impact:

- giam memory traffic
- giam latency preprocessing cho trajectory dai

Status:

- [done] Hoan thanh 2026-03-28 trong `chunk_streamer.cpp`
- PDB split path da chuyen sang frame init khong deep-copy template molecule
- regression smoke da chay tu `test/` va khong thay doi output quan sat duoc

Target:

- `chunk_streamer.cpp`

3.2 PCA scratch reuse

- dua `dev` thanh scratch buffer tai su dung thay vi allocate moi frame
- xem lai loops de giu contiguous access va tranh branch khong can thiet

Status:

- [done] Hoan thanh 2026-03-28 trong `titan_pca.cpp`
- da tai su dung ca `frame_mat` va `dev` ma khong doi semantics weighted /
  unweighted PCA
- regression `example/3` da khop output hash-identical

Expected impact:

- giam allocator pressure
- cai thien throughput o trajectory dai

Target:

- `titan_pca.cpp`

3.3 GIST reduction

- thay full per-thread grid clone bang tiled reduction hoac blocked merge
- neu can, giam bo nho thread-local theo slab / tile thay vi toan bo grid

Status:

- [done] Hoan thanh 2026-03-28 trong `titan_gist.cpp`
- da thay full-grid clone bang sparse per-thread accumulator + touched-voxel
  merge
- da doi chieu `nt=1` / `nt=4` ma khong doi output khoa hoc

Expected impact:

- giam peak RAM khi tang so thread
- giam serial tail reduction

Target:

- `titan_gist.cpp`

3.4 Defensive numerical guards

- dat guard sat kernel cho screened electrostatics
- thong nhat xu ly input bat thuong thanh error / warning ro rang

Status:

- [done] Hoan thanh 2026-03-28
- da them gate input trong `arguments.cpp`
- da them near-zero / invalid distance guard trong scalar va AVX screening path
- regression smoke invalid-input da fail-fast ro rang tu `test/`

Expected impact:

- tranh path rare bi NaN / Inf lan rong

Target:

- `energy.cpp`
- `energy_intrin.cpp`

Exit criteria Phase 1:

- khong doi output khoa hoc ngoai sai so cho phep
- benchmark medium giam thoi gian va peak RAM o it nhat 2-3 hotspot

================================================================================
4. Phase 2: Data Layout + Shared Backend Refactor (2-4 tuan)
================================================================================

Muc tieu:

- loai bo debt kien truc can tro moi toi uu hoa lon hon

4.1 Flat trajectory backend

- thay cac `vector<vector<double>>` trong DCCM / allostery bang flat row-major
  buffers
- chuan hoa layout cho per-frame va per-atom access pattern
- uu tien alignment 32-byte de giu duoc AVX2 path sach

Status:

- [done] Slice 1 da vao `dccm_network.cpp`
- `calculate_and_export()` da dung flat/aligned row-major cho covariance path
- `calculate_lspmd()` da dung flat row-major cho transpose buffers va bo
  allocation scratch theo tung pair
- [done] Slice 2 da flatten `calculate_allostery_information_optimized()`
  sang flat/aligned row-major cho projection bins
- [done] Slice 3 da gom legacy allostery vao shared flat kernel va adaptive
  bin storage
- [done] Slice 4 da bo `get_ca_traj_aos()` va noi PCA thang vao `traj_x/y/z`
- phase 4.1 coi nhu da xong cho quick-win flat backend

Target:

- `dccm_network.cpp`
- cac helper trajectory / network lien quan

4.2 Shared sparse graph backend

- tao backend CSR / COO dung chung cho:
  - DCCM network
  - LJ contact graph
  - TPT flux graph
  - cac module network trong tuong lai
- edge pruning nen xay ra som thay vi tao full dense roi moi nen

Status:

- [done] Slice 1 da them `titan_sparse_graph.h` voi `CSRTopology` va
  `WeightedCSRGraph`
- [done] TPT da migrate sang backend chung ma van giu hash-identical cho
  regression TPT core
- [done] LJ da migrate phan adjacency/unwrapping sang backend chung; output
  `lj_features_0.dat`, `lj_summary_0.dat`, `lj_timeseries_0.dat`, `msm_*.dat`,
  va `tica_*` khop voi baseline hien co
- [done] Slice 2 da migrate DCCM weighted graph / Brandes path sang backend
  chung; `lspmd_*.csv` khop hash voi baseline
- [done] Slice 3 da fix race trong `titan_lj.cpp` va lam on dinh hash gate cho
  `lj_cluster_rho_r_0.dat`
- [done] Slice 4 da ra not Brandes scratch / reduction trong
  `dccm_network.cpp` ma van giu hash gate `lspmd_*.csv`
- [done] `shared neighborhood backend` da co slice 1; neu lam tiep thi la mo
  rong reuse sang them module contact / structure khac

Target:

- module moi: `titan_sparse_graph.*` hoac backend tuong duong

4.3 Shared neighborhood backend

- dua cell-list / neighbor search thanh backend dung chung thay vi de logic nay
  phan tan
- uu tien tai su dung pattern dang tot trong `titan_lj.cpp`

Status:

- [done] Slice 1 da them `titan_neighbor_grid.h` voi `HalfShellGrid`
- [done] LJ `build_contact_graph()` da migrate sang backend chung va van giu
  semantics MIC / padded-bbox theo mode periodic
- [done] LLPS `build_contacts()` da migrate sang backend chung
- [done] DCCM `build_lspmd_candidate_pairs()` da dung chung half-shell grid cho
  pair pruning
- [done] helper da co co `periodic` mode rieng de tranh wrap nham o non-PBC
- [done] Repo-local LJ smoke `test/shared_neighborhood_lj_smoke.exe` da
  validate half-shell grid + MIC / padded-bbox semantics ma khong can bat ky
  `D:` data nao
- [done] Repo-local LLPS smoke `test/llps_neighbor_grid_smoke.exe` da pass;
  rerun `test/shared_neighbor_llps_post2` va `test/shared_neighbor_llps_post3`
  khop hash voi nhau
- [note] `example/5/llps_contact_map.dat` va `example/5/llps_density_rna.dx`
  khac voi rerun hien tai nen can refresh baseline truoc khi dung lam hash gate

Target:

- `titan_lj.cpp`
- future interaction / contact / structure modules

Exit criteria Phase 2:

- DCCM / network path khong con phu thuoc chinh vao `vector<vector<T>>`
- co backend sparse graph tai su dung duoc o it nhat 2 module

================================================================================
5. Phase 3: Algorithmic Scaling (3-6 tuan)
================================================================================

Muc tieu:

- sua cac diem nghen dense bang thay doi thuat toan, khong chi micro-opt

5.1 PCA scalable

- bo sung mode truncated eigen / randomized SVD khi chi can top components
- can nhac block covariance accumulation hoac direct low-rank solver

Status:

- [done] Slice 1: truncated spectral solver cho PCA K lon
  - them `-pca-solver=<full|partial>`, `-pca-topk=<N>`, va
    `-pca-max-k-full=<N>` trong `arguments.h`, `parse_arguments.cpp`,
    `arguments.cpp`, va `help.cpp`
  - `titan_pca.cpp` da them duong `partial(subspace)` theo huong subspace-iteration
    + Ritz spectrum, co fallback ve `full` neu khong hoi tu
  - runtime log da co `PCA Solver`; output `pca_fel.csv` / `pca_fes_2d.dat`
    van giu format cu

Target:

- `titan_pca.cpp`

5.2 Sparse TPT

- [done] bo persistent dense `net_flux_`; flux graph duoc nen truc tiep vao
  `edges_`
- [done] dominant path extraction da chay tren edge list / sparse residual graph
- [done] chuyen committor solve sang sparse linear system
- [done] bo persistent dense `reverse_tpm_`; backward solve build truc tiep tu
  `tpm_` + `stationary_dist_`
- [done] tranh recompute gross flux / rescan dense pair trong flux compression
- [done] bo persistent dense `tpm_`; transition graph da duoc normalize sang CSR
  forward / reverse va flux worker da duyet adjacency sparse
- [done] toi uu reverse-edge lookup / adjacency queries trong TPT residual path,
  tan dung backend CSR chung da co thay vi mo rong them dense helper rieng
  - `titan_tpt.cpp` da them `reverse_edge_lookup_` de tra cuu canh nguoc O(1)
    theo edge-index, bo scan adjacency lap lai
  - duong flux worker giu TitanExecutor parallel path, khong doi output file TPT

Target:

- `titan_tpt.cpp`

5.3 DCCM / allostery pruning

- cho phep residue-distance pruning, kNN pruning, hoac contact-prior pruning
- tach transfer entropy thanh stage tuy chon vi day la path dat nhat
- nghien cuu estimator kNN / KDE cho mot so metric thay vi histogram moi cap

Status:

- [done] Slice 1: pruning opt-in + TE optional
  - `dccm_network.cpp` da dung `build_lspmd_candidate_pairs()` de sinh
    candidate-pair list cho allostery khi bat `-allostery-cut`
  - `engine_postproc.cpp` da log `Allostery Pair Filter`,
    `Allostery Candidate Pairs`, va `Allostery TE Stage`
  - default path khong thay doi output; regression
    `test/dccm_allostery_pruning_smoke` khop hash voi
    `test/dccm_allostery_smoke`
  - LSP-MD khong bi anh huong; regression
    `test/allostery_pruning_lspmd_guard_post2` khop hash voi
    `test/shared_sparse_graph_dccm_lspmd_post1`
  - validation input `test/allostery_pruning_invalid_cut` da fail-fast dung
    thong diep neu `-allostery-cut < 0`

- [done] Slice 2: ranking-aware pruning
  - `dccm_network.cpp` da them ranking-aware `kNN` pruning theo frame-voted
    nearest-neighbor votes + `min_dist2` tie-break, va hop nhat an toan voi
    cutoff-pruning khi cung bat
  - `engine_postproc.cpp` / `help.cpp` da doi log va CLI wording sang
    `frame-voted kNN` de phan anh dung heuristic hien tai
  - `parse_arguments.cpp` da fail-fast neu `-allostery-cut-mode` khong nam
    trong `avg|first`
  - regression `test/allostery_pruning_knn_only_post1`,
    `test/allostery_pruning_knn_only_post2`, va
    `test/allostery_pruning_knn_only_post3` cho thay `913 / 81810` candidate
    pairs, output `network_*.csv` hash-identical qua rerun, va TE header-only
    dung nhu mong doi khi bat `-allostery-no-te`
  - regression `test/allostery_pruning_cut_knn_union_post1` giu
    `9342 / 81810` candidate pairs, khong nho hon cutoff-only baseline
    `test/dccm_allostery_pruned_cut12`
  - validation input `test/allostery_pruning_invalid_knn` va
    `test/allostery_pruning_invalid_cut_mode` deu fail-fast dung thong diep
  - default path va LSP-MD guard van on dinh: `test/dccm_allostery_pruning_smoke`
    khop hash voi `test/dccm_allostery_smoke`, va
    `test/allostery_pruning_lspmd_guard_post3` khop hash voi
    `test/shared_sparse_graph_dccm_lspmd_post1`

- [done] Slice 3: estimator / prior refinement
  - da them `-allostery-knn-prior=<frame|contact>` de tach `frame-voted`
    nearest-neighbor heuristic khoi true contact-frequency prior
  - da them `-allostery-knn-contact-cut=<A>`; contact prior hien xep hang theo
    tan suat contact tren anchor `Cbeta/N`, tie-break bang `min_dist2`, va
    fallback ve `Calpha` neu cache anchor chua san sang
  - da them `-allostery-prune-window=<avg|first>`; `-allostery-cut-mode`
    duoc giu lam alias tuong thich nguoc de tach ro y nghia `pruning window`
    khoi distance cutoff
  - regression `test/allostery_default_slice3_post1..post2`,
    `test/allostery_knn_frame_slice3_post1..post2`, va
    `test/allostery_knn_contact_slice3_post1..post2` deu deterministic;
    `frame-voted` giu `913 / 81810` candidate pairs, con `contact-frequency`
    cho `999 / 81810`
  - `test/allostery_pruning_lspmd_guard_post4` van khop hash voi
    `example/2/lspmd_*.csv`
  - validation input `test/allostery_invalid_knn_prior_slice3` va
    `test/allostery_invalid_contact_cut_slice3` deu fail-fast dung thong diep

- [done] Slice 4: estimator upgrades
  - da them `-allostery-estimator=<hist|knn|kde>` va
    `-allostery-kde-bandwidth=<val>` trong parser/validation/help
  - allostery kernel trong `dccm_network.cpp` da co estimator hook cho MI/GC:
    `hist` (legacy), `kde` (smooth histogram), va `knn` (distance-weighted adjustment)
  - `engine_postproc.cpp` da log `Allostery Estimator` de theo doi mode runtime,
    giu output `raw_mi_edges.csv`, `network_gc_edges.csv`, `network_te_edges.csv`

Target:

- `dccm_network.cpp`

5.4 MSM spectral cost


Status:

- [done] ra soat ITS scan va eigensolve per lag trong `titan_msm.cpp`
  - tach helper dung chung de build TPM theo lag (`count -> symmetrize -> normalize`)
  - ITS worker chay song song qua `TitanExecutor`, precompute boundary mot lan,
    va giam luong state tam giu theo lag
  - them log profiler `ITS Matrix+Spectral`, `ITS Total`, va `ITS Solver Usage`

- [done] bo sung che do spectral cho K lon (giu tuong thich nguoc)
  - them `-msm-its-solver=<full|schur|partial>`
  - them `-msm-its-topk=<N>` va `-msm-its-max-k-full=<N>`
  - co fallback ve full solver neu solver duoc chon khong hoi tu
  - `help.cpp`, `parse_arguments.cpp`, `arguments.cpp` da cap nhat validation

- [done] giu nguyen output ITS hien tai
  - van ghi `msm_its_scan.dat` va `msm_its_tpm/tpm_lag*.dat`
  - khong doi ten file/format cot de tranh vo workflow post-process cu

- [done] toi uu cho truong hop K rat lon
  - che do `partial` da doi sang subspace-iteration + Ritz spectrum (khong goi
    full eigensolver tren ma tran `K x K`), giam them chi phi spectral khi
    `K` lon trong ITS scan
Target:

- `titan_msm.cpp`

Exit criteria Phase 3:

- medium/large dataset scale tot hon ro rang
- memory growth giam theo dung muc tieu doi voi PCA / TPT / DCCM

================================================================================
6. Phase 4: Toolchain + Runtime Portability
================================================================================

Muc tieu:

- tranh phu thuoc qua muc vao may build cu the

Cong viec:

- giu `-march=native` cho local benchmarking
- them build portable hoac runtime CPU dispatch cho:
  - scalar
  - AVX2
  - AVX2 + FMA
- gioi han code moi trong subset C++23 ma `g++13/libstdc++13` ho tro on dinh
- khong xem `g++13` la target toolchain bat buoc cho local benchmark; no la moc
  tuong thich de tranh dung nham tinh nang thu vien / ngon ngu qua moi
- [done] Slice 1: harden build/profile plumbing ma khong doi logic tinh toan
  - `Makefile` dung object/dependency tach theo `build/native/` va
    `build/portable/` de tranh stale object contamination
  - `titan_build_info.h` + startup log in ro `Build Profile`, `C++23 Floor`,
    `Compiled SIMD`, `Runtime SIMD`, `Runtime Supports Build SIMD`,
    va `std::jthread`
  - native/portable regression trong `test/phase4_*` hash-identical cho
    allostery (`raw_mi_edges.csv`, `network_gc_edges.csv`, `network_te_edges.csv`)
    va LSP-MD (`lspmd_edges.csv`, `lspmd_degree_centrality.csv`,
    `lspmd_betweenness_centrality.csv`)
- [done] Slice 2: centralize runtime dispatch cho cac kernel nong nhat,
  giu scalar fallback khong doi logic va so sanh `native` vs `portable`
  - them `titan_simd_dispatch.h` de gom shared helper cho reduction,
    dot-product, squared-distance, add-row, va 4-distance hotspots
  - wire runtime dispatch/log metadata qua `titan_build_info.h` va startup log
  - native/portable regression duoi `test/phase4_*` duoc dung de xac nhan
    output khoa hoc khong doi
  - portable profile van la `AVX2 + FMA + generic tuning`; true scalar-baseline
    portability can cleanup AVX assumption rong hon va la follow-up rieng
- [done] Slice 3: migrate classic IFP / H-bond AVX screening sang shared
  runtime-dispatch helper
  - them shared `point-vs-4` cutoff helpers trong `titan_simd_dispatch.h`
    cho contiguous va gathered coordinate paths
  - `titan_ifp.cpp` bo AVX block rieng va dung helper shared cho ligand atom
    screening
  - `h_bond.cpp` bo AVX2 prefilter rieng va dung helper gathered shared cho
    hydrogen-bond group pruning
  - regression duoi `test/phase4_slice3_ifp_*` hash-identical cho
    `ifp_occupancy.csv`, `out_hydrogen_bonds*.dat`,
    `out_list_res_HB_*.dat`, `stat_average_HB.dat`, va
    `stat_hydrogenbonds.dat`
  - ket qua khop ca theo truc `baseline native truoc-sua -> native sau-sua`
    va `native -> portable`
- [done] Slice 4: migrate LLPS autocorrelation popcount sang shared
  runtime-dispatch helper
  - `titan_llps.cpp` bo local AVX2 popcount block trong `finalize()` va goi
    helper XOR-overlap popcount chung tu `titan_simd_dispatch.h`
  - regression trong `test/phase4_slice4_llps_baseline_native`,
    `test/phase4_slice4_llps_native`, va
    `test/phase4_slice4_llps_portable` hash-identical cho
    `llps_summary.dat`, `llps_timeseries.dat`, `llps_degree_dist.dat`,
    `llps_dbscan_timeseries.csv`, `llps_dbscan_droplets.pdb`,
    `llps_contact_map.dat`, `llps_density_protein.dx`,
    `llps_density_rna.dx`, va `polymer_scaling.dat`
  - `llps_autocorr.dat` tiep tuc vang mat o ca 3 run vi fixture chua co du
    post-percolation frames cho `C(t)`, nen absence nay duoc xem la gate on
    dinh cua slice
- [done] Slice 5: migrate DSSP CA-distance prefilter sang shared
  runtime-dispatch helper
  - `titan_dssp.cpp` bo local AVX distance-register block trong
    `compute_hbond_topology()` va dung
    `titan_simd_dispatch::point_vs_four_cutoff_mask(...)`
  - logic hbond energy `hbond_energy_scalar()` va residue exclusion
    (`donor_acceptor_excluded`) giu nguyen; chi doi lop prefilter 4-lane
  - regression trong `test/phase4_slice5_dssp_baseline_native`,
    `test/phase4_slice5_dssp_native`, va
    `test/phase4_slice5_dssp_portable` hash-identical cho
    `dssp_mol1.dat` va `dssp_mol1_content.dat`
- [done] Slice 6: migrate GBSA Shrake-Rupley SASA burial filter sang shared
  runtime-dispatch helper
  - `gbsa_intrin.cpp` bo local AVX2 burial-mask block trong
    `compute_sasa_shrake_rupley()` va dung
    `titan_simd_dispatch::point_vs_four_cutoff_mask(...)`
  - logic khoa hoc giu nguyen: burial test van la strict `< Rj^2`, SASA
    accumulation va downstream GBSA math khong doi
  - regression trong `test/phase4_slice6_gbsa_baseline_native`,
    `test/phase4_slice6_gbsa_native`, va
    `test/phase4_slice6_gbsa_portable` hash-identical cho
    `gbsa_summary.dat`, `gbsa_timeseries_avg.dat`,
    `traj_0_gbsa_snapshots.dat`, `gbsa_decomposition_ligand.dat`, va
    `gbsa_decomposition_receptor.dat`
- [done] Slice 7: migrate GBSA effective Born-radii helper sang shared
  runtime-dispatch helper
  - `gbsa_intrin.cpp` bo local AVX horizontal-sum + register loop trong
    `cal_effective_born_radii_avx2()` va dung
    `titan_simd_dispatch::gbsa_inverse_r4_volume_sum(...)`
  - `titan_simd_dispatch.h` da them helper Born-radii dung chung voi 3 mode:
    scalar, AVX2, va AVX2+FMA; AVX2/FMA chi thay doi cach thuc thi, khong doi
    cong thuc `rho^3 / 3 / r^4`, gate `r2 > 1e-6`, hay downstream Born-radius
    formulas
  - regression trong `test/phase4_slice7_gbsa_baseline_native`,
    `test/phase4_slice7_gbsa_native`, va
    `test/phase4_slice7_gbsa_portable` hash-identical cho
    `gbsa_summary.dat`, `gbsa_timeseries_avg.dat`,
    `traj_0_gbsa_snapshots.dat`, `gbsa_decomposition_ligand.dat`, va
    `gbsa_decomposition_receptor.dat`
- [done] Slice 8: migrate GBSA LJ energy helper sang shared
  runtime-dispatch helper
  - `gbsa_intrin.cpp` bo local AVX register math lap lai trong
    `cal_energy_gbsa_avx2()` cho duong OPLS/AMBER-style LJ mixing va dung
    `titan_simd_dispatch::gbsa_lj_mixed_sigma_epsilon_r6_fourlanes(...)`
  - `titan_simd_dispatch.h` da them helper LJ dung chung voi 2 mode:
    scalar va AVX2; AVX2 chi thay doi cach thuc thi, khong doi cong thuc
    `4 * epsilon * (sr12 - sr6)`, sigma mixing, epsilon mixing, hay
    `r6` floor `1e-12`
  - regression trong `test/phase4_slice8_gbsa_baseline_native`,
    `test/phase4_slice8_gbsa_native`, va
    `test/phase4_slice8_gbsa_portable` hash-identical cho
    `gbsa_summary.dat`, `gbsa_timeseries_avg.dat`,
    `traj_0_gbsa_snapshots.dat`, `gbsa_decomposition_ligand.dat`, va
    `gbsa_decomposition_receptor.dat`
- [done] Slice 9: migrate AMBER energy helper path sang shared
  runtime-dispatch helper
  - `energy_intrin.cpp` va `energy_intrin_nofilter.cpp` bo local AVX register
    math lap lai trong filtered + no-filter AMBER LJ helper paths va dung
    `titan_simd_dispatch::lj_mixed_sigma_epsilon_distance_fourlanes(...)`
  - `titan_simd_dispatch.h` da them helper LJ dung chung voi 2 mode:
    scalar va AVX2; AVX2 chi thay doi cach thuc thi, khong doi sigma mixing
    arithmetic, epsilon mixing `sqrt(16 * eps_i * eps_j)`, hay screening
    semantics voi `safe_d`
  - regression trong `test/phase4_slice9_energy_baseline_native`,
    `test/phase4_slice9_energy_native`,
    `test/phase4_slice9_energy_portable`,
    `test/phase4_slice9_energy_screening_baseline_native`,
    `test/phase4_slice9_energy_screening_native`, va
    `test/phase4_slice9_energy_screening_portable` deu hash-identical cho
    toan bo output files tru `arg.in` va `log.txt`
- [done] Slice 10: migrate AMBER/CHARMM force helper path sang shared
  runtime-dispatch helper
  - `force_intrin.cpp` bo local AVX register math lap lai trong
    `force_energy_add_intrin(...)` cho nhanh AMBER/CHARMM LJ force helper va
    dung `titan_simd_dispatch::lj_mixed_sigma_epsilon_distance_c6_fourlanes(...)`
  - `titan_simd_dispatch.h` da them helper LJ dung chung tra ve ca `c6` va
    energy voi 2 mode: scalar va AVX2; AVX2 chi thay doi cach thuc thi,
    khong doi sigma mixing arithmetic, epsilon mixing `sqrt(eps_i * eps_j)`,
    hay `c6`-based self-interaction guard
  - regression trong `test/phase4_slice10_force_baseline_native`,
    `test/phase4_slice10_force_native`, va
    `test/phase4_slice10_force_portable` deu hash-identical cho toan bo
    output files tru `arg.in`, `log.txt`, va `run_stdout.txt`
  - `titan_rf_model.cpp` cung da duoc khoa compile `-O0` bang `Makefile`
    cho moi profile thay vi dua vao source pragma
- [done] Slice 11: migrate legacy OPLS force helper path sang shared
  runtime-dispatch helper
  - `force_intrin.cpp` bo local AVX register math lap lai trong nhanh OPLS
    cua `force_energy_add_intrin(...)` va dung
    `titan_simd_dispatch::lj_sigma_product_distance_c6_fourlanes(...)`
  - `titan_simd_dispatch.h` da them helper LJ dung chung tra ve ca `c6` va
    energy voi 2 mode: scalar va AVX2; helper nay giu nguyen legacy OPLS AVX
    lane math order dua tren sigma-product lanes, khong doi `c6`-based
    self-interaction guard
  - regression trong `test/phase4_slice11_force_opls_baseline_native`,
    `test/phase4_slice11_force_opls_native`, va
    `test/phase4_slice11_force_opls_portable` deu hash-identical cho toan bo
    output files tru `arg.in`, `log.txt`, va `run_stdout.txt`
- [done] Slice 12: dong audit AVX local con sot va chot Phase 4
  - `dccm_network.cpp` bo helper `compute_4dist_avx` local va dead SIMD
    fallback code; `calculate_lspmd()` dung truc tiep
    `titan_simd_dispatch::compute_4dist(...)` va shared dot-product helper
  - `MOL.h` bo local `_mm256` loop trong `reduce_forces()` va dung
    `titan_simd_dispatch::add_row_inplace(...)` cho phan prefix vectorizable
  - `find_pull_direction.cpp` bo local AVX translate loop va dung
    `titan_simd_dispatch::translate_xyz(...)`
  - regression trong `test/phase4_slice12_native_lspmd_example2`,
    `test/phase4_slice12_portable_lspmd_example2`,
    `test/phase4_slice12_force_opls_native`,
    `test/phase4_slice12_force_opls_portable`,
    `test/phase4_slice12_pull_native`,
    `test/phase4_slice12_pull_native_rerun`, va
    `test/phase4_slice12_pull_portable`
  - hash gates:
    `lspmd_edges.csv`, `lspmd_degree_centrality.csv`,
    `lspmd_betweenness_centrality.csv`; OPLS `cal_force.txt`,
    `cal_force.xvg`, `total_force.dat`, `total_energy.dat`,
    `electrostatic_force.dat`, `vanderwaals_force.dat`; pull-direction
    `DE_out_vector.dat`, `DE_vector.pdb`, `best_conf.pdb`, va aligned snapshot
  - ket qua sau audit:
    AVX local con lai trong `energy_intrin*.cpp`, `force_intrin.cpp`,
    `gbsa_intrin.cpp`, `titan_gist.cpp`, va `cal_sub_groups.cpp` duoc giu lai
    co chu dich vi do la full kernels hoac microkernels dac thu; `avx_math.h`,
    `BoxMatrix.h`, `functions.h`, va `simd_compat.h` la support headers
  - `DE_log.dat` / `DE_report.dat` cua pull-direction khong duoc dung lam hash
    gate khoa hoc vi co chen thong tin report/runtime theo profile

Target:

- `Makefile`
- cac kernel compute chinh

Exit criteria:

- Titan build duoc o che do portable
- benchmark native van la baseline nhanh nhat tren may local

================================================================================
7. Phase 5: Module moi nen them sau khi co backend tot
================================================================================

Cac module duoi day kha thi va co gia tri cao, nhung nen lam sau khi Phase 2-3
da tao duoc nen tang sparse / flat-buffer tot.

7.1 `titan_mbar`

Status:

- `Shipped core + standalone outputs`
  - `5.1a core + frame weights + thermo scan`
  - `5.1b weighted TICA FES`
  - `5.1c weighted PCA FES`
  - `5.1d NPT enthalpy-like observable path`
- `Shipped equilibrium integrations`
  - `5.1e per-frame temperature / demux-aware MBAR`
  - current regression artifacts `per_frame`, `demux`, va `NPT` da ghi
    `MBAR Converged: yes` trong native/portable outputs hien tai
  - `5.1f MBAR hardening + fail-closed gate`
  - `5.1g shared equilibrium weight-source abstraction`
  - `5.1h optional default weighted TICA/PCA selector`
  - `5.1i MBAR-weighted MSM occupancy summaries`
  - `5.1j TPT context-only overlays from MBAR populations`
  - `5.1k MBAR-weighted UMAP postproc (phase 1)`
  - `5.1l MBAR-weighted DCCM/allostery/LSP-MD context (phase 1)`
- `Production boundaries`
  - standalone core + optional equilibrium postproc source; khong overwrite hay
    thay the DHAM
  - interface production cho integration la
    `-eq-weight-source <auto|uniform|dham|mbar>`
  - `auto` van giu DHAM la default neu co weights hop le; `mbar` chi kich hoat
    khi user opt-in ro rang
  - khong doi kinetics hay thay `TPM/pi` cua `MSM/TPT/DHAM`; MBAR hien chua la
    kinetic reweighting engine production

Pham vi v1 da ship:

- MBAR cho `temperature ladder` REMD/HREMD
- xuat `mbar_frame_weights.dat`, `mbar_summary.dat`,
  `mbar_window_audit.dat`
- neu `-mbar-scan-points > 1` thi xuat `mbar_cv.dat`
- reweight FES doc lap sang `target temperature`:
  `mbar_tica_fes_2d.dat`, `mbar_tica_cluster_representatives.csv`,
  `mbar_pca_fes_2d.dat`, `mbar_pca_cluster_representatives.csv`
- co ca hai mode thermodynamic-state assignment:
  `fixed` theo ladder-per-file va `per_frame` theo demux hoac nearest-ladder
  mapping tu XVG temperature column
- `nvt` dung observable energy; `npt` dung observable enthalpy-like tu XVG

Ung dung trong Titan:

- dien giai ensemble muc tieu cho REMD/HREMD o nhiet do dich
- weighted basin ranking va representative-frame selection tren TICA/PCA
- heat-capacity / melting / phase-transition diagnostics cho protein va LJ
- xu ly nhanh workflow reweighting nghiem tuc hon so voi histogram FES thuan

Tich hop vao Titan hien tai:

- tai su dung shape input cua DHAM:
  `energy_list.txt`, va `temps.txt` chi con la optional ladder override
- `-mbar-domain=tica` di qua metadata cua TICA/MSM
- `-mbar-domain=pca` di qua metadata PCA va `pca_fel.csv`
- `-mbar-domain=both` xuat ca hai khong doi logic output cu
- `window-match` va `dt-energy` dung cung semantics voi DHAM
- MBAR solver duoc tach khoi standalone export de co the cap frame weights som
  hon cho tang postproc equilibrium
- `-eq-weight-source auto` van giu `DHAM -> PCA/TICA` weighted path lam mac dinh
- `-eq-weight-source mbar` la interface production duy nhat de opt-in MBAR vao
  default `TICA/PCA` FES, `msm_equilibrium_summary.dat`, va
  `tpt_equilibrium_overlay.dat`
- consumer equilibrium hien tai cua MBAR la `TICA/PCA/MSM/TPT`; tat ca deu di
  theo policy explicit opt-in va khong co implicit promotion tu `auto` sang
  `mbar`

Luu y thiet ke:

- v1 chi ho tro temp-ladder, chua lam generic `u_kn`
- raw GROMACS T-REMD replica outputs voi `T(t)` thay doi theo frame da duoc ho
  tro qua `-mbar-temp-mode per_frame`
- trong `per_frame`, TITAN co the tu infer ladder nhiet do tu replica input;
  neu co `-mbar-demux-map` thi demux la route authoritative de dung state
  assignment va suy ra ladder trung binh theo state
- `-mbar-demux-map` la duong authoritative cho state assignment; neu khong co,
  TITAN fallback sang nearest-ladder mapping tu XVG temperature column
- MBAR trong repo hien tai duoc xem nhu `equilibrium reweighting engine`, khong
  phai kinetic reweighting engine
- `WHAM` khong con la deliverable chinh o Phase 5; neu can thi de lai
  nhu note compatibility ve sau
- NPT pressure-shift MBAR tong quat tu `E + P*V` van la mot extension sau;
  `5.1d` hien tai tap trung vao enthalpy-like observable path on dinh cho
  workflow GROMACS NPT

Delivered hardening + current equilibrium-safe slices:

- [done] `5.1f MBAR hardening + fail-closed gate`
  - downstream `mbar_*` weighted FES da chuyen sang fail-closed theo
    `converged=true`
  - benchmark harness da co `mbar_fixed_tica` va `mbar_per_frame_tica`
  - regression native/portable explicit da cover `fixed`, `per_frame`,
    `demux`, va `NPT`
- [done] `5.1g shared equilibrium weight-source abstraction`
  - da tach abstraction chung `auto | uniform | DHAM | MBAR` o tang
    postproc equilibrium
  - duoc dung lai cho default FES va equilibrium-only summaries; khong di vao
    kinetics
- [done] `5.1h optional default weighted TICA/PCA selector`
  - da co `-eq-weight-source <auto|uniform|dham|mbar>`
  - `auto` giu behavior cu; chi khi opt-in `mbar` moi feed vao
    `tica_fes_2d.dat` / `pca_fes_2d.dat`
- [done] `5.1i MBAR-weighted MSM occupancy summaries`
  - da them `msm_equilibrium_summary.dat` cho microstate/macrostate populations
    tu MBAR weights
  - khong rebuild TPM, khong doi ITS, khong doi committor inputs
- [done] `5.1j TPT context-only overlays from MBAR populations`
  - da them `tpt_equilibrium_overlay.dat` cho basin/source/sink populations va
    per-state equilibrium context
  - `TPT Model Source` va kinetics van giu nguyen; MBAR khong thay `TPM/pi`
- [done] `5.1k MBAR-weighted UMAP postproc (phase 1)`
  - da them `umap_fes_2d.dat`, `umap_cluster_representatives.csv`, va
    `umap_equilibrium_summary.dat` khi `-eq-weight-source mbar`
  - embedding fit van giu unweighted; MBAR chi anh huong equilibrium sidecars
    va khong doi output mac dinh neu user khong opt-in
- [done] `5.1l MBAR-weighted DCCM/allostery/LSP-MD context (phase 1)`
  - da them `dccm_equilibrium_context.dat`,
    `allostery_equilibrium_ranking.csv`, va
    `lspmd_equilibrium_ranking.csv`
  - core correlation / network construction / LSP-MD estimator van giu nguyen;
    MBAR chi them contextual overlay va ranking
- [done] `5.1m Cross-module MBAR validation matrix (phase 1)`
  - benchmark harness da cover regression cho `fixed`, `per_frame`, `demux`,
    `NPT`, selector `PCA`, `UMAP`, `MSM`, `TPT`, `DCCM/allostery`, `LSP-MD`
  - da them guard cases cho `auto`, `uniform`, explicit `dham`, va requested
    `mbar` nhung khong resolve duoc
  - da them cross-case hash guard de khoa invariant `auto == uniform` tren mot
    so output dai dien khi khong co implicit promotion sang `mbar`

Follow-up sau `5.1m` phase 1:

- mo rong validation matrix tren fixture lon hon / data REMD thuc hon neu can
- co the them hash/content gates cho nhieu output hon neu sau nay can freeze
  regression o muc chi tiet hon

Khong plan trong dot nay:

- khong lam `MBAR -> DHAM`
- khong reuse `-eq-weight-source` cho kinetics
- khong xem `MBAR -> MSM/TPT` kinetics nhu production deliverable trong wave nay

MBAR kinetics research track:

- exploratory / experimental only; tach khoi production roadmap equilibrium
- huong nay chi duoc mo bang interface experimental rieng trong tuong lai,
  khong dung chung contract voi `-eq-weight-source`
- `MBAR -> MSM kinetic reweighting` la research-only cho toi khi co validation
  synthetic/analytic ro rang
- `MBAR -> TPT core kinetics` la research-only; khong duoc mo ta nhu thay the
  `TPM/pi` core trong wave hien tai
- khong muc nao trong research track duoc danh dau production-ready neu chua co
  gate validation rieng

Acceptance cho wave tich hop tiep theo:

- roadmap phai giu dung ground truth hien tai:
  - MBAR da co standalone weighted `TICA/PCA` exports
  - default `DHAM -> PCA/TICA` path van la mac dinh khi `-eq-weight-source auto`
  - `-eq-weight-source mbar` chi feed MBAR vao equilibrium FES/summaries
  - `TPT` core van chi nhan built MSM hoac `DHAM TPM/pi`; khong nhan MBAR
- truth sync:
  - current `per_frame`, `demux`, va `NPT` MBAR artifacts phai duoc ghi dung la
    `Converged: yes`
  - roadmap khong con mo ta MBAR la standalone-only
- cho `5.1f`:
  - tat ca case `fixed` native/portable converge
  - case `per_frame/demux/NPT` hoac converge, hoac fail closed va khong xuat
    weighted downstream outputs
  - benchmark harness co it nhat 1 case MBAR `fixed` va 1 case MBAR
    `per_frame`
- cho `5.1h-5.1j`:
  - da giu nguyen hanh vi o default outputs neu user khong chon MBAR explicit
  - cac tich hop MBAR moi chi anh huong den equilibrium summaries/FES, khong doi
    kinetics
  - smoke guard da xac nhan `auto` khong sinh `msm_equilibrium_summary.dat` hay
    `tpt_equilibrium_overlay.dat`
- equilibrium postproc invariants:
  - `TICA/PCA/MSM/TPT` hien da la MBAR consumer o muc equilibrium-only
  - `UMAP` va `DCCM/allostery/LSP-MD` da ship phase 1 voi explicit opt-in
    invariant
  - default outputs duoi `auto` phai giu nguyen hanh vi cu
- cho `5.1k-5.1m`:
  - benchmark harness da co explicit regression cho `UMAP`, `TPT overlay`,
    `DCCM/allostery/LSP-MD`, `auto`, `uniform`, `dham`, va unresolved-`mbar`
  - cross-case hash checks da xac nhan mot nhom output dai dien la identical
    giua `auto` va `uniform` khi khong co `DHAM` de implicit promotion
- research-track invariants:
  - moi huong `MBAR -> MSM/TPT kinetics` phai co validation gate rieng
  - khong muc nao trong research track duoc danh dau `done` hay production-ready
    neu chua co synthetic/analytic validation
- validation matrix:
  - moi consumer moi cua shared equilibrium weights phai co regression rieng
  - production acceptance va research acceptance phai duoc tach ro

7.2 `titan_hbond_kinetics`

- hydrogen bond
- salt bridge
- water bridge
- residence time / intermittent correlation

Status:

- [done] Public `hbond kinetics` module da co trong workspace
  - da them interface `-hbk` voi `kind = hbond|salt|water`, `domain =
    inter|intra|both`, `maxlag`, va cutoff/angle controls; parser accept ca
    `-opt=value` va `-opt value`, trong khi help/fixtures canonicalize ro
    `arg.in =` va `CLI space`
  - analyzer moi la streaming time-correlation core cho `continuous` /
    `intermittent C(t)`, residence statistics, va lifetime histogram; khong
    thay doi contract cua legacy `-hb`, `IFP`, hay `global_water_bridges_*`
  - primitive detection da duoc tach ra helper dung chung de reuse donor /
    acceptor / hydrogen lookup, salt typing, va water-neighbor logic ma khong
    copy lai detector
  - output moi gom `hbond_kinetics_summary.dat`,
    `hbond_kinetics_{hbond,salt,water}_pairs.csv`,
    `hbond_kinetics_{hbond,salt,water}_ct.dat`, va
    `hbond_kinetics_{hbond,salt,water}_lifetimes.dat`
  - smoke C++ `test/titan_hbond_kinetics_smoke.cpp` da pass; fixture local va
    benchmark harness cases cho `-hbk` da co trong workspace hien tai

7.3 `titan_validation`

- Chapman-Kolmogorov test
- VAMP-2 / score-based model selection
- bootstrap uncertainty cho MSM / TPT

Status:

- [done] Public `validation` module da co trong workspace
  - da them interface `-validation` voi CK multiples, VAMP candidate grid,
    `-validation-auto-pick`, va bootstrap seed/sample controls
  - preselect phase da chen vao sau khi TICA solve xong va truoc khi build MSM
    chinh; `auto-pick` co the route nguoc selected `(components, lag, clusters)`
    vao pipeline `MSM -> DHAM -> TPT`
  - postproc validation da xuat `validation_summary.dat`,
    `validation_ck_summary.dat`, `validation_ck_states.dat`,
    `validation_vamp2_ranking.csv`, `validation_selected_model.dat`,
    va bootstrap sidecars cho `MSM/TPT`
  - kinetic `tpt_kinetics.dat` va legacy `tpt_bootstrap_*` da route sang shared
    validation core; DHAM-backed TPT bootstrap da co path thuc thi
  - build `native` va `portable` da pass; smoke C++ `test/titan_validation_smoke.cpp`
    da pass trong workspace hien tai

7.4 `titan_corr_knn`

- generalized correlation
- mutual information
- transfer entropy bang estimator tot hon histogram dense

Status:

- [done] Public `corr-knn` module + exact kNN estimators da co trong workspace
  - da them interface `-corr-knn` voi `k`, `tau`, pair-pruning, `no-te`, va raw
    outputs `corr_knn_summary.dat`, `corr_knn_mi_pairs.csv`,
    `corr_knn_gc_pairs.csv`, `corr_knn_te_pairs.csv`
  - `MI` dung KSG type-1 tren fluctuation vector 3D theo residue node; `TE` dung
    conditional-MI kieu Frenzel-Pompe; `GC` duoc suy ra tu `MI`
  - `-allostery-estimator=knn` da route sang core moi de giu compatibility, trong
    khi `hist` va `kde` van giu path cu
  - smoke C++ va benchmark local native/portable cho `corr_knn`,
    `allostery hist/kde/knn`, va reuse path da pass trong workspace hien tai

7.5 `titan_sq`

- structure factor `S(q)`
- streaming RDF / pair correlation
- phu hop LJ / LLPS / condensed-phase workflows

Status:

- [done] Streaming `RDF/S(q)` analyzer da co trong workspace
  - da them `-sq`, channel `all|trace|backbone|sticker`, partials theo channel,
    auto-PBC policy, va output `sq_summary.dat`, `sq_rdf_*`,
    `sq_structure_factor_*`
  - core tinh toan la streaming pair-correlation; `S(q)` duoc suy ra o cuoi run tu
    `RDF` bang isotropic transform co Lorch window
  - orthorhombic periodic path dung neighbor-grid; triclinic/dodecahedron fallback
    sang dense pair loop + `MIC`; non-PBC chi warn va skip
  - smoke/test fixture/harness local cho `-sq` da pass trong workspace hien tai,
    gom case nho native/portable va case `example/5` cho `backbone/sticker`

7.6 `titan_sparse_graph`

- graph analytics dung chung
- community detection
- edge filtering
- bottleneck ranking

Status:

- [done] Public unified sparse-graph analytics layer da co trong workspace
  - da them interface `-graph`, `-graph-source`, `-graph-community`,
    `-graph-filter-topk`, `-graph-filter-component`, voi canonical syntax
    `arg.in =` / `CLI space`
  - core moi da consume graph noi bo tu `DCCM`, `LSP-MD`, thresholded
    `allostery GC/TE`, va `TPT`; export `graph_summary.dat`,
    `graph_<source>_edges.csv`, `graph_<source>_nodes.csv`,
    `graph_<source>_communities.csv`, `graph_<source>_community_summary.csv`,
    `graph_<source>_supergraph.csv`, `graph_<source>_bottleneck_edges.csv`, va
    `graph_<source>_bottleneck_nodes.csv`
  - compatibility/fan-out da co cho legacy community outputs cua DCCM,
    `lspmd_*` centrality files, va `tpt_bottleneck_*`; `-community=louvain`
    da route vao graph core cho source `dccm`
  - build `native` va `portable` da pass; smoke
    `test/titan_sparse_graph_smoke.cpp` va benchmark local native/portable cho
    `graph_dccm`, `graph_lspmd`, `graph_allostery_gc`, `graph_allostery_te`,
    `graph_tpt`, va `graph_auto` da pass trong workspace hien tai

7.7 `titan_awh_analysis`

Status:

- [done] Nen tang AWH RC-binned IFP da co san
  - `titan_awh_ifp.cpp` da doc explicit RC `.xvg` per walker, canh hang theo
    thoi gian frame, chia RC thanh cac bin deu nhau, va xuat:
    `awh_rc_bin_summary.dat`, `awh_ifp_by_bin.csv`,
    `awh_ifp_heatmap.dat`, `awh_rc_frame_map.dat`, va
    `awh_rep_bin_*.pdb`
  - workflow nay uu tien `structural/state mapping` tren biased ensemble, khong
    lay TICA/MSM lam truc chinh

- [done] AWH-1 provenance + overlap diagnostics
  - them `awh_walker_diagnostics.dat` va `awh_rc_overlap_diagnostics.dat`
    de report `matched/unmatched frames`, RC/time range, occupied bins theo
    walker, dominant walker theo bin, va tinh trang `single_walker` /
    `multiwalker` overlap
  - `awh_rc_bin_summary.dat` da bo sung header diagnostics:
    `occupied_bins`, `empty_bins`, `single_walker_bins`,
    `multiwalker_overlap_bins`
  - runtime warning da ro rang hon khi co qua nhieu frame khong match,
    bin trong, hoac overlap multiwalker kem
  - da verify bang smoke local duoi `test/awh_inputs_local`,
    `test/awh_ifp_only_local_native`, `test/awh_ifp_local_native`, va
    `test/awh_ifp_local_portable`; classic `ifp_occupancy.csv` giu hash-identical
    khi bat them `-awh-ifp`, va cac output AWH moi giu hash-identical giua
    native / portable

- [done] AWH-2 richer RC-binned structural mapping
  - them `awh_bin_structural_summary.dat` de tong hop structural load theo
    RC bin: `active_walkers`, `dominant_walker`, representative frame, tong
    residue-pair count, tong contact-event load, va breakdown theo 6 loai
    IFP (`HBond`, `SaltBridge`, `PiPi`, `CationPi`, `Hydrophobic`, `Halogen`)
  - them `awh_bin_top_contacts.csv` de rut gon top contact theo moi loai
    interaction trong tung RC bin, uu tien cach doc co che nhanh hon so voi
    quet toan bo `awh_ifp_by_bin.csv`
  - `awh_rc_frame_map.dat` da duoc enrich them cot `rc_left/center/right`,
    `bin_frames`, `active_walkers`, `dominant_walker`, `dominant_fraction`,
    va `overlap_status`; cac cot cu duoc giu nguyen o dau file
  - da verify bang smoke local:
    `ifp_occupancy.csv` van hash-identical giua baseline IFP va AWH path;
    cac file AWH cu giu nguyen hash ngoai `awh_rc_frame_map.dat` doi co chu
    dich, va toan bo file AWH moi + file enrich deu hash-identical giua
    native / portable

- [done] AWH-3 walker-consistency / stability
  - them `awh_walker_stability_by_bin.dat` de classify moi RC bin thanh:
    `empty`, `single_walker`, `multiwalker_balanced`,
    `multiwalker_imbalanced`; dong thoi report `effective_walkers`,
    `dominant_fraction`, va fractional support cua tung walker
  - them `awh_walker_pair_overlap.dat` de tong hop overlap theo cap walker:
    `shared_bins`, `bin_jaccard`, `shared_frames`, `frame_overlap_fraction`,
    va so bin `balanced_shared` / `imbalanced_shared`
  - `awh_rc_bin_summary.dat` da bo sung them header diagnostics:
    `multiwalker_balanced_bins` va `multiwalker_imbalanced_bins`
  - da verify bang 2 smoke local:
    fixture cu cho case `no overlap`, va fixture consistency moi
    `test/awh_ifp_consistency_native_post1`,
    `test/awh_ifp_consistency_portable_post1` cho thay ro ca
    `multiwalker_imbalanced`, `multiwalker_balanced`, va `single_walker`
  - `ifp_occupancy.csv` van hash-identical giua baseline IFP va AWH path; cac
    output AWH giu hash-identical giua native / portable

- [done] AWH-4 optional exploratory context
  - da them provenance `global_frame_index -> RC bin` qua
    `awh_frame_assignments.dat`
  - da them `awh_exploratory_context.dat` de overlay PCA/UMAP/MSM/TPT nhu lop
    context phu tren biased ensemble
  - giu wording `exploratory only`; khong dien giai unbiased kinetics hay
    equilibrium populations tu overlay nay

- workflow phan tich cho GROMACS AWH biased ensemble dua tren explicit RC `.xvg` ingest, RC binning, va snapshot selection theo reaction coordinate
- uu tien structural/state mapping va IFP dọc theo reaction coordinate de mo ta co che unbinding, thay vi xem TICA/MSM la truc chinh
- PCA/TICA/UMAP/MSM/TPT neu co dung thi chi la lop context phu; MSM/TPT tren AWH van chi nen dung o muc exploratory cho state/pathway mapping tuong doi, khong dien giai equilibrium populations hay absolute kinetics neu chua co reweighting/validation rieng
- xem tai lieu chi tiet trong `TITAN_AWH_ANALYSIS_PLAN.md`

Luu y:

- GPU offload chua nen uu tien truoc khi flat buffer + sparse backend da on
- module moi nen tai su dung backend chung, tranh copy them trajectory

================================================================================
8. Thu tu thuc thi khuyen nghi
================================================================================

Thu tu de xuat:

1. Phase 0
2. Phase 1
3. Phase 2
4. Phase 3
5. Phase 4 song song cuoi Phase 2 / dau Phase 3
6. Phase 5 sau khi backend chung da on

Neu can chia sprint ngan, nen bat dau bang 4 work items:

1. defensive numerical guards
2. chunk streamer copy removal
3. benchmark harness + regression guard
4. flat trajectory backend cho DCCM / network

Da xong trong sprint gan day:

1. benchmark harness + regression guard
2. defensive numerical guards
3. chunk streamer copy removal
4. PCA scratch reuse
5. GIST memory reduction

Sau do moi sang:

1. sparse TPT backend cho `tpm_`
2. shared sparse graph backend
3. portable build profile

================================================================================
9. Nhung diem can giu lai
================================================================================

Nhung phan dang di dung huong va nen duoc tai su dung:

- `titan_lj.cpp` da co pattern cell-list hop ly
- `titan_rms.cpp` kha gon va de toi uu tiep
- `titan_msm.cpp` da tranh duoc resize ngay tho bang doubling growth
- dense pair cache da co guard, can giu y tuong guard nay khi them module moi

================================================================================
10. Dinh nghia thanh cong
================================================================================

Roadmap nay duoc xem la thanh cong khi:

- Titan scale tot hon tren medium / large datasets
- peak RAM giam o cac module PCA / TPT / DCCM / GIST / streaming
- cac kernel AVX2 nhanh hon nhung van co fallback an toan
- module moi co the xay tren backend chung thay vi tu mang theo debt rieng
