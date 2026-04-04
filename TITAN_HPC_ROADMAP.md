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
    `AVX2 + FMA + generic tuning`; full runtime dispatch van la phase sau
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
  - allostery da co residue-distance / contact-prior pruning opt-in qua
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
- [next] Phase 3 / 5.3 `DCCM / allostery pruning` slice 2
  - them kNN pruning / contact-prior ranking de tranh cat het long-range pair
    co tin hieu, thay vi chi dung cutoff hinh hoc thuan tuy
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
- [next] toi uu reverse-edge lookup / chia se sparse graph backend thay vi moi
  module tu quan ly CSR rieng

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

- [next] Slice 3: estimator / prior refinement
  - can nhac true contact-frequency prior hoac estimator kNN/KDE cho mot so
    metric de giam bias cua histogram theo tung cap residue
  - tach ro `pruning window` khoi ten `cut-mode` neu CLI bat dau mo rong them
    nhieu prior khac nhau

Target:

- `dccm_network.cpp`

5.4 MSM spectral cost

- ra soat ITS scan va eigensolve per lag
- neu K lon, can nhac partial spectrum / Schur / reuse intermediate structures

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

7.1 `titan_mbar` hoac `titan_wham`

- mo rong tu `titan_dham`
- phu hop cho free-energy / reweighting workflow nghiem tuc hon

7.2 `titan_hbond_kinetics`

- hydrogen bond
- salt bridge
- water bridge
- residence time / intermittent correlation

7.3 `titan_validation`

- Chapman-Kolmogorov test
- VAMP-2 / score-based model selection
- bootstrap uncertainty cho MSM / TPT

7.4 `titan_corr_knn`

- generalized correlation
- mutual information
- transfer entropy bang estimator tot hon histogram dense

7.5 `titan_sq`

- structure factor `S(q)`
- streaming RDF / pair correlation
- phu hop LJ / LLPS / condensed-phase workflows

7.6 `titan_sparse_graph`

- graph analytics dung chung
- community detection
- edge filtering
- bottleneck ranking

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
