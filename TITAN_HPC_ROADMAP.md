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
  - portable: khong phu thuoc `-march=native`
  - native: giu toi uu local

Deliverables:

- benchmark harness co output wall time + peak RAM
- bo input regression nho de test lai sau moi refactor
- tai lieu build profile ro rang

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

Target:

- `chunk_streamer.cpp`

3.2 PCA scratch reuse

- dua `dev` thanh scratch buffer tai su dung thay vi allocate moi frame
- xem lai loops de giu contiguous access va tranh branch khong can thiet

Expected impact:

- giam allocator pressure
- cai thien throughput o trajectory dai

Target:

- `titan_pca.cpp`

3.3 GIST reduction

- thay full per-thread grid clone bang tiled reduction hoac blocked merge
- neu can, giam bo nho thread-local theo slab / tile thay vi toan bo grid

Expected impact:

- giam peak RAM khi tang so thread
- giam serial tail reduction

Target:

- `titan_gist.cpp`

3.4 Defensive numerical guards

- dat guard sat kernel cho screened electrostatics
- thong nhat xu ly input bat thuong thanh error / warning ro rang

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

Target:

- module moi: `titan_sparse_graph.*` hoac backend tuong duong

4.3 Shared neighborhood backend

- dua cell-list / neighbor search thanh backend dung chung thay vi de logic nay
  phan tan
- uu tien tai su dung pattern dang tot trong `titan_lj.cpp`

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

- chuyen committor solve sang sparse linear system
- bo `gross_flux_` / `net_flux_` dense full matrix khi state graph da thua
- dominant path extraction nen chay tren edge list / sparse residual graph

Target:

- `titan_tpt.cpp`

5.3 DCCM / allostery pruning

- cho phep residue-distance pruning, kNN pruning, hoac contact-prior pruning
- tach transfer entropy thanh stage tuy chon vi day la path dat nhat
- nghien cuu estimator kNN / KDE cho mot so metric thay vi histogram moi cap

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
- khong gia dinh tat ca tinh nang C++23 deu on dinh tren libstdc++ cua GCC 13

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

1. benchmark harness + regression guard
2. chunk streamer copy removal
3. PCA scratch reuse
4. GIST memory reduction

Sau do moi sang:

1. flat trajectory backend cho DCCM
2. sparse TPT
3. shared sparse graph backend

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

