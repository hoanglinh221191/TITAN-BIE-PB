# Báo Cáo Phân Tích Cải Thiện Hiệu Năng Cấu Trúc Dữ Liệu

Dựa trên việc phân tích các file `.h` trong thư mục `src`, dưới đây là các cấu trúc dữ liệu có thể được tối ưu hóa để cải thiện hiệu năng (đặc biệt là cache locality và thời gian truy xuất), cùng với các vị trí cụ thể trong mã nguồn:

## 1. Vector Lồng Nhau (`std::vector<std::vector<T>>`)

**Vấn đề:** Việc sử dụng `std::vector` chứa `std::vector` gây phân mảnh bộ nhớ (memory fragmentation). Mỗi vector bên trong sẽ cấp phát bộ nhớ động tại một vùng heap khác nhau, dẫn đến hiện tượng "cache miss" liên tục khi duyệt qua các phần tử, làm giảm đáng kể hiệu suất đọc/ghi bộ nhớ.

**Giải pháp:** 
- Chuyển sang sử dụng mảng 1 chiều (flattened 1D array) `std::vector<T>` và tính toán chỉ số (index math): `index = row * num_cols + col`.
- Sử dụng các lớp ma trận (Matrix classes) hoặc `std::mdspan` (nếu dùng C++23) để bọc mảng 1 chiều này nhằm giữ nguyên giao diện 2D dễ đọc.

**Các vị trí phát hiện:**
- `gbsa_intrin.h` (dòng 59): `std::vector<std::vector<GBSA_energy_snap>> trajectories;`
- `nucleic_profiler.h` (dòng 60, 61): `worker_base_hits_`, `worker_backbone_hits_`
- `titan_dssp.h` (dòng 60, 139, 141): `trajectory_fractions_`, `trajectory_frame_profiles_`
- `titan_llps.h` (dòng 54, 245, 254): `contact_freq`, `macro_chains_`, `adj_ring_`
- `titan_llps_cluster.h` (dòng 766, 767): `cluster_pro`, `cluster_rna`
- `titan_membrane.h` (dòng 66): `std::vector<std::vector<std::size_t>> tails;`
- `titan_rotational_diffusion.h` (dòng 158): `std::vector<std::vector<Point3>> bodies_;`
- `titan_umap.h`, `titan_umap_knn.h`: Các hàm tính toán KNN.
- `titan_validation_core.h`: `StateSegments`, `implied_timescale_samples`, v.v.
- `tabi/tabi_types.h`, `tabi/tabi_fmm_m2l_prototype.h`

## 2. Lạm Dụng `std::map` 

**Vấn đề:** `std::map` trong C++ thường được cài đặt dưới dạng cây Đỏ-Đen (Red-Black tree). Nó cung cấp độ phức tạp $O(\log n)$ cho các thao tác thêm, xóa, sửa, nhưng lại có tính liền kề bộ nhớ (cache locality) rất kém do mỗi node được cấp phát động ngẫu nhiên trên heap. 

**Giải pháp:**
- Nếu không cần duyệt dữ liệu theo thứ tự (không cần ordered), hãy chuyển sang `std::unordered_map` (Hash map), giúp giảm độ phức tạp xuống $O(1)$.
- Đối với dữ liệu nhỏ hoặc ít khi thêm/xóa mà chỉ đọc, có thể dùng các flat maps như mảng đã sắp xếp `std::vector<std::pair<K, V>>` kết hợp `std::lower_bound` hoặc thư viện như `absl::flat_hash_map`.
- *Lưu ý:* Với `std::tuple` dùng làm key (ví dụ trong `titan_dihedral.h`), bạn sẽ cần định nghĩa thêm hàm băm (hash function) tùy chỉnh khi dùng `unordered_map`.

**Các vị trí phát hiện:**
- `titan_dihedral.h` (dòng 165-169): Sử dụng `std::map` với key là `std::tuple`. Đây có thể là một nút thắt cổ chai lớn (bottleneck) nếu kích thước map lớn.
- `titan_hbond_kinetics.h` (dòng 71): `std::map<std::size_t, std::uint64_t> lifetime_histogram_;` (Histogram thường không cần map dạng cây).
- `titan_membrane.h` (dòng 266-297): Rất nhiều map lưu trữ số đếm (`protein_contact_counts_`, `lipid_contact_counts_`, `order_accumulators_`, ...). Chuỗi `std::string` làm key cũng chậm hơn key là số nguyên.
- `titan_performance_observer.h` (dòng 24, 75): `std::map<std::string, double> counters;`
- `titan_s2.h` (dòng 218, 301, 302): `atom_by_name`, `aggregate_results_`.
- `tabi/analyzer_impl.h`: `pbsa_residue_stats_`, `pbsa_residue_traj_stats_`.

## 3. Lạm Dụng `std::set`

**Vấn đề:** Tương tự như `std::map`, `std::set` là một cây tìm kiếm nhị phân với tính cache locality rất thấp.

**Giải pháp:**
- Dùng `std::unordered_set`.
- Nếu quá trình chỉ thêm một lần và sau đó tìm kiếm (build-once, read-many), việc dùng `std::vector` + `std::sort` + `std::unique` kết hợp `std::binary_search` mang lại hiệu năng cao hơn rất nhiều.

**Các vị trí phát hiện:**
- `titan_membrane.h` (dòng 180, 206, 217, 281): `std::set<std::size_t> near_lipids`, `std::set<std::string> near_lipid_ids_`.

## 4. `std::list` Từ Thư Viện Bên Ngoài (Eigen)

**Vấn đề:** `std::list` là danh sách liên kết đôi (doubly-linked list), hoàn toàn không thân thiện với cache do các node nằm rời rạc trên heap.
**Phát hiện:** `std::list` có xuất hiện trong thư viện Eigen (`src/eigen/Eigen/src/IterativeLinearSolvers/IncompleteCholesky.h`). Vì đây là thư viện bên thứ ba (third-party), nên không nhất thiết phải sửa mã nguồn của nó trừ khi bạn có ý định dùng phiên bản tối ưu riêng.

## 5. Cải Tiến & Giảm Thiểu Template (Template Bloat)

**Vấn đề:** Template quá phức tạp làm chậm thời gian biên dịch (compile time) và có thể gây code bloat (mã nguồn phình to) nếu instance hóa quá nhiều.

**Giải pháp đề xuất chung:**
1. **Pimpl Idiom / Tách Interface:** Đưa các phần code không phụ thuộc vào template parameter ra các class base (non-template class) hoặc file `.cpp` để giảm thiểu lượng code phải sinh ra tại compile-time ở file header.
2. **Hạn chế Inline logic quá lớn:** Các phương thức template quá dài nên được cố gắng hạn chế nếu không yêu cầu inlining khắt khe.
3. **Type Erasure:** Cân nhắc thay thế template bằng tính đa hình tĩnh (Polymorphism / Interfaces) hoặc Type Erasure (như `std::function`, `std::any` nếu không thuộc các vùng code đòi hỏi hiệu suất quá cao (hot path)), để giảm số dòng code biên dịch lại khi thay đổi kiểu.

## 6. TICA Cho Hệ Nhiều Replica / Nhiều Snapshot

**Tình huống thực tế:** Với hệ 48 replica, MBAR báo 174720 aligned samples. Job `13_mbar_equilibrium_reweighting` ghi nhận `TICA Total Frames = 174720`, `TICA Feature Dim = 4656` với 97 C-alpha, `TICA Backend = exact_dual`, và `TICA + MSM Analysis = 35235.761 s` trong khi MBAR chỉ khoảng 255 s, UMAP khoảng 206 s. TICA/FES/UMAP/MSM downstream có thể phải xử lý toàn bộ số frame này. Nếu dùng feature mặc định, TITAN dựng vector khoảng cách cặp C-alpha với số chiều:

`D = n_ca * (n_ca - 1) / 2`

Vì vậy chi phí dựng feature là `O(n_frames * n_ca^2)`, nhưng với job đã đo, phần streaming accumulation không phải nút thắt chính; nút thắt nằm ở post-processing `exact_dual`, nơi TITAN giữ/cache ma trận feature `n_frames x D` rồi chạy SVD/projection-cache rebuild. `-nt=8` không tự động làm bước SVD/projection-cache này nhanh tuyến tính.

**Các vị trí phát hiện:**
- `impl_engine.cpp` (khoảng dòng 2126-2154): mặc định chọn feature TICA là pairwise C-alpha distances nếu không dùng `-feature-source` hẹp hơn.
- `impl_engine.cpp` (khoảng dòng 4373-4381): mỗi frame dựng toàn bộ vector khoảng cách cặp C-alpha bằng vòng lặp lồng nhau.
- `titan_kinetics.cpp` (`build_exact_dual_basis_`): backend exact-dual chạy thin SVD trên feature cache rất lớn; đây là ứng viên nút thắt chính trong job 48 replica đã đo.
- `titan_kinetics.cpp` (các hàm `build_randomized_basis_`, `solve_reduced_from_basis_`, `project_all_frames_to_basis_`): backend randomized đọc lại feature spool qua nhiều pass và nhân ma trận theo row, nhưng tránh SVD full tall matrix kiểu exact-dual.
- `titan_kinetics.h`: `frame_cache_` / spool đang lưu full feature vector theo frame cho backend large-system.
- `impl_engine.cpp` (khối ước lượng backend): job log ước lượng `TICA Estimated Frames = 316` trong khi thực tế là 174720 frame; cần sửa estimator cho input-list/replica runs để auto backend không chọn/khuyến khích sai.

**Cải tiến đề xuất:**
1. **Ưu tiên P0 cho case 48 replica:** tránh `exact_dual` cho ma trận cỡ `174720 x 4656`; dùng `randomized` hoặc backend streaming-sketch làm mặc định/khuyến nghị khi `n_frames * D` vượt ngưỡng.
2. **Sửa estimator backend:** log hiện ước lượng 316 frame nhưng thực tế 174720 frame. Backend auto và memory estimate phải dựa trên số frame sau `-b/-e/-dt` nhân với số trajectory/replica hữu dụng.
3. **Tách projection-cache theo mức cần thiết:** downstream MSM/UMAP/FES chỉ cần TICA-space projections, không phải lúc nào cũng cần giữ full raw feature history. Nên lưu/cache projection hoặc compressed basis rows khi đã đủ để replay.
4. **Giảm ghi/đọc full feature matrix khi dùng `randomized`:** cân nhắc backend streaming-sketch, trong đó pairwise distances được nhân trực tiếp vào random sketch theo block, chỉ lưu representation đã giảm chiều thay vì spool full `n_frames x D` khi không cần replay full feature.
5. **Thêm instrumentation:** log riêng thời gian `exact_dual_svd`, `feature build`, `feature cache/spool write/read`, `randomized sketch`, `reduced solve`, `projection replay`, kèm `n_frames`, `n_ca`, `D`, backend và `tica_rank`.
6. **Tùy chọn chạy nhanh trước khi refactor code:** với validation/FES sơ bộ, dùng `-tica-backend randomized`, giảm `-tica-rank` xuống 16-32 nếu chấp nhận giảm rank, hoặc dùng `-feature-source dihedral`, `-feature-source qnative`, `-feature-source cmap` thay vì mặc định pairwise C-alpha.

**Ghi chú MBAR/TICA weights:** Log kiểu `mapped weights sum = 0` là vấn đề mapping/metadata trọng số riêng, không phải nguyên nhân trực tiếp làm TICA chậm. Trong job 13, `mbar_summary.dat` ghi `MBAR Usable Samples = 174720` nhưng `MBAR Mapped Global Frames = 3640`; `mbar_frame_weights.dat` có nhiều dòng đầu `global_frame_index = -1`. Điều này cho thấy mapping global frame weights chưa bao phủ toàn bộ TICA frame space và cần sửa/diagnose riêng ở metadata chọn cho `mbar-domain=both`.

## Tổng Kết

Các bước để thực hiện tối ưu hiệu năng dựa trên bản báo cáo này:
1. Với hệ lớn/nhiều replica, ưu tiên trước **TICA pairwise-feature path**: đo và tối ưu feature build, randomized/spool passes, projection replay, và giảm full feature matrix khi có thể.
2. Viết một lớp (class) tiện ích đại diện cho **Ma Trận 2D phẳng (Flat 2D Matrix)** và thay thế dần các `std::vector<std::vector<T>>` ở hot path, không thay đại trà khi chưa đo.
3. Thay các đoạn sử dụng **`std::map` bằng `std::unordered_map` hoặc dense/flat vector**, đặc biệt là cho những biểu đồ đếm (counters) và ánh xạ dữ liệu (như chuỗi `string` sang số liệu).
4. Sử dụng **dense vector / flat index** hoặc Hash Maps với custom hasher cho `std::tuple` thay cho việc dùng `std::map<std::tuple, ...>`.
5. Xem xét lại chiến lược lưu trữ dữ liệu chuỗi (String keys) trong map để tối giản chi phí băm. Sử dụng index thay thế cho string (ví dụ ánh xạ tên nguyên tử sang số nguyên rồi dùng array/vector map) nếu có thể.

## 7. TICA Timing Instrumentation And Projection Cache Update (2026-06-17)

Implemented first-stage instrumentation for the high-frame-count TICA path:

- `TitanKinetics::TimingStats` now separates feature build, cache write/read, spool write/read, exact-dual center/copy, exact-dual SVD, randomized sketch, randomized QR, randomized compression, randomized SVD, reduced-solve replay, reduced eigensolve, and projection replay.
- `run_tica_msm_postproc` logs these phase timings and frame/byte counters after TICA/MSM post-processing, so future job 13/14 reruns can identify whether the bottleneck is SVD, replay, or I/O.
- `impl_engine.cpp` reports feature-build time from the trajectory loop separately from `accumulate_frame`, projection, and MSM metadata bookkeeping.
- `TitanKinetics::project_all_frames_to_sink` streams solved TICA projections row-by-row into a caller-provided sink.
- `TitanMSM::project_pending_features` now fills `proj_mat_` directly from that sink instead of first materializing a second full projection matrix and then copying it.
- After the reduced projection cache is built, exact-dual raw feature rows and MSM pending raw features are released. This does not remove the need for raw features during exact-dual SVD, but it reduces memory pressure for downstream MSM/FES/UMAP exports.

Remaining larger refactor:

- The current randomized backend still writes and rereads the full `n_frames x D` feature spool for sketch, compression, reduced solve, and projection replay.
- A true streaming-sketch backend should remain opt-in until measured. The next target is to compute reduced sketch rows during feature construction or in larger blocks, so TITAN can avoid keeping/replaying the full raw feature matrix when downstream only needs low-dimensional TICA projections.

## 8. TICA P0 Auto Backend Policy Update (2026-06-29)

Implemented the P0 guard for large high-dimensional TICA runs:

- `src/titan_tica_backend_policy.h` centralizes the auto-backend decision used by `impl_engine.cpp`.
- Auto mode now computes `raw_feature_values = estimated_frames * feature_dim`.
- If `raw_feature_values >= 250000000`, auto mode selects `randomized` before considering RAM-fit `exact_dense` or `exact_dual`.
- Explicit user requests such as `-tica-backend exact_dual` are still honored, but the log reports that the large-matrix limit was exceeded.
- The 48-replica case `174720 x 4656 = 813496320` therefore auto-selects `randomized`.
- `impl_engine.cpp` logs `TICA Raw Feature Values`, `TICA Large-Matrix Limit`, and `TICA Large-Matrix Policy` so reruns of job 13/14 show exactly why the backend was selected.

Next priority after this update is to rerun the 48-replica job with the new binary and inspect the timing instrumentation. If randomized is still dominated by spool read/write or projection replay, the next structural refactor remains the streaming-sketch/reduced-spool backend from section 7.
