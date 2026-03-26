# TitanDHAM — Hướng dẫn Fix Toàn diện

Tài liệu này chứa **code patch cụ thể** cho từng bug được nhận diện,
theo thứ tự ưu tiên từ nghiêm trọng nhất.

---

## Fix 1 (NGHIÊM TRỌNG) — TICA Reference Frame bị reset giữa các replica

### Vấn đề
`tica_ref_initialized` được khai báo **bên trong** `while(true)` outer loop
(dòng ~820 `impl_engine.cpp`). Mỗi lần vòng lặp ngoài chạy qua một replica mới,
biến này bị reset về `false` → frame đầu tiên của replica đó trở thành
reference mới → mỗi replica được align vào một hệ tọa độ khác nhau →
toàn bộ không gian TICA mất tính nhất quán across replicas.

### Fix — impl_engine.cpp

**Bước 1:** Tìm đoạn khai báo workspace TICA (khoảng dòng 806-825):

```cpp
// Fix 5: hoist TICA workspace OUTSIDE streaming loops
// (pairwise distance vector + Kabsch mobile/reference matrices — zero heap
// per frame)
const int tica_n_ca =
    Arguments.msm.flag_tica && tica_analyzer
        ? static_cast<int>(dccm_analyzer.get_resnrs().size())
        : 0;
// ... (các khai báo khác)
bool tica_ref_initialized = false;   // <-- ĐÂY LÀ VẤN ĐỀ
```

**Bước 2:** Di chuyển TOÀN BỘ khối khai báo workspace TICA ra ngoài
`while(true)` outer loop, đặt ngay sau khi `tica_analyzer` được khởi tạo
(khoảng dòng 696, sau khối `if (Arguments.msm.flag_tica && !tica_analyzer)`).

```cpp
// === ĐẶT SAU KHỐI KHỞI TẠO tica_analyzer (dòng ~696) ===
// TICA workspace — khai báo MỘT LẦN, dùng chung cho TẤT CẢ replica.
// Quan trọng: tica_ref_initialized = false chỉ một lần duy nhất.
// Tất cả replica sẽ được align vào reference frame của REPLICA ĐẦU TIÊN.
int tica_n_ca_global = 0;
int tica_feat_dim_global = 0;
Eigen::VectorXd tica_features_global;
Eigen::MatrixXd tica_mobile_global;
Eigen::MatrixXd tica_ref_centered_global;   // Reference dùng chung
Eigen::Matrix3d tica_R_global;
bool tica_ref_initialized_global = false;   // Chỉ false một lần duy nhất

// Khởi tạo kích thước sau khi tica_analyzer được tạo
if (Arguments.msm.flag_tica && tica_analyzer) {
    tica_n_ca_global = static_cast<int>(dccm_analyzer.get_resnrs().size());
    tica_feat_dim_global = tica_analyzer->feature_dim();
    if (tica_n_ca_global > 0) {
        tica_features_global.resize(tica_feat_dim_global);
        tica_mobile_global.resize(tica_n_ca_global, 3);
        tica_ref_centered_global.resize(tica_n_ca_global, 3);
    }
}
```

**Bước 3:** Bên trong outer loop, **xóa** đoạn khai báo workspace TICA cũ
(các dòng khai báo `tica_n_ca`, `tica_feat_dim`, `tica_features`,
`tica_mobile`, `tica_ref_centered`, `tica_R`, `tica_ref_initialized`).
Thay bằng reference đến biến global:

```cpp
// Bên trong while(true) outer loop, thay thế khai báo cũ bằng:
// (không khai báo lại — dùng biến _global đã khai báo ở trên)
const int& tica_n_ca = tica_n_ca_global;
const int& tica_feat_dim = tica_feat_dim_global;
auto& tica_features = tica_features_global;
auto& tica_mobile = tica_mobile_global;
auto& tica_ref_centered = tica_ref_centered_global;
auto& tica_R = tica_R_global;
bool& tica_ref_initialized = tica_ref_initialized_global;
```

**Bước 4:** Tại ranh giới trajectory (dòng ~1993), giữ nguyên
`reset_lag_buffer()` nhưng **KHÔNG** reset `tica_ref_initialized`:

```cpp
// Reset TICA lag buffer at trajectory boundary (giữ nguyên)
if (Arguments.msm.flag_tica && tica_analyzer) {
    tica_analyzer->reset_lag_buffer();
    // KHÔNG reset tica_ref_initialized ở đây.
    // Tất cả replica dùng chung reference frame của replica đầu tiên.
}
```

---

## Fix 2 (NGHIÊM TRỌNG) — Energy-DTraj Alignment (dt synchronization)

### Vấn đề
File `.xvg` từ `gmx energy` và file XTC trajectory có thể có output frequency
khác nhau. Ví dụ: XTC ghi mỗi 10 ps, energy ghi mỗi 2 ps → 5 energy frames
ứng với 1 trajectory frame → nếu map 1:1 theo index thì sai hoàn toàn.

### Fix — Cách tiếp cận đúng

Thay vì dùng index, **map theo timestamp**. Thêm field vào `DhamArgs`:

```cpp
// Trong arguments.h, struct DhamArgs:
struct DhamArgs {
    bool flag_dham{false};
    std::string temps_file{};
    std::string ener_list_file{};
    double target_temp{300.0};
    int    max_iter{5000};
    double tolerance{1e-6};
    double dt_xtc_ps{10.0};    // Output frequency của XTC (ps) — từ mdp ncoordout
    double dt_energy_ps{2.0};  // Output frequency của energy (ps) — từ mdp nstxout-energy
    // Nếu dt_xtc_ps / dt_energy_ps = N nguyên → lấy mỗi N energy frame 1 lần.
    // Nếu không chia hết → dùng nearest-neighbor matching theo thời gian.
};
```

**Hàm parse .xvg với timestamp matching** (đặt trong `titan_dham.cpp`):

```cpp
// Đọc file .xvg GROMACS, trả về vector<pair<double,double>> = {time_ps, energy_kJ}
// Bỏ qua dòng bắt đầu bằng '#' hoặc '@'
static std::vector<std::pair<double,double>> parse_xvg(const std::string& path) {
    std::vector<std::pair<double,double>> result;
    std::ifstream f(path);
    if (!f.is_open()) {
        TitanUI::print_error("DHAM: Cannot open energy file: " + path);
        return result;
    }
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#' || line[0] == '@') continue;
        std::istringstream ss(line);
        double t, e;
        if (ss >> t >> e)
            result.emplace_back(t, e);
    }
    return result;
}

// Map energy values vào dtraj timestamps
// dtraj_times: vector thời gian của từng frame trong dtraj (ps)
// xvg_data: output của parse_xvg()
// Trả về vector<double> cùng độ dài với dtraj_times
static std::vector<double> align_energy_to_dtraj(
    const std::vector<double>& dtraj_times,
    const std::vector<std::pair<double,double>>& xvg_data)
{
    std::vector<double> aligned(dtraj_times.size(),
                                std::numeric_limits<double>::quiet_NaN());
    if (xvg_data.empty()) return aligned;

    size_t j = 0; // pointer vào xvg_data
    for (size_t i = 0; i < dtraj_times.size(); ++i) {
        double t_target = dtraj_times[i];
        // Tìm entry xvg gần nhất về thời gian
        while (j + 1 < xvg_data.size() &&
               std::abs(xvg_data[j+1].first - t_target) <
               std::abs(xvg_data[j].first  - t_target))
            ++j;
        // Chỉ gán nếu sai lệch < 0.5 * dt_energy (tránh gán nhầm)
        if (std::abs(xvg_data[j].first - t_target) < 0.5 * 
            (xvg_data.size() > 1 ? xvg_data[1].first - xvg_data[0].first : 2.0))
            aligned[i] = xvg_data[j].second;
    }

    // Kiểm tra NaN: nếu > 10% frame bị NaN, cảnh báo
    size_t nan_count = std::count_if(aligned.begin(), aligned.end(),
        [](double v){ return std::isnan(v); });
    if (nan_count > aligned.size() / 10) {
        TitanUI::print_warning("DHAM: " + std::to_string(nan_count) +
            "/" + std::to_string(aligned.size()) +
            " frames could not be matched to energy values. "
            "Check --dham-dt-xtc and --dham-dt-energy.");
    }
    return aligned;
}
```

**Lưu timestamp khi accumulate frame** trong `impl_engine.cpp`,
trong khối TICA accumulation:

```cpp
// Trong inner loop, sau tica_analyzer->accumulate_frame(tica_features):
if (Arguments.dham.flag_dham && msm_analyzer) {
    // Lưu timestamp của frame này để sau mapping với energy
    msm_analyzer->record_frame_time(
        pro_in.m_mol[chunk_snap_idx].time_of_snap,  // ps
        input_name - Arguments.b_input              // replica index
    );
}
```

---

## Fix 3 — Replica ID Tracking trong TitanMSM

### Vấn đề
`add_projected_frame()` không lưu thông tin replica. Sau K-means, không biết
frame nào thuộc replica nào để xây dựng $C^k_{ij}$ riêng biệt.

### Fix — titan_msm.h: thêm overload backward-compatible

```cpp
// Trong class TitanMSM (titan_msm.h):

// Thêm vào private members:
std::vector<int>    frame_replica_ids_;  // replica_id tương ứng với proj_mat_ row
std::vector<double> frame_timestamps_;   // timestamp (ps) của từng frame

// Public interface mới (backward compatible — replica_id=-1 = mode đơn replica):
void add_projected_frame(const Eigen::VectorXd& proj_x, int replica_id = -1);
void record_frame_time(double time_ps, int replica_id = -1);

// Accessor cho DHAM postproc:
const std::vector<int>&    get_frame_replica_ids() const { return frame_replica_ids_; }
const std::vector<double>& get_frame_timestamps()  const { return frame_timestamps_; }

// Xây dựng per-replica discrete trajectory sau K-means:
// Trả về map: replica_id → vector<int> dtraj của replica đó
std::unordered_map<int, std::vector<int>> get_per_replica_dtrajs() const;
std::unordered_map<int, std::vector<double>> get_per_replica_timestamps() const;
```

### Fix — titan_msm.cpp: implementation

```cpp
void TitanMSM::add_projected_frame(const Eigen::VectorXd& proj_x,
                                    int replica_id) {
    // Grow proj_mat_ với doubling strategy (giữ nguyên logic cũ)
    if (n_proj_frames_ >= proj_mat_.rows()) {
        int64_t new_rows = std::max(int64_t(1024), proj_mat_.rows() * 2);
        proj_mat_.conservativeResize(new_rows, n_comp_);
        frame_replica_ids_.reserve(static_cast<size_t>(new_rows));
        frame_timestamps_.reserve(static_cast<size_t>(new_rows));
    }
    proj_mat_.row(n_proj_frames_++) = proj_x.transpose();
    frame_replica_ids_.push_back(replica_id);
    frame_timestamps_.push_back(0.0);  // overwrite bởi record_frame_time
}

void TitanMSM::record_frame_time(double time_ps, int /*replica_id*/) {
    if (!frame_timestamps_.empty())
        frame_timestamps_.back() = time_ps;
}

std::unordered_map<int, std::vector<int>>
TitanMSM::get_per_replica_dtrajs() const {
    std::unordered_map<int, std::vector<int>> result;
    for (size_t i = 0; i < discrete_traj_.size(); ++i) {
        int rid = (i < frame_replica_ids_.size()) ? frame_replica_ids_[i] : -1;
        result[rid].push_back(discrete_traj_[i]);
    }
    return result;
}

std::unordered_map<int, std::vector<double>>
TitanMSM::get_per_replica_timestamps() const {
    std::unordered_map<int, std::vector<double>> result;
    for (size_t i = 0; i < frame_timestamps_.size(); ++i) {
        int rid = (i < frame_replica_ids_.size()) ? frame_replica_ids_[i] : -1;
        result[rid].push_back(frame_timestamps_[i]);
    }
    return result;
}
```

### Fix — impl_engine.cpp: truyền replica_id khi accumulate

```cpp
// Trong inner loop (dòng ~1721), thay thế:
//   msm_analyzer->add_projected_frame(proj);
// bằng:
if (Arguments.msm.flag_msm && msm_analyzer && tica_analyzer->is_solved()) {
    Eigen::VectorXd proj = tica_analyzer->project_frame(tica_features);
    int current_replica = input_name - Arguments.b_input;  // 0-indexed
    msm_analyzer->add_projected_frame(proj, current_replica);
    msm_analyzer->record_frame_time(
        pro_in.m_mol[chunk_snap_idx].time_of_snap,
        current_replica);
}
```

---

## Fix 4 — Validation: số replica vs số nhiệt độ

### Vấn đề
Nếu `temps_file` có 48 nhiệt độ nhưng `input.txt` chỉ có 30 file, hoặc ngược lại,
toàn bộ DHAM sẽ tính sai im lặng.

### Fix — arguments.cpp: thêm vào `post_process()`

```cpp
// Thêm vào cuối post_process() trong arguments.cpp,
// sau phần validate listfile hiện tại:
if (dham.flag_dham) {
    // Đọc và đếm nhiệt độ
    if (dham.temps_file.empty()) {
        TitanUI::print_error("DHAM: --dham-temps is required.");
        exit(EXIT_FAILURE);
    }
    std::ifstream tf(dham.temps_file);
    if (!tf.is_open()) {
        TitanUI::print_error("DHAM: Cannot open temps file: " + dham.temps_file);
        exit(EXIT_FAILURE);
    }
    int n_temps = 0;
    double tmp;
    while (tf >> tmp) ++n_temps;
    tf.close();

    // Đếm số replica từ input list
    int n_replicas_from_input = (flag_listfile)
        ? [&]() {
            std::ifstream lf(list_filename);
            std::string ln;
            int cnt = 0;
            while (std::getline(lf, ln)) if (!ln.empty()) ++cnt;
            return cnt;
          }()
        : (e_input - b_input + 1);

    if (n_temps != n_replicas_from_input) {
        TitanUI::print_error(
            "DHAM: Number of temperatures in " + dham.temps_file +
            " (" + std::to_string(n_temps) + ") does not match "
            "number of replica trajectories (" +
            std::to_string(n_replicas_from_input) + "). Aborting.");
        exit(EXIT_FAILURE);
    }

    // Đọc và đếm energy files
    if (!dham.ener_list_file.empty()) {
        std::ifstream ef(dham.ener_list_file);
        if (!ef.is_open()) {
            TitanUI::print_error("DHAM: Cannot open energy list file: " +
                                 dham.ener_list_file);
            exit(EXIT_FAILURE);
        }
        int n_ener_files = 0;
        std::string eline;
        while (std::getline(ef, eline)) if (!eline.empty()) ++n_ener_files;
        ef.close();
        if (n_ener_files != n_temps) {
            TitanUI::print_error(
                "DHAM: Number of energy files (" +
                std::to_string(n_ener_files) +
                ") must match number of temperatures (" +
                std::to_string(n_temps) + ").");
            exit(EXIT_FAILURE);
        }
    }
    TitanUI::print_info("[DHAM] Validated: " + std::to_string(n_temps) +
                        " replicas, temps file OK, energy list OK.");
}
```

---

## Fix 5 — LSE zero-guard trong TitanDHAM::solve()

### Vấn đề
`c_out_(i,k)` có thể bằng 0 nếu microstate `i` không được thăm ở replica `k`.
Gọi `.log()` trên 0 → `-inf` → sau `exp()` → 0 (underflow) hoặc NaN.

### Fix — titan_dham.cpp: trong vòng lặp tự tương thích

```cpp
// Trong hàm solve(), trước vòng lặp LSE:
// Thay vì:
//   Eigen::ArrayXXd log_terms = f_k_ - u_ik_ + c_out_.log();   // BUG nếu c_out_=0

// Dùng:
// Thêm tiny epsilon vào c_out_ trước khi log để tránh log(0)
// Giá trị 1e-300 đủ nhỏ để không ảnh hưởng kết quả vật lý
// (count=0 sẽ đóng góp exp(-∞) ≈ 0 vào tổng, đúng về mặt toán học)
constexpr double LOG_ZERO_GUARD = 1e-300;

Eigen::ArrayXXd log_c = (c_out_.array() + LOG_ZERO_GUARD).log(); // [n_states x n_replicas]

// Tính log-terms: mỗi hàng i, mỗi cột k:
// log_term(i,k) = f_k(k) - u_ik(i,k) + log(c_out(i,k))
Eigen::ArrayXXd log_terms = 
    f_k_.transpose().array().replicate(n_states_, 1)  // broadcast f_k theo hàng
    - u_ik_.array()
    + log_c;

// LSE theo từng hàng (sum over replicas k for each state i):
Eigen::VectorXd new_weights(n_states_);
for (int i = 0; i < n_states_; ++i) {
    double M = log_terms.row(i).maxCoeff();
    if (!std::isfinite(M)) {
        // State này không được thăm ở BẤT KỲ replica nào → weight = 0
        new_weights(i) = 0.0;
        continue;
    }
    double sum_exp = (log_terms.row(i) - M).exp().sum();
    new_weights(i) = std::exp(M) * sum_exp;
}

// Normalize để tránh overflow tích lũy qua các iteration
double total = new_weights.sum();
if (total > 1e-300)
    new_weights /= total;
```

---

## Tóm tắt thứ tự implement

| Thứ tự | Fix | File | Mức độ ảnh hưởng |
|--------|-----|------|-----------------|
| 1 | TICA reference alignment | `impl_engine.cpp` | 🔴 Data correctness |
| 2 | Energy-DTraj timestamp matching | `titan_dham.cpp` (mới) | 🔴 Data correctness |
| 3 | Replica ID trong TitanMSM | `titan_msm.h/.cpp` | 🟠 Required for DHAM |
| 4 | Validation n_replicas == n_temps | `arguments.cpp` | 🟡 Safety |
| 5 | LSE zero-guard | `titan_dham.cpp` (mới) | 🟡 Numerical stability |

Fix 1 có thể implement **ngay bây giờ** mà không cần viết bất kỳ code DHAM nào —
nó sẽ cải thiện tính chính xác của TICA/MSM đơn replica lẫn multi-replica.
Fix 3 là backward-compatible (default `replica_id = -1`), có thể merge song song.
Fix 2, 4, 5 chỉ cần khi bắt đầu viết `titan_dham.cpp`.
