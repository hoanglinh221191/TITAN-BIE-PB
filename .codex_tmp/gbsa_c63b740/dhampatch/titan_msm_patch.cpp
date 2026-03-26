// ============================================================================
// titan_msm.h — PATCH: additions needed for DHAM integration
//
// Chỉ chứa những phần CẦN THÊM vào titan_msm.h hiện tại.
// Tìm class TitanMSM và chèn các dòng được đánh dấu [ADD].
// ============================================================================

// ─── TRONG private members, THÊM VÀO SAU các member hiện có ─────────────────
//
//   [ADD] Replica tracking for DHAM
//   std::vector<int>    frame_replica_ids_;   // replica_id per proj_mat_ row
//   std::vector<double> frame_timestamps_;    // time_ps per proj_mat_ row


// ─── TRONG public interface, THÊM VÀO SAU add_projected_frame() ──────────────
//
//   Thay prototype hiện tại:
//     void add_projected_frame(const Eigen::VectorXd& proj_x);
//
//   Bằng overload backward-compatible:
//     void add_projected_frame(const Eigen::VectorXd& proj_x,
//                              int replica_id = -1);
//
//   [ADD] thêm:
//     void record_frame_time(double time_ps, int replica_id = -1);
//
//   [ADD] Accessors for DHAM postproc:
//     const std::vector<int>&    get_frame_replica_ids() const;
//     const std::vector<double>& get_frame_timestamps()  const;
//
//   [ADD] Per-replica dtraj extraction (called after train_kmeans):
//     // Returns map: replica_id → discrete trajectory of that replica
//     std::unordered_map<int, std::vector<int>>
//         get_per_replica_dtrajs() const;
//
//     // Returns map: replica_id → timestamps (ps) of that replica
//     std::unordered_map<int, std::vector<double>>
//         get_per_replica_timestamps() const;

// ─── TRONG titan_msm.cpp, SỬA add_projected_frame() ─────────────────────────
//
// void TitanMSM::add_projected_frame(const Eigen::VectorXd& proj_x,
//                                     int replica_id) {
//     if (n_proj_frames_ >= proj_mat_.rows()) {
//         int64_t new_rows = std::max(int64_t(1024), proj_mat_.rows() * 2);
//         proj_mat_.conservativeResize(new_rows, n_comp_);
//         // [ADD] reserve vectors to avoid repeated realloc
//         frame_replica_ids_.reserve(static_cast<size_t>(new_rows));
//         frame_timestamps_.reserve(static_cast<size_t>(new_rows));
//     }
//     proj_mat_.row(n_proj_frames_++) = proj_x.transpose();
//     // [ADD]
//     frame_replica_ids_.push_back(replica_id);
//     frame_timestamps_.push_back(0.0);   // overwritten by record_frame_time
// }
//
// // [ADD] new method
// void TitanMSM::record_frame_time(double time_ps, int /*replica_id*/) {
//     if (!frame_timestamps_.empty())
//         frame_timestamps_.back() = time_ps;
// }
//
// // [ADD] new method
// std::unordered_map<int, std::vector<int>>
// TitanMSM::get_per_replica_dtrajs() const {
//     std::unordered_map<int, std::vector<int>> result;
//     for (size_t i = 0; i < discrete_traj_.size(); ++i) {
//         int rid = (i < frame_replica_ids_.size())
//                   ? frame_replica_ids_[i] : -1;
//         result[rid].push_back(discrete_traj_[i]);
//     }
//     return result;
// }
//
// // [ADD] new method
// std::unordered_map<int, std::vector<double>>
// TitanMSM::get_per_replica_timestamps() const {
//     std::unordered_map<int, std::vector<double>> result;
//     for (size_t i = 0; i < frame_timestamps_.size(); ++i) {
//         int rid = (i < frame_replica_ids_.size())
//                   ? frame_replica_ids_[i] : -1;
//         result[rid].push_back(frame_timestamps_[i]);
//     }
//     return result;
// }
//
// // [ADD] Thêm vào clear_projected():
// void TitanMSM::clear_projected() {
//     n_proj_frames_ = 0;
//     traj_break_indices_.clear();
//     discrete_traj_.clear();
//     clustered_ = false;
//     msm_built_  = false;
//     eigen_solved_ = false;
//     frame_replica_ids_.clear();   // [ADD]
//     frame_timestamps_.clear();    // [ADD]
// }

// ─── TRONG arguments.h, THÊM DhamArgs vào struct arguments ──────────────────
//
// struct DhamArgs {
//     bool   flag_dham       {false};
//     std::string temps_file {};       // file chứa danh sách nhiệt độ (K)
//     std::string ener_list_file {};   // file chứa danh sách path .xvg
//     double target_temp     {300.0};  // nhiệt độ mục tiêu (K)
//     int    max_iter        {5000};
//     double tolerance       {1e-6};
//     double dt_xtc_ps       {10.0};   // output freq của XTC (ps)
//     double dt_energy_ps    {2.0};    // output freq của energy .xvg (ps)
// };
//
// // Thêm vào struct arguments public members (sau IfpArgs ifp):
// DhamArgs dham;

// ─── TRONG impl_engine.cpp — Fix 1: TICA reference alignment ─────────────────
//
// Tìm đoạn (khoảng dòng 806):
//   // Fix 5: hoist TICA workspace OUTSIDE streaming loops
//   const int tica_n_ca = ...
//   ...
//   bool tica_ref_initialized = false;   ← XÓA DÒNG NÀY
//
// Di chuyển TOÀN BỘ khối này ra TRƯỚC while(true) outer loop (dòng 348).
// Thêm suffix _global vào tên biến để tránh nhầm lẫn:
//
//   // === ĐẶT TRƯỚC while(true) — dòng ~346 ===
//   int tica_n_ca_global = 0, tica_feat_dim_global = 0;
//   Eigen::VectorXd tica_features_global;
//   Eigen::MatrixXd tica_mobile_global, tica_ref_centered_global;
//   Eigen::Matrix3d tica_R_global = Eigen::Matrix3d::Identity();
//   bool tica_ref_initialized_global = false;   // chỉ false DUY NHẤT một lần
//
//   if (Arguments.msm.flag_tica && tica_analyzer) {
//       tica_n_ca_global = static_cast<int>(dccm_analyzer.get_resnrs().size());
//       tica_feat_dim_global = tica_analyzer->feature_dim();
//       if (tica_n_ca_global > 0) {
//           tica_features_global.resize(tica_feat_dim_global);
//           tica_mobile_global.resize(tica_n_ca_global, 3);
//           tica_ref_centered_global.resize(tica_n_ca_global, 3);
//       }
//   }
//
// Bên trong while(true) outer loop, thay khai báo cũ bằng alias:
//
//   const int& tica_n_ca       = tica_n_ca_global;
//   const int& tica_feat_dim   = tica_feat_dim_global;
//   auto& tica_features        = tica_features_global;
//   auto& tica_mobile          = tica_mobile_global;
//   auto& tica_ref_centered    = tica_ref_centered_global;
//   auto& tica_R               = tica_R_global;
//   bool& tica_ref_initialized = tica_ref_initialized_global;
//
// Tại ranh giới trajectory (dòng ~1993), GIỮ NGUYÊN reset_lag_buffer(),
// KHÔNG reset tica_ref_initialized:
//
//   if (Arguments.msm.flag_tica && tica_analyzer) {
//       tica_analyzer->reset_lag_buffer();
//       // KHÔNG thêm: tica_ref_initialized_global = false;
//   }

// ─── TRONG impl_engine.cpp — Fix 3: truyền replica_id khi accumulate ─────────
//
// Tìm dòng ~1725:
//   msm_analyzer->add_projected_frame(proj);
//
// Thay bằng:
//   int current_replica = input_name - Arguments.b_input;
//   msm_analyzer->add_projected_frame(proj, current_replica);
//   msm_analyzer->record_frame_time(
//       pro_in.m_mol[chunk_snap_idx].time_of_snap, current_replica);
