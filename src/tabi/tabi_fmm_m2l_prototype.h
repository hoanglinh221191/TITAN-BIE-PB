#ifndef TITAN_TABI_FMM_M2L_PROTOTYPE_H
#define TITAN_TABI_FMM_M2L_PROTOTYPE_H

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

class TitanExecutor;

namespace titan::tabi::fmm_m2l {

enum class KernelKind {
  Laplace,
  Yukawa,
};

struct PrototypeConfig {
  KernelKind kernel = KernelKind::Laplace;
  double kappa = 0.0;
  int max_leaf_size = 64;
  int expansion_order = 4;
  double separation_theta = 0.6;
  int threads = 1;
};

struct PrototypeMetrics {
  std::size_t sources = 0;
  std::size_t targets = 0;
  std::size_t cells = 0;
  std::size_t near_pairs = 0;
  std::size_t far_pairs = 0;
  std::size_t m2l_operators = 0;
  double setup_seconds = 0.0;
  double apply_seconds = 0.0;
  double direct_seconds = 0.0;
  double rms_abs_error = 0.0;
  double max_abs_error = 0.0;
  double rel_rms_error = 0.0;
  bool finite_outputs = true;
};

struct JcpChannelPrototypeMetrics {
  std::size_t panels = 0;
  std::size_t channels = 0;
  std::size_t near_pairs = 0;
  std::size_t far_pairs = 0;
  double setup_seconds = 0.0;
  double apply_seconds = 0.0;
  double direct_seconds = 0.0;
  double aggregate_rel_rms_error = 0.0;
  double max_channel_rel_rms_error = 0.0;
  double max_abs_error = 0.0;
  bool finite_outputs = true;
};

struct JcpK1BridgeMetrics {
  std::size_t panels = 0;
  std::size_t near_pairs = 0;
  std::size_t far_pairs = 0;
  double setup_seconds = 0.0;
  double apply_seconds = 0.0;
  double direct_seconds = 0.0;
  double rel_rms_error = 0.0;
  double rms_abs_error = 0.0;
  double max_abs_error = 0.0;
  bool finite_outputs = true;
};

struct JcpDerivativeBridgeMetrics {
  std::size_t panels = 0;
  std::size_t channels = 0;
  std::size_t near_pairs = 0;
  std::size_t far_pairs = 0;
  double setup_seconds = 0.0;
  double apply_seconds = 0.0;
  double direct_seconds = 0.0;
  double aggregate_rel_rms_error = 0.0;
  double max_channel_rel_rms_error = 0.0;
  double k2_rel_rms_error = 0.0;
  double k3_rel_rms_error = 0.0;
  double k4_rel_rms_error = 0.0;
  double k2_reference_rms = 0.0;
  double k3_reference_rms = 0.0;
  double k4_reference_rms = 0.0;
  double max_abs_error = 0.0;
  bool finite_outputs = true;
};

struct JcpDirectParityMetrics {
  std::size_t panels = 0;
  double phi_rel_rms_error = 0.0;
  double phin_rel_rms_error = 0.0;
  double aggregate_rel_rms_error = 0.0;
  double max_abs_error = 0.0;
  bool finite_outputs = true;
};

struct JcpOfflineOperatorResult {
  std::size_t panels = 0;
  std::size_t near_pairs = 0;
  std::size_t far_pairs = 0;
  double setup_seconds = 0.0;
  double apply_seconds = 0.0;
  std::vector<double> phi;
  std::vector<double> phin;
  bool finite_outputs = true;
};

struct JcpAnalyticOperatorResult {
  std::size_t panels = 0;
  std::size_t near_pairs = 0;
  std::size_t far_pairs = 0;
  double setup_seconds = 0.0;
  double apply_seconds = 0.0;
  std::size_t context_nodes = 0;
  std::size_t context_work_items = 0;
  std::vector<double> phi;
  std::vector<double> phin;
  bool finite_outputs = true;
};

struct JcpGuardedGmresMetrics {
  std::size_t panels = 0;
  std::size_t near_pairs = 0;
  std::size_t far_pairs = 0;
  int direct_iterations = 0;
  int adapter_iterations = 0;
  double direct_initial_residual = 0.0;
  double direct_final_residual = 0.0;
  double adapter_initial_residual = 0.0;
  double adapter_final_residual = 0.0;
  double direct_true_residual = 0.0;
  double solution_rel_error = 0.0;
  double direct_matvec_ms = 0.0;
  double adapter_matvec_ms = 0.0;
  double speedup_vs_direct_matvec = 0.0;
  double adapter_setup_ms = 0.0;
  double adapter_apply_ms = 0.0;
  bool direct_converged = false;
  bool adapter_converged = false;
  bool adapter_residual_decreased = false;
  bool finite_outputs = true;
};

struct PrototypePointCloud {
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> z;
  std::vector<double> q;
  std::vector<std::size_t> panel_id;
};

struct PrototypePanel {
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
  double nx = 0.0;
  double ny = 0.0;
  double nz = 1.0;
  double area = 1.0;
  double scalar_density = 1.0;
};

struct FmmM2lOverlapPanelGroup {
  std::vector<std::size_t> owned_panels;
  std::vector<std::size_t> factor_panels;
  std::size_t halo_count{0};
};

class FmmM2lAnalyticContext {
public:
  FmmM2lAnalyticContext(const PrototypeConfig &config,
                        std::vector<PrototypePanel> panels, double eps_ratio,
                        double kappa);
  ~FmmM2lAnalyticContext();

  FmmM2lAnalyticContext(FmmM2lAnalyticContext &&) noexcept;
  FmmM2lAnalyticContext &
  operator=(FmmM2lAnalyticContext &&) noexcept;
  FmmM2lAnalyticContext(const FmmM2lAnalyticContext &) = delete;
  FmmM2lAnalyticContext &
  operator=(const FmmM2lAnalyticContext &) = delete;

  [[nodiscard]] bool matches(const PrototypeConfig &config,
                             std::size_t panel_count, double eps_ratio,
                             double kappa) const noexcept;
  [[nodiscard]] std::size_t panel_count() const noexcept;
  [[nodiscard]] std::size_t context_nodes() const noexcept;
  [[nodiscard]] std::size_t context_work_items() const noexcept;
  [[nodiscard]] double setup_seconds() const noexcept;
  [[nodiscard]] std::vector<std::vector<std::size_t>>
  block_panel_ids(std::size_t max_block_size) const;
  [[nodiscard]] std::vector<FmmM2lOverlapPanelGroup>
  overlap_panel_groups(std::size_t owned_block_size,
                       std::size_t overlap_panel_cap) const;
  [[nodiscard]] std::vector<std::vector<std::size_t>>
  near_panel_graph(std::size_t neighbor_cap) const;

  JcpAnalyticOperatorResult apply(const std::vector<double> &phi,
                                  const std::vector<double> &phin,
                                  TitanExecutor &executor) const;

private:
  struct Impl;
  std::unique_ptr<Impl> impl_;
};

const char *kernel_kind_name(KernelKind kernel) noexcept;

PrototypeMetrics run_synthetic_probe(const PrototypeConfig &config,
                                     std::size_t panel_count, unsigned seed);

PrototypePointCloud make_mini_protein_like_point_cloud(std::size_t panel_count,
                                                       unsigned seed);

std::vector<PrototypePanel>
make_tabi_panel_like_fixture(std::size_t panel_count, unsigned seed);

PrototypePointCloud point_cloud_from_panels(const std::vector<PrototypePanel> &panels);

PrototypeMetrics run_point_cloud_probe(const PrototypeConfig &config,
                                       const PrototypePointCloud &points);

JcpChannelPrototypeMetrics
run_jcp_channel_probe(const PrototypeConfig &config,
                      const std::vector<PrototypePanel> &panels);

std::string classify_scope_d_decision(const JcpChannelPrototypeMetrics &metrics);

JcpK1BridgeMetrics
run_jcp_k1_bridge_probe(const PrototypeConfig &config,
                        const std::vector<PrototypePanel> &panels);

std::string classify_scope_e_decision(const JcpK1BridgeMetrics &metrics);

JcpDerivativeBridgeMetrics
run_jcp_derivative_bridge_probe(const PrototypeConfig &config,
                                const std::vector<PrototypePanel> &panels);

std::string
classify_scope_f_decision(const JcpDerivativeBridgeMetrics &metrics);

JcpDirectParityMetrics compare_jcp_direct_reference(
    const std::vector<PrototypePanel> &panels,
    const std::vector<double> &production_phi,
    const std::vector<double> &production_phin, double eps_ratio,
    double kappa);

JcpOfflineOperatorResult evaluate_jcp_offline_operator(
    const PrototypeConfig &config, const std::vector<PrototypePanel> &panels,
    double eps_ratio, double kappa);

JcpAnalyticOperatorResult evaluate_jcp_offline_operator_analytic(
    const PrototypeConfig &config, const std::vector<PrototypePanel> &panels,
    double eps_ratio, double kappa);

JcpAnalyticOperatorResult evaluate_jcp_offline_operator_analytic(
    const PrototypeConfig &config, const std::vector<PrototypePanel> &panels,
    const std::vector<double> &phi, const std::vector<double> &phin,
    double eps_ratio, double kappa);

JcpAnalyticOperatorResult evaluate_jcp_offline_operator_analytic(
    const PrototypeConfig &config, const std::vector<PrototypePanel> &panels,
    const std::vector<double> &phi, const std::vector<double> &phin,
    double eps_ratio, double kappa, TitanExecutor &executor);

JcpGuardedGmresMetrics run_jcp_fmm_m2l_guarded_gmres_probe(
    const PrototypeConfig &config, const std::vector<PrototypePanel> &panels,
    double eps_ratio, double kappa, int restart, int max_iterations,
    double tolerance);

std::string classify_scope_a_decision(const PrototypeMetrics &metrics,
                                      const PrototypeConfig &config);

} // namespace titan::tabi::fmm_m2l

#endif // TITAN_TABI_FMM_M2L_PROTOTYPE_H
