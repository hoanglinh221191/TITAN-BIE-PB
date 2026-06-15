#ifndef TITAN_TABI_FMM_M2L_PRECONDITIONER_H
#define TITAN_TABI_FMM_M2L_PRECONDITIONER_H

#include "tabi_fmm_m2l_prototype.h"
#include "tabi_preconditioner.h"
#include "tabi_types.h"

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

class TitanExecutor;

namespace TitanPBE::Tabi {

struct JcpSparseIluPreconditioner {
  struct Block2 {
    double a00{0.0};
    double a01{0.0};
    double a10{0.0};
    double a11{0.0};
  };

  std::size_t panel_count{0};
  bool coupled{true};
  std::vector<double> inv_diag;
  std::vector<std::size_t> row_ptr;
  std::vector<std::size_t> col_index;
  std::vector<Block2> values;
  std::vector<Block2> diag_inverse;
  std::size_t nnz_blocks{0};
  double avg_degree{0.0};
  std::size_t max_degree{0};
  std::size_t valid_blocks{0};
  std::size_t valid_panel_count{0};
  std::size_t ill_conditioned_blocks{0};
  std::string quality_status{"off"};

  void apply(const std::vector<double> &in, std::vector<double> &out) const;
  [[nodiscard]] std::size_t estimated_bytes() const noexcept;
};

std::size_t estimate_fmm_m2l_block_jacobi_bytes(
    const titan::tabi::fmm_m2l::FmmM2lAnalyticContext &context,
    std::size_t diagonal_size, bool coupled, std::size_t block_panel_cap);

std::shared_ptr<JcpBlockPreconditioner>
build_fmm_m2l_block_jacobi_preconditioner(
    const std::vector<Panel> &panels,
    const titan::tabi::fmm_m2l::FmmM2lAnalyticContext &context, bool coupled,
    double eps_ratio, double kappa, double lhs_phi, double lhs_phin,
    const std::vector<double> &diagonal, std::size_t block_panel_cap,
    TitanExecutor &executor, const std::string &near_integration = "centroid",
    const std::string &quad_buffer = "off");

std::size_t estimate_fmm_m2l_overlap_schwarz_bytes(
    const titan::tabi::fmm_m2l::FmmM2lAnalyticContext &context,
    std::size_t diagonal_size, bool coupled, std::size_t owned_panel_cap,
    std::size_t overlap_panel_cap);

std::shared_ptr<JcpBlockPreconditioner>
build_fmm_m2l_overlap_schwarz_preconditioner(
    const std::vector<Panel> &panels,
    const titan::tabi::fmm_m2l::FmmM2lAnalyticContext &context, bool coupled,
    double eps_ratio, double kappa, double lhs_phi, double lhs_phin,
    const std::vector<double> &diagonal, std::size_t owned_panel_cap,
    std::size_t overlap_panel_cap, TitanExecutor &executor,
    const std::string &near_integration = "centroid",
    const std::string &quad_buffer = "off");

std::size_t estimate_fmm_m2l_nearfield_ilu0_bytes(
    const titan::tabi::fmm_m2l::FmmM2lAnalyticContext &context,
    std::size_t diagonal_size, bool coupled, std::size_t neighbor_cap);

std::shared_ptr<JcpSparseIluPreconditioner>
build_fmm_m2l_nearfield_ilu0_preconditioner(
    const std::vector<Panel> &panels,
    const titan::tabi::fmm_m2l::FmmM2lAnalyticContext &context, bool coupled,
    double eps_ratio, double kappa, double lhs_phi, double lhs_phin,
    const std::vector<double> &diagonal, std::size_t neighbor_cap,
    double drop_tol, TitanExecutor &executor,
    const std::string &near_integration = "centroid",
    const std::string &quad_buffer = "off");

} // namespace TitanPBE::Tabi

#endif // TITAN_TABI_FMM_M2L_PRECONDITIONER_H
