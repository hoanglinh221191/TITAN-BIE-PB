#ifndef TITAN_TABI_FMM_M2L_LOCAL_H
#define TITAN_TABI_FMM_M2L_LOCAL_H

#include "tabi_fmm_m2l_cartesian.h"

#include <array>
#include <cstddef>

namespace titan::tabi::fmm_m2l {

inline constexpr int kCartesianLocalSourceOrder = 4;
inline constexpr int kCartesianLocalTargetOrder = 4;
inline constexpr int kCartesianLocalKernelOrder = 8;
inline constexpr std::size_t kCartesianLocalTermMax = 35;
inline constexpr std::size_t kCartesianKernelJetTermCount = 165;
inline constexpr std::size_t kCartesianJcpKernelJetCount = 8;

enum class JcpLocalField : std::size_t {
  K1 = 0,
  K2 = 1,
  K3Potential = 2,
  K4Potential = 3,
};

struct CartesianJcpLocalP4 {
  std::array<std::array<double, kCartesianLocalTermMax>, 4> fields{};
  bool finite_outputs{true};
};

struct CartesianJcpM2lOperatorP4 {
  std::array<std::array<double, kCartesianKernelJetTermCount>,
             kCartesianJcpKernelJetCount>
      jets{};
  bool finite_outputs{true};
};

std::size_t cartesian_local_kernel_term_count() noexcept;
std::size_t cartesian_local_term_count() noexcept;

CartesianJcpLocalP4 translate_cartesian_jcp_m2l_p4(
    const CartesianJcpMoments &moments, double source_center_x,
    double source_center_y, double source_center_z, double target_center_x,
    double target_center_y, double target_center_z, double eps_ratio,
    double kappa);

CartesianJcpM2lOperatorP4 build_cartesian_jcp_m2l_operator_p4(
    double source_center_x, double source_center_y, double source_center_z,
    double target_center_x, double target_center_y, double target_center_z,
    double eps_ratio, double kappa);

CartesianJcpLocalP4 apply_cartesian_jcp_m2l_operator_p4(
    const CartesianJcpMoments &moments,
    const CartesianJcpM2lOperatorP4 &operator_p4);

void accumulate_cartesian_jcp_local(
    CartesianJcpLocalP4 &destination,
    const CartesianJcpLocalP4 &source) noexcept;

double evaluate_cartesian_jcp_local_scalar(
    const CartesianJcpLocalP4 &local, JcpLocalField field, double dx,
    double dy, double dz) noexcept;

std::array<double, 4> evaluate_cartesian_jcp_l2p(
    const PrototypePanel &target, double target_center_x,
    double target_center_y, double target_center_z,
    const CartesianJcpLocalP4 &local) noexcept;

std::array<double, 4> evaluate_cartesian_jcp_l2p_order(
    const PrototypePanel &target, double target_center_x,
    double target_center_y, double target_center_z,
    const CartesianJcpLocalP4 &local, int target_order) noexcept;

} // namespace titan::tabi::fmm_m2l

#endif // TITAN_TABI_FMM_M2L_LOCAL_H
