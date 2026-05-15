#ifndef TITAN_PBE_TABI_NEAR_H
#define TITAN_PBE_TABI_NEAR_H

#include "tabi_types.h"

#include <cstddef>
#include <string>
#include <vector>

namespace TitanPBE::Tabi {

struct JcpNearAccum {
  double phi{0.0};
  double phin{0.0};
};

struct JcpNearPairKernel {
  double k1{0.0};
  double k2{0.0};
  double k3{0.0};
  double k4{0.0};
  bool self_pair{false};
  bool one_ring{false};
  bool geometric_buffer{false};
  bool used_duffy{false};
  bool used_centroid{false};
  bool valid{false};
  std::string policy{"off"};
};

std::string resolve_near_integration_mode(const std::string &requested,
                                          const std::vector<Panel> &panels);

std::string resolve_quad_buffer_mode(const std::string &requested,
                                     const std::string &near_integration,
                                     const std::vector<Panel> &panels);

struct GeometricNearDecision {
  bool active{false};
  bool adaptive{false};
  double gap{0.0};
  double gap_ratio{0.0};
  double center_ratio{0.0};
  double max_edge{0.0};
};

GeometricNearDecision classify_geometric_near_buffer_pair(
    const Panel &target, const Panel &source,
    const std::string &quad_buffer) noexcept;

JcpNearPairKernel evaluate_jcp_near_pair_kernel(
    const Panel &target, std::size_t target_index, const Panel &source,
    std::size_t source_index, double eps_ratio, double kappa, bool coupled,
    const std::string &near_integration, const std::string &quad_buffer,
    TreecodeStats *stats = nullptr);

JcpNearAccum accumulate_jcp_near_ordered(
    const Panel &target, std::size_t target_index,
    const std::vector<Panel> &panels, const std::vector<double> &unknowns,
    double eps_ratio, double kappa, bool coupled,
    const std::vector<std::size_t> &order, std::size_t begin,
    std::size_t count, const std::string &near_integration,
    const std::string &quad_buffer,
    TreecodeStats *stats = nullptr);

} // namespace TitanPBE::Tabi

#endif // TITAN_PBE_TABI_NEAR_H
