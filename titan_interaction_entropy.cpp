#include "titan_interaction_entropy.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <vector>

namespace TitanEntropy {
namespace {

bool compute_minus_tds(std::span<const double> values,
                       const InteractionEntropyOptions &options,
                       double &dimensionless, double &minus_tds,
                       std::size_t &filtered_samples) {
  if (values.empty() || options.rt_energy <= 0.0 ||
      !std::isfinite(options.rt_energy))
    return false;

  double avg = std::accumulate(values.begin(), values.end(), 0.0) /
               static_cast<double>(values.size());
  double sigma = 0.0;
  for (double e : values)
    sigma += (e - avg) * (e - avg);
  sigma = std::sqrt(sigma / static_cast<double>(values.size()));

  std::vector<double> scaled;
  scaled.reserve(values.size());
  for (double e : values) {
    if (sigma <= 1.0e-14 || options.filter_sigma <= 0.0 ||
        std::abs(e - avg) <= options.filter_sigma * sigma)
      scaled.push_back((e - avg) / options.rt_energy);
  }

  filtered_samples = scaled.size();
  if (scaled.empty())
    return false;

  double max_scaled = scaled[0];
  for (double e : scaled)
    max_scaled = std::max(max_scaled, e);

  double log_sum = 0.0;
  for (double e : scaled)
    log_sum += std::exp(e - max_scaled);

  dimensionless =
      max_scaled + std::log(log_sum) -
      std::log(static_cast<double>(scaled.size()));
  minus_tds = options.rt_energy * dimensionless;
  return true;
}

double estimate_lag1_ess(std::span<const double> values) {
  if (values.size() <= 1)
    return static_cast<double>(values.size());

  double avg = std::accumulate(values.begin(), values.end(), 0.0) /
               static_cast<double>(values.size());
  double sigma = 0.0;
  for (double e : values)
    sigma += (e - avg) * (e - avg);
  sigma = std::sqrt(sigma / static_cast<double>(values.size()));
  if (sigma <= 1.0e-14)
    return static_cast<double>(values.size());

  double lag_num = 0.0;
  for (std::size_t i = 1; i < values.size(); ++i)
    lag_num += (values[i - 1] - avg) * (values[i] - avg);
  const double rho =
      std::clamp(lag_num /
                     (static_cast<double>(values.size() - 1) * sigma * sigma),
                 -0.99, 0.99);
  return std::clamp(static_cast<double>(values.size()) * (1.0 - rho) /
                        (1.0 + rho),
                    1.0, static_cast<double>(values.size()));
}

} // namespace

InteractionEntropyResult
compute_interaction_entropy(std::span<const double> total_energy,
                            const InteractionEntropyOptions &options) {
  InteractionEntropyResult out;
  out.samples = total_energy.size();
  if (options.rt_energy <= 0.0 || !std::isfinite(options.rt_energy)) {
    out.status = "invalid_rt";
    return out;
  }
  if (total_energy.size() < options.min_samples) {
    out.status = "insufficient_samples";
    return out;
  }

  out.ess = estimate_lag1_ess(total_energy);
  if (!compute_minus_tds(total_energy, options, out.dimensionless,
                         out.minus_tds_energy, out.filtered_samples)) {
    out.status = "no_filtered_samples";
    return out;
  }
  out.status = "ok";

  if (options.block_size == 0)
    return out;

  out.block_size = options.block_size;
  std::vector<double> block_values;
  for (std::size_t begin = 0; begin < total_energy.size();
       begin += options.block_size) {
    const std::size_t end =
        std::min(total_energy.size(), begin + options.block_size);
    std::span<const double> block(total_energy.data() + begin, end - begin);
    if (block.size() < options.min_samples)
      continue;
    double block_dim = 0.0;
    double block_minus_tds = 0.0;
    std::size_t block_filtered = 0;
    if (compute_minus_tds(block, options, block_dim, block_minus_tds,
                          block_filtered))
      block_values.push_back(block_minus_tds);
  }

  out.block_count = block_values.size();
  if (!block_values.empty()) {
    out.block_mean_energy =
        std::accumulate(block_values.begin(), block_values.end(), 0.0) /
        static_cast<double>(block_values.size());
    double block_var = 0.0;
    for (double v : block_values)
      block_var += (v - out.block_mean_energy) * (v - out.block_mean_energy);
    if (block_values.size() > 1)
      out.block_std_energy =
          std::sqrt(block_var / static_cast<double>(block_values.size() - 1));
    out.block_drift_energy = block_values.back() - block_values.front();
  }
  return out;
}

InteractionEntropyResult
compute_interaction_entropy(std::span<const double> elec_energy,
                            std::span<const double> vdw_energy,
                            const InteractionEntropyOptions &options) {
  const std::size_t n = std::min(elec_energy.size(), vdw_energy.size());
  std::vector<double> total;
  total.reserve(n);
  for (std::size_t i = 0; i < n; ++i)
    total.push_back(elec_energy[i] + vdw_energy[i]);
  return compute_interaction_entropy(std::span<const double>(total.data(), n),
                                     options);
}

} // namespace TitanEntropy
