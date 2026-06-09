#include "titan_uncertainty.h"

#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

namespace {

bool near(double a, double b, double tol = 1.0e-12) {
  return std::fabs(a - b) <= tol;
}

void require(bool ok, const char *message) {
  if (!ok) {
    std::cerr << "FAIL: " << message << "\n";
    std::exit(1);
  }
}

} // namespace

int main() {
  using namespace titan_uncertainty;

  const auto small = partition_blocks(9, "auto");
  require(small.empty(), "frames < 10 should be insufficient for auto blocks");

  const auto mid = partition_blocks(25, "auto");
  require(mid.size() == 5, "25 frames should use 5 auto blocks");
  require(mid.front().begin == 0 && mid.back().end == 25,
          "auto blocks should cover the full range");

  const auto explicit_blocks = partition_blocks(10, "3");
  require(explicit_blocks.size() == 3, "explicit block count should parse");
  require(explicit_blocks[0].begin == 0 && explicit_blocks[0].end == 4,
          "block partition should distribute remainder deterministically");
  require(explicit_blocks[2].begin == 7 && explicit_blocks[2].end == 10,
          "last block range should be deterministic");

  const Stats stats = summarize_values({1.0, 2.0, 3.0, 4.0});
  require(stats.n == 4, "stats should count finite values");
  require(near(stats.mean, 2.5), "mean should match");
  require(near(stats.sem, std::sqrt((5.0 / 3.0) / 4.0)),
          "SEM should use sample variance");
  require(std::isfinite(stats.ci95_low) && std::isfinite(stats.ci95_high),
          "95% CI should be finite for n >= 2");

  const Stats one = summarize_values({42.0});
  require(one.n == 1 && std::isnan(one.sem),
          "one-sample SEM should be nan");

  const auto ci_a = bootstrap_ci95({1.0, 2.0, 3.0, 4.0}, 80, 7);
  const auto ci_b = bootstrap_ci95({1.0, 2.0, 3.0, 4.0}, 80, 7);
  require(near(ci_a.first, ci_b.first) && near(ci_a.second, ci_b.second),
          "bootstrap CI should be deterministic for the same seed");
  require(ci_a.first <= ci_a.second, "bootstrap CI should be ordered");

  const auto boot_stat_a =
      bootstrap_statistic({1.0, 2.0, 3.0, 4.0}, 80, 7);
  const auto boot_stat_b =
      bootstrap_statistic({1.0, 2.0, 3.0, 4.0}, 80, 7);
  require(near(boot_stat_a.ci95_low, boot_stat_b.ci95_low) &&
              near(boot_stat_a.ci95_high, boot_stat_b.ci95_high),
          "bootstrap statistic should be deterministic");

  std::vector<XYPoint> linear_points;
  for (int i = 1; i <= 20; ++i)
    linear_points.push_back(XYPoint{static_cast<double>(i),
                                    1.5 + 2.0 * static_cast<double>(i)});
  const auto linear_fit =
      linear_fit_with_blocks(linear_points, 1.0, 20.0, "5", 80, 11);
  require(linear_fit.point.ok, "linear block fit should be ok");
  require(near(linear_fit.point.slope, 2.0, 1.0e-12),
          "linear fit should recover slope");
  require(near(linear_fit.point.intercept, 1.5, 1.0e-12),
          "linear fit should recover intercept");
  require(linear_fit.slope_stats.n == 5 &&
              near(linear_fit.slope_stats.mean, 2.0),
          "linear block slopes should summarize");
  require(std::isfinite(linear_fit.slope_bootstrap_ci95_low) &&
              std::isfinite(linear_fit.slope_bootstrap_ci95_high),
          "linear fit bootstrap CI should be finite");

  std::vector<XYPoint> log_points;
  for (int i = 1; i <= 20; ++i) {
    const double t = static_cast<double>(i);
    log_points.push_back(XYPoint{t, 3.0 * std::pow(t, 1.25)});
  }
  const auto log_fit = loglog_fit_with_blocks(log_points, 1.0, 20.0, "5", 80, 13);
  require(log_fit.point.ok, "log-log block fit should be ok");
  require(near(log_fit.point.slope, 1.25, 1.0e-12),
          "log-log fit should recover alpha");

  const auto insufficient =
      linear_fit_with_blocks({XYPoint{1.0, 2.0}, XYPoint{2.0, 4.0}},
                             1.0, 2.0, "auto", 20, 3);
  require(!insufficient.point.ok &&
              insufficient.point.status == "insufficient_points",
          "too few fit rows should be insufficient");

  bool rejected = false;
  try {
    (void)partition_blocks(50, "0");
  } catch (const std::exception &) {
    rejected = true;
  }
  require(rejected, "invalid explicit block count should reject");

  std::cout << "PASS titan_uncertainty_smoke\n";
  return 0;
}
