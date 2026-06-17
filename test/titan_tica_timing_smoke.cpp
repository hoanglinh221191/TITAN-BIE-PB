#include "titan_kinetics.h"

#include <Eigen/Dense>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

namespace {

void require_true(bool condition, const std::string &message) {
  if (!condition) {
    std::cerr << "FAIL titan_tica_timing_smoke: " << message << "\n";
    std::exit(1);
  }
}

Eigen::VectorXd make_frame(int i) {
  const double t = 0.25 * static_cast<double>(i + 1);
  Eigen::VectorXd x(5);
  x << t, t * t, std::sin(t), std::cos(0.7 * t),
      static_cast<double>((i % 3) - 1);
  return x;
}

void add_frames(TitanKinetics &tica) {
  for (int i = 0; i < 14; ++i) {
    tica.accumulate_frame(make_frame(i));
  }
}

void assert_exact_dual_timing_and_sink() {
  TitanKinetics tica(5, 1, 2, 0.1, TicaBackend::ExactDual, 4);
  add_frames(tica);
  tica.solve_tica();
  require_true(tica.is_solved(), "exact_dual did not solve");

  const auto before_projection = tica.timing_stats();
  require_true(before_projection.cache_write_frames == tica.num_frames(),
               "exact_dual cache write frame counter mismatch");
  require_true(before_projection.exact_dual_svd_calls == 1,
               "exact_dual SVD call counter was not recorded");
  require_true(before_projection.reduced_solve_replay_frames ==
                   tica.num_frames(),
               "exact_dual reduced-solve replay counter mismatch");

  int64_t projected_rows = 0;
  double projection_checksum = 0.0;
  tica.project_all_frames_to_sink(
      [&](int64_t frame_index, const Eigen::Ref<const Eigen::VectorXd> &row) {
        require_true(frame_index == projected_rows,
                     "projection sink frame index is not sequential");
        require_true(row.size() == tica.n_components(),
                     "projection sink row has wrong component count");
        projection_checksum += row.sum();
        ++projected_rows;
      });
  require_true(projected_rows == tica.num_frames(),
               "projection sink did not receive every exact_dual frame");
  require_true(std::isfinite(projection_checksum),
               "projection sink checksum is not finite");

  const auto after_projection = tica.timing_stats();
  require_true(after_projection.projection_replay_frames == tica.num_frames(),
               "projection replay frame counter mismatch");
  require_true(after_projection.projection_replay_calls == 1,
               "projection replay call counter mismatch");

  tica.release_raw_feature_cache();
  require_true(tica.raw_feature_cache_rows() == 0,
               "exact_dual raw feature cache was not released");
}

void assert_randomized_timing_and_sink() {
  TitanKinetics tica(5, 1, 2, 0.1, TicaBackend::Randomized, 3);
  add_frames(tica);
  tica.solve_tica();
  require_true(tica.is_solved(), "randomized did not solve");

  const auto stats = tica.timing_stats();
  require_true(stats.spool_write_frames == tica.num_frames(),
               "randomized spool write frame counter mismatch");
  require_true(stats.randomized_sketch_read_frames == tica.num_frames(),
               "randomized sketch read counter mismatch");
  require_true(stats.randomized_compression_read_frames == tica.num_frames(),
               "randomized compression read counter mismatch");
  require_true(stats.reduced_solve_replay_frames == tica.num_frames(),
               "randomized reduced-solve replay counter mismatch");
  require_true(stats.randomized_svd_calls == 1,
               "randomized SVD call counter was not recorded");

  int64_t projected_rows = 0;
  tica.project_all_frames_to_sink(
      [&](int64_t, const Eigen::Ref<const Eigen::VectorXd> &row) {
        require_true(row.size() == tica.n_components(),
                     "randomized projection sink row has wrong component count");
        ++projected_rows;
      });
  require_true(projected_rows == tica.num_frames(),
               "projection sink did not receive every randomized frame");
}

} // namespace

int main() {
  assert_exact_dual_timing_and_sink();
  assert_randomized_timing_and_sink();
  std::cout << "PASS titan_tica_timing_smoke\n";
  return 0;
}
