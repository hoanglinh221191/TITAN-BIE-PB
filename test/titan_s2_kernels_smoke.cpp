#include "titan_s2_kernels.h"

#include "titan_executor.h"

#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

namespace s2 = titan_s2_kernel;

void require(bool condition, const std::string &message) {
  if (!condition) {
    std::cerr << "FAIL titan_s2_kernels_smoke: " << message << '\n';
    std::exit(1);
  }
}

void require_close(double actual, double expected, const std::string &message,
                   double tolerance = 1.0e-9) {
  if (!std::isfinite(actual) || std::fabs(actual - expected) > tolerance) {
    std::cerr << "FAIL titan_s2_kernels_smoke: " << message
              << " actual=" << actual << " expected=" << expected << '\n';
    std::exit(1);
  }
}

template <typename Fn>
void require_invalid_argument(Fn &&fn, const std::string &message) {
  bool rejected = false;
  try {
    fn();
  } catch (const std::invalid_argument &) {
    rejected = true;
  }
  require(rejected, message);
}

template <typename Fn>
void require_out_of_range(Fn &&fn, const std::string &message) {
  bool rejected = false;
  try {
    fn();
  } catch (const std::out_of_range &) {
    rejected = true;
  }
  require(rejected, message);
}

box make_box(double side = 10.0) {
  box frame_box;
  frame_box.set_a(side);
  frame_box.set_b(side);
  frame_box.set_c(side);
  frame_box.set_alpha(90.0);
  frame_box.set_beta(90.0);
  frame_box.set_gamma(90.0);
  return frame_box;
}

s2::BondVectorsSoA make_vectors(std::size_t count) {
  s2::BondVectorsSoA vectors;
  vectors.x.reserve(count);
  vectors.y.reserve(count);
  vectors.z.reserve(count);
  vectors.valid.reserve(count);
  for (std::size_t i = 0; i < count; ++i) {
    const double angle = 0.17 * static_cast<double>(i);
    vectors.x.push_back(std::cos(angle));
    vectors.y.push_back(std::sin(angle));
    vectors.z.push_back(0.0);
    vectors.valid.push_back(1);
  }
  return vectors;
}

void run_rigid_vector_contract() {
  s2::BondVectorsSoA vectors;
  for (std::size_t i = 0; i < 8; ++i) {
    vectors.x.push_back(1.0);
    vectors.y.push_back(0.0);
    vectors.z.push_back(0.0);
    vectors.valid.push_back(1);
  }
  std::vector<double> time_ps(vectors.size(), 0.0);
  for (std::size_t i = 0; i < time_ps.size(); ++i)
    time_ps[i] = static_cast<double>(i);

  TitanExecutor executor(4);
  const auto rows = s2::compute_correlation(vectors, time_ps, 3, executor);
  require(rows.size() == 3, "rigid vector emits requested lags");
  for (const auto &row : rows) {
    require_close(row.correlation, 1.0, "rigid vector has P2 correlation one");
    require(row.samples == vectors.size() - row.lag_index,
            "rigid vector sample count matches lag");
  }
}

void run_synthetic_fit_contract() {
  constexpr double expected_s2 = 0.35;
  constexpr double expected_tau_ps = 4.5;
  std::vector<s2::CorrelationRow> rows;
  for (std::size_t lag = 1; lag <= 16; ++lag) {
    const double time_ps = static_cast<double>(lag);
    rows.push_back(s2::CorrelationRow{
        lag, time_ps,
        expected_s2 + (1.0 - expected_s2) * std::exp(-time_ps / expected_tau_ps),
        1000 - 7 * lag});
  }

  const s2::FitResult fit = s2::fit_single_exp_plateau(rows);
  require(fit.ok, "synthetic exponential fit succeeds");
  require(!fit.used_tail_fallback, "synthetic exponential uses scan fit");
  require_close(fit.raw_plateau, expected_s2,
                "synthetic exponential recovers raw plateau", 1.0e-6);
  require_close(fit.s2, expected_s2, "synthetic exponential recovers S2",
                1.0e-6);
  require_close(fit.tau_ps, expected_tau_ps,
                "synthetic exponential recovers tau", 1.0e-5);
}

void run_raw_clamp_and_short_tail_contract() {
  const std::vector<s2::CorrelationRow> rows{
      {1, 1.0, 1.2, 2}, {2, 2.0, 1.4, 3}, {3, 3.0, 1.6, 5}};
  const s2::FitResult fit = s2::fit_single_exp_plateau(rows);
  require(fit.ok, "short tail fallback succeeds");
  require(fit.used_tail_fallback, "fewer than five lags use tail fallback");
  require(fit.tail_rows == 3, "fallback uses at least three tail rows");
  require_close(fit.raw_plateau, 1.46, "raw plateau remains inspectable");
  require_close(fit.s2, 1.0, "reported S2 clamps raw plateau");
  require(fit.tau_ps > 0.0, "fallback tau remains positive");
}

void run_mic_and_raw_contract() {
  const std::vector<double> from_x{9.5};
  const std::vector<double> from_y{0.0};
  const std::vector<double> from_z{0.0};
  const std::vector<double> to_x{0.5};
  const std::vector<double> to_y{0.0};
  const std::vector<double> to_z{0.0};
  const box valid_box = make_box();
  box invalid_box;
  s2::BondVectorsSoA mic;
  s2::BondVectorsSoA raw;

  s2::build_bond_vectors(from_x, from_y, from_z, to_x, to_y, to_z,
                         &valid_box, mic);
  s2::build_bond_vectors(from_x, from_y, from_z, to_x, to_y, to_z,
                         &invalid_box, raw);
  require(mic.size() == 1 && raw.size() == 1,
          "bond vector builder emits one row");
  require_close(mic.x[0], 1.0, "valid box applies MIC");
  require_close(raw.x[0], -1.0,
                "intentional invalid-box policy falls back to raw delta");
}

void run_mismatched_builder_reject_contract() {
  const std::vector<double> one{0.0};
  const std::vector<double> two{0.0, 1.0};
  s2::BondVectorsSoA vectors;
  require_invalid_argument([&]() {
    s2::build_bond_vectors(one, one, one, two, one, one, nullptr, vectors);
  }, "bond vector builder rejects mismatched coordinate spans");
}

void run_zero_sample_fit_reject_contract() {
  require(std::isnan(s2::CorrelationRow{}.correlation),
          "default correlation row is invalid until samples exist");
  s2::BondVectorsSoA invalid_vectors;
  for (std::size_t i = 0; i < 4; ++i)
    invalid_vectors.push(0.0, 0.0, 0.0, false);
  const s2::CorrelationRow scalar =
      s2::compute_correlation_lag_scalar(invalid_vectors, 1);
  const s2::CorrelationRow optimized =
      s2::compute_correlation_lag_optimized(invalid_vectors, 1);
  require(scalar.samples == 0 && optimized.samples == 0,
          "invalid bonds produce zero correlation samples");
  require(std::isnan(scalar.correlation) && std::isnan(optimized.correlation),
          "zero-sample lag correlation is nan");

  const std::vector<s2::CorrelationRow> rows{
      {1, 1.0, std::numeric_limits<double>::quiet_NaN(), 0},
      {2, 2.0, 0.7, 5},
      {3, 3.0, 0.6, 5}};
  require(!s2::fit_single_exp_plateau(rows).ok,
          "fewer than three usable lag rows cannot produce S2");
}

void run_zero_sample_fallback_skip_contract() {
  const std::vector<s2::CorrelationRow> rows{
      {1, 1.0, 99.0, 0}, {2, 2.0, 0.2, 2}, {3, 3.0, 0.4, 3},
      {4, 4.0, 0.6, 5}};
  const s2::FitResult fit = s2::fit_single_exp_plateau(rows);
  require(fit.ok && fit.used_tail_fallback,
          "three usable lag rows use tail fallback");
  require(fit.tail_rows == 3, "tail fallback counts usable rows only");
  require_close(fit.raw_plateau, 0.46,
                "tail fallback ignores zero-sample rows");
}

void run_invalid_fit_time_order_reject_contract() {
  const std::vector<s2::CorrelationRow> rows{
      {1, 1.0, 0.9, 10}, {2, 2.0, 0.8, 10}, {3, 2.0, 0.7, 10},
      {4, 4.0, 0.6, 10}, {5, 5.0, 0.5, 10}};
  const s2::FitResult fit = s2::fit_single_exp_plateau(rows);
  require(!fit.ok, "fit rejects non-increasing usable lag times");
  require(!fit.used_tail_fallback,
          "invalid lag ordering fails closed instead of hiding misuse");
}

void run_kernel_input_reject_contract() {
  s2::BondVectorsSoA malformed = make_vectors(4);
  malformed.y.pop_back();
  const std::vector<double> time_ps{0.0, 1.0, 2.0, 3.0};
  TitanExecutor executor(2);
  require_invalid_argument(
      [&]() { static_cast<void>(s2::compute_correlation_lag_scalar(malformed, 1)); },
      "scalar lag rejects malformed SoA");
  require_invalid_argument(
      [&]() { static_cast<void>(s2::compute_correlation_lag_optimized(malformed, 1)); },
      "optimized lag rejects malformed SoA");
  require_invalid_argument(
      [&]() { static_cast<void>(s2::compute_correlation_scalar(malformed, time_ps, 2)); },
      "scalar correlation rejects malformed SoA");
  require_invalid_argument(
      [&]() { static_cast<void>(s2::compute_correlation(malformed, time_ps, 2, executor)); },
      "parallel correlation rejects malformed SoA");
  const std::array<double, 9> identity{
      1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
  require_invalid_argument(
      [&]() { s2::rotate_bond_vectors_inplace_scalar(malformed, identity); },
      "scalar rotation rejects malformed SoA");
  require_invalid_argument(
      [&]() { s2::rotate_bond_vectors_inplace_optimized(malformed, identity); },
      "optimized rotation rejects malformed SoA");

  const s2::BondVectorsSoA vectors = make_vectors(4);
  const std::vector<double> repeated_time{0.0, 1.0, 1.0, 3.0};
  const std::vector<double> short_time{0.0, 1.0, 2.0};
  require_invalid_argument(
      [&]() { static_cast<void>(s2::compute_correlation_scalar(vectors, repeated_time, 2)); },
      "scalar correlation rejects invalid time axis");
  require_invalid_argument(
      [&]() { static_cast<void>(s2::compute_correlation(vectors, repeated_time, 2, executor)); },
      "parallel correlation rejects invalid time axis");
  require_invalid_argument(
      [&]() { static_cast<void>(s2::compute_correlation_scalar(vectors, short_time, 2)); },
      "scalar correlation rejects mismatched time axis");
  require_invalid_argument(
      [&]() { static_cast<void>(s2::compute_correlation(vectors, short_time, 2, executor)); },
      "parallel correlation rejects mismatched time axis");
  require_out_of_range(
      [&]() { static_cast<void>(s2::compute_correlation_lag_scalar(vectors, 0)); },
      "scalar lag rejects zero");
  require_out_of_range(
      [&]() { static_cast<void>(s2::compute_correlation_lag_optimized(vectors, 0)); },
      "optimized lag rejects zero");
  require_out_of_range(
      [&]() { static_cast<void>(s2::compute_correlation_lag_scalar(vectors, vectors.size())); },
      "scalar lag rejects upper bound");
  require_out_of_range(
      [&]() { static_cast<void>(s2::compute_correlation_lag_optimized(vectors, vectors.size())); },
      "optimized lag rejects upper bound");
  require(s2::compute_correlation_scalar(vectors, time_ps, 0).empty(),
          "scalar maxlag zero remains valid empty");
  require(s2::compute_correlation(vectors, time_ps, 0, executor).empty(),
          "parallel maxlag zero remains valid empty");
}

void require_same_vectors(const s2::BondVectorsSoA &actual,
                          const s2::BondVectorsSoA &expected,
                          const std::string &message) {
  require(actual.size() == expected.size(), message + " size");
  for (std::size_t i = 0; i < actual.size(); ++i) {
    require(actual.valid[i] == expected.valid[i],
            message + " valid " + std::to_string(i));
    require_close(actual.x[i], expected.x[i],
                  message + " x " + std::to_string(i), 1.0e-12);
    require_close(actual.y[i], expected.y[i],
                  message + " y " + std::to_string(i), 1.0e-12);
    require_close(actual.z[i], expected.z[i],
                  message + " z " + std::to_string(i), 1.0e-12);
  }
}

void run_builder_scalar_optimized_tail_contract() {
  const box frame_box = make_box();
  box invalid_box;
  for (const std::size_t tail : {1U, 2U, 3U}) {
    const std::size_t count = 4 + tail;
    std::vector<double> from_x(count), from_y(count), from_z(count);
    std::vector<double> to_x(count), to_y(count), to_z(count);
    for (std::size_t i = 0; i < count; ++i) {
      from_x[i] = i % 2 == 0 ? 9.4 : 1.0 + static_cast<double>(i);
      from_y[i] = 0.3 * static_cast<double>(i);
      from_z[i] = -0.2 * static_cast<double>(i);
      to_x[i] = i % 2 == 0 ? 0.6 : 2.5 + static_cast<double>(i);
      to_y[i] = 1.1 + 0.1 * static_cast<double>(i);
      to_z[i] = 0.7 - 0.05 * static_cast<double>(i);
    }
    s2::BondVectorsSoA scalar_mic;
    s2::BondVectorsSoA optimized_mic;
    s2::BondVectorsSoA scalar_raw;
    s2::BondVectorsSoA optimized_raw;
    s2::build_bond_vectors_scalar(from_x, from_y, from_z, to_x, to_y, to_z,
                                  &frame_box, scalar_mic);
    s2::build_bond_vectors_optimized(from_x, from_y, from_z, to_x, to_y, to_z,
                                     &frame_box, optimized_mic);
    s2::build_bond_vectors_scalar(from_x, from_y, from_z, to_x, to_y, to_z,
                                  &invalid_box, scalar_raw);
    s2::build_bond_vectors_optimized(from_x, from_y, from_z, to_x, to_y, to_z,
                                     nullptr, optimized_raw);
    require_same_vectors(optimized_mic, scalar_mic,
                         "SIMD MIC builder parity tail " + std::to_string(tail));
    require_same_vectors(
        optimized_raw, scalar_raw,
        "SIMD intentional raw-policy builder parity tail " +
            std::to_string(tail));
  }
}

void run_half_box_mic_rounding_parity_contract() {
  const box frame_box = make_box();
  const std::vector<double> from_x(4, 0.0);
  const std::vector<double> from_y(4, 0.0);
  const std::vector<double> from_z(4, 0.0);
  const std::vector<double> to_x{5.0, -5.0, 15.0, -15.0};
  const std::vector<double> to_y(4, 0.0);
  const std::vector<double> to_z(4, 0.0);
  s2::BondVectorsSoA scalar;
  s2::BondVectorsSoA optimized;
  s2::build_bond_vectors_scalar(from_x, from_y, from_z, to_x, to_y, to_z,
                                &frame_box, scalar);
  s2::build_bond_vectors_optimized(from_x, from_y, from_z, to_x, to_y, to_z,
                                   &frame_box, optimized);
  require_same_vectors(optimized, scalar,
                       "SIMD MIC half-away-from-zero rounding parity");
}

void run_rotate_scalar_optimized_tail_contract() {
  const std::array<double, 9> rotate_z_90{
      0.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0};
  for (const std::size_t tail : {1U, 2U, 3U}) {
    s2::BondVectorsSoA scalar = make_vectors(4 + tail);
    s2::BondVectorsSoA optimized = scalar;
    s2::BondVectorsSoA dispatched = scalar;
    s2::rotate_bond_vectors_inplace_scalar(scalar, rotate_z_90);
    s2::rotate_bond_vectors_inplace_optimized(optimized, rotate_z_90);
    s2::rotate_bond_vectors_inplace(dispatched, rotate_z_90);
    require_same_vectors(optimized, scalar,
                         "SIMD rotation parity tail " + std::to_string(tail));
    require_same_vectors(dispatched, scalar,
                         "runtime rotation parity tail " + std::to_string(tail));
  }
}

void run_scalar_optimized_tail_contract() {
  for (const std::size_t tail : {1U, 2U, 3U}) {
    const s2::BondVectorsSoA vectors = make_vectors(8 + tail);
    const std::size_t lag = 4;
    const s2::CorrelationRow scalar =
        s2::compute_correlation_lag_scalar(vectors, lag);
    const s2::CorrelationRow optimized =
        s2::compute_correlation_lag_optimized(vectors, lag);
    require(scalar.samples == optimized.samples,
            "scalar/optimized tail sample counts match");
    require_close(optimized.correlation, scalar.correlation,
                  "scalar/optimized P2 tail " + std::to_string(tail),
                  1.0e-12);
  }
}

void run_worker_determinism_contract() {
  const s2::BondVectorsSoA vectors = make_vectors(24);
  std::vector<double> time_ps(vectors.size(), 0.0);
  for (std::size_t i = 0; i < time_ps.size(); ++i)
    time_ps[i] = 0.5 * static_cast<double>(i);
  TitanExecutor one_worker(1);
  TitanExecutor four_workers(4);

  const auto one = s2::compute_correlation(vectors, time_ps, 7, one_worker);
  const auto four = s2::compute_correlation(vectors, time_ps, 7, four_workers);
  require(one.size() == four.size(), "worker result sizes match");
  for (std::size_t i = 0; i < one.size(); ++i) {
    require(one[i].lag_index == four[i].lag_index,
            "worker lag order is deterministic");
    require(one[i].samples == four[i].samples,
            "worker sample count is deterministic");
    require_close(one[i].lag_time_ps, four[i].lag_time_ps,
                  "worker lag time is deterministic");
    require_close(one[i].correlation, four[i].correlation,
                  "worker correlation is deterministic");
  }
}

} // namespace

int main() {
  run_rigid_vector_contract();
  run_synthetic_fit_contract();
  run_raw_clamp_and_short_tail_contract();
  run_mic_and_raw_contract();
  run_mismatched_builder_reject_contract();
  run_zero_sample_fit_reject_contract();
  run_zero_sample_fallback_skip_contract();
  run_invalid_fit_time_order_reject_contract();
  run_kernel_input_reject_contract();
  run_builder_scalar_optimized_tail_contract();
  run_half_box_mic_rounding_parity_contract();
  run_rotate_scalar_optimized_tail_contract();
  run_scalar_optimized_tail_contract();
  run_worker_determinism_contract();
  std::cout << "PASS titan_s2_kernels_smoke\n";
  return 0;
}
