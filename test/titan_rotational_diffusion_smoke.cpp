#include "titan_rotational_diffusion.h"

#include "titan_executor.h"

#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

namespace rotational = titan_rotational_diffusion;

constexpr double kPi = 3.14159265358979323846;

void require(bool condition, const std::string &message) {
  if (!condition) {
    std::cerr << "FAIL titan_rotational_diffusion_smoke: " << message << '\n';
    std::exit(1);
  }
}

void require_close(double actual, double expected, const std::string &message,
                   double tolerance = 1.0e-9) {
  if (!std::isfinite(actual) || std::fabs(actual - expected) > tolerance) {
    std::cerr << "FAIL titan_rotational_diffusion_smoke: " << message
              << " actual=" << actual << " expected=" << expected << '\n';
    std::exit(1);
  }
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

rotational::Point3 rotate_z_90(const rotational::Point3 &point) {
  return rotational::Point3{-point.y, point.x, point.z};
}

void require_same_covariance(const rotational::Covariance3 &actual,
                             const rotational::Covariance3 &expected,
                             const std::string &message) {
  for (std::size_t i = 0; i < actual.values.size(); ++i)
    require_close(actual.values[i], expected.values[i],
                  message + " lane " + std::to_string(i));
}

void run_rigid_triad_rotation_contract() {
  const std::array<rotational::Point3, 3> reference{
      rotational::Point3{0.0, 0.0, 0.0},
      rotational::Point3{2.0, 0.0, 0.0},
      rotational::Point3{0.0, 1.0, 0.0}};
  std::array<rotational::Point3, 3> rotated{};
  for (std::size_t i = 0; i < reference.size(); ++i)
    rotated[i] = rotate_z_90(reference[i]);

  const rotational::OrientationFit fit =
      rotational::fit_orientation_from_triad(reference, rotated, nullptr);
  require(fit.valid, "rigid triad rotation produces a stable orientation");
  require_close(rotational::angular_distance_rad(
                    rotational::Quaternion::identity(), fit.orientation),
                0.5 * kPi, "rigid triad recovers 90 degree rotation");
}

void run_hemisphere_contract() {
  const rotational::Quaternion q =
      rotational::normalized(rotational::Quaternion{0.5, 0.5, 0.5, 0.5});
  const rotational::Quaternion minus_q{-q.w, -q.x, -q.y, -q.z};
  const rotational::Quaternion canonical =
      rotational::canonicalize_hemisphere(minus_q, q);

  require_close(canonical.w, q.w, "q/-q canonicalization preserves w");
  require_close(canonical.x, q.x, "q/-q canonicalization preserves x");
  require_close(canonical.y, q.y, "q/-q canonicalization preserves y");
  require_close(canonical.z, q.z, "q/-q canonicalization preserves z");
  require_close(rotational::angular_distance_rad(q, minus_q), 0.0,
                "q/-q encode the same orientation");
}

void run_translation_and_pbc_invariance_contract() {
  const std::array<rotational::Point3, 3> reference{
      rotational::Point3{9.0, 9.0, 5.0},
      rotational::Point3{1.0, 9.0, 5.0},
      rotational::Point3{9.0, 1.0, 5.0}};
  const std::array<rotational::Point3, 3> translated{
      rotational::Point3{109.0, -41.0, 12.0},
      rotational::Point3{109.0, -39.0, 12.0},
      rotational::Point3{107.0, -41.0, 12.0}};
  const std::array<rotational::Point3, 3> wrapped{
      rotational::Point3{9.0, 9.0, 5.0},
      rotational::Point3{9.0, 1.0, 5.0},
      rotational::Point3{7.0, 9.0, 5.0}};
  const box frame_box = make_box();

  const rotational::OrientationFit translated_fit =
      rotational::fit_orientation_from_triad(reference, translated, &frame_box);
  const rotational::OrientationFit wrapped_fit =
      rotational::fit_orientation_from_triad(reference, wrapped, &frame_box);
  require(translated_fit.valid && wrapped_fit.valid,
          "translated and wrapped triads remain stable");
  require_close(rotational::angular_distance_rad(translated_fit.orientation,
                                                 wrapped_fit.orientation),
                0.0, "translation and PBC wrapping preserve orientation");
  require_close(rotational::angular_distance_rad(
                    rotational::Quaternion::identity(),
                    wrapped_fit.orientation),
                0.5 * kPi, "wrapped triad recovers 90 degree rotation");
}

void run_collinear_reject_contract() {
  const std::array<rotational::Point3, 3> reference{
      rotational::Point3{0.0, 0.0, 0.0},
      rotational::Point3{1.0, 0.0, 0.0},
      rotational::Point3{2.0, 0.0, 0.0}};
  const std::array<rotational::Point3, 3> frame = reference;

  require(!rotational::fit_orientation_from_triad(reference, frame, nullptr)
               .valid,
          "collinear triad fails closed");
}

void run_weighted_degeneracy_reject_contract() {
  const std::array<rotational::Point3, 3> reference{
      rotational::Point3{0.0, 0.0, 0.0},
      rotational::Point3{1.0, 0.0, 0.0},
      rotational::Point3{0.0, 1.0, 0.0}};
  const std::array<rotational::Point3, 3> frame = reference;
  const std::array<double, 3> masses{12.0, 12.0, 0.0};

  require(!rotational::fit_orientation(reference, frame, masses).valid,
          "zero-mass atoms cannot rescue a degenerate weighted body frame");
}

void run_invalid_box_reject_contract() {
  const std::array<rotational::Point3, 3> reference{
      rotational::Point3{0.0, 0.0, 0.0},
      rotational::Point3{1.0, 0.0, 0.0},
      rotational::Point3{0.0, 1.0, 0.0}};
  box invalid = make_box();
  invalid.set_gamma(0.0);
  require(!rotational::fit_orientation_from_triad(reference, reference,
                                                   &invalid)
               .valid,
          "degenerate periodic box fails closed before MIC reconstruction");
}

void run_covariance_tail_contract() {
  const std::vector<rotational::Point3> reference{
      {1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 10.0},
      {11.0, 13.0, 17.0}, {19.0, 23.0, 29.0}, {31.0, 37.0, 41.0},
      {43.0, 47.0, 53.0}};
  const std::vector<rotational::Point3> frame{
      {-2.0, 1.0, 4.0}, {-5.0, 4.0, 7.0}, {-8.0, 7.0, 11.0},
      {-13.0, 11.0, 18.0}, {-23.0, 19.0, 30.0}, {-37.0, 31.0, 42.0},
      {-47.0, 43.0, 54.0}};

  for (const std::size_t tail : {1U, 2U, 3U}) {
    const std::size_t count = 4 + tail;
    const rotational::Covariance3 scalar =
        rotational::covariance_scalar(reference.data(), frame.data(), count);
    const rotational::Covariance3 optimized =
        rotational::covariance_optimized(reference.data(), frame.data(), count);
    require_same_covariance(optimized, scalar,
                            "optimized covariance matches scalar tail " +
                                std::to_string(tail));
  }
}

void run_worker_determinism_contract() {
  const double s = std::sqrt(0.5);
  const std::vector<rotational::Quaternion> orientations{
      rotational::Quaternion::identity(),
      rotational::Quaternion{s, 0.0, 0.0, s},
      rotational::Quaternion{0.0, 0.0, 0.0, 1.0},
      rotational::Quaternion{-s, 0.0, 0.0, s}};
  const std::vector<double> time_ps{0.0, 1.0, 2.0, 3.0};
  TitanExecutor one_worker(1);
  TitanExecutor four_workers(4);

  const auto one = rotational::compute_lag_results(
      orientations, time_ps, 3, one_worker);
  const auto four = rotational::compute_lag_results(
      orientations, time_ps, 3, four_workers);
  const auto scalar =
      rotational::compute_lag_results_scalar(orientations, time_ps, 3);
  require(one.size() == 3 && four.size() == one.size(),
          "rotational lag kernel emits bounded rows");
  for (std::size_t i = 0; i < one.size(); ++i) {
    require(one[i].lag_index == four[i].lag_index,
            "worker count preserves lag order");
    require(one[i].samples == four[i].samples,
            "worker count preserves sample counts");
    require_close(one[i].lag_time_ps, four[i].lag_time_ps,
                  "worker count preserves lag times");
    require_close(one[i].rotational_msd_rad2, four[i].rotational_msd_rad2,
                  "worker count preserves rotational MSD");
    require_close(one[i].orientation_correlation,
                  four[i].orientation_correlation,
                  "worker count preserves orientation correlation");
    require_close(four[i].rotational_msd_rad2, scalar[i].rotational_msd_rad2,
                  "AVX2/scalar rotational MSD matches");
    require_close(four[i].c1, scalar[i].c1, "AVX2/scalar C1 matches");
    require_close(four[i].c2, scalar[i].c2, "AVX2/scalar C2 matches");
    require_same_covariance(four[i].covariance, scalar[i].covariance,
                            "AVX2/scalar quaternion covariance matches");
  }
  require_close(one[0].rotational_msd_rad2, 0.25 * kPi * kPi,
                "lag one rotational MSD handles hemisphere continuity");
}

void run_avx_quaternion_tail_contract() {
  std::vector<rotational::Quaternion> orientations;
  std::vector<double> time_ps;
  for (std::size_t frame = 0; frame < 8; ++frame) {
    const double angle = static_cast<double>(frame) * kPi / 12.0;
    orientations.push_back(
        rotational::Quaternion{std::cos(0.5 * angle), 0.0, 0.0,
                               std::sin(0.5 * angle)});
    time_ps.push_back(static_cast<double>(frame));
  }
  TitanExecutor four_workers(4);
  const auto optimized = rotational::compute_lag_results(
      orientations, time_ps, 3, four_workers);
  const auto scalar =
      rotational::compute_lag_results_scalar(orientations, time_ps, 3);
  require(optimized.size() == 3 && scalar.size() == optimized.size(),
          "AVX2 quaternion lag kernel emits rows for tail 3/2/1");
  for (std::size_t row = 0; row < optimized.size(); ++row) {
    require(optimized[row].samples == 7 - row,
            "AVX2 quaternion lag kernel keeps block plus scalar tail samples");
    require_close(optimized[row].rotational_msd_rad2,
                  scalar[row].rotational_msd_rad2,
                  "AVX2/scalar tail rotational MSD matches");
    require_close(optimized[row].c1, scalar[row].c1,
                  "AVX2/scalar tail C1 matches");
    require_close(optimized[row].c2, scalar[row].c2,
                  "AVX2/scalar tail C2 matches");
    require_same_covariance(optimized[row].covariance, scalar[row].covariance,
                            "AVX2/scalar quaternion covariance tail matches");
  }
}

void run_anisotropic_tensor_fit_contract() {
  const std::array<double, 3> expected_d{0.01, 0.02, 0.04};
  const std::array<double, 9> identity_axes{1.0, 0.0, 0.0, 0.0, 1.0,
                                            0.0, 0.0, 0.0, 1.0};
  std::vector<rotational::LagResult> rows;
  for (std::size_t lag = 1; lag <= 10; ++lag) {
    rotational::LagResult row;
    row.lag_index = lag;
    row.lag_time_ps = static_cast<double>(lag);
    row.samples = 1000;
    row.covariance = rotational::favro_covariance_model(
        row.lag_time_ps, expected_d, identity_axes);
    rows.push_back(row);
  }
  const rotational::TensorFit fit = rotational::fit_tensor(rows, "1,10");
  require(fit.ok, "Favro tensor fit converges for exact anisotropic model");
  for (std::size_t axis = 0; axis < expected_d.size(); ++axis)
    require_close(fit.d_per_ps[axis], expected_d[axis],
                  "Favro tensor fit recovers anisotropic diffusion axis " +
                      std::to_string(axis),
                      1.0e-6);
}

void run_impossible_tensor_fit_reject_contract() {
  std::vector<rotational::LagResult> rows;
  for (std::size_t lag = 1; lag <= 5; ++lag) {
    rotational::LagResult row;
    row.lag_index = lag;
    row.lag_time_ps = static_cast<double>(lag);
    row.samples = 1000;
    row.covariance.values = {0.9, 0.9, 0.9, 0.0, 0.0, 0.0};
    rows.push_back(row);
  }
  const rotational::TensorFit fit = rotational::fit_tensor(rows, "1,5");
  require(!fit.ok, "Favro fit rejects an impossible covariance series");
}

void run_invalid_timestamp_reject_contract() {
  const std::vector<rotational::Quaternion> orientations{
      rotational::Quaternion::identity(), rotational::Quaternion::identity(),
      rotational::Quaternion::identity()};
  const std::vector<double> time_ps{0.0, 1.0, 1.0};
  TitanExecutor executor(2);
  require(rotational::compute_lag_results(orientations, time_ps, 2, executor)
              .empty(),
          "rotational lag kernel rejects non-increasing timestamps");
}

} // namespace

int main() {
  run_rigid_triad_rotation_contract();
  run_hemisphere_contract();
  run_translation_and_pbc_invariance_contract();
  run_collinear_reject_contract();
  run_weighted_degeneracy_reject_contract();
  run_invalid_box_reject_contract();
  run_covariance_tail_contract();
  run_worker_determinism_contract();
  run_avx_quaternion_tail_contract();
  run_anisotropic_tensor_fit_contract();
  run_impossible_tensor_fit_reject_contract();
  run_invalid_timestamp_reject_contract();
  std::cout << "PASS titan_rotational_diffusion_smoke\n";
  return 0;
}
