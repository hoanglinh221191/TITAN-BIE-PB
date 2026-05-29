#include "titan_contact_kernels.h"

#include "BoxMatrix.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

namespace {

box make_box() {
  box frame_box;
  frame_box.set_a(10.0);
  frame_box.set_b(10.0);
  frame_box.set_c(10.0);
  frame_box.set_alpha(90.0);
  frame_box.set_beta(90.0);
  frame_box.set_gamma(90.0);
  return frame_box;
}

void require(bool condition, const std::string &message) {
  if (!condition) {
    std::cerr << "FAIL titan_contact_kernels_smoke: " << message << '\n';
    std::exit(1);
  }
}

void require_close(double actual, double expected, const std::string &message) {
  if (std::fabs(actual - expected) > 1.0e-9) {
    std::cerr << "FAIL titan_contact_kernels_smoke: " << message
              << " actual=" << actual << " expected=" << expected << '\n';
    std::exit(1);
  }
}

void sort_candidates(std::vector<titan_contact_kernel::ContactCandidate> &items) {
  std::sort(items.begin(), items.end(), [](const auto &lhs, const auto &rhs) {
    if (lhs.row != rhs.row)
      return lhs.row < rhs.row;
    return lhs.col < rhs.col;
  });
}

void require_same_candidates(
    std::vector<titan_contact_kernel::ContactCandidate> scalar,
    std::vector<titan_contact_kernel::ContactCandidate> optimized) {
  sort_candidates(scalar);
  sort_candidates(optimized);
  require(scalar.size() == optimized.size(),
          "scalar and optimized candidate counts differ");
  for (std::size_t i = 0; i < scalar.size(); ++i) {
    require(scalar[i].row == optimized[i].row &&
                scalar[i].col == optimized[i].col,
            "scalar and optimized candidate indices differ");
    require_close(optimized[i].distance_a, scalar[i].distance_a,
                  "optimized distance matches scalar");
  }
}

void run_interface_tail_and_mic_contract() {
  const BoxMatrix box_matrix(make_box());
  titan_contact_kernel::PointSoA rows;
  rows.x = {1.0, 5.0};
  rows.y = {1.0, 5.0};
  rows.z = {1.0, 5.0};

  titan_contact_kernel::PointSoA cols;
  cols.x = {9.0, 1.5, 5.5, 8.1, 7.7, 0.2, 4.0};
  cols.y = {1.0, 1.0, 5.0, 8.1, 7.7, 1.0, 4.0};
  cols.z = {1.0, 1.0, 5.0, 8.1, 7.7, 1.0, 4.0};

  std::vector<titan_contact_kernel::ContactCandidate> scalar;
  std::vector<titan_contact_kernel::ContactCandidate> optimized;
  titan_contact_kernel::scan_interface_contacts_scalar(
      rows, cols, box_matrix, 2.25 * 2.25, 0, rows.size(), scalar);
  titan_contact_kernel::scan_interface_contacts_optimized(
      rows, cols, box_matrix, 2.25 * 2.25, 0, rows.size(), optimized);

  require_same_candidates(scalar, optimized);
  require(scalar.size() == 5, "MIC/tail scan finds expected contacts");
  require_close(scalar.front().distance_a, 2.0,
                "orthorhombic MIC wraps row 0 to col 0");
}

void run_mol1_upper_triangle_contract() {
  const BoxMatrix box_matrix(make_box());
  titan_contact_kernel::PointSoA points;
  points.x = {1.0, 2.0, 9.8, 5.0, 5.9};
  points.y = {1.0, 1.0, 1.0, 5.0, 5.0};
  points.z = {1.0, 1.0, 1.0, 5.0, 5.0};

  std::vector<titan_contact_kernel::ContactCandidate> scalar;
  std::vector<titan_contact_kernel::ContactCandidate> optimized;
  titan_contact_kernel::scan_mol1_contacts_scalar(
      points, box_matrix, 1.5 * 1.5, 0, points.size(), scalar);
  titan_contact_kernel::scan_mol1_contacts_optimized(
      points, box_matrix, 1.5 * 1.5, 0, points.size(), optimized);

  require_same_candidates(scalar, optimized);
  require(scalar.size() == 3, "mol1 upper-triangle scan deduplicates i/j pairs");
  for (const auto &candidate : scalar)
    require(candidate.row < candidate.col, "mol1 candidate is upper-triangle");
}

void run_native_distance_range_contract() {
  const BoxMatrix box_matrix(make_box());
  titan_contact_kernel::PointSoA rows;
  rows.x = {1.0, 5.0, 9.5};
  rows.y = {1.0, 5.0, 1.0};
  rows.z = {1.0, 5.0, 1.0};

  titan_contact_kernel::PointSoA cols;
  cols.x = {9.0, 5.5, 0.4, 8.0, 1.5};
  cols.y = {1.0, 5.0, 1.0, 1.0, 1.0};
  cols.z = {1.0, 5.0, 1.0, 1.0, 1.0};

  const std::vector<titan_contact_kernel::NativePairIndex> pairs{
      {0, 0}, {1, 1}, {2, 2}, {0, 4}, {2, 3}};
  std::vector<double> scalar(pairs.size(), 0.0);
  std::vector<double> optimized(pairs.size(), 0.0);

  titan_contact_kernel::native_distances_scalar(rows, cols, pairs, 0,
                                                pairs.size(), box_matrix,
                                                scalar);
  titan_contact_kernel::native_distances_optimized(rows, cols, pairs, 0,
                                                   pairs.size(), box_matrix,
                                                   optimized);
  for (std::size_t i = 0; i < pairs.size(); ++i)
    require_close(optimized[i], scalar[i],
                  "native optimized distance matches scalar");
  require_close(scalar[0], 2.0, "native distance uses MIC wrapping");
}

void run_optional_benchmark() {
  if (std::getenv("TITAN_CONTACT_KERNEL_BENCH") == nullptr)
    return;

  const BoxMatrix box_matrix(make_box());
  titan_contact_kernel::PointSoA rows;
  titan_contact_kernel::PointSoA cols;
  const std::size_t n = 1000;
  rows.x.reserve(n);
  rows.y.reserve(n);
  rows.z.reserve(n);
  cols.x.reserve(n);
  cols.y.reserve(n);
  cols.z.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    rows.x.push_back(static_cast<double>(i % 100) * 0.1);
    rows.y.push_back(static_cast<double>((i * 3) % 100) * 0.1);
    rows.z.push_back(static_cast<double>((i * 7) % 100) * 0.1);
    cols.x.push_back(static_cast<double>((i * 11) % 100) * 0.1);
    cols.y.push_back(static_cast<double>((i * 13) % 100) * 0.1);
    cols.z.push_back(static_cast<double>((i * 17) % 100) * 0.1);
  }

  std::vector<titan_contact_kernel::ContactCandidate> scalar;
  std::vector<titan_contact_kernel::ContactCandidate> optimized;
  const auto scalar_start = std::chrono::steady_clock::now();
  titan_contact_kernel::scan_interface_contacts_scalar(
      rows, cols, box_matrix, 1.5 * 1.5, 0, rows.size(), scalar);
  const auto scalar_end = std::chrono::steady_clock::now();
  const auto optimized_start = std::chrono::steady_clock::now();
  titan_contact_kernel::scan_interface_contacts_optimized(
      rows, cols, box_matrix, 1.5 * 1.5, 0, rows.size(), optimized);
  const auto optimized_end = std::chrono::steady_clock::now();
  require_same_candidates(scalar, optimized);

  const auto scalar_us = std::chrono::duration_cast<std::chrono::microseconds>(
                             scalar_end - scalar_start)
                             .count();
  const auto optimized_us =
      std::chrono::duration_cast<std::chrono::microseconds>(optimized_end -
                                                            optimized_start)
          .count();
  std::cout << "BENCH titan_contact_kernels scalar_us=" << scalar_us
            << " optimized_us=" << optimized_us
            << " contacts=" << optimized.size() << '\n';
}

} // namespace

int main() {
  run_interface_tail_and_mic_contract();
  run_mol1_upper_triangle_contract();
  run_native_distance_range_contract();
  run_optional_benchmark();
  std::cout << "PASS titan_contact_kernels_smoke\n";
  return 0;
}
