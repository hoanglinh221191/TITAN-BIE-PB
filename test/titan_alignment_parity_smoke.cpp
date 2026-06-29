#include <Eigen/Dense>

#include "titan_alignment.h"

#include <array>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

void require(bool ok, const std::string &message) {
  if (!ok) {
    std::cerr << "FAIL titan_alignment_parity_smoke: " << message << '\n';
    std::exit(1);
  }
}

Eigen::Matrix3d eigen_reference_kabsch(const Eigen::MatrixXd &ref_centered,
                                       const Eigen::MatrixXd &mobile_centered) {
  Eigen::Matrix3d h = Eigen::Matrix3d::Zero();
  for (Eigen::Index r = 0; r < ref_centered.rows(); ++r)
    for (Eigen::Index i = 0; i < 3; ++i)
      for (Eigen::Index j = 0; j < 3; ++j)
        h(i, j) += ref_centered(r, i) * mobile_centered(r, j);
  Eigen::JacobiSVD<Eigen::Matrix3d> svd(h, Eigen::ComputeFullU |
                                               Eigen::ComputeFullV);
  Eigen::Matrix3d u = svd.matrixU();
  Eigen::Matrix3d v = svd.matrixV();
  if ((v * u.transpose()).determinant() < 0.0)
    v.col(2) *= -1.0;
  return v * u.transpose();
}

std::vector<double> centered(const std::vector<std::array<double, 3>> &points) {
  std::array<double, 3> c{0.0, 0.0, 0.0};
  for (const auto &p : points) {
    c[0] += p[0];
    c[1] += p[1];
    c[2] += p[2];
  }
  for (double &x : c)
    x /= static_cast<double>(points.size());
  std::vector<double> out(points.size() * 3);
  for (std::size_t r = 0; r < points.size(); ++r)
    for (std::size_t col = 0; col < 3; ++col)
      out[r * 3 + col] = points[r][col] - c[col];
  return out;
}

std::vector<double> rotated_translated(
    const std::vector<std::array<double, 3>> &points,
    const titan::algebra::Mat3 &row_rotation,
    const std::array<double, 3> &translation) {
  std::vector<double> out(points.size() * 3);
  const titan::algebra::Mat3 inv = titan::algebra::transpose(row_rotation);
  for (std::size_t r = 0; r < points.size(); ++r) {
    const titan::algebra::Vec3 p{points[r][0], points[r][1], points[r][2]};
    const titan::algebra::Vec3 q = titan::algebra::mul_row(p, inv);
    out[r * 3 + 0] = q.x + translation[0];
    out[r * 3 + 1] = q.y + translation[1];
    out[r * 3 + 2] = q.z + translation[2];
  }
  return out;
}

Eigen::MatrixXd to_eigen(const std::vector<double> &row_major,
                         std::size_t rows) {
  Eigen::MatrixXd out(static_cast<Eigen::Index>(rows), 3);
  for (std::size_t r = 0; r < rows; ++r)
    for (std::size_t c = 0; c < 3; ++c)
      out(static_cast<Eigen::Index>(r), static_cast<Eigen::Index>(c)) =
          row_major[r * 3 + c];
  return out;
}

double rmsd(const std::vector<double> &a, const std::vector<double> &b) {
  double sum = 0.0;
  for (std::size_t i = 0; i < a.size(); ++i) {
    const double d = a[i] - b[i];
    sum += d * d;
  }
  return std::sqrt(sum / static_cast<double>(a.size() / 3));
}

void run_case(const std::string &name,
              const std::vector<std::array<double, 3>> &ref_points,
              const titan::algebra::Mat3 &row_rotation,
              double tolerance) {
  const std::size_t rows = ref_points.size();
  std::vector<double> ref = centered(ref_points);
  std::vector<double> mobile =
      rotated_translated(ref_points, row_rotation, {3.5, -2.0, 1.25});
  std::vector<double> mobile_ref = mobile;

  titan::algebra::Matrix3ConstView ref_view{ref.data(), rows, 3, 1};
  titan::algebra::Matrix3View mobile_view{mobile.data(), rows, 3, 1};
  titan::algebra::Mat3 r{};
  TitanAlignment::align_frame_inplace(ref_view, mobile_view, r);

  Eigen::MatrixXd ref_eigen = to_eigen(ref, rows);
  Eigen::MatrixXd mob_eigen = to_eigen(mobile_ref, rows);
  Eigen::RowVector3d centroid = mob_eigen.colwise().mean();
  mob_eigen.rowwise() -= centroid;
  const Eigen::Matrix3d r_ref = eigen_reference_kabsch(ref_eigen, mob_eigen);
  mob_eigen = mob_eigen * r_ref;

  std::vector<double> expected(rows * 3);
  for (std::size_t row = 0; row < rows; ++row)
    for (std::size_t col = 0; col < 3; ++col)
      expected[row * 3 + col] =
          mob_eigen(static_cast<Eigen::Index>(row),
                    static_cast<Eigen::Index>(col));

  require(std::abs(titan::algebra::determinant(r) - 1.0) <= 1.0e-10,
          name + ": determinant must be +1");
  require(rmsd(mobile, expected) <= tolerance,
          name + ": aligned coordinates differ from Eigen reference");
  require(rmsd(mobile, ref) <= 1.0e-8 ||
              rmsd(expected, ref) > 1.0e-8,
          name + ": production alignment unexpectedly worse than reference");
}

titan::algebra::Mat3 rz90() {
  titan::algebra::Mat3 r = titan::algebra::Mat3::identity();
  r(0, 0) = 0.0;
  r(0, 1) = -1.0;
  r(1, 0) = 1.0;
  r(1, 1) = 0.0;
  return r;
}

titan::algebra::Mat3 mixed_rotation() {
  titan::algebra::Mat3 r{};
  r(0, 0) = 0.36;
  r(0, 1) = -0.48;
  r(0, 2) = 0.80;
  r(1, 0) = 0.80;
  r(1, 1) = 0.60;
  r(1, 2) = 0.00;
  r(2, 0) = -0.48;
  r(2, 1) = 0.64;
  r(2, 2) = 0.60;
  return r;
}

} // namespace

int main() {
  run_case("identity",
           {{-1.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0},
            {0.0, -1.0, 0.0}},
           titan::algebra::Mat3::identity(), 1.0e-10);
  run_case("rz90",
           {{-1.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0},
            {0.0, -1.0, 0.0}},
           rz90(), 1.0e-10);
  run_case("full_rank",
           {{-1.2, 0.4, 2.1}, {0.7, -0.3, 1.5}, {2.4, 1.1, -0.9},
            {-0.8, -1.7, -1.2}, {1.6, 2.2, 0.3}, {-2.1, 0.9, -0.4}},
           mixed_rotation(), 1.0e-8);
  run_case("near_collinear",
           {{-2.0, 0.0, 0.0}, {-1.0, 1.0e-7, 0.0},
            {0.0, -1.0e-7, 2.0e-7}, {1.0, 2.0e-7, -1.0e-7},
            {2.0, -2.0e-7, 1.0e-7}},
           mixed_rotation(), 1.0e-7);
  std::cout << "PASS titan_alignment_parity_smoke\n";
  return 0;
}
