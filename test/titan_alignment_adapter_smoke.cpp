#include <Eigen/Dense>

#include "titan_alignment.h"

#include <cassert>
#include <cmath>

namespace {

void near(double actual, double expected, double eps = 1.0e-10) {
  assert(std::abs(actual - expected) <= eps);
}

void fill_column_major_fixture(Eigen::MatrixXd &ref, Eigen::MatrixXd &mobile) {
  ref.resize(4, 3);
  mobile.resize(4, 3);
  ref << -1.0, 0.0, 0.0,
          1.0, 0.0, 0.0,
          0.0, 1.0, 0.0,
          0.0,-1.0, 0.0;
  mobile << 0.0, 1.0, 0.0,
            0.0,-1.0, 0.0,
           -1.0, 0.0, 0.0,
            1.0, 0.0, 0.0;
}

template <typename Matrix>
void check_aligned(const Matrix &actual, const Matrix &expected) {
  for (typename Matrix::Index r = 0; r < expected.rows(); ++r)
    for (typename Matrix::Index c = 0; c < expected.cols(); ++c)
      near(actual(r, c), expected(r, c), 1.0e-8);
}

} // namespace

int main() {
  static_assert(TitanAlignment::detail::StridedMatrixLike<Eigen::MatrixXd>);
  static_assert(!TitanAlignment::detail::StridedMatrixLike<int>);

  Eigen::MatrixXd ref;
  Eigen::MatrixXd mobile;
  fill_column_major_fixture(ref, mobile);

  Eigen::Matrix3d R_out;
  TitanAlignment::align_frame_inplace(ref, mobile, R_out);
  near(R_out.determinant(), 1.0, 1.0e-8);
  check_aligned(mobile, ref);

  Eigen::MatrixXd mobile_for_kabsch;
  fill_column_major_fixture(ref, mobile_for_kabsch);
  Eigen::Matrix3d kabsch_R =
      TitanAlignment::kabsch_rotation(ref, mobile_for_kabsch);
  near(kabsch_R.determinant(), 1.0, 1.0e-8);

  Eigen::MatrixXd mobile_for_result;
  fill_column_major_fixture(ref, mobile_for_result);
  const auto result2 = TitanAlignment::align_frame(ref, mobile_for_result);
  Eigen::Matrix3d result_R = result2.R;
  Eigen::RowVector3d mob_centroid = result2.mob_centroid;
  near(result_R.determinant(), 1.0, 1.0e-8);
  near(mob_centroid(0), 0.0, 1.0e-8);
  near(mob_centroid(1), 0.0, 1.0e-8);
  near(mob_centroid(2), 0.0, 1.0e-8);
  check_aligned(mobile_for_result, ref);

  using RowMajorMatrix = Eigen::Matrix<double, Eigen::Dynamic, 3, Eigen::RowMajor>;
  RowMajorMatrix ref_rm(4, 3);
  RowMajorMatrix mobile_rm(4, 3);
  ref_rm << -1.0, 0.0, 0.0,
             1.0, 0.0, 0.0,
             0.0, 1.0, 0.0,
             0.0,-1.0, 0.0;
  mobile_rm << 0.0, 1.0, 0.0,
               0.0,-1.0, 0.0,
              -1.0, 0.0, 0.0,
               1.0, 0.0, 0.0;

  const TitanAlignment::AlignResult result =
      TitanAlignment::align_frame(ref_rm, mobile_rm);
  near(titan::algebra::determinant(result.R), 1.0, 1.0e-8);
  check_aligned(mobile_rm, ref_rm);

  return 0;
}
