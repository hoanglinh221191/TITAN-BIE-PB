#include "titan_alignment.h"

#include <Eigen/Dense>
#include <Eigen/SVD>

#include <cassert>

namespace {

titan::algebra::Mat3 to_mat3(const Eigen::Matrix3d &m) {
  titan::algebra::Mat3 out{};
  for (std::size_t r = 0; r < 3; ++r)
    for (std::size_t c = 0; c < 3; ++c)
      out(r, c) = m(static_cast<Eigen::Index>(r), static_cast<Eigen::Index>(c));
  return out;
}

} // namespace

namespace TitanAlignment {

titan::algebra::Mat3
kabsch_rotation(titan::algebra::Matrix3ConstView ref_centered,
                titan::algebra::Matrix3ConstView mobile_centered) {
  assert(ref_centered.rows == mobile_centered.rows);

  Eigen::Matrix3d H = Eigen::Matrix3d::Zero();
  for (std::size_t r = 0; r < ref_centered.rows; ++r) {
    for (std::size_t i = 0; i < 3; ++i)
      for (std::size_t j = 0; j < 3; ++j)
        H(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) +=
            ref_centered(r, i) * mobile_centered(r, j);
  }

  Eigen::JacobiSVD<Eigen::Matrix3d> svd(H, Eigen::ComputeFullU |
                                               Eigen::ComputeFullV);
  Eigen::Matrix3d U = svd.matrixU();
  Eigen::Matrix3d V = svd.matrixV();
  if ((V * U.transpose()).determinant() < 0.0)
    V.col(2) *= -1.0;
  return to_mat3(V * U.transpose());
}

AlignResult align_frame(titan::algebra::Matrix3ConstView ref_centered,
                        titan::algebra::Matrix3View mobile) {
  assert(ref_centered.rows == mobile.rows);

  AlignResult result{};

  for (std::size_t r = 0; r < mobile.rows; ++r) {
    result.mob_centroid.x += mobile(r, 0);
    result.mob_centroid.y += mobile(r, 1);
    result.mob_centroid.z += mobile(r, 2);
  }
  result.mob_centroid = result.mob_centroid / static_cast<double>(mobile.rows);

  for (std::size_t r = 0; r < mobile.rows; ++r) {
    mobile(r, 0) -= result.mob_centroid.x;
    mobile(r, 1) -= result.mob_centroid.y;
    mobile(r, 2) -= result.mob_centroid.z;
  }

  result.R = kabsch_rotation(
      ref_centered,
      titan::algebra::Matrix3ConstView{mobile.data, mobile.rows,
                                       mobile.row_stride,
                                       mobile.col_stride});

  for (std::size_t r = 0; r < mobile.rows; ++r) {
    const double x = mobile(r, 0);
    const double y = mobile(r, 1);
    const double z = mobile(r, 2);
    mobile(r, 0) = x * result.R(0, 0) + y * result.R(1, 0) +
                   z * result.R(2, 0);
    mobile(r, 1) = x * result.R(0, 1) + y * result.R(1, 1) +
                   z * result.R(2, 1);
    mobile(r, 2) = x * result.R(0, 2) + y * result.R(1, 2) +
                   z * result.R(2, 2);
  }

  return result;
}

void align_frame_inplace(titan::algebra::Matrix3ConstView ref_centered,
                         titan::algebra::Matrix3View mobile,
                         titan::algebra::Mat3 &R_out) {
  const AlignResult result = align_frame(ref_centered, mobile);
  R_out = result.R;
}

} // namespace TitanAlignment
