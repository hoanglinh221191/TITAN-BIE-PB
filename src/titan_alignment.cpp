#include "titan_alignment.h"

#include <Eigen/Dense>
#include <Eigen/SVD>

namespace {

Eigen::MatrixXd to_eigen(titan::algebra::Matrix3ConstView view) {
  Eigen::MatrixXd out(view.rows, 3);
  for (std::size_t r = 0; r < view.rows; ++r)
    for (std::size_t c = 0; c < 3; ++c)
      out(static_cast<Eigen::Index>(r), static_cast<Eigen::Index>(c)) =
          view(r, c);
  return out;
}

titan::algebra::Mat3 to_mat3(const Eigen::Matrix3d &m) {
  titan::algebra::Mat3 out{};
  for (std::size_t r = 0; r < 3; ++r)
    for (std::size_t c = 0; c < 3; ++c)
      out(r, c) = m(static_cast<Eigen::Index>(r), static_cast<Eigen::Index>(c));
  return out;
}

Eigen::Matrix3d to_eigen(const titan::algebra::Mat3 &m) {
  Eigen::Matrix3d out;
  for (std::size_t r = 0; r < 3; ++r)
    for (std::size_t c = 0; c < 3; ++c)
      out(static_cast<Eigen::Index>(r), static_cast<Eigen::Index>(c)) = m(r, c);
  return out;
}

} // namespace

namespace TitanAlignment {

titan::algebra::Mat3
kabsch_rotation(titan::algebra::Matrix3ConstView ref_centered,
                titan::algebra::Matrix3ConstView mobile_centered) {
  const Eigen::MatrixXd ref = to_eigen(ref_centered);
  const Eigen::MatrixXd mobile = to_eigen(mobile_centered);
  Eigen::Matrix3d H = ref.transpose() * mobile;
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
  AlignResult result;

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
                                       mobile.row_stride});

  const Eigen::Matrix3d R = to_eigen(result.R);
  for (std::size_t r = 0; r < mobile.rows; ++r) {
    Eigen::RowVector3d row(mobile(r, 0), mobile(r, 1), mobile(r, 2));
    row *= R;
    mobile(r, 0) = row(0);
    mobile(r, 1) = row(1);
    mobile(r, 2) = row(2);
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
