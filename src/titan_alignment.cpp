// ============================================================================
// TitanAlignment — Kabsch rotation implementation
// ============================================================================
#include "titan_alignment.h"
#include <Eigen/Dense>
#include <Eigen/SVD>

namespace TitanAlignment {

// -----------------------------------------------------------------------
// kabsch_rotation — core SVD-based rotation
// -----------------------------------------------------------------------
Eigen::Matrix3d kabsch_rotation(const Eigen::MatrixXd &ref_centered,
                                const Eigen::MatrixXd &mobile_centered) {
  // H = ref^T * mobile  (3×3 cross-covariance)
  Eigen::Matrix3d H = ref_centered.transpose() * mobile_centered;

  Eigen::JacobiSVD<Eigen::Matrix3d> svd(H, Eigen::ComputeFullU |
                                               Eigen::ComputeFullV);

  Eigen::Matrix3d U = svd.matrixU();
  Eigen::Matrix3d V = svd.matrixV();

  // Reflection check: ensure a proper rotation (det = +1)
  if ((V * U.transpose()).determinant() < 0.0)
    V.col(2) *= -1.0;

  // R: mobile * R ≈ ref
  return V * U.transpose();
}

// -----------------------------------------------------------------------
// align_frame — batch path (may allocate)
// -----------------------------------------------------------------------
AlignResult align_frame(const Eigen::MatrixXd &ref_centered,
                        Eigen::MatrixXd &mobile) {
  AlignResult result;

  // Center mobile
  result.mob_centroid = mobile.colwise().mean();
  mobile.rowwise() -= result.mob_centroid;

  result.R = kabsch_rotation(ref_centered, mobile);
  mobile *= result.R; // rotate in-place: (N×3)*(3×3) — no temp if eval-ed

  return result;
}

// -----------------------------------------------------------------------
// align_frame_inplace — streaming path, ZERO heap allocation
// -----------------------------------------------------------------------
void align_frame_inplace(const Eigen::MatrixXd &ref_centered,
                         Eigen::MatrixXd &mobile, Eigen::Matrix3d &R_out) {
  // Center in-place using Eigen expression (no temporary vector allocated)
  mobile.rowwise() -= mobile.colwise().mean().eval();

  R_out = kabsch_rotation(ref_centered, mobile);
  mobile *= R_out; // (N×3) × (3×3): Eigen avoids temp with aliasing rules
}

} // namespace TitanAlignment
