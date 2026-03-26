// ============================================================================
// TitanAlignment — Kabsch RMSD alignment shared by TICA, DCCM, PCA
//
// Design rules:
//   - kabsch_rotation()    : pure function, no side-effects
//   - align_frame()        : batch (allocates internally, for PCA post-loop)
//   - align_frame_inplace(): streaming (ZERO heap alloc, for TICA hot loop)
// ============================================================================
#ifndef TITAN_ALIGNMENT_H
#define TITAN_ALIGNMENT_H

#include <Eigen/Dense>

namespace TitanAlignment {

// -----------------------------------------------------------------------
// Result from a full align_frame() call (batch path)
// -----------------------------------------------------------------------
struct AlignResult {
  Eigen::Matrix3d R;               // rotation matrix
  Eigen::RowVector3d ref_centroid; // COG of reference (pre-centered)
  Eigen::RowVector3d mob_centroid; // COG of mobile (before centering)
};

// -----------------------------------------------------------------------
// Core Kabsch — both inputs MUST already be mean-centered.
// Returns R such that:  mobile_centered * R  ≈  ref_centered
// Handles reflection (determinant check).
// -----------------------------------------------------------------------
Eigen::Matrix3d kabsch_rotation(const Eigen::MatrixXd &ref_centered,   // (N×3)
                                const Eigen::MatrixXd &mobile_centered // (N×3)
);

// -----------------------------------------------------------------------
// Batch align — centers mobile, computes R, applies in-place.
// mobile is overwritten with aligned+centered coordinates.
// For PCA (called outside hot loops). Does allocate internally.
// -----------------------------------------------------------------------
AlignResult
align_frame(const Eigen::MatrixXd &ref_centered, // (N×3) pre-centered reference
            Eigen::MatrixXd &mobile // (N×3) in/out: centered + rotated
);

// -----------------------------------------------------------------------
// Streaming align — ZERO heap allocation.
// Pre-allocate mobile (N×3) OUTSIDE the hot loop and reuse each call.
//   mobile : in/out — raw coords → aligned+centered (overwritten in-place)
//   R_out  : pre-allocated Matrix3d written with the rotation (reused)
// -----------------------------------------------------------------------
void align_frame_inplace(
    const Eigen::MatrixXd &ref_centered, // (N×3) fixed reference
    Eigen::MatrixXd &mobile,             // (N×3) in/out
    Eigen::Matrix3d &R_out               // (3×3) pre-allocated, written here
);

} // namespace TitanAlignment
#endif // TITAN_ALIGNMENT_H
