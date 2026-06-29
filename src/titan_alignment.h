#ifndef TITAN_ALIGNMENT_H
#define TITAN_ALIGNMENT_H

#include "algebra/mat3.h"
#include "algebra/matrix_view.h"
#include "algebra/vec3.h"

namespace TitanAlignment {

struct AlignResult {
  titan::algebra::Mat3 R;
  titan::algebra::Vec3 ref_centroid;
  titan::algebra::Vec3 mob_centroid;
};

titan::algebra::Mat3
kabsch_rotation(titan::algebra::Matrix3ConstView ref_centered,
                titan::algebra::Matrix3ConstView mobile_centered);

AlignResult align_frame(titan::algebra::Matrix3ConstView ref_centered,
                        titan::algebra::Matrix3View mobile);

void align_frame_inplace(titan::algebra::Matrix3ConstView ref_centered,
                         titan::algebra::Matrix3View mobile,
                         titan::algebra::Mat3 &R_out);

} // namespace TitanAlignment

#endif // TITAN_ALIGNMENT_H
