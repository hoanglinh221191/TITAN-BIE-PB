#include "titan_alignment.h"

#include <cassert>
#include <cmath>
#include <vector>

namespace {
void near(double actual, double expected, double eps = 1.0e-10) {
  assert(std::abs(actual - expected) <= eps);
}
} // namespace

int main() {
  std::vector<double> ref = {
      -1.0, 0.0, 0.0,
       1.0, 0.0, 0.0,
       0.0, 1.0, 0.0,
       0.0,-1.0, 0.0,
  };
  std::vector<double> mobile = {
       0.0, 1.0, 0.0,
       0.0,-1.0, 0.0,
      -1.0, 0.0, 0.0,
       1.0, 0.0, 0.0,
  };

  titan::algebra::Matrix3ConstView ref_view{ref.data(), 4, 3};
  titan::algebra::Matrix3View mobile_view{mobile.data(), 4, 3};
  titan::algebra::Mat3 R{};
  TitanAlignment::align_frame_inplace(ref_view, mobile_view, R);

  near(titan::algebra::determinant(R), 1.0, 1.0e-8);
  for (std::size_t i = 0; i < ref.size(); ++i)
    near(mobile[i], ref[i], 1.0e-8);

  return 0;
}
