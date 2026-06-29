#include "algebra/mat3.h"
#include "algebra/matrix3_buffer.h"
#include "algebra/matrix_view.h"
#include "algebra/vec3.h"

#include <cassert>
#include <cmath>
#include <vector>

namespace {
void near(double actual, double expected, double eps = 1.0e-12) {
  assert(std::abs(actual - expected) <= eps);
}
} // namespace

int main() {
  using titan::algebra::Mat3;
  using titan::algebra::Matrix3Buffer;
  using titan::algebra::Matrix3ConstView;
  using titan::algebra::Matrix3View;
  using titan::algebra::Vec3;

  const Vec3 x{1.0, 0.0, 0.0};
  const Vec3 y{0.0, 1.0, 0.0};
  const Vec3 z = titan::algebra::cross(x, y);
  near(z.x, 0.0);
  near(z.y, 0.0);
  near(z.z, 1.0);
  near(titan::algebra::dot(z, x), 0.0);
  near(titan::algebra::norm(Vec3{3.0, 4.0, 0.0}), 5.0);

  const Mat3 identity = Mat3::identity();
  near(titan::algebra::determinant(identity), 1.0);
  const Vec3 moved = titan::algebra::mul(identity, Vec3{2.0, 3.0, 4.0});
  near(moved.x, 2.0);
  near(moved.y, 3.0);
  near(moved.z, 4.0);

  Mat3 non_identity{};
  non_identity(0, 0) = 1.0;
  non_identity(0, 1) = 2.0;
  non_identity(0, 2) = 3.0;
  non_identity(1, 0) = 0.0;
  non_identity(1, 1) = 1.0;
  non_identity(1, 2) = 4.0;
  non_identity(2, 0) = 5.0;
  non_identity(2, 1) = 6.0;
  non_identity(2, 2) = 0.0;
  near(titan::algebra::determinant(non_identity), 1.0);
  const Vec3 mixed = titan::algebra::mul(non_identity, Vec3{1.0, 2.0, 3.0});
  near(mixed.x, 14.0);
  near(mixed.y, 14.0);
  near(mixed.z, 17.0);
  const Mat3 transposed = titan::algebra::transpose(non_identity);
  near(transposed(0, 2), 5.0);
  near(transposed(2, 1), 4.0);

  std::vector<double> rows = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
  Matrix3View view{rows.data(), 2, 3};
  view(1, 2) = 7.0;
  Matrix3ConstView const_view{rows.data(), 2, 3};
  near(const_view(0, 1), 2.0);
  near(const_view(1, 2), 7.0);

  std::vector<double> columns = {1.0, 4.0, 2.0, 5.0, 3.0, 6.0};
  Matrix3ConstView column_major_view{columns.data(), 2, 1, 2};
  near(column_major_view(0, 0), 1.0);
  near(column_major_view(1, 0), 4.0);
  near(column_major_view(0, 2), 3.0);
  near(column_major_view(1, 2), 6.0);

  Matrix3Buffer buffer;
  buffer.resize(2);
  assert(buffer.rows() == 2);
  for (std::size_t i = 0; i < 6; ++i)
    near(buffer.data()[i], 0.0);
  buffer(0, 0) = 1.0;
  buffer(0, 1) = 2.0;
  buffer(0, 2) = 3.0;
  buffer(1, 0) = 4.0;
  buffer(1, 1) = 5.0;
  buffer(1, 2) = 6.0;
  Matrix3View buffer_view = buffer.view();
  buffer_view(1, 2) = 7.0;
  Matrix3ConstView buffer_const_view = buffer.const_view();
  near(buffer.data()[0], 1.0);
  near(buffer.data()[1], 2.0);
  near(buffer.data()[2], 3.0);
  near(buffer.data()[3], 4.0);
  near(buffer.data()[4], 5.0);
  near(buffer.data()[5], 7.0);
  near(buffer_const_view(1, 2), 7.0);
  Matrix3Buffer copied = buffer;
  copied(0, 0) = 8.0;
  near(buffer(0, 0), 1.0);
  near(copied(0, 0), 8.0);

  return 0;
}
