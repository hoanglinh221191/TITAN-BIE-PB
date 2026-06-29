#pragma once

#include "algebra/vec3.h"

#include <array>
#include <cstddef>

namespace titan::algebra {

struct Mat3 {
  std::array<double, 9> v{};

  [[nodiscard]] constexpr double operator()(std::size_t r,
                                            std::size_t c) const noexcept {
    return v[r * 3 + c];
  }

  [[nodiscard]] constexpr double &operator()(std::size_t r,
                                             std::size_t c) noexcept {
    return v[r * 3 + c];
  }

  [[nodiscard]] static constexpr Mat3 identity() noexcept {
    Mat3 m{};
    m(0, 0) = 1.0;
    m(1, 1) = 1.0;
    m(2, 2) = 1.0;
    return m;
  }
};

[[nodiscard]] constexpr Mat3 transpose(const Mat3 &m) noexcept {
  Mat3 out{};
  for (std::size_t r = 0; r < 3; ++r)
    for (std::size_t c = 0; c < 3; ++c)
      out(r, c) = m(c, r);
  return out;
}

[[nodiscard]] constexpr double determinant(const Mat3 &m) noexcept {
  return m(0, 0) * (m(1, 1) * m(2, 2) - m(1, 2) * m(2, 1)) -
         m(0, 1) * (m(1, 0) * m(2, 2) - m(1, 2) * m(2, 0)) +
         m(0, 2) * (m(1, 0) * m(2, 1) - m(1, 1) * m(2, 0));
}

[[nodiscard]] constexpr Vec3 mul(const Mat3 &m, Vec3 x) noexcept {
  return {m(0, 0) * x.x + m(0, 1) * x.y + m(0, 2) * x.z,
          m(1, 0) * x.x + m(1, 1) * x.y + m(1, 2) * x.z,
          m(2, 0) * x.x + m(2, 1) * x.y + m(2, 2) * x.z};
}

} // namespace titan::algebra
