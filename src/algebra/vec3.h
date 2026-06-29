#pragma once

#include <cassert>
#include <cmath>

namespace titan::algebra {

struct Vec3 {
  double x{0.0};
  double y{0.0};
  double z{0.0};

  constexpr Vec3() noexcept = default;
  constexpr Vec3(double x_in, double y_in, double z_in) noexcept
      : x(x_in), y(y_in), z(z_in) {}

  [[nodiscard]] constexpr double operator[](int idx) const noexcept {
    assert(idx >= 0 && idx < 3);
    return idx == 0 ? x : (idx == 1 ? y : z);
  }

  [[nodiscard]] constexpr double &operator[](int idx) noexcept {
    assert(idx >= 0 && idx < 3);
    return idx == 0 ? x : (idx == 1 ? y : z);
  }
};

[[nodiscard]] constexpr Vec3 operator+(Vec3 a, Vec3 b) noexcept {
  return {a.x + b.x, a.y + b.y, a.z + b.z};
}

[[nodiscard]] constexpr Vec3 operator-(Vec3 a, Vec3 b) noexcept {
  return {a.x - b.x, a.y - b.y, a.z - b.z};
}

[[nodiscard]] constexpr Vec3 operator*(Vec3 a, double s) noexcept {
  return {a.x * s, a.y * s, a.z * s};
}

[[nodiscard]] constexpr Vec3 operator*(double s, Vec3 a) noexcept {
  return a * s;
}

[[nodiscard]] constexpr Vec3 operator/(Vec3 a, double s) noexcept {
  return {a.x / s, a.y / s, a.z / s};
}

constexpr Vec3 &operator+=(Vec3 &a, Vec3 b) noexcept {
  a.x += b.x;
  a.y += b.y;
  a.z += b.z;
  return a;
}

[[nodiscard]] constexpr double dot(Vec3 a, Vec3 b) noexcept {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

[[nodiscard]] constexpr Vec3 cross(Vec3 a, Vec3 b) noexcept {
  return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z,
          a.x * b.y - a.y * b.x};
}

[[nodiscard]] constexpr double norm2(Vec3 a) noexcept { return dot(a, a); }

[[nodiscard]] inline double norm(Vec3 a) noexcept {
  return std::sqrt(norm2(a));
}

[[nodiscard]] inline Vec3 normalized_or_zero(Vec3 a,
                                             double eps = 1.0e-12) noexcept {
  const double n = norm(a);
  return n > eps ? (a / n) : Vec3{};
}

} // namespace titan::algebra
