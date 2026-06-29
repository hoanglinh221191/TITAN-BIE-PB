#pragma once

#include <cstddef>

namespace titan::algebra {

struct Matrix3ConstView {
  const double *data{nullptr};
  std::size_t rows{0};
  std::ptrdiff_t row_stride{3};

  [[nodiscard]] constexpr double operator()(std::size_t r,
                                            std::size_t c) const noexcept {
    return data[r * row_stride + static_cast<std::ptrdiff_t>(c)];
  }
};

struct Matrix3View {
  double *data{nullptr};
  std::size_t rows{0};
  std::ptrdiff_t row_stride{3};

  [[nodiscard]] constexpr double operator()(std::size_t r,
                                            std::size_t c) const noexcept {
    return data[r * row_stride + static_cast<std::ptrdiff_t>(c)];
  }

  [[nodiscard]] constexpr double &operator()(std::size_t r,
                                             std::size_t c) noexcept {
    return data[r * row_stride + static_cast<std::ptrdiff_t>(c)];
  }
};

} // namespace titan::algebra
