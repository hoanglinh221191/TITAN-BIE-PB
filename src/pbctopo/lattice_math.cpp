#include "lattice_math.h"

#include "../BoxMatrix.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <tuple>

namespace titan_pbctopo {
namespace {

int nearest_integer(double value) {
  return static_cast<int>(std::floor(value + 0.5));
}

double row_norm2(const double matrix[3][3], int row) {
  return matrix[row][0] * matrix[row][0] +
         matrix[row][1] * matrix[row][1] +
         matrix[row][2] * matrix[row][2];
}

double row_dot(const double matrix[3][3], int lhs, int rhs) {
  return matrix[lhs][0] * matrix[rhs][0] +
         matrix[lhs][1] * matrix[rhs][1] +
         matrix[lhs][2] * matrix[rhs][2];
}

int determinant3(const int matrix[3][3]) {
  return matrix[0][0] *
             (matrix[1][1] * matrix[2][2] -
              matrix[1][2] * matrix[2][1]) -
         matrix[0][1] *
             (matrix[1][0] * matrix[2][2] -
              matrix[1][2] * matrix[2][0]) +
         matrix[0][2] *
             (matrix[1][0] * matrix[2][1] -
              matrix[1][1] * matrix[2][0]);
}

int cofactor3(const int matrix[3][3], int row, int col) {
  int values[4]{};
  int next = 0;
  for (int r = 0; r < 3; ++r) {
    if (r == row)
      continue;
    for (int c = 0; c < 3; ++c) {
      if (c != col)
        values[next++] = matrix[r][c];
    }
  }
  const int minor = values[0] * values[3] - values[1] * values[2];
  return ((row + col) & 1) == 0 ? minor : -minor;
}

} // namespace

PbctopoLatticeMetric::PbctopoLatticeMetric(const box &Box) {
  const BoxMatrix source(Box);
  for (int row = 0; row < 3; ++row) {
    for (int col = 0; col < 3; ++col) {
      matrix_[row][col] = source.m[row][col];
      inverse_[row][col] = source.inv[row][col];
    }
  }
  initialize_properties();
}

PbctopoLatticeMetric::PbctopoLatticeMetric(
    const double matrix[3][3], const double inverse[3][3]) {
  for (int row = 0; row < 3; ++row) {
    for (int col = 0; col < 3; ++col) {
      matrix_[row][col] = matrix[row][col];
      inverse_[row][col] = inverse[row][col];
    }
  }
  initialize_properties();
}

void PbctopoLatticeMetric::initialize_properties() {
  for (int row = 0; row < 3; ++row) {
    for (int col = 0; col < 3; ++col) {
      reduced_matrix_[row][col] = matrix_[row][col];
      reduced_to_original_[row][col] = row == col ? 1 : 0;
    }
  }

  for (int iteration = 0; iteration < 48; ++iteration) {
    bool changed = false;
    for (int lhs = 0; lhs < 2; ++lhs) {
      for (int rhs = lhs + 1; rhs < 3; ++rhs) {
        const double lhs_norm = row_norm2(reduced_matrix_, lhs);
        const double rhs_norm = row_norm2(reduced_matrix_, rhs);
        if (rhs_norm + 1.0e-13 * std::max(1.0, lhs_norm) >= lhs_norm)
          continue;
        for (int axis = 0; axis < 3; ++axis) {
          std::swap(reduced_matrix_[lhs][axis],
                    reduced_matrix_[rhs][axis]);
          std::swap(reduced_to_original_[lhs][axis],
                    reduced_to_original_[rhs][axis]);
        }
        changed = true;
      }
    }
    for (int row = 1; row < 3; ++row) {
      for (int base = row - 1; base >= 0; --base) {
        const double denominator = row_norm2(reduced_matrix_, base);
        if (!std::isfinite(denominator) || denominator <= 1.0e-24)
          continue;
        const long long multiplier = std::llround(
            row_dot(reduced_matrix_, row, base) / denominator);
        if (multiplier == 0 || std::abs(multiplier) > 1000000)
          continue;
        for (int axis = 0; axis < 3; ++axis) {
          reduced_matrix_[row][axis] -=
              static_cast<double>(multiplier) *
              reduced_matrix_[base][axis];
          reduced_to_original_[row][axis] -=
              static_cast<int>(multiplier) *
              reduced_to_original_[base][axis];
        }
        changed = true;
      }
    }
    if (!changed)
      break;
  }

  const int determinant = determinant3(reduced_to_original_);
  if (std::abs(determinant) != 1) {
    for (int row = 0; row < 3; ++row) {
      for (int col = 0; col < 3; ++col) {
        reduced_matrix_[row][col] = matrix_[row][col];
        reduced_to_original_[row][col] = row == col ? 1 : 0;
      }
    }
    basis_reduced_ = false;
  } else {
    basis_reduced_ = false;
    for (int row = 0; row < 3; ++row) {
      for (int col = 0; col < 3; ++col) {
        original_to_reduced_[row][col] =
            cofactor3(reduced_to_original_, col, row) / determinant;
        basis_reduced_ =
            basis_reduced_ ||
            reduced_to_original_[row][col] != (row == col ? 1 : 0);
      }
    }
  }
  if (std::abs(determinant) != 1) {
    for (int row = 0; row < 3; ++row)
      for (int col = 0; col < 3; ++col)
        original_to_reduced_[row][col] = row == col ? 1 : 0;
  }
  for (int row = 0; row < 3; ++row) {
    for (int col = 0; col < 3; ++col) {
      reduced_inverse_[row][col] = 0.0;
      for (int inner = 0; inner < 3; ++inner) {
        reduced_inverse_[row][col] +=
            inverse_[row][inner] *
            static_cast<double>(original_to_reduced_[inner][col]);
      }
    }
  }

  orthogonal_basis_ = true;
  for (int lhs = 0; lhs < 3; ++lhs) {
    const double lhs_norm2 = row_norm2(reduced_matrix_, lhs);
    for (int rhs = lhs + 1; rhs < 3; ++rhs) {
      const double rhs_norm2 = row_norm2(reduced_matrix_, rhs);
      const double dot = row_dot(reduced_matrix_, lhs, rhs);
      const double scale = std::sqrt(lhs_norm2 * rhs_norm2);
      if (!std::isfinite(scale) || scale <= 0.0 ||
          std::fabs(dot) > 1.0e-12 * scale) {
        orthogonal_basis_ = false;
        return;
      }
    }
  }
}

std::array<double, 3>
PbctopoLatticeMetric::scaled_to_cartesian(double sx, double sy,
                                           double sz) const {
  return {matrix_[0][0] * sx + matrix_[1][0] * sy + matrix_[2][0] * sz,
          matrix_[0][1] * sx + matrix_[1][1] * sy + matrix_[2][1] * sz,
          matrix_[0][2] * sx + matrix_[1][2] * sy + matrix_[2][2] * sz};
}

std::array<double, 3>
PbctopoLatticeMetric::cartesian_to_scaled(double x, double y, double z) const {
  return {inverse_[0][0] * x + inverse_[1][0] * y + inverse_[2][0] * z,
          inverse_[0][1] * x + inverse_[1][1] * y + inverse_[2][1] * z,
          inverse_[0][2] * x + inverse_[1][2] * y + inverse_[2][2] * z};
}

PbctopoNearestImage PbctopoLatticeMetric::nearest_image_scaled(
    double dsx, double dsy, double dsz) const {
  PbctopoNearestImage best;
  best.d2 = std::numeric_limits<double>::infinity();
  if (!std::isfinite(dsx) || !std::isfinite(dsy) || !std::isfinite(dsz))
    return best;

  const double source_scaled[3]{dsx, dsy, dsz};
  double reduced_scaled[3]{};
  for (int reduced_axis = 0; reduced_axis < 3; ++reduced_axis) {
    for (int source_axis = 0; source_axis < 3; ++source_axis) {
      reduced_scaled[reduced_axis] +=
          static_cast<double>(original_to_reduced_[source_axis][reduced_axis]) *
          source_scaled[source_axis];
    }
  }

  const int seed_x = nearest_integer(reduced_scaled[0]);
  const int seed_y = nearest_integer(reduced_scaled[1]);
  const int seed_z = nearest_integer(reduced_scaled[2]);
  auto consider = [&](int nx, int ny, int nz) {
    const double reduced_residual[3]{
        reduced_scaled[0] - static_cast<double>(nx),
        reduced_scaled[1] - static_cast<double>(ny),
        reduced_scaled[2] - static_cast<double>(nz)};
    const std::array<double, 3> cart{
        reduced_matrix_[0][0] * reduced_residual[0] +
            reduced_matrix_[1][0] * reduced_residual[1] +
            reduced_matrix_[2][0] * reduced_residual[2],
        reduced_matrix_[0][1] * reduced_residual[0] +
            reduced_matrix_[1][1] * reduced_residual[1] +
            reduced_matrix_[2][1] * reduced_residual[2],
        reduced_matrix_[0][2] * reduced_residual[0] +
            reduced_matrix_[1][2] * reduced_residual[1] +
            reduced_matrix_[2][2] * reduced_residual[2]};
    const double d2 = cart[0] * cart[0] + cart[1] * cart[1] +
                      cart[2] * cart[2];
    const double tolerance =
        1.0e-13 * std::max({1.0, std::fabs(d2), std::fabs(best.d2)});
    const auto tie_key = std::tuple{cart[0], cart[1], cart[2]};
    const auto best_key = std::tuple{best.x, best.y, best.z};
    if (!best.valid || d2 + tolerance < best.d2 ||
        (std::fabs(d2 - best.d2) <= tolerance && tie_key < best_key)) {
      const int reduced_image[3]{-nx, -ny, -nz};
      int original_image[3]{};
      for (int original_axis = 0; original_axis < 3; ++original_axis) {
        for (int reduced_axis = 0; reduced_axis < 3; ++reduced_axis) {
          original_image[original_axis] +=
              reduced_to_original_[reduced_axis][original_axis] *
              reduced_image[reduced_axis];
        }
      }
      best = {dsx + original_image[0],
              dsy + original_image[1],
              dsz + original_image[2],
              cart[0],
              cart[1],
              cart[2],
              d2,
              original_image[0],
              original_image[1],
              original_image[2],
              true};
    }
  };

  // For a diagonal Gram matrix the closest-vector problem separates into
  // three independent nearest-integer problems.
  if (orthogonal_basis_) {
    consider(seed_x, seed_y, seed_z);
    return best;
  }

  for (int dx = -1; dx <= 1; ++dx)
    for (int dy = -1; dy <= 1; ++dy)
      for (int dz = -1; dz <= 1; ++dz)
        consider(seed_x + dx, seed_y + dy, seed_z + dz);

  const double radius = std::sqrt(std::max(0.0, best.d2));
  double bounds[3]{};
  for (int axis = 0; axis < 3; ++axis) {
    const double inverse_column_norm =
        std::sqrt(reduced_inverse_[0][axis] * reduced_inverse_[0][axis] +
                  reduced_inverse_[1][axis] * reduced_inverse_[1][axis] +
                  reduced_inverse_[2][axis] * reduced_inverse_[2][axis]);
    bounds[axis] = radius * inverse_column_norm + 1.0e-12;
  }
  const int min_x =
      static_cast<int>(std::ceil(reduced_scaled[0] - bounds[0]));
  const int max_x =
      static_cast<int>(std::floor(reduced_scaled[0] + bounds[0]));
  const int min_y =
      static_cast<int>(std::ceil(reduced_scaled[1] - bounds[1]));
  const int max_y =
      static_cast<int>(std::floor(reduced_scaled[1] + bounds[1]));
  const int min_z =
      static_cast<int>(std::ceil(reduced_scaled[2] - bounds[2]));
  const int max_z =
      static_cast<int>(std::floor(reduced_scaled[2] + bounds[2]));
  for (int nx = min_x; nx <= max_x; ++nx)
    for (int ny = min_y; ny <= max_y; ++ny)
      for (int nz = min_z; nz <= max_z; ++nz)
        consider(nx, ny, nz);
  return best;
}

PbctopoNearestImage PbctopoLatticeMetric::nearest_image_cartesian(
    double dx, double dy, double dz) const {
  const auto scaled = cartesian_to_scaled(dx, dy, dz);
  return nearest_image_scaled(scaled[0], scaled[1], scaled[2]);
}

double PbctopoLatticeMetric::shortest_basis_length() const {
  double shortest = std::numeric_limits<double>::infinity();
  for (int row = 0; row < 3; ++row) {
    const double length =
        std::sqrt(reduced_matrix_[row][0] * reduced_matrix_[row][0] +
                  reduced_matrix_[row][1] * reduced_matrix_[row][1] +
                  reduced_matrix_[row][2] * reduced_matrix_[row][2]);
    shortest = std::min(shortest, length);
  }
  return shortest;
}

} // namespace titan_pbctopo
