#ifndef TITAN_PBCTOPO_LATTICE_MATH_H
#define TITAN_PBCTOPO_LATTICE_MATH_H

#include "../box_type.h"

#include <array>

namespace titan_pbctopo {

struct PbctopoNearestImage {
  double sx = 0.0;
  double sy = 0.0;
  double sz = 0.0;
  double x = 0.0;
  double y = 0.0;
  double z = 0.0;
  double d2 = 0.0;
  int image_x = 0;
  int image_y = 0;
  int image_z = 0;
  bool valid = false;
};

class PbctopoLatticeMetric {
public:
  explicit PbctopoLatticeMetric(const box &Box);
  PbctopoLatticeMetric(const double matrix[3][3],
                       const double inverse[3][3]);

  std::array<double, 3> scaled_to_cartesian(double sx, double sy,
                                             double sz) const;
  std::array<double, 3> cartesian_to_scaled(double x, double y,
                                             double z) const;
  PbctopoNearestImage nearest_image_scaled(double dsx, double dsy,
                                            double dsz) const;
  PbctopoNearestImage nearest_image_cartesian(double dx, double dy,
                                               double dz) const;
  double shortest_basis_length() const;
  bool basis_reduced() const noexcept { return basis_reduced_; }

private:
  void initialize_properties();

  double matrix_[3][3]{};
  double inverse_[3][3]{};
  double reduced_matrix_[3][3]{};
  double reduced_inverse_[3][3]{};
  int reduced_to_original_[3][3]{};
  int original_to_reduced_[3][3]{};
  bool orthogonal_basis_ = false;
  bool basis_reduced_ = false;
};

} // namespace titan_pbctopo

#endif
