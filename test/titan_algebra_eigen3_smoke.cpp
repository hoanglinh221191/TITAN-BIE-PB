#include "algebra/symmetric_eigen.h"

#include <array>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

namespace {

using titan::algebra::Mat3;
using titan::algebra::Vec3;

void require(bool ok, const std::string &message) {
  if (!ok) {
    std::cerr << "FAIL titan_algebra_eigen3_smoke: " << message << '\n';
    std::exit(1);
  }
}

void near(double actual, double expected, const std::string &message,
          double eps = 1.0e-10) {
  if (std::abs(actual - expected) > eps) {
    std::cerr << "FAIL titan_algebra_eigen3_smoke: " << message
              << " actual=" << actual << " expected=" << expected << '\n';
    std::exit(1);
  }
}

Mat3 symmetric(double xx, double xy, double xz, double yy, double yz,
               double zz) {
  Mat3 out{};
  out(0, 0) = xx;
  out(0, 1) = xy;
  out(0, 2) = xz;
  out(1, 0) = xy;
  out(1, 1) = yy;
  out(1, 2) = yz;
  out(2, 0) = xz;
  out(2, 1) = yz;
  out(2, 2) = zz;
  return out;
}

Vec3 column(const Mat3 &m, std::size_t idx) {
  return {m(0, idx), m(1, idx), m(2, idx)};
}

void check_eigen_system(const Mat3 &a,
                        const titan::algebra::SymmetricEigen3 &es,
                        double residual_eps) {
  require(es.values[0] <= es.values[1] + 1.0e-12 &&
              es.values[1] <= es.values[2] + 1.0e-12,
          "eigenvalues must be sorted ascending");
  for (std::size_t i = 0; i < 3; ++i) {
    const Vec3 v = column(es.vectors, i);
    near(titan::algebra::norm(v), 1.0, "eigenvector norm", 1.0e-9);
    const Vec3 av = titan::algebra::mul(a, v);
    const Vec3 lv = v * es.values[i];
    require(titan::algebra::norm(av - lv) <= residual_eps,
            "A*v residual too large");
    for (std::size_t j = i + 1; j < 3; ++j)
      near(titan::algebra::dot(v, column(es.vectors, j)), 0.0,
           "eigenvectors must be orthogonal", 1.0e-8);
  }
}

} // namespace

int main() {
  {
    const Mat3 a = symmetric(3.0, 0.0, 0.0, 1.0, 0.0, 2.0);
    const auto es = titan::algebra::symmetric_eigen3(a);
    near(es.values[0], 1.0, "diagonal lambda 0");
    near(es.values[1], 2.0, "diagonal lambda 1");
    near(es.values[2], 3.0, "diagonal lambda 2");
    check_eigen_system(a, es, 1.0e-12);
  }
  {
    const Mat3 a = symmetric(4.0, 1.2, -0.4, 2.5, 0.7, 1.3);
    const auto es = titan::algebra::symmetric_eigen3(a);
    check_eigen_system(a, es, 1.0e-10);
  }
  {
    const Mat3 a = symmetric(1.0, 1.0e-8, 0.0, 1.0 + 2.0e-8, 1.0e-8,
                             1.0 + 4.0e-8);
    const auto es = titan::algebra::symmetric_eigen3(a);
    check_eigen_system(a, es, 1.0e-12);
  }
  {
    const std::array<double, 16> a = {
        1.0, 0.0, 0.0, 0.0,
        0.0, 2.0, 0.0, 0.0,
        0.0, 0.0, 3.0, 0.0,
        0.0, 0.0, 0.0, 4.0,
    };
    const auto dominant = titan::algebra::dominant_eigenvector_symmetric4(a);
    near(std::abs(dominant[3]), 1.0, "dominant diagonal 4x4 component");
    near(dominant[0], 0.0, "dominant diagonal 4x4 x", 1.0e-12);
    near(dominant[1], 0.0, "dominant diagonal 4x4 y", 1.0e-12);
    near(dominant[2], 0.0, "dominant diagonal 4x4 z", 1.0e-12);
  }
  std::cout << "PASS titan_algebra_eigen3_smoke\n";
  return 0;
}
