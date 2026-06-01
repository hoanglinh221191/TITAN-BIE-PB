#include "titan_dihedral_kernels.h"

#include "BoxMatrix.h"
#include "titan_executor.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>

namespace {

void require(bool condition, const std::string &message) {
  if (!condition) {
    std::cerr << "FAIL titan_dihedral_kernels_smoke: " << message << '\n';
    std::exit(1);
  }
}

void require_close(double actual, double expected, const std::string &message,
                   double tolerance = 1.0e-9) {
  if (std::isnan(expected)) {
    require(std::isnan(actual), message);
    return;
  }
  if (std::fabs(actual - expected) > tolerance) {
    std::cerr << "FAIL titan_dihedral_kernels_smoke: " << message
              << " actual=" << actual << " expected=" << expected << '\n';
    std::exit(1);
  }
}

box make_box() {
  box frame_box;
  frame_box.set_a(10.0);
  frame_box.set_b(10.0);
  frame_box.set_c(10.0);
  frame_box.set_alpha(90.0);
  frame_box.set_beta(90.0);
  frame_box.set_gamma(90.0);
  return frame_box;
}

MOL make_frame() {
  MOL frame(make_box(), 0.0, "AMBER99SB");
  frame.pos.x = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 4.9};
  frame.pos.y = {0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 0.0};
  frame.pos.z = {0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0};
  return frame;
}

void run_canonical_contract() {
  using titan_dihedral_kernel::compute_dihedral_deg_from_points;
  require_close(compute_dihedral_deg_from_points(
                    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                    1.0, 1.0),
                90.0, "canonical +90");
  require_close(compute_dihedral_deg_from_points(
                    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                    1.0, -1.0),
                -90.0, "canonical -90");
  require_close(compute_dihedral_deg_from_points(
                    1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0,
                    1.0, 0.0),
                180.0, "canonical trans");
  require(std::isnan(compute_dihedral_deg_from_points(
              0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0,
              1.0)),
          "degenerate torsion returns nan");
}

void run_scalar_optimized_match() {
  titan_dihedral_kernel::TorsionPlanSoA plan;
  plan.push(0, 1, 2, 3, 0, titan_dihedral_kernel::TorsionKind::Phi);
  plan.push(0, 1, 2, 6, 0, titan_dihedral_kernel::TorsionKind::Psi);
  plan.push(4, 1, 2, 3, 1, titan_dihedral_kernel::TorsionKind::Omega);
  plan.push(0, 1, 2, 99, 2, titan_dihedral_kernel::TorsionKind::Chi1);
  plan.push(7, 1, 2, 3, 3, titan_dihedral_kernel::TorsionKind::Chi2);

  MOL frame = make_frame();
  titan_dihedral_kernel::TorsionFrameResult scalar;
  titan_dihedral_kernel::TorsionFrameResult optimized;
  BoxMatrix matrix(frame.Box);
  titan_dihedral_kernel::compute_torsions_scalar(plan, frame, true, &matrix,
                                                 0, plan.size(),
                                                 scalar.angle_deg);
  titan_dihedral_kernel::compute_torsions_optimized(plan, frame, true, &matrix,
                                                    0, plan.size(),
                                                    optimized.angle_deg);
  require(scalar.angle_deg.size() == optimized.angle_deg.size(),
          "scalar and optimized result sizes match");
  for (std::size_t i = 0; i < scalar.angle_deg.size(); ++i) {
    if (std::isnan(scalar.angle_deg[i]))
      require(std::isnan(optimized.angle_deg[i]), "nan lane matches");
    else
      require_close(optimized.angle_deg[i], scalar.angle_deg[i],
                    "optimized lane matches scalar");
  }
}

void run_executor_order_contract() {
  titan_dihedral_kernel::TorsionPlanSoA plan;
  for (int i = 0; i < 17; ++i) {
    plan.push(0, 1, 2, (i % 2) == 0 ? 3 : 6, 0,
              titan_dihedral_kernel::TorsionKind::Phi);
  }
  MOL frame = make_frame();
  BoxMatrix matrix(frame.Box);
  TitanExecutor one_worker(1);
  TitanExecutor four_workers(4);
  titan_dihedral_kernel::TorsionFrameResult one;
  titan_dihedral_kernel::TorsionFrameResult four;
  titan_dihedral_kernel::compute_torsions_parallel(plan, frame, true, &matrix,
                                                   one_worker, one);
  titan_dihedral_kernel::compute_torsions_parallel(plan, frame, true, &matrix,
                                                   four_workers, four);
  require(one.angle_deg.size() == four.angle_deg.size(),
          "parallel result sizes match");
  for (std::size_t i = 0; i < one.angle_deg.size(); ++i) {
    require_close(four.angle_deg[i], one.angle_deg[i],
                  "parallel output order is stable");
  }
}

} // namespace

int main() {
  run_canonical_contract();
  run_scalar_optimized_match();
  run_executor_order_contract();
  std::cout << "PASS titan_dihedral_kernels_smoke\n";
  return 0;
}
