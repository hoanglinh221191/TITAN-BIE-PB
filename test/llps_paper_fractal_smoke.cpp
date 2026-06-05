#define TITAN_LLPS_TESTING 1

#include "titan_llps_cluster.h"

#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>

static void fail(const char *message) {
  std::cerr << "FAIL llps_paper_fractal_smoke: " << message << "\n";
  std::exit(1);
}

int main() {
  std::vector<std::array<double, 3>> points;
  for (int x = 0; x < 4; ++x) {
    for (int y = 0; y < 4; ++y) {
      for (int z = 0; z < 4; ++z) {
        points.push_back({0.25 * x, 0.25 * y, 0.25 * z});
      }
    }
  }
  const auto result = LLPSPaperFractalAccumulator::debug_fit_points(points);
  if (!result.valid)
    fail("compact cubic dense phase should be valid");
  if (result.scale_count_mass_radius < 3)
    fail("mass-radius scale count too small");
  if (result.scale_count_box_count < 3)
    fail("box-count scale count too small");
  if (result.df_mass_radius < 2.0 || result.df_mass_radius > 3.3)
    fail("mass-radius df outside compact 3D range");
  if (result.df_box_count < 2.0 || result.df_box_count > 3.3)
    fail("box-count df outside compact 3D range");
  if (result.fit_r2_mass_radius < 0.85)
    fail("mass-radius fit r2 too low");
  if (result.fit_r2_box_count < 0.85)
    fail("box-count fit r2 too low");
  std::cout << "PASS llps_paper_fractal_smoke\n";
  return 0;
}
