#include "titan_tica_estimator.h"

#include <cstdlib>
#include <iostream>
#include <string>

namespace {

void require(bool condition, const std::string &message) {
  if (!condition) {
    std::cerr << "FAIL titan_tica_estimator_smoke: " << message << '\n';
    std::exit(1);
  }
}

} // namespace

int main() {
  TitanTicaFrameEstimateConfig remd48;
  remd48.begin_time_ps = 100000.0;
  remd48.end_time_ps = 827850.0;
  remd48.dt_ps = 200.0;
  remd48.begin_set = true;
  remd48.end_set = true;
  remd48.trajectory_count = 48;
  remd48.number_of_files = 48;
  remd48.fallback_estimated_frames = 316.0;

  const auto estimated = estimate_tica_selected_frames(remd48);
  require(estimated.frames == 174720,
          "48-replica selected-frame estimate should be 174720");
  require(estimated.source == "time-window",
          "48-replica estimate should report time-window source");

  TitanTicaFrameEstimateConfig fallback;
  fallback.dt_ps = 0.0;
  fallback.end_set = false;
  fallback.trajectory_count = 0;
  fallback.number_of_files = 0;
  fallback.fallback_estimated_frames = -12.0;

  const auto fallback_estimated = estimate_tica_selected_frames(fallback);
  require(fallback_estimated.frames == 1,
          "fallback estimate should never be less than one frame");
  require(fallback_estimated.source == "fallback",
          "fallback estimate should report fallback source");

  std::cout << "PASS titan_tica_estimator_smoke\n";
  return 0;
}
