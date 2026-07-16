#include "component_shift_domain.h"

#include <algorithm>
#include <cmath>

namespace titan_pbctopo {
namespace {

constexpr int kLegacyImageCoeff[27][3] = {
    {1, 0, 0},   {-1, 0, 0},  {0, 1, 0},   {0, -1, 0},  {0, 0, 1},
    {0, 0, -1},  {1, 1, 0},   {1, -1, 0},  {-1, 1, 0},  {-1, -1, 0},
    {1, 0, 1},   {1, 0, -1},  {-1, 0, 1},  {-1, 0, -1}, {0, 1, 1},
    {0, -1, 1},  {0, 1, -1},  {0, -1, -1}, {1, 1, 1},   {-1, 1, 1},
    {1, -1, 1},  {1, 1, -1},  {-1, -1, 1}, {-1, 1, -1}, {1, -1, -1},
    {-1, -1, -1}, {0, 0, 0}};

} // namespace

std::vector<PbctopoShift> build_component_image_shifts(box Box,
                                                        int shell_radius) {
  double v1[3]{}, v2[3]{}, v3[3]{};
  Box.get_v1(v1);
  Box.get_v2(v2);
  Box.get_v3(v3);

  shell_radius = std::clamp(shell_radius, 1, 3);
  const std::size_t width = static_cast<std::size_t>(2 * shell_radius + 1);
  std::vector<PbctopoShift> shifts;
  shifts.reserve(width * width * width);
  auto append = [&](int ix, int iy, int iz) {
    const std::size_t index = shifts.size();
    shifts.push_back({index,
                      ix,
                      iy,
                      iz,
                      ix * v1[0] + iy * v2[0] + iz * v3[0],
                      ix * v1[1] + iy * v2[1] + iz * v3[1],
                      ix * v1[2] + iy * v2[2] + iz * v3[2]});
  };

  for (const auto &coeff : kLegacyImageCoeff)
    append(coeff[0], coeff[1], coeff[2]);
  for (int shell = 2; shell <= shell_radius; ++shell) {
    for (int ix = -shell; ix <= shell; ++ix) {
      for (int iy = -shell; iy <= shell; ++iy) {
        for (int iz = -shell; iz <= shell; ++iz) {
          if (std::max({std::abs(ix), std::abs(iy), std::abs(iz)}) == shell)
            append(ix, iy, iz);
        }
      }
    }
  }

  std::stable_sort(shifts.begin(), shifts.end(),
                   [](const PbctopoShift &lhs, const PbctopoShift &rhs) {
                     const double lhs_norm2 =
                         lhs.tx * lhs.tx + lhs.ty * lhs.ty + lhs.tz * lhs.tz;
                     const double rhs_norm2 =
                         rhs.tx * rhs.tx + rhs.ty * rhs.ty + rhs.tz * rhs.tz;
                     if (std::abs(lhs_norm2 - rhs_norm2) > 1.0e-12)
                       return lhs_norm2 < rhs_norm2;
                     return lhs.shift_index < rhs.shift_index;
                   });
  return shifts;
}

} // namespace titan_pbctopo
