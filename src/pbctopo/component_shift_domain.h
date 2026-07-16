#ifndef TITAN_PBCTOPO_COMPONENT_SHIFT_DOMAIN_H
#define TITAN_PBCTOPO_COMPONENT_SHIFT_DOMAIN_H

#include "../box_type.h"
#include "types.h"

#include <vector>

namespace titan_pbctopo {

std::vector<PbctopoShift> build_component_image_shifts(
    box Box, int shell_radius = 1);

} // namespace titan_pbctopo

#endif
