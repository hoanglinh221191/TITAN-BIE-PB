#ifndef TITAN_PBCTOPO_VIRTUAL_SEED_H
#define TITAN_PBCTOPO_VIRTUAL_SEED_H

namespace titan_pbctopo {

struct PbctopoVirtualSeedCenter {
  double cx = 1.5;
  double cy = 1.5;
  double cz = 1.5;
};

double pbctopo_normalize_virtual_center(double center);
PbctopoVirtualSeedCenter
pbctopo_center_delta_mic(PbctopoVirtualSeedCenter latest,
                         PbctopoVirtualSeedCenter previous);
PbctopoVirtualSeedCenter
extrapolate_virtual_seed_center(PbctopoVirtualSeedCenter previous,
                                PbctopoVirtualSeedCenter latest);

} // namespace titan_pbctopo

#endif // TITAN_PBCTOPO_VIRTUAL_SEED_H
