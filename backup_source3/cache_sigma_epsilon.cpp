#include "cache_sigma_epsilon.h"
#include <cmath>

int cache_sigma_epsilon::set_size(size_t MOL_size, const arguments &Arguments,
                                  size_t mol_1_size, size_t mol_2_size) {
  cutoff_dis = Arguments.distance_energy_cutoff;
  force_field_name = Arguments.force_f;

  if ((Arguments.flag_cal_force == true && MOL_size >= 200000000) ||
      (MOL_size >= 400000000) || Arguments.flag_cache == 2) {
    return 1;
  }

  elec_ener.assign(MOL_size, 0.0);
  vdw_ener.assign(MOL_size, 0.0);

  if (Arguments.flag_cal_force == true) {
    elec_force_x.assign(MOL_size, 0.0);
    elec_force_y.assign(MOL_size, 0.0);
    elec_force_z.assign(MOL_size, 0.0);
    vdw_force_x.assign(MOL_size, 0.0);
    vdw_force_y.assign(MOL_size, 0.0);
    vdw_force_z.assign(MOL_size, 0.0);
  }

  return 0;
}
