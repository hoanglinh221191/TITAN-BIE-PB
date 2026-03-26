#ifndef CACHE_SIGMA_EPSILON_H
#define CACHE_SIGMA_EPSILON_H

#include "arguments.h"
#include "c6c12.h"
#include "miscellaneous.h"
#include <string>
#include <vector>

// Force standard includes

struct cache_energy {
  double elec_ener{0.0};
  double vdw_ener{0.0};

  cache_energy() = default;
  cache_energy(double e, double v) : elec_ener(e), vdw_ener(v) {}

  void set_energy(double e, double v) {
    elec_ener = e;
    vdw_ener = v;
  }

  cache_energy &operator+=(const cache_energy &other) {
    elec_ener += other.elec_ener;
    vdw_ener += other.vdw_ener;
    return *this;
  }
};

#include "force.h"

struct cache_force {
  force forces[2];

  cache_force() = default;
  cache_force(double ex, double ey, double ez, double vx, double vy,
              double vz) {
    set_force(ex, ey, ez, vx, vy, vz);
  }
  void set_force(const force &f1, const force &f2) {
    forces[0] = f1;
    forces[1] = f2;
  }

  void set_force(double e_x, double e_y, double e_z, double v_x, double v_y,
                 double v_z) {
    forces[0].Fx = e_x;
    forces[0].Fy = e_y;
    forces[0].Fz = e_z;
    forces[1].Fx = v_x;
    forces[1].Fy = v_y;
    forces[1].Fz = v_z;
  }

  cache_force &operator+=(const cache_force &other) {
    forces[0] += other.forces[0];
    forces[1] += other.forces[1];
    return *this;
  }
};

// OPT: Structure of Arrays (SoA) layout to maximize memory bandwidth and
// auto-vectorization
struct cache_sigma_epsilon {
  // 1D Arrays for Energy
  std::vector<double> elec_ener;
  std::vector<double> vdw_ener;

  // 1D Arrays for Force
  std::vector<double> elec_force_x;
  std::vector<double> elec_force_y;
  std::vector<double> elec_force_z;
  std::vector<double> vdw_force_x;
  std::vector<double> vdw_force_y;
  std::vector<double> vdw_force_z;

  std::vector<double> charge;
  std::string force_field_name;
  double cutoff_dis{-1.1};

  cache_sigma_epsilon() = default;

  cache_sigma_epsilon(double dis_cutoff)
      : cutoff_dis{dis_cutoff}, force_field_name{"amber99sb"} {}

  cache_sigma_epsilon(size_t MOL_size) : force_field_name{"amber99sb"} {
    elec_ener.assign(MOL_size, 0.0);
    vdw_ener.assign(MOL_size, 0.0);
    elec_force_x.assign(MOL_size, 0.0);
    elec_force_y.assign(MOL_size, 0.0);
    elec_force_z.assign(MOL_size, 0.0);
    vdw_force_x.assign(MOL_size, 0.0);
    vdw_force_y.assign(MOL_size, 0.0);
    vdw_force_z.assign(MOL_size, 0.0);
  }

  cache_sigma_epsilon(size_t MOL_size, const arguments &Arguments)
      : force_field_name{"amber99sb"} {
    elec_ener.assign(MOL_size, 0.0);
    vdw_ener.assign(MOL_size, 0.0);
    if (Arguments.flag_cal_force == 1) {
      elec_force_x.assign(MOL_size, 0.0);
      elec_force_y.assign(MOL_size, 0.0);
      elec_force_z.assign(MOL_size, 0.0);
      vdw_force_x.assign(MOL_size, 0.0);
      vdw_force_y.assign(MOL_size, 0.0);
      vdw_force_z.assign(MOL_size, 0.0);
    }
  }

  int set_size(size_t MOL_size, const arguments &Arguments, size_t mol_1_size,
               size_t mol_2_size);

  void set_elec_vdw(size_t pos, const double *elec, const double *vdw,
                    size_t sz) {
#pragma omp simd
    for (size_t i = 0; i < sz; ++i) {
      elec_ener[pos + i] = elec[i];
      vdw_ener[pos + i] = vdw[i];
    }
  }

  void set_elec_vdw(size_t pos, double elec, double vdw) {
    elec_ener[pos] = elec;
    vdw_ener[pos] = vdw;
  }

  cache_energy get_energy(size_t pos) const {
    return {elec_ener[pos], vdw_ener[pos]};
  }

  cache_force get_force(size_t pos) const {
    cache_force f;
    f.set_force(elec_force_x[pos], elec_force_y[pos], elec_force_z[pos],
                vdw_force_x[pos], vdw_force_y[pos], vdw_force_z[pos]);
    return f;
  }

  void set_force(size_t pos, double ex, double ey, double ez, double vx,
                 double vy, double vz) {
    elec_force_x[pos] = ex;
    elec_force_y[pos] = ey;
    elec_force_z[pos] = ez;
    vdw_force_x[pos] = vx;
    vdw_force_y[pos] = vy;
    vdw_force_z[pos] = vz;
  }
};
#endif // !CACHE_SIGMA_EPSILON_H