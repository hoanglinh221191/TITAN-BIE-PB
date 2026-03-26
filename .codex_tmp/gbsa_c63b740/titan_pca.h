#ifndef TITAN_PCA_H
#define TITAN_PCA_H

#include "arguments.h"
#include "info.h"
#include "titan_ui.h"
#include <Eigen/Dense>
#include <fstream>
#include <string>
#include <vector>

namespace TitanPCA {

// Structure representing a single grid bin on the 2D Free Energy Surface
struct FES_Grid {
  double pc1;
  double pc2;
  int count;
  double free_energy;
  int sink_id; // -1 means unassigned, else indexes a local minimum
};

// Structure representing a macrostate cluster defined by a prominent density
// sink
struct ClusterNode {
  int id;
  double min_pc1;
  double min_pc2;
  double min_energy;
  int representative_frame_idx;
  double min_euclidean_distance;
  int population;
  double percentage;
};

// Main API isolated from impl_engine.cpp handling Covariance, FES Histograms &
// Peak Clustering
void execute_pca_fes_analysis(
    const arguments &args,
    const std::vector<std::vector<std::array<double, 3>>> &ca_traj,
    const std::vector<std::string> &ca_resname,
    const std::vector<int> &ca_resnrs, const std::vector<std::string> &ca_chain,
    const std::vector<int> &ca_is_lig, info &info_out, std::fstream &log);

} // namespace TitanPCA

#endif
