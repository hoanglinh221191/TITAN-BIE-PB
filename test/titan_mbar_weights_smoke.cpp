#include "titan_mbar_weights.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

namespace {

void require(bool condition, const std::string &message) {
  if (!condition) {
    std::cerr << "FAIL titan_mbar_weights_smoke: " << message << '\n';
    std::exit(1);
  }
}

double sum_weights(const std::vector<double> &weights) {
  double sum = 0.0;
  for (double value : weights)
    sum += value;
  return sum;
}

TitanMBAR::FrameWeightRow weight_row(int replica, std::size_t local_frame,
                                     double time_ps, double weight,
                                     int64_t global_frame_index) {
  TitanMBAR::FrameWeightRow row;
  row.replica_id = replica;
  row.local_frame_index = local_frame;
  row.time_ps = time_ps;
  row.weight = weight;
  row.sampled_temperature_k = 300.0;
  row.sampled_state_index = replica;
  row.global_frame_index = global_frame_index;
  return row;
}

} // namespace

int main() {
  TitanMBAR::FrameMetadata tica_metadata;
  tica_metadata.replica_ids = {0, 0, 0, 1, 1, 1};
  tica_metadata.timestamps_ps = {100.0, 200.0, 300.0, 100.0, 200.0, 300.0};
  tica_metadata.global_indices = {0, 1, 2, 3, 4, 5};

  const std::vector<TitanMBAR::FrameWeightRow> dense_weights = {
      weight_row(0, 0, 100.0, 0.10, 0),
      weight_row(0, 1, 200.0, 0.20, 1),
      weight_row(0, 2, 300.0, 0.30, 2),
      weight_row(1, 0, 100.0, 0.40, 3),
      weight_row(1, 1, 200.0, 0.50, 4),
      weight_row(1, 2, 300.0, 0.60, 5),
  };

  const auto tica_remap =
      TitanMBAR::remap_frame_weights_to_metadata(dense_weights, tica_metadata);
  require(tica_remap.mapped_global_frames == 6,
          "dense TICA metadata should map every frame");
  require(tica_remap.global_frame_weights.size() == 6,
          "dense TICA remap should allocate through max global index");
  require(std::fabs(sum_weights(tica_remap.global_frame_weights) - 2.10) < 1e-12,
          "dense TICA remap should preserve total weight");

  TitanMBAR::FrameMetadata pca_metadata;
  pca_metadata.replica_ids = {0, 1};
  pca_metadata.timestamps_ps = {200.0, 300.0};
  pca_metadata.global_indices = {10, 11};

  const auto pca_remap =
      TitanMBAR::remap_frame_weights_to_metadata(dense_weights, pca_metadata);
  require(pca_remap.mapped_global_frames == 2,
          "sparse PCA metadata should map only matching replica/time rows");
  require(pca_remap.skipped_missing_metadata == 4,
          "sparse PCA metadata should report skipped missing rows");
  require(pca_remap.global_frame_weights.size() == 12,
          "sparse PCA remap should allocate through target max global index");
  require(std::fabs(pca_remap.global_frame_weights[10] - 0.20) < 1e-12,
          "PCA remap should place replica 0 time 200 ps at global 10");
  require(std::fabs(pca_remap.global_frame_weights[11] - 0.60) < 1e-12,
          "PCA remap should place replica 1 time 300 ps at global 11");
  require(sum_weights(pca_remap.global_frame_weights) > 0.0,
          "sparse PCA remap should have positive mapped weight sum");

  TitanMBAR::FrameMetadata invalid_metadata;
  invalid_metadata.replica_ids = {0, 0};
  invalid_metadata.timestamps_ps = {100.0, 200.0};
  invalid_metadata.global_indices = {-1, 3};

  const auto invalid_remap =
      TitanMBAR::remap_frame_weights_to_metadata(dense_weights, invalid_metadata);
  require(invalid_remap.mapped_global_frames == 1,
          "invalid target global index should be skipped");
  require(invalid_remap.skipped_invalid_global_index == 1,
          "invalid target global index should be counted");
  require(invalid_remap.skipped_missing_metadata == 4,
          "missing target metadata should be counted separately");

  std::cout << "PASS titan_mbar_weights_smoke\n";
  return 0;
}
