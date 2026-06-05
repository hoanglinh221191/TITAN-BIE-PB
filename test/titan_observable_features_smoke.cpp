#include "titan_observable_features.h"

#include <Eigen/Dense>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

namespace fs = std::filesystem;

void require(bool condition, const std::string &message) {
  if (!condition) {
    std::cerr << "FAIL titan_observable_features_smoke: " << message << '\n';
    std::exit(1);
  }
}

std::string read_file(const fs::path &path) {
  std::ifstream in(path, std::ios::binary);
  require(in.good(), "could not open " + path.string());
  return std::string(std::istreambuf_iterator<char>(in),
                     std::istreambuf_iterator<char>());
}

} // namespace

int main() {
  TitanObservableFeatureRouter router;
  router.configure(true, true, "mixed");
  std::vector<FeatureColumn> columns;
  columns.push_back(
      {"dccm", "ca_pair_A_1_B_2", "A", false, "A:1:ALA-B:2:LEU", "mol1"});
  columns.push_back({"qnative", "Q", "fraction", false, "", "mol1"});
  columns.push_back({"dihedral", "A_2_LEU_phi_sin", "sin", true,
                     "A:2:LEU", "mol1"});
  router.begin_trajectory(columns, "mixed");
  require(router.active(), "router is active");
  require(router.export_enabled(), "export is enabled");
  require(router.consume_tica(), "TICA consumption is enabled");
  require(router.feature_dim() == 3, "feature dimension matches metadata");

  Eigen::VectorXd finite(3);
  finite << 1.0, 0.5, -0.25;
  require(router.accept_frame(0, 0.0, finite), "finite frame accepted");

  Eigen::VectorXd nonfinite(3);
  nonfinite << 1.0, std::numeric_limits<double>::quiet_NaN(), 0.0;
  require(!router.accept_frame(1, 1.0, nonfinite),
          "non-finite frame rejected");

  Eigen::VectorXd wrong_dim(2);
  wrong_dim << 1.0, 2.0;
  require(!router.accept_frame(2, 2.0, wrong_dim),
          "dimension mismatch rejected");

  const fs::path out_dir =
      fs::temp_directory_path() / "titan_observable_features_smoke";
  fs::create_directories(out_dir);
  const fs::path matrix = out_dir / "features_by_frame.csv";
  const fs::path metadata = out_dir / "features_metadata.csv";
  const fs::path summary = out_dir / "feature_export_summary.dat";
  std::fstream log(out_dir / "feature_export.log",
                   std::ios::out | std::ios::trunc);
  router.export_outputs(matrix.string(), metadata.string(), summary.string(),
                        log);

  const std::string matrix_text = read_file(matrix);
  const std::string metadata_text = read_file(metadata);
  const std::string summary_text = read_file(summary);
  require(matrix_text.find("frame,time_ps,dccm:ca_pair_A_1_B_2,qnative:Q,"
                           "dihedral:A_2_LEU_phi_sin") != std::string::npos,
          "feature matrix header is stable");
  require(metadata_text.find("column,source,name,unit,periodic,residue,role") !=
              std::string::npos,
          "metadata header is present");
  require(summary_text.find("frames_accepted: 1") != std::string::npos,
          "summary records accepted frame count");
  require(summary_text.find("frames_skipped_nonfinite: 1") != std::string::npos,
          "summary records non-finite skipped frame count");
  require(summary_text.find("frames_skipped_dimension: 1") != std::string::npos,
          "summary records dimension skipped frame count");
  require(summary_text.find("tica_consumed: true") != std::string::npos,
          "summary records TICA consumption status");

  std::cout << "PASS titan_observable_features_smoke\n";
  return 0;
}
