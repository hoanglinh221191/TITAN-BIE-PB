#define TITAN_LLPS_TESTING 1

#include "arguments.h"
#include "box_type.h"
#include "info.h"
#include "titan_llps_cluster.h"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

static void fail(const char *message) {
  std::cerr << "FAIL llps_paper_fractal_integration_smoke: " << message
            << "\n";
  std::exit(1);
}

std::string arguments::set_Buildnumber() { return "test"; }

std::string arguments::get_build_number() const { return "test"; }

static LLPSNode make_node(int id, double x, double y, double z, int type = 1) {
  LLPSNode node;
  node.id = id;
  node.type = type;
  node.chain_idx = id;
  node.x = x;
  node.y = y;
  node.z = z;
  node.cluster_id = 1;
  node.is_core = true;
  return node;
}

static bool file_exists(const std::string &path) {
  std::ifstream in(path);
  return in.good();
}

static std::string read_text(const std::string &path) {
  std::ifstream in(path);
  std::ostringstream buffer;
  buffer << in.rdbuf();
  return buffer.str();
}

int main() {
  LLPSPaperFractalConfig config;
  config.enabled = true;
  config.min_scale_nm = 0.1;
  config.max_scale_fraction = 1.0;

  std::vector<LLPSNode> nodes;
  nodes.push_back(make_node(0, 1.0, 50.0, 50.0));
  nodes.push_back(make_node(1, 1.2, 50.4, 50.1, 2));
  nodes.push_back(make_node(2, 0.8, 49.8, 50.3));
  nodes.push_back(make_node(3, 1.4, 50.1, 49.7, 2));
  nodes.push_back(make_node(4, 99.0, 50.2, 49.9));
  nodes.push_back(make_node(5, 98.8, 49.7, 50.2, 2));
  nodes.push_back(make_node(6, 99.2, 50.3, 50.4));
  nodes.push_back(make_node(7, 98.6, 49.9, 49.6, 2));
  std::vector<bool> dense(nodes.size(), true);
  const box cubic_box(100.0, 100.0, 100.0, 90.0, 90.0, 90.0);

  LLPSPaperFractalAccumulator accumulator;
  accumulator.configure(config);
  accumulator.process_frame(1.25, 7, nodes, dense, 1, cubic_box);
  if (accumulator.frames().size() != 1)
    fail("expected one PBC seam frame");
  const auto &frame = accumulator.frames()[0];
  if (!frame.valid)
    fail("PBC seam frame should be valid");
  if (frame.dense_nodes != 8)
    fail("expected eight dense nodes");
  if (frame.dense_protein_nodes != 4 || frame.dense_rna_nodes != 4)
    fail("dense protein/RNA counts changed");
  if (frame.dense_radius_nm >= 0.5)
    fail("dense radius should use MIC around the periodic seam");
  if (frame.mean_local_valency <= 0.0)
    fail("mixed dense frame should compute local valency");
  if (frame.protein_rna_neighbor_fraction <= 0.0 ||
      frame.protein_rna_neighbor_fraction > 1.0)
    fail("mixed dense frame should compute protein-RNA neighbor fraction");

  LLPSPaperFractalAccumulator disabled;
  config.enabled = false;
  disabled.configure(config);
  disabled.process_frame(1.25, 7, nodes, dense, 1, cubic_box);
  if (!disabled.frames().empty())
    fail("disabled config should not emit rows");

  LLPSPaperFractalAccumulator too_few;
  config.enabled = true;
  too_few.configure(config);
  std::vector<LLPSNode> few_nodes(nodes.begin(), nodes.begin() + 7);
  std::vector<bool> few_dense(few_nodes.size(), true);
  too_few.process_frame(2.0, 8, few_nodes, few_dense, 1, cubic_box);
  if (too_few.frames().size() != 1)
    fail("expected one invalid row for too few dense nodes");
  if (too_few.frames()[0].valid)
    fail("too few dense nodes should be invalid");
  if (too_few.frames()[0].invalid_reason != "dense_nodes_lt_8")
    fail("too few dense nodes should preserve invalid reason");

  accumulator.process_frame(2.0, 8, few_nodes, few_dense, 1, cubic_box);

  arguments args;
  info output_info(args);
  std::fstream log("llps_paper_fractal_integration_smoke.log",
                   std::ios::out | std::ios::trunc);
  accumulator.export_outputs(output_info, log, 0);
  log.close();

  const std::string timeseries =
      output_info.get_llps_paper_fractal_timeseries_filename();
  const std::string scales =
      output_info.get_llps_paper_fractal_scales_filename();
  const std::string summary =
      output_info.get_llps_paper_fractal_summary_filename();
  const std::string readme =
      output_info.get_llps_paper_fractal_readme_filename();
  if (!file_exists(timeseries))
    fail("missing fractal timeseries export");
  if (!file_exists(scales))
    fail("missing fractal scales export");
  if (!file_exists(summary))
    fail("missing fractal summary export");
  if (!file_exists(readme))
    fail("missing fractal readme export");
  const std::string timeseries_text = read_text(timeseries);
  if (timeseries_text.find("df_mass_radius") == std::string::npos ||
      timeseries_text.find("df_box_count") == std::string::npos)
    fail("timeseries header missing fractal dimensions");
  const std::string summary_text = read_text(summary);
  if (summary_text.find(
          "interpretation: morphology_descriptor_not_phase_boundary") ==
      std::string::npos)
    fail("summary missing morphology interpretation warning");

  std::cout << "PASS llps_paper_fractal_integration_smoke\n";
  return 0;
}
