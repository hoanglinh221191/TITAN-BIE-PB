#define TITAN_LLPS_TESTING 1

#include "box_type.h"
#include "titan_llps_cluster.h"

#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

static void fail(const char *message) {
  std::cerr << "FAIL llps_paper_fractal_integration_smoke: " << message
            << "\n";
  std::exit(1);
}

static LLPSNode make_node(int id, double x, double y, double z) {
  LLPSNode node;
  node.id = id;
  node.type = 1;
  node.chain_idx = id;
  node.x = x;
  node.y = y;
  node.z = z;
  node.cluster_id = 1;
  node.is_core = true;
  return node;
}

int main() {
  LLPSPaperFractalConfig config;
  config.enabled = true;
  config.min_scale_nm = 0.1;
  config.max_scale_fraction = 1.0;

  std::vector<LLPSNode> nodes;
  nodes.push_back(make_node(0, 1.0, 50.0, 50.0));
  nodes.push_back(make_node(1, 1.2, 50.4, 50.1));
  nodes.push_back(make_node(2, 0.8, 49.8, 50.3));
  nodes.push_back(make_node(3, 1.4, 50.1, 49.7));
  nodes.push_back(make_node(4, 99.0, 50.2, 49.9));
  nodes.push_back(make_node(5, 98.8, 49.7, 50.2));
  nodes.push_back(make_node(6, 99.2, 50.3, 50.4));
  nodes.push_back(make_node(7, 98.6, 49.9, 49.6));
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
  if (frame.dense_protein_nodes != 8 || frame.dense_rna_nodes != 0)
    fail("dense protein/RNA counts changed");
  if (frame.dense_radius_nm >= 0.5)
    fail("dense radius should use MIC around the periodic seam");

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

  std::cout << "PASS llps_paper_fractal_integration_smoke\n";
  return 0;
}
