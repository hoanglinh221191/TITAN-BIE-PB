#include "titan_executor.h"
#include "titan_sparse_graph.h"
#include "titan_sparse_graph_module.h"

#include <cstdlib>
#include <iostream>
#include <tuple>
#include <vector>

namespace {

using titan_sparse_graph_module::AnalysisConfig;
using titan_sparse_graph_module::SourceGraph;
using titan_sparse_graph_module::SourceKind;

[[noreturn]] void fail(const char *message) {
  std::cerr << "[titan_sparse_graph_smoke] " << message << "\n";
  std::exit(1);
}

SourceGraph make_graph(SourceKind source, bool directed, int node_count) {
  SourceGraph graph;
  graph.source = source;
  graph.directed = directed;
  graph.nodes.resize(static_cast<std::size_t>(node_count));
  for (int i = 0; i < node_count; ++i) {
    auto &node = graph.nodes[static_cast<std::size_t>(i)];
    node.node_id = "node_" + std::to_string(i);
    node.label = node.node_id;
    node.type = directed ? "Microstate" : "Residue";
    node.resname = "GLY";
    node.resnr = i + 1;
    node.chain = "A";
    node.state = directed ? i : -1;
  }
  return graph;
}

std::size_t argmax(const std::vector<double> &values) {
  if (values.empty())
    return 0;
  std::size_t best = 0;
  for (std::size_t i = 1; i < values.size(); ++i) {
    if (values[i] > values[best])
      best = i;
  }
  return best;
}

} // namespace

int main() {
  {
    const auto graph = titan_sparse::WeightedCSRGraph<double>::from_directed_edges(
        4, std::vector<std::tuple<int, int, double>>{{0, 1, 1.0},
                                                     {0, 2, 2.0},
                                                     {2, 3, 3.0}});
    if (graph.node_count() != 4 || graph.values().size() != 3) {
      fail("directed CSR builder returned unexpected graph size");
    }
    if (graph.row_begin(0) != 0 || graph.row_end(0) != 2 ||
        graph.row_begin(1) != 2 || graph.row_end(1) != 2 ||
        graph.row_begin(2) != 2 || graph.row_end(2) != 3) {
      fail("directed CSR builder row offsets mismatch");
    }
  }

  TitanExecutor executor(2);

  {
    SourceGraph graph = make_graph(SourceKind::AllosteryTe, true, 5);
    graph.edges = {
        {0, 1, 4.0, 4.0},
        {1, 2, 3.0, 3.0},
        {3, 4, 2.0, 2.0},
    };
    AnalysisConfig cfg;
    cfg.filter_component = "largest";
    const auto result = titan_sparse_graph_module::analyze_graph(graph, cfg, executor);
    if (result.filtered_graph.nodes.size() != 3 ||
        result.filtered_graph.edges.size() != 2 || result.component_count != 1) {
      fail("largest weak-component filtering should keep the 3-node directed component");
    }
  }

  {
    SourceGraph graph = make_graph(SourceKind::Dccm, false, 4);
    graph.edges = {
        {0, 1, 1.0, 1.0},
        {1, 2, 1.0, 1.0},
        {2, 3, 1.0, 1.0},
    };
    AnalysisConfig cfg;
    const auto result = titan_sparse_graph_module::analyze_graph(graph, cfg, executor);
    if (result.edge_betweenness.size() != 3) {
      fail("chain graph should produce 3 edge betweenness values");
    }
    const std::size_t max_edge = argmax(result.edge_betweenness);
    if (max_edge != 1) {
      fail("middle edge of an undirected chain should have the highest edge betweenness");
    }
  }

  {
    SourceGraph graph = make_graph(SourceKind::Dccm, false, 6);
    graph.edges = {
        {0, 1, 1.0, 1.0}, {0, 2, 1.0, 1.0}, {1, 2, 1.0, 1.0},
        {3, 4, 1.0, 1.0}, {3, 5, 1.0, 1.0}, {4, 5, 1.0, 1.0},
        {2, 3, 0.01, 0.01},
    };
    AnalysisConfig cfg;
    cfg.community_method = "louvain";
    const auto result = titan_sparse_graph_module::analyze_graph(graph, cfg, executor);
    if (!result.community_available ||
        result.community_membership.size() != graph.nodes.size()) {
      fail("Louvain community detection should run on undirected DCCM graphs");
    }
    if (result.community_membership[0] != result.community_membership[1] ||
        result.community_membership[0] != result.community_membership[2] ||
        result.community_membership[3] != result.community_membership[4] ||
        result.community_membership[3] != result.community_membership[5] ||
        result.community_membership[0] == result.community_membership[3]) {
      fail("two-clique graph should separate into two communities");
    }
  }

  {
    SourceGraph graph = make_graph(SourceKind::Tpt, true, 4);
    graph.edges = {
        {0, 1, 5.0, 5.0},
        {1, 2, 4.0, 4.0},
        {2, 3, 5.0, 5.0},
        {0, 2, 0.2, 0.2},
    };
    AnalysisConfig cfg;
    const auto result = titan_sparse_graph_module::analyze_graph(graph, cfg, executor);
    if (result.edge_betweenness.size() != 4 ||
        result.node_betweenness.size() != 4) {
      fail("directed TPT-like graph should export node/edge bottleneck metrics");
    }
    const std::size_t top_edge = argmax(result.edge_betweenness);
    if (top_edge != 1) {
      fail("middle reactive edge should dominate directed bottleneck ranking");
    }
  }

  std::cout << "titan_sparse_graph_smoke PASS\n";
  return 0;
}
