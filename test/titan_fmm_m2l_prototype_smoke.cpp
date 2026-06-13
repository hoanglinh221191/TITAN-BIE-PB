#include "tabi/tabi_fmm_m2l_prototype.h"

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

using titan::tabi::fmm_m2l::KernelKind;
using titan::tabi::fmm_m2l::JcpChannelPrototypeMetrics;
using titan::tabi::fmm_m2l::JcpDerivativeBridgeMetrics;
using titan::tabi::fmm_m2l::JcpK1BridgeMetrics;
using titan::tabi::fmm_m2l::PrototypeConfig;
using titan::tabi::fmm_m2l::PrototypeMetrics;
using titan::tabi::fmm_m2l::classify_scope_a_decision;
using titan::tabi::fmm_m2l::classify_scope_d_decision;
using titan::tabi::fmm_m2l::classify_scope_e_decision;
using titan::tabi::fmm_m2l::classify_scope_f_decision;
using titan::tabi::fmm_m2l::kernel_kind_name;
using titan::tabi::fmm_m2l::make_mini_protein_like_point_cloud;
using titan::tabi::fmm_m2l::make_tabi_panel_like_fixture;
using titan::tabi::fmm_m2l::point_cloud_from_panels;
using titan::tabi::fmm_m2l::run_jcp_channel_probe;
using titan::tabi::fmm_m2l::run_jcp_derivative_bridge_probe;
using titan::tabi::fmm_m2l::run_jcp_k1_bridge_probe;
using titan::tabi::fmm_m2l::run_point_cloud_probe;
using titan::tabi::fmm_m2l::run_synthetic_probe;

void fail(const std::string &message) {
  std::cerr << "FAIL titan_fmm_m2l_prototype_smoke: " << message << '\n';
  std::exit(1);
}

void expect(bool condition, const std::string &message) {
  if (!condition) {
    fail(message);
  }
}

void expect_close(double lhs, double rhs, double tol, const std::string &message) {
  if (std::abs(lhs - rhs) > tol) {
    fail(message + " lhs=" + std::to_string(lhs) + " rhs=" + std::to_string(rhs));
  }
}

KernelKind parse_kernel(const std::string &value) {
  if (value == "laplace") {
    return KernelKind::Laplace;
  }
  if (value == "yukawa") {
    return KernelKind::Yukawa;
  }
  throw std::invalid_argument("bad kernel: " + value);
}

void print_csv_row(std::size_t panel_count, const PrototypeConfig &config,
                   const PrototypeMetrics &metrics) {
  const double speedup =
      metrics.apply_seconds > 0.0 ? metrics.direct_seconds / metrics.apply_seconds : 0.0;
  std::cout << panel_count << ',' << kernel_kind_name(config.kernel) << ','
            << config.expansion_order << ',' << config.separation_theta << ','
            << (metrics.finite_outputs ? "ok" : "nonfinite") << ','
            << std::setprecision(12) << metrics.setup_seconds << ','
            << metrics.apply_seconds << ',' << metrics.direct_seconds << ','
            << metrics.near_pairs << ',' << metrics.far_pairs << ','
            << metrics.m2l_operators << ',' << metrics.rms_abs_error << ','
            << metrics.max_abs_error << ',' << metrics.rel_rms_error << ','
            << speedup << ',' << classify_scope_a_decision(metrics, config) << '\n';
}

int run_probe_cli(int argc, char **argv) {
  PrototypeConfig config;
  std::size_t panels = 256;
  unsigned seed = 7;
  bool header = false;
  std::string fixture = "synthetic";
  for (int i = 2; i < argc; ++i) {
    const std::string arg = argv[i];
    auto read_value = [&](const std::string &prefix) -> std::string {
      if (arg.rfind(prefix, 0) != 0) {
        return {};
      }
      return arg.substr(prefix.size());
    };
    if (arg == "--csv-header") {
      header = true;
    } else if (const std::string v = read_value("--panels="); !v.empty()) {
      panels = static_cast<std::size_t>(std::stoull(v));
    } else if (const std::string v = read_value("--kernel="); !v.empty()) {
      config.kernel = parse_kernel(v);
    } else if (const std::string v = read_value("--fixture="); !v.empty()) {
      fixture = v;
    } else if (const std::string v = read_value("--kappa="); !v.empty()) {
      config.kappa = std::stod(v);
    } else if (const std::string v = read_value("--leaf="); !v.empty()) {
      config.max_leaf_size = std::stoi(v);
    } else if (const std::string v = read_value("--order="); !v.empty()) {
      config.expansion_order = std::stoi(v);
    } else if (const std::string v = read_value("--theta="); !v.empty()) {
      config.separation_theta = std::stod(v);
    } else if (const std::string v = read_value("--threads="); !v.empty()) {
      config.threads = std::stoi(v);
    } else if (const std::string v = read_value("--seed="); !v.empty()) {
      seed = static_cast<unsigned>(std::stoul(v));
    } else {
      throw std::invalid_argument("unknown argument: " + arg);
    }
  }
  if (header) {
    std::cout << "panel_count,kernel,order,theta,status,setup_seconds,apply_seconds,"
                 "direct_seconds,near_pairs,far_pairs,m2l_operators,rms_abs_error,"
                 "max_abs_error,rel_rms_error,speedup_vs_direct,decision_class\n";
  }
  PrototypeMetrics metrics;
  if (fixture == "synthetic") {
    metrics = run_synthetic_probe(config, panels, seed);
  } else if (fixture == "mini-protein-like") {
    metrics = run_point_cloud_probe(config,
                                    make_mini_protein_like_point_cloud(panels, seed));
  } else if (fixture == "tabi-panel-like") {
    metrics = run_point_cloud_probe(
        config, point_cloud_from_panels(make_tabi_panel_like_fixture(panels, seed)));
  } else {
    throw std::invalid_argument("unknown fixture: " + fixture);
  }
  print_csv_row(panels, config, metrics);
  return 0;
}

int run_jcp_channel_probe_cli(int argc, char **argv) {
  PrototypeConfig config;
  std::size_t panels = 1024;
  unsigned seed = 47;
  bool header = false;
  for (int i = 2; i < argc; ++i) {
    const std::string arg = argv[i];
    auto read_value = [&](const std::string &prefix) -> std::string {
      if (arg.rfind(prefix, 0) != 0) {
        return {};
      }
      return arg.substr(prefix.size());
    };
    if (arg == "--csv-header") {
      header = true;
    } else if (const std::string v = read_value("--panels="); !v.empty()) {
      panels = static_cast<std::size_t>(std::stoull(v));
    } else if (const std::string v = read_value("--kernel="); !v.empty()) {
      config.kernel = parse_kernel(v);
    } else if (const std::string v = read_value("--kappa="); !v.empty()) {
      config.kappa = std::stod(v);
    } else if (const std::string v = read_value("--leaf="); !v.empty()) {
      config.max_leaf_size = std::stoi(v);
    } else if (const std::string v = read_value("--order="); !v.empty()) {
      config.expansion_order = std::stoi(v);
    } else if (const std::string v = read_value("--theta="); !v.empty()) {
      config.separation_theta = std::stod(v);
    } else if (const std::string v = read_value("--threads="); !v.empty()) {
      config.threads = std::stoi(v);
    } else if (const std::string v = read_value("--seed="); !v.empty()) {
      seed = static_cast<unsigned>(std::stoul(v));
    } else {
      throw std::invalid_argument("unknown argument: " + arg);
    }
  }
  if (header) {
    std::cout << "panel_count,kernel,order,theta,status,setup_seconds,apply_seconds,"
                 "direct_seconds,channels,aggregate_rel_rms_error,"
                 "max_channel_rel_rms_error,max_abs_error,near_pairs,far_pairs,"
                 "speedup_vs_direct,"
                 "scope_d_gate\n";
  }
  const JcpChannelPrototypeMetrics metrics =
      run_jcp_channel_probe(config, make_tabi_panel_like_fixture(panels, seed));
  const double speedup =
      metrics.apply_seconds > 0.0 ? metrics.direct_seconds / metrics.apply_seconds : 0.0;
  std::cout << panels << ',' << kernel_kind_name(config.kernel) << ','
            << config.expansion_order << ',' << config.separation_theta << ','
            << (metrics.finite_outputs ? "ok" : "nonfinite") << ','
            << std::setprecision(12) << metrics.setup_seconds << ','
            << metrics.apply_seconds << ',' << metrics.direct_seconds << ','
            << metrics.channels << ',' << metrics.aggregate_rel_rms_error << ','
            << metrics.max_channel_rel_rms_error << ',' << metrics.max_abs_error << ','
            << metrics.near_pairs << ',' << metrics.far_pairs << ',' << speedup << ','
            << classify_scope_d_decision(metrics) << '\n';
  return 0;
}

int run_jcp_k1_bridge_probe_cli(int argc, char **argv) {
  PrototypeConfig config;
  std::size_t panels = 1024;
  unsigned seed = 61;
  bool header = false;
  for (int i = 2; i < argc; ++i) {
    const std::string arg = argv[i];
    auto read_value = [&](const std::string &prefix) -> std::string {
      if (arg.rfind(prefix, 0) != 0) {
        return {};
      }
      return arg.substr(prefix.size());
    };
    if (arg == "--csv-header") {
      header = true;
    } else if (const std::string v = read_value("--panels="); !v.empty()) {
      panels = static_cast<std::size_t>(std::stoull(v));
    } else if (const std::string v = read_value("--kappa="); !v.empty()) {
      config.kappa = std::stod(v);
    } else if (const std::string v = read_value("--leaf="); !v.empty()) {
      config.max_leaf_size = std::stoi(v);
    } else if (const std::string v = read_value("--order="); !v.empty()) {
      config.expansion_order = std::stoi(v);
    } else if (const std::string v = read_value("--theta="); !v.empty()) {
      config.separation_theta = std::stod(v);
    } else if (const std::string v = read_value("--threads="); !v.empty()) {
      config.threads = std::stoi(v);
    } else if (const std::string v = read_value("--seed="); !v.empty()) {
      seed = static_cast<unsigned>(std::stoul(v));
    } else {
      throw std::invalid_argument("unknown argument: " + arg);
    }
  }
  if (header) {
    std::cout << "panel_count,order,theta,status,setup_seconds,apply_seconds,"
                 "direct_seconds,rel_rms_error,rms_abs_error,max_abs_error,"
                 "near_pairs,far_pairs,speedup_vs_direct,scope_e_gate\n";
  }
  const JcpK1BridgeMetrics metrics =
      run_jcp_k1_bridge_probe(config, make_tabi_panel_like_fixture(panels, seed));
  const double speedup =
      metrics.apply_seconds > 0.0 ? metrics.direct_seconds / metrics.apply_seconds : 0.0;
  std::cout << panels << ',' << config.expansion_order << ','
            << config.separation_theta << ','
            << (metrics.finite_outputs ? "ok" : "nonfinite") << ','
            << std::setprecision(12) << metrics.setup_seconds << ','
            << metrics.apply_seconds << ',' << metrics.direct_seconds << ','
            << metrics.rel_rms_error << ',' << metrics.rms_abs_error << ','
            << metrics.max_abs_error << ',' << metrics.near_pairs << ','
            << metrics.far_pairs << ',' << speedup << ','
            << classify_scope_e_decision(metrics) << '\n';
  return 0;
}

int run_jcp_derivative_bridge_probe_cli(int argc, char **argv) {
  PrototypeConfig config;
  std::size_t panels = 1024;
  unsigned seed = 73;
  bool header = false;
  config.kappa = 0.2;
  for (int i = 2; i < argc; ++i) {
    const std::string arg = argv[i];
    auto read_value = [&](const std::string &prefix) -> std::string {
      if (arg.rfind(prefix, 0) != 0) {
        return {};
      }
      return arg.substr(prefix.size());
    };
    if (arg == "--csv-header") {
      header = true;
    } else if (const std::string v = read_value("--panels="); !v.empty()) {
      panels = static_cast<std::size_t>(std::stoull(v));
    } else if (const std::string v = read_value("--kappa="); !v.empty()) {
      config.kappa = std::stod(v);
    } else if (const std::string v = read_value("--leaf="); !v.empty()) {
      config.max_leaf_size = std::stoi(v);
    } else if (const std::string v = read_value("--order="); !v.empty()) {
      config.expansion_order = std::stoi(v);
    } else if (const std::string v = read_value("--theta="); !v.empty()) {
      config.separation_theta = std::stod(v);
    } else if (const std::string v = read_value("--threads="); !v.empty()) {
      config.threads = std::stoi(v);
    } else if (const std::string v = read_value("--seed="); !v.empty()) {
      seed = static_cast<unsigned>(std::stoul(v));
    } else {
      throw std::invalid_argument("unknown argument: " + arg);
    }
  }
  if (header) {
    std::cout << "panel_count,order,theta,status,setup_seconds,apply_seconds,"
                 "direct_seconds,channels,aggregate_rel_rms_error,"
                 "max_channel_rel_rms_error,k2_rel_rms_error,k3_rel_rms_error,"
                 "k4_rel_rms_error,max_abs_error,near_pairs,far_pairs,"
                 "speedup_vs_direct,scope_f_gate\n";
  }
  const JcpDerivativeBridgeMetrics metrics = run_jcp_derivative_bridge_probe(
      config, make_tabi_panel_like_fixture(panels, seed));
  const double speedup =
      metrics.apply_seconds > 0.0 ? metrics.direct_seconds / metrics.apply_seconds : 0.0;
  std::cout << panels << ',' << config.expansion_order << ','
            << config.separation_theta << ','
            << (metrics.finite_outputs ? "ok" : "nonfinite") << ','
            << std::setprecision(12) << metrics.setup_seconds << ','
            << metrics.apply_seconds << ',' << metrics.direct_seconds << ','
            << metrics.channels << ',' << metrics.aggregate_rel_rms_error << ','
            << metrics.max_channel_rel_rms_error << ',' << metrics.k2_rel_rms_error
            << ',' << metrics.k3_rel_rms_error << ',' << metrics.k4_rel_rms_error
            << ',' << metrics.max_abs_error << ',' << metrics.near_pairs << ','
            << metrics.far_pairs << ',' << speedup << ','
            << classify_scope_f_decision(metrics) << '\n';
  return 0;
}

void run_smoke_tests() {
  PrototypeConfig direct_cfg;
  direct_cfg.separation_theta = 0.0;
  direct_cfg.max_leaf_size = 16;
  direct_cfg.expansion_order = 4;
  direct_cfg.threads = 1;
  const PrototypeMetrics direct = run_synthetic_probe(direct_cfg, 96, 11);
  expect(direct.finite_outputs, "direct-only outputs must be finite");
  expect(direct.far_pairs == 0, "theta=0 must keep far pairs at zero");
  expect(direct.near_pairs == 96ull * 96ull, "theta=0 must cover all pairs as near");
  expect(direct.rel_rms_error <= 1.0e-14, "direct-only prototype must match reference");

  PrototypeConfig mixed_cfg;
  mixed_cfg.separation_theta = 0.35;
  mixed_cfg.max_leaf_size = 12;
  mixed_cfg.expansion_order = 2;
  mixed_cfg.threads = 2;
  const PrototypeMetrics mixed = run_synthetic_probe(mixed_cfg, 160, 17);
  expect(mixed.finite_outputs, "mixed near/far outputs must be finite");
  expect(mixed.near_pairs > 0, "mixed fixture must have near work");
  expect(mixed.far_pairs > 0, "mixed fixture must have far work");
  expect(mixed.m2l_operators > 0, "mixed fixture must create cached operator keys");
  expect(mixed.rel_rms_error < 0.35, "low-order M2L approximation is unexpectedly poor");

  PrototypeConfig yukawa_cfg = mixed_cfg;
  yukawa_cfg.kernel = KernelKind::Yukawa;
  yukawa_cfg.kappa = 0.2;
  const PrototypeMetrics yukawa = run_synthetic_probe(yukawa_cfg, 160, 17);
  expect(yukawa.finite_outputs, "Yukawa outputs must be finite");
  expect(yukawa.far_pairs > 0, "Yukawa fixture must have far work");
  expect(yukawa.rel_rms_error < 0.35, "low-order Yukawa approximation is unexpectedly poor");

  PrototypeConfig nt1 = mixed_cfg;
  nt1.threads = 1;
  PrototypeConfig nt4 = mixed_cfg;
  nt4.threads = 4;
  const PrototypeMetrics one = run_synthetic_probe(nt1, 192, 23);
  const PrototypeMetrics four = run_synthetic_probe(nt4, 192, 23);
  expect(one.near_pairs == four.near_pairs, "near-pair count must be deterministic");
  expect(one.far_pairs == four.far_pairs, "far-pair count must be deterministic");
  expect(one.m2l_operators == four.m2l_operators, "operator count must be deterministic");
  expect_close(one.rel_rms_error, four.rel_rms_error, 1.0e-13,
               "nt=1 and nt=4 relative error must match");

  PrototypeConfig mini_cfg;
  mini_cfg.separation_theta = 0.55;
  mini_cfg.max_leaf_size = 32;
  mini_cfg.expansion_order = 2;
  mini_cfg.threads = 4;
  const auto mini_points = make_mini_protein_like_point_cloud(512, 31);
  const PrototypeMetrics mini = run_point_cloud_probe(mini_cfg, mini_points);
  expect(mini.finite_outputs, "mini-protein-like outputs must be finite");
  expect(mini.near_pairs > 0, "mini-protein-like fixture must have near work");
  expect(mini.far_pairs > 0, "mini-protein-like fixture must have far work");
  expect(mini.rel_rms_error < 0.05,
         "mini-protein-like low-order approximation is unexpectedly poor");

  PrototypeConfig panel_cfg;
  panel_cfg.separation_theta = 0.5;
  panel_cfg.max_leaf_size = 32;
  panel_cfg.expansion_order = 2;
  panel_cfg.threads = 4;
  const auto panel_points =
      point_cloud_from_panels(make_tabi_panel_like_fixture(768, 41));
  const PrototypeMetrics panel = run_point_cloud_probe(panel_cfg, panel_points);
  expect(panel.finite_outputs, "TABI-panel-like outputs must be finite");
  expect(panel.near_pairs > 0, "TABI-panel-like fixture must have near work");
  expect(panel.far_pairs > 0, "TABI-panel-like fixture must have far work");
  expect(panel.rel_rms_error < 0.05,
         "TABI-panel-like low-order approximation is unexpectedly poor");

  const JcpChannelPrototypeMetrics jcp =
      [&]() {
        PrototypeConfig channel_cfg = panel_cfg;
        channel_cfg.separation_theta = 0.8;
        return run_jcp_channel_probe(channel_cfg, make_tabi_panel_like_fixture(768, 53));
      }();
  expect(jcp.finite_outputs, "JCP-channel-shaped outputs must be finite");
  expect(jcp.channels == 4, "JCP-channel-shaped probe must report four channels");
  expect(jcp.near_pairs > 0, "JCP-channel-shaped probe must report near work");
  expect(jcp.far_pairs > 0, "JCP-channel-shaped probe must report far work");
  expect(jcp.aggregate_rel_rms_error < 0.05,
         "JCP-channel-shaped aggregate approximation is unexpectedly poor");
  expect(jcp.max_channel_rel_rms_error < 0.1,
         "JCP-channel-shaped channel approximation is unexpectedly poor");

  PrototypeConfig k1_cfg = panel_cfg;
  k1_cfg.separation_theta = 0.7;
  k1_cfg.kappa = 0.2;
  const JcpK1BridgeMetrics k1 =
      run_jcp_k1_bridge_probe(k1_cfg, make_tabi_panel_like_fixture(768, 67));
  expect(k1.finite_outputs, "JCP K1 bridge outputs must be finite");
  expect(k1.far_pairs > 0, "JCP K1 bridge must report far work");
  expect(k1.rel_rms_error < 0.05, "JCP K1 bridge approximation is unexpectedly poor");

  PrototypeConfig deriv_cfg = panel_cfg;
  deriv_cfg.separation_theta = 0.7;
  deriv_cfg.kappa = 0.2;
  const JcpDerivativeBridgeMetrics deriv = run_jcp_derivative_bridge_probe(
      deriv_cfg, make_tabi_panel_like_fixture(384, 79));
  expect(deriv.finite_outputs, "JCP derivative bridge outputs must be finite");
  expect(deriv.channels == 3, "JCP derivative bridge must report three channels");
  expect(deriv.far_pairs > 0, "JCP derivative bridge must report far work");
  expect(deriv.aggregate_rel_rms_error < 0.2,
         "JCP derivative bridge aggregate approximation is unexpectedly poor");
}

} // namespace

int main(int argc, char **argv) {
  try {
    if (argc > 1 && std::string(argv[1]) == "--probe") {
      return run_probe_cli(argc, argv);
    }
    if (argc > 1 && std::string(argv[1]) == "--jcp-channel-probe") {
      return run_jcp_channel_probe_cli(argc, argv);
    }
    if (argc > 1 && std::string(argv[1]) == "--jcp-k1-bridge-probe") {
      return run_jcp_k1_bridge_probe_cli(argc, argv);
    }
    if (argc > 1 && std::string(argv[1]) == "--jcp-derivative-bridge-probe") {
      return run_jcp_derivative_bridge_probe_cli(argc, argv);
    }
    run_smoke_tests();
    std::cout << "PASS titan_fmm_m2l_prototype_smoke\n";
    return 0;
  } catch (const std::exception &ex) {
    fail(ex.what());
  }
}
