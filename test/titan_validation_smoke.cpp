#include "titan_executor.h"
#include "titan_tpt.h"
#include "titan_validation_core.h"

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <vector>

namespace {

using TitanValidation::ChapmanKolmogorovSummary;
using TitanValidation::DhamRebuildContext;
using TitanValidation::ModelConfig;
using TitanValidation::StateSegments;

[[noreturn]] void fail(const char *message) {
  std::cerr << "[titan_validation_smoke] " << message << "\n";
  std::exit(1);
}

double next_unit(std::uint64_t &state) {
  state = state * 6364136223846793005ULL + 1442695040888963407ULL;
  const std::uint64_t mantissa = (state >> 11) & ((1ULL << 53) - 1ULL);
  return static_cast<double>(mantissa) / static_cast<double>(1ULL << 53);
}

StateSegments make_markov_segments(double stay_prob, int n_segments, int seg_len,
                                   std::uint64_t seed) {
  StateSegments segments;
  segments.reserve(static_cast<std::size_t>(n_segments));
  int state = 0;
  for (int seg = 0; seg < n_segments; ++seg) {
    std::vector<int> segment;
    segment.reserve(static_cast<std::size_t>(seg_len));
    for (int frame = 0; frame < seg_len; ++frame) {
      if (next_unit(seed) > stay_prob)
        state = 1 - state;
      segment.push_back(state);
    }
    segments.push_back(std::move(segment));
  }
  return segments;
}

StateSegments make_nonmarkov_segments(int n_segments, int seg_len) {
  StateSegments segments;
  segments.reserve(static_cast<std::size_t>(n_segments));
  for (int seg = 0; seg < n_segments; ++seg) {
    std::vector<int> segment;
    segment.reserve(static_cast<std::size_t>(seg_len));
    for (int frame = 0; frame < seg_len; ++frame) {
      segment.push_back((frame % 4 >= 2) ? 1 : 0); // 0,0,1,1 repeat
    }
    segments.push_back(std::move(segment));
  }
  return segments;
}

TitanMSM make_three_state_base_msm() {
  TitanMSM msm(3, 2, 0.1);
  std::uint64_t noise = 0x123456789abcdef0ULL;
  int state = 0;
  for (int seg = 0; seg < 6; ++seg) {
    for (int frame = 0; frame < 60; ++frame) {
      const double u = next_unit(noise);
      if (u > 0.92) {
        if (state == 0) {
          state = 1;
        } else if (state == 2) {
          state = 1;
        } else {
          state = next_unit(noise) > 0.5 ? 0 : 2;
        }
      }
      const double center = state == 0 ? -3.0 : (state == 1 ? 0.0 : 3.0);
      Eigen::VectorXd proj(2);
      proj[0] = center + 0.08 * (next_unit(noise) - 0.5);
      proj[1] = 0.5 * center + 0.08 * (next_unit(noise) - 0.5);
      msm.add_projected_frame(proj, 0, seg * 1000 + frame, 1.0);
      msm.record_frame_time(static_cast<double>(seg * 100 + frame), 0);
    }
    msm.add_trajectory_break();
  }
  return msm;
}

TitanTPT make_reference_tpt(const Eigen::MatrixXd &tpm, const Eigen::VectorXd &pi) {
  TitanTPT tpt(tpm, pi);
  tpt.compute({0}, {1}, 4, 0.0, nullptr, 1);
  return tpt;
}

} // namespace

int main() {
  TitanExecutor executor(2);

  {
    const auto markov = make_markov_segments(0.94, 8, 120, 0xdeadbeefULL);
    const auto nonmarkov = make_nonmarkov_segments(8, 120);
    const ChapmanKolmogorovSummary ck_markov =
        TitanValidation::compute_ck_for_segments(2, 1, markov, false, {2},
                                                 "markov");
    const ChapmanKolmogorovSummary ck_nonmarkov =
        TitanValidation::compute_ck_for_segments(2, 1, nonmarkov, false, {2},
                                                 "nonmarkov");
    if (!ck_markov.valid || ck_markov.metrics.empty())
      fail("Markov CK summary should be valid");
    if (!ck_nonmarkov.valid || ck_nonmarkov.metrics.empty())
      fail("Non-Markov CK summary should be valid");
    if (!(ck_markov.metrics.front().weighted_frobenius <
          ck_nonmarkov.metrics.front().weighted_frobenius)) {
      fail("CK Frobenius error should be smaller for the Markov chain");
    }
  }

  {
    TitanMSM base = make_three_state_base_msm();
    const std::vector<ModelConfig> candidates = {
        {2, 1, 2},
        {2, 1, 3},
    };
    auto sweep = TitanValidation::evaluate_vamp_candidates(
        base, 2, candidates, 3, 80, 42u, executor);
    if (!sweep.valid)
      fail("VAMP sweep should find a valid candidate");
    if (sweep.selected.n_clusters != 3)
      fail("Three-cluster candidate should outrank the underclustered model");

    auto sweep_repeat = TitanValidation::evaluate_vamp_candidates(
        base, 2, candidates, 3, 80, 42u, executor);
    if (sweep_repeat.selected.n_clusters != sweep.selected.n_clusters)
      fail("VAMP selection should be deterministic");
  }

  {
    const auto segments = make_markov_segments(0.93, 8, 120, 0xabcdefULL);
    TitanValidation::MsmBootstrapAggregate bootstrap;
    std::string error;
    if (!TitanValidation::run_msm_bootstrap(segments, 2, 1, 0.1, {2, 3}, 1, 8,
                                            42u, 2, executor, bootstrap,
                                            error)) {
      fail("MSM bootstrap should succeed on synthetic segments");
    }
    if (bootstrap.successful_samples <= 0)
      fail("MSM bootstrap should record successful samples");

    Eigen::MatrixXd tpm;
    Eigen::VectorXd pi;
    std::int64_t transitions = 0;
    if (!TitanValidation::build_markov_model_from_segments(
            2, 1, segments, false, tpm, pi, transitions, error)) {
      fail("Reference TPT model construction failed");
    }
    TitanTPT reference_tpt = make_reference_tpt(tpm, pi);
    const auto kinetic_diag =
        TitanValidation::compute_tpt_kinetic_diagnostics_from_segments(
            segments, 2, 1, 0.1, reference_tpt);
    if (kinetic_diag.observed_a_to_b_events <= 0 ||
        kinetic_diag.observed_b_to_a_events <= 0) {
      fail("Synthetic kinetic diagnostics should observe both directions");
    }

    TitanValidation::TptBootstrapAggregate tpt_bootstrap;
    if (!TitanValidation::run_tpt_bootstrap(
            segments, 2, 1, 0.1, false, {0}, {1}, 4, 0.0, 2, 8, 7u, executor,
            tpt_bootstrap, error)) {
      fail("MSM-backed TPT bootstrap should succeed");
    }
    if (tpt_bootstrap.successful_samples <= 0)
      fail("MSM-backed TPT bootstrap should record successful samples");
  }

  {
    DhamRebuildContext context;
    context.available = true;
    context.n_states = 2;
    context.lag_frames = 1;
    context.n_replicas = 2;
    context.temps_k = {300.0, 320.0};
    context.target_temp_k = 310.0;
    context.max_iter = 256;
    context.tolerance = 1e-6;
    context.ensemble = "nvt";
    context.per_replica_dtrajs[0] = make_markov_segments(0.90, 1, 120, 11u)[0];
    context.per_replica_dtrajs[1] = make_markov_segments(0.88, 1, 120, 29u)[0];
    for (int replica = 0; replica < 2; ++replica) {
      std::vector<double> times;
      std::vector<double> observable;
      times.reserve(120);
      observable.reserve(120);
      for (int frame = 0; frame < 120; ++frame) {
        times.push_back(static_cast<double>(frame));
        const int state = context.per_replica_dtrajs[replica][frame];
        observable.push_back(state == 0 ? -5.0 + 0.1 * replica
                                        : -2.0 + 0.1 * replica);
      }
      context.per_replica_timestamps_ps[replica] = std::move(times);
      context.per_replica_observables[replica] = std::move(observable);
    }

    Eigen::MatrixXd dham_tpm;
    Eigen::VectorXd dham_pi;
    std::int64_t transitions = 0;
    std::string error;
    if (!TitanValidation::build_dham_model_from_context(
            context, 1, dham_tpm, dham_pi, transitions, error)) {
      fail("DHAM rebuild context should produce a valid model");
    }
    if (transitions <= 0 || dham_tpm.rows() != 2 || dham_pi.size() != 2)
      fail("DHAM rebuild output dimensions mismatch");

    TitanValidation::TptBootstrapAggregate dham_bootstrap;
    if (!TitanValidation::run_dham_tpt_bootstrap(context, 4, 0.0, {0}, {1}, 4,
                                                 9u, 2, executor,
                                                 dham_bootstrap, error)) {
      fail("DHAM-backed TPT bootstrap should succeed");
    }
    if (dham_bootstrap.successful_samples <= 0)
      fail("DHAM-backed TPT bootstrap should record successful samples");
  }

  std::cout << "titan_validation_smoke PASS\n";
  return 0;
}

