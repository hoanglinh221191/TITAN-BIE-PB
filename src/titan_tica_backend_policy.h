#ifndef TITAN_TICA_BACKEND_POLICY_H
#define TITAN_TICA_BACKEND_POLICY_H

#include "titan_kinetics.h"

#include <algorithm>
#include <cmath>

struct TitanTicaBackendPolicyConfig {
  TicaBackend requested_backend{TicaBackend::Auto};
  long double available_ram_bytes{0.0L};
  long double ram_budget_fraction{0.30L};
  long double dense_memory_estimate{0.0L};
  long double dual_memory_estimate{0.0L};
  long double estimated_frames{0.0L};
  int feature_dim{0};
  bool need_projection_cache{false};
  long double large_matrix_value_limit{250000000.0L};
};

struct TitanTicaBackendPolicyDecision {
  TicaBackend selected_backend{TicaBackend::Auto};
  long double ram_budget{0.0L};
  long double raw_feature_values{0.0L};
  long double large_matrix_value_limit{250000000.0L};
  bool large_matrix_limit_exceeded{false};
  bool selected_by_large_matrix_limit{false};
  bool hard_memory_cap_exceeded{false};
  bool explicit_exact_dual_overridden{false};
};

inline long double titan_tica_positive_or_zero(long double value) {
  return std::isfinite(value) && value > 0.0L ? value : 0.0L;
}

inline TitanTicaBackendPolicyDecision
select_tica_backend_for_policy(const TitanTicaBackendPolicyConfig &cfg) {
  TitanTicaBackendPolicyDecision decision;
  decision.selected_backend = cfg.requested_backend;
  decision.ram_budget =
      cfg.available_ram_bytes > 0.0L
          ? cfg.available_ram_bytes * std::max(0.0L, cfg.ram_budget_fraction)
          : 0.0L;
  decision.raw_feature_values =
      titan_tica_positive_or_zero(cfg.estimated_frames) *
      static_cast<long double>(std::max(0, cfg.feature_dim));
  decision.large_matrix_value_limit =
      titan_tica_positive_or_zero(cfg.large_matrix_value_limit);
  decision.large_matrix_limit_exceeded =
      decision.large_matrix_value_limit > 0.0L &&
      decision.raw_feature_values >= decision.large_matrix_value_limit;
  decision.hard_memory_cap_exceeded =
      decision.ram_budget > 0.0L &&
      titan_tica_positive_or_zero(cfg.dual_memory_estimate) >
          decision.ram_budget;

  if (cfg.requested_backend == TicaBackend::Auto) {
    if (decision.large_matrix_limit_exceeded) {
      decision.selected_backend = TicaBackend::Randomized;
      decision.selected_by_large_matrix_limit = true;
    } else if (cfg.available_ram_bytes > 0.0L &&
               cfg.dense_memory_estimate <= decision.ram_budget) {
      decision.selected_backend = TicaBackend::ExactDense;
    } else if (cfg.available_ram_bytes <= 0.0L ||
               cfg.dual_memory_estimate <= decision.ram_budget) {
      decision.selected_backend = TicaBackend::ExactDual;
    } else {
      decision.selected_backend = TicaBackend::Randomized;
    }
  } else if (cfg.requested_backend == TicaBackend::ExactDual &&
             (decision.hard_memory_cap_exceeded ||
              (cfg.available_ram_bytes <= 0.0L &&
               decision.large_matrix_limit_exceeded))) {
    decision.selected_backend = TicaBackend::Randomized;
    decision.explicit_exact_dual_overridden = true;
  }

  // Projection-cache workflows need a backend that can replay per-frame
  // coordinates. Never let this exact_dense -> exact_dual rewrite bypass the
  // same hard cap applied to an explicit exact_dual request.
  if (cfg.need_projection_cache &&
      decision.selected_backend == TicaBackend::ExactDense) {
    decision.selected_backend = decision.hard_memory_cap_exceeded
                                    ? TicaBackend::Randomized
                                    : TicaBackend::ExactDual;
  }

  return decision;
}

#endif // TITAN_TICA_BACKEND_POLICY_H
