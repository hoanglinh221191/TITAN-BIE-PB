#include "titan_tica_backend_policy.h"

#include <cmath>
#include <cstdint>
#include <iostream>
#include <string>

namespace {

void require_true(bool ok, const std::string &message) {
  if (!ok) {
    std::cerr << "FAIL: " << message << "\n";
    std::exit(1);
  }
}

long double gib(long double value) {
  return value * 1024.0L * 1024.0L * 1024.0L;
}

void large_pairwise_case_prefers_randomized_even_when_ram_is_available() {
  TitanTicaBackendPolicyConfig cfg;
  cfg.requested_backend = TicaBackend::Auto;
  cfg.available_ram_bytes = gib(256.0L);
  cfg.ram_budget_fraction = 0.30L;
  cfg.estimated_frames = 174720.0L;
  cfg.feature_dim = 4656;
  cfg.dense_memory_estimate =
      6.0L * cfg.feature_dim * cfg.feature_dim * sizeof(double);
  cfg.dual_memory_estimate =
      3.0L * cfg.estimated_frames * cfg.feature_dim * sizeof(double);
  cfg.need_projection_cache = true;

  const auto decision = select_tica_backend_for_policy(cfg);

  require_true(decision.selected_backend == TicaBackend::Randomized,
               "48-replica pairwise case should auto-select randomized");
  require_true(decision.large_matrix_limit_exceeded,
               "large feature matrix limit should be reported");
  require_true(decision.raw_feature_values ==
                   static_cast<long double>(174720) *
                       static_cast<long double>(4656),
               "raw feature value count should match N*D");
}

void small_case_keeps_existing_exact_dense_behavior() {
  TitanTicaBackendPolicyConfig cfg;
  cfg.requested_backend = TicaBackend::Auto;
  cfg.available_ram_bytes = gib(32.0L);
  cfg.ram_budget_fraction = 0.30L;
  cfg.estimated_frames = 1000.0L;
  cfg.feature_dim = 128;
  cfg.dense_memory_estimate =
      6.0L * cfg.feature_dim * cfg.feature_dim * sizeof(double);
  cfg.dual_memory_estimate =
      3.0L * cfg.estimated_frames * cfg.feature_dim * sizeof(double);
  cfg.need_projection_cache = false;

  const auto decision = select_tica_backend_for_policy(cfg);

  require_true(decision.selected_backend == TicaBackend::ExactDense,
               "small dense-fit case should keep exact_dense");
  require_true(!decision.large_matrix_limit_exceeded,
               "small case should not report large matrix limit");
}

void projection_cache_still_rewrites_small_exact_dense_to_exact_dual() {
  TitanTicaBackendPolicyConfig cfg;
  cfg.requested_backend = TicaBackend::Auto;
  cfg.available_ram_bytes = gib(32.0L);
  cfg.ram_budget_fraction = 0.30L;
  cfg.estimated_frames = 1000.0L;
  cfg.feature_dim = 128;
  cfg.dense_memory_estimate =
      6.0L * cfg.feature_dim * cfg.feature_dim * sizeof(double);
  cfg.dual_memory_estimate =
      3.0L * cfg.estimated_frames * cfg.feature_dim * sizeof(double);
  cfg.need_projection_cache = true;

  const auto decision = select_tica_backend_for_policy(cfg);

  require_true(decision.selected_backend == TicaBackend::ExactDual,
               "small projection-cache case should still use exact_dual");
}

void explicit_backend_requests_are_honored() {
  TitanTicaBackendPolicyConfig cfg;
  cfg.requested_backend = TicaBackend::ExactDual;
  cfg.available_ram_bytes = gib(256.0L);
  cfg.ram_budget_fraction = 0.30L;
  cfg.estimated_frames = 174720.0L;
  cfg.feature_dim = 4656;
  cfg.dense_memory_estimate =
      6.0L * cfg.feature_dim * cfg.feature_dim * sizeof(double);
  cfg.dual_memory_estimate =
      3.0L * cfg.estimated_frames * cfg.feature_dim * sizeof(double);
  cfg.need_projection_cache = true;

  const auto decision = select_tica_backend_for_policy(cfg);

  require_true(decision.selected_backend == TicaBackend::ExactDual,
               "explicit exact_dual request should be honored");
  require_true(decision.large_matrix_limit_exceeded,
               "large explicit request should still report diagnostic flag");
}

} // namespace

int main() {
  large_pairwise_case_prefers_randomized_even_when_ram_is_available();
  small_case_keeps_existing_exact_dense_behavior();
  projection_cache_still_rewrites_small_exact_dense_to_exact_dual();
  explicit_backend_requests_are_honored();
  std::cout << "PASS titan_tica_backend_policy_smoke\n";
  return 0;
}
