#include "titan_hbond_kinetics.h"
#include "titan_interaction_primitives.h"

#include <cstdlib>
#include <iostream>

namespace {

using titan_hbond_kinetics::StreamingKinetics;

[[noreturn]] void fail(const char *message) {
  std::cerr << "[titan_hbond_kinetics_smoke] " << message << "\n";
  std::exit(1);
}

} // namespace

int main() {
  if (!titan_interaction::is_hbond_donor("OG", "SER"))
    fail("SER OG should be recognized as an H-bond donor");
  if (!titan_interaction::is_hbond_acceptor("OD1", "ASP"))
    fail("ASP OD1 should be recognized as an H-bond acceptor");
  if (!titan_interaction::is_hbond_donor("N6", "A"))
    fail("A:N6 should be recognized as a nucleic-base donor");
  if (!titan_interaction::is_hbond_acceptor("O4", "U"))
    fail("U:O4 should be recognized as a nucleic-base acceptor");
  if (!titan_interaction::is_acidic("OE1", "GLU"))
    fail("GLU OE1 should be recognized as acidic");
  if (!titan_interaction::is_basic("NZ", "LYS"))
    fail("LYS NZ should be recognized as basic");

  StreamingKinetics kinetics(2);
  kinetics.process_active_ids({0});
  kinetics.process_active_ids({0});
  kinetics.process_active_ids({});
  kinetics.process_active_ids({0, 1, 1});
  kinetics.finalize();

  const auto &origins = kinetics.origin_counts();
  const auto &intermittent = kinetics.intermittent_hits();
  const auto &continuous = kinetics.continuous_hits();
  if (origins.size() != 3 || intermittent.size() != 3 || continuous.size() != 3)
    fail("maxlag=2 should expose lag bins 0..2");
  if (origins[0] != 4 || intermittent[0] != 4 || continuous[0] != 4)
    fail("lag-0 normalization counts are incorrect");
  if (origins[1] != 2 || intermittent[1] != 1 || continuous[1] != 1)
    fail("lag-1 continuous/intermittent counts mismatch");
  if (origins[2] != 2 || intermittent[2] != 1 || continuous[2] != 0)
    fail("lag-2 continuous/intermittent counts mismatch");

  const auto &stats = kinetics.key_stats();
  if (stats.size() < 2)
    fail("two distinct primitive ids should remain separate in the accumulator");
  if (stats[0].total_active_frames != 3 || stats[0].event_count != 2 ||
      stats[0].total_lifetime_frames != 3 || stats[0].max_lifetime_frames != 2) {
    fail("primary-key lifetime statistics mismatch");
  }
  if (stats[1].total_active_frames != 1 || stats[1].event_count != 1 ||
      stats[1].total_lifetime_frames != 1 || stats[1].max_lifetime_frames != 1) {
    fail("secondary-key lifetime statistics mismatch");
  }

  const auto &hist = kinetics.lifetime_histogram();
  auto one_it = hist.find(1);
  auto two_it = hist.find(2);
  if (one_it == hist.end() || one_it->second != 2 || two_it == hist.end() ||
      two_it->second != 1) {
    fail("lifetime histogram contents mismatch");
  }

  std::cout << "titan_hbond_kinetics_smoke PASS\n";
  return 0;
}
