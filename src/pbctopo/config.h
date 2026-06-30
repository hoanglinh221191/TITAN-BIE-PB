#ifndef TITAN_PBCTOPO_CONFIG_H
#define TITAN_PBCTOPO_CONFIG_H

struct arguments;

namespace titan_pbctopo {

enum class AnchorBackend {
  Auto,
  Prim,
  SpatialAll,
};

constexpr int kDefaultExactComponentLimit = 6;
constexpr int kMaxExactComponentLimit = 16;
constexpr int kDefaultBeamWidth = 12;
constexpr int kMaxBeamWidth = 4096;

struct RuntimeCutoffs {
  double bond_cutoff = 2.15;
  double clash_cutoff = 1.05;
  double pair_contact_cutoff = 4.5;
  double continuity_cutoff = 20.0;
};

RuntimeCutoffs resolve_runtime_cutoffs(const arguments &Arguments);
AnchorBackend resolve_anchor_backend(const arguments &Arguments);
const char *anchor_backend_name(AnchorBackend backend);

} // namespace titan_pbctopo

#endif
