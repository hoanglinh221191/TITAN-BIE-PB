#include "pull_direction_pathway.h"

#include "MOL.h"
#include "functions.h"
#include "pathway_de_variants.h"
#include "titan_ui.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <span>
#include <sstream>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

constexpr double kPathTol = 1.0e-8;
constexpr std::size_t kNearestResidueEdges = 3;
constexpr double kOutwardStopTol = 1.0e-6;
constexpr double kMinProgressSeedGapNorm = 0.12;
constexpr double kMinPatternStep = 0.05;
constexpr std::size_t kPathPatternMaxSweeps = 32;
constexpr double kPathwayProgressDedupEps = 1.0e-6;
constexpr double kPathwayDirectMinCellNorm = 1.0e-4;

using Point3 = std::array<double, 3>;

struct Interval {
  double begin{0.0};
  double end{0.0};
};

struct PathNode {
  std::size_t chainnum{0};
  std::size_t respos{0};
  std::string chain{" "};
  std::string res_name{"UNK"};
  Point3 com{0.0, 0.0, 0.0};
  std::size_t first_local_atom{0};
};

struct PathSegment {
  double s0{0.0};
  double length{0.0};
  Point3 a{0.0, 0.0, 0.0};
  Point3 d{0.0, 0.0, 0.0};
};

struct PathContext {
  std::size_t group_index{0};
  std::string group_name;
  std::vector<PathNode> all_nodes;
  std::vector<PathNode> principal_nodes;
  std::vector<PathSegment> segments;
  std::vector<Interval> valid_intervals;
  std::size_t pruned_component_nodes{0};
  std::size_t pruned_branch_nodes{0};
  bool used_geometric_edges{false};
  double total_length{0.0};
};

struct PathResult {
  bool ok{false};
  std::string error;
  PathContext path;
};

struct CandidateResult {
  cache_energy energy{};
  double raw_score{0.0};
  double progress{0.0};
  double angle_deg{0.0};
  force vector{};
  bool valid{false};
};

struct PathProjection {
  double progress{0.0};
  double distance{0.0};
  Point3 point{0.0, 0.0, 0.0};
};

struct PathwayState {
  PathProjection current_projection{};
  std::vector<Interval> forward_intervals;
  PathwayMode mode{PathwayMode::Advance};
};

struct PathwayWorkerContext {
  MOL scratch{};
  cache_sigma_epsilon eval_cache{};
  std::size_t eval_quota_remaining{0};
  std::size_t evaluations_done{0};
};

struct PathwayEvalBudget {
  std::size_t max_evals{0};
  std::size_t block_size{8};
  std::atomic_size_t reserved{0};
};

struct PathwayDirectStats {
  bool enabled{false};
  bool used_fallback_budget{false};
  bool stopped_by_eval_budget{false};
  bool stopped_by_local_threshold{false};
  bool stopped_by_stagnation{false};
  std::size_t effective_eval_budget{0};
  std::size_t expansions{0};
  std::size_t batches{0};
  std::size_t max_batch_cells{0};
  std::size_t cells_remaining{0};
  std::size_t evaluated_candidates{0};
  std::size_t best_expansion{0};
  std::size_t phase_evaluations_used{0};
  double best_cell_size{0.0};
  double best_half_width_norm{0.0};
  double wall_seconds{0.0};
};

struct PathwayDirectCell {
  double center_norm{0.5};
  double half_width_norm{0.5};
  PathwayCandidateState state{};
  std::size_t depth{0};
  std::size_t id{0};
};

struct PathwayDirectDomain {
  std::vector<Interval> intervals;
  std::vector<double> prefix_lengths;
  double total_valid_length{0.0};
};

struct PathwayReportRow {
  PathwayCandidateState representative{};
  std::size_t representative_index{0};
  std::size_t count{0};
};

struct ResidueKey {
  std::size_t chainnum{0};
  std::size_t respos{0};
  bool operator==(const ResidueKey &rhs) const noexcept {
    return chainnum == rhs.chainnum && respos == rhs.respos;
  }
};

struct ResidueKeyHash {
  std::size_t operator()(const ResidueKey &key) const noexcept {
    return (key.chainnum << 32U) ^ key.respos;
  }
};

struct ResidueAcc {
  PathNode node;
  Point3 weighted{0.0, 0.0, 0.0};
  double mass{0.0};
};

class DisjointSet {
public:
  explicit DisjointSet(std::size_t n) : parent(n), rank(n, 0) {
    std::iota(parent.begin(), parent.end(), 0);
  }
  std::size_t find(std::size_t x) {
    if (parent[x] != x)
      parent[x] = find(parent[x]);
    return parent[x];
  }
  bool unite(std::size_t a, std::size_t b) {
    a = find(a);
    b = find(b);
    if (a == b)
      return false;
    if (rank[a] < rank[b])
      std::swap(a, b);
    parent[b] = a;
    if (rank[a] == rank[b])
      ++rank[a];
    return true;
  }

private:
  std::vector<std::size_t> parent;
  std::vector<int> rank;
};

Point3 sub_pt(const Point3 &lhs, const Point3 &rhs) {
  return {lhs[0] - rhs[0], lhs[1] - rhs[1], lhs[2] - rhs[2]};
}

Point3 add_pt(const Point3 &lhs, const Point3 &rhs) {
  return {lhs[0] + rhs[0], lhs[1] + rhs[1], lhs[2] + rhs[2]};
}

Point3 scale_pt(const Point3 &p, double scale) {
  return {p[0] * scale, p[1] * scale, p[2] * scale};
}

double dot_pt(const Point3 &lhs, const Point3 &rhs) {
  return lhs[0] * rhs[0] + lhs[1] * rhs[1] + lhs[2] * rhs[2];
}

double norm_sq(const Point3 &p) { return dot_pt(p, p); }

double dist_sq(const Point3 &lhs, const Point3 &rhs) {
  return norm_sq(sub_pt(lhs, rhs));
}

double dist_pt(const Point3 &lhs, const Point3 &rhs) {
  return std::sqrt(dist_sq(lhs, rhs));
}

Point3 mol_com(const MOL &mol, const Multi_MOL_parameters &parameters) {
  Point3 com{0.0, 0.0, 0.0};
  double mass_sum = 0.0;
  for (std::size_t i = 0; i < mol.size(); ++i) {
    const double mass = parameters[i].get_mass();
    com[0] += mol.pos.x[i] * mass;
    com[1] += mol.pos.y[i] * mass;
    com[2] += mol.pos.z[i] * mass;
    mass_sum += mass;
  }
  if (mass_sum <= 0.0)
    return {0.0, 0.0, 0.0};
  return scale_pt(com, 1.0 / mass_sum);
}

force force_from_points(const Point3 &from, const Point3 &to) {
  return {to[0] - from[0], to[1] - from[1], to[2] - from[2]};
}

double angle_to_previous_vector(const force &vector, const force &old_vector) {
  const double old_mod = old_vector.modul();
  const double new_mod = vector.modul();
  if (old_mod < 1.0e-12 || new_mod < 1.0e-12)
    return 0.0;
  force new_scaled{vector};
  force old_scaled{old_vector};
  new_scaled *= old_mod;
  old_scaled *= new_mod;
  return 360.0 / 3.14159265358979323846 *
         std::atan2(force(new_scaled - old_scaled).modul(),
                    force(new_scaled + old_scaled).modul());
}

std::size_t clamp_pattern_tail_budget(std::size_t total_budget) {
  if (total_budget == 0)
    return 0;
  return std::clamp<std::size_t>(total_budget / 8, 48, 160);
}

std::string residue_label(const PathNode &node) {
  std::ostringstream oss;
  oss << node.res_name << node.respos << " (chain " << node.chain << ", id "
      << node.chainnum << ")";
  return oss.str();
}

double base_score(const cache_energy &trial, const cache_energy &reference,
                  int flag_score_function) {
  if (flag_score_function == 1) {
    return trial.elec_ener + trial.vdw_ener - reference.elec_ener -
           reference.vdw_ener +
           2.0 * std::fabs(trial.elec_ener + trial.vdw_ener);
  }
  return -1.0 * std::exp(1.6863 * (-trial.elec_ener - trial.vdw_ener +
                                   reference.elec_ener + reference.vdw_ener));
}

cache_sigma_epsilon make_eval_cache(const cache_sigma_epsilon &source) {
  cache_sigma_epsilon scratch;
  scratch.charge = source.charge;
  scratch.force_field_name = source.force_field_name;
  scratch.cutoff_dis = source.cutoff_dis;
  scratch.charmm_pair_sigma = source.charmm_pair_sigma;
  scratch.charmm_pair_epsilon = source.charmm_pair_epsilon;
  scratch.charmm_type_count = source.charmm_type_count;
  scratch.charmm_pair_sparse_keys = source.charmm_pair_sparse_keys;
  scratch.charmm_pair_sparse_sigma = source.charmm_pair_sparse_sigma;
  scratch.charmm_pair_sparse_epsilon = source.charmm_pair_sparse_epsilon;
  scratch.martini_pair_sigma = source.martini_pair_sigma;
  scratch.martini_pair_epsilon = source.martini_pair_epsilon;
  scratch.martini_type_count = source.martini_type_count;
  return scratch;
}

void eval_energy_model(const MOL &rec, const MOL &lg,
                       const std::vector<mol_info> &info_rec,
                       const std::vector<mol_info> &info_lg,
                       cache_sigma_epsilon &CachedData, cache_energy &Energy,
                       const arguments &Arguments,
                       const Multi_MOL_parameters &rec_parameters,
                       const Multi_MOL_parameters &lg_parameters) {
  if (Arguments.force_f.find("CHARMM") != std::string::npos ||
      Arguments.force_f.find("AMBER") != std::string::npos) {
#if defined(__AVX__)
    cal_energy_amber_intrin(rec, lg, info_rec, info_lg, CachedData, 0,
                            lg.size(), Energy, 0, Arguments, rec_parameters,
                            lg_parameters);
#else
    cal_energy_amber(rec, lg, info_rec, info_lg, CachedData, 0, lg.size(),
                     Energy, 0, Arguments, rec_parameters, lg_parameters);
#endif
  } else if (Arguments.force_f.find("OPLSAA") != std::string::npos) {
#if defined(__AVX__)
    cal_energy_opls_intrin(rec, lg, info_rec, info_lg, CachedData, 0, lg.size(),
                           Energy, 0, Arguments, rec_parameters, lg_parameters);
#else
    cal_energy_opls(rec, lg, info_rec, info_lg, CachedData, 0, lg.size(),
                    Energy, 0, Arguments, rec_parameters, lg_parameters);
#endif
  } else if (Arguments.force_f.find("MARTINI") != std::string::npos) {
    cal_energy_martini(rec, lg, info_rec, info_lg, CachedData, 0, lg.size(),
                       Energy, 0, Arguments, rec_parameters, lg_parameters);
  } else {
    cal_energy_gromos(rec, lg, info_rec, info_lg, CachedData, 0, lg.size(),
                      Energy, 0, Arguments, rec_parameters, lg_parameters);
  }
}

cache_energy reference_energy(const Multi_MOL &rec, const Multi_MOL &lg,
                              cache_sigma_epsilon &CachedData,
                              const arguments &Arguments,
                              TitanExecutor &executor) {
  const std::size_t worker_count = std::max<std::size_t>(1, Arguments.nt);
  const std::size_t slots =
      std::max<std::size_t>(1, std::min(executor.capacity(), worker_count));
  std::vector<cache_energy> parts(slots);
  const std::size_t lg_size = lg.m_mol[0].size();
  executor.parallel_for(
      lg_size, worker_count, std::max<std::size_t>(1, lg_size / slots),
      [&](std::size_t worker_id, std::size_t begin, std::size_t end) {
        if (Arguments.force_f.find("CHARMM") != std::string::npos ||
            Arguments.force_f.find("AMBER") != std::string::npos) {
#if defined(__AVX__)
          cal_energy_amber_intrin(rec.m_mol[0], lg.m_mol[0], rec.info, lg.info,
                                  CachedData, begin, end, parts[worker_id], 0,
                                  Arguments, rec.parameters, lg.parameters);
#else
          cal_energy_amber(rec.m_mol[0], lg.m_mol[0], rec.info, lg.info,
                           CachedData, begin, end, parts[worker_id], 0,
                           Arguments, rec.parameters, lg.parameters);
#endif
        } else if (Arguments.force_f.find("OPLSAA") != std::string::npos) {
#if defined(__AVX__)
          cal_energy_opls_intrin(rec.m_mol[0], lg.m_mol[0], rec.info, lg.info,
                                 CachedData, begin, end, parts[worker_id], 0,
                                 Arguments, rec.parameters, lg.parameters);
#else
          cal_energy_opls(rec.m_mol[0], lg.m_mol[0], rec.info, lg.info,
                          CachedData, begin, end, parts[worker_id], 0,
                          Arguments, rec.parameters, lg.parameters);
#endif
        } else if (Arguments.force_f.find("MARTINI") != std::string::npos) {
          cal_energy_martini(rec.m_mol[0], lg.m_mol[0], rec.info, lg.info,
                             CachedData, begin, end, parts[worker_id], 0,
                             Arguments, rec.parameters, lg.parameters);
        } else {
          cal_energy_gromos(rec.m_mol[0], lg.m_mol[0], rec.info, lg.info,
                            CachedData, begin, end, parts[worker_id], 0,
                            Arguments, rec.parameters, lg.parameters);
        }
      });
  cache_energy total;
  for (const auto &part : parts)
    total += part;
  return total;
}

void add_edge(std::map<std::pair<std::size_t, std::size_t>, double> &edges,
              std::size_t lhs, std::size_t rhs, double weight) {
  if (lhs == rhs)
    return;
  if (lhs > rhs)
    std::swap(lhs, rhs);
  const auto key = std::make_pair(lhs, rhs);
  const auto it = edges.find(key);
  if (it == edges.end())
    edges.emplace(key, weight);
  else
    it->second = std::min(it->second, weight);
}

std::vector<std::vector<std::pair<std::size_t, double>>>
build_adj(std::size_t n,
          const std::map<std::pair<std::size_t, std::size_t>, double> &edges) {
  std::vector<std::vector<std::pair<std::size_t, double>>> adj(n);
  for (const auto &[edge, weight] : edges) {
    adj[edge.first].push_back({edge.second, weight});
    adj[edge.second].push_back({edge.first, weight});
  }
  return adj;
}

std::vector<std::size_t>
component_nodes(std::size_t start,
                const std::vector<std::vector<std::pair<std::size_t, double>>> &adj) {
  std::vector<std::size_t> nodes;
  if (start >= adj.size())
    return nodes;
  std::vector<unsigned char> seen(adj.size(), 0);
  std::vector<std::size_t> stack{start};
  seen[start] = 1;
  while (!stack.empty()) {
    const std::size_t node = stack.back();
    stack.pop_back();
    nodes.push_back(node);
    for (const auto &[next, weight] : adj[node]) {
      (void)weight;
      if (seen[next] != 0)
        continue;
      seen[next] = 1;
      stack.push_back(next);
    }
  }
  return nodes;
}

std::vector<Interval> merge_intervals(std::vector<Interval> intervals) {
  if (intervals.empty())
    return intervals;
  std::sort(intervals.begin(), intervals.end(),
            [](const Interval &lhs, const Interval &rhs) {
              if (lhs.begin != rhs.begin)
                return lhs.begin < rhs.begin;
              return lhs.end < rhs.end;
            });
  std::vector<Interval> merged;
  merged.push_back(intervals.front());
  for (std::size_t i = 1; i < intervals.size(); ++i) {
    if (intervals[i].begin <= merged.back().end + 1.0e-9)
      merged.back().end = std::max(merged.back().end, intervals[i].end);
    else
      merged.push_back(intervals[i]);
  }
  return merged;
}

std::vector<Interval> within_radius(const PathSegment &segment,
                                    const Point3 &center, double radius_sq) {
  std::vector<Interval> out;
  const double a = norm_sq(segment.d);
  if (a <= kPathTol)
    return out;
  const Point3 m = sub_pt(segment.a, center);
  const double b = 2.0 * dot_pt(m, segment.d);
  const double c = norm_sq(m) - radius_sq;
  const double disc = b * b - 4.0 * a * c;
  if (disc < 0.0) {
    if (c <= 0.0)
      out.push_back({0.0, 1.0});
    return out;
  }
  const double sqrt_disc = std::sqrt(std::max(0.0, disc));
  const double t1 = (-b - sqrt_disc) / (2.0 * a);
  const double t2 = (-b + sqrt_disc) / (2.0 * a);
  const double begin = std::clamp(std::min(t1, t2), 0.0, 1.0);
  const double end = std::clamp(std::max(t1, t2), 0.0, 1.0);
  if (begin <= end + 1.0e-12)
    out.push_back({begin, end});
  return out;
}

std::vector<Interval> subtract_interval(const Interval &outer,
                                        const Interval &inner) {
  if (inner.end <= outer.begin || inner.begin >= outer.end)
    return {outer};
  std::vector<Interval> out;
  if (inner.begin > outer.begin + 1.0e-9)
    out.push_back({outer.begin, std::min(inner.begin, outer.end)});
  if (inner.end < outer.end - 1.0e-9)
    out.push_back({std::max(inner.end, outer.begin), outer.end});
  return out;
}

std::vector<Interval> shell_intervals(const PathSegment &segment,
                                      const Point3 &center, double min_radius,
                                      double max_radius) {
  auto outer = within_radius(segment, center, max_radius * max_radius);
  if (outer.empty())
    return outer;
  if (min_radius <= 0.0)
    return outer;
  const double min_eff = std::max(0.0, min_radius - 1.0e-9);
  const auto inner = within_radius(segment, center, min_eff * min_eff);
  if (inner.empty())
    return outer;
  std::vector<Interval> out;
  for (const auto &o : outer) {
    std::vector<Interval> pending{o};
    for (const auto &i : inner) {
      std::vector<Interval> next;
      for (const auto &candidate : pending) {
        const auto pieces = subtract_interval(candidate, i);
        next.insert(next.end(), pieces.begin(), pieces.end());
      }
      pending = std::move(next);
    }
    out.insert(out.end(), pending.begin(), pending.end());
  }
  return out;
}

PathResult build_path(const Multi_MOL &rec, const Multi_MOL &lg,
                      groups_cal &groups, const arguments &Arguments) {
  PathResult result;
  result.path.group_index = groups.group_pathway_index();
  result.path.group_name = groups.group_name(result.path.group_index);

  std::unordered_map<std::size_t, std::size_t> rec_local_from_global;
  const std::size_t group1 = groups.group_1_index();
  rec_local_from_global.reserve(groups.groups_atom_number(group1));
  for (std::size_t i = 0; i < groups.groups_atom_number(group1); ++i)
    rec_local_from_global.emplace(groups.atom_index(group1, i), i);

  std::unordered_map<ResidueKey, std::size_t, ResidueKeyHash> seen;
  std::vector<ResidueAcc> residues;
  residues.reserve(groups.groups_atom_number(result.path.group_index));
  for (std::size_t i = 0; i < groups.groups_atom_number(result.path.group_index); ++i) {
    const std::size_t global_atom = groups.atom_index(result.path.group_index, i);
    const auto local_it = rec_local_from_global.find(global_atom);
    if (local_it == rec_local_from_global.end()) {
      result.error = "Pathway atom is not present in routed g1.";
      return result;
    }
    const std::size_t local_atom = local_it->second;
    const mol_info &atom_info = rec.info[local_atom];
    const ResidueKey key{atom_info.get_chainnum(), atom_info.get_respos()};
    auto [it, inserted] = seen.emplace(key, residues.size());
    if (inserted) {
      ResidueAcc acc;
      acc.node.chainnum = atom_info.get_chainnum();
      acc.node.respos = atom_info.get_respos();
      acc.node.chain = atom_info.get_chain();
      acc.node.res_name = atom_info.get_original_res_name();
      acc.node.first_local_atom = local_atom;
      residues.push_back(std::move(acc));
      it->second = residues.size() - 1;
    }
    ResidueAcc &acc = residues[it->second];
    const double mass = rec.parameters[local_atom].get_mass();
    acc.mass += mass;
    acc.weighted[0] += rec.m_mol[0].pos.x[local_atom] * mass;
    acc.weighted[1] += rec.m_mol[0].pos.y[local_atom] * mass;
    acc.weighted[2] += rec.m_mol[0].pos.z[local_atom] * mass;
  }
  for (auto &acc : residues) {
    if (acc.mass <= 0.0)
      continue;
    acc.node.com = scale_pt(acc.weighted, 1.0 / acc.mass);
    result.path.all_nodes.push_back(acc.node);
  }
  if (result.path.all_nodes.size() < 2) {
    result.error = "Pathway group collapses to fewer than 2 residue COM nodes.";
    return result;
  }

  const Point3 lg_start = mol_com(lg.m_mol[0], lg.parameters);
  std::size_t start = 0;
  double best_d2 = std::numeric_limits<double>::max();
  for (std::size_t i = 0; i < result.path.all_nodes.size(); ++i) {
    const double d2 = dist_sq(result.path.all_nodes[i].com, lg_start);
    if (d2 < best_d2) {
      best_d2 = d2;
      start = i;
    }
  }

  std::map<std::pair<std::size_t, std::size_t>, double> edges;
  std::map<std::size_t, std::vector<std::size_t>> by_chain;
  for (std::size_t i = 0; i < result.path.all_nodes.size(); ++i)
    by_chain[result.path.all_nodes[i].chainnum].push_back(i);
  for (auto &[chain, nodes] : by_chain) {
    (void)chain;
    std::sort(nodes.begin(), nodes.end(), [&](std::size_t lhs, std::size_t rhs) {
      if (result.path.all_nodes[lhs].respos != result.path.all_nodes[rhs].respos)
        return result.path.all_nodes[lhs].respos <
               result.path.all_nodes[rhs].respos;
      return lhs < rhs;
    });
    for (std::size_t i = 1; i < nodes.size(); ++i)
      add_edge(edges, nodes[i - 1], nodes[i],
               dist_pt(result.path.all_nodes[nodes[i - 1]].com,
                       result.path.all_nodes[nodes[i]].com));
  }

  auto adj = build_adj(result.path.all_nodes.size(), edges);
  auto component = component_nodes(start, adj);
  if (component.size() < result.path.all_nodes.size()) {
    result.path.used_geometric_edges = true;
    for (std::size_t i = 0; i < result.path.all_nodes.size(); ++i) {
      std::vector<std::pair<double, std::size_t>> nearest;
      for (std::size_t j = 0; j < result.path.all_nodes.size(); ++j) {
        if (i == j)
          continue;
        nearest.push_back({dist_sq(result.path.all_nodes[i].com,
                                   result.path.all_nodes[j].com),
                           j});
      }
      std::sort(nearest.begin(), nearest.end());
      for (std::size_t k = 0; k < std::min(kNearestResidueEdges, nearest.size());
           ++k) {
        const double distance = std::sqrt(nearest[k].first);
        if (distance > Arguments.de.dis_along_path)
          continue;
        add_edge(edges, i, nearest[k].second, distance);
      }
    }
    adj = build_adj(result.path.all_nodes.size(), edges);
    component = component_nodes(start, adj);
  }
  if (component.size() < 2) {
    result.error = "Pathway component near g2 has fewer than 2 nodes.";
    return result;
  }
  result.path.pruned_component_nodes =
      result.path.all_nodes.size() - component.size();

  std::unordered_map<std::size_t, std::size_t> local;
  for (std::size_t i = 0; i < component.size(); ++i)
    local.emplace(component[i], i);

  struct EdgeRec { std::size_t a{}, b{}; double w{}; };
  std::vector<EdgeRec> mst_input;
  for (const auto &[edge, weight] : edges) {
    const auto a = local.find(edge.first);
    const auto b = local.find(edge.second);
    if (a == local.end() || b == local.end())
      continue;
    mst_input.push_back({a->second, b->second, weight});
  }
  std::sort(mst_input.begin(), mst_input.end(),
            [](const EdgeRec &lhs, const EdgeRec &rhs) {
              if (lhs.w != rhs.w)
                return lhs.w < rhs.w;
              if (lhs.a != rhs.a)
                return lhs.a < rhs.a;
              return lhs.b < rhs.b;
            });
  DisjointSet dsu(component.size());
  std::vector<std::vector<std::pair<std::size_t, double>>> tree(component.size());
  for (const auto &edge : mst_input) {
    if (!dsu.unite(edge.a, edge.b))
      continue;
    tree[edge.a].push_back({edge.b, edge.w});
    tree[edge.b].push_back({edge.a, edge.w});
  }

  const std::size_t local_start = local[start];
  std::vector<int> parent(component.size(), -1);
  std::vector<double> geodesic(component.size(), 0.0);
  std::vector<std::size_t> stack{local_start};
  parent[local_start] = static_cast<int>(local_start);
  while (!stack.empty()) {
    const std::size_t node = stack.back();
    stack.pop_back();
    for (const auto &[next, weight] : tree[node]) {
      if (parent[next] != -1)
        continue;
      parent[next] = static_cast<int>(node);
      geodesic[next] = geodesic[node] + weight;
      stack.push_back(next);
    }
  }

  std::size_t local_terminal = local_start;
  for (std::size_t node = 0; node < component.size(); ++node) {
    const auto &candidate = result.path.all_nodes[component[node]];
    const auto &best = result.path.all_nodes[component[local_terminal]];
    const double cand_d2 = dist_sq(candidate.com, lg_start);
    const double best_terminal_d2 = dist_sq(best.com, lg_start);
    if (geodesic[node] > geodesic[local_terminal] + 1.0e-9) {
      local_terminal = node;
      continue;
    }
    if (std::fabs(geodesic[node] - geodesic[local_terminal]) <= 1.0e-9) {
      if (cand_d2 > best_terminal_d2 + 1.0e-9) {
        local_terminal = node;
        continue;
      }
      if (std::fabs(cand_d2 - best_terminal_d2) <= 1.0e-9 &&
          (candidate.chainnum < best.chainnum ||
           (candidate.chainnum == best.chainnum &&
            candidate.respos < best.respos))) {
        local_terminal = node;
      }
    }
  }

  std::vector<std::size_t> principal;
  for (std::size_t cursor = local_terminal;; cursor = static_cast<std::size_t>(parent[cursor])) {
    principal.push_back(component[cursor]);
    if (cursor == local_start)
      break;
  }
  std::reverse(principal.begin(), principal.end());
  for (std::size_t global_id : principal) {
    if (!result.path.principal_nodes.empty() &&
        dist_pt(result.path.principal_nodes.back().com,
                result.path.all_nodes[global_id].com) <= kPathTol)
      continue;
    result.path.principal_nodes.push_back(result.path.all_nodes[global_id]);
  }
  if (result.path.principal_nodes.size() < 2) {
    result.error = "Principal pathway degenerates to fewer than 2 unique nodes.";
    return result;
  }
  result.path.pruned_branch_nodes =
      component.size() - result.path.principal_nodes.size();

  double s = 0.0;
  for (std::size_t i = 1; i < result.path.principal_nodes.size(); ++i) {
    PathSegment segment;
    segment.a = result.path.principal_nodes[i - 1].com;
    segment.d = sub_pt(result.path.principal_nodes[i].com, segment.a);
    segment.length = std::sqrt(norm_sq(segment.d));
    if (segment.length <= kPathTol)
      continue;
    segment.s0 = s;
    s += segment.length;
    result.path.segments.push_back(segment);
  }
  result.path.total_length = s;
  if (result.path.segments.empty() || result.path.total_length <= kPathTol) {
    result.error = "Principal pathway has zero usable arclength.";
    return result;
  }

  std::vector<Interval> valid;
  for (const auto &segment : result.path.segments) {
    auto shell = shell_intervals(segment, lg_start,
                                 Arguments.de.de_min_scan_radius,
                                 Arguments.de.de_scan_radius);
    for (const auto &interval : shell) {
      valid.push_back(
          {segment.s0 + interval.begin * segment.length,
           segment.s0 + interval.end * segment.length});
    }
  }
  result.path.valid_intervals = merge_intervals(std::move(valid));
  if (result.path.valid_intervals.empty()) {
    result.error =
        "Forward pathway does not intersect the requested scan-radius window.";
    return result;
  }

  result.ok = true;
  return result;
}

Point3 point_on_path(const PathContext &path, double s) {
  if (s <= 0.0)
    return path.principal_nodes.front().com;
  if (s >= path.total_length)
    return path.principal_nodes.back().com;
  for (const auto &segment : path.segments) {
    if (s <= segment.s0 + segment.length + 1.0e-9) {
      const double t = (s - segment.s0) / segment.length;
      return add_pt(segment.a, scale_pt(segment.d, t));
    }
  }
  return path.principal_nodes.back().com;
}

double clamp_to_intervals(double value, const std::vector<Interval> &intervals) {
  if (intervals.empty())
    return value;
  if (value <= intervals.front().begin)
    return intervals.front().begin;
  if (value >= intervals.back().end)
    return intervals.back().end;
  for (const auto &interval : intervals) {
    if (value >= interval.begin && value <= interval.end)
      return value;
  }
  double best = intervals.front().begin;
  double best_gap = std::numeric_limits<double>::max();
  for (const auto &interval : intervals) {
    const double left = std::fabs(value - interval.begin);
    if (left < best_gap) {
      best_gap = left;
      best = interval.begin;
    }
    const double right = std::fabs(value - interval.end);
    if (right < best_gap) {
      best_gap = right;
      best = interval.end;
    }
  }
  return best;
}

double sample_interval(std::mt19937 &mt, const std::vector<Interval> &intervals) {
  double total = 0.0;
  for (const auto &interval : intervals)
    total += interval.end - interval.begin;
  std::uniform_real_distribution<double> pick(0.0, total);
  double marker = pick(mt);
  for (const auto &interval : intervals) {
    const double len = interval.end - interval.begin;
    if (marker <= len)
      return interval.begin + marker;
    marker -= len;
  }
  return intervals.back().end;
}

PathProjection project_point_to_path(const PathContext &path, const Point3 &point,
                                     double min_progress = 0.0) {
  PathProjection best;
  best.distance = std::numeric_limits<double>::max();
  best.progress = min_progress;
  best.point = point_on_path(path, min_progress);

  for (const auto &segment : path.segments) {
    const double segment_begin = segment.s0;
    const double segment_end = segment.s0 + segment.length;
    if (segment_end + 1.0e-9 < min_progress)
      continue;

    const double begin_t =
        min_progress <= segment_begin ? 0.0
        : min_progress >= segment_end ? 1.0
                                      : (min_progress - segment_begin) /
                                            segment.length;
    const Point3 offset = sub_pt(point, segment.a);
    const double denom = norm_sq(segment.d);
    double t = denom <= kPathTol ? begin_t : dot_pt(offset, segment.d) / denom;
    t = std::clamp(t, begin_t, 1.0);
    const Point3 projected = add_pt(segment.a, scale_pt(segment.d, t));
    const double distance = dist_pt(point, projected);
    const double progress = segment.s0 + t * segment.length;
    if (distance + 1.0e-9 < best.distance ||
        (std::fabs(distance - best.distance) <= 1.0e-9 &&
         progress < best.progress)) {
      best.distance = distance;
      best.progress = progress;
      best.point = projected;
    }
  }

  if (!std::isfinite(best.distance))
    best.distance = 0.0;
  return best;
}

std::vector<Interval> build_forward_intervals(const PathContext &path,
                                              const Point3 &center,
                                              double current_progress,
                                              double min_radius,
                                              double max_radius) {
  std::vector<Interval> valid;
  for (const auto &segment : path.segments) {
    const double segment_end = segment.s0 + segment.length;
    if (segment_end + 1.0e-9 < current_progress)
      continue;

    auto shell = shell_intervals(segment, center, min_radius, max_radius);
    for (const auto &interval : shell) {
      double begin = segment.s0 + interval.begin * segment.length;
      double end = segment.s0 + interval.end * segment.length;
      begin = std::max(begin, current_progress);
      end = std::min(end, segment_end);
      if (begin <= end + 1.0e-9)
        valid.push_back({begin, end});
    }
  }
  return merge_intervals(std::move(valid));
}

std::string pathway_mode_label(PathwayMode mode) {
  return mode == PathwayMode::Recovery ? "Recovery" : "Advance";
}

bool old_vector_points_outward(const PathContext &path, const Point3 &start_com,
                               const PathProjection &current_projection,
                               const force &old_vector) {
  const Point3 future = add_pt(
      start_com, {old_vector.Fx, old_vector.Fy, old_vector.Fz});
  const PathProjection future_projection =
      project_point_to_path(path, future, current_projection.progress);
  return future_projection.distance >
         current_projection.distance + kOutwardStopTol;
}

std::size_t total_pathway_evaluations_used(
    const std::vector<PathwayWorkerContext> &workers) {
  std::size_t total = 0;
  for (const auto &worker : workers)
    total += worker.evaluations_done;
  return total;
}

bool reserve_eval_budget(PathwayEvalBudget &budget,
                         PathwayWorkerContext &worker) {
  if (budget.max_evals == 0) {
    ++worker.evaluations_done;
    return true;
  }
  if (worker.eval_quota_remaining != 0) {
    --worker.eval_quota_remaining;
    ++worker.evaluations_done;
    return true;
  }

  std::size_t current = budget.reserved.load(std::memory_order_relaxed);
  while (current < budget.max_evals) {
    const std::size_t grant =
        std::min(budget.block_size, budget.max_evals - current);
    if (budget.reserved.compare_exchange_weak(current, current + grant,
                                              std::memory_order_relaxed,
                                              std::memory_order_relaxed)) {
      worker.eval_quota_remaining = grant - 1;
      ++worker.evaluations_done;
      return true;
    }
  }
  return false;
}

std::size_t pathway_parallel_grain(std::size_t total_items,
                                   std::size_t requested_workers,
                                   std::size_t ligand_atom_count) {
  const std::size_t workers = std::max<std::size_t>(1, requested_workers);
  if (total_items <= workers)
    return 1;
  const std::size_t min_items_per_worker =
      ligand_atom_count >= 2048 ? 4
      : ligand_atom_count >= 1024 ? 6
      : ligand_atom_count >= 256  ? 8
                                : 12;
  return std::max<std::size_t>(
      min_items_per_worker,
      (total_items + workers * 2 - 1) / (workers * 2));
}

std::size_t pathway_eval_budget_block_size(std::size_t population_size,
                                           std::size_t worker_slots) {
  const std::size_t workers = std::max<std::size_t>(1, worker_slots);
  const std::size_t per_worker =
      std::max<std::size_t>(1, (population_size + workers - 1) / workers);
  return std::clamp<std::size_t>(std::max<std::size_t>(8, per_worker / 2), 8,
                                 64);
}

bool pathway_candidate_improves(const PathwayCandidateState &lhs,
                                const PathwayCandidateState &rhs,
                                PathwayMode mode) {
  return pathway_candidate_better(lhs, 0, rhs, 0, mode);
}

std::vector<PathwayReportRow> build_grouped_report_rows(
    const std::vector<PathwayCandidateState> &population, PathwayMode mode) {
  std::vector<std::pair<double, std::size_t>> by_progress;
  by_progress.reserve(population.size());
  for (std::size_t i = 0; i < population.size(); ++i) {
    if (population[i].valid)
      by_progress.push_back({population[i].progress, i});
  }
  std::sort(by_progress.begin(), by_progress.end(),
            [](const auto &lhs, const auto &rhs) {
              if (std::fabs(lhs.first - rhs.first) > kPathwayProgressDedupEps)
                return lhs.first < rhs.first;
              return lhs.second < rhs.second;
            });

  std::vector<PathwayReportRow> rows;
  for (const auto &[progress, index] : by_progress) {
    (void)progress;
    if (rows.empty() ||
        std::fabs(population[index].progress -
                  rows.back().representative.progress) > kPathwayProgressDedupEps) {
      PathwayReportRow row;
      row.representative = population[index];
      row.representative_index = index;
      row.count = 1;
      rows.push_back(row);
      continue;
    }
    auto &row = rows.back();
    ++row.count;
    if (pathway_candidate_better(population[index], index, row.representative,
                                 row.representative_index, mode)) {
      row.representative = population[index];
      row.representative_index = index;
    }
  }

  std::stable_sort(rows.begin(), rows.end(),
                   [&](const PathwayReportRow &lhs,
                       const PathwayReportRow &rhs) {
                     return pathway_candidate_better(
                         lhs.representative, lhs.representative_index,
                         rhs.representative, rhs.representative_index, mode);
                   });
  return rows;
}

std::string variant_pattern_section_title(PullDeVariant variant) {
  switch (variant) {
  case PullDeVariant::NlShadeLbc:
    return "NL-SHADE-LBC Pattern Refinement";
  case PullDeVariant::Jso:
    return "JSO Pattern Refinement";
  default:
    return "Pathway Pattern Refinement";
  }
}

PathwayState prepare_pathway_state(const PathContext &path,
                                   const Point3 &start_com,
                                   const arguments &Arguments) {
  PathwayState state;
  state.current_projection = project_point_to_path(path, start_com);
  state.forward_intervals = build_forward_intervals(
      path, start_com, state.current_projection.progress,
      Arguments.de.de_min_scan_radius, Arguments.de.de_scan_radius);
  state.mode = state.current_projection.distance > Arguments.de.dis_along_path
                   ? PathwayMode::Recovery
                   : PathwayMode::Advance;
  return state;
}

PathwayCandidateState evaluate_candidate(
    const PathContext &path, const PathwayState &state, double progress,
    const Point3 &start_com, const Multi_MOL &rec, const Multi_MOL &lg,
    PathwayWorkerContext &worker, const cache_energy &energy_ref,
    const arguments &Arguments, const de_arg &DE_arguments,
    PathwayEvalBudget &eval_budget) {
  PathwayCandidateState result;
  result.progress = progress;
  result.delta_progress = std::max(
      0.0, progress - state.current_projection.progress);
  const Point3 target = point_on_path(path, progress);
  result.vector = force_from_points(start_com, target);
  result.pull_distance = result.vector.modul();
  result.angle_deg =
      angle_to_previous_vector(result.vector, DE_arguments.old_vector);
  if (DE_arguments.flag_old_vec &&
      result.angle_deg > DE_arguments.restrain_angle) {
    return result;
  }
  if (!reserve_eval_budget(eval_budget, worker)) {
    return result;
  }
  worker.scratch = lg.m_mol[0];
  worker.scratch.translate(result.vector.Fx, result.vector.Fy, result.vector.Fz);
  eval_energy_model(rec.m_mol[0], worker.scratch, rec.info, lg.info,
                    worker.eval_cache, result.energy, Arguments, rec.parameters,
                    lg.parameters);
  result.raw_score =
      base_score(result.energy, energy_ref, Arguments.de.de_flag_score_function);
  result.energy_key = pathway_energy_key(result.raw_score);
  result.valid = std::isfinite(result.raw_score) &&
                 std::isfinite(result.energy_key) &&
                 std::isfinite(result.pull_distance);
  return result;
}

PullDeBatchEvaluation evaluate_candidate_batch(
    std::span<const double> proposals, std::span<PathwayCandidateState> results,
    const PathContext &path, const PathwayState &state, const Point3 &start_com,
    const Multi_MOL &rec, const Multi_MOL &lg,
    std::vector<PathwayWorkerContext> &workers, const cache_energy &energy_ref,
    const arguments &Arguments, const de_arg &DE_arguments,
    PathwayEvalBudget &eval_budget, TitanExecutor &executor,
    std::size_t requested_workers) {
  PullDeBatchEvaluation batch;
  if (proposals.empty())
    return batch;

  std::vector<std::pair<double, std::size_t>> ordered;
  ordered.reserve(proposals.size());
  for (std::size_t i = 0; i < proposals.size(); ++i)
    ordered.push_back({proposals[i], i});
  std::sort(ordered.begin(), ordered.end(),
            [](const auto &lhs, const auto &rhs) {
              if (std::fabs(lhs.first - rhs.first) > kPathwayProgressDedupEps)
                return lhs.first < rhs.first;
              return lhs.second < rhs.second;
            });

  std::vector<double> unique_progress;
  unique_progress.reserve(ordered.size());
  std::vector<std::size_t> proposal_to_unique(proposals.size(), 0);
  for (const auto &[progress, index] : ordered) {
    if (unique_progress.empty() ||
        std::fabs(progress - unique_progress.back()) > kPathwayProgressDedupEps) {
      unique_progress.push_back(progress);
    }
    proposal_to_unique[index] = unique_progress.size() - 1;
  }
  std::vector<PathwayCandidateState> unique_results(unique_progress.size());

  const std::size_t before = total_pathway_evaluations_used(workers);
  executor.parallel_for(
      unique_progress.size(), requested_workers,
      pathway_parallel_grain(unique_progress.size(), requested_workers,
                             lg.m_mol[0].size()),
      [&](std::size_t worker_id, std::size_t begin, std::size_t end) {
        auto &worker = workers[worker_id];
        for (std::size_t i = begin; i < end; ++i) {
          unique_results[i] =
              evaluate_candidate(path, state, unique_progress[i], start_com, rec,
                                 lg, worker, energy_ref, Arguments,
                                 DE_arguments, eval_budget);
        }
      });
  const std::size_t after = total_pathway_evaluations_used(workers);
  for (std::size_t i = 0; i < proposals.size(); ++i)
    results[i] = unique_results[proposal_to_unique[i]];
  batch.evaluations_used = after - before;
  batch.budget_exhausted =
      (eval_budget.max_evals != 0 && after >= eval_budget.max_evals);
  return batch;
}

PathwayPatternStats refine_pathway_candidate_pattern(
    PathwayCandidateState &best_candidate, const PathwayProgressDomain &domain,
    const PathwayBatchEvaluator &evaluator, double initial_step,
    PathwayMode mode) {
  PathwayPatternStats stats;
  stats.enabled = true;
  stats.initial_score = best_candidate.raw_score;
  const auto started_at = std::chrono::steady_clock::now();
  double step = std::max(initial_step, kMinPatternStep);

  while (stats.sweeps < kPathPatternMaxSweeps && step > kMinPatternStep) {
    std::array<double, 2> proposals{
        domain.sanitize(best_candidate.progress + step),
        domain.sanitize(best_candidate.progress - step)};
    std::array<PathwayCandidateState, 2> trials{};
    const PullDeBatchEvaluation eval_info = evaluator(
        std::span<const double>(proposals), std::span<PathwayCandidateState>(trials));
    stats.evaluations_used += eval_info.evaluations_used;
    if (eval_info.budget_exhausted)
      stats.stopped_by_eval_budget = true;
    ++stats.sweeps;

    std::size_t best_trial = trials.size();
    for (std::size_t i = 0; i < trials.size(); ++i) {
      if (!trials[i].valid)
        continue;
      if (pathway_candidate_improves(trials[i], best_candidate, mode)) {
        best_trial = i;
        break;
      }
    }

    if (best_trial < trials.size()) {
      best_candidate = trials[best_trial];
      stats.improved = true;
      ++stats.accepted_moves;
      continue;
    }

    if (stats.stopped_by_eval_budget)
      break;
    step *= 0.5;
  }

  stats.final_score = best_candidate.raw_score;
  stats.wall_seconds =
      std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                    started_at)
          .count();
  return stats;
}

PathwayPatternStats refine_pathway_variant_multi_seed_pattern(
    PathwayCandidateState &best_candidate,
    const std::vector<PathwayCandidateState> &candidate_pool,
    const PathwayProgressDomain &domain, const PathwayBatchEvaluator &evaluator,
    const PullDeVariantStats &variant_stats, double total_length,
    PathwayMode mode) {
  PathwayPatternStats aggregate;
  aggregate.enabled = true;
  aggregate.initial_score = best_candidate.raw_score;
  const auto started_at = std::chrono::steady_clock::now();
  const double seed_gap =
      std::max(0.1, total_length * kMinProgressSeedGapNorm);
  const auto seed_indices =
      select_pathway_pattern_seed_indices(candidate_pool, 3, seed_gap, mode);
  aggregate.seeds_attempted = seed_indices.empty() ? 1 : seed_indices.size();
  aggregate.distinct_seed_basins = aggregate.seeds_attempted;
  const double initial_step = pathway_variant_pattern_initial_step(
      variant_stats, std::max(0.1, domain.max_progress - domain.min_progress));

  if (seed_indices.empty()) {
    PathwayCandidateState local_best = best_candidate;
    PathwayPatternStats local_stats = refine_pathway_candidate_pattern(
        local_best, domain, evaluator, initial_step, mode);
    aggregate.sweeps = local_stats.sweeps;
    aggregate.accepted_moves = local_stats.accepted_moves;
    aggregate.evaluations_used = local_stats.evaluations_used;
    aggregate.stopped_by_eval_budget = local_stats.stopped_by_eval_budget;
    if (pathway_candidate_improves(local_best, best_candidate, mode)) {
      best_candidate = local_best;
      aggregate.improved = true;
    }
  } else {
    for (const std::size_t seed_index : seed_indices) {
      PathwayCandidateState local_best = candidate_pool[seed_index];
      PathwayPatternStats local_stats = refine_pathway_candidate_pattern(
          local_best, domain, evaluator, initial_step, mode);
      aggregate.sweeps += local_stats.sweeps;
      aggregate.accepted_moves += local_stats.accepted_moves;
      aggregate.evaluations_used += local_stats.evaluations_used;
      if (local_stats.improved &&
          pathway_candidate_improves(local_best, best_candidate, mode)) {
        best_candidate = local_best;
        aggregate.improved = true;
      }
      if (local_stats.stopped_by_eval_budget) {
        aggregate.stopped_by_eval_budget = true;
        break;
      }
    }
  }

  aggregate.final_score = best_candidate.raw_score;
  aggregate.wall_seconds =
      std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                    started_at)
          .count();
  return aggregate;
}

std::size_t effective_pathway_direct_eval_budget(const de_arg &DE_arguments) {
  if (DE_arguments.pull_max_evals != 0)
    return DE_arguments.pull_max_evals;
  return std::max<std::size_t>(DE_arguments.pop_num * DE_arguments.ite_num,
                               1024);
}

double pathway_direct_local_threshold_size(const PathwayDirectDomain &domain) {
  const double span = std::max(1.0e-6, domain.total_valid_length);
  return std::max(kPathwayDirectMinCellNorm, 0.08 / span);
}

PathwayDirectDomain build_pathway_direct_domain(
    const std::vector<Interval> &intervals) {
  PathwayDirectDomain domain;
  domain.intervals = intervals;
  domain.prefix_lengths.reserve(intervals.size());
  double prefix = 0.0;
  for (const auto &interval : domain.intervals) {
    domain.prefix_lengths.push_back(prefix);
    prefix += std::max(0.0, interval.end - interval.begin);
  }
  domain.total_valid_length = prefix;
  return domain;
}

double pathway_direct_normalized_to_progress(const PathwayDirectDomain &domain,
                                             double normalized) {
  if (domain.intervals.empty())
    return 0.0;
  normalized = std::clamp(normalized, 0.0, 1.0);
  if (domain.total_valid_length <= kPathTol)
    return domain.intervals.front().begin;
  double cursor = normalized * domain.total_valid_length;
  for (std::size_t i = 0; i < domain.intervals.size(); ++i) {
    const double length =
        std::max(0.0, domain.intervals[i].end - domain.intervals[i].begin);
    if (cursor <= length + kPathTol || i + 1 == domain.intervals.size())
      return domain.intervals[i].begin + std::min(cursor, length);
    cursor -= length;
  }
  return domain.intervals.back().end;
}

double pathway_direct_cell_size_metric(const PathwayDirectCell &cell) {
  return cell.half_width_norm;
}

bool pathway_direct_cell_preferred(const PathwayDirectCell &lhs,
                                   const PathwayDirectCell &rhs,
                                   PathwayMode mode) {
  const double lhs_size = pathway_direct_cell_size_metric(lhs);
  const double rhs_size = pathway_direct_cell_size_metric(rhs);
  if (lhs.state.valid != rhs.state.valid)
    return lhs.state.valid;
  if (lhs.state.valid && rhs.state.valid &&
      !pathway_candidate_equivalent(lhs.state, rhs.state, mode)) {
    return pathway_candidate_better(lhs.state, lhs.id, rhs.state, rhs.id, mode);
  }
  if (std::fabs(lhs_size - rhs_size) > 1.0e-12)
    return lhs_size > rhs_size;
  if (lhs.depth != rhs.depth)
    return lhs.depth < rhs.depth;
  return lhs.id < rhs.id;
}

std::size_t select_pathway_direct_cell_index(
    const std::vector<PathwayDirectCell> &frontier, PathwayMode mode) {
  double largest_size = 0.0;
  for (const auto &cell : frontier)
    largest_size = std::max(largest_size, pathway_direct_cell_size_metric(cell));
  const double size_cutoff = largest_size * 0.5;

  std::size_t chosen = frontier.size();
  for (std::size_t i = 0; i < frontier.size(); ++i) {
    if (pathway_direct_cell_size_metric(frontier[i]) + 1.0e-12 < size_cutoff)
      continue;
    if (chosen == frontier.size() ||
        pathway_direct_cell_preferred(frontier[i], frontier[chosen], mode)) {
      chosen = i;
    }
  }
  if (chosen != frontier.size())
    return chosen;

  chosen = 0;
  for (std::size_t i = 1; i < frontier.size(); ++i) {
    if (pathway_direct_cell_preferred(frontier[i], frontier[chosen], mode))
      chosen = i;
  }
  return chosen;
}

bool pathway_direct_cell_dominates(const PathwayDirectCell &lhs,
                                   const PathwayDirectCell &rhs,
                                   PathwayMode mode) {
  const double lhs_size = pathway_direct_cell_size_metric(lhs);
  const double rhs_size = pathway_direct_cell_size_metric(rhs);
  const bool no_worse_size = lhs_size + 1.0e-12 >= rhs_size;
  const bool lhs_better =
      pathway_candidate_better(lhs.state, lhs.id, rhs.state, rhs.id, mode);
  const bool rhs_better =
      pathway_candidate_better(rhs.state, rhs.id, lhs.state, lhs.id, mode);
  const bool no_worse_score = lhs_better || !rhs_better;
  const bool strictly_better = lhs_size > rhs_size + 1.0e-12 || lhs_better;
  return no_worse_size && no_worse_score && strictly_better;
}

std::vector<std::size_t> select_pathway_direct_cell_batch_indices(
    const std::vector<PathwayDirectCell> &frontier, std::size_t max_cells,
    PathwayMode mode) {
  std::vector<std::size_t> pareto;
  pareto.reserve(frontier.size());
  for (std::size_t i = 0; i < frontier.size(); ++i) {
    bool dominated = false;
    for (std::size_t j = 0; j < frontier.size(); ++j) {
      if (i == j)
        continue;
      if (pathway_direct_cell_dominates(frontier[j], frontier[i], mode)) {
        dominated = true;
        break;
      }
    }
    if (!dominated)
      pareto.push_back(i);
  }
  if (pareto.empty())
    return {select_pathway_direct_cell_index(frontier, mode)};

  std::sort(pareto.begin(), pareto.end(), [&](std::size_t lhs, std::size_t rhs) {
    if (pathway_candidate_equivalent(frontier[lhs].state, frontier[rhs].state,
                                     mode)) {
      const double lhs_size = pathway_direct_cell_size_metric(frontier[lhs]);
      const double rhs_size = pathway_direct_cell_size_metric(frontier[rhs]);
      if (std::fabs(lhs_size - rhs_size) > 1.0e-12)
        return lhs_size > rhs_size;
      return frontier[lhs].id < frontier[rhs].id;
    }
    return pathway_candidate_better(frontier[lhs].state, frontier[lhs].id,
                                    frontier[rhs].state, frontier[rhs].id,
                                    mode);
  });

  if (pareto.size() > max_cells)
    pareto.resize(max_cells);
  std::sort(pareto.begin(), pareto.end());
  return pareto;
}

PathwayDirectStats run_pathway_direct_search(
    PathwayCandidateState &best_candidate,
    std::vector<PathwayCandidateState> &ranked_candidates, std::size_t &best_stage,
    bool &global_budget_stop, const PathwayDirectDomain &direct_domain,
    const PathwayBatchEvaluator &evaluator, std::size_t direct_phase_budget,
    PathwayMode mode, TitanExecutor &executor, std::size_t requested_workers,
    std::size_t ligand_atom_count) {
  PathwayDirectStats stats;
  stats.enabled = true;
  stats.effective_eval_budget = direct_phase_budget;
  const auto started_at = std::chrono::steady_clock::now();
  const double local_threshold = pathway_direct_local_threshold_size(direct_domain);

  PathwayDirectCell root;
  std::array<double, 1> root_progress{
      pathway_direct_normalized_to_progress(direct_domain, 0.5)};
  std::array<PathwayCandidateState, 1> root_eval{};
  PullDeBatchEvaluation eval_info =
      evaluator(std::span<const double>(root_progress),
                std::span<PathwayCandidateState>(root_eval));
  root.center_norm = 0.5;
  root.state = root_eval[0];
  std::size_t phase_evals = eval_info.evaluations_used;
  if (eval_info.budget_exhausted) {
    stats.stopped_by_eval_budget = true;
    global_budget_stop = true;
  }

  std::vector<PathwayDirectCell> frontier;
  frontier.push_back(root);
  if (root.state.valid) {
    best_candidate = root.state;
    ranked_candidates.push_back(root.state);
    stats.best_cell_size = pathway_direct_cell_size_metric(root);
  }

  std::size_t next_cell_id = 1;
  const std::size_t max_batch_capacity = std::max<std::size_t>(1, requested_workers * 2);
  std::size_t expansions_since_best = 0;

  struct SplitPlan {
    std::size_t frontier_index{0};
    PathwayDirectCell middle{};
    std::array<PathwayDirectCell, 2> children{};
  };
  std::vector<SplitPlan> plans(max_batch_capacity);
  std::vector<double> child_norm(max_batch_capacity * 2, 0.0);
  std::vector<PathwayCandidateState> child_states(max_batch_capacity * 2);

  while (!frontier.empty()) {
    if (direct_phase_budget != 0 && phase_evals >= direct_phase_budget) {
      stats.stopped_by_eval_budget = true;
      global_budget_stop = true;
      break;
    }

    const std::size_t chosen_preview =
        select_pathway_direct_cell_index(frontier, mode);
    const double preview_size =
        pathway_direct_cell_size_metric(frontier[chosen_preview]);
    if (preview_size <= local_threshold) {
      stats.stopped_by_local_threshold = true;
      break;
    }
    if (expansions_since_best >= 32) {
      stats.stopped_by_stagnation = true;
      break;
    }

    const std::size_t remaining_evals =
        direct_phase_budget > phase_evals ? direct_phase_budget - phase_evals : 0;
    std::size_t batch_cells =
        std::min(max_batch_capacity,
                 std::min(frontier.size(),
                          std::max<std::size_t>(1, remaining_evals / 2)));
    if (batch_cells == 0)
      break;

    const auto selected_indices =
        select_pathway_direct_cell_batch_indices(frontier, batch_cells, mode);
    if (selected_indices.empty())
      break;

    std::size_t active_plans = 0;
    for (const std::size_t frontier_index : selected_indices) {
      PathwayDirectCell chosen = frontier[frontier_index];
      const double split_half = chosen.half_width_norm / 3.0;
      if (split_half < kPathwayDirectMinCellNorm)
        continue;

      SplitPlan plan;
      plan.frontier_index = frontier_index;
      plan.middle = chosen;
      plan.middle.depth = chosen.depth + 1;
      plan.middle.half_width_norm = split_half;
      plan.middle.id = next_cell_id++;
      for (std::size_t child = 0; child < 2; ++child) {
        plan.children[child] = plan.middle;
        plan.children[child].center_norm =
            chosen.center_norm + (child == 0 ? -split_half : split_half);
        plan.children[child].center_norm =
            std::clamp(plan.children[child].center_norm, 0.0, 1.0);
        plan.children[child].id = next_cell_id++;
        child_norm[active_plans * 2 + child] = pathway_direct_normalized_to_progress(
            direct_domain, plan.children[child].center_norm);
      }
      plans[active_plans++] = plan;
    }
    if (active_plans == 0)
      break;

    stats.max_batch_cells = std::max(stats.max_batch_cells, active_plans);
    eval_info = evaluator(
        std::span<const double>(child_norm.data(), active_plans * 2),
        std::span<PathwayCandidateState>(child_states.data(), active_plans * 2));
    phase_evals += eval_info.evaluations_used;
    if (eval_info.budget_exhausted) {
      stats.stopped_by_eval_budget = true;
      global_budget_stop = true;
    }

    for (std::size_t plan_index = 0; plan_index < active_plans; ++plan_index) {
      frontier[plans[plan_index].frontier_index] = plans[plan_index].middle;
      for (std::size_t child_index = 0; child_index < 2; ++child_index) {
        PathwayDirectCell child = plans[plan_index].children[child_index];
        child.state = child_states[plan_index * 2 + child_index];
        frontier.push_back(child);
        if (!child.state.valid)
          continue;
        ranked_candidates.push_back(child.state);
        if (!best_candidate.valid ||
            pathway_candidate_better(child.state, child.id, best_candidate,
                                     best_stage, mode)) {
          best_candidate = child.state;
          best_stage = stats.expansions + plan_index + 1;
          stats.best_expansion = best_stage;
          stats.best_cell_size = pathway_direct_cell_size_metric(child);
          stats.best_half_width_norm = child.half_width_norm;
          expansions_since_best = 0;
        }
      }
    }

    expansions_since_best += active_plans;
    stats.expansions += active_plans;
    ++stats.batches;
    if (stats.stopped_by_eval_budget)
      break;
  }

  stats.cells_remaining = frontier.size();
  stats.evaluated_candidates = ranked_candidates.size();
  stats.phase_evaluations_used = phase_evals;
  if (best_candidate.valid && stats.best_expansion == 0)
    stats.best_expansion = best_stage;
  stats.wall_seconds =
      std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                    started_at)
          .count();
  return stats;
}

void write_pathway_status_file(const info &output_info, const std::string &status,
                               const std::string &reason, double path_distance,
                               bool used_old_vector, const de_arg &DE_arguments) {
  std::fstream status_file(output_info.get_pull_direction_status_filename(),
                           std::ios::out);
  status_file << std::fixed << std::setprecision(4);
  status_file << "STATUS=" << status << "\n";
  status_file << "REASON=" << reason << "\n";
  status_file << "PATH_DISTANCE=" << path_distance << "\n";
  status_file << "USED_OLD_VECTOR=" << (used_old_vector ? 1 : 0) << "\n";
  if (DE_arguments.snapshot_metadata_valid)
    status_file << "SNAPSHOT_TIME_PS=" << DE_arguments.snapshot_time_ps << "\n";
}

} // namespace

int run_pathway_pull_direction_legacy_unused(
    Multi_MOL &rec, Multi_MOL &lg, cache_sigma_epsilon &CachedData,
    de_arg &DE_arguments, arguments &Arguments, info &output_info,
    TitanExecutor &executor, groups_cal &Groups_from_index_file) {
#if 0
  if (rec.m_mol[0].mindist_mol(lg.m_mol[0], Arguments, DE_arguments) == true) {
    TitanUI::print_info("Automatic termination is on!");
    return 1;
  }

  const auto built = build_path(rec, lg, Groups_from_index_file, Arguments);
  if (!built.ok) {
    TitanUI::print_error("Pull-direction pathway: " + built.error);
    return 1;
  }
  if (built.path.pruned_component_nodes > 0) {
    TitanUI::print_warning(
        "Pathway auto-order ignored " +
        std::to_string(built.path.pruned_component_nodes) +
        " residue node(s) outside the component nearest g2.");
  }
  if (built.path.pruned_branch_nodes > 0) {
    TitanUI::print_warning(
        "Pathway auto-order pruned " +
        std::to_string(built.path.pruned_branch_nodes) +
        " branch residue node(s) outside the principal path.");
  }

  const Point3 start_com = mol_com(lg.m_mol[0], lg.parameters);
  const cache_energy energy_ref =
      reference_energy(rec, lg, CachedData, Arguments, executor);
  const std::size_t worker_count = std::max<std::size_t>(1, Arguments.nt);

  std::vector<double> progress(DE_arguments.pop_num, 0.0);
  std::vector<CandidateResult> candidates(DE_arguments.pop_num);
  {
    std::random_device rd;
    std::mt19937 mt(rd());
    for (std::size_t i = 0; i < progress.size(); ++i)
      progress[i] = sample_interval(mt, built.path.valid_intervals);
  }

  executor.parallel_for(
      progress.size(), worker_count, 8,
      [&](std::size_t worker_id, std::size_t begin, std::size_t end) {
        (void)worker_id;
        cache_sigma_epsilon eval_cache = make_eval_cache(CachedData);
        MOL scratch{lg.m_mol[0]};
        for (std::size_t i = begin; i < end; ++i)
          candidates[i] = evaluate_candidate(built.path, progress[i], start_com,
                                             rec, lg, eval_cache, energy_ref,
                                             Arguments, scratch);
      });

  double prev_best = std::numeric_limits<double>::max();
  int stagnation = 0;
  for (std::size_t iter = 0; iter < DE_arguments.ite_num; ++iter) {
    TitanUI::show_progress("Pulling Direction", iter, DE_arguments.ite_num,
                           "(Pathway Evolution)");
    const std::vector<double> old_progress = progress;
    executor.parallel_for(
        progress.size(), worker_count, 8,
        [&](std::size_t worker_id, std::size_t begin, std::size_t end) {
          cache_sigma_epsilon eval_cache = make_eval_cache(CachedData);
          MOL scratch{lg.m_mol[0]};
          std::random_device rd;
          std::mt19937 mt(static_cast<unsigned>(
              rd() ^ (0x9e3779b9U + worker_id * 1315423911U +
                      static_cast<unsigned>(iter * 2654435761U))));
          std::uniform_int_distribution<std::size_t> choose(0, progress.size() - 1);
          std::uniform_real_distribution<double> prob(0.0, 1.0);
          for (std::size_t i = begin; i < end; ++i) {
            std::size_t a = choose(mt), b = choose(mt), c = choose(mt);
            while (a == b || a == c || b == c || a == i || b == i || c == i) {
              a = choose(mt);
              b = choose(mt);
              c = choose(mt);
            }
            double trial = old_progress[i];
            if (DE_arguments.algorithm == 0) {
              if (prob(mt) < DE_arguments.cross_prob || i == choose(mt))
                trial = old_progress[a] +
                        DE_arguments.diff_weight *
                            (old_progress[b] - old_progress[c]);
            } else if (prob(mt) < DE_arguments.cross_prob) {
              trial = old_progress[a] +
                      DE_arguments.F * (old_progress[b] - old_progress[c]);
            } else {
              trial = old_progress[a] +
                      0.5 * (DE_arguments.F + 1.0) *
                          (old_progress[b] + old_progress[c] -
                           2.0 * old_progress[a]);
            }
            trial = std::clamp(trial, 0.0, built.path.total_length);
            trial = clamp_to_intervals(trial, built.path.valid_intervals);
            const auto candidate = evaluate_candidate(
                built.path, trial, start_com, rec, lg, eval_cache, energy_ref,
                Arguments, scratch);
            if (candidate.total_score < candidates[i].total_score) {
              progress[i] = trial;
              candidates[i] = candidate;
            }
          }
        });

    double iter_best = std::numeric_limits<double>::max();
    for (const auto &candidate : candidates)
      iter_best = std::min(iter_best, candidate.total_score);
    if (iter_best < prev_best - 1.0e-8) {
      prev_best = iter_best;
      stagnation = 0;
    } else {
      ++stagnation;
    }
    if (stagnation >= DE_arguments.stagnation_limit)
      break;
  }

  std::size_t best_index = 0;
  for (std::size_t i = 1; i < candidates.size(); ++i) {
    if (candidates[i].total_score < candidates[best_index].total_score)
      best_index = i;
  }
  const Point3 best_com = add_pt(
      start_com,
      {candidates[best_index].vector.Fx, candidates[best_index].vector.Fy,
       candidates[best_index].vector.Fz});
  const double best_angle =
      angle_to_previous_vector(candidates[best_index].vector,
                               DE_arguments.old_vector);

  std::vector<std::pair<double, std::size_t>> ranked(candidates.size());
  for (std::size_t i = 0; i < candidates.size(); ++i)
    ranked[i] = {candidates[i].total_score, i};
  std::sort(ranked.begin(), ranked.end());

  std::fstream de_output(output_info.get_DE_output_filename(), std::ios::out);
  std::fstream de_log(output_info.get_DE_log_filename(),
                      std::ios::out | std::ios::app);
  std::fstream de_report("DE_report.dat", std::ios::out);
  de_report << std::fixed << std::setprecision(4);
  de_report << "TITAN Pathway Pull Report\n\n";
  de_report << "g1 reference: " << Arguments.group_1 << "\n";
  de_report << "g2 pulled: " << Arguments.group_2 << "\n";
  de_report << "pathway group: " << built.path.group_name << " ["
            << built.path.group_index << "]\n";
  de_report << "pathway rescue radius: " << Arguments.de.dis_along_path
            << " A\n";
  de_report << "raw nodes: " << built.path.all_nodes.size() << "\n";
  de_report << "principal nodes: " << built.path.principal_nodes.size() << "\n";
  de_report << "pruned disconnected: " << built.path.pruned_component_nodes
            << "\n";
  de_report << "pruned branches: " << built.path.pruned_branch_nodes << "\n";
  de_report << "used geometric edges: "
            << (built.path.used_geometric_edges ? "YES" : "NO") << "\n";
  if (built.path.pruned_component_nodes > 0)
    de_report << "warning: disconnected pathway nodes were ignored because they "
                 "are outside the start component near g2.\n";
  if (built.path.pruned_branch_nodes > 0)
    de_report << "warning: branch pathway nodes were pruned to keep a single "
                 "principal path.\n";
  de_report << "start residue: " << residue_label(built.path.principal_nodes.front())
            << "\n";
  de_report << "terminal residue: "
            << residue_label(built.path.principal_nodes.back()) << "\n";
  de_report << "arc length: " << built.path.total_length << " A\n";
  de_report << "best score: " << candidates[best_index].total_score << "\n";
  de_report << "best raw score: " << candidates[best_index].raw_score << "\n";
  de_report << "best progress: " << candidates[best_index].progress << " A\n";
  de_report << "pull vector: (" << candidates[best_index].vector.Fx << ", "
            << candidates[best_index].vector.Fy << ", "
            << candidates[best_index].vector.Fz << ")\n";
  de_report << "pull distance: " << candidates[best_index].vector.modul()
            << " A\n";
  de_report << "angle to prev vector: " << best_angle << " deg\n\n";
  de_report << "Top candidates\n";
  de_report << std::setw(8) << "Rank" << std::setw(8) << "Index"
            << std::setw(14) << "Score" << std::setw(14) << "RawScore"
            << std::setw(14) << "Progress" << std::setw(12) << "Pull(A)"
            << std::setw(14) << "Elec" << std::setw(14) << "VdW" << "\n";
  for (std::size_t rank = 0; rank < std::min<std::size_t>(50, ranked.size());
       ++rank) {
    const auto idx = ranked[rank].second;
    de_report << std::setw(8) << rank + 1 << std::setw(8) << idx
              << std::setw(14) << candidates[idx].total_score << std::setw(14)
              << candidates[idx].raw_score << std::setw(14)
              << candidates[idx].progress << std::setw(12)
              << candidates[idx].vector.modul() << std::setw(14)
              << candidates[idx].energy.elec_ener << std::setw(14)
              << candidates[idx].energy.vdw_ener << "\n";
  }

  de_output << std::fixed << std::setprecision(4);
  de_output << candidates[best_index].vector.Fx << "  "
            << candidates[best_index].vector.Fy << "  "
            << candidates[best_index].vector.Fz << "  Angle: " << best_angle
            << " -- Modul: " << candidates[best_index].vector.modul()
            << " Angstrom\n";
  de_output << "# Pull g2 toward the best reachable point on the auto-ordered "
               "pathway of g1\n";
  de_output << "# Start residue: " << residue_label(built.path.principal_nodes.front())
            << "\n";
  de_output << "# Terminal residue: "
            << residue_label(built.path.principal_nodes.back()) << "\n";
  de_output << "# Lig COM: (" << start_com[0] << ", " << start_com[1] << ", "
            << start_com[2] << ") -> (" << best_com[0] << ", " << best_com[1]
            << ", " << best_com[2] << ")\n";
  de_log << std::fixed << std::setprecision(4);
  de_log << candidates[best_index].vector.Fx << "  "
         << candidates[best_index].vector.Fy << "  "
         << candidates[best_index].vector.Fz << "  Angle: " << best_angle
         << " -- Modul: " << candidates[best_index].vector.modul()
         << " Angstrom  Score: " << candidates[best_index].total_score
         << "  Progress: " << candidates[best_index].progress
         << "\n--------------------------\n";

  std::fstream best_conf("best_conf.pdb", std::ios::out);
  MOL best_candidate{lg.m_mol[0]};
  best_candidate.translate(candidates[best_index].vector.Fx,
                           candidates[best_index].vector.Fy,
                           candidates[best_index].vector.Fz);
  best_conf.precision(3);
  for (std::size_t i = 0; i < best_candidate.size(); ++i) {
    best_conf << "ATOM" << std::setw(7) << i + 1 << std::setw(5)
              << lg.info[i].get_atom() << std::setw(4)
              << lg.info[i].get_original_res_name() << " A" << std::setw(4) << 1
              << std::setw(12) << std::fixed << best_candidate.pos.x[i]
              << std::setw(8) << best_candidate.pos.y[i] << std::setw(8)
              << best_candidate.pos.z[i] << std::setw(7) << 1.00
              << std::setw(6) << 0.00 << '\n';
  }

  std::fstream vector_pdb(output_info.get_DE_pdb_vector_filename(),
                          std::ios::out | std::ios::app);
  vector_pdb.precision(3);
  vector_pdb << "ATOM" << std::setw(7) << 1 << std::setw(5) << "CA"
             << std::setw(4) << "ALA"
             << " A" << std::setw(4) << 1 << std::setw(12) << std::fixed
             << start_com[0] << std::setw(8) << start_com[1] << std::setw(8)
             << start_com[2] << " " << std::setw(7) << 1.00 << " "
             << std::setw(6) << 0.00 << '\n';
  vector_pdb << "ATOM" << std::setw(7) << 2 << std::setw(5) << "CB"
             << std::setw(4) << "ALA"
             << " A" << std::setw(4) << 1 << std::setw(12) << std::fixed
             << best_com[0] << std::setw(8) << best_com[1] << std::setw(8)
             << best_com[2] << " " << std::setw(7) << 1.00 << " "
             << std::setw(6) << 0.00 << '\n';
  vector_pdb << "TER\n";
  return 0;
#endif
  (void)rec;
  (void)lg;
  (void)CachedData;
  (void)DE_arguments;
  (void)Arguments;
  (void)output_info;
  (void)executor;
  (void)Groups_from_index_file;
  return 1;
}

int run_pathway_pull_direction(Multi_MOL &rec, Multi_MOL &lg,
                               cache_sigma_epsilon &CachedData,
                               de_arg &DE_arguments, arguments &Arguments,
                               info &output_info, TitanExecutor &executor,
                               groups_cal &Groups_from_index_file) {
  if (rec.m_mol[0].mindist_mol(lg.m_mol[0], Arguments, DE_arguments) == true) {
    TitanUI::print_info("Automatic termination is on!");
    return 1;
  }

  const auto built = build_path(rec, lg, Groups_from_index_file, Arguments);
  if (!built.ok) {
    TitanUI::print_error("Pull-direction pathway: " + built.error);
    return 1;
  }
  if (built.path.pruned_component_nodes > 0) {
    TitanUI::print_warning(
        "Pathway auto-order ignored " +
        std::to_string(built.path.pruned_component_nodes) +
        " residue node(s) outside the component nearest g2.");
  }
  if (built.path.pruned_branch_nodes > 0) {
    TitanUI::print_warning(
        "Pathway auto-order pruned " +
        std::to_string(built.path.pruned_branch_nodes) +
        " branch residue node(s) outside the principal path.");
  }

  const Point3 start_com = mol_com(lg.m_mol[0], lg.parameters);
  const PathwayState pathway_state =
      prepare_pathway_state(built.path, start_com, Arguments);
  const bool off_path =
      pathway_state.current_projection.distance > Arguments.de.dis_along_path;

  std::string stop_reason;
  if (DE_arguments.flag_old_vec && off_path &&
      old_vector_points_outward(built.path, start_com,
                                pathway_state.current_projection,
                                DE_arguments.old_vector)) {
    stop_reason = "OLD_VECTOR_POINTS_AWAY_FROM_PATHWAY";
  } else if (pathway_state.forward_intervals.empty()) {
    stop_reason = "NO_REACHABLE_FORWARD_TARGET";
  }

  const auto started_at = std::chrono::steady_clock::now();
  const cache_energy energy_ref =
      reference_energy(rec, lg, CachedData, Arguments, executor);
  const std::size_t requested_workers = std::max<std::size_t>(1, Arguments.nt);
  const std::size_t worker_slots =
      std::max<std::size_t>(1, std::min(executor.capacity(), requested_workers));
  std::vector<PathwayWorkerContext> workers(worker_slots);
  for (auto &worker : workers) {
    worker.scratch = lg.m_mol[0];
    worker.eval_cache = make_eval_cache(CachedData);
  }

  PathwayEvalBudget eval_budget;
  eval_budget.block_size =
      pathway_eval_budget_block_size(DE_arguments.pop_num, worker_slots);

  const PullOptimizerKind optimizer_kind =
      classify_pull_optimizer(DE_arguments.pull_opt);
  const bool use_direct_family =
      optimizer_is_direct_family(DE_arguments.pull_opt);
  const bool use_variant_de =
      (optimizer_kind == PullOptimizerKind::DeShade ||
       optimizer_kind == PullOptimizerKind::DeLShade ||
       optimizer_kind == PullOptimizerKind::DeJso ||
       optimizer_kind == PullOptimizerKind::DeJsoPattern ||
       optimizer_kind == PullOptimizerKind::DeNlShadeLbc ||
       optimizer_kind == PullOptimizerKind::DeNlShadeLbcPattern);
  const bool use_variant_pattern =
      (optimizer_kind == PullOptimizerKind::DeJsoPattern ||
       optimizer_kind == PullOptimizerKind::DeNlShadeLbcPattern);
  const bool use_de_pattern = (optimizer_kind == PullOptimizerKind::DePattern);
  const bool use_direct_pattern =
      (optimizer_kind == PullOptimizerKind::DirectPattern);
  const std::size_t variant_total_budget =
      use_variant_de
          ? (DE_arguments.pull_max_evals != 0
                 ? DE_arguments.pull_max_evals
                 : DE_arguments.pop_num * (DE_arguments.ite_num + 1))
          : DE_arguments.pull_max_evals;
  const std::size_t direct_total_budget =
      use_direct_family ? effective_pathway_direct_eval_budget(DE_arguments) : 0;
  std::size_t variant_phase_budget = variant_total_budget;
  std::size_t direct_phase_budget = direct_total_budget;
  std::size_t pattern_reserved_budget = 0;
  if (use_variant_pattern && variant_total_budget != 0) {
    pattern_reserved_budget = clamp_pattern_tail_budget(variant_total_budget);
    if (variant_total_budget > pattern_reserved_budget)
      variant_phase_budget = variant_total_budget - pattern_reserved_budget;
  }
  if (use_direct_pattern && direct_total_budget != 0) {
    pattern_reserved_budget = clamp_pattern_tail_budget(direct_total_budget);
    if (direct_total_budget > pattern_reserved_budget)
      direct_phase_budget = direct_total_budget - pattern_reserved_budget;
  }
  eval_budget.max_evals = use_variant_de
                              ? variant_phase_budget
                              : use_direct_family ? direct_phase_budget
                                                  : DE_arguments.pull_max_evals;

  auto sanitize_progress = [&](double progress) {
    progress = std::clamp(progress, pathway_state.forward_intervals.front().begin,
                          pathway_state.forward_intervals.back().end);
    return clamp_to_intervals(progress, pathway_state.forward_intervals);
  };
  const PathwayProgressDomain domain{
      pathway_state.forward_intervals.empty()
          ? 0.0
          : pathway_state.forward_intervals.front().begin,
      pathway_state.forward_intervals.empty()
          ? 0.0
          : pathway_state.forward_intervals.back().end,
      [&](std::mt19937 &mt) {
        return sample_interval(mt, pathway_state.forward_intervals);
      },
      sanitize_progress};
  const PathwayDirectDomain direct_domain =
      build_pathway_direct_domain(pathway_state.forward_intervals);

  const auto batch_evaluator =
      [&](std::span<const double> proposals,
          std::span<PathwayCandidateState> evaluated_candidates) {
        std::vector<double> sanitized(proposals.size(), 0.0);
        for (std::size_t i = 0; i < proposals.size(); ++i)
          sanitized[i] = domain.sanitize(proposals[i]);
        return evaluate_candidate_batch(
            std::span<const double>(sanitized), evaluated_candidates, built.path,
            pathway_state, start_com, rec, lg, workers, energy_ref, Arguments,
            DE_arguments, eval_budget, executor, requested_workers);
      };

  std::vector<PathwayCandidateState> report_candidates;
  PathwayCandidateState best_candidate;
  PullDeVariantStats de_variant_stats;
  PathwayLegacyDeStats legacy_stats;
  PathwayDirectStats direct_stats;
  PathwayPatternStats pattern_stats;
  std::size_t best_stage = 0;
  bool stopped_by_eval_budget = false;

  if (stop_reason.empty()) {
    bool ok = false;
    if (use_variant_de) {
      PullDeVariantOptions variant_options;
      variant_options.variant = de_variant_from_optimizer(optimizer_kind);
      variant_options.seed = DE_arguments.pull_seed;
      variant_options.initial_population = DE_arguments.pop_num;
      variant_options.minimum_population = DE_arguments.min_pop_num;
      variant_options.max_generations = DE_arguments.ite_num;
      variant_options.planned_evaluations = variant_total_budget;
      variant_options.memory_size = DE_arguments.memory_size;
      variant_options.p_best_rate = DE_arguments.p_best_rate;
      variant_options.archive_scale = DE_arguments.archive_scale;
      variant_options.stagnation_limit = DE_arguments.stagnation_limit;
      ok = run_pathway_de_variant(variant_options, domain, pathway_state.mode,
                                  batch_evaluator, report_candidates,
                                  best_candidate, de_variant_stats);
      if (use_variant_pattern && ok) {
        eval_budget.max_evals = variant_total_budget;
        pattern_stats = refine_pathway_variant_multi_seed_pattern(
            best_candidate, report_candidates, domain, batch_evaluator,
            de_variant_stats, built.path.total_length, pathway_state.mode);
        if (pattern_stats.improved)
          report_candidates.push_back(best_candidate);
      }
      stopped_by_eval_budget = de_variant_stats.stopped_by_eval_budget ||
                               pattern_stats.stopped_by_eval_budget;
      best_stage = de_variant_stats.best_generation;
    } else if (use_direct_family) {
      direct_stats = run_pathway_direct_search(
          best_candidate, report_candidates, best_stage, stopped_by_eval_budget,
          direct_domain, batch_evaluator, direct_phase_budget, pathway_state.mode,
          executor, requested_workers, lg.m_mol[0].size());
      direct_stats.used_fallback_budget = (DE_arguments.pull_max_evals == 0);
      ok = best_candidate.valid;
      if (use_direct_pattern && ok) {
        eval_budget.max_evals = direct_total_budget;
        const double initial_step = std::max(
            0.1, 0.25 * std::max(0.1, domain.max_progress - domain.min_progress));
        pattern_stats = refine_pathway_candidate_pattern(
            best_candidate, domain, batch_evaluator, initial_step,
            pathway_state.mode);
        pattern_stats.seeds_attempted = 1;
        pattern_stats.distinct_seed_basins = 1;
        if (pattern_stats.improved)
          report_candidates.push_back(best_candidate);
      }
      stopped_by_eval_budget = direct_stats.stopped_by_eval_budget ||
                               pattern_stats.stopped_by_eval_budget;
    } else {
      ok = run_pathway_legacy_de(DE_arguments, domain, pathway_state.mode,
                                 batch_evaluator, report_candidates,
                                 best_candidate, legacy_stats);
      if (use_de_pattern && ok) {
        const double initial_step = std::max(
            0.1, 0.25 * std::max(0.1, domain.max_progress - domain.min_progress));
        pattern_stats = refine_pathway_candidate_pattern(
            best_candidate, domain, batch_evaluator, initial_step,
            pathway_state.mode);
        pattern_stats.seeds_attempted = 1;
        pattern_stats.distinct_seed_basins = 1;
        if (pattern_stats.improved)
          report_candidates.push_back(best_candidate);
      }
      stopped_by_eval_budget = legacy_stats.stopped_by_eval_budget ||
                               pattern_stats.stopped_by_eval_budget;
      best_stage = legacy_stats.best_generation;
    }
    if (!ok || !best_candidate.valid)
      stop_reason = "NO_VALID_FORWARD_DIRECTION";
  }

  const double wall_seconds =
      std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                    started_at)
          .count();
  const std::size_t total_evaluations = total_pathway_evaluations_used(workers);
  const double evals_per_second =
      wall_seconds > 0.0 ? static_cast<double>(total_evaluations) / wall_seconds
                         : 0.0;

  if (!stop_reason.empty()) {
    write_pathway_status_file(output_info, "STOP", stop_reason,
                              pathway_state.current_projection.distance,
                              DE_arguments.flag_old_vec, DE_arguments);
    std::fstream de_log(output_info.get_DE_log_filename(),
                        std::ios::out | std::ios::app);
    std::fstream de_report("DE_report.dat", std::ios::out);
    de_report << std::fixed << std::setprecision(4);
    de_report << "TITAN Pathway Pull Report\n\n";
    de_report << "status: STOP\n";
    de_report << "reason: " << stop_reason << "\n";
    de_report << "optimizer: " << DE_arguments.pull_opt << "\n";
    de_report << "pathway group: " << built.path.group_name << " ["
              << built.path.group_index << "]\n";
    de_report << "pathway radius: " << Arguments.de.dis_along_path << " A\n";
    de_report << "current mode: " << pathway_mode_label(pathway_state.mode)
              << "\n";
    de_report << "current progress: " << pathway_state.current_projection.progress
              << " A\n";
    de_report << "current path distance: "
              << pathway_state.current_projection.distance << " A\n";
    de_report << "terminal residue: "
              << residue_label(built.path.principal_nodes.back()) << "\n";
    de_report << "total evaluations: " << total_evaluations << "\n";
    de_report << "wall time: " << wall_seconds << " s\n";
    de_report << "eval/s: " << evals_per_second << "\n";
    de_log << "STOP  Reason: " << stop_reason << "  PathDist: "
           << pathway_state.current_projection.distance
           << "  Progress: " << pathway_state.current_projection.progress
           << "\n--------------------------\n";
    TitanUI::print_info("Pathway pull-direction stopped: " + stop_reason);
    return 0;
  }

  const Point3 best_com =
      add_pt(start_com, {best_candidate.vector.Fx, best_candidate.vector.Fy,
                         best_candidate.vector.Fz});
  const double best_angle = best_candidate.angle_deg;
  const auto report_rows =
      build_grouped_report_rows(report_candidates, pathway_state.mode);

  write_pathway_status_file(output_info, "CONTINUE", "OK",
                            pathway_state.current_projection.distance,
                            DE_arguments.flag_old_vec, DE_arguments);

  std::fstream de_output(output_info.get_DE_output_filename(), std::ios::out);
  std::fstream de_log(output_info.get_DE_log_filename(),
                      std::ios::out | std::ios::app);
  std::fstream de_report("DE_report.dat", std::ios::out);
  std::fstream best_conf("best_conf.pdb", std::ios::out);
  std::fstream vector_pdb(output_info.get_DE_pdb_vector_filename(),
                          std::ios::out | std::ios::app);

  de_report << std::fixed << std::setprecision(4);
  de_report << "==============================================================="
               "=================\n";
  de_report << "  TITAN Pathway Pull Report\n";
  de_report << "==============================================================="
               "=================\n\n";
  de_report << "  Parameters\n";
  de_report << "  ----------\n";
  de_report << "  Optimizer          : " << DE_arguments.pull_opt << "\n";
  de_report << "  Optimizer family   : "
            << (use_direct_family ? "Direct" : "DE") << "\n";
  de_report << "  Variant            : "
            << (use_variant_de
                    ? std::string(pull_de_variant_name(de_variant_stats.variant))
                    : DE_arguments.pull_opt)
            << "\n";
  de_report << "  Pathway group      : " << built.path.group_name << " ["
            << built.path.group_index << "]\n";
  de_report << "  Pathway radius     : " << Arguments.de.dis_along_path << " A\n";
  de_report << "  Current mode       : "
            << pathway_mode_label(pathway_state.mode) << "\n";
  de_report << "  Current progress   : "
            << pathway_state.current_projection.progress << " A\n";
  de_report << "  Current path dist  : "
            << pathway_state.current_projection.distance << " A\n";
  de_report << "  Arc length         : " << built.path.total_length << " A\n";
  de_report << "  Start residue      : "
            << residue_label(built.path.principal_nodes.front()) << "\n";
  de_report << "  Terminal residue   : "
            << residue_label(built.path.principal_nodes.back()) << "\n";
  de_report << "  Scan radius        : " << DE_arguments.scan_radius << " A\n";
  de_report << "  Min scan radius    : " << DE_arguments.min_scan_radius
            << " A\n";
  de_report << "  Restrain angle     : " << DE_arguments.restrain_angle
            << " deg\n";
  de_report << "  Seed               : " << DE_arguments.pull_seed << "\n\n";

  de_report << "  Benchmark\n";
  de_report << "  ---------\n";
  de_report << "  Total evaluations  : " << total_evaluations << "\n";
  de_report << "  Wall time          : " << wall_seconds << " s\n";
  de_report << "  Eval/s             : " << evals_per_second << "\n";
  de_report << "  Best stage         : " << best_stage << "\n";
  de_report << "  Eval budget stop   : "
            << (stopped_by_eval_budget ? "YES" : "NO") << "\n\n";

  if (use_direct_family) {
    de_report << "  DIRECT Stats\n";
    de_report << "  ------------\n";
    de_report << "  Effective budget   : " << direct_stats.effective_eval_budget
              << "\n";
    de_report << "  Used fallback budg.: "
              << (direct_stats.used_fallback_budget ? "YES" : "NO") << "\n";
    de_report << "  Expansions         : " << direct_stats.expansions << "\n";
    de_report << "  Batches            : " << direct_stats.batches << "\n";
    de_report << "  Max batch cells    : " << direct_stats.max_batch_cells
              << "\n";
    de_report << "  Cells remaining    : " << direct_stats.cells_remaining
              << "\n";
    de_report << "  Best expansion     : " << direct_stats.best_expansion
              << "\n";
    de_report << "  Best cell size     : " << direct_stats.best_cell_size
              << "\n";
    de_report << "  Local threshold    : "
              << pathway_direct_local_threshold_size(direct_domain) << "\n";
    de_report << "  Local stop         : "
              << (direct_stats.stopped_by_local_threshold ? "YES" : "NO")
              << "\n";
    de_report << "  Stagnation stop    : "
              << (direct_stats.stopped_by_stagnation ? "YES" : "NO")
              << "\n\n";
  } else if (use_variant_de) {
    de_report << "  DE Variant Stats\n";
    de_report << "  ----------------\n";
    de_report << "  Reduction mode   : "
              << (de_variant_stats.nonlinear_population_reduction ? "nonlinear"
                                                                  : "linear/fixed")
              << "\n";
    de_report << "  Initial pop       : " << de_variant_stats.initial_population
              << "\n";
    de_report << "  Final pop         : " << de_variant_stats.final_population
              << "\n";
    de_report << "  Reductions        : "
              << de_variant_stats.population_reductions << "\n";
    de_report << "  Generations done  : "
              << de_variant_stats.generations_completed << "\n";
    de_report << "  Best generation   : " << de_variant_stats.best_generation
              << "\n";
    de_report << "  Final p-best      : " << de_variant_stats.final_p_best_rate
              << "\n";
    if (de_variant_stats.selective_pressure > 0.0) {
      de_report << "  Selective press.  : "
                << de_variant_stats.selective_pressure << "\n";
    }
    if (de_variant_stats.nonlinear_population_reduction) {
      de_report << "  Bias power F      : " << de_variant_stats.bias_power_f
                << "\n";
      de_report << "  Bias power CR     : " << de_variant_stats.bias_power_cr
                << "\n";
    }
    de_report << "  Stagnation stop   : "
              << (de_variant_stats.stopped_by_stagnation ? "YES" : "NO")
              << "\n\n";
  }

  if (pattern_stats.enabled) {
    const std::string title = use_direct_pattern
                                  ? std::string("DIRECT Pattern Refinement")
                              : use_variant_de
                                  ? variant_pattern_section_title(
                                        de_variant_stats.variant)
                                  : std::string("Pattern Refinement");
    de_report << "  " << title << "\n";
    de_report << "  " << std::string(title.size(), '-') << "\n";
    if (use_variant_pattern) {
      de_report << "  Variant budget     : " << variant_phase_budget << "\n";
      de_report << "  Pattern reserve    : " << pattern_reserved_budget << "\n";
    }
    if (use_direct_pattern) {
      de_report << "  Variant budget     : " << direct_phase_budget << "\n";
      de_report << "  Pattern reserve    : " << pattern_reserved_budget << "\n";
    }
    de_report << "  Raw score before   : " << pattern_stats.initial_score << "\n";
    de_report << "  Raw score after    : " << pattern_stats.final_score << "\n";
    de_report << "  Improved           : "
              << (pattern_stats.improved ? "YES" : "NO") << "\n";
    de_report << "  Sweeps             : " << pattern_stats.sweeps << "\n";
    de_report << "  Accepted moves     : " << pattern_stats.accepted_moves
              << "\n";
    de_report << "  Extra evals        : " << pattern_stats.evaluations_used
              << "\n";
    de_report << "  Seeds tried        : " << pattern_stats.seeds_attempted
              << "\n";
    de_report << "  Distinct basins    : "
              << pattern_stats.distinct_seed_basins << "\n\n";
  }

  de_report << "  Best Candidate\n";
  de_report << "  --------------\n";
  de_report << "  Best raw score     : " << best_candidate.raw_score << "\n";
  de_report << "  Best energy key    : " << best_candidate.energy_key << "\n";
  de_report << "  Best progress      : " << best_candidate.progress << " A\n";
  de_report << "  delta_progress     : " << best_candidate.delta_progress
            << " A\n";
  de_report << "  Pull vector        : (" << best_candidate.vector.Fx << ", "
            << best_candidate.vector.Fy << ", " << best_candidate.vector.Fz
            << ")\n";
  de_report << "  Pull distance      : " << best_candidate.pull_distance
            << " A\n";
  de_report << "  Angle to old       : " << best_angle << " deg\n\n";

  de_report << "  Ranked Candidate Table\n";
  de_report << "  ----------------------\n";
  de_report << std::setw(8) << "Rank" << std::setw(8) << "Index"
            << std::setw(8) << "Count" << std::setw(14) << "RawScore"
            << std::setw(14) << "EKey" << std::setw(14) << "Progress"
            << std::setw(14) << "dProg" << std::setw(12) << "Pull(A)"
            << std::setw(12) << "Angle" << std::setw(18) << "Elec"
            << std::setw(18) << "VdW" << "\n";
  for (std::size_t rank = 0;
       rank < std::min<std::size_t>(50, report_rows.size());
       ++rank) {
    const auto &row = report_rows[rank];
    const auto &candidate = row.representative;
    de_report << std::setw(8) << rank + 1 << std::setw(8)
              << row.representative_index << std::setw(8) << row.count
              << std::scientific << std::setw(14) << candidate.raw_score
              << std::setw(14) << candidate.energy_key << std::fixed
              << std::setw(14) << candidate.progress << std::setw(14)
              << candidate.delta_progress << std::setw(12)
              << candidate.pull_distance << std::setw(12)
              << candidate.angle_deg
              << std::scientific << std::setw(18)
              << candidate.energy.elec_ener << std::setw(18)
              << candidate.energy.vdw_ener << std::fixed << "\n";
  }

  de_output << std::fixed << std::setprecision(4);
  de_output << best_candidate.vector.Fx << "  " << best_candidate.vector.Fy
            << "  " << best_candidate.vector.Fz << "  Angle: " << best_angle
            << " -- Modul: " << best_candidate.vector.modul() << " Angstrom\n";
  de_output << "# Pull g2 toward the best reachable point on the forward tail "
               "of the auto-ordered pathway of g1\n";
  de_output << "# Current progress: " << pathway_state.current_projection.progress
            << " A\n";
  de_output << "# Terminal residue: "
            << residue_label(built.path.principal_nodes.back()) << "\n";
  de_output << "# Lig COM: (" << start_com[0] << ", " << start_com[1] << ", "
            << start_com[2] << ") -> (" << best_com[0] << ", " << best_com[1]
            << ", " << best_com[2] << ")\n";

  de_log << std::fixed << std::setprecision(4);
  de_log << best_candidate.vector.Fx << "  " << best_candidate.vector.Fy
         << "  " << best_candidate.vector.Fz << "  Angle: " << best_angle
         << " -- Modul: " << best_candidate.vector.modul()
         << " Angstrom  RawScore: " << best_candidate.raw_score
         << "  EnergyKey: " << best_candidate.energy_key
         << "  Progress: " << best_candidate.progress
         << "\n--------------------------\n";

  MOL best_conf_mol{lg.m_mol[0]};
  best_conf_mol.translate(best_candidate.vector.Fx, best_candidate.vector.Fy,
                          best_candidate.vector.Fz);
  best_conf.precision(3);
  for (std::size_t i = 0; i < best_conf_mol.size(); ++i) {
    best_conf << "ATOM" << std::setw(7) << i + 1 << std::setw(5)
              << lg.info[i].get_atom() << std::setw(4)
              << lg.info[i].get_original_res_name() << " A" << std::setw(4) << 1
              << std::setw(12) << std::fixed << best_conf_mol.pos.x[i]
              << std::setw(8) << best_conf_mol.pos.y[i] << std::setw(8)
              << best_conf_mol.pos.z[i] << std::setw(7) << 1.00
              << std::setw(6) << 0.00 << '\n';
  }

  vector_pdb.precision(3);
  vector_pdb << "ATOM" << std::setw(7) << 1 << std::setw(5) << "CA"
             << std::setw(4) << "ALA"
             << " A" << std::setw(4) << 1 << std::setw(12) << std::fixed
             << start_com[0] << std::setw(8) << start_com[1] << std::setw(8)
             << start_com[2] << " " << std::setw(7) << 1.00 << " "
             << std::setw(6) << 0.00 << '\n';
  vector_pdb << "ATOM" << std::setw(7) << 2 << std::setw(5) << "CB"
             << std::setw(4) << "ALA"
             << " A" << std::setw(4) << 1 << std::setw(12) << std::fixed
             << best_com[0] << std::setw(8) << best_com[1] << std::setw(8)
             << best_com[2] << " " << std::setw(7) << 1.00 << " "
             << std::setw(6) << 0.00 << '\n';
  vector_pdb << "TER\n";
  return 0;
}
