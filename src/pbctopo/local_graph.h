#ifndef TITAN_PBCTOPO_LOCAL_GRAPH_H
#define TITAN_PBCTOPO_LOCAL_GRAPH_H

#include "types.h"
#include "../box_type.h"

#include <cstddef>
#include <cstdint>
#include <string>
#include <unordered_set>
#include <vector>

class MOL;
class Multi_MOL;
class mol_info;

namespace titan_pbctopo {

struct PbctopoRange {
  std::size_t begin = 0;
  std::size_t end = 0;
  std::size_t source_index = 0;
  bool topo_backed = false;
};

struct PbctopoSystemStaticData {
  std::size_t atom_count = 0;
  std::size_t heavy_count = 0;
  std::vector<char> heavy_mask;
  std::vector<PbctopoRange> ranges;
  std::vector<PbctopoContinuityPair> continuity_pairs_local;
  std::vector<std::size_t> continuity_anchor_local;
};

struct PbctopoTrajectoryContext {
  std::size_t total_atoms_hint = 0;
  std::size_t heavy_atoms_hint = 0;
  std::size_t local_edges_hint = 0;
  std::size_t projection_edges_hint = 0;
  std::size_t component_hint = 0;
  std::size_t plane_hint = 0;
  std::size_t virtual_state_hint = 512;
  std::size_t coarse_seed_hint = 32;
  std::size_t coarse_candidate_hint = 16;
  std::size_t continuity_pair_hint = 0;
  PbctopoSystemStaticData rec_static;
  PbctopoSystemStaticData lig_static;
  std::vector<PbctopoContinuityPair> continuity_pairs;
  std::vector<std::size_t> continuity_anchor_indices;
};

std::vector<PbctopoRange> build_ranges_for_system(const Multi_MOL &system,
                                                  const MOL &mol);

void build_heavy_atom_subset(const std::vector<PbctopoAtom> &atoms,
                             std::vector<PbctopoAtom> &heavy_atoms);

PbctopoSystemStaticData
build_pbctopo_system_static_data(const Multi_MOL &system, const MOL &mol);

PbctopoTrajectoryContext
build_pbctopo_trajectory_context(const Multi_MOL &mol1_system,
                                 const Multi_MOL &mol2_system,
                                 std::size_t frame_count);

void build_projection_edge_subset(
    const std::vector<PbctopoAtom> &base_atoms,
    const std::vector<PbctopoAtom> &projection_atoms,
    const std::vector<PbctopoLocalEdge> &local_edges,
    std::vector<PbctopoLocalEdge> &projection_edges,
    std::unordered_set<std::uint64_t> &projection_edge_set,
    std::vector<std::size_t> &projection_index);

void build_projection_continuity_subset(
    const std::vector<PbctopoAtom> &base_atoms,
    const std::vector<PbctopoAtom> &projection_atoms,
    const std::vector<PbctopoContinuityPair> &continuity_pairs,
    std::vector<PbctopoContinuityPair> &projection_continuity_pairs,
    std::vector<std::size_t> &projection_index);

void build_projection_anchor_subset(
    const std::vector<PbctopoAtom> &base_atoms,
    const std::vector<std::size_t> &anchor_indices,
    std::vector<std::size_t> &projection_anchor_indices,
    const std::vector<std::size_t> &projection_index);

void append_atoms_for_system(const MOL &mol,
                             const std::vector<mol_info> *info,
                             const PbctopoSystemStaticData *static_data,
                             int owner,
                             std::vector<PbctopoAtom> &atoms);

void build_local_neighborhood_graph(
    const std::vector<PbctopoAtom> &atoms, box Box,
    std::vector<PbctopoLocalEdge> &edges,
    std::unordered_set<std::uint64_t> &edge_set, double bond_cutoff);

void append_topology_bond_edges_for_system(
    const Multi_MOL &system, int owner, const std::vector<PbctopoAtom> &atoms,
    box Box, std::vector<PbctopoLocalEdge> &edges,
    std::unordered_set<std::uint64_t> &edge_set, PbctopoTimingStats &timing);

std::string pbctopo_topology_source_label(const PbctopoTimingStats &timing);

} // namespace titan_pbctopo

#endif
