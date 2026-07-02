#include "local_graph.h"

#include "../BoxMatrix.h"
#include "../Multi_MOL.h"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <functional>
#include <limits>
#include <string>
#include <unordered_map>

namespace titan_pbctopo {
namespace {

double wrap_unit(double value) {
  value -= std::floor(value);
  if (value >= 1.0)
    value -= 1.0;
  if (value < 0.0)
    value += 1.0;
  return value;
}

std::uint64_t atom_pair_key(std::size_t a, std::size_t b) {
  if (a > b)
    std::swap(a, b);
  return (static_cast<std::uint64_t>(a) << 32) | static_cast<std::uint64_t>(b);
}

double cart_dist_sq_from_scaled_delta(box Box, double dsx, double dsy,
                                      double dsz) {
  double v1[3]{}, v2[3]{}, v3[3]{};
  Box.get_v1(v1);
  Box.get_v2(v2);
  Box.get_v3(v3);
  const double dx = v1[0] * dsx + v2[0] * dsy + v3[0] * dsz;
  const double dy = v1[1] * dsx + v2[1] * dsy + v3[1] * dsz;
  const double dz = v1[2] * dsx + v2[2] * dsy + v3[2] * dsz;
  return dx * dx + dy * dy + dz * dz;
}

void scaled_to_cartesian(box Box, double sx, double sy, double sz, double &x,
                         double &y, double &z) {
  double v1[3]{}, v2[3]{}, v3[3]{};
  Box.get_v1(v1);
  Box.get_v2(v2);
  Box.get_v3(v3);
  x = sx * v1[0] + sy * v2[0] + sz * v3[0];
  y = sx * v1[1] + sy * v2[1] + sz * v3[1];
  z = sx * v1[2] + sy * v2[2] + sz * v3[2];
}

std::vector<PbctopoRange>
build_ranges_from_index(const std::vector<::index> &Index,
                        std::size_t atom_count) {
  std::vector<PbctopoRange> ranges;
  ranges.reserve(Index.size());
  for (std::size_t k = 0; k + 1 < Index.size(); ++k) {
    const std::size_t begin = Index[k].get_res_pos_e(0) - 1;
    const std::size_t end = Index[k + 1].get_res_pos_e(0) - 1;
    if (begin >= atom_count || end > atom_count || begin >= end)
      continue;
    ranges.push_back({begin, end, ranges.size(), true});
  }
  return ranges;
}

std::vector<PbctopoRange>
build_ranges_from_info(const std::vector<mol_info> &info,
                       std::size_t atom_count) {
  std::vector<PbctopoRange> ranges;
  if (info.size() != atom_count || atom_count == 0)
    return ranges;

  std::size_t begin = 0;
  std::size_t source_index = 0;
  for (std::size_t i = 1; i < atom_count; ++i) {
    if (info[i].get_chainnum() != info[i - 1].get_chainnum() ||
        info[i].get_chain() != info[i - 1].get_chain()) {
      if (begin < i)
        ranges.push_back({begin, i, source_index++, true});
      begin = i;
    }
  }
  if (begin < atom_count)
    ranges.push_back({begin, atom_count, source_index, true});
  return ranges;
}

std::vector<PbctopoRange>
build_ranges_from_topology(const MolTopology &topology,
                           std::size_t atom_count) {
  std::vector<PbctopoRange> ranges;
  const auto &molecule_ranges = topology.molecule_ranges();
  ranges.reserve(molecule_ranges.size());
  for (const MolTopologyMoleculeRange &range : molecule_ranges) {
    if (range.begin >= range.end || range.end > atom_count)
      continue;
    ranges.push_back({range.begin, range.end, ranges.size(), true});
  }
  return ranges;
}

char first_alpha_upper(const std::string &text) {
  for (unsigned char ch : text) {
    if (std::isalpha(ch))
      return static_cast<char>(std::toupper(ch));
  }
  return '\0';
}

bool detect_heavy_atom(const std::vector<mol_info> *info,
                       std::size_t local_atom_index) {
  if (info == nullptr || local_atom_index >= info->size())
    return true;
  const auto &atom_info = (*info)[local_atom_index];
  const std::string atom_name = atom_info.get_atom();
  const std::string &atom_type = atom_info.get_atom_type();
  const char atom_symbol = first_alpha_upper(atom_name);
  const char type_symbol = first_alpha_upper(atom_type);
  return atom_symbol != 'H' && type_symbol != 'H';
}

struct PbctopoResidueChainKey {
  int owner = 0;
  std::uint64_t chain_key = 0;

  bool operator==(const PbctopoResidueChainKey &other) const {
    return owner == other.owner && chain_key == other.chain_key;
  }
};

struct PbctopoResidueChainKeyHash {
  std::size_t operator()(const PbctopoResidueChainKey &key) const {
    const std::uint64_t owner_bits =
        static_cast<std::uint64_t>(static_cast<std::uint32_t>(key.owner));
    const std::uint64_t mixed =
        key.chain_key ^ (owner_bits * 0x9e3779b97f4a7c15ull);
    return static_cast<std::size_t>(mixed ^ (mixed >> 32));
  }
};

std::vector<int>
build_topology_residue_sequence_indices(const std::vector<PbctopoAtom> &atoms) {
  std::vector<int> residue_sequence_indices(
      atoms.size(), std::numeric_limits<int>::min());
  struct ChainState {
    int last_residue_id = std::numeric_limits<int>::min();
    int current_rank = -1;
    bool initialized = false;
  };
  std::unordered_map<PbctopoResidueChainKey, ChainState,
                     PbctopoResidueChainKeyHash>
      chain_states;
  chain_states.reserve(std::max<std::size_t>(8, atoms.size() / 64));

  for (std::size_t atom_idx = 0; atom_idx < atoms.size(); ++atom_idx) {
    const auto &atom = atoms[atom_idx];
    if (!atom.has_topology_key)
      continue;
    const PbctopoResidueChainKey key{atom.owner, atom.chain_key};
    auto &state = chain_states[key];
    if (!state.initialized) {
      state.initialized = true;
      state.last_residue_id = atom.residue_id;
      state.current_rank = 0;
    } else if (atom.residue_id != state.last_residue_id) {
      state.last_residue_id = atom.residue_id;
      ++state.current_rank;
    }
    residue_sequence_indices[atom_idx] = state.current_rank;
  }

  return residue_sequence_indices;
}

std::vector<PbctopoContinuityPair> build_local_continuity_pairs(
    const Multi_MOL &system, const PbctopoSystemStaticData &static_data) {
  std::vector<PbctopoContinuityPair> continuity_pairs;
  continuity_pairs.reserve(static_data.atom_count);
  if (system.info.size() == static_data.atom_count &&
      static_data.atom_count > 0) {
    auto atom_priority = [&](std::size_t atom_idx) {
      const std::string atom_name = system.info[atom_idx].get_atom();
      if (atom_name == "CA")
        return 0;
      if (atom_name == "C4'" || atom_name == "C4*")
        return 1;
      if (atom_name == "BB")
        return 2;
      if (atom_idx < static_data.heavy_mask.size() &&
          static_data.heavy_mask[atom_idx] != 0)
        return 3;
      return 4;
    };

    auto append_range_pairs = [&](std::size_t begin, std::size_t end) {
      if (begin >= end || begin >= static_data.atom_count)
        return;
      end = std::min(end, static_data.atom_count);

      std::size_t prev_anchor = std::numeric_limits<std::size_t>::max();
      std::string prev_chain;
      std::size_t prev_chainnum = std::numeric_limits<std::size_t>::max();
      std::size_t i = begin;
      while (i < end) {
        const auto &start_info = system.info[i];
        const std::size_t chainnum = start_info.get_chainnum();
        const std::string chain = start_info.get_chain();
        const std::size_t resnr = start_info.get_resnr();
        std::size_t j = i + 1;
        while (j < end && system.info[j].get_chainnum() == chainnum &&
               system.info[j].get_chain() == chain &&
               system.info[j].get_resnr() == resnr) {
          ++j;
        }

        std::size_t anchor = std::numeric_limits<std::size_t>::max();
        int best_priority = std::numeric_limits<int>::max();
        for (std::size_t atom_idx = i; atom_idx < j; ++atom_idx) {
          const int priority = atom_priority(atom_idx);
          if (priority < best_priority) {
            best_priority = priority;
            anchor = atom_idx;
          }
        }

        if (anchor != std::numeric_limits<std::size_t>::max()) {
          if (prev_anchor != std::numeric_limits<std::size_t>::max() &&
              chainnum == prev_chainnum && chain == prev_chain) {
            continuity_pairs.push_back({prev_anchor, anchor});
          }
          prev_anchor = anchor;
          prev_chainnum = chainnum;
          prev_chain = chain;
        }
        i = j;
      }
    };

    if (static_data.ranges.empty()) {
      append_range_pairs(0, static_data.atom_count);
    } else {
      for (const auto &range : static_data.ranges)
        append_range_pairs(range.begin, range.end);
    }

    if (!continuity_pairs.empty())
      return continuity_pairs;
  }

  for (const auto &range : static_data.ranges) {
    std::size_t prev = std::numeric_limits<std::size_t>::max();
    for (std::size_t atom_idx = range.begin; atom_idx < range.end; ++atom_idx) {
      if (atom_idx < static_data.heavy_mask.size() &&
          static_data.heavy_mask[atom_idx] == 0) {
        continue;
      }
      if (prev != std::numeric_limits<std::size_t>::max())
        continuity_pairs.push_back({prev, atom_idx});
      prev = atom_idx;
    }
  }
  return continuity_pairs;
}

void append_global_continuity_pairs(
    const std::vector<PbctopoContinuityPair> &local_pairs, std::size_t offset,
    std::vector<PbctopoContinuityPair> &global_pairs) {
  global_pairs.reserve(global_pairs.size() + local_pairs.size());
  for (const auto &pair : local_pairs)
    global_pairs.push_back({pair.a + offset, pair.b + offset});
}

void build_anchor_index_list(
    std::size_t atom_count,
    const std::vector<PbctopoContinuityPair> &continuity_pairs,
    std::vector<std::size_t> &anchor_indices) {
  anchor_indices.clear();
  if (atom_count == 0 || continuity_pairs.empty())
    return;
  std::vector<char> seen(atom_count, 0);
  for (const auto &pair : continuity_pairs) {
    if (pair.a < atom_count && seen[pair.a] == 0) {
      seen[pair.a] = 1;
      anchor_indices.push_back(pair.a);
    }
    if (pair.b < atom_count && seen[pair.b] == 0) {
      seen[pair.b] = 1;
      anchor_indices.push_back(pair.b);
    }
  }
}

void append_global_anchor_indices(
    const std::vector<std::size_t> &local_indices, std::size_t offset,
    std::vector<std::size_t> &global_indices) {
  global_indices.reserve(global_indices.size() + local_indices.size());
  for (std::size_t local_idx : local_indices)
    global_indices.push_back(local_idx + offset);
}

} // namespace

std::vector<PbctopoRange> build_ranges_for_system(const Multi_MOL &system,
                                                  const MOL &mol) {
  std::vector<PbctopoRange> ranges =
      build_ranges_from_index(system.Index, mol.pos.x.size());
  if (!ranges.empty())
    return ranges;

  ranges = build_ranges_from_info(system.info, mol.pos.x.size());
  if (!ranges.empty())
    return ranges;

  if (!mol.pos.x.empty())
    ranges.push_back({0, mol.pos.x.size(), 0, false});
  return ranges;
}

void build_heavy_atom_subset(const std::vector<PbctopoAtom> &atoms,
                             std::vector<PbctopoAtom> &heavy_atoms) {
  heavy_atoms.clear();
  heavy_atoms.reserve(atoms.size());
  for (const auto &atom : atoms) {
    if (atom.is_heavy)
      heavy_atoms.push_back(atom);
  }
  if (heavy_atoms.empty())
    heavy_atoms = atoms;
}

PbctopoSystemStaticData
build_pbctopo_system_static_data(const Multi_MOL &system, const MOL &mol) {
  PbctopoSystemStaticData data;
  data.atom_count = mol.pos.x.size();
  data.ranges = build_ranges_for_system(system, mol);
  data.heavy_mask.assign(data.atom_count, 1);
  data.heavy_count = 0;
  for (std::size_t atom_idx = 0; atom_idx < data.atom_count; ++atom_idx) {
    const bool is_heavy = detect_heavy_atom(&system.info, atom_idx);
    data.heavy_mask[atom_idx] = is_heavy ? 1 : 0;
    data.heavy_count += is_heavy ? 1u : 0u;
  }
  if (data.heavy_count == 0 && data.atom_count > 0) {
    std::fill(data.heavy_mask.begin(), data.heavy_mask.end(), 1);
    data.heavy_count = data.atom_count;
  }
  PbctopoSystemStaticData continuity_data = data;
  continuity_data.ranges =
      build_ranges_from_topology(system.topology, mol.pos.x.size());
  if (continuity_data.ranges.empty())
    continuity_data.ranges =
        build_ranges_from_info(system.info, mol.pos.x.size());
  if (continuity_data.ranges.empty())
    continuity_data.ranges = data.ranges;
  data.continuity_pairs_local =
      build_local_continuity_pairs(system, continuity_data);
  build_anchor_index_list(data.atom_count, data.continuity_pairs_local,
                          data.continuity_anchor_local);
  return data;
}

PbctopoTrajectoryContext
build_pbctopo_trajectory_context(const Multi_MOL &mol1_system,
                                 const Multi_MOL &mol2_system,
                                 std::size_t frame_count) {
  PbctopoTrajectoryContext context;
  if (frame_count == 0 || mol1_system.m_mol.empty())
    return context;

  const MOL &rec = mol1_system.m_mol.front();
  context.rec_static = build_pbctopo_system_static_data(mol1_system, rec);
  std::size_t lig_atoms = 0;
  std::size_t lig_heavy = 0;
  if (!mol2_system.m_mol.empty() && !mol2_system.m_mol.front().pos.x.empty()) {
    const MOL &lig = mol2_system.m_mol.front();
    context.lig_static = build_pbctopo_system_static_data(mol2_system, lig);
    lig_atoms = context.lig_static.atom_count;
    lig_heavy = context.lig_static.heavy_count;
  }

  context.total_atoms_hint =
      std::max<std::size_t>(1, context.rec_static.atom_count + lig_atoms);
  context.heavy_atoms_hint =
      std::max<std::size_t>(1, context.rec_static.heavy_count + lig_heavy);
  context.local_edges_hint =
      std::max<std::size_t>(1024, context.total_atoms_hint * 8);
  context.projection_edges_hint =
      std::max<std::size_t>(1024, context.heavy_atoms_hint * 8);
  append_global_continuity_pairs(context.rec_static.continuity_pairs_local, 0,
                                 context.continuity_pairs);
  append_global_continuity_pairs(context.lig_static.continuity_pairs_local,
                                 context.rec_static.atom_count,
                                 context.continuity_pairs);
  append_global_anchor_indices(context.rec_static.continuity_anchor_local, 0,
                               context.continuity_anchor_indices);
  append_global_anchor_indices(context.lig_static.continuity_anchor_local,
                               context.rec_static.atom_count,
                               context.continuity_anchor_indices);
  context.continuity_pair_hint =
      std::max<std::size_t>(64, context.continuity_pairs.size() + 1);
  context.component_hint =
      std::max<std::size_t>(8, mol1_system.Index.size() + mol2_system.Index.size());
  context.plane_hint = std::max<std::size_t>(32, context.heavy_atoms_hint * 2 + 4);
  context.virtual_state_hint = 512;
  context.coarse_seed_hint = 32;
  context.coarse_candidate_hint = 16;
  return context;
}

void build_projection_edge_subset(
    const std::vector<PbctopoAtom> &base_atoms,
    const std::vector<PbctopoAtom> &projection_atoms,
    const std::vector<PbctopoLocalEdge> &local_edges,
    std::vector<PbctopoLocalEdge> &projection_edges,
    std::unordered_set<std::uint64_t> &projection_edge_set,
    std::vector<std::size_t> &projection_index) {
  projection_edges.clear();
  projection_edge_set.clear();
  if (projection_atoms.size() == base_atoms.size()) {
    projection_edges = local_edges;
    projection_edge_set.reserve(local_edges.size() * 2 + 1);
    for (const auto &edge : projection_edges)
      projection_edge_set.insert(atom_pair_key(edge.a, edge.b));
    return;
  }

  const std::size_t invalid_index = std::numeric_limits<std::size_t>::max();
  projection_index.assign(base_atoms.size(), invalid_index);
  for (std::size_t idx = 0; idx < projection_atoms.size(); ++idx)
    projection_index[projection_atoms[idx].global_atom_id] = idx;

  projection_edges.reserve(local_edges.size());
  projection_edge_set.reserve(local_edges.size() + 1);
  for (const auto &edge : local_edges) {
    const std::size_t mapped_a =
        projection_index[base_atoms[edge.a].global_atom_id];
    const std::size_t mapped_b =
        projection_index[base_atoms[edge.b].global_atom_id];
    if (mapped_a == invalid_index || mapped_b == invalid_index)
      continue;
    if (!projection_edge_set.insert(atom_pair_key(mapped_a, mapped_b)).second)
      continue;
    projection_edges.push_back({mapped_a, mapped_b, edge.ref_d2,
                                edge.inv_ref_d2, edge.break_cutoff_d2});
  }
}

void build_projection_continuity_subset(
    const std::vector<PbctopoAtom> &base_atoms,
    const std::vector<PbctopoAtom> &projection_atoms,
    const std::vector<PbctopoContinuityPair> &continuity_pairs,
    std::vector<PbctopoContinuityPair> &projection_continuity_pairs,
    std::vector<std::size_t> &projection_index) {
  projection_continuity_pairs.clear();
  if (projection_atoms.size() == base_atoms.size()) {
    projection_continuity_pairs = continuity_pairs;
    return;
  }
  const std::size_t invalid_index = std::numeric_limits<std::size_t>::max();
  if (projection_index.size() != base_atoms.size())
    projection_index.assign(base_atoms.size(), invalid_index);
  for (const auto &pair : continuity_pairs) {
    const std::size_t mapped_a =
        projection_index[base_atoms[pair.a].global_atom_id];
    const std::size_t mapped_b =
        projection_index[base_atoms[pair.b].global_atom_id];
    if (mapped_a == invalid_index || mapped_b == invalid_index)
      continue;
    projection_continuity_pairs.push_back({mapped_a, mapped_b});
  }
}

void build_projection_anchor_subset(
    const std::vector<PbctopoAtom> &base_atoms,
    const std::vector<std::size_t> &anchor_indices,
    std::vector<std::size_t> &projection_anchor_indices,
    const std::vector<std::size_t> &projection_index) {
  projection_anchor_indices.clear();
  if (anchor_indices.empty())
    return;
  if (projection_index.size() != base_atoms.size()) {
    projection_anchor_indices = anchor_indices;
    return;
  }
  const std::size_t invalid_index = std::numeric_limits<std::size_t>::max();
  projection_anchor_indices.reserve(anchor_indices.size());
  for (std::size_t anchor_idx : anchor_indices) {
    if (anchor_idx >= base_atoms.size())
      continue;
    const std::size_t mapped =
        projection_index[base_atoms[anchor_idx].global_atom_id];
    if (mapped == invalid_index)
      continue;
    projection_anchor_indices.push_back(mapped);
  }
}

void append_atoms_for_system(const MOL &mol,
                             const std::vector<mol_info> *info,
                             const PbctopoSystemStaticData *static_data,
                             int owner,
                             std::vector<PbctopoAtom> &atoms) {
  for (std::size_t local = 0; local < mol.pos.x.size(); ++local) {
    double r[3]{mol.pos.x[local], mol.pos.y[local], mol.pos.z[local]};
    double s[3]{};
    mol.Box.to_scaled(r, s);
    s[0] = wrap_unit(s[0]);
    s[1] = wrap_unit(s[1]);
    s[2] = wrap_unit(s[2]);
    scaled_to_cartesian(mol.Box, s[0], s[1], s[2], r[0], r[1], r[2]);

    PbctopoAtom atom;
    atom.global_atom_id = atoms.size();
    atom.fragment_id = 0;
    atom.owner = owner;
    atom.local_atom_index = local;
    if (info != nullptr && local < info->size()) {
      const auto &atom_info = (*info)[local];
      const std::uint64_t chainnum =
          static_cast<std::uint64_t>(atom_info.get_chainnum() & 0xffffffffu);
      const std::uint64_t chainhash =
          static_cast<std::uint64_t>(std::hash<std::string>{}(atom_info.get_chain()));
      atom.chain_key = (chainnum << 32) ^ (chainhash & 0xffffffffu);
      atom.residue_id = static_cast<int>(atom_info.get_resnr());
      atom.has_topology_key = true;
    }
    atom.is_heavy =
        (static_data != nullptr && local < static_data->heavy_mask.size())
            ? (static_data->heavy_mask[local] != 0)
            : detect_heavy_atom(info, local);
    atom.x = r[0];
    atom.y = r[1];
    atom.z = r[2];
    atom.sx = s[0];
    atom.sy = s[1];
    atom.sz = s[2];
    atoms.push_back(atom);
  }
}

void build_local_neighborhood_graph(
    const std::vector<PbctopoAtom> &atoms, box Box,
    std::vector<PbctopoLocalEdge> &edges,
    std::unordered_set<std::uint64_t> &edge_set, double bond_cutoff) {
  edges.clear();
  edge_set.clear();
  if (atoms.size() < 2)
    return;
  const std::vector<int> residue_sequence_indices =
      build_topology_residue_sequence_indices(atoms);

  const double kBondCutoff = bond_cutoff;
  const double kBondCutoff2 = bond_cutoff * bond_cutoff;
  edges.reserve(atoms.size() * 2);
  edge_set.reserve(atoms.size() * 2);
  const BoxMatrix box_matrix(Box);
  const double sx_bound = std::min(
      0.5, kBondCutoff *
                   std::sqrt(box_matrix.inv[0][0] * box_matrix.inv[0][0] +
                             box_matrix.inv[0][1] * box_matrix.inv[0][1] +
                             box_matrix.inv[0][2] * box_matrix.inv[0][2]) +
               1.0e-9);
  const double sy_bound = std::min(
      0.5, kBondCutoff *
                   std::sqrt(box_matrix.inv[1][0] * box_matrix.inv[1][0] +
                             box_matrix.inv[1][1] * box_matrix.inv[1][1] +
                             box_matrix.inv[1][2] * box_matrix.inv[1][2]) +
               1.0e-9);
  const double sz_bound = std::min(
      0.5, kBondCutoff *
                   std::sqrt(box_matrix.inv[2][0] * box_matrix.inv[2][0] +
                             box_matrix.inv[2][1] * box_matrix.inv[2][1] +
                             box_matrix.inv[2][2] * box_matrix.inv[2][2]) +
               1.0e-9);
  const int nx =
      std::max(1, static_cast<int>(std::floor(1.0 / std::max(sx_bound, 1.0e-6))));
  const int ny =
      std::max(1, static_cast<int>(std::floor(1.0 / std::max(sy_bound, 1.0e-6))));
  const int nz =
      std::max(1, static_cast<int>(std::floor(1.0 / std::max(sz_bound, 1.0e-6))));
  const double cell_sx = 1.0 / static_cast<double>(nx);
  const double cell_sy = 1.0 / static_cast<double>(ny);
  const double cell_sz = 1.0 / static_cast<double>(nz);
  auto wrap_index = [](int value, int mod) {
    int wrapped = value % mod;
    if (wrapped < 0)
      wrapped += mod;
    return wrapped;
  };

  std::unordered_map<ClashCellKey, std::vector<std::size_t>, ClashCellKeyHash>
      bond_grid;
  bond_grid.reserve(std::max<std::size_t>(64, atoms.size() / 2));
  for (std::size_t i = 0; i < atoms.size(); ++i) {
    const auto &ai = atoms[i];
    const int cx = std::min(nx - 1, static_cast<int>(std::floor(ai.sx / cell_sx)));
    const int cy = std::min(ny - 1, static_cast<int>(std::floor(ai.sy / cell_sy)));
    const int cz = std::min(nz - 1, static_cast<int>(std::floor(ai.sz / cell_sz)));
    for (int dx = -1; dx <= 1; ++dx) {
      for (int dy = -1; dy <= 1; ++dy) {
        for (int dz = -1; dz <= 1; ++dz) {
          const ClashCellKey key{wrap_index(cx + dx, nx),
                                 wrap_index(cy + dy, ny),
                                 wrap_index(cz + dz, nz)};
          const auto it = bond_grid.find(key);
          if (it == bond_grid.end())
            continue;
          for (std::size_t j : it->second) {
            const auto &aj = atoms[j];
            if (ai.owner != aj.owner)
              continue;
            if (ai.has_topology_key && aj.has_topology_key) {
              if (ai.chain_key != aj.chain_key)
                continue;
              const int ai_residue_rank = residue_sequence_indices[i];
              const int aj_residue_rank = residue_sequence_indices[j];
              if (ai_residue_rank != std::numeric_limits<int>::min() &&
                  aj_residue_rank != std::numeric_limits<int>::min()) {
                if (std::abs(ai_residue_rank - aj_residue_rank) > 1)
                  continue;
              } else if (std::abs(ai.residue_id - aj.residue_id) > 1) {
                continue;
              }
            }
            double dsx = aj.sx - ai.sx;
            double dsy = aj.sy - ai.sy;
            double dsz = aj.sz - ai.sz;
            dsx -= std::round(dsx);
            dsy -= std::round(dsy);
            dsz -= std::round(dsz);
            const double d2 = cart_dist_sq_from_scaled_delta(Box, dsx, dsy, dsz);
            if (d2 > kBondCutoff2)
              continue;
            const std::uint64_t pair_key = atom_pair_key(i, j);
            if (!edge_set.insert(pair_key).second)
              continue;
            const double ref_d2 = std::max(d2, 1.0e-6);
            edges.push_back(
                {j, i, d2, 1.0 / ref_d2, std::max(9.0, ref_d2 * 4.0)});
          }
        }
      }
    }
    bond_grid[{cx, cy, cz}].push_back(i);
  }
}

void append_topology_bond_edges_for_system(
    const Multi_MOL &system, int owner, const std::vector<PbctopoAtom> &atoms,
    box Box, std::vector<PbctopoLocalEdge> &edges,
    std::unordered_set<std::uint64_t> &edge_set, PbctopoTimingStats &timing) {
  if (!system.topology.has_bonds())
    return;

  const auto &topology_bonds = system.topology.bonds();
  timing.pbctopo_topology_bonds_available += topology_bonds.size();
  if (system.info.empty())
    return;

  std::vector<std::size_t> local_to_atom(
      system.info.size(), std::numeric_limits<std::size_t>::max());
  for (std::size_t atom_idx = 0; atom_idx < atoms.size(); ++atom_idx) {
    const auto &atom = atoms[atom_idx];
    if (atom.owner != owner)
      continue;
    if (atom.local_atom_index >= local_to_atom.size())
      continue;
    local_to_atom[atom.local_atom_index] = atom_idx;
  }

  for (const MolTopologyBond &bond : topology_bonds) {
    if (bond.a >= local_to_atom.size() || bond.b >= local_to_atom.size()) {
      ++timing.pbctopo_topology_bond_edges_skipped;
      continue;
    }
    const std::size_t a = local_to_atom[bond.a];
    const std::size_t b = local_to_atom[bond.b];
    if (a == std::numeric_limits<std::size_t>::max() ||
        b == std::numeric_limits<std::size_t>::max() || a == b) {
      ++timing.pbctopo_topology_bond_edges_skipped;
      continue;
    }

    const std::uint64_t pair_key = atom_pair_key(a, b);
    if (!edge_set.insert(pair_key).second) {
      ++timing.pbctopo_topology_bond_edges_skipped;
      continue;
    }

    double dsx = atoms[b].sx - atoms[a].sx;
    double dsy = atoms[b].sy - atoms[a].sy;
    double dsz = atoms[b].sz - atoms[a].sz;
    dsx -= std::round(dsx);
    dsy -= std::round(dsy);
    dsz -= std::round(dsz);
    const double d2 = cart_dist_sq_from_scaled_delta(Box, dsx, dsy, dsz);
    const double ref_d2 = std::max(d2, 1.0e-6);
    edges.push_back({a, b, d2, 1.0 / ref_d2, std::max(9.0, ref_d2 * 4.0)});
    ++timing.pbctopo_topology_bond_edges_added;
  }
}

std::string pbctopo_topology_source_label(const PbctopoTimingStats &timing) {
  return timing.pbctopo_topology_bonds_available > 0
             ? "gromacs_topology+geometry_cutoff"
             : "geometry_cutoff";
}

} // namespace titan_pbctopo
