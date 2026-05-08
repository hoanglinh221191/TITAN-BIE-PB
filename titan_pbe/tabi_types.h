#ifndef TITAN_PBE_TABI_TYPES_H
#define TITAN_PBE_TABI_TYPES_H

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

class TitanExecutor;

namespace TitanPBE::Tabi {

struct Atom {
  double x{0.0};
  double y{0.0};
  double z{0.0};
  double q{0.0};
  double radius{0.0};
};

struct Panel {
  double x{0.0};
  double y{0.0};
  double z{0.0};
  double nx{0.0};
  double ny{0.0};
  double nz{0.0};
  double area{0.0};
  double diagonal{1.0};
  double solid_angle{0.0};
  std::size_t atom_index{0};
};

struct Vertex {
  double x{0.0};
  double y{0.0};
  double z{0.0};
  double nx{0.0};
  double ny{0.0};
  double nz{0.0};
  std::size_t atom_index{0};
};

struct Triangle {
  std::size_t v0{0};
  std::size_t v1{0};
  std::size_t v2{0};
  std::size_t atom_index{0};
  double cx{0.0};
  double cy{0.0};
  double cz{0.0};
  double nx{0.0};
  double ny{0.0};
  double nz{0.0};
  double area{0.0};
};

struct SurfaceMesh {
  std::vector<Vertex> vertices;
  std::vector<Triangle> triangles;
  std::vector<Panel> panels;
  std::size_t exposed_candidates{0};
  std::size_t buried_candidates{0};
  std::size_t solid_angle_fallbacks{0};
  std::size_t nonmanifold_edges{0};
  double avg_valence{0.0};
  double diag_min{1.0};
  double diag_max{1.0};
  std::size_t cell_count{0};
  std::size_t max_neighbor_count{0};
  double mean_neighbor_count{0.0};
  double neighbor_build_ms{0.0};
  double mesh_build_ms{0.0};
  double solid_angle_ms{0.0};
};

struct Options {
  double eps_in{1.0};
  double eps_out{80.0};
  double salt_m{0.15};
  double temperature_k{298.15};
  int max_iter{10000};
  double residual_tol{1.0e-6};
  std::string sum_mode{"auto"};
  double probe_radius{0.0};
  std::string surface_mode{"vdw"};
  std::string surface_geometry{"triangle"};
  std::string diagonal_mode{"solid-angle"};
  int taylor_order{3};
  std::string formulation{"jcp2013"};
  std::string energy_projection{"boundary"};
};

struct Result {
  std::string status{"ok"};
  std::string reaction_energy_source{"tabi_boundary_integral"};
  std::string convergence_reason{"not_started"};
  std::string residual_history{};
  std::size_t panel_count{0};
  std::size_t memory_bytes{0};
  std::size_t peak_memory_bytes{0};
  double surface_area_a2{0.0};
  double surface_build_ms{0.0};
  double surface_neighbor_ms{0.0};
  double surface_mesh_ms{0.0};
  double solid_angle_ms{0.0};
  std::size_t exposed_points{0};
  std::size_t buried_points{0};
  double total_charge{0.0};
  int iterations{0};
  double residual{0.0};
  bool converged{false};
  std::string sum_mode_effective{"direct"};
  double tree_build_ms{0.0};
  std::size_t tree_nodes{0};
  std::size_t tree_leaf_count{0};
  std::size_t near_interactions{0};
  std::size_t far_interactions{0};
  std::string matvec_mode{"direct"};
  double rhs_ms{0.0};
  double gmres_ms{0.0};
  double matvec_ms_total{0.0};
  double matvec_ms_avg{0.0};
  std::size_t matvec_calls{0};
  double orthogonalization_ms{0.0};
  double preconditioner_ms{0.0};
  std::string surface_geometry{"triangle"};
  std::size_t mesh_vertices{0};
  std::size_t mesh_triangles{0};
  std::size_t boundary_nodes{0};
  double avg_valence{0.0};
  std::size_t nonmanifold_edges{0};
  std::size_t solid_angle_fallbacks{0};
  double diag_min{1.0};
  double diag_max{1.0};
  int taylor_order{3};
  double interaction_build_ms{0.0};
  double tree_eval_ms{0.0};
  double tree_moment_ms{0.0};
  double tree_traversal_ms{0.0};
  double near_eval_ms{0.0};
  double far_eval_ms{0.0};
  double energy_kcal_mol{0.0};
  double energy_kj_mol{0.0};
  double energy_kt{0.0};
  std::string formulation{"jcp2013"};
  std::string energy_projection{"boundary"};
  double rhs_scale{1.0};
  double dielectric_contrast{0.0};
  double operator_contrast{0.0};
  std::size_t unknown_count{0};
  std::string bie_system{"off"};
  double lhs_diag_phi{0.0};
  double lhs_diag_phin{0.0};
  std::string energy_equation{"off"};
  std::string kernel_convention{"off"};
  double source_norm_s1{0.0};
  double source_norm_s2{0.0};
  double phi_norm{0.0};
  double phin_norm{0.0};
  double kappa_inv_a{0.0};
  double coulomb_factor_kcal_mol{0.0};
  double boundary_energy_kcal_mol{0.0};
  double reference_energy_kcal_mol{0.0};
  double selected_energy_scale{1.0};
  bool fallback_used{false};
};

struct TreecodeStats {
  double tree_build_ms{0.0};
  std::size_t tree_nodes{0};
  std::size_t tree_leaf_count{0};
  std::size_t near_interactions{0};
  std::size_t far_interactions{0};
  std::string matvec_mode{"direct"};
  int taylor_order{3};
  double interaction_build_ms{0.0};
  double tree_eval_ms{0.0};
  double tree_moment_ms{0.0};
  double tree_traversal_ms{0.0};
  double near_eval_ms{0.0};
  double far_eval_ms{0.0};
  std::size_t matvec_calls{0};
};

struct TreeNode {
  double cx{0.0};
  double cy{0.0};
  double cz{0.0};
  double radius{0.0};
  double area_sum{0.0};
  std::size_t first_child{0};
  std::size_t child_count{0};
  std::size_t first_panel{0};
  std::size_t panel_count{0};
};

struct TreeWorkspace {
  std::vector<std::size_t> order;
  std::vector<std::uint64_t> morton;
  std::vector<TreeNode> nodes;
  std::size_t root_index{0};
  std::size_t leaf_count{0};
  double build_ms{0.0};
  int taylor_order{3};
};

std::vector<Panel> build_surface(const std::vector<Atom> &atoms,
                                 const Options &options,
                                 Result *result = nullptr,
                                 TitanExecutor *executor = nullptr);

void apply_direct_operator(const std::vector<Panel> &panels,
                           const std::vector<double> &density,
                           double contrast, double kappa,
                           TitanExecutor &executor,
                           std::vector<double> &out,
                           TreecodeStats *stats = nullptr);

TreeWorkspace build_tree_workspace(const std::vector<Panel> &panels,
                                   TitanExecutor &executor);

void apply_treecode_operator(const std::vector<Panel> &panels,
                             const std::vector<double> &density,
                             double contrast, double kappa,
                             const TreeWorkspace &tree,
                             TitanExecutor &executor,
                             std::vector<double> &out,
                             TreecodeStats *stats = nullptr);

Result solve_linear_pb(const std::vector<Atom> &atoms, const Options &options,
                       TitanExecutor &executor);

} // namespace TitanPBE::Tabi

#endif // TITAN_PBE_TABI_TYPES_H
