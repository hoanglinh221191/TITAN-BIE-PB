#ifndef TITAN_PBE_H
#define TITAN_PBE_H

#include "MOL.h"
#include "MOL_para.h"
#include "arguments.h"

#include <cstddef>
#include <cstdint>
#include <fstream>
#include <map>
#include <string>
#include <vector>

class TitanExecutor;
class EnergyDecompositionProvider;

namespace TitanPBE {

struct FrameResult {
  std::size_t trajectory_id{0};
  std::size_t frame_index{0};
  std::size_t global_index{0};
  double time_ps{0.0};
  std::size_t atom_count{0};
  double total_charge{0.0};
  std::size_t nx{0};
  std::size_t ny{0};
  std::size_t nz{0};
  std::size_t grid_points{0};
  std::size_t memory_bytes{0};
  double grid_origin_x_A{0.0};
  double grid_origin_y_A{0.0};
  double grid_origin_z_A{0.0};
  double grid_center_x_A{0.0};
  double grid_center_y_A{0.0};
  double grid_center_z_A{0.0};
  double grid_spacing_A{0.0};
  std::size_t solute_points{0};
  std::size_t solvent_points{0};
  bool focus_enabled{false};
  int focus_levels{1};
  double focus_initial_scale{0.0};
  double focus_initial_margin{0.0};
  std::string solver_mode{"sor"};
  int mg_levels{1};
  int mg_cycles{0};
  int mg_pre{0};
  int mg_post{0};
  int mg_coarse_iter{0};
  int solvent_iterations{0};
  int vacuum_iterations{0};
  int solvent_iterations_final{0};
  int vacuum_iterations_final{0};
  double solvent_rms{0.0};
  double vacuum_rms{0.0};
  double solvent_max_dphi{0.0};
  double vacuum_max_dphi{0.0};
  bool solvent_converged{false};
  bool vacuum_converged{false};
  std::string surface_effective{"vdw"};
  std::string ion_access{"dielectric"};
  double ion_radius_A{2.0};
  std::size_t ion_accessible_points{0};
  std::size_t stern_excluded_points{0};
  bool face_refine{false};
  std::string pbe_model{"linear"};
  std::string pbe_model_effective{"linear"};
  std::string linear_seed_solver{"none"};
  std::string pbe_biomodel{"pbe"};
  std::string dielectric_model{"two-dielectric"};
  std::string nl_solver{"picard"};
  std::string nl_inner_solver{"sor"};
  std::string nl_mg_mode{"off"};
  double gaussian_exp{2.0};
  double gaussian_sigma{1.0};
  double density_cutoff{0.759};
  double surface_cutoff{0.0};
  double gap_diel{22.0};
  std::string radius_convention{"native"};
  std::string focus_parent_grid{};
  std::string focus_child_grid{};
  std::string probe_status{"off"};
  int newton_inner_iterations{0};
  double newton_step_norm{0.0};
  double newton_residual_before{0.0};
  double newton_residual_after{0.0};
  std::string newton_line_search_status{"off"};
  double newton_correction_residual{0.0};
  double newton_correction_max_residual{0.0};
  double newton_accepted_lambda{0.0};
  int newton_rejected_steps{0};
  double newton_jac_diag_min{0.0};
  double newton_jac_diag_max{0.0};
  std::size_t newton_screen_points{0};
  std::size_t newton_active_nonlinear_points{0};
  std::string newton_stop_reason{"off"};
  int newton_inner_levels{1};
  int newton_inner_cycles{0};
  double newton_inner_residual_before{0.0};
  double newton_inner_residual_after{0.0};
  std::string newton_inner_effective{"off"};
  double gaussian_density_sum{0.0};
  std::size_t gaussian_cutoff_points{0};
  double gaussian_eps_min{0.0};
  double gaussian_eps_max{0.0};
  double gaussian_eps_mean{0.0};
  double gaussian_charge_norm_error{0.0};
  std::string boundary_mode{"coulombic"};
  std::string charge_assignment{"trilinear"};
  double rhs_scale{0.0};
  double charge_assigned{0.0};
  double charge_moment_error{0.0};
  double boundary_rms{0.0};
  double boundary_max{0.0};
  double linear_residual{0.0};
  int nl_iterations{0};
  int nl_blocks{0};
  double nl_relax{0.5};
  double nl_max_dphi{0.0};
  double nl_residual{0.0};
  std::size_t nl_clipped_points{0};
  int nl_continuation_steps{0};
  int nl_backtracks{0};
  std::string nl_status{"linear"};
  double energy_kcal_mol{0.0};
  double energy_kj_mol{0.0};
  double energy_kt{0.0};
  double solvent_phase_kcal_mol{0.0};
  double solvent_phase_kj_mol{0.0};
  double solvent_phase_kt{0.0};
  double vacuum_phase_kcal_mol{0.0};
  double vacuum_phase_kj_mol{0.0};
  double vacuum_phase_kt{0.0};
  double atom_corrected_kcal_mol{0.0};
  double atom_corrected_kj_mol{0.0};
  double atom_corrected_kt{0.0};
  double grid_corrected_kcal_mol{0.0};
  double grid_corrected_kj_mol{0.0};
  double grid_corrected_kt{0.0};
  double iscm_water_phase_kcal_mol{0.0};
  double iscm_water_phase_kj_mol{0.0};
  double iscm_water_phase_kt{0.0};
  double iscm_vacuum_phase_kcal_mol{0.0};
  double iscm_vacuum_phase_kj_mol{0.0};
  double iscm_vacuum_phase_kt{0.0};
  double iscm_corrected_kcal_mol{0.0};
  double iscm_corrected_kj_mol{0.0};
  double iscm_corrected_kt{0.0};
  std::size_t iscm_boundary_points{0};
  std::size_t iscm_charged_boundary_points{0};
  double iscm_induced_charge_sum{0.0};
  double iscm_induced_charge_l1{0.0};
  double iscm_induced_charge_l2{0.0};
  double iscm_induced_charge_max{0.0};
  std::string iscm_status{"off"};
  std::size_t delphi_face_crossings{0};
  std::size_t delphi_eps_face_matches{0};
  std::size_t delphi_boundary_grid_points{0};
  std::size_t delphi_charged_boundary_points{0};
  std::string delphi_iscm_status{"off"};
  std::string vdwms_enabled{"off"};
  std::size_t vdwms_boundary_initial{0};
  std::size_t vdwms_boundary_final{0};
  std::size_t vdwms_added{0};
  std::size_t vdwms_removed{0};
  std::size_t vdwms_exposed_points{0};
  double vdwms_midpoint_rms{0.0};
  double coulomb_energy_kcal_mol{0.0};
  double coulomb_energy_kj_mol{0.0};
  double coulomb_energy_kt{0.0};
  std::string reaction_energy_source{"atom"};
  bool dgsolv_enabled{false};
  std::string dgsolv_status{"off"};
  double complex_energy_kcal_mol{0.0};
  double complex_energy_kj_mol{0.0};
  double mol1_energy_kcal_mol{0.0};
  double mol1_energy_kj_mol{0.0};
  double mol2_energy_kcal_mol{0.0};
  double mol2_energy_kj_mol{0.0};
  double deltaG_solv_kcal_mol{0.0};
  double deltaG_solv_kj_mol{0.0};
  double deltaG_solv_kt{0.0};
  bool pbsa_enabled{false};
  std::string pbsa_status{"off"};
  double pbsa_elec_kj_mol{0.0};
  double pbsa_vdw_kj_mol{0.0};
  double pbsa_mm_kj_mol{0.0};
  double pbsa_delta_sasa_a2{0.0};
  double pbsa_sasa_complex_a2{0.0};
  double pbsa_sasa_mol1_a2{0.0};
  double pbsa_sasa_mol2_a2{0.0};
  double pbsa_nonpolar_kj_mol{0.0};
  double pbsa_delta_solv_kj_mol{0.0};
  double pbsa_delta_solv_kcal_mol{0.0};
  double pbsa_no_entropy_kj_mol{0.0};
  double pbsa_no_entropy_kcal_mol{0.0};
  std::size_t pbsa_receptor_atom_count{0};
  std::size_t pbsa_ligand_atom_count{0};
  std::string pbsa_complex_status{"off"};
  std::string pbsa_receptor_status{"off"};
  std::string pbsa_ligand_status{"off"};
  int pbsa_receptor_solvent_iterations{0};
  int pbsa_receptor_vacuum_iterations{0};
  int pbsa_ligand_solvent_iterations{0};
  int pbsa_ligand_vacuum_iterations{0};
  int pbsa_receptor_mg_cycles{0};
  int pbsa_ligand_mg_cycles{0};
  double pbsa_receptor_solvent_rms{0.0};
  double pbsa_receptor_vacuum_rms{0.0};
  double pbsa_ligand_solvent_rms{0.0};
  double pbsa_ligand_vacuum_rms{0.0};
  double pbsa_receptor_solvent_max_dphi{0.0};
  double pbsa_receptor_vacuum_max_dphi{0.0};
  double pbsa_ligand_solvent_max_dphi{0.0};
  double pbsa_ligand_vacuum_max_dphi{0.0};
  double pbsa_time_complex_ms{0.0};
  double pbsa_time_receptor_ms{0.0};
  double pbsa_time_ligand_ms{0.0};
  double pbsa_time_components_ms{0.0};
  double pbsa_time_residue_ms{0.0};
  double pbsa_time_frame_ms{0.0};
  std::string status{"ok"};
};

class Analyzer {
public:
  Analyzer() = default;

  void configure(const arguments &args);
  [[nodiscard]] bool enabled() const noexcept { return args_.flag_pbe; }

  void begin_trajectory(std::size_t trajectory_id,
                        const std::string &source_label);
  void end_trajectory(std::fstream &log);

  FrameResult consume_frame(const MOL &mol1, const Multi_MOL_parameters &params1,
                            const MOL *mol2,
                            const Multi_MOL_parameters *params2,
                            std::size_t trajectory_id,
                            std::size_t frame_index,
                            std::size_t global_index,
                            TitanExecutor &executor, std::fstream &log,
                            const EnergyDecompositionProvider *energy_provider =
                                nullptr);

  void finalize(std::fstream &log);

private:
  struct Atom {
    double x{0.0};
    double y{0.0};
    double z{0.0};
    double q{0.0};
    double radius{0.0};
    std::string element;
    std::string atom_name;
    std::string resname;
    int resnr{0};
    std::string chain;
  };

  struct GridShape {
    double ox{0.0};
    double oy{0.0};
    double oz{0.0};
    double h{0.5};
    std::size_t nx{0};
    std::size_t ny{0};
    std::size_t nz{0};
    std::size_t n{0};
    std::size_t memory_bytes{0};
  };

  struct SolveStats {
    int iterations{0};
    double rms{0.0};
    double max_dphi{0.0};
    bool converged{false};
  };

  struct NewtonCorrectionStats {
    SolveStats sweeps{};
    double residual_rms{0.0};
    double residual_max{0.0};
    double jac_diag_min{0.0};
    double jac_diag_max{0.0};
    std::size_t screen_points{0};
    std::size_t active_nonlinear_points{0};
    int inner_levels{1};
    int inner_cycles{0};
    double inner_residual_before{0.0};
    double inner_residual_after{0.0};
    std::string inner_effective{"sor"};
    bool residual_decreased{false};
  };

  struct ChargeDiagnostics {
    double assigned_charge{0.0};
    double moment_error{0.0};
  };

  struct BoundaryDiagnostics {
    double rms{0.0};
    double max_abs{0.0};
  };

  struct MapSamplePoint {
    std::string kind{};
    std::size_t rank{0};
    std::size_t id{0};
    double hint{0.0};
  };

  struct ParityAudit {
    double eps_sum_min{0.0};
    double eps_sum_max{0.0};
    double eps_sum_mean{0.0};
    double midpoint_eps_sum_min{0.0};
    double midpoint_eps_sum_max{0.0};
    double midpoint_eps_sum_mean{0.0};
    double midpoint_eps_sum_interior_mean{0.0};
    std::size_t screen_points{0};
    double screen_sum{0.0};
    double screen_rms{0.0};
    double screen_max{0.0};
    double source_sum{0.0};
    double source_rms{0.0};
    double source_max{0.0};
    double residual_max_abs{0.0};
    double derivative_max_abs{0.0};
    double derivative_max_rel{0.0};
    std::size_t derivative_samples{0};
  };

  struct ReactionSolve {
    GridShape shape{};
    SolveStats solvent{};
    SolveStats vacuum{};
    int solvent_total_iterations{0};
    int vacuum_total_iterations{0};
    int focus_levels{1};
    double focus_initial_scale{0.0};
    double focus_initial_margin{0.0};
    std::string solver_mode{"sor"};
    int mg_levels{1};
    int mg_cycles{0};
    int mg_pre{0};
    int mg_post{0};
    int mg_coarse_iter{0};
    std::size_t solute_points{0};
    std::size_t solvent_points{0};
    std::size_t ion_accessible_points{0};
    std::size_t stern_excluded_points{0};
    std::size_t peak_memory_bytes{0};
    std::string pbe_model_effective{"linear"};
    std::string linear_seed_solver{"none"};
    std::string nl_solver{"picard"};
    std::string nl_inner_solver{"sor"};
    std::string nl_mg_mode{"off"};
    std::string focus_parent_grid{};
    std::string focus_child_grid{};
    std::string probe_status{"off"};
    int newton_inner_iterations{0};
    double newton_step_norm{0.0};
    double newton_residual_before{0.0};
    double newton_residual_after{0.0};
    std::string newton_line_search_status{"off"};
    double newton_correction_residual{0.0};
    double newton_correction_max_residual{0.0};
    double newton_accepted_lambda{0.0};
    int newton_rejected_steps{0};
    double newton_jac_diag_min{0.0};
    double newton_jac_diag_max{0.0};
    std::size_t newton_screen_points{0};
    std::size_t newton_active_nonlinear_points{0};
    std::string newton_stop_reason{"off"};
    int newton_inner_levels{1};
    int newton_inner_cycles{0};
    double newton_inner_residual_before{0.0};
    double newton_inner_residual_after{0.0};
    std::string newton_inner_effective{"off"};
    double gaussian_density_sum{0.0};
    std::size_t gaussian_cutoff_points{0};
    double gaussian_eps_min{0.0};
    double gaussian_eps_max{0.0};
    double gaussian_eps_mean{0.0};
    double gaussian_charge_norm_error{0.0};
    double rhs_scale{0.0};
    double charge_assigned{0.0};
    double charge_moment_error{0.0};
    double boundary_rms{0.0};
    double boundary_max{0.0};
    double linear_residual{0.0};
    int nl_iterations{0};
    int nl_blocks{0};
    double nl_relax{0.5};
    double nl_max_dphi{0.0};
    double nl_residual{0.0};
    std::size_t nl_clipped_points{0};
    int nl_continuation_steps{0};
    int nl_backtracks{0};
    std::string nl_status{"linear"};
    double energy_kcal_mol{0.0};
    double energy_kj_mol{0.0};
    double energy_kt{0.0};
    double solvent_phase_kcal_mol{0.0};
    double solvent_phase_kj_mol{0.0};
    double solvent_phase_kt{0.0};
    double vacuum_phase_kcal_mol{0.0};
    double vacuum_phase_kj_mol{0.0};
    double vacuum_phase_kt{0.0};
    double atom_corrected_kcal_mol{0.0};
    double atom_corrected_kj_mol{0.0};
    double atom_corrected_kt{0.0};
    double grid_corrected_kcal_mol{0.0};
    double grid_corrected_kj_mol{0.0};
    double grid_corrected_kt{0.0};
    double iscm_water_phase_kcal_mol{0.0};
    double iscm_water_phase_kj_mol{0.0};
    double iscm_water_phase_kt{0.0};
    double iscm_vacuum_phase_kcal_mol{0.0};
    double iscm_vacuum_phase_kj_mol{0.0};
    double iscm_vacuum_phase_kt{0.0};
    double iscm_corrected_kcal_mol{0.0};
    double iscm_corrected_kj_mol{0.0};
    double iscm_corrected_kt{0.0};
    std::size_t iscm_boundary_points{0};
    std::size_t iscm_charged_boundary_points{0};
    double iscm_induced_charge_sum{0.0};
    double iscm_induced_charge_l1{0.0};
    double iscm_induced_charge_l2{0.0};
    double iscm_induced_charge_max{0.0};
    std::string iscm_status{"off"};
    std::size_t delphi_face_crossings{0};
    std::size_t delphi_eps_face_matches{0};
    std::size_t delphi_boundary_grid_points{0};
    std::size_t delphi_charged_boundary_points{0};
    std::string delphi_iscm_status{"off"};
    std::string vdwms_enabled{"off"};
    std::size_t vdwms_boundary_initial{0};
    std::size_t vdwms_boundary_final{0};
    std::size_t vdwms_added{0};
    std::size_t vdwms_removed{0};
    std::size_t vdwms_exposed_points{0};
    double vdwms_midpoint_rms{0.0};
    double coulomb_energy_kcal_mol{0.0};
    double coulomb_energy_kj_mol{0.0};
    double coulomb_energy_kt{0.0};
    std::string reaction_energy_source{"atom"};
    ParityAudit parity{};
    std::vector<double> atom_reaction_kj_mol;
    bool converged{false};
  };

  struct FieldSolve {
    GridShape shape{};
    SolveStats stats{};
    int total_iterations{0};
    int levels{1};
    double initial_scale{0.0};
    double initial_margin{0.0};
    std::string solver_mode{"sor"};
    int mg_levels{1};
    int mg_cycles{0};
    int mg_pre{0};
    int mg_post{0};
    int mg_coarse_iter{0};
    std::size_t solute_points{0};
    std::size_t solvent_points{0};
    std::size_t ion_accessible_points{0};
    std::size_t stern_excluded_points{0};
    std::size_t peak_memory_bytes{0};
    std::string pbe_model_effective{"linear"};
    std::string linear_seed_solver{"none"};
    std::string nl_solver{"picard"};
    std::string nl_inner_solver{"sor"};
    std::string nl_mg_mode{"off"};
    std::string focus_parent_grid{};
    std::string focus_child_grid{};
    std::string probe_status{"off"};
    int newton_inner_iterations{0};
    double newton_step_norm{0.0};
    double newton_residual_before{0.0};
    double newton_residual_after{0.0};
    std::string newton_line_search_status{"off"};
    double newton_correction_residual{0.0};
    double newton_correction_max_residual{0.0};
    double newton_accepted_lambda{0.0};
    int newton_rejected_steps{0};
    double newton_jac_diag_min{0.0};
    double newton_jac_diag_max{0.0};
    std::size_t newton_screen_points{0};
    std::size_t newton_active_nonlinear_points{0};
    std::string newton_stop_reason{"off"};
    int newton_inner_levels{1};
    int newton_inner_cycles{0};
    double newton_inner_residual_before{0.0};
    double newton_inner_residual_after{0.0};
    std::string newton_inner_effective{"off"};
    double gaussian_density_sum{0.0};
    std::size_t gaussian_cutoff_points{0};
    double gaussian_eps_min{0.0};
    double gaussian_eps_max{0.0};
    double gaussian_eps_mean{0.0};
    double gaussian_charge_norm_error{0.0};
    double rhs_scale{0.0};
    double charge_assigned{0.0};
    double charge_moment_error{0.0};
    double boundary_rms{0.0};
    double boundary_max{0.0};
    double linear_residual{0.0};
    int nl_iterations{0};
    int nl_blocks{0};
    double nl_relax{0.5};
    double nl_max_dphi{0.0};
    double nl_residual{0.0};
    std::size_t nl_clipped_points{0};
    int nl_continuation_steps{0};
    int nl_backtracks{0};
    std::string nl_status{"linear"};
    ParityAudit parity{};
    std::vector<double> atom_phi;
    AlignedDoubleVector phi_grid;
    AlignedDoubleVector rho_grid;
    AlignedDoubleVector eps_sum_grid;
    AlignedDoubleVector eps_x_grid;
    AlignedDoubleVector eps_y_grid;
    AlignedDoubleVector eps_z_grid;
    AlignedDoubleVector screen_grid;
    std::vector<unsigned char> inside_grid;
    std::string vdwms_enabled{"off"};
    std::size_t vdwms_boundary_initial{0};
    std::size_t vdwms_boundary_final{0};
    std::size_t vdwms_added{0};
    std::size_t vdwms_removed{0};
    std::size_t vdwms_exposed_points{0};
    double vdwms_midpoint_rms{0.0};
  };

  struct ForensicAudit {
    std::size_t face_cross_x{0};
    std::size_t face_cross_y{0};
    std::size_t face_cross_z{0};
    std::size_t eps_face_count_x{0};
    std::size_t eps_face_count_y{0};
    std::size_t eps_face_count_z{0};
    std::size_t internal_face_count_x{0};
    std::size_t internal_face_count_y{0};
    std::size_t internal_face_count_z{0};
    double eps_face_min{0.0};
    double eps_face_max{0.0};
    double eps_face_mean{0.0};
    double eps_face_mean_x{0.0};
    double eps_face_mean_y{0.0};
    double eps_face_mean_z{0.0};
    double vdw_volume_estimate_A3{0.0};
    double rho_sum{0.0};
    double rho_l1{0.0};
    double rho_l2{0.0};
    double rho_mx{0.0};
    double rho_my{0.0};
    double rho_mz{0.0};
    double denom_min{0.0};
    double denom_max{0.0};
    double denom_mean{0.0};
    double screen_sum{0.0};
    double screen_max{0.0};
  };

  struct VdwmsAudit {
    std::string enabled{"off"};
    std::size_t boundary_initial{0};
    std::size_t boundary_final{0};
    std::size_t added{0};
    std::size_t removed{0};
    std::size_t exposed_points{0};
    double midpoint_rms{0.0};
  };

  struct DelphiVdwMediaMap {
    std::vector<std::int32_t> face_media;
    std::vector<std::int32_t> face_owner;
    VdwmsAudit audit{};
  };

  struct IscmSample {
    std::string kind{};
    std::size_t rank{0};
    std::size_t i{0};
    std::size_t j{0};
    std::size_t k{0};
    double x{0.0};
    double y{0.0};
    double z{0.0};
    double phi_center{0.0};
    double phi_neighbor_average{0.0};
    double charge_adjustment{0.0};
    double induced_charge{0.0};
    double atom_potential_sum{0.0};
    double energy_contribution_kt{0.0};
  };

  struct IscmAudit {
    std::size_t boundary_points{0};
    std::size_t charged_boundary_points{0};
    double induced_charge_sum{0.0};
    double induced_charge_l1{0.0};
    double induced_charge_l2{0.0};
    double induced_charge_max{0.0};
    double water_phase_kt{0.0};
    double water_phase_kj_mol{0.0};
    double water_phase_kcal_mol{0.0};
    double vacuum_phase_kt{0.0};
    double vacuum_phase_kj_mol{0.0};
    double vacuum_phase_kcal_mol{0.0};
    double corrected_kt{0.0};
    double corrected_kj_mol{0.0};
    double corrected_kcal_mol{0.0};
    std::size_t delphi_face_crossings{0};
    std::size_t delphi_eps_face_matches{0};
    std::string status{"off"};
    std::vector<IscmSample> samples;
  };

  struct PbsaComponents {
    double elec_kj_mol{0.0};
    double vdw_kj_mol{0.0};
    double mm_kj_mol{0.0};
    double sasa_complex_a2{0.0};
    double sasa_mol1_a2{0.0};
    double sasa_mol2_a2{0.0};
    double delta_sasa_a2{0.0};
    double nonpolar_kj_mol{0.0};
    double delta_solv_kj_mol{0.0};
    double no_entropy_kj_mol{0.0};
    std::string status{"ok"};
  };

  struct ScalarStats {
    std::size_t count{0};
    double mean{0.0};
    double m2{0.0};
    double min{0.0};
    double max{0.0};
  };

  struct PbsaEntropyResult {
    std::size_t samples{0};
    std::size_t filtered_samples{0};
    double ess{0.0};
    double minus_tds_kj_mol{0.0};
    std::size_t block_size{0};
    std::size_t block_count{0};
    double block_mean_kj_mol{0.0};
    double block_std_kj_mol{0.0};
    double block_drift_kj_mol{0.0};
    std::string status{"disabled"};
  };

  struct ResidueKey {
    std::string molecule{"unknown"};
    std::size_t chainnum{0};
    std::string chain{" "};
    std::size_t resnr{0};
    std::size_t respos{0};
    std::string resname{"UND"};
  };

  struct ResidueContribution {
    ResidueKey key{};
    double elec_kj_mol{0.0};
    double vdw_kj_mol{0.0};
    double mm_kj_mol{0.0};
    double polar_kj_mol{0.0};
    double nonpolar_kj_mol{0.0};
    double total_kj_mol{0.0};
  };

  struct ResidueFrameRow : ResidueContribution {
    std::size_t trajectory_id{0};
    std::size_t frame_index{0};
    std::size_t global_index{0};
    double time_ps{0.0};
  };

  struct ResidueStats {
    ResidueKey key{};
    std::size_t frames{0};
    ScalarStats elec{};
    ScalarStats vdw{};
    ScalarStats mm{};
    ScalarStats polar{};
    ScalarStats nonpolar{};
    ScalarStats total{};
  };

  struct MultigridLevel {
    GridShape shape{};
    AlignedDoubleVector phi;
    AlignedDoubleVector rho;
    AlignedDoubleVector rhs;
    AlignedDoubleVector eps_x;
    AlignedDoubleVector eps_y;
    AlignedDoubleVector eps_z;
    AlignedDoubleVector denom;
    AlignedDoubleVector screen;
    std::vector<unsigned char> inside;
    std::vector<unsigned char> ion_accessible;
    ChargeDiagnostics charge_diag{};
    BoundaryDiagnostics boundary_diag{};
  };

  struct NewtonCorrectionLevel {
    GridShape shape{};
    AlignedDoubleVector delta;
    AlignedDoubleVector rhs;
    AlignedDoubleVector eps_x;
    AlignedDoubleVector eps_y;
    AlignedDoubleVector eps_z;
    AlignedDoubleVector diag;
    AlignedDoubleVector screen;
    AlignedDoubleVector base_phi;
    double max_correction{1.0};
  };

  void collect_atoms(const MOL &mol, const Multi_MOL_parameters &params,
                     std::vector<Atom> &atoms, const char *label) const;
  void collect_test_pqr_atoms(std::vector<Atom> &atoms) const;
  [[nodiscard]] std::string inferred_element_symbol(const Atom &atom) const;
  [[nodiscard]] double effective_surface_radius(const Atom &atom) const;
  [[nodiscard]] std::string radius_media_source(const Atom &atom) const;
  [[nodiscard]] std::string effective_radius_convention() const;
  [[nodiscard]] bool uses_gaussian_model() const;
  [[nodiscard]] double gaussian_density_at(const std::vector<Atom> &atoms,
                                           double x, double y,
                                           double z) const;
  [[nodiscard]] double gaussian_dielectric(double density, double eps_in,
                                           double eps_out) const;
  [[nodiscard]] std::vector<Atom>
  adjusted_surface_atoms(const std::vector<Atom> &atoms) const;
  GridShape build_shape(const std::vector<Atom> &atoms) const;
  GridShape build_shape(const std::vector<Atom> &atoms, double scale,
                        double margin) const;
  void ensure_workspace(const GridShape &shape);
  void classify_surface(const GridShape &shape, const std::vector<Atom> &atoms,
                        TitanExecutor &executor);
  void classify_ion_accessibility(const GridShape &shape,
                                  const std::vector<Atom> &atoms,
                                  TitanExecutor &executor);
  ChargeDiagnostics assign_charges(const GridShape &shape,
                                   const std::vector<Atom> &atoms);
  [[nodiscard]] bool oracle_replay_applies(const GridShape &shape) const;
  void apply_oracle_surface_map(const GridShape &shape);
  [[nodiscard]] ChargeDiagnostics
  assign_oracle_charge_grid(const GridShape &shape,
                            const std::vector<Atom> &atoms);
  void apply_oracle_midpoint_epsilon_map(const GridShape &shape);
  [[nodiscard]] DelphiVdwMediaMap
  build_delphi_vdw_media_map(const GridShape &shape,
                             const std::vector<Atom> &surface_atoms,
                             TitanExecutor &executor) const;
  void apply_delphi_vdw_media_map(const GridShape &shape,
                                  const DelphiVdwMediaMap &media_map,
                                  double eps_in, double eps_out);
  [[nodiscard]] BoundaryDiagnostics
  apply_oracle_boundary_map(const GridShape &shape);
  void build_coefficients(const GridShape &shape,
                          const std::vector<Atom> &atoms, double eps_in,
                          double eps_out, double kappa,
                          TitanExecutor &executor);
  BoundaryDiagnostics set_boundary(const GridShape &shape,
                                   const std::vector<Atom> &atoms,
                                   double eps_out, double kappa);
  SolveStats solve_sor(const GridShape &shape, TitanExecutor &executor);
  SolveStats solve_nonlinear_sor_block(const GridShape &shape,
                                       const AlignedDoubleVector &source,
                                       int sweeps, double relax,
                                       double screen_scale,
                                       TitanExecutor &executor,
                                       std::size_t &clipped_points);
  double compute_linear_residual(const GridShape &shape) const;
  double compute_nonlinear_residual(const GridShape &shape,
                                    const AlignedDoubleVector &source,
                                    double screen_scale) const;
  std::size_t update_nonlinear_source(const GridShape &shape,
                                      AlignedDoubleVector &source,
                                      double alpha, double clip) const;
  double compute_nonlinear_residual_from_phi(
      const GridShape &shape, const AlignedDoubleVector &field,
      const AlignedDoubleVector &source, double screen_scale) const;
  double compute_nonlinear_residual_vector(
      const GridShape &shape, const AlignedDoubleVector &field,
      const AlignedDoubleVector &source, double screen_scale,
      AlignedDoubleVector &residual, TitanExecutor &executor) const;
  NewtonCorrectionStats
  solve_newton_correction_sor(const GridShape &shape,
                              const AlignedDoubleVector &residual,
                              double residual_rms, double alpha, double clip,
                              double screen_scale, AlignedDoubleVector &delta,
                              int max_sweeps, TitanExecutor &executor) const;
  NewtonCorrectionStats
  solve_newton_correction_vcycle(const GridShape &shape,
                                 const AlignedDoubleVector &residual,
                                 double residual_rms, double alpha,
                                 double clip, double screen_scale,
                                 AlignedDoubleVector &delta,
                                 int max_cycles) const;
  NewtonCorrectionStats
  solve_newton_correction(const GridShape &shape,
                          const AlignedDoubleVector &residual,
                          double residual_rms, double alpha, double clip,
                          double screen_scale, AlignedDoubleVector &delta,
                          int max_iterations, const std::string &inner_solver,
                          TitanExecutor &executor) const;
  FieldSolve solve_newton_nonlinear_field(
      const std::vector<Atom> &atoms, double eps_in, double eps_out,
      double kappa, TitanExecutor &executor, const FrameResult *map_result);
  double grid_reaction_energy_kcal(const GridShape &shape,
                                   const AlignedDoubleVector &solvent_phi,
                                   const AlignedDoubleVector &vacuum_phi) const;
  double grid_phase_energy_kcal(const GridShape &shape,
                                const AlignedDoubleVector &field) const;
  double atom_phase_energy_kcal(const std::vector<Atom> &atoms,
                                const std::vector<double> &atom_phi) const;
  double analytic_coulomb_energy_kcal(const std::vector<Atom> &atoms,
                                      double eps) const;
  FieldSolve solve_focused_field(const std::vector<Atom> &atoms, double eps_in,
                                 double eps_out, double kappa,
                                 TitanExecutor &executor,
                                 const FrameResult *map_result);
  FieldSolve solve_vcycle_field(const std::vector<Atom> &atoms, double eps_in,
                                double eps_out, double kappa,
                                TitanExecutor &executor,
                                const FrameResult *map_result);
  FieldSolve solve_linear_field(const std::vector<Atom> &atoms, double eps_in,
                                double eps_out, double kappa,
                                TitanExecutor &executor,
                                const FrameResult *map_result);
  FieldSolve solve_nonlinear_field(const std::vector<Atom> &atoms,
                                   double eps_in, double eps_out, double kappa,
                                   TitanExecutor &executor,
                                   const FrameResult *map_result);
  FieldSolve solve_physical_field(const std::vector<Atom> &atoms,
                                  double eps_in, double eps_out, double kappa,
                                  TitanExecutor &executor,
                                  const FrameResult *map_result);
  void prepare_multigrid_level(MultigridLevel &level, const GridShape &shape,
                               const std::vector<Atom> &atoms, double eps_in,
                               double eps_out, double kappa,
                               bool physical_rhs, bool physical_boundary,
                               TitanExecutor &executor);
  void copy_level_to_workspace(const MultigridLevel &level);
  SolveStats smooth_multigrid_level(MultigridLevel &level, int iterations,
                                    bool stop_on_convergence) const;
  double compute_multigrid_residual(const MultigridLevel &level,
                                    AlignedDoubleVector &residual) const;
  void restrict_residual_to_coarse(const MultigridLevel &fine,
                                   const AlignedDoubleVector &fine_residual,
                                   MultigridLevel &coarse) const;
  void prolongate_add_correction(const MultigridLevel &coarse,
                                 MultigridLevel &fine) const;
  void run_multigrid_vcycle(std::vector<MultigridLevel> &levels,
                            std::size_t level_index);
  void prepare_newton_correction_hierarchy(
      const GridShape &shape, const AlignedDoubleVector &residual,
      double alpha, double clip, double screen_scale,
      std::vector<NewtonCorrectionLevel> &levels) const;
  SolveStats smooth_newton_correction_level(NewtonCorrectionLevel &level,
                                            int iterations,
                                            bool stop_on_convergence) const;
  double compute_newton_correction_level_residual(
      const NewtonCorrectionLevel &level,
      AlignedDoubleVector &residual) const;
  void restrict_newton_correction_residual(
      const NewtonCorrectionLevel &fine,
      const AlignedDoubleVector &fine_residual,
      NewtonCorrectionLevel &coarse) const;
  void prolongate_add_newton_correction(
      const NewtonCorrectionLevel &coarse, NewtonCorrectionLevel &fine) const;
  void run_newton_correction_vcycle(
      std::vector<NewtonCorrectionLevel> &levels,
      std::size_t level_index) const;
  GridShape coarsen_shape(const GridShape &fine) const;
  [[nodiscard]] int effective_vcycle_levels(const GridShape &fine_shape) const;
  [[nodiscard]] std::size_t
  multigrid_memory_bytes(const std::vector<MultigridLevel> &levels) const;
  ReactionSolve solve_reaction_field(const std::vector<Atom> &atoms,
                                     TitanExecutor &executor,
                                     const FrameResult *map_result);
  PbsaComponents
  compute_pbsa_components(const MOL &mol1, const Multi_MOL_parameters &params1,
                          const MOL &mol2,
                          const Multi_MOL_parameters &params2,
                          const std::vector<Atom> &mol1_atoms,
                          const std::vector<Atom> &mol2_atoms,
                          const std::vector<Atom> &complex_atoms,
                          double deltaG_solv_kj_mol,
                          TitanExecutor &executor,
                          const EnergyDecompositionProvider *energy_provider)
      const;
  [[nodiscard]] double compute_sasa_total(const std::vector<Atom> &atoms,
                                          TitanExecutor &executor) const;
  [[nodiscard]] std::vector<double>
  compute_sasa_per_atom(const std::vector<Atom> &atoms,
                        TitanExecutor &executor) const;
  [[nodiscard]] std::vector<ResidueFrameRow> compute_pbsa_residue_rows(
      const MOL &mol1, const MOL &mol2, const std::vector<Atom> &mol1_atoms,
      const std::vector<Atom> &mol2_atoms,
      const std::vector<Atom> &complex_atoms,
      const ReactionSolve &complex_solve, const ReactionSolve &mol1_solve,
      const ReactionSolve &mol2_solve, const FrameResult &result,
      TitanExecutor &executor,
      const EnergyDecompositionProvider *energy_provider) const;
  void accumulate_pbsa_residue_rows(
      const FrameResult &result, const std::vector<ResidueFrameRow> &rows);
  std::vector<double> interpolate_atoms(const GridShape &shape,
                                        const std::vector<Atom> &atoms) const;
  bool interpolate_field_value(const GridShape &shape,
                               const AlignedDoubleVector &field, double x,
                               double y, double z, double &value) const;
  void seed_from_focus_source(const GridShape &shape,
                              const GridShape &source_shape,
                              const AlignedDoubleVector &source_phi,
                              TitanExecutor &executor);
  void seed_from_focus_source(const GridShape &shape,
                              const GridShape &source_shape,
                              const AlignedDoubleVector &source_phi,
                              TitanExecutor &executor, bool boundary_only);
  [[nodiscard]] int effective_focus_levels() const noexcept;
  [[nodiscard]] std::string effective_solver_mode() const;
  [[nodiscard]] double effective_nl_relax() const noexcept;
  [[nodiscard]] double nonlinear_alpha() const noexcept;
  [[nodiscard]] double rhs_scale(const GridShape &shape) const noexcept;
  [[nodiscard]] bool nonlinear_dimensionless_mode() const noexcept;
  [[nodiscard]] double nonlinear_rhs_scale(const GridShape &shape,
                                           double alpha) const noexcept;
  [[nodiscard]] double nonlinear_argument_value(double field_value,
                                                double alpha) const noexcept;
  [[nodiscard]] double nonlinear_source_value(double phi, double alpha,
                                              double clip) const noexcept;
  [[nodiscard]] double nonlinear_source_derivative(double phi, double alpha,
                                                   double clip) const noexcept;
  void scale_nonlinear_field_to_solve_units(AlignedDoubleVector &field,
                                            double alpha) const;
  void scale_nonlinear_field_to_physical_units(AlignedDoubleVector &field,
                                               double alpha) const;
  [[nodiscard]] ParityAudit compute_parity_audit(
      const GridShape &shape, const AlignedDoubleVector *source,
      const AlignedDoubleVector *residual) const;
  [[nodiscard]] ForensicAudit compute_forensic_audit(
      const GridShape &shape) const;
  [[nodiscard]] IscmAudit compute_iscm_audit(
      const GridShape &shape, const std::vector<Atom> &atoms,
      const AlignedDoubleVector &solvent_phi,
      const AlignedDoubleVector &solvent_screen,
      const std::vector<unsigned char> &solvent_inside,
      const AlignedDoubleVector &solvent_eps_x,
      const AlignedDoubleVector &solvent_eps_y,
      const AlignedDoubleVector &solvent_eps_z,
      double grid_vacuum_phase_kcal_mol) const;
  void copy_iscm_audit_to_solve(const IscmAudit &audit,
                                ReactionSolve &solve) const;
  void update_newton_derivative_audit(const GridShape &shape,
                                      const AlignedDoubleVector &source,
                                      const AlignedDoubleVector &delta,
                                      double alpha, double clip,
                                      double screen_scale,
                                      ParityAudit &audit) const;
  void write_energy_header_if_needed();
  void append_energy_row(const FrameResult &result);
  void maybe_write_parity_stage(const FrameResult &result,
                                const ReactionSolve &solve) const;
  [[nodiscard]] bool barnase_forensic_enabled() const;
  void maybe_write_barnase_forensic_audits(
      const FrameResult &result, const std::vector<Atom> &atoms,
      const FieldSolve &solvent, const FieldSolve &vacuum,
      const ReactionSolve &solve, const ForensicAudit &solvent_audit) const;
  void maybe_write_iscm_surface_charge_samples(
      const FrameResult &result, const IscmAudit &audit) const;
  void maybe_write_atom_phi_audit(const FrameResult &result,
                                  const std::vector<Atom> &atoms,
                                  const FieldSolve &solvent,
                                  const FieldSolve &vacuum,
                                  const ReactionSolve &solve) const;
  void maybe_write_linear_forensic_points(
      const FrameResult &result, const std::vector<Atom> &atoms,
      const FieldSolve &solvent, const FieldSolve &vacuum,
      const ReactionSolve &solve) const;
  void maybe_write_parity_points(const FrameResult &result,
                                 const GridShape &shape,
                                 const AlignedDoubleVector &source,
                                 const AlignedDoubleVector &residual,
                                 const AlignedDoubleVector &delta,
                                 double alpha, double clip,
                                 double screen_scale) const;
  [[nodiscard]] bool parity_map_dump_enabled() const noexcept;
  [[nodiscard]] std::vector<MapSamplePoint>
  collect_map_sample_points(const GridShape &shape,
                            const std::vector<Atom> &atoms,
                            const ReactionSolve &solve) const;
  void maybe_write_parity_map_dumps(const FrameResult &result,
                                    const std::vector<Atom> &atoms,
                                    const FieldSolve &solvent,
                                    const FieldSolve &vacuum,
                                    const ReactionSolve &solve) const;
  void accumulate_pbsa(const FrameResult &result);
  [[nodiscard]] PbsaEntropyResult compute_pbsa_entropy() const;
  [[nodiscard]] PbsaEntropyResult
  compute_pbsa_entropy_for_values(const std::vector<double> &values) const;
  void write_pbsa_energy_csv(const PbsaEntropyResult &entropy) const;
  void write_pbsa_decomposition(const PbsaEntropyResult &entropy) const;
  void write_pbsa_residue_outputs() const;
  void write_pbsa_by_trajectory() const;
  void write_pbsa_ie_diagnostics(const PbsaEntropyResult &global_entropy) const;
  void write_pbsa_performance() const;
  void write_pbsa_summary(std::fstream &log);
  static void add_scalar(ScalarStats &stats, double value);
  static double scalar_stddev(const ScalarStats &stats);
  void maybe_write_maps(const GridShape &shape, const FrameResult &result,
                        bool write_phi_now);
  [[nodiscard]] std::string maybe_write_probe_sites(
      const FrameResult &result, const std::vector<Atom> &atoms,
      const ReactionSolve &solve);
  void write_cube(const std::string &filename, const GridShape &shape,
                  const std::string &field_name, bool epsilon_field,
                  bool ion_field = false) const;
  [[nodiscard]] bool should_write_map() const;
  [[nodiscard]] std::string map_suffix(const FrameResult &result) const;

  PbeArgs args_{};
  std::string current_source_label_{};
  std::size_t current_trajectory_id_{0};

  std::ofstream energy_out_;
  bool energy_header_written_{false};
  mutable bool parity_stage_header_written_{false};
  mutable bool parity_points_header_written_{false};
  mutable bool forensic_inputs_header_written_{false};
  mutable bool forensic_surface_header_written_{false};
  mutable bool forensic_operator_header_written_{false};
  mutable bool radius_media_audit_header_written_{false};
  mutable bool atom_phi_audit_header_written_{false};
  mutable bool iscm_surface_samples_header_written_{false};
  mutable bool field_samples_header_written_{false};
  mutable bool boundary_samples_header_written_{false};
  mutable bool charge_grid_audit_header_written_{false};
  bool wrote_first_map_{false};

  std::vector<Atom> atoms_;
  AlignedDoubleVector phi_;
  AlignedDoubleVector rho_;
  AlignedDoubleVector eps_x_;
  AlignedDoubleVector eps_y_;
  AlignedDoubleVector eps_z_;
  AlignedDoubleVector denom_;
  AlignedDoubleVector screen_;
  std::vector<unsigned char> inside_;
  std::vector<unsigned char> ion_accessible_;
  double last_gaussian_density_sum_{0.0};
  std::size_t last_gaussian_cutoff_points_{0};
  double last_gaussian_eps_min_{0.0};
  double last_gaussian_eps_max_{0.0};
  double last_gaussian_eps_mean_{0.0};
  VdwmsAudit last_vdwms_audit_{};

  std::size_t frames_total_{0};
  std::size_t frames_converged_{0};
  std::size_t frames_failed_{0};
  std::size_t max_grid_points_seen_{0};
  std::size_t max_memory_bytes_seen_{0};
  std::size_t max_solute_points_seen_{0};
  std::size_t max_solvent_points_seen_{0};
  std::size_t max_ion_accessible_points_seen_{0};
  std::size_t max_stern_excluded_points_seen_{0};
  int max_focus_levels_seen_{1};
  int max_mg_levels_seen_{1};
  int max_mg_cycles_seen_{0};
  std::size_t nonlinear_frames_{0};
  std::size_t nonlinear_converged_{0};
  std::size_t nonlinear_failed_{0};
  int max_nl_blocks_seen_{0};
  int max_nl_iterations_seen_{0};
  std::size_t max_nl_clipped_points_seen_{0};
  double max_nl_dphi_seen_{0.0};
  double max_nl_residual_seen_{0.0};
  int max_newton_inner_iterations_seen_{0};
  double max_newton_step_norm_seen_{0.0};
  double max_newton_residual_before_seen_{0.0};
  double max_newton_residual_after_seen_{0.0};
  double max_newton_correction_residual_seen_{0.0};
  double max_newton_correction_max_residual_seen_{0.0};
  int max_newton_rejected_steps_seen_{0};
  int max_newton_inner_levels_seen_{1};
  int max_newton_inner_cycles_seen_{0};
  double max_newton_inner_residual_before_seen_{0.0};
  double max_newton_inner_residual_after_seen_{0.0};
  double max_linear_residual_seen_{0.0};
  double max_boundary_rms_seen_{0.0};
  double max_boundary_max_seen_{0.0};
  double mg_cycle_mean_{0.0};
  double energy_mean_{0.0};
  double energy_m2_{0.0};
  double energy_min_{0.0};
  double energy_max_{0.0};
  double iter_mean_{0.0};
  std::size_t dgsolv_frames_{0};
  double dgsolv_mean_{0.0};
  double dgsolv_m2_{0.0};
  double dgsolv_min_{0.0};
  double dgsolv_max_{0.0};
  std::size_t pbsa_frames_{0};
  std::size_t pbsa_failed_{0};
  ScalarStats pbsa_elec_stats_{};
  ScalarStats pbsa_vdw_stats_{};
  ScalarStats pbsa_mm_stats_{};
  ScalarStats pbsa_delta_sasa_stats_{};
  ScalarStats pbsa_nonpolar_stats_{};
  ScalarStats pbsa_polar_stats_{};
  ScalarStats pbsa_no_entropy_stats_{};
  std::vector<double> pbsa_ie_mm_energies_kj_;
  std::vector<FrameResult> pbsa_rows_;
  std::vector<ResidueFrameRow> pbsa_residue_frame_rows_;
  std::map<std::string, ResidueStats> pbsa_residue_stats_;
  std::map<std::string, ResidueStats> pbsa_residue_traj_stats_;
  double pbsa_residue_closure_max_abs_kj_{0.0};
  double pbsa_residue_closure_elec_max_abs_kj_{0.0};
  double pbsa_residue_closure_vdw_max_abs_kj_{0.0};
  double pbsa_residue_closure_polar_max_abs_kj_{0.0};
  double pbsa_residue_closure_nonpolar_max_abs_kj_{0.0};
};

} // namespace TitanPBE

#endif
