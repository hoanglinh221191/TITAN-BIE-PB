// titan_gist.h — Grid Inhomogeneous Solvation Theory (GIST) Water Analysis
// Computes per-voxel density, enthalpy (E_sw), and entropy for solvation
// thermodynamics. Outputs OpenDX maps for VMD/PyMOL visualization.
#ifndef TITAN_GIST_H
#define TITAN_GIST_H

#include "MOL.h"
#include "Multi_MOL.h"
#include "info.h"

#include <Eigen/Dense>
#include <array>
#include <string>
#include <thread>
#include <vector>

// ============================================================================
// GistVoxel — one cell of the 3D grid (cache-line aligned)
// ============================================================================
struct alignas(64) GistVoxel {
  uint32_t count{0};       // water visit counter (serial per-frame)
  uint32_t pad_{0};        // align doubles to 8-byte boundary
  double energy_sw{0.0};   // Σ E_sw (kJ/mol)
  double energy_ww{0.0};   // Σ E_ww (kJ/mol)  [Phase 2]
  double dipole_x{0.0};    // Σ d̂_x (unit vector)
  double dipole_y{0.0};    // Σ d̂_y
  double dipole_z{0.0};    // Σ d̂_z
  // Total: 4+4 + 5*8 = 48 bytes → padded to 64 by alignas
};

// ============================================================================
// TitanGIST — main analysis class
// ============================================================================
class TitanGIST {
public:
  TitanGIST() = default;

  // Initialize grid from ligand reference frame + solute parameters
  void init(const MOL &lig_ref, double padding, double resolution,
            double temperature, double cutoff, size_t n_threads,
            const std::string &water_model,
            const Multi_MOL_parameters &rec_params,
            const std::vector<mol_info> &rec_info,
            const Multi_MOL_parameters &lig_params,
            const std::vector<mol_info> &lig_info);

  // Per-frame: rotate water into ref frame, bin, compute E_sw
  void process_frame(const MOL &pro, const MOL &lig, const MOL &water,
                     const std::vector<mol_info> &pro_info,
                     const Eigen::MatrixXd &ref_ca_centered,
                     const Eigen::Matrix3d &R, bool is_frame_zero,
                     bool tica_available);

  // Post-processing: compute thermodynamics and export .dx files
  void export_dx_maps(const std::string &prefix, int n_frames,
                      int min_count) const;

  int n_voxels() const { return nx_ * ny_ * nz_; }

  // BUG2 fix: set DCCM's raw Cα indices for TICA centroid
  void set_tica_ca_indices(const std::vector<int> &dccm_ca_raw_indices) {
    tica_ca_indices_.clear();
    tica_ca_indices_.reserve(dccm_ca_raw_indices.size());
    for (int idx : dccm_ca_raw_indices)
      tica_ca_indices_.push_back(static_cast<size_t>(idx));
  }

private:
  // Grid
  std::vector<GistVoxel> grid_;
  int nx_{0}, ny_{0}, nz_{0};
  double voxel_size_{0.5};
  Eigen::Vector3d origin_{0, 0, 0};

  // Water model constants
  double q_O_{0}, q_H_{0}, sigma_O_{0}, eps_O_{0};
  double rho_bulk_{0}, E_ww_bulk_{0}, temperature_{300.0};

  // Coulomb/LJ cutoff
  double cut2_{81.0};  // cutoff² (default 9² = 81 Å²)

  // Thread count
  size_t nt_{1};

  // Merged solute (rec + lig) heavy atom arrays — contiguous for AVX2
  // Layout: [0, n_rec_heavy_) = receptor, [n_rec_heavy_, n_solute_) = ligand
  // All arrays padded to multiple of 4 (pad values = 0) for clean AVX2 loop
  std::vector<size_t> rec_heavy_idx_;   // original atom indices in rec
  std::vector<size_t> lig_heavy_idx_;   // original atom indices in lig
  size_t n_rec_heavy_{0};
  size_t n_solute_{0};                  // total length (padded)
  size_t n_solute_real_{0};             // unpadded count

  std::vector<double> solute_x_, solute_y_, solute_z_;     // coords
  std::vector<double> solute_charge_;                       // partial charge
  std::vector<double> solute_sigma_half_;  // precomputed: 0.5*(σ_O + σ_i)
  std::vector<double> solute_eps_mixed_;   // precomputed: sqrt(ε_O * ε_i)

  // Kabsch reference (built from Cα if TICA unavailable)
  Eigen::MatrixXd self_ref_centered_;
  Eigen::Vector3d ref_com_{0, 0, 0};
  bool ref_built_{false};
  bool ref_frame_locked_{false};    // BUG1: prevent re-init on multi-file
  std::vector<size_t> ca_indices_;  // Cα atom indices (all CA/BB/BB1)
  std::vector<size_t> tica_ca_indices_; // BUG2: DCCM's Cα subset for TICA

  // Per-frame scratch buffer
  std::vector<std::array<double, 3>> wat_buf_;

  // Helpers
  inline int voxel_idx(double x, double y, double z) const noexcept;

  double compute_esw(double ox, double oy, double oz, double h1x, double h1y,
                     double h1z, double h2x, double h2y, double h2z) const;

  void write_dx(const std::string &path, const std::vector<double> &data,
                const std::string &quantity_name) const;
};

#endif // TITAN_GIST_H
