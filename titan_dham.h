// ============================================================================
// TitanDHAM — Dynamic Histogram Analysis Method
// Reweights REMD/H-REMD MSM to a target temperature using energy information.
//
// Algorithm: Kong & Brooks, J. Chem. Phys. 105(6), 2414 (1996)
//            dHAM formulation per Rosta & Hummer, JCTC 2015
//
// Dependencies: Eigen/Dense, TitanUI, parse_xvg.h
// ============================================================================
#pragma once

#include "parse_xvg.h"
#include "titan_ui.h"
#include <Eigen/Dense>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// ============================================================================
// Constants
// ============================================================================
static constexpr double DHAM_kB_kJ = 0.008314463;  // kJ/(mol·K)
static constexpr double DHAM_LOG_ZERO_GUARD = 1e-300;

// Row-major layout for AVX2-friendly row scans in iterate_once()
using MatrixXdRM = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                                 Eigen::RowMajor>;

// ============================================================================
// TitanDHAM
// ============================================================================
class TitanDHAM {
public:
    // -------------------------------------------------------------------------
    // Constructor
    //   n_states:   number of microstates (K-means clusters)
    //   n_replicas: number of REMD replicas
    //   temps:      temperature of each replica (K), length = n_replicas
    //   target_temp: temperature to reweight to (K)
    // -------------------------------------------------------------------------
    TitanDHAM(int n_states, int n_replicas,
               const std::vector<double>& temps,
               double target_temp = 300.0);

    // -------------------------------------------------------------------------
    // Data loading
    // -------------------------------------------------------------------------

    // Accumulate one transition i→j for replica k (called after K-means)
    void accumulate_counts(int replica_k, int state_i, int state_j);

    // Set the mean potential energy (kJ/mol) of microstate i at replica k
    void set_state_energy(int replica_k, int state_i, double mean_energy_kJ);

    // -------------------------------------------------------------------------
    // Solve self-consistent equations (WHAM-like iteration)
    //   max_iter:  maximum number of iterations
    //   tol:       convergence tolerance on free energy shifts f_k
    // Returns true if converged.
    // -------------------------------------------------------------------------
    bool solve(int max_iter = 5000, double tol = 1e-6);

    // -------------------------------------------------------------------------
    // Results
    // -------------------------------------------------------------------------

    // Stationary distribution at target temperature (length = n_states)
    const Eigen::VectorXd& get_stationary_distribution() const { return p_i_; }

    // Unbiased transition matrix at target temperature (n_states × n_states)
    Eigen::MatrixXd get_unbiased_transition_matrix() const;

    // Free energy per microstate at target temperature (kJ/mol, relative to min)
    Eigen::VectorXd get_free_energy_profile() const;

    // Per-replica free energy shifts f_k (dimensionless)
    const Eigen::VectorXd& get_free_energy_shifts() const { return f_k_; }

    // Compute per-frame weights W_t = p_i / N_i for each replica trajectory.
    std::unordered_map<int, std::vector<double>>
    get_frame_weights(
        const std::unordered_map<int, std::vector<int>>& per_replica_dtrajs) const;

    // Flatten per-replica frame weights in replica order 0..n_replicas-1.
    std::vector<double>
    get_flattened_frame_weights(
        const std::unordered_map<int, std::vector<int>>& per_replica_dtrajs,
        int n_replicas) const;

    // Export results to files
    // eu_par: energy unit divisor (1.0 = kJ/mol, 4.184 = kcal/mol)
    void write_results(const std::string& pi_file,
                       const std::string& fes_file,
                       const std::string& tpm_file,
                       std::fstream& log,
                       double eu_par = 1.0) const;

    // -------------------------------------------------------------------------
    // Static helpers — public so postproc code can call them
    // -------------------------------------------------------------------------

    // Parse GROMACS .xvg file → {time_ps, energy_kJ} pairs
    // Uses high-perf TitanIO::parse_xvg_time_value under the hood
    static std::vector<std::pair<double,double>>
    parse_xvg(const std::string& path);

    // Map energy values to dtraj timestamps via nearest-neighbour matching.
    // Returns vector<double> of same length as dtraj_times_ps.
    // NaN is inserted where no match is found within 0.5*dt_energy_ps.
    static std::vector<double>
    align_energy_to_dtraj(const std::vector<double>& dtraj_times_ps,
                          const std::vector<std::pair<double,double>>& xvg_data,
                          double dt_energy_ps);

    // Compute per-state mean energy from dtraj + aligned energy array.
    // Returns matrix [n_states × n_replicas].
    // Frame-state pairs where energy is NaN are excluded from average.
    static Eigen::MatrixXd
    compute_state_energies(
        const std::unordered_map<int, std::vector<int>>&    per_replica_dtrajs,
        const std::unordered_map<int, std::vector<double>>& per_replica_energies,
        int n_states,
        int n_replicas,
        int num_threads = 1);

    // Read list of file paths from a text file (one path per line)
    static std::vector<std::string>
    read_file_list(const std::string& list_path);

    // Read temperatures from a text file (one temperature per line, in K)
    static std::vector<double>
    read_temps(const std::string& temps_path);

private:
    int n_states_;
    int n_replicas_;
    double beta_target_;               // 1 / (kB * T_target)
    std::vector<double> betas_;        // 1 / (kB * T_k) for each replica

    // Sparse count maps — one per replica:
    // C_k_sparse_[k][i * n_states_ + j] = transitions i→j in replica k
    // This avoids storing dense n_states x n_states matrices for every replica.
    std::vector<std::unordered_map<std::uint64_t, double>> C_k_sparse_;

    // c_out_(i,k) = sum_j C_k(i,j) = row sums (outgoing counts)
    // Maintained online during accumulate_counts() and rebuilt from sparse maps
    // before solve() as a consistency check.
    MatrixXdRM  c_out_;              // [n_states × n_replicas]  RowMajor

    // u_ik_(i,k) = beta_k * U_k(i)  — dimensionless reduced potential
    MatrixXdRM  u_ik_;               // [n_states × n_replicas]  RowMajor

    // Pre-computed log(c_out + guard) — invariant during solver iterations
    MatrixXdRM  log_c_;              // [n_states × n_replicas]  RowMajor

    // Optimized variables
    Eigen::VectorXd f_k_;             // free energy shift per replica
    Eigen::VectorXd p_i_;             // stationary distribution at target T

    bool solved_ = false;

    // Rebuild c_out_ from C_k_ row sums
    void rebuild_c_out();

    // One iteration of DHAM self-consistent equations using LSE trick
    double iterate_once();
};

// ============================================================================
// Inline implementation of I/O helpers
// ============================================================================

inline std::vector<std::pair<double,double>>
TitanDHAM::parse_xvg(const std::string& path) {
    auto result = TitanIO::parse_xvg_time_value(path, 0, 1);
    if (!result) {
        TitanUI::print_error("DHAM: " + result.error());
        return {};
    }
    return std::move(result.value());
}

inline std::vector<double>
TitanDHAM::align_energy_to_dtraj(
    const std::vector<double>& dtraj_times_ps,
    const std::vector<std::pair<double,double>>& xvg_data,
    double dt_energy_ps)
{
    const double nan_val = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> aligned(dtraj_times_ps.size(), nan_val);
    if (xvg_data.empty()) return aligned;

    // Compute actual dt_energy from data if not provided or zero
    double dt_e = dt_energy_ps;
    if (dt_e <= 0.0 && xvg_data.size() > 1)
        dt_e = xvg_data[1].first - xvg_data[0].first;
    if (dt_e <= 0.0) dt_e = 2.0; // fallback: 2 ps
    const double match_tol =
        0.5 * dt_e +
        1e-9 * std::max(1.0, std::abs(dt_e));

    size_t j = 0;
    for (size_t i = 0; i < dtraj_times_ps.size(); ++i) {
        double t_target = dtraj_times_ps[i];
        // Fix B: Reset j if target time goes backwards (trajectory restart)
        if (j > 0 && xvg_data[j].first > t_target + match_tol)
            j = 0;
        while (j + 1 < xvg_data.size() &&
               std::abs(xvg_data[j+1].first - t_target) <
               std::abs(xvg_data[j].first   - t_target))
            ++j;
        if (std::abs(xvg_data[j].first - t_target) <= match_tol)
            aligned[i] = xvg_data[j].second;
    }

    size_t nan_count = 0;
    for (double v : aligned) if (std::isnan(v)) ++nan_count;
    if (nan_count > aligned.size() / 10) {
        TitanUI::print_warning(
            "DHAM: " + std::to_string(nan_count) + "/" +
            std::to_string(aligned.size()) +
            " trajectory frames could not be matched to energy values.");
    }
    return aligned;
}

inline std::vector<std::string>
TitanDHAM::read_file_list(const std::string& list_path) {
    std::vector<std::string> result;
    std::ifstream f(list_path);
    if (!f.is_open()) {
        TitanUI::print_error("DHAM: Cannot open file list: " + list_path);
        return result;
    }
    std::string line;
    while (std::getline(f, line)) {
        size_t start = line.find_first_not_of(" \t\r\n");
        if (start == std::string::npos) continue;
        size_t end = line.find_last_not_of(" \t\r\n");
        result.push_back(line.substr(start, end - start + 1));
    }
    return result;
}

inline std::vector<double>
TitanDHAM::read_temps(const std::string& temps_path) {
    std::vector<double> temps;
    std::ifstream f(temps_path);
    if (!f.is_open()) {
        TitanUI::print_error("DHAM: Cannot open temps file: " + temps_path);
        return temps;
    }
    double t;
    while (f >> t) {
        if (t <= 0.0) {
            TitanUI::print_error("DHAM: Invalid temperature <= 0: " + std::to_string(t));
            continue;
        }
        temps.push_back(t);
    }
    return temps;
}
