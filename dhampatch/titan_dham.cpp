// ============================================================================
// TitanDHAM — Dynamic Histogram Analysis Method (Implementation)
// ============================================================================
#include "titan_dham.h"
#include <algorithm>
#include <iomanip>
#include <numeric>
#include <stdexcept>

// ============================================================================
// Constructor
// ============================================================================
TitanDHAM::TitanDHAM(int n_states, int n_replicas,
                       const std::vector<double>& temps,
                       double target_temp)
    : n_states_(n_states), n_replicas_(n_replicas)
{
    if (static_cast<int>(temps.size()) != n_replicas) {
        TitanUI::print_error("TitanDHAM: temps.size() (" +
            std::to_string(temps.size()) + ") != n_replicas (" +
            std::to_string(n_replicas) + ").");
        throw std::invalid_argument("TitanDHAM: temps/n_replicas mismatch.");
    }

    beta_target_ = 1.0 / (DHAM_kB_kJ * target_temp);

    betas_.resize(n_replicas);
    for (int k = 0; k < n_replicas; ++k)
        betas_[k] = 1.0 / (DHAM_kB_kJ * temps[k]);

    // Pre-allocate count matrices
    C_k_.resize(n_replicas,
                Eigen::MatrixXd::Zero(n_states, n_states));

    c_out_  = Eigen::MatrixXd::Zero(n_states, n_replicas);
    u_ik_   = Eigen::MatrixXd::Zero(n_states, n_replicas);

    // Initialize f_k to zero (unconstrained — first replica anchored during solve)
    f_k_ = Eigen::VectorXd::Zero(n_replicas);
    p_i_ = Eigen::VectorXd::Constant(n_states, 1.0 / n_states);
}

// ============================================================================
// accumulate_counts
// ============================================================================
void TitanDHAM::accumulate_counts(int replica_k, int state_i, int state_j) {
    if (replica_k < 0 || replica_k >= n_replicas_ ||
        state_i  < 0 || state_i  >= n_states_   ||
        state_j  < 0 || state_j  >= n_states_) return;
    C_k_[static_cast<size_t>(replica_k)](state_i, state_j) += 1.0;
}

// ============================================================================
// set_state_energy
// ============================================================================
void TitanDHAM::set_state_energy(int replica_k, int state_i,
                                   double mean_energy_kJ) {
    if (replica_k < 0 || replica_k >= n_replicas_ ||
        state_i  < 0 || state_i  >= n_states_) return;
    // u_ik(i,k) = beta_k * U_k(i)  (dimensionless reduced potential)
    u_ik_(state_i, static_cast<int>(replica_k)) =
        betas_[static_cast<size_t>(replica_k)] * mean_energy_kJ;
}

// ============================================================================
// rebuild_c_out — row sums of each C_k
// ============================================================================
void TitanDHAM::rebuild_c_out() {
    for (int k = 0; k < n_replicas_; ++k)
        c_out_.col(k) = C_k_[static_cast<size_t>(k)].rowwise().sum();
}

// ============================================================================
// iterate_once — one DHAM self-consistent step using LSE trick
// Returns max |Δf_k| for convergence check.
// ============================================================================
double TitanDHAM::iterate_once() {
    // ------------------------------------------------------------
    // Step 1: Update p_i (stationary distribution at target T)
    //
    // p_i(i) ∝ Σ_k [ c_out(i,k) / Σ_k' c_out(i,k') * exp(f_k - u_ik(i,k)) ]
    //
    // Using LSE: for each state i, compute
    //   log_terms(i,k) = f_k(k) - u_ik(i,k) + log(c_out(i,k) + guard)
    // then
    //   log_denom(i) = LSE over k of log_terms(i,k)
    //   p_i(i) ∝ exp(log_denom(i) + (beta_target - beta_k) * U_k(i))
    //
    // The physically correct target-T reweighting adds
    //   (beta_target - beta_k) * U_k(i) per replica, then averages.
    // ------------------------------------------------------------

    // log_c(i,k) = log(c_out(i,k) + guard)
    Eigen::ArrayXXd log_c =
        (c_out_.array() + DHAM_LOG_ZERO_GUARD).log();  // [n_states × n_replicas]

    // Broadcast f_k_ as a row vector, replicate across states
    Eigen::ArrayXXd fk_mat =
        f_k_.transpose().array().replicate(n_states_, 1);  // [n_states × n_replicas]

    // log_terms(i,k) = f_k(k) - u_ik(i,k) + log_c(i,k)
    Eigen::ArrayXXd log_terms = fk_mat - u_ik_.array() + log_c;

    // LSE over replicas (axis=1) for each state i
    Eigen::VectorXd log_sum_exp(n_states_);
    for (int i = 0; i < n_states_; ++i) {
        double M = log_terms.row(i).maxCoeff();
        if (!std::isfinite(M)) {
            // State i completely unvisited across all replicas → p_i = 0
            log_sum_exp(i) = -std::numeric_limits<double>::infinity();
        } else {
            double sum_e = (log_terms.row(i) - M).exp().sum();
            log_sum_exp(i) = M + std::log(sum_e + DHAM_LOG_ZERO_GUARD);
        }
    }

    // p_i ∝ exp(log_sum_exp)  (normalize after)
    // Shift by max to avoid numerical overflow in exp
    double lse_max = log_sum_exp.maxCoeff();
    if (!std::isfinite(lse_max)) {
        TitanUI::print_error("TitanDHAM::iterate_once: all states unvisited. "
                             "Check count matrices.");
        return 0.0;
    }
    Eigen::VectorXd p_new = (log_sum_exp.array() - lse_max).exp();

    // Apply target-temperature reweighting:
    // multiply each state by exp(-beta_target * <U_i>)
    // where <U_i> = mean over replicas weighted by c_out
    // Approximation: use u_ik_.rowwise().mean() / beta_mean as average U(i)
    // Exact approach: weight by c_out(i,k)
    for (int i = 0; i < n_states_; ++i) {
        double c_total = c_out_.row(i).sum();
        if (c_total < 1e-30) { p_new(i) = 0.0; continue; }
        // Weighted mean energy at target beta:
        // exp(sum_k c_out(i,k)/c_total * (beta_target - beta_k) * U_k(i))
        double log_reweight = 0.0;
        for (int k = 0; k < n_replicas_; ++k) {
            double w = c_out_(i, k) / c_total;
            if (w < 1e-15) continue;
            // u_ik_(i,k) = beta_k * U_k(i)  → U_k(i) = u_ik_(i,k) / beta_k
            double U_ki = (std::abs(betas_[static_cast<size_t>(k)]) > 1e-30)
                          ? u_ik_(i, k) / betas_[static_cast<size_t>(k)]
                          : 0.0;
            log_reweight += w * (beta_target_ - betas_[static_cast<size_t>(k)]) * U_ki;
        }
        p_new(i) *= std::exp(log_reweight);
    }

    // Normalize p_i
    double p_sum = p_new.sum();
    if (p_sum > 1e-300)
        p_new /= p_sum;
    else {
        TitanUI::print_warning(
            "TitanDHAM::iterate_once: p_i sum is ~0. Keeping previous p_i.");
        return 0.0;
    }
    p_i_ = p_new;

    // ------------------------------------------------------------
    // Step 2: Update f_k (free energy shifts)
    //
    // f_k(k) = -log( Σ_i p_i(i) * exp(-u_ik(i,k)) / Z )
    //        = -log( Σ_i p_i(i) * exp(-u_ik(i,k)) )  + const
    //
    // LSE: for each replica k,
    //   log_terms2(i,k) = log(p_i(i) + guard) - u_ik(i,k)
    //   f_k_new(k) = -LSE over i of log_terms2
    // Anchor: f_k[0] = 0  (arbitrary gauge)
    // ------------------------------------------------------------
    Eigen::ArrayXXd log_pi_mat =
        (p_i_.array() + DHAM_LOG_ZERO_GUARD).log().replicate(1, n_replicas_);
    Eigen::ArrayXXd log_terms2 = log_pi_mat - u_ik_.array();

    Eigen::VectorXd f_k_new(n_replicas_);
    for (int k = 0; k < n_replicas_; ++k) {
        double M2 = log_terms2.col(k).maxCoeff();
        if (!std::isfinite(M2)) {
            f_k_new(k) = 0.0;
        } else {
            double sum_e2 = (log_terms2.col(k) - M2).exp().sum();
            f_k_new(k) = -(M2 + std::log(sum_e2 + DHAM_LOG_ZERO_GUARD));
        }
    }

    // Anchor replica 0 to f=0
    double anchor = f_k_new(0);
    f_k_new.array() -= anchor;

    double max_delta = (f_k_new - f_k_).cwiseAbs().maxCoeff();
    f_k_ = f_k_new;
    return max_delta;
}

// ============================================================================
// solve
// ============================================================================
bool TitanDHAM::solve(int max_iter, double tol) {
    TitanUI::print_module_header("DHAM Self-Consistent Solver");

    // Validate that at least some counts exist
    double total_counts = 0.0;
    for (int k = 0; k < n_replicas_; ++k)
        total_counts += C_k_[static_cast<size_t>(k)].sum();
    if (total_counts < 1.0) {
        TitanUI::print_error(
            "TitanDHAM::solve: No counts accumulated. "
            "Call accumulate_counts() first.");
        return false;
    }

    // Build c_out from count matrices
    rebuild_c_out();

    TitanUI::print_info("  States: " + std::to_string(n_states_) +
                        "  Replicas: " + std::to_string(n_replicas_) +
                        "  Total counts: " + std::to_string(static_cast<long>(total_counts)));

    // Initialize f_k to 0, p_i to uniform
    f_k_.setZero();
    p_i_ = Eigen::VectorXd::Constant(n_states_, 1.0 / n_states_);

    int iter = 0;
    double delta = std::numeric_limits<double>::max();

    for (iter = 0; iter < max_iter; ++iter) {
        delta = iterate_once();
        if (iter % 500 == 0) {
            TitanUI::print_info("  DHAM iter " + std::to_string(iter) +
                                "  max |Δf_k| = " + std::to_string(delta));
        }
        if (delta < tol) break;
    }

    solved_ = true;
    if (delta < tol) {
        TitanUI::print_info("  DHAM converged in " + std::to_string(iter) +
                            " iterations (|Δf_k| = " + std::to_string(delta) + ").");
        return true;
    } else {
        TitanUI::print_warning(
            "  DHAM did not converge in " + std::to_string(max_iter) +
            " iterations (final |Δf_k| = " + std::to_string(delta) +
            "). Results may be inaccurate. Consider increasing max_iter.");
        return false;
    }
}

// ============================================================================
// get_unbiased_transition_matrix
// ============================================================================
Eigen::MatrixXd TitanDHAM::get_unbiased_transition_matrix() const {
    if (!solved_) {
        TitanUI::print_error("TitanDHAM: solve() must be called first.");
        return Eigen::MatrixXd::Zero(n_states_, n_states_);
    }

    // Reweight count matrices:
    // C_unbiased(i,j) = Σ_k C_k(i,j) * w_k(i)
    // where w_k(i) = exp(f_k - u_ik(i,k) - f_ref(i))
    // f_ref(i) = LSE over k for normalization
    Eigen::MatrixXd C_unbiased = Eigen::MatrixXd::Zero(n_states_, n_states_);

    for (int i = 0; i < n_states_; ++i) {
        // log_w(k) = f_k(k) - u_ik(i,k) for replica k, unnormalized
        Eigen::VectorXd log_w(n_replicas_);
        for (int k = 0; k < n_replicas_; ++k)
            log_w(k) = f_k_(k) - u_ik_(i, k);

        double M = log_w.maxCoeff();
        if (!std::isfinite(M)) continue;
        Eigen::VectorXd w = (log_w.array() - M).exp();
        double w_sum = w.sum();
        if (w_sum < 1e-300) continue;
        w /= w_sum;

        for (int k = 0; k < n_replicas_; ++k) {
            if (w(k) < 1e-15) continue;
            C_unbiased.row(i) += w(k) * C_k_[static_cast<size_t>(k)].row(i);
        }
    }

    // Symmetrize and row-normalize → TPM
    Eigen::MatrixXd C_sym = 0.5 * (C_unbiased + C_unbiased.transpose());
    Eigen::MatrixXd T = Eigen::MatrixXd::Zero(n_states_, n_states_);
    for (int i = 0; i < n_states_; ++i) {
        double rs = C_sym.row(i).sum();
        if (rs < 1e-10)
            T.row(i) = Eigen::RowVectorXd::Constant(n_states_, 1.0 / n_states_);
        else
            T.row(i) = C_sym.row(i) / rs;
    }
    return T;
}

// ============================================================================
// get_free_energy_profile
// ============================================================================
Eigen::VectorXd TitanDHAM::get_free_energy_profile() const {
    if (!solved_) return Eigen::VectorXd::Zero(n_states_);

    // G(i) = -kT * log(p_i(i)),  shifted so min = 0
    const double kT = 1.0 / beta_target_;  // kJ/mol
    Eigen::VectorXd G(n_states_);

    // Find min non-zero p_i for reference
    double p_max = p_i_.maxCoeff();
    if (p_max < 1e-300) return Eigen::VectorXd::Zero(n_states_);

    for (int i = 0; i < n_states_; ++i) {
        double p = p_i_(i);
        if (p < 1e-300)
            G(i) = std::numeric_limits<double>::infinity();
        else
            G(i) = -kT * std::log(p / p_max);  // relative to most populated state
    }
    return G;
}

// ============================================================================
// write_results
// ============================================================================
void TitanDHAM::write_results(const std::string& pi_file,
                                const std::string& fes_file,
                                const std::string& tpm_file,
                                std::fstream& log) const {
    if (!solved_) {
        TitanUI::print_warning("TitanDHAM::write_results: solve() was not called.");
        return;
    }

    // Stationary distribution
    {
        std::ofstream f(pi_file);
        if (f.is_open()) {
            f << "# DHAM Reweighted Stationary Distribution at Target Temperature\n";
            f << "# state_id  p_i\n";
            for (int i = 0; i < n_states_; ++i)
                f << i << "  " << std::fixed << std::setprecision(8) << p_i_(i) << "\n";
            TitanUI::print_info("DHAM: stationary distribution written to " + pi_file);
            TitanUI::log_kv(log, "DHAM Pi File", pi_file);
        }
    }

    // Free energy profile
    {
        Eigen::VectorXd G = get_free_energy_profile();
        std::ofstream f(fes_file);
        if (f.is_open()) {
            f << "# DHAM Reweighted Free Energy Profile (kJ/mol, min=0)\n";
            f << "# state_id  G_kJ_per_mol\n";
            for (int i = 0; i < n_states_; ++i) {
                f << i << "  " << std::fixed << std::setprecision(4);
                if (std::isinf(G(i)))
                    f << "inf";
                else
                    f << G(i);
                f << "\n";
            }
            TitanUI::print_info("DHAM: free energy profile written to " + fes_file);
            TitanUI::log_kv(log, "DHAM FES File", fes_file);
        }
    }

    // Unbiased transition matrix
    {
        Eigen::MatrixXd T = get_unbiased_transition_matrix();
        std::ofstream f(tpm_file);
        if (f.is_open()) {
            f << "# DHAM Unbiased Transition Matrix at Target Temperature\n";
            f << "# " << n_states_ << "x" << n_states_ << "\n";
            f << std::fixed << std::setprecision(8);
            for (int i = 0; i < n_states_; ++i) {
                for (int j = 0; j < n_states_; ++j) {
                    if (j > 0) f << " ";
                    f << T(i, j);
                }
                f << "\n";
            }
            TitanUI::print_info("DHAM: unbiased TPM written to " + tpm_file);
            TitanUI::log_kv(log, "DHAM TPM File", tpm_file);
        }
    }
}

// ============================================================================
// compute_state_energies (static)
// ============================================================================
Eigen::MatrixXd TitanDHAM::compute_state_energies(
    const std::unordered_map<int, std::vector<int>>&    per_replica_dtrajs,
    const std::unordered_map<int, std::vector<double>>& per_replica_energies,
    int n_states, int n_replicas)
{
    // result(i,k) = mean energy of state i in replica k  (kJ/mol)
    // Initialized to 0; states with no visits keep 0
    Eigen::MatrixXd result = Eigen::MatrixXd::Zero(n_states, n_replicas);
    Eigen::MatrixXi counts  = Eigen::MatrixXi::Zero(n_states, n_replicas);

    for (int k = 0; k < n_replicas; ++k) {
        auto it_d = per_replica_dtrajs.find(k);
        auto it_e = per_replica_energies.find(k);
        if (it_d == per_replica_dtrajs.end() ||
            it_e == per_replica_energies.end()) continue;

        const auto& dtraj   = it_d->second;
        const auto& energies = it_e->second;

        size_t n_frames = std::min(dtraj.size(), energies.size());
        for (size_t t = 0; t < n_frames; ++t) {
            double e = energies[t];
            if (std::isnan(e)) continue;  // Skip unmatched frames
            int s = dtraj[t];
            if (s < 0 || s >= n_states) continue;
            result(s, k) += e;
            counts(s, k) += 1;
        }
    }

    // Divide by counts; states with 0 visits: warn
    for (int i = 0; i < n_states; ++i) {
        for (int k = 0; k < n_replicas; ++k) {
            if (counts(i, k) > 0)
                result(i, k) /= static_cast<double>(counts(i, k));
            // else: stays 0 — caller should call set_state_energy only for visited states
        }
    }

    return result;
}
