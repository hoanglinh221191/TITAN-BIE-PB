// ============================================================================
// TitanDHAM — Dynamic Histogram Analysis Method (Implementation)
// ============================================================================
#include "titan_dham.h"
#include <algorithm>
#include <iomanip>
#include <numeric>
#include <stdexcept>
#include <thread>

namespace {

std::uint64_t pack_transition_key(int state_i, int state_j, int n_states) noexcept {
    return static_cast<std::uint64_t>(state_i) *
               static_cast<std::uint64_t>(n_states) +
           static_cast<std::uint64_t>(state_j);
}

int unpack_transition_source(std::uint64_t key, int n_states) noexcept {
    return static_cast<int>(key / static_cast<std::uint64_t>(n_states));
}

int unpack_transition_dest(std::uint64_t key, int n_states) noexcept {
    return static_cast<int>(key % static_cast<std::uint64_t>(n_states));
}

} // namespace

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

    C_k_sparse_.resize(n_replicas);

    c_out_  = MatrixXdRM::Zero(n_states, n_replicas);
    u_ik_   = MatrixXdRM::Constant(
        n_states, n_replicas, std::numeric_limits<double>::infinity());
    log_c_  = MatrixXdRM::Zero(n_states, n_replicas);

    // Initialize f_k to zero (first replica anchored during solve)
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

    const std::uint64_t key =
        pack_transition_key(state_i, state_j, n_states_);
    auto& counts = C_k_sparse_[static_cast<size_t>(replica_k)];
    counts[key] += 1.0;
    c_out_(state_i, replica_k) += 1.0;
    solved_ = false;
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
    solved_ = false;
}

// ============================================================================
// rebuild_c_out — row sums of sparse count maps
// ============================================================================
void TitanDHAM::rebuild_c_out() {
    c_out_.setZero();
    for (int k = 0; k < n_replicas_; ++k) {
        for (const auto& entry : C_k_sparse_[static_cast<size_t>(k)]) {
            const int state_i = unpack_transition_source(entry.first, n_states_);
            if (state_i >= 0 && state_i < n_states_)
                c_out_(state_i, k) += entry.second;
        }
    }
}

// ============================================================================
// iterate_once — one DHAM self-consistent step using LSE trick
// Returns max |Δf_k| for convergence check.
// ============================================================================
double TitanDHAM::iterate_once() {
    // ------------------------------------------------------------
    // Step 1: Update p_i (stationary distribution at target T)
    // ------------------------------------------------------------

    Eigen::VectorXd log_p_new(n_states_);
    for (int i = 0; i < n_states_; ++i) {
        double max_value = -std::numeric_limits<double>::infinity();
        for (int k = 0; k < n_replicas_; ++k) {
            if (!std::isfinite(u_ik_(i, k)) || c_out_(i, k) <= 0.0)
                continue;
            const double value = f_k_(k) - u_ik_(i, k) + log_c_(i, k);
            if (value > max_value)
                max_value = value;
        }

        if (!std::isfinite(max_value)) {
            log_p_new(i) = -std::numeric_limits<double>::infinity();
        } else {
            double sum_e = 0.0;
            for (int k = 0; k < n_replicas_; ++k) {
                if (!std::isfinite(u_ik_(i, k)) || c_out_(i, k) <= 0.0)
                    continue;
                const double value = f_k_(k) - u_ik_(i, k) + log_c_(i, k);
                sum_e += std::exp(value - max_value);
            }
            log_p_new(i) = max_value + std::log(sum_e + DHAM_LOG_ZERO_GUARD);
        }
    }


    // Apply target-temperature reweighting in log-space.
    for (int i = 0; i < n_states_; ++i) {
        const double c_total = c_out_.row(i).sum();
        if (c_total < 1e-30) {
            log_p_new(i) = -std::numeric_limits<double>::infinity();
            continue;
        }
        double log_reweight = 0.0;
        for (int k = 0; k < n_replicas_; ++k) {
            const double w = c_out_(i, k) / c_total;
            if (w < 1e-15 || !std::isfinite(u_ik_(i, k)))
                continue;
            const double beta_k = betas_[static_cast<size_t>(k)];
            const double U_ki = (std::abs(beta_k) > 1e-30)
                                ? u_ik_(i, k) / beta_k
                                : 0.0;
            log_reweight += w * (beta_target_ - beta_k) * U_ki;
        }
        log_p_new(i) += log_reweight;
    }

    // p_i proportional to exp(log_p_new) - shift by max to avoid overflow.
    const double max_log_p = log_p_new.maxCoeff();
    if (!std::isfinite(max_log_p)) {
        TitanUI::print_error("TitanDHAM::iterate_once: all states unvisited.");
        return std::numeric_limits<double>::infinity();
    }
    Eigen::VectorXd p_new = (log_p_new.array() - max_log_p).exp();

    // Normalize p_i
    double p_sum = p_new.sum();
    if (p_sum > 1e-300)
        p_new /= p_sum;
    else {
        TitanUI::print_warning(
            "TitanDHAM::iterate_once: p_i sum is ~0. Aborting solve.");
        return std::numeric_limits<double>::infinity();
    }
    p_i_ = p_new;

    // ------------------------------------------------------------
    // Step 2: Update f_k (free energy shifts)
    // f_k(k) = -log( Σ_i p_i(i) * exp(-u_ik(i,k)) )
    // Anchor: f_k[0] = 0
    // ------------------------------------------------------------
    const Eigen::VectorXd log_pi =
        (p_i_.array() + DHAM_LOG_ZERO_GUARD).log().matrix();
    Eigen::VectorXd f_k_new(n_replicas_);
    for (int k = 0; k < n_replicas_; ++k) {
        double max_value = -std::numeric_limits<double>::infinity();
        for (int i = 0; i < n_states_; ++i) {
            if (!std::isfinite(u_ik_(i, k)) || p_i_(i) <= 0.0)
                continue;
            const double value = log_pi(i) - u_ik_(i, k);
            if (value > max_value)
                max_value = value;
        }

        if (!std::isfinite(max_value)) {
            TitanUI::print_error(
                "TitanDHAM::iterate_once: encountered non-finite free-energy update.");
            return std::numeric_limits<double>::infinity();
        } else {
            double sum_e2 = 0.0;
            for (int i = 0; i < n_states_; ++i) {
                if (!std::isfinite(u_ik_(i, k)) || p_i_(i) <= 0.0)
                    continue;
                const double value = log_pi(i) - u_ik_(i, k);
                sum_e2 += std::exp(value - max_value);
            }
            f_k_new(k) = -(max_value + std::log(sum_e2 + DHAM_LOG_ZERO_GUARD));
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
    solved_ = false;

    // Rebuild c_out from sparse count maps as a consistency check.
    rebuild_c_out();
    const double total_counts = c_out_.sum();
    if (total_counts < 1.0) {
        TitanUI::print_error(
            "TitanDHAM::solve: No counts accumulated. "
            "Call accumulate_counts() first.");
        return false;
    }

    // Opt 1: Pre-compute log(c_out + guard) once — invariant during iterations
    log_c_ = (c_out_.array() + DHAM_LOG_ZERO_GUARD).log();

    TitanUI::print_info("  States: " + std::to_string(n_states_) +
                        "  Replicas: " + std::to_string(n_replicas_) +
                        "  Total counts: " + std::to_string(static_cast<long>(total_counts)));

    // Initialize f_k to 0, p_i to uniform
    f_k_.setZero();
    p_i_ = Eigen::VectorXd::Constant(n_states_, 1.0 / n_states_);

    int iter = 0;
    double delta = std::numeric_limits<double>::max();
    bool numerical_failure = false;

    for (iter = 0; iter < max_iter; ++iter) {
        delta = iterate_once();
        if (iter % 500 == 0) {
            TitanUI::print_info("  DHAM iter " + std::to_string(iter) +
                                "  max |Δf_k| = " + std::to_string(delta));
        }
        if (!std::isfinite(delta)) {
            numerical_failure = true;
            break;
        }
        if (delta < tol) break;
    }

    if (numerical_failure) {
        TitanUI::print_error(
            "  DHAM aborted due to numerical instability. Results were not accepted.");
        return false;
    }

    if (delta < tol) {
        solved_ = true;
        TitanUI::print_info("  DHAM converged in " + std::to_string(iter) +
                            " iterations (|Δf_k| = " + std::to_string(delta) + ").");
        return true;
    } else {
        TitanUI::print_warning(
            "  DHAM did not converge in " + std::to_string(max_iter) +
            " iterations (final |Δf_k| = " + std::to_string(delta) +
            "). Results may be inaccurate.");
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

    Eigen::MatrixXd C_unbiased = Eigen::MatrixXd::Zero(n_states_, n_states_);
    MatrixXdRM row_weights = MatrixXdRM::Zero(n_states_, n_replicas_);
    Eigen::VectorXd log_w(n_replicas_);

    for (int i = 0; i < n_states_; ++i) {
        for (int k = 0; k < n_replicas_; ++k)
            log_w(k) = f_k_(k) - u_ik_(i, k);

        double M = log_w.maxCoeff();
        if (!std::isfinite(M)) continue;
        Eigen::VectorXd w = (log_w.array() - M).exp();
        double w_sum = w.sum();
        if (w_sum < 1e-300) continue;
        w /= w_sum;
        row_weights.row(i) = w.transpose();
    }

    for (int k = 0; k < n_replicas_; ++k) {
        for (const auto& entry : C_k_sparse_[static_cast<size_t>(k)]) {
            const int state_i = unpack_transition_source(entry.first, n_states_);
            const double weight = row_weights(state_i, k);
            if (weight < 1e-15)
                continue;
            const int state_j = unpack_transition_dest(entry.first, n_states_);
            C_unbiased(state_i, state_j) += weight * entry.second;
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

    const double kT = 1.0 / beta_target_;  // kJ/mol
    Eigen::VectorXd G(n_states_);

    double p_max = p_i_.maxCoeff();
    if (p_max < 1e-300) return Eigen::VectorXd::Zero(n_states_);

    for (int i = 0; i < n_states_; ++i) {
        double p = p_i_(i);
        if (p < 1e-300)
            G(i) = std::numeric_limits<double>::infinity();
        else
            G(i) = -kT * std::log(p / p_max);
    }
    return G;
}

// ============================================================================
// get_frame_weights
// ============================================================================
std::unordered_map<int, std::vector<double>>
TitanDHAM::get_frame_weights(
    const std::unordered_map<int, std::vector<int>>& per_replica_dtrajs) const {
    if (!solved_) {
        TitanUI::print_error(
            "TitanDHAM::get_frame_weights: solve() must run first.");
        return {};
    }

    std::vector<double> N_i(static_cast<size_t>(n_states_), 0.0);
    for (const auto& entry : per_replica_dtrajs) {
        for (int s : entry.second) {
            if (s >= 0 && s < n_states_)
                N_i[static_cast<size_t>(s)] += 1.0;
        }
    }

    std::unordered_map<int, std::vector<double>> result;
    result.reserve(per_replica_dtrajs.size());
    for (const auto& entry : per_replica_dtrajs) {
        std::vector<double> weights(entry.second.size(), 0.0);
        for (size_t t = 0; t < entry.second.size(); ++t) {
            const int s = entry.second[t];
            if (s < 0 || s >= n_states_)
                continue;
            const double Ni = N_i[static_cast<size_t>(s)];
            if (Ni > 0.0)
                weights[t] = p_i_(s) / Ni;
        }
        result.emplace(entry.first, std::move(weights));
    }
    return result;
}

// ============================================================================
// get_flattened_frame_weights
// ============================================================================
std::vector<double> TitanDHAM::get_flattened_frame_weights(
    const std::unordered_map<int, std::vector<int>>& per_replica_dtrajs,
    int n_replicas) const {
    auto per_replica_weights = get_frame_weights(per_replica_dtrajs);

    size_t total = 0;
    for (int k = 0; k < n_replicas; ++k) {
        auto it = per_replica_dtrajs.find(k);
        if (it != per_replica_dtrajs.end())
            total += it->second.size();
    }

    std::vector<double> flat;
    flat.reserve(total);
    for (int k = 0; k < n_replicas; ++k) {
        auto it = per_replica_weights.find(k);
        if (it == per_replica_weights.end())
            continue;
        flat.insert(flat.end(), it->second.begin(), it->second.end());
    }
    return flat;
}

// ============================================================================
// write_results
// ============================================================================
void TitanDHAM::write_results(const std::string& pi_file,
                                 const std::string& fes_file,
                                const std::string& tpm_file,
                                std::fstream& log,
                                double eu_par) const {
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
        std::string unit_label = (eu_par > 1.1) ? "kcal/mol" : "kJ/mol";
        std::ofstream f(fes_file);
        if (f.is_open()) {
            f << "# DHAM Reweighted Free Energy Profile (" << unit_label << ", min=0)\n";
            f << "# state_id  G_" << unit_label << "\n";
            for (int i = 0; i < n_states_; ++i) {
                f << i << "  " << std::fixed << std::setprecision(4);
                if (std::isinf(G(i)))
                    f << "inf";
                else
                    f << G(i) / eu_par;
                f << "\n";
            }
            TitanUI::print_info("DHAM: free energy profile (" + unit_label +
                                ") written to " + fes_file);
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
    int n_states, int n_replicas, int num_threads)
{
    Eigen::MatrixXd result = Eigen::MatrixXd::Zero(n_states, n_replicas);
    Eigen::MatrixXi counts  = Eigen::MatrixXi::Zero(n_states, n_replicas);

    // Opt 4: Each thread processes disjoint columns — zero locks needed
    auto worker = [&](int k_start, int k_end) {
        for (int k = k_start; k < k_end; ++k) {
            auto it_d = per_replica_dtrajs.find(k);
            auto it_e = per_replica_energies.find(k);
            if (it_d == per_replica_dtrajs.end() ||
                it_e == per_replica_energies.end()) continue;

            const auto& dtraj   = it_d->second;
            const auto& energies = it_e->second;

            size_t n_frames = std::min(dtraj.size(), energies.size());
            for (size_t t = 0; t < n_frames; ++t) {
                double e = energies[t];
                if (std::isnan(e)) continue;
                int s = dtraj[t];
                if (s < 0 || s >= n_states) continue;
                result(s, k) += e;
                counts(s, k) += 1;
            }
        }
    };

    int n_use = std::min(num_threads, n_replicas);
    if (n_use <= 1) {
        worker(0, n_replicas);
    } else {
        std::vector<std::jthread> threads;
        threads.reserve(static_cast<size_t>(n_use));
        int chunk = (n_replicas + n_use - 1) / n_use;
        for (int t = 0; t < n_use; ++t) {
            int ks = t * chunk;
            int ke = std::min(ks + chunk, n_replicas);
            if (ks >= ke) continue;
            threads.emplace_back(worker, ks, ke);
        }
    }

    for (int i = 0; i < n_states; ++i) {
        for (int k = 0; k < n_replicas; ++k) {
            if (counts(i, k) > 0)
                result(i, k) /= static_cast<double>(counts(i, k));
            else
                result(i, k) = std::numeric_limits<double>::quiet_NaN();
        }
    }

    return result;
}
