#include "postproc/postproc_common.h"

using namespace postproc_detail;

// ============================================================================
// 6. DHAM post-processing
// ============================================================================
void run_dham_postproc(const arguments &args,
                       std::unique_ptr<TitanMSM> &msm_analyzer,
                       TitanExecutor &executor, info &info_output,
                       std::fstream &log,
                       std::vector<double> &out_frame_weights,
                       Eigen::MatrixXd *out_reweighted_tpm,
                       Eigen::VectorXd *out_reweighted_pi,
                       bool *out_converged,
                       TitanValidation::DhamRebuildContext
                           *out_validation_dham_context) {
  out_frame_weights.clear();
  if (out_reweighted_tpm)
    out_reweighted_tpm->resize(0, 0);
  if (out_reweighted_pi)
    out_reweighted_pi->resize(0);
  if (out_converged)
    *out_converged = false;
  if (out_validation_dham_context)
    *out_validation_dham_context = {};
  if (!args.dham.flag_dham || !msm_analyzer)
    return;

  auto start_dham = std::chrono::steady_clock::now();
  TitanUI::print_module_header("DHAM Post-Processing");
  const bool use_npt_observable = (args.dham.ensemble == "npt");
  if (args.lj.flag_lj) {
    TitanUI::log_kv(log, "DHAM Feature Domain", "LJ phase-2 feature space");
  } else {
    TitanUI::log_kv(log, "DHAM Feature Domain",
                    "protein network -> TICA -> MSM feature space");
  }
  TitanUI::log_kv(log, "DHAM Ensemble", args.dham.ensemble);
  if (use_npt_observable) {
    TitanUI::log_kv(log, "DHAM Thermodynamic Observable",
                    "enthalpy-like reweighting H = column1 + P*V(box)");
    TitanUI::log_kv(log, "DHAM Target Pressure (bar)",
                    std::to_string(args.dham.target_pressure_bar));
    TitanUI::print_info(
        "DHAM: NPT mode enabled. Using H = E_col1 + P*V with P = " +
        std::to_string(args.dham.target_pressure_bar) +
        " bar and V from trajectory box metadata.");
    TitanUI::print_warning(
        "DHAM NPT: for best thermodynamic consistency, the first numeric "
        "column in each .xvg should be configurational/potential energy, not "
        "a mixed observable.");
  } else {
    TitanUI::log_kv(log, "DHAM Thermodynamic Observable",
                    "energy-based NVT/REMD reweighting");
  }

  // 1. Read temperatures
  if (args.dham.temps_file.empty()) {
    TitanUI::print_error("DHAM: -dham-temps not set. Skipping DHAM.");
    return;
  }
  auto temps = TitanDHAM::read_temps(args.dham.temps_file);
  if (temps.empty()) {
    TitanUI::print_error("DHAM: No temperatures read. Skipping.");
    return;
  }

  // 2. Read energy file list
  const bool segment_audit = !args.dham.segment_map_file.empty();
  if (args.dham.ener_list_file.empty() && !segment_audit) {
    TitanUI::print_error("DHAM: -dham-ener not set. Skipping DHAM.");
    return;
  }
  auto ener_paths = segment_audit ? std::vector<std::string>{} : TitanDHAM::read_file_list(args.dham.ener_list_file);
  if (ener_paths.empty() && !segment_audit) {
    TitanUI::print_error("DHAM: No energy files read. Skipping.");
    return;
  }

  const int n_replicas = static_cast<int>(temps.size());
  if (!segment_audit && static_cast<int>(ener_paths.size()) != n_replicas) {
    TitanUI::print_error(
        "DHAM: #temps (" + std::to_string(n_replicas) +
        ") != #energy_files (" + std::to_string(ener_paths.size()) + ").");
    return;
  }

  // Validate n_replicas if explicitly set
  if (args.dham.n_replicas > 0 && args.dham.n_replicas != n_replicas) {
    TitanUI::print_warning(
        "DHAM: --dham-n-replicas=" + std::to_string(args.dham.n_replicas) +
        " != temps file count (" + std::to_string(n_replicas) +
        "). Using temps file count.");
  }

  TitanUI::print_info("DHAM: " + std::to_string(n_replicas) + " replicas, " +
                       "target T = " + std::to_string(args.dham.target_temp) + " K");
  // Fix C: Document the 1-file-per-replica assumption
  TitanUI::print_warning(
      "DHAM: Each input trajectory file is treated as one replica. "
      "Concatenation does not establish physical continuity; REMD inputs "
      "need audited segment/exchange metadata.");

  // 3. Extract per-replica data from MSM
  if (!msm_analyzer->is_clustered()) {
    TitanUI::print_error("DHAM: MSM K-means has not run. "
        "Ensure -msm flag is set and enough frames are collected.");
    return;
  }
  if (!segment_audit && !msm_analyzer->is_built()) {
    TitanUI::print_error(
        "DHAM: MSM transition model has not been built successfully.");
    return;
  }

  const int dham_lag = msm_analyzer->msm_lag_time();
  if (dham_lag <= 0) {
    TitanUI::print_error("DHAM: invalid MSM lag time. Aborting DHAM.");
    return;
  }

  auto per_replica_dtrajs     = msm_analyzer->get_per_replica_dtrajs();
  auto per_replica_timestamps = msm_analyzer->get_per_replica_timestamps();
  auto per_replica_box_volumes = use_npt_observable
      ? msm_analyzer->get_per_replica_box_volumes_nm3()
      : std::unordered_map<int, std::vector<double>>{};

  if (per_replica_dtrajs.empty()) {
    TitanUI::print_error("DHAM: No per-replica discrete trajectories. "
                          "Was replica_id passed to add_projected_frame?");
    return;
  }

  if (segment_audit) {
    try {
      const auto metadata = titan_dham_segments::Metadata::read(
          args.dham.segment_map_file, n_replicas);
      std::vector<titan_dham_segments::CountAudit> audits;
      for (int k = 0; k < n_replicas; ++k) {
        const auto &rep = metadata.replicas[k];
        if (std::abs(rep.temperature - temps[k]) > 1e-6 ||
            (use_npt_observable && std::abs(rep.pressure - args.dham.target_pressure_bar) > 1e-6))
          throw std::invalid_argument("DHAM segment metadata ensemble mismatch");
        const auto &states = per_replica_dtrajs.at(k);
        const auto &times = per_replica_timestamps.at(k);
        if (use_npt_observable) {
          const auto &volumes = per_replica_box_volumes.at(k);
          if (volumes.size() != times.size())
            throw std::invalid_argument("DHAM segment volume dimensions differ");
          for (size_t i = 0; i < times.size(); ++i) {
            const double expected = metadata.match(k, times[i]).volume;
            if (!std::isfinite(volumes[i]) ||
                std::abs(volumes[i] - expected) > std::max(1e-4, 2e-6 * expected))
              throw std::invalid_argument("DHAM segment trajectory volume mismatch");
          }
        }
        audits.push_back(titan_dham_segments::audit_counts(
            metadata, k, states, times, dham_lag, msm_analyzer->n_clusters()));
      }
      std::ofstream counts("dham_segment_counts.csv"), decisions("dham_segment_decisions.csv");
      if (!counts || !decisions) throw std::runtime_error("Cannot write DHAM segment audit");
      counts << "replica,state_from,state_to,candidate_count,physical_rates_validated\n";
      decisions << "replica,decision,windows\n";
      size_t eligible = 0;
      for (int k = 0; k < n_replicas; ++k) {
        for (const auto &[pair, n] : audits[k].counts)
          counts << k << ',' << pair.first << ',' << pair.second << ',' << n << ",false\n";
        for (const auto &[decision, n] : audits[k].decisions) {
          decisions << k << ',' << titan_dham_segments::name(decision) << ',' << n << '\n';
          if (decision == titan_dham_segments::Decision::eligible) eligible += n;
        }
      }
      counts.flush();decisions.flush();
      if (!counts || !decisions) throw std::runtime_error("Incomplete DHAM segment audit output");
      TitanUI::log_kv(log, "DHAM Segment Audit Eligible", std::to_string(eligible));
      TitanUI::log_kv(log, "DHAM Mode", "segment/exchange count audit; estimator not run");
      TitanUI::print_info("DHAM segment/exchange candidate counts written. Exact endpoint energies required; no nearest substitution.");
      TitanUI::print_warning("DHAM segment audit does not correct exchange censoring. No FES, rates, weights or bootstrap validation are produced by this mode.");
    } catch (const std::exception &error) {
      TitanUI::print_error(std::string("DHAM segment audit rejected: ") + error.what());
    }
    return;
  }

  // 4. Parse + align energy .xvg files with optional per-replica
  // window matching against the common trajectory/energy overlap.
  std::unordered_map<int, std::vector<double>> per_replica_energies;
  std::unordered_map<int, std::vector<std::pair<double, double>>> parsed_xvg;
  std::mutex parsed_xvg_mutex;
  const int n_io_threads =
      std::min(static_cast<int>(executor.capacity()), n_replicas);
  const int chunk = (n_replicas + n_io_threads - 1) / n_io_threads;
  executor.parallel_do(static_cast<size_t>(n_io_threads), [&](size_t t_sz) {
    const int t = static_cast<int>(t_sz);
    const int ks = t * chunk;
    const int ke = std::min(ks + chunk, n_replicas);
    if (ks >= ke)
      return;
    for (int k = ks; k < ke; ++k) {
      auto xvg_data = TitanDHAM::parse_xvg(ener_paths[static_cast<size_t>(k)]);
      if (!xvg_data.empty())
        xvg_data = filter_postproc_energy_pairs(xvg_data, args);
      std::lock_guard<std::mutex> lock(parsed_xvg_mutex);
      parsed_xvg[k] = std::move(xvg_data);
    }
  });
  TitanUI::log_kv(log, "DHAM Window Matching", args.dham.window_match);
  if (args.dham.dt_energy_ps > 0.0) {
    TitanUI::print_info("DHAM: Using user-specified dt_energy = " +
                        std::to_string(args.dham.dt_energy_ps) + " ps.");
  } else {
    TitanUI::print_info(
        "DHAM: Auto-detecting dt_energy independently for each energy file.");
  }
  if (args.dham.window_match == "off") {
    TitanUI::print_info(
        "DHAM: Window matching disabled. TITAN will align against the full "
        "trajectory time axis.");
  } else if (args.dham.window_match == "strict") {
    TitanUI::print_info(
        "DHAM: Strict window matching enabled. TITAN will abort if any "
        "replica has a trajectory/energy mismatch.");
  } else {
    TitanUI::print_info(
        "DHAM: Auto window matching enabled. TITAN will trim each replica to "
        "the common trajectory/energy overlap before DHAM.");
  }

  auto subset_by_indices = [](const auto &src,
                              const std::vector<std::size_t> &indices) {
    using value_type = typename std::decay_t<decltype(src)>::value_type;
    std::vector<value_type> out;
    out.reserve(indices.size());
    for (std::size_t idx : indices) {
      if (idx < src.size())
        out.push_back(src[idx]);
    }
    return out;
  };
  auto format_ps = [](double value) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(3) << value;
    return oss.str();
  };

  std::size_t window_total_original_frames = 0;
  std::size_t window_total_kept_frames = 0;
  std::size_t window_total_dropped_frames = 0;
  std::size_t window_total_unmatched_frames = 0;
  std::size_t window_trimmed_replicas = 0;
  std::size_t window_replicas_without_overlap = 0;
  std::size_t window_missing_energy_files = 0;
  std::size_t window_truncated_start = 0;
  std::size_t window_truncated_end = 0;
  std::size_t window_irregular_dt = 0;
  std::size_t window_sparse_energy = 0;
  bool strict_window_violation = false;
  std::vector<DhamWindowAuditRow> window_audit_rows;
  window_audit_rows.reserve(static_cast<std::size_t>(n_replicas));

  for (int k = 0; k < n_replicas; ++k) {
    auto it_xvg = parsed_xvg.find(k);
    auto it_dtraj = per_replica_dtrajs.find(k);
    auto it_ts = per_replica_timestamps.find(k);
    DhamWindowAuditRow audit;
    audit.replica_id = k;
    audit.traj_present =
        (it_dtraj != per_replica_dtrajs.end() && it_ts != per_replica_timestamps.end());
    audit.energy_present =
        (it_xvg != parsed_xvg.end() && !it_xvg->second.empty());
    if (audit.traj_present) {
      audit.traj_frames = it_dtraj->second.size();
      if (!it_ts->second.empty()) {
        audit.traj_min_ps = it_ts->second.front();
        audit.traj_max_ps = it_ts->second.back();
      }
    }
    if (it_xvg != parsed_xvg.end()) {
      audit.energy_points = it_xvg->second.size();
      if (!it_xvg->second.empty()) {
        audit.energy_min_ps = it_xvg->second.front().first;
        audit.energy_max_ps = it_xvg->second.back().first;
      }
    }
    if (it_xvg == parsed_xvg.end() || it_xvg->second.empty() ||
        it_dtraj == per_replica_dtrajs.end() ||
        it_ts == per_replica_timestamps.end()) {
      if (!audit.energy_present)
        ++window_missing_energy_files;
      audit.status = !audit.traj_present ? "missing_trajectory" :
                     (!audit.energy_present ? "missing_or_empty_energy" :
                                              "incomplete_inputs");
      window_audit_rows.push_back(std::move(audit));
      continue;
    }

    const double dt_energy_ps =
        (args.dham.dt_energy_ps > 0.0) ? args.dham.dt_energy_ps : 0.0;
    if (args.dham.window_match == "off") {
      const auto matched = TitanDHAM::window_match_energy_to_dtraj(
          it_ts->second, it_xvg->second, dt_energy_ps);
      audit.inferred_dt_energy_ps = matched.inferred_dt_energy_ps;
      audit.overlap_found = matched.overlap_found;
      audit.overlap_start_ps = matched.overlap_start_ps;
      audit.overlap_end_ps = matched.overlap_end_ps;
      audit.kept_frames = matched.kept_frames;
      audit.dropped_outside_overlap = matched.dropped_outside_overlap;
      audit.unmatched_in_overlap = matched.unmatched_in_overlap;
      audit.sparse_energy = (audit.energy_points < 2);
      audit.irregular_energy_dt =
          detect_irregular_energy_dt(it_xvg->second, matched.inferred_dt_energy_ps);
      per_replica_energies[k] = TitanDHAM::align_energy_to_dtraj(
          it_ts->second, it_xvg->second, dt_energy_ps);
      audit.status = "legacy_off";
      if (audit.irregular_energy_dt)
        ++window_irregular_dt;
      if (audit.sparse_energy)
        ++window_sparse_energy;
      window_audit_rows.push_back(std::move(audit));
      continue;
    }

    const auto matched = TitanDHAM::window_match_energy_to_dtraj(
        it_ts->second, it_xvg->second, dt_energy_ps);
    audit.inferred_dt_energy_ps = matched.inferred_dt_energy_ps;
    audit.overlap_found = matched.overlap_found;
    audit.overlap_start_ps = matched.overlap_start_ps;
    audit.overlap_end_ps = matched.overlap_end_ps;
    audit.kept_frames = matched.kept_frames;
    audit.dropped_outside_overlap = matched.dropped_outside_overlap;
    audit.unmatched_in_overlap = matched.unmatched_in_overlap;
    const double coverage_tol =
        std::max(1e-6, 2.0 * std::max(1e-6, matched.inferred_dt_energy_ps));
    audit.truncated_start =
        std::isfinite(matched.traj_min_ps) && std::isfinite(matched.energy_min_ps) &&
        matched.energy_min_ps > matched.traj_min_ps + coverage_tol;
    audit.truncated_end =
        std::isfinite(matched.traj_max_ps) && std::isfinite(matched.energy_max_ps) &&
        matched.energy_max_ps < matched.traj_max_ps - coverage_tol;
    audit.sparse_energy = (audit.energy_points < 2);
    audit.irregular_energy_dt =
        detect_irregular_energy_dt(it_xvg->second, matched.inferred_dt_energy_ps);
    window_total_original_frames += matched.original_frames;
    window_total_kept_frames += matched.kept_frames;
    window_total_dropped_frames += matched.dropped_outside_overlap;
    window_total_unmatched_frames += matched.unmatched_in_overlap;
    if (audit.truncated_start)
      ++window_truncated_start;
    if (audit.truncated_end)
      ++window_truncated_end;
    if (audit.irregular_energy_dt)
      ++window_irregular_dt;
    if (audit.sparse_energy)
      ++window_sparse_energy;

    const bool no_usable_overlap =
        !matched.overlap_found || matched.kept_frames == 0;
    const bool has_window_loss =
        matched.dropped_outside_overlap > 0 || matched.unmatched_in_overlap > 0;

    if (no_usable_overlap) {
      ++window_replicas_without_overlap;
      const std::string msg =
          "DHAM: Replica " + std::to_string(k) +
          " has no usable overlap between trajectory [" +
          format_ps(matched.traj_min_ps) + ", " + format_ps(matched.traj_max_ps) +
          "] ps and energy [" + format_ps(matched.energy_min_ps) + ", " +
          format_ps(matched.energy_max_ps) + "] ps.";
      if (args.dham.window_match == "strict") {
        TitanUI::print_warning(msg);
        strict_window_violation = true;
      } else {
        TitanUI::print_warning(msg + " Replica will be ignored by DHAM.");
      }
      audit.status = "no_overlap";
      window_audit_rows.push_back(std::move(audit));
      per_replica_dtrajs.erase(k);
      per_replica_timestamps.erase(k);
      per_replica_energies.erase(k);
      if (use_npt_observable)
        per_replica_box_volumes.erase(k);
      continue;
    }

    if (matched.dropped_outside_overlap > 0)
      ++window_trimmed_replicas;
    if (args.dham.window_match == "strict" && has_window_loss) {
      TitanUI::print_warning(
          "DHAM: Replica " + std::to_string(k) +
          " violated strict window matching. kept=" +
          std::to_string(matched.kept_frames) + "/" +
          std::to_string(matched.original_frames) + ", dropped_outside=" +
          std::to_string(matched.dropped_outside_overlap) +
          ", unmatched_inside=" +
          std::to_string(matched.unmatched_in_overlap) + ".");
      strict_window_violation = true;
      audit.status = "strict_rejected";
      window_audit_rows.push_back(std::move(audit));
      continue;
    }

    if (has_window_loss) {
      TitanUI::print_info(
          "DHAM: Replica " + std::to_string(k) + " window-matched to [" +
          format_ps(matched.overlap_start_ps) + ", " +
          format_ps(matched.overlap_end_ps) + "] ps; kept " +
          std::to_string(matched.kept_frames) + "/" +
          std::to_string(matched.original_frames) + " frames, dropped " +
          std::to_string(matched.dropped_outside_overlap) + ", unmatched " +
          std::to_string(matched.unmatched_in_overlap) + ".");
    }

    if (audit.truncated_start && audit.truncated_end) {
      audit.status = "truncated_both";
    } else if (audit.truncated_start) {
      audit.status = "truncated_start";
    } else if (audit.truncated_end) {
      audit.status = "truncated_end";
    } else if (audit.irregular_energy_dt) {
      audit.status = "irregular_dt";
    } else if (audit.sparse_energy) {
      audit.status = "sparse_energy";
    } else if (has_window_loss) {
      audit.status = "trimmed_to_overlap";
    } else {
      audit.status = "ok";
    }

    per_replica_dtrajs[k] =
        subset_by_indices(it_dtraj->second, matched.kept_indices);
    per_replica_timestamps[k] = matched.kept_dtraj_times_ps;
    per_replica_energies[k] = matched.aligned_values;
    if (use_npt_observable) {
      auto it_volume = per_replica_box_volumes.find(k);
      if (it_volume != per_replica_box_volumes.end()) {
        per_replica_box_volumes[k] =
            subset_by_indices(it_volume->second, matched.kept_indices);
      }
    }
    window_audit_rows.push_back(std::move(audit));
  }

  if (args.dham.window_match != "off") {
    TitanUI::log_kv(log, "DHAM Window Frames Kept",
                    std::to_string(window_total_kept_frames) + "/" +
                        std::to_string(window_total_original_frames));
    TitanUI::log_kv(log, "DHAM Window Frames Dropped",
                    std::to_string(window_total_dropped_frames));
    TitanUI::log_kv(log, "DHAM Window Frames Unmatched",
                    std::to_string(window_total_unmatched_frames));
    TitanUI::log_kv(log, "DHAM Window Replicas Trimmed",
                    std::to_string(window_trimmed_replicas));
    TitanUI::log_kv(log, "DHAM Window Replicas No Overlap",
                    std::to_string(window_replicas_without_overlap));
    if (window_total_dropped_frames > 0 || window_total_unmatched_frames > 0) {
      TitanUI::print_info(
        "DHAM: Window matching kept " +
          std::to_string(window_total_kept_frames) + "/" +
          std::to_string(window_total_original_frames) +
          " trajectory frames across replicas.");
    }
  }
  TitanUI::log_kv(log, "DHAM Window Energy Missing",
                  std::to_string(window_missing_energy_files));
  TitanUI::log_kv(log, "DHAM Window Truncated Start",
                  std::to_string(window_truncated_start));
  TitanUI::log_kv(log, "DHAM Window Truncated End",
                  std::to_string(window_truncated_end));
  TitanUI::log_kv(log, "DHAM Window Irregular dt",
                  std::to_string(window_irregular_dt));
  TitanUI::log_kv(log, "DHAM Window Sparse Energy",
                  std::to_string(window_sparse_energy));
  write_dham_window_audit_file(info_output.get_dham_window_audit_filename(),
                               window_audit_rows);
  TitanUI::print_info("DHAM: Exported window audit to " +
                      info_output.get_dham_window_audit_filename());
  if (strict_window_violation) {
    TitanUI::log_kv(log, "DHAM Converged", "no");
    TitanUI::log_kv(log, "DHAM Abort Reason",
                    "strict window matching rejected one or more replicas");
    TitanUI::print_warning(
        "DHAM: strict window matching failed. Skipping DHAM reweighting.");
    return;
  }
  std::size_t post_match_candidate_transitions = 0;
  for (const auto &[replica_id, dtraj] : per_replica_dtrajs) {
    (void)replica_id;
    if (dtraj.size() > static_cast<std::size_t>(dham_lag)) {
      post_match_candidate_transitions +=
          dtraj.size() - static_cast<std::size_t>(dham_lag);
    }
  }
  TitanUI::log_kv(log, "DHAM Window Candidate Transitions",
                  std::to_string(post_match_candidate_transitions));
  if (post_match_candidate_transitions == 0) {
    TitanUI::log_kv(log, "DHAM Converged", "no");
    TitanUI::log_kv(log, "DHAM Abort Reason",
                    "zero lag-valid transitions after window matching");
    TitanUI::print_warning(
        "DHAM: Window-matched data contain no lag-valid transitions for the "
        "current MSM lag. Merge matching trajectory/energy windows, reduce "
        "-msm-lag, or disable automatic trimming with -dham-window-match off.");
    return;
  }

  std::unordered_map<int, std::vector<double>> per_replica_observables;
  if (use_npt_observable) {
    const double pv_factor =
        args.dham.target_pressure_bar * kBarNm3ToKJPerMol;
    size_t replicas_missing_volume = 0;
    size_t volume_size_mismatches = 0;
    size_t frames_missing_valid_volume = 0;
    size_t frames_combined = 0;

    for (int k = 0; k < n_replicas; ++k) {
      auto it_energy = per_replica_energies.find(k);
      if (it_energy == per_replica_energies.end())
        continue;

      auto it_volume = per_replica_box_volumes.find(k);
      if (it_volume == per_replica_box_volumes.end()) {
        ++replicas_missing_volume;
        continue;
      }

      const auto &energies = it_energy->second;
      const auto &volumes_nm3 = it_volume->second;
      if (volumes_nm3.size() != energies.size()) {
        ++volume_size_mismatches;
      }

      std::vector<double> observable(
          energies.size(), std::numeric_limits<double>::quiet_NaN());
      const size_t n_frames = std::min(energies.size(), volumes_nm3.size());
      for (size_t t = 0; t < n_frames; ++t) {
        const double e = energies[t];
        const double v = volumes_nm3[t];
        if (!std::isfinite(e) || !std::isfinite(v) || v < 0.0) {
          ++frames_missing_valid_volume;
          continue;
        }
        observable[t] = e + pv_factor * v;
        ++frames_combined;
      }
      per_replica_observables[k] = std::move(observable);
    }

    if (replicas_missing_volume > 0) {
      TitanUI::print_error(
          "DHAM NPT: Missing per-frame box volumes for " +
          std::to_string(replicas_missing_volume) + "/" +
          std::to_string(n_replicas) +
          " replicas. Re-run with the current TITAN build so TICA/MSM caches "
          "trajectory volume metadata.");
      return;
    }
    if (frames_combined == 0) {
      TitanUI::print_error(
          "DHAM NPT: No valid E(t) + P*V(t) samples were constructed from the "
          "aligned energy traces and trajectory box metadata.");
      return;
    }
    if (volume_size_mismatches > 0) {
      TitanUI::print_warning(
          "DHAM NPT: " + std::to_string(volume_size_mismatches) +
          " replica volume series had a different length than the aligned "
          "energy series; truncation to the common prefix was applied.");
    }
    if (frames_missing_valid_volume > 0) {
      TitanUI::print_warning(
          "DHAM NPT: " + std::to_string(frames_missing_valid_volume) +
          " aligned frames had invalid energy/volume data and were excluded "
          "from the enthalpy-like observable.");
    }
    TitanUI::print_info(
        "DHAM NPT: Constructed " + std::to_string(frames_combined) +
        " frame-level H = E_col1 + P*V samples.");
  } else {
    per_replica_observables = per_replica_energies;
  }
  if (per_replica_observables.empty()) {
    TitanUI::log_kv(log, "DHAM Converged", "no");
    TitanUI::log_kv(log, "DHAM Abort Reason",
                    "no usable aligned observables remained after window matching");
    TitanUI::print_warning(
        "DHAM: No usable aligned thermodynamic observables remain after "
        "window matching.");
    return;
  }

  if (out_validation_dham_context) {
    out_validation_dham_context->available = true;
    out_validation_dham_context->n_states = msm_analyzer->n_clusters();
    out_validation_dham_context->lag_frames = dham_lag;
    out_validation_dham_context->n_replicas = n_replicas;
    out_validation_dham_context->temps_k = temps;
    out_validation_dham_context->target_temp_k = args.dham.target_temp;
    out_validation_dham_context->max_iter = args.dham.max_iter;
    out_validation_dham_context->tolerance = args.dham.tolerance;
    out_validation_dham_context->ensemble = args.dham.ensemble;
    out_validation_dham_context->per_replica_dtrajs = per_replica_dtrajs;
    out_validation_dham_context->per_replica_timestamps_ps =
        per_replica_timestamps;
    out_validation_dham_context->per_replica_observables =
        per_replica_observables;
  }

  // 5. Determine number of microstates from MSM
  const int n_states = msm_analyzer->n_clusters();

  // 6. Compute per-state mean thermodynamic observables (Opt 4: multi-threaded)
  auto state_energies = TitanDHAM::compute_state_energies(
      per_replica_dtrajs, per_replica_observables, n_states, n_replicas,
      static_cast<int>(executor.capacity()));

  // 7. Create TitanDHAM and load data
  TitanDHAM dham(n_states, n_replicas, temps, args.dham.target_temp);
  std::vector<std::vector<unsigned char>> valid_energy(
      static_cast<size_t>(n_replicas),
      std::vector<unsigned char>(static_cast<size_t>(n_states), 0));
  size_t valid_state_energy_pairs = 0;

  // Load state energies (skip NaN = unvisited state/replica pairs)
  for (int k = 0; k < n_replicas; ++k) {
    for (int i = 0; i < n_states; ++i) {
      double e = state_energies(i, k);
      if (std::isfinite(e)) {
        valid_energy[static_cast<size_t>(k)][static_cast<size_t>(i)] = 1;
        ++valid_state_energy_pairs;
        dham.set_state_energy(k, i, e);
      }
    }
  }
  const size_t total_state_energy_pairs =
      static_cast<size_t>(n_states) * static_cast<size_t>(n_replicas);
  const size_t missing_state_energy_pairs =
      total_state_energy_pairs - valid_state_energy_pairs;
  if (missing_state_energy_pairs > 0) {
    TitanUI::print_warning(
        "DHAM: Missing aligned thermodynamic observables for " +
        std::to_string(missing_state_energy_pairs) + "/" +
        std::to_string(total_state_energy_pairs) +
        " state/replica pairs. Transitions touching those states will be skipped.");
  }

  // Build count matrices from per-replica dtrajs (use same lag as MSM)
  // Skip transitions that cross trajectory breaks or states with missing energy.
  TitanUI::print_info("DHAM: Using MSM lag = " + std::to_string(dham_lag) +
                        " frames for count matrix construction.");
  size_t total_counts = 0;
  size_t skipped_breaks = 0;
  size_t skipped_missing_energy = 0;
  size_t skipped_invalid_states = 0;
  for (int k = 0; k < n_replicas; ++k) {
    auto it = per_replica_dtrajs.find(k);
    if (it == per_replica_dtrajs.end()) continue;
    const auto &dtraj = it->second;
    auto it_ts = per_replica_timestamps.find(k);
    const bool have_times = (it_ts != per_replica_timestamps.end() &&
                              it_ts->second.size() == dtraj.size());
    // Estimate expected gap = lag * dt_frame
    double dt_frame = 1.0, expected_gap = dham_lag;
    if (have_times && it_ts->second.size() > 1) {
      dt_frame = it_ts->second[1] - it_ts->second[0];
      expected_gap = dham_lag * dt_frame;
    }
    for (size_t t = 0; t + static_cast<size_t>(dham_lag) < dtraj.size(); ++t) {
      if (have_times) {
        double actual_gap = it_ts->second[t + static_cast<size_t>(dham_lag)]
                          - it_ts->second[t];
        if (std::abs(actual_gap - expected_gap) > 0.1 * std::abs(expected_gap)) {
          ++skipped_breaks;
          continue; // trajectory break — skip phantom transition
        }
      }
      const int state_i = dtraj[t];
      const int state_j = dtraj[t + static_cast<size_t>(dham_lag)];
      if (state_i < 0 || state_i >= n_states || state_j < 0 || state_j >= n_states) {
        ++skipped_invalid_states;
        continue;
      }
      if (!valid_energy[static_cast<size_t>(k)][static_cast<size_t>(state_i)] ||
          !valid_energy[static_cast<size_t>(k)][static_cast<size_t>(state_j)]) {
        ++skipped_missing_energy;
        continue;
      }
      dham.accumulate_counts(k, state_i, state_j);
      ++total_counts;
    }
  }
  if (skipped_breaks > 0) {
    TitanUI::print_warning("DHAM: Skipped " + std::to_string(skipped_breaks) +
        " transitions across trajectory breaks (kept " +
        std::to_string(total_counts) + ").");
  }
  if (skipped_missing_energy > 0) {
    TitanUI::print_warning("DHAM: Skipped " +
        std::to_string(skipped_missing_energy) +
        " transitions involving states without aligned thermodynamic observables.");
  }
  if (skipped_invalid_states > 0) {
    TitanUI::print_warning("DHAM: Skipped " +
        std::to_string(skipped_invalid_states) +
        " transitions with invalid microstate indices.");
  }

  // 8. Solve DHAM
  bool converged = dham.solve(args.dham.max_iter, args.dham.tolerance);
  TitanUI::log_kv(log, "DHAM Converged", converged ? "yes" : "no");
  if (out_converged)
    *out_converged = converged;
  if (converged) {
    if (out_reweighted_pi)
      *out_reweighted_pi = dham.get_stationary_distribution();
    if (out_reweighted_tpm)
      *out_reweighted_tpm = dham.get_unbiased_transition_matrix();
    auto per_replica_weights = dham.get_frame_weights(per_replica_dtrajs);
    const auto &frame_replica_ids = msm_analyzer->get_frame_replica_ids();
    const auto &frame_global_indices = msm_analyzer->get_frame_global_indices();

    if (frame_replica_ids.size() != frame_global_indices.size()) {
      TitanUI::print_warning(
          "DHAM: projected frame metadata mismatch: replica_ids=" +
          std::to_string(frame_replica_ids.size()) + ", global_indices=" +
          std::to_string(frame_global_indices.size()) +
          ". Some frame weights may be dropped.");
    }

    int64_t max_global_index = -1;
    for (int64_t global_idx : frame_global_indices) {
      if (global_idx > max_global_index)
        max_global_index = global_idx;
    }
    if (max_global_index >= 0) {
      out_frame_weights.assign(static_cast<size_t>(max_global_index) + 1, 0.0);
    }

    std::unordered_map<int, size_t> per_replica_offsets;
    per_replica_offsets.reserve(per_replica_weights.size());
    std::vector<unsigned char> assigned(out_frame_weights.size(), 0);
    const size_t tracked_projected_frames =
        std::min(frame_replica_ids.size(), frame_global_indices.size());
    size_t mapped_frames = 0;
    size_t dropped_missing_replica_weights = 0;
    size_t dropped_invalid_global_index = 0;
    size_t duplicate_global_indices = 0;

    for (size_t frame_idx = 0; frame_idx < tracked_projected_frames; ++frame_idx) {
      const int replica_id = frame_replica_ids[frame_idx];
      auto it_weights = per_replica_weights.find(replica_id);
      size_t &replica_offset = per_replica_offsets[replica_id];
      if (it_weights == per_replica_weights.end() ||
          replica_offset >= it_weights->second.size()) {
        ++dropped_missing_replica_weights;
        continue;
      }

      const int64_t global_idx = frame_global_indices[frame_idx];
      if (global_idx < 0 ||
          static_cast<size_t>(global_idx) >= out_frame_weights.size()) {
        ++dropped_invalid_global_index;
        ++replica_offset;
        continue;
      }

      const double w = it_weights->second[replica_offset++];
      const size_t global_slot = static_cast<size_t>(global_idx);
      if (assigned[global_slot]) {
        ++duplicate_global_indices;
      }
      out_frame_weights[global_slot] = w;
      assigned[global_slot] = 1;
      ++mapped_frames;
    }

    double weight_sum = 0.0;
    for (double w : out_frame_weights)
      weight_sum += w;
    TitanUI::print_info(
        "DHAM: Computed " + std::to_string(mapped_frames) +
        " global frame weights for reweighted FES.");
    TitanUI::print_info("DHAM: Weight sum = " + std::to_string(weight_sum) +
                        " (expected ~1.0)");
    if (dropped_missing_replica_weights > 0) {
      TitanUI::print_warning(
          "DHAM: Dropped " +
          std::to_string(dropped_missing_replica_weights) +
          " projected frames while rebuilding global weights.");
    }
    if (dropped_invalid_global_index > 0) {
      TitanUI::print_warning(
          "DHAM: Dropped " + std::to_string(dropped_invalid_global_index) +
          " projected frames with invalid global indices.");
    }
    if (duplicate_global_indices > 0) {
      TitanUI::print_warning(
          "DHAM: Encountered " + std::to_string(duplicate_global_indices) +
          " duplicate global frame indices while assigning weights.");
    }

    // 9. Export results
    dham.write_results(info_output.get_dham_pi_filename(),
                       info_output.get_dham_fes_filename(),
                       info_output.get_dham_tpm_filename(), log,
                       args.eu_par);
  } else {
    TitanUI::print_warning(
        "DHAM: solve() did not converge - results were not written and "
        "PCA/TICA FES will use uniform weights. Increase --dham-max-iter "
        "or check the aligned energy data.");
  }

  auto time_dham = std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::steady_clock::now() - start_dham);
  TitanUI::log_time(log, "DHAM Analysis", time_dham.count() / 1000.0);
}
