#include "cas_analysis.h"
#include "gbsa_intrin.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <thread>

CASAnalyzer::CASAnalyzer() : global_total_frames(0), has_reference(false) {}

void CASAnalyzer::init(Multi_MOL &pro_in, Multi_MOL &lig,
                       const std::string &mol_2_type) {
  auto init_cas_data = [&](Multi_MOL &mol,
                           std::unordered_map<int, CAS_Result> &cas_data,
                           std::vector<int> &unique) {
    unique.clear();
    for (const auto &info : mol.info) {
      int res_id = info.get_resnr();
      if (cas_data.find(res_id) == cas_data.end()) {
        unique.emplace_back(res_id);
        CAS_Result new_res;
        new_res.res_id = res_id;
        new_res.res_name = info.get_res_name();
        cas_data[res_id] = new_res;
      }
    }
  };
  if (cas_receptor_data.empty())
    init_cas_data(pro_in, cas_receptor_data, unique_res_rec);
  if (cas_ligand_data.empty() &&
      mol_2_type.find("protein") != std::string::npos)
    init_cas_data(lig, cas_ligand_data, unique_res_lig);
}

void CASAnalyzer::process_frame(Multi_MOL &pro_in, Multi_MOL &lig,
                                size_t num_snap, size_t num_threads,
                                const std::vector<d_min> &dis_min,
                                cache_sigma_epsilon &CachedData,
                                const arguments &Arguments, info &info_output,
                                const AlignedDoubleVector &rec_free_radii,
                                const AlignedDoubleVector &lg_free_radii,
                                const AlignedDoubleVector &rec_complex_radii,
                                const AlignedDoubleVector &lg_complex_radii,
                                const GBSA_energy_snap &frame_gbsa_total) {
  global_total_frames++;
  if (!has_reference) {
    ref_pro_frame = pro_in.m_mol[num_snap];
    if (lig.m_mol.size() > 0) {
      ref_lig_frame = lig.m_mol[(lig.m_mol.size() > 1) ? num_snap : 0];
    }
    ref_pro_info = pro_in.info;
    ref_lig_info = lig.info;
    has_reference = true;
  }

  // ── Parallelized perform_cas_scan — Residue Shard Strategy ──
  auto perform_cas_scan = [&](Multi_MOL &mol,
                              std::unordered_map<int, CAS_Result> &cas_data,
                              const std::vector<int> &unique_res,
                              const std::string &type) {
    if (unique_res.empty())
      return;

    const size_t n_workers = num_threads;
    const size_t n_res = unique_res.size();

    // Partial result maps — one per worker, merged after join.
    std::vector<std::unordered_map<int, CAS_Result>> partial_maps(n_workers);
    for (size_t w = 0; w < n_workers; ++w) {
      for (int res_id : unique_res) {
        CAS_Result r;
        r.res_id = cas_data[res_id].res_id;
        r.res_name = cas_data[res_id].res_name;
        partial_maps[w][res_id] = r;
      }
    }

    // Shard unique_res across workers: worker w → [begin, end)
    auto get_shard = [&](size_t w) -> std::pair<size_t, size_t> {
      size_t begin = n_res * w / n_workers;
      size_t end = n_res * (w + 1) / n_workers;
      return {begin, end};
    };

    size_t active_lg_snap = (lig.m_mol.size() > 1) ? num_snap : 0;
    const MOL &active_lg_mol = lig.m_mol[active_lg_snap];
    size_t rec_size_ = pro_in.m_mol[num_snap].size();
    size_t active_lg_size_ = active_lg_mol.size();

    // Launch n_workers jthreads — each owns a deep copy of mol frame + params.
    std::vector<std::jthread> workers;
    workers.reserve(n_workers);

    for (size_t w = 0; w < n_workers; ++w) {
      auto [shard_begin, shard_end] = get_shard(w);
      if (shard_begin >= shard_end)
        continue; // empty shard

      // Deep copy — each worker mutates its own copy, never touches shared mol.
      MOL worker_mol_frame = mol.m_mol[num_snap];
      Multi_MOL_parameters worker_mol_params = mol.parameters;
      size_t m_size = worker_mol_frame.size();

      workers.emplace_back([&, w, shard_begin, shard_end,
                            worker_mol_frame = std::move(worker_mol_frame),
                            worker_mol_params = std::move(worker_mol_params),
                            m_size]() mutable {
        // Per-worker local radii buffers (no class member → no race)
        AlignedDoubleVector local_mut_free_radii(m_size, 0.0);
        AlignedDoubleVector local_mut_complex_radii(m_size, 0.0);

        auto &partial = partial_maps[w]; // only this worker writes here

        for (size_t ri = shard_begin; ri < shard_end; ++ri) {
          int res_id = unique_res[ri];

          // --- Find min distance (read-only dis_min) ---
          double min_d = -1.0;
          for (size_t t = 0; t < num_threads; ++t) {
            for (const auto &d_b : dis_min[t].base) {
              if (type == "receptor" && (int)d_b.resnr_rec == res_id) {
                if (min_d < 0 || d_b.d < min_d)
                  min_d = d_b.d;
              } else if (type == "ligand" && (int)d_b.resnr_ab == res_id) {
                if (min_d < 0 || d_b.d < min_d)
                  min_d = d_b.d;
              }
            }
          }
          if (min_d < 0)
            continue;

          // --- Mutate worker's private copy ---
          worker_mol_frame.mutate_to_alanine(res_id, mol.info,
                                             worker_mol_params);

          // --- Born radii on worker's private frame ---
          std::fill(local_mut_free_radii.begin(), local_mut_free_radii.end(),
                    0.0);
          std::fill(local_mut_complex_radii.begin(),
                    local_mut_complex_radii.end(), 0.0);

          cal_effective_born_radii_avx2(
              worker_mol_frame, worker_mol_params.radius, local_mut_free_radii,
              Arguments.igb, 0, m_size, nullptr, nullptr);

          if (type == "receptor" && active_lg_size_ > 0) {
            cal_effective_born_radii_avx2(
                worker_mol_frame, worker_mol_params.radius,
                local_mut_complex_radii, Arguments.igb, 0, m_size,
                &active_lg_mol, &lig.parameters.radius);
          } else if (type == "ligand" && rec_size_ > 0) {
            cal_effective_born_radii_avx2(
                worker_mol_frame, worker_mol_params.radius,
                local_mut_complex_radii, Arguments.igb, 0, m_size,
                &pro_in.m_mol[num_snap], &pro_in.parameters.radius);
          }

          // --- GBSA with 1 thread (worker IS the thread) ---
          GBSA_energy_snap mut_gbsa_total;
          GBSA_Aggregator local_agg(Arguments, info_output);
          info_ener_res local_infoener;

          if (type == "receptor") {
            cal_energy_gbsa_avx2(
                worker_mol_frame, active_lg_mol, mol.info, lig.info, CachedData,
                0, m_size, 0, active_lg_size_, mut_gbsa_total, num_snap,
                Arguments, worker_mol_params, lig.parameters, mol.lcpo_p1,
                mol.lcpo_p2, lig.lcpo_p1, lig.lcpo_p2, local_mut_free_radii,
                lg_free_radii, local_mut_complex_radii, lg_complex_radii,
                local_agg, &local_infoener);
          } else {
            cal_energy_gbsa_avx2(
                pro_in.m_mol[num_snap], worker_mol_frame, pro_in.info, mol.info,
                CachedData, 0, rec_size_, 0, m_size, mut_gbsa_total, num_snap,
                Arguments, pro_in.parameters, worker_mol_params, pro_in.lcpo_p1,
                pro_in.lcpo_p2, mol.lcpo_p1, mol.lcpo_p2, rec_free_radii,
                local_mut_free_radii, rec_complex_radii,
                local_mut_complex_radii, local_agg, &local_infoener);
          }

          // --- Accumulate into private partial map (no lock needed) ---
          auto &r = partial[res_id];
          r.ddG_vdw += (mut_gbsa_total.E_vdw - frame_gbsa_total.E_vdw);
          r.ddG_elec += (mut_gbsa_total.E_elec - frame_gbsa_total.E_elec);
          r.ddG_polar += (mut_gbsa_total.G_polar - frame_gbsa_total.G_polar);
          r.ddG_nonpolar +=
              (mut_gbsa_total.G_nonpolar - frame_gbsa_total.G_nonpolar);
          r.ddG_total +=
              (mut_gbsa_total.total_energy - frame_gbsa_total.total_energy);
          r.distance_sum += min_d;
          r.contact_frames++;

          // --- Restore worker's private copy ---
          worker_mol_frame.restore_mutation(res_id, mol.info,
                                            worker_mol_params);
        }
      }); // end jthread lambda
    } // end for w

    // jthread RAII: destructor calls request_stop() + join() automatically.
    // Explicitly join for clarity and to guarantee merge happens after.
    for (auto &jt : workers) {
      if (jt.joinable())
        jt.join();
    }

    // --- Merge partial maps into main cas_data (single-threaded) ---
    for (size_t w = 0; w < n_workers; ++w) {
      for (const auto &[res_id, partial_res] : partial_maps[w]) {
        auto &target = cas_data[res_id];
        target.ddG_vdw += partial_res.ddG_vdw;
        target.ddG_elec += partial_res.ddG_elec;
        target.ddG_polar += partial_res.ddG_polar;
        target.ddG_nonpolar += partial_res.ddG_nonpolar;
        target.ddG_total += partial_res.ddG_total;
        target.distance_sum += partial_res.distance_sum;
        target.contact_frames += partial_res.contact_frames;
      }
    }
  };

  perform_cas_scan(pro_in, cas_receptor_data, unique_res_rec, "receptor");
  if (Arguments.mol_2_type.find("protein") != std::string::npos) {
    perform_cas_scan(lig, cas_ligand_data, unique_res_lig, "ligand");
  }
}

void CASAnalyzer::export_results(const std::string &rec_file,
                                 const std::string &lig_file, double eu_par) {
  if (global_total_frames == 0)
    return;
  auto write_data = [&](const std::unordered_map<int, CAS_Result> &cas_data,
                        const std::string &file) {
    if (cas_data.empty())
      return;
    std::ofstream out(file);
    out << "Res_ID Res_Name Occupancy_Fraction Thermo_ddG Bound_ddG "
           "Avg_Distance\n";
    for (const auto &pair : cas_data) {
      const CAS_Result &res = pair.second;
      if (res.contact_frames == 0)
        continue;
      double occupancy = (double)res.contact_frames / global_total_frames;
      if (occupancy < 0.10)
        continue;
      double thermo_ddG = (res.ddG_total / global_total_frames) * eu_par;
      double bound_ddG = (res.ddG_total / res.contact_frames) * eu_par;
      double avg_distance = res.distance_sum / res.contact_frames;
      out << res.res_id << " " << res.res_name << " " << std::fixed
          << std::setprecision(4) << occupancy << " " << thermo_ddG << " "
          << bound_ddG << " " << avg_distance << "\n";
    }
  };
  write_data(cas_receptor_data, rec_file);
  write_data(cas_ligand_data, lig_file);
}

void CASAnalyzer::export_hotspot_pdb(const std::string &out_pdb_filename,
                                     double eu_par) {
  if (!has_reference || global_total_frames == 0)
    return;
  FILE *out_file = fopen(out_pdb_filename.c_str(), "w");
  if (!out_file)
    return;

  int atom_serial = 1;

  auto write_atoms = [&](const MOL &mol, const std::vector<mol_info> &info,
                         const std::unordered_map<int, CAS_Result> &cas_data,
                         size_t chain_offset) {
    for (size_t i = 0; i < mol.size(); ++i) {
      int res_id = info[i].get_resnr();
      double b_factor = 0.00;
      auto it = cas_data.find(res_id);
      if (it != cas_data.end()) {
        const CAS_Result &res = it->second;
        double occupancy = (double)res.contact_frames / global_total_frames;
        if (occupancy >= 0.10) {
          b_factor = (res.ddG_total / global_total_frames) * eu_par;
        }
      }

      size_t cid = info[i].get_chainnum() + chain_offset;
      char chain_char = (cid < 26) ? ('A' + cid) : ('a' + (cid - 26));

      fprintf(out_file,
              "ATOM  %5d %-4.4s %-4.4s%c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
              atom_serial++, info[i].get_atom().c_str(),
              info[i].get_res_name().c_str(), chain_char, info[i].get_resnr(),
              mol.pos.x[i], mol.pos.y[i], mol.pos.z[i], 1.00, b_factor);
    }
  };

  size_t max_pro_chainid = 0;
  for (const auto &mi : ref_pro_info) {
    if (mi.get_chainnum() > max_pro_chainid) {
      max_pro_chainid = mi.get_chainnum();
    }
  }

  write_atoms(ref_pro_frame, ref_pro_info, cas_receptor_data, 0);
  if (ref_lig_frame.size() > 0) {
    fprintf(out_file, "TER\n");
    write_atoms(ref_lig_frame, ref_lig_info, cas_ligand_data,
                max_pro_chainid + 1);
  }
  fprintf(out_file, "END\n");
  fclose(out_file);
}
