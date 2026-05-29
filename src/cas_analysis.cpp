#include "cas_analysis.h"
#include "gbsa_intrin.h"
#include "titan_thread_compat.h"
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <thread>

CASAnalyzer::CASAnalyzer() : global_total_frames(0), has_reference(false) {}

void CASAnalyzer::init(Multi_MOL &pro_in, Multi_MOL &lig,
                       const std::string &mol_2_type) {
  auto is_retained_alanine_atom = [](const std::string &atom_name) {
    return atom_name == "N" || atom_name == "H" || atom_name == "CA" ||
           atom_name == "HA" || atom_name == "C" || atom_name == "O" ||
           atom_name == "CB" || atom_name == "HB1" || atom_name == "HB2" ||
           atom_name == "HB3";
  };

  auto init_cas_data = [&](Multi_MOL &mol,
                           std::unordered_map<int, CAS_Result> &cas_data,
                           std::unordered_map<int, std::vector<CASMutationAtom>>
                               &mutation_plans,
                           std::vector<int> &unique) {
    unique.clear();
    mutation_plans.clear();

    for (const auto &info : mol.info) {
      const int res_id = info.get_resnr();
      if (cas_data.find(res_id) == cas_data.end()) {
        unique.emplace_back(res_id);
        CAS_Result new_res;
        new_res.res_id = res_id;
        new_res.res_name = info.get_res_name();
        cas_data[res_id] = new_res;
      }
    }

    mutation_plans.reserve(unique.size());
    for (size_t atom_idx = 0; atom_idx < mol.info.size(); ++atom_idx) {
      const auto &atom_info = mol.info[atom_idx];
      if (is_retained_alanine_atom(atom_info.get_atom())) {
        continue;
      }

      auto &plan = mutation_plans[atom_info.get_resnr()];
      plan.push_back(CASMutationAtom{
          atom_idx,
          mol.parameters[atom_idx].get_charge(),
          mol.parameters[atom_idx].get_epsilon(),
          mol.parameters[atom_idx].get_sigma(),
          mol.parameters[atom_idx].get_radius()});
    }
  };

  if (cas_receptor_data.empty()) {
    init_cas_data(pro_in, cas_receptor_data, receptor_mutation_plans,
                  unique_res_rec);
  }
  if (cas_ligand_data.empty() &&
      mol_2_type.find("protein") != std::string::npos) {
    init_cas_data(lig, cas_ligand_data, ligand_mutation_plans, unique_res_lig);
  }
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
  (void)info_output;

  global_total_frames++;
  if (!has_reference) {
    ref_pro_frame = pro_in.m_mol[num_snap];
    if (!lig.m_mol.empty()) {
      ref_lig_frame = lig.m_mol[(lig.m_mol.size() > 1) ? num_snap : 0];
    }
    ref_pro_info = pro_in.info;
    ref_lig_info = lig.info;
    has_reference = true;
  }

  size_t total_contacts = 0;
  for (const auto &thread_contacts : dis_min) {
    total_contacts += thread_contacts.base.size();
  }

  std::unordered_map<int, double> rec_min_dist;
  std::unordered_map<int, double> lig_min_dist;
  rec_min_dist.reserve(total_contacts);
  lig_min_dist.reserve(total_contacts);

  auto update_min_distance = [](std::unordered_map<int, double> &min_map,
                                int res_id, double distance) {
    const auto [it, inserted] = min_map.emplace(res_id, distance);
    if (!inserted && distance < it->second) {
      it->second = distance;
    }
  };

  for (const auto &thread_contacts : dis_min) {
    for (const auto &contact : thread_contacts.base) {
      update_min_distance(rec_min_dist, static_cast<int>(contact.resnr_rec),
                          contact.d);
      update_min_distance(lig_min_dist, static_cast<int>(contact.resnr_ab),
                          contact.d);
    }
  }

  auto apply_mutation_plan = [](const std::vector<CASMutationAtom> &plan,
                                Multi_MOL_parameters &parameters) {
    for (const auto &atom : plan) {
      parameters[atom.atom_index].set_charge(0.0);
      parameters[atom.atom_index].set_epsilon_sigma_radius(0.0, 0.0, 0.001);
    }
  };

  auto restore_mutation_plan = [](const std::vector<CASMutationAtom> &plan,
                                  Multi_MOL_parameters &parameters) {
    for (const auto &atom : plan) {
      parameters[atom.atom_index].set_charge(atom.charge);
      parameters[atom.atom_index].set_epsilon_sigma_radius(
          atom.epsilon, atom.sigma, atom.radius);
    }
  };

  auto perform_cas_scan = [&](Multi_MOL &mol,
                              std::unordered_map<int, CAS_Result> &cas_data,
                              const std::unordered_map<int, double>
                                  &active_contact_distances,
                              const std::unordered_map<
                                  int, std::vector<CASMutationAtom>>
                                  &mutation_plans,
                              bool mutate_receptor) {
    if (active_contact_distances.empty()) {
      return;
    }

    std::vector<int> active_residues;
    std::vector<double> active_min_distances;
    active_residues.reserve(active_contact_distances.size());
    active_min_distances.reserve(active_contact_distances.size());
    for (const auto &[res_id, min_distance] : active_contact_distances) {
      active_residues.push_back(res_id);
      active_min_distances.push_back(min_distance);
    }

    const size_t n_res = active_residues.size();
    const size_t n_workers =
        (n_res == 0) ? size_t{0} : std::min(num_threads, n_res);

    auto get_shard = [&](size_t worker_idx) -> std::pair<size_t, size_t> {
      const size_t begin = n_res * worker_idx / n_workers;
      const size_t end = n_res * (worker_idx + 1) / n_workers;
      return {begin, end};
    };

    const size_t active_lg_snap = (lig.m_mol.size() > 1) ? num_snap : 0;
    const MOL &active_lg_mol = lig.m_mol[active_lg_snap];
    const size_t rec_size_ = pro_in.m_mol[num_snap].size();
    const size_t active_lg_size_ = active_lg_mol.size();
    const MOL &worker_frame = mol.m_mol[num_snap];

    std::vector<CAS_Result> frame_results(n_res);
    for (size_t idx = 0; idx < n_res; ++idx) {
      auto &result = frame_results[idx];
      result.res_id = active_residues[idx];
      result.res_name = cas_data[result.res_id].res_name;
    }

    std::vector<titan::jthread> workers;
    workers.reserve(n_workers);

    for (size_t worker_idx = 0; worker_idx < n_workers; ++worker_idx) {
      const auto shard = get_shard(worker_idx);
      const size_t shard_begin = shard.first;
      const size_t shard_end = shard.second;
      if (shard_begin >= shard_end) {
        continue;
      }

      Multi_MOL_parameters worker_mol_params = mol.parameters;
      const size_t m_size = worker_frame.size();

      workers.emplace_back([&, shard_begin, shard_end,
                            worker_mol_params = std::move(worker_mol_params),
                            m_size]() mutable {
        AlignedDoubleVector local_mut_free_radii(m_size, 0.0);
        AlignedDoubleVector local_mut_complex_radii(m_size, 0.0);

        for (size_t ri = shard_begin; ri < shard_end; ++ri) {
          const int res_id = active_residues[ri];
          const double min_d = active_min_distances[ri];
          auto &result = frame_results[ri];
          result.distance_sum += min_d;
          result.contact_frames++;

          const auto plan_it = mutation_plans.find(res_id);
          if (plan_it == mutation_plans.end() || plan_it->second.empty()) {
            continue;
          }

          apply_mutation_plan(plan_it->second, worker_mol_params);

          std::fill(local_mut_free_radii.begin(), local_mut_free_radii.end(),
                    0.0);
          std::fill(local_mut_complex_radii.begin(),
                    local_mut_complex_radii.end(), 0.0);

          cal_effective_born_radii_avx2(
              worker_frame, worker_mol_params.radius, local_mut_free_radii,
              Arguments.igb, 0, m_size, nullptr, nullptr);

          if (mutate_receptor && active_lg_size_ > 0) {
            cal_effective_born_radii_avx2(
                worker_frame, worker_mol_params.radius,
                local_mut_complex_radii, Arguments.igb, 0, m_size,
                &active_lg_mol, &lig.parameters.radius);
          } else if (!mutate_receptor && rec_size_ > 0) {
            cal_effective_born_radii_avx2(
                worker_frame, worker_mol_params.radius,
                local_mut_complex_radii, Arguments.igb, 0, m_size,
                &pro_in.m_mol[num_snap], &pro_in.parameters.radius);
          }

          GBSA_energy_snap mut_gbsa_total;
          if (mutate_receptor) {
            cal_energy_gbsa_avx2(
                worker_frame, active_lg_mol, mol.info, lig.info, CachedData, 0,
                m_size, 0, active_lg_size_, mut_gbsa_total, num_snap, Arguments,
                worker_mol_params, lig.parameters, mol.lcpo_p1, mol.lcpo_p2,
                lig.lcpo_p1, lig.lcpo_p2, local_mut_free_radii, lg_free_radii,
                local_mut_complex_radii, lg_complex_radii, nullptr, nullptr,
                nullptr, false);
          } else {
            cal_energy_gbsa_avx2(
                pro_in.m_mol[num_snap], worker_frame, pro_in.info, mol.info,
                CachedData, 0, rec_size_, 0, m_size, mut_gbsa_total, num_snap,
                Arguments, pro_in.parameters, worker_mol_params,
                pro_in.lcpo_p1, pro_in.lcpo_p2, mol.lcpo_p1, mol.lcpo_p2,
                rec_free_radii, local_mut_free_radii, rec_complex_radii,
                local_mut_complex_radii, nullptr, nullptr, nullptr, false);
          }

          restore_mutation_plan(plan_it->second, worker_mol_params);

          result.ddG_vdw += (mut_gbsa_total.E_vdw - frame_gbsa_total.E_vdw);
          result.ddG_elec += (mut_gbsa_total.E_elec - frame_gbsa_total.E_elec);
          result.ddG_polar +=
              (mut_gbsa_total.G_polar - frame_gbsa_total.G_polar);
          result.ddG_nonpolar +=
              (mut_gbsa_total.G_nonpolar - frame_gbsa_total.G_nonpolar);
          result.ddG_total +=
              (mut_gbsa_total.total_energy - frame_gbsa_total.total_energy);
        }
      });
    }

    for (auto &jt : workers) {
      if (jt.joinable()) {
        jt.join();
      }
    }

    for (const auto &result : frame_results) {
      auto &target = cas_data[result.res_id];
      target.ddG_vdw += result.ddG_vdw;
      target.ddG_elec += result.ddG_elec;
      target.ddG_polar += result.ddG_polar;
      target.ddG_nonpolar += result.ddG_nonpolar;
      target.ddG_total += result.ddG_total;
      target.distance_sum += result.distance_sum;
      target.contact_frames += result.contact_frames;
    }
  };

  perform_cas_scan(pro_in, cas_receptor_data, rec_min_dist,
                   receptor_mutation_plans, true);
  if (Arguments.mol_2_type.find("protein") != std::string::npos) {
    perform_cas_scan(lig, cas_ligand_data, lig_min_dist,
                     ligand_mutation_plans, false);
  }
}

void CASAnalyzer::export_results(const std::string &rec_file,
                                 const std::string &lig_file, double eu_par) {
  if (global_total_frames == 0) {
    return;
  }

  auto write_data = [&](const std::unordered_map<int, CAS_Result> &cas_data,
                        const std::string &file) {
    if (cas_data.empty()) {
      return;
    }

    std::ofstream out(file);
    out << "Res_ID Res_Name Occupancy_Fraction Thermo_ddG Bound_ddG "
           "Avg_Distance\n";
    for (const auto &pair : cas_data) {
      const CAS_Result &res = pair.second;
      if (res.contact_frames == 0) {
        continue;
      }

      const double occupancy =
          static_cast<double>(res.contact_frames) / global_total_frames;
      if (occupancy < 0.10) {
        continue;
      }

      const double thermo_ddG =
          (res.ddG_total / global_total_frames) * eu_par;
      const double bound_ddG = (res.ddG_total / res.contact_frames) * eu_par;
      const double avg_distance = res.distance_sum / res.contact_frames;
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
  if (!has_reference || global_total_frames == 0) {
    return;
  }

  FILE *out_file = fopen(out_pdb_filename.c_str(), "w");
  if (!out_file) {
    return;
  }

  int atom_serial = 1;

  auto write_atoms = [&](const MOL &mol, const std::vector<mol_info> &info,
                         const std::unordered_map<int, CAS_Result> &cas_data,
                         size_t chain_offset) {
    for (size_t i = 0; i < mol.size(); ++i) {
      const int res_id = info[i].get_resnr();
      double b_factor = 0.00;
      const auto it = cas_data.find(res_id);
      if (it != cas_data.end()) {
        const CAS_Result &res = it->second;
        const double occupancy =
            static_cast<double>(res.contact_frames) / global_total_frames;
        if (occupancy >= 0.10) {
          b_factor = (res.ddG_total / global_total_frames) * eu_par;
        }
      }

      const size_t cid = info[i].get_chainnum() + chain_offset;
      const char chain_char = (cid < 26) ? ('A' + cid) : ('a' + (cid - 26));

      fprintf(out_file,
              "ATOM  %5d %-4.4s %-4.4s%c%4zu    %8.3f%8.3f%8.3f%6.2f%6.2f\n",
              atom_serial++, info[i].get_atom().c_str(),
              info[i].get_res_name().c_str(), chain_char,
              static_cast<size_t>(info[i].get_resnr()), mol.pos.x[i],
              mol.pos.y[i], mol.pos.z[i], 1.00, b_factor);
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
