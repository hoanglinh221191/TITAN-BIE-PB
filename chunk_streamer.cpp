#include "chunk_streamer.h"
#include "titan_ui.h"
#include <algorithm>
#include <filesystem>
#include <immintrin.h> // AVX2 intrinsics
#include <iostream>
#include <thread>

// XDR headers included HERE (not in .h) to avoid `matrix` typedef clash with
// Eigen
extern "C" {
#include "xdrfile.h"
#include "xdrfile_trr.h"
#include "xdrfile_xtc.h"
}

namespace {

void initialize_stream_frame(MOL &frame, const Multi_MOL &source_topology,
                             const box &frame_box, double time_ps,
                             size_t atom_count) {
  frame.force_f = source_topology.force_f;
  frame.mol_type = std::string(source_topology.get_mol_type());
  frame.Box = frame_box;
  frame.time_of_snap = time_ps;
  frame.flag_set_time = true;
  frame.compute_forces = false;
  frame.pos.resize(atom_count);
}

void initialize_stream_water_frame(MOL &frame, std::string_view force_field,
                                   const box &frame_box, double time_ps,
                                   size_t atom_count) {
  frame.force_f = std::string(force_field);
  frame.mol_type = "water";
  frame.Box = frame_box;
  frame.time_of_snap = time_ps;
  frame.flag_set_time = true;
  frame.compute_forces = false;
  frame.pos.resize(atom_count);
}

} // namespace

ChunkStreamer::ChunkStreamer(arguments &args, info &output_info,
                             groups_cal &groups, std::string_view force_field,
                             size_t num_threads)
    : args_ref(args), info_ref(output_info), groups_ref(groups),
      force_f(force_field), nt(num_threads) {}

ChunkStreamer::~ChunkStreamer() {
  // XDR cleanup
  if (xdr_handle_) {
    xdrfile_close(xdr_handle_);
    xdr_handle_ = nullptr;
  }
  if (xdr_pos_buf_) {
    delete[] xdr_pos_buf_;
    xdr_pos_buf_ = nullptr;
  }
  // PDB cleanup
  if (stream_pro.is_open())
    stream_pro.close();
  if (stream_lig.is_open())
    stream_lig.close();
}

std::expected<void, std::string> ChunkStreamer::open_streams(
    const std::string &name_pro, const std::string &name_lig,
    const Multi_MOL &init_pro, const Multi_MOL &init_lig) {

  // Store topology references for later use in parallel_split_coordinates
  init_pro_ = &init_pro;
  init_lig_ = &init_lig;

  // ── XDR PATH (XTC/TRR) ──────────────────────────────────────────────
  if (args_ref.flag_xdr) {
    namespace fs = std::filesystem;
    const fs::path traj_path = fs::absolute(fs::path(args_ref.traj_file));
    const std::string tf = traj_path.string();
    std::error_code ec;
    const bool use_lig_group = !args_ref.flag_nocal;

    if (!fs::exists(traj_path, ec) || ec) {
      return std::unexpected("XDR: trajectory file not found: " + tf);
    }

    const auto traj_size = fs::file_size(traj_path, ec);
    if (ec) {
      return std::unexpected("XDR: cannot stat trajectory file: " + tf);
    }
    if (traj_size == 0) {
      return std::unexpected("XDR: trajectory file is empty: " + tf);
    }

    // 1. Get natoms from XDR file header
    int natoms_check = 0;
    int rc;
    if (args_ref.xdr_type == "xtc")
      rc = read_xtc_natoms(const_cast<char *>(tf.c_str()), &natoms_check);
    else
      rc = read_trr_natoms(const_cast<char *>(tf.c_str()), &natoms_check);

    if (rc != exdrOK)
      return std::unexpected("XDR: Cannot read natoms from " + tf +
                             " (size=" + std::to_string(traj_size) +
                             " bytes, rc=" + std::to_string(rc) + ")");

    // 2. Validate natoms: XDR must have >= atoms referenced by index groups
    size_t topology_natoms =
        groups_ref.groups_atom_number(groups_ref.group_1_index());
    if (use_lig_group)
      topology_natoms += groups_ref.groups_atom_number(groups_ref.group_2_index());
    if (args_ref.n_wat >= 0)
      topology_natoms += groups_ref.groups_atom_number(args_ref.n_wat);

    if (natoms_check < static_cast<int>(topology_natoms))
      return std::unexpected("XDR: natoms too small — XDR has " +
                             std::to_string(natoms_check) +
                             " but index groups reference " +
                             std::to_string(topology_natoms) + " atoms");
    if (natoms_check != static_cast<int>(topology_natoms)) {
      TitanUI::print_info("XDR will stream " + std::to_string(natoms_check) +
                          " atoms/frame and subset " +
                          std::to_string(topology_natoms) +
                          " selected atoms on-the-fly.");
    }
    xdr_natoms_ = natoms_check;

    auto cache_group_indices = [&](size_t group_index,
                                   std::vector<size_t> &out_indices,
                                   size_t &out_size)
        -> std::expected<void, std::string> {
      out_size = groups_ref.groups_atom_number(group_index);
      out_indices.clear();
      out_indices.reserve(out_size);
      for (size_t i = 0; i < out_size; ++i) {
        const size_t atom_index_1based =
            groups_ref.atom_index(group_index, i);
        if (atom_index_1based == 0) {
          return std::unexpected("XDR: invalid atom index in selected group.");
        }
        out_indices.push_back(atom_index_1based - 1);
      }
      return {};
    };

    auto pro_cache =
        cache_group_indices(groups_ref.group_1_index(), pro_atom_indices0_,
                            pro_size_);
    if (!pro_cache)
      return std::unexpected(pro_cache.error());

    if (use_lig_group) {
      auto lig_cache =
          cache_group_indices(groups_ref.group_2_index(), lig_atom_indices0_,
                              lig_size_);
      if (!lig_cache)
        return std::unexpected(lig_cache.error());
    } else {
      lig_atom_indices0_.clear();
      lig_size_ = 0;
    }

    if (args_ref.n_wat >= 0) {
      auto wat_cache =
          cache_group_indices(args_ref.n_wat, wat_atom_indices0_, wat_size_);
      if (!wat_cache)
        return std::unexpected(wat_cache.error());
    } else {
      wat_atom_indices0_.clear();
      wat_size_ = 0;
    }

    size_t max_selected_atom_index = 0;
    auto update_max_index = [&](const std::vector<size_t> &indices) {
      if (!indices.empty()) {
        max_selected_atom_index =
            std::max(max_selected_atom_index,
                     *std::max_element(indices.begin(), indices.end()));
      }
    };
    update_max_index(pro_atom_indices0_);
    if (use_lig_group)
      update_max_index(lig_atom_indices0_);
    update_max_index(wat_atom_indices0_);
    if (natoms_check <= static_cast<int>(max_selected_atom_index)) {
      return std::unexpected(
          "XDR: natoms too small for selected atom indices - XDR has " +
          std::to_string(natoms_check) + " atoms but needs at least " +
          std::to_string(max_selected_atom_index + 1));
    }

    // 3. Open XDR stream
    xdr_handle_ = xdrfile_open(tf.c_str(), "r");
    if (!xdr_handle_)
      return std::unexpected("XDR: Cannot open " + tf);

    // 4. Pre-allocate position buffer ONCE
    xdr_pos_buf_ = new rvec[xdr_natoms_];
    xdr_eof_ = false;
    xdr_next_target_time_ = -1e18;

    return {}; // success — don't fall through to PDB
  }

  // ── PDB PATH (unchanged) ────────────────────────────────────────────
  if (args_ref.flag_1_pdb) {
    stream_pro.open(args_ref.pdb_name, std::ios::in);
    if (!stream_pro.is_open())
      return std::unexpected("Failed to open stream: " + args_ref.pdb_name);

    if (frame0_end_pos_pro >= 0) {
      stream_pro.seekg(frame0_end_pos_pro);
    } else {
      // Fallback: parse-skip if position was not saved
      Multi_MOL dummy_mol;
      auto sk =
          dummy_mol.read_pdb_chunk(stream_pro, 1, args_ref, info_ref, force_f);
      if (!sk)
        return std::unexpected("Failed to skip Frame 0 in stream");
    }
  } else {
    stream_pro.open(name_pro, std::ios::in);
    stream_lig.open(name_lig, std::ios::in);

    if (!stream_pro.is_open())
      return std::unexpected("Failed to open stream: " + name_pro);
    if (!stream_lig.is_open())
      return std::unexpected("Failed to open stream: " + name_lig);

    // Skip frame 0
    if (frame0_end_pos_pro >= 0) {
      stream_pro.seekg(frame0_end_pos_pro);
    } else {
      Multi_MOL dummy_mol;
      dummy_mol.read_pdb_chunk(stream_pro, 1, args_ref, info_ref, force_f);
    }

    if (frame0_end_pos_lig >= 0) {
      stream_lig.seekg(frame0_end_pos_lig);
    } else {
      Multi_MOL dummy_mol;
      dummy_mol.read_pdb_chunk(stream_lig, 1, args_ref, info_ref, force_f);
    }
  }
  return {};
}

// ── FIX 1: parallel_split_coordinates with topology copy ──────────────
void ChunkStreamer::parallel_split_coordinates(
    const Multi_MOL &source_chunk, const Multi_MOL &init_pro,
    const Multi_MOL &init_lig, Multi_MOL &pro_out, Multi_MOL &lig_out,
    Multi_MOL &wat_out, size_t num_frames) {
  // Resize outputs
  pro_out.m_mol.resize(num_frames);
  lig_out.m_mol.resize(num_frames);
  if (args_ref.n_wat >= 0)
    wat_out.m_mol.resize(num_frames);

  // Pre-calculate target sizes to avoid redundant lookups
  size_t pro_size = groups_ref.groups_atom_number(groups_ref.group_1_index());
  size_t lig_size =
      args_ref.flag_nocal ? 0
                          : groups_ref.groups_atom_number(groups_ref.group_2_index());
  size_t wat_size =
      (args_ref.n_wat >= 0) ? groups_ref.groups_atom_number(args_ref.n_wat) : 0;

  // Lambda to process a range of frames
  auto process_frames = [&](size_t start_f, size_t end_f) {
    for (size_t f = start_f; f < end_f; ++f) {
      // ── FIX: Full structural copy from Frame 0 (topology: charge,
      //    type_of_atom, res_name, atom_name, etc.) then overwrite coords.
      pro_out.m_mol[f] = init_pro.m_mol[0];
      lig_out.m_mol[f] = init_lig.m_mol[0];

      // Overwrite box + time from this frame's data
      pro_out.m_mol[f].Box = source_chunk.m_mol[f].Box;
      pro_out.m_mol[f].time_of_snap = source_chunk.m_mol[f].time_of_snap;
      lig_out.m_mol[f].Box = source_chunk.m_mol[f].Box;
      lig_out.m_mol[f].time_of_snap = source_chunk.m_mol[f].time_of_snap;

      // Overwrite ONLY coordinates from source_chunk
      pro_out.m_mol[f].pos.x.resize(pro_size);
      pro_out.m_mol[f].pos.y.resize(pro_size);
      pro_out.m_mol[f].pos.z.resize(pro_size);

      for (size_t j = 0; j < pro_size; ++j) {
        size_t atom_idx =
            groups_ref.atom_index(groups_ref.group_1_index(), j) - 1;
        if (atom_idx < source_chunk.m_mol[f].pos.x.size()) {
          pro_out.m_mol[f].pos.x[j] = source_chunk.m_mol[f].pos.x[atom_idx];
          pro_out.m_mol[f].pos.y[j] = source_chunk.m_mol[f].pos.y[atom_idx];
          pro_out.m_mol[f].pos.z[j] = source_chunk.m_mol[f].pos.z[atom_idx];
        }
      }

      // 2. Ligand: overwrite coords only
      lig_out.m_mol[f].pos.x.resize(lig_size);
      lig_out.m_mol[f].pos.y.resize(lig_size);
      lig_out.m_mol[f].pos.z.resize(lig_size);

      for (size_t j = 0; j < lig_size; ++j) {
        size_t atom_idx =
            groups_ref.atom_index(groups_ref.group_2_index(), j) - 1;
        if (atom_idx < source_chunk.m_mol[f].pos.x.size()) {
          lig_out.m_mol[f].pos.x[j] = source_chunk.m_mol[f].pos.x[atom_idx];
          lig_out.m_mol[f].pos.y[j] = source_chunk.m_mol[f].pos.y[atom_idx];
          lig_out.m_mol[f].pos.z[j] = source_chunk.m_mol[f].pos.z[atom_idx];
        }
      }

      // 3. Setup Water (no topology copy needed for water)
      if (wat_size > 0) {
        wat_out.m_mol[f].Box = source_chunk.m_mol[f].Box;
        wat_out.m_mol[f].time_of_snap = source_chunk.m_mol[f].time_of_snap;
        wat_out.m_mol[f].pos.x.resize(wat_size);
        wat_out.m_mol[f].pos.y.resize(wat_size);
        wat_out.m_mol[f].pos.z.resize(wat_size);

        for (size_t j = 0; j < wat_size; ++j) {
          size_t atom_idx = groups_ref.atom_index(args_ref.n_wat, j) - 1;
          if (atom_idx < source_chunk.m_mol[f].pos.x.size()) {
            wat_out.m_mol[f].pos.x[j] = source_chunk.m_mol[f].pos.x[atom_idx];
            wat_out.m_mol[f].pos.y[j] = source_chunk.m_mol[f].pos.y[atom_idx];
            wat_out.m_mol[f].pos.z[j] = source_chunk.m_mol[f].pos.z[atom_idx];
          }
        }
      }
    }
  };

  // Spin up threads
  std::vector<std::jthread> workers;
  size_t chunk_per_thread = (num_frames + nt - 1) / nt;
  for (size_t t = 0; t < nt; ++t) {
    size_t start = t * chunk_per_thread;
    size_t end = std::min(start + chunk_per_thread, num_frames);
    if (start < end)
      workers.emplace_back(process_frames, start, end);
  }
}

std::expected<size_t, std::string>
ChunkStreamer::read_next_chunk(size_t chunk_size, Multi_MOL &pro_in,
                               Multi_MOL &lig, Multi_MOL &water_traj) {
  pro_in.m_mol.clear();
  lig.m_mol.clear();
  if (args_ref.n_wat >= 0)
    water_traj.m_mol.clear();

  // ── XDR PATH (XTC/TRR) ──────────────────────────────────────────────
  if (args_ref.flag_xdr) {
    if (xdr_eof_ || !xdr_handle_)
      return 0;

    constexpr double NM_TO_ANG = 10.0;

    const size_t pro_size = pro_size_;
    const size_t lig_size = lig_size_;
    const size_t wat_size = wat_size_;

    // Reserve frame slots only. Per-frame topology is shared from frame 0.
    pro_in.m_mol.reserve(chunk_size);
    lig.m_mol.reserve(chunk_size);
    if (wat_size > 0)
      water_traj.m_mol.reserve(chunk_size);

    size_t frames_read = 0;
    size_t skip_count = 0;

    // Read frames until chunk is filled or EOF/end-time
    while (frames_read < chunk_size) {
      int step;
      float time_ps, lambda = 0.0f;
      int rc;

      if (args_ref.xdr_type == "xtc") {
        rc = read_xtc(xdr_handle_, xdr_natoms_, &step, &time_ps, xdr_box_,
                      xdr_pos_buf_, &xdr_prec_);
      } else {
        int has_prop = 0;
        rc = read_trr(xdr_handle_, xdr_natoms_, &step, &time_ps, &lambda,
                      xdr_box_, xdr_pos_buf_, nullptr, nullptr, &has_prop);
      }

      if (rc != exdrOK) {
        xdr_eof_ = true;
        break;
      }

      double time_val = static_cast<double>(time_ps);

      // ── A. END FILTER — hard stop ──
      if (args_ref.flag_end_set && time_val > args_ref.end_time) {
        xdr_eof_ = true;
        break;
      }

      // ── B. BEGIN FILTER — skip frames before -b ──
      if (args_ref.flag_begin_set && time_val < args_ref.begin_time) {
        ++skip_count;
        if (skip_count % 50 == 0) {
          TitanUI::show_continuous_progress(
              "Skipping", skip_count,
              std::to_string(static_cast<int>(time_val)) + " ps (before -b=" +
                  std::to_string(static_cast<int>(args_ref.begin_time)) + ")");
        }
        continue;
      }
      if (skip_count > 0) {
        std::cout << "\n";
        TitanUI::print_info("[TITAN] Skipped " + std::to_string(skip_count) +
                            " frames. Starting analysis at t=" +
                            std::to_string(static_cast<int>(time_val)) + " ps");
        skip_count = 0;
      }

      // ── C. DT FILTER — stride ──
      if (args_ref.dt > 0.0) {
        if (xdr_next_target_time_ < -1e17) {
          // Initialize on first accepted-by-begin frame
          xdr_next_target_time_ =
              args_ref.flag_begin_set ? args_ref.begin_time : time_val;
        }
        if (time_val < xdr_next_target_time_ - 1e-4) {
          continue; // skip this frame (doesn't hit target time)
        }
        // Advance target: keep stepping until strictly past time_val
        while (xdr_next_target_time_ <= time_val + 1e-4)
          xdr_next_target_time_ += args_ref.dt;
      }

      // ── Frame passed all filters — store at frames_read index ──

      // ── Compute box from GROMACS matrix ──
      double v1x = static_cast<double>(xdr_box_[0][0]) * NM_TO_ANG;
      double v1y = static_cast<double>(xdr_box_[0][1]) * NM_TO_ANG;
      double v1z = static_cast<double>(xdr_box_[0][2]) * NM_TO_ANG;
      double v2x = static_cast<double>(xdr_box_[1][0]) * NM_TO_ANG;
      double v2y = static_cast<double>(xdr_box_[1][1]) * NM_TO_ANG;
      double v2z = static_cast<double>(xdr_box_[1][2]) * NM_TO_ANG;
      double v3x = static_cast<double>(xdr_box_[2][0]) * NM_TO_ANG;
      double v3y = static_cast<double>(xdr_box_[2][1]) * NM_TO_ANG;
      double v3z = static_cast<double>(xdr_box_[2][2]) * NM_TO_ANG;

      double la = std::sqrt(v1x * v1x + v1y * v1y + v1z * v1z);
      double lb = std::sqrt(v2x * v2x + v2y * v2y + v2z * v2z);
      double lc = std::sqrt(v3x * v3x + v3y * v3y + v3z * v3z);

      constexpr double RAD_TO_DEG = 180.0 / 3.14159265358979;
      double cos_alpha = (lb > 0 && lc > 0)
                             ? (v2x * v3x + v2y * v3y + v2z * v3z) / (lb * lc)
                             : 0.0;
      double cos_beta = (la > 0 && lc > 0)
                            ? (v1x * v3x + v1y * v3y + v1z * v3z) / (la * lc)
                            : 0.0;
      double cos_gamma = (la > 0 && lb > 0)
                             ? (v1x * v2x + v1y * v2y + v1z * v2z) / (la * lb)
                             : 0.0;
      cos_alpha = std::max(-1.0, std::min(1.0, cos_alpha));
      cos_beta = std::max(-1.0, std::min(1.0, cos_beta));
      cos_gamma = std::max(-1.0, std::min(1.0, cos_gamma));

      box frame_box;
      frame_box.set_a(la);
      frame_box.set_b(lb);
      frame_box.set_c(lc);
      frame_box.set_alpha(std::acos(cos_alpha) * RAD_TO_DEG);
      frame_box.set_beta(std::acos(cos_beta) * RAD_TO_DEG);
      frame_box.set_gamma(std::acos(cos_gamma) * RAD_TO_DEG);
      frame_box.setv();

      // ── Protein: topology copy + extract coords from rvec ──
      pro_in.m_mol.emplace_back();
      MOL &pro_frame = pro_in.m_mol.back();
      initialize_stream_frame(pro_frame, *init_pro_, frame_box, time_val,
                              pro_size);
      for (size_t j = 0; j < pro_size; ++j) {
        const size_t ai = pro_atom_indices0_[j];
        pro_frame.pos.x[j] = static_cast<double>(xdr_pos_buf_[ai][0]) * NM_TO_ANG;
        pro_frame.pos.y[j] = static_cast<double>(xdr_pos_buf_[ai][1]) * NM_TO_ANG;
        pro_frame.pos.z[j] = static_cast<double>(xdr_pos_buf_[ai][2]) * NM_TO_ANG;
      }

      // ── Ligand: topology copy + extract coords from rvec ──
      lig.m_mol.emplace_back();
      MOL &lig_frame = lig.m_mol.back();
      initialize_stream_frame(lig_frame, *init_lig_, frame_box, time_val,
                              lig_size);
      for (size_t j = 0; j < lig_size; ++j) {
        const size_t ai = lig_atom_indices0_[j];
        lig_frame.pos.x[j] = static_cast<double>(xdr_pos_buf_[ai][0]) * NM_TO_ANG;
        lig_frame.pos.y[j] = static_cast<double>(xdr_pos_buf_[ai][1]) * NM_TO_ANG;
        lig_frame.pos.z[j] = static_cast<double>(xdr_pos_buf_[ai][2]) * NM_TO_ANG;
      }

      // ── Water: extract coords only (no topology copy needed) ──
      if (wat_size > 0) {
        water_traj.m_mol.emplace_back();
        MOL &water_frame = water_traj.m_mol.back();
        initialize_stream_water_frame(water_frame, force_f, frame_box, time_val,
                                      wat_size);
        for (size_t j = 0; j < wat_size; ++j) {
          const size_t ai = wat_atom_indices0_[j];
          water_frame.pos.x[j] =
              static_cast<double>(xdr_pos_buf_[ai][0]) * NM_TO_ANG;
          water_frame.pos.y[j] =
              static_cast<double>(xdr_pos_buf_[ai][1]) * NM_TO_ANG;
          water_frame.pos.z[j] =
              static_cast<double>(xdr_pos_buf_[ai][2]) * NM_TO_ANG;
        }
      }

      ++frames_read;
    }

    if (frames_read == 0)
      return 0;

    return frames_read;
  }

  // ── PDB PATH (unchanged) ────────────────────────────────────────────
  if (args_ref.flag_1_pdb == false) {
    size_t pro_frames = 0;
    size_t lig_frames = 0;

    // REQUIREMENT 2: Read only if stream is still alive and hasn't failed
    if (!stream_pro.eof() && !stream_pro.fail()) {
      auto expected_pro = pro_in.read_pdb_chunk(stream_pro, chunk_size,
                                                args_ref, info_ref, force_f);
      if (!expected_pro)
        return std::unexpected("Receptor read failed: " + expected_pro.error());
      pro_frames = expected_pro.value();
      if (pro_frames > 0)
        last_pro_frame = pro_in.m_mol.back(); // Track last known frame
    }

    if (!stream_lig.eof() && !stream_lig.fail()) {
      auto expected_lig = lig.read_pdb_chunk(stream_lig, chunk_size, args_ref,
                                             info_ref, force_f);
      if (!expected_lig)
        return std::unexpected("Ligand read failed: " + expected_lig.error());
      lig_frames = expected_lig.value();
      if (lig_frames > 0)
        last_lig_frame = lig.m_mol.back(); // Track last known frame
    }

    size_t actual_frames = std::max(pro_frames, lig_frames);

    // Stop if BOTH trajectories are entirely exhausted
    if (actual_frames == 0)
      return 0;

    // Sync shorter trajectory by holding its last known frame
    if (pro_frames < actual_frames) {
      if (pro_frames == 0)
        TitanUI::print_warning(
            "Receptor trajectory ended early. Freezing at last frame.");
      pro_in.m_mol.resize(actual_frames, last_pro_frame);
    }
    if (lig_frames < actual_frames) {
      if (lig_frames == 0)
        TitanUI::print_warning(
            "Ligand trajectory ended early. Freezing at last frame.");
      lig.m_mol.resize(actual_frames, last_lig_frame);
    }

    return actual_frames;
  } else {
    Multi_MOL chunk_pdb;
    auto expected_pro = chunk_pdb.read_pdb_chunk(stream_pro, chunk_size,
                                                 args_ref, info_ref, force_f);

    if (!expected_pro)
      return std::unexpected("Chunk read failed: " + expected_pro.error());
    if (expected_pro.value() == 0)
      return 0; // EOF

    size_t frames = expected_pro.value();
    // PDB chunks already contain topology data, pass stored init references
    parallel_split_coordinates(chunk_pdb, *init_pro_, *init_lig_, pro_in, lig,
                               water_traj, frames);

    return frames;
  }
}
