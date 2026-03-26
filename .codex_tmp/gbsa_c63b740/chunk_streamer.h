#ifndef CHUNK_STREAMER_H
#define CHUNK_STREAMER_H

#include "Multi_MOL.h"
#include "arguments.h"
#include "groups.h"
#include "info.h"
#include <cmath> // acos, sqrt for triclinic box conversion
#include <expected>
#include <fstream>
#include <string>
#include <string_view>
#include <vector>

// Forward-declare XDR types to avoid including xdrfile.h
// (which defines `typedef float matrix[3][3]` conflicting with Eigen)
struct XDRFILE;

class ChunkStreamer {
public:
  ChunkStreamer(arguments &args, info &output_info, groups_cal &groups,
                std::string_view force_f, size_t num_threads);
  ~ChunkStreamer();

  void set_frame0_positions(std::streampos pos_pro, std::streampos pos_lig) {
    frame0_end_pos_pro = pos_pro;
    frame0_end_pos_lig = pos_lig;
  }

  std::expected<void, std::string> open_streams(const std::string &name_pro,
                                                const std::string &name_lig,
                                                const Multi_MOL &init_pro,
                                                const Multi_MOL &init_lig);

  std::expected<size_t, std::string> read_next_chunk(size_t chunk_size,
                                                     Multi_MOL &pro_in,
                                                     Multi_MOL &lig,
                                                     Multi_MOL &water_traj);

private:
  arguments &args_ref;
  info &info_ref;
  groups_cal &groups_ref;
  std::string force_f;
  size_t nt;

  std::ifstream stream_pro;
  std::ifstream stream_lig;

  std::streampos frame0_end_pos_pro = -1;
  std::streampos frame0_end_pos_lig = -1;

  MOL last_pro_frame;
  MOL last_lig_frame;

  void parallel_split_coordinates(const Multi_MOL &source_chunk,
                                  const Multi_MOL &init_pro,
                                  const Multi_MOL &init_lig, Multi_MOL &pro_out,
                                  Multi_MOL &lig_out, Multi_MOL &wat_out,
                                  size_t num_frames);

  // ── XDR state ──
  // Using opaque pointer + float arrays to avoid xdrfile.h in header
  XDRFILE *xdr_handle_ = nullptr;
  int xdr_natoms_ = 0;
  float (*xdr_pos_buf_)[3] = nullptr; // rvec = float[3], allocated as array
  float xdr_box_[3][3]{};             // matrix = float[3][3]
  float xdr_prec_ = 1000.0f;
  bool xdr_eof_ = false;
  size_t xdr_file_idx_ = 0;

  const Multi_MOL *init_pro_ = nullptr;
  const Multi_MOL *init_lig_ = nullptr;

  // ── XDR time filter state (persists across chunks) ──
  double xdr_next_target_time_ = -1e18; // for -dt stride
};

#endif
