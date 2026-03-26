#include "Multi_MOL.h"
#include "functions.h"
#include "titan_ui.h"
#include <algorithm>
#include <charconv>
#include <fstream>
#include <iostream>
#include <string>
#include <string_view>
using std::string_view;
#include <expected>
#include <system_error>

// -----------------------------------------------------------------------------
// 1. HIGH-PERFORMANCE HELPER FUNCTIONS
// -----------------------------------------------------------------------------

double fast_parse_double(std::string_view data) {
  const char *begin = data.data();
  const char *end = begin + data.size();
  while (begin < end && *begin == ' ') {
    ++begin;
  }
  double value;
  auto [ptr, ec] = std::from_chars(begin, end, value);
  if (ec == std::errc())
    return value;
  return 0.0;
}

size_t fast_parse_u64(std::string_view data) {
  const char *begin = data.data();
  const char *end = begin + data.size();
  while (begin < end && *begin == ' ') {
    ++begin;
  }
  size_t value;
  auto [ptr, ec] = std::from_chars(begin, end, value);
  if (ec == std::errc())
    return value;
  return 0;
}

// HÃ m trim khoáº£ng tráº¯ng siÃªu tá»‘c Ä‘á»™ tráº£
// vá» string_view, KHÃ”NG cáº¥p phÃ¡t bá»™ nhá»›
// Ä‘á»™ng
constexpr std::string_view trim_spaces(std::string_view sv) {
  sv.remove_prefix(std::min(sv.find_first_not_of(" \t\r\n"), sv.size()));
  sv.remove_suffix(sv.empty() ? 0
                              : sv.size() - sv.find_last_not_of(" \t\r\n") - 1);
  return sv;
}

// -----------------------------------------------------------------------------
// 2. CORE PARSING LOGIC
// -----------------------------------------------------------------------------

int check_and_settime_new(double &time_of_snap, std::string_view protein_1,
                          std::string_view protein_2, arguments &Arguments) {
  if (protein_1.length() >= 5 &&
      protein_1.find("TITLE") != std::string_view::npos) {
    size_t found = protein_1.find('=');
    size_t found_2 = protein_1.find("step");
    size_t time_length = std::min(protein_1.length(), found_2);

    if (found != std::string_view::npos && protein_1.length() > found + 2) {
      time_of_snap = fast_parse_double(
          protein_1.substr(found + 1, time_length - found - 1));
      return 0;
    } else {
      found = protein_2.find("t=");
      found_2 = protein_2.find("step");
      time_length = std::min(protein_2.length(), found_2);
      if (found != std::string_view::npos) {
        time_of_snap = fast_parse_double(
            protein_2.substr(found + 2, time_length - found - 2));
      }
      return 0;
    }
  }
  return 1;
}

// FIX: ÄÃ£ bá»• sung tham sá»‘ 'force_f' bá»‹ thiáº¿u
// Ä‘á»ƒ hÃ m emplace() cÃ³ thá»ƒ gá»i Ä‘Æ°á»£c
int Multi_MOL::helper_parse_pdb(arguments &Arguments, box &Box_temp,
                                bool &flag_set_time, double &time_of_snap_temp,
                                std::string_view protein_1, size_t &I,
                                bool &flag_set_box, ::info &output_info,
                                std::string_view protein_2,
                                std::string_view input_name,
                                std::string_view force_f) {

  if (!flag_set_time) {
    time_of_snap_temp = static_cast<double>(I);
    if (check_and_settime_new(time_of_snap_temp, protein_1, protein_2,
                              Arguments) == 0) {
      flag_set_time = true;
    }
  }

  // Tá»‘i Æ°u C++20: DÃ¹ng starts_with thay vÃ¬ substr() == "..."
  if (!flag_set_box && protein_1.length() >= 54 &&
      protein_1.starts_with("CRYST1")) {
    Box_temp = box(fast_parse_double(protein_1.substr(6, 9)),
                   fast_parse_double(protein_1.substr(15, 9)),
                   fast_parse_double(protein_1.substr(24, 9)),
                   fast_parse_double(protein_1.substr(33, 7)),
                   fast_parse_double(protein_1.substr(40, 7)),
                   fast_parse_double(protein_1.substr(47, 7)));
    flag_set_box = true;
  }

  if (protein_1.length() >= 54 && protein_1.starts_with("ATOM")) {
    // FAST EXIT: If we passed the end time, signal EOF to stop reading the file
    // entirely
    if (Arguments.flag_end_set && time_of_snap_temp > Arguments.end_time)
      return 1;

    if (!Arguments.flag_begin_set && I == 0)
      output_info.set_current_begin_time(time_of_snap_temp);

    if (check_in_time(Arguments, time_of_snap_temp)) {

      // STRICT -dt LOGIC
      bool accept_dt = true;
      if (Arguments.dt > 0) {
        if (this->next_target_time == -1e9) {
          this->next_target_time = Arguments.flag_begin_set
                                       ? Arguments.begin_time
                                       : time_of_snap_temp;
        }
        if (time_of_snap_temp < this->next_target_time - 1e-4)
          accept_dt = false;
      }

      if (accept_dt) {
        std::string_view atom = trim_spaces(protein_1.substr(12, 4));
        std::string_view resname = trim_spaces(protein_1.substr(16, 5));
        try {
          this->emplace(I, Box_temp, time_of_snap_temp, force_f,
                        mol_info(fast_parse_u64(protein_1.substr(22, 4)),
                                 "NULL", std::string(atom),
                                 std::string(protein_1.substr(21, 1)),
                                 std::string(resname)),
                        [&] {
                          MOL m;
                          m.pos.x.emplace_back(
                              fast_parse_double(protein_1.substr(30, 8)));
                          m.pos.y.emplace_back(
                              fast_parse_double(protein_1.substr(38, 8)));
                          m.pos.z.emplace_back(
                              fast_parse_double(protein_1.substr(46, 8)));
                          return m;
                        }());
        } catch (...) {
          return 1;
        }
      } else {
        if (I < this->m_mol.size()) {
          this->m_mol[I].pos.x.clear();
          this->m_mol[I].pos.y.clear();
          this->m_mol[I].pos.z.clear();
        }
        flag_set_box = false;
      }
    } else {
      flag_set_box = false;
    }
  }

  if (protein_1.starts_with("ENDMDL") || protein_1.starts_with("END")) {
    // ONLY finalize frame if atoms were actually pushed
    if (I < this->m_mol.size() && !this->m_mol[I].pos.x.empty()) {
      // Missing atoms check
      size_t current_atoms = this->m_mol[I].pos.x.size();
      size_t expected_atoms = current_atoms;
      if (I > 0 && (I - 1) < this->m_mol.size()) {
        expected_atoms = this->m_mol[I - 1].pos.x.size();
      }

      if (current_atoms != expected_atoms) {
        TitanUI::print_warning("Incomplete frame detected at " +
                               std::to_string(time_of_snap_temp) +
                               " ps. Skipping...");
        this->m_mol[I].pos.x.clear();
        this->m_mol[I].pos.y.clear();
        this->m_mol[I].pos.z.clear();
      } else {
        check_and_pop_size_2last_snaps(I);
        if (!flag_set_box && I != 0 && this->m_mol.size() > I)
          this->m_mol[I].Box = this->m_mol[I - 1].Box;
        else
          flag_set_box = false;

        ++I; // Successfully committed frame

        // Advance Target Time for next accepted frame
        if (Arguments.dt > 0) {
          while (this->next_target_time <= time_of_snap_temp + 1e-4) {
            this->next_target_time += Arguments.dt;
          }
        }
      }
    }
    if (I < this->m_mol.size()) {
      this->m_mol[I].pos.x.clear();
      this->m_mol[I].pos.y.clear();
      this->m_mol[I].pos.z.clear();
    }
    flag_set_time = false;
  }
  return 0;
}

// -----------------------------------------------------------------------------
// 3. FILE PARSING LOOP
// -----------------------------------------------------------------------------

int Multi_MOL::parse_pdb(std::string_view input_name, std::string_view force_f,
                         arguments &Arguments, ::info &output_info,
                         size_t max_frames) {
  // KhÃ´ng cáº§n sync_with_stdio(false) vÃ¬ Ä‘ang dÃ¹ng ifstream, khÃ´ng dÃ¹ng
  // std::cin
  std::ifstream input_pdb(std::string(input_name), std::ios::binary);
  if (!input_pdb.is_open()) {
    std::cerr << "Can not open " << input_name << " file!\n";
    return 1;
  }

  double time_of_snap_temp{};
  box Box_temp;
  bool flag_set_time{false}, flag_set_box{false};
  size_t I{}, limit_snap{10};

  const size_t BUFFER_SIZE = 1024 * 1024; // 1 MB buffer
  std::vector<char> buffer(BUFFER_SIZE);
  std::string remainder;
  remainder.reserve(256);

  while (input_pdb) {
    input_pdb.read(buffer.data(), BUFFER_SIZE);
    std::streamsize bytes_read = input_pdb.gcount();
    if (bytes_read <= 0)
      break;

    std::string_view chunk(buffer.data(), static_cast<size_t>(bytes_read));
    size_t start = 0;
    size_t end = chunk.find('\n');

    if (max_frames > 0 && I >= max_frames) {
      break;
    }

    while (end != std::string_view::npos) {
      if (max_frames > 0 && I >= max_frames) {
        break;
      }

      std::string_view line_view = chunk.substr(start, end - start);

      if (!line_view.empty() && line_view.back() == '\r') {
        line_view.remove_suffix(1);
      }

      if (!remainder.empty()) {
        remainder.append(line_view);
        if (!remainder.empty()) {
          if (1 == helper_parse_pdb(Arguments, Box_temp, flag_set_time,
                                    time_of_snap_temp, remainder, I,
                                    flag_set_box, output_info, remainder,
                                    input_name, force_f)) {
            return 0;
          }
        }
        remainder.clear();
      } else {
        if (!line_view.empty()) {
          // FIX: XÃ³a chuá»—i rÃ¡c "line)\n"; vÃ  truyá»n thÃªm
          // force_f
          if (1 == helper_parse_pdb(Arguments, Box_temp, flag_set_time,
                                    time_of_snap_temp, line_view, I,
                                    flag_set_box, output_info, line_view,
                                    input_name, force_f)) {
            return 0;
          }
        }
      }
      start = end + 1;
      end = chunk.find('\n', start);
    }

    if (start < static_cast<size_t>(bytes_read)) {
      remainder.assign(chunk.substr(start));
    } else {
      remainder.clear();
    }
  }

  if (!remainder.empty() && !(max_frames > 0 && I >= max_frames)) {
    helper_parse_pdb(Arguments, Box_temp, flag_set_time, time_of_snap_temp,
                     remainder, I, flag_set_box, output_info, remainder,
                     input_name, force_f);
  }

  if (this->m_mol.size() >= 2) {
    output_info.set_current_dt(
        this->m_mol[this->m_mol.size() - 1].time_of_snap -
        this->m_mol[this->m_mol.size() - 2].time_of_snap);
  }

  if (!Arguments.flag_end_set && !this->m_mol.empty()) {
    output_info.set_current_end_time(
        this->m_mol[this->m_mol.size() - 1].time_of_snap);
  }

  if (Arguments.flag_begin_set && this->m_mol.empty() &&
      time_of_snap_temp < Arguments.begin_time) {
    std::cout << "WARNING!!! Your begin time is " << Arguments.begin_time
              << " ps while time of the last snapshot is smaller than that!\n";
  }
  return 0;
}

std::expected<size_t, std::string>
Multi_MOL::read_pdb_chunk(std::ifstream &input_pdb, size_t chunk_size,
                          arguments &Arguments, ::info &output_info,
                          std::string_view force_f) {
  if (!input_pdb.is_open())
    return std::unexpected("Input PDB stream is closed.");
  if (input_pdb.eof())
    return 0; // Normal EOF

  size_t frames_read = 0;
  double time_of_snap_temp = 0.0;
  box Box_temp;
  bool flag_set_time = false, flag_set_box = false;
  size_t I = 0;

  const size_t BUFFER_SIZE = 1024 * 1024;
  std::vector<char> buf(BUFFER_SIZE);
  std::string remainder;
  remainder.reserve(512);

  while (frames_read < chunk_size && input_pdb) {
    input_pdb.read(buf.data(), BUFFER_SIZE);
    std::streamsize bytes_read = input_pdb.gcount();
    if (bytes_read <= 0)
      break;

    std::string_view chunk_view(buf.data(), bytes_read);
    size_t start = 0;
    size_t end = chunk_view.find('\n');

    while (end != std::string_view::npos) {
      if (frames_read >= chunk_size) {
        break;
      }

      std::string_view line_view = chunk_view.substr(start, end - start);
      if (!line_view.empty() && line_view.back() == '\r') {
        line_view.remove_suffix(1);
      }

      size_t old_I = I;

      if (!remainder.empty()) {
        remainder.append(line_view);
        if (!remainder.empty()) {
          if (1 == helper_parse_pdb(Arguments, Box_temp, flag_set_time,
                                    time_of_snap_temp, remainder, I,
                                    flag_set_box, output_info, remainder,
                                    "stream", force_f)) {
            if (Arguments.flag_end_set &&
                time_of_snap_temp > Arguments.end_time) {
              this->m_mol.resize(I);
              input_pdb.seekg(0, std::ios_base::end);
              input_pdb.setstate(std::ios::eofbit);
              return frames_read;
            }
            this->m_mol.resize(I);
            return std::unexpected("Format error at line: " + remainder);
          }
        }
        remainder.clear();
      } else {
        if (!line_view.empty()) {
          if (1 == helper_parse_pdb(Arguments, Box_temp, flag_set_time,
                                    time_of_snap_temp, line_view, I,
                                    flag_set_box, output_info, line_view,
                                    "stream", force_f)) {
            if (Arguments.flag_end_set &&
                time_of_snap_temp > Arguments.end_time) {
              this->m_mol.resize(I);
              input_pdb.seekg(0, std::ios_base::end);
              input_pdb.setstate(std::ios::eofbit);
              return frames_read;
            }
            this->m_mol.resize(I);
            return std::unexpected("Format error at line: " +
                                   std::string(line_view));
          }
        }
      }

      if (I > old_I) {
        // ── A. END FILTER — hard stop (same as XDR path) ──
        if (Arguments.flag_end_set && time_of_snap_temp > Arguments.end_time) {
          // Exclude this over-limit frame
          --I;
          this->m_mol.resize(I);
          remainder.clear();
          // Signal EOF so streaming loop stops entirely
          input_pdb.seekg(0, std::ios_base::end);
          input_pdb.setstate(std::ios::eofbit);
          return frames_read;
        }

        // ── B. BEGIN FILTER — skip frames before -b ──
        if (Arguments.flag_begin_set &&
            time_of_snap_temp < Arguments.begin_time) {
          // Discard this frame: roll back I, don't increment frames_read
          --I;
          this->m_mol.resize(I);
          goto next_line; // skip to next line without counting this frame
        }

        // ── C. DT FILTER — stride (skip frames between dt intervals) ──
        if (Arguments.dt > 0.0) {
          static thread_local double pdb_next_target_time = -1e30;
          if (pdb_next_target_time < -1e17) {
            // Initialize on first accepted-by-begin frame
            pdb_next_target_time = Arguments.flag_begin_set
                                       ? Arguments.begin_time
                                       : time_of_snap_temp;
          }
          if (time_of_snap_temp < pdb_next_target_time - 1e-4) {
            // Frame doesn't hit target time — discard it
            --I;
            this->m_mol.resize(I);
            goto next_line;
          }
          // Advance target: keep stepping until strictly past time
          while (pdb_next_target_time <= time_of_snap_temp + 1e-4)
            pdb_next_target_time += Arguments.dt;
        }

        frames_read++;
      }

    next_line:
      start = end + 1;
      end = chunk_view.find('\n', start);
    }

    if (frames_read >= chunk_size) {
      std::streamoff rewind_bytes = static_cast<std::streamoff>(bytes_read) -
                                    static_cast<std::streamoff>(start);
      input_pdb.seekg(-rewind_bytes, std::ios_base::cur);
      remainder.clear();
      break;
    }

    if (start < static_cast<size_t>(bytes_read)) {
      remainder.assign(chunk_view.substr(start));
    } else {
      remainder.clear();
    }
  }

  if (!remainder.empty() && frames_read < chunk_size) {
    size_t old_I = I;
    if (1 == helper_parse_pdb(Arguments, Box_temp, flag_set_time,
                              time_of_snap_temp, remainder, I, flag_set_box,
                              output_info, remainder, "stream", force_f)) {
    }
    if (I > old_I) {
      // Apply same time filters as main loop
      bool accept = true;
      if (Arguments.flag_end_set && time_of_snap_temp > Arguments.end_time) {
        --I;
        accept = false;
      }
      if (accept && Arguments.flag_begin_set &&
          time_of_snap_temp < Arguments.begin_time) {
        --I;
        accept = false;
      }
      if (accept)
        frames_read++;
    }
  }

  // Erase any allocated ghost frames that were skipped by time filters
  this->m_mol.resize(I); // Always resize to exactly I, even if I == 0

  if (this->parameters.charge.empty() && !this->info.empty()) {
    this->gen_vector_parameters();
  }

  return frames_read;
}
