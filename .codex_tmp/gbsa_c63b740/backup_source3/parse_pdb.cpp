#include "Multi_MOL.h"
#include "functions.h"
#include <algorithm>
#include <charconv>
#include <fstream>
#include <iostream>
#include <string>
#include <string_view>
using std::string_view;
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

// Hàm trim khoảng trắng siêu tốc độ trả về string_view, KHÔNG cấp phát bộ nhớ
// động
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

// FIX: Đã bổ sung tham số 'force_f' bị thiếu để hàm emplace() có thể gọi được
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

  // Tối ưu C++20: Dùng starts_with thay vì substr() == "..."
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
    if (time_of_snap_temp > Arguments.end_time)
      return 1;

    if (!Arguments.flag_begin_set && I == 0) {
      output_info.set_current_begin_time(time_of_snap_temp);
    }

    if (check_in_time(Arguments, time_of_snap_temp)) {
      // FIX: Loại bỏ std::string allocation, dùng string_view kết hợp
      // trim_spaces siêu nhẹ
      std::string_view atom = trim_spaces(protein_1.substr(12, 4));
      std::string_view resname = trim_spaces(protein_1.substr(16, 5));

      try {
        // Tùy thuộc vào constructor của mol_info, bạn có thể truyền thẳng
        // string_view. Ở đây tôi ép kiểu tạm sang std::string để đảm bảo code
        // bạn biên dịch được.
        this->emplace(I, Box_temp, time_of_snap_temp, force_f,
                      mol_info(fast_parse_u64(protein_1.substr(22, 4)), "NULL",
                               std::string(atom),
                               std::string(protein_1.substr(21, 1)),
                               std::string(resname)),
                      [&]{ MOL m; m.pos.x.push_back(fast_parse_double(protein_1.substr(30, 8))); m.pos.y.push_back(fast_parse_double(protein_1.substr(38, 8))); m.pos.z.push_back(fast_parse_double(protein_1.substr(46, 8))); return m; }());
      } catch (const std::exception &e) {
        std::cerr << "Exception parsing ATOM line: " << e.what()
                  << "\nLine: " << protein_1 << "\n";
        return 1;
      }
    } else {
      flag_set_box = false;
    }
  }

  if (protein_1.starts_with("ENDMDL") || protein_1.starts_with("END")) {
    if (check_in_time(Arguments, time_of_snap_temp)) {
      if (this->m_mol.empty()) {
        std::cerr << "Error occurred in " << input_name
                  << " file!\nCannot parse atom coordinates of molecule\n";
        return 1;
      }

      check_and_pop_size_2last_snaps(I);

      if (!flag_set_box && I != 0 && this->m_mol.size() > I) {
        this->m_mol[I].Box = this->m_mol[I - 1].Box;
      } else {
        flag_set_box = false;
      }
      ++I;
    }
    flag_set_time = false;
  }
  return 0;
}

// -----------------------------------------------------------------------------
// 3. FILE PARSING LOOP
// -----------------------------------------------------------------------------

int Multi_MOL::parse_pdb(std::string_view input_name,
                         std::string_view input_topo, arguments &Arguments,
                         ::info &output_info) {
  // Không cần sync_with_stdio(false) vì đang dùng ifstream, không dùng std::cin
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

    while (end != std::string_view::npos) {
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
          // FIX: Xóa chuỗi rác "line)\n"; và truyền thêm force_f
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

  if (!remainder.empty()) {
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