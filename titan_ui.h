#ifndef TITAN_UI_H
#define TITAN_UI_H

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

namespace TitanUI {

inline void show_progress(const std::string &module_name, size_t current,
                          size_t total, const std::string &extra_info = "") {
  if (total == 0)
    return;

  // current ranges from 0 to total - 1. Add 1 to make it reach 100% at the last
  // step.
  float progress = static_cast<float>(current + 1) / total;
  int barWidth = 40;
  std::cout << "\r\033[1;36m[TITAN]\033[0m \033[1;33m" << std::setw(15)
            << std::left << module_name << "\033[0m [";

  int pos = barWidth * progress;
  for (int i = 0; i < barWidth; ++i) {
    if (i < pos)
      std::cout << "\033[1;32m=\033[0m";
    else if (i == pos)
      std::cout << "\033[1;32m>\033[0m";
    else
      std::cout << " ";
  }

  std::cout << "] " << std::setw(3) << std::right << int(progress * 100.0)
            << " % " << extra_info << "   " << std::flush;

  if (current + 1 >= total) {
    std::cout << std::endl;
  }
}

inline void show_continuous_progress(const std::string &module_name,
                                     size_t current_frame,
                                     const std::string &extra_info = "") {
  std::cout << "\r\033[1;36m[TITAN]\033[0m \033[1;33m" << std::setw(15)
            << std::left << module_name << "\033[0m Frame: \033[1;32m"
            << std::setw(6) << std::right << current_frame << "\033[0m "
            << extra_info << "   " << std::flush;
}

inline void print_module_header(const std::string &header_name) {
  std::cout << "\n\033[1;35m>>> " << header_name << " <<<\033[0m\n";
}

inline void print_logo() {
  std::cout
      << "\n\033[1;36m"
      << "  ████████╗██╗████████╗ █████╗ ███╗   ██╗\n"
      << "  ╚══██╔══╝██║╚══██╔══╝██╔══██╗████╗  ██║\n"
      << "     ██║   ██║   ██║   ███████║██╔██╗ ██║\n"
      << "     ██║   ██║   ██║   ██╔══██║██║╚██╗██║\n"
      << "     ██║   ██║   ██║   ██║  ██║██║ ╚████║\n"
      << "     ╚═╝   ╚═╝   ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═══╝\n"
      << "\033[0m"
      << "  \033[1;33mThermodynamic & Interaction Trajectory ANalyzer\033[0m\n";
  std::cout
      << "  "
         "\033[1;32m===============================================\033[0m\n\n";
}

inline void log_section(std::fstream &log_file, const std::string &title) {
  if (!log_file.is_open())
    return;
  log_file << "\n=============================================================="
              "==================\n";
  log_file << "[ " << title << " ]\n";
  log_file << "================================================================"
              "================\n";
}

inline void log_kv(std::fstream &log_file, const std::string &key,
                   const std::string &value) {
  if (!log_file.is_open())
    return;
  log_file << std::left << std::setw(30) << key << " : " << value << "\n";
}

inline void log_message(std::fstream &log_file, const std::string &message) {
  if (!log_file.is_open())
    return;
  log_file << message << "\n";
}

inline void log_time(std::fstream &log_file, const std::string &task,
                     double time_sec) {
  if (!log_file.is_open())
    return;
  log_file << "  - " << std::left << std::setw(40) << task << " : "
           << std::fixed << std::setprecision(3) << std::right << std::setw(10)
           << time_sec << " s\n";
}

inline void print_error(const std::string &error_msg) {
  std::cerr << "\n\033[1;31m[ERROR]\033[0m " << error_msg << "\n";
}

inline void print_warning(const std::string &warning_msg) {
  std::cout << "\n\033[1;33m[WARNING]\033[0m " << warning_msg << "\n";
}

inline void print_info(const std::string &info_msg) {
  std::cout << info_msg << "\n";
}

} // namespace TitanUI

#endif // TITAN_UI_H
