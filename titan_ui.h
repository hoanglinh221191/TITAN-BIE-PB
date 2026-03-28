#ifndef TITAN_UI_H
#define TITAN_UI_H

#include <fstream>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <sstream>
#include <string>

namespace TitanUI {

inline std::string clamp_display_text(const std::string &text,
                                      std::size_t max_width) {
  if (text.size() <= max_width)
    return text;
  if (max_width <= 3)
    return text.substr(0, max_width);
  return text.substr(0, max_width - 3) + "...";
}

inline void print_runtime_module_summary(const std::string &requested,
                                         const std::string &active,
                                         const std::string &skipped = {},
                                         const std::string &streaming = {}) {
  if (requested.empty() && active.empty() && skipped.empty() &&
      streaming.empty())
    return;

  std::cout << "\n\033[1;36m[TITAN Runtime]\033[0m\n";
  if (!requested.empty())
    std::cout << "  \033[1;33mRequested\033[0m : " << requested << "\n";
  if (!active.empty())
    std::cout << "  \033[1;32mActive   \033[0m : " << active << "\n";
  if (!streaming.empty())
    std::cout << "  \033[1;36mStreaming\033[0m : " << streaming << "\n";
  if (!skipped.empty())
    std::cout << "  \033[1;90mSkipped  \033[0m : " << skipped << "\n";
  std::cout << "\n";
}

inline void write_logo_art(std::ostream &out) {
  out << "                                                  /\\\\ \n"
      << "                                                 //\\\\\\\\ \n"
      << "                                                ||--|| \n"
      << "  ████████╗██╗████████╗ █████╗ ███╗   ██╗       ||  || \n"
      << "  ╚══██╔══╝██║╚══██╔══╝██╔══██╗████╗  ██║      _||  ||_ \n"
      << "     ██║   ██║   ██║   ███████║██╔██╗ ██║     /\\\\||  ||/\\\\ \n"
      << "     ██║   ██║   ██║   ██╔══██║██║╚██╗██║    // ||__|| \\\\\\\\ \n"
      << "     ██║   ██║   ██║   ██║  ██║██║ ╚████║   |/  |____|  \\\\| \n"
      << "     ╚═╝   ╚═╝   ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═══╝       / || \\\\ \n"
      << "                                               /  ||  \\\\ \n"
      << "                                              |/  \\\\/  \\\\| \n"
      << "                                                  \\\\/ \n";
}

inline void write_logo_art_log(std::ostream &out) {
  out << "                                                     /\\\\\n"
      << "                                                    /  \\\\\n"
      << "                                                   /++++\\\\\n"
      << "                                                  /  ()  \\\\\n"
      << "  _______ ___ _______   ___   _   _               |      |\n"
      << " |__   __|_ _|__   __| / _ \\ | \\ | |              |  ||  |\n"
      << "    | |   | |   | |   / /_\\ \\|  \\| |             /|__||__|\\\n"
      << "    | |   | |   | |   |  _  || . ` |            /_/_/  \\_\\_\\\n"
      << "    | |  _| |_  | |   | | | || |\\  |               /_/\\/\\_\\\n"
      << "    |_| |_____| |_|   |_| |_||_| \\_|                /_/  \\_\\\n";
}

inline void show_progress(const std::string &module_name, size_t current,
                          size_t total, const std::string &status_text = "",
                          const std::string &module_text = "") {
  if (total == 0)
    return;

  // current ranges from 0 to total - 1. Add 1 to make it reach 100% at the last
  // step.
  float progress = static_cast<float>(current + 1) / total;
  int barWidth = 28;
  const std::string display_name = clamp_display_text(module_name, 18);
  const std::string display_status = clamp_display_text(status_text, 36);
  const std::string display_modules = clamp_display_text(module_text, 42);
  std::cout << "\r\033[1;36m[TITAN]\033[0m \033[1;33m" << std::setw(18)
            << std::left << display_name << "\033[0m [";

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
            << "%";
  if (!display_status.empty())
    std::cout << " | " << display_status;
  if (!display_modules.empty())
    std::cout << " | Active: " << display_modules;
  std::cout << "\033[K" << std::flush;

  if (current + 1 >= total) {
    std::cout << std::endl;
  }
}

inline void show_continuous_progress(const std::string &module_name,
                                     size_t current_frame,
                                     const std::string &status_text = "",
                                     const std::string &module_text = "") {
  const std::string display_name = clamp_display_text(module_name, 18);
  const std::string display_status = clamp_display_text(status_text, 44);
  const std::string display_modules = clamp_display_text(module_text, 40);
  std::ostringstream oss;
  oss << "Frame " << current_frame;
  if (!display_status.empty())
    oss << " | " << display_status;

  std::cout << "\r\033[1;36m[TITAN]\033[0m \033[1;33m" << std::setw(18)
            << std::left << display_name << "\033[0m "
            << "\033[1;32m" << oss.str() << "\033[0m";
  if (!display_modules.empty())
    std::cout << " | Active: " << display_modules;
  std::cout << "\033[K" << std::flush;
}

inline void print_module_header(const std::string &header_name) {
  std::cout << "\n\033[1;35m>>> " << header_name << " <<<\033[0m\n";
}

inline void print_logo() {
  std::cout << "\n\033[1;36m";
  write_logo_art(std::cout);
  std::cout << "\033[0m";
  std::cout
      << "  \033[1;33mThermodynamic & Interaction Trajectory ANalyzer\033[0m\n";
  std::cout
      << "  "
         "\033[1;32m===============================================\033[0m\n\n";
}

inline void log_logo(std::fstream &log_file) {
  if (!log_file.is_open())
    return;
  log_file << "\n";
  write_logo_art_log(log_file);
  log_file << "  Thermodynamic & Interaction Trajectory ANalyzer\n";
  log_file << "  ===============================================\n\n";
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
