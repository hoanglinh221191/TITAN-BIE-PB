#include "titan_performance_observer.h"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <thread>

namespace fs = std::filesystem;

namespace {

void require_true(bool value, const std::string &message) {
  if (!value)
    throw std::runtime_error(message);
}

std::string read_text(const fs::path &path) {
  std::ifstream in(path);
  if (!in)
    throw std::runtime_error("failed to open " + path.string());
  return std::string(std::istreambuf_iterator<char>(in),
                     std::istreambuf_iterator<char>());
}

bool contains(const std::string &text, const std::string &needle) {
  return text.find(needle) != std::string::npos;
}

void remove_if_exists(const fs::path &path) {
  std::error_code ec;
  fs::remove(path, ec);
}

} // namespace

int main() {
  const fs::path root = fs::temp_directory_path() / "titan_perf_observer_smoke";
  fs::create_directories(root);
  const fs::path summary = root / "summary.dat";
  const fs::path modules = root / "modules.csv";
  const fs::path frames = root / "frames.csv";
  const fs::path log_path = root / "log.txt";
  remove_if_exists(summary);
  remove_if_exists(modules);
  remove_if_exists(frames);
  remove_if_exists(log_path);

  {
    std::fstream log(log_path, std::ios::out | std::ios::trunc);
    TitanPerformanceObserver observer;
    observer.configure("off");
    observer.add_frame_elapsed("CMap", 0, 0.0, 1.0);
    observer.add_counter("CMap", "contacts", 2.0);
    observer.export_outputs(summary.string(), modules.string(), frames.string(),
                            log);
    require_true(!fs::exists(summary) && !fs::exists(modules) &&
                     !fs::exists(frames),
                 "off mode exported files");
  }

  {
    std::fstream log(log_path, std::ios::out | std::ios::trunc);
    TitanPerformanceObserver observer;
    observer.configure("detail");
    observer.add_frame_elapsed("QNative", 2, 20.0, 3.0);
    observer.add_frame_elapsed("CMap", 1, 10.0, 2.0);
    observer.add_frame_elapsed("CMap", 2, 20.0, 4.0, "compute", false);
    observer.add_counter("CMap", "pairs", 8.0);
    observer.add_counter("CMap", "pairs", 2.0);
    observer.add_warning("CMap");
    observer.add_export_elapsed("QNative", 0.5);
    observer.begin_module_frame("Dihedral", 3, 30.0);
    std::this_thread::sleep_for(std::chrono::milliseconds(1));
    observer.end_module_frame("Dihedral");
    observer.export_outputs(summary.string(), modules.string(), frames.string(),
                            log);
  }

  require_true(fs::exists(summary), "missing summary");
  require_true(fs::exists(modules), "missing modules csv");
  require_true(fs::exists(frames), "missing detail frames csv");

  const std::string summary_text = read_text(summary);
  const std::string modules_text = read_text(modules);
  const std::string frames_text = read_text(frames);
  require_true(contains(summary_text, "mode: detail"),
               "summary does not record detail mode");
  require_true(contains(summary_text, "module_count: 3"),
               "summary does not count modules");
  require_true(contains(modules_text, "CMap,2,1,"),
               "CMap counters did not aggregate skipped frame correctly");
  require_true(contains(modules_text, "pairs=10.000000"),
               "counter aggregation failed");
  require_true(contains(modules_text, "QNative,1,1,"),
               "QNative module row missing");
  require_true(contains(modules_text, "Dihedral,1,1,"),
               "begin/end module row double counted or missing");

  const std::size_t cmap_pos = frames_text.find("CMap,1,");
  const std::size_t dihedral_pos = frames_text.find("Dihedral,3,");
  const std::size_t qnative_pos = frames_text.find("QNative,2,");
  require_true(cmap_pos != std::string::npos &&
                   dihedral_pos != std::string::npos &&
                   qnative_pos != std::string::npos,
               "detail rows missing expected modules");
  require_true(cmap_pos < dihedral_pos && dihedral_pos < qnative_pos,
               "detail rows are not deterministically sorted by module");

  std::cout << "PASS titan_performance_observer_smoke\n";
  return 0;
}
