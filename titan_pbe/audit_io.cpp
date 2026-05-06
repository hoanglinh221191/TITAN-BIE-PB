#include "analyzer_impl.h"

#include "../titan_ui.h"

#include <cstddef>
#include <iomanip>
#include <sstream>

namespace TitanPBE {
namespace {

std::size_t audit_idx3(std::size_t i, std::size_t j, std::size_t k,
                       std::size_t nx, std::size_t ny) noexcept {
  return (k * ny + j) * nx + i;
}

std::string audit_format_suffix(std::size_t value, int width) {
  std::ostringstream oss;
  oss << std::setw(width) << std::setfill('0') << value;
  return oss.str();
}

} // namespace

bool AnalyzerImpl::should_write_map() const {
  return args_.map != "off" && args_.map_frames != "none";
}

std::string AnalyzerImpl::map_suffix(const FrameResult &result) const {
  if (args_.map_frames == "last")
    return "_last";
  return "_traj" + audit_format_suffix(result.trajectory_id, 4) + "_frame" +
         audit_format_suffix(result.frame_index, 6);
}

void AnalyzerImpl::maybe_write_maps(const GridShape &shape,
                                    const FrameResult &result,
                                    bool write_phi_now) {
  if (!should_write_map())
    return;
  if (args_.map_frames == "first" && wrote_first_map_)
    return;

  const std::string suffix = map_suffix(result);
  if (write_phi_now && (args_.map == "phi" || args_.map == "all")) {
    write_cube("pbe_phi" + suffix + ".cube", shape, "phi", false);
  }
  if (!write_phi_now && (args_.map == "epsilon" || args_.map == "all")) {
    write_cube("pbe_epsilon" + suffix + ".cube", shape, "epsilon", true);
  }
  if (!write_phi_now && (args_.map == "ion" || args_.map == "all")) {
    write_cube("pbe_ion" + suffix + ".cube", shape, "ion_accessible", false,
               true);
  }
  const bool completed_first =
      (write_phi_now && (args_.map == "phi" || args_.map == "all")) ||
      (!write_phi_now && (args_.map == "epsilon" || args_.map == "ion"));
  if (args_.map_frames == "first" && completed_first)
    wrote_first_map_ = true;
}

void AnalyzerImpl::write_cube(const std::string &filename,
                              const GridShape &shape,
                              const std::string &field_name,
                              bool epsilon_field, bool ion_field) const {
  std::ofstream out(filename, std::ios::out | std::ios::trunc);
  if (!out.is_open()) {
    TitanUI::print_warning("PBE: failed to write " + filename);
    return;
  }
  out << "TITAN PBE " << field_name << "\n";
  out << "units Angstrom; values are TITAN internal PBE units\n";
  out << std::setw(5) << 0 << " " << std::setw(12) << std::fixed
      << std::setprecision(6) << shape.ox << " " << std::setw(12) << shape.oy
      << " " << std::setw(12) << shape.oz << "\n";
  out << std::setw(5) << shape.nx << " " << std::setw(12) << shape.h << " "
      << std::setw(12) << 0.0 << " " << std::setw(12) << 0.0 << "\n";
  out << std::setw(5) << shape.ny << " " << std::setw(12) << 0.0 << " "
      << std::setw(12) << shape.h << " " << std::setw(12) << 0.0 << "\n";
  out << std::setw(5) << shape.nz << " " << std::setw(12) << 0.0 << " "
      << std::setw(12) << 0.0 << " " << std::setw(12) << shape.h << "\n";

  out << std::scientific << std::setprecision(6);
  int per_line = 0;
  for (std::size_t k = 0; k < shape.nz; ++k) {
    for (std::size_t j = 0; j < shape.ny; ++j) {
      for (std::size_t i = 0; i < shape.nx; ++i) {
        const std::size_t id = audit_idx3(i, j, k, shape.nx, shape.ny);
        const double value =
            ion_field
                ? (ion_accessible_[id] ? 1.0 : 0.0)
                : (epsilon_field ? (inside_[id] ? args_.int_diel
                                                : args_.ext_diel)
                                 : phi_[id]);
        out << std::setw(13) << value;
        if (++per_line == 6) {
          out << "\n";
          per_line = 0;
        }
      }
    }
  }
  if (per_line != 0)
    out << "\n";
}

} // namespace TitanPBE
