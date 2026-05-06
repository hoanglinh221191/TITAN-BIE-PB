#ifndef TITAN_PBE_H
#define TITAN_PBE_H

#include "MOL.h"
#include "MOL_para.h"
#include "arguments.h"
#include "titan_pbe/types.h"

#include <cstddef>
#include <fstream>
#include <memory>
#include <string>

class TitanExecutor;
class EnergyDecompositionProvider;

namespace TitanPBE {

class AnalyzerImpl;

class Analyzer {
public:
  Analyzer();
  ~Analyzer();

  Analyzer(Analyzer &&) noexcept;
  Analyzer &operator=(Analyzer &&) noexcept;
  Analyzer(const Analyzer &) = delete;
  Analyzer &operator=(const Analyzer &) = delete;

  void configure(const arguments &args);
  [[nodiscard]] bool enabled() const noexcept;

  void begin_trajectory(std::size_t trajectory_id,
                        const std::string &source_label);
  void end_trajectory(std::fstream &log);

  FrameResult consume_frame(const MOL &mol1, const Multi_MOL_parameters &params1,
                            const MOL *mol2,
                            const Multi_MOL_parameters *params2,
                            std::size_t trajectory_id,
                            std::size_t frame_index,
                            std::size_t global_index,
                            TitanExecutor &executor, std::fstream &log,
                            const EnergyDecompositionProvider *energy_provider =
                                nullptr);

  void finalize(std::fstream &log);

private:
  std::unique_ptr<AnalyzerImpl> impl_;
};

} // namespace TitanPBE

#endif