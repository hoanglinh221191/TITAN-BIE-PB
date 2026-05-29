#include "titan_energy_decomposition.h"

#include "arguments.h"
#include "c6c12.h"
#include "charmm_nonbonded.h"
#include "martini_nonbonded.h"

#include <algorithm>
#include <cctype>
#include <cmath>

namespace {

std::string lowercase_copy(std::string value) {
  std::transform(value.begin(), value.end(), value.begin(),
                 [](unsigned char c) {
                   return static_cast<char>(std::tolower(c));
                 });
  return value;
}

martini::ElectrostaticConfig
martini_config_for_cutoff(const arguments &args, double cutoff_nm,
                          bool cutoff_active) {
  auto cfg = martini::electrostatic_config(args);
  if (cutoff_active && cutoff_nm > 0.0) {
    cfg.cutoff_nm = cutoff_nm;
    cfg.k_rf = 0.5 / (cfg.cutoff_nm * cfg.cutoff_nm * cfg.cutoff_nm);
    cfg.c_rf = 1.5 / cfg.cutoff_nm;
  }
  return cfg;
}

} // namespace

bool titan_energy_force_field_contains(const std::string &ff,
                                       const char *needle) {
  return lowercase_copy(ff).find(lowercase_copy(needle)) != std::string::npos;
}

double titan_energy_mic_delta(double d, double box_len) {
  if (box_len > 0.0 && std::isfinite(box_len))
    d -= std::nearbyint(d / box_len) * box_len;
  return d;
}

double titan_resolve_sg_cutoff_nm(const arguments &args) {
  if (args.sg_cutoff == "none")
    return -1.0;
  if (args.sg_cutoff == "value" && args.sg_cutoff_value_nm > 0.0)
    return args.sg_cutoff_value_nm;
  if (args.distance_energy_cutoff > 0.0)
    return args.distance_energy_cutoff;
  if (titan_energy_force_field_contains(args.force_f, "martini"))
    return martini::kDefaultCutoffNm;
  return -1.0;
}

void EnergyDecompositionProvider::configure(
    const Multi_MOL &rec, const Multi_MOL &lig,
    const cache_sigma_epsilon &cached_data, const arguments &args,
    double cutoff_nm, bool cutoff_active) {
  rec_ = &rec;
  lig_ = &lig;
  cached_data_ = &cached_data;
  args_ = &args;
  cutoff_nm_ = cutoff_nm;
  cutoff_active_ = cutoff_active && cutoff_nm > 0.0;
  force_field_ = cached_data.force_field_name.empty() ? args.force_f
                                                      : cached_data.force_field_name;
}

bool EnergyDecompositionProvider::ready() const noexcept {
  return rec_ != nullptr && lig_ != nullptr && cached_data_ != nullptr &&
         args_ != nullptr && has_force_parameters();
}

bool EnergyDecompositionProvider::has_force_parameters() const noexcept {
  return rec_ != nullptr && lig_ != nullptr && !rec_->parameters.empty() &&
         !lig_->parameters.empty();
}

EnergyPairDecomposition EnergyDecompositionProvider::pair_energy(
    const MOL &rec_frame, const MOL &lig_frame, std::size_t rec_pos,
    std::size_t lig_pos) const {
  if (!ready() || rec_ == nullptr || lig_ == nullptr)
    return {};
  return pair_energy_generic(rec_frame, *rec_, rec_pos, lig_frame, *lig_,
                             lig_pos);
}

EnergyPairDecomposition EnergyDecompositionProvider::rec_pair_energy(
    const MOL &rec_frame, std::size_t rec_pos_a, std::size_t rec_pos_b) const {
  if (cached_data_ == nullptr || args_ == nullptr || rec_ == nullptr)
    return {};
  return pair_energy_generic(rec_frame, *rec_, rec_pos_a, rec_frame, *rec_,
                             rec_pos_b);
}

EnergyPairDecomposition EnergyDecompositionProvider::lig_pair_energy(
    const MOL &lig_frame, std::size_t lig_pos_a, std::size_t lig_pos_b) const {
  if (cached_data_ == nullptr || args_ == nullptr || lig_ == nullptr)
    return {};
  return pair_energy_generic(lig_frame, *lig_, lig_pos_a, lig_frame, *lig_,
                             lig_pos_b);
}

const mol_info *EnergyDecompositionProvider::rec_info(
    std::size_t pos) const noexcept {
  if (rec_ == nullptr || pos >= rec_->info.size())
    return nullptr;
  return &rec_->info[pos];
}

const mol_info *EnergyDecompositionProvider::lig_info(
    std::size_t pos) const noexcept {
  if (lig_ == nullptr || pos >= lig_->info.size())
    return nullptr;
  return &lig_->info[pos];
}

EnergyPairDecomposition EnergyDecompositionProvider::pair_energy_generic(
    const MOL &frame_a, const Multi_MOL &mol_a, std::size_t pos_a,
    const MOL &frame_b, const Multi_MOL &mol_b, std::size_t pos_b) const {
  EnergyPairDecomposition out;
  if (cached_data_ == nullptr || args_ == nullptr || pos_a >= frame_a.size() ||
      pos_b >= frame_b.size() || pos_a >= mol_a.parameters.size() ||
      pos_b >= mol_b.parameters.size())
    return out;

  double dx = frame_a.pos.x[pos_a] - frame_b.pos.x[pos_b];
  double dy = frame_a.pos.y[pos_a] - frame_b.pos.y[pos_b];
  double dz = frame_a.pos.z[pos_a] - frame_b.pos.z[pos_b];
  dx = titan_energy_mic_delta(dx, frame_a.Box.get_a());
  dy = titan_energy_mic_delta(dy, frame_a.Box.get_b());
  dz = titan_energy_mic_delta(dz, frame_a.Box.get_c());
  const double d_nm = 0.1 * std::sqrt(dx * dx + dy * dy + dz * dz);
  if (d_nm <= 1e-10)
    return out;
  if (cutoff_active_ && cutoff_nm_ > 0.0 && d_nm > cutoff_nm_)
    return out;

  if (titan_energy_force_field_contains(force_field_, "martini")) {
    double sigma_ij{}, epsilon_ij{};
    martini::pair_params(*cached_data_, mol_a.parameters, mol_b.parameters, pos_a,
                         pos_b, sigma_ij, epsilon_ij);
    if (martini::excluded_pair(*args_, mol_a.info, mol_b.info, pos_a, pos_b,
                               d_nm, sigma_ij))
      return out;
    const auto cfg =
        martini_config_for_cutoff(*args_, cutoff_nm_, cutoff_active_);
    out.elec = martini::electro_energy(cfg, mol_a.parameters.charge[pos_a],
                                       mol_b.parameters.charge[pos_b], d_nm);
    out.vdw = martini::lj_energy(sigma_ij, epsilon_ij, d_nm);
    return out;
  }

  const double qi = mol_a.parameters.charge[pos_a];
  const double qj = mol_b.parameters.charge[pos_b];
  out.elec = 138.935457839 * qi * qj /
             (std::max(1e-12, args_->epsilon_in) * d_nm);

  if (titan_energy_force_field_contains(force_field_, "charmm")) {
    double sigma_ij{}, epsilon_ij{};
    charmm::pair_params(*cached_data_, mol_a.parameters, mol_b.parameters, pos_a,
                        pos_b, sigma_ij, epsilon_ij);
    if (args_->self_int &&
        ((sigma_ij > 0.0 && d_nm < sigma_ij) ||
         (args_->flag_interchain_energy &&
          mol_b.info[pos_b].get_chain() == mol_a.info[pos_a].get_chain())))
      return {};
    out.vdw = charmm::lj_energy(sigma_ij, epsilon_ij, d_nm);
  } else {
    double sigma_ij =
        0.5 * (mol_b.parameters.sigma[pos_b] + mol_a.parameters.sigma[pos_a]);
    const double epsilon_ij = std::sqrt(std::max(
        0.0, mol_b.parameters.epsilon[pos_b] * mol_a.parameters.epsilon[pos_a]));
    if (titan_energy_force_field_contains(force_field_, "opls"))
      sigma_ij = std::sqrt(std::max(
          0.0, mol_b.parameters.sigma[pos_b] * mol_a.parameters.sigma[pos_a]));
    double sr6 = 0.0;
    if (mol_b.parameters.sigma[pos_b] >= 0.0 &&
        mol_a.parameters.sigma[pos_a] >= 0.0)
      sr6 = c6(sigma_ij / d_nm);
    if (args_->self_int &&
        (sr6 > 1.0 ||
         (args_->flag_interchain_energy &&
          mol_b.info[pos_b].get_chain() == mol_a.info[pos_a].get_chain())))
      return {};
    out.vdw = epsilon_ij * (square(sr6) - sr6);
  }
  return out;
}
