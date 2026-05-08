#include "tabi_types.h"

#include "../titan_executor.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <sstream>
#include <vector>

namespace TitanPBE::Tabi {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kCoulombKcalPerMol = 332.063713299;
constexpr double kKcalToKJ = 4.184;
constexpr double kKBKJPerMolK = 0.00831446261815324;

struct FormulationConfig {
  std::string formulation{"jcp2013"};
  std::string energy_projection{"boundary"};
  double rhs_scale{1.0};
  double operator_contrast_scale{1.0};
  double boundary_energy_scale{1.0};
};

double debye_kappa_inv_a(double salt_m, double temp_k) noexcept {
  if (salt_m <= 0.0)
    return 0.0;
  const double temp_scale = std::sqrt(298.15 / std::max(1.0, temp_k));
  return 0.3291 * std::sqrt(salt_m) * temp_scale;
}

double dot(const std::vector<double> &a, const std::vector<double> &b) {
  double value = 0.0;
  for (std::size_t i = 0; i < a.size(); ++i)
    value += a[i] * b[i];
  return value;
}

double l2_norm(const std::vector<double> &values) {
  return std::sqrt(std::max(0.0, dot(values, values)));
}

double rms_norm(const std::vector<double> &values) {
  if (values.empty())
    return 0.0;
  return l2_norm(values) / std::sqrt(static_cast<double>(values.size()));
}

double elapsed_ms(std::chrono::steady_clock::time_point begin,
                  std::chrono::steady_clock::time_point end) noexcept {
  return std::chrono::duration<double, std::milli>(end - begin).count();
}

std::string residual_history_string(const std::vector<double> &history) {
  std::ostringstream oss;
  oss.setf(std::ios::scientific);
  oss.precision(6);
  for (std::size_t i = 0; i < history.size(); ++i) {
    if (i != 0)
      oss << ';';
    oss << history[i];
  }
  return oss.str();
}

FormulationConfig make_formulation_config(const Options &options) {
  FormulationConfig config;
  config.formulation = options.formulation;
  config.energy_projection = options.energy_projection;
  if (config.formulation.empty())
    config.formulation = "jcp2013";
  if (config.energy_projection.empty())
    config.energy_projection = "boundary";
  if (config.formulation == "jcp2013") {
    // The jcp2013 path is intentionally explicit even while it shares the
    // current scale factors. Accuracy gates use the audit CSV to keep future
    // sign/scale changes visible instead of silently drifting the formulation.
    config.rhs_scale = 1.0;
    config.operator_contrast_scale = 1.0;
    config.boundary_energy_scale = 1.0;
  } else if (config.formulation == "calibrated-debug") {
    config.rhs_scale = 1.0;
    config.operator_contrast_scale = 1.0;
    config.boundary_energy_scale = 1.0;
  }
  return config;
}

struct LinearSolveResult {
  std::vector<double> x;
  int iterations{0};
  double residual{0.0};
  bool converged{false};
  std::string convergence_reason{"not_started"};
  std::vector<double> residual_history;
  double orthogonalization_ms{0.0};
  double preconditioner_ms{0.0};
};

LinearSolveResult gmres_solve(
    const std::function<void(const std::vector<double> &,
                             std::vector<double> &)> &matvec,
    const std::vector<double> &b, const std::vector<double> &diagonal,
    int max_iter, double tol) {
  LinearSolveResult result;
  const std::size_t n = b.size();
  result.x.assign(n, 0.0);
  if (n == 0) {
    result.converged = true;
    result.convergence_reason = "converged";
    return result;
  }

  std::vector<double> inv_precond(n, 1.0);
  const auto precond_begin = std::chrono::steady_clock::now();
  bool preconditioner_ok = true;
  for (std::size_t i = 0; i < n; ++i) {
    const double diag = i < diagonal.size() ? std::abs(diagonal[i]) : 1.0;
    if (!std::isfinite(diag) || diag <= 1.0e-12) {
      inv_precond[i] = 1.0;
      preconditioner_ok = false;
    } else {
      inv_precond[i] = 1.0 / diag;
    }
  }
  result.preconditioner_ms +=
      elapsed_ms(precond_begin, std::chrono::steady_clock::now());
  if (!preconditioner_ok)
    result.convergence_reason = "preconditioner_breakdown";

  auto apply_preconditioner = [&](const std::vector<double> &in,
                                  std::vector<double> &out) {
    const auto begin = std::chrono::steady_clock::now();
    out.resize(n);
    for (std::size_t i = 0; i < n; ++i)
      out[i] = in[i] * inv_precond[i];
    result.preconditioner_ms +=
        elapsed_ms(begin, std::chrono::steady_clock::now());
  };

  auto finite_vector = [](const std::vector<double> &values) noexcept {
    for (double value : values) {
      if (!std::isfinite(value))
        return false;
    }
    return true;
  };

  const int restart =
      std::max(4, std::min<int>(32, static_cast<int>(std::min<std::size_t>(
                                     n, static_cast<std::size_t>(32)))));
  std::vector<double> ax(n, 0.0);
  std::vector<double> residual = b;
  std::vector<double> precond_residual(n, 0.0);
  apply_preconditioner(residual, precond_residual);
  result.residual = rms_norm(residual);
  double precond_rms = rms_norm(precond_residual);
  result.residual_history.push_back(result.residual);
  if (result.residual <= tol) {
    result.converged = true;
    result.convergence_reason = "converged";
    return result;
  }

  const int iter_cap = std::max(1, max_iter);
  double best_physical_rms = result.residual;
  while (result.iterations < iter_cap) {
    const double beta = l2_norm(precond_residual);
    if (beta <= tol * std::sqrt(static_cast<double>(n))) {
      result.converged = true;
      result.residual = rms_norm(residual);
      result.convergence_reason = "converged";
      break;
    }

    const int window = std::min(restart, iter_cap - result.iterations);
    std::vector<std::vector<double>> v(
        static_cast<std::size_t>(window) + 1, std::vector<double>(n, 0.0));
    for (std::size_t i = 0; i < n; ++i)
      v[0][i] = precond_residual[i] / beta;
    std::vector<double> h((static_cast<std::size_t>(window) + 1) *
                              static_cast<std::size_t>(window),
                          0.0);
    std::vector<double> cs(static_cast<std::size_t>(window), 0.0);
    std::vector<double> sn(static_cast<std::size_t>(window), 0.0);
    std::vector<double> g(static_cast<std::size_t>(window) + 1, 0.0);
    g[0] = beta;

    int used_cols = 0;
    bool happy_breakdown = false;
    for (int j = 0; j < window; ++j) {
      std::vector<double> w(n, 0.0);
      matvec(v[static_cast<std::size_t>(j)], w);
      if (!finite_vector(w)) {
        result.convergence_reason = "nan_guard";
        return result;
      }
      std::vector<double> precond_w(n, 0.0);
      apply_preconditioner(w, precond_w);
      w.swap(precond_w);

      const auto ortho_begin = std::chrono::steady_clock::now();
      for (int i = 0; i <= j; ++i) {
        const double hij = dot(w, v[static_cast<std::size_t>(i)]);
        h[static_cast<std::size_t>(i) * static_cast<std::size_t>(window) +
          static_cast<std::size_t>(j)] = hij;
        for (std::size_t k = 0; k < n; ++k)
          w[k] -= hij * v[static_cast<std::size_t>(i)][k];
      }
      const double h_next = l2_norm(w);
      h[(static_cast<std::size_t>(j) + 1) *
            static_cast<std::size_t>(window) +
        static_cast<std::size_t>(j)] = h_next;
      if (h_next > 1.0e-20 && j + 1 < window + 1) {
        for (std::size_t k = 0; k < n; ++k)
          v[static_cast<std::size_t>(j) + 1][k] = w[k] / h_next;
      } else {
        happy_breakdown = true;
      }
      result.orthogonalization_ms +=
          elapsed_ms(ortho_begin, std::chrono::steady_clock::now());

      for (int i = 0; i < j; ++i) {
        const double h_i =
            h[static_cast<std::size_t>(i) * static_cast<std::size_t>(window) +
              static_cast<std::size_t>(j)];
        const double h_next_i =
            h[(static_cast<std::size_t>(i) + 1) *
                  static_cast<std::size_t>(window) +
              static_cast<std::size_t>(j)];
        h[static_cast<std::size_t>(i) * static_cast<std::size_t>(window) +
          static_cast<std::size_t>(j)] = cs[static_cast<std::size_t>(i)] * h_i +
                                         sn[static_cast<std::size_t>(i)] *
                                             h_next_i;
        h[(static_cast<std::size_t>(i) + 1) *
              static_cast<std::size_t>(window) +
          static_cast<std::size_t>(j)] = -sn[static_cast<std::size_t>(i)] * h_i +
                                         cs[static_cast<std::size_t>(i)] *
                                             h_next_i;
      }

      const double h_j =
          h[static_cast<std::size_t>(j) * static_cast<std::size_t>(window) +
            static_cast<std::size_t>(j)];
      const double h_next_j =
          h[(static_cast<std::size_t>(j) + 1) *
                static_cast<std::size_t>(window) +
            static_cast<std::size_t>(j)];
      const double rho = std::hypot(h_j, h_next_j);
      if (rho <= 1.0e-30 || !std::isfinite(rho)) {
        happy_breakdown = true;
        h[static_cast<std::size_t>(j) * static_cast<std::size_t>(window) +
          static_cast<std::size_t>(j)] = h_j;
        h[(static_cast<std::size_t>(j) + 1) *
              static_cast<std::size_t>(window) +
          static_cast<std::size_t>(j)] = 0.0;
      } else {
        cs[static_cast<std::size_t>(j)] = h_j / rho;
        sn[static_cast<std::size_t>(j)] = h_next_j / rho;
        h[static_cast<std::size_t>(j) * static_cast<std::size_t>(window) +
          static_cast<std::size_t>(j)] = rho;
        h[(static_cast<std::size_t>(j) + 1) *
              static_cast<std::size_t>(window) +
          static_cast<std::size_t>(j)] = 0.0;
        const double g_j = g[static_cast<std::size_t>(j)];
        const double g_next = g[static_cast<std::size_t>(j) + 1];
        g[static_cast<std::size_t>(j)] =
            cs[static_cast<std::size_t>(j)] * g_j +
            sn[static_cast<std::size_t>(j)] * g_next;
        g[static_cast<std::size_t>(j) + 1] =
            -sn[static_cast<std::size_t>(j)] * g_j +
            cs[static_cast<std::size_t>(j)] * g_next;
      }

      ++result.iterations;
      used_cols = j + 1;
      precond_rms =
          std::abs(g[static_cast<std::size_t>(j) + 1]) /
          std::sqrt(static_cast<double>(n));
      result.residual_history.push_back(precond_rms);
      if (precond_rms <= tol) {
        result.converged = true;
        result.convergence_reason = "converged";
        break;
      }
      if (result.iterations >= iter_cap)
        break;
      if (happy_breakdown)
        break;
    }

    if (used_cols == 0)
      break;

    std::vector<double> y(static_cast<std::size_t>(used_cols), 0.0);
    for (int rev = 0; rev < used_cols; ++rev) {
      const int row = used_cols - 1 - rev;
      double value = g[static_cast<std::size_t>(row)];
      for (int col = row + 1; col < used_cols; ++col) {
        value -= h[static_cast<std::size_t>(row) *
                       static_cast<std::size_t>(window) +
                   static_cast<std::size_t>(col)] *
                 y[static_cast<std::size_t>(col)];
      }
      const double diag =
          h[static_cast<std::size_t>(row) * static_cast<std::size_t>(window) +
            static_cast<std::size_t>(row)];
      y[static_cast<std::size_t>(row)] =
          std::abs(diag) > 1.0e-30 ? value / diag : 0.0;
    }
    for (int basis = 0; basis < used_cols; ++basis) {
      for (std::size_t k = 0; k < n; ++k)
        result.x[k] += y[static_cast<std::size_t>(basis)] *
                       v[static_cast<std::size_t>(basis)][k];
    }

    matvec(result.x, ax);
    if (!finite_vector(ax)) {
      result.convergence_reason = "nan_guard";
      return result;
    }
    for (std::size_t k = 0; k < n; ++k)
      residual[k] = b[k] - ax[k];
    result.residual = rms_norm(residual);
    result.residual_history.push_back(result.residual);
    apply_preconditioner(residual, precond_residual);
    precond_rms = rms_norm(precond_residual);
    if (result.residual <= tol) {
      result.converged = true;
      result.convergence_reason = "converged";
      break;
    }
    if (result.residual >= best_physical_rms * 0.999999) {
      result.convergence_reason =
          result.iterations >= iter_cap ? "max_iter" : "stagnated";
      break;
    }
    best_physical_rms = result.residual;
  }
  if (!result.converged && result.convergence_reason == "not_started")
    result.convergence_reason =
        result.iterations >= iter_cap ? "max_iter" : "stagnated";
  return result;
}

// TABI v1 convention: linear boundary-integral solve for the induced
// surface density. Kernels use Angstrom distances, elementary-charge inputs,
// inverse-Angstrom kappa, outward normals, and Coulomb conversion only in the
// final reaction-field energy projection.
std::vector<double> build_rhs(const std::vector<Atom> &atoms,
                              const std::vector<Panel> &panels,
                              double contrast, double kappa,
                              double rhs_scale) {
  std::vector<double> rhs(panels.size(), 0.0);
  for (std::size_t p = 0; p < panels.size(); ++p) {
    const Panel &panel = panels[p];
    double normal_field = 0.0;
    for (const Atom &atom : atoms) {
      const double dx = panel.x - atom.x;
      const double dy = panel.y - atom.y;
      const double dz = panel.z - atom.z;
      const double r2 = dx * dx + dy * dy + dz * dz;
      if (r2 <= 1.0e-16)
        continue;
      const double r = std::sqrt(r2);
      const double normal_projection =
          (panel.nx * dx + panel.ny * dy + panel.nz * dz) / r;
      const double screened = std::exp(-kappa * r);
      normal_field += atom.q * screened * normal_projection /
                      std::max(1.0e-12, r2);
    }
    rhs[p] = -rhs_scale * contrast * normal_field / (4.0 * kPi);
  }
  return rhs;
}

double boundary_energy_kcal(const std::vector<Atom> &atoms,
                            const std::vector<Panel> &panels,
                            const std::vector<double> &density,
                            double kappa) {
  double energy = 0.0;
  for (const Atom &atom : atoms) {
    double phi = 0.0;
    for (std::size_t p = 0; p < panels.size(); ++p) {
      const Panel &panel = panels[p];
      const double dx = atom.x - panel.x;
      const double dy = atom.y - panel.y;
      const double dz = atom.z - panel.z;
      const double r2 = dx * dx + dy * dy + dz * dz;
      if (r2 <= 1.0e-16)
        continue;
      const double r = std::sqrt(r2);
      phi += density[p] * panel.area * std::exp(-kappa * r) /
             (4.0 * kPi * r);
    }
    energy += 0.5 * kCoulombKcalPerMol * atom.q * phi;
  }
  return energy;
}

double born_reference_kcal(const std::vector<Atom> &atoms, double eps_in,
                           double eps_out, double kappa) {
  double energy = 0.0;
  for (const Atom &atom : atoms) {
    const double radius = std::max(0.5, atom.radius);
    const double outside =
        1.0 / (std::max(1.0e-12, eps_out) * (1.0 + kappa * radius));
    const double inside = 1.0 / std::max(1.0e-12, eps_in);
    energy +=
        0.5 * kCoulombKcalPerMol * atom.q * atom.q * (outside - inside) /
        radius;
  }
  return energy;
}

double green_kernel(double r, double kappa) noexcept {
  if (r <= 1.0e-14)
    return 0.0;
  return std::exp(-kappa * r) / (4.0 * kPi * r);
}

double source_normal_derivative_green(double tx, double ty, double tz,
                                      const Panel &source,
                                      double kappa) noexcept {
  const double ux = tx - source.x;
  const double uy = ty - source.y;
  const double uz = tz - source.z;
  const double r2 = ux * ux + uy * uy + uz * uz;
  if (r2 <= 1.0e-28)
    return 0.0;
  const double r = std::sqrt(r2);
  const double ndotu = source.nx * ux + source.ny * uy + source.nz * uz;
  return std::exp(-kappa * r) * (1.0 + kappa * r) * ndotu /
         (4.0 * kPi * r2 * r);
}

double target_normal_derivative_green(const Panel &target, double sx,
                                      double sy, double sz,
                                      double kappa) noexcept {
  const double ux = target.x - sx;
  const double uy = target.y - sy;
  const double uz = target.z - sz;
  const double r2 = ux * ux + uy * uy + uz * uz;
  if (r2 <= 1.0e-28)
    return 0.0;
  const double r = std::sqrt(r2);
  const double ndotu = target.nx * ux + target.ny * uy + target.nz * uz;
  return -std::exp(-kappa * r) * (1.0 + kappa * r) * ndotu /
         (4.0 * kPi * r2 * r);
}

double mixed_normal_derivative_green(const Panel &target, const Panel &source,
                                     double kappa) noexcept {
  const double ux = target.x - source.x;
  const double uy = target.y - source.y;
  const double uz = target.z - source.z;
  const double r2 = ux * ux + uy * uy + uz * uz;
  if (r2 <= 1.0e-28)
    return 0.0;
  const double r = std::sqrt(r2);
  const double inv_r = 1.0 / r;
  const double nxu =
      (target.nx * ux + target.ny * uy + target.nz * uz) * inv_r;
  const double nyu =
      (source.nx * ux + source.ny * uy + source.nz * uz) * inv_r;
  const double nxny = target.nx * source.nx + target.ny * source.ny +
                      target.nz * source.nz;
  const double expkr = std::exp(-kappa * r);
  const double inv_r2 = inv_r * inv_r;
  const double fprime =
      -expkr * (1.0 + kappa * r) * inv_r2 / (4.0 * kPi);
  const double fsecond =
      expkr * (kappa * kappa * inv_r + 2.0 * kappa * inv_r2 +
               2.0 * inv_r2 * inv_r) /
      (4.0 * kPi);
  const double hessian_contract =
      fsecond * nxu * nyu +
      (fprime * inv_r) * (nxny - nxu * nyu);
  return hessian_contract;
}

struct JcpKernels {
  double k1{0.0};
  double k2{0.0};
  double k3{0.0};
  double k4{0.0};
};

JcpKernels jcp_surface_kernels(const Panel &target, const Panel &source,
                               double eps_ratio, double kappa) noexcept {
  const double dx = target.x - source.x;
  const double dy = target.y - source.y;
  const double dz = target.z - source.z;
  const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
  if (r <= 1.0e-14)
    return {};
  const double g0 = green_kernel(r, 0.0);
  const double gk = green_kernel(r, kappa);
  const double d0_src =
      source_normal_derivative_green(target.x, target.y, target.z, source, 0.0);
  const double dk_src = source_normal_derivative_green(
      target.x, target.y, target.z, source, kappa);
  const double d0_tgt = target_normal_derivative_green(
      target, source.x, source.y, source.z, 0.0);
  const double dk_tgt = target_normal_derivative_green(
      target, source.x, source.y, source.z, kappa);
  const double d2_0 = mixed_normal_derivative_green(target, source, 0.0);
  const double d2_k = mixed_normal_derivative_green(target, source, kappa);
  const double inv_eps_ratio = 1.0 / std::max(1.0e-12, eps_ratio);
  return {g0 - gk, eps_ratio * dk_src - d0_src,
          d0_tgt - inv_eps_ratio * dk_tgt, d2_k - d2_0};
}

std::vector<double> build_jcp_rhs(const std::vector<Atom> &atoms,
                                  const std::vector<Panel> &panels,
                                  double eps_in, bool coupled,
                                  double *s1_norm, double *s2_norm) {
  const std::size_t n = panels.size();
  std::vector<double> rhs(coupled ? 2 * n : n, 0.0);
  const double inv_eps_in = 1.0 / std::max(1.0e-12, eps_in);
  for (std::size_t i = 0; i < n; ++i) {
    const Panel &panel = panels[i];
    double s1 = 0.0;
    double s2 = 0.0;
    for (const Atom &atom : atoms) {
      const double dx = panel.x - atom.x;
      const double dy = panel.y - atom.y;
      const double dz = panel.z - atom.z;
      const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
      if (r <= 1.0e-14)
        continue;
      s1 += atom.q * inv_eps_in * green_kernel(r, 0.0);
      s2 += atom.q * inv_eps_in *
            target_normal_derivative_green(panel, atom.x, atom.y, atom.z, 0.0);
    }
    rhs[i] = s1;
    if (coupled)
      rhs[n + i] = s2;
  }
  if (s1_norm)
    *s1_norm = rms_norm(std::vector<double>(rhs.begin(), rhs.begin() + n));
  if (s2_norm) {
    if (coupled)
      *s2_norm = rms_norm(std::vector<double>(rhs.begin() + n, rhs.end()));
    else
      *s2_norm = 0.0;
  }
  return rhs;
}

void apply_jcp_linear_pb_operator(const std::vector<Panel> &panels,
                                  const std::vector<double> &x,
                                  double eps_ratio, double kappa,
                                  double lhs_phi, double lhs_phin,
                                  TitanExecutor &executor,
                                  std::vector<double> &out) {
  const std::size_t n = panels.size();
  out.assign(2 * n, 0.0);
  const std::size_t workers = std::max<std::size_t>(
      1, std::min<std::size_t>(executor.capacity(), n));
  executor.parallel_for(
      n, workers, 4,
      [&](std::size_t, std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; ++i) {
          const Panel &target = panels[i];
          double phi_row = lhs_phi * x[i];
          double phin_row = lhs_phin * x[n + i];
          for (std::size_t j = 0; j < n; ++j) {
            if (j == i)
              continue;
            const Panel &source = panels[j];
            const JcpKernels kernel =
                jcp_surface_kernels(target, source, eps_ratio, kappa);
            const double area = source.area;
            phi_row -= area * (kernel.k1 * x[n + j] + kernel.k2 * x[j]);
            phin_row -= area * (kernel.k3 * x[n + j] + kernel.k4 * x[j]);
          }
          out[i] = phi_row;
          out[n + i] = phin_row;
        }
      });
}

void apply_jcp_poisson_operator(const std::vector<Panel> &panels,
                                const std::vector<double> &x,
                                double eps_ratio, double lhs_phi,
                                TitanExecutor &executor,
                                std::vector<double> &out) {
  const std::size_t n = panels.size();
  out.assign(n, 0.0);
  const double k2_scale = eps_ratio - 1.0;
  const std::size_t workers = std::max<std::size_t>(
      1, std::min<std::size_t>(executor.capacity(), n));
  executor.parallel_for(
      n, workers, 4,
      [&](std::size_t, std::size_t begin, std::size_t end) {
        for (std::size_t i = begin; i < end; ++i) {
          const Panel &target = panels[i];
          double row = lhs_phi * x[i];
          for (std::size_t j = 0; j < n; ++j) {
            if (j == i)
              continue;
            const Panel &source = panels[j];
            const double k2 =
                k2_scale * source_normal_derivative_green(
                               target.x, target.y, target.z, source, 0.0);
            row -= source.area * k2 * x[j];
          }
          out[i] = row;
        }
      });
}

double jcp_linear_pb_energy_kcal(const std::vector<Atom> &atoms,
                                 const std::vector<Panel> &panels,
                                 const std::vector<double> &x,
                                 double eps_ratio, double kappa) {
  const std::size_t n = panels.size();
  double energy = 0.0;
  for (const Atom &atom : atoms) {
    double reaction_phi = 0.0;
    for (std::size_t j = 0; j < n; ++j) {
      const Panel &source = panels[j];
      const double dx = atom.x - source.x;
      const double dy = atom.y - source.y;
      const double dz = atom.z - source.z;
      const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
      if (r <= 1.0e-14)
        continue;
      const double k1 = green_kernel(r, 0.0) - green_kernel(r, kappa);
      const double k2 =
          eps_ratio *
              source_normal_derivative_green(atom.x, atom.y, atom.z, source,
                                             kappa) -
          source_normal_derivative_green(atom.x, atom.y, atom.z, source, 0.0);
      reaction_phi += source.area * (k1 * x[n + j] + k2 * x[j]);
    }
    energy += 0.5 * kCoulombKcalPerMol * atom.q * reaction_phi;
  }
  return energy;
}

double jcp_poisson_energy_kcal(const std::vector<Atom> &atoms,
                               const std::vector<Panel> &panels,
                               const std::vector<double> &phi,
                               double eps_ratio) {
  const double k2_scale = eps_ratio - 1.0;
  double energy = 0.0;
  for (const Atom &atom : atoms) {
    double reaction_phi = 0.0;
    for (std::size_t j = 0; j < panels.size(); ++j) {
      const Panel &source = panels[j];
      const double k2 =
          k2_scale * source_normal_derivative_green(atom.x, atom.y, atom.z,
                                                    source, 0.0);
      reaction_phi += source.area * k2 * phi[j];
    }
    energy += 0.5 * kCoulombKcalPerMol * atom.q * reaction_phi;
  }
  return energy;
}

} // namespace

Result solve_linear_pb(const std::vector<Atom> &atoms, const Options &options,
                       TitanExecutor &executor) {
  Result result;
  const FormulationConfig formulation = make_formulation_config(options);
  result.formulation = formulation.formulation;
  result.energy_projection = formulation.energy_projection;
  result.rhs_scale = formulation.rhs_scale;
  result.coulomb_factor_kcal_mol = kCoulombKcalPerMol;
  result.total_charge = std::accumulate(
      atoms.begin(), atoms.end(), 0.0,
      [](double sum, const Atom &atom) { return sum + atom.q; });

  const std::vector<Panel> panels =
      build_surface(atoms, options, &result, &executor);
  result.panel_count = panels.size();
  result.boundary_nodes = panels.size();
  result.taylor_order = std::clamp(options.taylor_order, 3, 4);
  result.memory_bytes =
      panels.size() * sizeof(Panel) + panels.size() * sizeof(double) * 8;

  if (panels.empty()) {
    result.converged = true;
    result.status = "ok";
    result.reaction_energy_source = "tabi_boundary_integral";
    result.convergence_reason = "converged";
    result.sum_mode_effective = "direct";
    return result;
  }

  const double kappa = debye_kappa_inv_a(options.salt_m,
                                         options.temperature_k);
  const double contrast =
      (options.eps_out - options.eps_in) /
      std::max(1.0e-12, options.eps_out + options.eps_in);
  const double operator_contrast = contrast * formulation.operator_contrast_scale;
  result.kappa_inv_a = kappa;
  result.dielectric_contrast = contrast;
  result.operator_contrast = operator_contrast;
  std::string mode = options.sum_mode;
  if (mode == "auto") {
    mode = formulation.formulation == "jcp2013"
               ? "direct"
               : (result.surface_geometry == "triangle"
               ? "direct"
               : (panels.size() > 1024 ? "treecode" : "direct"));
  }
  result.sum_mode_effective = mode;
  result.matvec_mode = mode;

  if (formulation.formulation == "jcp2013") {
    if (mode == "treecode") {
      throw std::runtime_error(
          "PBE: -pbe=tabi -pbe-tabi-formulation=jcp2013 does not support "
          "-pbe-tabi-sum=treecode yet; use auto/direct, or pin "
          "-pbe-tabi-formulation=titan-v1 for the legacy scalar treecode.");
    }
    result.matvec_mode = "direct";
    result.sum_mode_effective = "direct";
    result.taylor_order = 0;
    const double eps_ratio =
        std::max(1.0e-12, options.eps_out) /
        std::max(1.0e-12, options.eps_in);
    const double lhs_phi = 0.5 * (1.0 + eps_ratio);
    const double lhs_phin = 0.5 * (1.0 + 1.0 / eps_ratio);
    const bool coupled = kappa > 1.0e-14;
    result.operator_contrast = eps_ratio;
    result.unknown_count = coupled ? 2 * panels.size() : panels.size();
    result.bie_system = coupled ? "linear_pb_2x2" : "poisson";
    result.lhs_diag_phi = lhs_phi;
    result.lhs_diag_phin = coupled ? lhs_phin : 0.0;
    result.energy_equation = coupled ? "jcp2013_eq14" : "jcp2013_eq22";
    result.kernel_convention =
        "G0_Gk_source_normal_target_normal_outward_centroid_collocation";

    const auto rhs_begin = std::chrono::steady_clock::now();
    double s1_norm = 0.0;
    double s2_norm = 0.0;
    const std::vector<double> rhs = build_jcp_rhs(
        atoms, panels, options.eps_in, coupled, &s1_norm, &s2_norm);
    result.rhs_ms = elapsed_ms(rhs_begin, std::chrono::steady_clock::now());
    result.source_norm_s1 = s1_norm;
    result.source_norm_s2 = s2_norm;

    const int size_factor = coupled ? 8 : 4;
    const int max_iter =
        std::min(std::max(1, options.max_iter),
                 std::max(64, static_cast<int>(panels.size()) * size_factor));
    const double tol = std::max(1.0e-12, options.residual_tol);
    std::vector<double> diagonal(result.unknown_count, lhs_phi);
    if (coupled) {
      for (std::size_t i = panels.size(); i < diagonal.size(); ++i)
        diagonal[i] = lhs_phin;
    }

    double matvec_ms_total = 0.0;
    std::size_t matvec_calls = 0;
    auto matvec = [&](const std::vector<double> &x, std::vector<double> &out) {
      const auto begin = std::chrono::steady_clock::now();
      if (coupled) {
        apply_jcp_linear_pb_operator(panels, x, eps_ratio, kappa, lhs_phi,
                                     lhs_phin, executor, out);
      } else {
        apply_jcp_poisson_operator(panels, x, eps_ratio, lhs_phi, executor,
                                   out);
      }
      matvec_ms_total += elapsed_ms(begin, std::chrono::steady_clock::now());
      ++matvec_calls;
    };

    const auto gmres_begin = std::chrono::steady_clock::now();
    LinearSolveResult linear = gmres_solve(matvec, rhs, diagonal, max_iter, tol);
    result.gmres_ms = elapsed_ms(gmres_begin, std::chrono::steady_clock::now());
    result.iterations = linear.iterations;
    result.residual = linear.residual;
    result.converged = linear.converged;
    result.convergence_reason = linear.convergence_reason;
    result.residual_history = residual_history_string(linear.residual_history);
    result.matvec_ms_total = matvec_ms_total;
    result.matvec_calls = matvec_calls;
    result.matvec_ms_avg =
        matvec_calls > 0 ? matvec_ms_total / static_cast<double>(matvec_calls)
                         : 0.0;
    result.orthogonalization_ms = linear.orthogonalization_ms;
    result.preconditioner_ms = linear.preconditioner_ms;
    if (coupled) {
      result.phi_norm = rms_norm(std::vector<double>(
          linear.x.begin(), linear.x.begin() + panels.size()));
      result.phin_norm = rms_norm(std::vector<double>(
          linear.x.begin() + panels.size(), linear.x.end()));
    } else {
      result.phi_norm = rms_norm(linear.x);
      result.phin_norm = 0.0;
    }

    const double native_energy_scale = 4.0 * kPi;
    const double boundary_energy =
        native_energy_scale *
        (coupled ? jcp_linear_pb_energy_kcal(atoms, panels, linear.x,
                                             eps_ratio, kappa)
                 : jcp_poisson_energy_kcal(atoms, panels, linear.x,
                                           eps_ratio));
    const double born = born_reference_kcal(atoms, options.eps_in,
                                            options.eps_out, kappa);
    double energy = boundary_energy;
    std::string source = "tabi_boundary_integral";
    double selected_scale = 1.0;
    if (formulation.energy_projection == "born-check") {
      energy = born;
      source = "tabi_born_check";
      selected_scale = 0.0;
    } else if (formulation.energy_projection == "charge-potential-debug") {
      source = "tabi_charge_potential_debug";
    }
    result.boundary_energy_kcal_mol = boundary_energy;
    result.reference_energy_kcal_mol = born;
    result.selected_energy_scale = selected_scale;
    result.fallback_used = false;
    result.reaction_energy_source = source;
    if (!linear.converged || !std::isfinite(energy)) {
      result.status = "failed";
      if (!std::isfinite(energy))
        result.convergence_reason = "nan_guard";
    } else {
      result.status = "ok";
    }
    result.energy_kcal_mol = energy;
    result.energy_kj_mol = energy * kKcalToKJ;
    result.energy_kt =
        result.energy_kj_mol /
        (kKBKJPerMolK * std::max(1.0, options.temperature_k));
    result.peak_memory_bytes = result.memory_bytes;
    return result;
  }

  TreeWorkspace tree;
  TreecodeStats matvec_stats;
  matvec_stats.matvec_mode = mode;
  matvec_stats.taylor_order = mode == "treecode" ? result.taylor_order : 0;
  if (mode == "treecode") {
    tree = build_tree_workspace(panels, executor);
    tree.taylor_order = result.taylor_order;
    matvec_stats.tree_build_ms = tree.build_ms;
    matvec_stats.tree_nodes = tree.nodes.size();
    matvec_stats.tree_leaf_count = tree.leaf_count;
    result.tree_build_ms = tree.build_ms;
    result.tree_nodes = tree.nodes.size();
    result.tree_leaf_count = tree.leaf_count;
    result.memory_bytes += tree.nodes.size() * sizeof(TreeNode) +
                           tree.order.size() * sizeof(std::size_t) +
                           tree.morton.size() * sizeof(std::uint64_t);
  }

  const auto rhs_begin = std::chrono::steady_clock::now();
  const std::vector<double> rhs =
      build_rhs(atoms, panels, contrast, kappa, formulation.rhs_scale);
  result.rhs_ms = elapsed_ms(rhs_begin, std::chrono::steady_clock::now());
  const int max_iter =
      std::min(std::max(1, options.max_iter),
               std::max(32, static_cast<int>(panels.size()) * 4));
  const double tol = std::max(1.0e-12, options.residual_tol);
  double matvec_ms_total = 0.0;
  std::size_t matvec_calls = 0;
  auto matvec = [&](const std::vector<double> &x, std::vector<double> &out) {
    const auto begin = std::chrono::steady_clock::now();
    if (mode == "treecode")
      apply_treecode_operator(panels, x, operator_contrast, kappa, tree,
                              executor, out, &matvec_stats);
    else
      apply_direct_operator(panels, x, operator_contrast, kappa, executor, out,
                            &matvec_stats);
    matvec_ms_total += elapsed_ms(begin, std::chrono::steady_clock::now());
    ++matvec_calls;
  };
  const auto gmres_begin = std::chrono::steady_clock::now();
  std::vector<double> legacy_diagonal(panels.size(), 1.0);
  for (std::size_t i = 0; i < panels.size(); ++i)
    legacy_diagonal[i] = panels[i].diagonal;
  LinearSolveResult linear =
      gmres_solve(matvec, rhs, legacy_diagonal, max_iter, tol);
  result.gmres_ms = elapsed_ms(gmres_begin, std::chrono::steady_clock::now());
  result.iterations = linear.iterations;
  result.residual = linear.residual;
  result.converged = linear.converged;
  result.convergence_reason = linear.convergence_reason;
  result.residual_history = residual_history_string(linear.residual_history);
  result.matvec_ms_total = matvec_ms_total;
  result.matvec_calls = matvec_calls;
  result.matvec_ms_avg =
      matvec_calls > 0 ? matvec_ms_total / static_cast<double>(matvec_calls)
                       : 0.0;
  result.orthogonalization_ms = linear.orthogonalization_ms;
  result.preconditioner_ms = linear.preconditioner_ms;
  result.tree_build_ms = matvec_stats.tree_build_ms;
  result.tree_nodes = matvec_stats.tree_nodes;
  result.tree_leaf_count = matvec_stats.tree_leaf_count;
  result.near_interactions = matvec_stats.near_interactions;
  result.far_interactions = matvec_stats.far_interactions;
  result.matvec_mode = matvec_stats.matvec_mode;
  result.taylor_order =
      mode == "treecode" ? matvec_stats.taylor_order
                         : std::clamp(options.taylor_order, 3, 4);
  result.interaction_build_ms = matvec_stats.interaction_build_ms;
  result.tree_eval_ms = matvec_stats.tree_eval_ms;
  result.tree_moment_ms = matvec_stats.tree_moment_ms;
  result.tree_traversal_ms = matvec_stats.tree_traversal_ms;
  result.near_eval_ms = matvec_stats.near_eval_ms;
  result.far_eval_ms = matvec_stats.far_eval_ms;

  const double boundary_energy =
      boundary_energy_kcal(atoms, panels, linear.x, kappa) *
      formulation.boundary_energy_scale;
  const double born = born_reference_kcal(atoms, options.eps_in,
                                          options.eps_out, kappa);
  double energy = boundary_energy;
  double selected_scale = formulation.boundary_energy_scale;
  bool used_fallback = false;
  std::string source = "tabi_boundary_integral";
  if (formulation.energy_projection == "born-check") {
    energy = born;
    source = "tabi_born_check";
    selected_scale = 0.0;
  } else if (formulation.formulation == "calibrated-debug" &&
             std::isfinite(boundary_energy) &&
             std::abs(boundary_energy) > 1.0e-12 &&
             std::isfinite(born) && std::abs(born) > 1.0e-12) {
    selected_scale = born / boundary_energy;
    energy = boundary_energy * selected_scale;
    source = "tabi_calibrated_debug";
  } else if (formulation.energy_projection == "charge-potential-debug") {
    source = "tabi_charge_potential_debug";
  } else if (!std::isfinite(energy) || std::abs(energy) <= 1.0e-12 ||
             (std::abs(born) > 1.0e-12 && energy * born < 0.0)) {
    energy = born;
    used_fallback = true;
    source = "tabi_born_fallback";
  }
  result.boundary_energy_kcal_mol = boundary_energy;
  result.reference_energy_kcal_mol = born;
  result.selected_energy_scale = selected_scale;
  result.fallback_used = used_fallback;
  result.reaction_energy_source = source;
  if (!linear.converged && result.convergence_reason.empty())
    result.convergence_reason = "max_iter";
  if (!linear.converged || !std::isfinite(energy)) {
    result.status = "failed";
    if (!std::isfinite(energy))
      result.convergence_reason = "nan_guard";
  } else if (used_fallback) {
    result.status = "fallback";
  } else {
    result.status = "ok";
  }
  result.energy_kcal_mol = energy;
  result.energy_kj_mol = energy * kKcalToKJ;
  result.energy_kt =
      result.energy_kj_mol /
      (kKBKJPerMolK * std::max(1.0, options.temperature_k));
  result.peak_memory_bytes = result.memory_bytes;
  return result;
}

} // namespace TitanPBE::Tabi
