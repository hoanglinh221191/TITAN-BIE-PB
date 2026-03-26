#include "MOL.h"
#include "differential_evolution.h"
#include "functions.h"
#include "omp.h"
#include "titan_ui.h"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <random>
#include <string_view>


using std::string_view;
#include <limits>
#include <thread>

const double PI = 3.14159265;

void cartesian_to_spherical(double x, double y, double z, double &theta,
                            double &phi) {
  double r = std::sqrt(x * x + y * y + z * z);
  if (r < 1e-6) {
    theta = 0.0;
    phi = 0.0;
  } else {
    theta = std::acos(z / r);
    phi = std::atan2(y, x); // Kết quả từ -PI đến PI
    if (phi < 0)
      phi += 2 * 3.14159265358979323846; // Chuẩn hóa về [0, 2PI]
  }
}

double normalize_angle_2pi(double angle) {
  angle = fmod(angle, 2.0 * PI);
  if (angle < 0)
    angle += 2.0 * PI;
  return angle;
}

// Hàm kẹp góc Theta về [0, PI]
double clamp_theta(double theta) {
  if (theta < 0)
    return -theta; // Phản xạ
  if (theta > PI)
    return 2.0 * PI - theta; // Phản xạ bờ bên kia
  return theta;
}

bool restrain_dis_ligand_pathway(groups_cal &groups, Multi_MOL &rec,
                                 MOL &lg_candidate, arguments &Arguments,
                                 double &distance) {
  if (Arguments.de.flag_pathway == false)
    return true;
  double com_lg_candidate[3]{};
  lg_candidate.get_centroid(com_lg_candidate);

  bool flag_dis_min{false};

  for (size_t i{}; i < groups.groups_atom_number(groups.group_pathway_index());
       ++i) {
    if (rec.info[groups.atom_index(groups.group_pathway_index(), i) - 1]
            .get_atom_type()
            .substr(0, 1) != "H") {
      distance = rec.m_mol[0].squared_dis_between_atom_and_point(
          groups.atom_index(groups.group_pathway_index(), i) - 1,
          com_lg_candidate[0], com_lg_candidate[1], com_lg_candidate[2]);
      if (distance < Arguments.de.dis_along_path * Arguments.de.dis_along_path) {
        flag_dis_min = true;
        break;
      }
    }
  }
  return flag_dis_min;
}

double cal_min_dis_ligand_pathway(groups_cal &groups, Multi_MOL &rec,
                                  MOL &lg_candidate, arguments &Arguments) {
  if (Arguments.de.flag_pathway == false)
    return true;
  double com_lg_candidate[3]{};
  lg_candidate.get_centroid(com_lg_candidate);

  double mindis{10e10};

  for (size_t i{}; i < groups.groups_atom_number(groups.group_pathway_index());
       ++i) {
    if (rec.info[groups.atom_index(groups.group_pathway_index(), i) - 1]
            .get_atom_type()
            .substr(0, 1) != "H") {
      double dis{rec.m_mol[0].squared_dis_between_atom_and_point(
          groups.atom_index(groups.group_pathway_index(), i) - 1,
          com_lg_candidate[0], com_lg_candidate[1], com_lg_candidate[2])};
      if (dis < mindis) {
        mindis = dis;
      }
    }
  }
  return sqrt(mindis);
}

double
cal_min_dis_ligand_pathway_optimized(groups_cal &groups, Multi_MOL &rec,
                                     const std::span<double> com_lg_candidate,
                                     arguments &Arguments) {
  if (Arguments.de.flag_pathway == false)
    return 1.0;

  double mindis{10e10};

  size_t grp_size = groups.groups_atom_number(groups.group_pathway_index());
  for (size_t i{}; i < grp_size; ++i) {
    size_t atom_idx = groups.atom_index(groups.group_pathway_index(), i);
    if (rec.info[atom_idx - 1].get_atom_type().substr(0, 1) != "H") {
      double dis{rec.m_mol[0].squared_dis_between_atom_and_point(
          atom_idx - 1, com_lg_candidate[0], com_lg_candidate[1],
          com_lg_candidate[2])};
      if (dis < mindis) {
        mindis = dis;
      }
    }
  }
  return sqrt(mindis);
}

double cal_score_function(cache_energy &energy_trial, cache_energy &energy_ref,
                          int flag_function, double dist, double limit) {
  double base_score{};
  if (flag_function == 1) {
    base_score = energy_trial.elec_ener + energy_trial.vdw_ener -
                 energy_ref.elec_ener - energy_ref.vdw_ener +
                 2. * fabs(energy_trial.elec_ener + energy_trial.vdw_ener);
  } else {
    base_score =
        -1. * exp(1.6863 * (-energy_trial.elec_ener - energy_trial.vdw_ener +
                            energy_ref.elec_ener + energy_ref.vdw_ener));
  }
  double K_spring{1000.0};
  double penalty{dist > limit ? K_spring * (dist - limit) * (dist - limit)
                              : 0.};
  return base_score + penalty;
}

bool check_radius(de_arg &DE_arguments, force &tran_vector) {
  return (tran_vector.modul() < DE_arguments.scan_radius &&
          tran_vector.modul() > DE_arguments.min_scan_radius);
}

bool parse_angle(de_arg &DE_arguments, force &tran_vector, double &crit_angle) {

  if (!DE_arguments.flag_old_vec) {
    crit_angle = 0.0;
    return true;
  }

  // 1. Lấy hướng của vector mới (Unit Vector từ Spherical Coords)
  // tran_vector: Fx = theta, Fy = phi (Fz = rho không ảnh hưởng đến góc)
  double sin_theta = sin(tran_vector.Fx);
  double nx = sin_theta * cos(tran_vector.Fy);
  double ny = sin_theta * sin(tran_vector.Fy);
  double nz = cos(tran_vector.Fx);

  // 2. Lấy hướng của vector cũ (Cần chuẩn hóa nếu chưa)
  double old_mod = DE_arguments.old_vector.modul();
  if (old_mod < 1e-9)
    return true; // Tránh chia cho 0

  double ox = DE_arguments.old_vector.Fx / old_mod;
  double oy = DE_arguments.old_vector.Fy / old_mod;
  double oz = DE_arguments.old_vector.Fz / old_mod;

  // 3. Tính Cosine góc kẹp bằng Tích Vô Hướng (Dot Product)
  // cos(theta) = a_unit . b_unit
  double cos_val = nx * ox + ny * oy + nz * oz;

  // Kẹp giá trị trong [-1, 1] để tránh lỗi sai số máy tính
  if (cos_val > 1.0)
    cos_val = 1.0;
  if (cos_val < -1.0)
    cos_val = -1.0;

  // 4. Tính góc thực tế (để debug hoặc hiển thị nếu cần)
  // 180/PI = 57.2957795
  crit_angle = acos(cos_val) * 57.2957795;

  // 5. Kiểm tra điều kiện
  // Nếu góc < giới hạn thì trả về true
  return (crit_angle < DE_arguments.restrain_angle);
}

int parse_new_mol(MOL &mol_temp, double dX, double dY, double dZ,
                  arguments &Arguments, Multi_MOL &rec, Multi_MOL &lg_seed,
                  cache_sigma_epsilon &CachedData, cache_energy &e_cache) {
  mol_temp.translate(dX, dY, dZ);
  if (Arguments.force_f.find("CHARMM") != string::npos ||
      Arguments.force_f.find("AMBER") != string::npos ||
      Arguments.force_f.find("OPLSAA") != string::npos) {
    cal_energy_amber(rec.m_mol[0], mol_temp, rec.info, lg_seed.info, CachedData,
                     0, mol_temp.size(), e_cache, 0, Arguments, rec.parameters,
                     lg_seed.parameters);
  } else if (Arguments.force_f.find("MARTINI") != string::npos) {
    cal_energy_martini(rec.m_mol[0], mol_temp, rec.info, lg_seed.info,
                       CachedData, 0, mol_temp.size(), e_cache, 0, Arguments,
                       rec.parameters, lg_seed.parameters);
  } else {
    cal_energy_gromos(rec.m_mol[0], mol_temp, rec.info, lg_seed.info,
                      CachedData, 0, mol_temp.size(), e_cache, 0, Arguments,
                      rec.parameters, lg_seed.parameters);
  }
  return 0;
}

int Multi_MOL::gen_pop(MOL &mol_seed, de_arg &DE_arguments,
                       const vector<mol_info> &info, string_view mol_type,
                       vector<force> &tran_data) {
  const double PI_VAL = 3.14159265358979323846;

  std::random_device rd;
  std::mt19937 mt(rd());

  // CẢI TIẾN 1: Phân bố đều trên mặt cầu (Uniform Sphere Distribution)
  // Thay vì theta đều [0, PI], ta random cos(theta) trong [-1, 1] để tránh tụ
  // điểm ở cực
  std::uniform_real_distribution<double> dist_phi(0, 2.0 * PI_VAL);
  std::uniform_real_distribution<double> dist_cos_theta(-1.0, 1.0);
  std::uniform_real_distribution<double> dist_rho(DE_arguments.min_scan_radius,
                                                  DE_arguments.scan_radius);

  // Phân bố nhiễu Gaussian cho trường hợp fallback (cứu hộ khi bế tắc)
  std::normal_distribution<double> dist_noise(
      0.0, 0.1); // Sigma = 0.1 radian (~5.7 độ)

  // this->m_mol.clear(); // Đảm bảo sạch sẽ trước khi reserve
  this->m_mol.reserve(DE_arguments.pop_num);
  tran_data.reserve(DE_arguments.pop_num); // Reserve vector này luôn

  // Chuẩn bị vector tham chiếu (nếu có Restraint)
  double theta_ref = 0.0, phi_ref = 0.0;
  if (DE_arguments.flag_old_vec) {
    // Chuyển vector cũ sang tọa độ cầu để dùng làm "mỏ neo"
    cartesian_to_spherical(DE_arguments.old_vector.Fx,
                           DE_arguments.old_vector.Fy,
                           DE_arguments.old_vector.Fz, theta_ref, phi_ref);
  }

  for (size_t i{}; i < DE_arguments.pop_num; ++i) {
    double phi{}, theta{}, rho{};
    double dX{}, dY{}, dZ{};
    bool valid_sample = false;
    int attempts = 0;
    const int MAX_ATTEMPTS = 500; // Giới hạn số lần thử sinh ngẫu nhiên

    // Vòng lặp thử mẫu (Rejection Sampling Loop)
    while (attempts < MAX_ATTEMPTS) {
      phi = dist_phi(mt);
      double cos_theta = dist_cos_theta(mt);
      theta = std::acos(cos_theta);
      rho = dist_rho(mt);

      // Tính toán Cartesian
      dX = rho * std::sin(theta) * std::cos(phi);
      dY = rho * std::sin(theta) * std::sin(phi);
      dZ = rho * std::cos(theta);

      // Kiểm tra Restraint
      if (DE_arguments.flag_old_vec) {
        force v_new(dX, dY, dZ);
        // Tính góc giữa vector mới và vector cũ
        double dot = v_new.Fx * DE_arguments.old_vector.Fx +
                     v_new.Fy * DE_arguments.old_vector.Fy +
                     v_new.Fz * DE_arguments.old_vector.Fz;

        // Vector chênh lệch để tính cross product
        double modul1 = v_new.modul(); // Chính là rho
        double modul2 = DE_arguments.old_vector.modul();

        // Tránh chia cho 0
        if (modul1 < 1e-9 || modul2 < 1e-9) {
          attempts++;
          continue;
        }

        force v1(dX, dY, dZ);
        v1 *= modul2;
        force v2 = DE_arguments.old_vector;
        v2 *= modul1;

        double angle_deg = std::fabs(
            360.0 * std::atan2(force(v1 - v2).modul(), force(v1 + v2).modul()) /
            PI_VAL);

        if (angle_deg < DE_arguments.restrain_angle) {
          valid_sample = true;
          break; // Mẫu tốt, thoát vòng lặp thử
        }
      } else {
        valid_sample = true;
        break; // Không có restraint thì luôn đúng
      }
      attempts++;
    }

    // CƠ CHẾ FALLBACK: Nếu thử 500 lần vẫn không tìm ra góc thỏa mãn (do góc
    // quá hẹp) Ta cưỡng ép lấy vector cũ cộng thêm chút nhiễu (noise)
    if (!valid_sample && DE_arguments.flag_old_vec) {
      theta = theta_ref + dist_noise(mt);
      phi = phi_ref + dist_noise(mt);
      rho = dist_rho(mt);

      // Kẹp lại giá trị hợp lệ
      if (theta < 0)
        theta = -theta;
      if (theta > PI_VAL)
        theta = 2 * PI_VAL - theta;
      // Phi tự động wrap qua sin/cos nên không cần kẹp quá chặt

      // Tính lại Cartesian từ fallback
      dX = rho * std::sin(theta) * std::cos(phi);
      dY = rho * std::sin(theta) * std::sin(phi);
      dZ = rho * std::cos(theta);
    }

    // Lưu dữ liệu vector (Theta, Phi, Rho) dùng cho Mutation sau này
    tran_data.emplace_back(theta, phi, rho);

    // Tạo snapshot mới trong Multi_MOL
    MOL new_mol;
    new_mol.Box = mol_seed.Box;
    new_mol.set_time(mol_seed.get_time());
    new_mol.force_f = mol_seed.get_force_f();
    for (size_t j{}; j < mol_seed.size(); ++j) {
      new_mol.pos.x.emplace_back(mol_seed.pos.x[j] + dX);
      new_mol.pos.y.emplace_back(mol_seed.pos.y[j] + dY);
      new_mol.pos.z.emplace_back(mol_seed.pos.z[j] + dZ);
    }
    // Safe emplace back with native vector std::move instead of helper emplace
    this->m_mol.emplace_back(std::move(new_mol));
    this->m_mol.back().mol_type = string(mol_type);
  }

  // Safe info copy to array to prevent Segfault if empty
  this->info = info;

  return 0;
}

int cal_energy_part_DE(Multi_MOL &rec, Multi_MOL &lg_seed,
                       Multi_MOL &lg_candidates,
                       cache_sigma_epsilon &CachedData, de_arg &DE_arguments,
                       vector<cache_energy> &Energies, arguments &Arguments,
                       size_t begin, size_t end) {
  if (Arguments.force_f.find("CHARMM") != string::npos ||
      Arguments.force_f.find("AMBER") != string::npos ||
      Arguments.force_f.find("OPLSAA") != string::npos) {
    for (size_t i{begin}; i < end; ++i)
      cal_energy_amber(rec.m_mol[0], lg_candidates.m_mol[i], rec.info,
                       lg_seed.info, CachedData, 0, lg_seed.m_mol[0].size(),
                       Energies[i], i, Arguments, rec.parameters,
                       lg_seed.parameters);
  } else if (Arguments.force_f.find("MARTINI") != string::npos) {
    for (size_t i{begin}; i < end; ++i)
      cal_energy_martini(rec.m_mol[0], lg_candidates.m_mol[i], rec.info,
                         lg_seed.info, CachedData, 0, lg_seed.m_mol[0].size(),
                         Energies[i], i, Arguments, rec.parameters,
                         lg_seed.parameters);
  } else {
    for (size_t i{begin}; i < end; ++i)
      cal_energy_gromos(rec.m_mol[0], lg_candidates.m_mol[i], rec.info,
                        lg_seed.info, CachedData, 0, lg_seed.m_mol[0].size(),
                        Energies[i], i, Arguments, rec.parameters,
                        lg_seed.parameters);
  }
  return 0;
}

int impl_evo_1_ite(Multi_MOL &rec, Multi_MOL &lg_seed, Multi_MOL &lg_candidates,
                   cache_sigma_epsilon &CachedData, de_arg &DE_arguments,
                   vector<cache_energy> &Energies, arguments &Arguments,
                   size_t begin, size_t end, const vector<force> &tran_data,
                   vector<force> &out_tran, size_t thread_rank,
                   cache_energy &energy_ref, groups_cal &Groups_from_index_file,
                   int flag_score_function) {
  static thread_local std::random_device rd;
  static thread_local std::mt19937 mt(rd());
  std::uniform_int_distribution<size_t> dist_choose(0,
                                                    DE_arguments.pop_num - 1);
  std::uniform_real_distribution<double> dist_gencrit(0., 1.);

  for (size_t i{begin}; i < end; ++i) {
    size_t abc_pop[3]{};
    size_t temp_1{(dist_choose(mt) + thread_rank) % DE_arguments.pop_num},
        temp_2{dist_choose(mt)}, temp_3{dist_choose(mt)};

    while (temp_1 >= DE_arguments.pop_num || temp_1 == temp_2 ||
           temp_1 == temp_3 || temp_2 == temp_3 || temp_1 == i || temp_2 == i ||
           temp_3 == i) {
      temp_1 = dist_choose(mt);
      temp_2 = dist_choose(mt);
      temp_3 = dist_choose(mt);
    }
    abc_pop[0] = temp_1;
    abc_pop[1] = temp_2;
    abc_pop[2] = temp_3;
    size_t random_index{dist_choose(mt)};
    double crit{dist_gencrit(mt)};
    force tran_vector;
    MOL mol_temp{lg_candidates.m_mol[i]};
    cache_energy e_cache;
    double crit_angle{180.};
    bool is_valid_angle{}, is_valid_radius{};
    if (DE_arguments.algorithm == 0) {
      if (crit < DE_arguments.cross_prob || i == random_index) {
        tran_vector.Fx = tran_data[abc_pop[0]].Fx +
                         DE_arguments.diff_weight * (tran_data[abc_pop[1]].Fx -
                                                     tran_data[abc_pop[2]].Fx);
        tran_vector.Fy = tran_data[abc_pop[0]].Fy +
                         DE_arguments.diff_weight * (tran_data[abc_pop[1]].Fy -
                                                     tran_data[abc_pop[2]].Fy);
        tran_vector.Fz = tran_data[abc_pop[0]].Fz +
                         DE_arguments.diff_weight * (tran_data[abc_pop[1]].Fz -
                                                     tran_data[abc_pop[2]].Fz);
      }
      tran_vector.Fx = clamp_theta(tran_vector.Fx); // Theta
      tran_vector.Fy = normalize_angle_2pi(tran_vector.Fy);
      tran_vector.Fz = fabs(tran_vector.Fz);
      is_valid_angle = parse_angle(DE_arguments, tran_vector, crit_angle);
      is_valid_radius = check_radius(DE_arguments, tran_vector);
      if (!is_valid_angle || !is_valid_radius)
        continue;

      double dX{tran_vector.Fz * sin(tran_vector.Fx) * cos(tran_vector.Fy) -
                tran_data[i].Fz * sin(tran_data[i].Fx) * cos(tran_data[i].Fy)},
          dY{tran_vector.Fz * sin(tran_vector.Fx) * sin(tran_vector.Fy) -
             tran_data[i].Fz * sin(tran_data[i].Fx) * sin(tran_data[i].Fy)},
          dZ{tran_vector.Fz * cos(tran_vector.Fx) -
             tran_data[i].Fz * cos(tran_data[i].Fx)};

      parse_new_mol(mol_temp, dX, dY, dZ, Arguments, rec, lg_seed, CachedData,
                    e_cache);
    } else {
      if (crit < DE_arguments.cross_prob) {
        tran_vector.Fx = tran_data[abc_pop[0]].Fx +
                         DE_arguments.F * (tran_data[abc_pop[1]].Fx -
                                           tran_data[abc_pop[2]].Fx);
        tran_vector.Fy = tran_data[abc_pop[0]].Fy +
                         DE_arguments.F * (tran_data[abc_pop[1]].Fy -
                                           tran_data[abc_pop[2]].Fy);
        tran_vector.Fz = tran_data[abc_pop[0]].Fz +
                         DE_arguments.F * (tran_data[abc_pop[1]].Fz -
                                           tran_data[abc_pop[2]].Fz);
        tran_vector.Fx = clamp_theta(tran_vector.Fx); // Theta
        tran_vector.Fy = normalize_angle_2pi(tran_vector.Fy);
        tran_vector.Fz = fabs(tran_vector.Fz);

        double dX{tran_vector.Fz * sin(tran_vector.Fx) * cos(tran_vector.Fy) -
                  tran_data[i].Fz * sin(tran_data[i].Fx) *
                      cos(tran_data[i].Fy)},
            dY{tran_vector.Fz * sin(tran_vector.Fx) * sin(tran_vector.Fy) -
               tran_data[i].Fz * sin(tran_data[i].Fx) * sin(tran_data[i].Fy)},
            dZ{tran_vector.Fz * cos(tran_vector.Fx) -
               tran_data[i].Fz * cos(tran_data[i].Fx)};

        is_valid_angle = parse_angle(DE_arguments, tran_vector, crit_angle);
        is_valid_radius = check_radius(DE_arguments, tran_vector);
        if (!is_valid_angle || !is_valid_radius)
          continue;

        parse_new_mol(mol_temp, dX, dY, dZ, Arguments, rec, lg_seed, CachedData,
                      e_cache);
      } else {
        tran_vector.Fx =
            tran_data[abc_pop[0]].Fx +
            0.5 * (DE_arguments.F + 1.) *
                (tran_data[abc_pop[1]].Fx + tran_data[abc_pop[2]].Fx -
                 2. * tran_data[abc_pop[0]].Fx);
        tran_vector.Fy =
            tran_data[abc_pop[0]].Fy +
            0.5 * (DE_arguments.F + 1.) *
                (tran_data[abc_pop[1]].Fy + tran_data[abc_pop[2]].Fy -
                 2. * tran_data[abc_pop[0]].Fy);
        tran_vector.Fz =
            tran_data[abc_pop[0]].Fz +
            0.5 * (DE_arguments.F + 1.) *
                (tran_data[abc_pop[1]].Fz + tran_data[abc_pop[2]].Fz -
                 2. * tran_data[abc_pop[0]].Fz);
        tran_vector.Fx = clamp_theta(tran_vector.Fx); // Theta
        tran_vector.Fy = normalize_angle_2pi(tran_vector.Fy);
        tran_vector.Fz = fabs(tran_vector.Fz);

        double dX{tran_vector.Fz * sin(tran_vector.Fx) * cos(tran_vector.Fy) -
                  tran_data[i].Fz * sin(tran_data[i].Fx) *
                      cos(tran_data[i].Fy)},
            dY{tran_vector.Fz * sin(tran_vector.Fx) * sin(tran_vector.Fy) -
               tran_data[i].Fz * sin(tran_data[i].Fx) * sin(tran_data[i].Fy)},
            dZ{tran_vector.Fz * cos(tran_vector.Fx) -
               tran_data[i].Fz * cos(tran_data[i].Fx)};

        is_valid_angle = parse_angle(DE_arguments, tran_vector, crit_angle);
        is_valid_radius = check_radius(DE_arguments, tran_vector);
        if (!is_valid_angle || !is_valid_radius)
          continue;

        parse_new_mol(mol_temp, dX, dY, dZ, Arguments, rec, lg_seed, CachedData,
                      e_cache);
      }
    }

    double dis_candidate{cal_min_dis_ligand_pathway(Groups_from_index_file, rec,
                                                    mol_temp, Arguments)},
        dis_old_mem{cal_min_dis_ligand_pathway(
            Groups_from_index_file, rec, lg_candidates.m_mol[i], Arguments)};

    if (is_valid_angle &&
        is_valid_radius
        //&&
        //(restrain_dis_ligand_pathway(Groups_from_index_file,rec,mol_temp,Arguments)||
        // dis_candidate<=dis_old_mem)
        && (cal_score_function(e_cache, energy_ref, flag_score_function,
                               dis_candidate, Arguments.de.dis_along_path) <
            cal_score_function(Energies[i], energy_ref, flag_score_function,
                               dis_old_mem, Arguments.de.dis_along_path))) {
      out_tran[i] = tran_vector;
      Energies[i] = e_cache;
      lg_candidates.m_mol[i] = mol_temp;
    }
  }

  return 0;
}

int evo_1_ite(Multi_MOL &rec, Multi_MOL &lg_seed, Multi_MOL &lg_candidates,
              cache_sigma_epsilon &CachedData, de_arg &DE_arguments,
              vector<cache_energy> &Energies, arguments &Arguments,
              vector<force> &tran_data, cache_energy &energy_ref,
              groups_cal &Groups_from_index_file) {
  vector<force> out_tran(tran_data);
  {
    std::vector<std::jthread> workers(Arguments.nt);
    for (size_t i{}; i < workers.size(); ++i)
      workers[i] = std::jthread([&, i] {
        impl_evo_1_ite(
            rec, lg_seed, lg_candidates, CachedData, DE_arguments, Energies,
            Arguments, i * lg_candidates.m_mol.size() / workers.size(),
            (i + 1) * lg_candidates.m_mol.size() / workers.size(), tran_data,
            out_tran, i + 1, energy_ref, Groups_from_index_file,
            Arguments.de.de_flag_score_function);
      });
  } // jthreads auto-join here
  tran_data = out_tran;
  return 0;
}

int find_pull_direct(Multi_MOL &rec, Multi_MOL &lg,
                     cache_sigma_epsilon &CachedData, de_arg &DE_arguments,
                     arguments &Arguments, info &output_info,
                     groups_cal &Groups_from_index_file) {
  if (rec.m_mol[0].mindist_mol(lg.m_mol[0], Arguments, DE_arguments) == true) {
    cout << "Automatic termination is on!\n";
    return 1;
  }

  cache_energy energy_ref;
  const size_t n_threads = Arguments.nt;

  // Per-thread accumulators to avoid race condition on energy_ref
  {
    std::vector<cache_energy> energy_parts(n_threads);
    const size_t lg_size = lg.m_mol[0].size();
    {
      std::vector<std::jthread> workers(n_threads);
      for (size_t t = 0; t < n_threads; ++t)
        workers[t] = std::jthread([&, t] {
          size_t beg = t * lg_size / n_threads;
          size_t end = (t + 1) * lg_size / n_threads;
          if (Arguments.force_f.find("CHARMM") != string::npos ||
              Arguments.force_f.find("AMBER") != string::npos ||
              Arguments.force_f.find("OPLSAA") != string::npos) {
            cal_energy_amber(rec.m_mol[0], lg.m_mol[0], rec.info, lg.info,
                             CachedData, beg, end, energy_parts[t], 0,
                             Arguments, rec.parameters, lg.parameters);
          } else if (Arguments.force_f.find("MARTINI") != string::npos) {
            cal_energy_martini(rec.m_mol[0], lg.m_mol[0], rec.info, lg.info,
                               CachedData, beg, end, energy_parts[t], 0,
                               Arguments, rec.parameters, lg.parameters);
          } else {
            cal_energy_gromos(rec.m_mol[0], lg.m_mol[0], rec.info, lg.info,
                              CachedData, beg, end, energy_parts[t], 0,
                              Arguments, rec.parameters, lg.parameters);
          }
        });
    } // jthreads auto-join here
    // Merge (serial, no race)
    for (auto &ep : energy_parts) {
      energy_ref.elec_ener += ep.elec_ener;
      energy_ref.vdw_ener += ep.vdw_ener;
    }
  }
  TitanUI::print_info(
      "  Energy ref: elec=" + std::to_string(energy_ref.elec_ener) +
      "  vdw=" + std::to_string(energy_ref.vdw_ener));

  Multi_MOL lg_candidates;
  lg_candidates.m_mol.reserve(DE_arguments.pop_num);

  vector<force> tran_data;

  cout << "Generating population......";
  lg_candidates.gen_pop(lg.m_mol[0], DE_arguments, lg.info, lg.mol_type,
                        tran_data);
  cout << "done!\n";

  vector<double> Distances(DE_arguments.pop_num);
  double seed_centroid[3];
  lg.m_mol[0].get_centroid(seed_centroid);

#pragma omp parallel for
  for (int i = 0; i < DE_arguments.pop_num; ++i) {
    double theta = tran_data[i].Fx;
    double phi = tran_data[i].Fy;
    double rho = tran_data[i].Fz;

    double dX = rho * std::sin(theta) * std::cos(phi);
    double dY = rho * std::sin(theta) * std::sin(phi);
    double dZ = rho * std::cos(theta);

    double new_centroid[3];
    new_centroid[0] = seed_centroid[0] + dX;
    new_centroid[1] = seed_centroid[1] + dY;
    new_centroid[2] = seed_centroid[2] + dZ;

    Distances[i] = cal_min_dis_ligand_pathway_optimized(
        Groups_from_index_file, rec, new_centroid, Arguments);
  }

  vector<cache_energy> Energies(DE_arguments.pop_num);
  {
    std::vector<std::jthread> workers(n_threads);
    for (size_t i{}; i < n_threads; ++i)
      workers[i] = std::jthread([&, i] {
        cal_energy_part_DE(rec, lg, lg_candidates, CachedData, DE_arguments,
                           Energies, Arguments,
                           i * lg_candidates.m_mol.size() / n_threads,
                           (i + 1) * lg_candidates.m_mol.size() / n_threads);
      });
  } // jthreads auto-join here

  // --- DE Evolution Loop with early stopping & adaptive F ---
  double prev_best = std::numeric_limits<double>::max();
  int stagnation_count = 0;
  const int stag_limit = DE_arguments.stagnation_limit;

  for (size_t iter = 0; iter < DE_arguments.ite_num; ++iter) {
    TitanUI::show_progress("Pulling Direction", iter, DE_arguments.ite_num,
                           "(Evolution)");
    evo_1_ite(rec, lg, lg_candidates, CachedData, DE_arguments, Energies,
              Arguments, tran_data, energy_ref, Groups_from_index_file);

    // Track best score this iteration
    double iter_best = std::numeric_limits<double>::max();
    for (size_t j = 0; j < Energies.size(); ++j) {
      double s = cal_score_function(Energies[j], energy_ref,
                                    Arguments.de.de_flag_score_function, 0.,
                                    Arguments.de.dis_along_path);
      if (s < iter_best)
        iter_best = s;
    }

    if (iter_best < prev_best - 1e-8) {
      prev_best = iter_best;
      stagnation_count = 0;
    } else {
      ++stagnation_count;
    }

    // Adaptive F: increase when stagnating (explore), decrease when improving
    if (DE_arguments.adaptive_F) {
      if (stagnation_count > 5) {
        DE_arguments.diff_weight =
            std::min(DE_arguments.F_max, DE_arguments.diff_weight * 1.05);
      } else if (stagnation_count == 0) {
        DE_arguments.diff_weight =
            std::max(DE_arguments.F_min, DE_arguments.diff_weight * 0.98);
      }
    }

    // Early stopping
    if (stagnation_count >= stag_limit) {
      TitanUI::print_info("  DE converged early at iteration " +
                          std::to_string(iter) +
                          " (stagnation=" + std::to_string(stag_limit) + ")");
      break;
    }
  }

  // --- Open log files ---
  fstream de_output(output_info.get_DE_output_filename(), ios::out),
      de_log(output_info.get_DE_log_filename(), ios::out | ios::app),
      de_pdb_out(output_info.get_DE_pdb_vector_filename(), ios::out | ios::app);
  fstream de_report("DE_report.dat", ios::out),
      best_conf("best_conf.pdb", ios::out);

  // --- Find best candidate ---
  double min_e{cal_score_function(Energies[0], energy_ref,
                                  Arguments.de.de_flag_score_function, 0.,
                                  Arguments.de.dis_along_path)};
  size_t pos{};
  for (size_t i{1}; i < Energies.size(); ++i) {
    double s = cal_score_function(Energies[i], energy_ref,
                                  Arguments.de.de_flag_score_function, 0.,
                                  Arguments.de.dis_along_path);
    if (s < min_e) {
      pos = i;
      min_e = s;
    }
  }

  // --- Compute best pull vector ---
  double dX{tran_data[pos].Fz * sin(tran_data[pos].Fx) *
            cos(tran_data[pos].Fy)},
      dY{tran_data[pos].Fz * sin(tran_data[pos].Fx) * sin(tran_data[pos].Fy)},
      dZ{tran_data[pos].Fz * cos(tran_data[pos].Fx)};
  double modul_old{DE_arguments.old_vector.modul()},
      modul_new{tran_data[pos].Fz};
  force new_vec{dX, dY, dZ}, old_vec{DE_arguments.old_vector};
  new_vec *= modul_old;
  old_vec *= modul_new;
  double best_angle =
      360. / PI *
      atan2(force(new_vec - old_vec).modul(), force(new_vec + old_vec).modul());

  // --- Ligand COM (mass-weighted) ---
  double lg_com_x{}, lg_com_y{}, lg_com_z{}, total_mass{};
  for (size_t i{}; i < lg.m_mol[0].size(); ++i) {
    double m = lg.parameters[i].get_mass();
    lg_com_x += lg.m_mol[0].pos.x[i] * m;
    lg_com_y += lg.m_mol[0].pos.y[i] * m;
    lg_com_z += lg.m_mol[0].pos.z[i] * m;
    total_mass += m;
  }
  lg_com_x /= total_mass;
  lg_com_y /= total_mass;
  lg_com_z /= total_mass;

  // Best candidate COM
  double best_com_x{lg_com_x + dX}, best_com_y{lg_com_y + dY},
      best_com_z{lg_com_z + dZ};

  // --- Find nearest receptor heavy atom to ligand COM (pull origin) ---
  size_t nearest_rec_idx{};
  double nearest_rec_dist{1e20};
  for (size_t i{}; i < rec.m_mol[0].size(); ++i) {
    if (rec.info[i].get_atom_type().substr(0, 1) == "H")
      continue;
    double dx = rec.m_mol[0].pos.x[i] - lg_com_x;
    double dy = rec.m_mol[0].pos.y[i] - lg_com_y;
    double dz = rec.m_mol[0].pos.z[i] - lg_com_z;
    double d2 = dx * dx + dy * dy + dz * dz;
    if (d2 < nearest_rec_dist) {
      nearest_rec_dist = d2;
      nearest_rec_idx = i;
    }
  }
  nearest_rec_dist = sqrt(nearest_rec_dist);

  // --- Find nearest receptor heavy atom to best candidate COM (pull target)
  // ---
  size_t target_rec_idx{};
  double target_rec_dist{1e20};
  for (size_t i{}; i < rec.m_mol[0].size(); ++i) {
    if (rec.info[i].get_atom_type().substr(0, 1) == "H")
      continue;
    double dx = rec.m_mol[0].pos.x[i] - best_com_x;
    double dy = rec.m_mol[0].pos.y[i] - best_com_y;
    double dz = rec.m_mol[0].pos.z[i] - best_com_z;
    double d2 = dx * dx + dy * dy + dz * dz;
    if (d2 < target_rec_dist) {
      target_rec_dist = d2;
      target_rec_idx = i;
    }
  }
  target_rec_dist = sqrt(target_rec_dist);

  // --- DE_report.dat: Comprehensive report ---
  de_report << fixed;
  de_report.precision(4);
  de_report << "==============================================================="
               "=================\n";
  de_report << "  TITAN Differential Evolution — Pull Direction Report\n";
  de_report << "==============================================================="
               "=================\n\n";

  de_report << "  Parameters\n";
  de_report << "  ----------\n";
  de_report << "  Population size    : " << DE_arguments.pop_num << "\n";
  de_report << "  Max iterations     : " << DE_arguments.ite_num << "\n";
  de_report << "  Algorithm          : " << DE_arguments.algorithm << "\n";
  de_report << "  Diff weight (F)    : " << DE_arguments.diff_weight << "\n";
  de_report << "  Cross prob (CR)    : " << DE_arguments.cross_prob << "\n";
  de_report << "  Scan radius        : " << DE_arguments.scan_radius << " A\n";
  de_report << "  Min scan radius    : " << DE_arguments.min_scan_radius
            << " A\n";
  de_report << "  Restrain angle     : " << DE_arguments.restrain_angle
            << " deg\n";
  de_report << "  Adaptive F         : "
            << (DE_arguments.adaptive_F ? "ON" : "OFF") << "\n";
  de_report << "  Stagnation limit   : " << DE_arguments.stagnation_limit
            << "\n";
  de_report << "  Score function     : " << Arguments.de.de_flag_score_function
            << "\n";
  de_report << "  Force field        : " << Arguments.force_f << "\n\n";

  de_report << "  Reference Energy\n";
  de_report << "  ----------------\n";
  de_report << "  Elec  : " << energy_ref.elec_ener << " kcal/mol\n";
  de_report << "  VdW   : " << energy_ref.vdw_ener << " kcal/mol\n";
  de_report << "  Total : " << energy_ref.elec_ener + energy_ref.vdw_ener
            << " kcal/mol\n\n";

  de_report << "  Ligand COM (original) : (" << lg_com_x << ", " << lg_com_y
            << ", " << lg_com_z << ")\n";
  de_report << "  Ligand COM (best)     : (" << best_com_x << ", " << best_com_y
            << ", " << best_com_z << ")\n\n";

  de_report << "  Pull Direction Anatomy\n";
  de_report << "  ----------------------\n";
  de_report << "  Origin receptor atom  : "
            << rec.info[nearest_rec_idx].get_atom() << " "
            << rec.info[nearest_rec_idx].get_original_res_name()
            << rec.info[nearest_rec_idx].get_respos() << " (chain "
            << rec.info[nearest_rec_idx].get_chainnum()
            << ", dist=" << nearest_rec_dist << " A)\n";
  de_report << "  Target receptor atom  : "
            << rec.info[target_rec_idx].get_atom() << " "
            << rec.info[target_rec_idx].get_original_res_name()
            << rec.info[target_rec_idx].get_respos() << " (chain "
            << rec.info[target_rec_idx].get_chainnum()
            << ", dist=" << target_rec_dist << " A)\n";
  de_report << "  Pull vector (dX,dY,dZ): (" << dX << ", " << dY << ", " << dZ
            << ")\n";
  de_report << "  Pull distance         : " << sqrt(dX * dX + dY * dY + dZ * dZ)
            << " A\n";
  de_report << "  Angle to prev vector  : " << best_angle << " deg\n\n";

  de_report << "  Best Candidate\n";
  de_report << "  --------------\n";
  de_report << "  Index  : " << pos << "\n";
  de_report << "  Score  : " << min_e << "\n";
  de_report << "  Elec   : " << Energies[pos].elec_ener << " kcal/mol\n";
  de_report << "  VdW    : " << Energies[pos].vdw_ener << " kcal/mol\n";
  de_report << "  dElec  : " << Energies[pos].elec_ener - energy_ref.elec_ener
            << " kcal/mol\n";
  de_report << "  dVdW   : " << Energies[pos].vdw_ener - energy_ref.vdw_ener
            << " kcal/mol\n\n";

  // --- Per-candidate table (top 20 + all if small pop) ---
  de_report << "  Population Energy Table\n";
  de_report << "  -----------------------\n";
  de_report << setw(8) << "Rank" << setw(8) << "Index" << setw(14) << "Score"
            << setw(14) << "Elec" << setw(14) << "VdW" << setw(14) << "dElec"
            << setw(14) << "dVdW" << setw(10) << "Rho(A)" << setw(12)
            << "Angle(deg)" << "\n";
  de_report << string(108, '-') << "\n";

  // Sort candidates by score for ranked output
  std::vector<std::pair<double, size_t>> scored(Energies.size());
  for (size_t i{}; i < Energies.size(); ++i) {
    scored[i] = {cal_score_function(Energies[i], energy_ref,
                                    Arguments.de.de_flag_score_function, 0.,
                                    Arguments.de.dis_along_path),
                 i};
  }
  std::sort(scored.begin(), scored.end());

  size_t n_print = std::min(scored.size(), size_t(50));
  for (size_t rank{}; rank < n_print; ++rank) {
    auto [score, idx] = scored[rank];
    double rho = tran_data[idx].Fz;
    // Angle to old vector
    double cdX{rho * sin(tran_data[idx].Fx) * cos(tran_data[idx].Fy)},
        cdY{rho * sin(tran_data[idx].Fx) * sin(tran_data[idx].Fy)},
        cdZ{rho * cos(tran_data[idx].Fx)};
    double cmod_old{DE_arguments.old_vector.modul()};
    force cvec_new{cdX, cdY, cdZ}, cvec_old{DE_arguments.old_vector};
    cvec_new *= cmod_old;
    cvec_old *= rho;
    double cangle = 360. / PI *
                    atan2(force(cvec_new - cvec_old).modul(),
                          force(cvec_new + cvec_old).modul());
    de_report << setw(8) << rank + 1 << setw(8) << idx << setw(14) << score
              << setw(14) << Energies[idx].elec_ener << setw(14)
              << Energies[idx].vdw_ener << setw(14)
              << Energies[idx].elec_ener - energy_ref.elec_ener << setw(14)
              << Energies[idx].vdw_ener - energy_ref.vdw_ener << setw(10) << rho
              << setw(12) << cangle << "\n";
  }
  if (scored.size() > n_print)
    de_report << "  ... (" << scored.size() - n_print
              << " more candidates omitted)\n";
  de_report << "==============================================================="
               "=================\n";

  // --- DE_out_vector.dat: enriched with atom info ---
  de_output << fixed;
  de_output.precision(4);
  de_output << dX << "  " << dY << "  " << dZ << "  Angle: " << best_angle
            << " -- Modul: " << sqrt(dX * dX + dY * dY + dZ * dZ)
            << " Angstrom\n";
  de_output << "# Pull direction: from " << rec.info[nearest_rec_idx].get_atom()
            << " " << rec.info[nearest_rec_idx].get_original_res_name()
            << rec.info[nearest_rec_idx].get_respos() << " (chain "
            << rec.info[nearest_rec_idx].get_chainnum() << ") towards "
            << rec.info[target_rec_idx].get_atom() << " "
            << rec.info[target_rec_idx].get_original_res_name()
            << rec.info[target_rec_idx].get_respos() << " (chain "
            << rec.info[target_rec_idx].get_chainnum() << ")\n";
  de_output << "# Lig COM: (" << lg_com_x << ", " << lg_com_y << ", "
            << lg_com_z << ") -> (" << best_com_x << ", " << best_com_y << ", "
            << best_com_z << ")\n";
  de_output << "# Best score: " << min_e
            << "  elec: " << Energies[pos].elec_ener
            << "  vdw: " << Energies[pos].vdw_ener << "\n";

  // --- DE_log.dat: append ---
  de_log << fixed;
  de_log.precision(4);
  de_log << dX << "  " << dY << "  " << dZ << "  Angle: " << best_angle
         << " -- Modul: " << sqrt(dX * dX + dY * dY + dZ * dZ)
         << " Angstrom  Score: " << min_e << "\n--------------------------\n";

  // --- best_conf.pdb ---
  best_conf.precision(3);
  for (size_t i{}; i < lg_candidates.m_mol[pos].size(); ++i)
    best_conf << "ATOM" << setw(7) << i + 1 << setw(5) << lg.info[i].get_atom()
              << setw(4) << lg.info[i].get_original_res_name() << " A"
              << setw(4) << 1 << setw(12) << fixed
              << lg_candidates.m_mol[pos].pos.x[i] << setw(8) << fixed
              << lg_candidates.m_mol[pos].pos.y[i] << setw(8) << fixed
              << lg_candidates.m_mol[pos].pos.z[i] << setw(7) << 1.00 << setw(6)
              << 0.00 << '\n';

  // --- DE_vector.pdb ---
  de_pdb_out.precision(3);
  de_pdb_out << "ATOM" << setw(7) << 1 << setw(5) << "CA" << setw(4) << "ALA"
             << " A" << setw(4) << 1 << setw(12) << fixed << lg_com_x << setw(8)
             << fixed << lg_com_y << setw(8) << fixed << lg_com_z << " "
             << setw(7) << 1.00 << " " << setw(6) << 0.00 << '\n';
  de_pdb_out << "ATOM" << setw(7) << 2 << setw(5) << "CB" << setw(4) << "ALA"
             << " A" << setw(4) << 1 << setw(12) << fixed << lg_com_x + dX
             << setw(8) << fixed << lg_com_y + dY << setw(8) << fixed
             << lg_com_z + dZ << " " << setw(7) << 1.00 << " " << setw(6)
             << 0.00 << '\n';
  de_pdb_out << "TER\n";
  return 0;
}
