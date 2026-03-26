#include "MOL.h"
#include "differential_evolution.h"
#include "functions.h"
#include "omp.h"
#include <cmath>
#include <iomanip>
#include <random>
#include <string_view>
using std::string_view;
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
  if (Arguments.flag_pathway == false)
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
      if (distance < Arguments.dis_along_path * Arguments.dis_along_path) {
        flag_dis_min = true;
        break;
      }
    }
  }
  return flag_dis_min;
}

double cal_min_dis_ligand_pathway(groups_cal &groups, Multi_MOL &rec,
                                  MOL &lg_candidate, arguments &Arguments) {
  if (Arguments.flag_pathway == false)
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

double cal_min_dis_ligand_pathway_optimized(groups_cal &groups, Multi_MOL &rec,
                                            const std::span<double> com_lg_candidate,
                                            arguments &Arguments) {
  if (Arguments.flag_pathway == false)
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
      new_mol.pos.x.push_back(mol_seed.pos.x[j] + dX);
      new_mol.pos.y.push_back(mol_seed.pos.y[j] + dY);
      new_mol.pos.z.push_back(mol_seed.pos.z[j] + dZ);
    }
    emplace(i, mol_seed.Box, mol_seed.get_time(), mol_seed.get_force_f(),
            info[0], std::move(new_mol), string(mol_type));
  }
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
    for (size_t i{begin}; i < end; ++i) {
      cal_energy_amber(rec.m_mol[0], lg_candidates.m_mol[i], rec.info,
                       lg_seed.info, CachedData, 0, lg_seed.m_mol[0].size(),
                       Energies[i], i, Arguments, rec.parameters,
                       lg_seed.parameters);
    }
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
      if (false == is_valid_angle && false == is_valid_radius)
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
        if (false == is_valid_angle && false == is_valid_radius)
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
        if (false == is_valid_angle && false == is_valid_radius)
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
                               dis_candidate, Arguments.dis_along_path) <
            cal_score_function(Energies[i], energy_ref, flag_score_function,
                               dis_old_mem, Arguments.dis_along_path))) {
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
  vector<thread> threads(Arguments.nt);
  for (size_t i{}; i < threads.size(); ++i)
    threads[i] =
        thread(impl_evo_1_ite, ref(rec), ref(lg_seed), ref(lg_candidates),
               ref(CachedData), ref(DE_arguments), ref(Energies),
               ref(Arguments), i * lg_candidates.m_mol.size() / threads.size(),
               (i + 1) * lg_candidates.m_mol.size() / threads.size(),
               ref(tran_data), ref(out_tran), i + 1, ref(energy_ref),
               ref(Groups_from_index_file), Arguments.de_flag_score_function);
  for (auto &t : threads)
    t.join();
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
  vector<thread> threads(Arguments.nt);

  if (Arguments.force_f.find("CHARMM") != string::npos ||
      Arguments.force_f.find("AMBER") != string::npos ||
      Arguments.force_f.find("OPLSAA") != string::npos) {
    for (size_t i{}; i < threads.size(); ++i)
      threads[i] =
          thread(cal_energy_amber, ref(rec.m_mol[0]), ref(lg.m_mol[0]),
                 ref(rec.info), ref(lg.info), ref(CachedData),
                 i * lg.m_mol[0].size() / threads.size(),
                 (i + 1) * lg.m_mol[0].size() / threads.size(), ref(energy_ref),
                 0, ref(Arguments), ref(rec.parameters), ref(lg.parameters));
    for (auto &t : threads)
      t.join();
  }

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
    // Recompute dX, dY, dZ from tran_data[i] (theta, phi, rho)
    // tran_data stores: Fx=theta, Fy=phi, Fz=rho
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
  for (size_t i{}; i < threads.size(); ++i)
    threads[i] =
        thread(cal_energy_part_DE, ref(rec), ref(lg), ref(lg_candidates),
               ref(CachedData), ref(DE_arguments), ref(Energies),
               ref(Arguments), i * lg_candidates.m_mol.size() / threads.size(),
               (i + 1) * lg_candidates.m_mol.size() / threads.size());
  for (auto &t : threads)
    t.join();

  cout << "Evolution stage........";
  for (size_t i{}; i < DE_arguments.ite_num; ++i) {
    cout << "  " << i + 1;
    evo_1_ite(rec, lg, lg_candidates, CachedData, DE_arguments, Energies,
              Arguments, tran_data, energy_ref, Groups_from_index_file);
  }
  cout << '\n';

  fstream de_output(output_info.get_DE_output_filename(), ios::out),
      de_log(output_info.get_DE_log_filename(), ios::out | ios::app),
      de_pdb_out(output_info.get_DE_pdb_vector_filename(), ios::out | ios::app);
  fstream enery("out_energy_pop.dat", ios::out),
      best_conf("best_conf.pdb", ios::out);

  double min_e{cal_score_function(Energies[0], energy_ref,
                                  Arguments.de_flag_score_function, 0.,
                                  Arguments.dis_along_path)};
  size_t pos{};
  for (size_t i{1}; i < Energies.size(); ++i) {
    enery << i << ": " << Energies[i].elec_ener << "  " << Energies[i].vdw_ener
          << "  " << energy_ref.elec_ener << "  " << energy_ref.vdw_ener;

    double dX{tran_data[i].Fz * sin(tran_data[i].Fx) * cos(tran_data[i].Fy)},
        dY{tran_data[i].Fz * sin(tran_data[i].Fx) * sin(tran_data[i].Fy)},
        dZ{tran_data[i].Fz * cos(tran_data[i].Fx)};
    double modul_old{DE_arguments.old_vector.modul()},
        modul_new{tran_data[pos].Fz};
    force new_vec{dX, dY, dZ}, old_vec{DE_arguments.old_vector};
    new_vec *= modul_old;
    old_vec *= modul_new;
    enery << "  "
          << 360. / PI *
                 atan2(force(new_vec - old_vec).modul(),
                       force(new_vec + old_vec).modul())
          << '\n';

    if (cal_score_function(Energies[i], energy_ref,
                           Arguments.de_flag_score_function, 0.,
                           Arguments.dis_along_path) < min_e) {
      pos = i;
      min_e = cal_score_function(Energies[i], energy_ref,
                                 Arguments.de_flag_score_function, 0.,
                                 Arguments.dis_along_path);
    }
  }

  enery << "MIN ENERGY: " << min_e << ": elec: " << Energies[pos].elec_ener
        << "  vdw: " << Energies[pos].vdw_ener << '\n';

  double dX{tran_data[pos].Fz * sin(tran_data[pos].Fx) *
            cos(tran_data[pos].Fy)},
      dY{tran_data[pos].Fz * sin(tran_data[pos].Fx) * sin(tran_data[pos].Fy)},
      dZ{tran_data[pos].Fz * cos(tran_data[pos].Fx)};
  double modul_old{DE_arguments.old_vector.modul()},
      modul_new{tran_data[pos].Fz};
  force new_vec{dX, dY, dZ}, old_vec{DE_arguments.old_vector};
  new_vec *= modul_old;
  old_vec *= modul_new;

  de_output << dX << "  " << dY << "  " << dZ << "  Angle: "
            << 360. / PI *
                   atan2(force(new_vec - old_vec).modul(),
                         force(new_vec + old_vec).modul())
            << " -- Modul: " << new_vec.modul() << " Angstrom\n";
  de_log << dX << "  " << dY << "  " << dZ << "  Angle: "
         << 360. / PI *
                atan2(force(new_vec - old_vec).modul(),
                      force(new_vec + old_vec).modul())
         << " -- Modul: " << new_vec.modul()
         << " Angstrom\n--------------------------\n";
  double x{}, y{}, z{}, mass{};
  size_t atom{};
  for (size_t i{}; i < lg.m_mol[0].size(); ++i) {
    x += lg.m_mol[0].pos.x[i] * lg.parameters[atom].get_mass();
    y += lg.m_mol[0].pos.y[i] * lg.parameters[atom].get_mass();
    z += lg.m_mol[0].pos.z[i] * lg.parameters[atom].get_mass();
    mass += lg.parameters[atom].get_mass();
    ++atom;
  }
  x /= mass;
  y /= mass;
  z /= mass;
  best_conf.precision(3);
  for (size_t i{}; i < lg_candidates.m_mol[pos].size(); ++i)
    best_conf << "ATOM" << setw(7) << i + 1 << setw(5) << lg.info[i].get_atom()
              << setw(4) << lg.info[i].get_original_res_name() << " A"
              << setw(4) << 1 << setw(12) << fixed
              << lg_candidates.m_mol[pos].pos.x[i] << setw(8) << fixed
              << lg_candidates.m_mol[pos].pos.y[i] << setw(8) << fixed
              << lg_candidates.m_mol[pos].pos.z[i] << setw(7) << 1.00 << setw(6)
              << 0.00 << '\n';

  de_pdb_out.precision(3);
  de_pdb_out << "ATOM" << setw(7) << 1 << setw(5) << "CA" << setw(4) << "ALA"
             << " A" << setw(4) << 1 << setw(12) << fixed << x << setw(8)
             << fixed << y << setw(8) << fixed << z << " " << setw(7) << 1.00
             << " " << setw(6) << 0.00 << '\n';
  de_pdb_out << "ATOM" << setw(7) << 2 << setw(5) << "CB" << setw(4) << "ALA"
             << " A" << setw(4) << 1 << setw(12) << fixed << x + dX << setw(8)
             << fixed << y + dY << setw(8) << fixed << z + dZ << " " << setw(7)
             << 1.00 << " " << setw(6) << 0.00 << '\n';
  de_pdb_out << "TER\n";
  return 0;
}
