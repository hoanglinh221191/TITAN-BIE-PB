#ifndef MOL_H
#define MOL_H

#include "COM.h"
#include "MOL_info.h"
#include "MOL_para.h"
#include "amino_a.h"
#include "arguments.h"
#include "box_type.h"
#include "c6c12.h"
#include "differential_evolution.h"
#include "force.h"
#include "h_bond.h"
#include "index_v2.h"
#include "miscellaneous.h"
#include "uti_func.h"
#include <string_view>
using std::string_view;

#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <vector>

// Tránh dùng using namespace std; trong file header!
using std::pair;
using std::string;
using std::vector;

class MOL;

void solve_pbc(MOL &mol_1, MOL &mol_2, box &box_type);
class Multi_MOL;

class MOL {
private: // Các hàm helper nội bộ nên để private
  int cal_electric_field_1atom(
      const vector<mol_info> &info, const arguments &Arguments,
      const vector<Multi_MOL_parameters> &parameters,
      force &electric_field_vector,
      const vector<size_t> &list_heavy_atoms, // OPT: Const correctness
      size_t pos_o_water) const;              // OPT: Const correctness
  int helper_cal_electric_field(
      const vector<mol_info> &info, const arguments &Arguments,
      const vector<Multi_MOL_parameters> &parameters, vector<force> &elec_field,
      const vector<size_t> &list_heavy_atoms, // OPT: Const correctness
      const vector<size_t> &pos_o_waters,
      size_t begin,      // OPT: Const correctness
      size_t end) const; // OPT: Const correctness
  int cal_dipoles_specific_atoms(
      const vector<size_t> &list_atoms, // OPT: Const correctness
      const vector<Multi_MOL_parameters> &parameters,
      const vector<mol_info> &info, vector<force> &dipoles_out,
      std::string_view atoms_type)
      const; // OPT: Using std::string_view for zero-copy and const correctness
  int cal_dipoles_water_molecules(
      const vector<size_t> &list_atoms, // OPT: Const correctness
      const vector<Multi_MOL_parameters> &parameters,
      const vector<mol_info> &info,
      vector<force> &dipoles_out) const; // OPT: Const correctness

public:
  // --- TỐI ƯU HÓA BỐ TRÍ BỘ NHỚ (MEMORY ALIGNMENT) ---
  // Sắp xếp biến từ lớn xuống nhỏ để tránh memory padding

  // OPT: Memory Alignment - Reordered members largest to smallest
  vector<double> x;
  vector<double> y;
  vector<double> z;
  vector<com> Com;
  string force_f;
  string mol_type;
  h_bonds_wrapper list_h_bonds;
  box Box;

  // 2. Các biến nguyên thủy (Primitives) (8 bytes -> 1 byte)
  double time_of_snap;
  bool flag_set_time;

  // --- CONSTRUCTORS (Chuẩn Modern C++) ---

  // Default Constructor: Loại bỏ reserve hardcode để tránh bùng nổ RAM
  MOL() noexcept
      : force_f{"AMBER99SB"}, mol_type{"protein"}, time_of_snap{0.0},
        flag_set_time{false} {}

  MOL(const box &box_t, double time,
      std::string_view
          ff) // OPT: std::string_view for zero-copy string references
      : Box{box_t}, force_f{ff}, mol_type{"protein"}, time_of_snap{time},
        flag_set_time{true} {}

  // Constructor that previously took a single molecule (deprecated, keeping for
  // compilation compatibility, but it should ideally take x, y, z)
  MOL(double _x, double _y, double _z, const box &box_t, double time,
      std::string_view ff, std::string_view _mol_type = "protein")
      : Box{box_t}, force_f{ff}, mol_type{_mol_type}, time_of_snap{time},
        flag_set_time{true} {
    x.push_back(_x);
    y.push_back(_y);
    z.push_back(_z);
  }

  // OPT: RULE OF DEFAULT / ZERO: Hãy để trình biên dịch tự sinh ra Copy/Move
  // Constructor.
  MOL(const MOL &) = default;
  MOL(MOL &&) noexcept = default;
  MOL &operator=(const MOL &) = default;
  MOL &operator=(MOL &&) noexcept = default;
  ~MOL() = default;

  // --- CORE METHODS ---

  void clear_MOL() noexcept {
    // Cách clear memory vector triệt để nhất (Shrink to fit)
    x.clear();
    x.shrink_to_fit();
    y.clear();
    y.shrink_to_fit();
    z.clear();
    z.shrink_to_fit();
    Com.clear();
    Com.shrink_to_fit();
    list_h_bonds.dat.clear();
    list_h_bonds.dat.shrink_to_fit();
    flag_set_time = false;
  }

  size_t size() const {
#ifndef NDEBUG
    if (x.size() != y.size() || y.size() != z.size()) {
      throw std::runtime_error("MOL::size() mismatch in x, y, z dimensions");
    }
#endif
    return x.size();
  }

  template <typename T> void com_pushback(T &&_com) {
    Com.emplace_back(std::forward<T>(_com));
  }

  com &get_com(size_t i) noexcept { return Com[i]; }

  double get_time() const noexcept { return time_of_snap; }
  void set_time(double _time_of_snap) { time_of_snap = _time_of_snap; }

  string_view get_force_f() const noexcept {
    return force_f;
  } // OPT: Added noexcept
  string_view get_mol_type() const noexcept {
    return mol_type;
  } // OPT: Added const reference getter

  // --- MATH & DISTANCE INLINES (Tối ưu hóa vòng lặp) ---
  // Lưu ý: Chỉ inline những hàm toán học cực ngắn

  inline double dis_between2_atoms(size_t i, size_t j) const {
    double dx = x[i] - x[j];
    double dy = y[i] - y[j];
    double dz = z[i] - z[j];
    return std::sqrt(dx * dx + dy * dy + dz * dz);
  }

  inline double dis_between_atom_and_point(size_t i, double px, double py,
                                           double pz) const {
    double dx = x[i] - px;
    double dy = y[i] - py;
    double dz = z[i] - pz;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
  }

  inline double squared_dis_between_atom_and_point(size_t i, double px,
                                                   double py, double pz) const {
    double dx = x[i] - px;
    double dy = y[i] - py;
    double dz = z[i] - pz;
    return (dx * dx + dy * dy + dz * dz);
  }

  // --- CÁC HÀM XỬ LÝ KHÁC (Khai báo nguyên bản) ---

  int post_proc_hbonds(const MOL &mol, vector<mol_info> &info,
                       vector<mol_info> &counter_info);
  int imp_post_proc(const pair<char, size_t> &atom_dat, const MOL &mol,
                    vector<mol_info> &info, vector<mol_info> &counter_info);
  void take_com(vector<Multi_MOL_parameters> &parameters, size_t i, size_t end,
                vector<class index> &Index, vector<mol_info> &info);

  int create_pbc_images(arguments &Arguments, Multi_MOL &out_mol,
                        const vector<class index> &Index);
  int find_best_image_of_residues(Multi_MOL &pbc_images, size_t chain,
                                  vector<class index> &Index,
                                  vector<size_t> &rejected_images);
  // OPT: Const correctness
  void get_centroid(double *com) const;

  int find_images_suitable(
      size_t chain, size_t res, vector<class index> &Index,
      vector<size_t> &list_image) const; // OPT: Const correctness
  int impl_solv_pbc_clust(arguments &Arguments, vector<class index> &Index);
  double find_dis_res_point(size_t chain, size_t res, Multi_MOL &images,
                            vector<class index> &Index, size_t pos_image,
                            double *point) const; // OPT: Const correctness
  double find_min_dis_2res(size_t chain, size_t chain2, size_t res1,
                           size_t res2, Multi_MOL &images,
                           vector<class index> &Index,
                           size_t pos_image) const; // OPT: Const correctness
  int find_cog_chain(size_t ind, arguments &Arguments,
                     vector<class index> &Index, double *center_of_geometry);
  size_t find_chain_ref_pbc(arguments &Arguments, vector<class index> &Index);

  int cal_dis_charged_atoms(vector<double> &distances,
                            vector<double> &distances_negative_charges,
                            vector<double> &distances_positive_charges,
                            const vector<Multi_MOL_parameters> &parameters,
                            double *com) const; // OPT: Const correctness
  int find_surface_water(
      const vector<size_t> &list_waters,      // OPT: Const vector reference
      const vector<size_t> &list_heavy_atoms, // OPT: Const vector reference
      vector<size_t> &list_surface) const;    // OPT: Const correctness
  int find_sphere_water(
      const vector<size_t> &list_waters,      // OPT: Const vector reference
      const vector<size_t> &list_heavy_atoms, // OPT: Const vector reference
      vector<size_t> &list_surface,
      double d_crit) const; // OPT: Const correctness

  int cal_electric_field(
      const vector<mol_info> &info, const arguments &Arguments,
      const vector<Multi_MOL_parameters> &parameters,
      const vector<size_t> &list_heavy_atoms, // OPT: Const vector reference
      const vector<size_t> &pos_o_waters,     // OPT: Const vector reference
      vector<force> &elec_fields) const;      // OPT: Const correctness
  int cal_angle_dipoles_elecfield(
      const vector<mol_info> &info, const arguments &Arguments,
      vector<force> &elec_fields,
      const vector<Multi_MOL_parameters> &parameters,
      const vector<size_t> &list_heavy_atoms, // OPT: Const vector reference
      const vector<size_t> &pos_o_waters,     // OPT: Const vector reference
      vector<double> &cos_out) const;         // OPT: Const correctness

  double cal_minimum_distance_2mols(const MOL &other,
                                    const arguments &Arguments) const;
  double cal_minimum_distance_mol_atom(const MOL &other,
                                       const arguments &Arguments,
                                       size_t &index) const;
  int cal_min_max_mol(double &min_x, double &max_x, double &min_y,
                      double &max_y, double &min_z, double &max_z) const;

  void take_com_all_atoms(vector<Multi_MOL_parameters> &parameters, size_t I,
                          size_t end, vector<class index> &Index,
                          vector<mol_info> &info, multi_points &com_res);
  int translate(double x, double y, double z);
  bool mindist_mol(const MOL &lig, arguments &Arguments, de_arg &DE_arg) const;
};

#endif