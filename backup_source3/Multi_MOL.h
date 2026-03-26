#ifndef MULTI_MOL_H
#define MULTI_MOL_H

#include "MOL.h"
#include "MOL_info.h"
#include "MOL_para.h"
#include "amino_a.h"
#include "arguments.h"
#include "c6c12.h"
#include "index_v2.h"
#include "info.h"
#include <fstream>
#include <map>
#include <string>
#include <string_view>
#include <vector>

using std::fstream;
using std::map;
using std::string;
using std::string_view;
using std::vector;

struct topo_info_temp;

class Multi_MOL {
private:
  // --- PRIVATE METHODS ---
  // TỐI ƯU: Đã chuyển các tham số lớn (arguments, amino_acid) thành const
  // reference (const &) để tránh copy và bảo vệ dữ liệu.
  int find_best_image(size_t ind, int chain_ref, arguments &Arguments);
  int find_best_image_new(size_t ind, size_t chain_ref, arguments &Arguments,
                          vector<vector<size_t>> &list_chain_combi);
  int impl_solv_pbc_new(Multi_MOL &mol2, size_t begin, size_t end,
                        arguments &Arguments, const amino_acid &amine);
  int impl_solve_pbc_chain(size_t ind, arguments &Arguments,
                           vector<size_t> &rejected_images);
  int find_best_image_of_mol2(Multi_MOL &mol2, size_t i, arguments &Arguments);
  int aggressive_find_best_image_of_mol2(Multi_MOL &mol2, size_t ind,
                                         arguments &Arguments);
  int solv_pbc_envi(size_t ind, arguments &Arguments,
                    vector<size_t> &rejected_images, size_t k,
                    Multi_MOL &pbc_images, double *com_ref_chain);

  // TỐI ƯU: Truyền string bằng string_view siêu nhẹ
  void imp_post_proc(int begin, int end, string_view type,
                     const amino_acid &amine);
  int parse_pdb(std::string_view input_name,
                std::string_view force_f, // OPT: std::string_view zero-copy
                arguments &Arguments, ::info &output_info);

  int post_process_type_res(class index &in, size_t pos_res);
  int parse_type(const amino_acid &amine);
  string check_head_tail_res(class index &in, size_t pos_res);
  int post_process_coordinates(const amino_acid &amine);
  void post_process_coordinates_gromos(const amino_acid &amine);
  void post_process_coordinates_lig();
  void post_process_coordinates_charmm(const amino_acid &amine);
  void post_process_coordinates_oplsaa(const amino_acid &amine);
  void post_process_coordinates_user_define(string_view type);

  int process_atom_from_topo(fstream &input_topo_protein, string &protein_1,
                             string &lig_name, int &i,
                             vector<topo_info_temp> &info_topo);
  int post_process_atom_from_topo(vector<topo_info_temp> &info_topo);
  int helper_check_side_chain(const arguments &Arguments, size_t begin,
                              size_t end);

  int helper_topo_lg_gromacs_format(
      std::string_view input_topo, // OPT: std::string_view
      std::string_view force_f, para &parameters, fstream &input_topo_protein,
      string &ligname, vector<topo_info_temp> &info_topo);
  int helper_topo_lg_parm_format(
      std::string_view input_topo,                // OPT: std::string_view
      std::string_view force_f, para &parameters, // OPT: std::string_view
      fstream &input_topo_protein, string &lig_name);

  int impl_cal_distribution_distance_charged_atoms(
      vector<vector<double>> &distances,
      vector<vector<double>> &distances_negative_charges,
      vector<vector<double>> &distances_positive_charges, size_t begin,
      size_t end);
  int find_water(vector<size_t> &list_atoms, vector<size_t> &list_heavy_atoms);

  int helper_parse_pdb(
      arguments &Arguments, box &Box_temp, bool &flag_set_time,
      double &time_of_snap_temp, string_view input_pdb, size_t &I,
      bool &flag_set_box, ::info &output_info, string_view protein_1,
      std::string_view input_name, // OPT: std::string_view zero-copy
      std::string_view force_f);   // OPT: std::string_view
  int helper_parse_buffer(vector<string> &buffer_string, arguments &Arguments,
                          box &Box_temp, bool &flag_set_time,
                          double &time_of_snap_temp, size_t &I,
                          bool &flag_set_box, ::info &output_info,
                          std::string_view input_name,
                          size_t &count_snap, // OPT: std::string_view
                          bool &flag_tune_snap, size_t &limit_snap,
                          bool &flag_end);

public:
  // --- MEMORY LAYOUT TỐI ƯU (Từ lớn đến nhỏ) ---
  // OPT: Sắp xếp memory layout
  vector<MOL> m_mol;
  vector<Multi_MOL_parameters> parameters;
  vector<class index> Index;
  vector<mol_info> info;
  string mol_type{"protein"};
  string force_f{"AMBER99SB"};

  // --- CONSTRUCTORS (Modern C++) ---
  Multi_MOL() noexcept = default; // Compiler tự động gán mol_type và force_f

  Multi_MOL(string_view _force_f) : force_f{_force_f} {}

  // OPT: RULE OF DEFAULT: Cho phép compiler tự sinh Copy/Move
  Multi_MOL(const Multi_MOL &) = default;
  Multi_MOL(Multi_MOL &&) noexcept = default;
  Multi_MOL &operator=(const Multi_MOL &) = default;
  Multi_MOL &operator=(Multi_MOL &&) noexcept = default;
  ~Multi_MOL() = default;

  // --- CÁC HÀM GET/SET ---
  inline size_t size_index() const noexcept {
    return Index.size();
  } // Thêm const
  class index &get_index(size_t i) { return Index[i]; }
  const class index &get_index(size_t i) const {
    return Index[i];
  } // Hỗ trợ read-only

  void set_ff(string_view force_f_) { force_f = force_f_; }
  void set_mol_type(string_view temp) { mol_type = temp; }

  // TỐI ƯU: Trả về tham chiếu const thay vì copy chuỗi
  string_view get_mol_type() const noexcept { return mol_type; }

  // --- CÁC HÀM XỬ LÝ CHÍNH ---
  int gen_vector_parameters();
  int generate_index();

  int check_dna_or_rna(const amino_acid &amine);
  int post_process_nucleotide();
  int post_process_protein();
  int take_coordinate(std::string_view input_name,
                      std::string_view _force_f, // OPT: std::string_view
                      arguments &Arguments, const amino_acid &amine,
                      ::info &output_info, int flag_raw);

  void take_topo_lg(string &input_name, string &input_topo,
                    arguments &Arguments,
                    std::string_view _force_f, // OPT: std::string_view
                    para &parameter, const amino_acid &amine,
                    ::info &output_info,
                    std::string_view forced_type_mol); // OPT: std::string_view

  void topo_lg(std::string_view input_topo,
               std::string_view force_f, // OPT: std::string_view
               para &parameters);
  int topo_protein(
      std::string_view input_topo,
      std::string_view force_f, // OPT: std::string_view
      para &parameters,
      const vector<topo_info_temp>
          &info_topo); // OPT: const vector reference to fix by-value array trap

  void pbc_prob(Multi_MOL &lig, size_t i, size_t num_threads);
  int solve_pbc_int(size_t init, size_t nt);
  int solve_pbc_new(Multi_MOL &mol2, arguments &Arguments,
                    const amino_acid &amine);
  int solv_pbc_pbcn(Multi_MOL &mol2, arguments &Arguments);

  // TỐI ƯU: Bỏ const vô nghĩa ở std::move, truyền string theo const reference
  int emplace(size_t i, box &box_type, double time,
              std::string_view force_f, // OPT: string_view
              mol_info &&info_in, MOL &&mol_in);
  int emplace(size_t i, box &box_type, double time,
              std::string_view force_f, // OPT: string_view
              const mol_info &info_in, const MOL &mol_in,
              std::string_view _mol_type); // OPT: string_view
  int emplace(size_t i, box &box_type, double time,
              std::string_view force_f, // OPT: string_view
              const mol_info &info_in, MOL &&mol_in,
              std::string_view _mol_type); // OPT: string_view

  double find_min_dis(size_t ind, size_t ind2, size_t chain, size_t chain2,
                      Multi_MOL &images) const; // Thêm const
  double find_min_dis2(size_t ind, size_t ind2, size_t chain, size_t chain2,
                       Multi_MOL &images) const; // Thêm const

  int wrapper_impl_post_proc(arguments &Arguments, string_view type,
                             const amino_acid &amine);
  int add_coms(bool flag, size_t snap_index);
  int take_charge(amino_acid &amine, const para &parameter, size_t end,
                  size_t begin, std::string_view force_f,
                  int flag); // OPT: std::string_view
  int tran_multi_mol(double *TranVec);
  int find_box_from_min_max();
  int align_to_box();
  int com_mol(const size_t num_snap, double *com) const;
  int check_side_chain(const arguments &Arguments);
  int set_chainnum();
  int com_dis_2mol(Multi_MOL &rhs, arguments &Arguments,
                   const ::info &info_output, fstream &output_dis);
  int solve_pbc_cluster(arguments &Arguments, ::info &output_info);
  int init_chains_connect(arguments &Arguments, const amino_acid &amine,
                          map<size_t, size_t> &chains_connect);
  int solv_pbc_ref(arguments &Arguments, const amino_acid &amine,
                   map<size_t, size_t> &chains_connect, size_t ind,
                   size_t ref_ind);
  int impl_com_dis_2mol(Multi_MOL &rhs, const arguments &Arguments,
                        const ::info &info_output, fstream &output_dis,
                        size_t i);
  int cal_dipole_moment(const ::info &info_output, const arguments &Arguments,
                        std::string_view prefix_name); // OPT: std::string_view
  int impl_cal_dipole_moment(const arguments &Arguments, vector<force> &dipoles,
                             size_t begin, size_t end);
  int cal_distribution_distance_charged_atoms(const arguments &Arguments,
                                              const ::info &output_info);
  int cal_angle_dipole_water_elecField(const arguments &Arguments,
                                       ::info &output_info);
  bool check_and_pop_size_2last_snaps(
      size_t &I, std::string_view name = " "); // OPT: string_view zero copy

  double total_charge() const; // Thêm const
  int cal_dis_res_2_com(const arguments &Arguments, fstream &output);
  int cal_volumes(vector<double> &vol_out) const; // Thêm const
  int gen_pop(MOL &mol_seed, de_arg &DE_arguments, const vector<mol_info> &info,
              string_view mol_type, vector<force> &tran_data);
};

#endif // MULTI_MOL_H