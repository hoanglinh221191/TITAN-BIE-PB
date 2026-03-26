#ifndef MULTI_MOL_H
#define MULTI_MOL_H
#include"arguments.h"
#include"index_v2.h"
#include"amino_a.h"
#include"c6c12.h"
#include"MOL.h"
#include"MOL_info.h"
#include"MOL_para.h"
#include"info.h"
#include<map>
class Multi_MOL
{
	int find_best_image(size_t ind, int chain_ref, arguments& Arguments);
	int find_best_image_new(size_t ind, size_t chain_ref, arguments& Arguments, vector<vector<size_t>>& list_chain_combi);
	int impl_solv_pbc_new(Multi_MOL& mol2, size_t begin, size_t end, arguments& Arguments, const amino_acid& amine);
	int impl_solve_pbc_chain(size_t ind, arguments& Arguments, vector<size_t>& rejected_images);
	int find_best_image_of_mol2(Multi_MOL& mol2, size_t i, arguments& Arguments);
	int aggressive_find_best_image_of_mol2(Multi_MOL& mol2, size_t ind, arguments& Arguments);
	int solv_pbc_envi(size_t ind, arguments& Arguments, vector<size_t>& rejected_images, size_t k,
		Multi_MOL& pbc_images, double* com_ref_chain);
	void imp_post_proc(int begin, int end, string type, const amino_acid& amine);
	int parse_pdb(string input_name, string force_f, arguments& Arguments, ::info& output_info);
	int post_process_type_res(class index& in, size_t pos_res);
	int parse_type(const amino_acid& amine);
	string check_head_tail_res(class index& in, size_t pos_res);
	int post_process_coordinates(const amino_acid& amine);
	void post_process_coordinates_gromos(const amino_acid& amine);
	void post_process_coordinates_lig();
	void post_process_coordinates_charmm(const amino_acid& amine);
	void post_process_coordinates_oplsaa(const amino_acid& amine);
	void post_process_coordinates_user_define(string type);
	int process_atom_from_topo(fstream& input_topo_protein, string& protein_1, string& lig_name, int& i, vector<topo_info_temp>& info_topo);
	int post_process_atom_from_topo(vector<topo_info_temp>& info_topo);
	int helper_check_side_chain(arguments& Arguments, size_t begin, size_t end);
	int helper_topo_lg_gromacs_format(string input_topo, string force_f, para &parameters, fstream& input_topo_protein, string& ligname);
	int helper_topo_lg_parm_format(string input_topo, string force_f, para& parameters, fstream& input_topo_protein, string& lig_name);
	int impl_cal_distribution_distance_charged_atoms(vector<vector<double>>& distances, vector<vector<double>>& distances_negative_charges,
		vector<vector<double>>& distances_positive_charges, size_t begin, size_t end);
	int find_water(vector<size_t>& list_atoms, vector<size_t>& list_heavy_atoms);
	int helper_parse_pdb(arguments& Arguments, box& Box_temp, bool& flag_set_time, double& time_of_snap_temp, string_view& input_pdb, size_t& I,
		bool& flag_set_box, ::info& output_info, string& protein_1, string input_name);
	int helper_parse_buffer(vector<string>& buffer_string, arguments& Arguments, box& Box_temp, bool& flag_set_time, double& time_of_snap_temp, size_t& I,
		bool& flag_set_box, ::info& output_info, string input_name, size_t& count_snap, bool& flag_tune_snap, size_t& limit_snap, bool& flag_end);

public:
	vector<MOL> m_mol;
	vector<Multi_MOL_parameters> parameters;
	vector<class index> Index;
	vector<mol_info> info;
	string mol_type, force_f{ "AMBER99SB" };
	Multi_MOL() noexcept :mol_type{ "protein" }
	{
		m_mol.reserve(1000);
		Index.reserve(50);
		info.reserve(200000);
		parameters.reserve(200000);
	}
	Multi_MOL(string _force_f) :mol_type{ "protein" }, force_f{ _force_f }
	{
		m_mol.reserve(1000);
		Index.reserve(50);
		info.reserve(200000);
		parameters.reserve(200000);
	}

	Multi_MOL(const Multi_MOL&);
	Multi_MOL(Multi_MOL&& rhs) noexcept;
	Multi_MOL& operator=(const Multi_MOL& rhs);
	Multi_MOL& operator=(Multi_MOL&& rhs) noexcept;
	~Multi_MOL() = default;

	inline size_t size_index() noexcept
	{
		return Index.size();
	}
	class index& get_index(size_t i)
	{
		return Index[i];
	}
	void set_ff(string force_f_) { force_f = force_f_; }
	int gen_vector_parameters();
	int generate_index();
	void set_mol_type(string temp) { mol_type = temp; }
	string get_mol_type() { return mol_type; }
	int check_dna_or_rna(const amino_acid &amine);
	int post_process_nucleotide();
	int post_process_protein();
	int take_coordinate(string input_name, string _force_f, arguments& Arguments, const amino_acid& amine, ::info& output_info, int flag_raw);
	
	void take_topo_lg(string& input_name, string& input_topo, arguments& Arguments, string _force_f, para& parameter, const amino_acid& amine, 
		::info& output_info, string forced_type_mol);

	void topo_lg(string input_topo, string force_f, para &parameters);
	int topo_protein(string input_topo, string force_f, para& parameters, vector<topo_info_temp> info_topo);
	void pbc_prob(Multi_MOL& lig, size_t i, size_t num_threads);
	int solve_pbc_int(size_t init, size_t nt);
	int solve_pbc_new(Multi_MOL& mol2, arguments& Arguments, const amino_acid& amine);
	int emplace(size_t i, box& box_type, double time, string force_f, mol_info&& info_in, molecule&& mol_in);
	int emplace(size_t i, box& box_type, double time, string force_f, mol_info& info_in, molecule& mol_in, string _mol_type);
	int emplace(size_t i, box& box_type, double time, string force_f, mol_info& info_in, molecule&& mol_in, string _mol_type);
	double find_min_dis(size_t ind, size_t ind2, size_t chain, size_t chain2, Multi_MOL& images);
	double find_min_dis2(size_t ind, size_t ind2, size_t chain, size_t chain2, Multi_MOL& images);
	int wrapper_impl_post_proc(arguments& Arguments, string type, const amino_acid& amine);
	int add_coms(bool flag, size_t snap_index);
	int take_charge(amino_acid &amine, const para &parameter, size_t end, size_t begin, string force_f, int flag);
	int tran_multi_mol(double* TranVec);
	int find_box_from_min_max();
	int align_to_box();
	int com_mol(const size_t num_snap, double* com) const;
	int check_side_chain(arguments& Arguments);
	int set_chainnum();
	int com_dis_2mol(Multi_MOL& rhs, arguments& Arguments, const ::info& info_output, fstream& output_dis);
	int solve_pbc_cluster(arguments& Arguments, ::info& output_info);
	int init_chains_connect(arguments& Arguments, const amino_acid& amine, map<size_t, size_t>& chains_connect);
	int solv_pbc_ref(arguments& Arguments, const amino_acid& amine, map<size_t, size_t>& chains_connect, size_t ind, size_t ref_ind);
	int impl_com_dis_2mol(Multi_MOL& rhs, arguments& Arguments, const ::info& info_output, fstream& output_dis, size_t i);
	int cal_dipole_moment(const ::info& info_output, const arguments& Arguments, string prefix_name);
	int impl_cal_dipole_moment(const arguments& Arguments, vector<force>& dipoles, size_t begin, size_t end);
	int cal_distribution_distance_charged_atoms(const arguments& Arguments, const ::info& output_info);
	int cal_angle_dipole_water_elecField(const arguments& Arguments, ::info& output_info);
	bool check_and_pop_size_2last_snaps(size_t& I, string name = " ");
	double total_charge();
	int cal_dis_res_2_com(arguments& Arguments, fstream& output);
	int cal_volumes(vector<double>& vol_out);
	int gen_pop(MOL& mol_seed, de_arg& DE_arguments, vector<mol_info>& info, string mol_type, vector<force>& tran_data);
};
#endif // !MULTI_MOL_H

