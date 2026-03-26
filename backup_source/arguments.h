#ifndef ARGUMENTS_H
#define ARGUMENTS_H
#include<string>
#include<vector>
using namespace std;
struct arguments
{
private:
	string Buildnumber;
	string set_Buildnumber();
	bool parse_input_file_is_open();
public:
	bool cal_dis_sc_contact, opdb, file_par, new_al_1pdb, entropy, entropy_quasi, flag_index, flag_topo_1, flag_topo_2, flag_1_pdb, flag_listfile, flag_nocal
		, flag_pattern_input, flag_range_input, help, flag_pbc_ref, print_energy_res, flag_use_common_topo, flag_elec_screening;
	bool flag_all_folders, flag_separation_results, flag_sidechain_hbond, flag_cal_force, flag_use_CI, print_force_res,
		flag_begin_set, flag_end_set, flag_energy_same_force, flag_maxforce, flag_new_emap, self_int, flag_sidechain, 
		flag_print_sidechain, flag_all_respair;
	bool flag_cat_result_multi_infile, flag_sub_groups, solv_pbc_clust, pbc_envi, flag_cal_contact, pbc_aggressive, flag_cal_dis_res_2com,
		flag_set_type_1, flag_set_type_2, flag_cal_dipole_moment, flag_charged_atom_distance_dis, flag_angle_dipole_efield, flag_sep_opdb,
		heat_capacity_cp, flag_cal_vol, flag_interchain_energy, pbc_mem, flag_pathway, flag_de_auto_term, flag_find_pulling_direct, flag_cal_mindis;
	int flag_pqr, flag_all_side_chain_ct, begin_mol_2, solv_pbc, energy_map, force_map, argc, group_1, group_2,
		solv_pbc_int, solv_pbc_int2, solv_pbc_int3, flag_h_bond, group_entropy_com, group_entropy_rec, group_entropy_lg, k_neighbor, chain_ref, dt,
		b_input, e_input, step_input, flag_cache, com_dis_energy_2mol, com_dis_2mol, pbc_mol2_aggressive, de_flag_score_function;
	vector<int> sub_group_1, sub_group_2;
	vector<char> pbc_clust_axis;
	vector<string> pbc_envi_chain;
	size_t number_of_files, nt, number_charged_distances_bin, number_bin_angle_dipole_efield, number_mindis_bin;
	size_t de_pop_num, de_algorithm, de_ite, group_pathway;
	double d_standard, begin_time, end_time, angle_h, dis_h, limit_dis, eu_par, fu_par, dis_clust, d_crit_water, distance_energy_cutoff{ -1.1 },
		d_crit_minmax_mol{ 4. }, charged_distances_max, charged_distances_min, salt_concentration, epsilon_in;
	double de_diff_weight, de_cross_prob, de_F, de_scan_radius, de_min_scan_radius, de_restrain_angle, dis_along_path;
	string topo_name, topo_name_2, pdb_name, out_pdb_name, index_name, force_f, command, path_amino_acid, path_nonbonded, name_out_cal_dis, arg_filename, list_filename, emap_type, fmap_type, input_pattern, 
		emap_level, fmap_level, out_motif, eu, fu, ref_file, tu, input_file_type, mol_1_type, mol_2_type, water_type, ref_pdb_name;
	arguments(): Buildnumber{set_Buildnumber()}, flag_topo_1 { false }, flag_topo_2{ false }, cal_dis_sc_contact{ false }, opdb{ false }, file_par{ false }, new_al_1pdb{ false }, flag_pbc_ref{ false },
		entropy{ false }, entropy_quasi{ false }, flag_listfile{ false }, flag_pattern_input{ false }, flag_range_input{ false }, help{ false }, print_energy_res{false}
		, flag_index{ false }, flag_1_pdb{ false }, flag_all_folders{false},flag_separation_results{false}, flag_sidechain_hbond{false}, flag_cal_force{false}, flag_use_CI{false},
		print_force_res{ false }, flag_begin_set{ false }, flag_end_set{ false }, flag_energy_same_force{ false }, flag_maxforce{ false }, flag_new_emap{ false }, flag_cache{ 1 }, self_int{ false },
		flag_sidechain{ true }, flag_print_sidechain{ true }, flag_all_respair{ false }, flag_cat_result_multi_infile {
		false
	}, flag_sub_groups{ false }, solv_pbc_clust{ false }, pbc_envi{ false }, flag_sep_opdb{ false }, flag_cal_dis_res_2com{ false },
		flag_set_type_1{ false }, flag_set_type_2{ false }, flag_cal_contact{ true }, pbc_aggressive{ false }, flag_cal_dipole_moment{false}, flag_use_common_topo{false},
		flag_charged_atom_distance_dis{ false }, flag_angle_dipole_efield{ false }, heat_capacity_cp{false}, flag_cal_vol{false}, flag_interchain_energy{false}, pbc_mem{false},
		flag_elec_screening{false}, flag_cal_mindis{false},
		flag_pqr{ 2 }, flag_all_side_chain_ct{ 1 }, begin_mol_2{}, solv_pbc{ 0 }, solv_pbc_int{ 0 }, solv_pbc_int2{ 0 }, solv_pbc_int3{}, energy_map{ 2 }, force_map{ 2 }, argc{ 1 }, group_1{}, group_2{}, b_input{}, e_input{ 100000000 }, step_input{ 1 },
		group_entropy_com{ -1 }, group_entropy_rec{ -1 }, group_entropy_lg{ -1 }, nt{}, number_charged_distances_bin{400}, number_bin_angle_dipole_efield{200}, number_mindis_bin{100}
		, k_neighbor{ 4 }, chain_ref{ -1 }, dt{0}, input_pattern{}, number_of_files{1ull}, com_dis_2mol{0}, pbc_mol2_aggressive{ 0 },
		begin_time{ 0. }, end_time{ 100000000000. }, eu_par{ 1.0 }, fu_par{ 1.66 }, d_crit_water{3.5}, flag_nocal{ false }, limit_dis{ 10e10 }
		, charged_distances_max{ 50. }, charged_distances_min{ 0. }, salt_concentration{ 0.15 }, epsilon_in{1.00000000}
		, d_standard{ 6.5 }, topo_name{ "lig.itp" }, topo_name_2{ "no" }, pdb_name{ "complex.pdb" }, out_pdb_name{ "out_complex.pdb" }, ref_pdb_name{}, index_name{ "index.ndx" }, force_f{ "AMBER99SB" }, command{ "-" }, emap_type{ "residue-atom" }, emap_level{ "per-res" }, fmap_type{ "residue-atom" }, fmap_level{ "per-res" }, out_motif{ "*" },
		path_amino_acid{ "lib_ff/" }, path_nonbonded{ "lib_ff/" }, angle_h{ 135. }, dis_h{ 3.5 }, dis_clust{ 6.5 }, flag_h_bond{ 0 }, name_out_cal_dis{ "out_cal_dis_com_" }, arg_filename{ "arg.in" },
		list_filename{ "input.txt" }, eu{ "kJ/mol" }, fu{ "pN" }, ref_file{}, tu{ "ps" }, input_file_type{ ".pdb" }, water_type{"surface"}, com_dis_energy_2mol{ 4 }
		, mol_1_type{ "protein" }, mol_2_type{"lig"}, de_diff_weight{ 0.9 }, de_cross_prob{ 0.9 }, de_F{ 0.5 }, de_scan_radius{ 2. }, de_min_scan_radius{ 0.5 },
		de_restrain_angle{ 65. }, de_pop_num{ 6480 }, de_algorithm{}, de_ite{ 100 }, de_flag_score_function{ 1 }, group_pathway{},
		dis_along_path{ 4. }, flag_de_auto_term{ false }, flag_find_pulling_direct{ false }, flag_pathway{ false } {
	}
	int post_process();
	string get_build_number()const;
	int parse_arguments(int argc, char** argv);
	int parse_arguments_interactive(int argc);
	int parse_arg_file();
	bool check_topo(string NameofTopoFile, int Topo);
	int help_print();
	void set_number_of_files(size_t num) { number_of_files = num; }
	void setPdbname(string name) { pdb_name = name; }
	string get_pdbOutput() const;
	int recheck_and_change_pdb_input();
	string get_input_file_type()const { return input_file_type; }
	int add_input_file_type()
	{
		if (pdb_name.length() <= 4 || pdb_name.substr(pdb_name.length() - 4) != get_input_file_type())
		{
			pdb_name += get_input_file_type();
			return 0;
		}
		else return 1;
	}
};
#endif