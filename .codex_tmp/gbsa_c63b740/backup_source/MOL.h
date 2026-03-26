#ifndef MOL_H
#define MOL_H
#include"molecule.h"
#include"index_v2.h"
#include"COM.h"
#include"box_type.h"
#include"arguments.h"
#include"h_bond.h"
#include"c6c12.h"
#include"amino_a.h"
#include"MOL_info.h"
#include"MOL_para.h"
#include"uti_func.h"
#include"miscellaneous.h"
#include"force.h"
#include"differential_evolution.h"

void solve_pbc(vector<molecule>& mol_1, vector<molecule>& mol_2, box& box_type);
using namespace std;
class Multi_MOL;
class MOL
{
	int cal_electric_field_1atom(const vector<mol_info>& info, const arguments& Arguments, 
		const vector<Multi_MOL_parameters>& parameters, force& electric_field_vector, vector<size_t>& list_heavy_atoms, size_t pos_o_water);
	int helper_cal_electric_field(const vector<mol_info>& info, const arguments& Arguments,
		const vector<Multi_MOL_parameters>& parameters, vector<force>& elec_field, vector<size_t>& list_heavy_atoms, vector<size_t>& pos_o_waters, size_t begin, size_t end);
	int cal_dipoles_specific_atoms(vector<size_t>& list_atoms, const vector<Multi_MOL_parameters>& parameters, const vector<mol_info>& info,
		vector<force>& dipoles_out, string atoms_type);
	int cal_dipoles_water_molecules(vector<size_t>& list_atoms, const vector<Multi_MOL_parameters>& parameters, const vector<mol_info>& info,
		vector<force>& dipoles_out);
public:
	vector<molecule> Molecule;
	vector<com> Com;
	h_bonds_wrapper list_h_bonds;
	box Box;
	string force_f, mol_type;
	double time_of_snap;
	bool flag_set_time;
	MOL() noexcept :time_of_snap{ 0.0 }, flag_set_time{ false }, force_f{ "AMBER99SB" }, mol_type{"protein"}, list_h_bonds()
	{
		Molecule.reserve(200000);
		Com.reserve(2000);
	}
	MOL(box& box_t, double time, string ff) :time_of_snap{ time }, flag_set_time{ true }, force_f{ ff }, mol_type{ "protein" }, list_h_bonds()
	{
		Box = box_t;
		Molecule.reserve(200000);
		Com.reserve(2000);
	}
	MOL(molecule&& mol_in, box& box_t, double time, string ff) :Box{ box_t }, time_of_snap{ time }, flag_set_time{ true }, force_f{ ff }, mol_type{ "protein" }, list_h_bonds() { Molecule.emplace_back(std::move(mol_in)); }
	MOL(molecule& mol_in, box& box_t, double time, string ff) :Box{ box_t }, time_of_snap{ time }, flag_set_time{ true }, force_f{ ff }, mol_type{ "protein" }, list_h_bonds() { Molecule.emplace_back(mol_in); }
	MOL(molecule& mol_in, box& box_t, double time, string ff, string _mol_type) :Box{ box_t }, time_of_snap{ time }, flag_set_time{ true }, force_f{ ff }, mol_type{ _mol_type }, list_h_bonds() { Molecule.emplace_back(mol_in); }
	MOL(molecule&& mol_in, box& box_t, double time, string ff, string _mol_type) :Box{ box_t }, time_of_snap{ time }, flag_set_time{ true }, force_f{ ff }, mol_type{ _mol_type }, list_h_bonds() { Molecule.emplace_back(move(mol_in)); }
	MOL(const MOL& _mol) :time_of_snap{ _mol.time_of_snap }, force_f{_mol.force_f}, mol_type{_mol.mol_type}, Box{_mol.Box}, Molecule{_mol.Molecule}
		, Com{ _mol.Com }, list_h_bonds{ _mol.list_h_bonds }, flag_set_time{true} {}
	~MOL() = default;
	MOL& operator=(const MOL& rhs);
	MOL& operator=(MOL&& rhs) noexcept;

	inline void clear_MOL() noexcept
	{
		vector<molecule>().swap(Molecule);
		vector<com>().swap(Com);
		list_h_bonds.dat.clear();
		flag_set_time = false;
	}

	molecule& get_mol(size_t i)
	{
		if (i > Molecule.size() || Molecule.empty()) { EXIT_FAILURE; }
		return Molecule[i];
	}

	inline size_t get_mol_size() const noexcept
	{
		return Molecule.size();
	}

	template<typename T>
	inline void com_pushback(T&& _com)
	{
		Com.emplace_back(forward<T>(_com));
	}
	inline com& get_com(size_t i) noexcept
	{
		if (i > Com.size() || Com.empty()) { EXIT_FAILURE; }
		return Com[i];
	}

	inline double get_time() const noexcept
	{
		return time_of_snap;
	}
	template<typename T>
	inline void set_time(T&& _time_of_snap)
	{
		time_of_snap = std::forward<T>(_time_of_snap);
	}
	string get_force_f() { return force_f; }
	int post_proc_hbonds(const MOL& mol, vector<mol_info>& info, vector<mol_info>& counter_info);
	int imp_post_proc(const pair<char, size_t>& atom_dat, const MOL& mol, vector<mol_info>& info, vector<mol_info>& counter_info);
	void take_com(vector<Multi_MOL_parameters>& parameters, size_t i, size_t end, vector<class index>& Index, 
		vector<mol_info>& info);
	int create_pbc_images(arguments& Arguments, Multi_MOL& out_mol, const vector<class index>& Index);
	int find_best_image_of_residues(Multi_MOL& pbc_images, size_t chain, vector<class index>& Index, vector<size_t>& rejected_images);
	inline double dis_between2_atoms(size_t i, size_t j)
	{ return sqrt(square(Molecule[i].get_x() - Molecule[j].get_x()) + square(Molecule[i].get_y() - Molecule[j].get_y()) + square(Molecule[i].get_z() - Molecule[j].get_z())); }
	inline double dis_between_atom_and_point(size_t i, double x, double y, double z)
	{
		return sqrt(square(Molecule[i].get_x() - x) + square(Molecule[i].get_y() - y) + square(Molecule[i].get_z() - z));
	}
	inline double squared_dis_between_atom_and_point(size_t i, double x, double y, double z)
	{
		return (square(Molecule[i].get_x() - x) + square(Molecule[i].get_y() - y) + square(Molecule[i].get_z() - z));
	}
	void get_centroid(double* com);
	double find_min_dis_2res(size_t chain, size_t chain2, size_t res1, size_t res2, Multi_MOL& images, vector<class index>& Index, size_t pos_image);
	int find_images_suitable(size_t chain, size_t res, vector<class index>& Index, vector<size_t>& list_image);
	int impl_solv_pbc_clust(arguments& Arguments, vector<class index>& Index);
	double find_dis_res_point(size_t chain, size_t res, Multi_MOL& images, vector<class index>& Index, size_t pos_image, double* point);
	int find_cog_chain(size_t ind, arguments& Arguments, vector<class index>& Index, double* center_of_geometry);
	size_t find_chain_ref_pbc(arguments& Arguments, vector<class index>& Index);
	int cal_dis_charged_atoms(vector<double>& distances, vector<double>& distances_negative_charges, vector<double>& distances_positive_charges,
		const vector<Multi_MOL_parameters>& parameters, double* com);
	int find_surface_water(vector<size_t>& list_waters, vector<size_t>& list_other_atoms, vector<size_t>& list_surface);
	int find_sphere_water(vector<size_t>& list_waters, vector<size_t>& list_heavy_atoms, vector<size_t>& list_surface, double d_crit);
	int cal_electric_field(const vector<mol_info>& info, const arguments& Arguments,
		const vector<Multi_MOL_parameters>& parameters, vector<size_t>& list_heavy_atoms, vector<size_t>& pos_o_waters, vector<force>& elec_fields);
	int cal_angle_dipoles_elecfield(const vector<mol_info>& info, const arguments& Arguments, vector<force>& elec_fields,
		const vector<Multi_MOL_parameters>& parameters, vector<size_t>& list_heavy_atoms, vector<size_t>& pos_o_waters, vector<double>& cos_out);
	double cal_minimum_distance_2mols(MOL& other, const arguments& Arguments);
	double cal_minimum_distance_mol_atom(MOL& other, const arguments& Arguments, size_t& index);
	int cal_min_max_mol(double& min_x, double& max_x, double& min_y, double& max_y, double& min_z, double& max_z);
	void take_com_all_atoms(vector<Multi_MOL_parameters>& parameters, size_t I, size_t end, vector<class index>& Index, vector<mol_info>& info, multi_points& com_res);
	int translate(double x, double y, double z);
	bool mindist_mol(MOL& lig, arguments& Arguments, de_arg& DE_arg);
};

#endif
