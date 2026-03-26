#include"Multi_MOL.h"
#include<iostream>
#include<cmath>
#include<iomanip>
#include<fstream>
#include<thread>
#include<algorithm>
#include<ctime>
#include<functional>
#include"info.h"


int sup_cal_dis_res_com(vector<multi_points>& com_res, vector<double>& distances, vector<string>& name_res, bool flag_o);

void Multi_MOL::post_process_coordinates_user_define(string type)
{
	mol_type = type;
	generate_index();
}

Multi_MOL::Multi_MOL(const Multi_MOL& rhs) :mol_type{ rhs.mol_type }, force_f{ rhs.force_f }, Index{ rhs.Index }
, parameters{ rhs.parameters }, info{ rhs.info }, m_mol{ rhs.m_mol } {}
Multi_MOL::Multi_MOL(Multi_MOL&& rhs) noexcept : mol_type{ move(rhs.mol_type) }, force_f{ move(rhs.force_f) }, Index{ move(rhs.Index) }
, parameters{ move(rhs.parameters) }, info{ move(rhs.info) }, m_mol{ move(rhs.m_mol) } {}

Multi_MOL& Multi_MOL::operator=(const Multi_MOL& rhs)
{
	this->mol_type = rhs.mol_type;
	this->force_f = rhs.force_f;
	if (!rhs.parameters.empty())
		this->parameters = rhs.parameters;
	this->Index = rhs.Index;
	this->m_mol = rhs.m_mol;
	this->info = rhs.info;
	return *this;
}

Multi_MOL& Multi_MOL::operator=(Multi_MOL&& rhs) noexcept
{
	this->mol_type = move(rhs.mol_type);
	this->force_f = move(rhs.force_f);
	if (!rhs.parameters.empty())
		this->parameters = move(rhs.parameters);
	this->Index = move(rhs.Index);
	this->m_mol = move(rhs.m_mol);
	this->info = move(rhs.info);
	return *this;
}

string Multi_MOL::check_head_tail_res(index& in, size_t pos_res)
{
	if (in.get_type_base(pos_res) == "protein")
	{
		if (info[in.get_res_pos_e(pos_res) - 1].get_atom().find("N") != string::npos
			&& info[in.get_res_pos_e(pos_res)].get_atom().find("H") != string::npos
			&& info[in.get_res_pos_e(pos_res) + 1].get_atom().find("H") != string::npos
			&& info[in.get_res_pos_e(pos_res) + 2].get_atom().find("H") != string::npos)
			return "nterminal";
		size_t remain{ in.get_res_pos_e(pos_res) + 160 >= info.size() ? info.size() : in.get_res_pos_e(pos_res) + 160 }, pos{};
		for (size_t i = in.get_res_pos_e(pos_res) - 1; i < remain; ++i)
		{
			if (info[in.get_res_pos_e(pos_res) - 1].get_resnr() != info[i].get_resnr())
			{
				pos = i;
				break;
			}

		}
		if (pos == 0)pos = info.size();
		if (info[pos - 3].get_atom().find("C") != string::npos
			&& info[pos - 1].get_atom().find("O") != string::npos
			&& info[pos - 2].get_atom().find("O") != string::npos)
			return "cterminal";
		return "normal";
	}
	if (in.get_type_base(pos_res) == "dna" || in.get_type_base(pos_res) == "rna")
	{
		int flag_poo{}, flag_h3t{}, flag_h5t{};
		size_t remain{ in.get_res_pos_e(pos_res) + 50 >= info.size() ? info.size() : in.get_res_pos_e(pos_res) + 50 };
		for (size_t i = in.get_res_pos_e(pos_res) - 1; (i < remain) && (info[in.get_res_pos_e(pos_res) - 1].get_resnr() == info[i].get_resnr()); ++i)
		{
			if (info[i].get_atom() == "P" || info[i].get_atom() == "O1P" || info[i].get_atom() == "O2P")
				++flag_poo;
			if (info[i].get_atom() == "H3T")
				++flag_h3t;
			if (info[i].get_atom() == "H5T")
				++flag_h5t;
		}
		if (flag_h3t == 0 && 1 == flag_h5t)
			return "5terminal";
		if (1 == flag_h3t && 1 == flag_h5t)
			return "nsingle";
		if (3 == flag_poo && 0 == flag_h3t)
			return "normal";
		if (3 == flag_poo && 1 == flag_h3t)
			return "3terminal";
	}
	return "ligand";
}
int Multi_MOL::emplace(size_t i, box& box_type, double time, string force_f, mol_info&& info_in, molecule&& mol_in)
{
	this->force_f = force_f;
	if (i < m_mol.size())
		m_mol[i].Molecule.emplace_back(move(mol_in));
	else
	{
		m_mol.emplace_back(move(mol_in), box_type, time, force_f);
	}

	if (m_mol.size() <= 1)
		info.emplace_back(move(info_in));
	return 0;
}

int Multi_MOL::emplace(size_t i, box& box_type, double time, string force_f, mol_info& info_in, molecule& mol_in, string _mol_type)
{
	this->force_f = force_f;
	if (i < m_mol.size())
		m_mol[i].Molecule.emplace_back(mol_in);
	else
	{
		m_mol.emplace_back(mol_in, box_type, time, force_f, _mol_type);
	}

	if (m_mol.size() <= 1)
		info.emplace_back(info_in);
	return 0;
}

int Multi_MOL::emplace(size_t i, box& box_type, double time, string force_f, mol_info& info_in, molecule&& mol_in, string _mol_type)
{
	this->force_f = force_f;
	if (i < m_mol.size())
		m_mol[i].Molecule.emplace_back(move(mol_in));
	else
	{
		m_mol.emplace_back(move(mol_in), box_type, time, force_f, _mol_type);
	}

	if (m_mol.size() <= 1)
		info.emplace_back(info_in);
	return 0;
}

void Multi_MOL::pbc_prob(Multi_MOL& lig, size_t i, size_t num_threads)
{
	if (num_threads == 0)
		num_threads = 1;
	for (size_t j{ i * m_mol.size() / num_threads }; j < (i + 1)*m_mol.size() / num_threads; ++j)
		solve_pbc(m_mol[j].Molecule, lig.m_mol[j].Molecule, m_mol[j].Box);
}

void Multi_MOL::imp_post_proc(int begin, int end, string type, const amino_acid& amine)
{
	if (type == "AMBER")
	{
		post_process_coordinates(amine);
	}
	else if (type == "GROMOS")
	{
		post_process_coordinates_gromos(amine);
	}
	else if (type == "lig")
	{
		post_process_coordinates_lig();
	}
	else if (type == "CHARMM")
	{
		post_process_coordinates_charmm(amine);
	}
	else if (type == "protein-user" || type == "complex")
	{
		post_process_coordinates_user_define(type);
	}
	else // (type == "OPLSAA")
	{
		post_process_coordinates_oplsaa(amine);
	}
}
int Multi_MOL::wrapper_impl_post_proc(arguments& Arguments, string type, const amino_acid& amine)
{
	if (m_mol.empty())
		return 1;
	this->imp_post_proc(0, 1, type, amine);

	if (Arguments.flag_sidechain == true)
		this->check_side_chain(Arguments);
	return 0;
}
int Multi_MOL::add_coms(bool flag, size_t snap_index)
{
	if (flag == false)
	{
		for (size_t i{}; i < Index.size() - 1; ++i)
			for (size_t j{}; j < Index[i].get_size_v(); ++j)
			{
				m_mol[snap_index].Com.emplace_back(Index[i].get_res_x_e(j), Index[i].get_res_y_e(j), Index[i].get_res_z_e(j), Index[i].get_chain(),
					Index[i].get_res_name_e(j), Index[i].get_res_num_e(j), Index[i].get_res_pos_e(j), Index[i].get_chainid(), Index[i].get_res_pos_e(j));
			}
		return 0;
	}
	else
	{
		for (size_t j{}; j < Index[0].get_size_v(); ++j)
		{
			m_mol[snap_index].Com.emplace_back(Index[0].get_res_x_e(j), Index[0].get_res_y_e(j), Index[0].get_res_z_e(j), Index[0].get_chain()
				, Index[0].get_res_name_e(j), Index[0].get_res_num_e(j), Index[0].get_res_pos_e(j), Index[0].get_chainid(), Index[0].get_res_pos_e(j));
		}
		return 0;
	}
}
int Multi_MOL::gen_vector_parameters()
{
	if (m_mol.empty() || !parameters.empty())
		return 1;
	else
	{
		this->parameters = vector<Multi_MOL_parameters>(info.size());
		return 0;
	}
}

int Multi_MOL::com_mol(const size_t num_snap, std::span<double> com) const
{
	double mol_mass{};
	for (size_t i{}; i < m_mol[num_snap].get_mol_size(); ++i)
	{
		double mass{ parameters[i].get_mass() };
		com[0] += m_mol[num_snap].Molecule[i].get_x()*mass;
		com[1] += m_mol[num_snap].Molecule[i].get_y()*mass;
		com[2] += m_mol[num_snap].Molecule[i].get_z()*mass;
		mol_mass += mass;
	}
	for (size_t i{}; i < 3; ++i) com[i] /= mol_mass;
	return 0;
}

int Multi_MOL::check_side_chain(arguments& Arguments)
{
	if (Arguments.nt < 2 || Index.size()<5)
		this->helper_check_side_chain(Arguments, 0, Index.size() - 1);
	else
	{
		vector<thread> threads(Arguments.nt);
		for (size_t i{}; i < Arguments.nt; ++i)
			threads[i] = thread(&Multi_MOL::helper_check_side_chain, &(*this), ref(Arguments), i*(Index.size() - 1) / Arguments.nt, (i + 1)*(Index.size() - 1) / Arguments.nt);
		for (auto& t : threads)t.join();
	}
	return 0;
}

int Multi_MOL::helper_check_side_chain(arguments& Arguments, size_t begin, size_t end)
{
	for (size_t i{ begin }; i < end; ++i)
		for (size_t j{}; j < Index[i].get_size_v(); ++j)
			if (Arguments.flag_sidechain == true && Index[i].get_type_base(j) == "protein")
			{
				size_t end_pos{ j + 1 >= Index[i].get_size_v() ? Index[i + 1].get_res_pos_e(0) - 1 : Index[i].get_res_pos_e(j + 1) - 1 };
				for (size_t k{ Index[i].get_res_pos_e(j) - 1 }; k < end_pos; ++k)
					info[k].check_side_chain();
			}
	return 0;
}

int Multi_MOL::set_chainnum()
{
	size_t chainid{}, anum{};
	for (size_t i{}; i < this->Index.size() - 1; ++i)
	{
		for (size_t j{}; j < this->Index[i].get_size_v(); ++j)
		{
			size_t max_pos{ j != this->Index[i].get_size_v() - 1 ? this->Index[i].get_res_pos_e(j + 1) - 1 : this->Index[i + 1].get_res_pos_e(0) - 1 };
			for (size_t k = this->Index[i].get_res_pos_e(j) - 1; k < max_pos ; ++k)
			{
				if (anum >= this->info.size()) {
					cout << "overflow info! atom num: " << anum << "; info size: " << this->info.size() << "\n"; break;
				}
				this->info[anum].set_chainnum(chainid);
				this->info[anum].set_respos(this->Index[i].get_res_pos_e(j));
				this->Index[i].set_chainid(chainid);
				++anum;
			}
		}
		++chainid;
	}
	//cout << "anum: " << anum << " info size: " << info.size() << '\n';
	return 0;
}

int Multi_MOL::com_dis_2mol(Multi_MOL& rhs, arguments& Arguments, const ::info& info_output, fstream& output_dis)
{
	cout << "calculating COM distance......";
	for (size_t i{}; i < m_mol.size(); ++i)
	{
		if (i % 200 == 0)
			cout << i << " ";
		if (impl_com_dis_2mol(rhs, Arguments, info_output, output_dis, i) == 1)
			break;
	}
	return 0;
}

int Multi_MOL::impl_com_dis_2mol(Multi_MOL& rhs, arguments& Arguments, const ::info& info_output, fstream& output_dis, size_t i)
{
	double com_mol_1[3]{}, com_mol_2[3]{};

	com_mol(i, com_mol_1);

	if (i >= rhs.m_mol.size())
		return 1;
	else
		rhs.com_mol(i, com_mol_2);

	double dis{ sqrt(square(com_mol_1[0] - com_mol_2[0]) + square(com_mol_1[1] - com_mol_2[1]) + square(com_mol_1[2] - com_mol_2[2])) };
	if (Arguments.com_dis_2mol == 1)
		dis = fabs(com_mol_1[0] - com_mol_2[0]);
	else if (Arguments.com_dis_2mol == 2)
		dis = fabs(com_mol_1[1] - com_mol_2[1]);
	else if (Arguments.com_dis_2mol == 3)
		dis = fabs(com_mol_1[2] - com_mol_2[2]);
	else {}

	output_dis << setw(10) << left << i * info_output.get_current_dt() + info_output.get_current_begin_time()
		<< setw(12) << right << dis;

	double min_dist[3]{ 10e10,10e10,10e10 };

	for (size_t j{}; j < rhs.m_mol[i].Molecule.size(); ++j)
	{
		if (min_dist[0] > fabs(rhs.m_mol[i].Molecule[j].get_x() - com_mol_1[0]))
			min_dist[0] = fabs(rhs.m_mol[i].Molecule[j].get_x() - com_mol_1[0]);
		if (min_dist[1] > fabs(rhs.m_mol[i].Molecule[j].get_y() - com_mol_1[1]))
			min_dist[1] = fabs(rhs.m_mol[i].Molecule[j].get_y() - com_mol_1[1]);
		if (min_dist[2] > fabs(rhs.m_mol[i].Molecule[j].get_z() - com_mol_1[2]))
			min_dist[2] = fabs(rhs.m_mol[i].Molecule[j].get_z() - com_mol_1[2]);
	}

	output_dis << setw(12) << min_dist[0] << setw(12) << min_dist[1] << setw(12) << min_dist[2] << '\n';
	return 0;
}

int Multi_MOL::find_water(vector<size_t>& list_atoms_water, vector<size_t>& list_heavy_atoms)
{
	for (size_t i{}; i < Index.size() - 1; ++i)
	{
		for (size_t j{}; j < Index[i].get_size_v(); ++j)
		{
			size_t remain_index{ j + 1 == Index[i].get_size_v() ? Index[i + 1].get_res_pos_e(0) - 1 : Index[i].get_res_pos_e(j + 1) - 1 };
			if (Index[i].get_type_base(j) == "water")
			{
				for (size_t k{ Index[i].get_res_pos_e(j) - 1 }; k < remain_index; ++k)
					if ((this->force_f.find("OPLSAA")!=string::npos && this->info[k].get_atom_type()!= "opls_114") || this->info[k].get_atom_type()[0] != 'H')
						list_atoms_water.emplace_back(k);
			}
			
			if (Index[i].get_type_base(j) != "water" && Index[i].get_type_base(j) != "ion")
			{
				for (size_t k{ Index[i].get_res_pos_e(j) - 1 }; k < remain_index; ++k)
					if ((this->force_f.find("OPLSAA") != string::npos && this->info[k].get_atom_type() != "opls_114") || this->info[k].get_atom_type()[0] != 'H')
						list_heavy_atoms.emplace_back(k);
			}
		}
	}

	if (list_atoms_water.empty())
	{
		cout << "WARNING!! No water!\n";
	}

	return 0;
}

int Multi_MOL::cal_angle_dipole_water_elecField(const arguments& Arguments, ::info& output_info)
{
	vector<size_t> list_waters, list_heavy_atoms;

	this->find_water(list_waters, list_heavy_atoms);
	vector<double> cos_angle_dipoles_efield;

	for (size_t i{}; i < m_mol.size(); ++i)
	{
		vector<size_t> list_surface_waters;

		if (Arguments.water_type == "surface")
			m_mol[i].find_surface_water(list_waters, list_heavy_atoms, list_surface_waters);
		else
			m_mol[i].find_sphere_water(list_waters, list_heavy_atoms, list_surface_waters, Arguments.d_crit_water);

		vector<force> elec_fields(list_surface_waters.size());
		//cout <<"water: " << list_waters.size()<< " atoms. water sphere: " << list_surface_waters.size() << " atoms. Heavy: " << list_heavy_atoms.size() << " atoms\n";
		//cout <<"crit dis: "<< Arguments.d_crit_water << '\n';

		m_mol[i].cal_electric_field(info, Arguments, parameters, list_heavy_atoms, list_surface_waters, elec_fields);

		m_mol[i].cal_angle_dipoles_elecfield(info, Arguments, elec_fields, parameters, list_heavy_atoms, list_surface_waters, cos_angle_dipoles_efield);
	}

	double max{ -1000. }, min{ 1000. };
	for (size_t i{}; i < cos_angle_dipoles_efield.size(); ++i)
	{
		if (cos_angle_dipoles_efield[i] > max)max = cos_angle_dipoles_efield[i];
		if (cos_angle_dipoles_efield[i] < min)min = cos_angle_dipoles_efield[i];
	}
	double delta{ (max - min) / Arguments.number_bin_angle_dipole_efield };
	vector<size_t> histo(Arguments.number_bin_angle_dipole_efield);
	
	for (size_t j{}; j < cos_angle_dipoles_efield.size(); ++j)
	{
		size_t pos{ static_cast<size_t>((cos_angle_dipoles_efield[j] - min) / delta) };
		if (pos >= histo.size())
			pos = histo.size() - 1;
		++histo[pos];
	}
	
	fstream out_angle_dipole_efield(Arguments.pdb_name.substr(0, Arguments.pdb_name.length() - 4) + output_info.getDistributionAngleDipoleEfield(), ios::out);
	out_angle_dipole_efield << "#Cosine        Population (%)\n";
	for (size_t i{}; i < histo.size(); ++i)
		out_angle_dipole_efield << setw(15) << left << i * delta + min << " " << right << setw(15) << histo[i] *100./cos_angle_dipoles_efield.size()<< '\n';
	return 0;
}


bool Multi_MOL::check_and_pop_size_2last_snaps(size_t& I, string name )
{
	if (this->m_mol.size() >= 2
		&& this->m_mol[this->m_mol.size() - 1].Molecule.size() != this->m_mol[this->m_mol.size() - 2].Molecule.size())
	{
		cout << "WARNING!! Snap number" + name << I << " has different number of atoms (" <<
			this->m_mol[this->m_mol.size() - 1].Molecule.size() << " atoms) compared to snap " << I - 1
			<< " (" << this->m_mol[this->m_mol.size() - 2].Molecule.size() << " atoms)! Removed.\n";

		this->m_mol.pop_back();
		I = this->m_mol.size() - 1;
		return false;
	}
	return true;
}

double Multi_MOL::total_charge()
{
	double charge{};
	for (size_t i{}; i < parameters.size(); ++i)
		charge += parameters[i].get_charge();
	if (fabs(charge) < 10e-8)
		charge = 0.0000000;
	return charge;
}

int Multi_MOL::cal_dis_res_2_com(arguments& Arguments, fstream& output)
{
	size_t num_threads{ Arguments.nt };
	vector<thread> workers(num_threads);
	vector<vector<double>> distances;
	vector<string> name_res;
	for (size_t snap_num{}; snap_num < this->m_mol.size(); ++snap_num)
	{
		vector<double> dis_temp;
		vector<multi_points> pos_res(num_threads);
		size_t num_threads_rest(this->Index.size() - 1 > 0 ? this->Index.size() - 1 : 1);
		size_t ite_thread{};
		while (num_threads_rest > 0)
		{
			if (num_threads_rest <= num_threads)
			{
				for (size_t i{}; i < num_threads_rest; ++i)
				{
					workers[i] = thread(&MOL::take_com_all_atoms, this->m_mol[snap_num], ref(this->parameters), i + ite_thread, 
						this->Index[i + ite_thread + 1].get_res_pos_e(0), ref(this->Index), ref(this->info), ref(pos_res[i]));
				}
				for_each(workers.begin(), workers.begin() + num_threads_rest, mem_fn(&thread::join));
				num_threads_rest = 0ull;
			}
			else
			{
				for (size_t i{}; i < num_threads; ++i)
				{
					workers[i] = thread(&MOL::take_com_all_atoms, this->m_mol[snap_num], ref(this->parameters), i + ite_thread,
						this->Index[i + ite_thread + 1].get_res_pos_e(0), ref(this->Index), ref(this->info), ref(pos_res[i]));
				}
				for (auto& t : workers)t.join();
				num_threads_rest -= num_threads;
				ite_thread += num_threads;
			}
		}
		bool flag_o{ snap_num > 0 ? false : true };
		sup_cal_dis_res_com(pos_res, dis_temp, name_res, flag_o);
		distances.emplace_back(dis_temp);
	}

	double delta{ 0.2 }, min_dis{}, max_dis{ 100. };
	if (distances.empty())return 1;
	vector<pair<double, size_t>> tot_dis_res(distances[0].size()), sd_dis_res(distances[0].size());
	
	for (size_t i{}; i < distances[0].size(); ++i)
	{
		for (size_t j{}; j < distances.size(); ++j)
		{
			if (i >= tot_dis_res.size())continue;
			tot_dis_res[i].first += distances[j][i];
			++tot_dis_res[i].second;
		}
		tot_dis_res[i].first /= tot_dis_res[i].second;
		for (size_t j{}; j < distances.size(); ++j)
		{
			if (i >= tot_dis_res.size())continue;
			sd_dis_res[i].first += square(distances[j][i] - tot_dis_res[i].first);
		}
		output << setw(10) << left << name_res[i] << right << setw(12) << tot_dis_res[i].first << setw(12) << sqrt(sd_dis_res[i].first / distances.size()) << '\n';
	}
	return 0;
}

int Multi_MOL::cal_volumes(vector<double>& vol_out)
{
	for (auto& m : this->m_mol)
		vol_out.emplace_back(m.Box.volume());
	return 0;
}