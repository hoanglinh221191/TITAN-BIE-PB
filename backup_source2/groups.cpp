#include<fstream>
#include<vector>
#include<string>
#include<sstream>
#include<iostream>
#include<iomanip>
#include<cmath>
#include"Multi_MOL.h"
#include"groups.h"
#include"uti_func.h"
using namespace std;
int s_cal_energy(cache_sigma_epsilon& CachedData, const Multi_MOL &rec, const Multi_MOL &lg, size_t begin_mol_1, size_t end_mol_1, size_t begin_mol_2, size_t end_mol_2, size_t num_snap, double& E_E, double& E_V);
int groups_cal::index(arguments& Arguments)
{
	if (Arguments.flag_index == false)
		return 1;
	fstream index_file(Arguments.index_name, ios::in);
	group_cal group_temp;
	int flag_f{};
	if (!index_file.is_open())
	{
		cout << "index file is not found!\n";
		//Arguments.flag_index = false;
		return 1;
	}
	for (string temp{}; getline(index_file, temp);)
	{
		size_t found_1 = temp.find("[ ");
		size_t found_2 = temp.find(" ]");
		if (found_1 != string::npos && found_2 != string::npos)
		{
			++flag_f;
			if (flag_f >= 2)
			{
				groups.emplace_back(group_temp);
				group_temp.atom_index.clear();
			}
			group_temp.name = temp.substr(2, temp.length() - 4);
		}
		else
		{
			istringstream isstream_temp(temp);
			size_t atom{};
			while (isstream_temp >> atom)
			{
				group_temp.atom_index.emplace_back(atom);
			}
		}
	}
	index_file.close();
	groups.emplace_back(group_temp);
	group_temp.atom_index.clear();
	if (Arguments.group_1 == 0 && Arguments.group_2 == 0 && Arguments.flag_nocal==false)
		this->select_groups_to_cal_energy();
	else
	{
		this->group_1 = Arguments.group_1;
		this->group_2 = Arguments.group_2;

		if (this->group_1 >= this->group_num() || this->group_2 >= this->group_num())
			this->select_groups_to_cal_energy();
	}
	this->groups[group_1].check_conti();
	this->groups[group_2].check_conti();
	if (Arguments.entropy_quasi == true)
	{
		if (Arguments.group_entropy_com != -1)
		{
			this->group_entropy_com = Arguments.group_entropy_com;
			this->group_entropy_rec = Arguments.group_entropy_rec;
			this->group_entropy_lg = Arguments.group_entropy_lg;
		}
		else
			this->select_group_to_cal_entropy();
		this->groups[group_entropy_com].check_conti();
		this->groups[group_entropy_rec].check_conti();
		this->groups[group_entropy_lg].check_conti();
	}
	Arguments.flag_index = true;
	parse_index_flag = true;
	add_sub_groups(Arguments);
	
	this->groups[group_1].gen_setgroup();
	this->groups[group_2].gen_setgroup();

	return 0;
}
void groups_cal::print()
{
	fstream out("out_index_groups.dat", ios::out);
	for (size_t i = 0; i < groups.size(); ++i)
	{
		out << "[ " << groups[i].name << " ]\n";
		for (size_t j = 0; j < groups[i].atom_index.size(); ++j)
			out << groups[i].atom_index[j] << "  ";
		out << '\n';
	}
}
/*int groups_cal::calculate_energy_for_2_groups(const Multi_MOL &rec, const Multi_MOL &lg, cache_sigma_epsilon& CachedData, int begin, int end, size_t num_snap, int group_1, int group_2)
{
	double E_E = 0.0, E_V = 0.0;
	size_t begin_mol_1{}, begin_mol_2{}, end_mol_1{}, end_mol_2{};
	if (this->groups[group_1].atom_index[this->groups[group_1].atom_index.size() - 1] <= rec.m_mol[num_snap].Molecule.size())
	{
		begin_mol_1 = groups[group_1].atom_index[0];
		end_mol_1 = groups[group_1].atom_index[this->groups[group_1].atom_index.size() - 1];
	}
	if (this->groups[group_1].atom_index[0]> rec.m_mol[num_snap].Molecule.size()
		&& this->groups[group_1].atom_index[this->groups[group_1].atom_index.size() - 1] <= (rec.m_mol[num_snap].Molecule.size() + lg.m_mol[num_snap].Molecule.size()))
	{
		begin_mol_2 = groups[group_1].atom_index[0];
		end_mol_2 = groups[group_1].atom_index[this->groups[group_1].atom_index.size() - 1];
	}
	if (this->groups[group_2].atom_index[this->groups[group_2].atom_index.size() - 1] <= (rec.m_mol[num_snap].Molecule.size() + lg.m_mol[num_snap].Molecule.size())
		&& this->groups[group_2].atom_index[0] > rec.m_mol[num_snap].Molecule.size())
	{
		begin_mol_2 = groups[group_2].atom_index[0];
		end_mol_2 = groups[group_2].atom_index[this->groups[group_2].atom_index.size() - 1];
	}
	if (this->groups[group_2].atom_index[this->groups[group_2].atom_index.size() - 1] <= rec.m_mol[num_snap].Molecule.size())
	{
		begin_mol_1 = groups[group_2].atom_index[0];
		end_mol_1 = groups[group_2].atom_index[this->groups[group_2].atom_index.size() - 1];
	}
	if (this->groups[group_1].atom_index[this->groups[group_1].atom_index.size() - 1] <= rec.m_mol[num_snap].Molecule.size()
		&& this->groups[group_2].atom_index[this->groups[group_2].atom_index.size() - 1] <= rec.m_mol[num_snap].Molecule.size())
	{
		if (groups[group_1].atom_index[this->groups[group_1].atom_index.size() - 1]< groups[group_2].atom_index[0]) {
			begin_mol_1 = groups[group_1].atom_index[0];
			end_mol_1 = groups[group_1].atom_index[this->groups[group_1].atom_index.size() - 1];
			begin_mol_2 = groups[group_2].atom_index[0];
			end_mol_2 = groups[group_2].atom_index[this->groups[group_2].atom_index.size() - 1];
		}
		else if (groups[group_1].atom_index[0] > groups[group_2].atom_index[this->groups[group_2].atom_index.size() - 1])
		{
			begin_mol_2 = groups[group_1].atom_index[0];
			end_mol_2 = groups[group_1].atom_index[this->groups[group_1].atom_index.size() - 1];
			begin_mol_1 = groups[group_2].atom_index[0];
			end_mol_1 = groups[group_2].atom_index[this->groups[group_2].atom_index.size() - 1];
		}
		else
			return 1;
	}
	if (end_mol_1 == 1 || end_mol_2 == 1)
		return 1;
	s_cal_energy(CachedData, rec, lg, begin_mol_1, end_mol_1, begin_mol_2, end_mol_2, num_snap, E_E, E_V);
	this->E_elec.emplace_back(E_E);
	this->E_vdw.emplace_back(E_V);
	this->time_of_snap.emplace_back(rec.m_mol[num_snap].get_time());
	return 0;
}*/
int groups_cal::print_result_cal_to_file()
{
	string output_name = "energy_groups_" + to_string(group_1) + "_" + to_string(group_2) + ".dat";
	fstream output_file(output_name, ios::out);
	output_file << "@    title \"Energy interaction between groups " << group_1 << " and " << group_2 << '\n'
		<< "@    xaxis label \"Time (ps)\"\n"
		<< "@    yaxis label \"Energy (kJ/mol)\""
		<< "@ s0 legend \"Time\"\n"
		<< "@ s1 legend \"electrostatic\"\n"
		<< "@ s2 legend \"van der Waal\"\n"
		<< "@   ps  Elec_e(kJ / mol)   VdW_e(kJ / mol) Total_e(kJ / mol)\n";
	for (size_t i = 0; i < E_elec.size(); ++i)
	{
		output_file << setw(6) << time_of_snap[i] << setw(17) << E_elec[i] << setw(17) << E_vdw[i] << setw(17) << E_elec[i] + E_vdw[i] << '\n';
	}
	return 0;
}
int groups_cal::select_groups_to_cal_energy()
{
	cout << "Please choose 2 groups to calculate interaction energy:\n";
	for (int i = 0; i < 2; ++i)
	{
		for (int j = 0; j < groups.size(); ++j)
			cout << setw(4) << j << setw(25) << groups[j].name << ": " << setw(7) << groups[j].atom_index.size() << " atoms\n";
		if (i == 0)
			cin >> group_1;
		else
			cin >> group_2;
	}
	return 0;
}
int groups_cal::select_group_to_cal_entropy()
{
	cout << "Please choose 3 groups to calculate entropy:\n";
	bool test{ true };
	int count{ 0 }, in{};
	do
	{
		for (int j = 0; j < groups.size(); ++j)
			cout << setw(4) << j << setw(25) << groups[j].name << ": " << setw(7) << groups[j].atom_index.size() << " atoms\n";
		cin >> in;
		if (in < groups.size() && in>0)
		{
			++count;
			if (count == 1)group_entropy_com = in;
			else if (count == 2)group_entropy_rec = in;
			else group_entropy_lg = in;
		}
		else
			cout << "Your group is not in index file!\nPlease choose again!\n";
		if (count == 3)test = false;
	} while (test);
	return 0;
}

int s_cal_energy(cache_sigma_epsilon& CachedData, const Multi_MOL &rec, const Multi_MOL &lg, size_t begin_mol_1, size_t end_mol_1, size_t begin_mol_2, size_t end_mol_2, size_t num_snap, double& E_E, double& E_V)
{
	for (size_t j = begin_mol_1 - 1; j < end_mol_1 - 1; ++j)
		for (size_t i = begin_mol_2 - 1; i < end_mol_2 - 1; ++i)
		{
			double d{ 0.1*sqrt(square(rec.m_mol[num_snap].pos.x[i] - lg.m_mol[num_snap].pos.x[j]) 
				+ square(rec.m_mol[num_snap].pos.y[i] - lg.m_mol[num_snap].pos.y[j]) 
				+ square(rec.m_mol[num_snap].pos.z[i] - lg.m_mol[num_snap].pos.z[j])) };
			if (rec.m_mol[num_snap].force_f.find("GROMOS")!=string::npos)
				CachedData.parameters[rec.m_mol[0].Molecule.size()*j + i].set_elec_vdw(CachedData.charge[rec.m_mol[0].Molecule.size()*j + i] / d, 
					sqrt(rec.parameters[i].get_epsilon() * lg.parameters[j].get_epsilon())/ c12(d) - sqrt(rec.parameters[i].get_sigma()*lg.parameters[j].get_sigma()) / c6(d));
			else
				CachedData.parameters[rec.m_mol[0].Molecule.size()*j + i].set_elec_vdw(CachedData.charge[rec.m_mol[0].Molecule.size()*j + i] / d, 
					4.0 * sqrt(rec.parameters[i].get_epsilon()*lg.parameters[j].get_epsilon())*(c12(0.5*(rec.parameters[i].get_sigma() + lg.parameters[j].get_sigma()) / d)
						- c6(0.5 * (rec.parameters[i].get_sigma() + lg.parameters[j].get_sigma()) / d)));
			E_V += CachedData.parameters[rec.m_mol[0].Molecule.size()*j + i].vdw_ener;
			E_E += CachedData.parameters[rec.m_mol[0].Molecule.size()*j + i].elec_ener;
		}
	return 0;
}
int group_cal::check_conti()
{
	for (size_t i{1};i < atom_index.size();++i)
	{
		if (atom_index[i] != atom_index[i - 1] + 1)
		{
			this->is_continuous = false;
			return 1;
		}
	}
	this->is_continuous = true;//to clarify
	return 0;
}
bool group_cal::check_in(size_t ind_check)
{
	if (group_set.empty())
	{
		for (auto& atom : atom_index)
			if (atom == ind_check)
				return true;
		return false;
	}
	else
	{
		if (this->group_set.find(static_cast<size_t>(ind_check)) != group_set.end())
			return true;
		else
			return false;
	}
}
int groups_cal::create_group_entropy(const Multi_MOL& mol1, const Multi_MOL& mol2, const arguments& Arguments, string name_group)
{
	group_cal group_temp{ name_group }, group_temp_rec{ "rec_" + name_group }, group_temp_lg{"lg_"+name_group};
		vector<string> atoms_of_groups{};
		if (name_group == "calpha")
			atoms_of_groups.emplace_back("CA");
		if (name_group == "backbone")
		{
			atoms_of_groups.emplace_back("CA");
			atoms_of_groups.emplace_back("C");
			atoms_of_groups.emplace_back("N");
			atoms_of_groups.emplace_back("O");
		}
		for (size_t i{};i < mol1.m_mol[0].Molecule.size();++i)
		{
			if (Arguments.topo_name_2 == "no")
			{
				for (auto& atom : atoms_of_groups)
					if (mol1.info[i].get_atom() == atom)
					{
						group_temp.atom_index.emplace_back(i + 1);
						group_temp_rec.atom_index.emplace_back(i + 1);
					}
			}
			else
			{
				for (auto& atom : atoms_of_groups)
				{
					group_temp.atom_index.emplace_back(i + 1);
					group_temp_rec.atom_index.emplace_back(i + 1);
				}
			}
		}
		fstream test_topo(Arguments.topo_name, ios::in);
		bool test{ test_topo.is_open() };
		test_topo.close();
		for (size_t i{};i < mol2.m_mol[0].Molecule.size();++i)
		{
			if (test == false)
			{
				for (auto& atom : atoms_of_groups)
					if (mol2.info[i].get_atom() == atom)
					{
						group_temp.atom_index.emplace_back(i + 1 + mol1.m_mol[0].Molecule.size());
						group_temp_lg.atom_index.emplace_back(i + 1);
					}
			}
			else
			{
				for (auto& atom : atoms_of_groups)
				{
					group_temp.atom_index.emplace_back(i + 1 + mol1.m_mol[0].Molecule.size());
					group_temp_lg.atom_index.emplace_back(i + 1);
				}
			}
		}
		group_temp.check_conti();
		group_temp_rec.check_conti();
		group_temp_lg.check_conti();
		groups.emplace_back(group_temp);
		groups.emplace_back(group_temp_rec);
		groups.emplace_back(group_temp_lg);
		group_entropy_com= groups.size() - 3;
		group_entropy_rec = groups.size() - 2;
		group_entropy_lg = groups.size() - 1;
	return 0;
}

groups_cal::groups_cal(const vector<::index>& Index) :entropy{ false }, group_1{}, group_2{}, group_entropy_com{},
group_entropy_rec{}, group_entropy_lg{} {}

int groups_cal::creat_index_file(arguments& Arguments, Multi_MOL& mol_in)
{
	for (size_t i{}; i < mol_in.Index.size() - 1; ++i)
	{
		bool check_in{ false };
		for (size_t j{}; j < groups.size(); ++j)
		{
			if (mol_in.Index[i].get_type() == groups[j].name)
			{
				check_in = true; break;
			}
		}

		if (check_in == false)
			groups.emplace_back(mol_in.Index[i].get_type());
	}

	for (size_t i{}; i < groups.size(); ++i)
	{
		for (size_t j{}; j < mol_in.Index.size() - 1; ++j)
			if (mol_in.Index[j].get_type() == groups[i].name)
			{
				for (size_t k{ mol_in.Index[j].get_res_pos_e(0) }; k < mol_in.Index[j + 1].get_res_pos_e(0); ++k)
					groups[i].atom_index.emplace_back(k);
			}
	}

	fstream out_index_file(Arguments.index_name, ios::out);
	for (size_t i{}; i < groups.size(); ++i)
	{
		out_index_file << "[ " << groups[i].name << " ]\n";
		for (size_t j{}; j < groups[i].atom_index.size(); ++j)
		{
			out_index_file << setw(8) << groups[i].atom_index[j];
			if (to_string(groups[i].atom_index[j]).length() > 8) out_index_file << " ";
			if (j != 0 && j % 15 == 0) out_index_file << '\n';
		}
	}

	return 0;
}

int groups_cal::add_sub_group_1(arguments& Arguments)
{
	for (size_t i{}; i < Arguments.sub_group_1.size(); ++i)
	{
		sub_group_1.emplace_back(Arguments.sub_group_1[i], 0);
	}
	return 0;
}

int groups_cal::add_sub_group_2(arguments& Arguments)
{
	for (size_t i{}; i < Arguments.sub_group_2.size(); ++i)
	{
		sub_group_2.emplace_back(Arguments.sub_group_2[i], 0);
	}
	return 0;
}

int groups_cal::add_sub_groups(arguments& Arguments)
{
	if (Arguments.flag_sub_groups == false)return 1;
	add_sub_group_1(Arguments);
	add_sub_group_2(Arguments);
	check_sub_group();
	parse_subgroups();
	return 0;
}

int groups_cal::check_sub_group()
{
	for (size_t i{};i<sub_group_1.size();++i)
		if (sub_group_1[i].first < 0 || sub_group_1[i].first > groups.size())
		{
			cout << "Error occur! Sub group number " << sub_group_1[i].first << " is not contained in index file! ("<< groups.size() <<")\n";
			sub_group_1.erase(sub_group_1.begin() + i);
		}
	for (size_t i{};i<sub_group_2.size();++i)
		if (sub_group_2[i].first<0 || sub_group_2[i].first > groups.size())
		{
			cout << "Error occur! Sub group number " << sub_group_2[i].first << " is not contained in index file! (" << groups.size() << ")\n";
			sub_group_2.erase(sub_group_2.begin() + i);
		}
	return 0;
}

int helper_parse_subgroups(vector<pair<int, size_t> >& subgroups_index, groups_cal& groups, size_t group)
{
	for (size_t i{}; i < subgroups_index.size(); ++i)
	{
		if (groups.atom_index(subgroups_index[i].first, 0) - groups.atom_index(group, 0) >= 0)
		{
			subgroups_index[i].second = groups.atom_index(subgroups_index[i].first, 0) - groups.atom_index(group, 0);
		}
		else
		{
			cout << "Error occur! Sub group " << groups.group_name(subgroups_index[i].first) << " is not contained in group "
				<< groups.group_name(group) << "!\n";
			subgroups_index.erase(subgroups_index.begin() + i);
		}
	}
	return 0;
}

int groups_cal::parse_subgroups()
{
	helper_parse_subgroups(sub_group_1, *this, this->group_1_index());
	helper_parse_subgroups(sub_group_2, *this, this->group_2_index());
	return 0;
}

size_t groups_cal::add_group_pathway(arguments& Arguments)
{
	if (Arguments.flag_pathway == true)
	{
		group_pathway = Arguments.group_pathway;
		return 0;
	}
	else
		return 0;
}