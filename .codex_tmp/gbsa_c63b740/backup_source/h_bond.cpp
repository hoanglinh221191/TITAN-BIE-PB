#include"h_bond.h"
#include<cmath>
#include<thread>
//#include<mutex>
#include<fstream>
#include<iostream>
#include<iomanip>
#include"MOL.h"
#include"miscellaneous.h"
#include"functions.h"
using namespace std;
  double PI = 3.1415926535897;
//std::mutex emplace_mutex_1;
//std::mutex emplace_mutex_2;
double cal_angle_3_atoms(const molecule& mol1, const molecule& mol2, const molecule& mol3);
int print_hbond_per_snap(MOL& pro_in, MOL& lig, vector<mol_info>& info_rec, vector<mol_info>& info_lg, const arguments& Arguments);
inline bool check_sidechain_hbond( const arguments& Arguments, string atom)
{
	return (!Arguments.flag_sidechain_hbond || (Arguments.flag_sidechain_hbond && atom != "N" && atom != "O"));
}
int h_group::add_groups(const MOL& mol_in, const MOL& ref, size_t i, const size_t remain, char _name, double bond_length, vector<mol_info>& info)
{
	name = _name;
	index = i;
	size_t next{ i + remain > mol_in.Molecule.size() ? mol_in.Molecule.size() : i + remain };
	for (size_t j{ i };j < next;++j)
	{
		if (info[j].get_resnr() == info[i].get_resnr())
		{
			if ((info[j].get_atom_type()[0] == 'H' || (info[j].get_atom().length() >= 4 && info[j].get_atom()[1] == 'H') || info[j].get_atom()[0] == 'H')
				&& mol_in.Molecule[i].dis_2_atoms(mol_in.Molecule[j]) < bond_length)
				h_index.emplace_back(j);
		}
		else
			break;
	}
	return 0;
}
int h_groups::find_groups(const MOL& mol_in, const MOL& ref, vector<mol_info>& info, const arguments& Arguments)
{
	  size_t remain{ 120 };
	for (size_t i{};i < mol_in.Molecule.size();++i)
	{
		h_group group_temp;
		if ((info[i].get_atom_type()[0] == 'O' || info[i].get_atom()[0] == 'O') && check_sidechain_hbond(Arguments,info[i].get_atom()))
		{
			group_temp.add_groups(mol_in, ref, i, remain, 'O', 1.1, info);
			groups.emplace_back(group_temp);
		}
		else if ((info[i].get_atom_type()[0]=='N' || info[i].get_atom()[0] == 'N') && check_sidechain_hbond(Arguments, info[i].get_atom()))
		{
			group_temp.add_groups(mol_in, ref, i, remain, 'N', 1.5,info);
			groups.emplace_back(group_temp);
		}
		else if (info[i].get_atom_type()[0] == 'F' || info[i].get_atom()[0] == 'F')
		{
			group_temp.add_groups(mol_in, ref, i, remain, 'F', 1.1,info);
			groups.emplace_back(group_temp);
		}
		else continue;
	}
	return 0;
}
void cal_h_bonds(MOL& pro_in, MOL& lig, size_t begin, size_t end, const arguments& Arguments, h_groups& groups_1, h_groups& groups_2
, h_bonds_wrapper& dat_pro_in, h_bonds_wrapper& dat_lig)
{
	for (size_t i{ begin };i < end;++i)
	{
		for (size_t j{};j < groups_2.groups.size();++j)
		{
			double d{ pro_in.Molecule[groups_1.groups[i].index].dis_2_atoms(lig.Molecule[groups_2.groups[j].index]) };
			
			if ((!groups_1.groups[i].h_index.empty() || !groups_2.groups[j].h_index.empty()) 
				&&	d <= Arguments.dis_h && d> 0.01)
			{
				for (size_t k{};k < groups_1.groups[i].h_index.size();++k)
				{
					double angle{ cal_angle_3_atoms(pro_in.Molecule[groups_1.groups[i].index], pro_in.Molecule[groups_1.groups[i].h_index[k]], lig.Molecule[groups_2.groups[j].index]) };
					if (angle >= Arguments.angle_h)
					{
						//std::lock_guard<std::mutex> lock(emplace_mutex_1);
						dat_lig.dat.emplace_back(groups_1.groups[i].index, groups_2.groups[j].index, groups_1.groups[i].h_index[k], angle, d);
					}
				}
				for (size_t k{};k < groups_2.groups[j].h_index.size();++k)
				{
					double angle{ cal_angle_3_atoms(lig.Molecule[groups_2.groups[j].index], lig.Molecule[groups_2.groups[j].h_index[k]], pro_in.Molecule[groups_1.groups[i].index]) };
					if (angle >= Arguments.angle_h)
					{
						//std::lock_guard<std::mutex> lock(emplace_mutex_2);
						dat_pro_in.dat.emplace_back(groups_2.groups[j].index, groups_1.groups[i].index, groups_2.groups[j].h_index[k], angle, d);
					}
				}
			}
		}
	}
}
int find_h_bonds(MOL& pro_in, MOL& pro_in_ref, MOL& lig, MOL& lig_ref, h_groups& groups_1, h_groups& groups_2, const arguments& Arguments, size_t num_snap,
	vector<mol_info>& info_rec, vector<mol_info>& info_lg, vector<vector<size_t> >& Hbonds_vector)
{
	if (num_snap == 0)
	{
		if (Arguments.nt >= 2)
		{
			thread t1(&h_groups::find_groups, &groups_1, ref(pro_in), ref(pro_in_ref), ref(info_rec), ref(Arguments));
			thread t2(&h_groups::find_groups, &groups_2, ref(lig), ref(lig_ref), ref(info_lg), ref(Arguments));
			t1.join();
			t2.join();
		}
		else
		{
			groups_1.find_groups(pro_in, pro_in_ref, info_rec, Arguments);
			groups_2.find_groups(lig, lig_ref, info_lg, Arguments);
		}
	}

	vector<h_bonds_wrapper> h_bonds_pro{ Arguments.nt }, h_bonds_lig{ Arguments.nt };

	vector<thread> threads_h(Arguments.nt);
	for (size_t i{};i < Arguments.nt;++i)
	{
		threads_h[i] = thread(cal_h_bonds, ref(pro_in), ref(lig), i*groups_1.groups.size() / Arguments.nt, (i + 1)*groups_1.groups.size() / Arguments.nt
			, ref(Arguments), ref(groups_1), ref(groups_2), ref(h_bonds_pro[i]), ref(h_bonds_lig[i]));
	}
	for (auto& t : threads_h) { t.join(); }

	for (size_t i{}; i < Arguments.nt; ++i)
	{
		for (size_t j{}; j < h_bonds_pro[i].dat.size(); ++j)
			pro_in.list_h_bonds.dat.emplace_back(h_bonds_pro[i].dat[j]);
		for (size_t j{}; j < h_bonds_lig[i].dat.size(); ++j)
			lig.list_h_bonds.dat.emplace_back(h_bonds_lig[i].dat[j]);
	}

	//cout << "pos proc..";
	lig.post_proc_hbonds(lig, info_lg, info_rec);
	//cout << " lig done ";
	pro_in.post_proc_hbonds(pro_in, info_rec, info_lg);
	//cout << "beg.. "<<num_snap;
	print_hbond_per_snap(pro_in, lig, info_rec, info_lg, Arguments);
	//cout << " printed  ";
	if (num_snap >= Hbonds_vector.size())
	{
		Hbonds_vector.emplace_back(vector<size_t>(1,pro_in.list_h_bonds.dat.size() + lig.list_h_bonds.dat.size()));
	}
	else
	{
		Hbonds_vector[num_snap].emplace_back(pro_in.list_h_bonds.dat.size() + lig.list_h_bonds.dat.size());
	}
	//cout << "exported!\n";
	return 0;
}
inline double cal_angle_3_atoms(const molecule& mol1, const molecule& mol2, const molecule& mol3)
{
	double vector_1[3]{ mol1.get_x() - mol2.get_x(), mol1.get_y() - mol2.get_y(), mol1.get_z() - mol2.get_z() },
		vector_2[3]{ mol3.get_x() - mol2.get_x(),mol3.get_y() - mol2.get_y(),mol3.get_z() - mol2.get_z() };
	return fabs(acos((vector_1[0] * vector_2[0] + vector_1[1] * vector_2[1] + vector_1[2] * vector_2[2]) / (mol1.dis_2_atoms(mol2)*mol3.dis_2_atoms(mol2)))*180. / PI);
}
int MOL::post_proc_hbonds(const MOL& mol, vector<mol_info>& info, vector<mol_info>& counter_info)
{
	vector<pair<char, int> > max_hbonds;
	max_hbonds.emplace_back('O', 2);
	max_hbonds.emplace_back('N', 3);
	max_hbonds.emplace_back('H', 1);
	max_hbonds.emplace_back('F', 1);
	for (size_t i{};i < max_hbonds.size();++i)
		imp_post_proc(max_hbonds[i], mol, info, counter_info);
	return 0;
}
int MOL::imp_post_proc(const pair<char,size_t>& atom_dat, const MOL& mol, vector<mol_info>& info, vector<mol_info>& counter_info)
{
	vector<hbonds_pos> pos_summary;
	for (size_t i{}; i < list_h_bonds.dat.size();++i)
	{
		//if (info.size() <= list_h_bonds.dat[i].get_index_acceptor()) { cout << info.size() << " smaller than acceptor index " << list_h_bonds.dat[i].get_index_acceptor() << '\n'; continue; }
		string name{ info[list_h_bonds.dat[i].get_index_acceptor()].get_atom_type() };
		int count_appear{};
		size_t index{ list_h_bonds.dat[i].get_index_acceptor() };
		if (atom_dat.first == 'H')
		{
			//if (counter_info.size() <= list_h_bonds.dat[i].get_index_h()) { cout << counter_info.size() << " smaller than h index " << list_h_bonds.dat[i].get_index_h() << '\n'; continue; }
			name = counter_info[list_h_bonds.dat[i].get_index_h()].get_atom_type();
			index = list_h_bonds.dat[i].get_index_h();
		}
		if (name.find(atom_dat.first) != string::npos)
		{
			for (size_t j{}; j < pos_summary.size();++j)
			{
				if (index == pos_summary[j].index_mol)
				{
					pos_summary[j].index_rep.emplace_back(i);
					++count_appear;
				}
				else continue;
			}
			if (count_appear > 0) count_appear = 0;
			else pos_summary.emplace_back(index,i);
		}
	}
	for (size_t i{};i < pos_summary.size();++i)
	{
		if (pos_summary[i].index_rep.size() > atom_dat.second)
		{
			while (pos_summary[i].index_rep.size() != atom_dat.second)
			{
				double max_dis{ 0. };
				size_t pos_max_dis{}, pos_summary_del{};
				for (size_t j{};j < pos_summary[i].index_rep.size();++j)
				{
					if (max_dis < list_h_bonds.dat[pos_summary[i].index_rep[j]].get_dis())
					{
						max_dis = list_h_bonds.dat[pos_summary[i].index_rep[j]].get_dis();
						pos_max_dis = pos_summary[i].index_rep[j];
						pos_summary_del = j;
					}
				}
				//if (pos_summary_del >= pos_summary[i].index_rep.size())cout << " !!!pos_summary ";
				pos_summary[i].index_rep.erase(pos_summary[i].index_rep.begin() + pos_summary_del);
				//if (pos_max_dis >= list_h_bonds.dat.size())cout << " !!!list h bond.dat ";
				list_h_bonds.dat.erase(list_h_bonds.dat.begin() + pos_max_dis);
				
				for (size_t k{}; k < pos_summary.size(); ++k)
					for (size_t l{}; l < pos_summary[k].index_rep.size(); ++l)
						if (pos_summary[k].index_rep[l] > pos_max_dis) --pos_summary[k].index_rep[l];
			}
		}
	}
	return 0;
}
void h_groups::print_group(int m)
{
	fstream groups_out("out_hbonds_groups_"+to_string(m)+".dat", ios::out);
	for (size_t i{};i < groups.size();++i)
	{
		groups_out << groups[i].name << " " << groups[i].index << '\n';
		for (size_t j{};j < groups[i].h_index.size();++j)
			groups_out << groups[i].h_index[j] << "   ";
		groups_out << '\n';
	}
}

int stat_hbonds(vector<vector<size_t> >& Hbonds_vector, const info& infos, arguments& Arguments, vector<double>& time_of_snaps, vector<string>& list_name)
{
	if (Arguments.flag_h_bond == 0)return 1;
	fstream statHB(infos.get_stathb_name(), ios::out);
	statHB << "#Time (" << Arguments.tu << ")   Number of H-bonds   SD\n";
	for (size_t i{}; i < Hbonds_vector.size(); ++i)
	{
		double mean{}, sd{};
		for (size_t j{}; j < Hbonds_vector[i].size(); ++j)
		{
			mean += Hbonds_vector[i][j];
		}
		mean /= Hbonds_vector[i].size();
		for (size_t j{}; j < Hbonds_vector[i].size(); ++j)
			sd += square(Hbonds_vector[i][j] - mean);
		sd = sqrt(sd / Hbonds_vector[i].size());
		statHB << setw(8) << time_of_snaps[i] << setw(12) << mean << setw(12) << sd << '\n';
	}
	statHB.close();
	fstream stat_aveHB(infos.get_stat_ave_hb_name(), ios::out);
	for (size_t i{}; i < list_name.size(); ++i)
	{
		double mean{}, sd{};
		size_t Size{Hbonds_vector.size()};
		for (size_t j{}; j < Hbonds_vector.size(); ++j)
			if (Hbonds_vector[j].size() > i)
				mean += Hbonds_vector[j][i];
			else { Size = j; break; }
		mean /= Size;
		for (size_t j{}; j < Size; ++j)
				sd += square(Hbonds_vector[j][i] - mean);
		stat_aveHB << setw(15) << list_name[i] << setw(12) << mean << setw(12) << sqrt(sd / Size) << '\n';

	}
	stat_aveHB.close();
	return 0;
}

int print_hbond_per_snap(MOL& pro_in, MOL& lig, vector<mol_info>& info_rec, vector<mol_info>& info_lg, const arguments& Arguments)
{
	fstream out_h_bonds("out_hydrogen_bonds.dat", ios::app);
	for (size_t i{}; i < pro_in.list_h_bonds.dat.size(); ++i)
	{
		out_h_bonds << setw(4) << left << info_rec[pro_in.list_h_bonds.dat[i].get_index_acceptor()].get_atom() << setw(3) << right << info_rec[pro_in.list_h_bonds.dat[i].get_index_acceptor()].get_chain() << setw(5) << right << info_rec[pro_in.list_h_bonds.dat[i].get_index_acceptor()].get_resnr()
			<< setw(6) << right << info_lg[pro_in.list_h_bonds.dat[i].get_index_h()].get_atom() << setw(5) << right << info_lg[pro_in.list_h_bonds.dat[i].get_index_h()].get_resnr() << '\n';
	}
	for (size_t i{}; i < lig.list_h_bonds.dat.size(); ++i)
	{
		if (lig.list_h_bonds.dat[i].get_index_h() >= info_rec.size()) cout << "info_rec.size smaller than index h\n";
		if (lig.list_h_bonds.dat[i].get_index_acceptor() > info_lg.size())cout << "info_lg size smaller than index acc\n";

		out_h_bonds << setw(4) << left << info_lg[lig.list_h_bonds.dat[i].get_index_acceptor()].get_atom() << setw(3) << right << info_lg[lig.list_h_bonds.dat[i].get_index_acceptor()].get_chain() << setw(5) << right << info_lg[lig.list_h_bonds.dat[i].get_index_acceptor()].get_resnr()
			<< setw(6) << right << info_rec[lig.list_h_bonds.dat[i].get_index_h()].get_atom() << setw(5) << right << info_rec[lig.list_h_bonds.dat[i].get_index_h()].get_resnr() << '\n';
	}
	out_h_bonds.close();
	fstream out_h_bonds_rec("out_hydrogen_bonds_receptor.dat", ios::app);
	for (size_t i{}; i < pro_in.list_h_bonds.dat.size(); ++i)
	{
		out_h_bonds_rec << setw(4) << left << info_rec[pro_in.list_h_bonds.dat[i].get_index_acceptor()].get_atom() << setw(5) << right << info_rec[pro_in.list_h_bonds.dat[i].get_index_acceptor()].get_resnr()
			<< setw(6) << right << info_lg[pro_in.list_h_bonds.dat[i].get_index_h()].get_atom() << setw(5) << right << info_lg[pro_in.list_h_bonds.dat[i].get_index_h()].get_resnr() << '\n';
	}
	out_h_bonds_rec.close();
	fstream out_h_bonds_lig("out_hydrogen_bonds_ligand.dat", ios::app);
	for (size_t i{}; i < lig.list_h_bonds.dat.size(); ++i)
	{
		out_h_bonds_lig << setw(4) << left << info_lg[lig.list_h_bonds.dat[i].get_index_acceptor()].get_atom() << setw(5) << right << info_lg[lig.list_h_bonds.dat[i].get_index_acceptor()].get_resnr()
			<< setw(6) << right << info_rec[lig.list_h_bonds.dat[i].get_index_h()].get_atom() << setw(5) << right << info_rec[lig.list_h_bonds.dat[i].get_index_h()].get_resnr() << '\n';
	}
	out_h_bonds_lig.close();
	return 0;
}