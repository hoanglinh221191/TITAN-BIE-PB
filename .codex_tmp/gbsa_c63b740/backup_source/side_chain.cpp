#include"distance_v2.h"
#include"arguments.h"
#include"info.h"
#include"folder.h"
#include<fstream>
#include<vector>
#include<iomanip>
#include"stat_hb_and_sc_contact.h"
#include"Multi_MOL.h"
#include<unordered_map>
#include<iostream>

int traj_sc::add_sc(string _resname, size_t _resnum, size_t _sc_num, size_t _chainid, char _chain)
{
	res_sc sc_temp(_resname, _resnum, _sc_num, _chainid, _chain);
	if(this->map_residues_sc.try_emplace(sc_temp, residues_sc.size()).second==true)
	{
		residues_sc.emplace_back(_resname, _resnum, _sc_num, _chainid, _chain);
		
	}
	return 0;
}

int cal_average_sidechain_contact_time(vector<double>& time, vector<pair<double, double> >& sc_average_sd, vector<vector<size_t> >& out_sc_time, arguments& Arguments)
{
	if (out_sc_time.empty())return 1;

	for (size_t i{}; i < time.size(); ++i)
	{
		double average{}, sd{};
		size_t size{1};

		for (size_t j{}; j < out_sc_time.size(); ++j)
		{
			if (i < out_sc_time[j].size())
			{
				average += out_sc_time[j][i];
				++size;
			}
			else continue;
		}
		average /= size;

		for (size_t j{}; j < out_sc_time.size(); ++j)
		{
			if (i < out_sc_time[j].size())
				sd += square(out_sc_time[j][i] - average);
			else continue;
		}
		sd = sqrt(sd / size);
		if (Arguments.flag_use_CI == true)sd = 1.96*sd / sqrt(size);

		sc_average_sd.emplace_back(average, sd);
	}
	return 0;
}

int cal_sidechain_contact_time(arguments& Arguments, vector<double>& time, vector<vector<size_t> >& out_sc_time, info& info_output,
	folders& paths)
{
	vector<pair<double, double> > sc_average_sd;
	cal_average_sidechain_contact_time(time, sc_average_sd, out_sc_time, Arguments);

	fstream outfile_sc_time(info_output.getSCtimeout(), ios::out);
	string error_type{ Arguments.flag_use_CI == true ? "95%CI" : "Standard deviation" };
	outfile_sc_time << "#Time (" << Arguments.tu << ")   Average SC num    " << error_type << '\n';
	for (size_t i{}; i < time.size(); ++i)
	{
		outfile_sc_time << setw(8) << left << time[i] << right << setw(12) << sc_average_sd[i].first << setw(12) << sc_average_sd[i].second << '\n';
	}

	for (size_t i{}; i < paths.get_size(); ++i)
	{
		fstream out(paths.get_path(i) + "/out_sc_num.dat", ios::out);
		for (size_t j{}; j < out_sc_time[i].size(); ++j)
		{
			if (time.size() <= j)
				break;
			out << setw(8) << left << time[j] << right << setw(12) << out_sc_time[i][j] << '\n';
		}
	}
	return 0;
}

int add_num_sc(size_t scnum, int filenum, vector<vector<size_t> > & sc_num_time)
{
	if (filenum >= sc_num_time.size())
	{
		sc_num_time.emplace_back(vector<size_t>(1, scnum));
	}
	else
		sc_num_time[filenum].emplace_back(scnum);
	return 0;
}


void print_sidechain_contact(size_t num_threads, vector<d_min>& dis_min, arguments& Arguments, fstream& out_sc,
	vector<vector<size_t> >& out_sc_time, traj_sc& sc_rec_1traj, traj_sc& sc_lg_1traj, int filenum)
{
	size_t total_sc{};
	for (size_t i{}; i < num_threads; ++i)
	{
		total_sc += dis_min[i].base.size();
		for (size_t j{}; j < dis_min[i].base.size(); ++j)
		{
			if (Arguments.flag_topo_1 == true && Arguments.flag_topo_2 == false)
			{
				if(Arguments.flag_print_sidechain==true)
					out_sc << setw(3) << dis_min[i].base[j].atom_num << " " << setw(6) << dis_min[i].base[j].chain_rec << dis_min[i].base[j].resnr_rec << '\n';
				
				sc_lg_1traj.add_sc(dis_min[i].base[j].res_name_ab, dis_min[i].base[j].atom_num, 1,
					dis_min[i].base[j].chainid_ab, dis_min[i].base[j].chain_ab[0]);
				sc_rec_1traj.add_sc(dis_min[i].base[j].res_name_rec, dis_min[i].base[j].resnr_rec, 1,
					dis_min[i].base[j].chainid_rec, dis_min[i].base[j].chain_rec[0]);
			}
			else if (Arguments.flag_topo_1 == false && Arguments.flag_topo_2 == true)
			{
				if (Arguments.flag_print_sidechain == true)
					out_sc << setw(3) << dis_min[i].base[j].chain_ab << setw(3) << left << dis_min[i].base[j].resnr_ab << right << "       " << dis_min[i].base[j].atom_num_rec << '\n';
				
				sc_rec_1traj.add_sc(dis_min[i].base[j].res_name_rec, dis_min[i].base[j].atom_num_rec, 1,
					dis_min[i].base[j].chainid_rec, dis_min[i].base[j].chain_rec[0]);
				sc_lg_1traj.add_sc(dis_min[i].base[j].res_name_ab, dis_min[i].base[j].resnr_ab, 1,
					dis_min[i].base[j].chainid_ab, dis_min[i].base[j].chain_ab[0]);
			}
			else if (Arguments.flag_topo_1 == true && Arguments.flag_topo_2 == true && Arguments.flag_use_common_topo==false)
			{
				if (Arguments.flag_print_sidechain == true)
					out_sc << setw(3) << dis_min[i].base[j].atom_num << " " << setw(6) << dis_min[i].base[j].atom_num_rec << '\n';
				
				sc_rec_1traj.add_sc(dis_min[i].base[j].res_name_rec, dis_min[i].base[j].atom_num_rec, 1,
					dis_min[i].base[j].chainid_rec, dis_min[i].base[j].chain_rec[0]);
				sc_lg_1traj.add_sc(dis_min[i].base[j].res_name_ab, dis_min[i].base[j].atom_num, 1,
					dis_min[i].base[j].chainid_ab, dis_min[i].base[j].chain_ab[0]);
			}
			else
			{
				if (Arguments.flag_print_sidechain == true)
					out_sc << setw(3) << dis_min[i].base[j].chain_ab << setw(3) << left << dis_min[i].base[j].resnr_ab << right << "       " << dis_min[i].base[j].chain_rec << dis_min[i].base[j].resnr_rec << setw(12) << dis_min[i].base[j].d << '\n';

				sc_rec_1traj.add_sc(dis_min[i].base[j].res_name_rec, dis_min[i].base[j].resnr_rec, 1,
					dis_min[i].base[j].chainid_rec, dis_min[i].base[j].chain_rec[0]);
				sc_lg_1traj.add_sc(dis_min[i].base[j].res_name_ab, dis_min[i].base[j].resnr_ab, 1,
					dis_min[i].base[j].chainid_ab, dis_min[i].base[j].chain_ab[0]);
			}
		}
	}
	add_num_sc(total_sc, filenum, out_sc_time);
}

bool traj_sc::check_in(const res_sc& res_contact, size_t& pos)
{
	for (size_t i{}; i<residues_sc.size(); ++i)
	{
		if (residues_sc[i].equal_to(res_contact))
		{
			pos = i;
			return true;
		}
	}
	return false;
}

int traj_sc::filt_sc(vector<traj_sc>& sc_reslist)
{
	for (size_t i{}; i < sc_reslist.size(); ++i)
	{
		for (size_t j{}; j < sc_reslist[i].residues_sc.size(); ++j)
		{
			auto emplace_result = this->map_residues_sc.try_emplace(sc_reslist[i].residues_sc[j], residues_sc.size());
			if (emplace_result.second==false)
				this->residues_sc[emplace_result.first->second].sc_num += sc_reslist[i].residues_sc[j].sc_num;
			else
			{
				this->residues_sc.emplace_back(sc_reslist[i].residues_sc[j]);
			}
		}
		this->snap_num += sc_reslist[i].snap_num;
	}
	return 0;
}

int stat_sidechain_contact(vector<traj_sc>& sc_reslist, string output_filename)
{
	traj_sc stat_out_sc;
	stat_out_sc.filt_sc(sc_reslist);

	fstream output(output_filename, ios::out);
	output << "#Total snap for all trajectories: " << stat_out_sc.snap_num << " ; traj num: " << sc_reslist.size() << '\n';

	sc_reslist.clear();

	output << "#Chain  resnum   resname      average SC contacts\n";

	for (auto& sc : stat_out_sc.residues_sc)
	{
		output << sc.chain << setw(5) << sc.resnum << setw(5) << sc.resname << setw(15) << sc.sc_num * 100. / stat_out_sc.snap_num << '\n';
	}
	return 0;
}

int calc_contact_map(Multi_MOL& mol1, Multi_MOL& mol2, size_t num_threads, vector<d_min>& dis_min, arguments& Arguments, size_t num_snap
, vector<vector<map_pos>>& matrix, unordered_map<tuple<size_t, size_t, size_t, size_t>, pair<size_t, size_t>,tuple_hash>& pos_matrix)
{
	if (matrix.empty())
	{
		size_t num_res_1{}, num_res_2{};

		for (size_t i{}; i < mol1.Index.size() - 1; ++i)
			num_res_1 += mol1.Index[i].get_size_v();
		//std::cout << "group " << Arguments.group_1 << " res num:" << num_res_1 << '\n';
		for (size_t i{}; i < mol2.Index.size() - 1; ++i)
			num_res_2 += mol2.Index[i].get_size_v();
		//std::cout << "group " << Arguments.group_2 << " res num:" << num_res_2 << '\n';

		if (num_res_1 == 0 || num_res_2 == 0)
			return 1;

		matrix = vector<vector<map_pos>>(num_res_1);
		//std::cout << "matrix size " << matrix.size() << '\n';
		for (size_t i{}; i < matrix.size(); ++i)
		{
			matrix[i] = vector<map_pos>(num_res_2);
			//std::cout << matrix[i].size() << " ";
		}
		//std::cout << '\n';
		size_t pos1{};
		for (size_t i{}; i < mol1.Index.size() - 1; ++i)
			for (size_t j{}; j < mol1.Index[i].get_size_v(); ++j)
			{
				size_t pos2{};
				for (size_t k{}; k < mol2.Index.size() - 1; ++k)
					for (size_t l{}; l < mol2.Index[k].get_size_v(); ++l)
					{
						if (pos1 >= matrix.size() || pos2 >= matrix[pos1].size())
						{
							//std::cout << "overflow" << pos1<<"  "<<pos2<<"\n";
							break;
						}
						matrix[pos1][pos2].chain_1 = mol1.Index[i].get_chain();
						matrix[pos1][pos2].chain_2 = mol2.Index[k].get_chain();
						matrix[pos1][pos2].chain_1id = mol1.Index[i].get_chainid();
						matrix[pos1][pos2].chain_2id = mol2.Index[k].get_chainid();
						matrix[pos1][pos2].resnum_1 = mol1.Index[i].get_res_num_e(j);
						matrix[pos1][pos2].resnum_2 = mol2.Index[k].get_res_num_e(l);
						size_t res2{ mol2.Index[k].get_res_num_e(l) };
						if (mol2.get_mol_type() == "ligand")
						{
							matrix[pos1][pos2].resnum_2 = mol2.Index[k].get_res_pos_e(l);
							res2= mol2.Index[k].get_res_pos_e(l);
						}
						//cout << "mol1 chainid: " << mol1.Index[i].get_chainid() << " " << mol2.Index[k].get_chainid() << " "
						//	<< mol1.Index[i].get_res_num_e(j) << " " << mol2.Index[k].get_res_pos_e(l) << '\n';

						pos_matrix.try_emplace(tuple<size_t, size_t, size_t, size_t>( mol1.Index[i].get_chainid(), mol2.Index[k].get_chainid(), mol1.Index[i].get_res_num_e(j), res2), pair<size_t, size_t>(pos1, pos2));
						++pos2;
					}
				++pos1;
			}

	}

	for (size_t i{};i<dis_min.size();++i)
		for (size_t j{}; j < dis_min[i].base.size(); ++j)
		{
			tuple<size_t, size_t, size_t, size_t> matrix_pos{dis_min[i].base[j].chainid_rec, dis_min[i].base[j].chainid_ab ,dis_min[i].base[j].resnr_rec,dis_min[i].base[j].resnr_ab };

			//cout << dis_min[i].base[j].chainid_rec << dis_min[i].base[j].chainid_ab << dis_min[i].base[j].resnr_rec << dis_min[i].base[j].atom_num << '\n';

			if (mol2.get_mol_type() == "ligand")
				matrix_pos = tuple<size_t, size_t, size_t, size_t>{ dis_min[i].base[j].chainid_rec, dis_min[i].base[j].chainid_ab ,dis_min[i].base[j].resnr_rec,dis_min[i].base[j].atom_num };

			++matrix[pos_matrix[matrix_pos].first][pos_matrix[matrix_pos].second].cont_num;
		}
	return 0;
}

int print_cont_map(vector<vector<map_pos>>& matrix, fstream& out_contact_map, size_t num_snap)
{
	out_contact_map.precision(4);

	//cout << "matcontact size: " << matrix.size() << " \n";

	for (size_t i{}; i < matrix.size(); ++i)
	{
		for (size_t j{}; j < matrix[i].size(); ++j)
			out_contact_map << setw(9) << left << fixed << matrix[i][j].cont_num * 100. / num_snap;
		out_contact_map << '\n';
	}
	return 0;
}