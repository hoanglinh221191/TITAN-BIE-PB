#include"stat_ener_res_res.h"
#include<vector>
#include<iomanip>
#include<fstream>
#include<iostream>
#include <functional>
#include"Multi_MOL.h"
#include"info.h"
using namespace std;

/*bool info_ener_res::check_in(pair<info_ener, tuple<size_t, double, double>>& rhs)
{
	for (auto& ele : info_energy_respair)
	{
		if (ele.first.resnum_1 == rhs.first.resnum_1 && ele.first.resnum_2 == rhs.first.resnum_2
			&& ele.first.chain_1 == rhs.first.chain_1 && ele.first.chain_2 == rhs.first.chain_2
			&& ele.first.chain_1id == rhs.first.chain_1id && ele.first.chain_2id == rhs.first.chain_2id)
			return true;
		else
			continue;
	}
	return false;
}

bool info_ener_res::check_in(info_ener& rhs)
{
	for (auto& ele : info_energy_respair)
	{
		if (ele.first.resnum_1 == rhs.resnum_1 && ele.first.resnum_2 == rhs.resnum_2
			&& ele.first.chain_1 == rhs.chain_1 && ele.first.chain_2 == rhs.chain_2
			&& ele.first.chain_1id == rhs.chain_1id && ele.first.chain_2id == rhs.chain_2id)
			return true;
		else
			continue;
	}
	return false;
}*/

size_t info_ener_res::find(tuple<size_t, size_t, size_t, size_t>&& postemp, size_t begin, size_t end, size_t& result)
{
	if (position_map.find(postemp) != position_map.end())
	{
		result = position_map[postemp];
		return position_map[postemp];
	}
	else
	{
		result = info_energy_respair.size();
		return result;
	}
}

size_t info_ener_res::find(const info_ener& rhs, size_t begin, size_t end, size_t& result)
{
	tuple<size_t, size_t, size_t, size_t> postemp(rhs.chain_1id, rhs.chain_2id, rhs.resnum_1, rhs.resnum_2);
	if (position_map.find(postemp) != position_map.end())
	{
		result = position_map[postemp];
		return position_map[postemp];
	}
	else
	{
		result = info_energy_respair.size();
		return info_energy_respair.size();
	}
}

int info_ener_res::add_energy(pair<info_ener, tuple<size_t, double, double>>& rhs)
{
	size_t pos{ info_energy_respair.size() };
	tuple<size_t, size_t, size_t, size_t> temp(rhs.first.chain_1id, rhs.first.chain_2id, rhs.first.resnum_1, rhs.first.resnum_2);
	if(position_map.try_emplace(temp, pos).second==true)
	{
		info_energy_respair.emplace_back(rhs.second);
	}
	else
	{
		pos = position_map[temp];
		get<0>(info_energy_respair[pos]) += get<0>(rhs.second);
		get<1>(info_energy_respair[pos]) += get<1>(rhs.second);
		get<2>(info_energy_respair[pos]) += get<2>(rhs.second);
	}
		
	return 0;
}

int info_ener_res::add_energy(tuple<size_t, size_t, size_t, size_t>&& ref, double ee, double ev)
{
	size_t pos{ info_energy_respair.size()};
	if ( position_map.try_emplace(ref,pos).second ==false)
	{
		pos = position_map[ref];
		++get<0>(info_energy_respair[pos]);
		get<1>(info_energy_respair[pos]) += ee;
		get<2>(info_energy_respair[pos]) += ev;
	}
	else
	{
		info_energy_respair.emplace_back(tuple<size_t, double, double>(1u, ee, ev));
	}
	return 0;
}

int info_ener_res::add_energy(const tuple<size_t, size_t, size_t, size_t>& ref, double ee, double ev)
{
	size_t pos{ info_energy_respair.size() };
	if ( position_map.try_emplace(ref, pos).second == false)
	{
		pos = position_map[ref];
		++get<0>(info_energy_respair[pos]);
		get<1>(info_energy_respair[pos]) += ee;
		get<2>(info_energy_respair[pos]) += ev;
	}
	else
	{
		info_energy_respair.emplace_back(tuple<size_t, double, double>(1u, ee, ev));
	}
	return 0;
}

int info_ener_res::wrapper_impl_ener_maptype(Multi_MOL& mol1, Multi_MOL& mol2, info& output_info, const arguments& Arguments, size_t num_snap)
{
	vector<vector<double>> emap_elec, emap_vdw;
	info_ener_1mol_res mol1_ener_res, mol2_ener_res;

	this->cal_energy_maptype(mol1, mol2, emap_elec, emap_vdw, mol1_ener_res, mol2_ener_res, Arguments, num_snap);
	this->print_energy_maptype(mol1, mol2, emap_elec, emap_vdw, num_snap, output_info,Arguments);
	this->print_energy_pdbtype(mol1, mol2, mol1_ener_res, mol2_ener_res, num_snap, output_info,Arguments);

	return 0;
}

int info_ener_res::build_pos_map(Multi_MOL& mol1, Multi_MOL& mol2, vector<vector<size_t>>& pos_map, const arguments& Arguments)
{
	size_t pos1{};
	for (size_t i{}; i < mol1.Index.size() - 1; ++i)
		for (size_t j{}; j < mol1.Index[i].get_size_v(); ++j)
		{
			size_t pos2{};
			for (size_t k{}; k < mol2.Index.size() - 1; ++k)
				for (size_t l{}; l < mol2.Index[k].get_size_v(); ++l)
				{
					if (pos1 >= pos_map.size() || pos2 >= pos_map[pos1].size())
					{
						break;
					}

					;

					vector<size_t> find_result(Arguments.nt);
					find(tuple<size_t, size_t, size_t, size_t> ( mol1.Index[i].get_chainid(), mol2.Index[k].get_chainid(),
						mol1.Index[i].get_res_num_e(j), mol2.Index[k].get_res_num_e(l) ), 0, info_energy_respair.size(), ref(find_result[i]));

					size_t pos{ info_energy_respair.size() };
					for (auto& elec : find_result)
						if (elec < info_energy_respair.size())
							pos = elec;
					
					pos_map[pos1][pos2] = pos;

					++pos2;
				}
			++pos1;
		}
	return 0;
}

size_t info_ener_res::determine_pos(Multi_MOL& mol1, Multi_MOL& mol2, const arguments& Arguments,
	size_t i, size_t j, size_t l, size_t k)
{
	size_t pos{};
	find(tuple( mol1.Index[i].get_chainid(), mol2.Index[k].get_chainid(), mol1.Index[i].get_res_num_e(j), mol2.Index[k].get_res_num_e(l)), 
		0, info_energy_respair.size(), pos);
	return pos;
}

int info_ener_res::cal_energy_maptype(Multi_MOL& mol1, Multi_MOL& mol2, vector<vector<double>>& emap_elec, vector<vector<double>>& emap_vdw,
	info_ener_1mol_res& mol1_ener_res, info_ener_1mol_res& mol2_ener_res, const arguments& Arguments, size_t num_snap)
{
	if (info_energy_respair.empty())
	{
		cout << "info energy is empty\n";
		return 1;
	}
	cout << "calculating maptype....\n";
	size_t num_res_1{}, num_res_2{};

	for (size_t i{}; i < mol1.Index.size() - 1; ++i)
		num_res_1 += mol1.Index[i].get_size_v();
	std::cout << " res num:" << num_res_1 << '\n';
	for (size_t i{}; i < mol2.Index.size() - 1; ++i)
		num_res_2 += mol2.Index[i].get_size_v();
	std::cout << " res num:" << num_res_2 << '\n';

	if (num_res_1 == 0 || num_res_2 == 0)
	{
		cout << "num res of mol1: " << num_res_1 << '\n'
			<< "num res of mol2: " << num_res_2 << '\n';
		return 1;
	}


	emap_elec = vector<vector<double>>(num_res_1);
	emap_vdw = vector<vector<double>>(num_res_1);
	for (size_t i{}; i < emap_elec.size(); ++i)
	{
		emap_elec[i] = vector<double>(num_res_2,0.);
		emap_vdw[i] = vector<double>(num_res_2,0.);
	}

	cout << "res pair num: " << info_energy_respair.size() << '\n';

	//vector<thread> threads(Arguments.nt);

	size_t pos1{};
	for (size_t i{}; i < mol1.Index.size() - 1; ++i)
		for (size_t j{}; j < mol1.Index[i].get_size_v(); ++j)
		{
			size_t pos2{};
			for (size_t k{}; k < mol2.Index.size() - 1; ++k)
				for (size_t l{}; l < mol2.Index[k].get_size_v(); ++l)
				{
					if (pos1 >= emap_elec.size() || pos2 >= emap_elec[pos1].size())
					{
						//	std::cout << "overflow" << pos1<<"  "<<pos2<<" "<< emap_elec.size()<<"\n";
						break;
					}
					
					size_t pos{ determine_pos(mol1,mol2,Arguments,i,j,l,k) };

					if (pos != info_energy_respair.size())
					{
						emap_elec[pos1][pos2] += get<1>(info_energy_respair[pos]);
						emap_vdw[pos1][pos2] += get<2>(info_energy_respair[pos]);
						//cout << "\ncheck1 ";
						info_ener_1mol ener_1mol_temp(mol1.Index[i].get_chain(), mol1.Index[i].get_chainid(), mol1.Index[i].get_res_num_e(j));
						++get<0>(mol1_ener_res.info_energy_1mol_respair[ener_1mol_temp]);
						get<1>(mol1_ener_res.info_energy_1mol_respair[ener_1mol_temp]) += get<1>(info_energy_respair[pos]);
						get<2>(mol1_ener_res.info_energy_1mol_respair[ener_1mol_temp]) += get<2>(info_energy_respair[pos]);
						//cout << "check2 ";
						info_ener_1mol ener_1mol_temp2(mol2.Index[k].get_chain(), mol2.Index[k].get_chainid(), mol2.Index[k].get_res_num_e(l));
						++get<0>(mol2_ener_res.info_energy_1mol_respair[ener_1mol_temp2]);
						get<1>(mol2_ener_res.info_energy_1mol_respair[ener_1mol_temp2]) += get<1>(info_energy_respair[pos]);
						get<2>(mol2_ener_res.info_energy_1mol_respair[ener_1mol_temp2]) += get<2>(info_energy_respair[pos]);
						//cout << "check3 ";
					}
					//cout << "check ";
					++pos2;
				}
			++pos1;
		}

	cout << "build done\n";

	return 0;
}

int info_ener_res::print_energy_maptype(Multi_MOL& mol1, Multi_MOL& mol2, vector<vector<double>>& emap_elec, vector<vector<double>>& emap_vdw, 
	size_t num_snap, info& output_info, const arguments& Arguments)
{
	fstream output_elec("electrostatic" + output_info.getOuputEnerMaptypeAve(), ios::out), output_vdw("vanderWaals" + output_info.getOuputEnerMaptypeAve(), ios::out),
		output_total("total" + output_info.getOuputEnerMaptypeAve(), ios::out);

	output_elec.precision(4);
	output_vdw.precision(4);
	output_total.precision(4);
	for (size_t i{}; i < emap_elec.size(); ++i)
	{
		for (size_t j{}; j < emap_elec[i].size(); ++j)
		{
			output_elec << fixed << left << setw(9) << emap_elec[i][j] / num_snap / Arguments.eu_par<< " ";
			output_vdw << fixed << left << setw(9) << emap_vdw[i][j] / num_snap / Arguments.eu_par << " ";
			output_total << fixed << left << setw(9) << (emap_elec[i][j] + emap_vdw[i][j]) / num_snap / Arguments.eu_par << " ";
		}
		output_elec << '\n';
		output_vdw << '\n';
		output_total << '\n';
	}

	return 0;
}

int info_ener_res::print_energy_pdbtype(Multi_MOL& mol1, Multi_MOL& mol2, info_ener_1mol_res& mol1_ener_res, info_ener_1mol_res& mol2_ener_res,
	size_t num_snap, info& output_info, const arguments& Arguments)
{
	fstream output("total"+output_info.getOutputEnerPdbtypeAve(), ios::out);
	output.precision(3);

	size_t num{};
	for (size_t l{}; l < mol1.m_mol[0].Molecule.size(); ++l)
	{
		l > 99998 ? num = 99999 : num = l;
		size_t factor{ get<0>(mol1_ener_res.info_energy_1mol_respair[info_ener_1mol(mol1.info[l].get_chain(), mol1.info[l].get_chainnum(), mol1.info[l].get_resnr())]) };
		if (factor == 0 || Arguments.mol_1_type.find("protein")!=string::npos)
			factor = 1;
		string chain{ " " };
		if (!isdigit(mol1.info[l].get_chain()[0]))chain = mol1.info[l].get_chain();

		output << setw(4) << "ATOM" << setw(7) << num + 1 << setw(5) << mol1.info[l].get_atom() << " " << setw(4) << left << mol1.info[l].get_original_res_name() << right << setw(1) << chain
			<< setw(4) << mol1.info[l].get_resnr() << setw(12) << fixed << mol1.m_mol[0].Molecule[l].get_x() << setw(8) << fixed << mol1.m_mol[0].Molecule[l].get_y() << setw(8)
			<< mol1.m_mol[0].Molecule[l].get_z()  << setw(6) << "1.00" 
			<< setw(7) << (get<1>(mol1_ener_res.info_energy_1mol_respair[info_ener_1mol(mol1.info[l].get_chain(),mol1.info[l].get_chainnum(),mol1.info[l].get_resnr())]) / num_snap + get<2>(mol1_ener_res.info_energy_1mol_respair[info_ener_1mol(mol1.info[l].get_chain(), mol1.info[l].get_chainnum(), mol1.info[l].get_resnr())]) / num_snap)/Arguments.eu_par
			/ factor
			<< '\n';
	}

	for (size_t l{}; l < mol2.m_mol[0].Molecule.size(); ++l)
	{
		size_t factor{ 1 };
		if (Arguments.mol_2_type == "lig")
			factor = get<0>(mol2_ener_res.info_energy_1mol_respair[info_ener_1mol(mol2.info[l].get_chain(), mol2.info[l].get_chainnum(), mol2.info[l].get_resnr())])
			/get<0>(mol1_ener_res.info_energy_1mol_respair[info_ener_1mol(mol1.info[0].get_chain(), mol1.info[0].get_chainnum(), mol1.info[0].get_resnr())]);
		
		factor *= get<0>(mol2_ener_res.info_energy_1mol_respair[info_ener_1mol(mol2.info[l].get_chain(), mol2.info[l].get_chainnum(), mol2.info[l].get_resnr())]);
		
		if (factor == 0 || Arguments.mol_2_type.find("protein") != string::npos
			||get<0>(mol2_ener_res.info_energy_1mol_respair[info_ener_1mol(mol2.info[l].get_chain(), mol2.info[l].get_chainnum(), mol2.info[l].get_resnr())]) * factor <= 0)
			factor = 1;

		l > 99998 ? num = 99999 : num = l;
		string chain{ " " };
		if (!isdigit(mol2.info[l].get_chain()[0]))chain = mol2.info[l].get_chain();
		output << setw(4) << "ATOM" << setw(7) << num + 1 << setw(5) << mol2.info[l].get_atom() << " " << setw(4) << left << mol2.info[l].get_original_res_name() << right << setw(1) << chain
			<< setw(4) << mol2.info[l].get_resnr() << setw(12) << fixed << mol2.m_mol[0].Molecule[l].get_x() << setw(8) << fixed << mol2.m_mol[0].Molecule[l].get_y() << setw(8)
			<< mol2.m_mol[0].Molecule[l].get_z() << setw(6) << "1.00"
			<< setw(7) <<( get<1>(mol2_ener_res.info_energy_1mol_respair[info_ener_1mol(mol2.info[l].get_chain(), mol2.info[l].get_chainnum(), mol2.info[l].get_resnr())]) / num_snap + get<2>(mol2_ener_res.info_energy_1mol_respair[info_ener_1mol(mol2.info[l].get_chain(), mol2.info[l].get_chainnum(), mol2.info[l].get_resnr())]) / num_snap)/Arguments.eu_par
			/ factor
			<< '\n';
	}
	output.close();
	return 0;
}