#include"stat_hb_and_sc_contact.h"
#include"MOL.h"
#include<fstream>
#include<iomanip>

int traj_hb::add_hb(MOL& mol, vector<mol_info>& info_mol)
{
	for (size_t j{}; j < mol.list_h_bonds.dat.size(); ++j)
	{
		bool flag{ false };
		for (size_t i{}; i < this->residues_hb.size(); ++i)
			if (info_mol[mol.list_h_bonds.dat[j].get_index_acceptor()].get_chainnum() == this->residues_hb[i].chainid
				&& info_mol[mol.list_h_bonds.dat[j].get_index_acceptor()].get_resnr() == this->residues_hb[i].resnum)
			{
				++this->residues_hb[i].h_bond_num;
				flag = true;
				break;
			}

		if (flag == false)
		{
			this->residues_hb.emplace_back(info_mol[mol.list_h_bonds.dat[j].get_index_acceptor()].get_original_res_name()
				, info_mol[mol.list_h_bonds.dat[j].get_index_acceptor()].get_chain(), info_mol[mol.list_h_bonds.dat[j].get_index_acceptor()].get_resnr()
				, 1, info_mol[mol.list_h_bonds.dat[j].get_index_acceptor()].get_chainnum());
		}
	}
	return 0;
}

int add_hb_res(MOL& pro_in, MOL& lig, vector<mol_info>& info_rec, vector<mol_info>& info_lg,
	const arguments& Arguments, traj_hb& res_hb_rec_temp)
{
	res_hb_rec_temp.add_hb(pro_in,info_rec);

	for (size_t j{}; j < lig.list_h_bonds.dat.size(); ++j)
	{
		bool flag{ false };
		for (size_t i{}; i < res_hb_rec_temp.residues_hb.size(); ++i)
			if (info_rec[lig.list_h_bonds.dat[j].get_index_donor()].get_chainnum() == res_hb_rec_temp.residues_hb[i].chainid
				&& info_rec[lig.list_h_bonds.dat[j].get_index_donor()].get_resnr() == res_hb_rec_temp.residues_hb[i].resnum)
			{
				++res_hb_rec_temp.residues_hb[i].h_bond_num;
				flag = true;
				break;
			}

		if (flag == false)
		{
			res_hb_rec_temp.residues_hb.emplace_back(info_rec[lig.list_h_bonds.dat[j].get_index_donor()].get_original_res_name()
				, info_rec[lig.list_h_bonds.dat[j].get_index_donor()].get_chain(), info_rec[lig.list_h_bonds.dat[j].get_index_donor()].get_resnr()
				, 1, info_rec[lig.list_h_bonds.dat[j].get_index_donor()].get_chainnum());
		}
	}
	return 0;
}

bool traj_hb::check_in(const res_hb& res_hbond, size_t & pos)
{
	for (size_t i{};i<residues_hb.size();++i)
		if (this->residues_hb[i].equal_to(res_hbond) == true)
		{
			pos = i;
			return true;
		}
	return false;
}

int traj_hb::filt_hb(vector<traj_hb>& hb_reslist)
{
	for (size_t j{}; j < hb_reslist.size(); ++j)
	{
		for (size_t k{}; k < hb_reslist[j].residues_hb.size(); ++k)
		{
			size_t pos{};
			if (this->check_in(hb_reslist[j].residues_hb[k], pos) == true)
				this->residues_hb[pos].h_bond_num += hb_reslist[j].residues_hb[k].h_bond_num;
			else
				this->residues_hb.emplace_back(hb_reslist[j].residues_hb[k]);
		}
		this->snap_num += hb_reslist[j].snap_num;
	}
	return 0;
}

int stat_hb_res(vector<traj_hb>& hb_reslist, string output_filename)
{
	traj_hb stat_out_hb;

	stat_out_hb.filt_hb(hb_reslist);

	fstream output(output_filename, ios::out);
	output << "#Total snap for all trajectories: " << stat_out_hb.snap_num <<" ; traj num: "<< hb_reslist.size() << '\n';
	output << "#Chain  resnum   resname      average HBs\n";
	for (size_t i{}; i < stat_out_hb.residues_hb.size(); ++i)
	{
		output << stat_out_hb.residues_hb[i].chain << setw(5) << stat_out_hb.residues_hb[i].resnum << setw(5) << stat_out_hb.residues_hb[i].resname
			<< setw(15) << stat_out_hb.residues_hb[i].h_bond_num * 1. / stat_out_hb.snap_num << '\n';
	}
	return 0;
}