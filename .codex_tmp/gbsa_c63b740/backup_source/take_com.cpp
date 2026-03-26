#include"MOL.h"
//take center of mass of side chain of residues

void MOL::take_com(vector<Multi_MOL_parameters>& parameters, size_t I, size_t end, vector<index>& Index, vector<mol_info>& info)
{
	for (size_t i{}; i < Index[I].get_size_v() - 1; ++i)
	{
		double X{}, Y{}, Z{}, M{};
		for (size_t j = Index[I].get_res_pos_e(i)-1; j < Index[I].get_res_pos_e(i+1)-1; ++j)
		{
			if (info[j].is_side_chain() == false) continue;
			X += Molecule[j].get_x()*parameters[j].get_mass();
			Y += Molecule[j].get_y()*parameters[j].get_mass();
			Z += Molecule[j].get_z()*parameters[j].get_mass();
			M += parameters[j].get_mass();
		}
		Index[I].set_res_coor_e(i, X / M, Y / M, Z / M);
	}
	double X{}, Y{}, Z{}, M{};
	for (size_t j = Index[I].get_res_pos_e(Index[I].get_size_v()-1)-1; j < end -1; ++j)
	{
		if (info[j].is_side_chain() == false) continue;
		X += Molecule[j].get_x()*parameters[j].get_mass();
		Y += Molecule[j].get_y()*parameters[j].get_mass();
		Z += Molecule[j].get_z()*parameters[j].get_mass();
		M += parameters[j].get_mass();
	}
	Index[I].set_res_coor_e(Index[I].get_size_v() - 1, X / M, Y / M, Z / M);
}

// take center of geometry of residues
void MOL::take_com_all_atoms(vector<Multi_MOL_parameters>& parameters, size_t I, size_t end, vector<index>& Index, vector<mol_info>& info, 
	multi_points& com_res)
{
	for (size_t i{}; i < Index[I].get_size_v() - 1; ++i)
	{
		double X{}, Y{}, Z{}, M{};
		for (size_t j = Index[I].get_res_pos_e(i) - 1; j < Index[I].get_res_pos_e(i + 1) - 1; ++j)
		{
			X += Molecule[j].get_x() * parameters[j].get_mass();
			Y += Molecule[j].get_y() * parameters[j].get_mass();
			Z += Molecule[j].get_z() * parameters[j].get_mass();
			M += parameters[j].get_mass();
		}
		com_res.emplace(point(X / M, Y / M, Z / M), M, to_string(Index[I].get_res_num_e(i)));
	}
	double X{}, Y{}, Z{}, M{};
	for (size_t j = Index[I].get_res_pos_e(Index[I].get_size_v() - 1) - 1; j < end - 1; ++j)
	{
		X += Molecule[j].get_x() * parameters[j].get_mass();
		Y += Molecule[j].get_y() * parameters[j].get_mass();
		Z += Molecule[j].get_z() * parameters[j].get_mass();
		M += parameters[j].get_mass();
	}
	com_res.emplace(point(X / M, Y / M, Z / M), M, to_string(Index[I].get_res_num_e(Index[I].get_size_v() - 1)));
}