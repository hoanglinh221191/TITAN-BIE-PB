#include"MOL.h"
#include<fstream>
#include<iostream>
#include<cmath>
#include"functions.h"
using namespace std;

MOL& MOL::operator=(const MOL& rhs)
{
	this->Molecule = rhs.Molecule;
	this->Com = rhs.Com;
	this->time_of_snap = rhs.time_of_snap;
	this->Box = rhs.Box;
	this->force_f = rhs.force_f;
	this->mol_type = rhs.mol_type;
	this->list_h_bonds = rhs.list_h_bonds;
	return *this;
}
MOL& MOL::operator=(MOL&& rhs) noexcept
{
	Molecule = std::move(rhs.Molecule);
	Com = std::move(rhs.Com);
	time_of_snap = rhs.time_of_snap;
	Box = std::move(rhs.Box);
	list_h_bonds = std::move(rhs.list_h_bonds);
	force_f = rhs.force_f;
	mol_type = rhs.mol_type;
	return *this;
}

void MOL::get_centroid(std::span<double> com)
{
	com[0] = 0.;com[1] = 0.;com[2] = 0.;
	for (size_t i{};i < Molecule.size();++i)
	{
		com[0] += Molecule[i].get_x();
		com[1] += Molecule[i].get_y();
		com[2] += Molecule[i].get_z();
	}
	com[0] /= Molecule.size();com[1] /= Molecule.size();com[2] /= Molecule.size();
}

int MOL::cal_dis_charged_atoms(vector<double>& distances, vector<double>& distances_negative_charges, vector<double>& distances_positive_charges,
	const vector<Multi_MOL_parameters>& parameters, std::span<double> com)
{
	for (size_t i{}; i < Molecule.size(); ++i)
	{
		if (fabs(parameters[i].get_charge()) > 0.0000001)
		{
			double d_temp{ sqrt(square(Molecule[i].get_x() - com[0]) + square(Molecule[i].get_y() - com[1]) + square(Molecule[i].get_z() - com[2])) };
			
			distances.emplace_back(d_temp);

			if (parameters[i].get_charge() > 0.0000001)
				distances_positive_charges.emplace_back(d_temp);
			else
				distances_negative_charges.emplace_back(d_temp);
		}
	}
	return 0;
}

int MOL::cal_electric_field_1atom(const vector<mol_info>& info, const arguments& Arguments, const vector<Multi_MOL_parameters>& parameters,
	force& electric_field_vector, vector<size_t>& list_heavy_atoms, size_t pos_o_water)
{
	for (size_t i{}; i < list_heavy_atoms.size(); ++i)
	{
		auto& atom = list_heavy_atoms[i];
		if (atom == pos_o_water)continue;
		const double temp_charge{ 138.935457839* parameters[atom].get_charge() };
		const double d{ 0.1*sqrt(square(Molecule[atom].get_x() - Molecule[pos_o_water].get_x())
			+ square(Molecule[atom].get_y() - Molecule[pos_o_water].get_y())
			+ square(Molecule[atom].get_z() - Molecule[pos_o_water].get_z())) };
		const double par_temp{ temp_charge / pow(d,3) };

		electric_field_vector += force(0.1*(Molecule[pos_o_water].get_x() - Molecule[atom].get_x()) * par_temp,
			0.1*(Molecule[pos_o_water].get_y() - Molecule[atom].get_y()) * par_temp, 0.1*(Molecule[pos_o_water].get_z() - Molecule[atom].get_z()) * par_temp);
	}
	return 0;
}

int MOL::find_surface_water(vector<size_t>& list_waters, vector<size_t>& list_heavy_atoms, vector<size_t>& list_surface)
{
	for (size_t i{}; i < list_waters.size(); ++i)
	{
		bool flag_surface{ false };
		for (size_t j{}; j < list_heavy_atoms.size(); ++j)
		{
			if (dis_between2_atoms(list_waters[i], list_heavy_atoms[j]) <= 3.5000)
			{
				flag_surface = true; break;
			}
		}
		if(flag_surface==true)
			list_surface.emplace_back(list_waters[i]);
	}
	return 0;
}

int MOL::find_sphere_water(vector<size_t>& list_waters, vector<size_t>& list_heavy_atoms, vector<size_t>& list_surface, double d_crit)
{
	double x{}, y{}, z{};
	for (size_t i{}; i < list_heavy_atoms.size(); ++i)
	{
		x += Molecule[list_heavy_atoms[i]].get_x();
		y += Molecule[list_heavy_atoms[i]].get_y();
		z += Molecule[list_heavy_atoms[i]].get_z();
	}
	x /= list_heavy_atoms.size(); y /= list_heavy_atoms.size(); z /= list_heavy_atoms.size();

	for (size_t i{}; i < list_waters.size(); ++i)
	{
		if (dis_between_atom_and_point(list_waters[i], x, y, z) <= d_crit)
			list_surface.emplace_back(list_waters[i]);
	}
	return 0;
}

int MOL::cal_dipoles_specific_atoms(vector<size_t>& list_atoms, const vector<Multi_MOL_parameters>& parameters, const vector<mol_info>& info,
	vector<force>& dipoles_out, string atoms_type)
{
	if (atoms_type == "water")
		cal_dipoles_water_molecules(list_atoms, parameters, info, dipoles_out);
	else return 1;
	return 0;
}

int MOL::cal_dipoles_water_molecules(vector<size_t>& list_atoms, const vector<Multi_MOL_parameters>& parameters, const vector<mol_info>& info,
	vector<force>& dipoles_out)
{
	for (size_t i{}; i < list_atoms.size(); ++i)
	{
		auto& atom = list_atoms[i];
		if (atom + 2 < Molecule.size() && info[atom + 1].get_atom_type()[0] == 'H' && info[atom + 2].get_atom_type()[0] == 'H')
		{
			dipoles_out.emplace_back(force(Molecule[atom].get_x()*parameters[atom].get_charge() + Molecule[atom + 1].get_x()*parameters[atom + 1].get_charge() + Molecule[atom + 2].get_x()*parameters[atom + 2].get_charge(),
				Molecule[atom].get_y()*parameters[atom].get_charge() + Molecule[atom + 1].get_y()*parameters[atom + 1].get_charge() + Molecule[atom + 2].get_y()*parameters[atom + 2].get_charge(),
				Molecule[atom].get_z()*parameters[atom].get_charge() + Molecule[atom + 1].get_z()*parameters[atom + 1].get_charge() + Molecule[atom + 2].get_z()*parameters[atom + 2].get_charge()));
		}
		else
			continue;
	}
	return 0;
}

int MOL::helper_cal_electric_field(const vector<mol_info>& info, const arguments& Arguments,
	const vector<Multi_MOL_parameters>& parameters, vector<force>& elec_field, vector<size_t>& list_heavy_atoms, vector<size_t>& pos_o_waters, 
	size_t begin, size_t end)
{
	for (size_t i{ begin }; i < end; ++i)
	{
		if (i >= elec_field.size())elec_field.emplace_back(force());

		cal_electric_field_1atom(info, Arguments, parameters, elec_field[i], list_heavy_atoms, pos_o_waters[i]);
	}
	return 0;
}

int MOL::cal_electric_field(const vector<mol_info>& info, const arguments& Arguments,
	const vector<Multi_MOL_parameters>& parameters, vector<size_t>& list_heavy_atoms, vector<size_t>& pos_o_waters, vector<force>& elec_fields)
{
	vector<thread> threads_cal(Arguments.nt);

	for (size_t i{}; i < threads_cal.size(); ++i)
		threads_cal[i] = thread(&MOL::helper_cal_electric_field, this, ref(info), ref(Arguments), ref(parameters), ref(elec_fields)
			, ref(list_heavy_atoms), ref(pos_o_waters), i*pos_o_waters.size() / Arguments.nt, (i + 1)*pos_o_waters.size() / Arguments.nt);
	for (auto& t : threads_cal)
		t.join();
	return 0;
}

int MOL::cal_angle_dipoles_elecfield(const vector<mol_info>& info, const arguments& Arguments, vector<force>& elec_fields,
	const vector<Multi_MOL_parameters>& parameters, vector<size_t>& list_heavy_atoms, vector<size_t>& pos_o_waters, vector<double>& cos_out)
{
	vector<force> dipoles;
	cal_dipoles_specific_atoms(pos_o_waters, parameters, info, dipoles, "water");

	size_t SIZE{ dipoles.size() >= elec_fields.size() ? elec_fields.size() : dipoles.size() };

	for (size_t i{}; i < SIZE; ++i)
	{
		cos_out.emplace_back((dipoles[i].Fx*elec_fields[i].Fx + dipoles[i].Fy*elec_fields[i].Fy + dipoles[i].Fz*elec_fields[i].Fz) / (dipoles[i].modul()*elec_fields[i].modul()));
	}
	return 0;
}

int MOL::cal_min_max_mol(double& min_x, double& max_x, double& min_y, double& max_y, double& min_z, double& max_z)
{
	min_x = 10e13; min_y = 10e13; min_z = 10e13;
	max_x = -10e13; max_y = -10e13; max_z = -10e13;
	for (auto& atom : Molecule)
	{
		if (max_x < atom.get_x())max_x = atom.get_x();
		if (min_x > atom.get_x())min_x = atom.get_x();
		if (max_y < atom.get_y())max_y = atom.get_y();
		if (min_y > atom.get_y())min_y = atom.get_y();
		if (max_z < atom.get_z())max_z = atom.get_z();
		if (min_z > atom.get_z())min_z = atom.get_z();
	}
	return 0;
}

double MOL::cal_minimum_distance_2mols(MOL& other, const arguments& Arguments)
{
	double distance{ 10e12 };
	double min_x_this{}, min_x_other{}, max_x_this{}, max_x_other{}, min_y_this{}, min_y_other{}, max_y_this{}, max_y_other{}, min_z_this{}, min_z_other{}
	, max_z_this{}, max_z_other{};

	//this->cal_min_max_mol(min_x_this, max_x_this, min_y_this, max_y_this, min_z_this, max_z_this);
	//other.cal_min_max_mol(min_x_other, max_x_other, min_y_other, max_y_other, min_z_other, max_z_other);

	min_x_this -= Arguments.d_crit_minmax_mol; max_x_this += Arguments.d_crit_minmax_mol;
	min_y_this -= Arguments.d_crit_minmax_mol; max_y_this += Arguments.d_crit_minmax_mol;
	min_z_this -= Arguments.d_crit_minmax_mol; max_z_this += Arguments.d_crit_minmax_mol;

	min_x_other -= Arguments.d_crit_minmax_mol; max_x_other += Arguments.d_crit_minmax_mol;
	min_y_other -= Arguments.d_crit_minmax_mol; max_y_other += Arguments.d_crit_minmax_mol;
	min_z_other -= Arguments.d_crit_minmax_mol; max_z_other += Arguments.d_crit_minmax_mol;

	const double ddX{ 0.7*(max_x_this - min_x_this) }, ddY{ 0.7*(max_y_this - min_y_this) }, ddZ{ 0.7*(max_z_this - min_z_this) };

	if (((min_x_other >= min_x_this && min_x_other <= max_x_this) || (max_x_other >= min_x_this && max_x_other <= max_x_this))
		&& ((min_y_other >= min_y_this && min_y_other <= max_y_this) || (max_y_other >= min_y_this && max_y_other <= max_y_this))
		&& ((min_z_other >= min_z_this &&min_z_other <= max_z_this) || (max_z_other >= min_z_this && max_z_other <= max_z_this)))
	{
		for (auto& atom : Molecule)
			for (auto& atom2 : other.Molecule)
			{
				if (fabs(atom.get_x() - atom2.get_x()) > 0.7*Box.get_a() || fabs(atom.get_y() - atom2.get_y()) > 0.7*Box.get_b()
					|| fabs(atom.get_z() - atom2.get_z()) > 0.7*Box.get_c())
					continue;
				 double dis{ };
				 if (Arguments.pbc_mem == false)
					 dis = sqrt(square(atom.get_x() - atom2.get_x()) + square(atom.get_y() - atom2.get_y()) + square(atom.get_z() - atom2.get_z()));
				 else
					 dis = fabs(atom.get_z() - atom2.get_z())+fabs(atom.get_x()-atom2.get_x())+fabs(atom.get_y()-atom2.get_y());
				if (dis < distance)
					distance = dis;
			}
	}

	return distance;
}

double MOL::cal_minimum_distance_mol_atom(MOL& other, const arguments& Arguments, size_t& index)
{
	double distance{10e12};
	molecule atom_2{ other.Molecule[index] };
	double center[3]{};

	this->get_centroid(center);
	distance = sqrt(square(center[0] - atom_2.get_x()) + square(center[1] - atom_2.get_y()) + square(center[2] - atom_2.get_z()));
	return distance;
}

int MOL::translate(double x, double y, double z)
{
	for (auto& atom : Molecule)
	{
		atom.set_x(atom.get_x() + x);
		atom.set_y(atom.get_y() + y);
		atom.set_z(atom.get_z() + z);
	}
	return 0;
}

bool MOL::mindist_mol(MOL& lig, arguments& Arguments, de_arg& DE_arg)
{
	if (Arguments.flag_de_auto_term == true && DE_arg.flag_old_vec == true)
	{
		double mindist{ 10000. };
		for (auto& a : Molecule)
			for (auto& a_rhs : lig.Molecule)
			{
				double d{ a.dis_2_atoms(a_rhs) };
				if (d < mindist)
					mindist = d;
			}
		cout << mindist << '\n';
		if (mindist > Arguments.d_standard)
			return true;
		else
			return false;
	}
	else
		return false;
}