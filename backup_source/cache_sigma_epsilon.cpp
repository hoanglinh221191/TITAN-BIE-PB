#include"cache_sigma_epsilon.h"
#include<cmath>

cache_force& cache_force::operator=(const cache_force& rhs)
{
	forces = rhs.forces;
	return *this;
}
cache_force& cache_force::operator=(cache_force&& rhs)
{
	forces = move(rhs.forces);
	return *this;
}

cache_force& cache_force::operator+=(const cache_force& rhs)
{
	if (!rhs.forces.empty() && !this->forces.empty())
	{
		this->forces[0] += rhs.forces[0];
		this->forces[1] += rhs.forces[1];
	}
	return *this;
}

cache_force& cache_force::operator+=(cache_force&& rhs)
{
	if (!rhs.forces.empty() && !this->forces.empty())
	{
		this->forces[0] += move(rhs.forces[0]);
		this->forces[1] += move(rhs.forces[1]);
	}
	return *this;
}

cache_energy& cache_energy::operator=(const cache_energy& rhs)
{
	this->elec_ener = rhs.elec_ener;
	this->vdw_ener = rhs.vdw_ener;
	return *this;
}

cache_energy& cache_energy::operator=(cache_energy&& rhs)
{
	this->elec_ener = move(rhs.elec_ener);
	this->vdw_ener = move(rhs.vdw_ener);
	return *this;
}

cache_energy& cache_energy::operator+=(const cache_energy& rhs)
{
	this->elec_ener += rhs.elec_ener;
	this->vdw_ener += rhs.vdw_ener;

	return *this;
}
cache_energy& cache_energy::operator+=(cache_energy&& rhs)
{
	this->elec_ener += move(rhs.elec_ener);
	this->vdw_ener += move(rhs.vdw_ener); 
	return *this;
}

int cache_sigma_epsilon::set_size(size_t MOL_size, const arguments& Arguments, size_t mol_1_size, size_t mol_2_size)
{
	cutoff_dis = Arguments.distance_energy_cutoff;

	force_field_name = Arguments.force_f;

	if ((Arguments.flag_cal_force == true && MOL_size >= 200000000) || (MOL_size >= 400000000) || Arguments.flag_cache==2)
	{
		return 1;
	}

	parameters = vector<cache_energy>{ MOL_size, cache_energy{Arguments} };
	if (Arguments.flag_cal_force == true)
		force_cache=vector<cache_force>{ MOL_size, cache_force{Arguments} };

	return 0;
}

/*int cache_sigma_epsilon::precal_charges(const Multi_MOL& mol1, const Multi_MOL& mol2, size_t begin, size_t end)
{
	for (auto j = begin; j < end; ++j)
		for (size_t i{}; i < mol1.m_mol[0].Molecule.size(); ++i)
		{
			charge[mol1.m_mol[0].Molecule.size()*j + i] = 138.935457839*mol1.parameters[i].get_charge()*mol2.parameters[j].get_charge();
		}
	return 0;
}
/*int cache_sigma_epsilon::precal_parameters(const para& c6c12_para, string ff_name)
{
	for (size_t i{}; i<c6c12_para.base.size();++i)
		for (size_t j{}; j < c6c12_para.base.size(); ++j)
		{
			if (ff_name.find("GROMOS") == string::npos)
				precal.emplace(c6c12_para.base[i].atom_type + c6c12_para.base[j].atom_type, cached_value(0.5 *(c6c12_para.base[i].c6_sigma + c6c12_para.base[j].c6_sigma)
					, 4.0 * sqrt(c6c12_para.base[i].c12_epsilon*c6c12_para.base[j].c12_epsilon)));
			else
				precal.emplace(c6c12_para.base[i].atom_type + c6c12_para.base[j].atom_type, cached_value(sqrt(c6c12_para.base[i].c6_sigma * c6c12_para.base[j].c6_sigma),
					sqrt(c6c12_para.base[i].c12_epsilon * c6c12_para.base[j].c12_epsilon)));
		}
	return 0;
}*/