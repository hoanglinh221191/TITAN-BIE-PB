#include<string>
#include<thread>
#include<vector>
#include<fstream>
#include<iomanip>
#include<unordered_set>
#include<set>
#include <utility>
#include<chrono>
#include"functions.h"
#include"uti_func.h"
#include<cmath>
#include<cstdlib>
#include<algorithm>
#include"miscellaneous.h"
#include"force.h"
#include"cache_sigma_epsilon.h"
#include"new_emap.h"
#include"stat_ener_res_res.h"

#ifdef __AVX__
#include <immintrin.h>
#else
#warning AVX is not available.Code will not compile!
#endif

using namespace std;
int imp_cal_energy_lg(vector<data_out>& lg_ener_inter, const Multi_MOL& rec, const Multi_MOL& lg, size_t num_snap, size_t begin, size_t end, 
	const vector<d_min> & dmin, const cache_sigma_epsilon& CachedData, bool flag_all_mode, const vector<Multi_MOL_parameters>& rec_parameters, 
	const vector<Multi_MOL_parameters>& lg_parameters,
	unordered_map<pair<size_t, size_t>, pair<size_t, size_t>, pair_hash>& map_list_res, vector<pair<size_t, size_t>>& list_res, const arguments& Arguments);
int imp_cal_energy_rec(vector<data_out>& rec_ener_inter, info_ener_res& infoener_res, const Multi_MOL& rec, const Multi_MOL& lg, size_t num_snap, size_t begin, size_t end,
	const vector<d_min> & dmin, const cache_sigma_epsilon& CachedData, bool flag_all_mode, const vector<Multi_MOL_parameters>& rec_parameters, const vector<Multi_MOL_parameters>& lg_parameters, 
	unordered_map<pair<size_t, size_t>, pair<size_t, size_t>, pair_hash>& map_list_res, vector<pair<size_t, size_t>>& list_res, const arguments& Arguments);
void print_energy_map(size_t nt, vector<vector<data_out> >& data, fstream& energy_file, bool flag_all_mode, const arguments& Arguments,
	vector<energy_res>& Energy_res, size_t num_snap, info& info_output);
//bool check_in_energy_list(vector<data_out>& lg_ener_inter, string chain, size_t resnr);

inline int add_energy(const cache_sigma_epsilon & CachedData
	, size_t i, double & electro_energy, double & van_der_waals_energy,
	const vector<Multi_MOL_parameters>& rec_parameters, const vector<Multi_MOL_parameters>& lg_parameters,
	const MOL& rec, const MOL& lg, size_t snap, const arguments& Arguments, 
	const vector<mol_info>& info_rec, const vector<mol_info>& info_lg, size_t rec_pos = 0, size_t lg_pos = 0)
{
	if (i < CachedData.parameters.size())
	{
		electro_energy += CachedData.parameters[i].elec_ener;
		van_der_waals_energy += CachedData.parameters[i].vdw_ener;
	}
	else
	{
		if (rec_pos >= rec.Molecule.size() || lg_pos >= lg.Molecule.size())
			return 1;
		double dx{ rec.Molecule[rec_pos].get_x() - lg.Molecule[lg_pos].get_x() },
			dy{ rec.Molecule[rec_pos].get_y() - lg.Molecule[lg_pos].get_y() },
			dz{ rec.Molecule[rec_pos].get_z() - lg.Molecule[lg_pos].get_z() },
			charge{ 138.935457839*rec_parameters[rec_pos].get_charge()*lg_parameters[lg_pos].get_charge() / Arguments.epsilon_in },
			d{ 0.1*sqrt(square(dx) + square(dy) + square(dz)) };

		if (CachedData.cutoff_dis > 0. && d > CachedData.cutoff_dis) return 1;

		electro_energy += charge / d;

		if (CachedData.force_field_name.find("AMBER") != string::npos
			|| CachedData.force_field_name.find("CHARMM") != string::npos)
		{
			double temp_vdw{ sqrt(rec_parameters[rec_pos].get_epsilon()*lg_parameters[lg_pos].get_epsilon()) },
				c_6{ c6(0.5*(lg_parameters[lg_pos].get_sigma() + rec_parameters[rec_pos].get_sigma()) / d) }, c_12{ square(c_6) },
				sd{ square(1. / d) }, temp_par_vdw{ temp_vdw*(12. * c_12 * sd - 6. * c_6 * sd) };

			if (lg_parameters[lg_pos].get_sigma() < 0. || rec_parameters[rec_pos].get_sigma() < 0.)
				c_6 = 0.000;
			
			if (Arguments.self_int == true && (c_6 > 1. || (Arguments.flag_interchain_energy == true && info_lg[lg_pos].get_chain() == info_rec[rec_pos].get_chain())))
				return 1;
			
			van_der_waals_energy += temp_vdw * (c_12 - c_6);
		}
		else if (CachedData.force_field_name.find("OPLS") != string::npos)
		{
			double temp_vdw{ sqrt(rec_parameters[rec_pos].get_epsilon()*lg_parameters[lg_pos].get_epsilon()) },
				c_6{ c6(sqrt(lg_parameters[lg_pos].get_sigma() * rec_parameters[rec_pos].get_sigma()) / d) }, c_12{ square(c_6) },
				sd{ square(1. / d) }, temp_par_vdw{ temp_vdw*(12. * c_12 * sd - 6. * c_6 * sd) };

			if (lg_parameters[lg_pos].get_sigma() < 0. || rec_parameters[rec_pos].get_sigma() < 0.)
				c_6 = 0.000;

			if (Arguments.self_int == true && (c_6 > 1. || (Arguments.flag_interchain_energy == true && info_lg[lg_pos].get_chain() == info_rec[rec_pos].get_chain())))
				return 1;

			van_der_waals_energy += temp_vdw * (c_12 - c_6);
		}
		else
		{
			double sd{ 1. / square(d) }, temp_elec{ sd * charge/ Arguments.epsilon_in }, c_6{ c6(d) }, c_12{ square(c_6) };
			double epsilon{ sqrt(rec_parameters[rec_pos].get_epsilon()*lg_parameters[lg_pos].get_epsilon()) },
				sigma{ sqrt(rec_parameters[rec_pos].get_sigma()*lg_parameters[lg_pos].get_sigma()) },
				temp_vdw{ 12.* sd * epsilon / c_12 - 6.* sd *sigma / c_6 };

			van_der_waals_energy += epsilon / c_12 - sigma / c_6;
		}
	}

	return 0;
}

int helper_print_energy(size_t nt, vector<vector<data_out> > & data_lg, vector<vector<data_out> >& data_rec, fstream& lg_ener_inter, fstream& rec_ener_inter
	, bool flag_all_mode, const arguments& Arguments, vector<energy_res>& Energy_res_mol1, vector<energy_res>& Energy_res_mol2, 
	size_t num_snap, info& info_output)
{
	print_energy_map(nt, data_lg, lg_ener_inter, flag_all_mode, Arguments, Energy_res_mol2, num_snap, info_output);
	print_energy_map(nt, data_rec, rec_ener_inter, flag_all_mode, Arguments, Energy_res_mol1, num_snap,info_output);
	return 0;
}

int wrapper_print_energy(size_t nt, vector<vector<data_out> > & data_lg, vector<vector<data_out> >& data_rec
	, bool flag_all_mode, const arguments& Arguments, vector<energy_res>& Energy_res_mol1, vector<energy_res>& Energy_res_mol2, 
	size_t num_snap, info& info_output)
{
	fstream lg_ener_inter("ener_lg.txt", ios::out | ios::app);
	//lg_ener_inter << "Atom" << setw(12) << "Electro." << setw(12) << "Van der W." << '\n';
	fstream rec_ener_inter("ener_rec.txt", ios::out | ios::app);
	//rec_ener_inter << "Res" << setw(12) << "Electro." << setw(12) << "Van der W." << '\n';

	if (nt >= 2)
	{
		thread print(print_energy_map, nt, ref(data_lg), ref(lg_ener_inter), flag_all_mode, ref(Arguments), ref(Energy_res_mol2), num_snap, ref(info_output));
		thread print2(print_energy_map, nt, ref(data_rec), ref(rec_ener_inter), flag_all_mode, ref(Arguments), ref(Energy_res_mol1), num_snap, ref(info_output));
		print.join();
		print2.join();
	}
	else
	{
		helper_print_energy(nt, data_lg, data_rec, lg_ener_inter, rec_ener_inter, flag_all_mode, Arguments, 
			Energy_res_mol1, Energy_res_mol2, num_snap, info_output);
	}
	return 0;
}

int cal_energy_amber_screening(const MOL& rec, const MOL& lg, const vector<mol_info>& info_rec, const vector<mol_info>& info_lg, cache_sigma_epsilon& CachedData,
	size_t begin, size_t end, cache_energy& Energy, size_t num_snap, const arguments& Arguments
	, const vector<Multi_MOL_parameters>& rec_parameters, const vector<Multi_MOL_parameters>& lg_parameters)
{
	double epsilon_in{ Arguments.epsilon_in }, lambda_debye{ 3.07 / sqrt(Arguments.salt_concentration) };
	Energy.set_energy(0., 0.);

	for (size_t j{ begin }; j < end; ++j)
	{
		double temp_epsilon_lg{ 16. * lg_parameters[j].get_epsilon() }, temp_charge{ 138.935457839 * lg_parameters[j].get_charge() },
			temp_sigma{ lg_parameters[j].get_sigma() };
		size_t rec_size{ rec.Molecule.size() };
		for (size_t i{}; i < rec_size; ++i)
		{
			double d{ 0.1 * sqrt(square(rec.Molecule[i].get_x() - lg.Molecule[j].get_x())
				+ square(rec.Molecule[i].get_y() - lg.Molecule[j].get_y())
				+ square(rec.Molecule[i].get_z() - lg.Molecule[j].get_z())) },
				c_6{ c6(0.5 * (temp_sigma + rec_parameters[i].get_sigma()) / d) };

			if (Arguments.distance_energy_cutoff > 0. && d > Arguments.distance_energy_cutoff)
			{
				size_t pos{ rec.Molecule.size() * j + i };
				CachedData.parameters[pos].set_elec_vdw(0., 0.);
				continue;
			}

			if (Arguments.self_int == true && (d < 0.5 * (temp_sigma + rec_parameters[i].get_sigma())
				|| (Arguments.flag_interchain_energy == true && info_lg[j].get_chain() == info_rec[i].get_chain())))
			{
				continue;
			}

			if (CachedData.parameters.empty())
			{
				Energy.elec_ener += temp_charge * rec_parameters[i].get_charge() *exp(-d/lambda_debye) / d / epsilon_in;
				if (temp_sigma < 0. || rec_parameters[i].get_sigma() < 0.)
					Energy.vdw_ener += sqrt(temp_epsilon_lg * rec_parameters[i].get_epsilon()) * (square(c_6) - 0.000);
				else
					Energy.vdw_ener += sqrt(temp_epsilon_lg * rec_parameters[i].get_epsilon()) * (square(c_6) - c_6);
			}
			else
			{
				size_t pos{ rec.Molecule.size() * j + i };
				if (temp_sigma < 0. || rec_parameters[i].get_sigma() < 0.)
					CachedData.parameters[pos].set_elec_vdw(temp_charge * rec_parameters[i].get_charge() * exp(-d / lambda_debye) / d / epsilon_in, sqrt(rec_parameters[i].get_epsilon() * temp_epsilon_lg) * (square(c_6) - 0.000));
				else
					CachedData.parameters[pos].set_elec_vdw(temp_charge * rec_parameters[i].get_charge() * exp(-d / lambda_debye) / d / epsilon_in, sqrt(rec_parameters[i].get_epsilon() * temp_epsilon_lg) * (square(c_6) - c_6));
				Energy += CachedData.parameters[pos];
			}
		}
	}
	return 0;
}


void cal_energy_amber(const MOL &rec, const MOL &lg, const vector<mol_info>& info_rec, const vector<mol_info>& info_lg, cache_sigma_epsilon& CachedData, 
	size_t begin, size_t end, cache_energy& Energy, size_t num_snap, const arguments& Arguments
, const vector<Multi_MOL_parameters>& rec_parameters, const vector<Multi_MOL_parameters>& lg_parameters)
{
	Energy.set_energy(0., 0.);
	for (size_t j{ begin }; j < end; ++j)
	{
		const double pos_x{ lg.Molecule[j].get_x() }, pos_y{ lg.Molecule[j].get_y() }, pos_z{ lg.Molecule[j].get_z() };
		const size_t rec_size{ rec.Molecule.size() }, size_vector_rec{ rec.Molecule.size() / 4 };
		double elec_temp[4]{}, vdw_temp_d[4]{}, d[4]{};

		__m256d avx_temp_epsilon_lg = _mm256_set1_pd(16. * lg_parameters[j].get_epsilon());
		__m256d avx_temp_charge_lg = _mm256_set1_pd(138.935457839 * lg_parameters[j].get_charge() / Arguments.epsilon_in);
		__m256d avx_temp_sigma = _mm256_set1_pd(lg_parameters[j].get_sigma());
		__m256d avx_half = _mm256_set1_pd(0.5);

		for (size_t i{}; i <= (size_vector_rec - 1) * 4; i += 4)
		{
			//__m256d avx_array_result = _mm256_set1_pd(temp_charge);

			__m256d avx_array_d_rec = _mm256_set_pd(0.01 * (square(rec.Molecule[i].get_x() - pos_x) + square(rec.Molecule[i].get_y() - pos_y) + square(rec.Molecule[i].get_z() - pos_z)),
				0.01 * (square(rec.Molecule[i + 1].get_x() - pos_x) + square(rec.Molecule[i + 1].get_y() - pos_y) + square(rec.Molecule[i + 1].get_z() - pos_z)),
				0.01 * (square(rec.Molecule[i + 2].get_x() - pos_x) + square(rec.Molecule[i + 2].get_y() - pos_y) + square(rec.Molecule[i + 2].get_z() - pos_z)),
				0.01 * (square(rec.Molecule[i + 3].get_x() - pos_x) + square(rec.Molecule[i + 3].get_y() - pos_y) + square(rec.Molecule[i + 3].get_z() - pos_z)));
			
			avx_array_d_rec = _mm256_sqrt_pd(avx_array_d_rec);

			__m256d avx_array_charge_rec = _mm256_set_pd(rec_parameters[i].get_charge(), rec_parameters[i + 1].get_charge(),
				rec_parameters[i + 2].get_charge(), rec_parameters[i + 3].get_charge());
			__m256d avx_elec_result = _mm256_mul_pd(avx_temp_charge_lg, avx_array_charge_rec);
			avx_elec_result = _mm256_div_pd(avx_elec_result, avx_array_d_rec);

			__m256d vdw_temp = _mm256_set_pd( rec_parameters[i].get_epsilon(),  rec_parameters[i + 1].get_epsilon(),
				 rec_parameters[i + 2].get_epsilon(),  rec_parameters[i + 3].get_epsilon());
			vdw_temp = _mm256_mul_pd(avx_temp_epsilon_lg, vdw_temp);
			vdw_temp = _mm256_sqrt_pd(vdw_temp);

			__m256d vdw_temp2 = _mm256_set_pd(rec_parameters[i].get_sigma(), rec_parameters[i + 1].get_sigma(), rec_parameters[i + 2].get_sigma(), rec_parameters[i + 3].get_sigma());
			vdw_temp2 = _mm256_add_pd(avx_temp_sigma, vdw_temp2);

			vdw_temp2 = _mm256_mul_pd(vdw_temp2, avx_half);
			
			vdw_temp2 = _mm256_div_pd(vdw_temp2, avx_array_d_rec);
			vdw_temp2 = _mm256_mul_pd(vdw_temp2, vdw_temp2);
			__m256d vdw_temp3 = vdw_temp2;
			vdw_temp2 = _mm256_mul_pd(vdw_temp2, vdw_temp2);
			vdw_temp2 = _mm256_mul_pd(vdw_temp2, vdw_temp3);
			vdw_temp3 = _mm256_fmsub_pd(vdw_temp2, vdw_temp2, vdw_temp2);
			vdw_temp3 = _mm256_mul_pd(vdw_temp, vdw_temp3);

			
			if (CachedData.parameters.empty())
			{
				Energy.elec_ener += hsum_double_avx(avx_elec_result);
				Energy.vdw_ener += hsum_double_avx(vdw_temp3);
			}
			else
			{
				_mm256_storeu_pd(elec_temp, avx_elec_result);
				_mm256_storeu_pd(vdw_temp_d, vdw_temp3);

				if (Arguments.distance_energy_cutoff < 0.)
				{
					Energy.elec_ener += hsum_double_avx(avx_elec_result);
					Energy.vdw_ener += hsum_double_avx(vdw_temp3);
					CachedData.set_elec_vdw(rec.Molecule.size() * j + i, elec_temp, vdw_temp_d, 4);
				}
				else
				{
					_mm256_storeu_pd(d, avx_array_d_rec);
					for (size_t k{}; k < 4; ++k)
					{
						size_t pos{ rec.Molecule.size() * j + i + k };
						if (d[k] > Arguments.distance_energy_cutoff)
						{
							CachedData.parameters[pos].set_elec_vdw(0., 0.);
							Energy += CachedData.parameters[pos];
							const double temp_sigma{ lg_parameters[j].get_sigma() };
							if (Arguments.self_int == true && (d[k] < 0.5 * (temp_sigma + rec_parameters[i + k].get_sigma())
								|| (Arguments.flag_interchain_energy == true && info_lg[j].get_chain() == info_rec[i + k].get_chain())))
							{
								continue;
							}
							continue;
						}

						CachedData.parameters[pos].set_elec_vdw(elec_temp[k], vdw_temp_d[k]);

						Energy.elec_ener += elec_temp[k];
						Energy.vdw_ener += vdw_temp_d[k];
					}
				}
				
			}

		}

		if (size_vector_rec * 4 < rec_size)
			for (size_t i{ size_vector_rec * 4 }; i < rec_size; ++i)
			{
				const double temp_charge{ 138.935457839 * lg_parameters[j].get_charge() / Arguments.epsilon_in },
					temp_epsilon{ 16. * lg_parameters[j].get_epsilon() }, temp_sigma{ lg_parameters[j].get_sigma() };

				cal_one_ele(rec, lg, CachedData, Arguments, pos_x, pos_y, pos_z, temp_sigma, temp_charge, temp_epsilon, rec_parameters, lg_parameters, info_rec, info_lg, Energy, i, j);
			}
	}
}

void cal_energy_opls(const MOL &rec, const MOL &lg, const vector<mol_info>& info_rec, const vector<mol_info>& info_lg, cache_sigma_epsilon& CachedData,
	size_t begin, size_t end, cache_energy& Energy, size_t num_snap, const arguments& Arguments
	, const vector<Multi_MOL_parameters>& rec_parameters, const vector<Multi_MOL_parameters>& lg_parameters)
{
	Energy.set_energy(0., 0.);
	for (size_t j{ begin }; j < end; ++j)
	{
		double temp_epsilon_lg{ 16. *lg_parameters[j].get_epsilon() }, temp_charge{ 138.935457839* lg_parameters[j].get_charge() / Arguments.epsilon_in },
			temp_sigma{ lg_parameters[j].get_sigma() },
			pos_x{ lg.Molecule[j].get_x() }, pos_y{ lg.Molecule[j].get_y() }, pos_z{ lg.Molecule[j].get_z() };
		size_t rec_size{ rec.Molecule.size() };

		for (size_t i{}; i < rec_size; ++i)
		{
			double d{ 0.1*sqrt(square(rec.Molecule[i].get_x() - pos_x)
				+ square(rec.Molecule[i].get_y() - pos_y)
				+ square(rec.Molecule[i].get_z() - pos_z)) },
				c_6{ c6(sqrt(temp_sigma * rec_parameters[i].get_sigma()) / d) };
			
			if (Arguments.self_int == true && (c_6 > 1. || (Arguments.flag_interchain_energy == true && info_lg[j].get_chain() == info_rec[i].get_chain())))
			{
				continue;
			}

			if (Arguments.distance_energy_cutoff > 0. && d > Arguments.distance_energy_cutoff)
			{
				size_t pos{ rec.Molecule.size() * j + i };
				CachedData.parameters[pos].set_elec_vdw(0., 0.);
				continue;
			}

			if (CachedData.parameters.empty())
			{
				Energy.elec_ener += temp_charge * rec_parameters[i].get_charge() / d;
				Energy.vdw_ener += sqrt(temp_epsilon_lg * rec_parameters[i].get_epsilon())*(square(c_6) - c_6);
			}
			else
			{
				size_t pos{ rec.Molecule.size()*j + i };
				CachedData.parameters[pos].set_elec_vdw(temp_charge * rec_parameters[i].get_charge() / d, sqrt(rec_parameters[i].get_epsilon()*temp_epsilon_lg)*(square(c_6) - c_6));
				Energy += CachedData.parameters[pos];
			}
		}
	}
}

void cal_energy_gromos(const MOL &rec, const MOL &lg, const vector<mol_info>& info_rec, const vector<mol_info>& info_lg, cache_sigma_epsilon& CachedData,
	size_t begin, size_t end, cache_energy& Energy, size_t num_snap, const arguments& Arguments
	, const vector<Multi_MOL_parameters>& rec_parameters, const vector<Multi_MOL_parameters>& lg_parameters)
{
	Energy.set_energy(0., 0.);
	for (size_t j{ begin }; j < end; ++j)
	{
		double temp_charge{ 138.935457839* lg_parameters[j].get_charge() / Arguments.epsilon_in },
			temp_epsilon{ lg_parameters[j].get_epsilon() }, temp_sigma{ lg_parameters[j].get_sigma() },
			pos_x{ lg.Molecule[j].get_x() }, pos_y{ lg.Molecule[j].get_y() }, pos_z{ lg.Molecule[j].get_z() };
		size_t rec_size{ rec.Molecule.size() };
		for (size_t i{}; i < rec_size; ++i)
		{
			double d{ 0.1*sqrt(square(rec.Molecule[i].get_x() - pos_x)
				+ square(rec.Molecule[i].get_y() - pos_y)
				+ square(rec.Molecule[i].get_z() - pos_z)) }, c_6{ c6(d) };

			if (Arguments.distance_energy_cutoff > 0. && d > Arguments.distance_energy_cutoff) 
			{
				size_t pos{ rec.Molecule.size() * j + i };
				CachedData.parameters[pos].set_elec_vdw(0., 0.);
				continue;
			}

			if (Arguments.self_int == true && (i == j || (Arguments.flag_interchain_energy == true && info_lg[j].get_chain() == info_rec[i].get_chain())))
			{
				continue;
			}

			if (CachedData.parameters.empty())
			{
				Energy.elec_ener += temp_charge * rec_parameters[i].get_charge() / d;
				Energy.vdw_ener += sqrt(rec_parameters[i].get_epsilon()*temp_epsilon) / square(c_6) - sqrt(rec_parameters[i].get_sigma()*temp_sigma) / c_6;
			}
			else
			{
				size_t pos{ rec.Molecule.size()*j + i };
				CachedData.parameters[pos].set_elec_vdw(temp_charge * rec_parameters[i].get_charge() / d, sqrt(rec_parameters[i].get_epsilon()*temp_epsilon) / square(c_6) - sqrt(rec_parameters[i].get_sigma()*temp_sigma) / c_6);
				Energy += CachedData.parameters[pos];
			}
		}
	}
}


int imp_cal_energy_lg_atom(vector<data_out>& lg_ener_inter, const Multi_MOL& rec, const Multi_MOL& lg, size_t num_snap, size_t begin, 
	size_t end, const vector<d_min> & dmin, const cache_sigma_epsilon& CachedData,
	const vector<Multi_MOL_parameters>& rec_parameters, const vector<Multi_MOL_parameters>& lg_parameters, const arguments& Arguments)
{
	for (size_t i = begin; i < end; ++i)
	{
		double E_E{}, E_V{};
		for (size_t j{}; j < rec.m_mol[num_snap].Molecule.size(); ++j)
		{
			add_energy(CachedData, rec.m_mol[0].Molecule.size()*i + j, E_E, E_V,rec_parameters,lg_parameters,rec.m_mol[num_snap],lg.m_mol[num_snap],num_snap,
				Arguments, rec.info,lg.info,j,i);
		}
		lg_ener_inter.emplace_back(i + 1, " ", lg.info[i].get_original_res_name(), E_E, E_V);
	}
	return 0;
}

int imp_cal_energy_lg(vector<data_out>& lg_ener_inter, const Multi_MOL& rec, const Multi_MOL& lg, size_t num_snap, size_t begin, size_t end, 
	const vector<d_min> & dmin, const cache_sigma_epsilon& CachedData, bool flag_all_mode,
	const vector<Multi_MOL_parameters>& rec_parameters, const vector<Multi_MOL_parameters>& lg_parameters,
	unordered_map<pair<size_t, size_t>, pair<size_t, size_t>, pair_hash>& map_list_res, vector<pair<size_t,size_t>>& list_res, const arguments& Arguments)
{
	for (size_t ite{begin};ite<end;++ite)
	{
		double E_E{}, E_V{};
		pair<size_t, size_t> ptemp{ map_list_res[list_res[ite]] };
		size_t num = dmin[ptemp.first].base[ptemp.second].atom_num - 1 + 300 > lg.m_mol[num_snap].Molecule.size() ? lg.m_mol[num_snap].Molecule.size() : dmin[ptemp.first].base[ptemp.second].atom_num - 1 + 300;

		for (size_t j = dmin[ptemp.first].base[ptemp.second].atom_num - 1; j < num; ++j)
		{
			if (lg.info[j].get_resnr() == lg.info[dmin[ptemp.first].base[ptemp.second].atom_num - 1].get_resnr()
				&& lg.info[j].get_chain() == lg.info[dmin[ptemp.first].base[ptemp.second].atom_num - 1].get_chain())
			{
				if (flag_all_mode == true)
				{
					for (size_t k{}; k < rec.m_mol[num_snap].Molecule.size(); ++k)
						add_energy(CachedData, rec.m_mol[0].Molecule.size() * j + k, E_E, E_V, rec_parameters, lg_parameters, rec.m_mol[num_snap],
							lg.m_mol[num_snap], num_snap, Arguments, rec.info, lg.info, k, j);
				}
				else
				{
					size_t beginpos{ dmin[ptemp.first].base[ptemp.second].atom_num_rec - 1 >= 0 ? dmin[ptemp.first].base[ptemp.second].atom_num_rec - 1 : 0 };
					for (size_t k{ beginpos }; k < rec.m_mol[num_snap].Molecule.size(); ++k)
					{
						if (rec.info[k].get_resnr() == rec.info[dmin[ptemp.first].base[ptemp.second].atom_num_rec - 1].get_resnr()
							&& rec.info[k].get_chain() == rec.info[dmin[ptemp.first].base[ptemp.second].atom_num_rec - 1].get_chain())
						{
							add_energy(CachedData, rec.m_mol[0].Molecule.size() * j + k, E_E, E_V, rec_parameters, lg_parameters, rec.m_mol[num_snap],
								lg.m_mol[num_snap], num_snap, Arguments, rec.info, lg.info, k, j);
						}
						else break;
					}
				}
			}
			else
				break;
		}
		
		lg_ener_inter.emplace_back(dmin[ptemp.first].base[ptemp.second].resnr_ab, dmin[ptemp.first].base[ptemp.second].chain_ab, dmin[ptemp.first].base[ptemp.second].res_name_ab, E_E, E_V);

	}

	/*unordered_set<pair<size_t, size_t>, pair_hash> energy_list;

	for (size_t ii{ begin }; ii < end; ++ii)
	{
		for (size_t i{}; i < dmin[ii].base.size(); ++i)
		{
			if (energy_list.emplace(dmin[ii].base[i].chainid_ab, dmin[ii].base[i].resnr_ab).second == false)
				continue;

			double E_E{}, E_V{};
			size_t num = dmin[ii].base[i].atom_num - 1 + 300 > lg.m_mol[num_snap].Molecule.size() ? lg.m_mol[num_snap].Molecule.size() : dmin[ii].base[i].atom_num - 1 + 300;

			for (size_t j = dmin[ii].base[i].atom_num - 1; j < num; ++j)
			{
				if (lg.info[j].get_resnr() == lg.info[dmin[ii].base[i].atom_num - 1].get_resnr()
					&& lg.info[j].get_chain() == lg.info[dmin[ii].base[i].atom_num - 1].get_chain())
				{

					if (flag_all_mode == true)
					{
						for (size_t k{}; k < rec.m_mol[num_snap].Molecule.size(); ++k)
							add_energy(CachedData, rec.m_mol[0].Molecule.size() * j + k, E_E, E_V, rec_parameters, lg_parameters, rec.m_mol[num_snap],
								lg.m_mol[num_snap], num_snap, Arguments, rec.info, lg.info, k, j);
					}
					else
					{
						size_t beginpos{ dmin[ii].base[i].atom_num_rec - 1 >= 0 ? dmin[ii].base[i].atom_num_rec - 1 : 0 };
						for (size_t k{ beginpos }; k < rec.m_mol[num_snap].Molecule.size(); ++k)
						{
							if (rec.info[k].get_resnr() == rec.info[dmin[ii].base[i].atom_num_rec - 1].get_resnr()
								&& rec.info[k].get_chain() == rec.info[dmin[ii].base[i].atom_num_rec - 1].get_chain())
							{
								add_energy(CachedData, rec.m_mol[0].Molecule.size() * j + k, E_E, E_V, rec_parameters, lg_parameters, rec.m_mol[num_snap],
									lg.m_mol[num_snap], num_snap, Arguments, rec.info, lg.info, k, j);
							}
							else break;
						}
					}
				}
				else
					break;
			}
			
			lg_ener_inter.emplace_back(dmin[ii].base[i].resnr_ab, dmin[ii].base[i].chain_ab, dmin[ii].base[i].res_name_ab, E_E, E_V);
		}
	}
	*/
	return 0;
}

int imp_cal_energy_rec(vector<data_out>& rec_ener_inter, info_ener_res& infoener_res, const Multi_MOL& rec, const Multi_MOL& lg, size_t num_snap,
	size_t begin, size_t end, const vector<d_min> & dmin, const cache_sigma_epsilon& CachedData, bool flag_all_mode,
	const vector<Multi_MOL_parameters>& rec_parameters, const vector<Multi_MOL_parameters>& lg_parameters,
	unordered_map<pair<size_t, size_t>, pair<size_t, size_t>, pair_hash>& map_list_res, vector<pair<size_t, size_t>>& list_res, const arguments& Arguments)
{
	for (size_t ite{ begin }; ite < end; ++ite)
	{
		double E_E{}, E_V{};
		pair<size_t, size_t> ptemp{ map_list_res[list_res[ite]] };
		if (flag_all_mode == true)
		{
			auto rec_info_dmin{ rec.info[dmin[ptemp.first].base[ptemp.second].atom_num_rec - 1] };
			for (size_t k{}; k < lg.m_mol[num_snap].Molecule.size(); ++k)
			{
				for (size_t j{ dmin[ptemp.first].base[ptemp.second].atom_num_rec - 1 }; j < rec.m_mol[num_snap].Molecule.size(); ++j)
				{
					//cout << rec.info[dmin[ptemp.first].base[ptemp.second].atom_num_rec-1].get_respos() << "-" << rec.info[j].get_respos() << '\n';
					if (rec.info[j].get_resnr() == rec_info_dmin.get_resnr()
						&& rec.info[j].get_chainnum() == rec_info_dmin.get_chainnum())
						//&& rec.info[dmin[ptemp.first].base[ptemp.second].atom_num_rec - 1].get_respos()== rec.info[j].get_respos())
					{
						double ee{}, ev{};
						add_energy(CachedData, rec.m_mol[0].Molecule.size() * k + j, ee, ev, rec_parameters, lg_parameters, rec.m_mol[num_snap],
							lg.m_mol[num_snap], num_snap, Arguments, rec.info, lg.info, j, k);
						
						E_E += ee; E_V += ev;
						if (Arguments.energy_map == 4)
						{
							infoener_res.add_energy(tuple(rec.info[j].get_chainnum(), lg.info[k].get_chainnum(), rec.info[j].get_resnr(), lg.info[k].get_resnr()),
								ee, ev);
						}

					}
					else break;
				}
			}
		}
		else
		{
			size_t num = dmin[ptemp.first].base[ptemp.second].atom_num_rec - 1 + 300 > rec.info.size() ? rec.info.size() : dmin[ptemp.first].base[ptemp.second].atom_num_rec - 1 + 300;
			for (size_t j = dmin[ptemp.first].base[ptemp.second].atom_num_rec - 1; j < num; ++j)
			{
				auto& rec_info_dmin = rec.info[dmin[ptemp.first].base[ptemp.second].atom_num_rec - 1];
				auto& rec_info = rec.info[j];

				if (rec_info.get_resnr() == rec_info_dmin.get_resnr()
					&& rec_info.get_chain() == rec_info_dmin.get_chain())
				{
					size_t beginpos{ dmin[ptemp.first].base[ptemp.second].atom_num - 1 >= 0 ? dmin[ptemp.first].base[ptemp.second].atom_num - 1 : 0 };
					for (size_t k{ beginpos }; k < lg.m_mol[num_snap].Molecule.size(); ++k)
						if (lg.info[k].get_resnr() == lg.info[dmin[ptemp.first].base[ptemp.second].atom_num - 1].get_resnr()
							&& lg.info[k].get_chain() == lg.info[dmin[ptemp.first].base[ptemp.second].atom_num - 1].get_chain())
						{
							double ee{}, ev{};
							add_energy(CachedData, rec.m_mol[0].Molecule.size() * k + j, ee, ev, rec_parameters, lg_parameters, rec.m_mol[num_snap],
								lg.m_mol[num_snap], num_snap, Arguments, rec.info, lg.info, j, k);
							E_E += ee; E_V += ev;

							if (Arguments.energy_map == 4)
							{
								infoener_res.add_energy(tuple(rec_info.get_chainnum(), lg.info[k].get_chainnum(), rec_info.get_resnr(), lg.info[k].get_resnr()),
									ee, ev);
							}
						}
						else
							break;
				}
			}
		}
		rec_ener_inter.emplace_back(rec.info[dmin[ptemp.first].base[ptemp.second].atom_num_rec - 1].get_resnr()
			, rec.info[dmin[ptemp.first].base[ptemp.second].atom_num_rec - 1].get_chain(),
			rec.info[dmin[ptemp.first].base[ptemp.second].atom_num_rec - 1].get_res_name(), E_E, E_V);

	}

	/*unordered_set<pair<size_t, size_t>, pair_hash> energy_list;

	for (size_t ii{begin}; ii < end; ++ii)
	{
		for (size_t i{}; i < dmin[ii].base.size(); ++i)
		{
			if (energy_list.emplace(dmin[ii].base[i].chainid_rec, dmin[ii].base[i].resnr_rec).second == false)
				continue;

			double E_E{}, E_V{};

			if (flag_all_mode == true)
			{
				for (size_t k{}; k < lg.m_mol[num_snap].Molecule.size(); ++k)
				{
					size_t num = dmin[ii].base[i].atom_num_rec - 1 + 300 > rec.info.size() ? rec.info.size() : dmin[ii].base[i].atom_num_rec - 1 + 300;
					auto& rec_info_dmin = rec.info[dmin[ii].base[i].atom_num_rec - 1];
					for (size_t j = dmin[ii].base[i].atom_num_rec - 1; j < num; ++j)
					{
						auto& rec_info = rec.info[j];
						if (rec_info.get_resnr() == rec_info_dmin.get_resnr()
							&& rec_info.get_chain() == rec_info_dmin.get_chain())
						{
							double ee{}, ev{};
							add_energy(CachedData, rec.m_mol[0].Molecule.size() * k + j, ee, ev, rec_parameters, lg_parameters, rec.m_mol[num_snap],
								lg.m_mol[num_snap], num_snap, Arguments, rec.info, lg.info, j, k);
							E_E += ee; E_V += ev;

							if (Arguments.energy_map == 4)
							{
								infoener_res.add_energy(tuple(rec_info.get_chainnum(), lg.info[k].get_chainnum(), rec_info.get_resnr(), lg.info[k].get_resnr()),
									ee, ev);
							}
						}
						else break;
					}
				}
			}
			else
			{
				size_t num = dmin[ii].base[i].atom_num_rec - 1 + 300 > rec.info.size() ? rec.info.size() : dmin[ii].base[i].atom_num_rec - 1 + 300;
				for (size_t j = dmin[ii].base[i].atom_num_rec - 1; j < num; ++j)
				{
					auto& rec_info_dmin = rec.info[dmin[ii].base[i].atom_num_rec - 1];
					auto& rec_info = rec.info[j];

					if (rec_info.get_resnr() == rec_info_dmin.get_resnr()
						&& rec_info.get_chain() == rec_info_dmin.get_chain())
					{
						size_t beginpos{ dmin[ii].base[i].atom_num - 1 >= 0 ? dmin[ii].base[i].atom_num - 1 : 0 };
						for (size_t k{ beginpos }; k < lg.m_mol[num_snap].Molecule.size(); ++k)
							if (lg.info[k].get_resnr() == lg.info[dmin[ii].base[i].atom_num - 1].get_resnr()
								&& lg.info[k].get_chain() == lg.info[dmin[ii].base[i].atom_num - 1].get_chain())
							{
								double ee{}, ev{};
								add_energy(CachedData, rec.m_mol[0].Molecule.size() * k + j, ee, ev, rec_parameters, lg_parameters, rec.m_mol[num_snap],
									lg.m_mol[num_snap], num_snap, Arguments, rec.info, lg.info, j, k);
								E_E += ee; E_V += ev;

								if (Arguments.energy_map == 4)
								{
									infoener_res.add_energy(tuple(rec_info.get_chainnum(), lg.info[k].get_chainnum(), rec_info.get_resnr(), lg.info[k].get_resnr()),
										ee, ev);
								}
							}
							else
								break;
					}
				}
			}

			rec_ener_inter.emplace_back(dmin[ii].base[i].resnr_rec, dmin[ii].base[i].chain_rec, dmin[ii].base[i].res_name_rec, E_E, E_V);
		}
	}
	//cout << "infoener res size: " << infoener_res.info_energy_respair.size() << '\n';*/
	return 0;
}

int imp_cal_energy_rec_res(multi_snap_emap& rec_ener_inter, const Multi_MOL& rec, const Multi_MOL& lg, size_t num_snap, size_t begin, 
	size_t end, const vector<d_min> & dmin, const cache_sigma_epsilon& CachedData, const arguments& Arguments)
{
	unordered_set<tuple<size_t,size_t,size_t,size_t>,tuple_hash> dmin_noduplicate;
	for (size_t ii{ 0 }; ii < dmin.size(); ++ii)
	{
		multi_emap_per_res emap_list_temp{ num_snap };
		for (size_t i{}; i < dmin[ii].base.size(); ++i)
		{
			if (dmin_noduplicate.emplace(dmin[ii].base[i].chainid_ab, dmin[ii].base[i].chainid_rec, dmin[ii].base[i].resnr_ab, dmin[ii].base[i].resnr_rec).second == false)
				continue;
			double E_E{}, E_V{};
			size_t num = dmin[ii].base[i].atom_num_rec - 1 + 250 > rec.info.size() ? rec.info.size() : dmin[ii].base[i].atom_num_rec - 1 + 250;
			//if (Arguments.flag_topo_2 == true)
			//	num = dmin[ii].base[i].atom_num_rec;
			
			for (size_t j = dmin[ii].base[i].atom_num_rec - 1; j < num; ++j)
			{
				size_t num2{ dmin[ii].base[i].atom_num - 1 + 250 > lg.info.size() ? lg.info.size() : dmin[ii].base[i].atom_num - 1 + 250 };
				//if (Arguments.flag_topo_1 == true) 
					//num2 = dmin[ii].base[i].atom_num;
				
				for (size_t k = dmin[ii].base[i].atom_num - 1; k < num2; ++k)
					if (lg.info[k].get_resnr() == lg.info[dmin[ii].base[i].atom_num - 1].get_resnr()
						&& lg.info[k].get_chain() == lg.info[dmin[ii].base[i].atom_num - 1].get_chain()
						&& rec.info[j].get_resnr() == rec.info[dmin[ii].base[i].atom_num_rec - 1].get_resnr()
						&& rec.info[j].get_chain() == rec.info[dmin[ii].base[i].atom_num_rec - 1].get_chain())
					{
						add_energy(CachedData, rec.m_mol[0].Molecule.size()*k + j, E_E, E_V, rec.parameters, lg.parameters, rec.m_mol[num_snap],
							lg.m_mol[num_snap], num_snap, Arguments, rec.info, lg.info, j, k);
					}
			}
			size_t resnr_lg{ Arguments.flag_topo_1 == true ? dmin[ii].base[i].atom_num : lg.info[dmin[ii].base[i].atom_num - 1].get_resnr() };
			size_t resnr_rec{ Arguments.flag_topo_2 == true ? dmin[ii].base[i].atom_num_rec : rec.info[dmin[ii].base[i].atom_num_rec - 1].get_resnr() };
			
			emap_list_temp.ResList.emplace_back(resnr_rec, rec.info[dmin[ii].base[i].atom_num_rec - 1].get_original_res_name()
				, rec.info[dmin[ii].base[i].atom_num_rec - 1].get_chain(), data_out(resnr_lg, lg.info[dmin[ii].base[i].atom_num - 1].get_chain(), lg.info[dmin[ii].base[i].atom_num - 1].get_original_res_name(), E_E, E_V));
		}
		rec_ener_inter.emap_pdb.emplace_back(emap_list_temp);
	}
	return 0;
}

void energyMapOf2Pro(vector<thread>& threads, info_ener_res& infoener_res, const Multi_MOL& rec, const Multi_MOL& lg, size_t num_snap, const vector<d_min> & dmin,
	const cache_sigma_epsilon& CachedData, size_t nt, bool flag_all_mode, const arguments& Arguments, 
	vector<energy_res>& Energy_res_mol1, vector<energy_res>& Energy_res_mol2, info& info_output,
	const vector<Multi_MOL_parameters>& rec_parameters, const vector<Multi_MOL_parameters>& lg_parameters)
{
	auto enermaptimestart =std::chrono::steady_clock::now();
	vector<vector<data_out> > data_lg(nt), data_rec(nt);
	vector<info_ener_res> infoener_res_v(nt);

	unordered_map<pair<size_t, size_t>, pair<size_t, size_t>, pair_hash> map_list_res, map_list_res_mol1;
	vector<pair<size_t, size_t>> list_res_mol2,list_res_mol1;
	for (size_t ii{}; ii < dmin.size(); ++ii)
	{
		for (size_t i{}; i < dmin[ii].base.size(); ++i)
		{
			//cout << dmin[ii].base[i].chainid_ab << " " << dmin[ii].base[i].resnr_ab << '\n';
			if (true == map_list_res.try_emplace(pair<size_t, size_t>(dmin[ii].base[i].chainid_ab, dmin[ii].base[i].resnr_ab), pair<size_t, size_t>(ii, i)).second)
				list_res_mol2.emplace_back(dmin[ii].base[i].chainid_ab, dmin[ii].base[i].resnr_ab);
			if (true == map_list_res_mol1.try_emplace(pair<size_t, size_t>(dmin[ii].base[i].chainid_rec, dmin[ii].base[i].resnr_rec), pair<size_t, size_t>(ii, i)).second)
				list_res_mol1.emplace_back(dmin[ii].base[i].chainid_rec, dmin[ii].base[i].resnr_rec);
		}
	}
	//cout << "size_of map mol2: " << list_res_mol2.size() << '\n';
	//cout << "size of map mol1: " << list_res_mol1.size() << '\n';

	for (size_t i{}; i < nt; ++i)
		threads[i] = thread(imp_cal_energy_lg, ref(data_lg[i]), ref(rec), ref(lg), num_snap, static_cast<size_t>(i * list_res_mol2.size() / nt),
			static_cast<size_t>((i + 1)* list_res_mol2.size() / nt), ref(dmin), ref(CachedData), flag_all_mode, ref(rec_parameters), ref(lg_parameters), ref(map_list_res), ref(list_res_mol2), ref(Arguments));
	for (auto& t : threads)
		t.join();

	for (size_t i{}; i < nt; ++i)
		threads[i] = thread(imp_cal_energy_rec, ref(data_rec[i]), ref(infoener_res_v[i]), ref(rec), ref(lg), num_snap, static_cast<size_t>(i * list_res_mol1.size() / nt),
			static_cast<size_t>((i + 1)* list_res_mol1.size() / nt), ref(dmin), ref(CachedData), flag_all_mode, ref(rec_parameters), ref(lg_parameters), ref(map_list_res_mol1), ref(list_res_mol1), ref(Arguments));
	for (auto& t : threads)
		t.join();
	//cout << "time of cal map snap " << num_snap << " (ms): "
	//	<< std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - enermaptimestart).count() << " ";
	if(Arguments.energy_map == 4)
		for (auto& ii : infoener_res_v)
		{
			for (auto& ele : ii.position_map)
			{
				infoener_res.add_energy(ele.first, get<1>(ii.info_energy_respair[ele.second]), get<2>(ii.info_energy_respair[ele.second]));
			}
		}
	
	wrapper_print_energy(nt, data_lg, data_rec, flag_all_mode, Arguments, Energy_res_mol1, Energy_res_mol2, num_snap, info_output);

	//cout << "time of print map snap " << num_snap << " (ms): "
		//<< std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - enermaptimestart).count() << "\n";
}

void energyMapOfProLig(vector<thread>& threads, info_ener_res& infoener_res, const Multi_MOL& rec, const Multi_MOL& lg, size_t num_snap, const vector<d_min> & dmin,
	const cache_sigma_epsilon& CachedData, size_t nt, bool flag_all_mode, const arguments& Arguments, vector<energy_res>& Energy_res_mol1,
	vector<energy_res>& Energy_res_mol2, info& info_output,
	const vector<Multi_MOL_parameters>& rec_parameters, const vector<Multi_MOL_parameters>& lg_parameters)
{
	vector<vector<data_out> > data_lg(nt), data_rec(nt);
	vector<info_ener_res> infoener_res_v(nt);

	unordered_map<pair<size_t, size_t>, pair<size_t, size_t>, pair_hash> map_list_res, map_list_res_mol1;
	vector<pair<size_t, size_t>> list_res_mol2, list_res_mol1;
	for (size_t ii{}; ii < dmin.size(); ++ii)
	{
		for (size_t i{}; i < dmin[ii].base.size(); ++i)
		{
			if (true == map_list_res.try_emplace(pair<size_t, size_t>(dmin[ii].base[i].chainid_ab, dmin[ii].base[i].resnr_ab), pair<size_t, size_t>(ii, i)).second)
				list_res_mol2.emplace_back(dmin[ii].base[i].chainid_ab, dmin[ii].base[i].resnr_ab);
			if (true == map_list_res_mol1.try_emplace(pair<size_t, size_t>(dmin[ii].base[i].chainid_rec, dmin[ii].base[i].resnr_rec), pair<size_t, size_t>(ii, i)).second)
				list_res_mol1.emplace_back(dmin[ii].base[i].chainid_rec, dmin[ii].base[i].resnr_rec);
		}
	}


	for (size_t i{};i < nt;++i)
		threads[i] = thread(imp_cal_energy_lg_atom, ref(data_lg[i]), ref(rec), ref(lg), num_snap, static_cast<size_t>(i*lg.m_mol[num_snap].Molecule.size() / nt), 
			static_cast<size_t>((i + 1)*lg.m_mol[num_snap].Molecule.size() / nt), ref(dmin), ref(CachedData), ref(rec_parameters), ref(lg_parameters), ref(Arguments));
	for (auto& t : threads)
		t.join();


	for (size_t i{};i < nt;++i)
		threads[i] = thread(imp_cal_energy_rec, ref(data_rec[i]), ref(infoener_res_v[i]), ref(rec), ref(lg), num_snap, static_cast<size_t>(i * dmin.size() / nt),
			static_cast<size_t>((i + 1)*dmin.size() / nt), ref(dmin), ref(CachedData), flag_all_mode, ref(rec_parameters), ref(lg_parameters), ref(map_list_res_mol1), ref(list_res_mol1), ref(Arguments));
	for (auto& t : threads)
		t.join();

	for (auto& ii : infoener_res_v)
	{
		for (auto& ele : ii.position_map)
		{
			infoener_res.add_energy(ele.first, get<1>(ii.info_energy_respair[ele.second]), get<2>(ii.info_energy_respair[ele.second]));
		}
	}

	wrapper_print_energy(nt, data_lg, data_rec, flag_all_mode, Arguments, Energy_res_mol1, Energy_res_mol2, num_snap, info_output);
}


int energyMapOf2Ligs(const Multi_MOL& rec, const Multi_MOL& lg, size_t num_snap, const vector<d_min> & dmin, const cache_sigma_epsilon& CachedData, 
	size_t nt, const arguments& Arguments, vector<energy_res>& Energy_res_mol1, vector<energy_res>& Energy_res_mol2, info& info_output,
	const vector<Multi_MOL_parameters>& rec_parameters, const vector<Multi_MOL_parameters>& lg_parameters)
{
	vector<thread> ener(nt);
	vector<vector<data_out> > data_lg(nt), data_rec(nt);


	for (size_t i{};i < nt;++i)
		ener[i] = thread(imp_cal_energy_lg_atom, ref(data_lg[i]), ref(rec), ref(lg), num_snap, static_cast<size_t>(i*lg.m_mol[num_snap].Molecule.size() / nt), 
			static_cast<size_t>((i + 1)*lg.m_mol[num_snap].Molecule.size() / nt), ref(dmin), ref(CachedData), ref(rec_parameters), ref(lg_parameters), ref(Arguments));
	for (auto& t : ener)
		t.join();

	for (size_t i{};i < nt;++i)
		ener[i] = thread(imp_cal_energy_lg_atom, ref(data_rec[i]), ref(lg), ref(rec), num_snap, static_cast<size_t>(i*rec.m_mol[num_snap].Molecule.size() / nt), 
			static_cast<size_t>((i + 1)*rec.m_mol[num_snap].Molecule.size() / nt), ref(dmin), ref(CachedData), ref(rec_parameters), ref(lg_parameters), ref(Arguments));
	for (auto& t : ener)
		t.join();

	wrapper_print_energy(1, data_lg, data_rec, true, Arguments, Energy_res_mol1, Energy_res_mol2, num_snap, info_output);

	return 0;
}

int wrap_cal_energy_map(vector<thread>& threads, Multi_MOL& mol1, Multi_MOL& mol2, arguments& Arguments, size_t num_snap, const vector<d_min> & dmin, 
	const cache_sigma_epsilon& CachedData, vector<energy_res>& Energy_res_mol1, vector<energy_res>& Energy_res_mol2, info& info_output, 
	info_ener_res& infoener_res)
{
	//cout << Arguments.mol_1_type << "-" << Arguments.mol_2_type
	//	<<" "<< Arguments.energy_map << " " << Arguments.emap_type << '\n';
	if (Arguments.mol_1_type == "lig" && Arguments.mol_2_type == "lig"  
		|| (Arguments.energy_map== 3 && Arguments.emap_type=="atom-atom"))
	{
		energyMapOf2Ligs(mol1, mol2, num_snap, dmin, CachedData, Arguments.nt, Arguments, Energy_res_mol1, Energy_res_mol2, info_output, 
			mol1.parameters, mol2.parameters);
	}
	else if (Arguments.mol_1_type != "lig" && Arguments.mol_2_type == "lig"
		|| (Arguments.energy_map == 3 && Arguments.emap_type == "residue-atom"))
	{
		//cout << "energyMapOfProLig\n";
		energyMapOfProLig(threads, infoener_res, mol1, mol2, num_snap, dmin, CachedData, Arguments.nt, Arguments.emap_level=="all", Arguments,
			Energy_res_mol1, Energy_res_mol2, info_output, mol1.parameters, mol2.parameters);
	}
	else if (Arguments.mol_1_type == "lig" && Arguments.mol_2_type !="lig"
		|| (Arguments.energy_map == 3 && Arguments.emap_type == "atom-residue"))
	{
		//cout << "energyMapOfProLig 2\n";
		energyMapOfProLig(threads, infoener_res, mol2, mol1, num_snap, dmin, CachedData, Arguments.nt, Arguments.emap_level =="all", Arguments,
			Energy_res_mol1, Energy_res_mol2, info_output, mol2.parameters, mol1.parameters);
	}
	else
	{
		//cout << "energyMapOf2Pro" << '\n';
		energyMapOf2Pro(threads, infoener_res, mol1, mol2, num_snap, dmin, CachedData, Arguments.nt, Arguments.emap_level == "all", Arguments,
			Energy_res_mol1, Energy_res_mol2, info_output, mol1.parameters, mol2.parameters);
	}

	return 0;
}