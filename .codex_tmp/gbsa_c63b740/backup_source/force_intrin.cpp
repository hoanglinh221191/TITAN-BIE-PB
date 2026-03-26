#include<string>
#include<thread>
#include<vector>
#include<fstream>
#include<iomanip>
#include<unordered_set>
#include<set>
#include <utility>
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

int force_energy_add_intrin(const cache_sigma_epsilon& CachedData, size_t* i
	, double& electro_energy, double& van_der_waals_energy, force& force_elec, force& force_vdw,
	const vector<Multi_MOL_parameters>& rec_parameters, const vector<Multi_MOL_parameters>& lg_parameters,
	const Multi_MOL& rec, const Multi_MOL& lg, size_t snap, const arguments& Arguments, size_t rec_pos = 0, size_t lg_pos = 0)
{
	if (!CachedData.parameters.empty())
	{
		__m256d elec_temp = _mm256_set_pd(CachedData.parameters[i[0]].elec_ener, CachedData.parameters[i[1]].elec_ener,
			CachedData.parameters[i[2]].elec_ener, CachedData.parameters[i[3]].elec_ener);
		__m256d vdw_temp = _mm256_set_pd(CachedData.parameters[i[0]].vdw_ener, CachedData.parameters[i[1]].vdw_ener,
			CachedData.parameters[i[2]].vdw_ener, CachedData.parameters[i[3]].vdw_ener);
		electro_energy += hsum_double_avx(elec_temp);
		van_der_waals_energy += hsum_double_avx(vdw_temp);

		if (!CachedData.force_cache.empty())
		{
			for (size_t k{}; k < 4; ++k)
			{
				force_elec += CachedData.force_cache[i[k]].forces[0];
				force_vdw += CachedData.force_cache[i[k]].forces[1];
			}
		}
	}
	else
	{

		__m256d avx_charge = _mm256_set_pd(138.935457839 * rec_parameters[rec_pos].get_charge() * lg_parameters[lg_pos].get_charge() / Arguments.epsilon_in,
			138.935457839 * rec_parameters[rec_pos].get_charge() * lg_parameters[lg_pos + 1].get_charge() / Arguments.epsilon_in,
			138.935457839 * rec_parameters[rec_pos].get_charge() * lg_parameters[lg_pos + 2].get_charge() / Arguments.epsilon_in,
			138.935457839 * rec_parameters[rec_pos].get_charge() * lg_parameters[lg_pos + 3].get_charge() / Arguments.epsilon_in);
		__m256d avx_d = _mm256_set_pd(0.01 * (square(rec.m_mol[snap].Molecule[rec_pos].get_x() - lg.m_mol[snap].Molecule[lg_pos].get_x())
			+ square(rec.m_mol[snap].Molecule[rec_pos].get_y() - lg.m_mol[snap].Molecule[lg_pos].get_y())
			+ square(rec.m_mol[snap].Molecule[rec_pos].get_z() - lg.m_mol[snap].Molecule[lg_pos].get_z())),
			0.01 * (square(rec.m_mol[snap].Molecule[rec_pos].get_x() - lg.m_mol[snap].Molecule[lg_pos + 1].get_x())
				+ square(rec.m_mol[snap].Molecule[rec_pos].get_y() - lg.m_mol[snap].Molecule[lg_pos + 1].get_y())
				+ square(rec.m_mol[snap].Molecule[rec_pos].get_z() - lg.m_mol[snap].Molecule[lg_pos + 1].get_z())),
			0.01 * (square(rec.m_mol[snap].Molecule[rec_pos].get_x() - lg.m_mol[snap].Molecule[lg_pos + 2].get_x())
				+ square(rec.m_mol[snap].Molecule[rec_pos].get_y() - lg.m_mol[snap].Molecule[lg_pos + 2].get_y())
				+ square(rec.m_mol[snap].Molecule[rec_pos].get_z() - lg.m_mol[snap].Molecule[lg_pos + 2].get_z())),
			0.01 * (square(rec.m_mol[snap].Molecule[rec_pos].get_x() - lg.m_mol[snap].Molecule[lg_pos + 3].get_x())
				+ square(rec.m_mol[snap].Molecule[rec_pos].get_y() - lg.m_mol[snap].Molecule[lg_pos + 3].get_y())
				+ square(rec.m_mol[snap].Molecule[rec_pos].get_z() - lg.m_mol[snap].Molecule[lg_pos + 3].get_z())));

		avx_d = _mm256_sqrt_pd(avx_d);


		__m256d avx_electro_energy = _mm256_div_pd(avx_charge, avx_d);

		if (CachedData.force_field_name.find("AMBER") != string::npos
			|| CachedData.force_field_name.find("CHARMM") != string::npos)
		{
			__m256d vdw_temp = _mm256_set_pd(rec_parameters[rec_pos].get_epsilon() * lg_parameters[lg_pos].get_epsilon(),
				rec_parameters[rec_pos].get_epsilon() * lg_parameters[lg_pos + 1].get_epsilon(),
				rec_parameters[rec_pos].get_epsilon() * lg_parameters[lg_pos + 2].get_epsilon(),
				rec_parameters[rec_pos].get_epsilon() * lg_parameters[lg_pos + 3].get_epsilon());
			vdw_temp = _mm256_sqrt_pd(vdw_temp);

			__m256d vdw_temp2 = _mm256_set_pd(0.5 * (lg_parameters[lg_pos].get_sigma() + rec_parameters[rec_pos].get_sigma()),
				0.5 * (lg_parameters[lg_pos + 1].get_sigma() + rec_parameters[rec_pos].get_sigma()),
				0.5 * (lg_parameters[lg_pos + 2].get_sigma() + rec_parameters[rec_pos].get_sigma()),
				0.5 * (lg_parameters[lg_pos + 3].get_sigma() + rec_parameters[rec_pos].get_sigma()));
			vdw_temp2 = _mm256_div_pd(vdw_temp2, avx_d);

			vdw_temp2 = _mm256_mul_pd(vdw_temp2, vdw_temp2);
			__m256d vdw_temp3 = vdw_temp2;
			vdw_temp2 = _mm256_mul_pd(vdw_temp2, vdw_temp2);
			vdw_temp2 = _mm256_mul_pd(vdw_temp2, vdw_temp3);
			vdw_temp3 = _mm256_fmsub_pd(vdw_temp2, vdw_temp2, vdw_temp2);
			//vdw_temp3 = _mm256_mul_pd(vdw_temp2, vdw_temp2);
			//vdw_temp3 = _mm256_sub_pd(vdw_temp3, vdw_temp2);
			vdw_temp3 = _mm256_mul_pd(vdw_temp, vdw_temp3);


			double elec_temp[4], vdw_temp_d[4], vdw_temp_d2[4];
			_mm256_storeu_pd(elec_temp, avx_electro_energy);
			//_mm256_storeu_pd(d, avx_d);
			_mm256_storeu_pd(vdw_temp_d, vdw_temp3);
			_mm256_storeu_pd(vdw_temp_d2, vdw_temp2);

			for (size_t k{}; k < 4; ++k)
			{
				electro_energy += elec_temp[k];
				van_der_waals_energy += vdw_temp_d[k];
				if (Arguments.self_int == true
					&& (vdw_temp_d2[k] > 1. || (Arguments.flag_interchain_energy == true && lg.info[lg_pos].get_chain() == rec.info[rec_pos].get_chain())))
					continue;
			}


			if (Arguments.flag_cal_force == true)
			{
				__m256d dx = _mm256_set_pd(rec.m_mol[snap].Molecule[rec_pos].get_x() - lg.m_mol[snap].Molecule[lg_pos].get_x(),
					rec.m_mol[snap].Molecule[rec_pos].get_x() - lg.m_mol[snap].Molecule[lg_pos + 1].get_x(),
					rec.m_mol[snap].Molecule[rec_pos].get_x() - lg.m_mol[snap].Molecule[lg_pos + 2].get_x(),
					rec.m_mol[snap].Molecule[rec_pos].get_x() - lg.m_mol[snap].Molecule[lg_pos + 3].get_x());


				__m256d	dy = _mm256_set_pd(rec.m_mol[snap].Molecule[rec_pos].get_y() - lg.m_mol[snap].Molecule[lg_pos].get_y(),
					rec.m_mol[snap].Molecule[rec_pos].get_y() - lg.m_mol[snap].Molecule[lg_pos + 1].get_y(),
					rec.m_mol[snap].Molecule[rec_pos].get_y() - lg.m_mol[snap].Molecule[lg_pos + 2].get_y(),
					rec.m_mol[snap].Molecule[rec_pos].get_y() - lg.m_mol[snap].Molecule[lg_pos + 3].get_y());

				__m256d	dz = _mm256_set_pd(rec.m_mol[snap].Molecule[rec_pos].get_z() - lg.m_mol[snap].Molecule[lg_pos].get_z(),
					rec.m_mol[snap].Molecule[rec_pos].get_z() - lg.m_mol[snap].Molecule[lg_pos + 1].get_z(),
					rec.m_mol[snap].Molecule[rec_pos].get_z() - lg.m_mol[snap].Molecule[lg_pos + 2].get_z(),
					rec.m_mol[snap].Molecule[rec_pos].get_z() - lg.m_mol[snap].Molecule[lg_pos + 3].get_z());

				__m256d avx_dd = _mm256_mul_pd(avx_d, avx_d);
				__m256d	temp_elec = _mm256_div_pd(avx_charge, avx_dd);
				__m256d temp_par_vdw = _mm256_div_pd(vdw_temp, avx_dd);

				double temp_par_vdw_t[4], temp_par_elec[4], dx_t[4], dy_t[4], dz_t[4];
				_mm256_storeu_pd(temp_par_elec, temp_elec);
				_mm256_storeu_pd(temp_par_vdw_t, temp_par_vdw);
				_mm256_storeu_pd(dx_t, dx); _mm256_storeu_pd(dy_t, dy); _mm256_storeu_pd(dz_t, dz);
				for (size_t k{}; k < 4; ++k)
				{
					const double c_6{ c6(vdw_temp_d2[k]) };
					double temp_vdw_t{ temp_par_vdw_t[k] * (12. * square(c_6) - 6. * c_6) };

					force_elec.add(0.1 * temp_par_elec[k] * dx_t[k], 0.1 * temp_par_elec[k] * dy_t[k], 0.1 * temp_par_elec[k] * dz_t[k]);
					force_vdw.add(0.1 * temp_par_vdw_t[k] * dx_t[k], 0.1 * temp_par_vdw_t[k] * dy_t[k], 0.1 * temp_par_vdw_t[k] * dz_t[k]);
				}
			}
		}
		else// if (CachedData.force_field_name.find("OPLS") != string::npos)
		{
			__m256d vdw_temp = _mm256_set_pd(rec_parameters[rec_pos].get_epsilon() * lg_parameters[lg_pos].get_epsilon(),
				rec_parameters[rec_pos].get_epsilon() * lg_parameters[lg_pos + 1].get_epsilon(),
				rec_parameters[rec_pos].get_epsilon() * lg_parameters[lg_pos + 2].get_epsilon(),
				rec_parameters[rec_pos].get_epsilon() * lg_parameters[lg_pos + 3].get_epsilon());
			vdw_temp = _mm256_sqrt_pd(vdw_temp);

			double rec_sigma{ rec_parameters[rec_pos].get_sigma() > 0. ? rec_parameters[rec_pos].get_sigma() : 0. },
				lg_1_sigma{ lg_parameters[lg_pos].get_sigma() > 0. ? lg_parameters[lg_pos].get_sigma() : 0. },
				lg_2_sigma{ lg_parameters[lg_pos + 1].get_sigma() > 0. ? lg_parameters[lg_pos + 1].get_sigma() : 0. },
				lg_3_sigma{ lg_parameters[lg_pos + 2].get_sigma() > 0. ? lg_parameters[lg_pos + 2].get_sigma() : 0. },
				lg_4_sigma{ lg_parameters[lg_pos + 3].get_sigma() > 0. ? lg_parameters[lg_pos + 3].get_sigma() : 0. };

			__m256d vdw_temp2 = _mm256_set_pd((lg_1_sigma * rec_sigma),
				(lg_2_sigma * rec_sigma),
				(lg_3_sigma * rec_sigma),
				(lg_4_sigma * rec_sigma));

			vdw_temp2 = _mm256_div_pd(vdw_temp2, avx_d);
			vdw_temp2 = _mm256_sqrt_pd(vdw_temp2);

			vdw_temp2 = _mm256_mul_pd(vdw_temp2, vdw_temp2);
			__m256d vdw_temp3 = vdw_temp2;
			vdw_temp2 = _mm256_mul_pd(vdw_temp2, vdw_temp2);
			vdw_temp2 = _mm256_mul_pd(vdw_temp2, vdw_temp3); //c_6
			vdw_temp3 = _mm256_fmsub_pd(vdw_temp2, vdw_temp2, vdw_temp2);
			//vdw_temp3 = _mm256_mul_pd(vdw_temp2, vdw_temp2);
			//vdw_temp3 = _mm256_sub_pd(vdw_temp3, vdw_temp2);
			vdw_temp3 = _mm256_mul_pd(vdw_temp, vdw_temp3);

			double elec_temp[4], vdw_temp_d[4], vdw_temp_d2[4];
			_mm256_storeu_pd(elec_temp, avx_electro_energy);
			//_mm256_storeu_pd(d, avx_d);
			_mm256_storeu_pd(vdw_temp_d, vdw_temp3);
			_mm256_storeu_pd(vdw_temp_d2, vdw_temp2);

			for (size_t k{}; k < 4; ++k)
			{
				electro_energy += elec_temp[k];
				van_der_waals_energy += vdw_temp_d[k];
				if (Arguments.self_int == true
					&& (vdw_temp_d2[k] > 1. || (Arguments.flag_interchain_energy == true && lg.info[lg_pos].get_chain() == rec.info[rec_pos].get_chain())))
					continue;
			}


			if (Arguments.flag_cal_force == true)
			{
				__m256d dx = _mm256_set_pd(rec.m_mol[snap].Molecule[rec_pos].get_x() - lg.m_mol[snap].Molecule[lg_pos].get_x(),
					rec.m_mol[snap].Molecule[rec_pos].get_x() - lg.m_mol[snap].Molecule[lg_pos + 1].get_x(),
					rec.m_mol[snap].Molecule[rec_pos].get_x() - lg.m_mol[snap].Molecule[lg_pos + 2].get_x(),
					rec.m_mol[snap].Molecule[rec_pos].get_x() - lg.m_mol[snap].Molecule[lg_pos + 3].get_x());


				__m256d	dy = _mm256_set_pd(rec.m_mol[snap].Molecule[rec_pos].get_y() - lg.m_mol[snap].Molecule[lg_pos].get_y(),
					rec.m_mol[snap].Molecule[rec_pos].get_y() - lg.m_mol[snap].Molecule[lg_pos + 1].get_y(),
					rec.m_mol[snap].Molecule[rec_pos].get_y() - lg.m_mol[snap].Molecule[lg_pos + 2].get_y(),
					rec.m_mol[snap].Molecule[rec_pos].get_y() - lg.m_mol[snap].Molecule[lg_pos + 3].get_y());

				__m256d	dz = _mm256_set_pd(rec.m_mol[snap].Molecule[rec_pos].get_z() - lg.m_mol[snap].Molecule[lg_pos].get_z(),
					rec.m_mol[snap].Molecule[rec_pos].get_z() - lg.m_mol[snap].Molecule[lg_pos + 1].get_z(),
					rec.m_mol[snap].Molecule[rec_pos].get_z() - lg.m_mol[snap].Molecule[lg_pos + 2].get_z(),
					rec.m_mol[snap].Molecule[rec_pos].get_z() - lg.m_mol[snap].Molecule[lg_pos + 3].get_z());

				__m256d avx_dd = _mm256_mul_pd(avx_d, avx_d);
				__m256d	temp_elec = _mm256_div_pd(avx_charge, avx_dd);
				__m256d temp_par_vdw = _mm256_div_pd(vdw_temp, avx_dd);

				double temp_par_vdw_t[4], temp_par_elec[4], dx_t[4], dy_t[4], dz_t[4];
				_mm256_storeu_pd(temp_par_elec, temp_elec);
				_mm256_storeu_pd(temp_par_vdw_t, temp_par_vdw);
				_mm256_storeu_pd(dx_t, dx); _mm256_storeu_pd(dy_t, dy); _mm256_storeu_pd(dz_t, dz);


				for (size_t k{}; k < 4; ++k)
				{
					const double c_6{ c6(vdw_temp_d2[k]) };
					double temp_vdw_t{ temp_par_vdw_t[k] * (12. * square(c_6) - 6. * c_6) };

					force_elec.add(0.1 * temp_par_elec[k] * dx_t[k], 0.1 * temp_par_elec[k] * dy_t[k], 0.1 * temp_par_elec[k] * dz_t[k]);
					force_vdw.add(0.1 * temp_par_vdw_t[k] * dx_t[k], 0.1 * temp_par_vdw_t[k] * dy_t[k], 0.1 * temp_par_vdw_t[k] * dz_t[k]);
				}


				/*double temp_vdw{sqrt(rec_parameters[rec_pos].get_epsilon() * lg_parameters[lg_pos].get_epsilon())},
					c_6{ c6(sqrt(lg_parameters[lg_pos].get_sigma() * rec_parameters[rec_pos].get_sigma()) / d) }, c_12{ square(c_6) };

				if (lg_parameters[lg_pos].get_sigma() < 0. || rec_parameters[rec_pos].get_sigma() < 0.)
					c_6 = 0.000;

				if (Arguments.self_int == true && (c_6 > 1. || (Arguments.flag_interchain_energy == true && lg.info[lg_pos].get_chain() == rec.info[i].get_chain())))
					return 1;

				van_der_waals_energy += temp_vdw * (c_12 - c_6);

				if (Arguments.flag_cal_force == true)
				{
					const double dx{ rec.m_mol[snap].Molecule[rec_pos].get_x() - lg.m_mol[snap].Molecule[lg_pos].get_x() },
						dy{ rec.m_mol[snap].Molecule[rec_pos].get_y() - lg.m_mol[snap].Molecule[lg_pos].get_y() },
						dz{ rec.m_mol[snap].Molecule[rec_pos].get_z() - lg.m_mol[snap].Molecule[lg_pos].get_z() };
					double	sd{ square(1. / d) }, temp_par_vdw{ temp_vdw * (12. * c_12 * sd - 6. * c_6 * sd) },
						temp_elec{ charge * sd };
					force_elec.add(0.1 * temp_elec * dx, 0.1 * temp_elec * dy, 0.1 * temp_elec * dz);
					force_vdw.add(0.1 * temp_par_vdw * dx, 0.1 * temp_par_vdw * dy, 0.1 * temp_par_vdw * dz);
				}*/
			}
			/*else
			{
				double sd{ 1. / square(d) }, temp_elec{ sd * charge }, c_6{ c6(d) }, c_12{ square(c_6) };
				double epsilon{ sqrt(rec_parameters[rec_pos].get_epsilon() * lg_parameters[lg_pos].get_epsilon()) },
					sigma{ sqrt(rec_parameters[rec_pos].get_sigma() * lg_parameters[lg_pos].get_sigma()) },
					temp_vdw{ 12. * sd * epsilon / c_12 - 6. * sd * sigma / c_6 };

				van_der_waals_energy += epsilon / c_12 - sigma / c_6;

				if (Arguments.flag_cal_force == true)
				{
					const double dx{ rec.m_mol[snap].Molecule[rec_pos].get_x() - lg.m_mol[snap].Molecule[lg_pos].get_x() },
						dy{ rec.m_mol[snap].Molecule[rec_pos].get_y() - lg.m_mol[snap].Molecule[lg_pos].get_y() },
						dz{ rec.m_mol[snap].Molecule[rec_pos].get_z() - lg.m_mol[snap].Molecule[lg_pos].get_z() };
					force_elec.add(0.1 * temp_elec * dx, 0.1 * temp_elec * dy, 0.1 * temp_elec * dz);
					force_vdw.add(0.1 * temp_vdw * dx, 0.1 * temp_vdw * dy, 0.1 * temp_vdw * dz);
				}
			}*/
		}
	}
	return 0;
}