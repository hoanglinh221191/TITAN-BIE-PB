#include"miscellaneous.h"
#include"folder.h"
#include"uti_func.h"
#include"info.h"
#include<cmath>
#include<fstream>
#include<iomanip>

int emplace_energy_res(energy_res& input, vector<energy_res>& output)
{
	bool check_in{ false };
	size_t pos{};
	for (size_t i{}; i < output.size(); ++i)
	{
		if (input.resnr == output[i].resnr && input.chain == output[i].chain)
		{
			pos = i;
			check_in = true;
			break;
		}
	}

	if (check_in == true)
	{
		for (size_t i{}; i < input.base.size(); ++i)
			output[pos].base.emplace_back(input.base[i]);
	}
	else
	{
		output.emplace_back(input);
	}
	return 0;
}

int print_stat_energy_res(vector<energy_res> & average_energy_res, vector<energy_res>& sd_energy_res, info& info_output
, const arguments& Arguments, string name)
{
	fstream output_elec_energy_file("out_" + name +"_electrostatic_energy" + info_output.get_out_energy_res_motif(), ios::out),
		output_vdw_energy_file("out_" + name +"_vanderwaals_energy" + info_output.get_out_energy_res_motif(), ios::out),
		output_total_energy_file("out_" + name + "_total_energy" + info_output.get_out_energy_res_motif(), ios::out),
		output_elec_force_file, output_vdw_force_file, output_total_force_file;

	if (Arguments.print_force_res == true)
	{
		output_elec_force_file.open("out_"+name+"_electrostatic_force" + info_output.get_out_energy_res_motif(), ios::out);
		output_vdw_force_file.open("out_" + name +"_vanderwaals_force" + info_output.get_out_energy_res_motif(), ios::out);
		output_total_force_file.open("out_" + name +"_total_force" + info_output.get_out_energy_res_motif(), ios::out);
		output_elec_force_file << "# ChainResidResname  Average force (" << Arguments.eu << ")    SD (" << Arguments.eu << ")\n";
		output_total_force_file << "# ChainResidResname  Average force (" << Arguments.eu << ")    SD (" << Arguments.eu << ")\n";
		output_vdw_force_file << "# ChainResidResname  Average force (" << Arguments.eu << ")    SD (" << Arguments.eu << ")\n";
	}

	output_elec_energy_file << "# ChainResidResname  Average energy (" << Arguments.eu << ")    SD (" << Arguments.eu << ")\n";
	output_total_energy_file << "# ChainResidResname  Average energy (" << Arguments.eu << ")    SD (" << Arguments.eu << ")\n";
	output_vdw_energy_file << "# ChainResidResname  Average energy (" << Arguments.eu << ")    SD (" << Arguments.eu << ")\n";

	for (size_t i{}; i < average_energy_res.size(); ++i)
	{	
		for (size_t j{}; j < average_energy_res[i].base.size(); ++j)
		{
			output_elec_energy_file << average_energy_res[i].chain << setw(5)<< left << average_energy_res[i].resnr
				<< setw(4) << average_energy_res[i].resname << right
				<< setw(12) << setprecision(6) << average_energy_res[i].base[j].e_e / Arguments.eu_par
				<< setw(12) << setprecision(6) << sd_energy_res[i].base[j].e_e / Arguments.eu_par << '\n';

			output_vdw_energy_file<< average_energy_res[i].chain << setw(5) << left << average_energy_res[i].resnr
				<< setw(4) << average_energy_res[i].resname << right
				<< setw(12) << setprecision(6) << average_energy_res[i].base[j].e_v / Arguments.eu_par
				<< setw(12) << setprecision(6) << sd_energy_res[i].base[j].e_v / Arguments.eu_par << '\n';

			output_total_energy_file << average_energy_res[i].chain << setw(5) << left << average_energy_res[i].resnr
				<< setw(4) << average_energy_res[i].resname << right
				<< setw(12) << setprecision(6) << (average_energy_res[i].base[j].e_e + average_energy_res[i].base[j].e_v) / Arguments.eu_par
				<< setw(12) << setprecision(6) << (sd_energy_res[i].base[j].e_e + sd_energy_res[i].base[j].e_v) / Arguments.eu_par << '\n';

			if (Arguments.print_force_res == true)
			{
				output_elec_force_file << average_energy_res[i].chain << setw(5) << left << average_energy_res[i].resnr
					<< setw(4) << average_energy_res[i].resname << right
					<< setw(12) << setprecision(6) << force(average_energy_res[i].base[j].F_elec / Arguments.eu_par).modul()
					<< setw(12) << setprecision(6) << force(sd_energy_res[i].base[j].F_elec / Arguments.eu_par).modul() << '\n';

				output_vdw_force_file << average_energy_res[i].chain << setw(5) << left << average_energy_res[i].resnr
					<< setw(4) << average_energy_res[i].resname << right
					<< setw(12) << setprecision(6) << force(average_energy_res[i].base[j].F_vdw / Arguments.eu_par).modul()
					<< setw(12) << setprecision(6) << force(sd_energy_res[i].base[j].F_vdw / Arguments.eu_par).modul() << '\n';

				output_total_force_file << average_energy_res[i].chain << setw(5) << left << average_energy_res[i].resnr
					<< setw(4) << average_energy_res[i].resname << right
					<< setw(12) << setprecision(6) << force((average_energy_res[i].base[j].F_elec + average_energy_res[i].base[j].F_vdw) / Arguments.eu_par).modul()
					<< setw(12) << setprecision(6) << force((sd_energy_res[i].base[j].F_elec + sd_energy_res[i].base[j].F_vdw) / Arguments.eu_par).modul() << '\n';
			}
		}
	}
	return 0;
}

int stat_energy_res(vector<vector<energy_res>>& Energy_stat_res,
	int file_num, const arguments& Arguments, info& info_output, string mol_name)
{
	vector<energy_res> out_energy_res;

	for (size_t i{}; i < Energy_stat_res.size(); ++i)
	{
		for (size_t j{}; j < Energy_stat_res[i].size(); ++j)
		{
			emplace_energy_res(Energy_stat_res[i][j], out_energy_res);
		}
	}

	vector<energy_res> average_energy_res, sd_energy_res;

	for (size_t i{}; i < out_energy_res.size(); ++i)
	{
		double e_elec{}, e_vdw{};
		force f_elec{}, f_vdw{};
		for (size_t j{}; j < out_energy_res[i].base.size(); ++j)
		{
			e_elec += out_energy_res[i].base[j].e_e;
			e_vdw += out_energy_res[i].base[j].e_v;
			f_elec += out_energy_res[i].base[j].F_elec;
			f_vdw += out_energy_res[i].base[j].F_vdw;
		}

		e_elec /= out_energy_res[i].base.size();
		e_vdw /= out_energy_res[i].base.size();
		f_elec /= 1.* out_energy_res[i].base.size();
		f_vdw /= 1.* out_energy_res[i].base.size();

		average_energy_res.emplace_back(out_energy_res[i].resnr, out_energy_res[i].chain, out_energy_res[i].resname, e_elec, e_vdw, f_elec, f_vdw, 0.);

		double sd_e_elec{}, sd_e_vdw{};
		force sd_f_elec{}, sd_f_vdw{};
		for (size_t j{}; j < out_energy_res[i].base.size(); ++j)
		{
			sd_e_elec += square(out_energy_res[i].base[j].e_e - e_elec);
			sd_e_vdw += square(out_energy_res[i].base[j].e_v - e_vdw);
			sd_f_elec += force(square(out_energy_res[i].base[j].F_elec.Fx - f_elec.Fx), square(out_energy_res[i].base[j].F_elec.Fy - f_elec.Fy)
				, square(out_energy_res[i].base[j].F_elec.Fz - f_elec.Fz));
			sd_f_vdw += force(square(out_energy_res[i].base[j].F_vdw.Fx - f_vdw.Fx), square(out_energy_res[i].base[j].F_vdw.Fy - f_vdw.Fy)
				, square(out_energy_res[i].base[j].F_vdw.Fz - f_vdw.Fz));
		}

		sd_e_elec = sqrt(sd_e_elec / out_energy_res[i].base.size());
		sd_e_vdw = sqrt(sd_e_vdw / out_energy_res[i].base.size());
		sd_energy_res.emplace_back(out_energy_res[i].resnr, out_energy_res[i].chain, out_energy_res[i].resname, sd_e_elec, sd_e_vdw, sd_f_elec, sd_f_vdw, 0.);
	}

	print_stat_energy_res(average_energy_res, sd_energy_res, info_output, Arguments, mol_name);
	return 0;
}