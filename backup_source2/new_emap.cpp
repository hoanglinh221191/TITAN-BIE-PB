#include"new_emap.h"
#include"uti_func.h"
#include<cmath>
#include<fstream>
#include<iomanip>

int multi_snap_emap::create_res_list_uni(vector<pair_res>& list_uni_res)
{
	for (auto& emap : emap_pdb)
	{
		for (auto& res : emap.ResList)
		{
			bool checkin{ false };
			for (auto& uni_res : list_uni_res)
			{
				if (uni_res.resnr_mol1 == res.resnr && uni_res.resnr_mol2 == res.Res.resnr
					&& uni_res.chain_mol1 == res.chain && uni_res.chain_mol2 == res.Res.chain)
				{
					checkin = true; break;
				}
			}
			if (checkin == false)list_uni_res.emplace_back(res.resnr, res.Res.resnr, res.chain, res.Res.chain);
		}
	}
	return 0;
}

int multi_snap_emap::create_output_files(vector<pair_res> list_uni_res, vector<fstream>& emap_files, vector<fstream>& ave_emap_files,
	  info& info_out,   folders& paths, int input_name,   arguments& Arguments)
{
	for (auto& pair : list_uni_res)
	{
		emap_files.emplace_back(fstream(paths.get_path(input_name - Arguments.b_input) + "/" + info_out.get_emap() + pair.chain_mol1 + to_string(pair.resnr_mol1)
			+ "_" + pair.chain_mol2 + to_string(pair.resnr_mol2) + ".dat", ios::out));
		ave_emap_files.emplace_back(fstream(paths.get_path(input_name - Arguments.b_input) + "/" + info_out.get_emap() + "average_" + pair.chain_mol1 + to_string(pair.resnr_mol1)
			+ "_" + pair.chain_mol2 + to_string(pair.resnr_mol2) + ".dat", ios::out));
	}
	return 0;
}

int multi_snap_emap::ana_emap(  arguments& Arguments,   info& info_out,   folders& paths, int input_name)
{
	vector<pair_res> list_uni_res;
	create_res_list_uni(list_uni_res);

	vector<fstream> emap_files, ave_emap_files;
	create_output_files(list_uni_res, emap_files, ave_emap_files, info_out, paths, input_name, Arguments);

	for (size_t i{}; i < list_uni_res.size();++i)
	{
		emap_files[i] << left << "# Time" << "  Electrostatic (" << Arguments.eu << ")  van der Waals (" << Arguments.eu << ")\n";
		energy_res_base average_result, sd_result;
		size_t num{};
		for (auto& emap : emap_pdb)
		{
			for (auto& res : emap.ResList)
			{
				if (list_uni_res[i].resnr_mol1 == res.resnr && list_uni_res[i].resnr_mol2 == res.Res.resnr
					&& list_uni_res[i].chain_mol1 == res.chain && list_uni_res[i].chain_mol2 == res.Res.chain)
				{
					emap_files[i] << setw(12) << left << emap.time_of_snap * info_out.get_current_dt() + info_out.get_current_begin_time() << right << setw(15) << res.Res.E_E / Arguments.eu_par << setw(15) << res.Res.E_V / Arguments.eu_par << '\n';
					average_result.e_e += res.Res.E_E; average_result.e_v += res.Res.E_V;
					++num;
				}
			}
		}

		average_result.e_e /= (num*Arguments.eu_par);
		average_result.e_v /= (num*Arguments.eu_par);

		for (auto& emap : emap_pdb)
		{
			for (auto& res : emap.ResList)
			{
				if (list_uni_res[i].resnr_mol1 == res.resnr && list_uni_res[i].resnr_mol2 == res.Res.resnr
					&& list_uni_res[i].chain_mol1 == res.chain && list_uni_res[i].chain_mol2 == res.Res.chain)
				{
					sd_result.e_e += square(res.Res.E_E / Arguments.eu_par - average_result.e_e);
					sd_result.e_v += square(res.Res.E_V / Arguments.eu_par - average_result.e_v);
				}
			}
		}

		string sd{ "SD" };
		if (Arguments.flag_use_CI) {
			sd_result.e_e = 1.96*sqrt(sd_result.e_e) / num;
			sd_result.e_v = 1.96*sqrt(sd_result.e_v) / num;
			sd = "Ci (95%)";
		}
		else
		{
			sd_result.e_e = sqrt(sd_result.e_e / num);
			sd_result.e_v = sqrt(sd_result.e_v / num);
		}
		ave_emap_files[i] << setw(12) << "# Number of snaps: " << num << '\n'
			<< "Average electrostatic: " << average_result.e_e << " (" << Arguments.eu << ") " << sd << ": " << sd_result.e_e << '\n'
			<< "Average van der Waals: " << average_result.e_v << " (" << Arguments.eu << ") " << sd << ": " << sd_result.e_v << '\n';
	}
	return 0;
}