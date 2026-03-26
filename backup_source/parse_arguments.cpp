#include<string>
#include<vector>
#include<iostream>
#include<sstream>
#include<fstream>
#include<algorithm>
#include<fstream>
#include<cctype>
#include<regex>
#include"MOL.h"
#include"arguments.h"
#include"functions.h"
using namespace std;
int arguments::parse_arguments(int argc, char** argv)
{
	argc = argc;
	for (int i = 1; i < argc; ++i)
	{
		if (string(argv[i]) == "-cscd")
		{
			cal_dis_sc_contact = true;
		}
		if (string(argv[i]) == "-opdb" || string(argv[i]) == "-opdbsep")
		{
			opdb = true;
			if (string(argv[i]) == "-opdbsep") flag_separation_results = true;
			if (i + 1 < argc)
			{
				if (string(argv[i + 1]).find(".pdb"))
					out_pdb_name = string(argv[i + 1]);
			}
		}
		if (string(argv[i]) == "-pqr")
		{
			flag_pqr = 1;
		}
		if (string(argv[i]) == "-interchain-energy")
		{
			flag_interchain_energy = true;
		}
		if (string(argv[i]) == "-dipole-moment")
		{
			flag_cal_dipole_moment = true;
		}
		if (string(argv[i]) == "-act")
		{
			flag_all_side_chain_ct = 0;
		}
		if (string(argv[i]) == "-useForceOptionforEnergy")
		{
			flag_energy_same_force = true;
		}
		if (string(argv[i]) == "-h" || string(argv[i]) == "--help")
		{
			help = true;
			if (help_print() == 0)
				return 2;
		}
		if (string(argv[i]) == "-v" || string(argv[i]) == "--version")
		{
			//cout << "Build: " << get_build_number() << '\n';
			return 2;
		}
		if (string(argv[i]) == "-pbc")
		{
			solv_pbc = 1;
		}
		if (string(argv[i]) == "-cat-result")
		{
			flag_cat_result_multi_infile = true;
		}
		if (string(argv[i]) == "-pbci")
		{
			solv_pbc_int = 1;
			solv_pbc = 1;
			if (i + 1 < argc && isdigit(argv[i + 1][0]))
			{
				solv_pbc = 0;
				chain_ref = atoi(argv[i + 1]);
			}
			else
			{
				cout << "WARNING!! The option -pbci lacks reference chain!\n";
				solv_pbc_int = 0;
				solv_pbc = 0;
			}
		}
		if (string(argv[i]) == "-pbci2")
		{
			solv_pbc_int = 1;
			solv_pbc_int2 = 1;
			solv_pbc = 1;
			if (i + 1 < argc && isdigit(argv[i + 1][0]))
			{
				solv_pbc = 0;
				chain_ref = atoi(argv[i + 1]);
			}
			else
			{
				cout << "WARNING!! The option -pbci2 lacks reference chain!\n";
				solv_pbc_int = 0;
				solv_pbc_int2 = 0;
				solv_pbc = 0;
			}
		}
		if (string(argv[i]) == "-pbci3")
		{
			solv_pbc_int = 1;
			solv_pbc_int2 = 1;
			solv_pbc = 1;
			solv_pbc_int3 = 1;
			if (i + 1 < argc && isdigit(argv[i + 1][0]))
			{
				solv_pbc = 0;
				chain_ref = atoi(argv[i + 1]);
			}
			else
			{
				cout << "WARNING!! The option -pbci3 lacks reference chain!\n";
				solv_pbc_int = 0;
				solv_pbc_int2 = 0;
				solv_pbc_int3 = 0;
				solv_pbc = 0;
			}
		}
		if (string(argv[i]) == "-pbca")
		{
			solv_pbc_int = 1;
			solv_pbc = 1;
			pbc_aggressive = true;
			if (i + 1 < argc && isdigit(argv[i + 1][0]))
			{
				solv_pbc = 0;
				chain_ref = atoi(argv[i + 1]);
				if (chain_ref < 0)chain_ref = 0;
			}
			else
			{
				solv_pbc_int = 0;
				solv_pbc = 0;
			}
		}
		if (string(argv[i]) == "-pbcclust")
		{
			solv_pbc_clust = true;
			if (i + 1 < argc && string(argv[i + 1]).find("-") == string::npos)
			{
				string temp{ string(argv[i + 1]) };
				for (size_t i{}; i < temp.size(); ++i)
					pbc_clust_axis.emplace_back(temp[i]);
			}
			else
				solv_pbc_clust = false;
		}
		if (string(argv[i]) == "-1pdb")
		{
			flag_1_pdb = true;
		}
		if (string(argv[i]) == "-hb")
		{
			flag_h_bond = 1;
			if (argc > i + 1 && string(argv[i + 1]) == "side-chain")
				flag_sidechain_hbond = true;
		}
		if (string(argv[i]) == "-entropy")
		{
			entropy = true;
		}
		if (string(argv[i]) == "-nocal")
		{
			flag_nocal = true;
			flag_1_pdb = true;
		}
		if (string(argv[i]) == "-all-subfolders")
		{
			flag_all_folders = true;
			flag_separation_results = true;
		}
		if (string(argv[i]) == "-split")
		{
			flag_separation_results = true;
		}
		if (string(argv[i]) == "-cal-force")
		{
			flag_cal_force = true;
		}
		if (string(argv[i]) == "-cal-volume")
		{
			flag_cal_vol = true;
		}
		if (string(argv[i]) == "-cache")
		{
			flag_cache = 1;
			if (argc > i + 1 && string(argv[i + 1]) == "-no")
				flag_cache = 2;
		}
		if (string(argv[i]) == "-max-force")
		{
			flag_maxforce = true;
		}
		if (string(argv[i]) == "-new-emap")
		{
			flag_new_emap = true;
		}
		if (string(argv[i]) == "-energy-res")
		{
			print_energy_res = true;
			energy_map = 1;
		}
		if (string(argv[i]) == "-force-res")
		{
			print_force_res = true;
			force_map = 1;
		}
		if (string(argv[i]) == "-use-CI")
		{
			flag_use_CI = true;
		}
		if (string(argv[i]) == "-elec-screening")
		{
			flag_elec_screening = true;
		}
		if (string(argv[i]) == "-charged_atom_dis")
		{
			flag_charged_atom_distance_dis = true;
			if (i + 1 < argc && argv[i + 1][0] != '-' && isdigit(argv[i + 1][0]) != 0)
				number_charged_distances_bin = atoi(argv[i + 1]);
		}
		if (string(argv[i]) == "-charged_distances_min")
		{
			flag_charged_atom_distance_dis = true;
			if (i + 1 < argc && argv[i + 1][0] != '-' && isdigit(argv[i + 1][0]) != 0)
				charged_distances_min = atof(argv[i + 1]);
		}
		if (string(argv[i]) == "-charged_distances_max")
		{
			flag_charged_atom_distance_dis = true;
			if (i + 1 < argc && argv[i + 1][0] != '-' && isdigit(argv[i + 1][0]) != 0)
				charged_distances_max = atof(argv[i + 1]);
		}
		if (string(argv[i]) == "-resdis2com")
		{
			flag_cal_dis_res_2com = true;
		}
		if (string(argv[i]) == "-angle_dipole_efield")
		{
			flag_angle_dipole_efield = true;
			if (i + 1 < argc && argv[i + 1][0] != '-' && isdigit(argv[i + 1][0]) != 0)
				number_bin_angle_dipole_efield = atoi(argv[i + 1]);
		}
		if (string(argv[i]) == "-side_chain_contact"&&i + 1 < argc)
		{
			if (string(argv[i + 1]) == "full-res")
				flag_sidechain = false;
			else
				flag_sidechain = true;
		}
		if (string(argv[i]) == "-noprint-side_chain_contact" && i + 1 < argc)
		{
			flag_print_sidechain = false;
		}
		if (string(argv[i]) == "-mol-1-type" && i + 1 < argc)
		{
			if (string(argv[i + 1]) == "protein" || string(argv[i + 1]) == "complex" || string(argv[i + 1]) == "lig")
				mol_1_type = string(argv[i + 1]);
		}
		if (string(argv[i]) == "-mol-2-type" && i + 1 < argc)
		{
			if (string(argv[i + 1]) ==  "protein" || string(argv[i + 1]) == "protein-user" || string(argv[i + 1]) == "complex" || string(argv[i + 1]) == "lig")
				mol_2_type = string(argv[i + 1]);
		}
		if (string(argv[i]) == "-com_energy" && i + 1 < argc)
		{
			if (string(argv[i + 1]) == "x")
				com_dis_energy_2mol = 1;
			else if (string(argv[i + 1]) == "y")
				com_dis_energy_2mol = 2;
			else if (string(argv[i + 1]) == "z")
				com_dis_energy_2mol = 3;
			else if (string(argv[i + 1]) == "all")
				com_dis_energy_2mol = 0;
			else com_dis_energy_2mol = 4;
		}
		if (string(argv[i]) == "-b" && i + 1 < argc)
		{
			begin_time = atof(argv[i + 1]);
			flag_begin_set = true;
		}
		if (string(argv[i]) == "-e" && i + 1 < argc)
		{
			end_time = atof(argv[i + 1]);
			flag_end_set = true;
		}
		if (string(argv[i]) == "-dt"&&i + 1 < argc)
		{
			dt = atoi(argv[i + 1]);
		}
		if (string(argv[i]) == "-mol2" && i + 1 < argc)
		{
			begin_mol_2 = atoi(argv[i + 1]);
		}
		if (string(argv[i]) == "-topo" && i + 1 < argc)
		{
			topo_name = argv[i + 1];
			if (check_topo(topo_name, 1) == true)
				flag_topo_1 = true;
		}
		if (string(argv[i]) == "-topo2"&&i + 1 < argc)
		{
			topo_name_2 = argv[i + 1];
			if (check_topo(topo_name_2, 2) == true)
				flag_topo_2 = true;
		}
		if (string(argv[i]) == "-f" && i + 1 < argc)
		{
			pdb_name = argv[i + 1];
			flag_1_pdb = true;
		}
		if (string(argv[i]) == "-n" && i + 1 < argc)
		{
			index_name = argv[i + 1];
			flag_index = true;
		}
		if (string(argv[i]) == "-eu"&&i + 1 < argc)
		{
			eu = argv[i + 1];
		}
		if (string(argv[i]) == "-fu" && i + 1 < argc)
		{
			fu = argv[i + 1];
		}

		if (string(argv[i]) == "-ff" && i + 1 < argc)
		{
			force_f = argv[i + 1];
			transform(force_f.begin(), force_f.end(), force_f.begin(), ::toupper);
		}
		if (string(argv[i]) == "-epsilon_in" && i + 1 < argc)
		{
			epsilon_in = atof(argv[i + 1]);
			++argc;
		}
		if (string(argv[i]) == "-g" && i + 1 < argc)
		{
			auto found = string(argv[i + 1]).find("&");
			auto found_2 = string(argv[i + 1]).find("-");
			if (found != std::string::npos)
			{
				group_1 = stoi(string(argv[i + 1]).substr(0, found));
				group_2 = stoi(string(argv[i + 1]).substr(found + 1));
			}
			else if (found_2 != std::string::npos)
			{
				group_1 = stoi(string(argv[i + 1]).substr(0, found_2));
				group_2 = stoi(string(argv[i + 1]).substr(found_2 + 1));
			}
			else
			{
				cout << "Your argument -g is invalid!\n";
			}
			if (group_1 == group_2)self_int = true;
		}
		if (string(argv[i]) == "-sg"&&i + 1 < argc)
		{
			auto found = string(argv[i + 1]).find("-");
			if (found == string::npos)
			{
				cout << "Your argument -sg is invalid!\n"
					<< "You should use - to separate groups of two molecule!\n";
			}
			else
			{
				string groups_1{ string(argv[i + 1]).substr(0, string(argv[i + 1]).find("-")) },
					groups_2{ string(argv[i + 1]).substr(string(argv[i + 1]).find("-") + 1) };
				
				regex reg{ "(&?)([[:blank:]]*)(\\d*)([[:blank:]]*)(&*)" };
				
				sregex_iterator pos(groups_1.cbegin(), groups_1.cend(), reg), pos_2(groups_2.cbegin(),groups_2.cend(),reg), end;
				for (; pos != end; ++pos)
				{
					if (!pos->str(3).empty() && isdigit(pos->str(3)[0]))
					sub_group_1.emplace_back(stoi(pos->str(3)));
				}

				for (; pos_2 != end; ++pos_2)
				{
					if (!pos_2->str(3).empty() && isdigit(pos_2->str(3)[0]))
					sub_group_2.emplace_back(stoi(pos_2->str(3)));
				}
				flag_sub_groups = true;
			}
		}
		if (string(argv[i]) == "-emap")
		{
			energy_map = 1;
			if (i + 1 < argc
				&& (string(argv[i + 1]) == "residue-atom" || string(argv[i + 1]) == "residue-residue" || string(argv[i + 1]) == "atom-atom" || string(argv[i + 1]) == "atom-residue"))
				emap_type = string(argv[i + 1]);
		}
		if (string(argv[i]) == "-emap-level" && i + 1 < argc)
		{
			if (string(argv[1 + 1]) == "all" || string(argv[i + 1]) == "per-res")
				emap_level = string(argv[i + 1]);
		}
		if (string(argv[i]) == "-fmap")
		{
			force_map = 1;
			if (i + 1 < argc
				&& (string(argv[i + 1]) == "residue-atom" || string(argv[i + 1]) == "residue-residue" || string(argv[i + 1]) == "atom-atom" || string(argv[i + 1]) == "atom-residue"))
				fmap_type = string(argv[i + 1]);
		}
		if (string(argv[i]) == "-fmap-level" && i + 1 < argc)
		{
			if (string(argv[1 + 1]) == "all" || string(argv[i + 1]) == "per-res")
				fmap_level = string(argv[i + 1]);
		}
		if (string(argv[i]) == "-nt" && i + 1 < argc)
		{
			if (stoi(string(argv[i + 1])) < 0)nt = 0;
			else
				nt = stoi(string(argv[i + 1]));
		}
		if (string(argv[i]) == "-dsc" && i + 1 < argc)
		{
			d_standard = atof(argv[i + 1]);
		}
		if (string(argv[i]) == "-denfo" && i + 1 < argc)
		{
			distance_energy_cutoff = atof(argv[i + 1]);
			distance_energy_cutoff *= 0.1;
		}
		if (string(argv[i]) == "-path-a" && i + 1 < argc)
		{
			path_amino_acid = argv[i + 1];
			if (path_amino_acid[path_amino_acid.length() - 1] != '/')
				path_amino_acid += '/';
		}
		if (string(argv[i]) == "-path-n"&&i + 1 < argc)
		{
			path_nonbonded = argv[i + 1];
			if (path_nonbonded[path_nonbonded.length() - 1] != '/')
				path_nonbonded += '/';
		}
		if (string(argv[i]) == "-pattern-input" && i + 1 < argc)
		{
			flag_pattern_input = true;
			flag_1_pdb = true;
			if (string(argv[i + 1]).find("*") != string::npos)
			{
				input_pattern = string(argv[i + 1]);
				if (i + 2 < argc && argv[i + 2][0] != '-')
					b_input = atoi(argv[i + 2]);
				if (i + 3 < argc&&argv[i + 3][0] != '-')
					e_input = atoi(argv[i + 3]);
				if (i + 4 < argc&&argv[i + 4][0] != '-')
					step_input = atoi(argv[i + 4]);
				if (e_input != 10000000)
					flag_range_input = true;
			}
			else
				cout << "Can not parse -pattern-input option, it must contain * character!\n";
		}
		if (string(argv[i]) == "-ahb" &&i + 1 < argc)
		{
			angle_h = atof(argv[i + 1]);
		}
		if (string(argv[i]) == "-dhb" && i + 1 < argc)
		{
			dis_h = atof(argv[i + 1]);
		}
		if (string(argv[i]) == "-entropy-quasi" && i + 3 < argc)
		{
			entropy_quasi = true;
			if (isdigit(argv[i + 1][0]) && isdigit(argv[i + 2][0]) && isdigit(argv[i + 3][0]))
			{
				group_entropy_com = atoi(argv[i + 1]);
				group_entropy_lg = atoi(argv[i + 2]);
				group_entropy_rec = atoi(argv[i + 3]);
			}
			else
			{
				cout << "-entropy-quasi option is invalid\n";
				entropy_quasi = false;
			}
		}
		if (string(argv[i]) == "-a"&&i + 1 < argc)
		{
			arg_filename = string(argv[i + 1]);
			parse_arg_file();
			return 0;
		}
		if (string(argv[i]) == "-i" && i + 1 < argc)
		{
			list_filename = string(argv[i + 1]);
			fstream test_listFile(list_filename, ios::in);
			if (test_listFile.is_open())
				flag_listfile = true;
			else
			{
				flag_listfile = false;
				cout << "Cannot open " << list_filename << " file!\n";
				exit(EXIT_FAILURE);
			}
			test_listFile.close();
		}
		if (string(argv[i]) == "-water-type" && i + 1 < argc)
		{
			if(string(argv[i+1])=="surface" && string(argv[i+1])=="sphere")
				water_type = string(argv[i + 1]);
			else
			{
				cout << "WARNING!!! -water-type option is incorrect!\n";
			}
		}
		if (string(argv[i]) == "-dis-water" && i + 1 < argc)
		{
			d_crit_water = atof(argv[i + 1]);
		}
		if (string(argv[i]) == "-all-respair" && i + 1 < argc)
		{
			flag_all_respair = true;
		}
		if (string(argv[i]) == "-pull-direct")
		{
			flag_find_pulling_direct = true;
		}
		if (string(argv[i]) == "-de-diff-weight" && i + 1 < argc)
		{
			de_diff_weight = atof(argv[i + 1]);
		}
		if (string(argv[i]) == "-de-cross-prob" && i + 1 < argc)
		{
			de_cross_prob = atof(argv[i + 1]);
		}
		if (string(argv[i]) == "-de-F" && i + 1 < argc)
		{
			de_F = atof(argv[i + 1]);
		}
		if (string(argv[i]) == "-de-scan-radius" && i + 1 < argc)
		{
			de_scan_radius = atof(argv[i + 1]);
		}
		if (string(argv[i]) == "-de-min-scan-radius" && i + 1 < argc)
		{
			de_min_scan_radius = atof(argv[i + 1]);
		}
		if (string(argv[i]) == "-de-score-function" && i + 1 < argc)
		{
			de_flag_score_function = atoi(argv[i + 1]);
		}
		if (string(argv[i]) == "-de-pop-num" && i + 1 < argc)
		{
			de_pop_num = stoi(string(argv[i + 1]));
		}
		if (string(argv[i]) == "-de-ite-num" && i + 1 < argc)
		{
			de_ite = stoi(string(argv[i + 1]));
		}
		if (string(argv[i]) == "-de-algorithm" && i + 1 < argc)
		{
			de_algorithm = stoi(string(argv[i + 1]));
		}
		if (string(argv[i]) == "-de-angle" && i + 1 < argc)
		{
			de_restrain_angle = atof(argv[i + 1]);
		}
		if (string(argv[i]) == "-de-auto-term")
		{
			flag_de_auto_term = true;
		}
		command += argv[i];
		command += " ";
	}
	return 0;
}
int arguments::parse_arg_file()
{
	fstream in_file(arg_filename, ios::in);
	if (!in_file.is_open())
	{
		cout << "Parameter file " << arg_filename <<" is not found!\n";
		return 1;
	}
	file_par = true;
	for (string temp{}; getline(in_file,temp);)
	{
		temp = temp.substr(0, temp.find("#")); if (temp.length() == 0)continue;
		if (temp.length() > 1 && temp[temp.length() - 1] == '\r')
			temp = temp.substr(0, temp.length() - 1);

		if (temp.find("-cscd") != string::npos && temp.find("-f") == string::npos)
		{
			if (temp.find("=") == string::npos)
				cal_dis_sc_contact = true;
			else
				if (temp.substr(temp.find("=") + 1).find("false") || stoi(temp.substr(temp.find("=") + 1)) == 0)
					cal_dis_sc_contact = false;
				else if (temp.substr(temp.find("=") + 1).find("true") || stoi(temp.substr(temp.find("=") + 1)) > 0)
					cal_dis_sc_contact = true;
			++argc;
		}
		if (temp.find("-opdb") != string::npos && temp.find("-f") == string::npos)
		{
			if (temp.find("=") == string::npos)
				opdb = true;
			else
				if (temp.substr(temp.find("=") + 1).find("false") || stoi(temp.substr(temp.find("=") + 1)) == 0)
					opdb = false;
				else if (temp.substr(temp.find("=") + 1).find("true") || stoi(temp.substr(temp.find("=") + 1)) > 0)
					opdb = true;
				else
				{
					out_pdb_name = temp.substr(temp.find("=") + 1);
					out_pdb_name.erase(remove(out_pdb_name.begin(), out_pdb_name.end(), ' '), out_pdb_name.end());
				}
			if (temp.find("-opdbsep") != string::npos) flag_sep_opdb = true;
			++argc;
		}

		if (temp.find("-new_al") != string::npos && temp.find("-f") == string::npos)
		{
			new_al_1pdb = true;
			++argc;
		}
		if (temp.find("-cat-result") != string::npos && temp.find("-f") == string::npos)
		{
			flag_cat_result_multi_infile = true;
			++argc;
		}
		if (temp.find("-nocal") != string::npos && temp.find("-f") == string::npos)
		{
			flag_nocal = true;
			flag_1_pdb = true;
			++argc;
		}
		if (temp.find("-pqr") != string::npos && temp.find("-f") == string::npos)
		{
			flag_pqr = 1;
			++argc;
		}
		if (temp.find("-dipole-moment") != string::npos && temp.find("=") == string::npos && temp.find("-f") == string::npos)
		{
			flag_cal_dipole_moment = true;
			++argc;
		}
		if (temp.find("-act") != string::npos && temp.find("-f") == string::npos)
		{
			flag_all_side_chain_ct = 0;
			++argc;
		}
		if (temp.find("-interchain-energy") != string::npos && temp.find("-f") == string::npos)
		{
			flag_interchain_energy = true;
			self_int = true;
			++argc;
		}
		if (temp.find("-refpbc") != string::npos && temp.find("-f") == string::npos && temp.find("-refpbc")< temp.find("="))
		{
			solv_pbc_int = 1;
			chain_ref = 0;
			if (temp.find("=") != string::npos)
			{
				ref_pdb_name = temp.substr(temp.find("=") + 1);
				if (ref_pdb_name.length() <= 4 || ref_pdb_name.substr(ref_pdb_name.length() - 4) != ".pdb")
					ref_pdb_name += ".pdb";
				flag_pbc_ref = true;
			}
			++argc;
		}
		if (temp.find("-pbc") != string::npos && temp.find("-pbci")==string::npos && temp.find("-pbcclust")==string::npos
			&& temp.find("-pbc-envi")==string::npos && temp.find("-pbca")==string::npos && temp.find("-f") == string::npos)
		{
			solv_pbc = 1;
			++argc;
		}
		if (temp.find("-pbci") != string::npos && temp.find("-pbci") < temp.find("=") && temp.find("-f") == string::npos)
		{
			solv_pbc_int = 1;
			solv_pbc = 1;
			if (temp.find("-pbci2") != string::npos)
				solv_pbc_int2 = 1;
			if (temp.find("-pbci3") != string::npos)
				solv_pbc_int3 = 1;
			if (temp.find("-pbci4") != string::npos)
				pbc_mol2_aggressive = 1;
			if (temp.find("-pbci5") != string::npos)
				pbc_mol2_aggressive = 2;
			if (temp.find("=") != string::npos)
			{
				chain_ref = stoi(temp.substr(temp.find("=") + 1));
				if (chain_ref < 0)chain_ref = 0;
				solv_pbc = 0;
			}
			else
			{
				chain_ref = 0;
			}
			++argc;
		}
		if (temp.find("-pbca") != string::npos && temp.find("-f")==string::npos)
		{
			solv_pbc_int = 1;
			solv_pbc = 1;
			pbc_aggressive = true;
			pbc_mol2_aggressive = 2;
			if (temp.find("=") != string::npos)
			{
				chain_ref = stoi(temp.substr(temp.find("=") + 1));
				solv_pbc = 0;
			}
			++argc;
		}

		if (temp.find("-pbcm") != string::npos && temp.find("-f") == string::npos)
		{
			solv_pbc_int = 1;
			solv_pbc = 1;
			pbc_mem = true;
			if (temp.find("=") != string::npos)
			{
				chain_ref = stoi(temp.substr(temp.find("=") + 1));
				solv_pbc = 0;
			}
			++argc;
		}

		if (temp.find("-pbcclust") != string::npos && temp.find("-f") == string::npos)
		{
			solv_pbc_clust = true;
			if (temp.find("=") != string::npos)
			{
				istringstream is(temp.substr(temp.find("=") + 1));
				char t{};
				while (is >> t)
					pbc_clust_axis.emplace_back(t);
			}
			else
			{
				pbc_clust_axis = vector<char>{'x','y','z'};
			}
		}
		if (temp.find("-pbc-envi") != string::npos && temp.find("-f") == string::npos)
		{
			if (temp.find("=") != string::npos)
			{
				istringstream is(temp.substr(temp.find("=") + 1));
				string t{};
				while (is >> t)
				{
					if (t == "@")t = " ";
					pbc_envi_chain.emplace_back(t);
				}
				if (!pbc_envi_chain.empty()) pbc_envi = true;
			}
		}
		if (temp.find("-1pdb") != string::npos && temp.find("-f") == string::npos)
		{
			flag_1_pdb = true;
			++argc;
		}
		if (temp.find("-h") != string::npos && temp.find("-hb") == string::npos && temp.find("-f") == string::npos)
		{
			help = true;
			help_print();
			++argc;
		}
		if (temp.find("-v") != string::npos && temp.find("-f") == string::npos && temp.find("-volume")==string::npos)
		{
			cout << "Build: " << get_build_number() << '\n';
			++argc;
			return 2;
		}
		if (temp.find("-energy-res") != string::npos && temp.find("-f") == string::npos)
		{
			print_energy_res = true;
			energy_map = 1;
			++argc;
		}
		if (temp.find("-force-res") != string::npos && temp.find("-f") == string::npos)
		{
			print_force_res = true;
			force_map = 1;
			++argc;
		}
		if (temp.find("-all-subfolders") != string::npos && temp.find("-f") == string::npos)
		{
			flag_all_folders = true;
			flag_separation_results = true;
			++argc;
		}
		if (temp.find("-cal-force") != string::npos && temp.find("-f") > temp.find("-cal-force"))
		{
			flag_cal_force = true;
			++argc;
		}
		if (temp.find("-cache") != string::npos && temp.find("-f") == string::npos)
		{
			flag_cache = 1;
			size_t pos = temp.find("=");
			if (pos!= string::npos)
			if (temp.substr(pos + 1).find("no") != string::npos)
				flag_cache = 2;
			++argc;
		}
		if (temp.find("-use-CI") != string::npos && temp.find("-f") == string::npos)
		{
			flag_use_CI = true;
			++argc;
		}
		if (temp.find("-useForceOptionforEnergy") != string::npos && temp.find("-f") == string::npos)
		{
			flag_energy_same_force = true;
			++argc;
		}
		if (temp.find("-max-force") != string::npos && temp.find("-f") == string::npos)
		{
			flag_maxforce = true;
			++argc;
		}
		if (temp.find("-split") != string::npos && temp.find("-f") == string::npos)
		{
			flag_separation_results = true;
			if (temp.find("=") != string::npos)
				if (temp.substr(temp.find("=") + 1).find("false") != string::npos || temp.substr(temp.find("=") + 1).find("0") != string::npos)
					flag_separation_results = false;
			++argc;
		}
		if (temp.find("-entropy") != string::npos&& temp.find("-entropy-quasi") == string::npos
			&& temp.find("-f") == string::npos)
			entropy = true;
		if (temp.find("-side_chain_contact") != string::npos && temp.find("-noprint-side_chain_contact")==string::npos)
		{
			flag_sidechain = true;
			if (temp.find("=") != string::npos && temp.substr(temp.find("=")+1).find("full-res")!=string::npos)
			{
				flag_sidechain = false;
			}
			++argc;
		}
		if (temp.find("-noprint-side_chain_contact") != string::npos)
		{
			++argc;
			flag_print_sidechain = false;
		}
		if (temp.length()>=12 && temp.find("-com-energy") != string::npos && temp.find("=")>temp.find("-com-energy")
			&& temp.find("-f") == string::npos)
		{
			++argc;
			if (temp.find("=") == string::npos)
			{
				com_dis_energy_2mol = 4; continue;
			}

			string temp_axis{ temp.substr(temp.find("=") + 1) };
			temp_axis.erase(remove(temp_axis.begin(), temp_axis.end(), ' '), temp_axis.end());
			if (temp_axis == "x")
				com_dis_energy_2mol = 1;
			else if (temp_axis == "y")
				com_dis_energy_2mol = 2;
			else if (temp_axis == "z")
				com_dis_energy_2mol = 3;
			else if (temp_axis == "all")
				com_dis_energy_2mol = 0;
			else com_dis_energy_2mol = 4;
		}
		if (temp.length() >= 8 && temp.find("-com-dis") != string::npos 
			&& temp.find("=") > temp.find("-com-dis") && temp.find("-f") == string::npos)
		{
			++argc;

			if (temp.find("=") == string::npos)
			{
				com_dis_2mol = 4; continue;
			}

			string temp_axis{ temp.substr(temp.find("=") + 1) };
			temp_axis.erase(remove(temp_axis.begin(), temp_axis.end(), ' '), temp_axis.end());
			if (temp_axis == "x")
				com_dis_2mol = 1;
			else if (temp_axis == "y")
				com_dis_2mol = 2;
			else if (temp_axis == "z")
				com_dis_2mol = 3;
			else
				com_dis_2mol = 4;
		}
		if (temp.length() >= 7 && temp.find("-resdis2com") != string::npos && temp.find("-f")==string::npos)
		{
			++argc;
			flag_cal_dis_res_2com = true;
		}
		if (temp.length() >= 3 && temp.find("-b") != string::npos
			&& temp.find("-f") == string::npos)
		{
			begin_time = stod(temp.substr(temp.find("=") + 1));
			flag_begin_set = true;
			++argc;
		}
		if (temp.length() >= 3 && temp.find("-e") != string::npos && temp.substr(temp.find("-e") + 1).find("=") != string::npos && temp.find("-emap") == string::npos && temp.find("-new-emap") == string::npos
			&& temp.find("-entropy") == string::npos && temp.find("-emap-level") == string::npos && temp.find("-eu") == string::npos
			&& temp.find("-energy-res") == string::npos && temp.find("-pbc-envi")==string::npos && temp.find("-interchain-energy")==string::npos 
			&& temp.find("-elec")==string::npos && temp.find("-ep")==string::npos)
		{
			end_time = stod(temp.substr(temp.find("=") + 1));
			flag_end_set = true;
			++argc;
		}
		if (temp.length() >= 4 && temp.find("-dt") != string::npos
			&& temp.find("-f") == string::npos)
		{
			dt = stoi(temp.substr(temp.find("=") + 1));
			++argc;
		}
		if (temp.length() > 4 && temp.find("-eu") != string::npos)
		{
			eu = temp.substr(temp.find("=") + 1);
			++argc;
		}
		if (temp.length() > 5 && temp.find("-fu") != string::npos
			&& temp.find("-f")==string::npos)
		{
			fu = temp.substr(temp.find("=") + 1);
			++argc;
		}
		if (temp.length() > 10 && temp.find("-cal-volume") != string::npos && temp.find("-f") == string::npos)
		{
			flag_cal_vol = true;
			++argc;
		}
		if (temp.length() > 10 && temp.find("-cal-mindis") != string::npos && temp.find("-f") == string::npos)
		{
			flag_cal_mindis = true;
			++argc;
		}
		if (temp.length() >= 5 && temp.find("-mol2") != string::npos
			&& temp.find("-f") == string::npos)
		{
			begin_mol_2 = stoi(temp.substr(temp.find("=") + 1));
			++argc;
		}
		if (temp.length() > 13 && temp.find("-elec-screening") != string::npos && temp.find("-f") == string::npos)
		{
			flag_elec_screening = true;
			if (temp.find("=") != string::npos && temp.find("false") != string::npos && temp.find("false") > temp.find("="))
			{
				flag_elec_screening = false;
				++argc;
			}
		}
		if (temp.length() >= 11 && temp.find("-epsilon_in") != string::npos && temp.find("-f") == string::npos)
		{
			if (temp.find("=") != string::npos && temp.find("=") > temp.find("-epsilon_in"))
			{
				epsilon_in = stod(temp.substr(temp.find("=") + 1));
				++argc;
			}
		}
		if (temp.length() >= 17 && temp.find("-charged_atom_dis") != string::npos && temp.find("-f") == string::npos)
		{
			flag_charged_atom_distance_dis = true;
			++argc;
			if (temp.find("=") != string::npos)
				number_charged_distances_bin = stoi(temp.substr(temp.find("=") + 1));
		}
		if(temp.length() > 22 && temp.find("-charged_distances_max") != string::npos)
		{
			charged_distances_max = stod(temp.substr(temp.find("=") + 1));
			flag_charged_atom_distance_dis = true;
			++argc;
			
		}
		if (temp.length() > 22 && temp.find("-charged_distances_min") != string::npos)
		{
			charged_distances_min = stod(temp.substr(temp.find("=") + 1));
			flag_charged_atom_distance_dis = true;
			++argc;

		}
		if (temp.length() >= 20 && temp.find("-angle_dipole_efield") != string::npos && temp.find("-f") == string::npos)
		{
			flag_angle_dipole_efield = true;
			++argc;
			if (temp.find("=") != string::npos)
				number_bin_angle_dipole_efield = stoi(temp.substr(temp.find("=") + 1));
		}
		if (temp.length() >= 12 && temp.find("-mol-1-type") != string::npos)
		{
			size_t pos{ temp.find("=") + 1 };
			if (pos < temp.length())
			{
				if (temp.substr(pos).find("protein") != string::npos)
					mol_1_type = "protein";
				else if (temp.substr(pos).find("complex") != string::npos)
					mol_1_type = "complex";
				else if (temp.substr(pos).find("lig") != string::npos || temp.substr(pos).find("ligand") != string::npos)
					mol_1_type = "lig";
				else
					mol_1_type = "unknown";
				flag_set_type_1 = true;
			}
			
		}

		if (temp.length() >= 12 && temp.find("-mol-2-type") != string::npos && temp.find("-f") == string::npos)
		{
			size_t pos{ temp.find("=") + 1 };
			if (pos < temp.length())
			{
				if (temp.substr(pos).find("protein") != string::npos)
					mol_2_type = "protein";
				else if (temp.substr(pos).find("complex") != string::npos)
					mol_2_type = "complex";
				else if (temp.substr(pos).find("lig") != string::npos || temp.substr(pos).find("ligand") != string::npos)
					mol_2_type = "lig";
				else if (temp.substr(pos).find("protein-user") != string::npos)
					mol_2_type = "protein-user";
				else
					mol_2_type = "unknown";

				flag_set_type_2 = true;
			}

		}

		if (temp.find("-new-emap") != string::npos && temp.find("-f") == string::npos)
		{
			flag_new_emap = true;
			++argc;
		}
		if (temp.find("-emap") != string::npos && temp.find("-emap-level")==string::npos && temp.find("-new-emap")==string::npos
			&& temp.find("-f") == string::npos)
		{
			energy_map = 1;
			if (temp.length() >= 7 && temp.find("=") != string::npos && temp.substr(temp.find("=")).find("-") != string::npos)
			{
				energy_map = 3;
				string::size_type pos_e{ temp.find("=") }, pos_mi{ temp.substr(temp.find("=")).find("-") };
				if (pos_e < pos_mi
					&& (temp.substr(pos_e + 1, pos_mi -1)=="residue"|| temp.substr(pos_e + 1, pos_mi -1) == "atom")
					&&(temp.substr(pos_e + pos_mi + 1) =="residue"|| temp.substr(pos_e + pos_mi + 1)=="atom"))
				{
					emap_type = temp.substr(pos_e + 1);
				}
			}
			if (temp.length() >= 7 && temp.find("=") != string::npos && temp.substr(temp.find("=")).find("-") == string::npos
				&& temp.find("maptype")!=string::npos && temp.find("maptype")>temp.find("emap"))
			{
				energy_map = 4;
				emap_level = "all";
			}
			++argc;
		}
		if (temp.find("-emap-level") != string::npos && temp.find("-f") == string::npos)
		{
			energy_map = 1;
			if (temp.length() >= 12 && temp.find("=") != string::npos && temp.find("=") > temp.find("-emap-level"))
			{
				emap_level = temp.substr(temp.find("=") + 1);
				emap_level.erase(remove(emap_level.begin(), emap_level.end(), ' '), emap_level.end());
				if (emap_level != "per-res" && emap_level != "all")
				{
					emap_level = "per-res";
					cout << "Your emap-level parameter is invalid, it must be per-res or all. This options is changed to per-res\n";
				}
			}
			++argc;
		}
		if (temp.find("-fmap") != string::npos && temp.find("-fmap-level") == string::npos && temp.find("-f") == string::npos)
		{
			force_map = 1;
			if (temp.length() >= 7 && temp.find("=") != string::npos && temp.find("=") > temp.find("-fmap") 
				&& temp.substr(temp.find("=")).find("-") != string::npos)
			{
				force_map = 3;
				string::size_type pos_e{ temp.find("=") }, pos_mi{ temp.substr(temp.find("=")).find("-") };
				if (pos_e < pos_mi
					&& (temp.substr(pos_e + 1, pos_mi -1) == "residue" || temp.substr(pos_e + 1, pos_mi -1) == "atom")
					&& (temp.substr(pos_e + pos_mi + 1) == "residue" || temp.substr(pos_e + pos_mi + 1) == "atom"))
				{
					fmap_type = temp.substr(pos_e + 1);
				}
			}
			++argc;
		}
		if (temp.find("-fmap-level") != string::npos && temp.find("-f") == string::npos)
		{
			force_map = 1;
			if (temp.length() >= 12 && temp.find("=") != string::npos && temp.find("=")>temp.find("-fmap-level"))
			{
				fmap_level = temp.substr(temp.find("=") + 1);
				fmap_level.erase(remove(fmap_level.begin(), fmap_level.end(), ' '), fmap_level.end());
				if (fmap_level != "per-res" && fmap_level != "all")
				{
					fmap_level = "per-res";
					cout << "Your fmap-level parameter is invalid, it must be per-res or all\n";
				}
			}
			++argc;
		}
		if (temp.length() > 6 && temp.find("-topo") != string::npos && temp.find("-topo2")==string::npos && temp.find("-f") == string::npos)
		{
			string str1{ temp.substr(temp.find("=")+1) };
			str1.erase(remove(str1.begin(), str1.end(), ' '), str1.end());
			topo_name = str1;
			if (check_topo(topo_name, 1) == true)
				flag_topo_1 = true;
			++argc;
		}
		if (temp.length() > 7 && temp.find("-topo2") != string::npos)
		{
			string str1{ temp.substr(temp.find("=") + 1) };
			str1.erase(remove(str1.begin(), str1.end(), ' '), str1.end());
			topo_name_2 = str1;
			if (check_topo(topo_name_2, 2) == true)
				flag_topo_2 = true;
			++argc;
		}
		if (temp.length() > 12 && temp.find("-common_topo") != string::npos)
		{
			flag_use_common_topo = true;
			if (flag_topo_1 == false && flag_topo_2 == false)
			{
				if (temp.find("=") != string::npos && temp.find("=") > temp.find("-common_topo"))
				{
					string name{ temp.substr(temp.find("=") + 1) };
					name.erase(remove(name.begin(), name.end(), ' '), name.end());
					topo_name = name; topo_name_2 = name;
				}
				else
				{
					cout << "-common_topo option has no file name. Will use topol.top as file name\n";
					topo_name = "topol.top"; topo_name_2 = "topol.top";
				}
				flag_topo_1 = true; flag_topo_2 = true;
				mol_1_type = "protein-user"; mol_2_type = "protein-user";
			}
			++argc;
		}
		if (temp.length() >= 3 && temp.find("-f") != string::npos && temp.find("-ff") == string::npos
			&& temp.find("-force") == string::npos && temp.find("-fu") == string::npos
			&& temp.find("-fmap") == string::npos)
		{
			string str1{ temp.substr(temp.find("=") + 1) };
			str1.erase(remove(str1.begin(), str1.end(), ' '), str1.end());
			pdb_name = str1;
			flag_1_pdb = 1;
			++argc;
		}
		if (temp.length() >= 3 && temp.find("-n") != string::npos && temp.find("-nt") == string::npos && temp.find("path-n") == string::npos
			&& temp.find("-new_al") == string::npos&& temp.find("-nocal") == string::npos && temp.find("-new-emap")==string::npos
			&& temp.find("-noprint-side_chain_contact")==string::npos && temp.substr(temp.find("-e") + 1).find("=") !=string::npos
			&& temp.find("-num")==string::npos)
		{
			string str1{ temp.substr(temp.find("=") + 1) };
			str1.erase(remove(str1.begin(), str1.end(), ' '), str1.end());
			index_name = str1; flag_index = true;
			++argc;
		}
		if (temp.length() >= 4 && temp.find("-ff") != string::npos)
		{
			string str1{ temp.substr(temp.find("=")+1) };
			str1.erase(remove(str1.begin(), str1.end(), ' '), str1.end());
			force_f = str1;
			transform(force_f.begin(), force_f.end(), force_f.begin(), ::toupper);
			++argc;
		}
		if (temp.length() >= 5 && temp.find("-g") != string::npos && temp.find("-sg")==string::npos)
		{
			string str1{ temp.substr(temp.find("=")+1) };
			str1.erase(remove(str1.begin(), str1.end(), ' '), str1.end());
			if (str1.find("&") != string::npos) {
				group_1 = stoi(str1.substr(0, str1.find("&")));
				group_2 = stoi(str1.substr(str1.find("&") + 1));
			}
			else if (str1.find("-") != string::npos) {
				group_1 = stoi(str1.substr(0, str1.find("-")));
				group_2 = stoi(str1.substr(str1.find("-") + 1));
			}
			else
				cout << "Your parameter -g is invalid!\n";

			if (group_1 == group_2) self_int = true;
			++argc;
		}
		if (temp.length() >= 6 && temp.find("-sg") != string::npos && temp.find("=")!=string::npos)
		{
			temp = temp.substr(temp.find("=") + 1);
			const auto found = temp.find("-");
			if (found == string::npos)
			{
				cout << "Your argument -sg is invalid!\n"
					<< "You should use - to separate groups of two molecule!\n";
			}
			else
			{
				string groups_1{ temp.substr(0, temp.find("-")) },
					groups_2{ temp.substr(temp.find("-") + 1) };

				regex reg{ "(&?)([[:blank:]]*)(\\d*)([[:blank:]]*)(&*)" };

				sregex_iterator pos(groups_1.begin(), groups_1.end(), reg), pos_2(groups_2.cbegin(), groups_2.cend(), reg), end;
				
				for (; pos != end; ++pos)
				{
					if(!pos->str(3).empty() && isdigit(pos->str(3)[0]))
					sub_group_1.emplace_back(stoi(pos->str(3)));
				}

				for (; pos_2 != end; ++pos_2)
				{
					if(!pos_2->str(3).empty() && isdigit(pos_2->str(3)[0]))
					sub_group_2.emplace_back(stoi(pos_2->str(3)));
				}
			}
			flag_sub_groups = true;
			++argc;
			/*fstream test("test_sg.dat", ios::out);
			for (size_t i{}; i < sub_group_1.size(); ++i)test << sub_group_1[i] << " ";
			for (size_t i{}; i < sub_group_2.size(); ++i)test << sub_group_2[i] << " ";*/
		}
		if (temp.length() >= 4 && temp.find("-nt") != string::npos)
		{
			int nt_temp{ stoi(temp.substr(temp.find("=") + 1)) };
			if (nt_temp < 0)nt = 0;
			else
				nt = nt_temp;
			++argc;
		}
		if (temp.length() >= 5 && temp.find("-dsc") != string::npos)
		{
			d_standard = stod(temp.substr(temp.find("=") + 1));
			++argc;
		}
		
		if (temp.length() >= 8 && temp.find("-denfo") != string::npos && temp.find("-f") == string::npos)
		{
			distance_energy_cutoff = stod(temp.substr(temp.find("=") + 1));
			distance_energy_cutoff *= 0.1;
			++argc;
		}
		if (temp.length() >= 7 && temp.find("-path-a") != string::npos)
		{
			string str1{ temp.substr(temp.find("=") + 1) };
			str1.erase(remove(str1.begin(), str1.end(), ' '), str1.end());
			if (str1[str1.length() - 1] != '/')
				str1 += '/';
			path_amino_acid = str1;
			++argc;
		}
		if (temp.length() >= 7 && temp.find("-path-n") != string::npos)
		{
			string str1{ temp.substr(temp.find("=") + 1) };
			str1.erase(remove(str1.begin(), str1.end(), ' '), str1.end());
			if (str1[str1.length() - 1] != '/')
				str1 += '/';
			path_nonbonded = str1;
			++argc;
		}
		if (temp.length() >= 3 && temp.find("-hb") != string::npos)
		{
			flag_h_bond = 1;
			string::size_type pos{ temp.find("=") };
			if (pos != string::npos && temp.substr(pos + 1).find("side-chain") != string::npos)
				flag_sidechain_hbond = true;
			++argc;
		}
		if (temp.length() >= 5 && temp.find("-ahb") != string::npos)
		{
			angle_h = stod(temp.substr(temp.find("=") + 1));
			++argc;
		}
		if (temp.length() >= 5 && temp.find("-dhb") != string::npos)
		{
			dis_h = stod(temp.substr(temp.find("=") + 1));
			++argc;
		}
		if (temp.length() >= 8 && temp.find("-entropy-quasi") != string::npos)
		{
			istringstream is(temp.substr(temp.find("=") + 1));
			int g1{-1}, g2{-1}, g3{-1};
			if (!(is >> g1 >> g2 >> g3))
				cout << "Cannot parse entropy-quasi option!\n";
			if (g1 != -1 && g2 != -1 && g3 != -1)
			{
				entropy_quasi = true;
				group_entropy_com = g1;
				group_entropy_rec = g2;
				group_entropy_lg = g3;
			}
			else
			{
				cout << "Cannot parse entropy-quasi option!\n";
			}
			++argc;
		}
		if (temp.length() >= 16 && temp.find("-pattern-input") != string::npos && temp.find("=") !=string::npos)
		{
			if (temp.find("*") != string::npos)
			{
				flag_pattern_input = true;
				flag_1_pdb = true;
				istringstream is(temp.substr(temp.find("=") + 1));
				is >> input_pattern >> b_input >> e_input >> step_input;
				if (e_input != 100000000)
					flag_range_input = true;
				
				if (b_input > e_input)
				{
					int t{ b_input };
					b_input = e_input;
					e_input = t;
				}
				if (step_input <= 0)
				{
					cout << "Step input must be positive!\n";
					step_input = 1;
				}
			}
			else
				cout << "Can not parse pattern of input files!";
			++argc;
		}
		if (temp.length() >= temp.find("=") + 2 && temp.find("-i") != string::npos && (temp.find("-f")==string::npos || temp.find("-f") > temp.find("="))
			&& temp.find("=") > temp.find("-i") 
			&& temp.find("-input")==string::npos && temp.find("=")!=string::npos && temp.find("-ildn") == string::npos
			&& temp.find("-int")==string::npos && temp.find("-ite")==string::npos)
		{
			list_filename = temp.substr(temp.find("=") + 1);
			fstream test_listFile(list_filename, ios::in);
			if (test_listFile.is_open())
				flag_listfile = true;
			else
				flag_listfile = false;
			test_listFile.close();
			++argc;
		}
		if (temp.length() > temp.find("=") + 1 && temp.find("-f") == string::npos && temp.find("-water-type") != string::npos)
		{
			if (temp.find("=") == string::npos) { water_type = "surface"; continue; }
			water_type = temp.substr(temp.find("=") + 1);
			if (!(water_type == "surface" || water_type == "sphere"))
			{
				cout << "WARING!!! -water-type option is incorrect!\n";
				water_type = "surface";
			}
			++argc;
		}
		if (temp.length() > temp.find("=") + 1 && temp.find("-f") == string::npos && temp.find("-dis-water") != string::npos)
		{
			d_crit_water = stod(temp.substr(temp.find("=") + 1));
			++argc;
		}
		if (temp.length() > temp.find("=") + 1 && temp.find("-contact-map") != string::npos && temp.find("-f") == string::npos)
		{
			string str1 = temp.substr(temp.find("=") + 1);
			str1.erase(remove(str1.begin(), str1.end(), ' '), str1.end());
			if (str1 == "no")
				flag_cal_contact = false;
			++argc;
		}
		if (temp.length() > 12 && temp.find("-all-respair") != string::npos)
		{
			flag_all_respair = true;
			++argc;
		}

		//pulling

		if (temp.length() > 11 && temp.find("-pull-direct") != string::npos)
		{
			flag_find_pulling_direct = true;
			++argc;
		}
		if (temp.length() > 15 && temp.find("-de-diff-weight") != string::npos)
		{
			de_diff_weight = stod(temp.substr(temp.find("=") + 1));
			++argc;
		}
		if (temp.length() > 14 && temp.find("-de-cross-prob") != string::npos)
		{
			de_cross_prob = stod(temp.substr(temp.find("=") + 1));
			++argc;
		}
		if (temp.length() > 5 && temp.find("-de-F") != string::npos)
		{
			de_F = stod(temp.substr(temp.find("=") + 1));
			++argc;
		}
		if (temp.length() > 15 && temp.find("-de-scan-radius") != string::npos)
		{
			de_scan_radius = stod(temp.substr(temp.find("=") + 1));
			++argc;
		}
		if (temp.length() > 18 && temp.find("-de-min-scan-radius") != string::npos)
		{
			de_min_scan_radius = stod(temp.substr(temp.find("=") + 1));
			++argc;
		}
		if (temp.length() > 18 && temp.find("-de-score-function") != string::npos)
		{
			de_flag_score_function = stoi(temp.substr(temp.find("=") + 1));
			++argc;
		}
		if (temp.length() > 11 && temp.find("-de-pop-num") != string::npos)
		{
			de_pop_num = stoi(temp.substr(temp.find("=") + 1));
			++argc;
		}
		if (temp.length() > 11 && temp.find("-de-ite-num") != string::npos)
		{
			de_ite = stoi(temp.substr(temp.find("=") + 1));
			++argc;
		}
		if (temp.length() > 13 && temp.find("-de-algorithm") != string::npos)
		{
			de_algorithm = stoi(temp.substr(temp.find("=") + 1));
			++argc;
		}
		if (temp.length() > 9 && temp.find("-de-angle") != string::npos && temp.find("-de-angle") < temp.find("="))
		{
			de_restrain_angle = stod(temp.substr(temp.find("=") + 1));
			++argc;
		}
		if (temp.length() > 10 && temp.find("-de-auto-term") != string::npos)
		{
			flag_de_auto_term = true;
			++argc;
		}
	}
	in_file.close(); 
	return 0;
}
int arguments::parse_arguments_interactive(int argc_in)
{
	int flag;
	string force_f;
	do
	{
		cout << "Please choose the force field that you want to use:" << '\n';
		cout << "1. AMBER99SB" << '\n';
		cout << "2. AMBER99SB-ILDN" << '\n';
		cout << "3. AMBER99SB-BSC1" << '\n';
		cout << "4. AMBER99SB-DISP" << '\n';
		cout << "5. AMBER14SB" << '\n';
		cout << "6. GROMOS43A1" << '\n';
		cout << "7. GROMOS53A6" << '\n';
		cout << "8. GROMOS54A7" << '\n';
		cout << "9. CHARMM27" << '\n';
		cout << "10. CHARMM36" << '\n';
		cout << "11. OPLSAA" << '\n';
		getline(cin, force_f);
		istringstream st(force_f);
		st >> flag;
	} while (flag >7 || flag<0);
	if (flag == 1)
	{
		force_f = "AMBER99SB";

	}
	else if (flag == 2)
	{
		force_f = "AMBER99SB-ILDN";
	}
	else if (flag == 3)
	{
		force_f = "AMBER99SB-BSC1";
	}
	else if (flag == 4)
	{
		force_f = "AMBER99SB-DISP";
	}
	else if (flag == 5)
	{
		force_f = "AMBER14SB";
	}
	else if (flag == 6)
	{
		force_f = "GROMOS43A1";

	}
	else if (flag == 7)
	{
		force_f = "GROMOS53A6";
	}
	else if (flag == 8)
	{
		force_f = "GROMOS54A7";
	}
	else if (flag == 9)
	{
		force_f = "CHARMM27";
	}
	else if (flag == 10)
	{
		force_f = "CHARMM36";
	}
	else
		force_f = "OPLSAA";

	do
	{
		cout << "Do you want to create energy map? (yes = 1, no = 2, or custom = 3)";
		if (!(cin >> energy_map))
		{
			cin.clear(); cin.ignore(10000, '\n');
		}
	} while (energy_map != 1 && energy_map != 2 && energy_map != 3);
	if (energy_map == 3)
	{
		cout << "Please specify type of energy map calculation: (residue-atom | atom-atom | residue-residue | atom-residue)";
		string temp{};getline(cin, temp);
		string::size_type pos{ temp.find("-") };
		if (pos != string::npos
			&& (temp.substr(0, pos) == "residue" || temp.substr(0, pos) == "atom")
			&& (temp.substr(pos + 1) == "residue" || temp.substr(pos + 1) == "atom"))
		{
			emap_type = temp;
		}
		else
		{
			cout << "Cannot parse! will use default value (residue-atom)\n";
		}
	}

	do
	{
		cout << "Do you want to take only minimun sc-ct distance? (Y/N=1/2)";
		if (!(cin >> flag_all_side_chain_ct))
		{
			cin.clear(); cin.ignore(10000, '\n');
		}
	} while (flag_all_side_chain_ct != 1 && flag_all_side_chain_ct != 2);
	do
	{
		cout << "Do you want to creat .PQR file? (Y/N=1/2)";
		cin >> flag_pqr;
	} while (flag_pqr != 1 && flag_pqr != 2);

	cout << "Do you want to use topo file (no if you don't want to use topo): ";
	cin >> topo_name;
	if (check_topo(topo_name, 1) == true)
		flag_topo_1 = true;

	cout << "Number of threads: ";
	cin >> nt;
	return 0;
}
bool arguments::check_topo(string NameofTopoFile, int Topo)
{
	if (NameofTopoFile == "no")
		return false;
	else 
	{
		fstream input_topo_protein(NameofTopoFile, ios::in);
		bool topo_is_open{ input_topo_protein.is_open() };
		input_topo_protein.close();
		if (topo_is_open == true)
			return true;
		else
		{
			while (!topo_is_open)
			{
				cout << "Cannot open " << NameofTopoFile << ". Please type your lig topo file (no if you don't want to use topo): ";
				cin >> NameofTopoFile;
				cout << NameofTopoFile << "\n";
				topo_name = NameofTopoFile;
				if (NameofTopoFile == "no") 
					return false;
				else
				{
					input_topo_protein.open(NameofTopoFile, ios::in);
					if (input_topo_protein.is_open()) { input_topo_protein.close(); return true; }
					else
					{
						input_topo_protein.close();	continue;
					}
				}
			}
			return false;
		}
	}
	
}