#include"arguments.h"
#include<thread>
#include<algorithm>
#include<cctype>
#include<fstream>
#include<iostream>
#include<filesystem>
namespace fs = std::filesystem;

bool arguments::parse_input_file_is_open()
{
	bool check_file_open{ false };
	fstream test_file(pdb_name, ios::in);
	if (test_file.is_open() || flag_all_folders == true)check_file_open = true;
	test_file.close();
	return check_file_open;
}

int arguments::recheck_and_change_pdb_input()
{
	add_input_file_type();

	if (flag_1_pdb == true && parse_input_file_is_open() == false && flag_pattern_input == false && flag_listfile == false)
	{
		bool check_file_open{ false };
		do
		{
			pdb_name = string();
			cout << "Please type your input " << get_input_file_type() << " file:\n";
			cin >> pdb_name;
			add_input_file_type();
			cout << "Your input: " << pdb_name << '\n';
			if (pdb_name.length()<1)return 1;
			flag_1_pdb = true;
			check_file_open = parse_input_file_is_open();
		} while (check_file_open == false);
	}
	else
		return 0;
	return 0;
}

string arguments::set_Buildnumber() {
	return "20260126";
}
string arguments::get_build_number() const
{
	return Buildnumber;
}

string arguments::get_pdbOutput() const {
	fs::path p_temp = pdb_name;
	return "out_" + p_temp.stem().string() + ".pdb";
}

int arguments::post_process()
{
	recheck_and_change_pdb_input();

	if (pdb_name.length()>4 && pdb_name.substr(pdb_name.length() - 4) != get_input_file_type())
	{
		cout << "Your input file " << pdb_name << " is invalid, input file must be " << get_input_file_type() << " file!\n";
		exit(EXIT_FAILURE);
	}
	if (force_f != "AMBER99SB"&& force_f != "AMBER99SB-ILDN" && force_f != "AMBER99SB-BSC1" && force_f != "AMBER99SB-DISP" && force_f != "AMBER14SB"
		&& force_f != "CHARMM27" && force_f != "CHARMM36"
		&& force_f != "GROMOS43A1" && force_f != "GROMOS53A6" && force_f != "GROMOS54A7"
		&& force_f != "OPLSAA")
	{
		cout << "Force field " << force_f << " is not supported!\n"
			<< "This tool supports AMBER99SB, AMBER99SB-ILDN, AMBER99SB-BSC1, AMBER99SB-DISP, AMBER14SB, GROMOS43A1, GROMOS53A6, GROMO54A7, CHARMM27, CHARMM36 and OPLSAA\nPlease check your option.\n";
		exit(EXIT_FAILURE);
	}
	if (begin_time < 0. || end_time < 0.)
	{
		cout << "Your begin time and end time must not be negative!\n";
		exit(EXIT_FAILURE);
	}
	if (dt < 0)
	{
		cout << "Time interval must not be negative!\n";
		exit(EXIT_FAILURE);
	}
	if (group_1 < 0 || group_2 < 0)
	{
		cout << "Groups for calculating energy are invalid! ";
		cout << "Your index group(s) must not be negative!\n";
		exit(EXIT_FAILURE);
	}
	if (entropy_quasi == true && (group_entropy_com < 0 || group_entropy_rec < 0 || group_entropy_lg < 0))
	{
		cout << "Entropy-quasi option is invalid! ";
		cout << "Your entropy-quasi group(s) must not be negative!\n";
		exit(EXIT_FAILURE);
	}
	if (b_input < 0 || e_input < 0)
	{
		cout << "Begin/end input must be positive!\n";
		exit(EXIT_FAILURE);
	}
	if (b_input > e_input)
	{
		cout << "Begin input must be smaller than end!\n";
		exit(EXIT_FAILURE);
	}
	if (step_input <= 0)
	{
		cout << "Step input must be positive!\n";
		exit(EXIT_FAILURE);
	}
	if (nt == 0 || nt > std::thread::hardware_concurrency())
	{
		nt = std::thread::hardware_concurrency();
		cout << "The number of threads is invalid. Will use " << nt << " threads!\n";
	}

	if (flag_listfile == true)
	{
		fstream list(list_filename, ios::in);
		if (!list.is_open())
		{
			cout << "Cannot open " << list_filename << "!\n";
			flag_listfile = false;
			return 1;
		}
		list.close();
	}

	if (flag_energy_same_force == true)
	{
		emap_level = fmap_level;
		emap_type = fmap_type;
		print_energy_res = print_force_res;
		cout << "Note: Options for energy and force are same!\n";
	}

	if (solv_pbc_int == 1 && solv_pbc_clust == true)
	{
		cout << "WARNING!! Both intensive pbc solver and clustered pbc solver are turned on, will turn off clustered pbc solver!\n";
		solv_pbc_clust = false;
	}
	/*if (flag_cal_force == true && print_energy_res == false)
	{
		cout << "WARNING! -cal-force option is on, will activate -energy-res option\n";
		print_energy_res = true;
	}*/

	if (flag_maxforce == true && flag_cal_force == false)
	{
		cout << "Note: -max-force option is on, will activate -cal-force option!\n";
		flag_cal_force = true;
	}

	if (flag_topo_1 == false)
		mol_2_type = "protein";

	if (mol_2_type.find("protein") != string::npos && mol_1_type.find("protein") != string::npos)
		emap_type = "residue-residue";

	if (group_1 == group_2) self_int = true;

	string name{ " -f " + pdb_name };
	if (flag_1_pdb == false)
		name = "-f %rec_x" + get_input_file_type() + "%";
	else if (flag_pattern_input == true)
		name = "";
	else if (flag_listfile == true)
		name = " -i " + list_filename;
	else if (flag_all_folders == true)
		name = " -all-subfolders ";
	else {}

	if (command == "-")
		command = "-nt " + to_string(nt) + " -ff " + force_f + name + " -dsc " + to_string(d_standard) + " -act " + to_string(flag_all_side_chain_ct)
		+ " -mol2 " + to_string(begin_mol_2) + " -b " + to_string(begin_time) + " -e " + to_string(end_time) + " -dt " + to_string(dt) + " -g " + to_string(group_1) + "-" + to_string(group_2)
		+ " -entropy-quasi " + to_string(entropy_quasi) + " -entropy " + to_string(entropy) + " -pbc " + to_string(solv_pbc) + " -pbci " + to_string(solv_pbc_int);
	
	if (flag_h_bond == 1) { command += " -hb "; if (flag_sidechain_hbond == true) command += " side-chain "; }
	
	if (flag_listfile == true)
		command += " -i " + list_filename;
	if (flag_topo_1 == true)
		command += " -topo " + topo_name;
	if (flag_topo_2 == true)
		command += " -topo2 " + topo_name_2;
	if (flag_pattern_input == true)
		command += " -pattern-input " + input_pattern;
	if (flag_energy_same_force == true)
		command += " -useForceOptionforEnergy ";
	if (flag_cal_force == true)
		command += " -cal-force ";
	if (print_energy_res == true)
		command += " -energy-res ";
	if (print_force_res == true)
		command += " -force-res ";
	if (self_int == true)
		command += " -self-int ";
	if (flag_maxforce == true)
		command += " -max-force ";
	if (flag_new_emap == true)
		command += " -new-emap ";
	command += "-emap-type " + emap_type + " ";
	command += "-mol-1-type " + mol_1_type + " -mol-2-type " + mol_2_type + " ";

	if (flag_sub_groups == true)
	{
		command += " -sg ";
		for (size_t i{}; i < sub_group_1.size(); ++i)
			command += to_string(sub_group_1[i]) + " ";
		command += "- ";
		for (size_t i{}; i < sub_group_2.size(); ++i)
			command += to_string(sub_group_2[i]) + " ";
	}
	if (distance_energy_cutoff > 0.)
		command += " -denfo " + to_string(distance_energy_cutoff);

	transform(eu.begin(), eu.end(), eu.begin(), ::toupper);
	if (eu != "KJ/MOL" && eu != "KCAL/MOL")
	{
		cout << "Energy unit is invalid, only kcal/mol and kJ/mol are supported. Will use the default unit (kJ/mol)\n";
		eu = "kJ/mol";
	}
	else if (eu == "KJ/MOL")
		eu = "kJ/mol";
	else
	{
		eu = "kcal/mol";
		eu_par = 4.184;
	}

	transform(fu.begin(), fu.end(), fu.begin(), ::toupper);
	if (fu != "PN" && fu != "KJ/MOL/NM")
	{
		cout << "Force unit is invalid, only pN and kJ/mol/nm are supported. Will use the default unit (pN)\n";
		fu = "pN";
	}
	else if (fu == "PN")
	{
		fu = "pN";
	}
	else
	{
		fu = "kJ/mol/nm";
		fu_par = 1.;
	}

	transform(tu.begin(), tu.end(), tu.begin(), ::toupper);
	if (tu != "PS" && tu != "NS")
	{
		cout << "Time unit is invalid, only ps and ns are supported. Will use the default unit (ps)\n";
	}
	else if (tu == "PS") tu = "ps";
	else tu = "ns";

	if (tu == "ns") {
		begin_time *= 1000.; end_time *= 1000.; dt *= 1000;
	}

	size_t lnum{ command.length() / 120 };
	for (size_t i{}; i < lnum; ++i)
	{
		size_t end{ (i + 1) * 120 + i + 1 >= command.length()?command.length()-1: (i + 1) * 120 + i + 1 };
		command = command.substr(0, 120 * (i + 1) + i) + "\n" + command.substr(end);
	}

	if (flag_index == false)
	{
		fstream index_file("index.ndx", ios::in);
		if (index_file.is_open())flag_index = true;
	}

	if (flag_elec_screening == true)
	{
		command += " -elec-screening -epsilon " + to_string(epsilon_in) + " -salt_concentration " + to_string(salt_concentration);
	}
	return 0;
}