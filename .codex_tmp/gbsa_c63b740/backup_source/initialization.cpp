#include<iostream>
#include<fstream>
#include<string>
#include"arguments.h"
#include"MOL.h"
#include"groups.h"
#include"functions.h"

int helper_parse_groups(groups_cal& Groups_from_index_file, arguments& Arguments)
{
	return Groups_from_index_file.get_parse_index_flag() == false ? Groups_from_index_file.index(Arguments) : 0;
}

int parse_input_pdb(arguments& Arguments, fstream& log, string name_pro, string name_lig, string force_f, int input_name, vector<string>& input_line
	, bool topo_2_is_open, bool topo_is_open, Multi_MOL& pro_in, Multi_MOL& lig, groups_cal& Groups_from_index_file
, para& parameter, amino_acid& amine, info& output_info)
{
	int flag_groups{ helper_parse_groups(Groups_from_index_file, Arguments) };

	fstream test(name_pro, ios::in);
	if (!test.is_open())
	{
		if ((Arguments.flag_1_pdb == true || Arguments.argc == 1) && input_name > 0 
			&& Arguments.flag_listfile == false && Arguments.flag_pattern_input == false && Arguments.flag_all_folders == false)
			return 1;

		if (Arguments.flag_1_pdb == false)
		{
			cout << "Cannot open " << name_pro << '\n';
			log << "Cannot open " << name_pro << '\n';
			return 1;
		}
		else if ((Arguments.flag_listfile == true || Arguments.flag_all_folders == true ) && Arguments.flag_1_pdb == true && !input_line.empty() && input_name < input_line.size())
		{
			if (input_name == 0 )
			{
				if (Arguments.flag_listfile == true)
				{
					cout << "Will take .pdb files from " << Arguments.list_filename << " ...\n";
					log << "Will take .pdb files from " << Arguments.list_filename << " ...\n";
				}
				else if (Arguments.flag_all_folders == true)
				{
					cout << "Will take .pdb files from folder " << fs::current_path() << "\n";
					log << "Will take .pdb files from folder " << fs::current_path() << "\n";
				}
			}
			Arguments.pdb_name = input_line[input_name];
			if ((Arguments.pdb_name.length()>4 && Arguments.pdb_name.substr(Arguments.pdb_name.length() - 4) != ".pdb")
				|| Arguments.pdb_name.length()<=4)
				Arguments.pdb_name += ".pdb";

			cout << Arguments.pdb_name << "...\n";
			log << Arguments.pdb_name << "...\n";
		}
		else if (Arguments.flag_pattern_input == true && Arguments.flag_1_pdb == true)
		{
			if (input_name == Arguments.b_input)
			{
				cout << "Will take .pdb files from pattern " << Arguments.input_pattern << " ...\n";
				log << "Will take .pdb files from pattern " << Arguments.input_pattern << " ...\n";
			}

			if (Arguments.flag_range_input == true)
				Arguments.pdb_name = Arguments.input_pattern.substr(0, Arguments.input_pattern.find("*")) + input_line[input_name - Arguments.b_input]
				+ Arguments.input_pattern.substr(Arguments.input_pattern.find("*") + 1);
			else
				Arguments.pdb_name = Arguments.input_pattern.substr(0, Arguments.input_pattern.find("*")) + to_string(input_name)
				+ Arguments.input_pattern.substr(Arguments.input_pattern.find("*") + 1);

			Arguments.add_input_file_type();

			cout << Arguments.pdb_name << "...\n";
			log << Arguments.pdb_name << "...\n";
		}
		else
		{
			cout << "Analysing " << Arguments.pdb_name << " file....\n";
			log << "Analysing " << Arguments.pdb_name << " file....\n";
		}

		if (Arguments.flag_1_pdb == true)
		{
			log << "Using one file .pdb \nAnalysing "<< Arguments.pdb_name <<" file...\n";
			int flag_pre_proc{ 1 };
			if (topo_2_is_open == true)
				flag_pre_proc = 2;
			if (Arguments.flag_use_common_topo == true)
				flag_pre_proc = 3;
			if (Arguments.topo_name != "no" && (Arguments.file_par == false || topo_is_open))
			{
				const int r_value_pdb{ pre_proc_pdb(force_f, pro_in, lig, Groups_from_index_file, flag_pre_proc, Arguments, amine, output_info)};
				
				if (r_value_pdb == 1)
				{
					log << "Error occured! Your group is not contained in pdb file\n";
					exit(EXIT_FAILURE);
				}
				else if (r_value_pdb == 2)
					return 1;
				else {}

				if (flag_pre_proc < 3)
				{
					lig.topo_lg(Arguments.topo_name, force_f, parameter);
					if (flag_pre_proc == 2)
						pro_in.topo_lg(Arguments.topo_name_2, force_f, parameter);
				}
				else
				{
					vector<topo_info_temp> info_topo;
					int topo_parse_result_1{ lig.topo_protein(Arguments.topo_name, force_f, parameter, info_topo) },
						topo_parse_result_2{ pro_in.topo_protein(Arguments.topo_name_2, force_f, parameter, info_topo) };

					if (topo_parse_result_1 == 1 && topo_parse_result_2 == 1)
					{
						log << "Cannot parse topology files (" << Arguments.topo_name << " and " << Arguments.topo_name_2 << ")\n";
						return 1;
					}
				}
			}
			else
			{
				const int r_value_pdb{ pre_proc_pdb(force_f, pro_in, lig, Groups_from_index_file, 0, Arguments, amine, output_info) };
				
				if (r_value_pdb == 1)
				{
					log << "Error occured! Your group is not contained in pdb file\n";
					exit(EXIT_FAILURE);
				}
				else if (r_value_pdb == 2)
					return 1;
				else {}
			}
			if (flag_groups == 0)
			{
				log << "------------------------------------------------oOo-----------------------------------------------\n"
					<< "Interaction energy is calculated for group " << Groups_from_index_file.group_name(Groups_from_index_file.group_1_index()) << "[" << Groups_from_index_file.group_1_index() << "] "
					<< "(" << Groups_from_index_file.groups_atom_number(Groups_from_index_file.group_1_index()) <<" atoms) "
					<< " and group " << Groups_from_index_file.group_name(Groups_from_index_file.group_2_index()) << "[" << Groups_from_index_file.group_2_index() << "] "
					<< "("<< Groups_from_index_file.groups_atom_number(Groups_from_index_file.group_2_index()) << " atoms) "
					<< "\n-------------------------------------------------------------------------------------------------\n";
			}
			log << "Type of system: " << pro_in.get_mol_type() << " and " << lig.get_mol_type() << "\n";
		}
		else
			return 1;
	}
	else
	{
		cout << name_pro << " - " << name_lig << '\n';
		log << name_pro << " - " << name_lig << '\n' << "Analysing input files...\n";
		fstream test_2(name_lig, ios::in);
		if (!test_2.is_open()) { cout << "Cannot open " << name_lig << " file!\n"; log << "Cannot open " << name_lig << " file!\n"; }
		test_2.close();
		test.close();
		if (topo_is_open==true && topo_2_is_open == false && pro_in.m_mol.empty() && lig.m_mol.empty())
		{
			pro_in.take_coordinate(name_pro, force_f, Arguments, amine, output_info, 0);
			lig.take_topo_lg(name_lig, Arguments.topo_name, Arguments, force_f, parameter, amine, output_info, Arguments.mol_2_type);
		}
		else if (topo_is_open == true && topo_2_is_open == true && pro_in.m_mol.empty() && lig.m_mol.empty())
		{
			lig.take_topo_lg(name_lig, Arguments.topo_name, Arguments, force_f, parameter, amine, output_info, Arguments.mol_2_type);
			pro_in.take_topo_lg(name_pro, Arguments.topo_name_2, Arguments, force_f, parameter, amine, output_info, Arguments.mol_1_type);
		}
		else if (topo_is_open == false && topo_2_is_open == true && pro_in.m_mol.empty() && lig.m_mol.empty())
		{
			pro_in.take_topo_lg(name_pro, Arguments.topo_name_2, Arguments, force_f, parameter, amine, output_info, Arguments.mol_1_type);
			lig.take_coordinate(name_lig, force_f, Arguments, amine, output_info, 0);
		}
		else
		{
			if (pro_in.m_mol.empty() && lig.m_mol.empty())
			{
				pro_in.take_coordinate(name_pro, force_f, Arguments, amine, output_info, 0);
				lig.take_coordinate(name_lig, force_f, Arguments, amine, output_info, 0);
			}
		}
		if (pro_in.m_mol.empty() && lig.m_mol.empty())
			return 1;
		log << "Type of system: " << pro_in.get_mol_type() << " and " << lig.get_mol_type() << "\n";
		Arguments.setPdbname(name_lig.substr(name_lig.find("_")));
		output_info.setPdbName(Arguments);
	}
	return 0;
}