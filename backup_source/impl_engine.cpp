#include<string>
#include<thread>
#include<vector>
#include<fstream>
#include<iostream>
#include<iomanip>
#include"functions.h"
#include<cmath>
#include<cstdlib>
#include<algorithm>
#include<ctime>
#include<functional>
#include"differential_evolution.h"


#pragma warning(disable: 4996)
using namespace std;
using namespace chrono;
int engine(arguments& Arguments)
{
	Arguments.post_process(); 
	/************************************/
	folders paths;
	vector<string> list_paths;
	de_arg DE_arguments;
	if (Arguments.flag_find_pulling_direct == true)
		DE_arguments.parse_input(Arguments);
	paths.add_paths(Arguments, list_paths);
	paths.create_folders(Arguments); 

	int sc_flag{Arguments.flag_all_side_chain_ct};
	auto total_run_time{chrono::steady_clock::now()};
	chrono::milliseconds time_energymaptype{  };
	string force_f{ Arguments.force_f };

	if (Arguments.energy_map >0)
	{
		fstream lg_ener_inter("ener_lg.txt", ios::out);
		if(Arguments.mol_2_type.find("protein")!=string::npos)
			lg_ener_inter << "#Res" << setw(11) << "Electro." << setw(12) << "Van der W." << '\n';
		else
			lg_ener_inter << "#Atom/Res" << setw(11) << "Electro." << setw(12) << "Van der W." << '\n';
		fstream rec_ener_inter("ener_rec.txt", ios::out);
		rec_ener_inter << "#Res" << setw(12) << "Electro." << setw(12) << "Van der W." << '\n';
		lg_ener_inter.close();
		rec_ener_inter.close();
	}
	
	if (Arguments.force_map >0)
	{
		fstream lg_force_inter("force_lg.txt", ios::out);
		if (Arguments.mol_2_type.find("protein") != string::npos)
			lg_force_inter << "#Res" << setw(11) << "Electro." << setw(12) << "Van der W." << '\n';
		else
			lg_force_inter << "#Atom/Res" << setw(11) << "Electro." << setw(12) << "Van der W." << '\n';
		fstream rec_force_inter("force_rec.txt", ios::out);
		rec_force_inter << "#Res" << setw(12) << "Electro." << setw(12) << "Van der W." << '\n';
		lg_force_inter.close();
		rec_force_inter.close();
	}

	/*************************************/
	double d_standard{ Arguments.d_standard };
	int num_amine{};
	size_t num_threads{ Arguments.nt }, tot_snap_allfiles{};
	string  name_pro{}, name_lig{};
	para parameter;

	vector<double> E_v, E_e;
	vector<size_t> size_snap;
	vector<vector<double> > E_v_vector, E_e_vector, Volumes;
	vector<force> Force_e, Force_v;
	vector<vector<force> > Force_e_vector, Force_v_vector;
	vector<energy_res> Force_res_stat_mol1, Force_res_stat_mol2;
	vector<vector<energy_res>> Energy_res_stat_1, Energy_res_stat_2;
	vector<vector<energy_groups> >  EF_subgroups;
	vector<vector<map_pos>> average_contact_matrix;
	vector<mindistance_group> MIN_DIS;
	vector<vector<size_t> > Hbonds_vector, sidechain_contact_num;
	amino_acid amine;
	vector<double> time_of_snap;
	vector<string> input_line, cscd_out_name;
	all_cscd cscd_out_dis;
	info_ener_res infoener_res;
	bool topo_2_is_open{ Arguments.flag_topo_2 }, topo_is_open(Arguments.flag_topo_1), status{true};
	
	groups_cal Groups_from_index_file;

	info info_output{Arguments};
	vector<traj_hb> res_hb_rec, res_hb_lg;
	vector<traj_sc> res_sc_rec, res_sc_lg;

	fstream out_sc(info_output.getSCout(), ios::out);
	/************************************/
	fstream cal_out;
	cal_out.open("cal.txt", ios::out);
	/*************************************/
	fstream log;
	log.open("log.txt", ios::out);
	log.close();
	log.open("log.txt", ios::out | ios::app);
	/*************************************/
	/******************************************/
	
	if (Arguments.com_dis_energy_2mol != 4)
	{
		for (size_t i{}; i < paths.get_size(); ++i)
		{
			fstream output_tot_dis(paths.get_path(i) + "/" + info_output.getDisCom2Mol(), ios::out);
			fstream output_elec_dis(paths.get_path(i) + "/electrostatic_" + info_output.getDisCom2Mol(), ios::out);
			fstream output_vdw_dis(paths.get_path(i) + "/vanderWaals_" + info_output.getDisCom2Mol(), ios::out);

			output_elec_dis.close();
			output_tot_dis.close();
			output_vdw_dis.close();
		}
	}

	if (Arguments.flag_h_bond == 1)
	{
		fstream out_h_bonds("out_hydrogen_bonds.dat", ios::out); out_h_bonds.close();
		fstream out_h_bonds_rec("out_hydrogen_bonds_receptor.dat", ios::out); 
		out_h_bonds_rec << "# atom of receptor | atom of ligand\n";
		out_h_bonds_rec.close();
		fstream out_h_bonds_lig("out_hydrogen_bonds_ligand.dat", ios::out);
		out_h_bonds_lig << "# atom of ligand | atom of receptor\n";
		out_h_bonds_lig.close();
	}
	log << "Software Build: " << info_output.getBuild() << '\n';

	time_t timeRun{ time(nullptr) };
	log << "Current time: " << put_time(localtime(&timeRun), "%c") << '\n';

	log << "Force field: "<< force_f << "  | Number of threads: " << num_threads << '\n';
	log << "** Your options: " << Arguments.command << '\n';
	if (Arguments.flag_energy_same_force == true)
		log << "Note: Options for energy and force are same!\n";
	log << "Parsing c6 and c12 and aminoacids info..." << '\n';

	take_para(parameter, force_f, Arguments);
	take_amino_acids(amine, num_amine, force_f, Arguments);
	take_ff_para(amine, parameter);

	
	int input_name{};
	gen_input_list(Arguments, input_line, paths, log, input_name);
	
	while (true)
	{
		auto time_onetraj= std::chrono::steady_clock::now();

		if (!input_line.empty())
		{
			if (Arguments.flag_all_folders == false)
			{
				name_lig = "lg_" + input_line[input_name - Arguments.b_input] + Arguments.get_input_file_type();
				name_pro = "rec_" + input_line[input_name - Arguments.b_input] + Arguments.get_input_file_type();
			}
			else
			{
				
				name_lig = fs::path(input_line[input_name - Arguments.b_input]).parent_path().string()+"/lg_" + fs::path(input_line[input_name - Arguments.b_input]).stem().string() + Arguments.get_input_file_type();
				name_pro = fs::path(input_line[input_name - Arguments.b_input]).parent_path().string() + "/rec_" + fs::path(input_line[input_name - Arguments.b_input]).stem().string() + Arguments.get_input_file_type();
			}
		}
		else
		{
			name_lig = "lg_" + to_string(input_name) + Arguments.get_input_file_type();
			name_pro = "rec_" + to_string(input_name) + Arguments.get_input_file_type();
		}

		Multi_MOL pro_in(Arguments.force_f), lig(Arguments.force_f);
		h_groups groups_1, groups_2;
		vector<energy_res> Energy_res_mol1, Energy_res_mol2, Force_res_mol1, Force_res_mol2;
		//vector<vector<com> > COM_vec_rec, COM_vec_lg;
		vector<energy_groups> Energy_groups_out;
		vector<fstream> out_cal_dis_files;
		traj_hb hb_rec_1traj, hb_lg_1traj;
		traj_sc sc_rec_1traj, sc_lg_1traj;
		multi_snap_emap rec_ener_inter;
		vector<vector<map_pos>> cont_matrix;
		unordered_map<tuple<size_t, size_t, size_t, size_t>, pair<size_t, size_t>, tuple_hash> pos_matrix;
		auto start_pdb{ std::chrono::steady_clock::now() };
		auto start = std::chrono::steady_clock::now();

		if (parse_input_pdb(Arguments, log, name_pro, name_lig, force_f, input_name, input_line, topo_2_is_open
			, topo_is_open, pro_in, lig, Groups_from_index_file, parameter, amine, info_output) == 1)
		{
			if (Arguments.flag_nocal == true)
			{
				cout << "Nocal mode is turned on!\n"; return 1;
			}
			else
			{
				if (input_name > Arguments.b_input)--input_name;
				status = false;
				if (Arguments.flag_pattern_input == true && Arguments.flag_range_input == false)
					log << "Pattern input ends at " << Arguments.pdb_name << " file!\n";
				else
					log << "Error occured when parsing " << Arguments.get_input_file_type() << "(s) file!\n";
				break;
			}
		}

		info_output.get_info_group(Groups_from_index_file);

		set_chainnum(pro_in, lig, Arguments);
		/*************************************/
		if (Arguments.opdb == true && Arguments.flag_sep_opdb == false)
		{
			fstream out_pdb(Arguments.get_pdbOutput(), ios::out);
			out_pdb.close();
		}
		/*************************************/
		double time_cscd{}, time_new_emap{};
		std::chrono::milliseconds time_cal{}, time_emap{}, time_fmap{}, time_subgroups_cal{}, time_solve_pbc{}, 
			time_sc{}, time_hb{}, time_other{};
		vector<thread> threads(num_threads);

		if (Arguments.flag_nocal == true)
		{
			impl_nocal(pro_in, lig, log, Arguments, info_output, threads, amine, parameter, paths, input_line, groups_1, groups_2, 
				Hbonds_vector, input_name, Volumes);
		
			if (check_end_cal(Arguments, input_line, input_name) == true)
			{
				paths.out_volume(time_of_snap, Volumes, Arguments);
				return 0;
			}
			continue;
		}


		if (Arguments.com_dis_2mol != 0)
		{
			fstream output_dis_2mol_file(paths.get_path(input_name - Arguments.b_input) + "/" + info_output.getDisCom2Mol(), ios::out);
			output_dis_2mol_file.close();
		}

		if (Arguments.flag_cal_vol == true)
		{
			vector<double> volumes_temp;
			pro_in.cal_volumes(volumes_temp);
			Volumes.emplace_back(volumes_temp);
		}

		if (topo_is_open == true && topo_2_is_open == true)
		{
			log << "Taking charges from topology file(s)...\n";
			if (num_threads > 1)
			{
				log << "Taking charges of 2 molecules in parallel..." << '\n';
				for (size_t i{}; i < num_threads; ++i)
				{
					threads[i] = thread(&Multi_MOL::take_charge, &pro_in, ref(amine), ref(parameter), (i + 1) * pro_in.m_mol[0].Molecule.size() / num_threads, i * pro_in.m_mol[0].Molecule.size() / num_threads, force_f, 1);

				}
				for (auto& t : threads)t.join();
				for (size_t i{}; i < num_threads; ++i)
				{
					threads[i] = thread(&Multi_MOL::take_charge, &lig, ref(amine), ref(parameter), (i + 1) * lig.m_mol[0].Molecule.size() / num_threads, i * lig.m_mol[0].Molecule.size() / num_threads, force_f, 1);

				}
				for (auto& t : threads)t.join();
			}
			else
			{
				log << "Taking charges of 2 molecules...\n";
				pro_in.take_charge(amine, parameter, pro_in.m_mol[0].Molecule.size(), 0ull, force_f, 1);
				lig.take_charge(amine, parameter, lig.m_mol[0].Molecule.size(), 0ull, force_f, 1);
			}
		}
		else if (topo_is_open == true && Arguments.flag_topo_1 == true && topo_2_is_open == false)
		{
			log << "Taking charges of rec..." << '\n';
			//if (topo_2_is_open == false)
			//{
			if (num_threads > 1ull)
			{
				vector<thread> threads_charge(num_threads);
				for (size_t i{}; i < num_threads; ++i)
				{
					threads_charge[i] = thread(&Multi_MOL::take_charge, &pro_in, ref(amine), ref(parameter), (i + 1) * pro_in.m_mol[0].Molecule.size() / num_threads, i * pro_in.m_mol[0].Molecule.size() / num_threads, force_f, 0);
				}
				for (auto& t : threads_charge)t.join();
			}
			else
				pro_in.take_charge(amine, parameter, pro_in.m_mol[0].Molecule.size(), 0ull, force_f, 0);
			//}
			//else
			//	pro_in.take_charge(amine, parameter, pro_in.m_mol[0].Molecule.size(), 0ull, force_f, 1);
			log << "Taking charges of lig...\n";
			lig.take_charge(amine, parameter, lig.m_mol[0].Molecule.size(), 0ull, force_f, 1);
		}
		else if (topo_is_open == false && topo_2_is_open == true)
		{
			if (num_threads > 1)
			{
				for (size_t i = 0; i < num_threads; ++i)
				{
					threads[i] = thread(&Multi_MOL::take_charge, &lig, ref(amine), ref(parameter), (i + 1) *lig.m_mol[0].Molecule.size() / num_threads, i *lig.m_mol[0].Molecule.size() / num_threads, force_f, 0);

				}
				for (auto& t : threads)t.join();

				pro_in.take_charge(amine, parameter, pro_in.m_mol[0].Molecule.size(), 0ull, force_f, 1);
			}
			else
			{
				pro_in.take_charge(amine, parameter, pro_in.m_mol[0].Molecule.size(), 0ull, force_f, 1);
				lig.take_charge(amine, parameter, lig.m_mol[0].Molecule.size(), 0ull, force_f, 0);
			}
		}
		else
		{
			if (num_threads > 1)
			{
				log << "Taking charges of 2 molecules in parallel..." << '\n';
				for (size_t i{}; i < num_threads; ++i)
				{
					threads[i] = thread(&Multi_MOL::take_charge, &pro_in, ref(amine), ref(parameter), (i + 1) *pro_in.m_mol[0].Molecule.size() / num_threads, i *pro_in.m_mol[0].Molecule.size() / num_threads, force_f, 0);

				}
				for (auto& t : threads)t.join();
				for (size_t i{}; i < num_threads; ++i)
				{
					threads[i] = thread(&Multi_MOL::take_charge, &lig, ref(amine), ref(parameter), (i + 1) *lig.m_mol[0].Molecule.size() / num_threads, i *lig.m_mol[0].Molecule.size() / num_threads, force_f, 0);

				}
				for (auto& t : threads)t.join();
			}
			else
			{
				log << "Taking charges of 2 molecules...\n";
				pro_in.take_charge(amine, parameter, pro_in.m_mol[0].Molecule.size(), 0ull, force_f, 0);
				lig.take_charge(amine, parameter, lig.m_mol[0].Molecule.size(), 0ull, force_f, 0);
			}
			/***********************************************************************/
		}

		auto time_pdb = std::chrono::duration_cast<milliseconds>(std::chrono::steady_clock::now() - start_pdb);
		log << "Box type: " << pro_in.m_mol[0].Box.get_box_type() << '\n';

		cache_sigma_epsilon CachedData{Arguments.distance_energy_cutoff};
		if (Arguments.flag_cache == 1
			|| Arguments.force_map == 1 || Arguments.force_map == 3 || Arguments.flag_new_emap == true ||
			Arguments.energy_map == 1 || Arguments.energy_map == 3 || Arguments.energy_map == 4)
		{
			if (CachedData.set_size(pro_in.m_mol[0].Molecule.size()*lig.m_mol[0].Molecule.size(), Arguments, pro_in.m_mol[0].Molecule.size(), lig.m_mol[0].Molecule.size()) == 1)
			{
				cout << "Big system, turned off Cache!\n";
				log << "Big system, turned off Cache!\n";
				
			}
		}

		/********************************************************************/
		log << "Solve pbc problem \n";

		process_time_of_infiles(pro_in, Arguments, info_output, time_of_snap);

		auto start_pbc{ chrono::steady_clock::now() };
		if (Arguments.solv_pbc_int == 1)
		{
			if (Arguments.chain_ref == -2)
			{
				for (size_t i{}; i < num_threads; ++i)
				{
					threads[i] = thread(&Multi_MOL::solve_pbc_int, &pro_in, i, num_threads);
				}
				for (auto& t : threads)t.join();
			}
			if (Arguments.chain_ref > -2)
			{
				log << "Will use new pbc algorithm!\n";
				cout << "Will use new pbc algorithm!\n";
				if (Arguments.solv_pbc_int2 == 1)
					log << "pbc algorithm with minimum distance criteria\n";
				if(Arguments.solv_pbc_int3==1)
					log<< "pbc algorithm with minimum total distance criteria\n";
				if (Arguments.pbc_aggressive == true)
				{
					log << "Using new algorithm to solv mol2 PBC.\n";
					cout << "Using new algorithm to solv mol2 PBC.\n";
				}
				if (Arguments.flag_pbc_ref == true)
				{
					log << "solve pbc using reference file\n";
					cout << "solve pbc using reference file\n";
				}
				pro_in.solve_pbc_new(lig, Arguments, amine);
				align_complex(pro_in, lig, 0, 0);
			}
		}
		if (Arguments.solv_pbc == 1)
		{
			for (size_t i{}; i < num_threads; ++i)
			{
				threads[i] = thread(&Multi_MOL::pbc_prob, &pro_in, ref(lig), i, num_threads);
			}
			for (auto& t : threads)t.join();
			align_complex(pro_in, lig, 0, 0);
		}
		time_solve_pbc += chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - start_pbc);
		/***********************************************************************/

		if (Arguments.entropy_quasi == true)
		{
			log << "Calculating entropy by quasi harmonic method ...\n";
			cout << "Calculating entropy by quasi harmonic method ...\n";
			cout << "Aligning system......\n";
			if (Arguments.solv_pbc == 0 && Arguments.chain_ref == -1) align_complex(pro_in, lig, 0, 0);
			log << "Entropy of groups " << Groups_from_index_file.group_name(Groups_from_index_file.group_entropy_com_index()) << " -- "
				<< Groups_from_index_file.group_name(Groups_from_index_file.group_entropy_lg_index()) << " -- "
				<< Groups_from_index_file.group_name(Groups_from_index_file.group_entropy_rec_index()) << '\n';
			cal_entropy(pro_in, lig, Groups_from_index_file, Arguments, paths, static_cast<size_t>(input_name), info_output);
		}
		/***********************************************************************/

		if (pro_in.m_mol.size() != lig.m_mol.size())
		{
			cout << "WARNING! 2 molecules have a different number of snaps.\n";
			cout << "molecule 1 has " << pro_in.m_mol.size() << " snaps and molecule 2 has " << lig.m_mol.size() << " snaps \n";
			log << "*********************************************************\n";
			log << "WARNING! 2 molecules have a different number of snaps. ";
			log << "molecule has " << pro_in.m_mol.size() << " snaps and molecule 2 has " << lig.m_mol.size() << " snaps \n";
		}
		else
		{
			log << "********************** \n";
			log << "Number of snaps: " << pro_in.m_mol.size() << '\n';
			log << "********************** \n";
		}
		size_t NUM{ pro_in.m_mol.size() > lig.m_mol.size() ? lig.m_mol.size() : pro_in.m_mol.size() };
		tot_snap_allfiles += NUM;
		log << "Charge of " << Groups_from_index_file.group_name(Groups_from_index_file.group_1_index()) << ": " << pro_in.total_charge() << "e\n"
			"Charge of " << Groups_from_index_file.group_name(Groups_from_index_file.group_2_index()) << ": " << lig.total_charge() << "e\n"
			<< "**********************\n";
		log << "Calculating....\n";
		cout << "Calculating....\n";

		//pulling
		if (Arguments.flag_find_pulling_direct == true)
		{
			cout << "Finding of pulling direction mode is on!\n";
			log << "Finding of pulling direction mode is on!\n";
			find_pull_direct(pro_in, lig, CachedData, DE_arguments, Arguments, info_output, Groups_from_index_file);
			return 0;
		}

		for (size_t num_snap{}; num_snap < NUM; ++num_snap)
		{
			for (size_t p{ 0 }; p <= 100; p += 20)
				if (num_snap == NUM / 100 * p)cout << p << "% ";
			if (num_snap == NUM-1)cout << '\n';

			vector<d_min> dis_min(num_threads);
			vector<cache_energy> Energy(num_threads);
			vector<cache_force> Force(num_threads);
			
			double E_V_{}, E_E_{};
			force f_e_{}, f_v_{};

			log << setw(10) << pro_in.m_mol[num_snap].time_of_snap;
			if ((num_snap +1 ) % 15 == 0)
				log << '\n';

			auto start_cal= std::chrono::steady_clock::now();
			if (Arguments.force_f.find("CHARMM") != string::npos || Arguments.force_f.find("AMBER") != string::npos)
			{
				for (size_t i{}; i < num_threads; ++i)
				{
					if (Arguments.flag_elec_screening == false)
					{
						if (Arguments.flag_cal_force == false)
						{
							if (!CachedData.parameters.empty())
							{
								threads[i] = thread(cal_energy_amber, ref(pro_in.m_mol[num_snap]), ref(lig.m_mol[num_snap]), ref(pro_in.info), ref(lig.info), ref(CachedData), lig.m_mol[0].Molecule.size() * i / num_threads,
									lig.m_mol[0].Molecule.size() * (i + 1) / num_threads, ref(Energy[i]), num_snap, ref(Arguments), ref(pro_in.parameters), ref(lig.parameters));
							}
							else
								threads[i] = thread(cal_energy_amber_intrin, ref(pro_in.m_mol[num_snap]), ref(lig.m_mol[num_snap]), ref(pro_in.info), ref(lig.info), ref(CachedData), lig.m_mol[0].Molecule.size() * i / num_threads,
									lig.m_mol[0].Molecule.size() * (i + 1) / num_threads, ref(Energy[i]), num_snap, ref(Arguments), ref(pro_in.parameters), ref(lig.parameters));
						}
						else
							threads[i] = thread(cal_energy_force_amber, ref(pro_in.m_mol[num_snap]), ref(lig.m_mol[num_snap]), ref(pro_in.info), ref(lig.info), ref(CachedData), lig.m_mol[0].Molecule.size() * i / num_threads,
								lig.m_mol[0].Molecule.size() * (i + 1) / num_threads, ref(Energy[i]), ref(Force[i]), num_snap, ref(Arguments), ref(pro_in.parameters), ref(lig.parameters));
					}
					else
					{
						if (Arguments.flag_cal_force == false)
							threads[i] = thread(cal_energy_amber_screening, ref(pro_in.m_mol[num_snap]), ref(lig.m_mol[num_snap]), ref(pro_in.info), ref(lig.info), ref(CachedData), lig.m_mol[0].Molecule.size() * i / num_threads,
								lig.m_mol[0].Molecule.size() * (i + 1) / num_threads, ref(Energy[i]), num_snap, ref(Arguments), ref(pro_in.parameters), ref(lig.parameters));
						else
							threads[i] = thread(cal_energy_force_amber, ref(pro_in.m_mol[num_snap]), ref(lig.m_mol[num_snap]), ref(pro_in.info), ref(lig.info), ref(CachedData), lig.m_mol[0].Molecule.size() * i / num_threads,
								lig.m_mol[0].Molecule.size() * (i + 1) / num_threads, ref(Energy[i]), ref(Force[i]), num_snap, ref(Arguments), ref(pro_in.parameters), ref(lig.parameters));

					}
				}
				for (auto& t : threads)t.join();
			}
			else if (Arguments.force_f.find("OPLSAA") != string::npos)
			{
				for (size_t i{}; i < num_threads; ++i)
				{
					if (Arguments.flag_cal_force == false)
						threads[i] = thread(cal_energy_opls, ref(pro_in.m_mol[num_snap]), ref(lig.m_mol[num_snap]), ref(pro_in.info), ref(lig.info), ref(CachedData), lig.m_mol[0].Molecule.size()*i / num_threads,
							lig.m_mol[0].Molecule.size()*(i + 1) / num_threads, ref(Energy[i]), num_snap, ref(Arguments), ref(pro_in.parameters), ref(lig.parameters));
					else
						threads[i] = thread(cal_energy_force_opls, ref(pro_in.m_mol[num_snap]), ref(lig.m_mol[num_snap]), ref(pro_in.info), ref(lig.info), ref(CachedData), lig.m_mol[0].Molecule.size()*i / num_threads,
							lig.m_mol[0].Molecule.size()*(i + 1) / num_threads, ref(Energy[i]), ref(Force[i]), num_snap, ref(Arguments), ref(pro_in.parameters), ref(lig.parameters));
				}
				for (auto& t : threads)t.join();
			}
			else
			{
				for (size_t i{}; i < num_threads; ++i)
				{
					if (Arguments.flag_cal_force == false)
						threads[i] = thread(cal_energy_gromos, ref(pro_in.m_mol[num_snap]), ref(lig.m_mol[num_snap]), ref(pro_in.info), ref(lig.info), ref(CachedData), lig.m_mol[0].Molecule.size()*i / num_threads, 
							lig.m_mol[0].Molecule.size()*(i + 1) / num_threads, ref(Energy[i]), num_snap, ref(Arguments), ref(pro_in.parameters), ref(lig.parameters));
					else
						threads[i] = thread(cal_energy_force_gromos, ref(pro_in.m_mol[num_snap]), ref(lig.m_mol[num_snap]), ref(pro_in.info), ref(lig.info), ref(CachedData), lig.m_mol[0].Molecule.size()*i / num_threads, 
							lig.m_mol[0].Molecule.size()*(i + 1) / num_threads, ref(Energy[i]), ref(Force[i]), num_snap, ref(Arguments), ref(pro_in.parameters), ref(lig.parameters));
				}
				for (auto& t : threads)t.join();
			}
			time_cal += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_cal);
			//cout << "done cal snap " << num_snap << " ";
			if (Arguments.self_int == true)
				remove_double_counting(CachedData, Energy);

			for (size_t i{}; i < Energy.size(); ++i)
			{
				E_E_ += Energy[i].elec_ener;
				E_V_ += Energy[i].vdw_ener;
				if (Arguments.flag_cal_force == true)
				{
					f_e_ += Force[i].forces[0];
					f_v_ += Force[i].forces[1];
				}
			}
			
			process_adding_energies_forces(Arguments, input_name, E_v, E_e, E_E_, E_V_, f_e_, f_v_, E_v_vector, E_e_vector, Force_e, Force_v,
				Force_e_vector, Force_v_vector, size_snap, num_snap);
			//cout << " done check1 ";
			auto start_sc{ chrono::steady_clock::now()};

			size_t num_threads_rest(pro_in.Index.size() - 1 > 0 ? pro_in.Index.size() - 1 : 1);
			size_t ite_thread{};
			while (num_threads_rest > 0)
			{
				if (num_threads_rest <= num_threads)
				{
					for (size_t i{}; i < num_threads_rest; ++i)
					{
						threads[i] = thread(&MOL::take_com, &pro_in.m_mol[num_snap], ref(pro_in.parameters), i + ite_thread, pro_in.Index[i + ite_thread + 1].get_res_pos_e(0), ref(pro_in.Index), ref(pro_in.info));
					}
					for_each(threads.begin(), threads.begin() + num_threads_rest, mem_fn(&thread::join));
					num_threads_rest = 0ull;
				}
				else
				{
					for (size_t i{}; i < num_threads; ++i)
					{
						threads[i] = thread(&MOL::take_com, &pro_in.m_mol[num_snap], ref(pro_in.parameters), i + ite_thread, pro_in.Index[i + ite_thread + 1].get_res_pos_e(0), ref(pro_in.Index), ref(pro_in.info));
					}
					for (auto& t : threads)t.join();
					num_threads_rest -= num_threads;
					ite_thread += num_threads;
				}
			}


			num_threads_rest = lig.Index.size() - 1 > 0 ? lig.Index.size() - 1 : 1;
			ite_thread = 0;
			while (num_threads_rest > 0)
			{
				if (num_threads_rest <= num_threads)
				{
					for (size_t i{}; i < num_threads_rest; ++i)
					{
						threads[i] = thread(&MOL::take_com, &lig.m_mol[num_snap], ref(lig.parameters), i + ite_thread, lig.Index[i + ite_thread + 1].get_res_pos_e(0), ref(lig.Index), ref(lig.info));
					}
					for_each(threads.begin(), threads.begin() + num_threads_rest, mem_fn(&thread::join));
					num_threads_rest = 0ull;
				}
				else
				{
					for (size_t i{}; i < num_threads; ++i)
					{
						threads[i] = thread(&MOL::take_com, &lig.m_mol[num_snap], ref(lig.parameters), i + ite_thread, lig.Index[i + ite_thread + 1].get_res_pos_e(0), ref(lig.Index), ref(lig.info));
					}
					for (auto& t : threads)t.join();
					num_threads_rest -= num_threads;
					ite_thread += num_threads;
				}
			}

			impl_add_coms_2_mols(lig, pro_in, Arguments, num_snap);
			/*****************************************************************************************************************/

			for (size_t i{}; i < num_threads; ++i)
			{
				threads[i] = thread(take_min, ref(pro_in.m_mol[num_snap].Com), ref(lig.m_mol[num_snap].Com), ref(dis_min[i]), 
					i *pro_in.m_mol[num_snap].Com.size() / num_threads, (i + 1) *pro_in.m_mol[num_snap].Com.size() / num_threads, 0, 
					lig.m_mol[num_snap].Com.size(), ref(Arguments), sc_flag, d_standard);
			}
			for (auto& t : threads)t.join();

			time_sc += chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - start_sc);
			//cout << "\nsnap " << num_snap << " mol1 consists of " << pro_in.m_mol[num_snap].Com.size() << " mol2 consists of " << lig.m_mol[num_snap].Com.size() << " residues\n";
			//cout << "\nsnap: "<<num_snap<< ": numof dismin " << dis_min.size() << " ";
			//for (auto& d : dis_min)
			//	cout << d.base.size() << ' '; cout << "mem size: " << Arguments.nt * (sizeof(std::vector<dmin_base>) + (sizeof(dmin_base) * dis_min[0].base.size()));
			//cout << '\n';
			print_sidechain_contact(num_threads, dis_min, Arguments, out_sc, sidechain_contact_num, 
				sc_rec_1traj, sc_lg_1traj, input_name - Arguments.b_input);
			//cout << "size of sc rec 1 traj: " << sc_rec_1traj.residues_sc.size() << " " << sc_rec_1traj.map_residues_sc.size()
			//	<< " " << sizeof(sc_rec_1traj.residues_sc[0]) << " " << sizeof(sc_rec_1traj.map_residues_sc) << '\n';
			if (Arguments.cal_dis_sc_contact == true)
			{
				clock_t start_cscd{ clock() };
				impl_cscd(pro_in, lig, Arguments, log, threads, num_snap, sc_flag, input_name, cscd_out_name, cscd_out_dis, paths, out_cal_dis_files);
				time_cscd += clock() - start_cscd;
			}

			if (Arguments.flag_pqr == 1)
			{
				write_pqr(pro_in, lig, input_line, input_name, num_snap, Arguments, paths);
			}
			write_pdb(pro_in, lig, Arguments, num_snap, paths, input_line, input_name);

			auto start_emap= std::chrono::steady_clock::now();
			if (Arguments.flag_energy_same_force == false && (Arguments.energy_map == 1 || Arguments.energy_map == 3 || Arguments.energy_map == 4))
			{
				wrap_cal_energy_map(threads, pro_in, lig, Arguments, num_snap, dis_min, CachedData, Energy_res_mol1, Energy_res_mol2, info_output, infoener_res);
			}
			time_emap += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_emap);
			//cout << "check emap..";
			//cout << "Ener resmol1 size: " << Energy_res_mol1.size() << " size 1 ele: " << sizeof(Energy_res_mol1[0]) << '\n';
			//cout << "Ener resmol2 size: " << Energy_res_mol2.size() << " size 1 ele: " << sizeof(Energy_res_mol2[0]) << '\n';
			auto start_fmap{ std::chrono::steady_clock::now() };
			if (Arguments.flag_cal_force == true && (Arguments.force_map == 1 || Arguments.force_map == 3))
			{
				wrap_cal_force_map(threads, pro_in, lig, Arguments, num_snap, dis_min, CachedData, Force_res_mol1, Force_res_mol2, info_output);
			}
			time_fmap += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_fmap);

			auto start_hb{ std::chrono::steady_clock::now() };
			
			if (Arguments.flag_h_bond == 1)
			{
				find_h_bonds(pro_in.m_mol[num_snap], pro_in.m_mol[0], lig.m_mol[num_snap], lig.m_mol[0], groups_1, groups_2, Arguments, num_snap, pro_in.info, lig.info, Hbonds_vector);
				
				if (Arguments.nt >= 2)
				{
					thread t_hb_1(add_hb_res, ref(pro_in.m_mol[num_snap]), ref(lig.m_mol[num_snap]), 
						ref(pro_in.info), ref(lig.info), Arguments, ref(hb_rec_1traj));
					thread t_hb_2(add_hb_res, ref(lig.m_mol[num_snap]), ref(pro_in.m_mol[num_snap]),
						ref(lig.info), ref(pro_in.info), Arguments, ref(hb_lg_1traj));
					t_hb_1.join();
					t_hb_2.join();
				}
				else
				{
					add_hb_res(pro_in.m_mol[num_snap], lig.m_mol[num_snap], pro_in.info, lig.info, Arguments, hb_rec_1traj);
					add_hb_res(lig.m_mol[num_snap], pro_in.m_mol[num_snap], lig.info, pro_in.info, Arguments, hb_lg_1traj);
				}
			}
			time_hb += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_hb);
			//cout << "check hbond..";
			clock_t start_new_emap{ clock() };
			if (Arguments.flag_new_emap == true)
			{
				imp_cal_energy_rec_res(rec_ener_inter, pro_in, lig, num_snap, 0, 0, dis_min, CachedData, Arguments);
			}
			time_new_emap += clock() - start_new_emap;
			//cout << "check new emap..";
			auto start_subgroups = std::chrono::steady_clock::now();
			if (Arguments.flag_sub_groups == true)
			{
				Energy_force_subgroups(Groups_from_index_file, pro_in, lig, CachedData, EF_subgroups, Energy_groups_out, num_snap, NUM, Arguments, threads);
				if (Arguments.flag_cal_mindis == true)
					mindis_subgroups(Groups_from_index_file, pro_in, lig, num_snap, NUM, Arguments, threads, MIN_DIS);
			}
			time_subgroups_cal += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_subgroups);
			//cout << "check sub groups..";
			auto start_dis2commol{ std::chrono::steady_clock::now() };
			if(Arguments.com_dis_energy_2mol != 4)
				dis_com_2mol(pro_in, lig, num_snap, NUM, Arguments, input_name, paths, E_E_, E_V_, info_output);
			if (Arguments.com_dis_2mol != 0)
			{
				fstream output_dis_2mol_file(paths.get_path(input_name - Arguments.b_input) + "/" + info_output.getDisCom2Mol(), ios::out|ios::app);
				pro_in.impl_com_dis_2mol(lig, Arguments, info_output, output_dis_2mol_file, num_snap);
			}
			//cout << "begin check contact ";
			if (Arguments.flag_cal_contact == true)
			{
				calc_contact_map(pro_in, lig, num_threads, dis_min, Arguments, num_snap, cont_matrix, pos_matrix);
				
			}
			//cout << "checked contact..";
			time_other += std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_dis2commol);

			if (num_snap != 0 || (NUM == 1 && Arguments.energy_map != 4)) {
				pro_in.m_mol[num_snap].clear_MOL();
				lig.m_mol[num_snap].clear_MOL();
			}
			
		}

		auto start_other{ std::chrono::steady_clock::now() };

		if (!out_cal_dis_files.empty()) {
			for (auto& file : out_cal_dis_files) if (file.is_open()) file.close();
		}
		if (Arguments.print_energy_res == true)
		{
			if (Arguments.flag_energy_same_force == true)
			{
				print_energy_res(Force_res_mol1, info_output.get_group_1_name(), paths, input_name - Arguments.b_input, Arguments);
				print_energy_res(Force_res_mol2, info_output.get_group_2_name(), paths, input_name - Arguments.b_input, Arguments);
			}
			else
			{
				print_energy_res(Energy_res_mol1, info_output.get_group_1_name(), paths, input_name - Arguments.b_input, Arguments);
				print_energy_res(Energy_res_mol2, info_output.get_group_2_name(), paths, input_name - Arguments.b_input, Arguments);
			}

			Energy_res_stat_1.emplace_back(Energy_res_mol1);
			Energy_res_stat_2.emplace_back(Energy_res_mol2);
		}

		if (Arguments.flag_cal_force == true && Arguments.print_force_res == true)
		{
			print_force_res(Force_res_mol1, info_output.get_group_1_name(), paths, input_name - Arguments.b_input, Arguments);
			print_force_res(Force_res_mol2, info_output.get_group_2_name(), paths, input_name - Arguments.b_input, Arguments);
		}

		if (Arguments.flag_maxforce == true)
		{
			find_maxforce_res(Force_res_mol1, info_output.get_group_1_name(), paths, input_name - Arguments.b_input, Arguments, Force_res_stat_mol1);
			find_maxforce_res(Force_res_mol2, info_output.get_group_2_name(), paths, input_name - Arguments.b_input, Arguments, Force_res_stat_mol2);
		}

		if (Arguments.flag_new_emap == true)
		{
			rec_ener_inter.ana_emap(Arguments, info_output, paths, input_name);
		}

		if (Arguments.flag_h_bond == 1)
		{
			hb_lg_1traj.snap_num = NUM; hb_rec_1traj.snap_num = NUM;
			res_hb_lg.emplace_back(hb_lg_1traj);
			res_hb_rec.emplace_back(hb_rec_1traj);
		}

		if (Arguments.flag_cal_contact == true)
		{
			sc_lg_1traj.snap_num = NUM; sc_rec_1traj.snap_num = NUM;

			sc_rec_1traj.clear_map(); sc_lg_1traj.clear_map();

			res_sc_lg.emplace_back(sc_lg_1traj);
			res_sc_rec.emplace_back(sc_rec_1traj);

			string cont_map_filename{ info_output.getContmap() + to_string(input_name) + ".dat" };
			fstream cont_map_file(cont_map_filename, ios::out);
			
			if (average_contact_matrix.empty())
			{
				//cout <<"cont_matrix size: " << cont_matrix.size() << '\n';
				average_contact_matrix = vector < vector<map_pos>>(cont_matrix.size());
				for (size_t a_i{}; a_i < average_contact_matrix.size(); ++a_i)
				{
					average_contact_matrix[a_i] = vector<map_pos>(cont_matrix[a_i].size(), map_pos());
					for (size_t a_j{}; a_j < average_contact_matrix[a_i].size(); ++a_j)
						average_contact_matrix[a_i][a_j] = cont_matrix[a_i][a_j];
				}
			}
			else
			{
				for (size_t a_i{}; a_i < average_contact_matrix.size(); ++a_i)
				{
					for (size_t a_j{}; a_j < average_contact_matrix[a_i].size(); ++a_j)
					{
						if(a_i<cont_matrix.size() && a_j<cont_matrix[a_i].size())
							average_contact_matrix[a_i][a_j].cont_num += cont_matrix[a_i][a_j].cont_num;
					}
				}
			}
			
			print_cont_map(cont_matrix, cont_map_file, NUM);
		}
		auto end_other{ std::chrono::steady_clock::now() };

		/**************************************************************************/
		log << "\n------------o0Time stat for one trajectory0o------------\n";
		log << "\nTotal calculation time: " << std::chrono::duration_cast<milliseconds>(std::chrono::steady_clock::now() - start).count() / 1000. << " s\n";
		log << "Time of parsing pdb= " <<  time_pdb.count() / 1000. << " s\n";
		log << "Time of solving pbc problem= " << time_solve_pbc.count() / 1000. << " s\n";
		log << "Time of energy calculation= " << time_cal.count()/1000. << " s\n";
		log << "Time of sidechain contact calculation= " << time_sc.count() / 1000. << " s\n"
			<< "Time of distance between side chains-side chains calculation= " << time_cscd / CLOCKS_PER_SEC << " s\n";
		log << "Time of HB calculation= " << time_hb.count() / 1000. << " s\n"
			<< "Time of emap calculation= " << time_emap.count()/1000. << " s\n"
			<< "Time of new emap calculation= " << time_new_emap / CLOCKS_PER_SEC << " s\n"
			<< "Time of fmap calculation= " << time_fmap.count() / 1000. << " s\n"
			<< "Time of contact map and dis2mol calculation= " << time_other.count() / 1000. << " s\n"
			<< "Time of subgroups calculation= " << time_subgroups_cal.count()/1000. << " s\n"
			<< "Other= " << std::chrono::duration_cast<milliseconds>(end_other - start_other).count() / 1000. << " s\n"
			<< "Time of one trajectory= " << std::chrono::duration_cast<milliseconds>(std::chrono::steady_clock::now() - time_onetraj).count() / 1000. << " s\n"
			<< "======================================================\n";

		if (check_end_cal(Arguments, input_line, input_name) == true)
		{
			if (Arguments.energy_map == 4)
			{
				auto start_enermaptype{ chrono::steady_clock::now()};
				log << "\nBuilding energy maptype...\n";
				infoener_res.wrapper_impl_ener_maptype(pro_in, lig, info_output, Arguments, tot_snap_allfiles);
				log << "Energy maptype done.\n";
				time_energymaptype += chrono::duration_cast<milliseconds>( chrono::steady_clock::now() - start_enermaptype);
			}

			break;
		}
	}
	out_sc.close();
	if (status == false && input_name == Arguments.b_input) { cout << "Error occured!\n"; return 1; }
	log<< "Analysing results .............\n";
	auto time_ana_other{ chrono::steady_clock::now() };
	/*if (Arguments.nt > 1)
	{
		thread ana1(out, ref(E_e), ref(E_v), ref(time_of_snap), ref(cal_out), ref(Arguments), input_name, ref(E_e_vector), ref(E_v_vector));
		thread ana2(&folders::out_to_folders, &paths, ref(E_e), ref(E_v), ref(time_of_snap), ref(E_e_vector), ref(E_v_vector), ref(Arguments));
		ana1.join();
		ana2.join();
	}
	else
	{*/
		out(E_e, E_v, time_of_snap, cal_out, Arguments, input_name, E_e_vector, E_v_vector, size_snap);
		paths.out_to_folders(E_e, E_v, time_of_snap, E_e_vector, E_v_vector, Arguments);

	//}

	if (Arguments.flag_cal_force == true)
	{
		out_force(Force_e, Force_v, time_of_snap, Arguments, input_name, Force_e_vector, Force_v_vector, size_snap);
		paths.out_force_to_folders(Force_e, Force_v, time_of_snap, Force_e_vector, Force_v_vector,Arguments);
	}

	if (Arguments.flag_maxforce == true)
	{
		stat_maxforce_res(Force_res_stat_mol1, info_output.get_group_1_name(), info_output, Arguments);
		stat_maxforce_res(Force_res_stat_mol2, info_output.get_group_2_name(), info_output, Arguments);
	}

	cal_sidechain_contact_time(Arguments, time_of_snap, sidechain_contact_num, info_output, paths);

	impl_entropy_duan_zhang(paths, E_e_vector, E_v_vector, time_of_snap, Arguments, E_e, E_v);

	print_cscd(Arguments, cscd_out_dis, cscd_out_name);

	paths.out_volume(time_of_snap, Volumes, Arguments);

	if (Arguments.flag_sub_groups == true)
	{
		post_process_Energy_force_subgroups(time_of_snap, EF_subgroups, paths, Arguments, info_output);
		if (Arguments.flag_cal_mindis == true)
			post_proc_mindis(Groups_from_index_file, Arguments, info_output, MIN_DIS);
	}

	stat_hbonds(Hbonds_vector, info_output, Arguments, time_of_snap,input_line);
	
	if (Arguments.flag_h_bond == 1)
	{
		if (num_threads > 1)
		{
			thread s_hbres_1(stat_hb_res, ref(res_hb_rec), info_output.get_res_list_hb_rec_filename());
			thread s_hbres_2(stat_hb_res, ref(res_hb_lg), info_output.get_res_list_hb_lg_filename());
			s_hbres_1.join();
			s_hbres_2.join();
		}
		else
		{
			stat_hb_res(res_hb_rec, info_output.get_res_list_hb_rec_filename());
			stat_hb_res(res_hb_lg, info_output.get_res_list_hb_lg_filename());
		}
	}

	if (Arguments.flag_cal_contact == true)
	{
		stat_sidechain_contact(res_sc_rec, info_output.get_res_list_sc_contact_rec_filename());
		stat_sidechain_contact(res_sc_lg, info_output.get_res_list_sc_contact_lg_filename());
	}

	if (Arguments.print_energy_res == true)
	{
		if (num_threads > 1)
		{
			thread s_ener_res_1(stat_energy_res, ref(Energy_res_stat_1), 1, ref(Arguments), ref(info_output), info_output.get_group_1_name());
			thread s_ener_res_2(stat_energy_res, ref(Energy_res_stat_2), 1, ref(Arguments), ref(info_output), info_output.get_group_2_name());
			s_ener_res_1.join();
			s_ener_res_2.join();
		}
		else
		{
			stat_energy_res(Energy_res_stat_1, 1, Arguments, info_output, info_output.get_group_1_name());
			stat_energy_res(Energy_res_stat_2, 1, Arguments, info_output, info_output.get_group_2_name());
		}
	}

	if (Arguments.flag_cal_contact == true)
	{
		fstream ave_cont_file(info_output.getOuputAverageContactMap(), ios::out);
		//cout << "total snaps: " << tot_snap_allfiles << '\n';
		print_cont_map(average_contact_matrix, ave_cont_file, tot_snap_allfiles);
	}


	cout << "All calculations are done!\nThank you for using this tool!:)\n";
	log << "Thank you for using this tool!:)\n";
	log << "Email: nguyenhoanglinh9@duytan.edu.vn";
	log << "\n======================================================\n"
		<< "Time of energy maptype calculation= " << time_energymaptype.count() / 1000. << "s\n"
		<< "Time of analyze other results= " << chrono::duration_cast<milliseconds>(chrono::steady_clock::now() - time_ana_other).count() / 1000. << "s\n";

	auto totalruntime{ chrono::duration_cast<milliseconds>(chrono::steady_clock::now() - total_run_time).count() / 1000.};
	log	<< "Total run time= " << totalruntime << " s = " << totalruntime / 60 << " minutes\n";
	log << "DONE!\n=====================================================================\n";
	log.close();
	return 0;
}
